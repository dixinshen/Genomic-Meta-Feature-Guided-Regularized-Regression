opt_alpha = function(theta, mean_lap, z, c, alpha_init, method = "NR", rate = 0.01, thresh = 1e-7) {
    if (method == "GD") {
        # initialize alpha
        alpha = alpha_init
        lambda = exp(z %*% alpha)
        v = (lambda^2*c^2 + 2*lambda*(1-c)) / 2
        diff = 1e5
        iter = 0
        rate = rate
        # gradient descent loop
        while (iter <= 1000 & diff > thresh) {
            alpha_old = alpha
            # gradient
            g = crossprod(z, (-1/v+mean_lap^2+theta)*((1-c)*lambda+(c*lambda)^2))
            # gradient step
            alpha = alpha_old - rate*g
            lambda = exp(z %*% alpha)
            v = (lambda^2*c^2 + 2*lambda*(1-c)) / 2
            diff = max(abs(alpha-alpha_old))
            iter = iter + 1
        }
        return(list(alpha=as.vector(alpha), iter=iter, rate=rate))
    }

    if (method == "NR") {
        alpha = alpha_init
        diff = 1e5
        iter = 0
        # Newton-Raphson loop
        while (iter <= 100 & diff > thresh) {
            alpha_old = alpha
            lambda = exp(z %*% alpha_old)
            v = (lambda^2*c^2 + 2*lambda*(1-c)) / 2
            # Gradient
            g = crossprod(z, (-1/v+mean_lap^2+theta)*((1-c)*lambda+(c*lambda)^2))
            # Hessian
            H = crossprod(z, as.vector((lambda/v)^2*(1-c+c^2*lambda)^2 -
                                           (-1/v+mean_lap^2+theta)*(1-c+2*c^2*lambda)) * z)
            # NR step
            alpha = alpha_old - solve(H, g)
            diff = max(abs(alpha-alpha_old))
            iter = iter + 1
        }
        return(as.vector(alpha))
    }

    if (method == "BFGS") {
        likelihood = function(alpha, z, c, theta, mean_lap) {
            lambda = exp(z %*% alpha)
            v = (lambda^2*c^2 + 2*lambda*(1-c)) / 2
            return(as.numeric(-sum(log(v)) + crossprod(mean_lap^2+theta, v)))
        }
        g = function(alpha, z, c, theta, mean_lap) {
            lambda = exp(z %*% alpha)
            v = (lambda^2*c^2 + 2*lambda*(1-c)) / 2
            return(crossprod(z, (-1/v+mean_lap^2+theta)*((1-c)*lambda+(c*lambda)^2)))
        }

        alpha = optim(alpha_init, likelihood, g, z=z, c=c, theta=theta, mean_lap=mean_lap, method = "BFGS")$par
        return(alpha)
    }
}

quad_approx = function(x, beta, Ck, Ri, d, n, m) {
    # quadratic approximation of Cox partial likelihood
    eta = drop(x %*% beta)
    exp_eta = exp(eta)
    sum_exp_eta_prime = 0
    sum_exp_eta = numeric()
    for (i in m:1) {
        sum_exp_eta[i] = sum_exp_eta_prime + sum( exp_eta[(n-Ri[i]+1):(n-Ri[i+1])] )
        sum_exp_eta_prime = sum_exp_eta[i]
    }
    u_prime = 0
    u2_prime = 0
    u = numeric()
    u2 = numeric()
    for (k in 1:n) {
        if (Ck[k+1] - Ck[k] == 0) {
            u[k] = u_prime
            u2[k] = u2_prime
        } else {
            u[k] = u_prime + d[Ck[k+1]] / sum_exp_eta[Ck[k+1]]
            u2[k] = u2_prime + d[Ck[k+1]] / (sum_exp_eta[Ck[k+1]])^2
        }
        u_prime = u[k]
        u2_prime = u2[k]
    }
    B = exp_eta*u - exp_eta*exp_eta*u2   # weights
    return(B)
}

initialize_cox = function(y) {
    ## initialize survival times
    ## input: survival outcome (in non-decreasing order times)
    ## output: Ri(risk set), Ck(valid time set), m(# of unique event time), d(# of events at each event time)

    n = NROW(y)
    t_event <- y[,1][y[,2]==1]   # event time
    D <- unique(t_event)   # unique event time
    m <- length(D)   # number of unique times
    d <- numeric()
    for(i in 1:m) {
        d[i] <- sum(t_event==D[i])
    }

    Ck_prime <- 0
    Ck <- numeric()
    Ri <- numeric()
    for (k in 1:n) {
        Ck[k] <- Ck_prime
        for (j in (Ck_prime+1):m) {
            if (D[j] <= y[k,1]) {
                Ck[k] <- Ck[k] + 1
                Ri[Ck[k]] <- n - k + 1
            } else {
                break
            }
            Ck_prime <- Ck[k]
        }
        if (Ck_prime==m) {break}
    }
    if (k < n) {Ck[(k+1):n] <- Ck_prime}
    Ri <- c(Ri, 0)
    Ck <- c(0, Ck)

    return(list(Ri=Ri, Ck=Ck, m=m, d=d))
}

glmtune = function(x, y, z, c = 0.5, family = "cox", alpha_init = rep(0,NCOL(z)), standardize = TRUE, intercept = TRUE,
                   z_intercept = TRUE, method = "BFGS", itermax_lap = 500, itermax_mm = 100,
                   thresh_lap = 1e-4, thresh_mm = 1e-3, woodbury = NCOL(x)>NROW(x)) {
    ## xtune logistic&cox module main function
    ## input: data x, y, z. c controls ridge, elastic net, lasso penalty type
    ## output: estimated beta, alpha

    if(z_intercept == TRUE) {z = cbind(1, z)}
    n = NROW(y)
    p = NCOL(x)
    q = NCOL(z)

    # if (standardize == TRUE) {
    #     xm = colMeans(x)
    #     x = sweep(x, 2, xm, "-")
    #     xs = drop(sqrt(crossprod(rep(1/n,n), x^2)))
    #     x = sweep(x, 2, xs, "/")
    # }

    if(family=="cox") {
        ordered = order(y[,1])
        y = y[ordered, ]
        x = as.matrix(x[ordered, ])
        likelihood_vals = initialize_cox(y)
    }

    # Lapalace approximation loop
    alpha = alpha_init
    iter_lap = 0
    dev_lap = 1e5
    while(iter_lap < itermax_lap & dev_lap > thresh_lap) {
        alpha_prev = alpha
        lambda = exp(z %*% alpha_prev)
        v = as.vector((lambda^2*c^2 + 2*lambda*(1-c)) / 2)
        # mean of Laplace approximation
        if(family=="binomial") {
            mean_lap = coef(glmnet(x, y, alpha=0, family="binomial", penalty.factor=v, standardize=F, intercept=intercept,
                            lambda=sum(v)/p/n))[-1]
            sigmoid = 1/(1+exp(-x%*%mean_lap))
            B = as.vector(sigmoid*(1-sigmoid))
        } else if(family=="cox") {
            mean_lap = tryCatch(drop(coef(glmnet(x, y, alpha=0, family="cox", penalty.factor=v, standardize=F,
                                                      lambda=sum(v)/p/n)))
                                , warning = function(e) {
                                    drop(coef(glmnet(x, y, alpha=0, family="cox", penalty.factor=v, standardize=F),s=sum(v)/p/n))}
                                )

            B = quad_approx(x, mean_lap, likelihood_vals$Ck, likelihood_vals$Ri,
                            likelihood_vals$d, n, likelihood_vals$m)
        }

        # Marjorization-minimization loop
        alpha_mm = alpha_prev
        iter_mm = 0
        dev_mm = 1e5
        while(iter_mm < itermax_mm & dev_mm > thresh_mm) {
            alpha_prev_mm = alpha_mm
            # majorize: calculate theta = inverse of var_lap
            if(!woodbury) {
                # var-cov matrix of Laplace approximation
                var_lap = crossprod(x, B*x) + diag(v)
                theta = diag(solve(var_lap))
            } else if(woodbury) {
                # use Woodbury matrix identity
                xvinv = sweep(x,2,v,"/")
                sigma_y = diag(ifelse(B==0, 0, 1/(B+1e-5))) + tcrossprod(x, xvinv)
                theta = 1/v - colSums(xvinv*solve(sigma_y,xvinv))
            }

            # minimize: solve argmin alpha
            alpha_mm = opt_alpha(theta, mean_lap, z, c, alpha_prev_mm, method=method)
            iter_mm = iter_mm + 1
            dev_mm = max(abs(alpha_mm - alpha_prev_mm))

            lambda = exp(z %*% alpha_mm)
            v = as.vector((lambda^2*c^2 + 2*lambda*(1-c)) / 2)
        }
        alpha = alpha_mm
        iter_lap = iter_lap + 1
        dev_lap = max(abs(alpha - alpha_prev))

        if (iter_lap%in%c(seq(10,50,10),100,500) | iter_lap%%1000 == 0) {
            cat("Laplace iteration ", iter_lap, ": difference of alpha in Laplace loop is ",
                dev_lap, "\n", sep="")
        }
    }
    cat("====================================================================================\n",
        "Final Lapalace approximation iteration is ", iter_lap, ": difference of alpha is ",
        dev_lap, "\n", sep="")

    # # unstandardize
    # if(standardize==TRUE) {
    #     x = sweep(x, 2, xs, "*")
    #     x = sweep(x, 2, xm, "+")
    # }

    # with estimated alpha, solve for beta
    pen_vec = exp(z%*%alpha)
    if(family=="binomial") {
        fit = glmnet(x, y, family="binomial", alpha=c, penalty.factor=pen_vec, standardize=standardize, intercept=intercept,
                     lambda=sum(pen_vec)/p/n)
        beta = coef(fit)
    } else if (family=="cox") {
        fit = glmnet(x, y, family="cox", alpha=c, penalty.factor=pen_vec, standardize=standardize, lambda=sum(pen_vec)/p/n)
        beta = coef(fit)

    }

    return(list(beta=beta, alpha=alpha, penalty.factor=pen_vec))
}
