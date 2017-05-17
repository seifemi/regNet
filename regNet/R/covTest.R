
#
#Modification of function covTest from package covTest. Rounding of results was removed. Internal function.
#
covTest <- function (fitobj, x, y, sigma.est = "full", status = NULL, maxp = min(nrow(x), ncol(x))) 
{
    s4 = substring(fitobj$call, 1, 4)[1]
    s7 = substring(fitobj$call, 1, 7)[1]
    s8 = substring(fitobj$call, 1, 8)[1]
    if (s4 == "lars" & s7 != "lars.en" & s8 != "lars.glm") {
        calltype = "lars"
        family = "gaussian"
    }
    if (s7 == "lars.en") {
        calltype = "lars.en"
        family = "gaussian"
    }
    if (s8 == "lars.glm") {
        calltype = "lars.glm"
        family = fitobj$family
    }
    if (family == "cox") 
        stop("Cox model not yet implemented")
    if (calltype == "lars") {
        if (fitobj$type != "LASSO") {
            stop("Call to Lars must use type='LASSO'")
        }
    }
    if (calltype == "lars") {
        type = "lar"
    }
    if (calltype == "lars.en") {
        type = "lars.en"
    }
    if (calltype == "lars.glm") {
        type = "lars.glm"
    }
    if (calltype == "lars.glm" & sigma.est == "full") {
        sigma.est = 1
        cat("glm model; sigma set to 1", fill = TRUE)
    }
    n = nrow(x)
    p = ncol(x)
    my = mean(y)
    lambda.min.ratio = ifelse(nrow(x) < ncol(x), 0.1, 1e-04)
    jlist = unlist(fitobj$act)
    if (type == "lar") 
        lamlist = c(fitobj$lambda, 0)
    if (type == "lars.en") 
        lamlist = c(fitobj$lambda, 0)
    if (type == "lars.glm") 
        lamlist = c(fitobj$lambda, 0)
    maxp.call = maxp
    maxp = length(jlist)
    maxp = min(maxp, which(lamlist == 0))
    maxp = min(maxp, maxp.call)
    jlist = jlist[1:maxp]
    cov0 = cov = sig = rep(NA, maxp)
    yy = y - my
    if (family == "binomial") {
        glmobj = glmnet(x, y, family = "binomial", standardize = fitobj$standardize, 
            lambda.min.ratio = lambda.min.ratio)
    }
    if (family == "cox") {
        glmobj = glmnet(x, Surv(y, status), family = "cox", standardize = fitobj$standardize, 
            lambda.min.ratio = lambda.min.ratio)
    }
    if (family == "cox") {
        junk = calcz02(x, y, status)
        sc = junk$sceta
        inf = junk$infeta
        miinf = misqrt((t(inf) + inf)/2)
        yy = t(sc) %*% miinf
    }
    for (j in 1:maxp) {
        if (jlist[j] > 0) {
            lambda = lamlist[j + 1]
            if (type == "lar") 
                yhat = predict(fitobj, x, s = lambda, type = "fit", 
                  mode = "lam")$fit
            if (type == "lars.en") 
                yhat = (1 + fitobj$lambda2) * predict.lars.en(fitobj, 
                  x, lambda)
            if (type == "lars.glm" & family == "binomial") {
                yhat = as.vector(predict(glmobj, x, type = "link", 
                  s = lambda/n))
            }
            if (type == "lars.glm" & family == "cox") {
                yhat = as.vector(predict(glmobj, x, type = "link", 
                  s = lambda/n))
            }
            cov[j] = sum(yy * yhat)
            if (j == 1) {
                cov0[j] = 0
            }
            if (j > 1) {
                tt0 = which(fitobj$beta[j, ] != 0)
                if (type == "lar") {
                  aa = update(fitobj, x = x[, tt0, drop = F])
                  yhat0 = predict(aa, x[, tt0], type = "fit", 
                    s = lambda, mode = "lam")$fit
                }
                if (type == "lars.en") {
                  aa = update(fitobj, x = x[, tt0, drop = F])
                  yhat0 = (1 + fitobj$lambda2) * predict.lars.en(aa, 
                    x[, tt0], lambda)
                }
                if (type == "lars.glm") {
                  if (family == "binomial") {
                    if (length(tt0) == 1) {
                      tt0 = c(tt0, tt0)
                    }
                    glmobj0 = glmnet(x[, tt0, drop = F], y, family = "binomial", 
                      standardize = fitobj$standardize, lambda.min.ratio = lambda.min.ratio)
                    yhat0 = as.vector(predict(glmobj0, x[, tt0, 
                      drop = F], type = "link", s = lambda/n))
                  }
                  if (family == "cox") {
                    if (length(tt0) == 1) {
                      tt0 = c(tt0, tt0)
                    }
                    glmobj0 = glmnet(x[, tt0, drop = F], Surv(y, 
                      status), family = "cox", standardize = fitobj$standardize, 
                      lambda.min.ratio = lambda.min.ratio)
                    yhat0 = as.vector(predict(glmobj0, x[, tt0, 
                      drop = F], type = "link", s = lambda/n))
                  }
                }
                cov0[j] = sum(yy * yhat0)
            }
        }
    }
    if (is.numeric((sigma.est))) {
        sigma = sigma.est
        sigma.type = "known"
        null.dist = "Exp(1)"
        if (sigma.est <= 0) {
            stop("sigma.est must be positive")
        }
    }
    if (sigma.est == "full") {
        if (nrow(x) < ncol(x) + 1) 
            stop("Number of observations must exceed number of variables,\nwhen sigma.est is `full; you need to specify a numeric value for sigma.est")
        sigma.type = "unknown"
        aaa = lsfit(x, y)
        sigma = sqrt(sum(aaa$res^2)/(n - p))
        np = n - p
        null.dist = paste("F(2,", as.character(np), ")", sep = "")
    }
    tt = ((cov - cov0)/sigma^2)
    if (sigma.type == "known") {
        out = cbind(jlist, tt, 1 - pexp(tt, 1))
        dimnames(out)[[2]] = c("Predictor_Number", "Drop_in_covariance", 
            "P-value")
    }
    if (sigma.type == "unknown") {
        out = cbind(jlist, tt, 1 - pf(tt, 2, n - p))
        dimnames(out)[[2]] = c("Predictor_Number", "Drop_in_covariance", 
            "P-value")
    }
    dimnames(out)[[1]] = rep("", nrow(out))
    #return(list(results = round(out, 4), sigma = round(sigma, 
    #    4), null.dist = null.dist))
    
    ##Return results without rounding
    return(list(results = out, sigma = sigma, null.dist = null.dist))
}
 
