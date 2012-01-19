# Time series diagnostics
# Source: Robert H. Shumway and David S. Stoffer, 2011
# URL: http://www.stat.pitt.edu/stoffer/tsa3/

# ACF and PACF
acf2 <- function(series,max.lag=NULL){
  num=length(series)
  if (num > 49 & is.null(max.lag)) max.lag=ceiling(10+sqrt(num))
  if (num < 50 & is.null(max.lag))  max.lag=floor(5*log10(num))
  if (max.lag > (num-1)) stop("Number of lags exceeds number of observations")
  ACF=acf(series, max.lag, plot=FALSE)$acf[-1]
  PACF=pacf(series, max.lag, plot=FALSE)$acf
  LAG=1:max.lag/frequency(series)
  minA=min(ACF)
  minP=min(PACF)
  U=2/sqrt(num)
  L=-U
  minu=min(minA,minP,L)-.01
  old.par <- par(no.readonly = TRUE)
  par(mfrow=c(2,1), mar = c(3,3,2,0.8),
    oma = c(1,1.2,1,1), mgp = c(1.5,0.6,0))
  plot(LAG, ACF, type="h",ylim=c(minu,1), 
    main=paste("Series: ",deparse(substitute(series))))
    abline(h=c(0,L,U), lty=c(1,2,2), col=c(1,4,4))
  plot(LAG, PACF, type="h",ylim=c(minu,1))
    abline(h=c(0,L,U), lty=c(1,2,2), col=c(1,4,4))
  on.exit(par(old.par))  
  ACF<-round(ACF,2); PACF<-round(PACF,2)    
  return(cbind(ACF, PACF)) 
  }

# Lag plots
lag.plot1 <- function(series,max.lag=1,corr=TRUE,smooth=TRUE){ 
   name1=paste(deparse(substitute(series)),"(t-",sep="")
   name2=paste(deparse(substitute(series)),"(t)",sep="")
   data1=as.ts(series)
   max.lag=as.integer(max.lag)
   prow=ceiling(sqrt(max.lag))
   pcol=ceiling(max.lag/prow)
   a=acf(series,max.lag,plot=FALSE)$acf[-1]
   old.par <- par(no.readonly = TRUE)
   par(mfrow=c(prow,pcol), mar=c(2.5, 4, 2.5, 1), cex.main=1.1, font.main=1)
  for(h in 1:max.lag){                       
   plot(lag(series,-h), data1, xy.labels=FALSE, main=paste(name1,h,")",sep=""), ylab=name2, xlab="") 
    if (smooth==TRUE) 
    lines(lowess(ts.intersect(lag(series,-h),series)[,1],
                 ts.intersect(lag(series,-h),series)[,2]), col="red")
    if (corr==TRUE)
    legend("topright", legend=round(a[h], digits=2), text.col ="blue", bg="white", x.intersp=0)
   on.exit(par(old.par))
   }
}

# GAMM diagnostics
# Source: Gavin Simpson, 2011
# URL: https://github.com/gavinsimpson/random_code

## Model Checking function
tsDiagGamm <- function(x, timevar, observed, f = 0.3, type = "normalized") {
    resi <- resid(x$lme, type = type)
    fits <- fitted(x$lme)
    on.exit(layout(1))
    layout(matrix(1:6, ncol = 3, byrow = TRUE))
    plot(resi ~ fits, ylab = "Normalized Residuals",
         xlab = "Fitted Values", main = "Fitted vs. Residuals")
    lines(lowess(x = fits, y = resi, f = f), col = "blue",
          lwd = 2)
    plot(resi ~ timevar, ylab = "Normalized Residuals",
         xlab = "Time", main = "Time series of residuals")
    lines(lowess(x = timevar, y = resi, f = f), col = "blue", lwd = 2)
    plot(observed ~ fits, ylab = "Observed",
         xlab = "Fitted Values", main = "Fitted vs. Observed",
         type = "n")
    abline(a = 0, b = 1, col = "red")
    points(observed ~ fits)
    lines(lowess(x = fits, y = observed, f = f), col = "blue",
          lwd = 2)
    hist(resi, freq = FALSE, xlab = "Normalized Residuals")
    qqnorm(resi)
    qqline(resi)
    acf(resi, main = "ACF of Residuals")
}

## Functions for derivatives of GAM models
Deriv <- function(mod, n = 200, eps = 1e-7, newdata) {
    if(isTRUE(inherits(mod, "list")))
        mod <- mod$gam
    m.terms <- attr(terms(mod), "term.labels")
    if(missing(newdata)) {
        newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                       function(x) seq(min(x), max(x), length = n))
        names(newD) <- m.terms
    } else {
        newD <- newdata
    }
    X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
    newD <- newD + eps
    X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
    Xp <- (X1 - X0) / eps
    Xp.r <- NROW(Xp)
    Xp.c <- NCOL(Xp)
    ## dims of bs
    bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
    # number of smooth terms
    t.labs <- attr(mod$terms, "term.labels")
    nt <- length(t.labs)
    ## list to hold the derivatives
    lD <- vector(mode = "list", length = nt)
    names(lD) <- t.labs
    for(i in seq_len(nt)) {
        Xi <- Xp * 0
        want <- grep(t.labs[i], colnames(X1))
        Xi[, want] <- Xp[, want]
        df <- Xi %*% coef(mod)
        df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
        lD[[i]] <- list(deriv = df, se.deriv = df.sd)
        ## Xi <- Xp * 0 ##matrix(0, nrow = Xp.r, ncol = Xp.c)
        ## J <- bs.dims[i]
        ## Xi[,(i-1) * J + 1:J + 1] <- Xp[,(i-1) * J + 1:J +1]
        ## df <- Xi %*% coef(mod)
        ## df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
        ## lD[[i]] <- list(deriv = df, se.deriv = df.sd)
    }
    class(lD) <- "Deriv"
    lD$gamModel <- mod
    lD$eps <- eps
    lD$eval <- newD - eps
    return(lD)
}

confint.Deriv <- function(object, term, alpha = 0.05, ...) {
    l <- length(object) - 3
    term.labs <- names(object[seq_len(l)])
    if(missing(term))
        term <- term.labs
    Term <- match(term, term.labs)
    ##term <- term[match(term, term.labs)]
    if(any(miss <- is.na(Term)))
        stop(paste("'term'", term[miss], "not a valid model term."))
    ## if(is.na(term))
    ##     stop("'term' not a valid model term.")
    res <- vector(mode = "list", length = length(term))
    names(res) <- term
    residual.df <- length(object$gamModel$y) - sum(object$gamModel$edf)
    tVal <- qt(1 - (alpha/2), residual.df)
    ## tVal <- qt(1 - (alpha/2), object$gamModel$df.residual)
    for(i in seq_along(term)) {
        upr <- object[[term[i]]]$deriv + tVal * object[[term[i]]]$se.deriv
        lwr <- object[[term[i]]]$deriv - tVal * object[[term[i]]]$se.deriv
        res[[term[i]]] <- list(upper = drop(upr), lower = drop(lwr))
    }
    res$alpha = alpha
    res
}

signifD <- function(x, d, upper, lower, eval = 0) {
    miss <- upper > eval & lower < eval
    incr <- decr <- x
    want <- d > eval
    incr[!want | miss] <- NA
    want <- d < eval
    decr[!want | miss] <- NA
    list(incr = incr, decr = decr)
}

plot.Deriv <- function(x, alpha = 0.05, polygon = TRUE,
                       sizer = FALSE, term, eval = 0, lwd = 3,
                       col = "lightgrey", border = col,
                       ylab, xlab, ...) {
    l <- length(x) - 3
    ## get terms and check specified (if any) are in model
    term.labs <- names(x[seq_len(l)])
    if(missing(term))
        term <- term.labs
    Term <- match(term, term.labs)
    if(any(miss <- is.na(Term)))
        stop(paste("'term'", term[miss], "not a valid model term."))
    if(all(is.na(Term)))
        stop("All terms in 'term' not found in model.")
    l <- sum(!miss)
    nplt <- n2mfrow(l)
    ## tVal <- qt(1 - (alpha/2), x$gamModel$df.residual)
    residual.df <- length(x$gamModel$y) - sum(x$gamModel$edf)
    tVal <- qt(1 - (alpha/2), residual.df)
    if(missing(ylab))
        ylab <- expression(italic(hat(f)*"'"*(x)))
    if(missing(xlab)) {
        xlab <- attr(terms(x$gamModel), "term.labels")[Term]
        names(xlab) <- xlab
    }
    layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
    CI <- confint(x, term = term, alpha = alpha)
    for(i in seq_along(term)) {
    ## for(i in seq_len(l)) {
        upr <- CI[[term[i]]]$upper
        lwr <- CI[[term[i]]]$lower
        ylim <- range(upr, lwr)
        plot(x$eval[,term[i]], x[[term[i]]]$deriv, type = "n",
             ylim = ylim, ylab = ylab, xlab = xlab[term[i]], ...)
        if(isTRUE(polygon)) {
            polygon(c(x$eval[,term[i]], rev(x$eval[,term[i]])),
                    c(upr, rev(lwr)), col = col, border = border)
        } else {
            lines(x$eval[,term[i]], upr, lty = "dashed")
            lines(x$eval[,term[i]], lwr, lty = "dashed")
        }
        abline(h = 0, ...)
        if(isTRUE(sizer)) {
            lines(x$eval[,term[i]], x[[term[i]]]$deriv, lwd = 1)
            S <- signifD(x[[term[i]]]$deriv, x[[term[i]]]$deriv, upr, lwr,
                         eval = eval)
            lines(x$eval[,term[i]], S$incr, lwd = lwd, col = "blue")
            lines(x$eval[,term[i]], S$decr, lwd = lwd, col = "red")
        } else {
            lines(x$eval[,term[i]], x[[term[i]]]$deriv, lwd = 2)
        }
    }
    layout(1)
    invisible(x)
}
