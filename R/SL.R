#' @title SuperLearner wrapper for cv.glmnet that works with the extra families
#' @name SL.glmnet.extra
#' @description
#' A wrapper of [glmnet::cv.glmnet()] similar to [SuperLearner::SL.glmnet()], except that \code{family} can be \code{binomial_extra}. [SuperLearner::SL.glmnet()] only passes the name of \code{family} and thus cannot pass the full customized families like [binomial_extra()].
#' @param Y see [SuperLearner::SL.glmnet()]
#' @param X see [SuperLearner::SL.glmnet()]
#' @param newX see [SuperLearner::SL.glmnet()]
#' @param family similar to [SuperLearner::SL.glmnet()]. Unlike [SuperLearner::SL.glmnet()], the entire family object rather than the name of \code{family} (e.g., "gaussian", "binomial") will be passed to [glmnet::cv.glmnet()]
#' @param obsWeights see [SuperLearner::SL.glmnet()]
#' @param id see [SuperLearner::SL.glmnet()]
#' @param alpha see [SuperLearner::SL.glmnet()]
#' @param nfolds see [SuperLearner::SL.glmnet()]
#' @param nlambda see [SuperLearner::SL.glmnet()]
#' @param useMin see [SuperLearner::SL.glmnet()]
#' @param loss see [SuperLearner::SL.glmnet()]
#' @param ... see [SuperLearner::SL.glmnet()]
#' @examples
#' if (interactive()){
#'     set.seed(321)
#'     expit <- binomial()$linkinv
#'     X <- matrix(rnorm(100 * 5), nrow = 100)
#'     y <- expit(1 + X[,1]) + rnorm(100)
#'     require(SuperLearner)
#'     SL.library<-list(c("SL.glmnet.extra","screen.glmnet.extra"))
#'     SuperLearner(y, data.frame(X), family=binomial_extra(),
#'                  SL.library = SL.library, cvControl = list(V = 2))
#' }
#' @export
SL.glmnet.extra<-function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10, nlambda = 100, useMin = TRUE, loss = "deviance", ...) {
    if(!requireNamespace("glmnet",quietly = FALSE)){
        stop("loading required package (glmnet) failed",call.=FALSE)
    }
    if (!is.matrix(X)) {
        X <- model.matrix(~-1 + ., X)
        newX <- model.matrix(~-1 + ., newX)
    }
    fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, lambda = NULL, type.measure = loss, nfolds = nfolds, family = family, alpha = alpha, nlambda = nlambda, ...)
    pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, "lambda.min", "lambda.1se"))
    fit <- list(object = fitCV, useMin = useMin)
    class(fit) <- "SL.glmnet"
    out <- list(pred = pred, fit = fit)
    return(out)
}


#' @title SuperLearner screener using cv.glmnet that works with the extra families
#' @name screen.glmnet.extra
#' @description
#' A wrapper of [glmnet::cv.glmnet()] similar to [SuperLearner::screen.glmnet()], except that \code{family} can be \code{binomial_extra}. [SuperLearner::screen.glmnet()] only passes the name of \code{family} (e.g., "gaussian", "binomial") and thus cannot pass the full customized families like [binomial_extra()].
#' @param Y see [SuperLearner::screen.glmnet()]
#' @param X see [SuperLearner::screen.glmnet()]
#' @param family similar to [SuperLearner::screen.glmnet()]. Unlike [SuperLearner::screen.glmnet()], the entire family object rather than the name of \code{family} (e.g., "gaussian", "binomial") will be passed to [glmnet::cv.glmnet()]
#' @param alpha see [SuperLearner::screen.glmnet()]
#' @param minscreen see [SuperLearner::screen.glmnet()]
#' @param nfolds see [SuperLearner::screen.glmnet()]
#' @param nlambda see [SuperLearner::screen.glmnet()]
#' @param ... see [SuperLearner::SL.glmnet()]
#' @examples
#' if (interactive()){
#'     set.seed(321)
#'     expit <- binomial()$linkinv
#'     X <- matrix(rnorm(100 * 5), nrow = 100)
#'     y <- expit(1 + X[,1]) + rnorm(100)
#'     require(SuperLearner)
#'     SL.library<-list(c("SL.glmnet.extra","screen.glmnet.extra"))
#'     SuperLearner(y, data.frame(X), family=binomial_extra(),
#'                  SL.library = SL.library, cvControl = list(V = 2))
#' }
#' @export
screen.glmnet.extra<-function(Y, X, family, alpha = 1, minscreen = 2, nfolds = 10, nlambda = 100, ...) {
    if(!requireNamespace("glmnet",quietly = FALSE)){
        stop("loading required package (glmnet) failed",call.=FALSE)
    }
    if (!is.matrix(X)) {
        X <- model.matrix(~-1 + ., X)
    }
    fitCV <- glmnet::cv.glmnet(x = X, y = Y, lambda = NULL, type.measure = "deviance", nfolds = nfolds, family = family, alpha = alpha, nlambda = nlambda)
    whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] != 0)
    if (sum(whichVariable) < minscreen) {
        warning("fewer than minscreen variables passed the glmnet screen, increased lambda to allow minscreen variables")
        sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta), 2, function(x) sum((x != 0)))
        newCut <- which.max(sumCoef >= minscreen)
        whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[, newCut] != 0)
    }
    return(whichVariable)
}
