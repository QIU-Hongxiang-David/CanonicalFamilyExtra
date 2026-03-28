#' @title SuperLearner wrapper for cv.glmnet that works with the extra families
#' @name SL.glmnet.extra
#' @description
#' A wrapper of \code{\link[glmnet:cv.glmnet]{cv.glmnet}} similar to \code{\link[SuperLearner:SL.glmnet]{SL.glmnet}}, except that \code{family} can be \code{binomial_extra}. \code{\link[SuperLearner:SL.glmnet]{SL.glmnet}} only passes the name of \code{family} and thus cannot pass the full customized families like \code{\link[binomial_extra]{binomial_extra}}.
#' @examples
#' set.seed(123)
#' expit <- binomial()$linkinv
#' X <- matrix(rnorm(100 * 5), nrow = 100)
#' y <- expit(1 + X[,1]) + rnorm(100)
#' require(SuperLearner)
#' SuperLearner(y, data.frame(X), family=binomial_extra(), SL.library = c("SL.glmnet.extra"), cvControl = list(V = 2))
#' 
SL.glmnet.extra<-function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10, nlambda = 100, useMin = TRUE, loss = "deviance", ...) {
    SuperLearner:::.SL.require("glmnet")
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
