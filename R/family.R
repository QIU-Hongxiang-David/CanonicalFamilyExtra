#' @importFrom stats binomial make.link poisson quasi

clip_interval<-function(x,lower=-Inf,upper=Inf){
    out<-pmax(pmin(x,upper),lower)
    if(is.matrix(x)){
        out<-matrix(out,nrow=nrow(x),ncol=ncol(x))
    }
    out
}


#' @title Family object for fitting binomial (e.g., logistic, probit) regression model with potentially unbounded continuous outcome
#' @name binomial_extra
#' @description
#' A family object for fitting binomial models (i.e., generalized linear models with range contained in the open unit interval \eqn{(0,1)}) with continuous outcomes that may fall outside the unit interval \eqn{[0,1]}. Also works with \code{\link[glmnet:glmnet]{glmnet}}.
#' @details
#' This family is useful, for example, when the estimand is a conditional probability function while the outcome is a transformed psudo-outcome so that the estimator is multiply robust, or estimating a regression function with known bounds while the outcome might not respect the known bounds. Naive approaches such as \code{glm(family=binomial())}, \code{glm(family=quasibinomial())}, \code{glm(family=gaussian(link="logit"))}, \code{glm(family=quasi(link="logit",variance="constant"))} etc. might not work appropriately in such cases.
#'
#' Particularly for logistic model, because of using the binomial working likelihood and its canonical link, the model fitting is a convex problem and does not depend on starting value.
#'
#' The output has \code{family="gaussian"} by default to be compatible with other learners in \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}, because when the outcome is continuous, other learners might not perform correctly with \code{family="binomial"}.
#'
#' @section Warning:
#'
#' This function tweaks basic family objects and might remove some safety features. USE WITH CARE!!!
#' @param link see \code{\link[stats:family]{family}}. Default to \code{"logit"}
#' @param variance see \code{\link[stats:family]{family}}. Default to \code{"mu(1-mu)"}, same as logistic regression.
#' @param family The family of the returned family object. Either \code{"gaussian"} or \code{"binomial"}. Default to \code{"gaussian"}. Seems not to matter for glm.
#' @returns a family object
#'
#' @examples
#' set.seed(123)
#' expit <- binomial()$linkinv
#'
#' # glm
#' x <- rnorm(100)
#' y <- expit(1 + x) + rnorm(100)
#' glm(y~x, family = binomial_extra()) # or family=binomial_extra, or family="binomial_extra"
#' # Errors or not so reliable
#' \dontrun{
#' glm(y~x, family = binomial())
#' glm(y~x, family = quasibinomial())
#' glm(y~x, family = gaussian(link = "logit"))
#' # setting starting value might work, but is non-convex and requires starting value
#' glm(y~x, family = gaussian(link = "logit"), start = c(-1,0))
#' glm(y~x, family = quasi(link = "logit", variance = "constant"))
#' }
#'
#' #glmnet
#' X <- matrix(rnorm(100 * 5), nrow = 100)
#' y <- expit(1 + X[,1]) + rnorm(100)
#' require(glmnet)
#' glmnet(X, y, family = binomial_extra())
#' # or family=binomial_extra; cannot use family="binomial_extra"
#' # Errors
#' \dontrun{
#' glmnet(X, y, family = binomial())
#' glmnet(X, y, family = gaussian(link = "logit"))
#' glmnet(X, y, family = quasi(link = "logit", variance = "constant"))
#' }
#'
#' # other links/variance for glm
#' x <- rnorm(100)
#' y <- expit(1 + x) + rnorm(100)
#' glm(y~x, family = binomial_extra("probit")) # probit regression
#' glm(y~x, family = binomial_extra(variance = "constant")) # least squares
#'
#' # within SuperLearner
#' X <- matrix(rnorm(100 * 3), nrow = 100)
#' y <- expit(1 + X[,1]) + rnorm(10)
#' require(SuperLearner)
#' SuperLearner(y, data.frame(X), family=binomial_extra(),
#'              SL.library = c("SL.glm", "SL.ipredbagg"), cvControl = list(V = 2))
#' # Error in SL.ipredbagg because of wrong family
#' \dontrun{
#' SuperLearner(y, data.frame(X), family = binomial_extra(family = "binomial"),
#'              SL.library = c("SL.glm", "SL.ipredbagg"), cvControl = list(V = 2))
#' }
#'
#'
#' @export
binomial_extra<-function(link="logit",variance="mu(1-mu)",family=c("gaussian","binomial")){
    linktemp <- substitute(link)
    if (!is.character(linktemp))
        linktemp <- deparse(linktemp)
    if (linktemp %in% c("logit", "probit", "cloglog", "identity", "inverse", "log", "1/mu^2", "sqrt")){
        stats <- make.link(linktemp)
    }else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    } else {
        stats <- link
        linktemp <- stats$name %||% deparse(linktemp)
    }

    family<-match.arg(family)

    variancetemp <- substitute(variance)
    if (!is.character(variancetemp))
        variancetemp <- deparse(variancetemp)

    eps<-.Machine$double.eps

    out<-do.call(quasi,list(link=linktemp,variance=variancetemp))
    out$family<-family
    out$linkfun<-function(mu){
        .Call(stats:::C_logit_link, clip_interval(mu,eps,1-eps))
    }
    out$linkinv<-function(eta){
        clip_interval(.Call(stats:::C_logit_linkinv, eta),eps,1-eps)
    }
    out$initialize<-expression({
        n <- rep.int(1, nobs)
        mustart <- pmax(0.001, pmin(0.999, y))
    })
    out$validmu<-binomial()$validmu
    if(is.character(variancetemp) && variancetemp=="mu(1-mu)"){
        out$dev.resids<-function(y,mu,wt){
            mu<-clip_interval(mu,eps,1-eps)
            -2*wt*(y*log(mu)+(1-y)*log1p(-mu))
        }
    }
    out
}







#' @title Family object for fitting log-linear model with potentially unbounded continuous outcome, particularly with Poisson working likelihood
#' @name poisson_extra
#' @description
#' A family object for fitting Poisson models (i.e., generalized linear models with range contained in the open unit interval \eqn{(0,\infty)}) with continuous outcomes that may not be non-negative integers. Also works with \code{\link[glmnet:glmnet]{glmnet}}.
#' @details
#' This family is useful, for example, when the estimand is a conditional probability function while the outcome is a transformed psudo-outcome so that the estimator is multiply robust, or estimating a positive regression function while the outcome might be negative or non-integers. Naive approaches such as \code{glm(family=poisson())}, \code{glm(family=quasipoisson())}, \code{glm(family=gaussian(link="log"))}, \code{glm(family=quasi(link="log",variance="constant"))} etc. might not work appropriately or reliably in such cases.
#'
#' Particularly for log-linear model, because of using the Poisson working likelihood and its canonical link, the model fitting is a convex problem and does not depend on starting value.
#'
#' The output has \code{family="gaussian"} by default to be compatible with other learners in \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}, because when the outcome is continuous, other learners might not perform correctly with \code{family="poisson"}.
#'
#' @section Warning:
#'
#' This function tweaks basic family objects and might remove some safety features. \code{dev.resids} of the family object should not be interpreted as the usual deviance residual for statistical inference, but is -2 times the working log likelihood and only for optimization and model fitting. USE WITH CARE!!!
#' @param link see \code{\link[stats:family]{family}}. Default to \code{"log"}
#' @param variance see \code{\link[stats:family]{family}}. Default to \code{"mu"}, same as Poisson regression.
#' @param family The family of the returned family object. Either \code{"gaussian"} or \code{"poisson"}. Default to \code{"gaussian"}. Seems not to matter for glm.
#' @returns a family object
#'
#' @examples
#' set.seed(123)
#'
#' #########
#' # negative outcomes
#' #########
#' # glm
#' x <- rnorm(100)
#' y <- exp(-1 + x) + rnorm(100)
#' glm(y~x, family = poisson_extra()) # or family=poisson_extra, or family="poisson_extra"
#' # Errors or not so reliable
#' \dontrun{
#' glm(y~x, family = poisson())
#' glm(y~x, family = quasipoisson())
#' glm(y~x, family = gaussian(link = "log"))
#' # setting starting value might work, but is non-convex and requires starting value
#' glm(y~x, family = gaussian(link = "log"), start = c(-1,0))
#' glm(y~x, family = quasi(link = "log", variance = "constant"))
#' # setting starting value might work, but is non-convex and requires starting value
#' glm(y~x, family = quasi(link = "log", variance = "constant"), start = c(-1,0))
#' }
#'
#' #glmnet
#' X <- matrix(rnorm(100 * 5), nrow = 100)
#' y <- exp(1 + X[,1]) + rnorm(100)
#' require(glmnet)
#' glmnet(X, y, family = poisson_extra())
#' # or family=poisson_extra; cannot use family="poisson_extra"
#'
#' # Errors
#' \dontrun{
#' glmnet(X, y, family = poisson())
#' glmnet(X, y, family = gaussian(link = "log"))
#' glmnet(X, y, family = quasi(link = "log", variance = "constant"))
#' }
#'
#' # within SuperLearner
#' X <- matrix(rnorm(100 * 3), nrow = 100)
#' y <- exp(1 + X[,1]) + rnorm(10)
#' require(SuperLearner)
#' SuperLearner(y, data.frame(X), family=poisson_extra(),
#'              SL.library = c("SL.glm", "SL.ipredbagg"), cvControl = list(V = 2))
#' # Error in SL.ipredbagg because of wrong family
#' \dontrun{
#' SuperLearner(y, data.frame(X), family = poisson_extra(family = "poisson"),
#'              SL.library = c("SL.glm", "SL.ipredbagg"), cvControl = list(V = 2))
#' }
#'
#'
#' ########
#' # positive non-integer outcomes
#' ########
#' x <- rnorm(100)
#' y <- rexp(100, exp(1 - x))
#' glm(y~x, family = poisson_extra())
#' # Errors or not so reliable
#' \dontrun{
#' glm(y~x, family = poisson())
#' glm(y~x, family = quasipoisson())
#' # might work but is non-convex, so may depend on starting value
#' glm(y~x, family = gaussian(link = "log"))
#' # might work but is non-convex, so may depend on starting value
#' glm(y~x, family = quasi(link = "log", variance = "constant"))
#' }
#'
#'
#' @export
poisson_extra<-function(link="log",variance="mu",family=c("gaussian","poisson")){
    linktemp <- substitute(link)
    if (!is.character(linktemp))
        linktemp <- deparse(linktemp)
    if (linktemp %in% c("log", "identity", "sqrt")){
        stats <- make.link(linktemp)
    }else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    } else {
        stats <- link
        linktemp <- stats$name %||% deparse(linktemp)
    }
    family<-match.arg(family)

    variancetemp <- substitute(variance)
    if (!is.character(variancetemp))
        variancetemp <- deparse(variancetemp)

    eps<-.Machine$double.eps

    out<-do.call(quasi,list(link=linktemp,variance=variancetemp))
    out$family<-family
    out$linkfun<-function(mu){
        log(clip_interval(mu,eps))
    }
    out$linkinv<-function(eta){
        clip_interval(exp(eta),eps)
    }
    out$initialize<-expression({
        n <- rep.int(1, nobs)
        mustart <- pmax(y, .0001)
    })
    out$validmu<-poisson()$validmu
    if(is.character(variancetemp) && variancetemp=="mu"){
        out$dev.resids<-function(y,mu,wt){
            mu<-clip_interval(mu,eps)
            2 * wt * (-y * log(mu) + mu)
        }
    }
    out
}
