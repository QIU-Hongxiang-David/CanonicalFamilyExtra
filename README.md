# CanonicalFamilyExtra

Extra Canonical Link Family Objects for Generalized Linear Models (R package)

To install from Github:
```{r}
devtools::install_github("Qiu-Hongxiang-David/CanonicalFamilyExtra")
```

To install from CRAN:
```{r}
install.packages("CanonicalFamilyExtra")
```

Example:
```{r}
library(CanonicalFamilyExtra)
x <- rnorm(100)
y <- 1/(1 + exp(-x)) + rnorm(100)

# fit logistic model
glm(y~x, family = binomial_extra())

# fit log-linear model
glm(y~x, family = poisson_extra())
```

More naive approaches appear not working.
