# Deviance Matrix Factorization

This R package computes _deviance matrix factorizations_ (DMFs) of data
matrices, a generalization of singular value decompositions to data entries
that follow exponential family distributions such as binomial and Poisson.

```r
library(dmf)
x <- as.matrix(iris[, -5])
lv <- dmf(x, family = quasipoisson(), rank = 1) |> dmf_center()
plot(lv$L[, 1], col = iris[, 5])
```

