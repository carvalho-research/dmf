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

## References

* Wang, Liang and Carvalho, Luis E (2023). ["Deviance matrix factorization"](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-17/issue-2/Deviance-matrix-factorization/10.1214/23-EJS2174.full). Electronic Journal of Statistics 17(2): 3762-3810. DOI: 10.1214/23-EJS2174. 

