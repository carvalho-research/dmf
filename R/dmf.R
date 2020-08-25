# [ Deviance Matrix Factorization ]

# TODO:
# `dmf(x ~ row(1) + column(1 + z), rank, offset, symmetric?, use_diagonal?)`
# `row(1) = tcrossprod(ones(n), mu)`, mu has dimension p, equivalent to center
# `column(1) = tcrossprod(mu, ones(p))`, mu has dimension n
# In either case L and V have to augmented with `ones` and `mu`
# TODO:
# - scree plot to define rank, LOO CV
# - multinomial: alpha_k + l_i * v_j or, with centering,
#     alpha_{k, j} + l_i * v_j
# - correlation with weights
# - merge normal and symmetric cases, add option for symmetric

# for non-negative matrix factorization
ident_poisson <- function (...) {
  fam <- quasipoisson(link = "identity", ...)
  #fam$initialize <- expression({
  #  if (any(y < 0))
  #    stop("negative values not allowed for the 'quasiPoisson' family")
  #  mustart <- y + epsilon * max(y)
  #})
  fam$variance <- fam$linkinv <- function (mu) pmax(mu, .Machine$double.eps)
  fam
}

zeros <- function (x) matrix(0, nrow = nrow(x), ncol = ncol(x))
symmetrize <- function (x) .5 * (x + t(x))
norm2 <- function (x) norm(as.matrix(x[!is.na(x)]), "2")
normalize <- function (x, margin = 2)
  sweep(x, margin, apply(x, margin, norm2), `/`)
gsym_solve <- function (A, b, tol = sqrt(.Machine$double.eps)) {
  ea <- eigen(A, symmetric = TRUE)
  V <- ea$vectors; d <- ea$values
  valid <- d > max(tol * d[1], 0)
  if (!all(valid)) {
    V <- V[, valid, drop = FALSE]; d <- d[valid]
  }
  V %*% sweep(crossprod(V, b), 1, d, `/`)
}

family_initialize <- function (x, weights, family = gaussian()) {
  nr <- nrow(x); nc <- ncol(x)
  etastart <- NULL; start <- NULL; mustart <- NULL
  y <- c(x); nobs <- length(y) #; weights <- rep(1, nobs)
  eval(family$initialize)
  matrix(mustart, nrow = nr, ncol = nc)
}


#' Center \code{dmf}.
#'
#' @param lv DMF structure to be centered.
#' @param x0 Intercept column.
#' @param reorder Reorder columns according to norm of Lambda?
#' @return Centered DMF structure with \code{center} attribute.
#' @export
dmf_center <- function (lv, x0 = rep(1, nrow(lv$L)), reorder = TRUE) {
  q0 <- qr(x0)
  Lc <- qr.resid(q0, lv$L); Vc <- lv$V
  if (reorder) {
    so <- order(apply(Lc, 2, norm2), decreasing = TRUE)
    Lc <- Lc[, so]; Vc <- Vc[, so]
  }
  sl <- svd(Lc) # re-orthogonalize
  Lc <- sweep(sl$u, 2, sl$d, `*`)
  Vc <- Vc %*% sl$v
  list(L = Lc, V = Vc, deviance = lv$deviance, family = lv$family,
       center = drop(tcrossprod(lv$V, qr.coef(q0, lv$L))))
}

#' Normalize sign of \code{dmf}: the first row of Lambda in \code{dmf} is non-negative.
#'
#' @param lv DMF structure to be sign-normalized.
#' @return Sign-normalized DMF structure.
#' @export
dmf_sign <- function (lv) {
  s <- sign(lv$L[1, ])
  lv$L <- sweep(lv$L, 2, s, `*`)
  lv$V <- sweep(lv$V, 2, s, `*`)
  lv
}

# NOTE: weight[i, j] == 0 means "don't use (i,j)"; is.na(x[i, j]) means "missing, estimate it"
#' Perform deviance matrix factorization.
#'
#' @param x Input matrix to be factorized.
#' @param family Family object to specify deviance loss.
#' @param rank Decomposition rank.
#' @param weights Entrywise weight.
#' @param offset Entrywise offset.
#' @param control Algorithm control parameters (see \code{glm.control}).
#' @return DMF structure.
#' @export
dmf <- function (x, family = gaussian(),
                 rank = ncol(x), weights = 1, offset = zeros(x),
                 control = glm.control(epsilon = 1e-6, maxit = 100)) {
  n <- nrow(x); p <- ncol(x)
  rank <- as.integer(rank)
  if (n < p) stop("fewer observations than predictors")
  if (rank <= 0) stop("invalid non-positive rank")
  if (rank > p) stop("rank is larger than number of predictors")
  if (length(weights) == 1)
    weights <- rep(weights, n * p)
  else {
    if (is.vector(weights)) {
      if (length(weights) != n)
        stop("inconsistent number of weights")
      else
        weights <- rep(weights, p)
    } else {
      if (nrow(weights) != n && ncol(weights) != p)
        stop("inconsistent number of weights")
    }
  }
  valid <- (weights > 0) & (!is.na(x))
  if (any(apply(!valid, 1, all)) || any(apply(!valid, 2, all)))
    stop("full row or column with zero weights or missing values")

  # [ initialize ]
  mu <- family_initialize(x, weights, family)
  eta <- family$linkfun(mu)
  eta[!valid] <- mu[!valid] <- 0
  L <- (eta - offset)[, 1:rank, drop = FALSE] # V = I_p
  V <- matrix(nrow = p, ncol = rank)
  for (it in 1:control$maxit) {
    mu_eta <- matrix(family$mu.eta(eta), nrow = n, ncol = p)
    var <- matrix(family$variance(mu), nrow = n, ncol = p)
    # W = mu_eta / var * (x - mu) * weights
    is_inf_mu <- is.infinite(mu)
    S <- mu_eta / var * mu_eta * weights
    S[is_inf_mu] <- 0
    Z <- eta - offset + (x - mu) / mu_eta # working residuals
    Z[is_inf_mu] <- eta[is_inf_mu] - offset[is_inf_mu]
    Z[!valid] <- 0

    L <- normalize(L)
    for (j in 1:p) {
      aj <- valid[, j]
      sj <- sqrt(S[, j]); Lj <- sweep(L, 1, sj, `*`)[aj, , drop = FALSE]
      V[j, ] <- gsym_solve(crossprod(Lj), crossprod(Lj, (Z[, j] * sj)[aj]))
    }

    V <- normalize(V)
    for (i in 1:n) {
      ai <- valid[i, ]
      si <- sqrt(S[i, ]); Vi <- sweep(V, 1, si, `*`)[ai, , drop = FALSE]
      L[i, ] <- gsym_solve(crossprod(Vi), crossprod(Vi, (Z[i, ] * si)[ai]))
    }

    eta <- tcrossprod(L, V) + offset
    mu <- family$linkinv(eta)
    dr <- family$dev.resids(x[valid], mu[valid], weights[valid])
    dr[is.na(dr) | is.nan(dr)] <- 0
    dev_new <- sum(dr)
    if (it > 1 && (dev - dev_new) / (dev + .1) < control$epsilon) break
    if (control$trace) message("[", it, "] dev = ", dev_new)
    dev <- dev_new
    L_old <- L; V_old <- V
  }
  # [ orthogonalize ]
  slv <- svd(tcrossprod(L_old, V_old))
  L <- sweep(slv$u[, 1:rank, drop = FALSE], 2, slv$d[1:rank], `*`)
  V <- slv$v[, 1:rank, drop = FALSE]

  list(L = L, V = V, deviance = dev, family = family)
}


# NOTE: factorization here is .5 * (L * V' + V * L')
#' Perform a symmetric deviance matrix factorization, (L * V' + V * L') / 2.
#'
#' @param x Input matrix to be factorized.
#' @param family Family object to specify deviance loss.
#' @param rank Decomposition rank.
#' @param weights Entrywise weight.
#' @param offset Entrywise offset.
#' @param control Algorithm control parameters (see \code{glm.control}).
#' @return DMF structure.
#' @export
dmf_symm <- function (x, family = gaussian(),
                      rank = ncol(x), weights = 1, offset = zeros(x),
                      maxit = 100, epsilon = 1e-6, trace = TRUE) {
  n <- nrow(x)
  if (n != ncol(x)) stop("matrix is not square")
  rank <- as.integer(rank)
  if (rank <= 0) stop("invalid non-positive rank")
  if (rank > n) stop("rank is larger than number of predictors")
  if (length(weights) == 1)
    weights <- rep(weights, n * n)
  else {
    if (is.vector(weights)) {
      if (length(weights) != n)
        stop("inconsistent number of weights")
      else
        weights <- rep(weights, n)
    } else {
      if (nrow(weights) != n && ncol(weights) != n)
        stop("inconsistent number of weights")
    }
  }
  valid <- (weights > 0) & (!is.na(x))
  if (any(apply(!valid, 1, all)) || any(apply(!valid, 2, all)))
    stop("full row or column with zero weights or missing values")

  # [ initialize ]
  mu <- family_initialize(x, weights, family)
  eta <- family$linkfun(mu)
  eta[!valid] <- mu[!valid] <- 0
  L <- (eta - offset)[, 1:rank, drop = FALSE] # V = I_n
  V <- matrix(nrow = n, ncol = rank)
  for (it in 1:maxit) {
    mu_eta <- matrix(family$mu.eta(eta), nrow = n, ncol = n)
    var <- matrix(family$variance(mu), nrow = n, ncol = n)
    W <- symmetrize(mu_eta / var * (x - mu) * weights)
    S <- symmetrize(mu_eta ^ 2 / var * weights)
    Z <- eta - offset + W / S # working residuals
    Z[!valid] <- 0

    L <- normalize(L)
    for (j in 1:n) {
      aj <- valid[, j]
      sj <- sqrt(S[, j]); Lj <- sweep(L, 1, sj, `*`)[aj, , drop = FALSE]
      V[j, ] <- gsym_solve(crossprod(Lj), crossprod(Lj, (Z[, j] * sj)[aj]))
    }

    V <- normalize(V)
    for (i in 1:n) {
      ai <- valid[i, ]
      si <- sqrt(S[i, ]); Vi <- sweep(V, 1, si, `*`)[ai, , drop = FALSE]
      L[i, ] <- gsym_solve(crossprod(Vi), crossprod(Vi, (Z[i, ] * si)[ai]))
    }

    eta <- .5 * (tcrossprod(L, V) + tcrossprod(V, L)) + offset
    mu <- family$linkinv(eta)
    dev_new <- sum(family$dev.resids(x, mu, weights))
    if (it > 1 && (dev - dev_new) / (dev + .1) < epsilon) break
    if (trace) message("[", it, "] dev = ", dev_new)
    dev <- dev_new
    L_old <- L; V_old <- V
  }
  # [ orthogonalize ]
  slv <- svd(tcrossprod(L_old, V_old))
  L <- sweep(slv$u[, 1:rank, drop = FALSE], 2, slv$d[1:rank], `*`)
  V <- slv$v[, 1:rank, drop = FALSE]

  list(L = L, V = V, deviance = dev, family = family)
}


# [ facilities ]
dmf_fitted <- function (lv, r = ncol(lv$L), subset = seq_len(nrow(lv$L)),
                        linear = FALSE) {
  eta <- tcrossprod(lv$L[subset, 1:r, drop = FALSE],
                    lv$V[, 1:r, drop = FALSE])
  if (!is.null(lv$center)) eta <- sweep(eta, 2, lv$center, `+`)
  if (linear) eta else lv$family$linkinv(eta)
}

dmf_residuals <- function (x, lv, r = ncol(lv$L), subset = seq_len(nrow(lv$L)),
                           pearson = FALSE) {
  eta <- tcrossprod(lv$L[subset, 1:r, drop = FALSE],
                    lv$V[, 1:r, drop = FALSE])
  if (!is.null(lv$center)) eta <- sweep(eta, 2, lv$center, `+`)
  mu <- lv$family$linkinv(eta)
  if (pearson) (x[subset,] - mu) / sqrt(lv$family$variance(mu)) else
  # FIXME: weights
    sign(x[subset,] - mu) * sqrt(lv$family$dev.resids(x[subset,], mu, 1))
}

dmf_predict <- function (lv, x, linear = FALSE) {
  offset <- if (!is.null(lv$center)) lv$center else rep(0, length(x))
  gf <- glm.fit(lv$V, x, family = lv$family, offset = offset)
  if (linear) gf$linear.predictor else gf$fitted
}

dmf_aic <- function (lv, scale = 2) {
  n <- nrow(lv$L); p <- nrow(lv$V); q <- ncol(lv$L) # == ncol(lv$V)
  if (is.na(scale)) scale <- log(n * p)
  scale * ((n + p) * q - q ^ 2) + lv$deviance
}

