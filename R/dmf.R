# [ Deviance Matrix Factorization ]
#' @importFrom utils tail
#' @importFrom stats coef cor cov gaussian glm.control glm.fit lm quantile quasipoisson rbinom rgamma rnorm rpois
#' @importFrom dplyr %>% group_by mutate
#' @importFrom MASS rnegbin
NULL
#> NULL
#'
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
  s <- sign(lv$V[1, ])
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
#' @param start List containing initial `L` and `V` matrices as named entries.
#' @param init_svd Should the initial `L` and `V` be computed from a truncated SVD?
#' @param control Algorithm control parameters (see \code{glm.control}).
#' @return DMF structure.
#' @export
dmf <- function (x, family = gaussian(),
                 rank = ncol(x), weights = 1, offset = zeros(x),
                 start = NULL, init_svd = TRUE,
                 control = glm.control(epsilon = 1e-6, maxit = 100)) {
  n <- nrow(x); p <- ncol(x)
  rank <- as.integer(rank)
  if (n < p) stop("fewer observations than predictors")
  if (rank <= 0) stop("invalid non-positive rank")
  if (rank > p) stop("rank is larger than number of predictors")
  if (length(weights) == 1)
    weights <- rep(weights, n * p)else {
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
  dim(weights) <- dim(x)
  valid <- (weights > 0) & (!is.na(x))
  if (any(apply(!valid, 1, all)) || any(apply(!valid, 2, all)))
    stop("full row or column with zero weights or missing values")


  if (is.null(start)) {
    # [ initialize ]
    mu <- family_initialize(x, weights, family)
    eta <- family$linkfun(mu)
    #[transformation from init can gets negative]
    # valid <- valid & family$valideta(eta) ?
    if (init_svd) {
      se <- RSpectra::svds(eta - offset, rank)
      L <- sweep(se$u, 2, se$d, `*`); V <- se$v
    } else {
      L <- (eta - offset)[, 1:rank, drop = FALSE]
      V <- matrix(0, nrow = p, ncol = rank); diag(V) <- 1
    }
  } else {
    L <- start$L; V <- start$V
    if (nrow(L) != n || ncol(L) != rank)
      stop("dimensions of L are inconsistent")
    if (nrow(V) != p || ncol(V) != rank)
      stop("dimensions of V are inconsistent")
    eta <- tcrossprod(L, V) + offset
    mu <- family$linkinv(eta)
  }


  # [ iterate ]
  for (it in 1:control$maxit) {
    mu_eta <- matrix(family$mu.eta(eta), nrow = n, ncol = p)
    var <- matrix(family$variance(mu), nrow = n, ncol = p)
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
    dr[is.na(dr) | is.nan(dr)] <- 0 #[infinite mu setting dr to 0]
    dev_new <- sum(dr)
    if (it > 1 && (dev - dev_new) / (dev + .1) < control$epsilon) break
    if (control$trace) message("[", it, "] dev = ", dev_new)
    dev <- dev_new
    L_old <- L; V_old <- V
  }
  # [ orthogonalize ]
  ql <- qr(L_old); qv <- qr(V_old)
  sr <- svd(tcrossprod(qr.R(ql), qr.R(qv)))
  L <- qr.Q(ql) %*% sweep(sr$u, 2, sr$d, `*`)
  V <- qr.Q(qv) %*% sr$v

  list(L = L, V = V, deviance = dev, family = family,
       prior.weights = weights, iter = it)
}


# NOTE: factorization here is .5 * (L * V' + V * L')
#' Perform a symmetric deviance matrix factorization, (L * V' + V * L') / 2.
#'
#' @param x Input matrix to be factorized.
#' @param family Family object to specify deviance loss.
#' @param rank Decomposition rank.
#' @param weights Entrywise weight.
#' @param offset Entrywise offset.
#' @param maxit integer giving the maximal number of IWLS iterations.
#' @param epsilon positive convergence tolerance `\epsilon`; the iterations converge.
#' @param trace logical indicating if output should be produced for each iteration.
#' @return DMF structure.
#'
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


# assume `y` and `mu` are n x k x p, `wt` is n x p
y_log_y <- function (y, mu) ifelse(y != 0., y * log(y / mu), 0.)
mn_deviance <- function (y, mu, wt) {
  n <- dim(y)[1]; k <- dim(y)[2]; p <- dim(y)[3]
  dr <- 0.
  for (i in seq_len(n)) {
    for (j in seq_len(p)) {
      yij <- y[i, , j]; muij <- mu[i, , j]
      dr <- dr + 2 * wt[i, j] * (sum(y_log_y(yij, muij)) +
                                 y_log_y(1 - sum(yij), 1 - sum(muij)))
    }
  }
  dr
}


# X_[i, ., j] ~ind MN(w_ij, pi_ij)
# log(pi_ijl / pi_ijk) = (L_l)_i' * V_j
# where each L_l is n x q and V is p x q
# For identifiability, L~ = vec_k(L) is nk x q and has orthogonal columns and
# V in Stiefel(p, q), that is, L~ has the (incomplete) generalized left vectors
# up to scale, L~ * V' = U * D * V~' is a SVD.

#' Perform a deviance matrix factorization, for a multinomial tensor.
#'
#' @param x Input tensor to be factorized, of dimensions (n x k) x p.
#' @param rank Decomposition rank.
#' @param control Algorithm control parameters (see \code{glm.control}).
#' @return DMF structure, where L is (n x k) x rank and V is p x rank.
#' @export
dmf_multinom <- function (x, rank = dim(x)[3],
                          control = glm.control(epsilon = 1e-6, maxit = 100)) {
  n <- dim(x)[1]; k <- dim(x)[2] - 1; p <- dim(x)[3]
  rank <- as.integer(rank)

  # [ initialize ]
  eta <- array(dim = c(n, k, p))
  mu <- array(dim = c(n, k, p))
  wt <- apply(x, c(1, 3), sum) # counts
  x <- x[, -(k + 1), , drop = FALSE] # remove last level as reference
  for (i in seq_len(n)) {
    for (j in seq_len(p)) {
      mu[i, , j] <- (x[i, , j] + .5) / (1 + wt[i, j] + .5 * k)
      x[i, , j] <- x[i, , j] / wt[i, j]
      eta[i, , j] <- log(mu[i, , j]) - log(1 - sum(mu[i, , j]))
    }
  }
  L <- eta[, , 1:rank, drop = FALSE] # V = I_p
  V <- matrix(nrow = p, ncol = rank)
  dev <- mn_deviance(x, mu, wt)

  # [ iterate ]
  for (it in 1:control$maxit) {
    L <- array(normalize(matrix(L, n * k, rank)), dim = dim(L))
    for (j in seq_len(p)) {
      z <- numeric(rank); H <- matrix(0., rank, rank)
      for (i in seq_len(n)) {
        xij <- matrix(L[i, ,], nrow = k, ncol = rank); yij <- x[i, , j]
        etaij <- eta[i, , j]; muij <- mu[i, , j]
        W <- diag(muij, nrow = k, ncol = k) - tcrossprod(muij)
        z <- z + wt[i, j] * crossprod(xij, yij - muij + drop(W %*% etaij))
        H <- H + wt[i, j] * crossprod(xij, W %*% xij)
      }
      V[j, ] <- gsym_solve(H, z)
    }

    V <- normalize(V)
    for (i in seq_len(n)) {
      z <- numeric(k * rank); H <- matrix(0., k * rank, k * rank)
      for (j in seq_len(p)) {
        xij <- kronecker(diag(k), V[j, , drop = FALSE]); yij <- x[i, , j]
        etaij <- eta[i, , j]; muij <- mu[i, , j]
        W <- diag(muij, nrow = k, ncol = k) - tcrossprod(muij)
        z <- z + wt[i, j] * crossprod(xij, yij - muij + drop(W %*% etaij))
        H <- H + wt[i, j] * crossprod(xij, W %*% xij)
      }
      L[i, ,] <- matrix(gsym_solve(H, z), nrow = k, byrow = TRUE)
    }

    for (i in seq_len(n)) {
      etai <- tcrossprod(L[i, ,], V)
      for (j in seq_len(p)) {
        muj <- exp(c(etai[, j], 0))
        mu[i, , j] <- muj[-(k + 1)] / sum(muj)
      }
      eta[i, , ] <- etai
    }
    dev_new <- mn_deviance(x, mu, wt)
    if (it > 1 && (dev - dev_new) / (dev + .1) < control$epsilon) break
    if (control$trace) message("[", it, "] dev = ", dev_new)
    dev <- dev_new
    L_old <- L; V_old <- V
  }
  # [ orthogonalize ]
  slv <- svd(tcrossprod(matrix(L_old, n * k, rank), V_old)) # generalized svd
  L <- array(sweep(slv$u[, 1:rank, drop = FALSE], 2, slv$d[1:rank], `*`), dim = dim(L))
  V <- slv$v[, 1:rank, drop = FALSE]

  list(L = L, V = V, deviance = dev)
}



# [ facilities ]
svd_rank1update <- function (sx, l, v) {
  q <- ncol(sx$u) # == ncol(sx$v)
  nl <- norm2(l); nv <- norm2(v)
  delta <- nl * nv; l <- l / nl; v <- v / nv
  ul <- l - sx$u %*% crossprod(sx$u, l); nl <- norm2(ul)
  uv <- v - sx$v %*% crossprod(sx$v, v); nv <- norm2(uv)

  B <- delta * tcrossprod(c(crossprod(sx$u, l), nl),
                          c(crossprod(sx$v, v), nv))
  diag(B)[1:q] <- diag(B)[1:q] + sx$d[1:q]
  sb <- svd(B)
  list(d = sb$d,
       u = cbind(sx$u, ul / nl) %*% sb$u,
       v = cbind(sx$v, uv / nv) %*% sb$v)
}

#' Increases the rank of a DMF.
#'
#' @param dx DMF structure for x of rank q.
#' @param x The matrix that was partially decomposed.
#' @param adjust_svd Should the DMF be orthogonalized after increasing rank?
#' @param maxit Maximum number of iterations for rank-one offset DMF.
#' @param ... Other parameters?
#' @return DMF structure for x of rank q + 1.
#' @export
dmf_rank1update <- function (dx, x, adjust_svd = FALSE, maxit = 100, refine = TRUE, ...) {
  dx1 <- dmf(x, rank = 1, offset = tcrossprod(dx$L, dx$V), family = dx$family,
             weights = dx$prior.weights,
             control = glm.control(maxit = maxit), ...)
  # [changed to vecorized operation for rank 1 extension]
  # dx1 <- dmf_extend(x, dx$family,
  #                   origin_rank = ncol(dx$L),
  #                   extend_rank = ncol(dx$L) + 1, offset = tcrossprod(dx$L, dx$V),
  #                   weights = dx$prior.weights, control = glm.control(maxit = maxit))

  if (adjust_svd) {
    D <- apply(dx$L, 2, norm2)
    ds <- list(d = D, u = sweep(dx$L, 2, D, `/`), v = dx$V) %>%
      svd_rank1update(dx1$L, dx1$V) # adjust
    dx1$L <- sweep(ds$u, 2, ds$d, `*`)
    dx1$V <- ds$v
  } else {
    dx1$L <- cbind(dx$L, dx1$L)
    dx1$V <- cbind(dx$V, dx1$V)
  }
  if(refine){
    return(dmf(x, rank = ncol(dx1$L), family = dx$family,
               weights = dx$prior.weights, start = dx1)) # refine
  }else{
    return(dx1)
  }

}



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



generate_Y <- function (glm_family, eta_star, phi = 1, glm_weights = 1) {
  y_len <- prod(dim(eta_star))
  mu <- glm_family$linkinv(eta_star)
  family <- glm_family$family
  if (grepl("Negative Binomial", family)) family <- "negative.binomial"
  y <- switch(family,
              gaussian = rnorm(y_len, mu, sqrt(phi)),
              poisson = rpois(y_len, mu),
              binomial = rbinom(y_len, glm_weights, mu),
              #Gamma = rgamma(y_len ,shape = phi, scale = mu / phi),
              Gamma = rgamma(y_len ,shape = 1/phi, scale = mu  * phi), # var(data) = phi * mu^2
              negative.binomial = rnegbin(y_len, mu = mu, theta = phi),
              stop("family `", glm_family$family, "` not recognized"))
  dim(y) <- dim(eta_star)
  y
}

#' Perform a generalized Hosmer–Lemeshow test , (L * V' + V * L') / 2.
#'
#' @param x Input matrix to be factorized.
#' @param lv DMF structure to be tested.
#' @param G number of groups for eta to be partitioned
#' @param weights Entrywise weight.
#' @param offset Entrywise offset.
#' @param dispersion If cutomized dispersion is provided
#' @export
family_test <- function (x, lv, G, weights = 1, offset = zeros(x), dispersion = NULL) {

  n <- nrow(x); p <- ncol(x)
  if (n < p) stop("fewer observations than predictors")
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

  # [ finite mu fit]
  mu_hat <- dmf_fitted(lv)
  is_inf_mu <- is.infinite(mu_hat)

  G <- as.integer(G)
  if (G > n*p) stop('fewer observations than groups of partition')

  #if (endsWith(lv$family$family,'binomial')) x <- x / weights

  if (!is.null(dispersion)){
    phi_hat <- dispersion
  }else{
    if (lv$family$family == "Gamma") {
      phi_hat = 1
      #phi_hat <- 1 / mean((x - mu_hat) ^ 2) # todo: check Gamma dispersion correct?
    } else if (lv$family$family == "gaussian") {
      phi_hat <- mean((x[!is_inf_mu] - mu_hat[!is_inf_mu]) ^ 2)
    } else if (startsWith(lv$family$family, "quasi")) {
      phi_hat <- mean((x[!is_inf_mu] - mu_hat[!is_inf_mu]) ^ 2)
    }
    else {
      phi_hat <- 1
    }
  }
  eta_hat <- data.frame(eta = c(lv$family$linkfun(mu_hat)))
  cut_quantle <- seq(0, 1, 1 / G)
  # TODO: remove magrittr dependency
  #     : warning on num cut less than G
  eta_g <- eta_hat %>% group_by() %>% mutate(G = cut(eta, unique(quantile(eta, cut_quantle)),
                                                    include.lowest = TRUE))
  eta_g$G <- as.numeric(as.factor(eta_g$G))
  num_cut <- length(unique(eta_g$G))
  num_remove <- G - num_cut


  # [ construct statistics ]
  resids <- matrix((x - mu_hat)  * weights) / sqrt(n * p)
  var <- matrix(lv$family$variance(mu_hat)* weights) / (n * p)
  # resids[is_inf_mu] = 0 # setting infinite mu fit to 0
  # inf_group <- eta_g$G[is_inf_mu]


  # [ Vectorized the for loop]
  valid_group <- seq_len(num_cut)
  # valid_group <- valid_group[valid_group!=inf_group]
  s_vec <- sapply(valid_group, function(x) sum(resids[eta_g$G==x]))
  d_diag <- sapply(valid_group, function(x) phi_hat * sum(var[eta_g$G==x]))

  norm_vec = c(s_vec/sqrt(d_diag), rep(num_remove,0))
  norm_vec[is.na(norm_vec)] = 0


  # [Chisq is always returned]
  return(list(norm_vec = norm_vec,
  chisq_stat = crossprod(norm_vec)))
}


#' Conduct factorization rank selection based upon ACT
#'
#' @param x Input matrix to be factorized.
#' @param family Family object to specify deviance loss.
#' @param weights Entrywise weight.
#' @param offset Entrywise offset.
#' @export
act_rank <- function(x, family = gaussian(), weights = 1, offset = zeros(x)) {
  n <- nrow(x); p <- ncol(x); d <- ncol(x)
  if (n < p) stop("fewer observations than predictors")
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

  invlink_Y <- family$linkfun(x / weights)
  inf_ily <- is.infinite(invlink_Y)
  invlink_Y[inf_ily & (invlink_Y > 0)] <- max(invlink_Y[!inf_ily])
  invlink_Y[inf_ily & (invlink_Y < 0)] <- min(invlink_Y[!inf_ily])
  if (any(is.na(invlink_Y))) {
    warning("NA produced in g(Y) transformation, imputing with mean")
    invlink_Y[is.na(invlink_Y)] <- mean(invlink_Y, na.rm =TRUE)
  }

  lambdas <- eigen(cor(invlink_Y))$value
  m <- rep(0, d - 1)
  for (j in 1:(d - 1)) {
    m[j] <- 1 / (d - j) * (sum(1/(lambdas[(j + 1):d] - lambdas[j])) +
                           1 / ((3 * lambdas[j] + lambdas[j + 1]) / 4 - lambdas[j]))
  }

  rho_j <- (d - 1:(d - 1)) / (n - 1)
  m1 <- rho_j * m - (1 - rho_j) / lambdas[1:(d-1)]
  lambdas_adjust <- -1 / m1
  sum(lambdas_adjust > 1 + sqrt(d / n))
}


#' Conduct factorization rank selection based upon maximum eigenvalue gap
#' Onatski(2010). Algorithm first fit dmf with rank = q_max then conduct eigen
#' gap search on the covariance of fitted eta
#' @param x Input matrix to be factorized.
#' @param family Family object to specify deviance loss.
#' @param q_max Maximum rank eta.
#' @param weights Entrywise weight.
#' @param offset Entrywise offset.
#' @param max_iter Maximum iteration.
#' @param fit_full Fit eta with full rank matrix? Always recommended if computation speed is permitted
#' @export
onatski_rank <- function(x, family, q_max, weights = 1, offset = zeros(x),
                            max_iter = 10, fit_full = FALSE) {
  n <- nrow(x); p <- ncol(x)
  if (p < q_max + 5 ) stop('decrease rmax')

  if(fit_full){fit_result <- dmf(x, family, rank= p, weights, offset)
  }else{fit_result <- dmf(x, family, rank= q_max, weights, offset)  }

  eta <- tcrossprod(fit_result$L, fit_result$V)
  cov_matrix <- cov(eta)
  lambda_hats <- eigen(cov_matrix)$value
  tol <- 1e3
  j_update <- q_max + 1
  while (tol>=1){
    j <- j_update
    y_reg <- lambda_hats[j:(j+4)]
    x_reg <- (j + seq(-1, 3, 1))^(2/3)
    delta_temp <- as.numeric(2* abs(coef(lm(y_reg ~ x_reg))))[2]
    flag <- which((-diff(lambda_hats)>= delta_temp))
    if(length(flag) == 0){q_hat = 0}else{q_hat <- tail(flag[flag<=q_max],1)}
    j_update <- q_hat + 1
    tol <- abs(j_update -j)
  }

  return (list(q_hat = q_hat, delta = delta_temp, L =fit_result$L,
               V = fit_result$V, deviance = fit_result$deviance, family = fit_result$family))
}

#' Conduct force-sign SVD given eta and desired rank q
#' @param eta Input matrix to conduct force-sign svd.
#' @param q designed rank of the force-sign svd.
#' @export
dmf_identify <- function(eta, q){

  svd_eta <- svd(eta, nu = q, nv = q)
  negative_flag <- c(1:q)[svd_eta$v[1,]<0]

  svd_eta$u[,negative_flag] <- matrix(-svd_eta$u[,negative_flag])
  svd_eta$v[,negative_flag] <- matrix(-svd_eta$v[,negative_flag])
  # [dmf identify l.v]
  V <- svd_eta$v
  L <- tcrossprod(svd_eta$u, diag(svd_eta$d[1:q]))
  return (list(L = L, V = V))
}


#' extend rank from q to an arbitrary rank `\tilde{q}` through offset = eta_q
#' @param dx dmf rank q fit result
#' @param x Input matrix to be factorized.
#' @param adjust_svd Todo:fixme.
#' @param refine boolen to indicate whether or not to refine rank q + 1 fit `\Lambda` V^T
#' @param control Algorithm control parameters (see \code{glm.control}).
#' @return Extended DMF structure.
#' @export
dmf_rank1update <- function (dx, x, adjust_svd = FALSE, refine = TRUE, control = glm.control(epsilon = 1e-6, maxit = 100), ...) {

  # [can changed to vectorized operation for rank 1 extension]
  dx1 <- dmf(x, rank = 1, offset = tcrossprod(dx$L, dx$V), family = dx$family,
             weights = dx$prior.weights,
             control = control, ...)

  if (adjust_svd) {
    D <- apply(dx$L, 2, norm2)
    ds <- list(d = D, u = sweep(dx$L, 2, D, `/`), v = dx$V) %>%
      svd_rank1update(dx1$L, dx1$V) # adjust
    dx1$L <- sweep(ds$u, 2, ds$d, `*`)
    dx1$V <- ds$v
  } else {
    dx1$L <- cbind(dx$L, dx1$L)
    dx1$V <- cbind(dx$V, dx1$V)
  }
  if(refine){
    return(dmf(x, rank = ncol(dx1$L), family = dx$family,
               weights = dx$prior.weights, start = dx1)) # refine
  }else{
    return(dx1)
  }

}

dmf_extend <- function (x, dx, extend_rank = rank + 1, adjust_svd = FALSE,
                        refine = FALSE, control = glm.control(epsilon = 1e-6, maxit = 100)
                        ){

  # [ some definition inherited from offset dx fit]
  n <- nrow(x); p <- ncol(x); origin_rank <- ncol(dx$L);
  family <- dx$family; weights <- dx$prior.weights
  offset <- tcrossprod(dx$L, dx$V)
  rank <- extend_rank -origin_rank

  # [ no need for checking since everything is inherited from dx]
  valid <- (weights > 0) & (!is.na(x))

  # [ initialize ]
  mu <- family_initialize(x, weights, family)
  eta <- family$linkfun(mu)

  #[transformation from init can gets negative]
  # valid <- valid & family$valideta(eta) ?

  eta[!valid] <- mu[!valid] <- 0
  L <- (eta - offset)[, 1:rank, drop = FALSE] # V = I_p
  V <- matrix(nrow = p, ncol = rank)

  # [ iterate ]
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

      # [efficient vectorization for extend rank 1]
      if (rank == 1){
          # [V step]
          left <- colSums(sweep(sqrt(S), 1, L, "*")^2)
          right <- colSums(sweep(S * Z, 1, L, "*"))
          V <- right/left
          # [L step]
          left <- rowSums(sweep(sqrt(S), 2, V, "*")^2)
          right <- rowSums(sweep(S * Z, 2, V, "*"))
          L <- right/left
      }else{
      # [iteration over n,p for extend rank d]
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
      }
      eta <- tcrossprod(L, V) + offset
      mu <- family$linkinv(eta)
      dr <- family$dev.resids(x[valid], mu[valid], weights[valid])
      dr[is.na(dr) | is.nan(dr)] <- 0
      dev_new <- sum(dr)
      # if (it > 1 && (dev - dev_new) / (dev + .1) < control$epsilon) break
      if (it > 2 && (dev - dev_new) / (dev + .1) < control$epsilon) break
      if (control$trace) message("[", it, "] dev = ", dev_new)
      dev <- dev_new
      # print(dev_new)
      L_old <- L; V_old <- V
  }

  # [identify the results]
  dx1  <- dx; dx1$deviance <- dev; dx1$iter <- it
  if (adjust_svd) {
    # TODO: [projection with dimension d]
    # etaq = dmf_identify(offset, origin_rank)
    # LL_diag = diag(crossprod(etaq$L))
    # L = L_old - etaq$L %*% (1/LL_diag * crossprod(etaq$L, L_old))
    # V = V_old - t(tcrossprod(crossprod(V_old, etaq$V), etaq$V))
    # list(L = etaq_identified$L, V = etaq_identified$V, deviance = dev, family = family)

    # [needs to study Luis's efficient projection]
    # D <- apply(dx$L, 2, norm2)
    # ds <- list(d = D, u = sweep(dx$L, 2, D, `/`), v = dx$V) |>
    #   svd_rank1update(dx1$L, dx1$V) # adjust
    # dx1$L <- sweep(ds$u, 2, ds$d, `*`)
    # dx1$V <- ds$v
  } else {
    L_old <- cbind(dx$L, L_old)
    V_old <- cbind(dx$V, V_old)
    ql <- qr(L_old); qv <- qr(V_old)
    sr <- svd(tcrossprod(qr.R(ql), qr.R(qv)))
    L <- qr.Q(ql) %*% sweep(sr$u, 2, sr$d, `*`)
    V <- qr.Q(qv) %*% sr$v
    dx1$L <- L; dx1$V <- V
  }
  if(refine){
    return(dmf(x, rank = ncol(dx1$L), family = dx$family, weights = dx$prior.weights, start = dx1)) # refine
  }else{return(dx1)}
}

