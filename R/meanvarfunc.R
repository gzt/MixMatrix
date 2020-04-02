

.sstep <- function(data, centers, U, V, weights) {
  dims <- dim(data)
  p <- dims[1]
  # q <- dims[2]
  n <- dims[3]


  zigmult <- rep(weights, each = p * p)
  swept_data <- sweep(data, c(1, 2), centers)

  Stmp <- xatx(swept_data, V)
  for (obs in 1:n) Stmp[, , obs] <- Stmp[, , obs] + U
  Smatrix <- cubeinv(Stmp)

  ss <- rowSums(Smatrix * zigmult, FALSE, 2)

  ssxtmp <- cubemult(Smatrix * zigmult, data)
  ssx <- rowSums(ssxtmp, FALSE, 2)

  ssxxtmp <- cubemult(data, ssxtmp)
  ssxx <- rowSums(ssxxtmp, FALSE, 2)
  ssd <- detsum(Smatrix, weights)

  list(ss = ss, ssx = ssx, ssxx = ssxx, ssd = ssd)
}


.means_function <- function(data, V = NULL, ss = NULL, ssx = NULL, weights,
                           row.mean = FALSE, col.mean = FALSE,
                           model = "normal", ...) {
  dims <- dim(data)
  p <- dims[1]
  q <- dims[2]
  n <- dims[3]

  sumzig <- sum(weights)
  newcenters <- matrix(0, nrow = p, ncol = q)
  if (model == "normal") {
    for (obs in 1:n) {
      newcenters <- newcenters + data[, , obs] * weights[obs]
    }

    newcenters <- newcenters / sumzig
    if (row.mean) {
      # make it so that the mean is constant within a row
      newcenters <- matrix(rowMeans(newcenters),
        nrow = dims[1], ncol = dims[2]
      )
    }
    if (col.mean) {
      # make it so that the mean is constant within a column
      newcenters <- matrix(colMeans(newcenters),
        nrow = dims[1],
        ncol = dims[2], byrow = TRUE
      )
    }
  } else {
    if (row.mean && col.mean) {
      # make it so that the mean is constant within a row
      scalarmu <- matrixtrace(ssx %*% solve(V) %*% ones(q, p)) /
        matrixtrace(ss %*% ones(p, q) %*% solve(V) %*% ones(q, p))
      newcenters <- scalarmu * ones(p, q)
    } else if (col.mean) {
      # make it so that the mean is constant within a column
      # ie mu = p x 1, times ones 1 x q
      newcenters <- ones(p, p) %*% ssx / sum(ss)
    } else if (row.mean) {
      # make it so that the mean is constant within a row
      # ie  ones p x 1 times mu = 1 x q
      newcenters <- solve(ss) %*% ssx %*%
        (solve(V) %*% ones(q, q)) / sum(solve(V))
    } else {
      newcenters <- solve(ss) %*% ssx
    }
  }

  newcenters
}





.col_vars <- function(data, center, df = 0, weights, ss, ssx, ssxx,
                     col.variance = "none", col_set_var = FALSE,
                     varflag = FALSE, ...) {
  n <- sum(weights)
  p <- dim(data)[1]
  q <- dim(data)[2]
  dfmult <- df + p + q - 1
  if (col.variance == "I") {
    new_v <- diag(q)
  } else if (col_set_var) {
    nLL <- function(theta) {
      vardetmat <- vardet(q, theta, TRUE, col.variance)
      varinvmat <- varinv(q, theta, TRUE, col.variance)
      # SXOX = rowSums(axbt(ssxtmp,varinvmat,data ), dims = 2)
      SXOX <- ssx %*% varinvmat %*% t(rowSums(data, dims = 2))
      return(-n * p * vardetmat +
        dfmult * matrixtrace(SXOX +
          ss %*% center %*% varinvmat %*% t(center) -
          ssx %*% varinvmat %*% t(center) -
          center %*% varinvmat %*% t(ssx)))
    }
    if (!isTRUE(sign(nLL(0.01)) * sign(nLL(.99)) <= 0)) {
      warning("Endpoints of derivative of likelihood do not have opposite
                 sign. Check variance specification.")
      rho_col <- 0
      varflag <- TRUE
    } else {
      fit0 <- stats::uniroot(nLL, c(0.01, .999), ...)
      rho_col <- fit0$root
    }
    new_v <- varmatgenerate(q, rho_col, col.variance)
  } else {
    new_v <- (dfmult / (n * p)) * (ssxx - t(ssx) %*% center -
      t(center) %*% ssx +
      t(center) %*% ss %*% center)
    if (col.variance == "cor") {
      new_v <- stats::cov2cor(new_v)
      if (!all(is.finite(new_v))) {
        varflag <- TRUE
        new_v <- diag(q)
      }
    } else {
      new_v <- new_v / new_v[1, 1]
    }
  }
  ## Fix V to have unit variance on first component
  list(V = new_v, varflag = varflag)
}


.row_vars <- function(data, center, df = 0, weights, ss, ssx, ssxx,
                     row.variance = "none", row_set_var = FALSE,
                     varflag = FALSE, ...) {
  n <- sum(weights)
  p <- dim(data)[1]
  q <- dim(data)[2]
  dfmult <- df + p + q - 1

  if (row.variance == "I") {
    new_u <- diag(p) * n * (df + p - 1) * p / matrixtrace(ss * dfmult)
  } else if (row_set_var) {
    nLL <- function(theta) {
      vardetmat <- vardet(p, theta, TRUE, row.variance)
      Sigma <- varmatgenerate(p, theta, row.variance)
      var <- n * (df + p - 1) * p / matrixtrace(Sigma %*% ss * dfmult)
      varderivative <- varderiv(p, theta, row.variance)
      return(var * dfmult * matrixtrace(varderivative %*% ss) +
        n * (df + p - 1) * vardetmat)
    }
    if (!isTRUE(sign(nLL(0.01)) * sign(nLL(.999)) <= 0)) {
      warning("Endpoints of derivative of likelihood do not have opposite sign.
               Check variance specification.")
      rho_row <- 0
      varflag <- TRUE
    } else {
      fit0 <- stats::uniroot(nLL, c(0.01, .998), ...)
      rho_row <- fit0$root
    }

    new_u <- varmatgenerate(p, rho_row, row.variance)
    var <- n * (df + p - 1) * p / matrixtrace(new_u %*% ss * dfmult)
    new_u <- var * new_u
  } else {
    new_uinv <- (dfmult / (n * (df + p - 1))) * ss
    new_u <- solve(new_uinv)
    if (row.variance == "cor") {
      vartmp <- exp(mean(log(diag(new_u)))) # should be pos def so no problems
      if (!is.finite(vartmp)) {
        vartmp <- 1
        varflag <- TRUE
        warning("Variance estimate for correlation matrix not
                 positive definite.")
      }
      new_u <- vartmp * stats::cov2cor(new_u)
      # this cute trick preserves the determinant of the matrix
    }
  }

  list(U = new_u, varflag = varflag)
}


.varparse <- function(varoption) {
  varflag <- FALSE
  varopt <- varoption

  if (grepl("^i", x = varoption, ignore.case = TRUE)) {
    varflag <- TRUE
    varopt <- "I"
  }

  if (grepl("^co", x = varoption, ignore.case = TRUE)) {
    varflag <- FALSE
    varopt <- "cor"
  }
  if (grepl("^ar", x = varoption, ignore.case = TRUE)) {
    varflag <- TRUE
    varopt <- "AR(1)"
  }
  if (grepl("^cs", x = varoption, ignore.case = TRUE)) {
    varflag <- TRUE
    varopt <- "CS"
  }

  list(varflag = varflag, varopt = varopt)
}
