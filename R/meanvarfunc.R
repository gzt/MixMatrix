

.SStep <- function(data, centers, U, V, weights) {
  dims <- dim(data)
  p <- dims[1]
  q <- dims[2]
  n <- dims[3]


  zigmult = rep(weights, each = p * p)
  swept.data <- sweep(data, c(1, 2), centers)

  Stmp = xatx(swept.data, V)
  for (obs in 1:n) Stmp[, , obs] = Stmp[, , obs] + U
  Smatrix = cubeinv(Stmp)

  SS = rowSums(Smatrix * zigmult, FALSE, 2)

  SSXtmp = cubemult(Smatrix * zigmult, data)
  SSX = rowSums(SSXtmp, FALSE, 2)

  SSXXtmp = cubemult(data, SSXtmp)
  SSXX = rowSums(SSXXtmp, FALSE, 2)
  SSD = detsum(Smatrix, weights)

  list(SS = SS, SSX = SSX, SSXX = SSXX, SSD = SSD)
}

# if 'normal', require only data and weights. if otherwise, require U, V, SS, SSX.
.MeansFunction <- function(data, V = NULL, SS = NULL, SSX = NULL, weights, row.mean = FALSE, col.mean = FALSE, model = "normal", ...) {
  dims <- dim(data)
  p <- dims[1]
  q <- dims[2]
  n <- dims[3]

  sumzig = sum(weights)
  newcenters = matrix(0, nrow = p, ncol = q)
  if (model == "normal") {
    for (obs in 1:n) {
      newcenters = newcenters + data[, , obs] * weights[obs]
    }

    newcenters = newcenters / sumzig
    if (row.mean) {
      # make it so that the mean is constant within a row
      newcenters <- matrix(rowMeans(newcenters), nrow = dims[1], ncol = dims[2])
    }
    if (col.mean) {
      # make it so that the mean is constant within a column
      newcenters <- matrix(colMeans(newcenters), nrow = dims[1], ncol = dims[2], byrow = TRUE)
    }
  } else {
    if (row.mean && col.mean) {
      # make it so that the mean is constant within a row
      scalarmu = matrixtrace(SSX %*% solve(V) %*% ones(q, p)) / matrixtrace(SS %*% ones(p, q) %*% solve(V) %*% ones(q, p))
      newcenters <- scalarmu * ones(p, q)
    } else if (col.mean) {
      # make it so that the mean is constant within a column
      # ie mu = p x 1, times ones 1 x q
      newcenters <- ones(p, p) %*% SSX / sum(SS)
    } else if (row.mean) {
      # make it so that the mean is constant within a row
      # ie  ones p x 1 times mu = 1 x q
      newcenters = solve(SS) %*% SSX %*% (solve(V) %*% ones(q, q)) / sum(solve(V))
    } else {
      newcenters = solve(SS) %*% SSX
    }
  }

  newcenters
}





.colVars <- function(data, center, df = 0, weights, SS, SSX, SSXX,
                     col.variance = "none", col.set.var = FALSE, varflag = FALSE, ...) {
  n = sum(weights)
  p = dim(data)[1]
  q = dim(data)[2]
  dfmult = df + p + q - 1
  if (col.variance == "I") {
    new.V = diag(q)
  } else if (col.set.var) {
    nLL <- function(theta) {
      vardetmat <- vardet(q, theta, TRUE, col.variance)
      varinvmat <- varinv(q, theta, TRUE, col.variance)
      # SXOX = rowSums(axbt(SSXtmp,varinvmat,data ), dims = 2)
      SXOX = SSX %*% varinvmat %*% t(rowSums(data, dims = 2)) # only need trace to be equal
      return(-n * p * vardetmat + dfmult * matrixtrace(SXOX + SS %*% center %*% varinvmat %*% t(center) - SSX %*% varinvmat %*% t(center) - center %*% varinvmat %*% t(SSX)))
    }
    if (!isTRUE(sign(nLL(0.01)) * sign(nLL(.99)) <= 0)) {
      warning("Endpoints of derivative of likelihood do not have opposite sign. Check variance specification.")
      rho.col = 0
      varflag = TRUE
    } else {
      fit0 <- stats::uniroot(nLL, c(0.01, .999), ...)
      rho.col <- fit0$root
    }
    new.V <- varmatgenerate(q, rho.col, col.variance)
  } else {
    new.V = (dfmult / (n * p)) * (SSXX - t(SSX) %*% center - t(center) %*% SSX + t(center) %*% SS %*% center)
    if (col.variance == "cor") {
      new.V = stats::cov2cor(new.V)
      if (!all(is.finite(new.V))) {
        varflag = TRUE
        new.V = diag(q)
      }
    } else new.V = new.V / new.V[1, 1]
  }
  # Fix V to have unit variance on first component
  list(V = new.V, varflag = varflag)
}


.rowVars <- function(data, center, df = 0, weights, SS, SSX, SSXX,
                     row.variance = "none", row.set.var = FALSE, varflag = FALSE, ...) {
  n = sum(weights)
  p = dim(data)[1]
  q = dim(data)[2]
  dfmult = df + p + q - 1

  if (row.variance == "I") {
    new.U = diag(p) * n * (df + p - 1) * p / matrixtrace(SS * dfmult)
  } else if (row.set.var) {
    nLL <- function(theta) {
      vardetmat <- vardet(p, theta, TRUE, row.variance)
      Sigma = varmatgenerate(p, theta, row.variance)
      var = n * (df + p - 1) * p / matrixtrace(Sigma %*% SS * dfmult)
      varderivative = varderiv(p, theta, row.variance)
      return(var * dfmult * matrixtrace(varderivative %*% SS) + n * (df + p - 1) * vardetmat)
    }
    if (!isTRUE(sign(nLL(0.01)) * sign(nLL(.999)) <= 0)) {
      warning("Endpoints of derivative of likelihood do not have opposite sign. Check variance specification.")
      rho.row = 0
      varflag = TRUE
    } else {
      fit0 <- stats::uniroot(nLL, c(0.01, .998), ...)
      rho.row <- fit0$root
    }

    new.U <- varmatgenerate(p, rho.row, row.variance)
    var = n * (df + p - 1) * p / matrixtrace(new.U %*% SS * dfmult)
    new.U = var * new.U
  } else {
    newUinv = (dfmult / (n * (df + p - 1))) * SS
    new.U = solve(newUinv)
    if (row.variance == "cor") {
      vartmp = exp(mean(log(diag(new.U)))) # should be pos def so no problems
      if (!is.finite(vartmp)) {
        vartmp = 1
        varflag = TRUE
        warning("Variance estimate for correlation matrix not positive definite.")
      }
      new.U = vartmp * stats::cov2cor(new.U)
      # this cute trick preserves the determinant of the matrix
    }
  }

  list(U = new.U, varflag = varflag)
}


.varparse <- function(varoption) {
  varflag = FALSE
  varopt = varoption

  if (grepl("^i", x = varoption, ignore.case = TRUE)) {
    varflag = TRUE
    varopt = "I"
  }

  if (grepl("^co", x = varoption, ignore.case = TRUE)) {
    varflag = FALSE
    varopt = "cor"
  }
  if (grepl("^ar", x = varoption, ignore.case = TRUE)) {
    varflag = TRUE
    varopt = "AR(1)"
  }
  if (grepl("^cs", x = varoption, ignore.case = TRUE)) {
    varflag = TRUE
    varopt = "CS"
  }

  list(varflag = varflag, varopt = varopt)
}
