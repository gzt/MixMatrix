library(matrixdist)

context("Testing input type integrity")


test_that("trying wrong type of input", {
  expect_error(rmatrixnorm(n = 1, mean = matrix(c("A", 1))),"numeric")
  expect_error(rmatrixt(n = 1, df = 1, mean = matrix(c("A", 1))),"numeric")
  expect_error(rmatrixinvt(n = 1, df = 1, mean = matrix(c("A", 1))),"numeric")
  expect_error(dmatrixnorm(x = matrix(c("A", 1))),"numeric")
  expect_error(dmatrixt(x = matrix(c("A", 1)), df = 1),"numeric")
  expect_error(dmatrixinvt(x = matrix(c("A", 1)), df = 1),"numeric")

  expect_error(rmatrixnorm(
    n = 1,
    mean = matrix(c(0, 0)),
    U = matrix(c("A", 0, 0, 1), nrow = 2)),
    "non-numeric", ignore.case = TRUE)
  expect_error(rmatrixt(
    n = 1,
    df = 1,
    mean = matrix(c(0, 0)),
    U = matrix(c("A", 0, 0, 1), nrow = 2)),
    "non-numeric", ignore.case = TRUE)
  expect_error(rmatrixinvt(
    n = 1,
    df = 1,
    mean = matrix(c(0, 0)),
    U = matrix(c("A", 0, 0, 1), nrow = 2)),
    "non-numeric", ignore.case = TRUE)
  expect_error(dmatrixnorm(x = matrix(c(0, 0)),
                           U = matrix(c("A", 0, 0, 1), nrow = 2)),
               "non-numeric", ignore.case = TRUE)
  expect_error(dmatrixt(
    x = matrix(c(0, 0)),
    df = 1,
    U = matrix(c("A", 0, 0, 1), nrow = 2)),
    "non-numeric", ignore.case = TRUE)
  expect_error(dmatrixinvt(
    x = mean(matrix(c(0, 0))),
    df = 1,
    U = matrix(c("A", 0, 0, 1), nrow = 2)),
    "non-numeric", ignore.case = TRUE)


  expect_error(rmatrixnorm(
    n = "A",
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    R = matrix(c(5, 1, 2, 0, 6, 1, -1, 2, 10), nrow = 3),
    list = FALSE),
    "non-numeric", ignore.case = TRUE)
  expect_error(rmatrixnorm(
    n = 10,
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = "A"),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    R = matrix(c(5, 1, 2, 0, 6, 1, -1, 2, 10), nrow = 3),
    list = FALSE),
    "non-numeric", ignore.case = TRUE)
  expect_error(rmatrixnorm(
    n = 10,
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c("A", 1, 0, .1), nrow = 2),
    R = matrix(c(5, 1, 2, 0, 6, 1, -1, 2, 10), nrow = 3),
    list = FALSE))
  expect_error(rmatrixnorm(
    n = 10,
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    R = matrix(c("A", 1, 2, 0, 6, 1, -1, 2, 10), nrow = 3),
    list = FALSE))
  expect_error(rmatrixnorm(
    n = 10,
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    R = matrix(c("A", 1, 2, 0, 6, 1, -1, 2, 10), nrow = 3),
    list = TRUE))
  expect_error(rmatrixnorm(
    n = 10,
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    U = matrix(c("A", 1, 0, .1), nrow = 2),
    R = matrix(c(5, 1, 2, 0, 6, 1, -1, 2, 10), nrow = 3),
    list = FALSE),
    "non-numeric", ignore.case = TRUE)

  expect_error(dmatrixnorm(
    matrix(c("A", 0, -100, 0, 25, -1000), nrow = 2),
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    log = TRUE),
    "non-numeric", ignore.case = TRUE)
  expect_error(dmatrixnorm(
    matrix(c(100, 0, -100, 0, 25, -1000), nrow = "A"),
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    log = TRUE),
    "non-numeric", ignore.case = TRUE)
  expect_error(dmatrixnorm(
    matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    mean = matrix(c("A", 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    log = TRUE))
  expect_error(dmatrixnorm(
    matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c("A", 1, 0, .1), nrow = 2),
    log = TRUE))
  expect_error(dmatrixnorm(
    matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    L = matrix(c(2, 1, 0, .1), nrow = "A"),
    log = TRUE))
  expect_error(dmatrixnorm(
    matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow =
                    2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    R = matrix(c("A", 1, 2, 0, 6, 1, -1, 2, 10), nrow =
                 3),
    log = TRUE))
  expect_error(dmatrixnorm(
    matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow =
                    2),
    U = matrix(c("A", 1, 0, .1), nrow = 2),
    log = TRUE),
    "non-numeric", ignore.case = TRUE)
  expect_error(dmatrixnorm(
    matrix(c(100, 0, -100, 0, 25, -1000), nrow = 2),
    mean = matrix(c(100, 0, -100, 0, 25, -1000), nrow =
                    2),
    L = matrix(c(2, 1, 0, .1), nrow = 2),
    V = matrix(c("A", 1, 2, 0, 6, 1, -1, 2, 10), nrow =
                 3),
    log = TRUE),
    "non-numeric", ignore.case = TRUE)


  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean = diag(5)), tol =
                              "Q"), "non-numeric", ignore.case = TRUE)
  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean = diag(5)), max.iter =
                              "Q"), "non-numeric", ignore.case = TRUE)
  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean = diag(5)),
                            U = matrix("Q", nrow = 5, ncol = 5)),
               "non-numeric", ignore.case = TRUE)
  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean = diag(5)),
                            V = matrix("Q", nrow = 5, ncol = 5)),
               "non-numeric", ignore.case = TRUE)

  expect_error(MLmatrixt(rmatrixnorm(n = 100, mean = diag(5)), tol =
                              "Q"), "non-numeric", ignore.case = TRUE)
  expect_error(MLmatrixt(rmatrixnorm(n = 100, mean = diag(5)), max.iter =
                              "Q"), "non-numeric", ignore.case = TRUE)
  expect_error(MLmatrixt(rmatrixnorm(n = 100, mean = diag(5)), df =
                           "Q"), "non-numeric", ignore.case = TRUE)
  expect_error(MLmatrixt(rmatrixnorm(n = 100, mean = diag(5)),
                            U = matrix("Q", nrow = 5, ncol = 5)),
               "non-numeric", ignore.case = TRUE)
  expect_error(MLmatrixt(rmatrixnorm(n = 100, mean = diag(5)),
                            V = matrix("Q", nrow = 5, ncol = 5)),
               "non-numeric", ignore.case = TRUE)

  # expect_error(toepgenerate(7, "A"))
  # expect_error(toepgenerate("A", .4))



})

context("Testing input dimension integrity")



test_that("Testing bad matrix dimension input", {
  A <- diag(3)
  B <- diag(4)

  expect_error(rmatrixnorm(
    n = 1,
    mean = A,
    L = B,
    R = A,
    list = FALSE),
    "Non-conforming")
  expect_error(rmatrixnorm(
    n = 1,
    mean = A,
    L = A,
    R = B,
    list = FALSE),
    "Non-conforming")
  expect_error(rmatrixnorm(
    n = 1,
    mean = A,
    U = B,
    R = A,
    list = FALSE),
    "Non-conforming")
  expect_error(rmatrixnorm(
    n = 1,
    mean = A,
    L = A,
    V = B,
    list = FALSE),
    "Non-conforming")
  expect_error(rmatrixnorm(
    n = 1,
    mean = A,
    U = A,
    V = B,
    list = FALSE),
    "Non-conforming")
  expect_error(rmatrixnorm(
    n = 1,
    mean = A,
    U = A,
    R = B,
    list = FALSE),
    "Non-conforming")

  expect_error(rmatrixt(
    n = 1,
    df = 1,
    mean = A,
    L = B,
    R = A,
    list = FALSE),
    "Non-conforming")
  expect_error(rmatrixt(
    n = 1,
    df = 1,
    mean = A,
    L = A,
    R = B,
    list = FALSE),
    "Non-conforming")
  expect_error(rmatrixt(
    n = 1,
    df = 1,
    mean = A,
    U = B,
    R = A,
    list = FALSE),
    "Non-conforming")
  expect_error(rmatrixt(
    n = 1,
    df = 1,
    mean = A,
    L = A,
    V = B,
    list = FALSE),
    "Non-conforming")
  expect_error(rmatrixt(
    n = 1,
    df = 1,
    mean = A,
    U = A,
    V = B,
    list = FALSE),
    "Non-conforming")
  expect_error(rmatrixt(
    n = 1,
    df = 1,
    mean = A,
    U = A,
    R = B,
    list = FALSE),
    "Non-conforming")

  expect_error(rmatrixinvt(
    n = 1,
    df = 1,
    mean = A,
    L = B,
    R = A,
    list = FALSE),
    "Non-conforming")

  expect_error(rmatrixinvt(
    n = 1,
    df = 1,
    mean = A,
    L = A,
    R = B,
    list = FALSE),
    "Non-conforming")
  expect_error(rmatrixinvt(
    n = 1,
    df = 1,
    mean = A,
    U = B,
    R = A,
    list = FALSE
  ),
  "Non-conforming")
  expect_error(rmatrixinvt(
    n = 1,
    df = 1,
    mean = A,
    L = A,
    V = B,
    list = FALSE
  ),
  "Non-conforming")
  expect_error(rmatrixinvt(
    n = 1,
    df = 1,
    mean = A,
    U = A,
    V = B,
    list = FALSE
  ),
  "Non-conforming")
  expect_error(rmatrixinvt(
    n = 1,
    df = 1,
    mean = A,
    U = A,
    R = B,
    list = FALSE
  ),
  "Non-conforming")

  expect_error(dmatrixnorm(
    x = A,
    mean = A,
    L = B,
    R = A
  ),
  "Non-conforming")
  expect_error(dmatrixnorm(
    x = A,
    mean = A,
    L = A,
    R = B
  ),
  "Non-conforming")
  expect_error(dmatrixnorm(
    x = A,
    mean = A,
    U = B,
    R = A
  ),
  "Non-conforming")
  expect_error(dmatrixnorm(
    x = A,
    mean = A,
    L = A,
    V = B
  ), "Non-conforming")
  expect_error(dmatrixnorm(
    x = A,
    mean = A,
    U = A,
    V = B
  ), "Non-conforming")
  expect_error(dmatrixnorm(
    x = A,
    mean = A,
    U = A,
    R = B
  ), "Non-conforming")

  expect_error(dmatrixt(
    x = A,
    df = 1,
    mean = A,
    L = B,
    R = A
  ), "Non-conforming")
  expect_error(dmatrixt(
    x = A,
    df = 1,
    mean = A,
    L = A,
    R = B
  ), "Non-conforming")
  expect_error(dmatrixt(
    x = A,
    df = 1,
    mean = A,
    U = B,
    R = A
  ), "Non-conforming")
  expect_error(dmatrixt(
    x = A,
    df = 1,
    mean = A,
    L = A,
    V = B
  ), "Non-conforming")
  expect_error(dmatrixt(
    x = A,
    df = 1,
    mean = A,
    U = A,
    V = B
  ), "Non-conforming")
  expect_error(dmatrixt(
    x = A,
    df = 1,
    mean = A,
    U = A,
    R = B
  ), "Non-conforming")

  expect_error(dmatrixinvt(
    x = A,
    df = 1,
    mean = A,
    L = B,
    R = A
  ), "Non-conforming")
  expect_error(dmatrixinvt(
    x = A,
    df = 1,
    mean = A,
    L = A,
    R = B
  ), "Non-conforming")
  expect_error(dmatrixinvt(
    x = A,
    df = 1,
    mean = A,
    U = B,
    R = A
  ), "Non-conforming")
  expect_error(dmatrixinvt(
    x = A,
    df = 1,
    mean = A,
    L = A,
    V = B
  ), "Non-conforming")
  expect_error(dmatrixinvt(
    x = A,
    df = 1,
    mean = A,
    U = A,
    V = B
  ), "Non-conforming")
  expect_error(dmatrixinvt(
    x = A,
    df = 1,
    mean = A,
    U = A,
    R = B
  ), "Non-conforming")

  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean = A), U = B, V = A))
  expect_error(MLmatrixnorm(rmatrixnorm(n = 100, mean = A), U = A, V = B))
  expect_error(MLmatrixt(rmatrixnorm(n = 100, mean = A), U = B, V = A))
  expect_error(MLmatrixt(rmatrixnorm(n = 100, mean = A), U = A, V = B))

})

test_that("Bad input to generators", {

  n = 0
  rho = .5
  expect_error(ARgenerate(n,rho), "greater")
  expect_error(CSgenerate(n,rho), "greater")
  rho = 1.2
  n = 2
  expect_error(ARgenerate(n,rho), "less")
  expect_error(CSgenerate(n,rho), "less")
  rho = -.5
  expect_warning(ARgenerate(n,rho), "greater")
  expect_warning(CSgenerate(n,rho), "greater")
  rho = .999
  expect_warning(ARgenerate(n,rho), "correlation")
  expect_warning(CSgenerate(n,rho), "correlation")
  rho = -2.5
  expect_error(ARgenerate(n,rho), "greater")
  expect_error(CSgenerate(n,rho), "greater")

})

test_that("Testing bad input to LDA/QDA", {
  A <- rmatrixnorm(5, mean = matrix(0, nrow = 2, ncol = 2))
  B <- rmatrixnorm(5, mean = matrix(1, nrow = 2, ncol = 2))
  C <- array(c(A,B), dim = c(2,2,10))
  D <- array(0, dim = c(2,2,5))
  E <- array(c(A,D), dim = c(2,2,10))

  groups <- c(rep(1,5),rep(2,5))
  groups.empty <- factor(rep("1",10), levels = c("1","2"))
  priors = c(.5,.5)
  expect_error(
    matrixlda(c(C), grouping = c(rep(1,40),rep(2,40)), prior = priors),
    "array"
  )
  expect_error(
    matrixlda(C, grouping = c(rep(1,120), rep(2,120)), prior = priors),
    "are different"
  )
  expect_error(
    matrixlda(C, grouping = groups, prior = c(.5,.4)),
    "invalid 'prior'"
  )
  expect_error(
    matrixlda(C, grouping = groups, prior = c(.4,.4,.2)),
    "incorrect length"

  )
  expect_warning(
    matrixlda(C, grouping = groups.empty, prior = priors),
    "empty"
  )
  expect_error(
    matrixlda(E, grouping = groups, prior = priors)
  )


  expect_error(
    matrixqda(c(C), grouping = c(rep(1,40),rep(2,40)), prior = priors),
    "array"
  )
  expect_error(
    matrixqda(C, grouping = c(rep(1,120), rep(2,120)), prior = priors),
    "are different"
  )
  expect_error(
    matrixqda(C, grouping = groups, prior = c(.5,.4)),
    "invalid 'prior'"
  )
  expect_error(
    matrixqda(C, grouping = groups, prior = c(.4,.4,.2)),
    "incorrect length"

  )
  expect_warning(
    matrixqda(C, grouping = groups.empty, prior = priors),
    "empty"
  )
  expect_error(
    matrixqda(E, grouping = groups, prior = priors)
  )

})


test_that("Out of bounds numeric input: ", {
  A <- diag(5)
  A[5, 5] <- 0
  expect_error(rmatrixt(
    0,
    1,
    matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = diag(5)
  ))
  expect_error(rmatrixt(
    1,
    -1,
    matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = diag(5)
  ))
  expect_error(rmatrixt(
    1,
    1,
    matrix(0, nrow = 5, ncol = 5),
    U = -diag(5),
    V = diag(5)
  ))
  expect_error(rmatrixt(
    1,
    1,
    matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = -diag(5)
  ))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), R = A))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5), U = A, V = diag(5)))
  expect_error(rmatrixt(1, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = A))

  expect_error(rmatrixinvt(
    0,
    1,
    matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = diag(5)
  ))
  expect_error(rmatrixinvt(
    1,
    -1,
    matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = diag(5)
  ))
  expect_error(rmatrixinvt(
    1,
    1,
    matrix(0, nrow = 5, ncol = 5),
    U = -diag(5),
    V = diag(5)
  ))
  expect_error(rmatrixinvt(
    1,
    1,
    matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = -diag(5)
  ))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), R = A))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5), U = A, V = diag(5)))
  expect_error(rmatrixinvt(1, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = A))

  expect_error(rmatrixnorm(
    0,
    matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = diag(5)
  ))
  expect_error(rmatrixnorm(
    1,
    matrix(0, nrow = 5, ncol = 5),
    U = -diag(5),
    V = diag(5)
  ))
  expect_error(rmatrixnorm(
    1,
    matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = -diag(5)
  ))
  expect_error(rmatrixnorm(1, matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(rmatrixnorm(1, matrix(0, nrow = 5, ncol = 5), U = diag(5), R = A))
  expect_error(rmatrixnorm(1, 1, matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(rmatrixnorm(1, 1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = A))

  expect_error(dmatrixt(
    1,
    matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = diag(5)
  ))

  expect_error(dmatrixinvt(
    1,
    matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = diag(5)
  ))

  expect_error(dmatrixt(
    df = -1,
    x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = diag(5)
  ))
  expect_error(dmatrixt(
    df = 1,
    x = matrix(0, nrow = 5, ncol = 5),
    U = -diag(5),
    V = diag(5)
  ))
  expect_error(dmatrixt(
    df = 1,
    x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = -diag(5)
  ))
  expect_error(dmatrixt(
    df = 1,
    x = matrix(0, nrow = 5, ncol = 5),
    L = A,
    V = diag(5)
  ))
  expect_error(dmatrixt(
    df = 1,
    x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    R = A
  ))
  expect_error(dmatrixt(
    df = 1,
    x = matrix(0, nrow = 5, ncol = 5),
    U = A,
    V = diag(5)
  ))
  expect_error(dmatrixt(
    df = 1,
    x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = A
  ))

  expect_error(dmatrixinvt(
    df = -1,
    x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = diag(5)
  ))
  expect_error(dmatrixinvt(
    df = 1,
    x = matrix(0, nrow = 5, ncol = 5),
    U = -diag(5),
    V = diag(5)
  ))
  expect_error(dmatrixinvt(
    df = 1,
    x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = -diag(5)
  ))
  expect_error(dmatrixinvt(
    df = 1,
    x = matrix(0, nrow = 5, ncol = 5),
    L = A,
    V = diag(5)
  ))
  expect_error(dmatrixinvt(
    df = 1,
    x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    R = A
  ))
  expect_error(dmatrixinvt(
    df = 1,
    x = matrix(0, nrow = 5, ncol = 5),
    U = A,
    V = diag(5)
  ))
  expect_error(dmatrixinvt(
    df = 1,
    x = matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = A
  ))

  expect_error(dmatrixnorm(
    matrix(0, nrow = 5, ncol = 5),
    U = -diag(5),
    V = diag(5)
  ))
  expect_error(dmatrixnorm(
    matrix(0, nrow = 5, ncol = 5),
    U = diag(5),
    V = -diag(5)
  ))
  expect_error(dmatrixnorm(matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(dmatrixnorm(matrix(0, nrow = 5, ncol = 5), U = diag(5), R = A))
  expect_error(dmatrixnorm(1, matrix(0, nrow = 5, ncol = 5), L = A, V = diag(5)))
  expect_error(dmatrixnorm(1, matrix(0, nrow = 5, ncol = 5), U = diag(5), V = A))



})

context("Testing input positive definite/invertible integrity")


test_that("Bad rank in covariance:", {
  A = matrix(0, nrow = 2, ncol = 2)

  expect_error(rmatrixnorm(
    2,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(rmatrixnorm(
    2,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_length(rmatrixnorm(
    2,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2),
    force = TRUE
  ), 8)
  expect_length(rmatrixnorm(
    2,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2),
    force = TRUE
  ), 8)
  expect_error(rmatrixnorm(
    2,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    U = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(rmatrixnorm(
    2,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    V = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))

  expect_error(rmatrixt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(rmatrixt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_length(rmatrixt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2),
    force = TRUE
  ), 8)
  expect_error(rmatrixt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    U = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(rmatrixt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    V = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))


  expect_error(rmatrixinvt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol = 2)
  ))
  expect_error(rmatrixinvt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol = 2)
  ))
  expect_error(rmatrixinvt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    U = matrix(c(1, 1, .5, .5), nrow = 2, ncol = 2)
  ))
  expect_error(rmatrixinvt(
    2,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    V = matrix(c(1, 1, .5, .5), nrow = 2, ncol = 2)
  ))

  expect_error(dmatrixnorm(
    A,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixnorm(
    A,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixnorm(
    A,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    U = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixnorm(
    A,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    V = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))

  expect_error(dmatrixinvt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixinvt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixinvt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    U = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixinvt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    V = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))

  expect_error(dmatrixt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    L = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    R = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    U = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))
  expect_error(dmatrixt(
    A,
    5,
    mean = matrix(c(0), nrow = 2, ncol = 2),
    V = matrix(c(1, 1, .5, .5), nrow = 2, ncol =
                 2)
  ))



})
