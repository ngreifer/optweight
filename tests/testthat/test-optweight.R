skip_on_cran()

test_that("optweight works, binary", {
  skip_if_not_installed("cobalt")

  rlang::local_options(optweight_solver_l2 = "osqp",
                       optweight_solver_l1 = "osqp",
                       optweight_solver_linf = "osqp",
                       optweight_solver_entropy = "scs",
                       optweight_solver_log = "scs")

  eps <- if (capabilities("long.double")) 5e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  bw <- runif(nrow(test_data))

  expect_no_condition({
    ow0 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow0$covs, ow0$treat, ow0$weights)) <= eps))

  expect_true(all(ow0$weights >= 1e-8))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     estimand = "ATE",
                     tols = 0,
                     norm = "l2",
                     min.w = 1e-8)
  })

  expect_equal(ow$weights,
               ow0$weights,
               tolerance = eps)

  # tols
  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                    tols = .02,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= .02 + eps))

  tols <- process_tols(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       tols = seq(.01, .06, by = .01))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    tols = tols,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= attr(tols, "internal.tols") + eps))

  expect_lt(rms_dev(ow$weights),
            rms_dev(ow0$weights))

  #norms
  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    norm = "linf",
                    eps = 1e-4)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_lt(rms_dev(ow0$weights),
            rms_dev(ow$weights))

  expect_gt(max_abs_dev(ow0$weights),
            max_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    norm = "linf",
                    tols = tols,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

  expect_gt(max_abs_dev(ow$weights),
            max_abs_dev(ow2$weights))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    norm = "l1")
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_lt(rms_dev(ow0$weights),
            rms_dev(ow$weights))

  expect_gt(mean_abs_dev(ow0$weights),
            mean_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     norm = "l1",
                     tols = tols,
                     std.binary = TRUE,
                     eps = 1e-5)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

  expect_gt(mean_abs_dev(ow$weights),
            mean_abs_dev(ow2$weights))

  if (rlang::is_installed("scs")) {
    expect_no_condition({
      ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      norm = "entropy")
    })

    expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_lt(rms_dev(ow0$weights),
              rms_dev(ow$weights))

    expect_gt(rel_ent(ow0$weights),
              rel_ent(ow$weights))

    expect_no_condition({
      ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       norm = "entropy",
                       tols = tols,
                       std.binary = TRUE)
    })

    expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_gt(rel_ent(ow$weights),
              rel_ent(ow2$weights))

    expect_no_condition({
      ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      norm = "log")
    })

    expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_lt(rms_dev(ow0$weights),
              rms_dev(ow$weights))

    expect_gt(sum(-log(ow0$weights)),
              sum(-log(ow$weights)))

    expect_no_condition({
      ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       norm = "log",
                       tols = tols,
                       std.binary = TRUE)
    })

    expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_gt(sum(-log(ow$weights)),
              sum(-log(ow2$weights)))
  }
})

test_that("optweight works, binary, s.weights", {
  skip_if_not_installed("cobalt")

  rlang::local_options(optweight_solver_l2 = "osqp",
                       optweight_solver_l1 = "osqp",
                       optweight_solver_linf = "osqp",
                       optweight_solver_entropy = "scs",
                       optweight_solver_log = "scs")

  eps <- if (capabilities("long.double")) 5e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  bw <- runif(nrow(test_data))

  expect_no_condition({
    ow0 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data, s.weights = "SW")
  })

  sw <- test_data$SW

  expect_true(all(abs(cobalt::col_w_smd(ow0$covs, ow0$treat, ow0$weights,
                                        s.weights = sw)) <= eps))

  expect_true(all(ow0$weights >= 1e-8))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    estimand = "ATE",
                    s.weights = "SW",
                    tols = 0,
                    norm = "l2",
                    min.w = 1e-8)
  })

  expect_equal(ow$weights,
               ow0$weights,
               tolerance = eps)

  # tols
  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    tols = .02,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights,
                                        s.weights = sw)) <= .02 + eps))

  tols <- process_tols(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       tols = seq(.01, .06, by = .01))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    tols = tols,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights,
                                        s.weights = sw)) <= attr(tols, "internal.tols") + eps))

  expect_lt(rms_dev(ow$weights, sw = sw),
            rms_dev(ow0$weights, sw = sw))

  #norms
  expect_error({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    norm = "linf")
  })

  expect_error({
    ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     s.weights = "SW",
                     norm = "linf",
                     tols = tols,
                     std.binary = TRUE)
  })


  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    norm = "l1",
                    eps = 1e-5)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights,
                                        s.weights = sw)) <= eps))

  expect_lt(rms_dev(ow0$weights, sw = sw),
            rms_dev(ow$weights, sw = sw))

  expect_gt(mean_abs_dev(ow0$weights, sw = sw),
            mean_abs_dev(ow$weights, sw = sw))

  expect_no_condition({
    ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     s.weights = "SW",
                     norm = "l1",
                     tols = tols,
                     std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights,
                                        s.weights = sw)) <= attr(tols, "internal.tols") + eps))

  expect_gt(mean_abs_dev(ow$weights, sw = sw),
            mean_abs_dev(ow2$weights, sw = sw))

  if (rlang::is_installed("scs")) {
    expect_no_condition({
      ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      s.weights = "SW",
                      norm = "entropy")
    })

    expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights,
                                          s.weights = sw)) <= eps))

    expect_lt(rms_dev(ow0$weights, sw = sw),
              rms_dev(ow$weights, sw = sw))

    expect_gt(rel_ent(ow0$weights, sw = sw),
              rel_ent(ow$weights, sw = sw))

    expect_no_condition({
      ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       s.weights = "SW",
                       norm = "entropy",
                       tols = tols,
                       std.binary = TRUE)
    })

    expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights,
                                          s.weights = sw)) <= attr(tols, "internal.tols") + eps))

    expect_gt(rel_ent(ow$weights, sw = sw),
              rel_ent(ow2$weights, sw = sw))

    expect_no_condition({
      ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      s.weights = "SW",
                      norm = "log")
    })

    expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights,
                                          s.weights = sw)) <= eps))

    expect_lt(rms_dev(ow0$weights, sw = sw),
              rms_dev(ow$weights, sw = sw))

    expect_gt(sum(-sw * log(ow0$weights)),
              sum(-sw * log(ow$weights)))

    expect_no_condition({
      ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       s.weights = "SW",
                       norm = "log",
                       tols = tols,
                       std.binary = TRUE)
    })

    expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights,
                                          s.weights = sw)) <= attr(tols, "internal.tols") + eps))

    expect_gt(sum(-sw * log(ow$weights)),
              sum(-sw * log(ow2$weights)))
  }
})

test_that("optweight works, binary, b.weights", {
  skip_if_not_installed("cobalt")

  rlang::local_options(optweight_solver_l2 = "osqp",
                       optweight_solver_l1 = "osqp",
                       optweight_solver_linf = "osqp",
                       optweight_solver_entropy = "scs",
                       optweight_solver_log = "scs")

  eps <- if (capabilities("long.double")) 5e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  bw <- runif(nrow(test_data))

  expect_no_condition({
    ow0 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow0$covs, ow0$treat, ow0$weights)) <= eps))

  expect_true(all(ow0$weights >= 1e-8))

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow0$weights))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    estimand = "ATE",
                    b.weights = bw,
                    tols = 0,
                    norm = "l2",
                    min.w = 1e-8)
  })

  expect_equal(ow$weights,
               ow0$weights,
               tolerance = eps)

  # tols
  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    tols = .02,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= .02 + eps))

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow$weights))

  tols <- process_tols(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       tols = seq(.01, .06, by = .01))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    tols = tols,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= attr(tols, "internal.tols") + eps))

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow0$weights, bw))

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow$weights))

  #norms
  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    norm = "linf",
                    eps = 2e-5)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow$weights, bw))

  expect_gt(max_abs_dev(ow0$weights, bw),
            max_abs_dev(ow$weights, bw))

  expect_lt(max_abs_dev(ow$weights, bw),
            max_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw,
                     norm = "linf",
                     tols = tols,
                     std.binary = TRUE,
                     eps = 1e-5)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

  expect_gt(max_abs_dev(ow$weights, bw),
            max_abs_dev(ow2$weights, bw))

  expect_lt(max_abs_dev(ow2$weights, bw),
            max_abs_dev(ow2$weights))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    norm = "l1",
                    eps = 1e-5)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow$weights, bw))

  expect_gt(mean_abs_dev(ow0$weights, bw),
            mean_abs_dev(ow$weights, bw))

  expect_lt(mean_abs_dev(ow$weights, bw),
            mean_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw,
                     norm = "l1",
                     tols = tols,
                     std.binary = TRUE,
                     eps = 5e-5)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

  expect_gt(mean_abs_dev(ow$weights, bw),
            mean_abs_dev(ow2$weights, bw))

  expect_lt(mean_abs_dev(ow2$weights, bw),
            mean_abs_dev(ow2$weights))

  if (rlang::is_installed("scs")) {
    expect_no_condition({
      ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      b.weights = bw,
                      norm = "entropy")
    })

    expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_lt(rms_dev(ow0$weights, bw),
              rms_dev(ow$weights, bw))

    expect_gt(rel_ent(ow0$weights, bw),
              rel_ent(ow$weights, bw))

    expect_lt(rel_ent(ow$weights, bw),
              rel_ent(ow$weights))

    expect_no_condition({
      ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       b.weights = bw,
                       norm = "entropy",
                       tols = tols,
                       std.binary = TRUE)
    })

    expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_gt(rel_ent(ow$weights, bw),
              rel_ent(ow2$weights, bw))

    expect_lt(rel_ent(ow2$weights, bw),
              rel_ent(ow2$weights))

    expect_no_condition({
      ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      b.weights = bw,
                      norm = "log")
    })

    expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_lt(rms_dev(ow0$weights, bw),
              rms_dev(ow$weights, bw))

    expect_gt(sum(-log(ow0$weights / bw)),
              sum(-log(ow$weights / bw)))

    expect_no_condition({
      ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       b.weights = bw,
                       norm = "log",
                       tols = tols,
                       std.binary = TRUE)
    })

    expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_gt(sum(-log(ow$weights / bw)),
              sum(-log(ow2$weights / bw)))
  }
})

test_that("optweight works, binary, s.weights + b.weights", {
  skip_if_not_installed("cobalt")

  rlang::local_options(optweight_solver_l2 = "osqp",
                       optweight_solver_l1 = "osqp",
                       optweight_solver_linf = "osqp",
                       optweight_solver_entropy = "scs",
                       optweight_solver_log = "scs")

  eps <- if (capabilities("long.double")) 5e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  bw <- runif(nrow(test_data))

  expect_no_condition({
    ow0 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data, s.weights = "SW",
                     b.weights = bw)
  })

  sw <- test_data$SW

  expect_true(all(abs(cobalt::col_w_smd(ow0$covs, ow0$treat, ow0$weights,
                                        s.weights = sw)) <= eps))

  expect_true(all(ow0$weights >= 1e-8))

  expect_lt(rms_dev(ow0$weights, bw, sw = sw),
            rms_dev(ow0$weights, sw = sw))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    estimand = "ATE",
                    s.weights = "SW",
                    b.weights = bw,
                    tols = 0,
                    norm = "l2",
                    min.w = 1e-8)
  })

  expect_equal(ow$weights,
               ow0$weights,
               tolerance = eps)

  # tols
  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    b.weights = bw,
                    tols = .02,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights,
                                        s.weights = sw)) <= .02 + eps))

  expect_lt(rms_dev(ow$weights, bw, sw = sw),
            rms_dev(ow$weights, sw = sw))

  tols <- process_tols(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       tols = seq(.01, .06, by = .01))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    b.weights = bw,
                    tols = tols,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights,
                                        s.weights = sw)) <= attr(tols, "internal.tols") + eps))

  expect_lt(rms_dev(ow$weights, bw, sw = sw),
            rms_dev(ow0$weights, bw, sw = sw))

  expect_lt(rms_dev(ow$weights, bw, sw = sw),
            rms_dev(ow$weights, sw = sw))

  #norms
  expect_error({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    b.weights = bw,
                    norm = "linf")
  })

  expect_error({
    ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     s.weights = "SW",
                     b.weights = bw,
                     norm = "linf",
                     tols = tols,
                     std.binary = TRUE)
  })


  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    b.weights = bw,
                    norm = "l1")
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights,
                                        s.weights = sw)) <= eps))

  expect_lt(rms_dev(ow0$weights, bw, sw = sw),
            rms_dev(ow$weights, bw, sw = sw))

  expect_gt(mean_abs_dev(ow0$weights, bw, sw = sw),
            mean_abs_dev(ow$weights, bw, sw = sw))

  expect_lt(mean_abs_dev(ow$weights, bw, sw = sw),
            mean_abs_dev(ow$weights, sw = sw))

  expect_no_condition({
    ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     s.weights = "SW",
                     b.weights = bw,
                     norm = "l1",
                     tols = tols,
                     std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights,
                                        s.weights = sw)) <= attr(tols, "internal.tols") + eps))

  expect_gt(mean_abs_dev(ow$weights, sw = sw),
            mean_abs_dev(ow2$weights, sw = sw))

  expect_lt(mean_abs_dev(ow2$weights, bw, sw = sw),
            mean_abs_dev(ow2$weights, sw = sw))

  if (rlang::is_installed("scs")) {
    expect_no_condition({
      ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      s.weights = "SW",
                      b.weights = bw,
                      norm = "entropy")
    })

    expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights,
                                          s.weights = sw)) <= eps))

    expect_lt(rms_dev(ow0$weights, bw, sw = sw),
              rms_dev(ow$weights, bw, sw = sw))

    expect_gt(rel_ent(ow0$weights, bw, sw = sw),
              rel_ent(ow$weights, bw, sw = sw))

    expect_lt(rel_ent(ow$weights, bw, sw = sw),
              rel_ent(ow$weights, sw = sw))

    expect_no_condition({
      ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       s.weights = "SW",
                       b.weights = bw,
                       norm = "entropy",
                       tols = tols,
                       std.binary = TRUE)
    })

    expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights,
                                          s.weights = sw)) <= attr(tols, "internal.tols") + eps))

    expect_gt(rel_ent(ow$weights, sw = sw),
              rel_ent(ow2$weights, sw = sw))

    expect_lt(rel_ent(ow2$weights, bw, sw = sw),
              rel_ent(ow2$weights, sw = sw))

    expect_no_condition({
      ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      s.weights = "SW",
                      b.weights = bw,
                      norm = "log")
    })

    expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights,
                                          s.weights = sw)) <= eps))

    expect_lt(rms_dev(ow0$weights, bw, sw = sw),
              rms_dev(ow$weights, bw, sw = sw))

    expect_gt(sum(-sw * log(ow0$weights / bw)),
              sum(-sw * log(ow$weights / bw)))

    expect_no_condition({
      ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       s.weights = "SW",
                       b.weights = bw,
                       norm = "log",
                       tols = tols,
                       std.binary = TRUE)
    })

    expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights,
                                          s.weights = sw)) <= attr(tols, "internal.tols") + eps))

    expect_gt(sum(-sw * log(ow$weights / bw)),
              sum(-sw * log(ow2$weights / bw)))
  }
})

test_that("optweight works, binary, b.weights, solver: highs", {
  skip_if_not_installed("cobalt")
  skip_if_not_installed("highs")

  rlang::local_options(optweight_solver_l2 = "highs",
                       optweight_solver_l1 = "highs",
                       optweight_solver_linf = "highs",
                       optweight_solver_entropy = "scs",
                       optweight_solver_log = "scs")

  eps <- if (capabilities("long.double")) 5e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  bw <- runif(nrow(test_data))

  expect_no_condition({
    ow0 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow0$covs, ow0$treat, ow0$weights)) <= eps))

  expect_true(all(ow0$weights >= 1e-8))

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow0$weights))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    estimand = "ATE",
                    b.weights = bw,
                    tols = 0,
                    norm = "l2",
                    min.w = 1e-8)
  })

  expect_equal(ow$weights,
               ow0$weights,
               tolerance = eps)

  # tols
  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    tols = .02,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= .02 + eps))

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow$weights))

  tols <- process_tols(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       tols = seq(.01, .06, by = .01))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    tols = tols,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= attr(tols, "internal.tols") + eps))

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow0$weights, bw))

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow$weights))

  #norms
  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    norm = "linf",
                    eps = 2e-5)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow$weights, bw))

  expect_gt(max_abs_dev(ow0$weights, bw),
            max_abs_dev(ow$weights, bw))

  expect_lt(max_abs_dev(ow$weights, bw),
            max_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw,
                     norm = "linf",
                     tols = tols,
                     std.binary = TRUE,
                     eps = 1e-5)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

  expect_gt(max_abs_dev(ow$weights, bw),
            max_abs_dev(ow2$weights, bw))

  expect_lt(max_abs_dev(ow2$weights, bw),
            max_abs_dev(ow2$weights))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    norm = "l1",
                    eps = 1e-5)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow$weights, bw))

  expect_gt(mean_abs_dev(ow0$weights, bw),
            mean_abs_dev(ow$weights, bw))

  expect_lt(mean_abs_dev(ow$weights, bw),
            mean_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw,
                     norm = "l1",
                     tols = tols,
                     std.binary = TRUE,
                     eps = 5e-5)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

  expect_gt(mean_abs_dev(ow$weights, bw),
            mean_abs_dev(ow2$weights, bw))

  expect_lt(mean_abs_dev(ow2$weights, bw),
            mean_abs_dev(ow2$weights))

  if (rlang::is_installed("scs")) {
    expect_no_condition({
      ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      b.weights = bw,
                      norm = "entropy")
    })

    expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_lt(rms_dev(ow0$weights, bw),
              rms_dev(ow$weights, bw))

    expect_gt(rel_ent(ow0$weights, bw),
              rel_ent(ow$weights, bw))

    expect_lt(rel_ent(ow$weights, bw),
              rel_ent(ow$weights))

    expect_no_condition({
      ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       b.weights = bw,
                       norm = "entropy",
                       tols = tols,
                       std.binary = TRUE)
    })

    expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_gt(rel_ent(ow$weights, bw),
              rel_ent(ow2$weights, bw))

    expect_lt(rel_ent(ow2$weights, bw),
              rel_ent(ow2$weights))

    expect_no_condition({
      ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      b.weights = bw,
                      norm = "log")
    })

    expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_lt(rms_dev(ow0$weights, bw),
              rms_dev(ow$weights, bw))

    expect_gt(sum(-log(ow0$weights / bw)),
              sum(-log(ow$weights / bw)))

    expect_no_condition({
      ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       b.weights = bw,
                       norm = "log",
                       tols = tols,
                       std.binary = TRUE)
    })

    expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_gt(sum(-log(ow$weights / bw)),
              sum(-log(ow2$weights / bw)))
  }
})

test_that("optweight works, binary, b.weights, solver: lpsolve", {
  skip_if_not_installed("cobalt")
  skip_if_not_installed("lpSolve")

  rlang::local_options(optweight_solver_l2 = "osqp",
                       optweight_solver_l1 = "lpsolve",
                       optweight_solver_linf = "lpsolve",
                       optweight_solver_entropy = "scs",
                       optweight_solver_log = "scs")

  eps <- if (capabilities("long.double")) 5e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  bw <- runif(nrow(test_data))

  expect_no_condition({
    ow0 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow0$covs, ow0$treat, ow0$weights)) <= eps))

  expect_true(all(ow0$weights >= 1e-8))

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow0$weights))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    estimand = "ATE",
                    b.weights = bw,
                    tols = 0,
                    norm = "l2",
                    min.w = 1e-8)
  })

  expect_equal(ow$weights,
               ow0$weights,
               tolerance = eps)

  # tols
  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    tols = .02,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= .02 + eps))

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow$weights))

  tols <- process_tols(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       tols = seq(.01, .06, by = .01))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    tols = tols,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= attr(tols, "internal.tols") + eps))

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow0$weights, bw))

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow$weights))

  #norms
  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    norm = "linf",
                    eps = 2e-5)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow$weights, bw))

  expect_gt(max_abs_dev(ow0$weights, bw),
            max_abs_dev(ow$weights, bw))

  expect_lt(max_abs_dev(ow$weights, bw),
            max_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw,
                     norm = "linf",
                     tols = tols,
                     std.binary = TRUE,
                     eps = 1e-4)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

  expect_gt(max_abs_dev(ow$weights, bw),
            max_abs_dev(ow2$weights, bw))

  expect_lt(max_abs_dev(ow2$weights, bw),
            max_abs_dev(ow2$weights))

  expect_no_condition({
    ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    norm = "l1",
                    eps = 1e-5)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow$weights, bw))

  expect_gt(mean_abs_dev(ow0$weights, bw),
            mean_abs_dev(ow$weights, bw))

  expect_lt(mean_abs_dev(ow$weights, bw),
            mean_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw,
                     norm = "l1",
                     tols = tols,
                     std.binary = TRUE,
                     eps = 5e-5)
  })

  expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

  expect_gt(mean_abs_dev(ow$weights, bw),
            mean_abs_dev(ow2$weights, bw))

  expect_lt(mean_abs_dev(ow2$weights, bw),
            mean_abs_dev(ow2$weights))

  if (rlang::is_installed("scs")) {
    expect_no_condition({
      ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      b.weights = bw,
                      norm = "entropy")
    })

    expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_lt(rms_dev(ow0$weights, bw),
              rms_dev(ow$weights, bw))

    expect_gt(rel_ent(ow0$weights, bw),
              rel_ent(ow$weights, bw))

    expect_lt(rel_ent(ow$weights, bw),
              rel_ent(ow$weights))

    expect_no_condition({
      ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       b.weights = bw,
                       norm = "entropy",
                       tols = tols,
                       std.binary = TRUE)
    })

    expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_gt(rel_ent(ow$weights, bw),
              rel_ent(ow2$weights, bw))

    expect_lt(rel_ent(ow2$weights, bw),
              rel_ent(ow2$weights))

    expect_no_condition({
      ow <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      b.weights = bw,
                      norm = "log")
    })

    expect_true(all(abs(cobalt::col_w_smd(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_lt(rms_dev(ow0$weights, bw),
              rms_dev(ow$weights, bw))

    expect_gt(sum(-log(ow0$weights / bw)),
              sum(-log(ow$weights / bw)))

    expect_no_condition({
      ow2 <- optweight(A ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       b.weights = bw,
                       norm = "log",
                       tols = tols,
                       std.binary = TRUE)
    })

    expect_true(all(abs(cobalt::col_w_smd(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_gt(sum(-log(ow$weights / bw)),
              sum(-log(ow2$weights / bw)))
  }
})

test_that("optweight works, continuous", {
  skip_if_not_installed("cobalt")

  rlang::local_options(optweight_solver_l2 = "osqp",
                       optweight_solver_l1 = "osqp",
                       optweight_solver_linf = "osqp",
                       optweight_solver_entropy = "scs",
                       optweight_solver_log = "scs")

  eps <- if (capabilities("long.double")) 5e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  bw <- runif(nrow(test_data))

  expect_no_condition({
    ow0 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow0$covs, ow0$treat, ow0$weights)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow0$covs, ow0$treat)),
               cobalt::col_w_mean(cbind(ow0$covs, ow0$treat), ow0$weights),
               tolerance = eps)

  expect_true(all(ow0$weights >= 1e-8))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    estimand = "ATE",
                    tols = 0,
                    norm = "l2",
                    min.w = 1e-8)
  })

  expect_equal(ow$weights,
               ow0$weights,
               tolerance = eps)

  # tols
  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    tols = .02)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= .02 + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  tols <- process_tols(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       tols = seq(.01, .06, by = .01))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    tols = tols,
                    std.binary = TRUE)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= attr(tols, "internal.tols") + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow$weights),
            rms_dev(ow0$weights))

  #norms
  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    norm = "linf",
                    eps = 1e-5)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow0$weights),
            rms_dev(ow$weights))

  expect_gt(max_abs_dev(ow0$weights),
            max_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     norm = "linf",
                     tols = tols,
                     eps = 1e-5)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
               cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
               tolerance = eps)

  expect_gt(max_abs_dev(ow$weights),
            max_abs_dev(ow2$weights))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    norm = "l1",
                    eps = 1e-5)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow0$weights),
            rms_dev(ow$weights))

  expect_gt(mean_abs_dev(ow0$weights),
            mean_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     norm = "l1",
                     tols = tols)
  })

  # Note: L1 fairly inaccurate with osqp
  expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + 2 * eps))

  expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
               cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
               tolerance = eps)

  expect_gt(mean_abs_dev(ow$weights),
            mean_abs_dev(ow2$weights))

  if (rlang::is_installed("scs")) {
    expect_no_condition({
      ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      norm = "entropy")
    })

    expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
                 cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
                 tolerance = eps)

    expect_lt(rms_dev(ow0$weights),
              rms_dev(ow$weights))

    expect_gt(rel_ent(ow0$weights),
              rel_ent(ow$weights))

    expect_no_condition({
      ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       norm = "entropy",
                       tols = tols)
    })

    expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
                 cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
                 tolerance = eps)

    expect_gt(rel_ent(ow$weights),
              rel_ent(ow2$weights))

    expect_no_condition({
      ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      norm = "log")
    })

    expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
                 cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
                 tolerance = eps)

    expect_lt(rms_dev(ow0$weights),
              rms_dev(ow$weights))

    expect_gt(sum(-log(ow0$weights)),
              sum(-log(ow$weights)))

    expect_no_condition({
      ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       norm = "log",
                       tols = tols)
    })

    expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
                 cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
                 tolerance = eps)

    expect_gt(sum(-log(ow$weights)),
              sum(-log(ow2$weights)))
  }
})

test_that("optweight works, continuous, s.weights", {
  skip_if_not_installed("cobalt")

  rlang::local_options(optweight_solver_l2 = "osqp",
                       optweight_solver_l1 = "osqp",
                       optweight_solver_linf = "osqp",
                       optweight_solver_entropy = "scs",
                       optweight_solver_log = "scs")

  eps <- if (capabilities("long.double")) 5e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  bw <- runif(nrow(test_data))

  expect_no_condition({
    ow0 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data, s.weights = "SW")
  })

  sw <- test_data$SW

  expect_true(all(abs(cobalt::col_w_corr(ow0$covs, ow0$treat, ow0$weights,
                                        s.weights = sw)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow0$covs, ow0$treat),
                                  s.weights = sw),
               cobalt::col_w_mean(cbind(ow0$covs, ow0$treat), ow0$weights,
                                  s.weights = sw),
               tolerance = eps)

  expect_true(all(ow0$weights >= 1e-8))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    estimand = "ATE",
                    s.weights = "SW",
                    tols = 0,
                    norm = "l2",
                    min.w = 1e-8)
  })

  expect_equal(ow$weights,
               ow0$weights,
               tolerance = eps)

  # tols
  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    tols = .02)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights,
                                        s.weights = sw)) <= .02 + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat),
                                  s.weights = sw),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights,
                                  s.weights = sw),
               tolerance = eps)

  tols <- process_tols(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       tols = seq(.01, .06, by = .01))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    tols = tols)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights,
                                        s.weights = sw)) <= attr(tols, "internal.tols") + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat),
                                  s.weights = sw),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights,
                                  s.weights = sw),
               tolerance = eps)

  expect_lt(rms_dev(ow$weights, sw = sw),
            rms_dev(ow0$weights, sw = sw))

  #norms
  expect_error({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    norm = "linf")
  })

  expect_error({
    ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     s.weights = "SW",
                     norm = "linf",
                     tols = tols)
  })

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    norm = "l1")
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights,
                                        s.weights = sw)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat),
                                  s.weights = sw),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights,
                                  s.weights = sw),
               tolerance = eps)

  expect_lt(rms_dev(ow0$weights, sw = sw),
            rms_dev(ow$weights, sw = sw))

  expect_gt(mean_abs_dev(ow0$weights, sw = sw),
            mean_abs_dev(ow$weights, sw = sw))

  expect_no_condition({
    ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     s.weights = "SW",
                     norm = "l1",
                     tols = tols)
  })

  # Note: osqp fairly inaccurate with L1
  expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights,
                                        s.weights = sw)) <= attr(tols, "internal.tols") + 2 * eps))

  expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat),
                                  s.weights = sw),
               cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights,
                                  s.weights = sw),
               tolerance = eps)

  expect_gt(mean_abs_dev(ow$weights, sw = sw),
            mean_abs_dev(ow2$weights, sw = sw))

  if (rlang::is_installed("scs")) {
    expect_no_condition({
      ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      s.weights = "SW",
                      norm = "entropy")
    })

    expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights,
                                          s.weights = sw)) <= eps))

    expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat),
                                    s.weights = sw),
                 cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights,
                                    s.weights = sw),
                 tolerance = eps)

    expect_lt(rms_dev(ow0$weights, sw = sw),
              rms_dev(ow$weights, sw = sw))

    expect_gt(rel_ent(ow0$weights, sw = sw),
              rel_ent(ow$weights, sw = sw))

    expect_no_condition({
      ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       s.weights = "SW",
                       norm = "entropy",
                       tols = tols)
    })

    expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights,
                                          s.weights = sw)) <= attr(tols, "internal.tols") + eps))

    expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat),
                                    s.weights = sw),
                 cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights,
                                    s.weights = sw),
                 tolerance = eps)

    expect_gt(rel_ent(ow$weights, sw = sw),
              rel_ent(ow2$weights, sw = sw))

    expect_no_condition({
      ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      s.weights = "SW",
                      norm = "log")
    })

    expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights,
                                          s.weights = sw)) <= eps))

    expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat),
                                    s.weights = sw),
                 cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights,
                                    s.weights = sw),
                 tolerance = eps)

    expect_lt(rms_dev(ow0$weights, sw = sw),
              rms_dev(ow$weights, sw = sw))

    expect_gt(sum(-sw * log(ow0$weights)),
              sum(-sw * log(ow$weights)))

    expect_no_condition({
      ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       s.weights = "SW",
                       norm = "log",
                       tols = tols)
    })

    expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights,
                                          s.weights = sw)) <= attr(tols, "internal.tols") + eps))

    expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat),
                                    s.weights = sw),
                 cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights,
                                    s.weights = sw),
                 tolerance = eps)

    expect_gt(sum(-sw * log(ow$weights)),
              sum(-sw * log(ow2$weights)))
  }
})

test_that("optweight works, continuous, b.weights", {
  skip_if_not_installed("cobalt")

  rlang::local_options(optweight_solver_l2 = "osqp",
                       optweight_solver_l1 = "osqp",
                       optweight_solver_linf = "osqp",
                       optweight_solver_entropy = "scs",
                       optweight_solver_log = "scs")

  eps <- if (capabilities("long.double")) 5e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  bw <- runif(nrow(test_data))

  expect_no_condition({
    ow0 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow0$covs, ow0$treat, ow0$weights)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow0$covs, ow0$treat)),
               cobalt::col_w_mean(cbind(ow0$covs, ow0$treat), ow0$weights),
               tolerance = eps)

  expect_true(all(ow0$weights >= 1e-8))

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow0$weights))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    estimand = "ATE",
                    b.weights = bw,
                    tols = 0,
                    norm = "l2",
                    min.w = 1e-8)
  })

  expect_equal(ow$weights,
               ow0$weights,
               tolerance = eps)

  # tols
  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    tols = .02)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= .02 + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow$weights))

  tols <- process_tols(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       tols = seq(.01, .06, by = .01))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    tols = tols)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= attr(tols, "internal.tols") + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow0$weights, bw))

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow$weights))

  #norms
  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    norm = "linf",
                    eps = 1e-5)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow$weights, bw))

  expect_gt(max_abs_dev(ow0$weights, bw),
            max_abs_dev(ow$weights, bw))

  expect_lt(max_abs_dev(ow$weights, bw),
            max_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw,
                     norm = "linf",
                     tols = tols)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
               cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
               tolerance = eps)

  expect_gt(max_abs_dev(ow$weights, bw),
            max_abs_dev(ow2$weights, bw))

  expect_lt(max_abs_dev(ow2$weights, bw),
            max_abs_dev(ow2$weights))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    norm = "l1")
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow$weights, bw))

  expect_gt(mean_abs_dev(ow0$weights, bw),
            mean_abs_dev(ow$weights, bw))

  expect_lt(mean_abs_dev(ow$weights, bw),
            mean_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw,
                     norm = "l1",
                     tols = tols)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + 2 * eps))

  expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
               cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
               tolerance = eps)

  expect_gt(mean_abs_dev(ow$weights, bw),
            mean_abs_dev(ow2$weights, bw))

  expect_lt(mean_abs_dev(ow2$weights, bw),
            mean_abs_dev(ow2$weights))

  if (rlang::is_installed("scs")) {
    expect_no_condition({
      ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      b.weights = bw,
                      norm = "entropy")
    })

    expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
                 cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
                 tolerance = eps)

    expect_lt(rms_dev(ow0$weights, bw),
              rms_dev(ow$weights, bw))

    expect_gt(rel_ent(ow0$weights, bw),
              rel_ent(ow$weights, bw))

    expect_lt(rel_ent(ow$weights, bw),
              rel_ent(ow$weights))

    expect_no_condition({
      ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       b.weights = bw,
                       norm = "entropy",
                       tols = tols)
    })

    expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
                 cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
                 tolerance = eps)

    expect_gt(rel_ent(ow$weights, bw),
              rel_ent(ow2$weights, bw))

    expect_lt(rel_ent(ow2$weights, bw),
              rel_ent(ow2$weights))

    expect_no_condition({
      ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      b.weights = bw,
                      norm = "log")
    })

    expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
                 cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
                 tolerance = eps)

    expect_lt(rms_dev(ow0$weights, bw),
              rms_dev(ow$weights, bw))

    expect_gt(sum(-log(ow0$weights / bw)),
              sum(-log(ow$weights / bw)))

    expect_no_condition({
      ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       b.weights = bw,
                       norm = "log",
                       tols = tols)
    })

    expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
                 cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
                 tolerance = eps)

    expect_gt(sum(-log(ow$weights / bw)),
              sum(-log(ow2$weights / bw)))
  }
})

test_that("optweight works, continuous, s.weights + b.weights", {
  skip_if_not_installed("cobalt")

  rlang::local_options(optweight_solver_l2 = "osqp",
                       optweight_solver_l1 = "osqp",
                       optweight_solver_linf = "osqp",
                       optweight_solver_entropy = "scs",
                       optweight_solver_log = "scs")

  eps <- if (capabilities("long.double")) 5e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  bw <- runif(nrow(test_data))

  expect_no_condition({
    ow0 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data, s.weights = "SW",
                     b.weights = bw)
  })

  sw <- test_data$SW

  expect_true(all(abs(cobalt::col_w_corr(ow0$covs, ow0$treat, ow0$weights,
                                        s.weights = sw)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow0$covs, ow0$treat),
                                  s.weights = sw),
               cobalt::col_w_mean(cbind(ow0$covs, ow0$treat), ow0$weights,
                                  s.weights = sw),
               tolerance = eps)

  expect_true(all(ow0$weights >= 1e-8))

  expect_lt(rms_dev(ow0$weights, bw, sw = sw),
            rms_dev(ow0$weights, sw = sw))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    estimand = "ATE",
                    s.weights = "SW",
                    b.weights = bw,
                    tols = 0,
                    norm = "l2",
                    min.w = 1e-8)
  })

  expect_equal(ow$weights,
               ow0$weights,
               tolerance = eps)

  # tols
  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    b.weights = bw,
                    tols = .02)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights,
                                        s.weights = sw)) <= .02 + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat),
                                  s.weights = sw),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights,
                                  s.weights = sw),
               tolerance = eps)

  expect_lt(rms_dev(ow$weights, bw, sw = sw),
            rms_dev(ow$weights, sw = sw))

  tols <- process_tols(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       tols = seq(.01, .06, by = .01))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    b.weights = bw,
                    tols = tols)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights,
                                        s.weights = sw)) <= attr(tols, "internal.tols") + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat),
                                  s.weights = sw),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights,
                                  s.weights = sw),
               tolerance = eps)

  expect_lt(rms_dev(ow$weights, bw, sw = sw),
            rms_dev(ow0$weights, bw, sw = sw))

  expect_lt(rms_dev(ow$weights, bw, sw = sw),
            rms_dev(ow$weights, sw = sw))

  #norms
  expect_error({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    b.weights = bw,
                    norm = "linf")
  })

  expect_error({
    ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     s.weights = "SW",
                     b.weights = bw,
                     norm = "linf",
                     tols = tols)
  })


  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    s.weights = "SW",
                    b.weights = bw,
                    norm = "l1")
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights,
                                        s.weights = sw)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat),
                                  s.weights = sw),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights,
                                  s.weights = sw),
               tolerance = eps)

  expect_lt(rms_dev(ow0$weights, bw, sw = sw),
            rms_dev(ow$weights, bw, sw = sw))

  expect_gt(mean_abs_dev(ow0$weights, bw, sw = sw),
            mean_abs_dev(ow$weights, bw, sw = sw))

  expect_lt(mean_abs_dev(ow$weights, bw, sw = sw),
            mean_abs_dev(ow$weights, sw = sw))

  expect_no_condition({
    ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     s.weights = "SW",
                     b.weights = bw,
                     norm = "l1",
                     tols = tols)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights,
                                        s.weights = sw)) <= attr(tols, "internal.tols") + 2 * eps))

  expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat),
                                  s.weights = sw),
               cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights,
                                  s.weights = sw),
               tolerance = eps)

  expect_gt(mean_abs_dev(ow$weights, sw = sw),
            mean_abs_dev(ow2$weights, sw = sw))

  expect_lt(mean_abs_dev(ow2$weights, bw, sw = sw),
            mean_abs_dev(ow2$weights, sw = sw))

  if (rlang::is_installed("scs")) {
    expect_no_condition({
      ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      s.weights = "SW",
                      b.weights = bw,
                      norm = "entropy")
    })

    expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights,
                                          s.weights = sw)) <= eps))

    expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat),
                                    s.weights = sw),
                 cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights,
                                    s.weights = sw),
                 tolerance = eps)

    expect_lt(rms_dev(ow0$weights, bw, sw = sw),
              rms_dev(ow$weights, bw, sw = sw))

    expect_gt(rel_ent(ow0$weights, bw, sw = sw),
              rel_ent(ow$weights, bw, sw = sw))

    expect_lt(rel_ent(ow$weights, bw, sw = sw),
              rel_ent(ow$weights, sw = sw))

    expect_no_condition({
      ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       s.weights = "SW",
                       b.weights = bw,
                       norm = "entropy",
                       tols = tols)
    })

    expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights,
                                          s.weights = sw)) <= attr(tols, "internal.tols") + eps))

    expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat),
                                    s.weights = sw),
                 cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights,
                                    s.weights = sw),
                 tolerance = eps)

    expect_gt(rel_ent(ow$weights, sw = sw),
              rel_ent(ow2$weights, sw = sw))

    expect_lt(rel_ent(ow2$weights, bw, sw = sw),
              rel_ent(ow2$weights, sw = sw))

    expect_no_condition({
      ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      s.weights = "SW",
                      b.weights = bw,
                      norm = "log")
    })

    expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights,
                                          s.weights = sw)) <= eps))

    expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat),
                                    s.weights = sw),
                 cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights,
                                    s.weights = sw),
                 tolerance = eps)

    expect_lt(rms_dev(ow0$weights, bw, sw = sw),
              rms_dev(ow$weights, bw, sw = sw))

    expect_gt(sum(-sw * log(ow0$weights / bw)),
              sum(-sw * log(ow$weights / bw)))

    expect_no_condition({
      ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       s.weights = "SW",
                       b.weights = bw,
                       norm = "log",
                       tols = tols)
    })

    expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights,
                                          s.weights = sw)) <= attr(tols, "internal.tols") + eps))

    expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat),
                                    s.weights = sw),
                 cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights,
                                    s.weights = sw),
                 tolerance = eps)

    expect_gt(sum(-sw * log(ow$weights / bw)),
              sum(-sw * log(ow2$weights / bw)))
  }
})

test_that("optweight works, continuous, b.weights, solver: highs", {
  skip_if_not_installed("cobalt")
  skip_if_not_installed("highs")

  rlang::local_options(optweight_solver_l2 = "highs",
                       optweight_solver_l1 = "highs",
                       optweight_solver_linf = "highs",
                       optweight_solver_entropy = "scs",
                       optweight_solver_log = "scs")

  eps <- if (capabilities("long.double")) 5e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  bw <- runif(nrow(test_data))

  expect_no_condition({
    ow0 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow0$covs, ow0$treat, ow0$weights)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow0$covs, ow0$treat)),
               cobalt::col_w_mean(cbind(ow0$covs, ow0$treat), ow0$weights),
               tolerance = eps)

  expect_true(all(ow0$weights >= 1e-8))

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow0$weights))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    estimand = "ATE",
                    b.weights = bw,
                    tols = 0,
                    norm = "l2",
                    min.w = 1e-8)
  })

  expect_equal(ow$weights,
               ow0$weights,
               tolerance = eps)

  # tols
  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    tols = .02)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= .02 + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow$weights))

  tols <- process_tols(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       tols = seq(.01, .06, by = .01))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    tols = tols)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= attr(tols, "internal.tols") + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow0$weights, bw))

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow$weights))

  #norms
  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    norm = "linf",
                    eps = 1e-5)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow$weights, bw))

  expect_gt(max_abs_dev(ow0$weights, bw),
            max_abs_dev(ow$weights, bw))

  expect_lt(max_abs_dev(ow$weights, bw),
            max_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw,
                     norm = "linf",
                     tols = tols)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
               cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
               tolerance = eps)

  expect_gt(max_abs_dev(ow$weights, bw),
            max_abs_dev(ow2$weights, bw))

  expect_lt(max_abs_dev(ow2$weights, bw),
            max_abs_dev(ow2$weights))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    norm = "l1")
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow$weights, bw))

  expect_gt(mean_abs_dev(ow0$weights, bw),
            mean_abs_dev(ow$weights, bw))

  expect_lt(mean_abs_dev(ow$weights, bw),
            mean_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw,
                     norm = "l1",
                     tols = tols)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + 2 * eps))

  expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
               cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
               tolerance = eps)

  expect_gt(mean_abs_dev(ow$weights, bw),
            mean_abs_dev(ow2$weights, bw))

  expect_lt(mean_abs_dev(ow2$weights, bw),
            mean_abs_dev(ow2$weights))

  if (rlang::is_installed("scs")) {
    expect_no_condition({
      ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      b.weights = bw,
                      norm = "entropy")
    })

    expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
                 cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
                 tolerance = eps)

    expect_lt(rms_dev(ow0$weights, bw),
              rms_dev(ow$weights, bw))

    expect_gt(rel_ent(ow0$weights, bw),
              rel_ent(ow$weights, bw))

    expect_lt(rel_ent(ow$weights, bw),
              rel_ent(ow$weights))

    expect_no_condition({
      ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       b.weights = bw,
                       norm = "entropy",
                       tols = tols)
    })

    expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
                 cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
                 tolerance = eps)

    expect_gt(rel_ent(ow$weights, bw),
              rel_ent(ow2$weights, bw))

    expect_lt(rel_ent(ow2$weights, bw),
              rel_ent(ow2$weights))

    expect_no_condition({
      ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      b.weights = bw,
                      norm = "log")
    })

    expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
                 cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
                 tolerance = eps)

    expect_lt(rms_dev(ow0$weights, bw),
              rms_dev(ow$weights, bw))

    expect_gt(sum(-log(ow0$weights / bw)),
              sum(-log(ow$weights / bw)))

    expect_no_condition({
      ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       b.weights = bw,
                       norm = "log",
                       tols = tols)
    })

    expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
                 cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
                 tolerance = eps)

    expect_gt(sum(-log(ow$weights / bw)),
              sum(-log(ow2$weights / bw)))
  }
})

test_that("optweight works, continuous, b.weights, solver: lpsolve", {
  skip_if_not_installed("cobalt")
  skip_if_not_installed("lpSolve")

  rlang::local_options(optweight_solver_l2 = "osqp",
                       optweight_solver_l1 = "lpsolve",
                       optweight_solver_linf = "lpsolve",
                       optweight_solver_entropy = "scs",
                       optweight_solver_log = "scs")

  eps <- if (capabilities("long.double")) 5e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123)
  bw <- runif(nrow(test_data))

  expect_no_condition({
    ow0 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow0$covs, ow0$treat, ow0$weights)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow0$covs, ow0$treat)),
               cobalt::col_w_mean(cbind(ow0$covs, ow0$treat), ow0$weights),
               tolerance = eps)

  expect_true(all(ow0$weights >= 1e-8))

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow0$weights))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    estimand = "ATE",
                    b.weights = bw,
                    tols = 0,
                    norm = "l2",
                    min.w = 1e-8)
  })

  expect_equal(ow$weights,
               ow0$weights,
               tolerance = eps)

  # tols
  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    tols = .02)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= .02 + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow$weights))

  tols <- process_tols(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       tols = seq(.01, .06, by = .01))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    tols = tols)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= attr(tols, "internal.tols") + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow0$weights, bw))

  expect_lt(rms_dev(ow$weights, bw),
            rms_dev(ow$weights))

  #norms
  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    norm = "linf",
                    eps = 1e-5)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow$weights, bw))

  expect_gt(max_abs_dev(ow0$weights, bw),
            max_abs_dev(ow$weights, bw))

  expect_lt(max_abs_dev(ow$weights, bw),
            max_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw,
                     norm = "linf",
                     tols = tols)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

  expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
               cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
               tolerance = eps)

  expect_gt(max_abs_dev(ow$weights, bw),
            max_abs_dev(ow2$weights, bw))

  expect_lt(max_abs_dev(ow2$weights, bw),
            max_abs_dev(ow2$weights))

  expect_no_condition({
    ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                    data = test_data,
                    b.weights = bw,
                    norm = "l1")
  })

  expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

  expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
               cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
               tolerance = eps)

  expect_lt(rms_dev(ow0$weights, bw),
            rms_dev(ow$weights, bw))

  expect_gt(mean_abs_dev(ow0$weights, bw),
            mean_abs_dev(ow$weights, bw))

  expect_lt(mean_abs_dev(ow$weights, bw),
            mean_abs_dev(ow$weights))

  expect_no_condition({
    ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                     data = test_data,
                     b.weights = bw,
                     norm = "l1",
                     tols = tols)
  })

  expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + 2 * eps))

  expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
               cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
               tolerance = eps)

  expect_gt(mean_abs_dev(ow$weights, bw),
            mean_abs_dev(ow2$weights, bw))

  expect_lt(mean_abs_dev(ow2$weights, bw),
            mean_abs_dev(ow2$weights))

  if (rlang::is_installed("scs")) {
    expect_no_condition({
      ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      b.weights = bw,
                      norm = "entropy")
    })

    expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
                 cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
                 tolerance = eps)

    expect_lt(rms_dev(ow0$weights, bw),
              rms_dev(ow$weights, bw))

    expect_gt(rel_ent(ow0$weights, bw),
              rel_ent(ow$weights, bw))

    expect_lt(rel_ent(ow$weights, bw),
              rel_ent(ow$weights))

    expect_no_condition({
      ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       b.weights = bw,
                       norm = "entropy",
                       tols = tols)
    })

    expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
                 cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
                 tolerance = eps)

    expect_gt(rel_ent(ow$weights, bw),
              rel_ent(ow2$weights, bw))

    expect_lt(rel_ent(ow2$weights, bw),
              rel_ent(ow2$weights))

    expect_no_condition({
      ow <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                      data = test_data,
                      b.weights = bw,
                      norm = "log")
    })

    expect_true(all(abs(cobalt::col_w_corr(ow$covs, ow$treat, ow$weights)) <= eps))

    expect_equal(cobalt::col_w_mean(cbind(ow$covs, ow$treat)),
                 cobalt::col_w_mean(cbind(ow$covs, ow$treat), ow$weights),
                 tolerance = eps)

    expect_lt(rms_dev(ow0$weights, bw),
              rms_dev(ow$weights, bw))

    expect_gt(sum(-log(ow0$weights / bw)),
              sum(-log(ow$weights / bw)))

    expect_no_condition({
      ow2 <- optweight(Ac ~ X1 + X2 + X3 + X4 + X5 + X6,
                       data = test_data,
                       b.weights = bw,
                       norm = "log",
                       tols = tols)
    })

    expect_true(all(abs(cobalt::col_w_corr(ow2$covs, ow2$treat, ow2$weights)) <= attr(tols, "internal.tols") + eps))

    expect_equal(cobalt::col_w_mean(cbind(ow2$covs, ow2$treat)),
                 cobalt::col_w_mean(cbind(ow2$covs, ow2$treat), ow2$weights),
                 tolerance = eps)

    expect_gt(sum(-log(ow$weights / bw)),
              sum(-log(ow2$weights / bw)))
  }
})
