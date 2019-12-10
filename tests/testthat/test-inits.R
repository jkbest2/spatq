context("Test initial value estimates")

init_ref <- list(beta_n = c(`X_nfactor(time)1` = 0.587256825835408,
                            `X_nfactor(time)2` = 0.588205033094543,
                            `X_nfactor(time)3` = 0.626845789494911,
                            `X_nfactor(time)4` = 0.529744815469454,
                            `X_nfactor(time)5` = 0.521775416950299,
                            `X_nfactor(time)6` = 0.583392497921106,
                            `X_nfactor(time)7` = 0.547056321267843,
                            `X_nfactor(time)8` = 0.532524114409541,
                            `X_nfactor(time)9` = 0.521260303845579,
                            `X_nfactor(time)10` = 0.564833675336079,
                            `X_nfactor(time)11` = 0.534444289890176,
                            `X_nfactor(time)12` = 0.576547219052125,
                            `X_nfactor(time)13` = 0.523163219182192,
                            `X_nfactor(time)14` = 0.536185564263769,
                            `X_nfactor(time)15` = 0.579262285977399),
                 beta_w = c(`X_wfactor(time)1` = -9.66137703953886,
                            `X_wfactor(time)2` = -9.59897106599928,
                            `X_wfactor(time)3` = -9.81339803797218,
                            `X_wfactor(time)4` = -9.62052759991199,
                            `X_wfactor(time)5` = -9.57110816901162,
                            `X_wfactor(time)6` = -9.80899973301897,
                            `X_wfactor(time)7` = -9.98472714599767,
                            `X_wfactor(time)8` = -9.97653221464713,
                            `X_wfactor(time)9` = -9.82775902612745,
                            `X_wfactor(time)10` = -9.92465282444752,
                            `X_wfactor(time)11` = -10.0126089055202,
                            `X_wfactor(time)12` = -10.0586928698776,
                            `X_wfactor(time)13` = -9.93390607066071,
                            `X_wfactor(time)14` = -9.91918236385518,
                            `X_wfactor(time)15` = -10.1222641902379 ),
                 lambda_n = c(R_n = 0.10077216742229),
                 lambda_w = c(R_w = 0.879102621182762))

init_ref <- list(beta_n = c(0.587256825835408, 0.588205033094543,
                            0.626845789494911, 0.529744815469454,
                            0.521775416950299, 0.583392497921106,
                            0.547056321267843, 0.532524114409541,
                            0.521260303845579, 0.564833675336079,
                            0.534444289890176, 0.576547219052125,
                            0.523163219182192, 0.536185564263769,
                            0.579262285977399),
                 beta_w = c(-9.66137703953886, -9.59897106599928,
                            -9.81339803797218, -9.62052759991199,
                            -9.57110816901162, -9.80899973301897,
                            -9.98472714599767, -9.97653221464713,
                            -9.82775902612745, -9.92465282444752,
                            -10.0126089055202, -10.0586928698776,
                            -9.93390607066071, -9.91918236385518,
                            -10.1222641902379),
                 lambda_n = 0.10077216742229,
                 lambda_w = 0.879102621182762)
init_ref_len <- vapply(init_ref, length, 1)

init <- init_fixef(obj_sim$env$data)
init_len <- vapply(init, length, 1)

parlen <- vapply(c("X_n", "X_w", "R_n", "R_w"),
                 function(nm) ncol(obj_sim$env$data[[nm]]), 1)

test_that("Initial values are correct length", {
  expect_equal(names(init), names(init_ref))
  expect_equal(init_len, init_ref_len)
  expect_equal(unname(init_len), unname(parlen))
})
