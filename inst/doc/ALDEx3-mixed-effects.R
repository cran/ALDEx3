## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = TRUE--------------------------------------------------------------
set.seed(42)
library(ALDEx3)
data("oral_mouthwash_data", package = "ALDEx3")

Y <- oral_mouthwash_data$counts
keep_names <- row.names(Y[((rowSums(Y == 0)) / ncol(Y)) <= 0.75, ])
other <- colSums(Y[((rowSums(Y == 0)) / ncol(Y)) > 0.75, ])
Y <- rbind(Y[keep_names, ], other = other)

meta <- oral_mouthwash_data$metadata

fit_blmm <- aldex(
  Y,
  ~ treat * timec + (1 | participant_id),
  data = meta,
  method = "blmm",
  n.cores = 1,
  nsample = 250,
  scale = clr.sm,
  gamma = 0
)

head(summary(fit_blmm))

## ----eval = TRUE--------------------------------------------------------------
set.seed(42)
data("oral_mouthwash_data", package = "ALDEx3")

Y0 <- oral_mouthwash_data$counts
keep_names <- row.names(Y0[((rowSums(Y0 == 0)) / ncol(Y0)) <= 0.75, ])
other <- colSums(Y0[((rowSums(Y0 == 0)) / ncol(Y0)) > 0.75, ])
Y <- rbind(Y0[keep_names, ], other = other)

# Keep a small subset to make the comparison easy to inspect in the vignette.
Y <- Y[c(seq_len(min(8, nrow(Y) - 1)), nrow(Y)), ]
meta <- oral_mouthwash_data$metadata

K <- 1

fit_args <- list(
  data = meta,
  n.cores = 1,
  nsample = 16,
  scale = clr.sm,
  gamma = 0
)

fit_model <- function(method) {
  do.call(
    aldex,
    c(
      list(Y, ~ treat * timec + (1 | participant_id)),
      fit_args,
      list(method = method)
    )
  )
}

set.seed(42)
fit_lme4 <- fit_model("lme4")

set.seed(42)
fit_blmm <- fit_model("blmm")

time_lme4 <- replicate(K, {
  set.seed(42)
  system.time({
    invisible(fit_model("lme4"))
  })[["elapsed"]]
})

time_blmm <- replicate(K, {
  set.seed(42)
  system.time({
    invisible(fit_model("blmm"))
  })[["elapsed"]]
})

s_lme4 <- summary(fit_lme4)
s_blmm <- summary(fit_blmm)

cmp <- merge(
  s_lme4,
  s_blmm,
  by = c("parameter", "entity"),
  suffixes = c(".lme4", ".blmm")
)

cmp$abs_est_diff <- abs(cmp$estimate.lme4 - cmp$estimate.blmm)
cmp$abs_se_diff <- abs(cmp$std.error.lme4 - cmp$std.error.blmm)

cmp_tbl <- cmp[order(cmp$abs_est_diff, decreasing = TRUE), ][1:8, c(
  "parameter", "entity",
  "estimate.lme4", "estimate.blmm", "abs_est_diff",
  "std.error.lme4", "std.error.blmm", "abs_se_diff"
)]

cat('<div style="overflow-x:auto; width:100%;">\n')
knitr::kable(
  cmp_tbl,
  format = "html",
  digits = 4,
  caption = "Largest absolute differences between exact lme4 and blmm on the shared-draw subset."
)
cat('\n</div>\n')

runtime_tbl <- data.frame(
  engine = c("lme4", "blmm"),
  mean_elapsed_sec = c(mean(time_lme4), mean(time_blmm)),
  median_elapsed_sec = c(median(time_lme4), median(time_blmm)),
  speedup_vs_lme4 = c(1, mean(time_lme4) / mean(time_blmm))
)

cat('\n\n')
cat('<div style="overflow-x:auto; width:100%;">\n')
knitr::kable(
  runtime_tbl,
  format = "html",
  digits = 3,
  caption = paste0(
    "Elapsed runtime on the shared-draw vignette subset averaged over K = ",
    K,
    " timing run(s). Times will vary by hardware."
  )
)
cat('\n</div>\n')

