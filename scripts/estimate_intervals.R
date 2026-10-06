# ============================================================
# Indirect reference intervals for complete blood count parameters
# in an Uzbek outpatient population - refineR analysis.
#
# Requires the refineR package from CRAN:
#     install.packages("refineR")
#
# refineR (Ammer et al., Sci Rep 2021;11:16023) models the
# non-pathological subpopulation as a Box-Cox transformed normal
# distribution and estimates its parameters by minimising a regularised
# cost function, making no assumption about the distribution of the
# pathological results.
#
# Outputs written:
#   reference_intervals.csv        - the reference intervals (Tables 1-3)
#   reference_interval_plots.pdf   - one diagnostic plot per estimate
#   session_info.txt               - session info for reproducibility
# ============================================================

if (!requireNamespace("refineR", quietly = TRUE)) {
  stop("refineR is not installed. Run:  install.packages(\"refineR\")")
}
library(refineR)
library(dplyr)

DATA_DIR   <- Sys.getenv("DATA_DIR",   "data")       # de-identified analytic datasets
RESULT_DIR <- Sys.getenv("RESULT_DIR", "results")    # analysis outputs
data <- read.csv(file.path(DATA_DIR, "cbc_validated.csv"),
                 stringsAsFactors = FALSE)
cat("Loaded", nrow(data), "encounters\n\n")

## ---- Partitions ----
## Sex partitioning is applied from age 15, where the WHO haemoglobin
## thresholds become sex-specific (120 g/L for both sexes at 12-14 years;
## 120 g/L for women and 130 g/L for men from 15). Children are not
## sex-partitioned because red cell parameters do not differ materially
## by sex before puberty.
data$partition <- with(data, ifelse(
  AGE < 15, "Children 1-14",
  ifelse(SEX == "Man", "Men >=15", "Women >=15")))

## ---- Parameters ----
## refineR's "BoxCox" model handles both symmetric and right-skewed
## analytes, so no manual log transformation is applied.
PARAMS <- c("HGB", "RBC", "HCT", "MCV", "MCH", "MCHC", "RDW_CV",
            "PLT", "MPV", "PCT", "PDW", "WBC",
            "NEU_ABS", "LYMPH_ABS", "MON_ABS", "EOS_ABS", "BAS_ABS")

## Minimum n per partition. refineR is reliable from a few hundred
## observations; every partition here is far larger.
MIN_N <- 500

results <- list()
pdf(file.path(RESULT_DIR, "reference_interval_plots.pdf"), width = 9, height = 6)

for (prt in unique(data$partition)) {
  sub <- data[data$partition == prt, ]
  for (p in PARAMS) {
    v <- suppressWarnings(as.numeric(sub[[p]]))
    v <- v[is.finite(v)]
    # Zeros are kept. refineR accepts them, and a zero basophil or eosinophil
    # count is a normal result, not a failed measurement; dropping zeros
    # would cut off the lower tail of these parameters.
    if (length(v) < MIN_N) {
      cat(sprintf("  [%s] %s: skipped, n=%d below minimum\n", prt, p, length(v)))
      next
    }

    cat(sprintf("Fitting [%s] %s  (n=%d) ...\n", prt, p, length(v)))
    # The confidence bounds are percentiles of the bootstrap replicates, so
    # with 200 replicates each 90% bound rests on 10 of them (5 with 100).
    est <- try(findRI(Data = v, model = "BoxCox", NBootstrap = 200), silent = TRUE)
    if (inherits(est, "try-error")) {
      cat(sprintf("  [%s] %s: findRI failed\n", prt, p))
      next
    }

    # 95% reference interval with a 90% confidence interval on each limit.
    # CIprop = 0.90 because CLSI EP28-A3c specifies 90% for reference limits
    # (the getRI default is 0.95). The limits themselves come from the fit to
    # the full data (pointEst = "fullDataEst", the default), not from the
    # bootstrap.
    ri <- try(getRI(est, RIperc = c(0.025, 0.975), CIprop = 0.90), silent = TRUE)
    if (inherits(ri, "try-error")) {
      cat(sprintf("  [%s] %s: getRI failed\n", prt, p))
      next
    }

    plot(est, title = paste0(prt, " - ", p))

    results[[length(results) + 1]] <- data.frame(
      partition   = prt,
      parameter   = p,
      n           = length(v),
      ri_lower    = round(ri$PointEst[1], 3),
      ri_upper    = round(ri$PointEst[2], 3),
      lower_ci_lo = round(ri$CILow[1], 3),
      lower_ci_hi = round(ri$CIHigh[1], 3),
      upper_ci_lo = round(ri$CILow[2], 3),
      upper_ci_hi = round(ri$CIHigh[2], 3),
      emp_p2.5    = round(quantile(v, 0.025), 3),
      emp_p97.5   = round(quantile(v, 0.975), 3),
      stringsAsFactors = FALSE
    )

    # The full run takes hours. Results are saved after every fit, so an
    # interruption loses one parameter, not the whole run. They go to a
    # separate file so the published CSV is never replaced by a partial one.
    write.csv(bind_rows(results),
              file.path(RESULT_DIR, "reference_intervals_partial.csv"),
              row.names = FALSE)
  }
}
dev.off()

ri_all <- bind_rows(results)
row.names(ri_all) <- NULL

## ---- Assessment of separation ----
## An interval that reproduces the empirical percentiles of the mixed
## distribution has failed to separate health from disease. The
## displacement of each limit from its empirical counterpart is expressed
## as a percentage of the empirical interval width.
emp_width <- ri_all$emp_p97.5 - ri_all$emp_p2.5
ri_all$disp_lower_pct <- round(100 * abs(ri_all$ri_lower - ri_all$emp_p2.5) / emp_width, 1)
ri_all$disp_upper_pct <- round(100 * abs(ri_all$ri_upper - ri_all$emp_p97.5) / emp_width, 1)

cat("\n=== INDIRECT REFERENCE INTERVALS (refineR) ===\n")
for (prt in unique(ri_all$partition)) {
  cat("---", prt, "---\n")
  print(ri_all[ri_all$partition == prt,
               c("parameter", "n", "ri_lower", "ri_upper",
                 "emp_p2.5", "emp_p97.5", "disp_lower_pct", "disp_upper_pct")],
        row.names = FALSE)
  cat("\n")
}

disp_max <- pmax(ri_all$disp_lower_pct, ri_all$disp_upper_pct)
cat(sprintf("Displacement of the more shifted limit: %.1f-%.1f%%, median %.1f%%, above 5%%: %d of %d\n\n",
            min(disp_max), max(disp_max), median(disp_max),
            sum(disp_max > 5), length(disp_max)))

write.csv(ri_all, file.path(RESULT_DIR, "reference_intervals.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(RESULT_DIR, "session_info.txt"))

cat("Saved: reference_intervals.csv\n")
cat("Saved: reference_interval_plots.pdf  (inspect these before using any interval)\n")
cat("Saved: session_info.txt\n")

# ============================================================
# After running, check each estimate before using it:
#
# 1. In reference_interval_plots.pdf, the estimated healthy distribution
#    should follow the main peak of the histogram and leave out the
#    pathological tail. Discard any estimate where it does not.
#
# 2. If disp_lower_pct and disp_upper_pct are both near zero, the method
#    has not separated health from disease for that parameter.
#
# 3. A wide confidence interval (lower_ci_lo..lower_ci_hi,
#    upper_ci_lo..upper_ci_hi) means the limit is poorly determined, even
#    if the estimate looks sensible.
# ============================================================
