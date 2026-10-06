# ============================================================
# Figures - indirect reference intervals.
#
# Both figures read only results/reference_intervals.csv, so this
# script runs without the patient-level dataset. Run estimate_intervals.R
# first.
#
#   Figure 1  the four red cell parameters with the largest correction
#   Figure 2  all 17 parameters, each interval scaled to its empirical range
# ============================================================


library(ggplot2)

library(dplyr)

DATA_DIR   <- Sys.getenv("DATA_DIR",   "data")
RESULT_DIR <- Sys.getenv("RESULT_DIR", "results")
FIG_DIR    <- Sys.getenv("FIG_DIR",    "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

## Plain, high-contrast style. Journals often print in greyscale, so colour
## is never the only thing that distinguishes the series.
theme_pub <- theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
    panel.border     = element_rect(linewidth = 0.4, colour = "grey35"),
    axis.text        = element_text(colour = "grey20"),
    axis.title       = element_text(face = "plain"),
    legend.key       = element_blank(),
    legend.background = element_blank(),
    strip.background = element_rect(fill = "grey93", colour = "grey35",
                                    linewidth = 0.4),
    plot.title       = element_text(face = "bold", size = 10, hjust = 0)
  )
theme_set(theme_pub)

save_fig <- function(plot, name, width_mm, height_mm) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot,
         width = width_mm, height = height_mm, units = "mm", dpi = 300)
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot,
         width = width_mm, height = height_mm, units = "mm")
  cat("wrote", name, ".png / .pdf\n")
}

# ============================================================
# Figure 1 - indirect reference intervals against the empirical
# percentiles, for the red cell parameters with the largest correction.
# Shows the asymmetric shift of the limits.
# ============================================================
ri <- read.csv(file.path(RESULT_DIR, "reference_intervals.csv"),
               stringsAsFactors = FALSE)
sel <- c("HGB", "HCT", "MCV", "MCH")
f4 <- ri %>% filter(parameter %in% sel)
f4$parameter <- factor(f4$parameter, levels = sel,
                       labels = c("Haemoglobin (g/L)", "Haematocrit (%)",
                                  "MCV (fL)", "MCH (pg)"))
f4$partition <- factor(f4$partition,
                       levels = c("Men >=15", "Women >=15", "Children 1-14"),
                       labels = c("Men \u226515", "Women \u226515",
                                  "Children 1-14"))

f4long <- rbind(
  transform(f4[, c("parameter", "partition")],
            lo = f4$emp_p2.5, hi = f4$emp_p97.5,
            type = "Empirical percentiles (mixed population)"),
  transform(f4[, c("parameter", "partition")],
            lo = f4$ri_lower, hi = f4$ri_upper,
            type = "refineR indirect reference interval")
)
f4long$type <- factor(f4long$type,
                      levels = c("Empirical percentiles (mixed population)",
                                 "refineR indirect reference interval"))
f4long$ypos <- ifelse(f4long$type == levels(f4long$type)[1], 0.18, -0.18)

p4 <- ggplot(f4long,
             aes(y = as.numeric(partition) + ypos, colour = type)) +
  geom_segment(aes(x = lo, xend = hi,
                   yend = as.numeric(partition) + ypos),
               linewidth = 1.5, lineend = "round") +
  facet_wrap(~ parameter, scales = "free_x", nrow = 2) +
  scale_y_continuous(breaks = 1:3, labels = levels(f4long$partition),
                     limits = c(0.5, 3.5)) +
  scale_colour_manual(values = c("grey65", "#B2182B"), name = NULL) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "top",
        legend.direction = "vertical",
        panel.grid.major.y = element_blank())

save_fig(p4, "figure1_reference_intervals", 183, 110)


## ---- Figure 2: all 17 parameters at a glance ----------------------
## Every parameter, with each interval expressed as a fraction of its
## empirical range so that parameters on different scales share one panel.
ri <- read.csv(file.path(RESULT_DIR, "reference_intervals.csv"),
               stringsAsFactors = FALSE)

ri2 <- ri %>%
  mutate(width = emp_p97.5 - emp_p2.5,
         rel_lo = (ri_lower - emp_p2.5) / width,
         rel_hi = (ri_upper - emp_p2.5) / width,
         partition = factor(partition,
                            levels = c("Men >=15", "Women >=15", "Children 1-14"),
                            labels = c("Men \u226515", "Women \u226515",
                                       "Children 1\u201314")))

ord <- c("HGB","RBC","HCT","MCV","MCH","MCHC","RDW_CV","PLT","MPV","PCT",
         "PDW","WBC","NEU_ABS","LYMPH_ABS","MON_ABS","EOS_ABS","BAS_ABS")
## Axis labels as in Tables 1-3, in place of the CSV column names.
lab <- c("Haemoglobin","Red cell count","Haematocrit","MCV","MCH","MCHC",
         "RDW-CV","Platelets","MPV","Plateletcrit","PDW","Leukocytes",
         "Neutrophils","Lymphocytes","Monocytes","Eosinophils","Basophils")
stopifnot(length(ord) == length(lab))
ri2$parameter <- factor(ri2$parameter, levels = rev(ord), labels = rev(lab))

p5b <- ggplot(ri2) +
  annotate("rect", xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf,
           fill = "grey92") +
  geom_segment(aes(x = rel_lo, xend = rel_hi,
                   y = parameter, yend = parameter),
               linewidth = 1.6, colour = "#B2182B", lineend = "round") +
  geom_vline(xintercept = c(0, 1), colour = "grey55", linewidth = 0.35) +
  facet_wrap(~ partition) +
  scale_x_continuous(breaks = c(0, 0.5, 1),
                     labels = c("empirical\nP2.5", "", "empirical\nP97.5")) +
  labs(x = NULL, y = NULL,
       title = NULL) +
  theme(axis.text.x = element_text(size = 6.5),
        axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"),
        panel.grid.major.x = element_blank())

save_fig(p5b, "figure2_all_parameters", 183, 105)

cat("\nFigures written to", FIG_DIR, "\n")
