# ============================================================
# Figure 3: derived adult reference intervals against published values.
#
# Reads only results/reference_intervals.csv and the published values below,
# so the figure reproduces without patient data. This study's intervals come
# from the results file (Tables 1-3); the published values match Table 4.
#
# Each interval is drawn as a segment between its two limits. The published
# sources give only the limits, so a boxplot would show spread that was
# never reported.
#
# The published values are the same as in compare_published.py and must be
# kept in step with it by hand. The check at the end of this file confirms
# that every plotted interval from this study matches the results file.
# ============================================================

suppressMessages({
  library(ggplot2)
  library(dplyr)
})

RESULT_DIR <- Sys.getenv("RESULT_DIR", "results")
FIG_DIR    <- Sys.getenv("FIG_DIR",    "figures")

ri <- read.csv(file.path(RESULT_DIR, "reference_intervals.csv"),
               stringsAsFactors = FALSE)

## ---- published values ----
## Ozarda 2017 (PMC5493180); sex-specific only for HGB, RBC, HCT.
## Dacie & Lewis Table 2.2, mean +/-2SD converted to the 95% range.
## HCT converted from L/L to %.
pub <- read.csv(text = "
parameter,partition,source,lo,hi
HGB,Men >=15,Ozarda 2017,131,175
HGB,Women >=15,Ozarda 2017,110,152
RBC,Men >=15,Ozarda 2017,4.43,6.07
RBC,Women >=15,Ozarda 2017,3.96,5.31
HCT,Men >=15,Ozarda 2017,39.2,52.2
HCT,Women >=15,Ozarda 2017,33.7,46.1
MCV,both,Ozarda 2017,77.2,95.7
MCH,both,Ozarda 2017,25.2,32.2
MCHC,both,Ozarda 2017,319,350
RDW_CV,both,Ozarda 2017,12.2,16.3
PLT,both,Ozarda 2017,152,383
MPV,both,Ozarda 2017,6.2,11.8
WBC,both,Ozarda 2017,4.39,11.59
NEU_ABS,both,Ozarda 2017,2.04,7.54
LYMPH_ABS,both,Ozarda 2017,1.21,3.77
MON_ABS,both,Ozarda 2017,0.26,0.94
EOS_ABS,both,Ozarda 2017,0.02,0.50
BAS_ABS,both,Ozarda 2017,0.01,0.12
HGB,Men >=15,Dacie and Lewis,130,170
HGB,Women >=15,Dacie and Lewis,120,150
RBC,Men >=15,Dacie and Lewis,4.5,5.5
RBC,Women >=15,Dacie and Lewis,3.8,4.8
HCT,Men >=15,Dacie and Lewis,40,50
HCT,Women >=15,Dacie and Lewis,36,46
MCV,both,Dacie and Lewis,83,101
MCH,both,Dacie and Lewis,27,32
MCHC,both,Dacie and Lewis,315,345
RDW_CV,both,Dacie and Lewis,11.6,14.0
PLT,both,Dacie and Lewis,150,410
WBC,both,Dacie and Lewis,4,10
NEU_ABS,both,Dacie and Lewis,2,7
LYMPH_ABS,both,Dacie and Lewis,1,3
MON_ABS,both,Dacie and Lewis,0.2,1.0
EOS_ABS,both,Dacie and Lewis,0.02,0.5
BAS_ABS,both,Dacie and Lewis,0.02,0.1
", stringsAsFactors = FALSE)

LAB <- c(HGB = "Haemoglobin", RBC = "Red cell count", HCT = "Haematocrit",
         MCV = "MCV", MCH = "MCH", MCHC = "MCHC", RDW_CV = "RDW-CV",
         PLT = "Platelets", MPV = "MPV", WBC = "Leukocytes",
         NEU_ABS = "Neutrophils", LYMPH_ABS = "Lymphocytes",
         MON_ABS = "Monocytes", EOS_ABS = "Eosinophils",
         BAS_ABS = "Basophils")
UNIT <- c(HGB = "g/L", RBC = "x10^12/L", HCT = "%", MCV = "fL", MCH = "pg",
          MCHC = "g/L", RDW_CV = "%", PLT = "x10^9/L", MPV = "fL",
          WBC = "x10^9/L", NEU_ABS = "x10^9/L", LYMPH_ABS = "x10^9/L",
          MON_ABS = "x10^9/L", EOS_ABS = "x10^9/L", BAS_ABS = "x10^9/L")

ADULTS <- c("Men >=15", "Women >=15")
PARAMS <- names(LAB)

## Published values given for both sexes are drawn against each partition.
pub_long <- do.call(rbind, lapply(ADULTS, function(p) {
  rows <- pub[pub$partition == p | pub$partition == "both", ]
  rows$partition <- p
  rows
}))

derived <- ri %>%
  filter(partition %in% ADULTS, parameter %in% PARAMS) %>%
  transmute(parameter, partition, source = "This study",
            lo = ri_lower, hi = ri_upper)

dat <- rbind(derived, pub_long[, names(derived)]) %>%
  mutate(
    label  = factor(LAB[parameter], levels = rev(LAB[PARAMS])),
    source = factor(source,
                    levels = c("This study", "Ozarda 2017", "Dacie and Lewis")),
    facet  = paste0(LAB[parameter], " (", UNIT[parameter], ")"),
    sex    = ifelse(partition == "Men >=15", "Men", "Women"))

## Each parameter has its own x scale, since their units differ widely.
dat$facet <- factor(dat$facet,
                    levels = paste0(LAB[PARAMS], " (", UNIT[PARAMS], ")"))

## Journals often print in greyscale, where three shades of grey are hard to
## tell apart, so each source also has its own line width. The sources appear
## in the same order in every panel.
p <- ggplot(dat, aes(y = sex, colour = source, linewidth = source,
                     group = source)) +
  geom_linerange(aes(xmin = lo, xmax = hi),
                 position = position_dodge(width = 0.7)) +
  geom_point(aes(x = lo), position = position_dodge(width = 0.7),
             size = 1.4, shape = 124, show.legend = FALSE) +
  geom_point(aes(x = hi), position = position_dodge(width = 0.7),
             size = 1.4, shape = 124, show.legend = FALSE) +
  facet_wrap(~ facet, scales = "free_x", ncol = 3) +
  scale_colour_manual(values = c("This study" = "black",
                                 "Ozarda 2017" = "grey40",
                                 "Dacie and Lewis" = "grey72")) +
  scale_linewidth_manual(values = c("This study" = 1.3,
                                    "Ozarda 2017" = 0.85,
                                    "Dacie and Lewis" = 0.6)) +
  ## Keep intervals clear of the panel border.
  scale_x_continuous(expand = expansion(mult = 0.06)) +
  guides(linewidth = "none") +
  labs(x = NULL, y = NULL, colour = NULL) +
  theme_bw(base_size = 8) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(size = 7),
        axis.text = element_text(size = 6.5))

ggsave(file.path(FIG_DIR, "figure3_published_comparison.png"), p,
       width = 183, height = 200, units = "mm", dpi = 300)
ggsave(file.path(FIG_DIR, "figure3_published_comparison.pdf"), p,
       width = 183, height = 200, units = "mm")

cat("wrote figure3_published_comparison .png / .pdf\n")

## ---- check: plotted intervals match the results file ----
chk <- read.csv(file.path(RESULT_DIR, "reference_intervals.csv"),
                stringsAsFactors = FALSE)
for (i in seq_len(nrow(derived))) {
  r <- chk[chk$partition == derived$partition[i] &
           chk$parameter == derived$parameter[i], ]
  stopifnot(nrow(r) == 1,
            isTRUE(all.equal(derived$lo[i], r$ri_lower)),
            isTRUE(all.equal(derived$hi[i], r$ri_upper)))
}
cat(sprintf("all %d derived intervals match the results file\n", nrow(derived)))
