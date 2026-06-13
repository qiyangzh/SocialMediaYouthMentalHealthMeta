########################################################################################################
# Social Media Use Meta-meta-analysis with Jungup
########################################################################################################
# Authors: Qiyang Zhang
# Contact: qiyang39@nus.edu.sg
# Created: 2025/05/14
# Revised: 2026/6/12

# This file analyzes the included studies in the School-based Mental Health Interventions systematic 
# review, including preparing the data for analysis and meta-regressions.

########################################################################################################
# Initial Set-up
########################################################################################################
# Clear workspace
rm(list=ls(all=TRUE))

# Load packages
test<-require(googledrive)   #all gs_XXX() functions for reading data from Google
if (test == FALSE) {
  install.packages("googledrive")
  require(googledrive)
}
test<-require(googlesheets4)   #all gs_XXX() functions for reading data from Google
if (test == FALSE) {
  install.packages("googlesheets4")
  require(googlesheets4)
}
test<-require(plyr)   #rename()
if (test == FALSE) {
  install.packages("plyr")
  require(plyr)
}
test<-require(metafor)   #escalc(); rma();
if (test == FALSE) {
  install.packages("metafor")
  require(metafor)
}
test<-require(robumeta)
if (test == FALSE) {
  install.packages("robumeta")
  require(robumeta)
}
test<-require(weightr) #selection modeling
if (test == FALSE) {
  install.packages("weightr")
  require(weightr)
}
test<-require(clubSandwich) #coeftest
if (test == FALSE) {
  install.packages("clubSandwich")
  require(clubSandwich)
}
test<-require(tableone)   #CreateTableOne()
if (test == FALSE) {
  install.packages("tableone")
  require(tableone)
}
test<-require(flextable)   
if (test == FALSE) {
  install.packages("flextable")
  require(flextable)
}
test<-require(officer)   
if (test == FALSE) {
  install.packages("officer")
  require(officer)
}
test<-require(tidyverse)   
if (test == FALSE) {
  install.packages("tidyverse")
  require(tidyverse)
}
test<-require(ggrepel)   
if (test == FALSE) {
  install.packages("ggrepel")
  require(ggrepel)
}
test<-require(meta)   
if (test == FALSE) {
  install.packages("meta")
  require(meta)
}
test<-require(dplyr)  
if (test == FALSE) {
  install.packages("dplyr")
  require(dplyr)
}
test<-require(ggplot2)  
if (test == FALSE) {
  install.packages("ggplot2")
  require(ggplot2)
}
test<-require(tibble)  
if (test == FALSE) {
  install.packages("tibble")
  require(tibble)
}
test<-require(patchwork)  
if (test == FALSE) {
  install.packages("patchwork")
  require(patchwork)
}
test<-require(stringr)  
if (test == FALSE) {
  install.packages("stringr")
  require(stringr)
}
test<-require(tidyr)  
if (test == FALSE) {
  install.packages("tidyr")
  require(tidyr)
}
test<-require(forcats)  
if (test == FALSE) {
  install.packages("forcats")
  require(forcats)
}
rm(test)

########################################################################################################
# Load data
########################################################################################################
# set up to load from Google
drive_auth(email = "zhangqiyang0329@gmail.com")
id <- drive_find(pattern = "Social Media Use Review", type = "spreadsheet")$id[1]

# load findings and studies
gs4_auth(email = "zhangqiyang0329@gmail.com")
findings <- read_sheet(id, sheet = "Findings", col_types = "c")
studies <- read_sheet(id, sheet = "Studies", col_types = "c")   # includes separate effect sizes for each finding from a study
# Load raw data
raw <- read_sheet(id, sheet = "Longitudinal_Studies", col_types = "c")
rm(id)

########################################################################################################
# Clean data
########################################################################################################
# remove any empty rows & columns
studies <- subset(studies, is.na(studies$Study)==FALSE)
findings <- subset(findings, is.na(findings$Study)==FALSE)

studies <- subset(studies, is.na(studies$Drop)==TRUE)
findings <- subset(findings, is.na(findings$Drop)==TRUE)

# merge dataframes
full <- merge(studies, findings, by = c("Study"), all = FALSE, suffixes = c(".s", ".f"))

# format to correct variable types
nums <- c("Sample size", "Effect.Size", "Mean.age", 
          "Female.percent", "CI_upper", "CI_lower",
          "Above16", "Above15", "FiftyFivePercentFemale")

full[nums] <- lapply(full[nums], as.numeric)
rm(nums)
###############################################################
#Create unique identifiers (ES, study, program)
###############################################################
full$ESId <- seq_len(nrow(full))
full$StudyID <- as.numeric(as.factor(full$Study))
summary(full$StudyID)

full$Effect.Size <- as.numeric(full$Effect.Size)

################################################################
# Calculate meta-analytic variables: Variances (Lipsey & Wilson, 2000, Eq. 3.23)
################################################################
#calculate standard errors
full <- full[!is.na(full$Effect.Size) & !is.na(full$CI_lower) & !is.na(full$CI_upper), ]
full <- full[full$CI_upper >= full$CI_lower, ]
full$se  <- (full$CI_upper - full$CI_lower) / (2 * 1.96)
full$var <- full$se^2
full <- full[!is.na(full$var) & full$var > 0, ]
###########################
#forest plot
###########################
#-----------------------------
# 1. Fit RVE null model
#-----------------------------
forest_model <- robu(
  formula = Effect.Size ~ 1,
  studynum = StudyID,
  data = full,
  var.eff.size = var,
  rho = 0.8,
  small = TRUE
)

# Extract overall pooled estimate and 95% CI
overall_es <- forest_model$reg_table$b.r[1]
overall_lb <- forest_model$reg_table$CI.L[1]
overall_ub <- forest_model$reg_table$CI.U[1]

#-----------------------------
# 2. Prepare ES-level plotting data
#-----------------------------
plot_dat <- full %>%
  arrange(StudyID)

# Show study label only on first ES from each study
plot_dat$study_label <- ifelse(
  duplicated(plot_dat$StudyID),
  "",
  plot_dat$Study
)

# Optional: create effect size labels within study
plot_dat$es_label <- paste0("ES ", seq_len(nrow(plot_dat)))

# Row positions: leave extra space at bottom for overall effect
k <- nrow(plot_dat)
rows <- seq(from = k + 2, to = 3, by = -1)

#-----------------------------
# 3. Draw ES-level forest plot
#-----------------------------
forest <- forest(
  x = plot_dat$Effect.Size,
  vi = plot_dat$var,
  slab = plot_dat$study_label,
  rows = rows,
  xlab = "Effect Size",
  alim = c(min(plot_dat$Effect.Size - 1.96 * sqrt(plot_dat$var), na.rm = TRUE),
           max(plot_dat$Effect.Size + 1.96 * sqrt(plot_dat$var), na.rm = TRUE)),
  cex = 1
)
# 
# # Optional: add horizontal separators between studies
# study_end_rows <- rows[c(which(!duplicated(plot_dat$StudyID))[-1] - 1)]
# abline(h = study_end_rows + 0.5, lty = "dotted", col = "gray80")

#-----------------------------
# 4. Add overall pooled effect at the bottom
#-----------------------------
addpoly(
  x = overall_es,
  ci.lb = overall_lb,
  ci.ub = overall_ub,
  row = 1,
  mlab = sprintf(
    "Overall effect (RVE) = %.2f [95%% CI %.2f, %.2f]",
    overall_es, overall_lb, overall_ub
  ),
  cex = 1
)
########################################
#meta-regression
########################################
#Null Model
V_list <- impute_covariance_matrix(vi=full$var, cluster=full$StudyID, r=0.8)

MVnull <- rma.mv(yi=Effect.Size,
                 V=V_list,
                 random=~1 | StudyID/ESId,
                 test="t",
                 data=full,
                 method="REML")
MVnull

# #t-test of each covariate#
MVnull.coef <- coef_test(MVnull, cluster=full$StudyID, vcov="CR2")
MVnull.coef

#moderator analysis
full_moderator <- full %>%
  dplyr::filter(Study != "Fassi et al. 2024")

full_moderator$Above16.c <- full_moderator$Above16 - mean(full_moderator$Above16)
full_moderator$Above15.c <- full_moderator$Above15 - mean(full_moderator$Above15)
full_moderator$FiftyFivePercentFemale.c <- full_moderator$FiftyFivePercentFemale - mean(full_moderator$FiftyFivePercentFemale)


terms <- c("Above16.c")
terms <- c("FiftyFivePercentFemale.c")
#terms <- c("Mean.age", "Female.percent")
#terms <- c("Above16.c", "FiftyFivePercentFemale.c")
formula <- reformulate(termlabels = c(terms))
formula
full_moderator$ESId <- seq_len(nrow(full_moderator))
full_moderator$StudyID <- as.numeric(as.factor(full_moderator$Study))
summary(full$StudyID)
V_list <- impute_covariance_matrix(vi=full_moderator$var, cluster=full_moderator$StudyID, r=0.8)

MVfull <- rma.mv(yi=Effect.Size,
                 V=V_list,
                 mods=formula,
                 random=~1 | StudyID/ESId,
                 test="t",
                 data=full_moderator,
                 method="REML")
MVfull

#t-test of each covariate#
MVfull.coef <- coef_test(MVfull, cluster=full_moderator$StudyID, vcov="CR2")
MVfull.coef

###sensitivity analysis1
full_sens1 <- full %>%
  dplyr::filter(Study != "Fassi et al. 2024")
V_listfull_sens1 <- impute_covariance_matrix(vi=full_sens1$var, cluster=full_sens1$StudyID, r=0.8)

MVnullfull_sens1 <- rma.mv(yi=Effect.Size,
                 V=V_listfull_sens1,
                 random=~1 | StudyID/ESId,
                 test="t",
                 data=full_sens1,
                 method="REML")
MVnullfull_sens1

# #t-test of each covariate#
MVnull.coeffull_sens1 <- coef_test(MVnullfull_sens1, cluster=full_sens1$StudyID, vcov="CR2")
MVnull.coeffull_sens1

###sensitivity analysis2
full_sens2 <- full %>%
  dplyr::filter(Study != "Marciano et al. 2022")
V_listfull_sens2 <- impute_covariance_matrix(vi=full_sens2$var, cluster=full_sens2$StudyID, r=0.8)

MVnullfull_sens2 <- rma.mv(yi=Effect.Size,
                           V=V_listfull_sens2,
                           random=~1 | StudyID/ESId,
                           test="t",
                           data=full_sens2,
                           method="REML")
MVnullfull_sens2

# #t-test of each covariate#
MVnull.coeffull_sens2 <- coef_test(MVnullfull_sens2, cluster=full_sens2$StudyID, vcov="CR2")
MVnull.coeffull_sens2
#################################################################################
# Heterogeneity
#################################################################################
# 95% prediction intervals
print(PI_upper <- MVnull$b[1] + (1.96*sqrt(MVnull$sigma2[1] + MVnull$sigma2[2])))
print(PI_lower <- MVnull$b[1] - (1.96*sqrt(MVnull$sigma2[1] + MVnull$sigma2[2])))

#################################################################################
# Marginal Means
#################################################################################
# re-run model for each moderator to get marginal means for each #
# set up table to store results
means <- data.frame(moderator = character(0), group = character(0), beta = numeric(0), SE = numeric(0),
                    tstat = numeric(0), df = numeric(0), p_Satt = numeric(0))
mods <- c("as.factor(Above16)", "as.factor(FiftyFivePercentFemale)")

for(i in 1:length(mods)){
  # i <- 1
  formula <- reformulate(termlabels = c(mods[i], terms, "-1"))   # Worth knowing - if you duplicate terms, it keeps the first one
  mod_means <- rma.mv(yi=Effect.Size, #effect size
                      V = V_list, #variance (tHIS IS WHAt CHANGES FROM HEmodel)
                      mods = formula, #ADD COVS HERE
                      random = ~1 | StudyID/ESId, #nesting structure
                      test= "t", #use t-tests
                      data=full_moderator, #define data
                      method="REML") #estimate variances using REML
  coef_mod_means <- as.data.frame(coef_test(mod_means,#estimation model above
                                            cluster=full_moderator$StudyID, #define cluster IDs
                                            vcov = "CR2")) #estimation method (CR2 is best)
  # limit to relevant rows (the means you are interested in)
  coef_mod_means$moderator <- gsub(x = mods[i], pattern = "as.factor", replacement = "")
  coef_mod_means$group <- rownames(coef_mod_means)
  rownames(coef_mod_means) <- c()
  coef_mod_means <- subset(coef_mod_means, substr(start = 1, stop = nchar(mods[i]), x = coef_mod_means$group)== mods[i])
  coef_mod_means$group <- substr(x = coef_mod_means$group, start = nchar(mods[i])+1, stop = nchar(coef_mod_means$group))
  means <- dplyr::bind_rows(means, coef_mod_means)
}
means
########################
#Output officer
########################
myreport<-read_docx()

MVnull.coef
str(MVnull.coef)
MVnull.coef$coef <- row.names(as.data.frame(MVnull.coef))
row.names(MVnull.coef) <- c()
MVnull.coef <- MVnull.coef[c("Coef", "beta", "SE", "tstat", "df_Satt", "p_Satt")]
MVnull.coef
str(MVnull.coef)

MVfull.coef$coef <- row.names(as.data.frame(MVfull.coef))
row.names(MVfull.coef) <- c()
MVfull.coef <- MVfull.coef[c("Coef", "beta", "SE", "tstat", "df_Satt", "p_Satt")]

# MetaRegression Table
model_null <- flextable(head(MVnull.coef, n=nrow(MVnull.coef)))
colkeys <- c("beta", "SE", "tstat", "df_Satt")
model_null <- colformat_double(model_null,  j = colkeys, digits = 2)
model_null <- colformat_double(model_null,  j = c("p_Satt"), digits = 3)
#model_null <- autofit(model_null)
model_null <- add_header_lines(model_null, values = c("Null Model"), top = FALSE)
model_null <- theme_vanilla(model_null)

myreport <- body_add_par(x = myreport, value = "Table 5: Model Results", style = "Normal")
myreport <- body_add_flextable(x = myreport, model_null)

model_full <- flextable(head(MVfull.coef, n=nrow(MVfull.coef)))
model_full <- colformat_double(model_full,  j = c("beta"), digits = 2)
model_full <- colformat_double(model_full,  j = c("p_Satt"), digits = 3)
#model_full <- autofit(model_full)
model_full <- delete_part(model_full, part = "header")
model_full <- add_header_lines(model_full, values = c("Meta-Regression"))
model_full <- theme_vanilla(model_full)

myreport <- body_add_flextable(x = myreport, model_full)
myreport <- body_add_par(x = myreport, value = "", style = "Normal")

# Marginal Means Table
marginalmeans <- flextable(head(means, n=nrow(means)))
colkeys <- c("moderator", "group", "SE", "tstat", "df")
marginalmeans <- colformat_double(marginalmeans,  j = colkeys, digits = 2)
marginalmeans <- colformat_double(marginalmeans,  j = c("p_Satt"), digits = 3)
rm(colkeys)
marginalmeans <- theme_vanilla(marginalmeans)
marginalmeans <- merge_v(marginalmeans, j = c("moderator"))
myreport <- body_add_par(x = myreport, value = "Table: Marginal Means", style = "Normal")
myreport <- body_add_flextable(x = myreport, marginalmeans)
# Write to word doc
file = paste("TableResults.docx", sep = "")
print(myreport, file)


#publication bias , selection modeling
full_y <- full$Effect.Size
full_v <- full$var
# two-sided p-values for Z-tests under H0
p <- 2 * pnorm(-abs(full_y / sqrt(full_v)))

# Inspect distribution across bins to pick viable steps
table(cut(p, breaks = c(0, 0.025, 0.05, 0.10, 0.50, 1), include.lowest = TRUE))
weightfunct(full_y, full_v, steps = c(0.5, 1))
weightfunct(full_y, full_v, steps = c(.025, .50, 1))

#Total sample size
sum(full$Sample.size)

#funnel plot
MVfull.modeling <- rma(
  yi = Effect.Size,
  vi = var,
  test = "t",
  data = full,
  slab = Study,
  method = "REML"
)

pdf("contour_funnel_plot.pdf", width = 7, height = 6)

metafor::funnel(
  MVfull.modeling,
  xlab = "Effect size",
  ylab = "Standard Error",
  main = "Contour-Enhanced Funnel Plot",
  refline = coef(MVfull.modeling)[1],
  level = c(90, 95, 99),
  shade = c("white", "gray85", "gray70"),
  back = "white",
  pch = 21,
  bg = "gray60",
  col = "black",
  cex = 0.9
)

legend(
  "topright",
  legend = c("p > .10", ".05 < p ≤ .10", ".01 < p ≤ .05", "p ≤ .01"),
  fill = c("white", "gray85", "gray70", "gray55"),
  bty = "n",
  cex = 0.8
)

dev.off()

# heatmap with Quality score column
# Ensure variables are factors/characters as needed
full <- full %>%
  mutate(
    Y = as.factor(Y),
    Study = as.character(Study)
  )

# Correct quality-score table
quality_df <- tribble(
  ~Study,                 ~Quality.score,
  "Fassi et al. 2024",    10,
  "Ferguson et al. 2025", 9,
  "Ivie et al. 2020",     10,
  "Liu et al. 2022",      11,
  "Marciano et al. 2022", 9,
  "Nan et al. 2024",      9,
  "Yin et al. 2019",      9
)

# Aggregate by Study × Outcome
by_study_outcome <- full %>%
  dplyr::group_by(Study, Y) %>%
  dplyr::summarise(
    Year = first(as.integer(Year.s)),
    Effect.Size = mean(Effect.Size, na.rm = TRUE),
    .groups = "drop"
  )

# Keep only studies that appear in the quality table
by_study_outcome <- by_study_outcome %>%
  filter(Study %in% quality_df$Study)

# Order studies from newer to older
study_levels <- quality_df %>%
  mutate(Year = as.integer(sub(".*(20[0-9]{2})$", "\\1", Study))) %>%
  arrange(desc(Year), Study) %>%
  pull(Study)

# Apply same order to both datasets
by_study_outcome <- by_study_outcome %>%
  mutate(Study = factor(Study, levels = study_levels))

quality_df <- quality_df %>%
  mutate(Study = factor(Study, levels = study_levels))

# Left panel: Quality score column
p_quality <- ggplot(quality_df, aes(x = 1, y = Study)) +
  geom_tile(fill = "grey95", color = "white", linewidth = 0.3) +
  geom_text(aes(label = Quality.score), size = 6) +
  scale_x_continuous(
    breaks = 1,
    labels = "Quality score",
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Meta-analysis") +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_blank(),
    plot.margin = margin(5.5, 0, 5.5, 5.5)
  )

# Main heatmap
p_heatmap <- ggplot(by_study_outcome, aes(x = Y, y = Study, fill = Effect.Size)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", Effect.Size)), size = 6) +
  scale_fill_gradient2(
    low = "#2B6CB0",
    mid = "#F2F2F2",
    high = "#C53030",
    midpoint = 0,
    name = "Effect Size"
  ) +
  labs(
    x = "Mental Health Problems",
    y = NULL,
    title = "Heatmap of Effect Sizes by Outcome Category and Meta-analysis"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 0)
  )

# Combine
p_quality + p_heatmap + plot_layout(widths = c(1.2, 6))

# =============================================================================
# CLEAN COMBINED FIGURE: FOREST PLOT + SMALLER HEATMAP USING PATCHWORK
# Paste this AFTER your forest plot code and heatmap data-preparation code
# Required objects already created:
# plot_dat, rows, overall_es, overall_lb, overall_ub,
# quality_df, by_study_outcome
# =============================================================================

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("ggplotify", quietly = TRUE)) install.packages("ggplotify")

library(patchwork)
library(ggplotify)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Convert metafor forest plot to a patchwork-compatible object
# -----------------------------------------------------------------------------

p_forest <- ggplotify::as.ggplot(function() {
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(
    mar = c(4, 4, 3, 2),
    cex = 0.75
  )
  
  metafor::forest(
    x = plot_dat$Effect.Size,
    vi = plot_dat$var,
    slab = plot_dat$study_label,
    rows = rows,
    xlab = "Effect Size",
    alim = c(
      min(plot_dat$Effect.Size - 1.96 * sqrt(plot_dat$var), na.rm = TRUE),
      max(plot_dat$Effect.Size + 1.96 * sqrt(plot_dat$var), na.rm = TRUE)
    ),
    xlim = c(-1.4, 1.2),
    cex = 0.75,
    main = "Forest Plot of Effect Sizes",
    header = c("Study", "Estimate [95% CI]")
  )
  
  metafor::addpoly(
    x = overall_es,
    ci.lb = overall_lb,
    ci.ub = overall_ub,
    row = 1,
    mlab = sprintf(
      "Overall effect = %.2f [95%% CI %.2f, %.2f]",
      overall_es, overall_lb, overall_ub
    ),
    cex = 0.75
  )
})

# -----------------------------------------------------------------------------
# 2. Recreate quality-score panel with smaller text
# -----------------------------------------------------------------------------

p_quality_small <- ggplot(quality_df, aes(x = 1, y = Study)) +
  geom_tile(fill = "grey95", color = "white", linewidth = 0.3) +
  geom_text(aes(label = Quality.score), size = 3.0) +
  scale_x_continuous(
    breaks = 1,
    labels = "Quality score",
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Meta-analysis"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    axis.ticks = element_blank(),
    plot.margin = margin(3, 0, 3, 3)
  )

# -----------------------------------------------------------------------------
# 3. Recreate heatmap with smaller text
# -----------------------------------------------------------------------------

p_heatmap_small <- ggplot(by_study_outcome, aes(x = Y, y = Study, fill = Effect.Size)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", Effect.Size)), size = 3.0) +
  scale_fill_gradient2(
    low = "#2B6CB0",
    mid = "#F2F2F2",
    high = "#C53030",
    midpoint = 0,
    name = "Effect Size"
  ) +
  labs(
    x = "Mental health problems",
    y = NULL,
    title = "Heatmap of Effect Sizes by Outcome Category"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(3, 3, 3, 0)
  )

# -----------------------------------------------------------------------------
# 4. Combine quality-score panel and heatmap
# -----------------------------------------------------------------------------

p_heatmap_combined_small <- p_quality_small + p_heatmap_small +
  patchwork::plot_layout(widths = c(1.1, 6))

# -----------------------------------------------------------------------------
# 5. Combine forest plot and heatmap
# -----------------------------------------------------------------------------

combined_forest_heatmap_small <- p_forest / p_heatmap_combined_small +
  patchwork::plot_layout(heights = c(1.2, 1.6)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1)
  )

# Show combined figure
combined_forest_heatmap_small

# -----------------------------------------------------------------------------
# 6. Save figure
# -----------------------------------------------------------------------------

ggsave(
  filename = "combined_forest_heatmap_small.png",
  plot = combined_forest_heatmap_small,
  width = 14,
  height = 12,
  dpi = 300
)

ggsave(
  filename = "combined_forest_heatmap_small.pdf",
  plot = combined_forest_heatmap_small,
  width = 14,
  height = 12
)
#4 figures together panel: Year of meta
# Prepare year counts
plot_df_year <- full %>%
  dplyr::distinct(Study, .keep_all = TRUE) %>%   
  dplyr::mutate(Year = as.numeric(Year.f)) %>%
  dplyr::filter(!is.na(Year)) %>%
  dplyr::count(Year, name = "Count") %>%
  complete(Year = 2015:2026, fill = list(Count = 0))

# Plot
p_a <- ggplot(plot_df_year, aes(x = Year, y = Count)) +
  geom_col(
    fill = "grey50",
    color = "grey50",
    width = 0.9
  ) +
  geom_text(
    aes(label = Count),
    vjust = -0.3,
    size = 4
  ) +
  scale_x_continuous(
    breaks = 2019:2025,
    limits = c(2018.5, 2025.5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, max(plot_df_year$Count, na.rm = TRUE) + 1, 1),
    limits = c(0, max(plot_df_year$Count, na.rm = TRUE) + 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Year meta-analysis was made available",
    y = "Count of meta-analyses",
    title = "Year of meta-analyses"
  ) +
  theme_gray(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 11),   # show year labels
    axis.ticks.x = element_line(),           # show x-axis ticks
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
p_a

#age included
plot_df_age <- full %>%
  distinct(Study, .keep_all = TRUE) %>%   # one row per meta-analysis
  dplyr::mutate(
    `Minimum age` = as.numeric(`Minimum age`),
    `Maximum age` = as.numeric(`Maximum age`),
    `Mean age` = as.numeric(`Mean.age`)
  ) %>%
  dplyr::filter(!is.na(`Minimum age`), !is.na(`Maximum age`)) %>%
  dplyr::arrange(`Minimum age`, `Maximum age`) %>%
  dplyr::mutate(row_id = row_number())

plot_df_age
p_b <- ggplot(plot_df_age, aes(y = Study)) +
  geom_segment(
    aes(
      x = `Minimum age`,
      xend = `Maximum age`,
      yend = Study
    ),
    linewidth = 0.4,
    color = "grey20"
  ) +
  geom_point(
    aes(x = `Minimum age`),
    size = 1,
    color = "grey20"
  ) +
  geom_point(
    aes(x = `Maximum age`),
    size = 1,
    color = "grey20"
  ) +
  geom_point(
    aes(x = `Mean age`),
    size = 2.2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.8,
    na.rm = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 26),
    breaks = seq(0, 26, by = 1),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Age",
    y = "Meta-analysis",
    title = "Age range and mean age of participants"
  ) +
  theme_gray(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 11),
    axis.ticks.x = element_line(),
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p_b

# Year range included
plot_df_year_range <- full %>%
  distinct(Study, .keep_all = TRUE) %>%   # one row per study/meta-analysis
  mutate(
    `Year lower` = as.numeric(year.lower),
    `Year upper` = as.numeric(year.upper)
  ) %>%
  filter(!is.na(`Year lower`), !is.na(`Year upper`)) %>%
  arrange(`Year lower`, `Year upper`) %>%
  mutate(
    Studies = factor(Study, levels = Study)
  )

p_c <- ggplot(plot_df_year_range, aes(y = Studies)) +
  geom_segment(
    aes(
      x = `Year lower`,
      xend = `Year upper`,
      yend = Studies
    ),
    linewidth = 0.4,
    color = "grey20"
  ) +
  geom_point(
    aes(x = `Year lower`),
    size = 1,
    color = "grey20"
  ) +
  geom_point(
    aes(x = `Year upper`),
    size = 1,
    color = "grey20"
  ) +
  scale_x_continuous(
    limits = c(
      min(plot_df_year_range$`Year lower`, na.rm = TRUE) - 1,
      max(plot_df_year_range$`Year upper`, na.rm = TRUE) + 1
    ),
    breaks = seq(
      min(plot_df_year_range$`Year lower`, na.rm = TRUE),
      max(plot_df_year_range$`Year upper`, na.rm = TRUE),
      by = 2
    ),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Year",
    y = "Meta-analysis",
    title = "Publication year range of included studies"
  ) +
  theme_gray(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 11),
    axis.ticks.x = element_line(),
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p_c

# Clean and count mental health outcomes
plot_df_outcomes <- full %>%
  dplyr::distinct(Study, .keep_all = TRUE) %>%   # one row per meta-analysis
  dplyr::select(Study, Outcomes) %>%
  dplyr::filter(!is.na(Outcomes), Outcomes != "") %>%
  dplyr::mutate(
    outcome_raw = stringr::str_to_lower(Outcomes),
    outcome_raw = stringr::str_replace_all(outcome_raw, "\n", " "),
    outcome_raw = stringr::str_replace_all(outcome_raw, " and ", ", "),
    outcome_raw = stringr::str_replace_all(outcome_raw, ";", ",")
  ) %>%
  tidyr::separate_rows(outcome_raw, sep = ",") %>%
  dplyr::mutate(
    outcome = stringr::str_trim(outcome_raw),
    
    outcome = dplyr::case_when(
      stringr::str_detect(outcome, "depress") ~ "Depression",
      stringr::str_detect(outcome, "anxiety") ~ "Anxiety",
      stringr::str_detect(outcome, "stress") ~ "Stress",
      stringr::str_detect(outcome, "suicidal") ~ "Suicidal ideation",
      stringr::str_detect(outcome, "loneliness") ~ "Loneliness",
      stringr::str_detect(outcome, "negative affect") ~ "Negative affect",
      stringr::str_detect(outcome, "ill-being|illbeing") ~ "Ill-being",
      TRUE ~ stringr::str_to_title(outcome)
    )
  ) %>%
  dplyr::filter(outcome != "") %>%
  dplyr::count(outcome, name = "n") %>%
  dplyr::arrange(desc(n), outcome) %>%
  dplyr::mutate(
    outcome = forcats::fct_reorder(outcome, n),
    Label = n
  )

# Check the counted data
plot_df_outcomes

# Horizontal bar plot
p_d <- ggplot(plot_df_outcomes, aes(x = n, y = outcome)) +
  geom_col(
    fill = "grey60",
    color = "black",
    linewidth = 0.25,
    width = 0.75
  ) +
  geom_text(
    aes(label = Label),
    hjust = -0.2,
    size = 4
  ) +
  scale_x_continuous(
    limits = c(0, max(plot_df_outcomes$n, na.rm = TRUE) + 1),
    breaks = seq(0, max(plot_df_outcomes$n, na.rm = TRUE) + 1, by = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Frequency",
    y = "Mental health outcome",
    title = "Frequency of mental health outcomes"
  ) +
  theme_gray(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p_d
#combine them together, / (p_e | p_f)
install.packages("patchwork")  # only if not installed
library(patchwork)
combined_plot <- (p_a | p_b) / (p_c | p_d) +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0, 1)
  )

combined_plot
# Meta-Regression Analysis: Youth Social Media Use & Mental Health
# Longitudinal Studies
# =============================================================================
# Dependencies: metafor, dplyr, ggplot2, readxl, stringr
# Install if needed:
# install.packages(c("metafor", "dplyr", "ggplot2", "readxl", "stringr", "clubSandwich"))

library(metafor)
library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)

# =============================================================================
# SECTION 1: DATA LOADING & CLEANING
# =============================================================================
# Rename columns for convenience (adjust if your column names differ)
# Key columns based on the spreadsheet structure:
#   Study ID / reference label
#   Reference (full citation)
#   Country
#   Data collection (follow-up duration)
#   Age.mean, Age.range
#   Sample size
#   Female number
#   X (predictor), Y (outcome)
#   Reverse Direction  (1 = mental health -> social media; blank/0 = social media -> mental health)
#   Effect size        (Cohen's d, already computed)
#   Correlation        (r, used to derive d where available)
#   Sample size (n)

colnames(raw) <- make.names(colnames(raw), unique = TRUE)

# Preview structure
cat("=== Column names ===\n")
print(colnames(raw))
cat("\n=== Dimensions ===\n")
cat(nrow(raw), "rows x", ncol(raw), "columns\n\n")


# =============================================================================
# SECTION 2: FILTER TO USABLE EFFECT SIZES
# =============================================================================

# The column "Effect.size" holds Cohen's d values.
# Keep only rows that:
#   (a) have Drop != 1 (not dropped by reviewers)
#   (b) have a numeric, non-missing effect size
#   (c) are longitudinal (Longitudinal. == "L" or column not blank)

df <- raw %>%
  # Remove explicitly dropped rows (Drop column == 1)
  filter(is.na(Drop) | Drop != 1) %>%
  # Keep rows with a numeric effect size
  mutate(
    d = suppressWarnings(as.numeric(Effect.size)),
    reverse = suppressWarnings(as.numeric(Reverse.Direction))
  ) %>%
  filter(!is.na(d)) %>%
  # Recode reverse direction: 1 = MH -> SM; 0 = SM -> MH
  mutate(reverse = ifelse(is.na(reverse) | reverse == 0, 0L, 1L))

cat("Rows with usable effect sizes:", nrow(df), "\n\n")
# =============================================================================
# SECTION 3: HANDLE REVERSE-DIRECTION EFFECTS
#
# Two principled options. Choose ONE based on your research question.
#
# OPTION A — Sign-flip: make ALL effects represent SM -> MH direction.
#             Reverse-direction effects (MH -> SM) get their sign flipped
#             so a positive d = more SM predicts worse MH (or vice versa).
#             Use when you want ONE pooled estimate covering both directions.
#
# OPTION B — Subset: analyse forward (SM -> MH) and reverse (MH -> SM)
#             separately. Use when the two directions are distinct RQs.
# =============================================================================

# --- OPTION A:  ---
df <- df %>%
  mutate(
    d_unified = d,
    direction_label = ifelse(reverse == 1,
                             "MH -> Social Media (reversed)",
                             "SM -> Mental Health (forward)")
  )

# Sanity check
cat("Direction breakdown:\n")
print(table(df$direction_label))
cat("\n")
# =============================================================================
# SECTION 4: COMPUTE SAMPLING VARIANCE
#
# Cohen's d variance ≈  (n1+n2)/(n1*n2) + d^2/(2*(n1+n2))
#
# When only total N is available (no group split), for a within-person or
# regression-based d from a single-group longitudinal study we use:
#   vi = 1/n + d^2/(2*n)     [single-group approximation]
# =============================================================================

df <- df %>%
  mutate(
    n_total = suppressWarnings(as.numeric(Sample.size)),
    n_female = suppressWarnings(as.numeric(Female.number)),
    age_mean = suppressWarnings(as.numeric(Age.mean))
  ) %>%
  mutate(
    # Use single-group variance approximation
    vi = ifelse(
      !is.na(n_total) & n_total > 0,
      (1 / n_total) + (d_unified^2 / (2 * n_total)),
      NA_real_
    ),
    sei = sqrt(vi)
  ) %>%
  filter(!is.na(vi))

cat("Rows after variance computation:", nrow(df), "\n\n")


# =============================================================================
# SECTION 5: DERIVE MODERATORS
# =============================================================================

df <- df %>%
  mutate(
    # --- Follow-up duration in months ---
    followup_raw = tolower(as.character(Data.collection)),
    followup_months = case_when(
      str_detect(followup_raw, "year")  ~
        as.numeric(str_extract(followup_raw, "[0-9]+\\.?[0-9]*")) * 12,
      str_detect(followup_raw, "month") ~
        as.numeric(str_extract(followup_raw, "[0-9]+\\.?[0-9]*")),
      str_detect(followup_raw, "week")  ~
        as.numeric(str_extract(followup_raw, "[0-9]+\\.?[0-9]*")) / 4.33,
      str_detect(followup_raw, "day")   ~
        as.numeric(str_extract(followup_raw, "[0-9]+\\.?[0-9]*")) / 30.4,
      TRUE ~ NA_real_
    ),
    
    # --- Proportion female ---
    prop_female = case_when(
      !is.na(n_female) & !is.na(n_total) & n_total > 0 ~ n_female / n_total,
      TRUE ~ NA_real_
    ),
    
    # --- Sample size (log-transformed for regression) ---
    log_n = log(n_total),
    
    # --- Region (broad) ---
    country_clean = str_trim(as.character(Country)),
    region = case_when(
      country_clean %in% c("United States", "USA", "US", "Canada") ~ "North America",
      country_clean %in% c("UK", "England", "Sweden", "Norway", "Finland",
                           "Netherlands", "Belgium", "Germany", "Italy",
                           "Iceland", "Ireland", "Spain", "Estonia",
                           "Hungary", "Lithuania") ~ "Europe",
      country_clean %in% c("China", "Japan", "South Korea") ~ "East Asia",
      country_clean %in% c("Australia") ~ "Oceania",
      country_clean %in% c("Latin America", "Vietnam") ~ "Other",
      TRUE ~ "Other"
    ),
    
    # --- Direction as factor for subgroup analysis ---
    direction_factor = factor(reverse, levels = c(0, 1),
                              labels = c("SM -> MH", "MH -> SM")),
    
    # --- Study label for forest plots ---
    study_label = df$Studies
  )

cat("Moderator summary:\n")
cat("  Follow-up months — range:", range(df$followup_months, na.rm = TRUE), "\n")
cat("  Age mean         — range:", range(df$age_mean, na.rm = TRUE), "\n")
cat("  Sample size      — range:", range(df$n_total, na.rm = TRUE), "\n")
cat("  Region:\n")
print(table(df$region))
cat("\n")


# =============================================================================
# SECTION 6: ASSIGN CLUSTER ID FOR MULTILEVEL / ROBUST SE
#
# Multiple effect sizes from the same study share a non-independence structure.
# We use "Reference" as the cluster variable.
# =============================================================================

df <- df %>%
  mutate(
    cluster_id = as.integer(factor(Reference)),
    es_id      = row_number()   # unique effect-size ID for level-2 in rma.mv
  )

cat("Number of unique studies (clusters):", n_distinct(df$cluster_id), "\n")
cat("Number of effect sizes:", nrow(df), "\n\n")
# =============================================================================
# SECTION 7: RANDOM-EFFECTS BASE MODEL (multilevel, 3-level)
#
# 3-level model: effect sizes nested within studies.
# struct = "CS" (compound symmetry) for within-study correlation.
# =============================================================================

cat("======================================================\n")
cat("MODEL 1: Three-level random-effects (base model)\n")
cat("======================================================\n")

res_base <- rma.mv(
  yi = d_unified,
  V  = vi,
  random = ~ 1 | cluster_id / es_id,   # study / effect-size
  data   = df,
  method = "REML"
)
print(summary(res_base))

# Variance partition
cat("\n--- Variance partition ---\n")
I2_total <- sum(res_base$sigma2) / (sum(res_base$sigma2) + res_base$tau2 + mean(df$vi))
cat(sprintf("Level-2 (within-study)  sigma2: %.4f\n", res_base$sigma2[2]))
cat(sprintf("Level-3 (between-study) sigma2: %.4f\n", res_base$sigma2[1]))
cat(sprintf("Overall heterogeneity I2 approx: %.1f%%\n", I2_total * 100))


# =============================================================================
# Model 2 and 3: Direction subgroup model
# =============================================================================
df_reverse <- df[df$reverse == 1, ]
res_direction <- metafor::rma.mv(
  yi = d_unified,
  V  = vi,
  random = ~ 1 | cluster_id / es_id,
  data   = df_reverse,
  method = "REML"
)
print(summary(res_direction))

# Robust SEs for direction sensitivity model
robust_res_direction <- clubSandwich::coef_test(
  res_direction,
  vcov = "CR2",
  cluster = df_reverse$cluster_id
)

cat("\n--- Robust SEs: reverse direction model \n")
print(robust_res_direction)

##forward direction model
df_forward <- df[df$reverse != 1, ]
forward_direction <- metafor::rma.mv(
  yi = d_unified,
  V  = vi,
  random = ~ 1 | cluster_id / es_id,
  data   = df_forward,
  method = "REML"
)
print(summary(forward_direction))

# Robust SEs for direction sensitivity model
robust_forward_direction <- clubSandwich::coef_test(
  forward_direction,
  vcov = "CR2",
  cluster = df_forward$cluster_id
)

cat("\n--- Robust SEs: forward direction model \n")
print(robust_forward_direction)
# =============================================================================
# REVISED FULL META-REGRESSION
# Moderators: direction, follow-up duration, mean age, proportion female
# Removed: region and log sample size
# =============================================================================

cat("\n======================================================\n")
cat("REVISED MODEL 3: Meta-regression without region and sample size\n")
cat("======================================================\n")

# Select rows with complete data for the retained moderators
df_complete_revised <- df %>%
  filter(
    !is.na(followup_months),
    !is.na(age_mean),
    !is.na(prop_female)
  ) %>%
  mutate(es_id = row_number())

cat("N (complete moderators):", nrow(df_complete_revised),
    "from", n_distinct(df_complete_revised$cluster_id), "studies\n\n")

res_full_revised <- rma.mv(
  yi = d_unified,
  V  = vi,
  mods = ~ direction_factor + followup_months + age_mean + prop_female,
  random = ~ 1 | cluster_id / es_id,
  data   = df_complete_revised,
  method = "REML"
)

summary(res_full_revised)

# =============================================================================
# SECTION 11: PARSIMONIOUS FOLLOW-UP MODELS
# =============================================================================

cat("\n======================================================\n")
cat("MODEL 4: Continuous moderators only (larger sample)\n")
cat("======================================================\n")

df_cont <- df %>%
  filter(!is.na(followup_months), !is.na(age_mean)) %>%
  mutate(es_id = row_number())

res_cont <- rma.mv(
  yi = d_unified,
  V  = vi,
  mods = ~ direction_factor + followup_months + age_mean,
  random = ~ 1 | cluster_id / es_id,
  data   = df_cont,
  method = "REML"
)
print(summary(res_cont))


# =============================================================================
# SECTION 12: PUBLICATION BIAS ASSESSMENT
# =============================================================================

cat("\n======================================================\n")
cat("PUBLICATION BIAS\n")
cat("======================================================\n")

# Egger's test via weighted regression of d on SE
egger <- lm(d_unified ~ sei, data = df, weights = 1 / vi)
cat("Egger's test (intercept ≠ 0 suggests asymmetry):\n")
print(summary(egger)$coefficients)

# Rank-correlation test (Begg & Mazumdar)
ranktest_res <- ranktest(res_base)
cat("\nBegg & Mazumdar rank correlation:\n")
print(ranktest_res)

# Trim-and-fill (on base model — note: less reliable with multilevel data)
# Use a simpler random-effects model for this
res_re_simple <- rma(yi = d_unified, vi = vi, data = df, method = "REML")
tf <- trimfill(res_re_simple)
cat("\nTrim-and-fill adjusted estimate:\n")
print(tf)


# =============================================================================
# SECTION 13: VISUALISATIONS
# =============================================================================

# ---------- 13A: Forest plot (forward-direction effects only) ----------
df_forward <- df %>%
  dplyr::filter(reverse == 0) %>%
  dplyr::arrange(d_unified) %>%
  dplyr::mutate(
    study_label_clean = stringr::str_replace_all(study_label, "\\s+", " "),
    study_label_clean = stringr::str_trunc(study_label_clean, width = 35)
  )

k <- nrow(df_forward)

# Larger row spacing
row_spacing <- 1.65
row_pos <- seq(from = k * row_spacing, by = -row_spacing, length.out = k)

# Save as PNG with larger physical size
png(
  filename = "forest_plot_forward_larger_font.png",
  width = 2600,
  height = max(4200, k * 95),
  res = 300
)

par(
  mar = c(5.5, 4.5, 4.5, 2.5),
  cex = 1
)

forward_forest <- metafor::forest(
  x        = df_forward$d_unified,
  vi       = df_forward$vi,
  slab     = df_forward$study_label_clean,
  rows     = row_pos,
  ylim     = c(0, max(row_pos) + 5),
  
  # Whole plot region: study labels + forest + estimate column
  xlim     = c(-2.0, 1.45),
  
  # Effect-size plotting region
  alim     = c(-1, 1),
  at       = c(-1, -0.5, 0, 0.5, 1),
  
  main     = "Forest Plot: Social Media Use \u2192 Mental Health",
  xlab     = "Cohen's d",
  header   = c("Study", "Estimate [95% CI]"),
  
  refline  = 0,
  col      = "black",
  cex      = 0.72,
  psize    = 0.9,
  efac     = 0.7,
  digits   = 2
)

dev.off()

cat("Saved: forest_plot_forward_larger_font.png\n")

# ---------- 13B: Forest plot (reverse-direction effects) ----------
  df_reverse <- df %>%
    filter(reverse == 1) %>%
    arrange(d_unified)
  
  if (nrow(df_reverse) > 0) {
    k <- nrow(df_reverse)
    row_spacing <- 1.5
    
    row_pos <- seq(from = k * row_spacing, by = -row_spacing, length.out = k)
    
    # png("forest_plot_reverse.png",
    #     width  = 750,                    # narrower canvas -> columns sit closer
    #     height = max(500, k * 35),
    #     res    = 120)
    svg("forest_plot_reverse.svg", width = 7.5, height = max(6, k * 0.35))
    
    par(mar = c(5, 4, 4, 2))
    
    reverse_forest<- forest(
      x        = df_reverse$d_unified,
      vi       = df_reverse$vi,
      slab     = df_reverse$study_label,
      rows     = row_pos,
      ylim     = c(0, max(row_pos) + 3),
      xlim     = c(-3, 2.5),               # total coordinate width (slab on left, numbers on right)
      alim     = c(-1, 1),               # range of the CI plotting region
      at       = c(-1, -0.5, 0, 0.5, 1), # x-axis tick marks inside alim
      main     = "Forest Plot: Mental Health \u2192 Social Media Use",
      xlab     = "Cohen's d (sign-flipped to SM\u2192MH direction)",
      col      = "black",
      refline  = 0,
      cex      = 0.75
    )
    
    dev.off()
  }

# ---------- 13C: Funnel plot ----------

# Funnel plot for simple random-effects model

pdf("funnel_plot_simple_RE.pdf", width = 7, height = 6)

metafor::funnel(
  res_re_simple,
  xlab    = "Cohen's d",
  ylab    = "Standard Error",
  main    = "Contour-Enhanced Funnel Plot",
  refline = coef(res_re_simple)[1],
  level   = c(90, 95, 99),
  shade   = c("white", "gray85", "gray70"),
  back    = "white",
  pch     = 21,
  bg      = "gray60",
  col     = "black",
  cex     = 0.9
)

legend(
  "topright",
  legend = c("p > .10", ".05 < p ≤ .10", ".01 < p ≤ .05", "p ≤ .01"),
  fill   = c("white", "gray85", "gray70", "gray55"),
  bty    = "n",
  cex    = 0.8
)

dev.off()

cat("Saved: funnel_plot_simple_RE.pdf\n")

# =============================================================================
# CREATE AND COMBINE TWO BUBBLE PLOTS USING PATCHWORK
# =============================================================================

# Load packages
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")

library(ggplot2)
library(dplyr)
library(patchwork)

# =============================================================================
# 1. Bubble plot: effect size vs. follow-up duration
# =============================================================================

df_plot <- df %>%
  dplyr::filter(!is.na(followup_months))

p_bubble <- ggplot(df_plot, aes(
  x = followup_months,
  y = d_unified,
  size = 1 / vi,
  colour = direction_label
)) +
  geom_point(alpha = 0.55) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_smooth(
    aes(weight = 1 / vi, group = direction_label),
    method = "lm",
    se = TRUE,
    linewidth = 0.8
  ) +
  scale_colour_manual(values = c("steelblue", "tomato")) +
  scale_size_continuous(range = c(1.5, 8), guide = "none") +
  labs(
    title = "Effect Size vs. Follow-Up Duration",
    x = "Follow-up duration (months)",
    y = "Cohen's d",
    colour = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# =============================================================================
# 2. Bubble plot: effect size vs. mean age
# =============================================================================

df_plot2 <- df %>%
  dplyr::filter(!is.na(age_mean))

p_age <- ggplot(df_plot2, aes(
  x = age_mean,
  y = d_unified,
  size = 1 / vi,
  colour = direction_label
)) +
  geom_point(alpha = 0.55) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_smooth(
    aes(weight = 1 / vi, group = direction_label),
    method = "lm",
    se = TRUE,
    linewidth = 0.8
  ) +
  scale_colour_manual(values = c("steelblue", "tomato")) +
  scale_size_continuous(range = c(1.5, 8), guide = "none") +
  labs(
    title = "Effect Size vs. Mean Age",
    x = "Mean age (years)",
    y = "Cohen's d",
    colour = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# =============================================================================
# 3. Combine the two plots
# =============================================================================

combined_bubble <- (p_bubble | p_age) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1)
  )

# Show combined figure
combined_bubble

# =============================================================================
# 4. Save combined figure
# =============================================================================

ggsave(
  filename = "combined_bubble_plots.png",
  plot = combined_bubble,
  width = 14,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "combined_bubble_plots.pdf",
  plot = combined_bubble,
  width = 14,
  height = 6
)

cat("Saved: combined_bubble_plots.png and combined_bubble_plots.pdf\n")
# =============================================================================
# SECTION 14: RESULTS SUMMARY TABLE
# =============================================================================

results_table <- data.frame(
  Model = c(
    "Base 3-level RE",
    "Direction subgroup: SM -> MH",
    "Direction subgroup: MH -> SM",
    "Revised meta-regression",
    "Trim-and-fill adjusted"
  ),
  k = c(
    nrow(df),
    sum(df$reverse == 0),
    sum(df$reverse == 1),
    nrow(df_complete_revised),
    nrow(df)
  ),
  Estimate_d = round(c(
    as.numeric(res_base$b[1]),
    as.numeric(res_direction$b[1]),
    as.numeric(res_direction$b[1]) + as.numeric(res_direction$b[2]),
    as.numeric(res_full_revised$b[1]),
    as.numeric(tf$b[1])
  ), 4),
  CI_lower = round(c(
    res_base$ci.lb[1],
    res_direction$ci.lb[1],
    NA,
    res_full_revised$ci.lb[1],
    tf$ci.lb
  ), 4),
  CI_upper = round(c(
    res_base$ci.ub[1],
    res_direction$ci.ub[1],
    NA,
    res_full_revised$ci.ub[1],
    tf$ci.ub
  ), 4),
  p_value = round(c(
    res_base$pval[1],
    res_direction$pval[1],
    NA,
    res_full_revised$pval[1],
    tf$pval
  ), 4)
)

write.csv(results_table, "meta_regression_results.csv", row.names = FALSE)

cat("\nSaved: meta_regression_results.csv\n")
cat("\n=== RESULTS SUMMARY ===\n")
print(results_table)


# =============================================================================
# SECTION 15: SENSITIVITY — REMOVE ZERO EFFECT SIZES
# (Zero d values were coded for non-significant unreported effects)
# =============================================================================

cat("\n======================================================\n")
cat("SENSITIVITY: Excluding zero-coded non-significant effects\n")
cat("======================================================\n")

df_nonzero <- df %>% filter(d_unified != 0) %>% mutate(es_id = row_number())
cat("N after removing zeros:", nrow(df_nonzero), "\n")

res_nonzero <- rma.mv(
  yi = d_unified,
  V  = vi,
  random = ~ 1 | cluster_id / es_id,
  data   = df_nonzero,
  method = "REML"
)
cat(sprintf(
  "Pooled d = %.3f [%.3f, %.3f], p = %.4f\n",
  res_nonzero$b[1], res_nonzero$ci.lb, res_nonzero$ci.ub, res_nonzero$pval[1]
))

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Output files created in your working directory:\n")
cat("  forest_plot_forward.png\n")
cat("  forest_plot_reverse.png\n")
cat("  funnel_plot.png\n")
cat("  bubble_followup.png\n")
cat("  bubble_age.png\n")
cat("  meta_regression_results.csv\n")

#Year
# Prepare year counts
plot_df_year <- df %>%
  distinct(Studies, .keep_all = TRUE) %>%   # one row per longitudinal study
  mutate(Year = as.numeric(Year)) %>%
  filter(!is.na(Year)) %>%
  count(Year, name = "Count") %>%
  complete(Year = 2015:2023, fill = list(Count = 0))

# Plot
p_a <- ggplot(plot_df_year, aes(x = Year, y = Count)) +
  geom_col(
    fill = "grey50",
    color = "grey50",
    width = 0.9
  ) +
  geom_text(
    aes(label = Count),
    vjust = -0.3,
    size = 4
  ) +
  scale_x_continuous(
    breaks = 2015:2023,
    limits = c(2014.5, 2023.5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, max(plot_df_year$Count, na.rm = TRUE) + 1, 1),
    limits = c(0, max(plot_df_year$Count, na.rm = TRUE) + 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Year longitudinal study was made available",
    y = "Number of longitudinal studies in that year",
    title = "Year of longitudinal studies"
  ) +
  theme_gray(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 11),   # show year labels
    axis.ticks.x = element_line(),           # show x-axis ticks
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
p_a
library(ggplot2)
library(dplyr)

# Age included
plot_df_age <- df %>%
  distinct(Studies, .keep_all = TRUE) %>%   # one row per study
  mutate(
    `Minimum age` = as.numeric(`Minimum.age`),
    `Maximum age` = as.numeric(`Maximum.age`),
    `Mean age` = as.numeric(Age.mean)
  ) %>%
  filter(!is.na(`Minimum age`), !is.na(`Maximum age`)) %>%
  arrange(`Minimum age`, `Maximum age`) %>%
  mutate(
    Studies = factor(Studies, levels = Studies)
  )

p_b <- ggplot(plot_df_age, aes(y = Studies)) +
  geom_segment(
    aes(
      x = `Minimum age`,
      xend = `Maximum age`,
      yend = Studies
    ),
    linewidth = 0.4,
    color = "grey20"
  ) +
  geom_point(
    aes(x = `Minimum age`),
    size = 1,
    color = "grey20"
  ) +
  geom_point(
    aes(x = `Maximum age`),
    size = 1,
    color = "grey20"
  ) +
  geom_point(
    aes(x = `Mean age`),
    size = 2.2,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.8,
    na.rm = TRUE
  ) +
  scale_x_continuous(
    limits = c(7, 21),
    breaks = seq(7, 21, by = 1),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Age range of included participants",
    y = "Study",
    title = "Age range and mean age of participants in each study"
  ) +
  theme_gray(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 11),
    axis.ticks.x = element_line(),
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p_b

#Violin plot of sample size
plot_df_sample <- df %>%
  distinct(Studies, .keep_all = TRUE) %>%   # one row per meta-analysis
  mutate(`Sample size` = as.numeric(`Sample.size`)) %>%
  filter(!is.na(`Sample size`))

p_c <- ggplot(plot_df_sample, aes(x = "Density", y = `Sample size`)) +
  geom_violin(
    trim = FALSE,
    fill = "white",
    color = "grey40",
    linewidth = 0.5,
    width = 0.9
  ) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    color = "grey30",
    linewidth = 0.5,
    outlier.shape = 16,
    outlier.size = 2.5,
    outlier.color = "grey35"
  ) +
  scale_y_continuous(
    limits = c(0, max(plot_df_sample$`Sample size`, na.rm = TRUE) * 1.05),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = "Density",
    y = "Sample size in longitudinal studies",
    title = "Violin and box plot of sample size in longitudinal studies"
  ) +
  theme_gray(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
p_c

library(ggplot2)
library(dplyr)

# Percentage of female participants
plot_df_female <- df %>%
  distinct(Studies, .keep_all = TRUE) %>%
  mutate(
    `Sample size` = as.numeric(`Sample.size`),
    `Female number` = as.numeric(`Female.number`),
    Female_percent = (`Female number` / `Sample size`) * 100
  ) %>%
  filter(!is.na(Female_percent)) %>%
  arrange(Female_percent) %>%
  mutate(
    Studies = factor(Studies, levels = Studies)
  )

p_d <- ggplot(plot_df_female, aes(x = Female_percent, y = Studies)) +
  geom_col(
    fill = "grey70",
    color = "black",
    linewidth = 0.25,
    width = 0.75
  ) +
  geom_text(
    aes(label = paste0(round(Female_percent, 1), "%")),
    hjust = -0.15,
    size = 3
  ) +
  scale_x_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, by = 10),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Percentage of female participants",
    y = "Study",
    title = "Percentage of female participants in each study"
  ) +
  theme_gray(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 8),
    axis.ticks.x = element_line(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p_d
#combine them together
install.packages("patchwork")  # only if not installed
library(patchwork)
combined_plot <- (p_a | p_b) / (p_c | p_d) +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0, 1)
  )

combined_plot

# Make sure direction variable is correctly coded
df <- df %>%
  mutate(
    direction_factor = factor(
      reverse,
      levels = c(0, 1),
      labels = c("SM -> MH", "MH -> SM")
    )
  )

# Check sample sizes by direction
table(df$direction_factor)

df_forward <- df %>%
  filter(direction_factor == "SM -> MH") %>%
  mutate(
    es_id = row_number(),
    cluster_id = as.integer(factor(Reference))
  )

res_forward <- rma.mv(
  yi = d_unified,
  V  = vi,
  random = ~ 1 | cluster_id / es_id,
  data = df_forward,
  method = "REML"
)

summary(res_forward)

# Cluster-robust SEs
robust_forward <- coef_test(
  res_forward,
  vcov = "CR2",
  cluster = df_forward$cluster_id
)

robust_forward

df_reverse <- df %>%
  filter(direction_factor == "MH -> SM") %>%
  mutate(
    es_id = row_number(),
    cluster_id = as.integer(factor(Reference))
  )

res_reverse <- rma.mv(
  yi = d_unified,
  V  = vi,
  random = ~ 1 | cluster_id / es_id,
  data = df_reverse,
  method = "REML"
)

summary(res_reverse)

# Cluster-robust SEs
robust_reverse <- coef_test(
  res_reverse,
  vcov = "CR2",
  cluster = df_reverse$cluster_id
)

robust_reverse

subgroup_results <- data.frame(
  Direction = c("SM -> MH", "MH -> SM"),
  k = c(nrow(df_forward), nrow(df_reverse)),
  Studies = c(
    dplyr::n_distinct(df_forward$cluster_id),
    dplyr::n_distinct(df_reverse$cluster_id)
  ),
  Estimate = c(coef(res_forward)[1], coef(res_reverse)[1]),
  SE = c(res_forward$se[1], res_reverse$se[1]),
  CI_lower = c(res_forward$ci.lb[1], res_reverse$ci.lb[1]),
  CI_upper = c(res_forward$ci.ub[1], res_reverse$ci.ub[1]),
  p_value = c(res_forward$pval[1], res_reverse$pval[1])
)

subgroup_results

subgroup_results_rounded <- subgroup_results %>%
  mutate(
    Estimate = round(Estimate, 3),
    SE = round(SE, 3),
    CI_lower = round(CI_lower, 3),
    CI_upper = round(CI_upper, 3),
    p_value = round(p_value, 3)
  )

subgroup_results_rounded