########################################################################################################
# Social Media Use Meta-meta-analysis with Jungup
########################################################################################################
# Authors: Qiyang Zhang
# Contact: qiyang39@nus.edu.sg
# Created: 2025/05/14
# Revised: 2026/3/7

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
#robumeta
library(robumeta)
library(metafor)
library(dplyr)

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
forest(
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
library(ggplot2)
library(dplyr)
library(patchwork)
library(tibble)

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


