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
full <- merge(studies, findings, by = c("Study"), all = TRUE, suffixes = c(".s", ".f"))

# format to correct variable types
nums <- c("Sample size", "Effect.Size", "Mean.age", "Female.percent",
          "Above16", "Above15", "FiftyFivePercentFemale")

full[nums] <- lapply(full[nums], as.numeric)
rm(nums)
full$Effect.Size[which(full$reverse==1)]<- full$Effect.Size * -1

###############################################################
#Create unique identifiers (ES, study, program)
###############################################################
full$ESId <- as.numeric(rownames(full))
full$StudyID <- as.numeric(as.factor(full$Study))
summary(full$StudyID)

full$Effect.Size <- as.numeric(full$Effect.Size)

################################################################
# Calculate meta-analytic variables: Variances (Lipsey & Wilson, 2000, Eq. 3.23)
################################################################
#calculate standard errors
full$se<-sqrt(((full$`Sample size`)/((full$`Sample size`/2)*(full$`Sample size`/2)))+((full$Effect.Size*full$Effect.Size)/(2*(full$`Sample size`))))

#calculate variance
full$var<-full$se*full$se

###########################
#forest plot
###########################
study_averages <- full %>%
  dplyr::group_by(StudyID) %>%
  dplyr::summarise(avg_effect_size = mean(Effect.Size, na.rm = TRUE),
            across(everything(), ~ first(.x)),
            .groups = "drop")

MVnull <- robu(formula = Effect.Size ~ 1, studynum = StudyID, data = study_averages, var.eff.size = var)
m.gen <- metagen(TE = Effect.Size,
                 seTE = se,
                 studlab = Study,
                 data = study_averages,
                 sm = "SMD",
                 fixed = FALSE,
                 random = TRUE,
                 method.tau = "REML",
                 method.random.ci = "HK")
summary(m.gen)
png(file = "forestplot.png", width = 2800, height = 3000, res = 300)

meta::forest(m.gen,
             sortvar = TE,
             prediction = TRUE,
             print.tau2 = FALSE,
             leftlabs = c("Author", "g", "SE"))
dev.off()
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

### Output prediction interval ###
# dat_clean <- data.frame(yi = full$Effect.Size, se_g = full$se)
# dat_clean <- na.omit(dat_clean)
# 
# yi_clean <- dat_clean$yi
# se_g_clean <- dat_clean$se_g
# install.packages("pimeta")
# library(pimeta)
# 
# pima_result <- pima(yi_clean, se_g_clean, method = "HK")  # Using the Hartung-Knapp method
# 
# print(pima_result)

#moderator analysis
full <- full %>%
  dplyr::filter(Study != "Fassi et al. 2024")

full$Above16.c <- full$Above16 - mean(full$Above16)
full$Above15.c <- full$Above15 - mean(full$Above15)
full$FiftyFivePercentFemale.c <- full$FiftyFivePercentFemale - mean(full$FiftyFivePercentFemale)

terms <- c("FiftyFivePercentFemale.c")
terms <- c("Above15.c")
formula <- reformulate(termlabels = c(terms))
formula
full$ESId <- as.numeric(rownames(full))
full$StudyID <- as.numeric(as.factor(full$Study))
summary(full$StudyID)
V_list <- impute_covariance_matrix(vi=full$var, cluster=full$StudyID, r=0.8)

MVfull <- rma.mv(yi=Effect.Size,
                 V=V_list,
                 mods=formula,
                 random=~1 | StudyID/ESId,
                 test="t",
                 data=full,
                 method="REML")
MVfull

#t-test of each covariate#
MVfull.coef <- coef_test(MVfull, cluster=full$StudyID, vcov="CR2")
MVfull.coef

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
                      data=full, #define data
                      method="REML") #estimate variances using REML
  coef_mod_means <- as.data.frame(coef_test(mod_means,#estimation model above
                                            cluster=full$StudyID, #define cluster IDs
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
MVfull.modeling <- rma(yi=Effect.Size,
                       vi=var,
                       test="t",
                       data=full,
                       slab = Study,
                       method="REML")
metafor::funnel(MVfull.modeling) 

sel <- selmodel(MVfull.modeling, type = "stepfun", steps = 0.025)
plot(sel, ylim=c(0,5))

pdf(file="countour_funnel_color.pdf")
# contour-enhanced funnel plot
metafor::funnel(MVfull.modeling, level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), refline=0, legend=FALSE)
dev.off()



# heatmap
library(ggplot2)
library(dplyr)
library(forcats)

# Ensure variables are factors
full <- full %>%
  mutate(
    Y = as.factor(Y),
    Study = as.factor(Study)
  )

# Aggregate by Study × Outcome
by_study_outcome <- full %>%
  group_by(Study, Y) %>%
  summarise(
    Year = first(as.integer(Year.s)),   # use your correct year variable
    Effect.Size = mean(Effect.Size, na.rm = TRUE),
    n_rows = n(),
    .groups = "drop"
  )

# Order studies from NEWER → OLDER
study_levels <- by_study_outcome %>%
  distinct(Study, Year) %>%
  arrange(desc(Year), Study) %>%
  pull(Study)

# Apply ordering
by_study_outcome <- by_study_outcome %>%
  mutate(
    Study = factor(Study, levels = study_levels)
  )

# Plot heatmap
ggplot(by_study_outcome, aes(x = Y, y = Study, fill = Effect.Size)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#C53030",
    mid = "#F2F2F2",
    high = "#2B6CB0",
    midpoint = 0,
    name = "Effect Size"
  ) +
  labs(
    x = "Mental Health Problems",
    y = "Meta-analysis",
    title = "Heatmap of Effect Sizes by Outcome Category and Meta-analysis"
  ) +
  geom_text(aes(label = sprintf("%.2f", Effect.Size)), size = 5) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


