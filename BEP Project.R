## Loading Packages ----

#install.packages("foreign")
library(foreign)
# package for running SEM
#install.packages("lavaan")
library(lavaan)
# package for plotting SEM
#install.packages('semPlot')
library(semPlot)
#install.packages('ggplot2')
library(ggplot2)
#install.packages('htmlTable')
library(htmlTable)
#install.packages('dplyr')
library(dplyr)
#install.packages('stringr')
library(stringr)

## Opening and renaming dataset ----

setwd('C:/Users/nimah/OneDrive/Bachelor of Psychology')
BDATA = read.spss("DataBoneMetastases.sav", to.data.frame=TRUE)
saveRDS(BDATA, "BDATA.rds")

# Only run this line after saving:
BDATA <- readRDS("BDATA.rds")

## Data inspection ----

View(BDATA)

print(colnames(BDATA), max = ncol(BDATA))
str(BDATA, list.len = ncol(BDATA))

options(max.print = 100000)
Filter(Negate(is.null), lapply(BDATA, levels))
Filter(Negate(is.null), lapply(BDATA, function(x) if(is.factor(x)) length(levels(x)) else NULL))

data.frame(index = 1:ncol(BDATA), name = colnames(BDATA))
BDATA <- BDATA[, 1:1734]

## Data inspection: removing 'onbekend' and valgevoel ----

BDATA <- BDATA %>%
  select(-contains("valgvoel"))

# Columns not to touch
protected_cols <- c("idnr", "condition", "age", "gender", "lastmeas", "numbermeas")

BDATA_cleaned <- BDATA 

# Loop through and only remove 'onbekend'
for (col_name in names(BDATA_cleaned)) {
  if (!(col_name %in% protected_cols) && is.factor(BDATA_cleaned[[col_name]])) {   # Only touch non-protected factor columns
    levs <- levels(BDATA_cleaned[[col_name]]) # Get the current levels
    target <- which(str_trim(levs) == "onbekend") # Find the level that is 'onbekend'
    if (length(target) > 0) { # It removes the level and converts the data to NA without touching other levels
      is.na(levels(BDATA_cleaned[[col_name]]))[target] <- TRUE
    }
  }
}

# Check if it worked...
Filter(Negate(is.null), lapply(BDATA_cleaned, levels))
Filter(Negate(is.null), lapply(BDATA_cleaned, function(x) if(is.factor(x)) length(levels(x)) else NULL))

## Data selection ----

# Baseline -> 1.00
# Immediate post-treatment: week 1 -> 2.00
# Stabilization post-treatment: week 4 -> 5.00

# Again, columns not to touch
protected_cols <- c("idnr", "condition", "age", "gender", "lastmeas", "numbermeas")

# Select protected columns + any column ending in 1.00, 2.00, or 5.00
BDATA_slim <- BDATA_cleaned %>%
  select(
    all_of(protected_cols),
    matches("\\.1\\.00$|\\.2\\.00$|\\.5\\.00$")
  )

# Check if it worked...
options(max.print = 100000)
Filter(Negate(is.null), lapply(BDATA_slim, levels))
Filter(Negate(is.null), lapply(BDATA_slim, function(x) if(is.factor(x)) length(levels(x)) else NULL))

# Check how each factor is represented within the columns...
lapply(BDATA_slim, table, useNA = "always")

saveRDS(BDATA_slim, "BDATA_slim.rds") # Saving to continue later...

setwd('C:/Users/nimah/OneDrive/Bachelor of Psychology')
BDATA_slim <- readRDS("BDATA_slim.rds")

## Data inspection: removing irrelevant items ----

data.frame(index = 1:ncol(BDATA_slim), name = colnames(BDATA_slim))
BDATA_slim <- BDATA_slim[, 1:132]

## Examining direction of categories and whether data is ready for CFA ----

protected_cols <- c("idnr", "condition", "age", "gender", "lastmeas", "numbermeas")

# Convert to numeric using the underlying factor levels
BDATA_numeric <- BDATA_slim %>%
  mutate(across(-all_of(protected_cols), ~as.numeric(.)))

# Identify which items to be reversed, so that the intensity-scores mean the same across all items (otherwise they could cancel each other out)
BDATA_numeric <- BDATA_numeric %>%
  mutate(across(112:ncol(.), ~ (5 - .)))

## CFA as a safety check for parceling ----

item_cfa_model1 <- '
  # Psychological distresss
  psych_dist =~ vangst.1.00 + vgespann.1.00 + vzenuwen.1.00 + vneersla.1.00 + vwanhoop.1.00 + vprikkel.1.00 + vpieker.1.00
  
  # Basic mobility
  bas_act =~ vlopenbi.1.00 + vlopenbu.1.00 + vtraplop.1.00 + vverzorg.1.00
  
  # Complex mobility
  comp_act =~ vklusjes.1.00 + vhuiswer.1.00 + vboodsch.1.00 
  
  # Somatic distress
  somatic_anx =~ vkortade.1.00 + vduizeli.1.00 + vtinteli.1.00 + vrilleri.1.00 + vhoofdpy.1.00
  
  # Musculoskeletal pain
  musculo =~ vspierpy.1.00 + vpijnrug.1.00 + vpijnbot.1.00 
  
  # Cutaneous and barrier abnormalities
  cutane =~ vpijnhui.1.00 + vbranoog.1.00  + vhaaruit.1.00 + vjeuk.1.00 + + vmondsli.1.00 + vmonddro.1.00
  
  # Gastrointestinal distress
  gast_dist =~ vmissely.1.00 + vbraken.1.00 + vmaagzuu.1.00 + vetennee.1.00 + vbuikpyn.1.00 + vverstop.1.00 + vdiarree.1.00
  
  # Vitality
  vita =~ vmoeheid.1.00 + vfutloos.1.00 + vslapelo.1.00 + vsexverm.1.00 + vconmoei.1.00
  
  # Social functioning
  social_func =~ vfamilie.1.00 + vsociale.1.00
'

fit_item_cfa1 <- cfa(item_cfa_model1, 
                     data = BDATA_numeric, 
                     ordered = TRUE,
                     estimator = "WLSMV",
                     missing = "pairwise",
                     std.lv = TRUE)

summary(fit_item_cfa1, fit.measures = TRUE, standardized = TRUE)

lavInspect(fit_item_cfa1, "cov.lv") # no problematic observation...

modindic <- modificationIndices(fit_item_cfa1)
modindic[order(modindic$mi,decreasing=TRUE),][1:20,] 


item_cfa_model2 <- '
  # Psychological distresss
  psych_dist =~ vangst.1.00 + vgespann.1.00 + vzenuwen.1.00 + vneersla.1.00 + vwanhoop.1.00 + vprikkel.1.00 + vpieker.1.00 + vsexverm.1.00 
  
  # Basic mobility
  bas_act =~ vlopenbi.1.00 + vlopenbu.1.00 + vtraplop.1.00 + vverzorg.1.00
  
  # Complex activity
  comp_act =~ vklusjes.1.00 + vhuiswer.1.00 + vboodsch.1.00 
  
  # Gastrointestinal distress
  gast_dist =~ vmissely.1.00 + vbraken.1.00 + vmaagzuu.1.00 + vetennee.1.00 + vdiarree.1.00 + vduizeli.1.00 
  
  # Pain Experience
  pain_exp =~ vspierpy.1.00 + vpijnrug.1.00 + vbuikpyn.1.00 + vverstop.1.00 + vrilleri.1.00 + vhoofdpy.1.00 + vpijnbot.1.00 + vpijnhui.1.00
  
  # Chemotherapy Related
  chemo_rel =~ vtinteli.1.00 + vmondsli.1.00 + vhaaruit.1.00 +  vjeuk.1.00 + vbranoog.1.00   
  
  # Fatigue
  vita =~ vmoeheid.1.00 + vfutloos.1.00 + vslapelo.1.00 + vconmoei.1.00 + vkortade.1.00 + vmonddro.1.00
  
  # Social functioning
  social_func =~ vfamilie.1.00 + vsociale.1.00
'

fit_item_cfa2 <- cfa(item_cfa_model2, 
                     data = BDATA_numeric, 
                     ordered = TRUE,
                     estimator = "WLSMV",
                     missing = "pairwise",
                     std.lv = TRUE)

summary(fit_item_cfa2, fit.measures = TRUE, standardized = TRUE)

anova(fit_item_cfa1, fit_item_cfa2) # deviations seem to be better...


## Parceling ----

# Parcelling Psychological distress at three time-points

# Time 1.00
BDATA_numeric$psych_1 <- rowMeans(BDATA_numeric[, c("vangst.1.00", "vgespann.1.00", "vzenuwen.1.00", 
                                                    "vneersla.1.00", "vwanhoop.1.00", "vprikkel.1.00", 
                                                    "vpieker.1.00")], na.rm = TRUE)
# Time 2.00
BDATA_numeric$psych_2 <- rowMeans(BDATA_numeric[, c("vangst.2.00", "vgespann.2.00", "vzenuwen.2.00", 
                                                    "vneersla.2.00", "vwanhoop.2.00", "vprikkel.2.00", 
                                                    "vpieker.2.00")], na.rm = TRUE)
# Time 5.00
BDATA_numeric$psych_5 <- rowMeans(BDATA_numeric[, c("vangst.5.00", "vgespann.5.00", "vzenuwen.5.00", 
                                                    "vneersla.5.00", "vwanhoop.5.00", "vprikkel.5.00", 
                                                    "vpieker.5.00")], na.rm = TRUE)

# Parcelling Somatic distress at three time-points 

# Time 1.00
BDATA_numeric$somatic_1 <- rowMeans(BDATA_numeric[, c("vkortade.1.00", "vduizeli.1.00", "vtinteli.1.00",
                                                      "vrilleri.1.00", "vhoofdpy.1.00")], na.rm = TRUE)
# Time 2.00
BDATA_numeric$somatic_2 <- rowMeans(BDATA_numeric[, c("vkortade.2.00", "vduizeli.2.00", "vtinteli.2.00",
                                                      "vrilleri.2.00", "vhoofdpy.2.00")], na.rm = TRUE)
# Time 5.00
BDATA_numeric$somatic_5 <- rowMeans(BDATA_numeric[, c("vkortade.5.00", "vduizeli.5.00", "vtinteli.5.00",
                                                      "vrilleri.5.00", "vhoofdpy.5.00")], na.rm = TRUE)

# Parcelling Gastrointestinal distress at three time-points

# Time 1.00
BDATA_numeric$gast_1 <- rowMeans(BDATA_numeric[, c("vmissely.1.00", "vbraken.1.00", "vmaagzuu.1.00", 
                                                   "vetennee.1.00", "vbuikpyn.1.00", "vverstop.1.00", 
                                                   "vdiarree.1.00")], na.rm = TRUE)
# Time 2.00
BDATA_numeric$gast_2 <- rowMeans(BDATA_numeric[, c("vmissely.2.00", "vbraken.2.00", "vmaagzuu.2.00", 
                                                   "vetennee.2.00", "vbuikpyn.2.00", "vverstop.2.00", 
                                                   "vdiarree.2.00")], na.rm = TRUE)
# Time 5.00
BDATA_numeric$gast_5 <- rowMeans(BDATA_numeric[, c("vmissely.5.00", "vbraken.5.00", "vmaagzuu.5.00", 
                                                   "vetennee.5.00", "vbuikpyn.5.00", "vverstop.5.00", 
                                                   "vdiarree.5.00")], na.rm = TRUE)

# Parcelling Musculoskeletal pain at three time-points

# Time 1.00
BDATA_numeric$musculo_1 <- rowMeans(BDATA_numeric[, c("vspierpy.1.00", "vpijnrug.1.00", 
                                                      "vpijnbot.1.00")], na.rm = TRUE)
# Time 2.00
BDATA_numeric$musculo_2 <- rowMeans(BDATA_numeric[, c("vspierpy.2.00", "vpijnrug.2.00", 
                                                      "vpijnbot.2.00")], na.rm = TRUE)
# Time 5.00
BDATA_numeric$musculo_5 <- rowMeans(BDATA_numeric[, c("vspierpy.5.00", "vpijnrug.5.00", 
                                                      "vpijnbot.5.00")], na.rm = TRUE)

# Parcelling Cutaneous and barrier abnormalities at three time-points

# Time 1.00
BDATA_numeric$cutane_1 <- rowMeans(BDATA_numeric[, c("vpijnhui.1.00", "vbranoog.1.00", "vmonddro.1.00",
                                                     "vhaaruit.1.00", "vjeuk.1.00", "vmondsli.1.00")], na.rm = TRUE)
# Time 2.00
BDATA_numeric$cutane_2 <- rowMeans(BDATA_numeric[, c("vpijnhui.2.00", "vbranoog.2.00", "vmonddro.2.00",
                                                     "vhaaruit.2.00", "vjeuk.2.00", "vmondsli.2.00")], na.rm = TRUE)
# Time 5.00
BDATA_numeric$cutane_5 <- rowMeans(BDATA_numeric[, c("vpijnhui.5.00", "vbranoog.5.00", "vmonddro.5.00",
                                                     "vhaaruit.5.00", "vjeuk.5.00", "vmondsli.5.00")], na.rm = TRUE)

# Parcelling Basic mobility at three time-points

# Time 1.00
BDATA_numeric$bas_act_1 <- rowMeans(BDATA_numeric[, c("vlopenbi.1.00", "vlopenbu.1.00", 
                                                      "vtraplop.1.00", "vverzorg.1.00")], na.rm = TRUE)
# Time 2.00
BDATA_numeric$bas_act_2 <- rowMeans(BDATA_numeric[, c("vlopenbi.2.00", "vlopenbu.2.00", 
                                                      "vtraplop.2.00", "vverzorg.2.00")], na.rm = TRUE)
# Time 5.00
BDATA_numeric$bas_act_5 <- rowMeans(BDATA_numeric[, c("vlopenbi.5.00", "vlopenbu.5.00", 
                                                      "vtraplop.5.00", "vverzorg.5.00")], na.rm = TRUE)

# Parcelling Complex mobility at three time-points

# Time 1.00
BDATA_numeric$comp_act_1 <- rowMeans(BDATA_numeric[, c("vklusjes.1.00", "vhuiswer.1.00", 
                                                       "vboodsch.1.00")], na.rm = TRUE)
# Time 2.00
BDATA_numeric$comp_act_2 <- rowMeans(BDATA_numeric[, c("vklusjes.2.00", "vhuiswer.2.00", 
                                                       "vboodsch.2.00")], na.rm = TRUE)
# Time 5.00
BDATA_numeric$comp_act_5 <- rowMeans(BDATA_numeric[, c("vklusjes.5.00", "vhuiswer.5.00", 
                                                       "vboodsch.5.00")], na.rm = TRUE)

# Parcelling Social functioning at three time-points

# Time 1.00
BDATA_numeric$social_1 <- rowMeans(BDATA_numeric[, c("vfamilie.1.00", "vsociale.1.00")], na.rm = TRUE)
# Time 2.00
BDATA_numeric$social_2 <- rowMeans(BDATA_numeric[, c("vfamilie.2.00", "vsociale.2.00")], na.rm = TRUE)
# Time 5.00
BDATA_numeric$social_5 <- rowMeans(BDATA_numeric[, c("vfamilie.5.00", "vsociale.5.00")], na.rm = TRUE)

# Parcelling Vitality at three time-points

# Time 1.00
BDATA_numeric$vita_1 <- rowMeans(BDATA_numeric[, c("vmoeheid.1.00", "vfutloos.1.00", "vslapelo.1.00", 
                                                   "vsexverm.1.00", "vconmoei.1.00")], na.rm = TRUE)
# Time 2.00
BDATA_numeric$vita_2 <- rowMeans(BDATA_numeric[, c("vmoeheid.2.00", "vfutloos.2.00", "vslapelo.2.00", 
                                                   "vsexverm.2.00", "vconmoei.2.00")], na.rm = TRUE)
# Time 5.00
BDATA_numeric$vita_5 <- rowMeans(BDATA_numeric[, c("vmoeheid.5.00", "vfutloos.5.00", "vslapelo.5.00", 
                                                   "vsexverm.5.00", "vconmoei.5.00")], na.rm = TRUE)

saveRDS(BDATA_numeric, "BDATA_numeric.rds") # This is our new final data from now on, ready to be used


## Examining the parcels (explore a bit) ----

ggplot(BDATA_numeric, aes(x = social_2)) +
  geom_histogram(binwidth = 0.2, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of SOCIAL (T1)", x = "Mean Score", y = "Count") 


## NEXT TIME START FROM HERE: Reloading everything ----

library(foreign)
library(lavaan)
library(semPlot)
library(ggplot2)
library(htmlTable)
library(dplyr)
library(stringr)

setwd('C:/Users/nimah/OneDrive/Bachelor of Psychology')
BDATA_numeric <- readRDS("BDATA_numeric.rds")

options(max.print = 100000)

# Strip ghost levels and keep only 8 and 24 (will run into problems otherwise)
BDATA_numeric <- BDATA_numeric[BDATA_numeric$condition %in% c(8, 24), ]
BDATA_numeric$condition <- factor(BDATA_numeric$condition)

## Running factor analysis ----

T1_MODEL = '
    # Pattern of Factor Loadings
    Psych_1 =~ L11*psych_1 + L21*social_1 + L31*vita_1
    Phys_1  =~ L41*somatic_1 + L51*gast_1 + L61*musculo_1 + L71*cutane_1 + L81*vita_1
    Func_1  =~ L91*bas_act_1 + L101*comp_act_1 + L111*social_1 
    
    # Common Factor Variances (UVI Identification)
    Psych_1 ~~ 1*Psych_1
    Phys_1  ~~ 1*Phys_1
    Func_1  ~~ 1*Func_1

    # Residual Factor Variances
    psych_1    ~~ E11*psych_1
    social_1   ~~ E22*social_1
    vita_1     ~~ E33*vita_1
    somatic_1  ~~ E44*somatic_1
    gast_1     ~~ E55*gast_1
    musculo_1  ~~ E66*musculo_1
    cutane_1   ~~ E77*cutane_1
    bas_act_1  ~~ E88*bas_act_1
    comp_act_1 ~~ E99*comp_act_1

    # Common Factor Means (Fixed to 0 for baseline)
    Psych_1 ~ 0*1
    Phys_1  ~ 0*1
    Func_1  ~ 0*1

    # Intercepts
    psych_1    ~ INT11*1
    social_1   ~ INT21*1
    vita_1     ~ INT31*1
    somatic_1  ~ INT41*1
    gast_1     ~ INT51*1
    musculo_1  ~ INT61*1
    cutane_1   ~ INT71*1
    bas_act_1  ~ INT81*1
    comp_act_1 ~ INT91*1
'    

# Run the model
T1_FIT <- cfa(T1_MODEL, 
              data = BDATA_numeric, 
              group = "condition",
              estimator = "MLR", 
              missing = "fiml",
              std.lv = TRUE)

# plot of the model we fitted
semPaths(T1_FIT)
semPaths(T1_FIT,whatLabels="par")

# View the full statistical output
summary(T1_FIT, fit.measures = TRUE, standardized = TRUE)

# modification indices
mod_T1 <- modificationIndices(T1_FIT)
head(mod_T1[order(mod_T1$mi, decreasing=TRUE), ], 10)


T2_MODEL = '
    # Pattern of Factor Loadings
    Psych_2 =~ L12*psych_2 + L22*social_2 + L32*vita_2
    Phys_2  =~ L42*somatic_2 + L52*gast_2 + L62*musculo_2 + L72*cutane_2 + L82*vita_2
    Func_2  =~ L92*bas_act_2 + L102*comp_act_2 + L112*social_2 

    # Common Factor Variances (UVI Identification)
    Psych_2 ~~ 1*Psych_2
    Phys_2  ~~ 1*Phys_2
    Func_2  ~~ 1*Func_2

    # Residual Factor Variances
    psych_2    ~~ E12*psych_2
    social_2   ~~ E22*social_2
    vita_2     ~~ E32*vita_2
    somatic_2  ~~ E42*somatic_2
    gast_2     ~~ E52*gast_2
    musculo_2  ~~ E62*musculo_2
    cutane_2   ~~ E72*cutane_2
    bas_act_2  ~~ E82*bas_act_2
    comp_act_2 ~~ E92*comp_act_2

    # Common Factor Means
    Psych_2 ~ 0*1
    Phys_2  ~ 0*1
    Func_2  ~ 0*1

    # Intercepts
    psych_2    ~ INT12*1
    social_2   ~ INT22*1
    vita_2     ~ INT32*1
    somatic_2  ~ INT42*1
    gast_2     ~ INT52*1
    musculo_2  ~ INT62*1
    cutane_2   ~ INT72*1
    bas_act_2  ~ INT82*1
    comp_act_2 ~ INT92*1

    # Latent Correlations
    Psych_2 ~~ Phys_2
    Psych_2 ~~ Func_2
    Phys_2  ~~ Func_2
'    

# Run the model
T2_FIT <- cfa(T2_MODEL, 
              data = BDATA_numeric, 
              group = "condition",
              estimator = "MLR", 
              missing = "fiml",
              std.lv = TRUE)

# plot of the model we fitted
semPaths(T2_FIT)
semPaths(T2_FIT,whatLabels="par")

# View the full statistical output
summary(T2_FIT, fit.measures = TRUE, standardized = TRUE)


T5_MODEL = '
    # Pattern of Factor Loadings
    Psych_5 =~ L15*psych_5 + L25*social_5 + L35*vita_5
    Phys_5  =~ L45*somatic_5 + L55*gast_5 + L65*musculo_5 + L75*cutane_5 + L85*vita_5
    Func_5  =~ L95*bas_act_5 + L105*comp_act_5 + L115*social_5

    # Common Factor Variances (UVI Identification)
    Psych_5 ~~ 1*Psych_5
    Phys_5  ~~ 1*Phys_5
    Func_5  ~~ 1*Func_5

    # Residual Factor Variances
    psych_5    ~~ E15*psych_5
    social_5   ~~ E25*social_5
    vita_5     ~~ E35*vita_5
    somatic_5  ~~ E45*somatic_5
    gast_5     ~~ E55*gast_5
    musculo_5  ~~ E65*musculo_5
    cutane_5   ~~ E75*cutane_5
    bas_act_5  ~~ E85*bas_act_5
    comp_act_5 ~~ E95*comp_act_5

    # Common Factor Means
    Psych_5 ~ 0*1
    Phys_5  ~ 0*1
    Func_5  ~ 0*1

    # Intercepts
    psych_5    ~ INT15*1
    social_5   ~ INT25*1
    vita_5     ~ INT35*1
    somatic_5  ~ INT45*1
    gast_5     ~ INT55*1
    musculo_5  ~ INT65*1
    cutane_5   ~ INT75*1
    bas_act_5  ~ INT85*1
    comp_act_5 ~ INT95*1

    # Latent Correlations
    Psych_5 ~~ Phys_5
    Psych_5 ~~ Func_5
    Phys_5  ~~ Func_5
'    

# Run the T5 model
T5_FIT <- cfa(T5_MODEL, 
              data = BDATA_numeric, 
              group = "condition",
              estimator = "MLR", 
              missing = "fiml",
              std.lv = TRUE)

# plot of the model we fitted
semPaths(T5_FIT)
semPaths(T5_FIT,whatLabels="par")

# View the full statistical output
summary(T5_FIT, fit.measures = TRUE, standardized = TRUE)

# anova(T1_FIT, T2_FIT)
# anova(T1_FIT, T5_FIT)
# anova(T2_FIT, T5_FIT)

## Structural Equation Modeling ----

LONG_MODEL_STEP1 = '
    # FACTOR LOADINGS (Freely Estimated)
    Psych_1 =~ psych_1 + social_1 + vita_1
    Psych_2 =~ psych_2 + social_2 + vita_2
    Psych_5 =~ psych_5 + social_5 + vita_5

    Phys_1  =~ somatic_1 + gast_1 + musculo_1 + cutane_1 + vita_1
    Phys_2  =~ somatic_2 + gast_2 + musculo_2 + cutane_2 + vita_2
    Phys_5  =~ somatic_5 + gast_5 + musculo_5 + cutane_5 + vita_5

    Func_1  =~ bas_act_1 + comp_act_1 + social_1
    Func_2  =~ bas_act_2 + comp_act_2 + social_2
    Func_5  =~ bas_act_5 + comp_act_5 + social_5

    # INTERCEPTS (Freely Estimated)
    psych_1 ~ 1;    psych_2 ~ 1;    psych_5 ~ 1
    social_1 ~ 1;   social_2 ~ 1;   social_5 ~ 1
    vita_1 ~ 1;     vita_2 ~ 1;     vita_5 ~ 1
    somatic_1 ~ 1;  somatic_2 ~ 1;  somatic_5 ~ 1
    gast_1 ~ 1;     gast_2 ~ 1;     gast_5 ~ 1
    musculo_1 ~ 1;  musculo_2 ~ 1;  musculo_5 ~ 1
    cutane_1 ~ 1;   cutane_2 ~ 1;   cutane_5 ~ 1
    bas_act_1 ~ 1;  bas_act_2 ~ 1;  bas_act_5 ~ 1
    comp_act_1 ~ 1; comp_act_2 ~ 1; comp_act_5 ~ 1

    # RESIDUAL VARIANCES (Freely Estimated)
    psych_1    ~~ psych_1;    psych_2    ~~ psych_2;    psych_5    ~~ psych_5
    social_1   ~~ social_1;   social_2   ~~ social_2;   social_5   ~~ social_5
    vita_1     ~~ vita_1;     vita_2     ~~ vita_2;     vita_5     ~~ vita_5
    somatic_1  ~~ somatic_1;  somatic_2  ~~ somatic_2;  somatic_5  ~~ somatic_5
    gast_1     ~~ gast_1;     gast_2     ~~ gast_2;     gast_5     ~~ gast_5
    musculo_1  ~~ musculo_1;  musculo_2  ~~ musculo_2;  musculo_5  ~~ musculo_5
    cutane_1   ~~ cutane_1;   cutane_2   ~~ cutane_2;   cutane_5   ~~ cutane_5
    bas_act_1  ~~ bas_act_1;  bas_act_2  ~~ bas_act_2;  bas_act_5  ~~ bas_act_5
    comp_act_1 ~~ comp_act_1; comp_act_2 ~~ comp_act_2; comp_act_5 ~~ comp_act_5

    # RESIDUAL COVARIANCES
      # T1-T2
    psych_1 ~~ psych_2;       social_1 ~~ social_2;     vita_1 ~~ vita_2
    somatic_1 ~~ somatic_2;   gast_1 ~~ gast_2;         musculo_1 ~~ musculo_2
    cutane_1 ~~ cutane_2;     bas_act_1 ~~ bas_act_2;   comp_act_1 ~~ comp_act_2
      # T1-T5
    psych_1 ~~ psych_5;       social_1 ~~ social_5;     vita_1 ~~ vita_5
    somatic_1 ~~ somatic_5;   gast_1 ~~ gast_5;         musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_5;     bas_act_1 ~~ bas_act_5;   comp_act_1 ~~ comp_act_5
      # T2-T5
    psych_2 ~~ psych_5;       social_2 ~~ social_5;     vita_2 ~~ vita_5
    somatic_2 ~~ somatic_5;   gast_2 ~~ gast_5;         musculo_2 ~~ musculo_5
    cutane_2 ~~ cutane_5;     bas_act_2 ~~ bas_act_5;   comp_act_2 ~~ comp_act_5
    
    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5

    # LATENT FACTOR MEANS (Fixed to 0 for Step 1 identification)
    Psych_1 ~ 0*1; Psych_2 ~ 0*1; Psych_5 ~ 0*1
    Phys_1  ~ 0*1; Phys_2  ~ 0*1; Phys_5  ~ 0*1
    Func_1  ~ 0*1; Func_2  ~ 0*1; Func_5  ~ 0*1
'

FIT_STEP1 <- cfa(LONG_MODEL_STEP1, 
                 data = BDATA_numeric, 
                 group = "condition",
                 estimator = "MLR", 
                 missing = "fiml",
                 std.lv=TRUE) #UVI

summary(FIT_STEP1, fit.measures = TRUE, standardized = TRUE)

# modindic <- modificationIndices(FIT_STEP1)
# modindic[order(modindic$mi,decreasing=TRUE),][1:10,] 


LONG_MODEL_STEP2 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_24)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_24)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_8, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;    psych_5 ~ c(INT8_1, INT24_1)*1
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;   social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;     vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, INT24_4)*1;  somatic_5 ~ c(INT8_4, INT24_4)*1
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;     gast_5 ~ c(INT8_5, INT24_5)*1
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(INT8_6, INT24_6)*1;  musculo_5 ~ c(INT8_6, INT24_6)*1
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, INT24_7)*1;   cutane_5 ~ c(INT8_7, INT24_7)*1
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;  bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1; comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_5, E24_5)*gast_2;          gast_5 ~~ c(E8_5, E24_5)*gast_5
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_6, E24_6)*musculo_5
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2 <- cfa(LONG_MODEL_STEP2, 
                 data = BDATA_numeric, 
                 group = "condition",
                 estimator = "MLR", 
                 missing = "fiml",
                 std.lv=TRUE)

summary(FIT_STEP2, fit.measures = TRUE, standardized = TRUE)

# semPaths(FIT_STEP2,whatLabels="par") # unintelligible...

# test difference between models
lavTestLRT(FIT_STEP1, FIT_STEP2)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2)$uni[order(lavTestScore(FIT_STEP2)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2))

## RELAXING INT musculo_5 (MF) -> uniform recalibration----

LONG_MODEL_STEP2.1 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_24)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_24)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_8, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;    psych_5 ~ c(INT8_1, INT24_1)*1
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;   social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;     vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, INT24_4)*1;  somatic_5 ~ c(INT8_4, INT24_4)*1
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;     gast_5 ~ c(INT8_5, INT24_5)*1
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(INT8_6, INT24_6)*1;  musculo_5 ~ c(INT8_6, S24_6_T5)*1 #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, INT24_7)*1;   cutane_5 ~ c(INT8_7, INT24_7)*1
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;  bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1; comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_5, E24_5)*gast_2;          gast_5 ~~ c(E8_5, E24_5)*gast_5
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_6, E24_6)*musculo_5
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.1 <- cfa(LONG_MODEL_STEP2.1, 
                 data = BDATA_numeric, 
                 group = "condition",
                 estimator = "MLR", 
                 missing = "fiml",
                 std.lv=TRUE)

summary(FIT_STEP2.1, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2, FIT_STEP2.1)
lavTestLRT(FIT_STEP1, FIT_STEP2.1)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.1)$uni[order(lavTestScore(FIT_STEP2.1)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.1))

## RELAXING INT musculo_2 (MF) -> uniform recalibration----

LONG_MODEL_STEP2.2 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_24)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_24)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_8, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;    psych_5 ~ c(INT8_1, INT24_1)*1
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;   social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;     vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, INT24_4)*1;  somatic_5 ~ c(INT8_4, INT24_4)*1
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;     gast_5 ~ c(INT8_5, INT24_5)*1
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(INT8_6, S24_6_T2)*1;  musculo_5 ~ c(INT8_6, S24_6_T5)*1 #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, INT24_7)*1;   cutane_5 ~ c(INT8_7, INT24_7)*1
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;  bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1; comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_5, E24_5)*gast_2;          gast_5 ~~ c(E8_5, E24_5)*gast_5
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_6, E24_6)*musculo_5
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.2 <- cfa(LONG_MODEL_STEP2.2, 
                   data = BDATA_numeric, 
                   group = "condition",
                   estimator = "MLR", 
                   missing = "fiml",
                   std.lv=TRUE)

summary(FIT_STEP2.2, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.1, FIT_STEP2.2)
lavTestLRT(FIT_STEP1, FIT_STEP2.2)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.2)$uni[order(lavTestScore(FIT_STEP2.2)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.2))

## RELAXING INT musculo_5 (SF) -> uniform recalibration----

LONG_MODEL_STEP2.3 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_24)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_24)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_8, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;    psych_5 ~ c(INT8_1, INT24_1)*1
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;   social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;     vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, INT24_4)*1;  somatic_5 ~ c(INT8_4, INT24_4)*1
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;     gast_5 ~ c(INT8_5, INT24_5)*1
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(INT8_6, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1 #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, INT24_7)*1;   cutane_5 ~ c(INT8_7, INT24_7)*1
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;  bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1; comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_5, E24_5)*gast_2;          gast_5 ~~ c(E8_5, E24_5)*gast_5
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_6, E24_6)*musculo_5
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.3 <- cfa(LONG_MODEL_STEP2.3, 
                   data = BDATA_numeric, 
                   group = "condition",
                   estimator = "MLR", 
                   missing = "fiml",
                   std.lv=TRUE)

summary(FIT_STEP2.3, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.2, FIT_STEP2.3)
lavTestLRT(FIT_STEP1, FIT_STEP2.3)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.3)$uni[order(lavTestScore(FIT_STEP2.3)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.3))

## RELAXING INT musculo_2 (SF) -> uniform recalibration----

LONG_MODEL_STEP2.4 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_24)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_24)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_8, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;    psych_5 ~ c(INT8_1, INT24_1)*1
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;   social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;     vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, INT24_4)*1;  somatic_5 ~ c(INT8_4, INT24_4)*1
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;     gast_5 ~ c(INT8_5, INT24_5)*1
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1 #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, INT24_7)*1;   cutane_5 ~ c(INT8_7, INT24_7)*1
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;  bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1; comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_5, E24_5)*gast_2;          gast_5 ~~ c(E8_5, E24_5)*gast_5
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_6, E24_6)*musculo_5
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.4 <- cfa(LONG_MODEL_STEP2.4, 
                   data = BDATA_numeric, 
                   group = "condition",
                   estimator = "MLR", 
                   missing = "fiml",
                   std.lv=TRUE)

summary(FIT_STEP2.4, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.3, FIT_STEP2.4)
lavTestLRT(FIT_STEP1, FIT_STEP2.4)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.4)$uni[order(lavTestScore(FIT_STEP2.4)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.4))

## RELAXING RES VAR gast_2 (SF) -> non-uniform recalibration----

LONG_MODEL_STEP2.5 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_24)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_24)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_8, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;    psych_5 ~ c(INT8_1, INT24_1)*1
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;   social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;     vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, INT24_4)*1;  somatic_5 ~ c(INT8_4, INT24_4)*1
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;     gast_5 ~ c(INT8_5, INT24_5)*1
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1 #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, INT24_7)*1;   cutane_5 ~ c(INT8_7, INT24_7)*1
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;  bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1; comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5         #freed
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_6, E24_6)*musculo_5
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.5 <- cfa(LONG_MODEL_STEP2.5, 
                   data = BDATA_numeric, 
                   group = "condition",
                   estimator = "MLR", 
                   missing = "fiml",
                   std.lv=TRUE)

summary(FIT_STEP2.5, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.4, FIT_STEP2.5)
lavTestLRT(FIT_STEP1, FIT_STEP2.5)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.5)$uni[order(lavTestScore(FIT_STEP2.5)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.5))

## RELAXING INT psych_5 (MF) -> uniform recalibration----

LONG_MODEL_STEP2.6 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_24)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_24)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_8, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;      psych_5 ~ c(INT8_1, S24_1_T5)*1    #freed
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;     social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;       vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, INT24_4)*1;    somatic_5 ~ c(INT8_4, INT24_4)*1
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;       gast_5 ~ c(INT8_5, INT24_5)*1
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1 #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, INT24_7)*1;     cutane_5 ~ c(INT8_7, INT24_7)*1
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;    bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1;   comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5         #freed
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_6, E24_6)*musculo_5
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.6 <- cfa(LONG_MODEL_STEP2.6, 
                   data = BDATA_numeric, 
                   group = "condition",
                   estimator = "MLR", 
                   missing = "fiml",
                   std.lv=TRUE)

summary(FIT_STEP2.6, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.5, FIT_STEP2.6)
lavTestLRT(FIT_STEP1, FIT_STEP2.6)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.6)$uni[order(lavTestScore(FIT_STEP2.6)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.6))

## RELAXING INT gast_5 (MF) -> uniform recalibration----

LONG_MODEL_STEP2.7 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_24)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_24)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_8, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;      psych_5 ~ c(INT8_1, S24_1_T5)*1    #freed
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;     social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;       vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, INT24_4)*1;    somatic_5 ~ c(INT8_4, INT24_4)*1
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;       gast_5 ~ c(INT8_5, S24_5_T5)*1     #freed
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1 #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, INT24_7)*1;     cutane_5 ~ c(INT8_7, INT24_7)*1
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;    bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1;   comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5         #freed
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_6, E24_6)*musculo_5
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.7 <- cfa(LONG_MODEL_STEP2.7, 
                   data = BDATA_numeric, 
                   group = "condition",
                   estimator = "MLR", 
                   missing = "fiml",
                   std.lv=TRUE)

summary(FIT_STEP2.7, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.6, FIT_STEP2.7)
lavTestLRT(FIT_STEP1, FIT_STEP2.7)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.7)$uni[order(lavTestScore(FIT_STEP2.7)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.7))


## RELAXING INT somatic_2 (MF) -> uniform recalibration----

LONG_MODEL_STEP2.8 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_24)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_24)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_8, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;      psych_5 ~ c(INT8_1, S24_1_T5)*1     #freed
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;     social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;       vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, S24_4_T2)*1;    somatic_5 ~ c(INT8_4, INT24_4)*1   #freed
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;       gast_5 ~ c(INT8_5, S24_5_T5)*1      #freed
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1  #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, INT24_7)*1;     cutane_5 ~ c(INT8_7, INT24_7)*1
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;    bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1;   comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5         #freed
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_6, E24_6)*musculo_5
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.8 <- cfa(LONG_MODEL_STEP2.8, 
                   data = BDATA_numeric, 
                   group = "condition",
                   estimator = "MLR", 
                   missing = "fiml",
                   std.lv=TRUE)

summary(FIT_STEP2.8, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.7, FIT_STEP2.8)
lavTestLRT(FIT_STEP1, FIT_STEP2.8)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.8)$uni[order(lavTestScore(FIT_STEP2.8)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.8))

## RELAXING INT cutane_2 (MF) -> uniform recalibration----

LONG_MODEL_STEP2.9 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_24)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_24)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_8, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;      psych_5 ~ c(INT8_1, S24_1_T5)*1     #freed
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;     social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;       vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, S24_4_T2)*1;   somatic_5 ~ c(INT8_4, INT24_4)*1    #freed
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;       gast_5 ~ c(INT8_5, S24_5_T5)*1      #freed
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1  #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, S24_7_T2)*1;    cutane_5 ~ c(INT8_7, INT24_7)*1     #freed
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;    bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1;   comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5         #freed
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_6, E24_6)*musculo_5
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.9 <- cfa(LONG_MODEL_STEP2.9, 
                   data = BDATA_numeric, 
                   group = "condition",
                   estimator = "MLR", 
                   missing = "fiml",
                   std.lv=TRUE)

summary(FIT_STEP2.9, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.8, FIT_STEP2.9)
lavTestLRT(FIT_STEP1, FIT_STEP2.9)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.9)$uni[order(lavTestScore(FIT_STEP2.9)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.9))


## RELAXING FACT LOAD bas_act_5 (SF) -> reprioritization----

LONG_MODEL_STEP2.10 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_24)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_24)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_S_8, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5     #freed

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;      psych_5 ~ c(INT8_1, S24_1_T5)*1     #freed
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;     social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;       vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, S24_4_T2)*1;   somatic_5 ~ c(INT8_4, INT24_4)*1    #freed
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;       gast_5 ~ c(INT8_5, S24_5_T5)*1      #freed
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1  #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, S24_7_T2)*1;    cutane_5 ~ c(INT8_7, INT24_7)*1     #freed
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;    bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1;   comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5         #freed
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_6, E24_6)*musculo_5
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.10 <- cfa(LONG_MODEL_STEP2.10, 
                   data = BDATA_numeric, 
                   group = "condition",
                   estimator = "MLR", 
                   missing = "fiml",
                   std.lv=TRUE)

summary(FIT_STEP2.10, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.9, FIT_STEP2.10)
lavTestLRT(FIT_STEP1, FIT_STEP2.10)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.10)$uni[order(lavTestScore(FIT_STEP2.10)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.10))


## RELAXING RES VAR musculo_5 (SF) -> non-uniform recalibration----

LONG_MODEL_STEP2.11 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_24)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_24)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_S_8, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5     #freed

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;      psych_5 ~ c(INT8_1, S24_1_T5)*1     #freed
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;     social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;       vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, S24_4_T2)*1;   somatic_5 ~ c(INT8_4, INT24_4)*1    #freed
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;       gast_5 ~ c(INT8_5, S24_5_T5)*1      #freed
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1  #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, S24_7_T2)*1;    cutane_5 ~ c(INT8_7, INT24_7)*1     #freed
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;    bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1;   comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5         #freed
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_S_6_T5, E24_6)*musculo_5   #freed
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.11 <- cfa(LONG_MODEL_STEP2.11, 
                    data = BDATA_numeric, 
                    group = "condition",
                    estimator = "MLR", 
                    missing = "fiml",
                    std.lv=TRUE)

summary(FIT_STEP2.11, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.10, FIT_STEP2.11)
lavTestLRT(FIT_STEP1, FIT_STEP2.11)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.11)$uni[order(lavTestScore(FIT_STEP2.11)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.11))

## RELAXING FACT LOAD gast_5 (MF) -> reprioritization----

LONG_MODEL_STEP2.12 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_24)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_S_24)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5  #freed
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_S_8, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5     #freed

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;      psych_5 ~ c(INT8_1, S24_1_T5)*1     #freed
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;     social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;       vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, S24_4_T2)*1;   somatic_5 ~ c(INT8_4, INT24_4)*1    #freed
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;       gast_5 ~ c(INT8_5, S24_5_T5)*1      #freed
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1  #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, S24_7_T2)*1;    cutane_5 ~ c(INT8_7, INT24_7)*1     #freed
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;    bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1;   comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5         #freed
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_S_6_T5, E24_6)*musculo_5   #freed
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.12 <- cfa(LONG_MODEL_STEP2.12, 
                    data = BDATA_numeric, 
                    group = "condition",
                    estimator = "MLR", 
                    missing = "fiml",
                    std.lv=TRUE)

summary(FIT_STEP2.12, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.11, FIT_STEP2.12)
lavTestLRT(FIT_STEP1, FIT_STEP2.12)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.12)$uni[order(lavTestScore(FIT_STEP2.12)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.12))


## RELAXING FACT LOAD gast_2 (MF) -> reprioritization----

LONG_MODEL_STEP2.13 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_24)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_S_24_T2)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2       #freed
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_S_24_T5)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5  #freed
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_S_8_T5, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5       #freed

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;      psych_5 ~ c(INT8_1, S24_1_T5)*1     #freed
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;     social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;       vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, S24_4_T2)*1;   somatic_5 ~ c(INT8_4, INT24_4)*1    #freed
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;       gast_5 ~ c(INT8_5, S24_5_T5)*1      #freed
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1  #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, S24_7_T2)*1;    cutane_5 ~ c(INT8_7, INT24_7)*1     #freed
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;    bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1;   comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5             #freed
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_S_6_T5, E24_6)*musculo_5  #freed
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.13 <- cfa(LONG_MODEL_STEP2.13, 
                    data = BDATA_numeric, 
                    group = "condition",
                    estimator = "MLR", 
                    missing = "fiml",
                    std.lv=TRUE)

summary(FIT_STEP2.13, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.12, FIT_STEP2.13)
lavTestLRT(FIT_STEP1, FIT_STEP2.13)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.13)$uni[order(lavTestScore(FIT_STEP2.13)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.13))

## RELAXING FACT LOAD social_5 (MF) -> reprioritization----

LONG_MODEL_STEP2.14 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_S_24_T5)*social_5 + c(L3_8, L3_24)*vita_5  #freed

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_S_24_T2)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2       #freed
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_S_24_T5)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5  #freed
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_S_8_T5, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5       #freed

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;      psych_5 ~ c(INT8_1, S24_1_T5)*1     #freed
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;     social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;       vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(INT8_4, S24_4_T2)*1;   somatic_5 ~ c(INT8_4, INT24_4)*1    #freed
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;       gast_5 ~ c(INT8_5, S24_5_T5)*1      #freed
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1  #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, S24_7_T2)*1;    cutane_5 ~ c(INT8_7, INT24_7)*1     #freed
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;    bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1;   comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5             #freed
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_S_6_T5, E24_6)*musculo_5  #freed
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.14 <- cfa(LONG_MODEL_STEP2.14, 
                    data = BDATA_numeric, 
                    group = "condition",
                    estimator = "MLR", 
                    missing = "fiml",
                    std.lv=TRUE)

summary(FIT_STEP2.14, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.13, FIT_STEP2.14)
lavTestLRT(FIT_STEP1, FIT_STEP2.14)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.14)$uni[order(lavTestScore(FIT_STEP2.14)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.14))


## RELAXING INT somatic_2 (SF) -> uniform calibration----

LONG_MODEL_STEP2.15 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_S_24_T5)*social_5 + c(L3_8, L3_24)*vita_5  #freed

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_S_24_T2)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2       #freed
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_S_24_T5)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5  #freed
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_S_8_T5, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5       #freed

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;      psych_5 ~ c(INT8_1, S24_1_T5)*1     #freed
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;     social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;       vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(S8_4_T2, S24_4_T2)*1;   somatic_5 ~ c(INT8_4, INT24_4)*1   #freed #freed
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;       gast_5 ~ c(INT8_5, S24_5_T5)*1      #freed
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1  #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, S24_7_T2)*1;    cutane_5 ~ c(INT8_7, INT24_7)*1     #freed
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;    bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1;   comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_3)*vita_2;          vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5             #freed
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_S_6_T5, E24_6)*musculo_5  #freed
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.15 <- cfa(LONG_MODEL_STEP2.15, 
                    data = BDATA_numeric, 
                    group = "condition",
                    estimator = "MLR", 
                    missing = "fiml",
                    std.lv=TRUE)

summary(FIT_STEP2.15, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.14, FIT_STEP2.15)
lavTestLRT(FIT_STEP1, FIT_STEP2.15)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.15)$uni[order(lavTestScore(FIT_STEP2.15)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.15))

## RELAXING RES VAR vita_2 (MF) -> non-uniform calibration----

LONG_MODEL_STEP2.16 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_S_24_T5)*social_5 + c(L3_8, L3_24)*vita_5  #freed

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_S_24_T2)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2       #freed
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_S_24_T5)*gast_5 + c(L6_8, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5  #freed
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_S_8_T5, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5       #freed

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;      psych_5 ~ c(INT8_1, S24_1_T5)*1     #freed
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;     social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;       vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(S8_4_T2, S24_4_T2)*1;  somatic_5 ~ c(INT8_4, INT24_4)*1   #freed #freed
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;       gast_5 ~ c(INT8_5, S24_5_T5)*1      #freed
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1  #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, S24_7_T2)*1;    cutane_5 ~ c(INT8_7, INT24_7)*1     #freed
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;    bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1;   comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_S_3_T2)*vita_2;     vita_5 ~~ c(E8_3, E24_3)*vita_5        #freed
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5             #freed
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_S_6_T5, E24_6)*musculo_5  #freed
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.16 <- cfa(LONG_MODEL_STEP2.16, 
                    data = BDATA_numeric, 
                    group = "condition",
                    estimator = "MLR", 
                    missing = "fiml",
                    std.lv=TRUE)

summary(FIT_STEP2.16, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.15, FIT_STEP2.16)
lavTestLRT(FIT_STEP1, FIT_STEP2.16)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.16)$uni[order(lavTestScore(FIT_STEP2.16)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.16))

## RELAXING FACT LOAD musculo_5 (SF) -> reprioritization----

LONG_MODEL_STEP2.17 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_24)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_S_24_T5)*social_5 + c(L3_8, L3_24)*vita_5  #freed

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_S_24_T2)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2         #freed
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_S_24_T5)*gast_5 + c(L6_S_8_T5, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5    #freed #freed
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_S_8_T5, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5       #freed

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;      psych_5 ~ c(INT8_1, S24_1_T5)*1     #freed
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;     social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;       vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(S8_4_T2, S24_4_T2)*1;  somatic_5 ~ c(INT8_4, INT24_4)*1   #freed #freed
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;       gast_5 ~ c(INT8_5, S24_5_T5)*1      #freed
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1  #freed #freed #freed #freed
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, S24_7_T2)*1;    cutane_5 ~ c(INT8_7, INT24_7)*1     #freed
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;    bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1;   comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_S_3_T2)*vita_2;     vita_5 ~~ c(E8_3, E24_3)*vita_5        #freed
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5             #freed
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_S_6_T5, E24_6)*musculo_5  #freed
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
'

FIT_STEP2.17 <- cfa(LONG_MODEL_STEP2.17, 
                    data = BDATA_numeric, 
                    group = "condition",
                    estimator = "MLR", 
                    missing = "fiml",
                    std.lv=TRUE)

summary(FIT_STEP2.17, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.16, FIT_STEP2.17)
lavTestLRT(FIT_STEP1, FIT_STEP2.17)

# find the specific response shift mechanism to relax
lavTestScore(FIT_STEP2.17)$uni[order(lavTestScore(FIT_STEP2.17)$uni$X2, decreasing = TRUE),]
View(parTable(FIT_STEP2.17))


## RELAXING FACT LOAD musculo_5 (SF) -> reprioritization----

LONG_MODEL_STEP2.18 = '
    # FACTOR LOADINGS (Locked)
    # c(label_Condition8, label_Condition24)
    Psych_1 =~ c(L1_8, L1_24)*psych_1 + c(L2_8, L2_24)*social_1 + c(L3_8, L3_24)*vita_1
    Psych_2 =~ c(L1_8, L1_S_24_T2)*psych_2 + c(L2_8, L2_24)*social_2 + c(L3_8, L3_24)*vita_2
    Psych_5 =~ c(L1_8, L1_24)*psych_5 + c(L2_8, L2_S_24_T5)*social_5 + c(L3_8, L3_24)*vita_5

    Phys_1  =~ c(L4_8, L4_24)*somatic_1 + c(L5_8, L5_24)*gast_1 + c(L6_8, L6_24)*musculo_1 + c(L7_8, L7_24)*cutane_1 + c(L8_8, L8_24)*vita_1
    Phys_2  =~ c(L4_8, L4_24)*somatic_2 + c(L5_8, L5_S_24_T2)*gast_2 + c(L6_8, L6_24)*musculo_2 + c(L7_8, L7_24)*cutane_2 + c(L8_8, L8_24)*vita_2
    Phys_5  =~ c(L4_8, L4_24)*somatic_5 + c(L5_8, L5_S_24_T5)*gast_5 + c(L6_S_8_T5, L6_24)*musculo_5 + c(L7_8, L7_24)*cutane_5 + c(L8_8, L8_24)*vita_5
    
    Func_1  =~ c(L9_8, L9_24)*bas_act_1 + c(L10_8, L10_24)*comp_act_1 + c(L11_8, L11_24)*social_1
    Func_2  =~ c(L9_8, L9_24)*bas_act_2 + c(L10_8, L10_24)*comp_act_2 + c(L11_8, L11_24)*social_2
    Func_5  =~ c(L9_S_8_T5, L9_24)*bas_act_5 + c(L10_8, L10_24)*comp_act_5 + c(L11_8, L11_24)*social_5

    # INTERCEPTS (Locked)
    psych_1 ~ c(INT8_1, INT24_1)*1;    psych_2 ~ c(INT8_1, INT24_1)*1;      psych_5 ~ c(INT8_1, S24_1_T5)*1
    social_1 ~ c(INT8_2, INT24_2)*1;   social_2 ~ c(INT8_2, INT24_2)*1;     social_5 ~ c(INT8_2, INT24_2)*1
    vita_1 ~ c(INT8_3, INT24_3)*1;     vita_2 ~ c(INT8_3, INT24_3)*1;       vita_5 ~ c(INT8_3, INT24_3)*1
    somatic_1 ~ c(INT8_4, INT24_4)*1;  somatic_2 ~ c(S8_4_T2, S24_4_T2)*1;  somatic_5 ~ c(INT8_4, INT24_4)*1
    gast_1 ~ c(INT8_5, INT24_5)*1;     gast_2 ~ c(INT8_5, INT24_5)*1;       gast_5 ~ c(INT8_5, S24_5_T5)*1
    musculo_1 ~ c(INT8_6, INT24_6)*1;  musculo_2 ~ c(S8_6_T2, S24_6_T2)*1;  musculo_5 ~ c(S8_6_T5, S24_6_T5)*1
    cutane_1 ~ c(INT8_7, INT24_7)*1;   cutane_2 ~ c(INT8_7, S24_7_T2)*1;    cutane_5 ~ c(INT8_7, INT24_7)*1
    bas_act_1 ~ c(INT8_8, INT24_8)*1;  bas_act_2 ~ c(INT8_8, INT24_8)*1;    bas_act_5 ~ c(INT8_8, INT24_8)*1
    comp_act_1 ~ c(INT8_9, INT24_9)*1; comp_act_2 ~ c(INT8_9, INT24_9)*1;   comp_act_5 ~ c(INT8_9, INT24_9)*1

    # RESIDUAL VARIANCES (Locked)
    psych_1 ~~ c(E8_1, E24_1)*psych_1;        psych_2 ~~ c(E8_1, E24_1)*psych_2;        psych_5 ~~ c(E8_1, E24_1)*psych_5
    social_1 ~~ c(E8_2, E24_2)*social_1;      social_2 ~~ c(E8_2, E24_2)*social_2;      social_5 ~~ c(E8_2, E24_2)*social_5
    vita_1 ~~ c(E8_3, E24_3)*vita_1;          vita_2 ~~ c(E8_3, E24_S_3_T2)*vita_2;     vita_5 ~~ c(E8_3, E24_3)*vita_5
    somatic_1 ~~ c(E8_4, E24_4)*somatic_1;    somatic_2 ~~ c(E8_4, E24_4)*somatic_2;    somatic_5 ~~ c(E8_4, E24_4)*somatic_5
    gast_1 ~~ c(E8_5, E24_5)*gast_1;          gast_2 ~~ c(E8_S_5_T2, E24_5)*gast_2;     gast_5 ~~ c(E8_5, E24_5)*gast_5
    musculo_1 ~~ c(E8_6, E24_6)*musculo_1;    musculo_2 ~~ c(E8_6, E24_6)*musculo_2;    musculo_5 ~~ c(E8_S_6_T5, E24_6)*musculo_5
    cutane_1 ~~ c(E8_7, E24_7)*cutane_1;      cutane_2 ~~ c(E8_7, E24_7)*cutane_2;      cutane_5 ~~ c(E8_7, E24_7)*cutane_5
    bas_act_1 ~~ c(E8_8, E24_8)*bas_act_1;    bas_act_2 ~~ c(E8_8, E24_8)*bas_act_2;    bas_act_5 ~~ c(E8_8, E24_8)*bas_act_5
    comp_act_1 ~~ c(E8_9, E24_9)*comp_act_1;  comp_act_2 ~~ c(E8_9, E24_9)*comp_act_2;  comp_act_5 ~~ c(E8_9, E24_9)*comp_act_5

    # RESIDUAL COVARIANCES
    psych_1 ~~ psych_2;         psych_2 ~~ psych_5;         psych_1 ~~ psych_5
    social_1 ~~ social_2;       social_2 ~~ social_5;       social_1 ~~ social_5
    vita_1 ~~ vita_2;           vita_2 ~~ vita_5;           vita_1 ~~ vita_5
    somatic_1 ~~ somatic_2;     somatic_2 ~~ somatic_5;     somatic_1 ~~ somatic_5
    gast_1 ~~ gast_2;           gast_2 ~~ gast_5;           gast_1 ~~ gast_5
    musculo_1 ~~ musculo_2;     musculo_2 ~~ musculo_5;     musculo_1 ~~ musculo_5
    cutane_1 ~~ cutane_2;       cutane_2 ~~ cutane_5;       cutane_1 ~~ cutane_5
    bas_act_1 ~~ bas_act_2;     bas_act_2 ~~ bas_act_5;     bas_act_1 ~~ bas_act_5
    comp_act_1 ~~ comp_act_2;   comp_act_2 ~~ comp_act_5;   comp_act_1 ~~ comp_act_5

    # LATENT FACTOR COVARIANCES
      # SAME FACTOR OVER TIME
    Psych_1 ~~ Psych_2; Phys_1 ~~ Phys_2; Func_1 ~~ Func_2
    Psych_2 ~~ Psych_5; Phys_2 ~~ Phys_5; Func_2 ~~ Func_5
    Psych_1 ~~ Psych_5; Phys_1 ~~ Phys_5; Func_1 ~~ Func_5
      # INTER-FACTOR RELATIONSHIPS
    Psych_1 ~~ Phys_1; Psych_1 ~~ Func_1; Phys_1 ~~ Func_1
    Psych_2 ~~ Phys_2; Psych_2 ~~ Func_2; Phys_2 ~~ Func_2
    Psych_5 ~~ Phys_5; Psych_5 ~~ Func_5; Phys_5 ~~ Func_5
      # CROSS-LAGGED T1-T2
    Psych_1 ~~ Phys_2; Psych_1 ~~ Func_2
    Phys_1  ~~ Psych_2; Phys_1  ~~ Func_2
    Func_1  ~~ Psych_2; Func_1  ~~ Phys_2
      # CROSS-LAGGED T2-T5
    Psych_2 ~~ Phys_5; Psych_2 ~~ Func_5
    Phys_2  ~~ Psych_5; Phys_2  ~~ Func_5
    Func_2  ~~ Psych_5; Func_2  ~~ Phys_5
      # CROSS-LAGGED T1-T5
    Psych_1 ~~ Phys_5; Psych_1 ~~ Func_5
    Phys_1  ~~ Psych_5; Phys_1  ~~ Func_5
    Func_1  ~~ Psych_5; Func_1  ~~ Phys_5
    
    # LATENT FACTOR MEANS
    Psych_1 ~ c(0,0)*1; Psych_2 ~ c(K8_Ps2, K24_Ps2)*1; Psych_5 ~ c(K8_Ps5, K24_Ps5)*1
    Phys_1  ~ c(0,0)*1; Phys_2  ~ c(K8_Ph2, K24_Ph2)*1; Phys_5  ~ c(K8_Ph5, K24_Ph5)*1
    Func_1  ~ c(0,0)*1; Func_2  ~ c(K8_Fu2, K24_Fu2)*1; Func_5  ~ c(K8_Fu5, K24_Fu5)*1
    
    # MEAN CALCULATIONS (Condition 8)
    M_psy_T1_8   := INT8_1
    M_psy_T2_8   := INT8_1 + L1_8 * K8_Ps2
    M_psy_T5_8   := INT8_1 + L1_8 * K8_Ps5
    M_som_T1_8   := INT8_4
    M_som_T2_8   := S8_4_T2 + L4_8 * K8_Ph2
    M_som_T5_8   := INT8_4 + L4_8 * K8_Ph5
    M_gas_T1_8   := INT8_5
    M_gas_T2_8   := INT8_5 + L5_8 * K8_Ph2
    M_gas_T5_8   := INT8_5 + L5_8 * K8_Ph5
    M_mus_T1_8   := INT8_6
    M_mus_T2_8   := S8_6_T2 + L6_8 * K8_Ph2
    M_mus_T5_8   := S8_6_T5 + L6_S_8_T5 * K8_Ph5
    M_cut_T1_8   := INT8_7
    M_cut_T2_8   := INT8_7 + L7_8 * K8_Ph2
    M_cut_T5_8   := INT8_7 + L7_8 * K8_Ph5
    M_bas_T1_8   := INT8_8
    M_bas_T2_8   := INT8_8 + L9_8 * K8_Fu2
    M_bas_T5_8   := INT8_8 + L9_S_8_T5 * K8_Fu5
    M_com_T1_8   := INT8_9
    M_com_T2_8   := INT8_9 + L10_8 * K8_Fu2
    M_com_T5_8   := INT8_9 + L10_8 * K8_Fu5
    M_vita_T1_8  := INT8_3
    M_vita_T2_8  := INT8_3 + (L3_8 * K8_Ps2) + (L8_8 * K8_Ph2)
    M_vita_T5_8  := INT8_3 + (L3_8 * K8_Ps5) + (L8_8 * K8_Ph5)
    M_soc_T1_8   := INT8_2
    M_soc_T2_8   := INT8_2 + (L2_8 * K8_Ps2) + (L11_8 * K8_Fu2)
    M_soc_T5_8   := INT8_2 + (L2_8 * K8_Ps5) + (L11_8 * K8_Fu5)

    # MEAN CALCULATIONS (Condition 24)
    M_psy_T1_24  := INT24_1
    M_psy_T2_24  := INT24_1 + L1_S_24_T2 * K24_Ps2
    M_psy_T5_24  := S24_1_T5 + L1_24 * K24_Ps5
    M_som_T1_24  := INT24_4
    M_som_T2_24  := S24_4_T2 + L4_24 * K24_Ph2
    M_som_T5_24  := INT24_4 + L4_24 * K24_Ph5
    M_gas_T1_24  := INT24_5
    M_gas_T2_24  := INT24_5 + L5_S_24_T2 * K24_Ph2
    M_gas_T5_24  := S24_5_T5 + L5_S_24_T5 * K24_Ph5
    M_mus_T1_24  := INT24_6
    M_mus_T2_24  := S24_6_T2 + L6_24 * K24_Ph2
    M_mus_T5_24  := S24_6_T5 + L6_24 * K24_Ph5
    M_cut_T1_24  := INT24_7
    M_cut_T2_24  := S24_7_T2 + L7_24 * K24_Ph2
    M_cut_T5_24  := INT24_7 + L7_24 * K24_Ph5
    M_bas_T1_24  := INT24_8
    M_bas_T2_24  := INT24_8 + L9_24 * K24_Fu2
    M_bas_T5_24  := INT24_8 + L9_24 * K24_Fu5
    M_com_T1_24  := INT24_9
    M_com_T2_24  := INT24_9 + L10_24 * K24_Fu2
    M_com_T5_24  := INT24_9 + L10_24 * K24_Fu5
    M_vita_T1_24 := INT24_3
    M_vita_T2_24 := INT24_3 + (L3_24 * K24_Ps2) + (L8_24 * K24_Ph2)
    M_vita_T5_24 := INT24_3 + (L3_24 * K24_Ps5) + (L8_24 * K24_Ph5)
    M_soc_T1_24  := INT24_2
    M_soc_T2_24  := INT24_2 + (L2_24 * K24_Ps2) + (L11_24 * K24_Fu2)
    M_soc_T5_24  := INT24_2 + (L2_S_24_T5 * K24_Ps5) + (L11_24 * K24_Fu5)
    
'

FIT_STEP2.18 <- cfa(LONG_MODEL_STEP2.18, 
                    data = BDATA_numeric, 
                    group = "condition",
                    estimator = "MLR", 
                    missing = "fiml",
                    std.lv=TRUE)

summary(FIT_STEP2.18, fit.measures = TRUE, standardized = TRUE)

# test difference between models
lavTestLRT(FIT_STEP2.17, FIT_STEP2.18)
lavTestLRT(FIT_STEP1, FIT_STEP2.18)