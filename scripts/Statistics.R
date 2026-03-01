# ----- 1. Explanation of this script ------------------------------------------
# This script focuses on the statistical analysis of the manual feeding offers,
# associated compounds, and polyp activity under continuos exposure
# Analysis will be conducted using LMER and GLMER


# ----- 2. Load in needed packages ---------------------------------------------
# compute lmer and glmer
library(lme4)

# compute posthoc test for glmer - glht
library(multcomp)

# to easily read in all data files
library(readr)

# to easily clean data
library(tidyverse)
library(rstatix)

# work with times and duration
library(lubridate)


# ----- 3. Read in needed data files -------------------------------------------
# table with all manual feeding data
Feeding <- read_rds("in/Feeding_data.rds")

Activity <- read_csv2("in/Polyp_activity_data.csv")


# ----- 4. Prepare data for analysis -------------------------------------------
## ---- 4.1. Manual feeding offers ---------------------------------------------
# Convert "yes" and "no" to binary (1 and 0) for reaction, ingestion, and egestion
Fed_num <- Feeding %>%
  mutate(reaction = ifelse(reaction == "yes", 1, 0),
         ingestion = ifelse(ingestion == "yes", 1, 0),
         incorporation = ifelse(egestion == "yes", 0, 1),
         egestion = ifelse(egestion == "yes", 1, 0))

Overall <- Fed_num
Overall$particle <- factor(Overall$particle, 
                           levels = c("Copepods", "Artemia_cysts",
                                      "PE", "PVC", "PS",
                                      "PA", "PP", "PET",
                                      "Sand"))
levels(Overall$particle)



## ---- 4.2. Continuos exposure ------------------------------------------------
Polyps <- Activity %>%
  pivot_longer(cols=c(7:15), # use previous columns as new entries - categories
               names_to='tp',  # assign new names for the columns
               values_to='activity') %>%  
  # add foodnote column: wich feeding modus was used?
  separate(tp, sep="_", into=c("ran", "week", "foodnote")) %>% 
  # rename treatments
  mutate(Treatment = case_when(Treatment == "MPsallphases" ~ "MP",
                               Treatment == "Control" ~ "control",
                               Treatment == "MPs+Food" ~ "MP+food")) %>%
  rename(treat = "Treatment",
         col = "Colony",
         spec = "Species") %>%
  # remove unnecessary columns
   dplyr::select(-ran, -Tank)

# level activity, important for visualisation 
Polyps$activity <- factor(Polyps$activity, 
                          levels = c("active", "moderate", "inactive"))
# recheck
levels(Polyps$activity)

# level treatment, important for visualisation 
Polyps$treat <- factor(Polyps$treat, 
                       levels = c("control", "MP", "MP+food"))
# recheck
levels(Polyps$treat)

# Add numeric ranking
# for PA = active   --> 1
# for PA = moderate --> 0.5
# for PA = inactive --> 0
Polyps <- Polyps %>% 
  mutate(ranks = case_when(activity == "active"~ "1",
                           activity == "moderate"~ "0.5",
                           activity == "inactive"~ "0"),
         ranks = as.numeric(ranks),
         week = as.numeric(week))



# ----- 5. Statistics ----------------------------------------------------------
## ---- 5.1. Q1: Differences between Polymers ----------------------------------
# without biotic influences
# Compare all polymer types against each other
### --- 5.1.1. Pocillopora verrucosa -------------------------------------------

# Reaction
Q1_Pve <- Overall %>%
  filter(spec %in% c("Pve"),
         biofilm %in% c("blank", "steril"),
         trigger %in% c("NoTrigger"),
         student %in% c ("Celine", "Verena", "Moritz", "Tim"))

Q1_Pve_react <- glmer(reaction ~ particle + (1|col),
                family = binomial, data = Q1_Pve)

summary(glht(Q1_Pve_react, linfct = mcp(particle = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Fit: glmer(formula = reaction ~ particle + (1 | col), 
#            data = Q1_Pve, family = binomial)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   Artemia_cysts - Copepods == 0  0.27673    0.10902   2.538 0.089102 .  
# Artemia_cysts - Copepods == 0  0.16097    0.15478   1.040 1.000000    
# PE - Copepods == 0            -1.10428    0.16957  -6.512 1.85e-09 ***
# PVC - Copepods == 0           -1.35284    0.17657  -7.662 5.34e-13 ***
# PS - Copepods == 0            -1.32813    0.17165  -7.738 3.06e-13 ***
# PA - Copepods == 0            -0.94419    0.16589  -5.692 2.90e-07 ***
# PP - Copepods == 0            -1.05496    0.16837  -6.266 8.90e-09 ***
# PET - Copepods == 0           -2.12748    0.21099 -10.083  < 2e-16 ***
# Sand - Copepods == 0          -1.75603    0.18262  -9.616  < 2e-16 ***
# PE - Artemia_cysts == 0       -1.26525    0.16941  -7.469 2.27e-12 ***
# PVC - Artemia_cysts == 0      -1.51380    0.17643  -8.580  < 2e-16 ***
# PS - Artemia_cysts == 0       -1.48909    0.17150  -8.683  < 2e-16 ***
# PA - Artemia_cysts == 0       -1.10516    0.16573  -6.669 6.71e-10 ***
# PP - Artemia_cysts == 0       -1.21593    0.16821  -7.229 1.32e-11 ***
# PET - Artemia_cysts == 0      -2.28845    0.21089 -10.852  < 2e-16 ***
# Sand - Artemia_cysts == 0     -1.91699    0.18251 -10.504  < 2e-16 ***
# PVC - PE == 0                 -0.24856    0.18876  -1.317 1.000000    
# PS - PE == 0                  -0.22385    0.18417  -1.215 1.000000    
# PA - PE == 0                   0.16008    0.17897   0.894 1.000000    
# PP - PE == 0                   0.04932    0.18122   0.272 1.000000    
# PET - PE == 0                 -1.02321    0.22113  -4.627 7.41e-05 ***
# Sand - PE == 0                -0.65175    0.19437  -3.353 0.011986 *  
# PS - PVC == 0                  0.02471    0.19054   0.130 1.000000    
# PA - PVC == 0                  0.40864    0.18554   2.202 0.386886    
# PP - PVC == 0                  0.29787    0.18770   1.587 1.000000    
# PET - PVC == 0                -0.77465    0.22643  -3.421 0.009978 ** 
# Sand - PVC == 0               -0.40319    0.20040  -2.012 0.486505    
# PA - PS == 0                   0.38393    0.18087   2.123 0.405347    
# PP - PS == 0                   0.27317    0.18309   1.492 1.000000    
# PET - PS == 0                 -0.79936    0.22263  -3.591 0.005609 ** 
# Sand - PS == 0                -0.42790    0.19609  -2.182 0.386886    
# PP - PA == 0                  -0.11077    0.17786  -0.623 1.000000    
# PET - PA == 0                 -1.18329    0.21841  -5.418 1.33e-06 ***
# Sand - PA == 0                -0.81183    0.19126  -4.245 0.000416 ***
# PET - PP == 0                 -1.07252    0.22023  -4.870 2.34e-05 ***
# Sand - PP == 0                -0.70106    0.19335  -3.626 0.005182 ** 
# Sand - PET == 0                0.37146    0.23110   1.607 1.000000   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)

# Calculate Odds Ratio 
exp(try$test$coefficients)
# OUTPUT: MP - Food 
#  5.006343 
log(exp(EZ_1$test$coefficients))
# 1.610706 



# Ingestion
Q1_Pve_inge <- glmer(ingestion ~ particle + (1|col),
                  family = binomial, data = Q1_Pve)

summary(glht(Q1_Pve_inge, linfct = mcp(particle = "Tukey")), 
        test = adjusted("holm"))
# WARNING: degenerate Hessian with 1 negative eigenvalues
# OUTPUT:
# Fit: glmer(formula = Ingestion ~ particle + (1 | col), 
#            data = Q1_Pve, family = binomial)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# Artemia_cysts - Copepods == 0     0.0772     0.2297   0.336 1.000000    
# PE - Copepods == 0               -0.8643     0.3228  -2.677 0.180627    
# PVC - Copepods == 0              -3.8277     1.0337  -3.703 0.006610 ** 
# PS - Copepods == 0               -2.0468     0.4668  -4.385 0.000407 ***
# PA - Copepods == 0              -25.4379 40460.0287  -0.001 1.000000    
# PP - Copepods == 0               -4.0393     1.0306  -3.919 0.002931 ** 
# PET - Copepods == 0             -23.7165 28723.3290  -0.001 1.000000    
# Sand - Copepods == 0             -3.7073     1.0355  -3.580 0.009957 ** 
# PE - Artemia_cysts == 0          -0.9415     0.3185  -2.956 0.084219 .  
# PVC - Artemia_cysts == 0         -3.9049     1.0328  -3.781 0.005003 ** 
# PS - Artemia_cysts == 0          -2.1240     0.4642  -4.576 0.000171 ***
# PA - Artemia_cysts == 0         -25.5151 40460.0287  -0.001 1.000000    
# PP - Artemia_cysts == 0          -4.1165     1.0295  -3.999 0.002166 ** 
# PET - Artemia_cysts == 0        -23.7937 28723.3290  -0.001 1.000000    
# Sand - Artemia_cysts == 0        -3.7845     1.0349  -3.657 0.007657 ** 
# PVC - PE == 0                    -2.9634     1.0565  -2.805 0.130859    
# PS - PE == 0                     -1.1824     0.5142  -2.299 0.494115    
# PA - PE == 0                    -24.5735 40460.0287  -0.001 1.000000    
# PP - PE == 0                     -3.1749     1.0536  -3.014 0.072308 .  
# PET - PE == 0                   -22.8522 28723.3290  -0.001 1.000000    
# Sand - PE == 0                   -2.8430     1.0583  -2.686 0.180627    
# PS - PVC == 0                     1.7809     1.1086   1.606 1.000000    
# PA - PVC == 0                   -21.6102 40460.0287  -0.001 1.000000    
# PP - PVC == 0                    -0.2116     1.4400  -0.147 1.000000    
# PET - PVC == 0                  -19.8888 28723.3290  -0.001 1.000000    
# Sand - PVC == 0                   0.1204     1.4435   0.083 1.000000    
# PA - PS == 0                    -23.3911 40460.0287  -0.001 1.000000    
# PP - PS == 0                     -1.9925     1.1059  -1.802 1.000000    
# PET - PS == 0                   -21.6698 28723.3290  -0.001 1.000000    
# Sand - PS == 0                   -1.6606     1.1102  -1.496 1.000000    
# PP - PA == 0                     21.3986 40460.0287   0.001 1.000000    
# PET - PA == 0                     1.7213 49618.9837   0.000 1.000000    
# Sand - PA == 0                   21.7305 40460.0287   0.001 1.000000    
# PET - PP == 0                   -19.6772 28723.3290  -0.001 1.000000    
# Sand - PP == 0                    0.3319     1.4414   0.230 1.000000    
# Sand - PET == 0                  20.0092 28723.3290   0.001 1.000000 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# Incorporation
# Release
Q1_Pve_inco <- glmer(incorporation ~ particle + (1|col),
                 family = binomial, data = Q1_Pve)

summary(glht(Q1_Pve_inco, linfct = mcp(particle = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = incorporation ~ particle + (1 | col), 
#            data = Q1, family = binomial)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
# Artemia_cysts - Copepods == 0  1.072e-01  1.016e+00   0.106     1.00
# PE - Copepods == 0            -1.827e+00  9.581e-01  -1.907     1.00
# PVC - Copepods == 0            3.953e+01  6.711e+07   0.000     1.00
# PS - Copepods == 0             3.804e+01  2.740e+07   0.000     1.00
# PP - Copepods == 0            -8.315e+01  6.711e+07   0.000     1.00
# Sand - Copepods == 0           2.649e+01  3.148e+06   0.000     1.00
# PE - Artemia_cysts == 0       -1.934e+00  9.573e-01  -2.021     0.91
# PVC - Artemia_cysts == 0       3.942e+01  6.711e+07   0.000     1.00
# PS - Artemia_cysts == 0        3.793e+01  2.740e+07   0.000     1.00
# PP - Artemia_cysts == 0       -8.326e+01  6.711e+07   0.000     1.00
# Sand - Artemia_cysts == 0      2.639e+01  3.148e+06   0.000     1.00
# PVC - PE == 0                  4.135e+01  6.711e+07   0.000     1.00
# PS - PE == 0                   3.987e+01  2.740e+07   0.000     1.00
# PP - PE == 0                  -8.132e+01  6.711e+07   0.000     1.00
# Sand - PE == 0                 2.832e+01  3.148e+06   0.000     1.00
# PS - PVC == 0                 -1.486e+00  7.249e+07   0.000     1.00
# PP - PVC == 0                 -1.227e+02  9.491e+07   0.000     1.00
# Sand - PVC == 0               -1.303e+01  6.718e+07   0.000     1.00
# PP - PS == 0                  -1.212e+02  7.249e+07   0.000     1.00
# Sand - PS == 0                -1.155e+01  2.758e+07   0.000     1.00
# Sand - PP == 0                 1.096e+02  6.718e+07   0.000     1.00
# (Adjusted p values reported -- holm method)




#### - S5.1.1. Time calculations -----------------------------------------------
Time_Q1_Pve <- Q1_Pve %>%
  # bring duration columns in time format hh:mm:ss
  mutate(reaction_dur = lubridate::hms(reaction_dur),
         ingestion_dur = lubridate::hms(ingestion_dur)) %>% 
  # bring durations in numeric to be able to calculate with
  mutate(reaction_dur_num = as.numeric(reaction_dur),
         ingestion_dur_num = as.numeric(ingestion_dur))

# calculate mean reaction time
mean_react_time_Q1_Pve <- Time_Q1_Pve %>%
  group_by(particle) %>%
  get_summary_stats(reaction_dur_num, type = "mean_sd")
# format mean and sd back into seconds
mean_react_time_Q1_Pve <- mean_react_time_Q1_Pve %>% 
  mutate(mean_react_dur = seconds_to_period(mean),
         sd_react_dur = seconds_to_period(sd))

# calculate mean ingestion time
mean_ing_time_Q1_Pve <- Time_Q1_Pve %>%
  group_by(particle) %>%
  get_summary_stats(ingestion_dur_num, type = "mean_sd")
# format mean and sd back into seconds
mean_ing_time_Q1_Pve <- mean_ing_time_Q1_Pve %>% 
  mutate(mean_ing_dur = seconds_to_period(mean),
         sd_ing_dur = seconds_to_period(sd))

# testing whether reaction time differs
Wilcox_Tab1 <- Time_Q1_Pve %>%
  wilcox_test(reaction_dur_num ~ particle, data= .)

# testing whether ingestion time differs
Time_Q1_Pve %>%
  wilcox_test(ingestion_dur_num ~ particle, data= .)
# not enough observations


### --- 5.1.2. Stylophora pistillata -------------------------------------------

# Reaction
Q1_Spi <- Overall %>%
  filter(spec %in% c("Spi"),
         biofilm %in% c("blank", "steril"),
         trigger %in% c("NoTrigger"),
         student %in% c ("Celine", "Verena", "Moritz", "Tim"))

Q1_Spi_react <- glmer(reaction ~ particle + (1|col),
                      family = binomial, data = Q1_Spi)

summary(glht(Q1_Spi_react, linfct = mcp(particle = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Fit: glmer(formula = reaction ~ particle + (1 | col), 
#            data = Q1_Spi, family = binomial)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# Artemia_cysts - Copepods == 0  0.431605   0.159702   2.703 0.082569 .  
# PE - Copepods == 0            -0.840669   0.176047  -4.775 3.41e-05 ***
# PVC - Copepods == 0           -2.457556   0.259901  -9.456  < 2e-16 ***
# PS - Copepods == 0            -1.965063   0.218331  -9.000  < 2e-16 ***
# PA - Copepods == 0            -1.488649   0.197509  -7.537 1.25e-12 ***
# PP - Copepods == 0            -1.063800   0.181706  -5.855 1.10e-07 ***
# PET - Copepods == 0           -2.072836   0.230994  -8.974  < 2e-16 ***
# Sand - Copepods == 0          -1.484988   0.187876  -7.904 7.19e-14 ***
# PE - Artemia_cysts == 0       -1.272555   0.174224  -7.304 6.98e-12 ***
# PVC - Artemia_cysts == 0      -2.889443   0.258950 -11.158  < 2e-16 ***
# PS - Artemia_cysts == 0       -2.396949   0.217110 -11.040  < 2e-16 ***
# PA - Artemia_cysts == 0       -1.920535   0.196091  -9.794  < 2e-16 ***
# PP - Artemia_cysts == 0       -1.495686   0.180045  -8.307  < 2e-16 ***
# PET - Artemia_cysts == 0      -2.504722   0.229865 -10.897  < 2e-16 ***
# Sand - Artemia_cysts == 0     -1.916874   0.186387 -10.284  < 2e-16 ***
# PVC - PE == 0                 -1.616887   0.268015  -6.033 3.87e-08 ***
# PS - PE == 0                  -1.124394   0.228036  -4.931 1.64e-05 ***
# PA - PE == 0                  -0.647980   0.208294  -3.111 0.024248 *  
# PP - PE == 0                  -0.223131   0.193561  -1.153 0.836369    
# PET - PE == 0                 -1.232167   0.240158  -5.131 6.29e-06 ***
# Sand - PE == 0                -0.644319   0.199186  -3.235 0.017044 *  
# PS - PVC == 0                  0.492493   0.296919   1.659 0.485903    
# PA - PVC == 0                  0.968907   0.282115   3.434 0.008907 ** 
# PP - PVC == 0                  1.393757   0.271561   5.132 6.29e-06 ***
# PET - PVC == 0                 0.384720   0.306290   1.256 0.836369    
# Sand - PVC == 0                0.972568   0.275462   3.531 0.006632 ** 
# PA - PS == 0                   0.476414   0.244527   1.948 0.350576    
# PP - PS == 0                   0.901263   0.232235   3.881 0.001770 ** 
# PET - PS == 0                 -0.107773   0.272082  -0.396 1.000000    
# Sand - PS == 0                 0.480075   0.236819   2.027 0.350576    
# PP - PA == 0                   0.424849   0.212897   1.996 0.350576    
# PET - PA == 0                 -0.584187   0.255848  -2.283 0.224108    
# Sand - PA == 0                 0.003661   0.217925   0.017 1.000000    
# PET - PP == 0                 -1.009036   0.244135  -4.133 0.000644 ***
# Sand - PP == 0                -0.421188   0.203996  -2.065 0.350576    
# Sand - PET == 0                0.587848   0.248489   2.366 0.197963 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# Ingestion
Q1_Spi_inge <- glmer(ingestion ~ particle + (1|col),
                     family = binomial, data = Q1_Spi)

summary(glht(Q1_Spi_inge, linfct = mcp(particle = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Fit: glmer(formula = Ingestion ~ particle + (1 | col), 
#            data = Q1_Spi, family = binomial)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# Artemia_cysts - Copepods == 0  4.857e-01  2.400e-01   2.023 0.860772    
# PE - Copepods == 0            -4.526e-01  2.877e-01  -1.573 1.000000    
# PVC - Copepods == 0           -2.468e+00  7.655e-01  -3.224 0.027860 *  
# PS - Copepods == 0            -2.573e+00  6.301e-01  -4.083 0.001244 ** 
# PA - Copepods == 0            -2.645e+00  5.510e-01  -4.800 4.75e-05 ***
# PP - Copepods == 0            -3.315e+00  6.161e-01  -5.380 2.46e-06 ***
# PET - Copepods == 0           -2.864e+01  2.728e+05   0.000 1.000000    
# Sand - Copepods == 0          -2.573e+00  5.004e-01  -5.142 8.71e-06 ***
# PE - Artemia_cysts == 0       -9.383e-01  2.821e-01  -3.326 0.021167 *  
# PVC - Artemia_cysts == 0      -2.953e+00  7.634e-01  -3.868 0.002847 ** 
# PS - Artemia_cysts == 0       -3.059e+00  6.276e-01  -4.873 3.40e-05 ***
# PA - Artemia_cysts == 0       -3.131e+00  5.482e-01  -5.712 3.80e-07 ***
# PP - Artemia_cysts == 0       -3.800e+00  6.135e-01  -6.194 2.11e-08 ***
# PET - Artemia_cysts == 0      -2.912e+01  2.728e+05   0.000 1.000000    
# Sand - Artemia_cysts == 0     -3.059e+00  4.972e-01  -6.151 2.69e-08 ***
# PVC - PE == 0                 -2.015e+00  7.797e-01  -2.584 0.204975    
# PS - PE == 0                  -2.120e+00  6.473e-01  -3.275 0.024272 *  
# PA - PE == 0                  -2.193e+00  5.706e-01  -3.842 0.003047 ** 
# PP - PE == 0                  -2.862e+00  6.337e-01  -4.517 0.000182 ***
# PET - PE == 0                 -2.818e+01  2.728e+05   0.000 1.000000    
# Sand - PE == 0                -2.120e+00  5.219e-01  -4.063 0.001310 ** 
# PS - PVC == 0                 -1.054e-01  9.603e-01  -0.110 1.000000    
# PA - PVC == 0                 -1.777e-01  9.104e-01  -0.195 1.000000    
# PP - PVC == 0                 -8.473e-01  9.512e-01  -0.891 1.000000    
# PET - PVC == 0                -2.617e+01  2.728e+05   0.000 1.000000    
# Sand - PVC == 0               -1.054e-01  8.807e-01  -0.120 1.000000    
# PA - PS == 0                  -7.232e-02  8.000e-01  -0.090 1.000000    
# PP - PS == 0                  -7.419e-01  8.461e-01  -0.877 1.000000    
# PET - PS == 0                 -2.606e+01  2.728e+05   0.000 1.000000    
# Sand - PS == 0                 4.788e-09  7.659e-01   0.000 1.000000    
# PP - PA == 0                  -6.696e-01  7.890e-01  -0.849 1.000000    
# PET - PA == 0                 -2.599e+01  2.728e+05   0.000 1.000000    
# Sand - PA == 0                 7.232e-02  7.023e-01   0.103 1.000000    
# PET - PP == 0                 -2.532e+01  2.728e+05   0.000 1.000000    
# Sand - PP == 0                 7.419e-01  7.545e-01   0.983 1.000000    
# Sand - PET == 0                2.606e+01  2.728e+05   0.000 1.000000   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# Incorporation
# Release
Q1_Spi_inco <- glmer(incorporation ~ particle + (1|col),
                     family = binomial, data = Q1_Spi)

summary(glht(Q1_Spi_inco, linfct = mcp(particle = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = incorporation ~ particle + (1 | col), 
#            data = Q1_Spi, family = binomial)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
# Artemia_cysts - Copepods == 0  4.199e+01  6.258e+06   0.000   1.0000  
# PE - Copepods == 0            -1.514e+00  1.242e+00  -1.219   1.0000  
# PVC - Copepods == 0            2.962e+01  1.632e+07   0.000   1.0000  
# PS - Copepods == 0             4.558e+01  3.875e+07   0.000   1.0000  
# PA - Copepods == 0            -4.317e+00  1.419e+00  -3.043   0.0656 .
# PP - Copepods == 0            -4.770e+01  3.875e+07   0.000   1.0000  
# Sand - Copepods == 0          -2.931e+00  1.504e+00  -1.948   1.0000  
# PE - Artemia_cysts == 0       -4.350e+01  6.258e+06   0.000   1.0000  
# PVC - Artemia_cysts == 0      -1.237e+01  1.748e+07   0.000   1.0000  
# PS - Artemia_cysts == 0        3.596e+00  3.925e+07   0.000   1.0000  
# PA - Artemia_cysts == 0       -4.630e+01  6.258e+06   0.000   1.0000  
# PP - Artemia_cysts == 0       -8.969e+01  3.925e+07   0.000   1.0000  
# Sand - Artemia_cysts == 0     -4.492e+01  6.258e+06   0.000   1.0000  
# PVC - PE == 0                  3.113e+01  1.632e+07   0.000   1.0000  
# PS - PE == 0                   4.710e+01  3.875e+07   0.000   1.0000  
# PA - PE == 0                  -2.803e+00  1.237e+00  -2.266   0.6329  
# PP - PE == 0                  -4.618e+01  3.875e+07   0.000   1.0000  
# Sand - PE == 0                -1.417e+00  1.334e+00  -1.062   1.0000  
# PS - PVC == 0                  1.596e+01  4.204e+07   0.000   1.0000  
# PA - PVC == 0                 -3.394e+01  1.632e+07   0.000   1.0000  
# PP - PVC == 0                 -7.732e+01  4.204e+07   0.000   1.0000  
# Sand - PVC == 0               -3.255e+01  1.632e+07   0.000   1.0000  
# PA - PS == 0                  -4.990e+01  3.875e+07   0.000   1.0000  
# PP - PS == 0                  -9.328e+01  5.479e+07   0.000   1.0000  
# Sand - PS == 0                -4.852e+01  3.875e+07   0.000   1.0000  
# PP - PA == 0                  -4.338e+01  3.875e+07   0.000   1.0000  
# Sand - PA == 0                 1.386e+00  1.500e+00   0.924   1.0000  
# Sand - PP == 0                 4.477e+01  3.875e+07   0.000   1.0000  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)



#### - S5.1.2. Time calculations -----------------------------------------------
Time_Q1_Spi <- Q1_Spi %>%
  # bring duration columns in time format hh:mm:ss
  mutate(reaction_dur = lubridate::hms(reaction_dur),
         ingestion_dur = lubridate::hms(ingestion_dur)) %>% 
  # bring durations in numeric to be able to calculate with
  mutate(reaction_dur_num = as.numeric(reaction_dur),
         ingestion_dur_num = as.numeric(ingestion_dur))

# calculate mean reaction time
mean_react_time_Q1_Spi <- Time_Q1_Spi %>%
  group_by(particle) %>%
  get_summary_stats(reaction_dur_num, type = "mean_sd")
# format mean and sd back into seconds
mean_react_time_Q1_Spi <- mean_react_time_Q1_Spi %>% 
  mutate(mean_react_dur = seconds_to_period(mean),
         sd_react_dur = seconds_to_period(sd))

# calculate mean ingestion time
mean_ing_time_Q1_Spi <- Time_Q1_Spi %>%
  group_by(particle) %>%
  get_summary_stats(ingestion_dur_num, type = "mean_sd")
# format mean and sd back into seconds
mean_ing_time_Q1_Spi <- mean_ing_time_Q1_Spi %>% 
  mutate(mean_ing_dur = seconds_to_period(mean),
         sd_ing_dur = seconds_to_period(sd))

# testing whether reaction time differs
Wilcox_Tab2 <- Time_Q1_Spi %>%
  wilcox_test(reaction_dur_num ~ particle, data= .)

# testing whether ingestion time differs
Time_Q1_Spi %>%
  wilcox_test(ingestion_dur_num ~ particle, data= .)
# not enough observations



## ---- 5.2. Q2: Differences between Shapes ------------------------------------
# without biotic influences
# Compare all PP and PA as fragments against as fibers
### --- 5.2.1. Pocillopora verrucosa -------------------------------------------
# PA
Q2_Pve_PA <- Overall %>%
  filter(spec %in% c("Pve"),
         particle %in% c("PA"),
         biofilm %in% c("blank", "steril"),
         trigger %in% c("NoTrigger"))

# PP
Q2_Pve_PP <- Overall %>%
  filter(spec %in% c("Pve"),
         particle %in% c("PP"),
         biofilm %in% c("blank", "steril"),
         trigger %in% c("NoTrigger"))

# Reaction to PA
Q2_Pve_PA_react <- glmer(reaction ~ shape + (1|col),
                      family = binomial, data = Q2_Pve_PA)

summary(glht(Q2_Pve_PA_react, linfct = mcp(shape = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = reaction ~ shape + (1 | col), data = Q2_Pve_PA, 
#            family = binomial)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# fragment - fiber == 0  -1.3817     0.6278  -2.201   0.0277 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# Reaction to PP
Q2_Pve_PP_react <- glmer(reaction ~ shape + (1|col),
                         family = binomial, data = Q2_Pve_PP)

summary(glht(Q2_Pve_PP_react, linfct = mcp(shape = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = reaction ~ shape + (1 | col), data = Q2_Pve_PP, 
# family = binomial)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# fragment - fiber == 0  -1.2297     0.5522  -2.227   0.0259 *
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# Ingestion of PP
Q2_Pve_PP_inge <- glmer(ingestion ~ shape + (1|col),
                     family = binomial, data = Q2_Pve_PP)

summary(glht(Q2_Pve_PP_inge, linfct = mcp(shape = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
#  Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = ingestion ~ shape + (1 | col), data = Q2_Pve_PP, 
#            family = binomial)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# fragment - fiber == 0     -30.28 9899617.58       0        1
# (Adjusted p values reported -- holm method)

# Ingestion of PA
Q2_Pve_PA_inge <- glmer(ingestion ~ shape + (1|col),
                        family = binomial, data = Q2_Pve_PA)

# Incorporation of PP
Q2_Pve_PP_incor <- glmer(incorporation ~ shape + (1|col),
                        family = binomial, data = Q2_Pve_PP)

# Incorporation of PA
# No data


#### - S5.2.1. Time calculations -----------------------------------------------
Time_Q2_Pve <- Overall %>%
  filter(spec %in% c("Pve"),
         particle %in% c("PP", "PA"),
         biofilm %in% c("blank", "steril"),
         trigger %in% c("NoTrigger")) %>%
  # bring duration columns in time format hh:mm:ss
  mutate(reaction_dur = lubridate::hms(reaction_dur),
         ingestion_dur = lubridate::hms(ingestion_dur)) %>% 
  # bring durations in numeric to be able to calculate with
  mutate(reaction_dur_num = as.numeric(reaction_dur),
         ingestion_dur_num = as.numeric(ingestion_dur))

# calculate mean reaction time
mean_react_time_Q2_Pve <- Time_Q2_Pve %>%
  group_by(particle, shape) %>%
  get_summary_stats(reaction_dur_num, type = "mean_sd")

# format mean and sd back into seconds
mean_react_time_Q2_Pve <- mean_react_time_Q2_Pve %>% 
  mutate(mean_react_dur = seconds_to_period(mean),
         sd_react_dur = seconds_to_period(sd))
mean_react_time_Q2_Pve 

# calculate mean ingestion time
mean_ing_time_Q2_Pve <- Time_Q2_Pve %>%
  group_by(particle, shape) %>%
  get_summary_stats(ingestion_dur_num, type = "mean_sd")
# format mean and sd back into seconds
mean_ing_time_Q2_Pve <- mean_ing_time_Q2_Pve %>% 
  mutate(mean_ing_dur = seconds_to_period(mean),
         sd_ing_dur = seconds_to_period(sd))

# testing whether reaction time differs
# PP
Time_Q2_Pve %>%
  filter(particle=="PP") %>% 
  wilcox_test(reaction_dur_num ~ shape, data= .)
# OUTPUT: # A tibble: 1 × 7
# .y.              group1 group2      n1    n2 statistic             p
# * <chr>            <chr>  <chr>    <int> <int>     <dbl>         <dbl>
#  1 reaction_dur_num fiber  fragment    84    13      21.5 0.00000000483

# testing whether ingestion time differs
Time_Q2_Pve %>%
  filter(particle=="PP") %>% 
  wilcox_test(ingestion_dur_num ~ shape, data= .)
# not enough observations

# PA
Time_Q2_Pve %>%
  filter(particle=="PA") %>% 
  wilcox_test(reaction_dur_num ~ shape, data= .)
# OUTPUT: # A tibble: 1 × 7
# .y.              group1 group2      n1    n2 statistic         p
# * <chr>            <chr>  <chr>    <int> <int>     <dbl>     <dbl>
#  1 reaction_dur_num fiber  fragment    88    11      104.5 0.0000182

# testing whether ingestion time differs
Time_Q2_Pve %>%
  filter(particle=="PA") %>% 
  wilcox_test(ingestion_dur_num ~ shape, data= .)
# not enough observations




### --- 5.2.2. Stylophora pistillata -------------------------------------------
# PA
Q2_Spi_PA <- Overall %>%
  filter(spec %in% c("Spi"),
         particle %in% c("PA"),
         biofilm %in% c("blank", "steril"),
         trigger %in% c("NoTrigger"))

# PP
Q2_Spi_PP <- Overall %>%
  filter(spec %in% c("Spi"),
         particle %in% c("PP"),
         biofilm %in% c("blank", "steril"),
         trigger %in% c("NoTrigger"))

# Reaction to PA
Q2_Spi_PA_react <- glmer(reaction ~ shape + (1|col),
                         family = binomial, data = Q2_Spi_PA)

summary(glht(Q2_Spi_PA_react, linfct = mcp(shape = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = reaction ~ shape + (1 | col), data = Q2_Spi_PA, 
#            family = binomial)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# fragment - fiber == 0  0.08717    0.65208   0.134    0.894
# (Adjusted p values reported -- holm method)


# Reaction to PP
Q2_Spi_PP_react <- glmer(reaction ~ shape + (1|col),
                         family = binomial, data = Q2_Spi_PP)

summary(glht(Q2_Spi_PP_react, linfct = mcp(shape = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = reaction ~ shape + (1 | col), data = Q2_Spi_PP, 
# family = binomial)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# fragment - fiber == 0  -0.7033     0.3043  -2.311   0.0208 *
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# Ingestion of PP
Q2_Spi_PP_inge <- glmer(ingestion ~ shape + (1|col),
                        family = binomial, data = Q2_Spi_PP)

summary(glht(Q2_Spi_PP_inge, linfct = mcp(shape = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
#  Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = ingestion ~ shape + (1 | col), data = Q2_Spi_PP, 
#            family = binomial)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# fragment - fiber == 0    -2.842e+01  1.733e+07       0        1
# (Adjusted p values reported -- holm method)

# Ingestion of PA
Q2_Spi_PA_inge <- glmer(ingestion ~ shape + (1|col),
                        family = binomial, data = Q2_Spi_PA)
summary(glht(Q2_Spi_PA_inge, linfct = mcp(shape = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = ingestion ~ shape + (1 | col), data = Q2_Spi_PA, 
#            family = binomial)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# fragment - fiber == 0     -29.86 5457439.89       0        1
# (Adjusted p values reported -- holm method)

# Incorporation of PP
Q2_Spi_PP_incor <- glmer(incorporation ~ shape + (1|col),
                         family = binomial, data = Q2_Spi_PP)

# Incorporation of PA
Q2_Spi_PA_incor <- glmer(incorporation ~ shape + (1|col),
                         family = binomial, data = Q2_Spi_PA)



#### - S5.2.2. Time calculations -----------------------------------------------
Time_Q2_Spi <- Overall %>%
  filter(spec %in% c("Spi"),
         particle %in% c("PP", "PA"),
         biofilm %in% c("blank", "steril"),
         trigger %in% c("NoTrigger")) %>%
  # bring duration columns in time format hh:mm:ss
  mutate(reaction_dur = lubridate::hms(reaction_dur),
         ingestion_dur = lubridate::hms(ingestion_dur)) %>% 
  # bring durations in numeric to be able to calculate with
  mutate(reaction_dur_num = as.numeric(reaction_dur),
         ingestion_dur_num = as.numeric(ingestion_dur))

# calculate mean reaction time
mean_react_time_Q2_Spi <- Time_Q2_Spi %>%
  group_by(particle, shape) %>%
  get_summary_stats(reaction_dur_num, type = "mean_sd")

# format mean and sd back into seconds
mean_react_time_Q2_Spi <- mean_react_time_Q2_Spi %>% 
  mutate(mean_react_dur = seconds_to_period(mean),
         sd_react_dur = seconds_to_period(sd))
mean_react_time_Q2_Spi

# calculate mean ingestion time
mean_ing_time_Q2_Spi <- Time_Q2_Spi %>%
  group_by(particle, shape) %>%
  get_summary_stats(ingestion_dur_num, type = "mean_sd")
# format mean and sd back into seconds
mean_ing_time_Q2_Spi <- mean_ing_time_Q2_Spi %>% 
  mutate(mean_ing_dur = seconds_to_period(mean),
         sd_ing_dur = seconds_to_period(sd))

# testing whether reaction time differs
# PP
Time_Q2_Spi %>%
  filter(particle=="PP") %>% 
  wilcox_test(reaction_dur_num ~ shape, data= .)
# OUTPUT: # A tibble: 1 × 7
# .y.              group1 group2      n1    n2 statistic             p
# * <chr>            <chr>  <chr>    <int> <int>     <dbl>         <dbl>
#  1 reaction_dur_num fiber  fragment     64    16      75.5 0.0000000773

# testing whether ingestion time differs
Time_Q2_Spi %>%
  filter(particle=="PP") %>% 
  wilcox_test(ingestion_dur_num ~ shape, data= .)
# not enough observations

# PA
Time_Q2_Spi %>%
  filter(particle=="PA") %>% 
  wilcox_test(reaction_dur_num ~ shape, data= .)
# OUTPUT: # A tibble: 1 × 7
# .y.              group1 group2      n1    n2 statistic         p
# * <chr>            <chr>  <chr>    <int> <int>     <dbl>     <dbl>
#  1 reaction_dur_num fiber  fragment    45    16        76 0.00000279

# testing whether ingestion time differs
Time_Q2_Spi %>%
  filter(particle=="PA") %>% 
  wilcox_test(ingestion_dur_num ~ shape, data= .)
# not enough observations



## ---- 5.3. Q3: Differences between biotic factors ----------------------------
# Compare all biotic factors against each other
# separate for fragments and fibers
### --- 5.3.1. Pocillopora verrucosa -------------------------------------------
# Fragment
Q3_Pve_Frag <- Overall %>%
  filter(spec %in% c("Pve"),
         shape %in% c("fragment"),
         student %in% c ("Celine", "Verena", "Moritz", "Tim")) %>%
  unite(biotics, c(biofilm, trigger), sep = "_", remove = FALSE)

# Fiber
Q3_Pve_Fib <- Overall %>%
  filter(spec %in% c("Pve"),
         shape %in% c("fiber"),
         student %in% c ("Celine", "Verena", "Moritz", "Tim")) %>%
  unite(biotics, c(biofilm, trigger), sep = "_", remove = FALSE)


# Reaction
# Fragment
Q3_Pve_Frag_react <- glmer(reaction ~ biotics + (1|col),
                      family = binomial, data = Q3_Pve_Frag)

summary(glht(Q3_Pve_Frag_react, linfct = mcp(biotics = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = reaction ~ biotics + (1 | col), data = Q3_Pve_Frag, 
           family = binomial)
Linear Hypotheses:
  Estimate Std. Error z value Pr(>|z|)    
  biofilm_Trigger - biofilm_NoTrigger == 0   -0.4113     0.1266  -3.247  0.00466 ** 
  steril_NoTrigger - biofilm_NoTrigger == 0   0.2334     0.1133   2.061  0.11805    
  steril_Trigger - biofilm_NoTrigger == 0    -0.1814     0.1222  -1.485  0.15753    
  steril_NoTrigger - biofilm_Trigger == 0     0.6447     0.1225   5.261 8.58e-07 ***
  steril_Trigger - biofilm_Trigger == 0       0.2299     0.1308   1.758  0.15753    
  steril_Trigger - steril_NoTrigger == 0     -0.4148     0.1179  -3.518  0.00217 ** 
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'

# Fiber
Q3_Pve_Fib_react <- glmer(reaction ~ biotics + (1|col),
                           family = binomial, data = Q3_Pve_Fib)

summary(glht(Q3_Pve_Fib_react, linfct = mcp(biotics = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = reaction ~ biotics + (1 | col), data = Q3_Pve_Fib, 
           family = binomial)
Linear Hypotheses:
  Estimate Std. Error z value Pr(>|z|)    
  biofilm_Trigger - biofilm_NoTrigger == 0   -0.4335     0.1324  -3.275 0.004229 ** 
  steril_NoTrigger - biofilm_NoTrigger == 0   0.3516     0.1193   2.947 0.009614 ** 
  steril_Trigger - biofilm_NoTrigger == 0     0.1089     0.1231   0.885 0.376047    
  steril_NoTrigger - biofilm_Trigger == 0     0.7850     0.1261   6.226 2.86e-09 ***
  steril_Trigger - biofilm_Trigger == 0       0.5424     0.1296   4.184 0.000143 ***
  steril_Trigger - steril_NoTrigger == 0     -0.2426     0.1162  -2.088 0.073666 .  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'


# Ingestion
# Fragment
Q3_Pve_Frag_ingest <- glmer(ingestion ~ biotics + (1|col),
                           family = binomial, data = Q3_Pve_Frag)

summary(glht(Q3_Pve_Frag_ingest, linfct = mcp(biotics = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ingestion ~ biotics + (1 | col), data = Q3_Pve_Frag, 
           family = binomial)
Linear Hypotheses:
  Estimate Std. Error z value Pr(>|z|)    
biofilm_Trigger - biofilm_NoTrigger == 0    0.8322     0.4787   1.738   0.3286  
steril_NoTrigger - biofilm_NoTrigger == 0   0.9492     0.4212   2.254   0.1211  
steril_Trigger - biofilm_NoTrigger == 0     1.1920     0.4401   2.708   0.0406 *
steril_NoTrigger - biofilm_Trigger == 0     0.1170     0.3794   0.308   1.0000  
steril_Trigger - biofilm_Trigger == 0       0.3598     0.3947   0.911   1.0000  
steril_Trigger - steril_NoTrigger == 0      0.2428     0.3290   0.738   1.0000 
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'

# Fiber
Q3_Pve_Fib_ingest <- glmer(ingestion ~ biotics + (1|col),
                          family = binomial, data = Q3_Pve_Fib)

summary(glht(Q3_Pve_Fib_ingest, linfct = mcp(biotics = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ingestion ~ biotics + (1 | col), data = Q3_Pve_Fib, 
           family = binomial)
Linear Hypotheses:
  Estimate Std. Error z value Pr(>|z|)    
biofilm_Trigger - biofilm_NoTrigger == 0   1.44861    0.67974   2.131   0.1323  
steril_NoTrigger - biofilm_NoTrigger == 0 -1.56085    1.16647  -1.338   0.5426  
steril_Trigger - biofilm_NoTrigger == 0   -1.47388    1.17284  -1.257   0.5426  
steril_NoTrigger - biofilm_Trigger == 0   -3.00946    1.05603  -2.850   0.0262 *
steril_Trigger - biofilm_Trigger == 0     -2.92249    1.05236  -2.777   0.0274 *
steril_Trigger - steril_NoTrigger == 0     0.08696    1.42495   0.061   0.9513   
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'


# Incorporation
# Fragment
Q3_Pve_Frag_incorp <- glmer(incorporation ~ biotics + (1|col),
                            family = binomial, data = Q3_Pve_Frag)

summary(glht(Q3_Pve_Frag_incorp, linfct = mcp(biotics = "Tukey")), 
        test = adjusted("holm"))
# NA

# Fiber
Q3_Pve_Fib_incorp <- glmer(incorporation ~ biotics + (1|col),
                           family = binomial, data = Q3_Pve_Fib)
# is singular

summary(glht(Q3_Pve_Fib_incorp, linfct = mcp(biotics = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ingestion ~ biotics + (1 | col), data = Q3_Pve_Fib, 
           family = binomial)
Linear Hypotheses:
  Estimate Std. Error z value Pr(>|z|)    
biofilm_Trigger - biofilm_NoTrigger == 0  -1.163e+00  1.351e+00  -0.861        1
steril_NoTrigger - biofilm_NoTrigger == 0 -3.327e+01  1.188e+07   0.000        1
steril_Trigger - biofilm_NoTrigger == 0   -4.462e+01  6.711e+07   0.000        1
steril_NoTrigger - biofilm_Trigger == 0   -3.211e+01  1.188e+07   0.000        1
steril_Trigger - biofilm_Trigger == 0     -4.345e+01  6.711e+07   0.000        1
steril_Trigger - steril_NoTrigger == 0    -1.134e+01  6.815e+07   0.000        1
(Adjusted p values reported -- holm method)'



#### - S5.3.1. Time calculations -----------------------------------------------
# Fragment
Time_Q3_Pve_Frag <- Q3_Pve_Frag %>%
  # bring duration columns in time format hh:mm:ss
  mutate(reaction_dur = lubridate::hms(reaction_dur),
         ingestion_dur = lubridate::hms(ingestion_dur)) %>% 
  # bring durations in numeric to be able to calculate with
  mutate(reaction_dur_num = as.numeric(reaction_dur),
         ingestion_dur_num = as.numeric(ingestion_dur))

# calculate mean reaction time
mean_react_time_Q3_Pve_Frag <- Time_Q3_Pve_Frag %>%
  group_by(shape, biotics) %>%
  get_summary_stats(reaction_dur_num, type = "mean_sd")

# format mean and sd back into seconds
mean_react_time_Q3_Pve_Frag <- mean_react_time_Q3_Pve_Frag %>% 
  mutate(mean_react_dur = seconds_to_period(mean),
         sd_react_dur = seconds_to_period(sd))
mean_react_time_Q3_Pve_Frag

# calculate mean ingestion time
mean_ing_time_Q3_Pve_Frag <- Time_Q3_Pve_Frag %>%
  group_by(shape, biotics) %>%
  get_summary_stats(ingestion_dur_num, type = "mean_sd")
# format mean and sd back into seconds
mean_ing_time_Q3_Pve_Frag <- mean_ing_time_Q3_Pve_Frag %>% 
  mutate(mean_ing_dur = seconds_to_period(mean),
         sd_ing_dur = seconds_to_period(sd))

# testing whether reaction time differs
Time_Q3_Pve_Frag %>%
  wilcox_test(reaction_dur_num ~ biotics, data= .)
# OUTPUT:
# A tibble: 6 × 9
# .y.              group1            group2              n1    n2 statistic     p p.adj p.adj.signif
# * <chr>            <chr>             <chr>            <int> <int>     <dbl> <dbl> <dbl> <chr>       
# 1 reaction_dur_num biofilm_NoTrigger biofilm_Trigger    166   117     9660. 0.939     1 ns          
# 2 reaction_dur_num biofilm_NoTrigger steril_NoTrigger   166   211    18534. 0.31      1 ns          
# 3 reaction_dur_num biofilm_NoTrigger steril_Trigger     166   144    12988. 0.168     1 ns          
# 4 reaction_dur_num biofilm_Trigger   steril_NoTrigger   117   211    13092. 0.342     1 ns          
# 5 reaction_dur_num biofilm_Trigger   steril_Trigger     117   144     9215  0.172     1 ns          
# 6 reaction_dur_num steril_NoTrigger  steril_Trigger     211   144    15700. 0.574     1 ns   


# testing whether ingestion time differs
Time_Q3_Pve_Frag  %>%
  wilcox_test(ingestion_dur_num ~ biotics, data= .)
# OUTPUT: A tibble: 6 × 9
# .y.               group1            group2              n1    n2 statistic     p p.adj p.adj.signif
# * <chr>             <chr>             <chr>            <int> <int>     <dbl> <dbl> <dbl> <chr>       
#   1 ingestion_dur_num biofilm_NoTrigger biofilm_Trigger      8    12      44   0.785     1 ns          
# 2 ingestion_dur_num biofilm_NoTrigger steril_NoTrigger     8    24      73.5 0.334     1 ns          
# 3 ingestion_dur_num biofilm_NoTrigger steril_Trigger       8    21      88.5 0.844     1 ns          
# 4 ingestion_dur_num biofilm_Trigger   steril_NoTrigger    12    24     134.  0.761     1 ns          
# 5 ingestion_dur_num biofilm_Trigger   steril_Trigger      12    21     144.  0.496     1 ns          
# 6 ingestion_dur_num steril_NoTrigger  steril_Trigger      24    21     302   0.256     1 ns          


# Fiber
Time_Q3_Pve_Fib <- Q3_Pve_Fib %>%
  # bring duration columns in time format hh:mm:ss
  mutate(reaction_dur = lubridate::hms(reaction_dur),
         ingestion_dur = lubridate::hms(ingestion_dur)) %>% 
  # bring durations in numeric to be able to calculate with
  mutate(reaction_dur_num = as.numeric(reaction_dur),
         ingestion_dur_num = as.numeric(ingestion_dur))

# calculate mean reaction time
mean_react_time_Q3_Pve_Fib <- Time_Q3_Pve_Fib %>%
  group_by(shape, biotics) %>%
  get_summary_stats(reaction_dur_num, type = "mean_sd")

# format mean and sd back into seconds
mean_react_time_Q3_Pve_Fib <- mean_react_time_Q3_Pve_Fib %>% 
  mutate(mean_react_dur = seconds_to_period(mean),
         sd_react_dur = seconds_to_period(sd))
mean_react_time_Q3_Pve_Fib

# calculate mean ingestion time
mean_ing_time_Q3_Pve_Fib <- Time_Q3_Pve_Fib %>%
  group_by(shape, biotics) %>%
  get_summary_stats(ingestion_dur_num, type = "mean_sd")
# format mean and sd back into seconds
mean_ing_time_Q3_Pve_Fib <- mean_ing_time_Q3_Pve_Fib %>% 
  mutate(mean_ing_dur = seconds_to_period(mean),
         sd_ing_dur = seconds_to_period(sd))

# testing whether reaction time differs
Wil <- Time_Q3_Pve_Fib %>%
  wilcox_test(reaction_dur_num ~ biotics, data= .)
# OUTPUT:
# A tibble: 6 × 9
# .y.              group1            group2              n1    n2 statistic     p p.adj p.adj.signif
# * <chr>            <chr>             <chr>            <int> <int>     <dbl> <dbl> <dbl> <chr>       
# 1 reaction_dur_num biofilm_NoTrigger biofilm_Trigger    144   115     7798. 0.411     1 ns          
# 2 reaction_dur_num biofilm_NoTrigger steril_NoTrigger   144   207    15296. 0.664     1 ns          
# 3 reaction_dur_num biofilm_NoTrigger steril_Trigger     144   169    11938. 0.767     1 ns          
# 4 reaction_dur_num biofilm_Trigger   steril_NoTrigger   115   207    12920  0.189     1 ns          
# 5 reaction_dur_num biofilm_Trigger   steril_Trigger     115   169    10079  0.586     1 ns          
# 6 reaction_dur_num steril_NoTrigger  steril_Trigger     207   169    16680  0.423     1 ns   


# testing whether ingestion time differs
Time_Q3_Pve_Fib  %>%
  wilcox_test(ingestion_dur_num ~ biotics, data= .)
# OUTPUT: A tibble: 6 × 9
# .y.               group1            group2              n1    n2 statistic     p p.adj p.adj.signif
# * <chr>             <chr>             <chr>            <int> <int>     <dbl> <dbl> <dbl> <chr>       
# 1 ingestion_dur_num biofilm_NoTrigger biofilm_Trigger      3    13      31.5 0.114 0.684 ns          
# 2 ingestion_dur_num biofilm_NoTrigger steril_NoTrigger     3     1       2   1     1     ns          
# 3 ingestion_dur_num biofilm_NoTrigger steril_Trigger       3     1       3   0.346 1     ns          
# 4 ingestion_dur_num biofilm_Trigger   steril_NoTrigger    13     1       6   1     1     ns          
# 5 ingestion_dur_num biofilm_Trigger   steril_Trigger      13     1       9.5 0.517 1     ns          
# 6 ingestion_dur_num steril_NoTrigger  steril_Trigger       1     1       1   1     1     ns            


### --- 5.3.2. Stylophora pistillata -------------------------------------------
# Fragment
Q3_Spi_Frag <- Overall %>%
  filter(spec %in% c("Spi"),
         shape %in% c("fragment"),
         student %in% c ("Celine", "Verena", "Moritz", "Tim")) %>%
  unite(biotics, c(biofilm, trigger), sep = "_", remove = FALSE)

# Fiber
Q3_Spi_Fib <- Overall %>%
  filter(spec %in% c("Spi"),
         shape %in% c("fiber"),
         student %in% c ("Celine", "Verena", "Moritz", "Tim")) %>%
  unite(biotics, c(biofilm, trigger), sep = "_", remove = FALSE)



# Reaction
# Fragment
Q3_Spi_Frag_react <- glmer(reaction ~ biotics + (1|col),
                           family = binomial, data = Q3_Spi_Frag)

summary(glht(Q3_Spi_Frag_react, linfct = mcp(biotics = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = reaction ~ biotics + (1 | col), data = Q3_Pve_Frag, 
           family = binomial)
Linear Hypotheses:
  Estimate Std. Error z value Pr(>|z|)    
biofilm_Trigger - biofilm_NoTrigger == 0   0.01137    0.15039   0.076    1.000
steril_NoTrigger - biofilm_NoTrigger == 0  0.26343    0.14320   1.840    0.395
steril_Trigger - biofilm_NoTrigger == 0    0.16917    0.14713   1.150    1.000
steril_NoTrigger - biofilm_Trigger == 0    0.25207    0.14287   1.764    0.395
steril_Trigger - biofilm_Trigger == 0      0.15781    0.14681   1.075    1.000
steril_Trigger - steril_NoTrigger == 0    -0.09426    0.13942  -0.676    1.000
(Adjusted p values reported -- holm method)'

# Fiber
Q3_Spi_Fib_react <- glmer(reaction ~ biotics + (1|col),
                          family = binomial, data = Q3_Spi_Fib)

summary(glht(Q3_Spi_Fib_react, linfct = mcp(biotics = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = reaction ~ biotics + (1 | col), data = Q3_Pve_Fib, 
           family = binomial)
Linear Hypotheses:
  Estimate Std. Error z value Pr(>|z|)    
biofilm_Trigger - biofilm_NoTrigger == 0    0.1214     0.1463   0.830   0.8427  
steril_NoTrigger - biofilm_NoTrigger == 0   0.4034     0.1425   2.831   0.0279 *
steril_Trigger - biofilm_NoTrigger == 0     0.2584     0.1456   1.775   0.3037  
steril_NoTrigger - biofilm_Trigger == 0     0.2820     0.1353   2.084   0.1856  
steril_Trigger - biofilm_Trigger == 0       0.1369     0.1385   0.989   0.8427  
steril_Trigger - steril_NoTrigger == 0     -0.1451     0.1345  -1.078   0.8427 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'


# Ingestion
# Fragment
Q3_Spi_Frag_ingest <- glmer(ingestion ~ biotics + (1|col),
                            family = binomial, data = Q3_Spi_Frag)

summary(glht(Q3_Spi_Frag_ingest, linfct = mcp(biotics = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ingestion ~ biotics + (1 | col), data = Q3_Pve_Frag, 
           family = binomial)
Linear Hypotheses:
  Estimate Std. Error z value Pr(>|z|)    
biofilm_Trigger - biofilm_NoTrigger == 0  -0.02187    0.35056  -0.062    1.000
steril_NoTrigger - biofilm_NoTrigger == 0  0.52996    0.31248   1.696    0.476
steril_Trigger - biofilm_NoTrigger == 0    0.10232    0.33488   0.306    1.000
steril_NoTrigger - biofilm_Trigger == 0    0.55183    0.31459   1.754    0.476
steril_Trigger - biofilm_Trigger == 0      0.12419    0.33357   0.372    1.000
steril_Trigger - steril_NoTrigger == 0    -0.42764    0.29726  -1.439    0.601
(Adjusted p values reported -- holm method)'

# Fiber
Q3_Spi_Fib_ingest <- glmer(ingestion ~ biotics + (1|col),
                           family = binomial, data = Q3_Spi_Fib)

summary(glht(Q3_Spi_Fib_ingest, linfct = mcp(biotics = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ingestion ~ biotics + (1 | col), data = Q3_Pve_Fib, 
           family = binomial)
Linear Hypotheses:
  Estimate Std. Error z value Pr(>|z|)    
biofilm_Trigger - biofilm_NoTrigger == 0  -0.37643    0.65171  -0.578        1
steril_NoTrigger - biofilm_NoTrigger == 0 -0.06265    0.63375  -0.099        1
steril_Trigger - biofilm_NoTrigger == 0   -0.92602    0.77003  -1.203        1
steril_NoTrigger - biofilm_Trigger == 0    0.31379    0.59639   0.526        1
steril_Trigger - biofilm_Trigger == 0     -0.54959    0.73854  -0.744        1
steril_Trigger - steril_NoTrigger == 0    -0.86337    0.72395  -1.193        1
(Adjusted p values reported -- holm method)'


# Incorporation
# Fragment
Q3_Spi_Frag_incorp <- glmer(incorporation ~ biotics + (1|col),
                            family = binomial, data = Q3_Spi_Frag)

summary(glht(Q3_Spi_Frag_incorp, linfct = mcp(biotics = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = incorporation ~ biotics + (1 | col), data = Q3_Spi_Frag, 
           family = binomial)
Linear Hypotheses:
  Estimate Std. Error z value Pr(>|z|)
biofilm_Trigger - biofilm_NoTrigger == 0    1.0468     1.2386   0.845        1
steril_NoTrigger - biofilm_NoTrigger == 0   1.0953     0.9805   1.117        1
steril_Trigger - biofilm_NoTrigger == 0     1.3581     1.2132   1.119        1
steril_NoTrigger - biofilm_Trigger == 0     0.0485     1.2753   0.038        1
steril_Trigger - biofilm_Trigger == 0       0.3113     1.4636   0.213        1
steril_Trigger - steril_NoTrigger == 0      0.2628     1.2636   0.208        1
(Adjusted p values reported -- holm method)'

# Fiber
Q3_Spi_Fib_incorp <- glmer(incorporation ~ biotics + (1|col),
                           family = binomial, data = Q3_Spi_Fib)
# is singular

summary(glht(Q3_Spi_Fib_incorp, linfct = mcp(biotics = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ingestion ~ biotics + (1 | col), data = Q3_Pve_Fib, 
           family = binomial)
Linear Hypotheses:
  Estimate Std. Error z value Pr(>|z|)    
biofilm_Trigger - biofilm_NoTrigger == 0   4.842e+01  4.559e+05   0.000        1
steril_NoTrigger - biofilm_NoTrigger == 0  1.965e+01  1.308e+04   0.002        1
steril_Trigger - biofilm_NoTrigger == 0    1.987e+01  1.308e+04   0.002        1
steril_NoTrigger - biofilm_Trigger == 0   -2.877e+01  4.558e+05   0.000        1
steril_Trigger - biofilm_Trigger == 0     -2.854e+01  4.558e+05   0.000        1
steril_Trigger - steril_NoTrigger == 0     2.231e-01  1.483e+00   0.150        1
(Adjusted p values reported -- holm method)'



#### - S5.3.2. Time calculations -----------------------------------------------
# Fragment
Time_Q3_Spi_Frag <- Q3_Spi_Frag %>%
  # bring duration columns in time format hh:mm:ss
  mutate(reaction_dur = lubridate::hms(reaction_dur),
         ingestion_dur = lubridate::hms(ingestion_dur)) %>% 
  # bring durations in numeric to be able to calculate with
  mutate(reaction_dur_num = as.numeric(reaction_dur),
         ingestion_dur_num = as.numeric(ingestion_dur))

# calculate mean reaction time
mean_react_time_Q3_Spi_Frag <- Time_Q3_Spi_Frag %>%
  group_by(shape, biotics) %>%
  get_summary_stats(reaction_dur_num, type = "mean_sd")

# format mean and sd back into seconds
mean_react_time_Q3_Spi_Frag <- mean_react_time_Q3_Spi_Frag %>% 
  mutate(mean_react_dur = seconds_to_period(mean),
         sd_react_dur = seconds_to_period(sd))
mean_react_time_Q3_Spi_Frag

# calculate mean ingestion time
mean_ing_time_Q3_Spi_Frag <- Time_Q3_Spi_Frag %>%
  group_by(shape, biotics) %>%
  get_summary_stats(ingestion_dur_num, type = "mean_sd")
# format mean and sd back into seconds
mean_ing_time_Q3_Spi_Frag <- mean_ing_time_Q3_Spi_Frag %>% 
  mutate(mean_ing_dur = seconds_to_period(mean),
         sd_ing_dur = seconds_to_period(sd))

# testing whether reaction time differs
Time_Q3_Spi_Frag %>%
  wilcox_test(reaction_dur_num ~ biotics, data= .)
# OUTPUT:
# A tibble: 6 × 9
# .y.              group1            group2              n1    n2 statistic     p p.adj p.adj.signif
# * <chr>            <chr>             <chr>            <int> <int>     <dbl> <dbl> <dbl> <chr>       
# 1 reaction_dur_num biofilm_NoTrigger biofilm_Trigger    100    98     3460. 0.000189     0.000945    ***         
# 2 reaction_dur_num biofilm_NoTrigger steril_NoTrigger   100   124     6884. 0.108        0.142       ns          
# 3 reaction_dur_num biofilm_NoTrigger steril_Trigger     100   115     4976  0.069        0.142       ns          
# 4 reaction_dur_num biofilm_Trigger   steril_NoTrigger    98   124     8502. 0.0000000477 0.000000286 ****        
# 5 reaction_dur_num biofilm_Trigger   steril_Trigger      98   115     6496. 0.047        0.142       ns          
# 6 reaction_dur_num steril_NoTrigger  steril_Trigger     124   115     5435  0.000488     0.002       ** 


# testing whether ingestion time differs
Time_Q3_Spi_Frag  %>%
  wilcox_test(ingestion_dur_num ~ biotics, data= .)
# OUTPUT: A tibble: 6 × 9
# .y.               group1            group2              n1    n2 statistic     p p.adj p.adj.signif
# * <chr>             <chr>             <chr>            <int> <int>     <dbl> <dbl> <dbl> <chr>       
# 1 ingestion_dur_num biofilm_NoTrigger biofilm_Trigger     21    21      212. 0.828 0.981 ns          
# 2 ingestion_dur_num biofilm_NoTrigger steril_NoTrigger    21    40      300. 0.068 0.338 ns          
# 3 ingestion_dur_num biofilm_NoTrigger steril_Trigger      21    25      307  0.327 0.981 ns          
# 4 ingestion_dur_num biofilm_Trigger   steril_NoTrigger    21    40      332  0.18  0.72  ns          
# 5 ingestion_dur_num biofilm_Trigger   steril_Trigger      21    25      304. 0.368 0.981 ns          
# 6 ingestion_dur_num steril_NoTrigger  steril_Trigger      40    25      700  0.007 0.04  *  


# Fiber
Time_Q3_Spi_Fib <- Q3_Spi_Fib %>%
  # bring duration columns in time format hh:mm:ss
  mutate(reaction_dur = lubridate::hms(reaction_dur),
         ingestion_dur = lubridate::hms(ingestion_dur)) %>% 
  # bring durations in numeric to be able to calculate with
  mutate(reaction_dur_num = as.numeric(reaction_dur),
         ingestion_dur_num = as.numeric(ingestion_dur))

# calculate mean reaction time
mean_react_time_Q3_Spi_Fib <- Time_Q3_Spi_Fib %>%
  group_by(shape, biotics) %>%
  get_summary_stats(reaction_dur_num, type = "mean_sd")

# format mean and sd back into seconds
mean_react_time_Q3_Spi_Fib <- mean_react_time_Q3_Spi_Fib %>% 
  mutate(mean_react_dur = seconds_to_period(mean),
         sd_react_dur = seconds_to_period(sd))
mean_react_time_Q3_Spi_Fib

# calculate mean ingestion time
mean_ing_time_Q3_Spi_Fib <- Time_Q3_Spi_Fib %>%
  group_by(shape, biotics) %>%
  get_summary_stats(ingestion_dur_num, type = "mean_sd")
# format mean and sd back into seconds
mean_ing_time_Q3_Spi_Fib <- mean_ing_time_Q3_Spi_Fib %>% 
  mutate(mean_ing_dur = seconds_to_period(mean),
         sd_ing_dur = seconds_to_period(sd))

# testing whether reaction time differs
Time_Q3_Spi_Fib %>%
  wilcox_test(reaction_dur_num ~ biotics, data= .)
# OUTPUT:
# A tibble: 6 × 9
# .y.              group1            group2              n1    n2 statistic     p p.adj p.adj.signif
# * <chr>            <chr>             <chr>            <int> <int>     <dbl> <dbl> <dbl> <chr>       
# 1 reaction_dur_num biofilm_NoTrigger biofilm_Trigger     96   115     4983  0.213 0.816 ns          
# 2 reaction_dur_num biofilm_NoTrigger steril_NoTrigger    96   136     6417  0.821 1     ns          
# 3 reaction_dur_num biofilm_NoTrigger steril_Trigger      96   120     5130. 0.157 0.816 ns          
# 4 reaction_dur_num biofilm_Trigger   steril_NoTrigger   115   136     8528. 0.204 0.816 ns          
# 5 reaction_dur_num biofilm_Trigger   steril_Trigger     115   120     6740. 0.754 1     ns          
# 6 reaction_dur_num steril_NoTrigger  steril_Trigger     136   120     7304. 0.136 0.816 ns  


# testing whether ingestion time differs
Time_Q3_Spi_Fib  %>%
  wilcox_test(ingestion_dur_num ~ biotics, data= .)
# OUTPUT: A tibble: 6 × 9
# .y.               group1            group2              n1    n2 statistic     p p.adj p.adj.signif
# * <chr>             <chr>             <chr>            <int> <int>     <dbl> <dbl> <dbl> <chr>       
# 1 ingestion_dur_num biofilm_NoTrigger biofilm_Trigger      5     6      16   0.925 1     ns          
# 2 ingestion_dur_num biofilm_NoTrigger steril_NoTrigger     5     5      11.5 0.916 1     ns          
# 3 ingestion_dur_num biofilm_NoTrigger steril_Trigger       5     3       2.5 0.177 0.87  ns          
# 4 ingestion_dur_num biofilm_Trigger   steril_NoTrigger     6     5      13   0.772 1     ns          
# 5 ingestion_dur_num biofilm_Trigger   steril_Trigger       6     3       2   0.086 0.519 ns          
# 6 ingestion_dur_num steril_NoTrigger  steril_Trigger       5     3       2.5 0.174 0.87  ns   


## --- 5.4. Q4: What makes the difference? -------------------------------------
### -- 5.4.1. Type -------------------------------------------------------------
# Compare MP to natural food
Q4 <- Overall %>%
  filter(student %in% c("Celine", "Verena", "Moritz", "Tim"),
         particle != "Sand") %>%
  mutate(type = case_when(particle == "Artemia_cysts" ~ "Food",
                          particle == "Copepods" ~ "Food",
                          particle == "PE" ~ "MP",
                          particle == "PS" ~ "MP",
                          particle == "PVC" ~ "MP",
                          particle == "PA" ~ "MP",
                          particle == "PP" ~ "MP",
                          particle == "PET" ~ "MP"))

Q4$type <- as.factor(Q4$type)

Q4$type <- factor(Q4$type, 
                  levels = c("MP", "Food"))

# Reaction
Q4_react_type <- glmer(reaction ~ type + (1|col) + (1|particle) + (1|shape) + (1|biofilm) + (1|trigger) + (1|spec),
                      family = binomial, data = Q4)

EZ_1 <- summary(glht(Q4_react_type, linfct = mcp(type = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = reaction ~ type +  (1|col) + (1|particle) + (1|shape) + (1|biofilm) + (1|trigger) + (1|spec), 
                    data = Q4, family = binomial)
Linear Hypotheses:
  Estimate Std. Error z value Pr(>|z|)    
MP - Food == 0  1.6107     0.3445   4.676 2.92e-06 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_1$test$coefficients)
# OUTPUT: MP - Food 
#  5.006343 
log(exp(EZ_1$test$coefficients))
# 1.610706 


# Ingestion
Q4_ingest_type <- glmer(ingestion ~ type + (1|col) + (1|particle) + (1|shape) + (1|biofilm) + (1|trigger) + (1|spec),
                  family = binomial, data = Q4)

EZ_2 <- summary(glht(Q4_ingest_type, linfct = mcp(type = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
' Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ingestion ~ type + (1 | col) + (1 | spec) + (1 | 
    shape) + (1 | biofilm) + (1 | trigger), data = Q4, family = binomial)
Linear Hypotheses:
               Estimate Std. Error z value Pr(>|z|)    
MP - Food == 0   3.038      1.046   2.905  0.00367 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_2$test$coefficients)
# OUTPUT: MP - Food 
#  20.85625 
log(exp(EZ_2$test$coefficients))
#  3.037654 

# Incorporation
Q4_incorp_type <- glmer(incorporation ~ type + (1|col) + (1|particle) + (1|shape) + (1|biofilm) + (1|trigger) + (1|spec),
                        family = binomial, data = Q4)
# is singular
EZ_3 <- summary(glht(Q4_incorp_type, linfct = mcp(type = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
' Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = incorporation ~ type + (1 | col) + (1 | spec) + 
    (1 | shape) + (1 | biofilm) + (1 | trigger), data = Q4, family = binomial)
Linear Hypotheses:
               Estimate Std. Error z value Pr(>|z|)   
MP - Food == 0   3.740      1.918   1.949   0.0512 .
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_3$test$coefficients)
# OUTPUT: MP - Food 
#  42.09281 
log(exp(EZ_3$test$coefficients))
#  3.739877 

### -- 5.4.1.2. Polymer --------------------------------------------------------
# compare polymers, leave out sand and food, therefore leave out type as random factor
# Reaction
Q4_P <- Q4 %>%
  filter(type != "Food")
Q4_react_polymer <- glmer(reaction ~ particle + (1|col) + (1|shape) + (1|biofilm) + (1|trigger) + (1|spec),
                       family = binomial, data = Q4_P)
EZ_P <- summary(glht(Q4_react_polymer, linfct = mcp(particle = "Tukey")), 
        test = adjusted("holm"))


# OUTPUT:
'Multiple Comparisons of Means: Tukey Contrasts


Fit: glmer(formula = reaction ~ particle + (1 | col) + (1 | type) + 
    (1 | shape) + (1 | biofilm) + (1 | trigger) + (1 | spec), 
    data = Q4, family = binomial)

Linear Hypotheses:
                              Estimate Std. Error z value Pr(>|z|)    
PVC - PE == 0  -0.43946    0.08258  -5.322 9.26e-07 ***
PS - PE == 0    0.01687    0.07564   0.223 0.823522    
PA - PE == 0   -0.07543    0.07679  -0.982 0.683939    
PP - PE == 0    0.31433    0.07233   4.345 8.95e-05 ***
PET - PE == 0  -0.81579    0.08984  -9.080  < 2e-16 ***
PS - PVC == 0   0.45633    0.08236   5.541 3.01e-07 ***
PA - PVC == 0   0.36403    0.08342   4.364 8.95e-05 ***
PP - PVC == 0   0.75378    0.07934   9.500  < 2e-16 ***
PET - PVC == 0 -0.37633    0.09556  -3.938 0.000328 ***
PA - PS == 0   -0.09229    0.07656  -1.206 0.683939    
PP - PS == 0    0.29746    0.07209   4.126 0.000184 ***
PET - PS == 0  -0.83266    0.08965  -9.288  < 2e-16 ***
PP - PA == 0    0.38975    0.07329   5.318 9.26e-07 ***
PET - PA == 0  -0.74036    0.09061  -8.171 2.44e-15 ***
PET - PP == 0  -1.13011    0.08687 -13.009  < 2e-16 ***'

# Calculate Odds Ratio 
exp(EZ_P$test$coefficients)
log(exp(EZ_P$test$coefficients))


# Ingestion
Q4_ingest_type <- glmer(ingestion ~ particle + (1|col) + (1|shape) + (1|biofilm) + (1|trigger) + (1|spec),
                        family = binomial, data = Q4_P)

EZ_P2 <- summary(glht(Q4_ingest_type, linfct = mcp(particle = "Tukey")), 
                test = adjusted("holm"))
# OUTPUT:
'  Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ingestion ~ particle + (1 | col) + (1 | type) + 
    (1 | shape) + (1 | biofilm) + (1 | trigger) + (1 | spec), 
    data = Q4, family = binomial)
Linear Hypotheses:
                              Estimate Std. Error z value Pr(>|z|)    
PVC - PE == 0   -1.9494     0.3196  -6.101 1.37e-08 ***
PS - PE == 0    -0.8278     0.1893  -4.372 0.000135 ***
PA - PE == 0    -1.5940     0.2326  -6.854 1.01e-10 ***
PP - PE == 0    -2.8886     0.3273  -8.826  < 2e-16 ***
PET - PE == 0  -16.4795    51.2001  -0.322 1.000000    
PS - PVC == 0    1.1216     0.3325   3.373 0.006685 ** 
PA - PVC == 0    0.3555     0.3570   0.996 1.000000    
PP - PVC == 0   -0.9391     0.4265  -2.202 0.193766    
PET - PVC == 0 -14.5301    51.2012  -0.284 1.000000    
PA - PS == 0    -0.7661     0.2510  -3.052 0.018185 *  
PP - PS == 0    -2.0607     0.3408  -6.047 1.77e-08 ***
PET - PS == 0  -15.6517    51.2003  -0.306 1.000000    
PP - PA == 0    -1.2946     0.3641  -3.556 0.003771 ** 
PET - PA == 0  -14.8855    51.2005  -0.291 1.000000    
PET - PP == 0  -13.5909    51.2012  -0.265 1.000000  '

# Calculate Odds Ratio 
exp(EZ_P2$test$coefficients)
log(exp(EZ_P2$test$coefficients))


# Incorporation
Q4_incorp_particle <- glmer(incorporation ~ particle + (1|col) + (1|shape) + (1|biofilm) + (1|trigger) + (1|spec),
                        family = binomial, data = Q4_P)

EZ_P3 <- summary(glht(Q4_incorp_particle, linfct = mcp(particle = "Tukey")), 
                test = adjusted("holm"))
# OUTPUT:
' Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = incorporation ~ particle + (1 | col) + (1 | type) + 
    (1 | shape) + (1 | biofilm) + (1 | trigger) + (1 | spec), 
    data = Q4, family = binomial)
Linear Hypotheses:
                                Estimate Std. Error    z value Pr(>|z|)    
PVC - PE == 0  1.189e+05  1.937e+07   0.006  1.00000   
PS - PE == 0  -9.454e-01  6.929e-01  -1.364  1.00000   
PA - PE == 0  -2.440e+00  6.985e-01  -3.493  0.00477 **
PP - PE == 0  -5.311e+04  2.023e+07  -0.003  1.00000   
PS - PVC == 0 -1.189e+05  1.937e+07  -0.006  1.00000   
PA - PVC == 0 -1.189e+05  1.937e+07  -0.006  1.00000   
PP - PVC == 0 -1.720e+05  2.801e+07  -0.006  1.00000   
PA - PS == 0  -1.495e+00  6.446e-01  -2.319  0.18357   
PP - PS == 0  -5.311e+04  2.023e+07  -0.003  1.00000   
PP - PA == 0  -5.311e+04  2.023e+07  -0.003  1.00000 '

# Calculate Odds Ratio 
exp(EZ_P3$test$coefficients)
log(exp(EZ_P3$test$coefficients))



### -- 5.4.2. Shape ------------------------------------------------------------
# Compare MP shapes: Fragment vs. Fiber
Q4_2 <- Q4 %>%
  filter(shape %in% c("fragment", "fiber")) 
# Reaction
Q4_react_shape <- glmer(reaction ~ shape + (1|col) + (1|particle) + (1|biofilm) + (1|trigger) + (1|spec),
                       family = binomial, data = Q4_2)

EZ_4 <- summary(glht(Q4_react_shape, linfct = mcp(shape = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = reaction ~ shape + (1 | col) + (1 | particle) + 
    (1 | biofilm) + (1 | trigger) + (1 | spec), data = Q4_2, 
    family = binomial)
Linear Hypotheses:
                      Estimate Std. Error z value Pr(>|z|)
fragment - fiber == 0   0.04861    0.31727   0.153    0.878
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_4$test$coefficients)
# OUTPUT: fragment - fiber 
# 1.049811
log(exp(EZ_4$test$coefficients))
#  0.04861052

# Ingestion
Q4_ingest_shape <- glmer(ingestion ~ shape + (1|col) + (1|particle) + (1|biofilm) + (1|trigger) + (1|spec),
                        family = binomial, data = Q4_2)

EZ_5 <- summary(glht(Q4_ingest_shape, linfct = mcp(shape = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
' Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ingestion ~ shape + (1 | col) + (1 | particle) + 
    (1 | biofilm) + (1 | trigger) + (1 | spec), data = Q4_2, 
    family = binomial)
Linear Hypotheses:
                      Estimate Std. Error z value Pr(>|z|)  
fragment - fiber == 0   2.0539     0.9152   2.244   0.0248 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_5$test$coefficients)
# OUTPUT: fragment - fiber 
# 7.798579
log(exp(EZ_5$test$coefficients))
# 2.053942 

# Incorporation
Q4_incorp_shape <- glmer(incorporation ~ shape + (1|col) + (1|particle) + (1|biofilm) + (1|trigger) + (1|spec),
                        family = binomial, data = Q4_2)

EZ_6 <- summary(glht(Q4_incorp_shape, linfct = mcp(shape = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
' Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = incorporation ~ shape + (1 | col) + (1 | particle) + 
    (1 | biofilm) + (1 | trigger) + (1 | spec), data = Q4_2, 
    family = binomial)
Linear Hypotheses:
                      Estimate Std. Error z value Pr(>|z|)  
fragment - fiber == 0    4.953      2.208   2.243   0.0249 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_6$test$coefficients)
# OUTPUT: fragment - fiber 
# 7.798579
log(exp(EZ_6$test$coefficients))
# 4.953191 


### -- 5.4.3. Biofilm ----------------------------------------------------------
# Compare unfouled MP vs biofouled MP
Q4_3 <- Q4 %>%
  filter(biofilm != "blank")
# Reaction
Q4_react_biofilm <- glmer(reaction ~ biofilm + (1|col) + (1|particle) + (1|shape) + (1|trigger) + (1|spec),
                        family = binomial, data = Q4_3)

EZ_7 <- summary(glht(Q4_react_biofilm, linfct = mcp(biofilm = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = reaction ~ biofilm + (1 | col) + (1 | particle) + 
    (1 | shape) + (1 | trigger) + (1 | spec), data = Q4_3, family = binomial)
Linear Hypotheses:
                      Estimate Std. Error z value Pr(>|z|)    
steril - biofilm == 0   0.28279    0.04606    6.14 8.24e-10 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_7$test$coefficients)
# OUTPUT: steril - biofilm 
# 1.32683 
log(exp(EZ_7$test$coefficients))
#  0.2827923 

# Ingestion
Q4_ingest_biofilm <- glmer(ingestion ~ biofilm + (1|col) + (1|particle) + (1|shape) + (1|trigger) + (1|spec),
                         family = binomial, data = Q4_3)

EZ_8 <- summary(glht(Q4_ingest_biofilm, linfct = mcp(biofilm = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
' Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ingestion ~ biofilm + (1 | col) + (1 | particle) + 
    (1 | shape) + (1 | trigger) + (1 | spec), data = Q4_3, family = binomial)
Linear Hypotheses:
                      Estimate Std. Error z value Pr(>|z|)
steril - biofilm == 0   0.06334    0.15536   0.408    0.684
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_8$test$coefficients)
# OUTPUT: steril - biofilm 
# 1.065386 
log(exp(EZ_8$test$coefficients))
# 0.06333673 

# Incorporation
Q4_incorp_biofilm <- glmer(incorporation ~ biofilm + (1|col) + (1|particle) + (1|shape) + (1|trigger) + (1|spec),
                         family = binomial, data = Q4_3)

EZ_9 <- summary(glht(Q4_incorp_biofilm, linfct = mcp(biofilm = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
' Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
FFit: glmer(formula = incorporation ~ biofilm + (1 | col) + (1 | particle) + 
    (1 | shape) + (1 | trigger) + (1 | spec), data = Q4_3, family = binomial)
Linear Hypotheses:
                      Estimate Std. Error z value Pr(>|z|)
steril - biofilm == 0  -0.2604     0.5311   -0.49    0.624
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_9$test$coefficients)
# OUTPUT: steril - biofilm 
# 0.7707687  
log(exp(EZ_9$test$coefficients))
#   -0.2603669 


### -- 5.4.4. Trigger ----------------------------------------------------------
# Compare response in absence vs presense of a feeding trigger
# Reaction
Q4_react_trigger <- glmer(reaction ~ trigger + (1|col) + (1|particle) + (1|shape) + (1|biofilm) + (1|spec),
                          family = binomial, data = Q4)

EZ_10 <- summary(glht(Q4_react_trigger, linfct = mcp(trigger = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = reaction ~ trigger + (1 | col) + (1 | particle) + 
    (1 | shape) + (1 | biofilm) + (1 | spec), data = Q4, family = binomial)
Linear Hypotheses:
                         Estimate Std. Error z value Pr(>|z|)    
NoTrigger - Trigger == 0  0.24941    0.03966   6.288 3.21e-10 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_10$test$coefficients)
# OUTPUT: NoTrigger - Trigger 
# 1.283271
log(exp(EZ_10$test$coefficients))
# 0.249412 

# Ingestion
Q4_ingest_trigger <- glmer(ingestion ~ trigger + (1|col) + (1|particle) + (1|shape) + (1|biofilm) + (1|spec),
                           family = binomial, data = Q4)
# is singular
EZ_11 <- summary(glht(Q4_ingest_trigger, linfct = mcp(trigger = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
' Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ingestion ~ trigger + (1 | col) + (1 | particle) + 
    (1 | shape) + (1 | biofilm) + (1 | spec), data = Q4, family = binomial)
Linear Hypotheses:
                         Estimate Std. Error z value Pr(>|z|)
NoTrigger - Trigger == 0 -0.07243    0.09603  -0.754    0.451
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_11$test$coefficients)
# OUTPUT: NoTrigger - Trigger 
# 0.930134 
log(exp(EZ_11$test$coefficients))
# -0.07242666 

# Incorporation
Q4_incorp_trigger <- glmer(incorporation ~ trigger + (1|col) + (1|particle) + (1|shape) + (1|biofilm) + (1|spec),
                           family = binomial, data = Q4)

EZ_12 <- summary(glht(Q4_incorp_trigger, linfct = mcp(trigger = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
' Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = incorporation ~ trigger + (1 | col) + (1 | particle) + 
    (1 | shape) + (1 | biofilm) + (1 | spec), data = Q4, family = binomial)
Linear Hypotheses:
                         Estimate Std. Error z value Pr(>|z|)  
NoTrigger - Trigger == 0  -1.0154     0.4434   -2.29    0.022 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_12$test$coefficients)
# OUTPUT: NoTrigger - Trigger 
# 0.3622619 
log(exp(EZ_12$test$coefficients))
# -1.015388 


### -- 5.4.5. Species ----------------------------------------------------------
# Compare coral species: P. verrucosa vs. S. pistillata
Q4$spec <- factor(Q4$spec, 
                  levels = c("Spi", "Pve"))

# Reaction
Q4_react_spec <- glmer(reaction ~ spec + (1|col) + (1|particle) + (1|shape) + (1|trigger) + (1|biofilm),
                          family = binomial, data = Q4)

EZ_13 <- summary(glht(Q4_react_spec, linfct = mcp(spec = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = reaction ~ spec + (1 | col) + (1 | particle) + 
    (1 | shape) + (1 | trigger) + (1 | biofilm), data = Q4, family = binomial)
Linear Hypotheses:
               Estimate Std. Error z value Pr(>|z|)    
Pve - Spi == 0 0.39598    0.03985   9.938   <2e-16 ***
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_13$test$coefficients)
# OUTPUT: Pve - Spi
# 1.485843 
log(exp(EZ_13$test$coefficients))
# 0.3959824 

# Ingestion
Q4_ingest_spec <- glmer(ingestion ~ spec + (1|col) + (1|particle) + (1|shape) + (1|trigger) + (1|biofilm),
                           family = binomial, data = Q4)

EZ_14 <- summary(glht(Q4_ingest_spec, linfct = mcp(spec = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
' Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ingestion ~ spec + (1 | col) + (1 | particle) + 
    (1 | shape) + (1 | trigger) + (1 | biofilm), data = Q4, family = binomial)
Linear Hypotheses:
               Estimate Std. Error z value Pr(>|z|)    
Pve - Spi == 0  -0.9298     0.1025   -9.076   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_14$test$coefficients)
# OUTPUT: Pve - Spi
# 0.3946294 
log(exp(EZ_14$test$coefficients))
# 0.9298091 

# Incorporation
Q4_incorp_spec <- glmer(incorporation ~ spec + (1|col) + (1|particle) + (1|shape) + (1|trigger) + (1|biofilm),
                           family = binomial, data = Q4)

EZ_15 <- summary(glht(Q4_incorp_spec, linfct = mcp(spec = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
' Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = incorporation ~ spec + (1 | col) + (1 | particle) + 
    (1 | shape) + (1 | trigger) + (1 | biofilm), data = Q4, family = binomial)
Linear Hypotheses:
               Estimate Std. Error z value Pr(>|z|)
Pve - Spi == 0   -0.1776     0.4655   -0.381    0.703
(Adjusted p values reported -- holm method)'

# Calculate Odds Ratio 
exp(EZ_15$test$coefficients)
# OUTPUT: Pve - Spi 
# 0.8373105 
log(exp(EZ_15$test$coefficients))
# -0.1775604 


## --- 5.5. Associated molecules  ----------------------------------------------
molecules <- read_csv2("in/Molecules.csv")
data_wide <- molecules %>%
  unite("Polymer_Status", Polymer, State, sep = "_") %>%
  pivot_longer(cols = starts_with("C"), 
               names_to = "molecule", 
               values_to = "measured") %>%
  pivot_wider(names_from = Polymer_Status, values_from = measured, values_fill = list(measured = 0))

upset_data <- data_wide %>%
  # column_to_rownames("molecule") %>%  # Die Zutatenspalte als Zeilenbezeichner setzen
  as.data.frame()  # Stellen sicher, dass es ein DataFrame ist

library(vegan)


molecule_test2 <- upset_data %>%
  select(-total)

zutaten_matrix2 <- upset_data[, -1]  

for (poly in colnames(molecule_test2)[3:ncol(molecule_test2)]) {  # Spalten 3 bis Ende (Zutaten)
  cat("\nChi²-Test für", poly, ":\n")
  print(chisq.test(table(molecule_test2$molecule, molecule_test2[[poly]])))
}

molecule_test <- molecules %>%
  unite("Polymer_Status", Polymer, State, sep = "_") %>%
  select(-total)

zutaten_matrix <- molecule_test[, -1]  

for (zutat in colnames(molecule_test)[3:ncol(molecule_test)]) {  # Spalten 3 bis Ende (Zutaten)
  cat("\nChi²-Test für", zutat, ":\n")
  print(chisq.test(table(molecule_test$Polymer_Status, molecule_test[[zutat]])))
}


for_test <- molecules %>%
  pivot_longer(cols = starts_with("C"), 
               names_to = "molecule", 
               values_to = "measured")

# For adding to Matrix
test <- glmer(C8H9NO2 ~ Polymer + (1|State), 
              family = "binomial", data = molecules)

# get summary of GLMER
summary(glht(test, linfct = mcp(Polymer = "Tukey")), 
        test = adjusted("holm"))


## --- 5.6. Q5: Polyp activity under continous exposure ------------------------
# Compare MP + control food levels to MP + high food levels Overall
Q5 <- Polyps %>%
  filter(treat %in% c("MP", "MP+food"))

### -- 5.6.1. Pocillopora verrucosa --------------------------------------------
# create a subset with data of Pve
Q5_Pve <- subset(Q5, spec == "Pve") 

# LMER didn't show a good fit, therefore GLMER is used
Q5_Pve_test <- glmer(ranks ~ treat + (1|col) + (1|Time) + (1|week) + (1|foodnote), 
                     family = "poisson", data = Q5_Pve)

# get summary of GLMER
summary(glht(Q5_Pve_test, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ranks ~ treat + (1 | col) + (1 | Time) + (1 | 
    week) + (1 | foodnote), data = Q5_Pve, family = "poisson")
Linear Hypotheses:
                  Estimate Std. Error z value Pr(>|z|)
MP+food - MP == 0  0.02828    0.08990   0.315    0.753
(Adjusted p values reported -- holm method)'


### -- 5.6.2. Stylophora pistillata --------------------------------------------
# create a subset with data of Pve
Q5_Spi <- subset(Q5, spec == "Spi")

# LMER didn't show a good fit, therefore GLMER is used
Q5_Spi_test <- glmer(ranks ~ treat + (1|col) + (1|Time) + (1|week) + (1|foodnote), 
                     family = "poisson", data = Q5_Spi)

# get summary of GLMER
summary(glht(Q5_Spi_test, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
'Simultaneous Tests for General Linear Hypotheses
Multiple Comparisons of Means: Tukey Contrasts
Fit: glmer(formula = ranks ~ treat + (1 | col) + (1 | Time) + (1 | 
    week) + (1 | foodnote), data = Q5_Spi, family = "poisson")
Linear Hypotheses:
                  Estimate Std. Error z value Pr(>|z|)
MP+food - MP == 0 -0.09679    0.08321  -1.163    0.245
(Adjusted p values reported -- holm method)'




# ----- 6. Write tables --------------------------------------------------------
## ---- 6.1. Tables of summary -------------------------------------------------
# create a table to summarize the mean response rates
react_sum1 <- Overall %>%
  group_by(spec, particle, shape) %>%
  summarise(ones = sum(reaction, na.rm = TRUE),   # count 1
            n    = n(), .groups = "drop")
react_sum <- Overall %>%
  group_by(spec, particle, shape, biofilm, trigger) %>%
  summarise(ones = sum(reaction, na.rm = TRUE),   # count 1
            n    = n(), .groups = "drop")

ingest_sum1 <- Overall %>%
  group_by(spec, particle, shape) %>%
  summarise(ones = sum(ingestion, na.rm = TRUE),   # count 1
            n    = n(), .groups = "drop")
ingest_sum <- Overall %>%
  group_by(spec, particle, shape, biofilm, trigger) %>%
  summarise(ones = sum(ingestion, na.rm = TRUE),   # count 1
    n    = n(), .groups = "drop")

incorp_sum1 <- Overall %>%
  group_by(spec, particle, shape) %>%
  summarise(ones = sum(incorporation, na.rm = TRUE),   # count 1
            n    = n(), .groups = "drop")
incorp_sum <- Overall %>%
  group_by(spec, particle, shape, biofilm, trigger) %>%
  summarise(ones = sum(incorporation, na.rm = TRUE),   # count 1
            n    = n(), .groups = "drop")


# create a table to summarize the mean polyp activity for P. verrucosa
Pve <- subset(Polyps, spec == "Pve")
Pve_sum <- Pve %>%
  group_by(week, ID, treat, col) %>%
  get_summary_stats(ranks, type = "mean")


# create a table to summarize the mean polyp activity for S. pistillata
Spi <- subset(Polyps, spec == "Spi")
Spi_sum <- Spi %>%
  group_by(tp, ID, conc, col) %>%
  get_summary_stats(ranks, type = "mean")




