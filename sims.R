###-----------------------------------------------------------------------------###
###-----------------------------------------------------------------------------###
###-----------------------------------------------------------------------------###
### Robust variance estimation in small meta-analysis with the SMD
###
### This script includes functions and simulations that were used to assess the 
### properties that are reported in the main article
###
###-----------------------------------------------------------------------------###
###-----------------------------------------------------------------------------###
###-----------------------------------------------------------------------------###

# Load libraries
library(tidyverse)
library(readxl)

##-----------------------------------------------##
## Sourced Scripts
##-----------------------------------------------##
source("df-HC2.R")
source("df-HC3.R")
source("funcs.R")


#############################################
# TAKES: n_values;  An excel file containing sample sizes used
#        delta;     true effect, hard-coded value 
#                   (results reported in the paper contain delta = 0, 0.5, 1)
#        A;         hard-coded value from Hartung, 1999
#        B;         hard-coded value from Hartung, 1999
#        kValue;    the number of studies, hard-coded value 
#                   (results reported in the paper contain kValue = 2, 4, 6, 8, 10)
# RETURNS: 4*kValue by 30 dataframe
#############################################

# Read in excel file containing n values 
# Readers interested in running the simulations need to change the current directory accordingly
n_values <- read_excel("n_values.xlsx", col_names = FALSE) 
delta <- 0.5
A <- 0.95
B <- 1.05
kValue <- nrow(n_values) # number of studies
tau_sq <- c(0, 0.5, 1, 2) # degree of heterogeneity

# Run 10,000 iterations

df_sims <- data.frame()
for (l in 1:10000){
    currentResults <- main()
    df_sims <- rbind(df_sims, currentResults) 
}

names(df_sims) <-
  c("kValue", "currentTau", "j",
    "df", "df_hc2_est", "df_hc3_est",
    "deltaHat_val", "tauHat_value",
    "varDL_deltaHat", "varHK_deltaHat", "qofBeta", "hc1", "hc2", "hc3",
    "is_in_CI_DL", "is_in_CI_HK", "is_in_CI_Hartung",
    "is_in_CI_hc1", "is_in_CI_hc2", "is_in_CI_hc3",
    "is_in_CI_hc2_newDOF_est", "is_in_CI_hc3_newDOF_est",
    "DL_CI_length", "HK_CI_length", "Hartung_CI_length",
    "hc1_CI_length", "hc2_CI_length", "hc3_CI_length",
    "hc2_newDOF_est_CI_length", "hc3_newDOF_est_CI_length")

head(df_sims)
