##########################################################################################################################################################################
########################################################################################################################################################################## 
# Developing a Discrete-Event Simulation of Early Cancer Diagnosis in Primary Care: A Tutorial Using Linked Health Records
# DES R CODE: MULTI-STATE MODELLING USING MSM PACKAGE
#
##########################################################################################################################################################################
##########################################################################################################################################################################

# Note: set up renv for final version

# Clear the workspace
rm(list = ls());


# PSA vs deterministic
Model_method <- 'PSA'
Model_method <- 'Deterministic'

# True random or fixed random
Model_random <- 'True_random'
Model_random <- 'Fixed_random'

# Stage method
Stage_shift <- FALSE
alpha <- 10

# Default standard error assumption for parameters with no reported uncertainty
default_se_prop <- 0.20



# Intervention
Is_intervention <- FALSE
Prob_intervention <- 1


########################
# 1. LOAD R PACKAGES
########################

library(simmer)
library(simmer.plot)
library(dplyr)
library(survival)
library(ggplot2)
library(MASS)
library(msm)
library(nnet)
library(SHELF)
library(tidyr)







