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


################################
# 2. EXTERNAL MODEL FUNCTIONS
################################

setwd("~/HEPD/Model Concep/CPRD/Murray/Final data/Control")
source("DES_model_functions2.R")



######################
# 3. LOAD DATA SETS
######################


# COSTS
Data_Costs <- read.csv("Data_Costs2.csv", header=TRUE)

# UTILITIES
Data_Utilities <- read.csv("Data_Utilities.csv", header=TRUE)



# FULL PATIENT INFORMATION

# Source - Linked dataset combining:
# - Clinical Practice Research Datalink (CPRD)
# - Hospital Episode Statistics (HES)
# - National Cancer Registry

# Contains patient-level characteristic data alongside event data including GP consultations, referrals, diagnoses

# Variables:
# - e_patid - Individual patient ID
# - Event - Current event state
# - next_event - next event state patient visits after current state
# - Time - Time which patient enters current state
# - AgeDiag - Age of patient at cancer diagnosis 
# - upci -
# - cms_burden - 
Data_Full_PatientsN <- read.csv("Data_Full_PatientsConN2.csv", header=TRUE)



# COSTS OF REFERRAL APPOINTMENTS based on treatment specialist codes

# Source - Derived from HES outpatient data and NHS Reference Costs 2024

# Contains patient-level referral information linked to:
# - referral priority (e.g. 2WW, non-2WW)
# - treatment specialty (tretspef, from HES)
# - cancer type (from cancer registry)
# - consultant-led outpatient appointment costs to associated tretspef code
PatientSpecCost <- read.csv("PatientSpecCost2.csv", header=TRUE)



# EMERGENCY COSTS

# Source - NHS Reference Costs 2024 - Emergency Care

# Contains national average unit costs for outpatient attendances, excluding those with no treatment 
# or dead on arrival to hospital
EmergencyCost <- read.csv("EmergencyCosts.csv", header=TRUE)



# OUTPATIENT APPOINTMENT COSTS

# Source - NHS Reference Costs 2024 - Outpatient Care

# Contains costs for outpatient appointments, excluding those relating to cancer or children
OutPatientAppCost <- read.csv("Outpatient App Costs.csv", header=TRUE) %>%
  mutate(Number.of.attendances = as.numeric(Number.of.attendances),
         National.average.unit.cost = as.numeric(National.average.unit.cost))



# OUTPATIENT PROCEDURE COSTS

# Source - NHS Reference Costs 2024 - Outpatient Procedures

# Contains costs for outpatient procedures, excluding those relating to cancer or children
OutPatientProCost <- read.csv("Outpatient Pro Costs.csv", header=TRUE) %>% filter(X == 'FALSE', Number.of.procedures != "*") %>%
  mutate(Number.of.procedures = as.numeric(Number.of.procedures),
         National.average.unit.cost = as.numeric(National.average.unit.cost))


# INPATIENT COSTS

# Source - NHS Reference Costs 2024 - Admitted Patient Care

# Contains costs for inpatient procedures, excluding those relating to cancer or children
InPatientCost <- read.csv("Inpatient Costs.csv", header=TRUE) %>% filter(X == 'FALSE', Number.of.finished.consultant.episodes != "*") %>%
  mutate(Number.of.finished.consultant.episodes = as.numeric(Number.of.finished.consultant.episodes),
         National.average.unit.cost = as.numeric(National.average.unit.cost))



# EXPERT ELICITATION

Q1LP <- read.csv("Q1LP.csv", header=TRUE)
Q2LP <- read.csv("Q2LP.csv", header=TRUE)
Q3LP <- read.csv("Q3LP.csv", header=TRUE)
Q4LP <- read.csv("Q4LP.csv", header=TRUE)


# RANDOM NUMBERS

Rand_set <- read.csv('Random_set.csv', header=TRUE)
U_samp <- read.csv('U_samp.csv', header=TRUE)
Patient_sample <- read.csv('Patient_sample.csv', header=TRUE)


# Cancer order
Cancer_order <- data.frame(
  breast = 1,
  colon = 2,
  lung = 3,
  prostate = 4,
  other = 5
)

###########################
# 4. COSTS AND UTILITIES 
###########################

# Discounting rate
Disc_rate <- Data_Costs %>% filter(Item == "Disc_rate") %>% pull(Cost)

# If before 2016, inflate to 2016, then inflate to 2024

# Inflation index from PSSRU 2024
nhs_index <- c(
  "2016" = 100,
  "2017" = 102.09,
  "2018" = 103.35,
  "2019" = 105.00,
  "2020" = 107.25,
  "2021" = 109.92,
  "2022" = 112.75,
  "2023" = 121.00,
  "2024" = 126.22
)

gdp_index <- c(
  "2001" = 100,
  "2002" = 102.1,
  "2003" = 104.55,
  "2004" = 107.37,
  "2005" = 110.49,
  "2006" = 113.91,
  "2007" = 116.19,
  "2008" = 120.03,
  "2009" = 122.43,
  "2010" = 124.02,
  "2011" = 127.12,
  "2012" = 129.02,
  "2013" = 131.86,
  "2014" = 133.97,
  "2015" = 134.91,
  "2016" = 137.20
)


# Apply to diagnosis and treatment costs
Treatment_costs <- Data_Costs %>% filter(Item == 'Treatment')
Treatment_costs$Year <- as.character(Treatment_costs$Year)

Treatment_costs$Inf_Costs <- mapply(
  inflate_cost,
  cost = Treatment_costs$Cost,
  from_year = Treatment_costs$Year
)

Treatment_costs$Inf_SE <- mapply(
  inflate_cost,
  cost = Treatment_costs$SE,
  from_year = Treatment_costs$Year
)

Diagnosis_costs <- Data_Costs %>% filter(Item == 'Diagnostic')

Diagnosis_costs$Year <- as.character(Diagnosis_costs$Year)

Diagnosis_costs$Inf_Costs <- mapply(
  inflate_cost,
  cost = Diagnosis_costs$Cost,
  from_year = Diagnosis_costs$Year
)


# Intervention cost
Cost_int <- Data_Costs %>% filter(Item == "Cost_int") %>% pull(Cost)



# FROM PSSRU 2024

# GP appointment (Additional 4.47 per minute)
GP_Cost <- Data_Costs %>% filter(Item == "GP_Cost") %>% pull(Cost)
Add_GP_cost <- Data_Costs %>% filter(Item == "Add_GP_cost") %>% pull(Cost)


# 2WW referral letter
TWWRef_cost <- Data_Costs %>% filter(Item == "TWWRef_cost") %>% pull(Cost)

# Non-2WW referral letter
N2WWRef_cost <- Data_Costs %>% filter(Item == "N2WWRef_cost") %>% pull(Cost)





# FROM NHS REFERENCE COSTS 2024

# 2WW appointment
#   Stored in PatientSpecCost
# Non-2WW appointment
#   Stored in PatientSpecCost


# Emergency costs # weighted average
# Emergency department costs
Emergency1 <- EmergencyCost %>% dplyr::select(Number.of.attendances, National.average.unit.cost) %>%
  na.omit() %>% mutate(Percent = Number.of.attendances/sum(Number.of.attendances), 
                       Cost = National.average.unit.cost * Percent)
Emerg_cost1 <- sum(Emergency1$Cost)

# Non-elective inpatient costs
Emergency2 <- InPatientCost %>% filter(Department.Description == "Non-Elective Inpatient - Long Stay" | Department.Description == "Non-Elective Inpatient - Short Stay") %>%
  dplyr::select(Number.of.finished.consultant.episodes, National.average.unit.cost) %>% na.omit() %>% 
  mutate(Percent = Number.of.finished.consultant.episodes/sum(Number.of.finished.consultant.episodes), 
         Cost = National.average.unit.cost * Percent)
Emerg_cost2 <- sum(Emergency2$Cost)

# Proportion of emergency department visits vs non-elective inpatient visits
# From data ~70% patients go through A&E
#           ~20% other including from GP
# Clinical advice
ProbEmerg1VSEmerg2 <- 0.9



# Out/inpatient costs # weighted average 
# Day case and elective inpatient costs
Inpatient1 <- InPatientCost %>% filter(Department.Description == "Daycase" | Department.Description == "Elective Inpatients") %>%
  dplyr::select(Number.of.finished.consultant.episodes, National.average.unit.cost) %>% na.omit() %>% 
  mutate(Percent = Number.of.finished.consultant.episodes/sum(Number.of.finished.consultant.episodes), 
         Cost = National.average.unit.cost * Percent)
Inpatient_cost1 <- sum(Inpatient1$Cost)

# Outpatient appointment costs
Outpatient1 <- OutPatientAppCost %>% dplyr::select(Number.of.attendances, National.average.unit.cost) %>% na.omit() %>% 
  mutate(Percent = Number.of.attendances/sum(Number.of.attendances), 
         Cost = National.average.unit.cost * Percent)
Outpatient_cost1 <- sum(Outpatient1$Cost)

# Outpatient procedure costs
Outpatient2 <- OutPatientProCost %>% dplyr::select(Number.of.procedures, National.average.unit.cost) %>% na.omit() %>% 
  mutate(Percent = Number.of.procedures/sum(Number.of.procedures), 
         Cost = National.average.unit.cost * Percent)
Outpatient_cost2 <- sum(Outpatient2$Cost)

# Proportion of outpatient vs inpatient diagnosis routes 
# Calculated from data
ProbOutVSIn <- 0.79

 
# Proportion of outpatient appointments vs procedures 
# Clinical advice
ProbOutAppVSOutPro <- 0.9





# Diagnosis cost - DUMMY FOR NOW
# Will be cancer type specific
Diagnosis_cost <- Data_Costs %>% filter(Item == "Diagnosis_cost") %>% pull(Cost)





# Utility
# Start from model
Utility_baseline_start <- 0.7863


# Life table for patients with no diagnosis
LifeTable <- read.csv('LifeTables.csv', header=TRUE)



##########################################################################################
# 5. LOAD REQUIRED MULTI-STATE MODELS AND DISTRIBUTIONS (or run DES_model_Additional.R)
##########################################################################################

# ADD IN REQUIRED INFO

# Model for the number of consecetive GP appointments attended
No_GP_distribution <- readRDS("No_GP_distribution.rds")
# Distribution for the number of subsequent referral appointments
Sub_app_distributions <- readRDS("Sub_app_distributions.rds")
# Model to predict cancer type
Model_cancer <- readRDS("Model_cancer.rds")
# Model to predict cancer staging
Model_stage <- readRDS("Model_stage.rds")

# Multi-state model for GP loop
fitGP <- readRDS("fitGP.rds")
# Multi-state model for 2WW referral pathways
fit2WW <- readRDS("fit2WW.rds")
# Multi-state model for non-2WW referral pathways
fitN2WW <- readRDS("fitN2WW.rds")

# Cancer distributions
Prob_type <- read.csv("Prob_type.csv", header=TRUE)

# Survival 
fit_surv <- readRDS("fit_survival.RDS")




#####################################
# Deterministic or PSA quantities
#####################################

# GP costs
GP_CostDP <- Fixed_PSA_Det(Mean = GP_Cost)
Add_GP_costDP <- Fixed_PSA_Det(Mean = Add_GP_cost)

# Referral letter costs
TWWRef_costDP <- Fixed_PSA_Det(Mean = TWWRef_cost)
N2WWRef_costDP <- Fixed_PSA_Det(Mean = N2WWRef_cost)

# Expert elicitation
PSANo <- 5555
if (Model_random == 'Fixed_random') {
  ExpertE_Q1 <- Expert_elic_DES(Density=Q1LP, PSA_no=10000, PSA_i=PSANo, U_samp=U_samp)
  ExpertE_Q2 <- Expert_elic_DES(Density=Q2LP, PSA_no=10000, PSA_i=PSANo, U_samp=U_samp)
  ExpertE_Q3 <- Expert_elic_DES(Density=Q3LP, PSA_no=10000, PSA_i=PSANo, U_samp=U_samp)
  ExpertE_Q4 <- Expert_elic_DES(Density=Q4LP, PSA_no=10000, PSA_i=PSANo, U_samp=U_samp)
} else {
  ExpertE_Q1 <- Expert_elic_DES(Density=Q1LP, PSA_no=10000, PSA_i=PSANo)
  ExpertE_Q2 <- Expert_elic_DES(Density=Q2LP, PSA_no=10000, PSA_i=PSANo)
  ExpertE_Q3 <- Expert_elic_DES(Density=Q3LP, PSA_no=10000, PSA_i=PSANo)
  ExpertE_Q4 <- Expert_elic_DES(Density=Q4LP, PSA_no=10000, PSA_i=PSANo)
}

# Intervention cost
Cost_intDP <- Fixed_PSA_Det(Mean = Cost_int)

# Referral appointment costs
PSC1 <- PatientSpecCost %>% dplyr::select(-e_patid) %>% distinct()
FirstDP1 <- data.frame(matrix(unlist(apply(as.matrix(PSC1$FirstCost), 1,  Fixed_PSA_Det)), ncol=2, byrow=TRUE))
PSC1$FirstDet <- FirstDP1$X1 
PSC1$FirstPSA <- FirstDP1$X2
FollowDP1 <- data.frame(matrix(unlist(apply(as.matrix(PSC1$FollUpCost), 1,  Fixed_PSA_Det)), ncol=2, byrow=TRUE))
PSC1$FollowDet <- FollowDP1$X1 
PSC1$FollowPSA <- FollowDP1$X2


PatientSpecCostDP <- PatientSpecCost %>% left_join(PSC1, by=c('priority','tretspef','Type','Consultant','FirstCost','FollUpCost'))


################################
# 6. DES PATHWAY TRAJECTORIES
################################



# TRAJECTORIES FROM DIAGNOSIS TO STAGING
# Calls make_stage_traj function for each diagnosis scenario

# GP to emergency diagnosis
traj_emerg_from_GP <- make_stage_traj(
  traj_name = "traj_emerg_from_GP",
  fit_model = fitGP,
  pathway = "EmergDiag",
  prefix = "GP_Emerg",
  Covariates_fun = function(Cov_pat) {
    data.frame(AgeDiag=Cov_pat$AgeDiag,
               upci=Cov_pat$upci,
               cms_burden=Cov_pat$cms_burden)},
  Diag_cost1 = ProbEmerg1VSEmerg2 * Emerg_cost1 + 
              (1-ProbEmerg1VSEmerg2) * Emerg_cost2,
  Diag_cost2 = NA,
  Treat_cost = Treatment_costs %>% dplyr::select(Cancer, Stage, Inf_Costs)
)

# 2WW referral to emergency diagnosis
traj_emerg_from_2WWRef <- make_stage_traj(
  traj_name = "traj_emerg_from_2WWRef",
  fit_model = fit2WW,
  pathway = "EmergDiag",
  prefix = "2WW_Emerg",
  Covariates_fun = function(Cov_pat) {
    data.frame(AgeDiag=Cov_pat$AgeDiag)},
  Diag_cost1 = ProbEmerg1VSEmerg2 * Emerg_cost1 + 
    (1-ProbEmerg1VSEmerg2) * Emerg_cost2,
  Diag_cost2 = NA,
  Treat_cost = Treatment_costs %>% dplyr::select(Cancer, Stage, Inf_Costs)
)

# 2WW appointment to emergency diagnosis
traj_emerg_from_2WWApp <- make_stage_traj(
  traj_name = "traj_emerg_from_2WWApp",
  fit_model = fit2WW,
  pathway = "EmergDiag",
  prefix = "2WW_Emerg",
  Covariates_fun = function(Cov_pat) {
    data.frame(AgeDiag=Cov_pat$AgeDiag)},
  Diag_cost1 = ProbEmerg1VSEmerg2 * Emerg_cost1 + 
    (1-ProbEmerg1VSEmerg2) * Emerg_cost2,
  Diag_cost2 = Diagnosis_costs %>% dplyr::select(Cancer, Stage, Distribution, Inf_Costs),
  Treat_cost = Treatment_costs %>% dplyr::select(Cancer, Stage, Inf_Costs)
)

# Non-2WW referral to emergency diagnosis
traj_emerg_from_N2WWRef <- make_stage_traj(
  traj_name = "traj_emerg_from_N2WWRef",
  fit_model = fitN2WW,
  pathway = "EmergDiag",
  prefix = "N2WW_Emerg",
  Covariates_fun = function(Cov_pat) {
    data.frame(AgeDiag=Cov_pat$AgeDiag)},
  Diag_cost1 = ProbEmerg1VSEmerg2 * Emerg_cost1 + 
    (1-ProbEmerg1VSEmerg2) * Emerg_cost2,
  Diag_cost2 = NA,
  Treat_cost = Treatment_costs %>% dplyr::select(Cancer, Stage, Inf_Costs)
)

# Non-2WW appointment to emergency diagnosis
traj_emerg_from_N2WWApp <- make_stage_traj(
  traj_name = "traj_emerg_from_N2WWApp",
  fit_model = fitN2WW,
  pathway = "EmergDiag",
  prefix = "N2WW_Emerg",
  Covariates_fun = function(Cov_pat) {
    data.frame(AgeDiag=Cov_pat$AgeDiag)},
  Diag_cost1 = ProbEmerg1VSEmerg2 * Emerg_cost1 + 
    (1-ProbEmerg1VSEmerg2) * Emerg_cost2,
  Diag_cost2 = Diagnosis_costs %>% dplyr::select(Cancer, Stage, Distribution, Inf_Costs),
  Treat_cost = Treatment_costs %>% dplyr::select(Cancer, Stage, Inf_Costs)
)

# GP to in/outpatient diagnosis
traj_InOut_from_GP <- make_stage_traj(
  traj_name = "traj_InOut_from_GP",
  fit_model = fitGP,
  pathway = "InOutDiag",
  prefix = "GP_InOut",
  Covariates_fun = function(Cov_pat) {
    data.frame(AgeDiag=Cov_pat$AgeDiag,
               upci=Cov_pat$upci,
               cms_burden=Cov_pat$cms_burden)},
  Diag_cost1 = ProbOutVSIn * (ProbOutAppVSOutPro * Outpatient_cost1 +
                               (1-ProbOutAppVSOutPro) * Outpatient_cost2) +
              (1-ProbOutVSIn) * Inpatient_cost1,
  Diag_cost2 = NA,
  Treat_cost = Treatment_costs %>% dplyr::select(Cancer, Stage, Inf_Costs)
)

# 2WW referral to in/outpatient diagnosis
traj_InOut_from_2WWRef <- make_stage_traj(
  traj_name = "traj_InOut_from_2WWRef",
  fit_model = fit2WW,
  pathway = "InOutDiag",
  prefix = "2WW_InOut",
  Covariates_fun = function(Cov_pat) {
    data.frame(AgeDiag=Cov_pat$AgeDiag)},
  Diag_cost1 = ProbOutVSIn * (ProbOutAppVSOutPro * Outpatient_cost1 +
                               (1-ProbOutAppVSOutPro) * Outpatient_cost2) +
              (1-ProbOutVSIn) * Inpatient_cost1,
  Diag_cost2 = NA,
  Treat_cost = Treatment_costs %>% dplyr::select(Cancer, Stage, Inf_Costs)
)

# 2WW appointment to in/outpatient diagnosis
traj_InOut_from_2WWApp <- make_stage_traj(
  traj_name = "traj_InOut_from_2WWApp",
  fit_model = fit2WW,
  pathway = "InOutDiag",
  prefix = "2WW_InOut",
  Covariates_fun = function(Cov_pat) {
    data.frame(AgeDiag=Cov_pat$AgeDiag)},
  Diag_cost1 = ProbOutVSIn * (ProbOutAppVSOutPro * Outpatient_cost1 +
                                (1-ProbOutAppVSOutPro) * Outpatient_cost2) +
    (1-ProbOutVSIn) * Inpatient_cost1,
  Diag_cost2 = Diagnosis_costs %>% dplyr::select(Cancer, Stage, Distribution, Inf_Costs),
  Treat_cost = Treatment_costs %>% dplyr::select(Cancer, Stage, Inf_Costs)
)

# Non-2WW referral to in/outpatient diagnosis
traj_InOut_from_N2WWRef <- make_stage_traj(
  traj_name = "traj_InOut_from_N2WWRef",
  fit_model = fitN2WW,
  pathway = "InOutDiag",
  prefix = "N2WW_InOut",
  Covariates_fun = function(Cov_pat) {
    data.frame(AgeDiag=Cov_pat$AgeDiag)},
  Diag_cost1 = ProbOutVSIn * (ProbOutAppVSOutPro * Outpatient_cost1 +
                               (1-ProbOutAppVSOutPro) * Outpatient_cost2) +
              (1-ProbOutVSIn) * Inpatient_cost1,
  Diag_cost2 = NA,
  Treat_cost = Treatment_costs %>% dplyr::select(Cancer, Stage, Inf_Costs)
)

# Non-2WW appointment to in/outpatient diagnosis
traj_InOut_from_N2WWApp <- make_stage_traj(
  traj_name = "traj_InOut_from_N2WWApp",
  fit_model = fitN2WW,
  pathway = "InOutDiag",
  prefix = "N2WW_InOut",
  Covariates_fun = function(Cov_pat) {
    data.frame(AgeDiag=Cov_pat$AgeDiag)},
  Diag_cost1 = ProbOutVSIn * (ProbOutAppVSOutPro * Outpatient_cost1 +
                                (1-ProbOutAppVSOutPro) * Outpatient_cost2) +
    (1-ProbOutVSIn) * Inpatient_cost1,
  Diag_cost2 = Diagnosis_costs %>% dplyr::select(Cancer, Stage, Distribution, Inf_Costs),
  Treat_cost = Treatment_costs %>% dplyr::select(Cancer, Stage, Inf_Costs)
)

# Diagnosis from 2WW appointment
traj_diag_from_2WWApp <- make_stage_traj(
  traj_name = "traj_diag_from_2WW",
  fit_model = fit2WW,
  pathway = "2WWDiag",
  prefix = "2WW_Diag",
  Covariates_fun = function(Cov_pat) {
    data.frame(AgeDiag=Cov_pat$AgeDiag)},
  Diag_cost1 = NA,
  Diag_cost2 = Diagnosis_costs %>% dplyr::select(Cancer, Stage, Distribution, Inf_Costs),
  Treat_cost = Treatment_costs %>% dplyr::select(Cancer, Stage, Inf_Costs)
)

# Diagnosis from non-2WW appointment
traj_diag_from_N2WWApp <- make_stage_traj(
  traj_name = "traj_diag_from_N2WW",
  fit_model = fitN2WW,
  pathway = "N2WWDiag",
  prefix = "N2WW_Diag",
  Covariates_fun = function(Cov_pat) {
    data.frame(AgeDiag=Cov_pat$AgeDiag)},
  Diag_cost1 = NA,
  Diag_cost2 = Diagnosis_costs %>% dplyr::select(Cancer, Stage, Distribution, Inf_Costs),
  Treat_cost = Treatment_costs %>% dplyr::select(Cancer, Stage, Inf_Costs)
)


#set_attribute(keys = "total_utility", values = function() {
# Data_Utilities %>% filter(Item == 'AnxDiagnosis') %>% pull(Utility)}, mod = "+") %>%



##################################
# 2WW APPOINTMENT TRAJECTORY


traj_2WW_attend <- trajectory(name = 'traj_2WW_attend') %>%
  
  # INITIALISE ATTRIBUTES

  # Sample time spent in 2WW appointment state and next branch
  set_attribute(keys = c("time_in_2WWApp","branch_2WWApp", "Prob_2WWAppNC"), value = function() {
    Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"), 
                                  Data=Data_Full_PatientsN)
    
    if (Model_random == "Fixed_random") {
      u_time <-  Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                     j == get_attribute(.env = simmodel, keys = "2WW_loop_visits")) %>% 
        dplyr::pull(TWWApp_time)
      u_branch <-  Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                       j == get_attribute(.env = simmodel, keys = "2WW_loop_visits")) %>% 
        dplyr::pull(TWWApp_branch)
    } else{
      u_time <- NULL
      u_branch <- runif(1)
    }
    
    Pred <- simulate_next(current_state="2WWApp", 
                          Fit=fit2WW, 
                          Covariates = data.frame(
                            AgeDiag=Cov_pat$AgeDiag),
                          u_time = u_time)

    state_order <- c("2WWDiag", "EmergDiag", "InOutDiag", "Consultation")
    Pred$Probs <- Pred$Probs[state_order]
    
    if (Model_random == "Fixed_random") {
      branch <- which(u_branch <= cumsum(Pred$Probs))[1]
    } else {
      branch <- sample(seq_along(state_order), 1, prob = Pred$Probs)
    }
    Time_branch <- Pred$Time
    
    c(Pred$Time, branch, Pred$Probs["Consultation"])
  }) %>%
  
  # Add referral appointment costs given number of subsequent appointments, cancer
  # type and treatment specialist associated with the cancer type
  set_attribute(keys = "total_cost", values = function() {
    Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"), 
                                  Data=Data_Full_PatientsN)
    
    if (Model_random == "Fixed_random") {
      u_spec <- Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                    j == get_attribute(.env = simmodel, keys = "2WW_loop_visits")) %>% 
        dplyr::pull(TWW_spec)
      u_apps <-Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                   j == get_attribute(.env = simmodel, keys = "2WW_loop_visits")) %>% 
        dplyr::pull(TWW_Apps)
    } else {
      u_spec <- NULL
      u_apps <- NULL
    }
    
    Cost_apps1 <- cost_ref_appts(Fit=Sub_app_distributions$fit_2WW, 
                                 Priority="2WWApp", 
                                 Cancer=names(Cancer_order[get_attribute(.env = simmodel, keys = 'Cancer_Type')]),
                                 Cost_data=PatientSpecCostDP,
                                 u_spec = u_spec,
                                 u_apps = u_apps)

    
    No2WWs <- Cost_apps1$NoApps 
    
    # Discount costs separately in case total time exceeds one year
    if (No2WWs <= 1) {
      TWW_cost_first <- Disc_cost(cost = Cost_apps1$CostFirst, time_total = simmer::now(.env = simmodel))
      TWW_cost_follow <- No2WWs * Disc_cost(cost = Cost_apps1$CostFollow, time_total = (simmer::now(.env = simmodel) + get_attribute(.env = simmodel, keys = "time_in_2WWApp")))
    } else {
      TWW_appt_times <- seq((simmer::now(.env = simmodel)), 
                            (simmer::now(.env = simmodel) + get_attribute(.env = simmodel, keys = "time_in_2WWApp")), 
                            length.out = No2WWs)
      TWW_cost_first <- Disc_cost(cost = Cost_apps1$CostFirst, time_total = TWW_appt_times[1])
      TWW_cost_follow <- sum(unlist(sapply(TWW_appt_times[-1], function(t) Disc_cost(cost = Cost_apps1$CostFollow, time_total = t))))
    }
    TWW_cost_total <- TWW_cost_first + TWW_cost_follow
    
    # Diagnostic costs
    Type <- names(Cancer_order[get_attribute(.env = simmodel, keys = 'Cancer_Type')])
    
    if (Type == 'breast') {
      Diag_costs1 <- Diagnosis_costs %>% filter(Cancer == Type, Stage == 'Appointment')
      Diag_costs <- sample(Diag_costs1$Inf_Costs, 1, prob=Diag_costs1$Distribution)
    } else {
      Diag_costs <- sum(Diagnosis_costs %>% filter(Cancer == Type, Stage == 'Appointment') %>% pull(Inf_Costs))
    }
    Diag_costs <- Disc_cost(cost = Diag_costs, time_total = simmer::now(.env = simmodel)) # Discounting
    
    as.numeric(TWW_cost_total + Diag_costs)
  }, mod = "+") %>%
  
  
  # Time delay in current state
  timeout_from_attribute(key = "time_in_2WWApp") %>%
  
  # Discounted QALYs
  set_attribute(keys = "total_QALYs", values = function() {
    Tutility <- get_attribute(.env = simmodel, keys = "total_utility") 
    Ttime <- get_attribute(.env = simmodel, keys = "time_in_2WWApp") 
    DQALYs <- Disc_QALYs(utility=Tutility, time_total=simmer::now(.env = simmodel), current_time=Ttime, discount_rate=Disc_rate)
    
    as.numeric(DQALYs)
  }, mod = "+") %>%
  
  # BRANCH 
  
  branch(option = function() get_attribute(.env = simmodel, keys = "branch_2WWApp"), continue = c(T,T,T,T),
         
         # 1. 2WW diagnosis
         traj_diag_from_2WWApp,
         
         # 2. Emergency diagnosis
         traj_emerg_from_2WWApp,
         
         # 3. In/outpatient diagnosis
         traj_InOut_from_2WWApp,
         
         # 4. Return to GP (rollback)
         trajectory(name = "2WW_GP") %>%
           set_attribute(keys = "total_cost", values = function() {
             # Diagnostic costs
             Type <- names(Cancer_order[get_attribute(.env = simmodel, keys = 'Cancer_Type')])
            
             Diag_costs1 <- Diagnosis_costs %>% filter(Cancer == Type, Stage == 'No diagnosis')
             Prob <- (Diag_costs1$Distribution/100)/get_attribute(.env = simmodel, keys = "Prob_2WWAppNC")
             Diag_costs <- sample(c(Diag_costs1$Inf_Costs,0), 1, prob=c(Prob,0)) 
             Diag_costs <- Disc_cost(cost = Diag_costs, time_total = simmer::now(.env = simmodel)) # Discounting
             
             as.numeric(Diag_costs)
           }, mod = "+") %>%
           simmer::rollback(target = "visit_start", times = Inf)
  )




#####################################
# NON-2WW APPOINTMENT TRAJECTORY


traj_N2WW_attend <- trajectory(name = 'traj_N2WW_attend') %>%
  
  # INITIALISE ATTRIBUTES

  # Predict cancer type referral is for (only for appointment costing, not diagnosis)
  set_attribute(keys = 'Cancer_Type', values = function() {
    Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"), 
                                  Data=Data_Full_PatientsN)
    Cancer_type_pred(Covariates=Cov_pat, Event='N2WWDiag')
  }) %>%
  
  # Sample time spent in non-2WW appointment state and next branch
  set_attribute(keys = c("time_in_N2WWApp","branch_N2WWApp","Prob_N2WWAppNC"), value = function() {
    Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"), 
                                  Data=Data_Full_PatientsN)
    
    if (Model_random == "Fixed_random") {
      u_time <-  Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                     j == get_attribute(.env = simmodel, keys = "N2WW_loop_visits")) %>% 
        dplyr::pull(N2WWApp_time)
      u_branch <-  Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                       j == get_attribute(.env = simmodel, keys = "N2WW_loop_visits")) %>% 
        dplyr::pull(N2WWApp_branch)
    } else{
      u_time <- NULL
      u_branch <- runif(1)
    }
    
    Pred <- simulate_next(current_state="N2WWApp", 
                          Fit=fitN2WW, 
                          Covariates = data.frame(
                            AgeDiag=Cov_pat$AgeDiag),
                          u_time = u_time)
    
    state_order <- c("N2WWDiag", "EmergDiag", "InOutDiag", "Consultation")
    Pred$Probs <- Pred$Probs[state_order]
    
    if (Model_random == "Fixed_random") {
      branch <- which(u_branch <= cumsum(Pred$Probs))[1]
    } else {
      branch <- sample(seq_along(state_order), 1, prob = Pred$Probs)
    }
    Time_branch <- Pred$Time
    
    c(Pred$Time, branch, Pred$Probs["Consultation"])
  }) %>%
  
  # Add referral appointment costs given number of subsequent appointments, cancer
  # type and treatment specialist associated with the cancer type
  set_attribute(keys = "total_cost", values = function() {
    Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"), 
                                  Data=Data_Full_PatientsN)
    
    if (Model_random == "Fixed_random") {
      u_spec <- Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                    j == get_attribute(.env = simmodel, keys = "N2WW_loop_visits")) %>% 
        dplyr::pull(N2WW_spec)
      u_apps <- Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                    j == get_attribute(.env = simmodel, keys = "N2WW_loop_visits")) %>% 
        dplyr::pull(N2WW_Apps)
    } else {
      u_spec <- NULL
      u_apps <- NULL
    }
    
    Cost_apps1 <- cost_ref_appts(Fit=Sub_app_distributions$fit_N2WW, 
                                 Priority="N2WWApp", 
                                 Cancer=names(Cancer_order[get_attribute(.env = simmodel, keys = 'Cancer_Type')]),
                                 Cost_data=PatientSpecCostDP,
                                 u_spec = u_spec,
                                 u_apps = u_apps)
    
    NoN2WWs <- Cost_apps1$NoApps 
    
    # Discount costs separately in case total time exceeds one year
    if (NoN2WWs <= 1) {
      N2WW_cost_first <- Disc_cost(cost = Cost_apps1$CostFirst, time_total = simmer::now(.env = simmodel))
      N2WW_cost_follow <- NoN2WWs * Disc_cost(cost = Cost_apps1$CostFollow, time_total = (simmer::now(.env = simmodel) + get_attribute(.env = simmodel, keys = "time_in_N2WWApp")))
    } else {
      N2WW_appt_times <- seq((simmer::now(.env = simmodel)), 
                            (simmer::now(.env = simmodel) + get_attribute(.env = simmodel, keys = "time_in_N2WWApp")), 
                            length.out = NoN2WWs)
      N2WW_cost_first <- Disc_cost(cost = Cost_apps1$CostFirst, time_total = N2WW_appt_times[1])
      N2WW_cost_follow <- sum(unlist(sapply(N2WW_appt_times[-1], function(t) Disc_cost(cost = Cost_apps1$CostFollow, time_total = t))))
    }
    N2WW_cost_total <- N2WW_cost_first + N2WW_cost_follow
    
    as.numeric(N2WW_cost_total)
  }, mod = "+") %>%
  
  # Time delay in current state
  timeout_from_attribute(key = "time_in_N2WWApp") %>%
  
  # Discounted QALYs
  set_attribute(keys = "total_QALYs", values = function() {
    Tutility <- get_attribute(.env = simmodel, keys = "total_utility") 
    Ttime <- get_attribute(.env = simmodel, keys = "time_in_N2WWApp") 
    DQALYs <- Disc_QALYs(utility=Tutility, time_total=simmer::now(.env = simmodel), current_time=Ttime, discount_rate=Disc_rate)
    
    as.numeric(DQALYs)
  }, mod = "+") %>%
  
  # BRANCH
  
  branch(option = function() get_attribute(.env = simmodel, keys = "branch_N2WWApp"), continue = c(T,T,T,T),
         
         # 1. Non2WW diagnosis
         traj_diag_from_N2WWApp,
         
         # 2. Emergency diagnosis
         traj_emerg_from_N2WWApp,
         
         # 3. In/outpatient diagnosis
         traj_InOut_from_N2WWApp,
         
         # 4. Return to GP (rollback)
         trajectory(name = "N2WW_GP") %>%
           
           set_attribute(keys = "total_cost", values = function() {
             # Diagnostic costs
             Type <- names(Cancer_order[get_attribute(.env = simmodel, keys = 'Cancer_Type')])
             
             Diag_costs1 <- Diagnosis_costs %>% filter(Cancer == Type, Stage == 'No diagnosis')
             Prob <- (Diag_costs1$Distribution/100)/get_attribute(.env = simmodel, keys = "Prob_N2WWAppNC")
             Diag_costs <- sample(c(Diag_costs1$Inf_Costs,0), 1, prob=c(Prob,0)) 
             Diag_costs <- Disc_cost(cost = Diag_costs, time_total = simmer::now(.env = simmodel)) # Discounting
             
             as.numeric(Diag_costs)
           }, mod = "+") %>%
           
           simmer::rollback(target = "visit_start", times = Inf)
  )



#################################
# 2WW REFERRAL TRAJECTORY


traj_2WW <- trajectory(name = 'traj_2WW') %>%
  
  # INITIALISE ATTRIBUTES
  
  set_attribute(keys = "2WW_loop_visits", values = 1, mod = "+") %>%
  
  # Add in cost for referral letter
  set_attribute(keys = c("total_cost"), values = function() {

    Ref_cost <- (if (Model_method == 'Deterministic') TWWRef_costDP$Det else TWWRef_costDP$PSA)
    Disc_cost(cost = Ref_cost, time_total = simmer::now(.env = simmodel)) # Discounting
    }, mod = "+") %>%
  
  # Sample time spent in 2WW referral state and next branch
  set_attribute(keys = c("time_in_2WWRef","branch_2WWRef"), value = function() {
    Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"), 
                                  Data=Data_Full_PatientsN)
    
    if (Model_random == "Fixed_random") {
      u_time <-  Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                     j == get_attribute(.env = simmodel, keys = "2WW_loop_visits")) %>% 
        dplyr::pull(TWWRef_time)
      u_branch <- Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                       j == get_attribute(.env = simmodel, keys = "2WW_loop_visits")) %>% 
        dplyr::pull(TWWRef_branch)
    } else{
      u_time <- NULL
      u_branch <- runif(1)
    }
    
    Pred <- simulate_next(current_state="2WWRef",
                          Fit=fit2WW,
                          Covariates = data.frame(
                            AgeDiag=Cov_pat$AgeDiag),
                          u_time = u_time)

    state_order <- c("2WWApp", "EmergDiag", "InOutDiag")
    Pred$Probs <- Pred$Probs[state_order]
    
    if (Model_random == "Fixed_random") {
      branch <- which(u_branch <= cumsum(Pred$Probs))[1]
    } else {
      branch <- sample(seq_along(state_order), 1, prob = Pred$Probs)
    }
    Time_branch <- Pred$Time

    c(Time_branch, branch)
  }) %>%
  
  # Predict cancer type referral is for (only for pre referral diagnostic costing, not diagnosis)
  set_attribute(keys = 'Cancer_Type', values = function() {
    Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"), 
                                  Data=Data_Full_PatientsN)
    Type <- Cancer_type_pred(Covariates=Cov_pat, Event='2WWDiag')
    if (Type == 5) {
      Type <- sample(1:4, 1, prob=Prob_type$Freq)
    }
    as.numeric(Type)
  }) %>%
  
  # Add in cost for pre referral diagnostic costing
  set_attribute(keys = c("total_cost"), values = function() {
    
    Type <- names(Cancer_order[get_attribute(.env = simmodel, keys = 'Cancer_Type')])
    if (Type == 'breast') {
      Diag_costs <- 0
    } else {
      Diag_costs <- Diagnosis_costs %>% filter(Cancer == Type, Stage == 'Referral') %>% pull(Inf_Costs)
      Diag_costs <- Disc_cost(cost = Diag_costs, time_total = simmer::now(.env = simmodel)) # Discounting
    }
    as.numeric(Diag_costs)
  }, mod = "+") %>%
  
  # Time delay in current state
  timeout_from_attribute(key = "time_in_2WWRef") %>%
  
  # Discounted QALYs
  set_attribute(keys = "total_QALYs", values = function() {
    Tutility <- get_attribute(.env = simmodel, keys = "total_utility") 
    Ttime <- get_attribute(.env = simmodel, keys = "time_in_2WWRef") 
    Disc_QALYs(utility=Tutility, time_total=simmer::now(.env = simmodel), current_time=Ttime, discount_rate=Disc_rate)
  }, mod = "+") %>%
  
  # BRANCH
  
  branch(option = function() get_attribute(.env = simmodel, keys = "branch_2WWRef"), continue = c(T,T,T),
         
         # 1. Attend 2WW appointment
         traj_2WW_attend,
         
         # 2. Emergency diagnosis
         traj_emerg_from_2WWRef,
         
         # 3. In/outpatient diagnosis
         traj_InOut_from_2WWRef
  )







#################################
# NON-2WW REFERRAL TRAJECTORY


traj_non2WW <- trajectory(name = 'traj_Non2WW') %>%
  
  # INITIALISE ATTRIBUTES
  
  set_attribute(keys = "N2WW_loop_visits", values = 1, mod = "+") %>%
  
  # Add in cost for referral letter
  set_attribute(keys = c("total_cost"), values = function() {

    Ref_cost <- (if (Model_method == 'Deterministic') N2WWRef_costDP$Det else N2WWRef_costDP$PSA)
    Disc_cost(cost = Ref_cost, time_total = simmer::now(.env = simmodel)) # Discounting
    }, mod = "+") %>%
  
  # Sample time spent in non-2WW referral state and next branch
  set_attribute(keys = c("time_in_N2WWRef","branch_N2WWRef"), value = function() {
    Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"), 
                                  Data=Data_Full_PatientsN)
    
    if (Model_random == "Fixed_random") {
      u_time <-  Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                     j == get_attribute(.env = simmodel, keys = "N2WW_loop_visits")) %>% 
        dplyr::pull(N2WWRef_time)
      u_branch <-  Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                       j == get_attribute(.env = simmodel, keys = "N2WW_loop_visits")) %>% 
        dplyr::pull(N2WWRef_branch)
    } else{
      u_time <- NULL
      u_branch <- runif(1)
    }
    
    Pred <- simulate_next(current_state="N2WWRef", 
                          Fit=fitN2WW, 
                          Covariates = data.frame(
                            AgeDiag=Cov_pat$AgeDiag),
                          u_time = u_time)
    
    state_order <- c("N2WWApp", "EmergDiag", "InOutDiag")
    Pred$Probs <- Pred$Probs[state_order]
    
    if (Model_random == "Fixed_random") {
      branch <- which(u_branch <= cumsum(Pred$Probs))[1]
    } else {
      branch <- sample(seq_along(state_order), 1, prob = Pred$Probs)
    }
    Time_branch <- Pred$Time

    c(Time_branch, branch)
  }) %>%
  
  # Predict cancer type referral is for (only for pre referral diagnostic costing, not diagnosis)
  set_attribute(keys = 'Cancer_Type', values = function() {
    Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"), 
                                  Data=Data_Full_PatientsN)
    Type <- Cancer_type_pred(Covariates=Cov_pat, Event='2WWDiag')
    if (Type == 5) {
      Type <- sample(1:4, 1, prob=Prob_type$Freq)
    }
    as.numeric(Type)
  }) %>%
  
  # Add in cost for pre referral diagnostic costing
  set_attribute(keys = c("total_cost"), values = function() {
    
    Type <- names(Cancer_order[get_attribute(.env = simmodel, keys = 'Cancer_Type')])
    if (Type == 'breast') {
      Diag_costs <- 0
    } else {
      Diag_costs <- Diagnosis_costs %>% filter(Cancer == Type, Stage == 'Referral') %>% pull(Inf_Costs)
      Diag_costs <- Disc_cost(cost = Diag_costs, time_total = simmer::now(.env = simmodel)) # Discounting
    }
    as.numeric(Diag_costs)
  }, mod = "+") %>%
  
  # Time delay in current state
  timeout_from_attribute(key = "time_in_N2WWRef") %>%
  
  # Discounted QALYs
  set_attribute(keys = "total_QALYs", values = function() {
    Tutility <- get_attribute(.env = simmodel, keys = "total_utility") 
    Ttime <- get_attribute(.env = simmodel, keys = "time_in_N2WWRef") 
    Disc_QALYs(utility=Tutility, time_total=simmer::now(.env = simmodel), current_time=Ttime, discount_rate=Disc_rate)
  }, mod = "+") %>%
  
  # BRANCH
  
  branch(option = function() get_attribute(.env = simmodel, keys = "branch_N2WWRef"), continue = c(T,T,T),
         
         # 1. Attend N2WW appointment
         traj_N2WW_attend,
         
         # 2. Emergency diagnosis
         traj_emerg_from_N2WWRef,
         
         # 3. In/outpatient diagnosis
         traj_InOut_from_N2WWRef
  )






################################
# INITIAL GP LOOP TRAJECTORY


Full_Model_traj <- trajectory(name = 'traj') %>%
  
  # INITIALISE ATTRIBUTES

  # Sample patient id at model entry
  set_attribute(keys = "Patient_id", values = function() {
    
    if (Model_random == "Fixed_random") {
      sim_name <- get_name(.env=simmodel)
      sim_index <- as.integer(sub("patient", "", sim_name)) + 1
      Samp_patient <- Patient_sample[sim_index]
    } else {
      Samp_patient <- Sample_patient(Patient_data = Data_Full_PatientsN)
    }
    as.numeric(Samp_patient)
    })%>%
  
  # Initialise cumulative outcomes
  # Total running cost and QALYs for patient
  set_attribute(keys = c("total_cost","total_utility","total_QALYs"), values = function() c(0,Utility_baseline_start,0)) %>%
  
  # Set utility for anxiety and/or depression
  set_attribute(keys = "total_utility", values = function() {
    Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"), 
                                  Data=Data_Full_PatientsN)
    if (Cov_pat$AorD == 'Anxiety') {
      Data_Utilities %>% filter(Item == 'Anxiety') %>% pull(Utility)
    } else if (Cov_pat$AorD == 'Depression') {
      Data_Utilities %>% filter(Item == 'Depression') %>% pull(Utility)
    } else {
      Data_Utilities %>% filter(Item == 'AandD') %>% pull(Utility)
    }
  }, mod = "*") %>%
  
  set_attribute(keys = "total_utility", values = Data_Utilities %>% filter(Item == 'AnxDiagnosis') %>% pull(Utility), mod = "+") %>%
  
  # Initialise state loop visit totals
  set_attribute(keys = c("GP_loop_visits", "2WW_loop_visits", "N2WW_loop_visits"), values = c(0,0,0)) %>%
  
  # Initialise GP return visit counter and rollback flag
  set_attribute(keys = "GP_loop_visits", values = 1, mod = "+", tag = "visit_start") %>%
  set_attribute(keys = "GP_n_visits", values = 0) %>%
  
    
  # Does the patient receive the intervention?
  set_attribute(keys = c("Intervention", "Int_cost"), values = function() {
    if (!Is_intervention) return(c(0, 0))
    
    if (Model_random == "Fixed_random") {
      u <- Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"),
                               j == get_attribute(.env = simmodel, keys = "GP_loop_visits")) %>% 
                        dplyr::pull(Intervention)
      IsInt <- as.integer(u < Prob_intervention)
    } else {
      IsInt <- rbinom(1, 1, Prob_intervention)
    }
    
    expert_val <- if (Model_method == 'Deterministic') ExpertE_Q1$DET else ExpertE_Q1$PSA 
    Add_GP_costDP2 <- if (Model_method == 'Deterministic') Add_GP_costDP$Det else Add_GP_costDP$PSA
    Cost_intDP2 <- if (Model_method == 'Deterministic') Cost_intDP$Det else Cost_intDP$PSA 
    Add_cost <- IsInt * ((expert_val * Add_GP_costDP2) + Cost_intDP2)
    #Add_cost <- Disc_cost(cost=Add_cost, time_total=now(.env = simmodel)) # Depends how intervention cost is formed
    
    c(IsInt, Add_cost)
  }) %>%
  
  # Sample time in recurring GP appointments and branch to next state
  set_attribute(keys = c("time_in_GP","branch_GP"), value = function() {
    Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"), 
                                  Data=Data_Full_PatientsN)
    
    if (Model_random == "Fixed_random") {
      u_time <-  Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                     j == get_attribute(.env = simmodel, keys = "GP_loop_visits")) %>% 
        dplyr::pull(GP_time)
      u_branch <-  Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                       j == get_attribute(.env = simmodel, keys = "GP_loop_visits")) %>% 
        dplyr::pull(GP_branch)
    } else{
      u_time <- NULL
      u_branch <- runif(1)
    }
    
    Pred <- simulate_next(
      current_state = 'Consultation', 
      Fit = fitGP, 
      Covariates = data.frame(
        AgeDiag=Cov_pat$AgeDiag,
        upci=Cov_pat$upci,
        cms_burden=Cov_pat$cms_burden),
      u_time = u_time)
    
    state_order <- c("2WWRef", "N2WWRef", "EmergDiag", "InOutDiag", "No_diagnosis")
    Pred$Probs <- Pred$Probs[state_order]

    if (get_attribute(.env = simmodel, keys = "Intervention") == 1) {
      if (Intervention_affect) {
        Add_2WWRef <- min(alpha/100, (1-Pred$Probs["2WWRef"]-Pred$Probs["No_diagnosis"]))

      } else {
        Add_2WWRef <- min((if (Model_method == 'Deterministic') ExpertE_Q2$DET else ExpertE_Q2$PSA)/100, (1-Pred$Probs["2WWRef"]-Pred$Probs["No_diagnosis"]))
        Add_EmergR <- min((if (Model_method == 'Deterministic') ExpertE_Q3$DET else ExpertE_Q3$PSA)/100, (1-Pred$Probs["EmergDiag"]-Pred$Probs["No_diagnosis"]))
        
        Pred$Probs["EmergDiag"] <- Pred$Probs["EmergDiag"] + Add_EmergR
        Scale_E <- (1-Pred$Probs["EmergDiag"]-Pred$Probs["No_diagnosis"]) / (Pred$Probs["N2WWRef"] + Pred$Probs["InOutDiag"] + Pred$Probs["2WWRef"])
        Pred$Probs["N2WWRef"] <- Pred$Probs["N2WWRef"] * Scale_E
        Pred$Probs["EmergDiag"] <- Pred$Probs["EmergDiag"] * Scale_E
        Pred$Probs["InOutDiag"] <- Pred$Probs["InOutDiag"] * Scale_E
      }
      
      Pred$Probs["2WWRef"] <- Pred$Probs["2WWRef"] + Add_2WWRef
      Scale_2 <- (1-Pred$Probs["2WWRef"]-Pred$Probs["No_diagnosis"]) / (Pred$Probs["N2WWRef"] + Pred$Probs["InOutDiag"] + Pred$Probs["EmergDiag"])
      Pred$Probs["N2WWRef"] <- Pred$Probs["N2WWRef"] * Scale_2
      Pred$Probs["EmergDiag"] <- Pred$Probs["EmergDiag"] * Scale_2
      Pred$Probs["InOutDiag"] <- Pred$Probs["InOutDiag"] * Scale_2
      

    }
    
    if (Model_random == "Fixed_random") {
      branch <- which(u_branch <= cumsum(Pred$Probs))[1]
    } else {
      branch <- sample(seq_along(state_order), 1, prob = Pred$Probs)
    }
    Time_branch <- Pred$Time
    
    c(Time_branch, branch)
  }) %>%
  

  # Time out for time spent attending gp appointments
  timeout_from_attribute(key = "time_in_GP") %>%
  
  # Add in number of GP appointments and total cost for attending
  set_attribute(keys = c("GP_n_visits", "total_cost"), values = function() {
    
    GPCost <- (if (Model_method == 'Deterministic') GP_CostDP$Det else GP_CostDP$PSA)
    
    if (Model_random == "Fixed_random") {
      u_count <- Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                      j == get_attribute(.env = simmodel, keys = "GP_loop_visits")) %>% 
        dplyr::pull(GP_NoApps)
    } else {
      u_count <- NULL
    }
    
    CostGPApps <- cost_gp_appts(t_sim = get_attribute(.env = simmodel, keys = "time_in_GP"),
                  Fit = No_GP_distribution, 
                  GPcost = GPCost,
                  u_count = u_count)
    
    NoGPs <- CostGPApps$NoGPs
    CostGPs1 <- CostGPApps$CostGPs
    CostGPs2 <- get_attribute(.env = simmodel, keys = "Int_cost") # Additional intervention cost
    
    CostGPs <- CostGPs1 + CostGPs2 
    
    # Discount costs separately in case total time exceeds one year
    gp_appt_times <- seq((simmer::now(.env = simmodel) - get_attribute(.env = simmodel, keys = "time_in_GP")), 
                         simmer::now(.env = simmodel), 
                         length.out = NoGPs)
    gp_cost_total <- sum(sapply(gp_appt_times, function(t) Disc_cost(cost = CostGPs, time_total = t)))

    c(NoGPs, gp_cost_total)
    }, mod = "+") %>%
  
  # Total QALYs for time in GP appointments
  # include discounting
  set_attribute(keys = "total_QALYs", values = function() {
    Tutility <- get_attribute(.env = simmodel, keys = "total_utility") 
    Ttime <- get_attribute(.env = simmodel, keys = "time_in_GP") 
    Disc_QALYs(utility=Tutility, time_total=simmer::now(.env = simmodel), current_time=Ttime, discount_rate=Disc_rate)
  }, mod = "+") %>%
  
  
  # BRANCH
  
  branch(option = function() get_attribute(.env = simmodel, keys = "branch_GP"), continue = c(T,T,T,T,T),
         
         # 1. 2WW referral
         traj_2WW,
         
         # 2. Non-2WW referral
         traj_non2WW,
         
         # 3. Emergency diagnosis
         traj_emerg_from_GP,

         # 4. In/outpatient diagnosis
         traj_InOut_from_GP,
         
         # 5. No diagnosis - pathway ends
         trajectory("No_diagnosis") %>%
           set_attribute(keys = "no_diagnosis", values = 1) %>%
           
           set_attribute(keys = "total_utility", values = -(Data_Utilities %>% filter(Item == 'AnxDiagnosis') %>% pull(Utility)), mod = "+") %>%
           # Discounted QALYs - LIFETIME
           set_attribute(keys = "total_QALYs", values = function() {
             Tutility <- get_attribute(.env = simmodel, keys = "total_utility") 
             Ttime <- simmer::now(.env = simmodel)
             Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"),
                                           Data=Data_Full_PatientsN)
             if (Cov_pat$gender == 1) {
               Surv_age <- LifeTable$Male_Surv
             } else {
               Surv_age <- LifeTable$Female_Surv
             }
             calc_no_cancer_qaly(age_start=Cov_pat$AgeDiag, endTime=Ttime, utility=Tutility, life_table=Surv_age)
           }, mod = "+")
  )



########################################
# DEFINE THE SIMULATION ENVIRONMENT
########################################

simmodel <- simmer(name = "simmodel") %>%
  add_generator(
    name_prefix = "patient",
    trajectory = Full_Model_traj,
    distribution = at(rep(x = 0, times = 10)),
    mon = 2)


# RUN THE SIMULATION
simmodel %>% reset() %>% run()


############
# PROBLEMS

# Data_Full_PatientsN$AorD[is.na(Data_Full_PatientsN$AorD)] <- 'both'
# Check everything
# Possible error in diagnosis costings



##################################
# Simplify results & validation

arrivals <- get_mon_arrivals(simmodel)
head(arrivals)
arrivals
atts3 <- get_mon_attributes(simmodel)
atts3

atts3 %>% filter(name == "patient1")

# Useful keys
atts3_results <- atts3 %>% filter(key %in% c('Patient_id','total_cost', 'no_diagnosis',
                                             'branch_GP','branch_2WWRef','branch_N2WWRef',"Diag_GP_Emerg","Diag_GP_InOut",
                                             'branch_2WWApp',"Diag_2WW_Emerg","Diag_2WW_InOut",
                                             'branch_N2WWApp',"Diag_N2WW_Emerg","Diag_N2WW_InOut",
                                             "Diag_2WW_Diag","Diag_N2WW_Diag",
                                             "Stage1","Stage2","Stage3","Stage4"))

# Just keep maximum cost for each patient
atts3_results <- atts3_results %>%
  group_by(name, key) %>%
  arrange(name, time) %>%
  filter(
    key != "total_cost" | value == max(value)
  ) %>%
  ungroup()

atts3_results %>% print(n=200)



sim_results <- atts3_results %>%
  count(key) %>%
  #group_by(key) %>%
  mutate(prop = n / sum(n))

data_results <- Data_Full_PatientsN %>%
  count(Event, next_event) %>%
  group_by(Event) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


# GP loop
sim_results %>% filter(key %in% c('branch_2WWRef','branch_N2WWRef',"Diag_GP_Emerg","Diag_GP_InOut")) %>% mutate(prop = n / sum(n))
data_results %>% filter(Event == "Consultation", next_event != "Consultation") %>% mutate(prop = n / sum(n))

  
# 2WW 
sim_results %>% filter(key %in% c('branch_2WWApp',"Diag_2WW_Emerg","Diag_2WW_InOut")) %>% mutate(prop = n / sum(n))
data_results %>% filter((Event == "2WWRef" | Event == "2WWApp") & !(next_event == "Consultation" | next_event == "2WWDiag")) %>% 
  dplyr::select(-Event) %>% group_by(next_event) %>%
  summarise(n = sum(n), .groups = "drop") %>% mutate(prop = n / sum(n))


# non-2WW 
sim_results %>% filter(key %in% c('branch_N2WWApp',"Diag_N2WW_Emerg","Diag_N2WW_InOut")) %>% mutate(prop = n / sum(n))
data_results %>% filter((Event == "N2WWRef" | Event == "N2WWApp") & !(next_event == "Consultation" | next_event == "N2WWDiag")) %>%   
  dplyr::select(-Event) %>% group_by(next_event) %>%
  summarise(n = sum(n), .groups = "drop") %>% mutate(prop = n / sum(n))


# Staging
sim_results %>% filter(key %in% c('Stage1','Stage2','Stage3','Stage4')) %>% mutate(prop = n / sum(n))
data_results %>% filter(next_event %in% c('Stage_1', 'Stage_2', 'Stage_3', 'Stage_4')) %>% 
  dplyr::select(-Event) %>% group_by(next_event) %>%
  summarise(n = sum(n), .groups = "drop") %>% mutate(prop = n / sum(n))




