##########################################################################################################################################################################
##########################################################################################################################################################################
#
# Developing a Discrete-Event Simulation of Early Cancer Diagnosis in Primary Care: A Tutorial Using Linked Health Records
# DES R CODE: FUNCTIONS FILE
#
##########################################################################################################################################################################
##########################################################################################################################################################################




# 1. INFLATE_COST

# Function to inflate given cost to given year
# @param cost Cost to inflate
# @param from_year Year of current cost
# @param to_year Year of requested inflation cost

# @return Inflated cost at given 'to_year'

inflate_cost <- function(cost, from_year, nhs_idx=nhs_index, gdp_idx=gdp_index, to_year = "2024") {
  
  if (as.numeric(from_year) < 2016) {
    Cost2016 <- cost * (gdp_idx["2016"] / gdp_idx[from_year])
    cost <- Cost2016
    from_year <- '2016'
  } 
  
  Cost2024 <- cost * (nhs_idx[to_year] / nhs_idx[from_year])
  
  return(Cost2024)
}





# 2. SIMULATE_NEXT

# Simulate the probabilities of transitioning to each of the next states and the 
# time to that state
#
# @param current_state Current state from Consultation, 2WWRef, N2WWRef, 2WWApp, N2WWApp, 2WWDiag, 
#                       N2WWDiag, EmergDiag, InOutDiag, Stage_1, Stage_2, Stage_3, Stage_4
# @param Fit A fitted 'msm' model object
# @param Covariates Named list of patient-level covariate values passed to qmatrix.msm()
#
# @return Named list containing simulated waiting time until next transition and
#         transition probabilities from current state to all possible states

simulate_next <- function(current_state, Fit, Covariates=list(Time=0), u_time = NULL) {
  
  # Extract transition intensity matrix
  Qhat <- qmatrix.msm(Fit, ci = "none", covariates=Covariates)
  
  # Filter matrix to current state
  rates <- Qhat[current_state, ]
  rates[current_state] <- 0
  
  # Rate of leaving current state
  lambda <- sum(rates)
  
  # In case of absorbing state
  if(lambda == 0){
    return(list(time = Inf, next_state = current_state))
  }
  
  # Sample waiting time - time until the next transition
  if (is.null(u_time)) {
    Time <- rexp(1, lambda)
  } else {
    Time <- -log(1 - u_time) / lambda
  }
  # Distribution over destination states, given that a transition happens, from the current state
  # Convert transition intensities to conditional jump probabilities:
  # P(next state = j | leaving current state i) = q_ij / (-q_ii)
  Probs = rates / lambda
  
  return(list(Time = Time, Probs = Probs))
}





# 3. SAMPLE_PATIENT

# Sample a patient from the observed patient data by first sampling a cancer type
# according to its observed frequency, then samples one patient at random from
# within that cancer type.

# @param patient_data A data frame containing patient covariates - at least the variables
# 'Type' and 'e_patid'

# @return Named list with the sampled cancer type and sampled patient id number -
# 'e_patid'

Sample_patient <- function(Patient_data) {
  
  # Restrict to patients with the sampled cancer type
  # Sample one patient from within the selected cancer type and return
  Samp_patient <- Patient_data %>% 
    slice_sample(n = 1)
  
  return(Patient_id=Samp_patient$e_patid)
}




# 4. COST_GP_APPTS

# Simulate number of GP appointments visited in a given simulated time 't_sim' 
# using a fitted model (negative binomial) 
# Associated cost is also returned

# @param t_sim Simulated time patient consistently visits GP
# @param fit A fitted negative binomial model object compatible with 'predict()'
# @param GPcost Cost per GP appointment

# @return A named list with simulated number of GP appointments and total cost of 
# simulated GP appointments

cost_gp_appts <- function(t_sim, Fit, GPcost, u_count=NULL) {
  # Predict expected number of GP appointments for simulated time
  Mu <- predict(Fit,
                newdata = data.frame(Time = t_sim),
                type = "response")
  # Simulate realised appointment count
  # Enforce a minimum of one GP appointment
  if (is.null(u_count)) {
    SampGP <- rnbinom(1, size = Fit$theta, mu = Mu)
  } else {
    u_count <- min(max(u_count, 1e-12), 1 - 1e-12)
    SampGP <- qnbinom(u_count, size = Fit$theta, mu = Mu)
  }
  
  SampGP <- max(SampGP, 1)
  return(list(NoGPs=SampGP, CostGPs=GPcost))
}





# 5. PATIENT_COVARIATES

# Extracts covariates for a single patient given the patient ID number

# @param ID Patient identifier - 'e_patid'
# @param Data Data frame containing patient-level records and required covariate columns

# @return One-row data frame containing specified covariates

patient_covariates <- function(ID, Data) {
  
  # Filter to requested patient and return selected covariates
  Data <- Data %>% filter(e_patid == ID) %>%
    filter(row_number() == 1)
  
  return(Data %>% dplyr::select(AgeDiag,gender,upci,imd,cms_burden,smoking_hist,AorD))
}




# 6. COST_REF_APPTS

# Simulates number of subsequent referral appointments given referral priority and 
# cancer type 
# Associated cost is given for first and all subsequent appointments
# Costs are sampled based on frequency of treatment specialists recorded for a given cancer type

# @param Fit A fitted negative binomial distribution object 
# @param Priority Either 2WW or non-2WW
# @param Cancer Cancer type
# @param Cost_data Data frame containing treatment specialists for cancer type and priority.
#                  First and follow-up appointment costs included


# @return A named list with simulated number of subsequent appointments and associated cost

cost_ref_appts <- function(Fit, Priority, Cancer, Cost_data, Model_meth=Model_method, u_spec=NULL, u_apps=NULL) {
  
  # Restrict cost data to the requested priority and cancer type
  Cost_data <- Cost_data %>% filter(priority == Priority, Type == Cancer)
  
  # Sample a specialty using observed subgroup proportions
  Spec_probs <- prop.table(table(Cost_data$tretspef))
  Spec_levels <- names(Spec_probs)
  
  if (is.null(u_spec)) {
    Samp_spec <- sample(Spec_levels, size = 1, prob = Spec_probs)
  } else {
    u_spec <- min(max(u_spec, 1e-12), 1 - 1e-12)
    Samp_spec <- Spec_levels[which(u_spec <= cumsum(Spec_probs))[1]]
  }
  
  # Retrieve the consultant unit cost for the sampled specialty
  if (Model_meth == 'Deterministic') {
    Cost_first <- Cost_data %>% filter(tretspef == Samp_spec) %>% dplyr::select(FirstDet) %>% distinct()
    Cost_followup <- Cost_data %>% filter(tretspef == Samp_spec) %>% dplyr::select(FollowDet) %>% distinct()
  } else {
    Cost_first <- Cost_data %>% filter(tretspef == Samp_spec) %>% dplyr::select(FirstPSA) %>% distinct()
    Cost_followup <- Cost_data %>% filter(tretspef == Samp_spec) %>% dplyr::select(FollowPSA) %>% distinct()
  }
  
  # Simulate subsequent appointment count
  # Add 1 for first appointment not included in fitted distribution
  if (is.null(u_apps)) {
    NoApps <- rnegbin(1, mu = Fit$estimate["mu"], theta = Fit$estimate["size"])
  } else {
    u_apps <- min(max(u_apps, 1e-12), 1 - 1e-12)
    NoApps <- qnbinom(u_apps, mu = Fit$estimate["mu"], size = Fit$estimate["size"])
  }
  
  return(list(NoApps=NoApps, CostFirst=Cost_first, CostFollow=Cost_followup))
}





# 

Cancer_type_pred <- function(Fit=Model_cancer, Covariates, Event, Order=Cancer_order) {
  
  newdata = data.frame(Event=Event, gender=Covariates$gender, 
                       AgeDiag=Covariates$AgeDiag, imd=Covariates$imd,
                       smoking_hist=Covariates$smoking_hist,
                       upci=Covariates$upci, cms_burden=Covariates$cms_burden)
  
  newdata$Event <- as.factor(newdata$Event)
  newdata$gender <- as.factor(newdata$gender)
  newdata$smoking_hist <- as.factor(newdata$smoking_hist)
  
  Pred = predict(Fit, newdata = newdata,
          type = "class")
  Cancer <- as.numeric(Cancer_order[Pred])
  return(Cancer)
}


Cancer_stage_pred <- function(Fit=Model_stage, Covariates, Event, Time, Cancer) {
  
  newdata = data.frame(Type=Cancer, Event=Event, Time=Time, gender=Covariates$gender, 
                       AgeDiag=Covariates$AgeDiag, imd=Covariates$imd,
                       smoking_hist=Covariates$smoking_hist,
                       cms_burden=Covariates$cms_burden)
  newdata$Type <- as.factor(newdata$Type)
  newdata$Event <- as.factor(newdata$Event)
  newdata$gender <- as.factor(newdata$gender)
  newdata$smoking_hist <- as.factor(newdata$smoking_hist)
  
  Pred = predict(Fit, newdata = newdata,
                 type = "probs")
  return(Pred)
}


Cancer_stage_REG <- function(Fit=Model_stage_REG, Covariates, Pathway, Cancer) {

  newdata = data.frame(Type=Cancer, 
                       gender=Covariates$gender, 
                       AgeDiag=Covariates$AgeDiag, 
                       upci=Covariates$upci,
                       cms_burden=Covariates$cms_burden,
                       final_route=Pathway)
  newdata$Type <- as.factor(newdata$Type)
  newdata$gender <- as.factor(newdata$gender)
  newdata$final_route <- as.factor(newdata$final_route)
  
  Pred = predict(Fit, newdata = newdata,
                 type = "probs")
  return(Pred)
}



# 7. MULTI-STATE SUB TRAJECTORIES - MAKE_STAGE_TRAJ

# Creates a trajectory for the cancer staging after patient has been diagnosed from
# any pathway

# @param traj_name Name of the trajectory
# @param fit_model A fitted 'msm' model object passed to simulate_next
# @param pathway Current state passed to simulate_next
# @param prefix Character used to construct attribute and branch names

# @return A simmer trajectory object


make_stage_traj <- function(traj_name, fit_model, pathway, prefix, Covariates_fun, Diag_cost1, Diag_cost2, Treat_cost) {
  
  trajectory(name = traj_name) %>%
    
    # Predict cancer type 
    set_attribute(keys = 'Cancer_Type_Diag', values = function() {
      if (pathway %in% c("2WWDiag", "N2WWDiag")) {
        Type = get_attribute(.env = simmodel, keys = 'Cancer_Type')
      } else {
        Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"), Data=Data_Full_PatientsN)
        Type = Cancer_type_pred(Covariates=Cov_pat, Event=pathway)
      } 
      
      if (Type == 5) {
        Type <- sample(1:4, 1, prob=Prob_type)
      }
      as.numeric(Type)
    }) %>%
    
    # Discount cost of current pathway
    set_attribute(keys = "total_cost", values = function() {
      DCost1 <- Disc_cost(cost = Diag_cost1, time_total = simmer::now(.env = simmodel))
      
      if (all(is.na(Diag_cost2))) {
        DCost2 <- 0
      } else {
        Type = get_attribute(.env = simmodel, keys = 'Cancer_Type')
        D_cost2 <- Diag_cost2 %>% filter(Cancer == Type, Stage == 'Diagnosis')
      
        if (Type == 'lung') {
          D_cost2 <- sample(D_cost2$Inf_Costs, 1, prob=as.numeric(D_cost2$Distribution)/100)
        } else {
          D_cost2 <- sum(apply(D_cost2, 1, function(x) sample(c(as.numeric(x[4]),0), 1, prob=c(as.numeric(x[3])/100,1-as.numeric(x[3])/100))))
        }
        DCost2 <- Disc_cost(cost = D_cost2, time_total = simmer::now(.env = simmodel))
      }
      DCost <- DCost1 + DCost2
      as.numeric(DCost)}, mod = "+") %>%
    
    # Set attribute to determine which stage
    set_attribute(keys = paste0("Diag_", prefix), value = function() {
      # Simulate transition probabilities from msm model
      Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"),
                                    Data=Data_Full_PatientsN)

      Probs <- Cancer_stage_pred(Covariates=Cov_pat, 
                                  Event=pathway, 
                                  Time=simmer::now(.env = simmodel), 
                                  Cancer=names(Cancer_order[as.numeric(get_attribute(.env = simmodel, keys = 'Cancer_Type_Diag'))]))

      
      state_order <- c(1,2,3,4)
      Probs <- as.numeric(Probs[state_order])
      
      if (Model_random == "Fixed_random") {
        u_branch <-  Rand_set %>% filter(i == get_attribute(.env = simmodel, keys = "Patient_id"), 
                                         j == 1) %>% 
          dplyr::pull(Stage_branch)
        Stage <- which(u_branch <= cumsum(Probs))[1]
      } else{
        Stage <- sample(seq_along(state_order), 1, prob = Probs)
      }
      
      as.numeric(Stage)
    }) %>%

    set_attribute(keys = "total_cost", values = function() {
      Type <- names(Cancer_order[get_attribute(.env = simmodel, keys = 'Cancer_Type_Diag')])
      StageDiag <- get_attribute(.env = simmodel, keys = paste0("Diag_", prefix))
      DCost <- Disc_cost(cost = Treat_cost %>% filter(Cancer==Type, Stage==StageDiag) %>% pull(Inf_Costs),
                time_total = simmer::now(.env = simmodel))
      
      as.numeric(DCost)}, mod = "+") %>%
    
    set_attribute(keys = "total_utility", values = Data_Utilities %>% filter(Item == 'AnxDiagnosis') %>% pull(Utility), mod = "-") %>%         
    
    # Utilities - ALL, discounted QALYs - BREAST STAGE 1,2,3 ONLY
    set_attribute(keys = c("total_utility","total_QALYs","BreastTimeLag"), values = function() {
      
      Current_utility <- get_attribute(.env = simmodel, keys = "total_utility")
      Current_QALYs <- get_attribute(.env = simmodel, keys = "total_QALYs")
      Type <- names(Cancer_order[get_attribute(.env = simmodel, keys = 'Cancer_Type_Diag')])
      StageDiag <- get_attribute(.env = simmodel, keys = paste0("Diag_", prefix))
      UtilityAll <- Data_Utilities %>% filter(Item == Type)
      
      if (Type == 'breast' & StageDiag < 4) {
        
        Q_length1 <- c(0.25, 0.375, 0.75)[StageDiag]
        Q_length2 <- c(4.75, 4.625, 4.25)[StageDiag]
        
        Utility1 <- UtilityAll %>% filter(Stage == StageDiag) %>% slice_head(n=1) %>% pull(Utility)
        Tutility1 <- Current_utility + Utility1
        Q_add1 <- Disc_QALYs(utility=Tutility1, time_total=(simmer::now(.env = simmodel)+Q_length1), current_time=Q_length1)
        
        Utility2 <- UtilityAll %>% filter(Stage == StageDiag) %>% slice_tail(n=1) %>% pull(Utility)
        Tutility2 <- Current_utility + Utility2
        Q_add2 <- Disc_QALYs(utility=Tutility2, time_total=(simmer::now(.env = simmodel)+Q_length2), current_time=Q_length2)
        
        Q_add <- Q_add1 + Q_add2
        Q_length <- (Q_length1 + Q_length2)*365
        Utility <- Current_utility
        
      } else {
        Q_add <- 0
        Q_length <- 0
        
        Utility <- UtilityAll %>% filter(Stage == StageDiag)
        
        if (Utility$MorA == 'M') {
          Utility <- Current_utility * (Utility %>% pull(Utility))
        } else {
          Utility <- Current_utility + (Utility %>% pull(Utility))
        }
      }
      
      QALYs <- Current_QALYs + Q_add
      
      c(Utility, QALYs, Q_length)
      
    }) %>%
    
    # Time delay in current state - BREAST ONLY
    timeout_from_attribute(key = "BreastTimeLag") %>%
    
    # Discounted QALYs - LIFETIME
    set_attribute(keys = "total_QALYs", values = function() {
      Tutility <- get_attribute(.env = simmodel, keys = "total_utility") 
      Ttime <- simmer::now(.env = simmodel)
      Type <- names(Cancer_order[get_attribute(.env = simmodel, keys = 'Cancer_Type_Diag')])
      StageDiag <- get_attribute(.env = simmodel, keys = paste0("Diag_", prefix))
      Cov_pat <- patient_covariates(ID = get_attribute(.env = simmodel, keys = "Patient_id"),
                                    Data=Data_Full_PatientsN)
      newdata <- data.frame(Type=Type, stage_cat=StageDiag, AgeDiag=Cov_pat$AgeDiag)
      calc_extrapolated_qaly(newdata=newdata, utility=Tutility, endTime=Ttime)
    }, mod = "+")
           
}




# 8. DISC_COST - Add in function for discounting (even though most pathways in the DES are within one year)

# Converts time from days to years and applies discounting to a cost value

# @param cost Cost to be discounted
# @param time_days Time in days at which the cost occurs
# @param discount_rate Annual discount rate = 0.035

# @return Discounted cost

Disc_cost <- function(cost, time_total, discount_rate=Disc_rate) {
  # If time hasn't gone over 1 year, do not discount
  # Convert time from days to years
  if (time_total < 365) {
    time_years <- 0
  } else {
    #time_years <- floor((time_total)/365) # If you prefer to discount in year steps
    time_years <- ((time_total)/365) # If you prefer constant discounting
  }
  
  # Return discounted cost
  Dis_Cost <- cost/((1 + discount_rate)^time_years)
  return(Dis_Cost)
}




# Discount QALYs

Disc_QALYs <- function(utility, time_total, current_time, discount_rate=Disc_rate, cycle_length=1/52) {
  # If time hasn't gone over 1 year, do not discount
  # Convert time from days to years
  
  # If you prefer to discount in year steps
  #Dis_QALY <- 0
  #Q_time <- current_time
  #for (i in 1:ceiling(current_time/365)) {
  #  time_years <- floor((time_total-current_time)/365) + (i-1)
  #  Dis_QALY <- Dis_QALY + (utility/((1 + discount_rate)^time_years) * (min(Q_time, 365)/365))
  #  Q_time <- Q_time - 365
  #}
  
  # If you prefer constant discounting
  times <- unique(c(seq((time_total-current_time)/365, time_total/365, by=cycle_length), time_total/365))
  diff_times <- diff(times)
  times[times < 1] <- 0
  discount <- 1 / ((1 + discount_rate) ^ (times))
  
  Dis_QALY <- sum(utility * discount[-1] * diff_times)
  
  return(Dis_QALY)
}





# 20% se function
Fixed_PSA_Det <- function(Mean, se_prop = default_se_prop) {
  Det <- Mean
  SD <- Mean * se_prop
  PSA <- rnorm(1, mean = Mean, sd = SD)
  return(list(Det=Det, PSA=PSA))
}


# Expert elicitation function


# sampling function
sample_fun <- function(X, CDF, N, U=NULL) {
  if (is.null(U)) {
    U <- matrix(runif(N), N, 1)
  } else {
    U <- U
  }
  Samps <- apply(U, 1, function(u) {X[which(CDF >= u)][1]})
  return(Samps)
}

Expert_elic_DES <- function(Density, PSA_no=100, PSA_i, U_samp=NULL) {
  
  x <- Density[, 1]
  f <- Density[, dim(Density)[2]]
  
  dens <- f / sum(f)
  
  #if (is_deterministic) {
  Output_DET <- weighted.mean(x, dens)
  #} else {
  
  # cumulative distribution
  cdf <- cumsum(dens)
  
  Samples <- sample_fun(X=x, CDF=cdf, N=PSA_no, U=U_samp)
  Output_PSA <- Samples[PSA_i]
  #}
  
  return(list(DET=Output_DET, PSA=Output_PSA))
}



# Survival QALYs

calc_extrapolated_qaly <- function(fit=fit_surv2, newdata, utility, endTime,
                                   horizon_years = 60,
                                   cycle_length = 1/12,
                                   discount_rate = 0.035) {
  
  endTime <- endTime/365
  
  times <- seq(endTime, horizon_years+endTime, by = cycle_length)
  times[times < 1] <- 0
  
  surv <- summary(
    fit,
    newdata = newdata,
    type = "survival",
    t = times,
    tidy = TRUE
  )
  
  s <- surv$est
  
  # area under survival curve using cycle approximation
  discount <- 1 / ((1 + discount_rate) ^ (times))
  
  qaly <- sum(utility * s * discount * cycle_length)
  
  return(qaly)
}


# Survival curve

build_survival_curve <- function(qx, start_age, horizon = 60) {
  
  surv <- numeric(horizon)
  surv[1] <- 1
  
  for (t in 2:horizon) {
    age <- start_age + t - 1
    surv[t] <- surv[t-1] * (1 - qx[t-1])
  }
  
  surv
}



# life_table has columns: age, qx
calc_no_cancer_qaly <- function(age_start, endTime, utility, life_table,
                                discount_rate = 0.035,
                                max_age = 100) {
  
  endTime <- endTime/365
  ages <- age_start:max_age
  qx <- life_table[-(1:(age_start-40))]
  
  # probability alive at start of each future year
  surv <- cumprod(c(1, 1 - qx[-length(qx)]))
  
  years <- 0:(length(ages) - 1) + endTime
  years[years < 1] <- 0
  discount <- 1 / ((1 + discount_rate) ^ years)
  
  qaly <- sum(utility * surv * discount)
  return(qaly)
}









