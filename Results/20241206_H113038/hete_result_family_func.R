###############################################################################
# Results for CP
###############################################################################
rm(list = ls()); gc()
options(scipen = 999)
library(data.table)
library(arrow)
library(lmtest)
library(sandwich)
library(dplyr)
library(ggplot2)
library(fixest)

###############################################################################
### PREAMBLE ##################################################################
###############################################################################
Disk <- "E"
setwd(paste0(Disk, ":/H113038/wllo/cp_siblings"))

source(paste0(Disk, ":/H113038/wllo/daughter/tool/arrow_helpers.R"))
source("scripts/functions_for_eq1_and_eq2.R")

Disk <- "E"
output.dir <- paste0(Disk, ":/H113038/wllo/cp_siblings/data/processed")
par.dir <- paste0(Disk, ":/H113038/parquet")
result.dir <- paste0(Disk, ":/H113038/wllo/cp_siblings/results")
CleanData.dir <- paste0(Disk, ":/H113038/CleanData/processed_data")
figure.dir <- paste0(result.dir, "/figure")
ref.dir <- paste0(Disk, ":/H113038/wllo/cp_siblings/reference")
YEARS <- 89L:110L

GET_PARENTS_DETAIL <- FALSE


### Simple Functions ##########################################################
print.time.diff <- function(now) { cat("\t",'time:',round(difftime(Sys.time(),now,units='mins'),2),"mins\n") }
start_timer <- function() { start_time <<- Sys.time() }
end_timer <- function() { end_time <<- Sys.time(); print(end_time - start_time) }
# 

###############################################################################
### MAIN
###############################################################################
# set upper and lower bound for birth_year
birth_yr_ub <- 2030
birth_yr_lb <- 1980

for (m in c(0, 1)) {
  
  # ---------- Marriage ----------
  marr <- open_dataset(paste0(output.dir, "/dt_w_marr.parquet")) %>%
    filter(birth_year >= birth_yr_lb & birth_year < birth_yr_ub & age >= 25 & age < 65 & male == m) %>% 
    collect() %>%
    as.data.table()
  
  
  control_mean_marr <-  marr[year == 110 & CP == 0, .(outcome = "marr",
                                                      mean = mean(marr, na.rm = TRUE),
                                                      sd = sd(marr, na.rm = TRUE),
                                                      count = .N)]
  
  # set models
  models <- c(
    "~ CP + patient_male | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
    "~ CP + patient_male | birth_year + patient_birth_year + year + father_edu + mother_edu + f_birth_year + m_birth_year"
  )
  
  marr_reg_list <- list()
  for (f in lapply(models, function(x) { paste("marr", x) })) {
    print(f)
    marr_reg <- feols(as.formula(f),
                        data = marr,
                        cluster = c("ID"))
    marr_reg_list <- c(marr_reg_list, list(marr_reg))
  }
  esttable(marr_reg_list)
  
  rm(marr); gc()
  
  # --------- Fertility ---------
  fertility <- open_dataset(paste0(output.dir, "/dt_w_fertility.parquet")) %>%
    filter(birth_year >= birth_yr_lb & birth_year < birth_yr_ub & age >= 25 & age < 65 & male == m) %>% 
    collect() %>%
    as.data.table()
  
  control_mean_fertility <-  fertility[year == 110 & CP == 0, .(outcome = "fertility",
                                                      mean = mean(fertility, na.rm = TRUE),
                                                      sd = sd(fertility, na.rm = TRUE),
                                                      count = .N)]
  
  # set models
  models <- c(
    "~ CP + patient_male | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
    "~ CP + patient_male | birth_year + patient_birth_year + year + father_edu + mother_edu + f_birth_year + m_birth_year"
  )
  
  fertility_reg_list <- list()
  for (f in lapply(models, function(x) { paste("fertility", x) })) {
    print(f)
    fertility_reg <- feols(as.formula(f),
                      data = fertility,
                      cluster = c("ID"))
    fertility_reg_list <- c(fertility_reg_list, list(fertility_reg))
  }
  esttable(fertility_reg_list)
  
  rm(fertility); gc()
  
  control_mean <- rbind(
    control_mean_marr,
    control_mean_fertility
  )
  # save the control mean
  fwrite(control_mean, paste0(result.dir, "/hete/control_mean_family_func_post_1980_male", m, ".csv"))
  ### Save socioeconomic outcomes
  fwrite(esttable(marr_reg_list,
                  fertility_reg_list), paste0(result.dir, "/hete/family_func_eq1_post_1980_male", m, ".csv"))

  rm(marr_reg_list,
     fertility_reg_list);gc()
  rm(control_mean,
     control_mean_marr,
     control_mean_fertility
  ); gc()
}

# WIPWIPWIPWIP ------------------------------------------------------------


for (high in c(0, 1)) {
  
  # ---------- Marriage ----------
  if (high == 0) {
    marr <- open_dataset(paste0(output.dir, "/dt_w_marr.parquet")) %>%
      filter(birth_year >= birth_yr_lb & birth_year < birth_yr_ub & age >= 25 & age < 65 & (f_bachelor == 0 & m_bachelor == 0)) %>% 
      collect() %>%
      as.data.table()
  } else {
    marr <- open_dataset(paste0(output.dir, "/dt_w_marr.parquet")) %>%
      filter(birth_year >= birth_yr_lb & birth_year < birth_yr_ub & age >= 25 & age < 65 & (f_bachelor == 1 | m_bachelor == 1)) %>% 
      collect() %>%
      as.data.table()
  }
  
  
  control_mean_marr <-  marr[year == 110 & CP == 0, .(outcome = "marr",
                                                      mean = mean(marr, na.rm = TRUE),
                                                      sd = sd(marr, na.rm = TRUE),
                                                      count = .N)]
  
  
  # set models
  models <- c(
    "~ CP + male + patient_male | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
    "~ CP + male + patient_male | birth_year + patient_birth_year + year + father_edu + mother_edu + f_birth_year + m_birth_year"
  )
  
  marr_reg_list <- list()
  for (f in lapply(models, function(x) { paste("marr", x) })) {
    print(f)
    marr_reg <- feols(as.formula(f),
                        data = marr,
                        cluster = c("ID"))
    marr_reg_list <- c(marr_reg_list, list(marr_reg))
  }
  esttable(marr_reg_list)
  
  
  rm(marr); gc()
  

  # ---------- Fertility ----------

  if (high == 0) {
    fertility <- open_dataset(paste0(output.dir, "/dt_w_fertility.parquet")) %>%
      filter(birth_year >= birth_yr_lb & birth_year < birth_yr_ub & age >= 25 & age < 65 & (f_bachelor == 0 & m_bachelor == 0)) %>% 
      collect() %>%
      as.data.table()
  } else {
    fertility <- open_dataset(paste0(output.dir, "/dt_w_fertility.parquet")) %>%
      filter(birth_year >= birth_yr_lb & birth_year < birth_yr_ub & age >= 25 & age < 65 & (f_bachelor == 1 | m_bachelor == 1)) %>% 
      collect() %>%
      as.data.table()
  }
  
  
  control_mean_fertility <-  fertility[year == 110 & CP == 0, .(outcome = "fertility",
                                                      mean = mean(fertility, na.rm = TRUE),
                                                      sd = sd(fertility, na.rm = TRUE),
                                                      count = .N)]
  
  
  # set models
  models <- c(
    "~ CP + male + patient_male | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
    "~ CP + male + patient_male | birth_year + patient_birth_year + year + father_edu + mother_edu + f_birth_year + m_birth_year"
  )
  
  fertility_reg_list <- list()
  for (f in lapply(models, function(x) { paste("fertility", x) })) {
    print(f)
    fertility_reg <- feols(as.formula(f),
                      data = fertility,
                      cluster = c("ID"))
    fertility_reg_list <- c(fertility_reg_list, list(fertility_reg))
  }
  esttable(fertility_reg_list)
  
  
  rm(fertility); gc()
  
  control_mean <- rbind(
    control_mean_marr,
    control_mean_fertility
  )
  # save the control mean
  fwrite(control_mean, paste0(result.dir, "/hete/control_mean_family_func_post_1980_p_edu", high, ".csv"))
  ### Save socioeconomic outcomes
  fwrite(esttable(marr_reg_list,
                  fertility_reg_list), paste0(result.dir, "/hete/family_func_eq1_post_1980_p_edu", high, ".csv"))
  
  rm(marr_reg_list,
     fertility_reg_list);gc()
  rm(control_mean,
     control_mean_marr,
     control_mean_fertility
  ); gc()
}