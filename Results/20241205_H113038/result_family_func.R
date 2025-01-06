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

# ###############################################################################
# ### FUNCTIONS #################################################################
# ###############################################################################

###############################################################################
### MAIN
###############################################################################
# set upper and lower bound for birth_year
birth_yr_ub <- 2030
birth_yr_lb <- 1980

# ---------- Marriage ----------
# load panel
marr <- open_dataset(paste0(output.dir, "/dt_w_marr.parquet")) %>%
  filter(birth_year >= birth_yr_lb & birth_year < birth_yr_ub & age >= 25 & age < 65) %>% 
  collect() %>%
  as.data.table()

head(marr)
marr[, gender_combo := "male-male"]
marr[male  == 1 & patient_male == 0, gender_combo := "male-female"]
marr[male  == 0 & patient_male == 1, gender_combo := "female-male"]
marr[male  == 0 & patient_male == 0, gender_combo := "female-female"]
marr[, gender_combo := factor(gender_combo, levels = c(
  "male-male",
  "male-female",
  "female-male",
  "female-female"
))]

# set models
Ys <- c(
  "marr"
)
models <- c(
  "~ CP | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
  "~ CP + male + patient_male | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
  "~ CP + male + patient_male | birth_year + patient_birth_year + year + father_edu + mother_edu + f_birth_year + m_birth_year"
)


marr_reg_list <- list()
control_mean_marr <- data.table(variable = c(), control_mean = c())
for (Y in Ys) {
  message(Y)
  formulas <- lapply(models, function(x) {paste(Y, x)})
  
  for (i in 1:length(formulas)) {
    f <- formulas[i]
    print(f)
    # find all variables we need in the formula
    vars <- gsub("[|+~]", "", f)
    vars <- strsplit(vars, "[, ]+")[[1]]
    print(vars)
    # subset rows with non-NA values
    # 原:AMT_filtered <- AMT[complete.cases(dt[, ..vars])]
    temp <- complete.cases(marr[, ..vars])
    marr_filtered <- marr[temp]  # 改
    
    control_mean_marr <- rbind(control_mean_marr,
                                data.table(paste0(Y, i),
                                           mean(marr_filtered[CP == 0, get(Y)], na.rm = TRUE),
                                           sd(marr_filtered[CP == 0, get(Y)], na.rm = TRUE)))
    
    reg <- feols(as.formula(f[[1]]),
                 data = marr_filtered,
                 cluster = c("ID"))
    marr_reg_list <- c(marr_reg_list, list(reg))
    rm(marr_filtered); gc()
  }
}
#esttable(marr_reg_list)
rm(marr); gc()

# ---------- Fertility ----------
# load panel
fertility <- open_dataset(paste0(output.dir, "/dt_w_fertility.parquet")) %>%
  filter(birth_year >= birth_yr_lb & birth_year < birth_yr_ub & age >= 25 & age < 65) %>% 
  collect() %>%
  as.data.table()

head(fertility)
fertility[, gender_combo := "male-male"]
fertility[male  == 1 & patient_male == 0, gender_combo := "male-female"]
fertility[male  == 0 & patient_male == 1, gender_combo := "female-male"]
fertility[male  == 0 & patient_male == 0, gender_combo := "female-female"]
fertility[, gender_combo := factor(gender_combo, levels = c(
  "male-male",
  "male-female",
  "female-male",
  "female-female"
))]

# set models
Ys <- c(
  "fertility"
)
models <- c(
  "~ CP | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
  "~ CP + male + patient_male | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
  "~ CP + male + patient_male | birth_year + patient_birth_year + year + father_edu + mother_edu + f_birth_year + m_birth_year"
)


fertility_reg_list <- list()
control_mean_fertility <- data.table(variable = c(), control_mean = c())
for (Y in Ys) {
  message(Y)
  formulas <- lapply(models, function(x) {paste(Y, x)})
  
  for (i in 1:length(formulas)) {
    f <- formulas[i]
    print(f)
    # find all variables we need in the formula
    vars <- gsub("[|+~]", "", f)
    vars <- strsplit(vars, "[, ]+")[[1]]
    print(vars)
    # subset rows with non-NA values
    # 原:AMT_filtered <- AMT[complete.cases(dt[, ..vars])]
    temp <- complete.cases(fertility[, ..vars])
    fertility_filtered <- fertility[temp]  # 改
    
    control_mean_fertility <- rbind(control_mean_fertility,
                               data.table(paste0(Y, i),
                                          mean(fertility_filtered[CP == 0, get(Y)], na.rm = TRUE),
                                          sd(fertility_filtered[CP == 0, get(Y)], na.rm = TRUE)))
    
    reg <- feols(as.formula(f[[1]]),
                 data = fertility_filtered,
                 cluster = c("ID"))
    fertility_reg_list <- c(fertility_reg_list, list(reg))
    rm(fertility_filtered); gc()
  }
}
#esttable(fertility_reg_list)
rm(fertility); gc()

### Save socioeconomic outcomes
fwrite(esttable(marr_reg_list,
                fertility_reg_list), file.path(result.dir, "/main_result/family_func_eq1_post_1980.csv"))


control_mean <- rbind(
  control_mean_marr,
  control_mean_fertility
)

# save the control mean
fwrite(control_mean, file.path(result.dir, "/main_result/control_mean_family_func_post_1980.csv"))
