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

# ---------- 1. Socioeconomic ----------
# load panel
AMT <- open_dataset(paste0(output.dir, "/dt_w_AMT.parquet")) %>%
  filter(birth_year >= birth_yr_lb & birth_year < birth_yr_ub & age >= 25 & age < 65) %>% 
  collect() %>%
  as.data.table()

head(AMT)
AMT[income > 0, log_income := log(income)]
AMT[, gender_combo := "male-male"]
AMT[male  == 1 & patient_male == 0, gender_combo := "male-female"]
AMT[male  == 0 & patient_male == 1, gender_combo := "female-male"]
AMT[male  == 0 & patient_male == 0, gender_combo := "female-female"]
AMT[, gender_combo := factor(gender_combo, levels = c(
  "male-male",
  "male-female",
  "female-male",
  "female-female"
))]

# set models
Ys <- c(
  "income",
  "log_income",
  "work"
)
models <- c(
  "~ CP | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
  "~ CP + male + patient_male | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
  "~ CP + male + patient_male | birth_year + patient_birth_year + year + father_edu + mother_edu + f_birth_year + m_birth_year"
)

# ---------- Income ----------
labor_reg_list <- list()
control_mean_labor <- data.table(variable = c(), control_mean = c())
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
    AMT_filtered <- AMT[complete.cases(AMT[, ..vars])]  # 改
    control_mean_labor <- rbind(control_mean_labor,
                                data.table(paste0(Y, i),
                                           mean(AMT_filtered[CP == 0, get(Y)], na.rm = TRUE),
                                           sd(AMT_filtered[CP == 0, get(Y)], na.rm = TRUE)))
    
    reg <- feols(as.formula(f[[1]]),
                 data = AMT_filtered,
                 cluster = c("ID"))
    labor_reg_list <- c(labor_reg_list, list(reg))
    rm(AMT_filtered); gc()
  }
}
esttable(labor_reg_list)
rm(AMT); gc()

# ---------- Education ----------

edu <- open_dataset(paste0(output.dir, "/dt_w_educ.parquet")) %>%
  filter(birth_year >= birth_yr_lb & birth_year < birth_yr_ub & age >= 25 & age < 65) %>% 
  collect() %>%
  as.data.table()

edu[, `:=`(high_school = ifelse(highest_edu >= 12, 1, 0),
           bachelor = ifelse(highest_edu >= 16, 1, 0))]
edu[, gender_combo := "male-male"]
edu[male  == 1 & patient_male == 0, gender_combo := "male-female"]
edu[male  == 0 & patient_male == 1, gender_combo := "female-male"]
edu[male  == 0 & patient_male == 0, gender_combo := "female-female"]
edu[, gender_combo := factor(gender_combo, levels = c(
  "male-male",
  "male-female",
  "female-male",
  "female-female"
))]
head(edu)

Ys <- c("highest_edu",
        "bachelor")

models <- c(
  "~ CP | birth_year + patient_birth_year + f_birth_year + m_birth_year",
  "~ CP + male + patient_male | birth_year + patient_birth_year + f_birth_year + m_birth_year",
  "~ CP + male + patient_male | birth_year + patient_birth_year + father_edu + mother_edu + f_birth_year + m_birth_year"
)
edu_reg_list <- list()
control_mean_edu <- data.table(variable = c(), control_mean = c())
for (Y in Ys) {
  message(Y)
  formulas <- lapply(models, function(x) {paste(Y, x)})
  
  for (i in 1:length(formulas)) {
    f <- formulas[i]
    print(f)
    # find all variables we need in the formula
    vars <- gsub("[|+~*]", "", f)
    vars <- strsplit(vars, "[, ]+")[[1]]
    print(vars)
    # subset rows with non-NA values
    # 原 edu_filtered <- edu[complete.cases(dt[, ..vars])]
    edu_filtered <- edu[complete.cases(edu[, ..vars])]  # 改
    control_mean_edu <- rbind(control_mean_edu,
                              data.table(paste0(Y, i),
                                         mean(edu_filtered[CP == 0, get(Y)], na.rm = TRUE),
                                         sd(edu_filtered[CP == 0, get(Y)], na.rm = TRUE)))
    
    reg <- feols(as.formula(f[[1]]),
                 data = edu_filtered,
                 vcov = "hetero")
    edu_reg_list <- c(edu_reg_list, list(reg))
    rm(edu_filtered); gc()
  }
}
esttable(edu_reg_list)
rm(edu); gc()

### Save socioeconomic outcomes
fwrite(esttable(labor_reg_list,
                edu_reg_list), file.path(result.dir, "main_result/socioecon_eq1_post_1980.csv"))

rm(labor_reg_list,
   edu_reg_list);gc()

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
  "~ CP + male + patient_male | birth_year + patient_birth_year + year + father_edu + mother_edu + f_birth_year + m_birth_year",
  "~ CP * male + patient_male | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
  "~ CP * male + patient_male | birth_year + patient_birth_year + year + father_edu + mother_edu + f_birth_year + m_birth_year"
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
    vars <- gsub("[|+~*]", "", f)
    vars <- strsplit(vars, "[, ]+")[[1]]
    print(vars)
    # subset rows with non-NA values
    # 原:AMT_filtered <- AMT[complete.cases(dt[, ..vars])]

    temp <- complete.cases(marr[, ..vars])
    marr_filtered <- marr[temp]  # 改
    rm(temp); gc()
    
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
esttable(marr_reg_list)
rm(marr); gc()

# ---------- Fertility ----------
# load panel
fert <- open_dataset(paste0(output.dir, "/dt_w_fert.parquet")) %>%
  filter(birth_year >= birth_yr_lb & birth_year < birth_yr_ub & age >= 25 & age < 65) %>% 
  collect() %>%
  as.data.table()

head(fert)
fert[, gender_combo := "male-male"]
fert[male  == 1 & patient_male == 0, gender_combo := "male-female"]
fert[male  == 0 & patient_male == 1, gender_combo := "female-male"]
fert[male  == 0 & patient_male == 0, gender_combo := "female-female"]
fert[, gender_combo := factor(gender_combo, levels = c(
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
  "~ CP + male + patient_male | birth_year + patient_birth_year + year + father_edu + mother_edu + f_birth_year + m_birth_year",
  "~ CP * male + patient_male | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
  "~ CP * male + patient_male | birth_year + patient_birth_year + year + father_edu + mother_edu + f_birth_year + m_birth_year"
)


fert_reg_list <- list()
control_mean_fert <- data.table(variable = c(), control_mean = c())
for (Y in Ys) {
  message(Y)
  formulas <- lapply(models, function(x) {paste(Y, x)})
  
  for (i in 1:length(formulas)) {
    f <- formulas[i]
    print(f)
    # find all variables we need in the formula
    vars <- gsub("[|+~*]", "", f)
    vars <- strsplit(vars, "[, ]+")[[1]]
    print(vars)
    # subset rows with non-NA values
    # 原:AMT_filtered <- AMT[complete.cases(dt[, ..vars])]
    
    temp <- complete.cases(fert[, ..vars])
    fert_filtered <- fert[temp]  # 改
    rm(temp); gc()
    
    control_mean_fert <- rbind(control_mean_fert,
                               data.table(paste0(Y, i),
                                          mean(fert_filtered[CP == 0, get(Y)], na.rm = TRUE),
                                          sd(fert_filtered[CP == 0, get(Y)], na.rm = TRUE)))
    
    reg <- feols(as.formula(f[[1]]),
                 data = fert_filtered,
                 cluster = c("ID"))
    fert_reg_list <- c(fert_reg_list, list(reg))
    rm(fert_filtered); gc()
  }
}
esttable(fert_reg_list)
rm(fert); gc()

### Save family function outcomes
#fwrite(esttable(marr_reg_list,
#                fertility_reg_list), file.path(result.dir, "/main_result/family_func_eq1_post_1980.csv"))
fwrite(esttable(marr_reg_list,
                fert_reg_list), file.path(result.dir, "/main_result/family_func_eq1_post_1980_int.csv"))


control_mean <- rbind(
  control_mean_marr,
  control_mean_fert
)

# save the control mean
fwrite(control_mean, file.path(result.dir, "/main_result/control_mean_family_func_post_1980.csv"))
rm(marr_reg_list, fert_reg_list); gc()

# ---------- Health ----------

opdte <- open_dataset(paste0(output.dir, "/dt_w_opdte.parquet")) %>% 
  filter(birth_year >= birth_yr_lb & birth_year < birth_yr_ub & age >= 25 & age < 65) %>% 
  collect() %>% 
  as.data.table()
head(opdte)
opdte[, WEIGHT_2 := WEIGHT * WEIGHT]
opdte[, WEEK_2 := WEEK * WEEK]
opdte[, age_diff := patient_birth_year - birth_year]
opdte[, eme_dummy := ifelse(eme == 0, 0, 1)]
opdte[, gender_combo := "male-male"]
opdte[male  == 1 & patient_male == 0, gender_combo := "male-female"]
opdte[male  == 0 & patient_male == 1, gender_combo := "female-male"]
opdte[male  == 0 & patient_male == 0, gender_combo := "female-female"]
opdte[, gender_combo := factor(gender_combo, levels = c(
  "male-male",
  "male-female",
  "female-male",
  "female-female"
))]

control_mean_health <- rbind(
  opdte[ CP == 0, .(outcome = "major", mean = mean(major, na.rm = TRUE), sd = sd(major), count = .N)],
  opdte[CP == 0, .(outcome = "chronic", mean = mean(chronic, na.rm = TRUE), sd = sd(chronic), count = .N)],
  opdte[CP == 0, .(outcome = "total_dot", mean = mean(total_dot, na.rm = TRUE), sd = sd(total_dot), count = .N)],
  opdte[CP == 0, .(outcome = "outpatient", mean = mean(outpatient_count, na.rm = TRUE), sd = sd(outpatient_count), count = .N)],
  opdte[CP == 0, .(outcome = "emergency", mean = mean(eme, na.rm = TRUE), sd = sd(eme, na.rm = TRUE), count = .N)],
  opdte[CP == 0, .(outcome = "eme_dummy", mean = mean(eme_dummy, na.rm = TRUE), sd = sd(eme, na.rm = TRUE), count = .N)],
  opdte[CP == 0, .(outcome = "no diagnosis", mean = mean(useless, na.rm = TRUE), sd = sd(useless, na.rm = TRUE), count = .N)]
)

#control_mean_edu <- control_mean_edu[, .(outcome, mean, sd)]
#control_mean_labor <- control_mean_labor[, .(outcome, mean, sd)]
#control_mean_health <- control_mean_health[, .(outcome, mean, sd)]
control_mean <- rbind(
  control_mean_labor,
  control_mean_edu,
  control_mean_health)

# save the control mean
fwrite(control_mean, file.path(result.dir, "main_result/control_mean_post_1980.csv"))

Ys <- c("outpatient_count",
        "total_dot",
        "eme_dummy",
        "major")

models <- c(
  "~ CP | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
  "~ CP + male + patient_male | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
  "~ CP + male + patient_male | birth_year + patient_birth_year + father_edu + mother_edu + year + f_birth_year + m_birth_year"
)


health_reg_list <- list()
control_mean_health <- data.table(variable = c(), control_mean = c())
for (Y in Ys) {
  message(Y)
  print(Y)
  formulas <- lapply(models, function(x) {paste(Y, x)})
  
  for (i in 1:length(formulas)) {
    f <- formulas[[i]]
    print(f)
    # find all variables we need in the formula
    vars <- gsub("[|+~]", "", f)
    vars <- strsplit(vars, "[, ]+")[[1]]
    print(vars)
    # subset rows with non-NA values
    # 原:opdte_filtered <- opdte[complete.cases(dt[, ..vars])]
    opdte_filtered <- opdte[complete.cases(opdte[, ..vars])]  # 改
    control_mean_health <- rbind(control_mean_health,
                                data.table(paste0(Y, i),
                                           mean(opdte_filtered[CP == 0, get(Y)], na.rm = TRUE),
                                           sd(opdte_filtered[CP == 0, get(Y)], na.rm = TRUE)))
    
    reg <- feols(as.formula(f),
                 data = opdte_filtered,
                 cluster = c("ID"))
    health_reg_list <- c(health_reg_list, list(reg))
    rm(opdte_filtered, health_filtered); gc()
  }
}
esttable(health_reg_list)
rm(opdte); gc()


#control_mean_health <- control_mean_health[, .(outcome, mean, sd)]

control_mean <- rbind(
  control_mean_labor,
  control_mean_edu,
  control_mean_health
)

# save the control mean
fwrite(control_mean, file.path(result.dir, "main_result/control_mean_post_1980.csv"))
rm(control_mean,
   control_mean_edu,
   control_mean_health,
   control_mean_labor
); gc()


fwrite(esttable(health_reg_list),
       file.path(result.dir, "main_result/health_eq1_post_1980.csv"))

lapply(models, function(x) {paste(Ys[1], x)})
for (f in formulas) {
  print(f)
}


# mental diseases -----------------------------------------------------------

opdte <- open_dataset(paste0(output.dir, "/dt_w_opdte_md.parquet")) %>% 
  filter(birth_year >= birth_yr_lb & birth_year < birth_yr_ub & age >= 25 & age < 65) %>% 
  collect() %>% 
  as.data.table()
head(opdte)
opdte[, WEIGHT_2 := WEIGHT * WEIGHT]
opdte[, WEEK_2 := WEEK * WEEK]
opdte[, age_diff := patient_birth_year - birth_year]
opdte[, eme_dummy := ifelse(eme == 0, 0, 1)]
opdte[, gender_combo := "male-male"]
opdte[male  == 1 & patient_male == 0, gender_combo := "male-female"]
opdte[male  == 0 & patient_male == 1, gender_combo := "female-male"]
opdte[male  == 0 & patient_male == 0, gender_combo := "female-female"]
opdte[, gender_combo := factor(gender_combo, levels = c(
  "male-male",
  "male-female",
  "female-male",
  "female-female"
))]

control_mean_health <- rbind(
  opdte[ CP == 0, .(outcome = "md", mean = mean(md, na.rm = TRUE), sd = sd(md), count = .N)],
  opdte[CP == 0, .(outcome = "md_visit", mean = mean(md_visit, na.rm = TRUE), sd = sd(md_visit), count = .N)],
  opdte[CP == 0, .(outcome = "psy", mean = mean(psy, na.rm = TRUE), sd = sd(psy), count = .N)],
  opdte[CP == 0, .(outcome = "psy_visit", mean = mean(psy_visit, na.rm = TRUE), sd = sd(psy_visit), count = .N)],
  opdte[CP == 0, .(outcome = "anx", mean = mean(anx, na.rm = TRUE), sd = sd(anx, na.rm = TRUE), count = .N)],
  opdte[CP == 0, .(outcome = "anx_visit", mean = mean(anx_visit, na.rm = TRUE), sd = sd(anx_visit, na.rm = TRUE), count = .N)],
  opdte[CP == 0, .(outcome = "dep", mean = mean(dep, na.rm = TRUE), sd = sd(dep, na.rm = TRUE), count = .N)],
  opdte[CP == 0, .(outcome = "dep_visit", mean = mean(dep_visit, na.rm = TRUE), sd = sd(dep_visit, na.rm = TRUE), count = .N)]
)

#control_mean_edu <- control_mean_edu[, .(outcome, mean, sd)]
#control_mean_labor <- control_mean_labor[, .(outcome, mean, sd)]
#control_mean_health <- control_mean_health[, .(outcome, mean, sd)]
#control_mean <-  control_mean_health

# save the control mean
#fwrite(control_mean, file.path(result.dir, "main_result/control_mean_md_post_1980.csv"))

Ys <- c("md",
        "md_visit",
        "psy",
        "psy_visit",
        "anx",
        "anx_visit",
        "dep",
        "dep_visit")

models <- c(
  "~ CP | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
  "~ CP + male + patient_male | birth_year + patient_birth_year + year + f_birth_year + m_birth_year",
  "~ CP + male + patient_male | birth_year + patient_birth_year + father_edu + mother_edu + year + f_birth_year + m_birth_year"
)


health_reg_list <- list()
control_mean_health <- data.table(variable = c(), control_mean = c())
for (Y in Ys) {
  message(Y)
  print(Y)
  formulas <- lapply(models, function(x) {paste(Y, x)})
  
  for (i in 1:length(formulas)) {
    f <- formulas[[i]]
    print(f)
    # find all variables we need in the formula
    vars <- gsub("[|+~]", "", f)
    vars <- strsplit(vars, "[, ]+")[[1]]
    print(vars)
    # subset rows with non-NA values
    # 原:opdte_filtered <- opdte[complete.cases(dt[, ..vars])]
    opdte_filtered <- opdte[complete.cases(opdte[, ..vars])]  # 改
    control_mean_health <- rbind(control_mean_health,
                                 data.table(paste0(Y, i),
                                            mean(opdte_filtered[CP == 0, get(Y)], na.rm = TRUE),
                                            sd(opdte_filtered[CP == 0, get(Y)], na.rm = TRUE)))
    
    reg <- feols(as.formula(f),
                 data = opdte_filtered,
                 cluster = c("ID"))
    health_reg_list <- c(health_reg_list, list(reg))
    rm(opdte_filtered, health_filtered); gc()
  }
}
esttable(health_reg_list)
rm(opdte); gc()


#control_mean_health <- control_mean_health[, .(outcome, mean, sd)]

control_mean <- control_mean_health

# save the control mean
fwrite(control_mean, file.path(result.dir, "main_result/control_mean_md_post_1980.csv"))

fwrite(esttable(health_reg_list),
       file.path(result.dir, "main_result/health_md_eq1_post_1980.csv"))
