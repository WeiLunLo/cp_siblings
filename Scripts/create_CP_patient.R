###############################################################################
# Find siblings relationship, identify CP patients, and
# construct a matched sample

# Pipeline:
# create_CP_patient >
# create_panel >
# result & hete_result

# Updates: 
###############################################################################
rm(list = ls()); gc()


library(data.table)
library(arrow)
library(dplyr)
library(ggplot2)
library(stringr)

# ---------- PREAMBLE ----------
# setwd("E:/H114028/wllo/cp_siblings")

source("code/utils/arrow_helpers.R")
source("code/utils/event_year_wage_helper.R")
source("code/utils/functions_for_eq1_and_eq2.R")

output.dir <- "E:/H114028/wllo/cp_siblings/data/processed"
par.dir <- "E:/H114028/data/parquet"
csv.dir <- "E:/H114028/data/csv"
result.dir <- "E:/H114028/wllo/cp_siblings/results"

CleanData.dir <- "E:/H114028/CleanData/processed_data"
YEARS <- 89L:110L

# get the icd codes for CP
all <- FALSE
if (all) {
  CPRelatedICD_dt <- fread("../reference/icd.csv") %>% 
    filter(ccs == 6.3) %>% 
    select(icd9, icd10) %>% 
    mutate(icd9 = str_remove(icd9, "\\."),
           icd10 = str_remove(icd10, "\\."))
  
  CPRelatedICD <- c(CPRelatedICD_dt$icd9, CPRelatedICD_dt$icd10)
} else {
  CPRelatedICD <- c(343, 3430, 3431, 3432, 3433, 3434, 3438, 3439, "G800", "G801", "G802", "G804", "G808", "G809")
}




### Simple Functions ##########################################################
print.time.diff <- function(now) { cat("\t",'time:',round(difftime(Sys.time(),now,units='mins'),2),"mins\n") }
start_timer <- function() { start_time <<- Sys.time() }
end_timer <- function() { end_time <<- Sys.time(); print(end_time - start_time) }


###############################################################################
### FUNCTIONS #################################################################
###############################################################################
load.year.data <- function(year) {
  ## This function read the OPDTE files for a year and look for
  ## CP outpatient records
  ## Then, it add up the dots by ID.
  ## Each row is an individual
  message(paste("Start Loading:", year))
  start_timer()
  data_list <- list()
  for (m in 1:12) {
    open_nhi <- open_nhi("OPDTE", year, m) %>% 
      filter((ICD9CM_1 %in% CPRelatedICD) | (ICD9CM_2 %in% CPRelatedICD) | (ICD9CM_3 %in% CPRelatedICD)) %>% 
      select(ID, T_DOT, ICD9CM_1, ICD9CM_2, ICD9CM_3) %>%
      collect() %>% 
      as.data.table()
    data_list <- c(data_list, list(open_nhi))
  }
  dtSumDot <- rbindlist(data_list)
  dtSumDot <- dtSumDot[, .(totalDot = sum(T_DOT), count = .N, ICD9CM_1, ICD9CM_2, ICD9CM_3), by = ID]
  dtSumDot[ , Year := year + 1911]
  end_timer()
  return(dtSumDot)
}

get_twins <- function(year) {
  message(paste("Start to load:", year))
  birth_record <- fread(paste0(csv.dir, "/H_BHP_BIRTH", year, ".csv"), select = c("ID_M", "DEL_NO"))
  birth_record[, Yr := year + 1911]
  twins <- birth_record[DEL_NO > 1, .(ID_M, Yr)]
  return(twins)
}
#get_birth_record <- function() {
#  mohw <- fread(paste0(csv.dir, "/H_MOHW_MCHD109.csv"))
#  mohw[, Yr := as.numeric(Yr)]
#  mohw[, `:=`(ID_M = ifelse(ID_M == "", NA, ID_M),
#              ID_F = ifelse(ID_F == "", NA, ID_F),
#              ID_C = ifelse(ID_C == "", NA, ID_C))]
#  return(mohw)
#}

# I changed year range beause we don't have 92 and 104 death file
add_death_info <- function(dt) {
  # death
  id_list <- dt$ID
  dt_death_all <- data.table()
  for (year in c(89:91, 93:103, 105:110)) {
    f <- paste0(par.dir, "/H_OST_DEATH", year, ".parquet")
    print(paste("Start to Load:", f))
    death_info <- open_dataset(f) %>% 
      filter(ID %in% id_list) %>% 
      select(ID, D_DATE) %>%  
      collect() %>% 
      as.data.table()
    dt_death_all <- rbind(death_info, dt_death_all, fill=TRUE)
  }
  dt <- merge(dt, dt_death_all, all.x = TRUE, by = "ID")
  return(dt)
}


###############################################################################
### MAIN ######################################################################
###############################################################################
## Identify twins =============================================================
# Identify twins's mother using H_BHP_BIRTH

prsn_info <- open_dataset("E:/H114028/CleanData/processed_data/PersInfo_death.parquet") %>%
  mutate(BIRTH_Y = round(ID_BIRTHDAY / 10000)) %>% 
  filter(BIRTH_Y %in% 2004:2020) %>% # New
  collect() %>% 
  setDT()

enrol_relation <- open_dataset("E:/H114028/CleanData/processed_data/enrol_relation.parquet") %>%
  filter(ID %in% prsn_info$ID) %>% 
  collect() %>% 
  setDT()
setnames(enrol_relation, "M_ID", "ID_M")

enrol_relation[prsn_info, on = .(ID), Yr := i.BIRTH_Y]

twins_M <- lapply((90: 110), get_twins) %>% rbindlist()
# load mohw (婦幼)
# mohw <- get_birth_record()
# merge twins_M and mohw by mother's id and year of children's birth
twins <- merge(twins_M, enrol_relation, by = c("ID_M", "Yr")) %>%
  select(ID_M, ID) %>%
  distinct() %>% 
  setDT()
setnames(twins, "ID", "ID_C")
twin_id <- twins$ID_C


## Create CP Patient ==========================================================
# 1. get CP from OPDTE data
dt_CP <- lapply(YEARS, load.year.data) %>% rbindlist(fill = TRUE)
# remove twins and get the very first diagnosis record
dt_CP <- dt_CP[!(ID %in% twin_id)]

# 2. get CP from catastrophic data
catas <- fread(file.path(csv.dir, "H_NHI_CATAS111.csv"))
# keep CP but exclude twins
catas <- catas[DISE_CODE %in% CPRelatedICD & !(ID %in% twin_id)]
catas[, Year := APPL_DATE %/% 10000]
catas <- catas[, .(ID, Year)]

# merge CP from both sources
#all_CP <- rbind(dt_CP, catas, fill = TRUE)
#rm(catas); gc()
all_CP <- dt_CP

# keep only their first record
dtCPOnset <- all_CP[ , .(onsetYr = min(Year)), by = ID]

# save the file
write_parquet(dt_CP, paste0(output.dir, "/dt_CP.parquet"))
write_parquet(dtCPOnset, paste0(output.dir, "/CP_Onset.parquet"))

## Merge with personal info and children-parents relation =====================
dtCPOnset <- open_dataset(paste0(output.dir, "/CP_Onset.parquet")) %>%
  collect() %>% 
  as.data.table() # 68479 samples; in H114028: 145239

# load personal information
pers_info <- open_dataset(paste0(CleanData.dir, "/PersInfo.parquet")) %>%
  filter(valid == TRUE & !(ID %in% twin_id)) %>% 
  select(ID, ID_BIRTHDAY, ID_S) %>% 
  collect() %>% 
  as.data.table()
# make gender numeric
pers_info[, ID_S := as.numeric(ID_S)]
# create birth year
pers_info[, birth_year := ID_BIRTHDAY %/% 10000]

# merge onset info with birth year and calculate onset age
dtCPOnset <- merge(dtCPOnset, pers_info, by = "ID") # 67800; in H114028: 143946
dtCPOnset[, onset_age := onsetYr - birth_year]
pers_info <- pers_info[!(ID %in% dtCPOnset[onset_age < 0]$ID)]


dtCPOnset <- dtCPOnset[onset_age >= 0] # 67209; in H114028: 143943
write_parquet(dtCPOnset, paste0(output.dir, "/CP_Onset_w_pers_info.parquet"))
dtCPOnset <- open_dataset(paste0(output.dir, "/CP_Onset_w_pers_info.parquet")) %>%
  collect() %>% 
  as.data.table() # 67209 samples; in H114028: 143943

# load children-parent linkage data
# relation_dt <- open_dataset(paste0(CleanData.dir, "/hh_structure.parquet")) %>% 
#   select(ID, F_ID, M_ID, valid_F, valid_M) %>% 
#     collect() %>%
#     as.data.table() # 15647726
relation_dt <- open_dataset("E:/H114028/CleanData/processed_data/enrol_relation.parquet") |> 
  collect() |> 
  setDT()
# setnames(relation_dt, "ID_M", "M_ID")
#remove invalid parents
# relation_dt[, `:=`(
#   F_ID = ifelse(valid_F, F_ID, NA),
#   M_ID = ifelse(valid_M, M_ID, NA)
#   )]
relation_dt <- relation_dt[!(is.na(F_ID) & is.na(M_ID))] # 15454012; in H114028: 14232913
# relation_dt[, `:=`(valid_F = NULL, valid_M = NULL)]

both_p <- relation_dt[!is.na(F_ID) & !is.na(M_ID), .(F_ID, M_ID)] %>% distinct() %>% setDT()

both_p[, count_f := .N, by = F_ID]
both_p[, count_m := .N, by = M_ID]
both_p <- both_p[count_f == 1 & count_m == 1]
both_p[, count_f := NULL]
both_p[, count_m := NULL]
valid_father <- both_p$F_ID
valid_mother <- both_p$M_ID

relation_dt <- relation_dt[F_ID %in% valid_father | M_ID %in% valid_mother] # 11383113; in H114028: 7179927
father_mother <- relation_dt[!is.na(F_ID) & !is.na(M_ID)] # 9884546, H114028: 5503519
only_father <- relation_dt[is.na(M_ID), .(ID, F_ID)] # 739242; in H114028: 867772
only_father <- merge(only_father, both_p, by = "F_ID") # 739242; in H114028: 867772
only_mother <- relation_dt[is.na(F_ID), .(ID, M_ID)] # 759325; in H114028: 808636
only_mother <- merge(only_mother, both_p, by = "M_ID") # 759325; in H114028: 808636

relation_dt <- rbind(father_mother[, .(ID, F_ID, M_ID)],
                     only_father[, .(ID, F_ID, M_ID)],
                     only_mother[, .(ID, F_ID, M_ID)],
                     use.names = T)

rm(both_p, father_mother, only_father, only_mother, valid_father, valid_mother); gc()

# Merge pers_info with relation
pers_info <- merge(pers_info, relation_dt, by = "ID") 
head(pers_info)

# # Remove those diagnosed as a CP patient after 40 from our sample
# pers_info <- pers_info[!(ID %in% older_CP$ID)] # 15398202

write_parquet(pers_info, paste0(output.dir, "/pers_info_relation_for_CP.parquet"))
pers_info <- open_dataset(paste0(output.dir, "/pers_info_relation_for_CP.parquet")) %>% 
  collect() %>% 
  as.data.table()

pers_info[, both_parents := ifelse(!is.na(F_ID) & !is.na(M_ID), TRUE, FALSE)]
rm(mohw, twins, twins_M, relation_dt); gc()

###############################################################################
## construct siblings relation, but only requires the second child be CP
###############################################################################
# try
extensive_sib <- merge(pers_info, pers_info[, .(patient = ID, patient_BD = ID_BIRTHDAY, patient_S = ID_S, M_ID, F_ID)],
                       by = c("F_ID", "M_ID"), allow.cartesian = TRUE) # 27540546; in H114028: 17295302
# extensive_sib[, `:=`(F_ID = NULL, M_ID = NULL)]
extensive_sib <- extensive_sib[ID_BIRTHDAY < patient_BD] # 8005732; in H114028: 499583

# `count` denotes the number of younger siblings one has
extensive_sib[, count := .N, by = ID]
# extensive_sib[, num_children := count + 1]
# when patient_count = n, it means that the patient is the (n+1)th child
extensive_sib[, patient_count := .N, by = patient]

# identify CP patient (younger sibling)
extensive_sib[, CP := 0]
extensive_sib[patient %in% dtCPOnset$ID, CP := 1]
# extensive_sib[, num_CP_younger := sum(CP), by = ID]

# keep only the 1st child as ID and second as patient
extensive_sib <- extensive_sib[patient_count == 1] 
# the first-born children should not have CP
extensive_sib <- extensive_sib[!(ID %in% dtCPOnset$ID)] # 4090009; in H114028: 2583287

older_ID <- extensive_sib$ID
younger_ID <- extensive_sib$patient
# The oldest sibling should not be someone else's younger sibling
extensive_sib <- extensive_sib[!(ID %in% younger_ID)] 
# The second sibling should not be some first child's older sibling
extensive_sib <- extensive_sib[!(patient %in% older_ID)] # 4090003; in H114028: 2582931

# check again if appears more than once (has twin younger brothers/sisters)
extensive_sib[, count := .N, by = ID]
extensive_sib <- extensive_sib[count == 1] # 4034253; in H114028: 2545676
extensive_sib[, count := NULL]
extensive_sib[, `:=`(patient_count = NULL)]

extensive_sib[, patient_birth_year := patient_BD %/% 10000]

extensive_sib[, birth_month := ID_BIRTHDAY %% 10000 %/% 100]
head(extensive_sib)

# add birth weight and weeks of pregnancy
BIRTH_dt <- data.table()
for (yr in 90:110) {
  message(yr)
  # read birth record
  BIRTH <- fread(paste0(csv.dir, "/H_BHP_BIRTH", yr, ".csv"), select = c("ID_M", "BIRTH_YM", "SEX", "WEEK", "WEIGHT"))
  BIRTH[, birth_year := BIRTH_YM %/% 100]
  BIRTH[, birth_month := BIRTH_YM %% 100]
  
  # merge with dt using mother's ID and first-borns' birthday ans gender as keys
  BIRTH <- merge(extensive_sib, BIRTH, by.x = c("M_ID", "birth_year", "birth_month", "ID_S"), by.y = c("ID_M", "birth_year", "birth_month", "SEX"))
  
  # bind every year
  BIRTH_dt <- rbind(BIRTH_dt, BIRTH[, .(ID, WEEK, WEIGHT)])
  rm(BIRTH); gc()
}

extensive_sib <- merge(extensive_sib, BIRTH_dt, all.x = TRUE, by = "ID")
rm(BIRTH_dt)
# add parents' education
father_dt <- extensive_sib[!is.na(F_ID), .(F_ID)]
mother_dt <- extensive_sib[!is.na(M_ID), .(M_ID)]

#father_educ <- open_dataset(paste0(CleanData.dir, "/highest_educ.parquet")) %>%
#  filter(ID %in% father_dt$F_ID) %>% 
#  collect() %>% 
#  as.data.table()
#setnames(father_educ, c("ID", "highest_edu"), c("F_ID", "father_edu"))

#mother_educ <- open_dataset(paste0(CleanData.dir, "/highest_educ.parquet")) %>%
#  filter(ID %in% mother_dt$M_ID) %>% 
#  collect() %>% 
#  as.data.table()
#setnames(mother_educ, c("ID", "highest_edu"), c("M_ID", "mother_edu"))

#father_dt <- merge(father_dt, father_educ, all.x = TRUE, by = "F_ID")
#mother_dt <- merge(mother_dt, mother_educ, all.x = TRUE, by = "M_ID")


extensive_sib <- merge(extensive_sib, father_dt, all.x = TRUE, by = "F_ID")
extensive_sib <- merge(extensive_sib, mother_dt, all.x = TRUE, by = "M_ID")

# Add birth city
#B_CITY_dt <- open_dataset(paste0(CleanData.dir, "/birth_city.parquet")) %>% 
#  select(ID, Birth_County, county_code) %>% collect() %>% as.data.table() %>% distinct()

#extensive_sib <- merge(extensive_sib, B_CITY_dt, all.x = TRUE, by = "ID")
#rm(B_CITY_dt); gc()

extensive_sib <- distinct(extensive_sib)

# somebody has two weights
extensive_sib[, count := .N, by = ID]
extensive_sib <- extensive_sib[count == 1]
extensive_sib[, count := NULL]

# generate some more variables

extensive_sib[, patient_birth_month := patient_BD %% 10000 %/% 100]
extensive_sib[, male := ifelse(ID_S == 1, 1, 0)]
extensive_sib[, patient_male := ifelse(patient_S == 1, 1, 0)]
extensive_sib[, `:=`(ID_S = NULL, patient_S = NULL)]
extensive_sib[, `:=`(father_count = NULL, mother_count = NULL, both_parents = NULL)]
write_parquet(extensive_sib, paste0(output.dir, "/sample.parquet"))
head(extensive_sib)
rm(older_ID, younger_ID, pers_info); gc()

# ---- Add Wage ----
extensive_sib <- open_dataset(paste0(output.dir, "/sample.parquet")) %>% 
  collect() %>% 
  as.data.table()
head(extensive_sib)

pers_info <- open_dataset(paste0(CleanData.dir, "/PersInfo_death.parquet")) %>%
  select(ID, D_DATE) %>%
  collect() %>%
  as.data.table()
head(pers_info)

extensive_sib <- merge(extensive_sib, pers_info, all.x = TRUE, by = "ID")
rm(pers_info); gc()
extensive_sib[, death_year := D_DATE %/% 10000]

father_dt <- extensive_sib[, .(ID = F_ID,
                               first_YM = 100*birth_year + birth_month,
                               second_YM = 100*patient_birth_year + patient_birth_month)]
mother_dt <- extensive_sib[, .(ID = M_ID,
                               first_YM = 100*birth_year + birth_month,
                               second_YM = 100*patient_birth_year + patient_birth_month)]

# event-year wage of father
income_raw_f <- get_event_year_wage(father_dt[, .(ID, first_YM)], -12, -1, "first_YM")
income_raw_f <- remove_blacklist(income_raw_f)
income_f <- income_raw_f[, .(f_first_event_income = sum(ID1_AMT)), by = ID]

income_raw_f <- get_event_year_wage(father_dt[, .(ID, second_YM)], -12, -1, "second_YM")
income_raw_f <- remove_blacklist(income_raw_f)
income_f2 <- income_raw_f[, .(f_second_event_income = sum(ID1_AMT)), by = ID]
head(income_f)

income_f <- merge(income_f, income_f2, all = TRUE, by = "ID")
extensive_sib <- merge(extensive_sib, income_f, all.x = TRUE, by.x = "F_ID", by.y = "ID")

# event-year wage of mother
income_raw_m <- get_event_year_wage(mother_dt[, .(ID, first_YM)], -12, -1, "first_YM")
income_raw_m <- remove_blacklist(income_raw_m)
income_m <- income_raw_m[, .(m_first_event_income = sum(ID1_AMT)), by = ID]

income_raw_m <- get_event_year_wage(mother_dt[, .(ID, second_YM)], -12, -1, "second_YM")
income_raw_m <- remove_blacklist(income_raw_m)
income_m2 <- income_raw_m[, .(m_second_event_income = sum(ID1_AMT)), by = ID]

head(income_m2)

income_m <- merge(income_m, income_m2, all = TRUE, by = "ID")
extensive_sib <- merge(extensive_sib, income_m, all.x = TRUE, by.x = "M_ID", by.y = "ID")

F_ID_list <- father_dt$ID
M_ID_list <- mother_dt$ID
parents_pers_info <- open_dataset(paste0(CleanData.dir, "/PersInfo_death.parquet")) %>%
  filter((ID %in% F_ID_list | ID %in% M_ID_list) & valid) %>% 
  select(ID, ID_S, ID_BIRTHDAY) %>% 
  collect() %>% 
  as.data.table()
parents_pers_info[, birth_year := ID_BIRTHDAY %/% 10000]

extensive_sib <- merge(extensive_sib, parents_pers_info[, .(F_ID = ID, f_birth_year = birth_year)],
                       all.x = TRUE, by = "F_ID")
extensive_sib <- merge(extensive_sib, parents_pers_info[, .(M_ID = ID, m_birth_year = birth_year)],
                       all.x = TRUE, by = "M_ID")

extensive_sib[is.na(f_first_event_income) & birth_year >= 2001, f_first_event_income := 0]
extensive_sib[is.na(m_first_event_income) & birth_year >= 2001, m_first_event_income := 0]
extensive_sib[is.na(f_second_event_income) & patient_birth_year >= 2001, f_second_event_income := 0]
extensive_sib[is.na(m_second_event_income) & patient_birth_year >= 2001, m_second_event_income := 0]
extensive_sib[, `:=`(f_first_event_work = ifelse(f_first_event_income == 0, 0, 1),
                     m_first_event_work = ifelse(m_first_event_income == 0, 0, 1),
                     f_second_event_work = ifelse(f_second_event_income == 0, 0, 1),
                     m_second_event_work = ifelse(m_second_event_income == 0, 0, 1))]

extensive_sib[, `:=`(family_first_event_income = f_first_event_income + m_first_event_income,
                     family_second_event_income = f_second_event_income + m_second_event_income,
                     both_first_event_work = f_first_event_work * m_first_event_work,
                     both_second_event_work = f_second_event_work * m_second_event_work,
                     either_first_event_work = ifelse(f_first_event_work == 1 | m_first_event_work == 1,
                                                      1, 0),
                     either_second_event_work = ifelse(f_second_event_work == 1 | m_second_event_work == 1,
                                                       1, 0))]
# extensive_sib[, `:=`(f_bachelor = ifelse(father_edu >= 16, 1, 0),
#           m_bachelor = ifelse(father_edu >= 16, 1, 0))]
head(extensive_sib)
write_parquet(extensive_sib, paste0(output.dir, "/sample.parquet"))

# ---------- Parental Health ----------
f_health_list <- list()
for (y in 89:110) {
  start_timer()
  message(y)
  opdte <- read_opdte(y, extensive_sib$F_ID)
  
  setnames(opdte,
           c("ID",
             "outpatient_count",
             "total_dot",
             "chronic",
             "major",
             "eme",
             "useless"),
           c("F_ID",
             "f_outpatient_count",
             "f_total_dot",
             "f_chronic",
             "f_major",
             "f_eme",
             "f_useless"))
  
  father_health <- merge(extensive_sib[birth_year == y + 1911, .(ID, F_ID)], opdte, all.x = TRUE, by = "F_ID")
  father_health[, `:=`(f_outpatient_count = ifelse(is.na(f_outpatient_count), 0, f_outpatient_count),
                       f_total_dot = ifelse(is.na(f_total_dot), 0, f_total_dot),
                       f_chronic = ifelse(is.na(f_chronic), 0, f_chronic),
                       f_major = ifelse(is.na(f_major), 0, f_major),
                       f_eme = ifelse(is.na(f_eme), 0, f_eme),
                       f_useless = ifelse(is.na(f_useless), 0, f_useless))]
  f_health_list <- c(f_health_list, list(father_health))
  end_timer()
}
father_health <- rbindlist(f_health_list)
head(father_health)
write_parquet(father_health, paste0(output.dir, "/father_health.parquet"))
rm(father_health, f_health_list, opdte); gc()

m_health_list <- list()
for (y in 89:110) {
  start_timer()
  message(y)
  opdte <- read_opdte(y, extensive_sib$M_ID)
  setnames(opdte,
           c("ID",
             "outpatient_count",
             "total_dot",
             "chronic",
             "major",
             "eme",
             "useless"),
           c("M_ID",
             "m_outpatient_count",
             "m_total_dot",
             "m_chronic",
             "m_major",
             "m_eme",
             "m_useless"))
  
  mother_health <- merge(extensive_sib[birth_year == y + 1911, .(ID, M_ID)], opdte, all.x = TRUE, by = "M_ID")
  mother_health[, `:=`(m_outpatient_count = ifelse(is.na(m_outpatient_count), 0, m_outpatient_count),
                       m_total_dot = ifelse(is.na(m_total_dot), 0, m_total_dot),
                       m_chronic = ifelse(is.na(m_chronic), 0, m_chronic),
                       m_major = ifelse(is.na(m_major), 0, m_major),
                       m_eme = ifelse(is.na(m_eme), 0, m_eme),
                       m_useless = ifelse(is.na(m_useless), 0, m_useless))]
  m_health_list <- c(m_health_list, list(mother_health))
  end_timer()
}
mother_health <- rbindlist(m_health_list)
head(mother_health)
write_parquet(mother_health, paste0(output.dir, "/mother_health.parquet"))
