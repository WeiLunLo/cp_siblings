### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Results for CP
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
rm(list = ls()); gc()

library(data.table)
library(arrow)
library(lmtest)
library(sandwich)
library(dplyr)
library(ggplot2)
library(fixest)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PREAMBLE ----------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

setwd("E:/H114028/wllo/cp_siblings")

source("code/utils/arrow_helpers.R")
source("code/utils/functions_for_eq1_and_eq2.R")


output.dir <- "data/processed"
par.dir <- "E:/H114028/data/parquet"
result.dir <- "results"
CleanData.dir <- "E:/H114028/CleanData/processed_data"
figure.dir <- paste0(result.dir, "/figure")
ref.dir <- "reference"
YEARS <- 89L:110L
arrow_path <- par.dir <- "E:/H114028/data/parquet"

GET_PARENTS_DETAIL <- FALSE

CPI <- fread(paste0(ref.dir, "/CPI.csv")) %>% as.data.table()

### ICD code
ICD_chronic <- fread(paste0(ref.dir,"/icd_chronic_self_made.csv"))
icd_chronic <- c(ICD_chronic$`ICD-9-CM`, ICD_chronic$`ICD-10-CM`)
icd_chronic <- str_replace_all(icd_chronic, "\\.", "")

ICD_major <- fread(paste0(ref.dir,"/icd_severe_self_made.csv"))
icd_major <- c(ICD_major$`ICD-9-CM`, ICD_major$`ICD-10-CM`)
icd_major <- str_replace_all(icd_major, "\\.", "")

icd_useless <- c("Z0389", "Z711")

rm(ICD_chronic, ICD_major)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## MAIN ----------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Read in relation dt for info on first link # TEST
relation_dt <- open_dataset("E:/H114028/CleanData/processed_data/enrol_relation.parquet") |> 
  select(ID, F_RELATION_YR, M_RELATION_YR) |> 
  collect() |> 
  setDT()

# dt <- read_weighted_siblings()
dt <- open_dataset(paste0(output.dir, "/sample.parquet")) %>%
  collect() %>%
  setDT()
head(dt)

# add patient death info
dt <- add_death_info_patient(dt, dt$patient)
dt[, patient_death_year := patient_D_DATE %/% 10000]

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## Income ----------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
dt_w_AMT <- lapply((89:110), merge_with_AMT) %>% rbindlist()
dt_w_AMT <- dt_w_AMT[age >= 25 & age < 65]
dt_w_AMT[, CP_male := CP*male]
dt_w_AMT[, CP_patient_male := CP*patient_male]
 
# we don't have education level information now
# dt_w_AMT[, `:=`(f_bachelor = ifelse(father_edu >= 16, 1, 0),
#                 m_bachelor = ifelse(mother_edu >= 16, 1, 0))]

# TEST and review; Change dt_w_AMT to dt if works
find_early_link <- function(dt_w_AMT, relation_dt) {
  dt_w_AMT <- merge(dt_w_AMT, replation_dt, by = "ID")
  dt_w_AMT[, F_identify_age := F_RELATION_YR - birth_year]
  dt_w_AMT[, M_identify_age := M_RELATION_YR - birth_year]
  dt_w_AMT[, `:=`(F_RELATION_YR = NULL,
                  M_RELATION_YR = NULL)]
  
  dt_w_AMT <- merge(dt_w_AMT, replation_dt[, .(patient = ID, F_RELATION_YR, M_RELATION_YR)], by = "patient")
  dt_w_AMT[, p_F_identify_age := F_RELATION_YR - patient_birth_year]
  dt_w_AMT[, p_M_identify_age := M_RELATION_YR - patient_birth_year]
  
  # Flag older siblings and patients (younger siblings) that found their parents early
  dt_w_AMT[, flag_early_link := as.numeric((F_identify_age <= 20 | is.na(F_ID)) &
                                             (M_identify_age <= 20 | is.na(M_ID)))]
  dt_w_AMT[, flag_p_early_link := as.numeric((p_F_identify_age <= 20 | is.na(F_ID)) &
                                               (p_M_identify_age <= 20 | is.na(M_ID)))]
  
  return(dt_w_AMT)
}

dt_w_AMT <- find_early_link(dt_w_AMT, relation_dt)


write_parquet(dt_w_AMT, paste0(output.dir, "/dt_w_AMT.parquet"))
 
rm(dt_w_AMT); gc()



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## Outpatient ----------
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
dt_w_opdte <- lapply((89:110), merge_with_opdte) %>% rbindlist()

dt_w_opdte <- dt_w_opdte[age > 0]
# No education level data
# dt_w_opdte[, `:=`(f_bachelor = ifelse(father_edu >= 16, 1, 0),
#                   m_bachelor = ifelse(mother_edu >= 16, 1, 0))]

dt_w_opdte <- find_early_link(dt_w_opdte, relation_dt) # TEST

write_parquet(dt_w_opdte, paste0(output.dir, "/dt_w_opdte.parquet"))

head(dt_w_opdte)

rm(dt_w_opdte); gc()


# patch -------------------------------------------------------------------

# living together
# before: use HHREG 91 - 110
# here we deal with multiple registry
# it takes 4 hours...

# now I am working on using outpatient/drug location or enroll instead (home_city)

# WIP ---------------------------------------------------------------------

source("E:/H114028/wllo/cp_siblings/code/result.R")

source("E:/H114028/wllo/cp_siblings/code/hete_result.R")