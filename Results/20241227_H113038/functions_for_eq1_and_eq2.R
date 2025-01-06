###############################################################################
# Functions For Result of EQ1 & EQ2
###############################################################################
library(data.table)
library(arrow)
library(lmtest)
library(sandwich)
library(dplyr)
library(ggplot2)
library(fixest)
library(stringr)

source("E:/H113038/arrow/arrow_20231023/arrow_helpers.R")

Disk <- "E"
ref.dir <- paste0(Disk, ":/H113038/yangc/ref")
par.dir <- paste0(Disk, ":/H113038/parquet")

### ICD code
ICD_chronic <- fread(paste0(ref.dir,"/icd_chronic.csv"))
icd_chronic <- c(ICD_chronic$`ICD-9-CM`, ICD_chronic$`ICD-10-CM`)
icd_chronic <- str_replace_all(icd_chronic, "\\.", "")

ICD_major <- fread(paste0(ref.dir,"/icd_major.csv"))
icd_major <- c(ICD_major$`ICD-9-CM`, ICD_major$`ICD-10-CM`)
icd_major <- str_replace_all(icd_major, "\\.", "")

icd_useless <- c("Z0389", "Z711")

rm(ICD_chronic, ICD_major)

icd_mental <- c(290:319, paste0("F0", 1:9), paste0("F", 10:99))
icd_dep <-  c("29383", "296", "3004", "311", "F063", paste0("F", 30:39))
icd_mental_excluded <-  c(290:292, "2930", "2931", "294", "310", 303:305,
                    paste0("F0", 1:5), "F07", "F09", "F482")
icd_anxio <-  c("29384", "300", "308", "30981", "F064", paste0("F4", 0:2),
              paste0("F43", 0:1), paste0("F48", 8:9))
icd_psych <-  c("295", "297", "298", "301","F20", "F25", 
              "F22", "F34", 'F60')

CPI <- fread("E:/H113038/yangc/ref/CPI.csv") %>% as.data.table()
### Simple Functions ##########################################################
print.time.diff <- function(now) { cat("\t",'time:',round(difftime(Sys.time(),now,units='mins'),2),"mins\n") }
start_timer <- function() { start_time <<- Sys.time() }
end_timer <- function() { end_time <<- Sys.time(); print(end_time - start_time) }


###############################################################################
### FUNCTIONS #################################################################
###############################################################################

add_age_group_variables <- function(dt, y) {
  print(y)
  dt[, age := (y + 1911) - birth_year]
  dt[, age_group := age %/% 5]
  dt[, age_group_chr := paste0(age_group*5, "~", age_group*5 + 4)]
  
  return(dt)
}

# ---------- AMT ----------
read_AMT <- function(y, l) {
  AMT <- open_dataset(paste0(CleanData.dir, "/annualFile_AMT", y, ".parquet")) %>% 
    filter(ID %in% l) %>% 
    collect() %>% 
    as.data.table()
  AMT[, year := y]
  cpi <- CPI[year == y]$CPI
  print(paste("CPI:", cpi))
  AMT[, FTAMT := (FTAMT/(cpi/100))]
  AMT[, PTAMT := (PTAMT/(cpi/100))]
  return(AMT)
}

merge_with_AMT <- function(y) {
  message(y)
  AMT <- read_AMT(y, dt$ID)
  dt_w_AMT <- merge(dt[is.na(D_DATE) | death_year > y + 1911], AMT, all.x = TRUE, by = "ID")
  dt_w_AMT[, income := FTAMT + PTAMT]
  dt_w_AMT[, income := ifelse(PRIMER_OCCUPATION %in% c("C", "E", "G", "H", "I", "J", "K"), 0, income)]
  dt_w_AMT[, income := ifelse(is.na(income), 0, income)]
  dt_w_AMT[, work := ifelse(income == 0, 0, 1)]
  dt_w_AMT[, year := y]
  dt_w_AMT[, age := year + 1911 - birth_year]
  
  return(dt_w_AMT)
}


# ---------- Health ----------
read_ipdte <- function(y, l){
  message(paste("read inpatient:", y))
  # Ipdte
  dtIpdte <- open_dataset(paste0(arrow_path, "/H_NHI_IPDTE", y, ".parquet")) %>%
    filter(ID %in% l) %>% 
    select(ID, MED_DOT) %>% 
    collect() %>% 
    as.data.table()
  dtIpdte <- dtIpdte[, .(ipdte_dot = sum(MED_DOT)), by = ID]
  
  return(dtIpdte)
}

read_druge <- function(y, l){
  message(paste("read drug:", y))
  dtDruge <- open_nhi("DRUGE", y) %>%
    filter(ID %in% l) %>% 
    select(ID, T_APPL_DOT) %>% 
    collect() %>% 
    as.data.table()
  
  dtDruge <- dtDruge[, .(drug_dot = sum(T_APPL_DOT)), by = ID]
  return(dtDruge)
}

read_opdte <- function(y, l) {
  message(paste("read outpatient:", y))
  
  dt_outpatient <- open_nhi("OPDTE", y) %>% 
    filter(!is.na(ID) & ID != "" & ID %in% l) %>% 
    select(ID, T_DOT, ICD9CM_1, ICD9CM_2, ICD9CM_3, HOSP_ID, CASE_TYPE) %>%
    collect() %>% 
    as.data.table()
  
  dt_outpatient[, chronic := ifelse((ICD9CM_1 %in% icd_chronic |
                                       ICD9CM_2 %in% icd_chronic | 
                                       ICD9CM_3 %in% icd_chronic),
                                    1, 0)]
  dt_outpatient[, major := ifelse((ICD9CM_1 %in% icd_major |
                                     ICD9CM_2 %in% icd_major | 
                                     ICD9CM_3 %in% icd_major),
                                  1, 0)]
  dt_outpatient[, eme := ifelse(CASE_TYPE == "02", 1, 0)]
  
  # examine eme source
  #dt_outpatient[eme==1, `:=` (eme_1 := str_sub(ICD9CM_1, 1, 1),
  #                            eme_2 := str_sub(ICD9CM_2, 1, 1),
  #                            eme_3 := str_sub(ICD9CM_3, 1, 1))]
  
  # mental diseases 
  
  # general mental diseases
  dt_outpatient[ , md:= 0]
  dt_outpatient[str_sub(ICD9CM_1, 1, 3) %in% icd_mental | 
                  str_sub(ICD9CM_2, 1, 3) %in% icd_mental |
                  str_sub(ICD9CM_3, 1, 3) %in% icd_mental, md:= 1]
  
  # Schizophrenia
  dt_outpatient[ , psy:= 0]
  dt_outpatient[str_sub(ICD9CM_1, 1, 3) %in% icd_psych | 
                str_sub(ICD9CM_2, 1, 3) %in% icd_psych |
                str_sub(ICD9CM_3, 1, 3) %in% icd_psych, psy:= 1]
  
  # Anxiety disorders
  dt_outpatient[ , anx:= 0]
  dt_outpatient[str_sub(ICD9CM_1, 1, 3) %in% icd_anxio | 
                str_sub(ICD9CM_2, 1, 3) %in% icd_anxio |
                str_sub(ICD9CM_3, 1, 3) %in% icd_anxio |
                str_sub(ICD9CM_1, 1, 4) %in% icd_anxio | 
                str_sub(ICD9CM_2, 1, 4) %in% icd_anxio |
                str_sub(ICD9CM_3, 1, 4) %in% icd_anxio |
                str_sub(ICD9CM_1, 1, 5) %in% icd_anxio | 
                str_sub(ICD9CM_2, 1, 5) %in% icd_anxio |
                str_sub(ICD9CM_3, 1, 5) %in% icd_anxio , anx:= 1] 
  # Mood disorders/Depression
  dt_outpatient[ , dep:= 0]
  dt_outpatient[str_sub(ICD9CM_1, 1, 3) %in% icd_dep | 
                str_sub(ICD9CM_2, 1, 3) %in% icd_dep |
                str_sub(ICD9CM_3, 1, 3) %in% icd_dep |
                str_sub(ICD9CM_1, 1, 4) %in% icd_dep | 
                str_sub(ICD9CM_2, 1, 4) %in% icd_dep |
                str_sub(ICD9CM_3, 1, 4) %in% icd_dep |
                str_sub(ICD9CM_1, 1, 5) %in% icd_dep | 
                str_sub(ICD9CM_2, 1, 5) %in% icd_dep |
                str_sub(ICD9CM_3, 1, 5) %in% icd_dep , dep:= 1] 
  
  dt_outpatient[(str_sub(ICD9CM_1, 1, 3) %in% icd_mental_excluded | 
                   str_sub(ICD9CM_2, 1, 3) %in% icd_mental_excluded |
                   str_sub(ICD9CM_3, 1, 3) %in% icd_mental_excluded) &
                  (str_sub(ICD9CM_1, 1, 4) %in% icd_mental_excluded |
                     str_sub(ICD9CM_2, 1, 4) %in% icd_mental_excluded |
                     str_sub(ICD9CM_3, 1, 4) %in% icd_mental_excluded),
                `:=` (md = 0, psy = 0, anx = 0, dep = 0)]
  
  dt_outpatient[, useless_icd := ifelse((ICD9CM_1 %in% icd_useless |
                                           ICD9CM_2 %in% icd_useless | 
                                           ICD9CM_3 %in% icd_useless),
                                        1, 0)]
  
  dt_outpatient[, `:=`(ICD9CM_1 = NULL,
                       ICD9CM_2 = NULL,
                       ICD9CM_3 = NULL)]
  
  dt_outpatient <- dt_outpatient[, .(outpatient_count = .N,
                                     total_dot = sum(T_DOT),
                                     chronic = max(chronic),
                                     major = max(major),
                                     eme = sum(eme),
                                     md = max(md),
                                     md_visit = sum(md),
                                     psy = max(psy),
                                     psy_visit = sum(psy),
                                     anx = max(anx),
                                     anx_visit = sum(anx),
                                     dep = max(dep),
                                     dep_visit = sum(dep),
                                     useless = sum(useless_icd)),
                                 by = ID]
  return(dt_outpatient)
}

merge_with_opdte <- function(y) {
  message(y)
  start_timer()
  opdte <- read_opdte(y, dt$ID)
  ipdte <- read_ipdte(y, dt$ID)
  druge <- read_druge(y, dt$ID)
  
  print(nrow(opdte))
  dt_w_health <- merge(dt[birth_year >= y & (is.na(D_DATE) | death_year > y + 1911)],
                       opdte, all.x = TRUE, by = "ID")
  dt_w_health <- merge(dt_w_health, ipdte, all.x = TRUE, by = "ID")
  dt_w_health <- merge(dt_w_health, druge, all.x = TRUE, by = "ID")
  
  dt_w_health[, `:=`(outpatient_count = ifelse(is.na(outpatient_count), 0, outpatient_count),
                    total_dot = ifelse(is.na(total_dot), 0, total_dot),
                    chronic = ifelse(is.na(chronic), 0, chronic),
                    major = ifelse(is.na(major), 0, major),
                    eme = ifelse(is.na(eme), 0, eme),
                    md = ifelse(is.na(md), 0, md),
                    psy = ifelse(is.na(psy), 0, psy),
                    anx = ifelse(is.na(anx), 0, anx),
                    dep = ifelse(is.na(dep), 0, dep),
                    md_visit = ifelse(is.na(md_visit), 0, md_visit),
                    psy_visit = ifelse(is.na(psy_visit), 0, psy_visit),
                    anx_visit = ifelse(is.na(anx_visit), 0, anx_visit),
                    dep_visit = ifelse(is.na(dep_visit), 0, dep_visit),
                    useless = ifelse(is.na(useless), 0, useless),
                    ipdte_dot = ifelse(is.na(ipdte_dot), 0, ipdte_dot),
                    drug_dot = ifelse(is.na(drug_dot), 0, drug_dot))]
  dt_w_health[, year := y]
  dt_w_health[, age := year + 1911 - birth_year]
  end_timer()
  
  return(dt_w_health)
}

merge_with_opdte_patient <- function(y) {
  message(y)
  start_timer()
  opdte <- read_opdte(y, dt$patient)
  ipdte <- read_ipdte(y, dt$patient)
  druge <- read_druge(y, dt$patient)
  
  print(nrow(opdte))
  dt_w_health <- merge(dt[patient_birth_year >= y & (is.na(patient_D_DATE) | patient_death_year > y + 1911)],
                       opdte, all.x = TRUE, by = c("patient" = "ID"))
  dt_w_health <- merge(dt_w_health, ipdte, all.x = TRUE, by = c("patient" = "ID"))
  dt_w_health <- merge(dt_w_health, druge, all.x = TRUE, by = c("patient" = "ID"))
  
  dt_w_health[, `:=`(outpatient_count = ifelse(is.na(outpatient_count), 0, outpatient_count),
                     total_dot = ifelse(is.na(total_dot), 0, total_dot),
                     chronic = ifelse(is.na(chronic), 0, chronic),
                     major = ifelse(is.na(major), 0, major),
                     eme = ifelse(is.na(eme), 0, eme),
                     useless = ifelse(is.na(useless), 0, useless),
                     ipdte_dot = ifelse(is.na(ipdte_dot), 0, ipdte_dot),
                     drug_dot = ifelse(is.na(drug_dot), 0, drug_dot))]
  dt_w_health[, year := y]
  dt_w_health[, age := year + 1911 - birth_year]
  end_timer()
  
  return(dt_w_health)
}

# wip ---------------------------------------------------------------------


add_death_info_patient <- function(data, l) {
  # death
  dt_death_all <- data.table()
  for (year in 89:110) {
    f <- paste0(par.dir, "/H_OST_DEATH", year, ".parquet")
    print(paste("Start to Load:", f))
    death_info <- open_dataset(f) %>% 
      filter(ID %in% l) %>% 
      select(ID, D_DATE) %>%
      rename(patient_D_DATE = D_DATE) %>%
      collect() %>% 
      as.data.table()
    dt_death_all <- rbind(death_info, dt_death_all, fill=TRUE)
  }
  y <- merge(data, dt_death_all, all.x = TRUE, 
             by.x = c("patient"), by.y = c("ID"))
  return(y)
}

# ==============================================================================
# ICD needed
# Anxiety disorders: 293.84, 300.0/10/2/3/5/89/9, 308, 309.81 [ICD-9]
#                    F06.4, F40-42, F43.0/1, F48.8/9 [ICD-10]
# Depression = c("29383", "296", "3004", "311", "F063", paste0("F", 30:39))
# Psych = c("295", "297", "298", "301","F20", "F25", 
#              "F22", "F34", 'F60')
# Other metall illness: 293.89/9, 299, 300.11-19/6/7/81/82, 301, 302, 306, 307,
#                       309.0/1/2/3/4/82/83/89/9, 312-319 [ICD-9]
#                       F06.1/8, F43.2/8/9, F44, F45, F48.1, F50-F99
#
# All mental illness: 290-319 except cognitive disorder(290, 293.0/1, 294, 310)
#                     and substance-related disorders(291-292, 303-305) [ICD-9]
#                     F01-F99 except cognitive disorder(F01-F05, F07, F09, F48.2)
#                     and substance-related disorders(F10-F19) [ICD-10]
# ==============================================================================