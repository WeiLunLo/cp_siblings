###############################################################################
# Results for CP
###############################################################################
rm(list = ls()); gc()

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

CPI <- fread(paste0(ref.dir, "/CPI.csv")) %>% as.data.table()

### ICD code
ICD_chronic <- fread(paste0(ref.dir,"/icd_chronic.csv"))
icd_chronic <- c(ICD_chronic$`ICD-9-CM`, ICD_chronic$`ICD-10-CM`)
icd_chronic <- str_replace_all(icd_chronic, "\\.", "")

ICD_major <- fread(paste0(ref.dir,"/icd_major.csv"))
icd_major <- c(ICD_major$`ICD-9-CM`, ICD_major$`ICD-10-CM`)
icd_major <- str_replace_all(icd_major, "\\.", "")

icd_useless <- c("Z0389", "Z711")

rm(ICD_chronic, ICD_major)


###############################################################################
### MAIN
###############################################################################
# dt <- read_weighted_siblings()
dt <- open_dataset(paste0(output.dir, "/sample.parquet")) %>%
  collect() %>%
  as.data.table()

head(dt)
dt[, death_year := D_DATE %/% 10000]

# add patient death info
dt <- add_death_info_patient(dt, dt$patient)
dt[, patient_death_year := patient_D_DATE %/% 10000]
############
## income
############
dt_w_AMT <- lapply((89:110), merge_with_AMT) %>% rbindlist()
dt_w_AMT <- dt_w_AMT[age >= 25 & age < 65]
dt_w_AMT[, CP_male := CP*male]
dt_w_AMT[, CP_patient_male := CP*patient_male]
dt_w_AMT[, `:=`(f_bachelor = ifelse(father_edu >= 16, 1, 0),
                m_bachelor = ifelse(mother_edu >= 16, 1, 0))]
write_parquet(dt_w_AMT, paste0(output.dir, "/dt_w_AMT.parquet"))

rm(dt_w_AMT); gc()

###############
## education
###############

highest_educ <- open_dataset(paste0(CleanData.dir, "/highest_educ.parquet")) %>% 
  filter(!is.na(ID) & ID != "" & (ID %in% dt$ID | ID %in% dt$patient)) %>% 
  collect() %>%
  as.data.table()

#marr <- open_dataset(paste0(CleanData.dir, "/marriage_long.parquet")) %>% 
#  filter(ID != "" & ID %in% dt$ID) %>% 
#  collect() %>% 
#  as.data.table()

dt_w_educ <- merge(dt, highest_educ, by = "ID")
dt_w_educ[, age := 110 + 1911 - birth_year]
dt_w_educ[, `:=`(f_bachelor = ifelse(father_edu >= 16, 1, 0),
                 m_bachelor = ifelse(mother_edu >= 16, 1, 0))]

# younger sibling's edu level
dt_w_educ <- merge(dt_w_educ, highest_educ, by.x = "patient", by.y = "ID", all.x = T)

dt_w_educ[, patient_age := 110 + 1911 - patient_birth_year]
dt_w_educ <- rename(dt_w_educ, 
                    highest_edu = highest_edu.x, patient_highest_edu = highest_edu.y)
dt_w_educ[, `:=`(bachelor = ifelse(highest_edu >= 16, 1, 0),
                 patient_bachelor = ifelse(is.na(patient_highest_edu),
                                           NA,
                                           ifelse(patient_highest_edu >= 16, 1, 0))
                 )]

write_parquet(dt_w_educ, paste0(output.dir, "/dt_w_educ.parquet"))

rm(dt_w_educ); gc()

###############
## marriage
###############

dt_marr <- open_dataset(paste0(CleanData.dir, "/marriage_long.parquet")) %>% 
  filter(ID != "" & ID %in% dt$ID) %>% 
  filter(year %in% 89:110) %>%
  collect() %>% 
  as.data.table()


merge_with_marr <- function(y) {
  message(y)
  marr <- merge(dt[is.na(D_DATE) | death_year > y + 1911], dt_marr[year==y], 
                     all.x = TRUE, by = "ID")
  marr <- marr[ , year:=y]
  marr[, marr := ifelse(MARR == 2, 1, 0)] # keep NA sample
  marr[, age := year + 1911 - birth_year]
  
  return(marr)
}

dt_w_marr <- lapply((89:110), merge_with_marr) %>% rbindlist()
dt_w_marr <- dt_w_marr[age>0]
dt_w_marr[, `:=`(f_bachelor = ifelse(father_edu >= 16, 1, 0),
                 m_bachelor = ifelse(mother_edu >= 16, 1, 0))]

write_parquet(dt_w_marr, paste0(output.dir, "/dt_w_marr.parquet"))

rm(dt_w_marr); gc()

###############
## outpatient
###############
dt_w_opdte <- lapply((89:110), merge_with_opdte) %>% rbindlist()

dt_w_opdte <- dt_w_opdte[age > 0]
dt_w_opdte[, `:=`(f_bachelor = ifelse(father_edu >= 16, 1, 0),
                  m_bachelor = ifelse(mother_edu >= 16, 1, 0))]

write_parquet(dt_w_opdte, paste0(output.dir, "/dt_w_opdte.parquet"))

head(dt_w_opdte)

rm(dt_w_opdte); gc()

###############
## fertility
###############

dt_fertility <- fread(paste0(parsed_path, "/H_MOHW_MCHD109.csv"))
# remove invalid records
dt_fertility <- dt_fertility[str_length(Yr)==4]
# keep only those born after 2004
dt_fertility <- dt_fertility[ , child_birth_year:= as.integer(Yr)][child_birth_year>=2004]

fert_father <- dt_fertility[!is.na(ID_F) & ID_F!="", list(fertility = .N), by = .(ID_F, child_birth_year)]
fert_mother <- dt_fertility[!is.na(ID_M) & ID_M!="", list(fertility = .N), by = .(ID_M, child_birth_year)]
fert_count <- rbind(rename(fert_father, ID = ID_F),
                    rename(fert_mother, ID = ID_M))

merge_with_fertility <- function(y) {
  message(y)
  
  fert <- merge(dt[is.na(D_DATE) | death_year > y + 1911], 
                fert_count[child_birth_year== y + 1911, .(ID, fertility)], 
                all.x = TRUE, by = "ID")
  fert[ , year:= y]
  fert[is.na(fertility), fertility := 0]
  fert[, age := year + 1911 - birth_year]
  
  return(fert)
}

dt_w_fert <- lapply((93:109), merge_with_fertility) %>% rbindlist()
dt_w_fert <- dt_w_fert[age > 0]
dt_w_fert[, `:=`(f_bachelor = ifelse(father_edu >= 16, 1, 0),
                 m_bachelor = ifelse(mother_edu >= 16, 1, 0))]

write_parquet(dt_w_fert, paste0(output.dir, "/dt_w_fert.parquet"))
rm(dt_w_fert); gc()

# patch -------------------------------------------------------------------

# living together
# use HHREG 91 - 110
# here we deal with multiple registry
# it takes 4 hours...
dt_living <- data.table()
for (yr in 91:110) {
  hhreg_ID <- open_dataset(paste0(output.dir, "/sample.parquet")) %>%
    select(ID, patient) %>%
    # patient_HOST
    merge(
      open_dataset(paste0(par.dir, "/chinese_hhreg/S_MOI_HHREG", yr, ".parquet")) %>%
        filter(!is.na(HOST)) %>%
        select(ID, HOST) %>%
        distinct(),
      by.x = c("ID"), by.y = c("ID")
    )
  
  hhreg_patient <- open_dataset(paste0(output.dir, "/sample.parquet")) %>%
    select(ID, patient) %>%
    # patient_HOST
    merge(
      open_dataset(paste0(par.dir, "/chinese_hhreg/S_MOI_HHREG", yr, ".parquet")) %>%
        filter(!is.na(HOST)) %>%
        select(ID, HOST) %>%
        distinct(),
      by.x = c("patient"), by.y = c("ID")
    ) %>%
    group_by(patient) %>%
    mutate(registry_n = n()) %>%
    collect() %>%
    as.data.table()
  
  # match multi-ID-registry to single-patient-registry in each loop
  
  dt_living_yr <- data.table()
  for (registry in unique(hhreg_patient$registry_n)){
    print(registry)
    print(Sys.time())
    dt_temp <- merge(
      hhreg_ID,
      hhreg_patient[, .SD[registry], by = .(patient)],
      by = c("ID", "patient", "HOST")
    ) %>%
      mutate(year = yr,
             living_together = T) %>%
      collect() %>%
      as.data.table()
    dt_living_yr <- rbind(dt_living_yr, dt_temp)
    print(Sys.time())
    rm(dt_temp)
  }
  dt_living <- rbind(dt_living, dt_living_yr)
  rm(hhreg_ID, hhreg_patient)
  gc()
}
write_parquet(dt_living, paste0(output.dir, "/dt_living.parquet"))

# append living status to dataset

# Labor
dt_w_AMT <- open_dataset(paste0(output.dir, "/dt_w_AMT.parquet")) %>%
  # append patient's death information
  left_join(open_dataset(paste0(CleanData.dir, "/persInfo.parquet")) %>%
              mutate(patient_death_year = D_DATE %/% 10000) %>%
              select(ID, patient_death_year),
            by = c("patient" = "ID")) %>%
  # append living status
  left_join(open_dataset(paste0(output.dir, "/dt_living.parquet")) %>%
              select(ID, patient, year, living_together), 
            by = c("ID", "patient", "year")) %>%
  collect() %>%
  as.data.table()

dt_w_AMT[year >=91 & year <=110 & is.na(living_together) &
               (is.na(patient_death_year) | patient_death_year > 1911 + year), 
             living_together:= FALSE]

write_parquet(dt_w_AMT, paste0(output.dir, "/dt_w_AMT.parquet"))
rm(dt_w_AMT); gc()

# WIP ---------------------------------------------------------------------


# Education

dt_w_educ <- open_dataset(paste0(output.dir, "/dt_w_educ.parquet")) %>%
  # append patient's death information
  left_join(open_dataset(paste0(CleanData.dir, "/persInfo.parquet")) %>%
              mutate(patient_death_year = D_DATE %/% 10000) %>%
              select(ID, patient_death_year),
            by = c("patient" = "ID")) %>%
  # append living status
  left_join(open_dataset(paste0(output.dir, "/dt_living.parquet")) %>%
              select(ID, patient, year, living_together), 
            by = c("ID", "patient", "year")) %>%
  collect() %>%
  as.data.table()

dt_w_educ[year >=91 & year <=110 & is.na(living_together) &
               (is.na(patient_death_year) | patient_death_year > 1911 + year), 
             living_together:= FALSE]

write_parquet(dt_w_educ, paste0(output.dir, "/dt_w_educ.parquet"))
rm(dt_w_educ); gc()

# WIP ---------------------------------------------------------------------


# Marriage
dt_w_marr <- open_dataset(paste0(output.dir, "/dt_w_marr.parquet")) %>%
  # append patient's death information
  left_join(open_dataset(paste0(CleanData.dir, "/persInfo.parquet")) %>%
              mutate(patient_death_year = D_DATE %/% 10000) %>%
              select(ID, patient_death_year),
            by = c("patient" = "ID")) %>%
  # append living status
  left_join(open_dataset(paste0(output.dir, "/dt_living.parquet")) %>%
              select(ID, patient, year, living_together), 
            by = c("ID", "patient", "year")) %>%
  collect() %>%
  as.data.table()

dt_w_marr[year >=91 & year <=110 & is.na(living_together) &
               (is.na(patient_death_year) | patient_death_year > 1911 + year), 
             living_together:= FALSE]

write_parquet(dt_w_marr, paste0(output.dir, "/dt_w_marr.parquet"))
rm(dt_w_marr); gc()

# Fertility

dt_w_fert <- open_dataset(paste0(output.dir, "/dt_w_fert.parquet")) %>%
  # append patient's death information
  left_join(open_dataset(paste0(CleanData.dir, "/persInfo.parquet")) %>%
              mutate(patient_death_year = D_DATE %/% 10000) %>%
              select(ID, patient_death_year),
            by = c("patient" = "ID")) %>%
  # append living status
  left_join(open_dataset(paste0(output.dir, "/dt_living.parquet")) %>%
              select(ID, patient, year, living_together), 
            by = c("ID", "patient", "year")) %>%
  collect() %>%
  as.data.table()

dt_w_fert[year >=91 & year <=110 & is.na(living_together) &
               (is.na(patient_death_year) | patient_death_year > 1911 + year), 
             living_together:= FALSE]

write_parquet(dt_w_fert, paste0(output.dir, "/dt_w_fert.parquet"))
rm(dt_fert); gc()

# Outpatient

dt_w_opdte <- open_dataset(paste0(output.dir, "/dt_w_opdte.parquet")) %>%
  # append patient's death information
  left_join(open_dataset(paste0(CleanData.dir, "/persInfo.parquet")) %>%
              mutate(patient_death_year = D_DATE %/% 10000) %>%
              select(ID, patient_death_year),
            by = c("patient" = "ID")) %>%
  # append living status
  left_join(open_dataset(paste0(output.dir, "/dt_living.parquet")) %>%
              select(ID, patient, year, living_together), 
            by = c("ID", "patient", "year")) %>%
  collect() %>%
  as.data.table()

dt_w_opdte[year >=91 & year <=110 & is.na(living_together) &
               (is.na(patient_death_year) | patient_death_year > 1911 + year), 
             living_together:= FALSE]

write_parquet(dt_w_opdte, paste0(output.dir, "/dt_w_opdte.parquet"))
rm(dt_w_opdte); gc()

# testing -----------------------------------------------------------------
###############
## outpatient New
###############
dt_w_opdte <- lapply((89:110), merge_with_opdte) %>% rbindlist()

dt_w_opdte <- dt_w_opdte[age > 0]
dt_w_opdte[, `:=`(f_bachelor = ifelse(father_edu >= 16, 1, 0),
                  m_bachelor = ifelse(mother_edu >= 16, 1, 0))]

write_parquet(dt_w_opdte, paste0(output.dir, "/dt_w_opdte_md.parquet"))

head(dt_w_opdte)

rm(dt_w_opdte); gc()