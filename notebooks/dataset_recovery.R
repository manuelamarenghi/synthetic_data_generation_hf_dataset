library(dplyr)
library(lubridate)
library(readxl)
library(fuzzyjoin)
library(tidyr)
library(data.table)

load("C:\\path\\HF_orig.RData")

dataset <- data %>%
  filter(tipo_prest %in% c(41, 42, 43)) %>%
  mutate(
    data_studio_out = as.numeric(as.POSIXct(data_studio_out, tz = "UTC")),
  ) %>%
  
  ############# comorbidity score ####################

data <- data %>%
  mutate(
    alcohol = ifelse(is.na(alcohol), 0, alcohol),
    tumor = ifelse(is.na(tumor), 0, tumor),
    coagulopathy = ifelse(is.na(coagulopathy), 0, coagulopathy),
    anemia = ifelse(is.na(anemia), 0, anemia),
    metastatic = ifelse(is.na(metastatic), 0, metastatic),
    dementia = ifelse(is.na(dementia), 0, dementia),
    renal = ifelse(is.na(renal), 0, renal),
    hemiplegia = ifelse(is.na(hemiplegia), 0, hemiplegia),
    wtloss = ifelse(is.na(wtloss), 0, wtloss),
    arrhythmia = ifelse(is.na(arrhythmia), 0, arrhythmia),
    pulmonarydz = ifelse(is.na(pulmonarydz), 0, pulmonarydz),
    electrolytes = ifelse(is.na(electrolytes), 0, electrolytes),
    compdiabetes = ifelse(is.na(compdiabetes), 0, compdiabetes),
    liver = ifelse(is.na(liver), 0, liver),
    pvd = ifelse(is.na(pvd), 0, pvd),
    psychosis = ifelse(is.na(psychosis), 0, psychosis),
    pulmcirc = ifelse(is.na(pulmcirc), 0, pulmcirc),
    hivaids = ifelse(is.na(hivaids), 0, hivaids),
    hypertension = ifelse(is.na(hypertension), 0, hypertension),
    chf = ifelse(is.na(chf), 0, chf),
    ICD = ifelse(is.na(ICD), 0, ICD),
    SHOCK = ifelse(is.na(SHOCK), 0, SHOCK),
    CABG = ifelse(is.na(CABG), 0, CABG),
    PTCA = ifelse(is.na(PTCA), 0, PTCA),
    dec_intra = ifelse(is.na(dec_intra), 0, dec_intra)
  )

comorb_vars <- c(
  "alcohol","tumor","coagulopathy","anemia","metastatic",
  "dementia","renal","hemiplegia","wtloss","arrhythmia",
  "pulmonarydz","electrolytes","compdiabetes","liver",
  "pvd","psychosis","pulmcirc","hivaids","hypertension","chf"
)


comorb_patient_level <- data %>%
  group_by(COD_REG) %>%
  summarise(
    across(
      all_of(comorb_vars),
      ~ as.integer(max(.x, na.rm = TRUE) >= 1)
    ),
    .groups = "drop"
  )

score_mapping <- data.frame(
  condition_name = c("metastatic","alcohol","tumor","hivaids", "psychosis", "liver", 
                     "wtloss", "dementia", "hemiplegia", "coagulopathy", 
                     "electrolytes", "renal", "chf", "anemia", "compdiabetes", 
                     "hypertension", "pulmcirc", "pulmonarydz", 
                     "pvd", "arrhythmia"),
  score_weight = c(18, 11, 10, 0, 8, 8, 6, 6, 5, 5, 4, 4, 4, 3, 2, 2, 2, 2, 1, 1)
)

setdiff(score_mapping$condition_name, colnames(comorb_patient_level))

weights <- setNames(
  score_mapping$score_weight,
  score_mapping$condition_name
)

comorb_patient_level <- comorb_patient_level %>%
  mutate(
    MCS_score = rowSums(
      across(
        all_of(names(weights)),
        ~ .x * weights[cur_column()]
      )
    )
  )

write.csv2(comorb_patient_level,
           "C:\\path\\HF local data\\msc.csv",
           row.names = FALSE,
           fileEncoding = "UTF-8")

dataset <- data %>%
  filter(tipo_prest %in% c(41, 42, 43)) %>%
  mutate(
    data_studio_out = as.numeric(as.POSIXct(data_studio_out, tz = "UTC")),
  ) %>%
  select(COD_REG, sesso, eta, gruppo, class_prest, dec_intra, qt_prest_Sum,
         ICD, SHOCK, CABG, PTCA) %>%
  group_by(COD_REG) %>%
  arrange(data_prest) %>%
  mutate(
    # This assigns the SAME time_index (e.g., 0, 5, 17, 20) to events sharing a date
    data_prest = as.numeric(as.Date(data_prest) - first(as.Date(data_prest))), 
    time_index = as.numeric(data_prest - first(data_prest)) 
  ) %>%
  ungroup() %>%
  as.data.table()

dataset <- dataset %>%
  left_join(
    comorb_patient_level %>% select(COD_REG, MCS_score),
    by = "COD_REG"
  )

length(unique(dataset$COD_REG))
#187514

dataset %>%
  count(COD_REG, name = "n_eventi") %>%
  filter(n_eventi > 2)

# 88937 pazienti con più di 2 eventi 

dataset %>%
  distinct(COD_REG, desc_studio_out) %>%
  count(desc_studio_out)

#  desc_studio_out     n
#        DECEDUTO 87256
#           PERSO  1444
#        TRONCATO 98814

dataset %>%
  distinct(COD_REG, MCS_score) %>%
  count(MCS_score)

#sesso     n
#     F 95253
#     M 92261

dataset %>%
  distinct(COD_REG, MCS_score) %>%
  mutate(
    fascia_MCS = cut(
      MCS_score,
      breaks = c(-Inf, 4, 9, 14, 20, Inf),
      labels = c("0","1", "2", "3", "4"),
      right = TRUE
    )
  ) %>%
  count(fascia_MCS)
#  fascia_MCS     n
#          0 35957
#          1 78614
#          2 36668
#          3 20081
#          4 16194

new <- dataset %>%
  mutate(
    sesso = as.numeric(factor(sesso)),
    gruppo = as.numeric(factor(gruppo)),
    class_prest = as.numeric(factor(class_prest))
  )

write.csv2(new,
           "C:\\path\\HF local data\\dataset_ricoveri_encode_index.csv",
           row.names = FALSE,
           fileEncoding = "UTF-8")

