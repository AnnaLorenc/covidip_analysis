library(tidyverse)
library(ggplot2)
library(dplyr)
library(data.table)
library(readr)
library(ggbeeswarm)
library(scales)
library(ggpubr)
library(ggsignif)

sample <-  read_csv(sample_datalocation)
patient <- read_csv(patient_datalocation)


combined_data <-
  left_join(sample, patient, by = c("patient_id"))

combined_data <- left_join(megatable, patient_data, by = "patient_id")

combined_data_1 <-   combined_data %>%
  mutate(
    admission_WHO_severity = case_when(
      class_dss %in% c("Low", "Moderate", "Severe") &
        admission_WHO_ordinal_scale %in% c(0, 1, 2) ~ "Low",
      class_dss %in% c("Low", "Moderate", "Severe") &
        admission_WHO_ordinal_scale %in% c(5, 6, 7, 8) ~ "Severe",
      class_dss %in% c("Low", "Moderate", "Severe") &
        admission_WHO_ordinal_scale %in% c(3, 4) ~ "Moderate",
      class_dss %in% c("seroneg_control", "seropos_control") ~ "Healthy"
    )
  ) %>%
  mutate(class_dshi = ifelse(patient_id %in% c("p028", "p049", "p053", "p071"), "HI", "COVID")) %>%
  mutate(
    trajectory_week = case_when(
      `Delta_WHO_COVIDBleed1_to_7days_or discharge_if_sooner` > 0 ~ "Worsen",
      `Delta_WHO_COVIDBleed1_to_7days_or discharge_if_sooner` <
        0 ~ "Improve",
      TRUE ~ "Stable"
    )
  ) %>%
  mutate(
    CRP_since = parse_number(covid_ip_bleed_01) - parse_number(date_COVIDIPbleed_CRP),
    albumin_since = parse_number(covid_ip_bleed_01) - parse_number(date_COVIDIPbleed_albumin),
    ferritin_since = parse_number(covid_ip_bleed_01) - parse_number(date_COVIDIPbleed_ferritin),
    FBC_lymph_since = parse_number(covid_ip_bleed_01) - parse_number(date_COVIDIPbleed_lymph),
    FBC_mono_since = parse_number(covid_ip_bleed_01) - parse_number(date_COVIDIPbleed_mono),
    FBC_neut_since = parse_number(covid_ip_bleed_01) - parse_number(date_COVIDIPbleed_neut),
    Ddimer_since = parse_number(covid_ip_bleed_01) - parse_number(date_COVIDIPbleed_Ddimer),
    ALT_since = parse_number(covid_ip_bleed_01) - parse_number(date_COVIDIPbleed_ALT)
  )

#ensure clinical parameters are only considered if within 3 days of 1st COVID bleed
combined_data_2 <- combined_data_1 %>%
  mutate(value_COVIDIPbleed_CRP = (ifelse((CRP_since <= 3 &
                                             CRP_since >= -3), value_COVIDIPbleed_CRP, NA
  ))) %>%
  mutate(value_COVIDIPbleed_albumin = (ifelse((albumin_since <= 3 &
                                                 albumin_since >= -3),
                                              value_COVIDIPbleed_albumin,
                                              NA
  ))) %>%
  mutate(value_COVIDIPbleed_ferritin = (ifelse((ferritin_since <= 3 &
                                                  ferritin_since >= -3),
                                               value_COVIDIPbleed_ferritin,
                                               NA
  ))) %>%
  mutate(value_COVIDIPbleed_lymph = (ifelse((FBC_lymph_since <= 3 &
                                               FBC_lymph_since >= -3),
                                            value_COVIDIPbleed_lymph,
                                            NA
  ))) %>%
  mutate(value_COVIDIPbleed_mono = (ifelse((FBC_mono_since <= 3 &
                                              FBC_mono_since >= -3),
                                           value_COVIDIPbleed_mono,
                                           NA
  ))) %>%
  mutate(value_COVIDIPbleed_neut = (ifelse((FBC_neut_since <= 3 &
                                              FBC_neut_since >= -3),
                                           value_COVIDIPbleed_neut,
                                           NA
  ))) %>%
  mutate(value_COVIDIPbleed_ALT = (ifelse((ALT_since <= 3 &
                                             ALT_since >= -3), value_COVIDIPbleed_ALT, NA
  )))
#for_models_2 <- subset(for_models_1,duplicated(patient_id) | duplicated(patient_id, fromLast=TRUE))

prognostic_param <- as.vector(
  c(
    "value_COVIDIPbleed_ferritin",
    "value_COVIDIPbleed_CRP",
    "value_COVIDIPbleed_lymph",
    "value_COVIDIPbleed_Ddimer",
    "value_COVIDIPbleed_ALT",
    "value_COVIDIPbleed_albumin",
    "Time_06/Cells_06/Singlets1_06/Singlets2_06/CD45p_06/Lymphocytes_06/CD3p_CD3n_06/T_06 | Count_back",
    "IP10_av_cyto_cyto",
    "IL6_av_cyto_cyto",
    "IL10_av_cyto_cyto",
    "Spike_IgM_Norm_sero",
    "RBD_IgM_Norm_sero",
    "Spike_IgG_Norm_sero",
    "RBD_IgG_Norm_sero",
    "peak_WHO_ordinal_scale"
  )
)

prognostic_readouts <-
  as.vector(c("length_stay_hospital_post_Covidbleed01",
              "trajectory_week"))

#choose first visit only, COVID only, only admitted as moderate or severe
combined_data_3 <- combined_data_2 %>%
  mutate(visit_number = str_sub(Clinical_sample, -1)) %>%
  filter(visit_number == 1) %>%
  filter((class_dss %in% c("Low", "Moderate", "Severe"))) %>%
  filter(admission_WHO_severity %in% c("Moderate", "Severe"))


#####Kruskal Wallis test for parameters of interest in Fig 6a

for_stats <-
  combined_data_3 %>% select(
    patient_id,
    trajectory_week,
    value_COVIDIPbleed_CRP,
    IL6_av_cyto_cyto,
    IL10_av_cyto_cyto,
    IP10_av_cyto_cyto
  ) %>%
  pivot_longer(value_COVIDIPbleed_CRP:IP10_av_cyto_cyto) %>%
  group_by(name)

#report Dunn multiple comparison test results, with and without p-value adjustment
tibble <- for_stats %>%
  rstatix::dunn_test(., value ~ trajectory_week, p.adjust.method = "bonferroni") %>%   #you can end here or reformat if the other version suits you better
  mutate(comparison = paste(group1, group2, sep = "-")) %>%
  select(name, p:comparison) %>%
  pivot_wider(
    id_cols = name,
    names_from = comparison,
    values_from = p:p.adj.signif
  )

write.csv(tibble, "trajectory_stats_20200710.csv")
