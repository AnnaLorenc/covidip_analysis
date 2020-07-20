## ----read data in-----------------------------------------------------------------------------------------------------------------

#Load in sample, patient and flow cell population scheme. Assumes the files in
#the /data subdir of the current directory. Most current versions of the files
#are on the website.
sample <-  read_csv(sample_datalocation)
patient <- read_csv(patient_datalocation)


#Join sample and patient information
pat_flow_sero_cyto_ifna <-
  left_join(sample, patient, by = c("patient_id"))

#Subset, reshape
fortests <-
  pat_flow_sero_cyto_ifna %>%
  select(
    grep(
      "Median|_back|freq$|_cyto|Ratio|sero$",
      colnames(pat_flow_sero_cyto_ifna),
      val = T
    ),
    age,
    sex,
    patient_id,
    Clinical_sample,
    starts_with("class_")
  ) %>%
  pivot_longer(
    names_to = "param",
    cols = grep(
      "Median|_back|freq$|_cyto|Ratio|sero$",
      colnames(pat_flow_sero_cyto_ifna),
      val = T
    ),
    values_to = "par"
  )

#Add grouping into controls/covid patients
fortests <- fortests %>%
  mutate(class_ds = case_when(
    class_dss %in% c("Low", "Moderate", "Severe") ~ "COVID",
    TRUE ~ class_dss
  ))
