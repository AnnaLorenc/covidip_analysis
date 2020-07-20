## ----variance---------------------------------------------------------------------------------------------------------------
###Identify parameters with different intra-individual variance in sick than in healthy (changing_with_without)
##Use only main covid-ip cohort vs controls
#Higher variance in sick: parameters chaning in the course of the disease. For these using weights in fiting models is counterproductive.
#Higher variance in healthy: focusing  of the immune system

for_variance <- pat_flow_sero_cyto_ifna %>%
  filter(class_dss != "LRTI_nonCovid")

for_variance <- for_variance %>%
  mutate(cohort = ifelse(
    class_dss %in% c("seroneg_control", "seropos_control"),
    "control",
    "COVID"
  ))


test_var <-
  compare_individual_variances (
    grep(
      "freq|_back|_cyto|Ratio|Median|sero$",
      colnames(for_variance),
      val = T
    ),
    for_variance,
    direction = "g"
  )

test_var <- test_var [sapply(test_var, nrow) > 0] %>%
  bind_rows


test_var_focusing <-
  compare_individual_variances (
    grep(
      "freq|_back|_cyto|Ratio|Median|sero",
      colnames(for_variance),
      val = T
    ),
    for_variance,
    direction = "l"
  )
test_var_focusing <-
  test_var_focusing [sapply(test_var_focusing, nrow) > 0] %>%
  bind_rows


#The same, but now using only seronegative conrtols. Low N.
test_var_sn <-
  compare_individual_variances (
    grep(
      "freq|_back|_cyto|Ratio|Median|sero",
      colnames(for_variance),
      val = T
    ),
    for_variance %>% filter(class_dss !=
                              "seropos_control"),
    direction = "g"
  )

test_var_sn <- test_var_sn [sapply(test_var_sn, nrow) > 0] %>%
  bind_rows
test_var_focusing_sn <-
  compare_individual_variances (
    grep(
      "freq|_back|_cyto|Ratio|Median|sero",
      colnames(for_variance),
      val = T
    ),
    for_variance %>% filter(class_dss !=
                              "seropos_control"),
    direction = "l"
  )
test_var_focusing_sn <-
  test_var_focusing_sn [sapply(test_var_focusing_sn, nrow) > 0] %>%
  bind_rows


changing_with_without <- full_join(
  test_var %>% filter(
    all < cutoff_variance_test |
      only3 < cutoff_variance_test | only2 < cutoff_variance_test
  ) %>%
    mutate(
      within_to_between = median_healthy_ind_sd / to_scale_by_within_healthy_sd,
      type = "sick_higher"
    ) %>%
    select(-c(sick_NA)),
  test_var_sn %>% filter(
    all < cutoff_variance_test |
      only3 < cutoff_variance_test | only2 < cutoff_variance_test
  ) %>%
    mutate(
      within_to_between = median_healthy_ind_sd / to_scale_by_within_healthy_sd,
      type = "sick_higher"
    ) %>%
    select(-c(sick_NA)),
  by = "par",
  suffix = c("_all", "_sn")
)


focusing_with_without <- full_join(
  test_var_focusing %>% filter(
    all < cutoff_variance_test |
      only3 < cutoff_variance_test | only2 < cutoff_variance_test
  ) %>%
    mutate(
      within_to_between = median_healthy_ind_sd / to_scale_by_within_healthy_sd,
      type = "focusing"
    ) %>%
    select(-c(sick_NA)),
  test_var_focusing_sn %>% filter(
    all < cutoff_variance_test |
      only3 < cutoff_variance_test | only2 < cutoff_variance_test
  ) %>%
    mutate(
      within_to_between = median_healthy_ind_sd / to_scale_by_within_healthy_sd,
      type = "focusing"
    ) %>%
    select(-c(sick_NA)),
  by = "par",
  suffix = c("_all", "_sn")
)
