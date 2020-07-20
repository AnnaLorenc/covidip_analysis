# Add results for parameters corrected for age/sex

results_corrected_log_all <-
  ct2 %>%
  filter(!is.null(test_par_logpar)) %>%
  select(test_par_logpar) %>%
  unnest(cols = c(test_par_logpar)) %>%
  filter(effect != "ran_pars") %>%
  select(param, group , term, estimate , std.error, p.value) %>%
  pivot_wider(
    id_cols = param,
    names_from = c(group, term),
    values_from = estimate:p.value
  ) %>%
  transmute(
    param,
    healthy_all_mean = `estimate_sick_(Intercept)`,
    LRTI_mean = `estimate_sick_LRTI_(Intercept)`,
    COVID_mean = healthy_all_mean + estimate_sick_sickCOVID,
    COVIDlow_mean = ifelse(
      exists("estimate_sick_sevL_(Intercept)"),
      `estimate_sick_sevL_(Intercept)`,
      NA_real_
    ),
    COVIDmoderate_mean = `estimate_sick_sevM_(Intercept)`,
    COVIDsevere_mean = `estimate_sick_CLMS_(Intercept)` + estimate_sick_CLMS_sick_CLMSSevere,
    healthy_all_vs_COVID.effsize = estimate_sick_sickCOVID / `std.error_sick_(Intercept)`,
    healthy_all_vs_COVID.pval = p.value_sick_sickCOVID,
    healthy_all_vs_LRTI.pval = p.value_sick_sickLRTI_nonCovid,
    LRTI_vs_COVID.pval = p.value_sick_LRTI_sick_LRTICOVID,
    COVIDmoderate_vs_healthy_all.pval = p.value_sick_CLMS_sick_CLMSModerate,
    COVIDlow_vs_healthy_all.pval = ifelse(
      exists("p.value_sick_CLMS_sick_CLMSLow"),
      p.value_sick_CLMS_sick_CLMSLow ,
      NA_real_
    ),
    COVIDsevere_vs_healthy_all.pval = p.value_sick_CLMS_sick_CLMSSevere,
    COVIDmoderate_vs_COVIDlow.pval = ifelse(
      exists("p.value_sick_sevL_sick_sevLModerate"),
      p.value_sick_sevL_sick_sevLModerate ,
      NA_real_
    ),
    COVIDsevere_vs_COVIDmoderate.pval = p.value_sick_sevM_sick_sevMSevere,
    COVIDsevere_vs_COVIDlow.pval = ifelse(
      exists("p.value_sick_sevL_sick_sevLSevere"),
      p.value_sick_sevL_sick_sevLSevere,
      NA_real_
    ),
    COVIDlow_vs_LRTI.pval = ifelse(
      exists("p.value_sick_sevL_sick_sevLnonCovid"),
      p.value_sick_sevL_sick_sevLnonCovid,
      NA_real_
    ),
    COVIDmoderate_vs_LRTI.pval = ifelse(
      exists("p.value_sick_sevM_sick_sevMnonCovid"),
      p.value_sick_sevM_sick_sevMnonCovid,
      NA_real_
    ),
    COVIDsevere_vs_LRTI.pval = p.value_sick_SLMS_sick_SLMSSevere
  ) %>%
  rename_with(.cols = !starts_with("param"),  ~ paste0(., "_all_log_corrected")) %>%
  mutate(param = gsub(param, pattern = "_logpar", rep = ""))

results_corrected_log_sn <-
  ct2 %>%
  filter(!is.null(test_par_logpar)) %>%
  select(test_par_logpar) %>%
  unnest(cols = c(test_par_logpar)) %>%
  filter(effect != "ran_pars") %>%
  select(param, group , term, estimate , std.error, p.value) %>%
  pivot_wider(
    id_cols = param,
    names_from = c(group, term),
    values_from = estimate:p.value
  ) %>%
  transmute(
    param,
    healthy_sn_mean = `estimate_sick_(Intercept)`,
    LRTI_mean = `estimate_sick_LRTI_(Intercept)`,
    COVID_mean = healthy_sn_mean + estimate_sick_sickCOVID,
    COVIDlow_mean = ifelse(
      exists("estimate_sick_sevL_(Intercept)"),
      `estimate_sick_sevL_(Intercept)`,
      NA_real_
    ),
    COVIDmoderate_mean = `estimate_sick_sevM_(Intercept)`,
    COVIDsevere_mean = `estimate_sick_CLMS_(Intercept)` + estimate_sick_CLMS_sick_CLMSSevere,
    healthy_sn_vs_COVID.effsize = estimate_sick_sickCOVID / `std.error_sick_(Intercept)`,
    healthy_sn_vs_COVID.pval = p.value_sick_sickCOVID,
    healthy_sn_vs_LRTI.pval = p.value_sick_sickLRTI_nonCovid,
    LRTI_vs_COVID.pval = p.value_sick_LRTI_sick_LRTICOVID,
    COVIDmoderate_vs_healthy_sn.pval = p.value_sick_CLMS_sick_CLMSModerate,
    COVIDlow_vs_healthy_sn.pval = ifelse(
      exists("p.value_sick_CLMS_sick_CLMSLow"),
      p.value_sick_CLMS_sick_CLMSLow ,
      NA_real_
    ),
    COVIDsevere_vs_healthy_sn.pval = p.value_sick_CLMS_sick_CLMSSevere,
    COVIDmoderate_vs_COVIDlow.pval = ifelse(
      exists("p.value_sick_sevL_sick_sevLModerate"),
      p.value_sick_sevL_sick_sevLModerate ,
      NA_real_
    ),
    COVIDsevere_vs_COVIDmoderate.pval = p.value_sick_sevM_sick_sevMSevere,
    COVIDsevere_vs_COVIDlow.pval = ifelse(
      exists("p.value_sick_sevL_sick_sevLSevere"),
      p.value_sick_sevL_sick_sevLSevere,
      NA_real_
    ),
    COVIDlow_vs_LRTI.pval = ifelse(
      exists("p.value_sick_sevL_sick_sevLnonCovid"),
      p.value_sick_sevL_sick_sevLnonCovid,
      NA_real_
    ),
    COVIDmoderate_vs_LRTI.pval = ifelse(
      exists("p.value_sick_sevM_sick_sevMnonCovid"),
      p.value_sick_sevM_sick_sevMnonCovid,
      NA_real_
    ),
    COVIDsevere_vs_LRTI.pval = p.value_sick_SLMS_sick_SLMSSevere
  ) %>%
  rename_with(.cols = !starts_with("param"),  ~ paste0(., "_sn_log_corrected")) %>%
  mutate(param = gsub(param, pattern = "_logpar", rep = ""))



results_corrected_all <-
  ct2 %>%
  filter(!is.null(test_par_par)) %>%
  select(test_par_par) %>%
  unnest(cols = c(test_par_par)) %>%
  filter(effect != "ran_pars") %>%
  select(param, group , term, estimate , std.error, p.value) %>%
  pivot_wider(
    id_cols = param,
    names_from = c(group, term),
    values_from = estimate:p.value
  ) %>%
  transmute(
    param,
    healthy_all_mean = `estimate_sick_(Intercept)`,
    LRTI_mean = `estimate_sick_LRTI_(Intercept)`,
    COVID_mean = healthy_all_mean + estimate_sick_sickCOVID,
    COVIDlow_mean = ifelse(
      exists("estimate_sick_sevL_(Intercept)"),
      `estimate_sick_sevL_(Intercept)`,
      NA_real_
    ),
    COVIDmoderate_mean = `estimate_sick_sevM_(Intercept)`,
    COVIDsevere_mean = `estimate_sick_CLMS_(Intercept)` + estimate_sick_CLMS_sick_CLMSSevere,
    healthy_all_vs_COVID.effsize = estimate_sick_sickCOVID / `std.error_sick_(Intercept)`,
    healthy_all_vs_COVID.pval = p.value_sick_sickCOVID,
    healthy_all_vs_LRTI.pval = p.value_sick_sickLRTI_nonCovid,
    LRTI_vs_COVID.pval = p.value_sick_LRTI_sick_LRTICOVID,
    COVIDmoderate_vs_healthy_all.pval = p.value_sick_CLMS_sick_CLMSModerate,
    COVIDlow_vs_healthy_all.pval = ifelse(
      exists("p.value_sick_CLMS_sick_CLMSLow"),
      p.value_sick_CLMS_sick_CLMSLow ,
      NA_real_
    ),
    COVIDsevere_vs_healthy_all.pval = p.value_sick_CLMS_sick_CLMSSevere,
    COVIDmoderate_vs_COVIDlow.pval = ifelse(
      exists("p.value_sick_sevL_sick_sevLModerate"),
      p.value_sick_sevL_sick_sevLModerate ,
      NA_real_
    ),
    COVIDsevere_vs_COVIDmoderate.pval = p.value_sick_sevM_sick_sevMSevere,
    COVIDsevere_vs_COVIDlow.pval = ifelse(
      exists("p.value_sick_sevL_sick_sevLSevere"),
      p.value_sick_sevL_sick_sevLSevere,
      NA_real_
    ),
    COVIDlow_vs_LRTI.pval = ifelse(
      exists("p.value_sick_sevL_sick_sevLnonCovid"),
      p.value_sick_sevL_sick_sevLnonCovid,
      NA_real_
    ),
    COVIDmoderate_vs_LRTI.pval = ifelse(
      exists("p.value_sick_sevM_sick_sevMnonCovid"),
      p.value_sick_sevM_sick_sevMnonCovid,
      NA_real_
    ),
    COVIDsevere_vs_LRTI.pval = p.value_sick_SLMS_sick_SLMSSevere
  ) %>%
  rename_with(.cols = !starts_with("param"),  ~ paste0(., "_all_corrected")) %>%
  mutate(param = gsub(param, pattern = "_par", rep = ""))


results_corrected_sn <-
  ct2 %>%
  filter(!is.null(test_par_par)) %>%
  select(test_par_par) %>%
  unnest(cols = c(test_par_par)) %>%
  filter(effect != "ran_pars") %>%
  select(param, group , term, estimate , std.error, p.value) %>%
  pivot_wider(
    id_cols = param,
    names_from = c(group, term),
    values_from = estimate:p.value
  ) %>%
  transmute(
    param,
    healthy_sn_mean = `estimate_sick_(Intercept)`,
    LRTI_mean = `estimate_sick_LRTI_(Intercept)`,
    COVID_mean = healthy_sn_mean + estimate_sick_sickCOVID,
    COVIDlow_mean = ifelse(
      exists("estimate_sick_sevL_(Intercept)"),
      `estimate_sick_sevL_(Intercept)`,
      NA_real_
    ),
    COVIDmoderate_mean = `estimate_sick_sevM_(Intercept)`,
    COVIDsevere_mean = `estimate_sick_CLMS_(Intercept)` + estimate_sick_CLMS_sick_CLMSSevere,
    healthy_sn_vs_COVID.effsize = estimate_sick_sickCOVID / `std.error_sick_(Intercept)`,
    healthy_sn_vs_COVID.pval = p.value_sick_sickCOVID,
    healthy_sn_vs_LRTI.pval = p.value_sick_sickLRTI_nonCovid,
    LRTI_vs_COVID.pval = p.value_sick_LRTI_sick_LRTICOVID,
    COVIDmoderate_vs_healthy_sn.pval = p.value_sick_CLMS_sick_CLMSModerate,
    COVIDlow_vs_healthy_sn.pval = ifelse(
      exists("p.value_sick_CLMS_sick_CLMSLow"),
      p.value_sick_CLMS_sick_CLMSLow ,
      NA_real_
    ),
    COVIDsevere_vs_healthy_sn.pval = p.value_sick_CLMS_sick_CLMSSevere,
    COVIDmoderate_vs_COVIDlow.pval = ifelse(
      exists("p.value_sick_sevL_sick_sevLModerate"),
      p.value_sick_sevL_sick_sevLModerate ,
      NA_real_
    ),
    COVIDsevere_vs_COVIDmoderate.pval = p.value_sick_sevM_sick_sevMSevere,
    COVIDsevere_vs_COVIDlow.pval = ifelse(
      exists("p.value_sick_sevL_sick_sevLSevere"),
      p.value_sick_sevL_sick_sevLSevere,
      NA_real_
    ),
    COVIDlow_vs_LRTI.pval = ifelse(
      exists("p.value_sick_sevL_sick_sevLnonCovid"),
      p.value_sick_sevL_sick_sevLnonCovid,
      NA_real_
    ),
    COVIDmoderate_vs_LRTI.pval = ifelse(
      exists("p.value_sick_sevM_sick_sevMnonCovid"),
      p.value_sick_sevM_sick_sevMnonCovid,
      NA_real_
    ),
    COVIDsevere_vs_LRTI.pval = p.value_sick_SLMS_sick_SLMSSevere
  ) %>%
  rename_with(.cols = !starts_with("param"),  ~ paste0(., "_sn_corrected")) %>%
  mutate(param = gsub(param, pattern = "_par", rep = ""))


results_corrected <-
  list(
    results_corrected_all,
    results_corrected_sn,
    results_corrected_log_all,
    results_corrected_log_sn
  ) %>%
  reduce(., .f = full_join)
