

#Prepare data: gather per parameter info:
#-for which contrasts there is enough data,
#-whether weights are to be used (based on intraindividual variance),
#-whether all controls might be used (based on seonegative-seropositive differences,
#-when difference among controls with pva<0.05 or impossible to compute)
#All variables are tested on normal and log10 scale

whole <- fortests %>%
  group_by(param) %>%
  nest() %>%
  mutate(
    test_ready = map(data, prepare_for_testing),
    designations_to_test = map(test_ready, check_which_params_might_be_tested),
    about_weights = !(param %in% changing_with_without$par),
    about_difference_in_seropos = param %in% diff_within_controls
  )


#whether sex/age correction is needed - this part relies on running model fit/access to its results
whole <- whole %>%
  mutate(
    par_sexage = !!resu_rb %>%
      filter(scale_arg == "par", par == param) %>% .$what_signif_0.01 %>%
      ifelse(length(.) == 0, "none", .),
    logpar_sexage = !!resu_rb %>%
      filter(scale_arg == "logpar", par == param) %>% .$what_signif_0.01 %>%
      ifelse(length(.) == 0, "none", .)
  )



#Test: use all controls
whole1 <- whole %>%
  mutate(test_par = pmap(
    list(
      whole_data = test_ready,
      scale_arg = "par",
      weights = about_weights,
      exclude_low = FALSE,
      designations_to_test = designations_to_test
    ),
    test_raw
  ))

whole2 <- whole1 %>%
  mutate(test_logpar = pmap(
    list(
      whole_data = test_ready,
      scale_arg = "logpar",
      weights = about_weights,
      exclude_low = FALSE,
      designations_to_test = designations_to_test
    ),
    test_raw
  ))

#Test: use seroneg controls only
whole3 <- whole2 %>%
  mutate(
    test_ready_sn = map(test_ready, function(x)
      filter(x, class_dss != "seropos_control", .preserve = FALSE)),
    designations_to_test_sn = map(test_ready_sn, check_which_params_might_be_tested)
  )

whole4 <- whole3 %>%
  mutate(test_par_sn = pmap(
    list(
      whole_data = test_ready_sn,
      scale_arg = "par",
      weights = about_weights,
      exclude_low = FALSE,
      designations_to_test = designations_to_test_sn
    ),
    test_raw
  ))

whole5 <- whole4 %>%
  mutate(test_logpar_sn = pmap(
    list(
      whole_data = test_ready_sn,
      scale_arg = "logpar",
      weights = about_weights,
      exclude_low = FALSE,
      designations_to_test = designations_to_test
    ),
    test_raw
  ))
