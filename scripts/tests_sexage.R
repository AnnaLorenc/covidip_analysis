## ----sex_age_correct--------------------------------------------------------------------------------------------------------
#Prepare parameters which need sex/age correction: correct all samples to be female/have 35yo
#This part requires age specification, which we cannot release due to anonymity reasons.


residuals_list <- list()
corrected_list <- list()

for (your_param in (resu_rb %>% filter(what_signif_0.01 != "none") %>% .$par %>%
                    unique())) {
  #  your_param="NK_01/Non_BT_01 |freq"
  print(your_param)
  #prepare whole_data
  whole_data <- pat_flow_sero_cyto_ifna %>%
    filter(!is.na(get(your_param))) %>%
    select(age, sex, patient_id, par = your_param, Clinical_sample, class_dss) %>%
    as_tibble()
  
  
  whole_data <- whole_data %>%
    mutate(par = ifelse(is.infinite(par), NA_real_, par)) %>%
    filter(!is.na(par))
  
  #ad hoc fix for a new parameter which is identical across samples
  if (nrow(whole_data) < 30 |
      your_param == "Time_06/Cells_06/Singlets1_06/Singlets2_06/CD45p_06/Not_lymphocytes_06/CD15p_06/Total_neutrophils_06/Mature_neutrophil_06/Mature_Immunosuppressive_Neutrophil_06 | Median (CD62L)") {
    warnings_weird_param[[your_param]] <- your_param
    next
  }
  
  
  whole_data <-   whole_data %>%
    mutate(sick = ifelse(
      class_dss == "LRTI_nonCovid",
      NA_character_,
      ifelse(
        class_dss %in% c("seropos_control", "seroneg_control"),
        "healthy",
        "sick"
      )
    ))
  
  minpar <- whole_data %>%  #changing to evade 0 issues in log10
    filter(par != 0) %>% .$par %>% min(., na.rm = T)
  
  whole_data <- whole_data %>%
    mutate(logpar = ifelse(par == 0, log10(0.001 * minpar), log10(par)))
  
  print(dim(whole_data))
  if (length(unique(whole_data$par)) < 3) {
    print(whole_data$par)
    next
  }
  
  model_desired <- resu_rb %>%
    filter(par == your_param) %>%
    unique()
  
  residuals_list[[your_param]] <- list()
  corrected_list[[your_param]] <- list()
  
  for (scale_arg in model_desired$scale_arg) {
    print(scale_arg)
    if (your_param %in% (
      resu_rb %>% dplyr:::filter(scale_arg == !!scale_arg, what_signif_0.01 !=
                                 "none") %>% .$par
    )) {
      print("correcting...")
      
      for_this_param <- whole_data %>%
        filter(class_status_detailed == "control")
      
      
      
      final_model_full_control <-
        run_a_model(
          model_desired,
          for_this_param,
          scale_arg = scale_arg,
          weights = TRUE,
          what_signif = model_desired %>% filter(scale_arg == !!scale_arg) %>% .$what_signif_0.01
        )
      #add to the list of residuals
      which_factors_significant <- final_model_full_control %>%
        terms(.) %>% labels(.)
      
      print(which_factors_significant)
      
      whole_data <- whole_data %>%
        filter_at(vars(sex, age), ~ !is.na(.))
      
      tryCatch({
        a <-
          prepare_pred_intervals(
            final_model_full_control,
            whole_data,
            scale_arg = scale_arg,
            add_to_real_data = TRUE
          )
        
        a <- a %>%
          mutate(resids = get(scale_arg) - fit)
        intercept <- final_model_full_control %>%
          summary() %>%
          coefficients() %>%
          .["(Intercept)", "Estimate"]
        
        if (grepl("age", which_factors_significant)) {
          slope <- final_model_full_control %>%
            summary() %>%
            coefficients() %>%
            .["age", "Estimate"]
        } else{
          slope = 0
        }
        
        print(intercept)
        residuals_list[[your_param]][[scale_arg]] <- a
        corrected_list[[your_param]][[scale_arg]] <-
          transmute(a, Clinical_sample,
                    cor = resids + intercept + 35 *
                      slope) #making everyone 35 yold
        colnames(corrected_list[[your_param]][[scale_arg]]) <-
          c("Clinical_sample", paste0(eval(your_param), "_", scale_arg))
        
        
      },
      
      error = function(cond) {
        message("impossible to correct")
        message("Here's the original error message:")
        message(cond)
        
        # Choose a return value in case of error
        return(NA)
      })
    } else{
      print("something wrong..")
    }
  }
}




corrected_list[lengths(corrected_list) == 0]

corrected <- corrected_list[lengths(corrected_list) > 0] %>%
  lapply(., function(x) {
    print(names(x))
    if (length(names(x)) > 1) {
      resu = full_join(x$logpar, x$par)
    } else{
      resu = x[[1]]
    }
    return(resu)
  }) %>%
  reduce(.x = ., .f = full_join)



## ----sex_age_corrected_data_tests-------------------------------------------------------------------------------------------
#Tests on data corrected for sex/age
#This part requires age specification, which we cannot release due to anonymity reasons.

corrected_for_tests <- full_join(
  pat_flow_sero_cyto_ifna %>%
    select(age, sex, patient_id, Clinical_sample, starts_with("class_")),
  corrected
) %>%
  pivot_longer(
    names_to = "param",
    cols = grep(
      "Median|_back|freq|_cyto|Ratio|sero",
      colnames(corrected),
      val = T
    ),
    values_to = "par"
  )

ct <- corrected_for_tests %>%
  group_by(param) %>%
  nest() %>%
  mutate(
    test_ready = map(data, prepare_for_testing, perform_log = FALSE),
    designations_to_test = map(test_ready, check_which_params_might_be_tested),
    about_weights = !(
      gsub(pattern = "_par|_logpar", rep = "", param) %in% changing_with_without$par
    ),
    about_difference_in_seropos = param %in% diff_within_controls
  )

ct <- ct %>%
  mutate(
    test_ready_sn = map(test_ready, function(x)
      filter(x, class_ds != "seropos_control", .preserve = FALSE)),
    designations_to_test_sn = map(test_ready_sn, check_which_params_might_be_tested)
  )


ct1 <- ct %>%
  mutate(test_par = pmap(
    list(
      whole_data = test_ready,
      scale_arg = "par",
      weights = about_weights,
      exclude_low = FALSE,
      designations_to_test = designations_to_test
    )  ,
    test_raw
  ),
  test_par_sn = pmap(
    list(
      whole_data = test_ready_sn,
      scale_arg = "par",
      weights = about_weights,
      exclude_low = FALSE,
      designations_to_test = designations_to_test_sn
    )  ,
    test_raw
  ))

ct1 <- ct1 %>%
  mutate(
    clean_param = gsub(param, pattern = "_par|_logpar", rep = ""),
    scale_arg = gsub(param, pattern = ".*_", rep = "")
  )

ct2 <- pivot_wider(
  ct1,
  id_cols = param,
  names_from = scale_arg,
  values_from = c(test_par, test_par_sn)
)



