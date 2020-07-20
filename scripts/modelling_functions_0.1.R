library(tidyverse)
library(ggplot2)
library(dplyr)
library(data.table)
library(readr)
library(ggbeeswarm)
library(lme4)
library(lmerTest)
library(effects)
library(broom)
library(broom.mixed)
library(rstatix)


run_a_model <- function(model_desired,
                        for_this_param,
                        scale_arg,
                        conf_level = .9,
                        weights = NULL,
                        what_signif = NULL,
                        column_with_test = "test",
                        model_specified_by_name = FALSE) {
  #model_desired is resu[[your_param]], with  cols:(column_with_test and column scale_arg) OR a name of a model
  
  if (!model_specified_by_name) {
    model_desired <- model_desired %>%
      filter(scale_arg == !!scale_arg) %>%
      select(column_with_test) %>%
      unlist(., use.names = FALSE)
    print(model_desired)
  }
  
  models_to_choose_from <- list(
    fit1_lmer =  paste0(scale_arg, "~ 1 + age+ (1|patient_id)"),
    fit2_lmer =  paste0(scale_arg, " ~ 1 + age+sex+ (1|patient_id)"),
    fit3_lmer =  paste0(scale_arg, " ~ 1 + sex+ (1|patient_id)"),
    fit1_lm = paste0(scale_arg, " ~ 1 + age"),
    fit2_lm = paste0(scale_arg, " ~ 1 + age + sex"),
    fit3_lm = paste0(scale_arg, " ~ 1 +  sex")
  )
  
  if (!is.null(weights)) {
    weights <- for_this_param %>% data.table() %>%
      .[, .(.N, par), by = patient_id] %>%
      .[, .(weights = 1 / N)] %>% unlist()
  }
  
  if (model_desired == "fit2_lmer") {
    model_desired = ifelse(
      what_signif == "age+sex",
      "fit2_lmer",
      ifelse(what_signif == "age", "fit1_lmer" , "fit3_lmer")
    )
  }
  if (model_desired == "fit2_lm") {
    model_desired = ifelse(
      what_signif == "age+sex",
      "fit2_lm",
      ifelse(what_signif == "age", "fit1_lm" , "fit3_lm")
    )
  }
  #prepare model# lme4 as intervals do not work with other
  if (grepl(pattern = "lmer", model_desired)) {
    this_param_model = substitute(expression(
      lme4:::lmer(
        formula(lmer_model_to_run),
        data = for_this_param,
        REML = TRUE,
        weights = weights
      )
    ),
    list(lmer_model_to_run = models_to_choose_from[[model_desired]]))
  } else{
    this_param_model = substitute(expression(lm(
      formula(lm_model_to_run),
      data = for_this_param,
      weights = weights
    )),
    list(lm_model_to_run = models_to_choose_from[[model_desired]]))
  }
  
  final_model = eval(eval(this_param_model))
  
  return(final_model)
}



prepare_pred_intervals <-
  function(model,
           newdata,
           scale_arg,
           add_to_real_data = TRUE) {
    #prepare prediction intervals dispatch
    if (class(model) %in% c("lmerMod", "lmerModLmerTest")) {
      prepare_pred_intervals_lmer(
        lmer_fit = model,
        newdata = newdata,
        scale_arg = scale_arg,
        add_to_real_data = add_to_real_data
      )
    } else{
      if ("lm" %in% class(model)) {
        prepare_pred_intervals_lm(lm_fit = model, newdata = newdata)
      }
    }
  }

prepare_pred_intervals_lm <- function(lm_fit, newdata) {
  #prepare prediction intervals when model used is lm
  newdat <-
    data.frame(newdata, predict(lm_fit, newdata, interval = "prediction")) %>%
    rename(c("plo" = "lwr", "phi" = "upr")) %>%
    left_join(.,
              newdata %>% select(Clinical_sample, patient_id, sex, age, eval(scale_arg)))
  return(newdat)
}

prepare_pred_intervals_lmer <-
  function(lmer_fit,
           newdata,
           scale_arg,
           add_to_real_data = TRUE) {
    #prepare prediction intervals when model used is lmer
    #add_to_real_data is about adding to the input columns so the table is wider
    if (add_to_real_data) {
      newdata[[paste0(scale_arg, "_real")]] <- newdata[[scale_arg]]

    }
    newdata <- newdata %>%
      filter(!is.na(sex), !is.na(age), !is.na(!!scale_arg)) %>%
      select(Clinical_sample,
             patient_id,
             sex,
             age,
             eval(scale_arg),
             eval(paste0(scale_arg, "_real")))
    
    newdata[[scale_arg]] <-
      lme4:::predict.merMod(lmer_fit, newdata, re.form = NA)
    
    mm <- model.matrix(terms(lmer_fit), model.frame(data = newdata))
    print(dim(newdata))
    print(dim(mm))
    pvar1 <- diag(mm %*% tcrossprod(vcov(lmer_fit), mm))
    cf = VarCorr(lmer_fit) %>% as.data.frame() %>% .[1, "vcov"]
    tvar1 <- pvar1 +  cf
    cmult <- 2 ## could use 1.96
    newdat <- data.frame(
      newdata,
      plo = newdata[[scale_arg]] - cmult * sqrt(pvar1),
      phi = newdata[[scale_arg]] + cmult * sqrt(pvar1),
      tlo = newdata[[scale_arg]] * sqrt(tvar1),
      thi = newdata[[scale_arg]] + cmult * sqrt(tvar1)
    ) %>% bind_cols
    
    newdat <- newdat %>%
      rename("fit" = !!scale_arg)
    
    if (add_to_real_data) {
      newdat <- newdat %>%
        rename(!!scale_arg := paste0(scale_arg, "_real"))
    }
    return(newdat)
  }


test_raw <- function(whole_data,
                     scale_arg,
                     weights = FALSE,
                     exclude_low = FALSE,
                     designations_to_test = c("sick", "class_WHO_severity", "sick_sevM")) {
  #Main testing function
  if (exclude_low) {
    whole_data <- whole_data %>% filter(class_WHO_severity != "Low")
  }
  
  if (sum(!is.na(whole_data[[scale_arg]])) > 3 &
      length(intersect(designations_to_test, colnames(whole_data))) > 0) {
    #this is to eliminate NAs in  logpar

    
    if (weights) {
      weights <- whole_data %>% data.table() %>%
        .[, .(.N, par), by = patient_id] %>%
        .[, .(weights = 1 / N)] %>% unlist()
      
 
    } else{
      weights = NULL
    }
    
    res = list()
    for (fac in intersect(designations_to_test, colnames(whole_data))) {
      a = NULL
      is_singular_a = FALSE
      if (length(whole_data$patient_id) > length(unique(whole_data$patient_id))) {
        if (levels(factor(whole_data[[fac]]))[1] %in% whole_data[[fac]]) {
          #checking whether the ref level is present in this subset
        
          a <-
            lmer(formula(paste0(
              scale_arg, "~ ", fac,  "+(1|patient_id)"
            )),
            data = whole_data,
            weights = weights)
          is_singular_a = isSingular(a)
        }
      }
      
      if (is.null(a) | is_singular_a) {
        if (levels(whole_data[[fac]])[1] %in% whole_data[[fac]])
          lm(formula(paste0(scale_arg,  "~ ", fac)),
             data = whole_data,
             weights = weights)
      }
      res[[fac]] <- a %>%
        tidy() %>%
        mutate(group = fac)
    }
    res = bind_rows(res)
  } else{
    res = NULL
  }
  return(res)
}




compare_individual_variances <-
  function(params, dataset, direction) {
    #Collect intraindividual variance estimate (as sd) from all individuals with 2 or 3 samples.
    #Compare distributions of these estimates in controls/patients
    alternative = direction
    sapply(params, function(par) {
      print(par)
      a = dataset %>% select(Clinical_sample, patient_id, cohort, par = eval(par)) %>%
        filter(!is.na(par)) %>%
        group_by(patient_id) %>%
        summarise(
          sd = sd(par, na.rm = T),
          cohort = unique(cohort),
          n = n()
        )
      
      
      result <- a %>%
        mutate(
          group = ifelse(cohort == "control", "control", "sick"),
          no_of_datapoints = factor(n, levels = c(1:3))
        ) %>%
        group_by(group, no_of_datapoints, .drop = FALSE) %>%
        summarise(n = n()) %>%
        filter(!is.na(group)) %>%
        pivot_wider(
          names_from = c(group, no_of_datapoints),
          values_from = n,
          values_fill = 0
        )
      
      result$all <- tryCatch({
        a %>% filter(!is.na(sd)) %>%
          mutate(group = cohort == "control") %>%
          wilcox.test(sd ~ group, data = ., alternative = alternative) %>%
          .$p.value
      },
      error = function(cond) {
        return(NA)
      })
      result$only3 <- tryCatch({
        a %>% filter(!is.na(sd), n == 3) %>%
          mutate(group = cohort == "control") %>%
          wilcox.test(sd ~ group, data = ., alternative = alternative) %>%
          .$p.value
      },
      error = function(cond) {
        return(NA)
      })
      
      result$only2 <- tryCatch({
        a %>% filter(!is.na(sd), n == 2) %>%
          mutate(group = cohort == "control") %>%
          wilcox.test(sd ~ group, data = ., alternative = alternative) %>%
          .$p.value
      },
      error = function(cond) {
        return(NA)
      })
      
      result$median_healthy_ind_sd <-  tryCatch(
        a %>%
          filter(cohort == "control") %>%
          summarise(sd = median(sd, na.rm = T)) %>% unlist(),
        error = function(cond) {
          return(NA)
        }
      )
      
      result$par = par
      result$to_scale_by_within_healthy_sd <-  tryCatch({
        dataset %>% select(Clinical_sample, patient_id, cohort, par = eval(par)) %>%
          filter(!is.na(par), cohort == "control") %>%
          .$par %>% sd(., na.rm = T)
      },
      error = function(cond) {
        return(NA)
      })
      return(result)
    })
    
  }



prepare_for_testing <- function(whole_data, perform_log = TRUE) {
  #Add columns with info for contrasts of interest, perform log etc.
  whole_data <- whole_data %>%
    filter(!is.na(par), !is.infinite(par)) %>%
    mutate(
      sick = ifelse(grepl("control", class_ds), "control", class_ds) %>% factor(),
      sick_LRTI = factor(sick, levels = c("LRTI_nonCovid", "COVID", "control")),
      #for control vs covid, control vs LRTI, LRTI vs covid
      sick_CLMS = ifelse(grepl("control", class_ds), "control", class_dss) %>%
        factor(),
      sick_SLMS = factor(
        sick_CLMS,
        levels = c("LRTI_nonCovid", "control", "Low", "Moderate", "Severe")
      ),
      sick_sevM = factor(
        gsub(pat = ".*_", rep = "", class_dss),
        levels = c("Moderate", "Low", "Severe", "control", "nonCovid")
      ),
      sick_sevL = factor(
        gsub(pat = ".*_", rep = "", class_dss),
        levels = c("Low", "Severe", "Moderate", "control", "nonCovid")
      )
    )
  if (perform_log) {
    minpar <- min(whole_data$par[whole_data$par > 0], na.rm = T)
    # print(minpar)
    
    whole_data <- whole_data %>%
      mutate(logpar = case_when(
        par == 0 ~ log10(0.001 * minpar),
        par < 0 ~ NA_real_,
        TRUE ~ log10(par)
      ))
  }
  return(whole_data)
}

extract_controls  <- function(whole_data, perform_log = FALSE) {
  #take control subset of the data
  whole_data <- whole_data %>%
    filter(!is.na(par), !is.infinite(par)) %>%
    filter(grepl("control", class_ds))
  
  if (perform_log) {
    minpar <- min(whole_data$par[whole_data$par > 0], na.rm = T)
    # print(minpar)
    
    whole_data <- whole_data %>%
      mutate(par = case_when(
        par == 0 ~ log10(0.001 * minpar),
        par < 0 ~ NA_real_,
        TRUE ~ log10(par)
      ))
  }
  return(whole_data)
  
}


test_ser_controls <- function(this_param_data) {
  lm(par ~ class_ds, data = this_param_data) %>%
#    summary() %>%
    tidy()
  
}

check_which_params_might_be_tested <- function(test_ready,
                                               designations_to_test = c("sick",
                                                                        "sick_LRTI",
                                                                        "sick_CLMS",
                                                                        "sick_SLMS",
                                                                        "sick_sevM" ,
                                                                        "sick_sevL")) {
  designations_to_test <-
    designations_to_test[sapply(designations_to_test, function(i)
      length(levels(test_ready[[i]])) <= length(unique(test_ready[[i]])))]
  if (length(unique(test_ready$par)) < 2) {
    designations_to_test <- c()
  }
  return(designations_to_test)
}

select_columns <-
  function(full_results_param,
           columns_pivotal,
           columns_all,
           cutoff_pval_snonly) {
    #Based on test results (sex/age dependence, parameter type, seropos vs seroneg), decide which pval/effsize is the one to use
    
    if (grepl("_cyto|_LIPS|Count_back", full_results_param$param)) {
      columns_pivotal <- grep("log", columns_pivotal, val = T)
      columns_all <- grep("log",  columns_all , val = T)
    } else{
      columns_pivotal <-
        grep("log", columns_pivotal, val = T, invert = TRUE)
      columns_all <- grep("log",  columns_all , val = T, invert = TRUE)
    }
    
    
    is_sexage_correction_important_col <-
      grep(pattern = "sexage",
           x = columns_pivotal,
           value = TRUE)
    should_snonly_beused_col <-
      grep("pvalseropos", columns_pivotal, value = TRUE)
    

    if (full_results_param[[is_sexage_correction_important_col]] != "none") {
      columns_all <-
        c(
          c("param", "about_weights", "pval_seropos_final"),
          grep("corrected", columns_all, value = TRUE)
        )
    } else{
      columns_all <-
        grep("corrected",
             columns_all,
             value = TRUE,
             invert = TRUE)
    }
    
    if (full_results_param[[should_snonly_beused_col]] < cutoff_pval_snonly |
        is.na(full_results_param[[should_snonly_beused_col]])) {
      columns_all <-
        grep("_all", columns_all, value = TRUE, inv = TRUE) #take sn columns only
    } else{
      columns_all <-
        grep("_sn", columns_all, value = TRUE, invert = TRUE)
    }  #take all columns only
    
    return(columns_all)
    
  }
