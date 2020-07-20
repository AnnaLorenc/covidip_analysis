## ----test_sexage------------------------------------------------------------------------------------------------------------
###Using controls only, test for significant influence of age/sex.
#This part requires age specification, which we cannot release due to anonymity reasons.
#Models fitted are: mixed model via lmer, with patient as random effect and lm, when there is not enough data points to estimate patient effects.
#As patients differ in the number of datapoints they have, they are weighted by the number of data points.



for_sex_age_models <- pat_flow_sero_cyto_ifna %>%
  filter(class_dss %in% c("seroneg_control", "seropos_control"))


resu <- list()

for (your_param in c(grep(
  "freq|_cyto|_sero$|Ratio|_back",
  colnames(for_sex_age_models),
  val = T
))) {
  relaxed_significance_cutoff <- 0.05 #to perform rough filtering
  print(your_param)
  

  for_this_param <- for_sex_age_models %>%
    data.table() %>%
    .[!is.na(get(your_param)) &
        !is.na(age) &
        !is.na(sex), .(age, sex, patient_id, par = get(your_param))]
  
  if (length(unique(for_this_param$sex)) < 2) {
    next
  }
  
  minpar = for_this_param %>%  #evade 0 issues in log10
    filter(par != 0) %>% .$par %>% min(., na.rm = T)
  for_this_param <- for_this_param %>%
    mutate(logpar = ifelse(par == 0, log10(0.001 * minpar), log10(par)))
  
  weights <- for_this_param %>% data.table() %>%
    .[, .(.N, par), by = patient_id] %>%
    .[, .(weights = 1 / N)] %>% unlist()
  
  if (length(unique(for_this_param$par)) > 3) {
    results = list()
    for (scale_arg in c("par", "logpar")) {
      try({
        print(scale_arg)
        formula_1_lm <- paste0(scale_arg, " ~ 1 + age")
        formula_2_lm <- paste0(scale_arg, " ~ 1 + age + sex")
        
        lm_fits <- list(
          fit1_lm = lm(formula(formula_1_lm), for_this_param, weights = weights),
          fit2_lm = lm(formula(formula_2_lm), for_this_param, weights = weights)
        )
        
        
        formula_0_lmer <- paste0(scale_arg, " ~1 +(1|patient_id)")
        formula_1_lmer <-
          paste0(scale_arg, "~ 1 + age+ (1|patient_id)")
        formula_2_lmer <-
          paste0(scale_arg, " ~ 1 + age+sex+ (1|patient_id)")
        formula_3_lmer <-
          paste0(scale_arg, " ~ 1 + sex+ (1|patient_id)")
        
        lmer_fits <- list()
        if (length(for_this_param$patient_id) > length(unique(for_this_param$patient_id))) {
          #lme4 here as lmerTest does not cooperate with tidy yet
          lmer_fits <- list(
            fit0_lmer = lmerTest:::lmer(
              formula(formula_0_lmer),
              data = for_this_param,
              REML = TRUE,
              weights = weights
            ),
            fit1_lmer = lmerTest:::lmer(
              formula(formula_1_lmer),
              data = for_this_param,
              REML = TRUE,
              weights = weights
            ),
            fit2_lmer = lmerTest:::lmer(
              formula(formula_2_lmer),
              data = for_this_param,
              REML = TRUE,
              weights = weights
            ),
            fit3_lmer = lmerTest:::lmer(
              formula(formula_3_lmer),
              data = for_this_param,
              REML = TRUE,
              weights = weights
            )
          )
        }
        
        all_fits <- c(lmer_fits, lm_fits)
        
        if (!(coef(all_fits[["fit1_lm"]]) %>% is.na() %>% any())) {
          if (length(lmer_fits) > 0) {
            model_comparison <- anova(
              lmer_fits[["fit0_lmer"]],
              lmer_fits[["fit1_lmer"]],
              lmer_fits[["fit2_lmer"]],
              lmer_fits[["fit3_lmer"]],
              lm_fits[["fit1_lm"]],
              lm_fits[["fit2_lm"]],
              lm_fits[["fit3_lm"]]
            ) %>%
              tidy() %>%
              mutate(term = c(names(all_fits))) %>%
              rename("test" = "term")
            #   print(model_comparison )
            model_comparison <- model_comparison %>%
              mutate(par = !!your_param,
                     scale_arg = !!scale_arg)
            model_comparison$isSingu = sapply(model_comparison$test, function(x)
              ifelse(
                grepl("lmer", x),
                isSingular(lmer_fits[[x]]),
                coef(lm_fits[[x]]) %>% is.na() %>% any()
              ))
            #   print(model_comparison)
            
          } else{
            model_comparison <-
              anova(lm_fits[["fit1_lm"]], lm_fits[["fit2_lm"]]) %>%
              tidy()
            
            model_comparison$test <- c("fit1_lm", "fit2_lm")
            model_comparison <- model_comparison %>%
              mutate(par = !!your_param,
                     scale_arg = !!scale_arg)
            model_comparison$isSingu = sapply(model_comparison$test, function(x)
              ifelse(
                grepl("lmer", x),
                isSingular(lmer_fits[[x]]),
                coef(lm_fits[[x]]) %>% is.na() %>% any()
              ))
            model_comparison$AIC <- sapply(lm_fits, AIC)
            
            #    print(model_comparison)
          }
        } else{
          if (length(lmer_fits) > 0) {
            model_comparison <- anova(lmer_fits[["fit0_lmer"]],
                                      lmer_fits[["fit1_lmer"]],
                                      lmer_fits[["fit2_lmer"]],
                                      lmer_fits[["fit3_lmer"]]) %>%
              tidy() %>%
              rename("test" = "term") %>%
              mutate(par = !!your_param,
                     scale_arg = !!scale_arg)
            
            model_comparison$isSingu = sapply(model_comparison$test, function(x) {
              print(x)
              ifelse(
                grepl("lmer", x),
                isSingular(lmer_fits[[x]]),
                coef(lm_fits[[x]]) %>% is.na() %>% any()
              )
            })
            #     print(model_comparison)
          }
          
        }
        
        
        #now pick the best model by AIC
        model_comparison_winner <- model_comparison %>%
          slice(which.min(AIC))
        
        #extract pval for age, sex
        best_model = model_comparison_winner %>% filter(p.value < relaxed_significance_cutoff) %>%
          .$test
        
        #here test and estimate as below only if significant (below relaxed_significance_cutoff)
        
        if (length(best_model) > 0) {
          res <-  all_fits[[best_model]]
          #    print(res)
          if (class(res) %in% c("merMod", "lmerModLmerTest")) {
            res <-
              broom.mixed:::tidy.merMod(
                x = res,
                effects = "fixed",
                conf.int = T,
                conf.level = .9
              )
          } else{
            res <- tidy(x = res,
                        conf.int = T,
                        conf.level = .9)
          }
          #   print( res)
          res <- res %>%
            mutate(
              par = !!your_param,
              best_model = !!best_model,
              scale_arg = !!scale_arg
            ) %>%
            pivot_wider(
              data = .,
              id_cols = c("par", "best_model", "scale_arg"),
              names_from = "term",
              values_from = setdiff(
                colnames(.),
                c("par", "best_model", "scale_arg", "term")
              )
            )
          
          results[[scale_arg]] <-
            bind_cols(model_comparison_winner, res)
        } else{
          NULL
        }
      })
    }
  }
  
  resu[[your_param]] <- bind_rows(results)
  
}

resu_rb <- lapply(resu, as_tibble) %>%
  bind_rows()


##----signif_sex_age---------------------------------------------------------------------------------------------------------
###Identify parameters for which sex/age are significant, identify which model
##was finally detected as best on desired significance level.
#This part requires detailed results from model fitting -  age specification, which we cannot release due to anonymity reasons.


resu_rb <- resu_rb %>%
  mutate(
    what_signif_0.01 = case_when(
      best_model %in% c("fit1_lm", "fit2_lm") ~
        case_when(
          p.value_age < cutoff_age_sex &
            p.value_sexM < cutoff_age_sex ~ "age+sex",
          p.value_age < cutoff_age_sex &
            (p.value_sexM > cutoff_age_sex | is.na(p.value_sexM)) ~ "age",
          (p.value_age > cutoff_age_sex |
             is.na(p.value_age)) & p.value_sexM < cutoff_age_sex ~ "sex",
          TRUE ~ "none"
        ),
      best_model == "fit1_lmer" &
        p.value_age < cutoff_age_sex ~ "age",
      best_model == "fit2_lmer" &
        p.value_age < cutoff_age_sex & p.value_sexM < cutoff ~ "age+sex",
      best_model == "fit2_lmer" &
        p.value_age < cutoff_age_sex & p.value_sexM > cutoff ~ "age",
      best_model == "fit2_lmer" &
        p.value_age > cutoff_age_sex & p.value_sexM < cutoff ~ "sex",
      best_model == "fit3_lmer" &
        p.value_sexM < cutoff_age_sex ~ "sex",
      TRUE ~ "none"
    )
  ) %>%
  select(1:12, 16:42, 47)

colnames(resu_rb) <- gsub("\\.\\.\\..*", "", colnames(resu_rb))
