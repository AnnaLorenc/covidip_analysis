#Combine results across tests for the final table

full_results <- list(
  whole5 %>%
    select(about_weights, par_sexage , logpar_sexage),
  results_seropos,
  results_all_par,
  results_all_logpar,
  results_sn_par,
  results_sn_logpar,
  results_corrected
) %>% reduce(., .f = full_join)


#Define and add _final set of columns, chosing values according to: pvals in
#sexage; seropos;  log or linear according to parameter type.
columns_pivotal <-
  c("par_sexage",
    "logpar_sexage",
    "pvalseropos_par",
    "pvalseropos_log")

columns_all  <- setdiff(colnames(full_results), columns_pivotal)

cols_for_extraction <- full_results %>%
  group_split() %>%
  lapply(., function(x)
    select_columns(
      full_results_param = x,
      columns_pivotal = columns_pivotal,
      columns_all = columns_all,
      cutoff_pval_snonly = cutoff_pval_snonly
    ))

names(cols_for_extraction) <- full_results %>%
  group_split() %>%
  sapply(., function(x)
    x$param)


too_look_for <- c(
  "param",
  "pval_seropos_final" ,
  "about_weights",
  "healthy_.*_mean",
  "COVID_mean",
  "COVIDlow_mean",
  "COVIDmoderate_mean",
  "COVIDsevere_mean",
  "LRTI_mean",
  "effsize",
  "healthy_.*vs_COVID.pval.*" ,
  "healthy_.*vs_LRTI.pval.*" ,
  "COVIDlow_vs_healthy.*pval.*",
  "COVIDmoderate_vs_healthy.*pval.*",
  "COVIDsevere_vs_healthy.*pval.*",
  "LRTI_vs_COVID",
  "COVIDlow_vs_LRTI",
  "COVIDmoderate_vs_LRTI",
  "COVIDsevere_vs_LRTI",
  "COVIDsevere_vs_COVIDlow",
  "COVIDsevere_vs_COVIDmoderate",
  "COVIDmoderate_vs_COVIDlow"
)

colnames_toget <- sapply(1:length(cols_for_extraction), function(i)
  sapply(too_look_for, function(x) {
    a = grep(x, cols_for_extraction[[i]], val = T)
    if (length(a) == 0) {
      a = x
    }
    return(a)
  })) %>% t()

colnames_toget <- data.frame(paramet = names(cols_for_extraction), colnames_toget)

final_selection <-
  lapply(1:nrow(colnames_toget), function(i)
    full_results %>% select(colnames_toget[i, -1] %>% unlist()) %>% filter(param ==
                                                                       colnames_toget$paramet[i])) %>%
  bind_rows()


colnames(final_selection) <-
  paste0(
    c(
      "param",
      "pval_seropos_final" ,
      "about_weights",
      "healthy_mean",
      "COVID_mean",
      "COVIDlow_mean",
      "COVIDmoderate_mean",
      "COVIDsevere_mean",
      "LRTI_mean",
      "healthy_vs_COVID.effsize",
      "healthy_vs_COVID.pval" ,
      "healthy_vs_LRTI.pval" ,
      "COVIDlow_vs_healthy.pval",
      "COVIDmoderate_vs_healthy.pval",
      "COVIDsevere_vs_healthy.pval",
      "LRTI_vs_COVID.pval",
      "COVIDlow_vs_LRTI.pval",
      "COVIDmoderate_vs_LRTI.pval",
      "COVIDsevere_vs_LRTI.pval",
      "COVIDsevere_vs_COVIDlow.pval",
      "COVIDsevere_vs_COVIDmoderate.pval",
      "COVIDmoderate_vs_COVIDlow.pval"
    ),
    "_final"
  )


#Main result table (also on the webpage)

full_results_better <- full_join(full_results,
                                 final_selection %>% select(-c(
                                   "pval_seropos_final_final" , "about_weights_final"
                                 )),
                                 by = c("param" = "param_final"))
