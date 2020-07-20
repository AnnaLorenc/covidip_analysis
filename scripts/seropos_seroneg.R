## ----seropos_seroneg--------------------------------------------------------------------------------------------------------
###Compare seropos_seroneg controls. Run tests on log/normal scale; use log10 for LIPS/cytokines and cell backcounts

seropos_seroneg_comparison <- fortests %>%
  group_by(param) %>%
  nest() %>%
  mutate(
    test_ready_controls = map(data, extract_controls, perform_log = FALSE),
    test_ready_controls_log = map(data, extract_controls, perform_log = TRUE),
    seropos_seroneg = map(test_ready_controls, test_ser_controls),
    seropos_seroneg_log = map(test_ready_controls_log, possibly(test_ser_controls, NULL)),
    pvalseropos_par = map(
      seropos_seroneg,
      ~ .x %>% filter(., term == "class_dsseropos_control") %>%
        select(p.value) %>% unlist(use.names =
                                     F)
    ),
    pvalseropos_log = map(
      seropos_seroneg_log,
      ~ ifelse(
        !is.null(.x),
        .x %>% filter(., term == "class_dsseropos_control") %>%
          select(p.value) %>% unlist(use.names =
                                       F),
        NA_real_
      )
    )
    
  ) %>%
  unnest(pvalseropos_par, pvalseropos_log)


seropos_seroneg_comparison <- seropos_seroneg_comparison %>%
  mutate(pval_seropos_final = ifelse(
    grepl("_cyto|LIPS$|Count_back", param),
    pvalseropos_log,
    pvalseropos_par
  ))

#Which parameters differ between seropos/seroneg?
diff_within_controls <- seropos_seroneg_comparison %>%
  filter(pval_seropos_final < cutoff_serpos_seroneg |
           is.na(pval_seropos_final)) %>%
  .$param


results_seropos <- seropos_seroneg_comparison%>%
  select(param, pvalseropos_par, pvalseropos_log, pval_seropos_final)



