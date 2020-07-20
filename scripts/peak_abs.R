#Compare peak Abs values per patient


ab_peak <-
  pat_flow_sero_cyto_ifna %>% select(
    patient_id,
    starts_with("class_ds") & !ends_with("sample"),
    ends_with("sero") &
      !contains("isSig") & !starts_with("N_")
  ) %>%
  group_by(patient_id) %>%
  mutate_at(vars(ends_with("sero")), ~ max(.x, na.rm = T)) %>%
  unique() %>%
  ungroup %>%
  filter(class_dss != "LRTI_nonCovid") %>%
  pivot_longer(cols = RBD_IgG_LIPS_sero:RBD_IgM_Norm_sero) %>%
  group_by(name)

ab_peak %>%
  rstatix::dunn_test(., value ~ class_dss) %>%  
  mutate(comparison = paste(group1, group2, sep = "-")) %>%
  select(name, p:comparison) %>%
  pivot_wider(
    id_cols = name,
    names_from = comparison,
    values_from = p:p.adj.signif
  )
