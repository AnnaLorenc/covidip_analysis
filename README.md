# covidip_analysis

This repo contains scripts to perform analyses reported in [https://www.medrxiv.org/content/10.1101/2020.06.08.20125112v1].  To replicate our analysis execute the workflow from the file master_analysis.R or its parts.
Script assumes thad in the directory data/data_from_webpage/ there are per-sample and per-patient results from https://www.immunophenotype.org/index.php/data/bulk-data-downloads/, so please place these two files there or change sample_datalocation and patient_datalocation accordingly.

For the privacy reasons we are unable to provide exact age of the patients. Hence, we make the code to check for sex/age dependency of immune parameters available, but you won't be able to run it with the real source data. Instead we provide necessary intermediate files for the parts of the pipeline which depend on age/sex correction. 
