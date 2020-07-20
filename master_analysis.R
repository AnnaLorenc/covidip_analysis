
# Scripts to perform analyses reported in
# https://www.medrxiv.org/content/10.1101/2020.06.08.20125112v1. 

# Script assumes that in the directory data/data_from_webpage/ there are
# per-sample and per-patient results from
# https://www.immunophenotype.org/index.php/data/bulk-data-downloads/, so please
# place these two files there or change sample_datalocation and
# patient_datalocation accordingly.
#
# For the privacy reasons we are unable to provide exact age of the patients.
# Hence, we make the code to check for sex/age dependency of immune parameters
# available, but you won't be able to run it with the real source data. Instead
# we provide necessary intermediate files for the parts of the pipeline which
# depend on age/sex correction.
##
##Contact 
##code: anna.lorenc@kcl.ac.uk
##project: covidip@kcl.ac.uk


source('scripts/modelling_functions_0.1.R')
dir.create(path = "outputs")

## Read data in. Assumes the data in the /data subdirectory of the current directory, change if necessary.
sample_datalocation <-  "data/data_from_webpage/2020-06-30flow_sero_cyto__ifnaexport.csv"
patient_datalocation <- "data/data_from_webpage/pat_upd_dates20-07-02.csv"

source('scripts/read_data_in.R', echo=TRUE)


######## Main statistical testing framework
## Test for differences between seropositive and seronegative controls.
## Resulting diff_within_controls contains parameters which are different at
## significance level cutoff_serpos_seroneg
cutoff_serpos_seroneg <- 0.05

source('scripts/seropos_seroneg.R', echo=TRUE)
write_csv(results_seropos, path = "outputs/results_seropos.csv")


###Using controls only, test for significant influence of age/sex. This part
#requires age specification, which we cannot release due to anonymity reasons.
#Do not run - results of this testing are provided.
#cutoff_age_sex <- 0.01
#source('scripts/age_sex_significance.R', echo=TRUE)

resu_rb <-read_csv("data/additional_data/resu_rb_results.csv")


###Identify parameters with different intra-individual variance in sick than in healthy (changing_with_without)
##To be conservative, use only main covid-ip cohort

cutoff_variance_test <- 0.01
source('scripts/test_variance.R', echo=TRUE)
write_csv(changing_with_without, path = "outputs/changing_with_without.csv")

##main statistical tests
#Prepare data: gather per parameter info:
# for which contrasts there is enough data,
# whether weights are to be used (based on intraindividual variance),
# whether all controls might be used (based on seonegative-seropositive differences,
# when difference among controls with pva<0.05 or impossible to compute)
#takes a while to compute so consider saving output as suggested
source('scripts/tests_main.R', echo=TRUE)
write_rds(whole5, path="outputs/whole5.rds")

##extract relevant pvals/effect sizes etc
source('scripts/combine_results_raw.R', echo=TRUE)

##Tests after correcting for age/sex for the relevant subset of parameters. 
#This part requires age specification, which we cannot release due to anonymity reasons. Do not run - results of this testing are provided.
#source('scripts/tests_sexage.R', echo=TRUE)
#source('scripts/combine_results_sexagecorrect.R', echo=TRUE)
results_corrected <- read_csv(file ="data/additional_data/results_corrected_withCF_2020-06-27.csv")
cutoff_pval_snonly <- 0.05 #when to use only seronegative controls for the final data
source('scripts/prepare_final_result_table.R', echo=TRUE)

write_csv(full_results_better, path = "outputs/full_results.csv")



######## Additional tests

####Compare peak antibody values
#### Fig 2a and FigS1d - peak antibody comparison
source("scripts/peak_abs.R")
write_csv(ab_peak, path='outputs/ab_peak.csv')

###Prognostic test
### Fig 6a-e
source("scripts/prognostic_figure.R")
write_csv(fig6a_tests, path='outputs/trajectory_stats.csv')
