install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)



#########################################################
##### EXPOSURE DATA #####################################
#########################################################

## USING DATA FROM OPEN GWAS RESOURCE 

ieugwasr::get_access_token()

Sys.setenv("OPENGWAS_JWT" = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJlbTg5MkBjYW0uYWMudWsiLCJpYXQiOjE3MTI3NjAxNTAsImV4cCI6MTcxMzk2OTc1MH0.Xe9gsjvlB91byBAaUC7bfpxbQxWDuEZZ2ZNwvWAHIyp8JIZ3PfA7oetwnwwQQiwsHftQv6vDvpJQgrawV_ip6FNhdBAJU9lh-DDWqw62RmHN2Vg9pgtpxgEpDGDNlP5wyMA6TywG9txWH7aJclJO5BEgipx8Gqh-qa-3nfmUDxJEcCUEfvjdemj1tK6rd2ZAclL7DaYrjONKpbMlz6lmWqyM212hJfs2KZyB0e5f_jx79Hn21WmMzqff6ILIFrT93cWyl1T0XPL-jgZ9q1ww-4fTQMLYVEqBLoktK_TI5fjo_NETjB8k7baIMRxudLfUSwPHx2PLbFyQba1pnzH6kg")

exposure_dat <- extract_instruments("ieu-a-2")

outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes = "ieu-a-7")

# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
res <- mr(dat)





## USING DATA FROM CSV OR TXT FILES 

# here using the example of a csv file which exists in the package 

bmi2_file <- system.file("extdata/bmi.csv", package = "TwoSampleMR")

# reading in this data 

bmi_exp_dat <- read_exposure_data(
  filename = bmi2_file,
  sep = ",",
  snp_col = "rsid",
  beta_col = "effect",
  se_col = "SE",
  effect_allele_col = "a1",
  other_allele_col = "a2",
  eaf_col = "a1_freq",
  pval_col = "p-value",
  units_col = "Units",
  gene_col = "Gene",
  samplesize_col = "n"
)

# checking it 

head(bmi_exp_dat)


# if the phenotype column is not provided we will assume that the 
# phenotype's name is simply exposure 
# this is entered in the exposure column and can be renamed 
# manually 

bmi_exp_dat$exposure <- "BMI"





## IF THE DATA ALREADY EXISTS IN A DATAFRAME IN R 

# can convert into the correct format using the format_data() function

random_df <- data.frame(
  SNP = c("rs1", "rs2"),
  beta = c(1, 2),
  se = c(1, 2),
  effect_allele = c("A", "T")
)
random_df


# this can be formatted like so

random_exp_dat <- format_data(random_df, type = "exposure")

random_exp_dat





## OBTAINING INSTRUMENTS FROM EXISTING CATALOGUES 

# a number of sources of instruments have already been curated and
# are available for use 
# they are provided as data objects in the MRInstruments package 
# to install these:

remotes::install_github("MRCIEU/MRInstruments")


# ACCESSING THE DATAFRAMES 
library(MRInstruments)
data(gwas_catalog)
head(gwas_catalog)


# to obtain instruments for body mass index using the Speliotes et al 2010 study 
bmi_gwas <-
  subset(gwas_catalog,
         grepl("Speliotes", Author) &
           Phenotype == "Body mass index")
bmi_exp_dat <- format_data(bmi_gwas)
bmi_exp_dat




## PROTEINS 

# independent top hits from GWASs on 47 protein
# levels in whole blood are stored in the proteomic_qtls 
# data object 

data(proteomic_qtls)
head(proteomic_qtls)


# to obtain instruments for the ApoH protein 

apoh_exp_dat <-
  format_proteomic_qtls(subset(proteomic_qtls, analyte == "ApoH"))

apoh_exp_dat



## GENE EXPRESSION LEVELS 

# independent top hits from GWASs on 32432 gene identifiers 
# and in 44 tissues are available from the GTEx study

data(gtex_eqtl)
head(gtex_eqtl)

# to obtain instruments for the IRAK1BP1 gene expression 
# levels in subcutaneous adipose tissue:

irak1bp1_exp_dat <-
  format_gtex_eqtl(subset(
    gtex_eqtl,
    gene_name == "IRAK1BP1" & tissue == "Adipose Subcutaneous"
  ))
irak1bp1_exp_dat




## IEU GWAS DATABASE 

# The IEU GWAS database contains the entire summary statistics for thousands of GWASs. 
# You can use this database to define the instruments for a particular exposure. 
# You can also use this database to obtain the effects for constructing polygenic 
# risk scores using different p-value thresholds

# to obtain a list and details of the available GWASs 
ao <- available_outcomes()
head(ao)


# To extract instruments for a particular trait using a particular study, for 
# example to obtain SNPs for body mass index using the Locke et al. 2015 GIANT
# study, you specify the study ID as follows:

bmi2014_exp_dat <- extract_instruments(outcomes = 'ieu-a-2')

# this returns a set of LD clumped SNPs that are GWAS 
# significant for BMI




# CLUMPING 

# For standard two sample MR it is important to ensure that the instruments for 
# the exposure are independent. Once instruments have been identified for an 
# exposure variable, the IEU GWAS database can be used to perform clumping.

# You can provide a list of SNP IDs, the SNPs will be extracted from 1000 
# genomes data, LD calculated between them, and amongst those SNPs that have 
# LD R-square above the specified threshold only the SNP with the lowest P-value
# will be retained. To do this, use the following command:

bmi_exp_dat <- clump_data(bmi2014_exp_dat)

# Note that for the instruments in the MRInstruments package the SNPs are already LD clumped.
# If a variant is dropped from your unclumped data it could be because it is absent from the reference panel




#################################################################
######### OUTCOME DATA #########################################
################################################################

ao <- available_outcomes()


## EXTRACTING PARTICULAR SNPs from particular studies 

# we need to identify the SNPs that influence the exposure and 
# then extract those snps from the outcome 

# so first we collect our exposure SNPs 
bmi_exp_dat <- extract_instruments(outcomes = 'ieu-a-2')


# now we find a suitable GWAS for coronary heart disease 
ao[grepl("heart disease", ao$trait), ]

# the cardiogram is the most recent one so we extract that 

chd_out_dat <- extract_outcome_data(
  snps = bmi_exp_dat$SNP,
  outcomes = 'ieu-a-7'
)


# the extract outcome data function is flexible
# the snps argument only requires an array of rsids 
# the outcomes argument can be a vector of outcomes eg.

chd_out_dat2 <- extract_outcome_data(
  snps = c("rs234", "rs17097147"), 
  outcomes = c('ieu-a-2', 'ieu-a-7')
)

# this extracts two snps from each of the outcomes 



## PROXIES 
# we add parameters to the outcome search to include proxies 

# proxies = TRUE or FALSE 
# rsq = numeric value of r2 
# palindromes = allow palindromic SNPs?
# maf threshold = if palindromes allowed what is the maximum minor allele frequency
# of palindromes allowed 




## USING LOCAL GWAS SUMMARY DATA 
# this may be data that is not present in the IEU database
# can still be used to perform analysis 

# suppose there was a GWAS summary file called "gwas_summary.csv" with 
# two million rows

outcome_dat <- read_outcome_data(
  snps = bmi_exp_dat$SNP,
  filename = "gwas_summary.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "effect",
  se_col = "SE",
  effect_allele_col = "a1",
  other_allele_col = "a2",
  eaf_col = "a1_freq",
  pval_col = "p-value",
  units_col = "Units",
  gene_col = "Gene",
  samplesize_col = "n"
)

# this returns an outcome dataframe with only the SNPs that 
# were requested 









###############################################################
################### HARMONISE DATA ############################
###############################################################

# The IEU GWAS database contains data that is already harmonised, meaning that 
# the non-effect allele is aligned to the human genome reference sequence 
# (build 37). It’s still recommended to harmonise, but in principle everything 
# should be on the forward strand and effect alleles always relating to the 
# same allele. Some discrepancies could arise if there are multi-allelic
# variants that are represented as different bi-allelic variants in different studies.


# to harmonise the data we do the following:
dat <- harmonise_data(
  exposure_dat = bmi_exp_dat, 
  outcome_dat = chd_out_dat
)



# running the MR 

mr(dat, method_list = c("mr_egger_regression", "mr_ivw"))
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
res_single <- mr_singlesnp(dat)
res_loo <- mr_leaveoneout(dat)


# PLOTS 
res <- mr(dat)
p1 <- mr_scatter_plot(res, dat)
p1
ggsave(p1[[1]], file = "filename.pdf", width = 7, height = 7)


# forest plot
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]



# leave one out 
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]


# funnel plots 
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p4[[1]]



