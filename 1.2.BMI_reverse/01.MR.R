library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(vroom)

#####################
# Reverse MR  (proteins to BMI)
#####################

wd <- "/scratch/richards/satoshi.yoshiji/09.proMR/14.BMI_noMHC_proxy/2.1.pQTL_to_BMI_opengwas/" 
setwd(wd)
output_dir <- paste0(wd, "output/")

# read the outcome GWAS
outcome_path <- '/scratch/richards/satoshi.yoshiji/09.proMR/01.exposure/BMI_opengwas/ukb-b-19953.tsv' # Use UKB BMI here
outcome_GWAS <- vroom(outcome_path)
outcome_GWAS <- outcome_GWAS %>% dplyr::rename(SNP = ID) 

args <- commandArgs(trailingOnly=TRUE) 
exp_path <- args[1]
protname <- args[2]

# example
# exp_path <- '/scratch/richards/satoshi.yoshiji/11.pQTL/01.pQTL_Tianyuan/decodeaptamer_sep/NPNT.6342_10.tsv'
# protname <- 'NPNT.6342_10'

#exposure
exp_dat <- read_exposure_data(filename = exp_path,
                             sep='\t',
                             snp_col='variant',
                             beta_col = 'beta_unadj',
                             se_col = 'se',
                             effect_allele_col = 'Amin',
                             other_allele_col = 'Amaj',
                             eaf_col = 'MAF',
                             pval_col = 'pval')

# read the outcome
formatted_outcome <- format_data(outcome_GWAS, snps = exp_dat$SNP, type="outcome", snp_col="SNP", beta_col="ES", se_col = "SE", eaf_col = "AF", 
                                 effect_allele_col = "ALT", other_allele_col = "REF", pval_col = "pval", chr_col = "seqnames", pos_col = "start")
formatted_outcome$id.outcome <- 'outcome'

# For proxy search, snappy v1.0 was used (https://gitlab.com/richards-lab/vince.forgetta/snappy/-/blob/master/snappy)

# harmonize
exp_dat_outcome <-harmonise_data(exposure_dat=exp_dat, outcome_dat=formatted_outcome)

# exclude rare variants
exp_dat_outcome %<>% filter((eaf.exposure > 0.001 & eaf.exposure < 0.999) & (eaf.outcome > 0.001 & eaf.outcome < 0.999))
exp_data_outcome_name <- paste0(output_dir,"harmonized/", protname, ".harmonized.txt")
write_tsv(exp_dat_outcome, file=exp_data_outcome_name)

# mr result
mr_results <- mr(exp_dat_outcome)
mr_name <- paste0(output_dir, "mr/", protname, ".mr.txt")
write.table(mr_results, file=mr_name, sep = '\t', quote = F)

# odds ratio
OR <- generate_odds_ratios(mr_results)
OR_name <- paste0(output_dir, "or/", protname, ".or.txt")
write_tsv(OR, file=OR_name)

# # scatter plot
# pdf_name <- paste0(output_dir, "pdf/", protname, ".scatter.pdf")
# pdf(pdf_name)
# mr_scatter_plot(mr_results, exp_dat_outcome)[[1]] + xlim(0,1)
# dev.off()

# horizontal pleiotropy
pleio_res <- mr_pleiotropy_test(exp_dat_outcome)
pleio_name <- paste0(output_dir,"pleio/", protname, ".pleio.txt")
write_tsv(pleio_res, file=pleio_name)

# hetero test
tryCatch({
  res_single <- mr_singlesnp(exp_dat_outcome)
  res_single_beta <- res_single$b
  res_single_se <- res_single$se
  res_Isq <- Isq(res_single_beta, res_single_se)
  hetero_res <- mr_heterogeneity(exp_dat_outcome)
  hetero_res$i_squared <- res_Isq
  hetero_name <- paste0(output_dir,"hetero/", protname, ".hetero.txt")
  write_tsv(hetero_res, file=hetero_name)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
