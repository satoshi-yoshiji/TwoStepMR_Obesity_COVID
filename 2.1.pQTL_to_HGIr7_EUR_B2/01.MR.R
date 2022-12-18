library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(vroom)

wd <- "/scratch/richards/satoshi.yoshiji/09.proMR/14.BMI_noMHC_proxy/3.pQTL_to_HGIr7_EUR_B2/" #!!! don't forget the slash (/)at the end of the full path
setwd(wd)

output_dir <- paste0(wd, "output/")
system(paste0('mkdir -p ', output_dir, 'harmonized/'))
system(paste0('mkdir -p ', output_dir, 'or/'))
system(paste0('mkdir -p ', output_dir, 'hetero/'))
system(paste0('mkdir -p ', output_dir, 'pleio/'))
system(paste0('mkdir -p ', output_dir, 'steiger/'))

# read the outcome GWAS
outcome_path <- '/scratch/richards/satoshi.yoshiji/05.COVID_proteome/01.COVID_GWAS/HGI_r7/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv'
outcome_GWAS <- vroom(outcome_path)
outcome_GWAS %<>% dplyr::mutate(SNP = rsid) # Note this column has to be changed beforehand
outcome_GWAS %<>% dplyr::rename(seqnames = `#CHR`, start = POS ) #IEUGWAS format
  
args <- commandArgs(trailingOnly=TRUE) #input UNZIPPED file
exp_path <- args[1]
protname <- args[2]

print(protname)

# example
# exp_path <- '/scratch/richards/satoshi.yoshiji/11.pQTL/03.cispQTL_decode2021_sep/NPNT.6342_10.tsv'
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
exp_dat$samplesize <- 35559

# read the outcome
formatted_outcome <- format_data(outcome_GWAS, snps = exp_dat$SNP, type="outcome", snp_col="SNP", beta_col="all_inv_var_meta_beta", se_col = "all_inv_var_meta_sebeta", eaf_col = "all_meta_AF", 
                             effect_allele_col = "ALT", other_allele_col = "REF", pval_col = "all_inv_var_meta_p", chr_col = "seqnames", pos_col = "start", ncase_col = 'all_inv_var_meta_cases', ncontrol = 'all_inv_var_meta_controls')
formatted_outcome$id.outcome <- 'outcome'
formatted_outcome %<>% mutate(samplesize = ncase + ncontrol)

# harmonize
exp_dat_outcome <-harmonise_data(exposure_dat=exp_dat, outcome_dat=formatted_outcome)
# exclude rare variants
exp_dat_outcome %<>% filter((eaf.exposure > 0.001 & eaf.exposure < 0.999) & (eaf.outcome > 0.001 & eaf.outcome < 0.999)) 

exp_data_outcome_name <- paste0(output_dir,"harmonized/", protname, ".harmonized.txt")
write_tsv(exp_dat_outcome, file=exp_data_outcome_name)

# mr result
mr_results <- mr(exp_dat_outcome)
#mr_name <- paste0(output_dir, "mr/", protname, ".mr.txt")
#write.table(mr_results, file=mr_name, sep = '\t', quote = F)

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
# to add isquared
res_single <- mr_singlesnp(exp_dat_outcome)
res_single_beta <- res_single$b
res_single_se <- res_single$se
res_Isq <- Isq(res_single_beta, res_single_se)
hetero_res <- mr_heterogeneity(exp_dat_outcome)
hetero_res$i_squared <- res_Isq
hetero_name <- paste0(output_dir,"hetero/", protname, ".hetero.txt")
write_tsv(hetero_res, file=hetero_name)

# steiger
steiger <- directionality_test(exp_dat_outcome)
steiger_name <- paste0(output_dir, "steiger/", protname, ".steiger.txt")
write_tsv(steiger, file=steiger_name)
print(paste0(protname, ': done'))
