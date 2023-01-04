library(TwoSampleMR)
library(ggplot2)

# working directory
wd <- "/scratch/richards/satoshi.yoshiji/github/TwostepMR_obesity_COVID/1.1.BMI_to_proteins/"
setwd(wd)
output_dir <- paste0(wd, "output/")

# exposure
exp_path <- '/scratch/richards/satoshi.yoshiji/09.proMR/01.exposure/BMI/BMI_noMHC_clumped.tsv'
exp_dat <- read_exposure_data(filename = exp_path,
                             sep='\t',
                             snp_col='SNP',
                             beta_col = 'beta.exposure',
                             se_col = 'se.exposure',
                             effect_allele_col = 'effect_allele.exposure',
                             other_allele_col = 'other_allele.exposure',
                             eaf_col = 'eaf.exposure',
                             pval_col = 'pval.exposure')

# outcome
args <- commandArgs(trailingOnly=TRUE) 
outcome_path <- args[1]
protname <- args[2]
# e.g.,
# outcome_path <- '/scratch/richards/public/decode_proteomics_2021/6342_10_NPNT_Nephronectin.txt.gz'
# protname <- '6342_10_NPNT_Nephronectin'

# For proxy search, you may use snappy v1.0 (https://gitlab.com/richards-lab/vince.forgetta/snappy/-/blob/master/snappy)

# standardize protein name
protname <- paste0(strsplit(protname,'_')[[1]][3], '.',  strsplit(protname,'_')[[1]][1], '_', strsplit(protname, '_')[[1]][2]) # e.g.m "NPNT.6342_10"

# import outcome GWAS
outcome_dat <-  read_outcome_data(
  snps= exp_dat$SNP,
  filename = outcome_path,
  sep='\t',
  snp_col='rsids',
  beta_col = 'Beta',
  se_col = 'SE',
  effect_allele_col = 'effectAllele',
  other_allele_col = 'otherAllele',
  eaf_col = 'ImpMAF',
  pval_col = 'Pval',
  samplesize_col = 'N'
)

# harmonize
exp_dat_outcome <- harmonise_data(exp_dat, outcome_dat, action = 2)
exp_data_outcome_name <- paste0(output_dir,"harmonized/", protname, ".harmonized.txt")
write.table(exp_dat_outcome, file=exp_data_outcome_name, sep = '\t', quote = F, row.names = F)

# mr 
mr_results <- mr(exp_dat_outcome, method_list = c('mr_ivw',
                                                  'mr_weighted_median',
                                                  'mr_weighted_mode',
                                                  'mr_egger_regression'))
mr_name <- paste0(output_dir, "mr/", protname, ".mr.txt")
write.table(mr_results, file=mr_name, sep = '\t', quote = F, row.names = F)

#odds ratio
OR <- generate_odds_ratios(mr_results)
OR_name <- paste0(output_dir, "or/", protname, ".or.txt")
write.table(OR, file=OR_name, sep = '\t', quote = F, row.names = F)

#scatter plot
# pdf_name <- paste0(output_dir, "pdf/", protname, ".scatter_narrow.pdf")
# pdf(pdf_name)
# mr_scatter_plot(mr_results, exp_dat_outcome)[[1]]
# dev.off()

#horizontal pleiotropy
pleio_res <- mr_pleiotropy_test(exp_dat_outcome)
pleio_name <- paste0(output_dir,"pleio/", protname, ".pleio.txt")
write.table(pleio_res, file=pleio_name, sep = '\t', quote = F, row.names = F)

#hetero test
hetero_res <- mr_heterogeneity(exp_dat_outcome)
hetero_res$isquared <- 100*(hetero_res$Q - hetero_res$Q_df)/hetero_res$Q  # I2 = 100%×(Q - df)/
hetero_name <- paste0(output_dir,"hetero/", protname, ".hetero.txt")
write.table(hetero_res, file=hetero_name, sep = '\t', quote = F, row.names = F)

# #forest
# res_single <- mr_singlesnp(exp_dat_outcome)
# forest <- mr_forest_plot(res_single)
# forest_name <- paste0(output_dir,"forest/", protname, ".forest.pdf")
# ggsave(forest_name, 
#        plot= forest[[1]],
#        width = 5,
#        height = 25
# )
