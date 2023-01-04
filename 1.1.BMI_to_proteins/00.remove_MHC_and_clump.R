library(tidyverse)
library(vroom)
library(TwoSampleMR)
library(ieugwasr)

################
# BMI GWAS
# Donwloaded from https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz
# Even though these SNPs are already conditionally independent, we took a conservative approach and performed conventional clumping
################
gwas <- vroom('/scratch/richards/satoshi.yoshiji/09.proMR/01.exposure/BMI/Meta-analysis_Locke_et_al+UKBiobank_2018_top_941_from_COJO_analysis_UPDATED.txt')

################
# remove MHC region
# chr6:28,477,797-33,448,354 (GRCh37)
# https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
#################

# select MHC region
mhc <- gwas %>%
dplyr::arrange(CHR, POS) %>%
dplyr::filter(CHR == 6) %>%
dplyr::filter(POS >= 28477797 & POS <= 33448354)

# create MHC region's snp list
mhcsnp <- mhc %>% dplyr::select(SNP)

# remove MHC's snp from the full GWAS
gwas_nomhc <- gwas %>% dplyr::filter(!(SNP %in% unlist(mhcsnp)))

#format exposure
exp_dat <- format_data(dat = gwas_nomhc,
                       type = 'exposure',
                       snp_col = 'SNP',
                       beta_col = 'BETA',
                       se_col = 'SE',
                       eaf_col = 'Freq_Tested_Allele_in_HRS',
                       effect_allele_col = 'Tested_Allele',
                       other_allele_col = 'Other_Allele',
                       pval_col = 'P',
                       samplesize_col = 'N')

##################
# clumping
##################

exp_clump <- exp_dat %>%
  mutate(rsid=exp_dat$SNP) %>%
  mutate(pval=exp_dat$pval.exposure)
exp_clump <- ld_clump(dat = exp_clump, plink_bin = "/home/richards/satoshi.yoshiji/local/bin/plink", bfile = "/home/richards/satoshi.yoshiji/scratch/database/1KG/EUR/EUR", clump_p = 5e-8) # clump_kb = 10000, clump_r2 = 0.001

# save results
write.table(exp_clump, file='BMI_noMHC_clumped.tsv', sep='\t', quote=F, row.names=F)
