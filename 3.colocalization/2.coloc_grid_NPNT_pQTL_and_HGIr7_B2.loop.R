library(tidyverse)
library(vroom)
library(coloc)
#library(locuscomparer)
library(magrittr)

setwd("/scratch/richards/satoshi.yoshiji/09.proMR/24.colocalization_grid/1.NPNT_HGIr7_B2/27combinations/")

########### pQTL ###############
# NPNT pQTL
# rs34712979 
# 4:105897896 (GRCh38); 4:106819053 (GRCh37)

trait1_gwas <- vroom('/scratch/richards/public/decode_proteomics_2021/6342_10_NPNT_Nephronectin.txt.gz')
trait1_sel <- trait1_gwas %>% filter(Chrom == 'chr4') %>% filter(Pos >= (105897896 - 500000)) %>% filter(Pos <= (105897896 + 500000))

# rename
trait1_sel <- trait1_sel %>% 
  dplyr::rename(
    position = Pos,
    beta = Beta,
    SE = SE,
    ALT = effectAllele,
    REF = otherAllele,
    MAF = ImpMAF,
    N = N,
    pvalues = Pval,
    rsid = rsids
  )
trait1_sel %<>% dplyr::select(position, beta, SE, ALT, REF, MAF, N, pvalues, rsid)

# give varbeta
trait1_sel$varbeta <- (trait1_sel$SE)^2

# limit rsid to one
trait1_sel$rsid <- unlist(lapply(strsplit(trait1_sel$rsid, ","), function(x) x[1]))

# remove variants w/o rsid
trait1_rsid <- trait1_sel %>%  filter(!is.na(rsid))

# remove duplicated snp
trait1_nodup <- trait1_rsid %>% group_by(rsid, pvalues) %>% dplyr::arrange(pvalues) %>% dplyr::filter(!duplicated(rsid)) %>% ungroup()

# rename to pqtl_nodup and add id column
trait1_nodup$chr <- 4 
trait1_nodup %<>% unite(col = 'id', chr, position, REF, ALT, sep=":", remove = F) 

########### COVID ###############

trait2_gwas <- vroom('/scratch/richards/satoshi.yoshiji/09.proMR/disease_GWAS/COVID/HGI_r7/COVID19_HGI_B2_ALL_eur_leave23andme_20220403.tsv')
trait2_sel <- trait2_gwas %>% filter(`#CHR` == '4') %>% filter(POS >= (105897896 - 500000)) %>% filter(POS <= (105897896 + 500000))

# rename
trait2_sel <- trait2_sel %>% 
  dplyr::rename(
    chr = `#CHR`,
    position = POS,
    beta = all_inv_var_meta_beta,
    SE = all_inv_var_meta_sebeta,
    ALT = ALT,
    REF = REF,
    MAF = all_meta_AF,
    pvalues = all_inv_var_meta_p,
    MAF = all_meta_AF,
    rsid = rsid
  ) 

trait2_sel %<>% dplyr::select(chr, position, beta, SE, ALT, REF, MAF, pvalues, rsid) # Note that this has rsid (GTEx does not)

# give varbeta
trait2_sel$varbeta <- (trait2_sel$SE)^2

# limit rsid to one
trait2_sel$rsid <- unlist(lapply(strsplit(trait2_sel$rsid, ","), function(x) x[1]))

# remove variants w/o rsid
trait2_rsid <- trait2_sel %>%  filter(!is.na(rsid))

# remove duplicated snp
trait2_nodup <- trait2_rsid %>% group_by(rsid, pvalues) %>% dplyr::arrange(pvalues) %>% dplyr::filter(!duplicated(rsid)) %>% ungroup()

# rename to pqtl_nodup and add id column
trait2_nodup$chr <- 4 
trait2_nodup %<>% unite(col = 'id', chr, position, REF, ALT, sep=":", remove = F) 

####################
# keep only common SNPs
####################
trait1_nodup_norsid <- trait1_nodup %>% filter(!is.na(rsid)) # all of pQTL variants must have rsid

common_id <- unique(intersect(trait1_nodup_norsid$id, trait2_nodup$id))

trait1_common <- trait1_nodup %>% filter(id %in% common_id) 
trait2_common <- trait2_nodup %>% filter(id %in% common_id) 

#trait2_common$chr <- 4
write.table(trait1_common , file = 'trait1_1Mb_around_rs34712979.tsv', row.names = F, sep = '\t', quote = F)

#covid_common$chr <- 4
write.table(trait2_common , file = 'trait2_1Mb_around_rs34712979.tsv', row.names = F, sep = '\t', quote = F)

####################
#make lists
####################
# make it a list
trait1_list <- trait1_common %>% dplyr::select(-rsid) %>% as.list()
trait1_list$type <- 'quant'
trait1_list$N <- 35559
trait1_list$MAF <- trait1_common$MAF
str(trait1_list)

max(trait1_list$MAF)

# check dataset
check_dataset(trait1_list)
plot_dataset(trait1_list)

# make it a list
trait2_list <- trait2_common %>% dplyr::select(-rsid) %>% as.list()
trait2_list$type <- 'cc'
trait2_list$N <- (32519+2062805)
trait2_list$cc <- 32519/(32519+2062805)
str(trait2_list)

# check dataset
check_dataset(trait2_list)
plot_dataset(trait2_list)

# colocalizatiaon
grid <- read_csv('grid.csv', col_names = F)
colnames(grid) <- c('number', 'p1', 'p2', 'p12')

all_coloc_res <- data.frame()

for(i in 1:nrow(grid)){
  priors <- grid[i,]
  number <- priors[1,1] %>% unlist()
  p1 <- priors[1,2] %>% unlist()
  p2 <- priors[1,3] %>% unlist()
  p12 <- priors[1,4] %>% unlist()
  
  coloc_res <- coloc.abf(
    dataset1 = trait1_list,
    dataset2 = trait2_list,
    p1 = p1, p2 = p2, p12 = p12,
    #p1 = 1e-04, p2 = 1e-04, p12 = 1e-05 # coloc default 
  ) 
  
  coloc_res_df <- data.frame(t(coloc_res$summary))
  coloc_res_df %<>% mutate(p1 = p1, p2 = p2, p12 = p12) %>% relocate(p1, p2, p12)
  all_coloc_res <- rbind(all_coloc_res, coloc_res_df)
  # write_tsv(coloc_res_df, file = paste0(number, 'coloc_res_df_p1_', p1, '_p2_', p2, '_p12_', p12, '.tsv'))
}
write_tsv(all_coloc_res, file = 'all_coloc_res.tsv')
