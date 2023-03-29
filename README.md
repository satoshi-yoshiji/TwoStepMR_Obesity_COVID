# TwostepMR_obesity_COVID  

### This repository contains codes for the following project:
---
## Proteome-wide Mendelian randomization implicates nephronectin as an actionable mediator of the effect of obesity on COVID-19 severity  
Satoshi Yoshiji, Guillaume Butler-Laporte, Tianyuan Lu, Julian Daniel Sunday Willett, Chen-Yang Su, Tomoko Nakanishi, David R. Morrison, Yiheng Chen, Kevin Liang, Michael Hultström, Yann Ilboudo, Zaman Afrasiabi, Shanshan Lan, Naomi Duggan, Chantal DeLuca, Mitra Vaezi, Chris Tselios, Xiaoqing Xue, Meriem Bouab, Fangyi Shi, Laetitia Laurent, Hans Markus Münter, Marc Afilalo, Jonathan Afilalo, Vincent Mooser, Nicholas J Timpson, Hugo Zeberg, Sirui Zhou, Vincenzo Forgetta, Yossi Farjoun, J. Brent Richards  
### Nature Metabolism (Nat Metab 5, 248–264 (2023). https://doi.org/10.1038/s42255-023-00742-w
---
#### Two-step MR approach:  
We identified circulating protein mediators of the effect of obesity on COVID-19 severity using a two-step MR approach. Each step is described below.

#### Step 1 MR: Identifying the effect of BMI on plasma protein levels  
・ Fist step is to generate instrumenta variables by running `00.remove_MHC_and_clump.R` in `1.1.BMI_to_proteins/`  
・ Then, run `01.MR.R` in parallel using `01.run_MR_in_parallel.sh`. The script performs MR analyses using PBS array jobs (You may use slurm array jobs instead). To run this script, You have to provide a nested list of pQTL paths (replace `/scratch/richards/satoshi.yoshiji/11.pQTL/decode_batch/listbatch.txt` with your paths).  
· pQTL used in step 1 MR can be found at https://www.decode.com/summarydata/  
⋅ `02.summarize_mr_results.sh` will collect the MR results.  
· For reverse MR, go to `1.2.BMI_reverse`. You can run reverse MR (i.e., from proteins to BMI) with `01.MR.R` in parallel using `01.run_MR_in_parallel.sh`.  
⋅ `02.summarize_mr_results.sh` will collect the reverse MR results.  
· `01.filter_mr_results.R` in the `1.3.collect_BMI_driven_proteins/` directory will curate the results for step 1 MR.  

#### Step 2 MR: Identifying the effect of BMI-driven proteins on COVID-19 severity
· The COVID-19 severity outcomes consists of two outcomes: critically ill COVID-19 (HGIr7_EUR_A2) and COVID-19 hospitalization (HGIr7_EUR_B2).  
⋅ In each directory, you can run `01.MR.R` in parallel using `01.run_MR_in_parallel.sh`. This performs MR analyses for each protein using PBS array jobs (You may use slurm array jobs instead). To run this script, You have to provide a nested list of pQTL path of BMI-driven proteins, identified in Step 1 MR (replace `/scratch/richards/satoshi.yoshiji/09.proMR/14.BMI_noMHC_proxy/2.2.proteins_for_step2_without_reverse/pQTL_batch/listbatch.txt` with your paths).  
· cis-pQTL used in step 2 MR can be found in Supplementary Table 2 of the deCODE study (Ferkingstad, E. et al. Nat Genet 2021).
⋅ `02.summarize_mr_results.sh` will collect the MR results.  

#### Colocalization: 
· The colocalization analyses evaluate whether cis-pQTL of the putatively causal proteins (NPNT and HSD17B14) and COVID-19 severity outcomes share a single causal variant.  
· `1.coloc_grid_NPNT_pQTL_and_HGIr7_A2.loop.R` performs colocalization analysis of cis-pQTL for NPNT with critically ill COVID-19  
· `2.coloc_grid_NPNT_pQTL_and_HGIr7_B2.loop.R` performs colocalization analysis of cis-pQTL for NPNT with COVID-19 hospitalization  
· `3.coloc_grid_HSD17B14_pQTL_and_HGIr7_A2.loop.R` performs colocalization analysis of cis-pQTL for HSD17B14 with critically ill COVID-19  
· `4.coloc_grid_HSD17B14_pQTL_and_HGIr7_B2.loop.R` performs colocalization analysis of cis-pQTL for HSD17B14 with COVID-19 hospitalization  
· Note that grid.csv provides parameters for prior probabilities for coloc.

For more information, please refer to our publication at Nature Metabolism (Nat Metab 5, 248–264 (2023). https://doi.org/10.1038/s42255-023-00742-w

