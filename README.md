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
・ Then, run `01.MR.R`

#### Step 2 MR: Identifying the effect of BMI-driven proteins on COVID-19 severity
· The COVID-19 severity outcomes consists of two outcomes: critically ill COVID-19 (HGIr7_EUR_A2) and COVID-19 hospitalization (HGIr7_EUR_B2).  
⋅ In each directory, you can run `01.MR.R`

#### Colocalization: 
· The colocalization analyses evaluate whether cis-pQTL of the putatively causal proteins (NPNT and HSD17B14) and COVID-19 severity outcomes share a single causal variant.  
· `1.coloc_grid_NPNT_pQTL_and_HGIr7_A2.loop.R` performs colocalization analysis of cis-pQTL for NPNT with critically ill COVID-19  
· `2.coloc_grid_NPNT_pQTL_and_HGIr7_B2.loop.R` performs colocalization analysis of cis-pQTL for NPNT with COVID-19 hospitalization  
· `3.coloc_grid_HSD17B14_pQTL_and_HGIr7_A2.loop.R` performs colocalization analysis of cis-pQTL for HSD17B14 with critically ill COVID-19  
· `4.coloc_grid_HSD17B14_pQTL_and_HGIr7_B2.loop.R` performs colocalization analysis of cis-pQTL for HSD17B14 with COVID-19 hospitalization  
· Note that grid.csv provides parameters for prior probabilities for coloc.

If you use these codes, please cite our publication (Nat Metab 5, 248–264 (2023). https://doi.org/10.1038/s42255-023-00742-w

