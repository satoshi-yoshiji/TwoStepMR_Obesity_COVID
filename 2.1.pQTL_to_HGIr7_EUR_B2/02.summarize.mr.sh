#!/bin/bash
#SBATCH --job-name=sum
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH -t 1:00:00
#SBATCH --account=rrg-brichard

#cd $SLURM_SUBMIT_DIR
rm *.tmp

#=====================
# the script starts here
#=====================
readlink -f output/or/* | sort -V > or_list.txt

for protein in `cat or_list.txt`
do

protein_name=`basename $protein ".or.txt"`
pleio=output/pleio/${protein_name}.pleio.txt
hetero=

################
# collect Wald rato results
################
wald_result=`grep "Wald" $protein`
awk -v prot="$protein_name" -v res="$wald_result" 'BEGIN {OFS="\t"; print prot, res}' >> summary_wald_or.tmp

################
# collect IVW results
################
ivw_result=`grep "Inverse" $protein`
awk -v prot="$protein_name" -v res="$ivw_result" 'BEGIN {OFS="\t"; print prot, res}' >> summary_ivw_or.tmp

################
# collect Egger results
################
egger_result=`grep "Egger" $protein`
awk -v prot="$protein_name" -v res="$egger_result" 'BEGIN {OFS="\t"; print prot, res}' >> summary_egger_or.tmp

################
# collect pleio results
################
pleio_result=`awk '(NR==2) {print $0}' output/pleio/${protein_name}.pleio.txt`
awk -v prot="$protein_name" -v res="$pleio_result" 'BEGIN {OFS="\t"; print prot, res}' >> summary_pleio_or.tmp

################
# collect hetero results
################
hetero_result=`grep "Inverse" output/hetero/${protein_name}.hetero.txt`
awk -v prot="$protein_name" -v res="$hetero_result" 'BEGIN {OFS="\t"; print prot, res}' >> summary_hetero_or.tmp

################
# combine Wald and IVW results
################

cat summary_wald_or.tmp summary_ivw_or.tmp summary_egger_or.tmp >> summary_results.tmp

################
# collect steiger results
################
steiger_result=`awk '(NR==2) {print $0}' output/steiger/${protein_name}.steiger.txt`
awk -v prot="$protein_name" -v res="$steiger_result" 'BEGIN {OFS="\t"; print prot, res}' >> steiger.tmp

done

awk 'BEGIN {print "protein\tid.exposure\tid.outcome\toutcome\texposure\tmethod\tnsnp\tb\tse\tpval\tlo_ci\tup_ci\tor\tor_lci95\tor_uci95"}1' summary_results.tmp > output/summary_results.txt

awk 'BEGIN {print "protein\tid.exposure\tid.outcome\toutcome\texposure\tegger_intercept\tse\tpval"}1' summary_pleio_or.tmp > output/summary_pleio.txt

awk 'BEGIN {print "protein\tid.exposure\tid.outcome\toutcome\texposure\tmethod\tQ\tQ_df\tQ_pval\tisquared"}1' summary_hetero_or.tmp > output/summary_hetero.txt

awk 'BEGIN {print "id.exposure\tid.outcome\texposure\toutcome\tsnp_r2.exposure\tsnp_r2.outcome\tcorrect_causal_direction\tsteiger_pval"}1' steiger.tmp > output/summary_steiger.txt

#awk 'BEGIN {print "protein\tid.exposure\tid.outcome\toutcome\texposure\tegger_intercept\tse\pval"}1' summary_or.tmp > summary_or.txt
rm *.tmp;

