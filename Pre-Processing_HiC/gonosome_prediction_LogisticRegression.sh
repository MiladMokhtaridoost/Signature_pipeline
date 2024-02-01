#!/bin/bash
#SBATCH -N 1
#SBATCH -c 38
#SBATCH --mem=100G
#SBATCH -t 24:00:00
#SBATCH --output=%x.e%j



#---------EDIT---------#
path=			# PATH/TO/Pre-Processing_HiC  ex. [/.../Signature-main/Pipeline/Pre-Processing_HiC]
list=			# full pathway and file name of a single column text file containing a list of your samples/datasets with unknown gonosomal sex
#----------------------#



#-------DONTEDIT-------#
module load R/4.2.1

PMR=$path/output
toscr=$path/scripts
train_f=$toscr/gonosome_train_data_f.txt
train_m=$toscr/gonosome_train_data_m.txt

# prepare data
echo "# Merging PMR files into correct format"
Rscript $toscr/gonosomes_PMR_script_merge.R $PMR $list
echo "Complete"

# run data through logistic regression model
echo "# Running Logistic Regression"
unknown_data=$PMR/gonosome_unknown_data.txt
Rscript $toscr/gonosome_prediction_LogisticRegression.R $PMR $train_f $train_m $unknown_data $list
echo "Complete"
#----------------------#
