#!/bin/bash
#SBATCH -N 1
#SBATCH -c 38
#SBATCH --mem=100G
#SBATCH -t 24:00:00
#SBATCH --output=%x.e%j



#--------------EDIT-------------------#
prefix=			# SAMPLE_NAME   ex. [Astrocyte_Spine]
inpath=			# PATH/TO/SAMFILES/FROM/BWA/RESULTS
toscr=			# PATH/TO/Pre-Processing_HiC  ex. [/.../Signature-main/Pipeline/Pre-Processing_HiC]
chrom_len_file=		# hg38.chrom.sizes from chosen genome assembly
chrList=		# text file containing a list of chromosomes in order from smallest to largest (length)
#-------------------------------------#



#------------DONTEDIT-----------------#
module load python/3.8.1
module load perl/5.28.0
i=0

outpath=$toscr/output
mkdir -p $outpath

cd $inpath


# produce PMR files #

for mapped_reads_file in $(ls *.sam); do

        #set index
        i=$(expr $i + 1)

        #set name variables
        outfile=$prefix.run$i.ProportionMappedReads.txt
        tmpfile=tmp.run$i

        #read in each line of chromosome file
        true > $outpath/$tmpfile.num.mapped.per.chrom.txt

        while IFS= read -r chrom
        do
        echo $chrom >> $outpath/$tmpfile.num.mapped.per.chrom.txt
        # grab total number of mapped reads per chrom, save to file
        grep -w "$chrom" $mapped_reads_file | wc -l >> $outpath/$tmpfile.num.mapped.per.chrom.txt
        done < "$chrList"

        # calculate mapped read coverage per chromosome
        python3 $toscr/scripts/gonosome_identity_v2.py $chrom_len_file $outpath/$tmpfile.num.mapped.per.chrom.txt > $outpath/$outfile

        #convert spaces to tabs in output file
        perl -p -i -e "s/ /\t/g" $outpath/$outfile

done

#-------------------------------------#
