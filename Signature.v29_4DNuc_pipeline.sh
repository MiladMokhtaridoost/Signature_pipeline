module load bedtools/2.24.0 
module load R/4.2.1 
module load python/2.7.9


# reading the configuration file
echo -e "Loading configuration ...\n"
typeclass=$(echo $1 | cut -f2 -d"=")
specie=$(echo $2 | cut -f2 -d"=")
tissue=$(echo $3 | cut -f2 -d"=")
res=$(echo $4 | cut -f2 -d"=")
name=$(echo $5 | cut -f2 -d"=")
out=$(echo $6 | cut -f2 -d"=")
pval=$(echo $7 | cut -f2 -d"=")
dir=$(echo $8 | cut -f2 -d"=")
in=$(echo $9 | cut -f2 -d"=")
prtT=0


echo "Starting Signature ..."

case $tissue in
*.txt)

echo "-----------------------------------------------------------"
echo "You have selected a list of cells/tissues to run Signature:"
	cat $tissue
echo "-----------------------------------------------------------"
echo "-----------------------------------------------------------"
echo "-----------------------------------------------------------"

prtT=1
	
	while IFS= read -r cell ; do
	case $typeclass in 

	#------------------#
	# CIS INTERACTIONS #
	#------------------#
	cis)
	echo "You have selected CIS interaction analysis ..."	
		if [ -d $out/$cell.cis_tmp ] ; then
			echo "Output directory $out/$cell.cis_tmp exists"
		else
			echo "Making directory $out/$cell.cis_tmp ... "
			mkdir $out/$cell.cis_tmp
		fi
		
		awk -F"[_/]" '{print $2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9}' $dir/$specie/$cell/cis.$res\_iced.sorted.txt | awk -F" " '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$5-$2"\t"$7}'  > $out/$cell.cis_tmp/tmp
		cd $out/$cell.cis_tmp
		awk -F"\t" '{print>$1".bed"}' $out/$cell.cis_tmp/tmp 
		rm $out/$cell.cis_tmp/*chrM*	#removing chrM
		find `pwd` -name \*.bed -printf '%f\n' | sort | uniq > $out/$cell.cis_tmp/tmp.$cell.list.cis.txt 
		cis_all=$out/$cell.cis_tmp/tmp.$cell.list.cis.txt

	echo "-----------------------------------------------------------"
	echo -e "\n## calculating cis z-scores for $cell ..."
		echo "- - - - - -"
		echo "Starting R"
		echo "- - - - - -"
			Rscript $in/cis_LWLR.R $out/$cell.cis_tmp $out/$cell.cis_tmp $in/LWLR_train_cis.R $res
			# cat all chrom files into one df with all zscores (even non-sig pvalues)
			cat $out/$cell.cis_tmp/chr*.txt > $out/$cell.cis_tmp/tmp.zscores.cis.LOESS.txt
			mv $out/$cell.cis_tmp/tmp.zscores.cis.LOESS.txt $out/$cell.cis_tmp/tmp.cis.$cell.$name.$res\_LOESS.all.txt

	echo "-----------------------------------------------------------"
	echo "removing z-scores with p-values less than $pval ..."
	echo "-----------------------------------------------------------"
		awk -F"\t" -v p=$pval '{if ($(NF) <= p) print $0;}' $out/$cell.cis_tmp/tmp.cis.$cell.$name.$res\_LOESS.all.txt > $out/$cell.cis_tmp/tmp.zscores.cut.cis.LOESS.txt

	mv $out/$cell.cis_tmp/tmp.zscores.cut.cis.LOESS.txt $out/$cell.cis_tmp/tmp.zscores.cis.LOESS.txt
	mv $out/$cell.cis_tmp/tmp.zscores.cis.LOESS.txt $out/$cell.cis_tmp/tmp.cis.$cell.$name.$res\_LOESS.txt

	echo "format sig zscore df"
		cat $out/$cell.cis_tmp/tmp.cis.$cell.$name.$res\_LOESS.txt > $out/$cell.cis_tmp/$name.cis_all.$cell.$res.zscores.LOESS.txt 
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$12}' $out/$cell.cis_tmp/$name.cis_all.$cell.$res.zscores.LOESS.txt > $out/tmp.$name.cis_all.$cell.$res.zscores.LOESS.mergeR.txt	
		sed -i "1i ID \t $cell" $out/tmp.$name.cis_all.$cell.$res.zscores.LOESS.mergeR.txt	
	echo "format all zscore df"
		cat $out/$cell.cis_tmp/tmp.cis.$cell.$name.$res\_LOESS.all.txt > $out/$cell.cis_tmp/$name.cis_all.$cell.$res.zscores.LOESS_zscoreall.txt 
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$12}' $out/$cell.cis_tmp/$name.cis_all.$cell.$res.zscores.LOESS_zscoreall.txt > $out/tmp.$name.cis_all.$cell.$res.zscores.LOESS.mergeR_zscoreall.txt	
		sed -i "1i ID \t $cell" $out/tmp.$name.cis_all.$cell.$res.zscores.LOESS.mergeR_zscoreall.txt	
	echo "format sig pvals df"
		cat $out/$cell.cis_tmp/tmp.cis.$cell.$name.$res\_LOESS.txt > $out/$cell.cis_tmp/$name.cis_all.$cell.$res.zscores.LOESS_pvalue.txt 
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$13}' $out/$cell.cis_tmp/$name.cis_all.$cell.$res.zscores.LOESS_pvalue.txt > $out/tmp.$name.cis_all.$cell.$res.zscores.LOESS.mergeR_pvalue.txt	
		sed -i "1i ID \t $cell" $out/tmp.$name.cis_all.$cell.$res.zscores.LOESS.mergeR_pvalue.txt	
	echo "format all pvals df"
		cat $out/$cell.cis_tmp/tmp.cis.$cell.$name.$res\_LOESS.all.txt > $out/$cell.cis_tmp/$name.cis_all.$cell.$res.zscores.LOESS_pvalueAll.txt 
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$13}' $out/$cell.cis_tmp/$name.cis_all.$cell.$res.zscores.LOESS_pvalueAll.txt > $out/tmp.$name.cis_all.$cell.$res.zscores.LOESS.mergeR_pvalueAll.txt	
		sed -i "1i ID \t $cell" $out/tmp.$name.cis_all.$cell.$res.zscores.LOESS.mergeR_pvalueAll.txt	
	
	;;
	
	
	#-----------------------#
	# pairwise INTERACTIONS #
	#-----------------------#
	transpairwise)
	echo "You have selected genome-wide TRANS pairwise analysis ..."
		if [ -d $out/$cell.trans.tmp ] ; then
			echo "Output directory $out/$cell.trans.tmp exists ..."
		else
			echo "Making directory $out/$cell.trans.tmp ..."
			mkdir $out/$cell.trans.tmp
		fi
		
		awk -F"[_/]" '{print $2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9}' $dir/$specie/$cell/trans.$res\_iced.sorted.txt | awk -F" " '{print $1$4".txt""\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$5-$2"\t"$7}' > $out/tmp ; mv $out/tmp $out/$cell.trans.tmp/
		cd $out/$cell.trans.tmp
		awk -F"\t" '{print>$1}' $out/$cell.trans.tmp/tmp
		rm $out/$cell.trans.tmp/*chrM*.txt	#removing chrM
		cat $out/$cell.trans.tmp/chr*chr*.txt > tmp.chrTransAll.txt
		find `pwd` -name chr\*.txt -printf '%f\n' | cut -d"-" -f1 | sort | uniq > $out/$cell.trans.tmp/$cell.list.trans.txt
		list=$out/$cell.trans.tmp/$cell.list.trans.txt		
	
	echo "-----------------------------------------------------------"
	echo -e "\n## calculating trans pairwise z-scores for $cell ..."
		echo "- - - - - -"
		echo "Starting R"
		echo "- - - - - -"
			while IFS= read -r file ; do
			Rscript $in/trans.pairwise_LWLR.R $file $out/$cell.trans.tmp $in/LWLR_train.R $res
			done < "$list"
			# cat all chrom files into one df with all zscores (even non-sig pvalues)
			cat $out/$cell.trans.tmp/*trans.LOESS*.txt > $out/tmp.trans.$name.$cell.$res\_zscore_all.txt	 
	
	echo "-----------------------------------------------------------"
	echo "keeping lowest z-score per interaction ..."    # removing lowest absolute value zscore (duplicated because of loess on both anchor and target)
	echo "-----------------------------------------------------------"
		Rscript $in/lowest_zscore_pairwise.R $out/tmp.trans.$name.$cell.$res\_zscore_all.txt  
	
	echo "-----------------------------------------------------------"
	echo "removing z-scores with p-values less than $pval ..."
	echo "-----------------------------------------------------------"
		awk -F"\t" -v p=$pval '{if ($14 <= p) print $0;}' $out/tmp.trans.$name.$cell.$res\_zscore_all.txt.tmp.zscores.trans.lowest.score.txt > $out/tmp.trans.$name.$cell.$res\_zscore.txt
	
	echo "format sig zscore df"
		awk -F"\t" '{print "A"$3"."$4"."$5"B"$6"."$7"."$8"\t"$13}' $out/tmp.trans.$name.$cell.$res\_zscore.txt > $out/tmp.$name.trans_all.$cell.$res.zscores.zscore.mergeR.txt
		sed -i "1i ID \t $cell" $out/tmp.$name.trans_all.$cell.$res.zscores.zscore.mergeR.txt
		lim=$(echo $(awk -F"\t" '{print $10}' $out/tmp.trans.$name.$cell.$res\_zscore.txt | awk '($1 >= 1)' | sort -n | awk -v p=$percentile '{all[NR] = $0} END{print all[int(NR*p)]}'))
		awk -v x=$lim '($13 >= x)' $out/tmp.trans.$name.$cell.$res\_zscore.txt | cat -n > $out/tmp.$cell.percentile.positive.bed
		awk -v z=$cell -F" " '{print $10"\t"z}' $out/tmp.trans.$name.$cell.$res\_zscore.txt > $out/tmp.$cell.violin.txt 
	echo "format all zscore df"
		cp $out/tmp.trans.$name.$cell.$res\_zscore_all.txt.tmp.zscores.trans.lowest.score.txt $out/tmp.trans.$name.$cell.$res\_zscoreall.txt
		awk -F"\t" '{print "A"$3"."$4"."$5"B"$6"."$7"."$8"\t"$13}' $out/tmp.trans.$name.$cell.$res\_zscore_all.txt.tmp.zscores.trans.lowest.score.txt > $out/tmp.$name.trans_all.$cell.$res.zscores.mergeR_zscoreall.txt
		sed -i "1i ID \t $cell" $out/tmp.$name.trans_all.$cell.$res.zscores.mergeR_zscoreall.txt
		lim=$(echo $(awk -F"\t" '{print $10}' $out/tmp.trans.$name.$cell.$res\_zscoreall.txt | awk '($1 >= 1)' | sort -n | awk -v p=$percentile '{all[NR] = $0} END{print all[int(NR*p)]}'))
		awk -v x=$lim '($13 >= x)' $out/tmp.trans.$name.$cell.$res\_zscoreall.txt | cat -n > $out/tmp.$cell.percentile.positive_zscoreall.bed
		awk -v z=$cell -F" " '{print $10"\t"z}' $out/tmp.trans.$name.$cell.$res\_zscoreall.txt > $out/tmp.$cell.zscoreall.violin.txt 
	echo "format sig pvals df"
		cp $out/tmp.trans.$name.$cell.$res\_zscore.txt $out/tmp.trans.$name.$cell.$res\_pvalue.txt
		awk -F"\t" '{print "A"$3"."$4"."$5"B"$6"."$7"."$8"\t"$14}' $out/tmp.trans.$name.$cell.$res\_zscore.txt > $out/tmp.$name.trans_all.$cell.$res.zscores.mergeR_pvalue.txt
		sed -i "1i ID \t $cell" $out/tmp.$name.trans_all.$cell.$res.zscores.mergeR_pvalue.txt
		lim=$(echo $(awk -F"\t" '{print $10}' $out/tmp.trans.$name.$cell.$res\_pvalue.txt | awk '($1 >= 1)' | sort -n | awk -v p=$percentile '{all[NR] = $0} END{print all[int(NR*p)]}'))
		awk -v x=$lim '($13 >= x)' $out/tmp.trans.$name.$cell.$res\_pvalue.txt | cat -n > $out/tmp.$cell.percentile.positive_pvalue.bed
		awk -v z=$cell -F" " '{print $10"\t"z}' $out/tmp.trans.$name.$cell.$res\_pvalue.txt > $out/tmp.$cell.pvalue.violin.txt
	echo "format all pvals df"		
		cp $out/tmp.trans.$name.$cell.$res\_zscore_all.txt.tmp.zscores.trans.lowest.score.txt $out/tmp.trans.$name.$cell.$res\_pvalueAll.txt
		awk -F"\t" '{print "A"$3"."$4"."$5"B"$6"."$7"."$8"\t"$14}' $out/tmp.trans.$name.$cell.$res\_zscore_all.txt.tmp.zscores.trans.lowest.score.txt > $out/tmp.$name.trans_all.$cell.$res.zscores.mergeR_pvalueAll.txt
		sed -i "1i ID \t $cell" $out/tmp.$name.trans_all.$cell.$res.zscores.mergeR_pvalueAll.txt
		lim=$(echo $(awk -F"\t" '{print $10}' $out/tmp.trans.$name.$cell.$res\_pvalueAll.txt | awk '($1 >= 1)' | sort -n | awk -v p=$percentile '{all[NR] = $0} END{print all[int(NR*p)]}'))
		awk -v x=$lim '($13 >= x)' $out/tmp.trans.$name.$cell.$res\_pvalueAll.txt | cat -n > $out/tmp.$cell.percentile.positive_pvalueAll.bed
		awk -v z=$cell -F" " '{print $10"\t"z}' $out/tmp.trans.$name.$cell.$res\_pvalueAll.txt > $out/tmp.$cell.pvalueAll.violin.txt
		
	echo "-----------------------------------------------------------"
	echo "cleaning up ..."
		rm -r $out/chr*
        echo "-----------------------------------------------------------"
	
	;;
	
	
	#-----------------------#
	#  1vsAll INTERACTIONS  #
	#-----------------------#
	trans1vsAll)
	echo -e "\nYou have selected genome-wide trans-1vsAll analysis"
	
        if [ -d $out ] ; then
            echo "Now processing: $cell"
        else
            mkdir $out
	    echo "Now processing: $cell"
        fi

	echo "-------------------------------------------"	
		echo -e "\n## preparing data ..."
		awk -F"[_/]" '{print $2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9}' $dir/$specie/$cell/trans.$res\_iced.sorted.txt | awk -F" " '{print $1$4".txt""\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$5-$2"\t"$7}' > $out/tmp
		cd $out
		awk -F"\t" '{print>$1}' $out/tmp
		rm $out/*chrM*.txt	#removing chrM
		echo "$out/tmp"
	
	echo -e "\n-------------------------------------------"
	echo -e "\n## calculating trans one vs all z-scores ..."
		echo "- - - - - -"
		echo "Starting R"
		echo "- - - - - -"
			Rscript $in/trans.1vsAll_LWLR.R $out $out $in/LWLR_train.R $res
		echo "- - - - - - - - - - - - - - - - - - - - - -"
		echo "Complete: trans.1vsAll_LWLR.R"
		echo "Aggregated files produced: $(find $out/*aggregated.txt | wc -l)"
		echo "- - - - - - - - - - - - - - - - - - - - - -"
		# cat all chrom files into one df with all zscores (even non-sig pvalues)
		cat $out/*aggregated.txt > $out/tmp.trans.$name.$cell.$res\_zscore.all.txt
	
	echo -e "\n-------------------------------------------"
	echo -e "\n## keeping lowest z-score per interaction ..."    # removing lowest absolute value zscore (duplicated because of loess on both anchor and target)
		echo "- - - - - -"
		echo "Starting R"
		echo "- - - - - -"
			Rscript $in/lowest_zscore_1vsAll.R $out/tmp.trans.$name.$cell.$res\_zscore.all.txt
		echo "- - - - - - - - - - - - - - - - -"
		echo "Complete: lowest_zscore_1vsAll.R"
		echo "- - - - - - - - - - - - - - - - -"
		
	echo -e "\n-------------------------------------------"
	echo -e "\n## removing z-scores with p-values less than $pval ..."

		awk -F"\t" -v p=$pval '{if ($12 <= p) print $0;}' $out/tmp.trans.$name.$cell.$res\_zscore.all.txt.tmp.zscores.trans.lowest.score.txt > $out/tmp.trans.$name.$cell.$res\_zscore.cut.txt
		mv $out/tmp.trans.$name.$cell.$res\_zscore.cut.txt $out/tmp.trans.$name.$cell.$res\_zscore.txt
	
	echo -e "\n-------------------------------------------"
	echo -e "\n## reformatting individual data frames ..."

	echo "> sig zscore df"
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$11}' $out/tmp.trans.$name.$cell.$res\_zscore.txt > $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR.txt
		sed -i "1i ID \t $cell" $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR.txt
	echo "> all zscore df"
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$11}' $out/tmp.trans.$name.$cell.$res\_zscore.all.txt.tmp.zscores.trans.lowest.score.txt > $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_zscoreall.txt
		sed -i "1i ID \t $cell" $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_zscoreall.txt

	echo "> sig pvals df"
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$12}' $out/tmp.trans.$name.$cell.$res\_zscore.txt > $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_pvalue.txt
		sed -i "1i ID \t $cell" $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_pvalue.txt
	echo "> all pvals df"
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$12}' $out/tmp.trans.$name.$cell.$res\_zscore.all.txt.tmp.zscores.trans.lowest.score.txt > $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_pvalueAll.txt
		sed -i "1i ID \t $cell" $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_pvalueAll.txt

	echo "> sig qvals df"
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$14}' $out/tmp.trans.$name.$cell.$res\_zscore.txt > $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_qvalue.txt
		sed -i "1i ID \t $cell" $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_qvalue.txt
	echo "> all qvals df"
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$14}' $out/tmp.trans.$name.$cell.$res\_zscore.all.txt.tmp.zscores.trans.lowest.score.txt > $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_qvalueAll.txt
		sed -i "1i ID \t $cell" $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_qvalueAll.txt
		
	echo -e "\n-------------------------------------------"
	echo -e "\n## cleaning up ..."
		mkdir -p $out/delete.$cell
		echo "$out/delete.$cell"
		mv $out/chr* $out/delete.$cell
	echo -e "\n-----------------------------------------------------------"
	echo -e "-----------------------------------------------------------"
	;;
	
	
	#----------------------------#
	#  1vsAll (NS) INTERACTIONS  #
	#----------------------------#
	trans-1vsAll-noSex)
		if [ -d $out ] ; then
			echo "Output directory $out exists ..."
		else
			echo "Making directory $out ..."
			mkdir $out
		fi
		
		awk -F"[_/]" '{print $2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9}' $dir/$specie/$cell/trans.$res\_iced.sorted.txt | awk -F" " '{print $1$4".txt""\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$5-$2"\t"$7}' > $out/tmp
		cd $out
		awk -F"\t" '{print>$1}' $out/tmp
		rm $out/*chrM*.txt	#removing chrM
		rm $out/*chrX*.txt	#removing chrX
		rm $out/*chrY*.txt	#removing chrY
		
	echo "-----------------------------------------------------------"
	echo -e "\n## calculating trans one vs all z-scores (sex chromosomes removed) for $cell ..."
		echo "- - - - - -"
		echo "Starting R"
		echo "- - - - - -"
			Rscript $in/1vsAll.loess.trans.R $out $out $in/LWLR_train.R $res
			# cat all chrom files into one df with all zscores (even non-sig pvalues)
			cat $out/*aggregated.txt > $out/tmp.trans.$name.$cell.$res\_zscore.all.txt	
			rm $out/*aggregated.txt

	echo "-----------------------------------------------------------"
	echo "keeping lowest z-score per interaction ..."    # removing lowest absolute value zscore (duplicated because of loess on both anchor and target)
	echo "-----------------------------------------------------------"
		Rscript $in/lowest_zscore_1vsAll.R $out/tmp.trans.$name.$cell.$res\_zscore.all.txt  

	echo "-----------------------------------------------------------"
	echo "removing z-scores with p-values less than $pval ..."
	echo "-----------------------------------------------------------"
		awk -F"\t" -v p=$pval '{if ($13 <= p) print $0;}' $out/tmp.trans.$name.$cell.$res\_zscore.all.txt.tmp.zscores.trans.lowest.score.txt > $out/tmp.trans.$name.$cell.$res\_zscore.cut.txt
		mv $out/tmp.trans.$name.$cell.$res\_zscore.cut.txt $out/tmp.trans.$name.$cell.$res\_zscore.txt
		
	Rscript $in/floess.trans.R $out/tmp
    
	echo "format sig zscore df"
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$12}' $out/tmp.trans.$name.$cell.$res\_zscore.txt > $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR.txt
		sed -i "1i ID \t $cell" $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR.txt
		awk -v z=$cell -F" " '{print $10"\t"z}' $out/tmp.trans.$name.$cell.$res\_zscore.txt > $out/tmp.$cell.violin.txt
	echo "format all zscore df"
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$12}' $out/tmp.trans.$name.$cell.$res\_zscore.all.txt.tmp.zscores.trans.lowest.score.txt > $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_zscoreall.txt
		sed -i "1i ID \t $cell" $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_zscoreall.txt
		awk -v z=$cell -F" " '{print $10"\t"z}' $out/tmp.trans.$name.$cell.$res\_zscore.all.txt > $out/tmp.$cell.zscoreall.violin.txt
	echo "format sig pvals df"
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$13}' $out/tmp.trans.$name.$cell.$res\_zscore.txt > $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_pvalue.txt
		sed -i "1i ID \t $cell" $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_pvalue.txt
		awk -v z=$cell -F" " '{print $10"\t"z}' $out/tmp.trans.$name.$cell.$res\_zscore.txt > $out/tmp.$cell.pvalue.violin.txt
	echo "format all pvals df"
		awk -F"\t" '{print "A"$2"."$3"."$4"B"$5"."$6"."$7"\t"$13}' $out/tmp.trans.$name.$cell.$res\_zscore.all.txt.tmp.zscores.trans.lowest.score.txt > $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_pvalueAll.txt
		sed -i "1i ID \t $cell" $out/$name.trans_all.$cell.$res.zscores.zscore.mergeR_pvalueAll.txt
		awk -v z=$cell -F" " '{print $10"\t"z}' $out/tmp.trans.$name.$cell.$res\_zscore.all.txt > $out/tmp.$cell.pvalueAll.violin.txt
	
	echo "-----------------------------------------------------------"
	echo "cleaning up ..."
		rm $out/chr*
	echo "-----------------------------------------------------------"	
	;;
	
	esac
	done < "$tissue"
	


echo -e "\nAll cells are finished processing"
echo "-----------------------------------------------------------"
cd $out


echo "Merging cells together ..."
	#file formatting for sig zscores
	$in/merger.awk $out/*.mergeR.txt > $out/tmp.totals.txt

	#file formatting for ALL zscores
	$in/merger.awk $out/*.mergeR_zscoreall.txt > $out/tmp.totals_zscoreall.txt

	#file formatting for sig pvalues
	$in/merger.awk $out/*.mergeR_pvalue.txt > $out/tmp.totals_pvalue.txt

	#file formatting for ALL pvalues
	$in/merger.awk $out/*.mergeR_pvalueAll.txt > $out/tmp.totals_pvalueAll.txt
	
	#file formatting for sig pvalues
	$in/merger.awk $out/*.mergeR_qvalue.txt > $out/tmp.totals_qvalue.txt

	#file formatting for ALL pvalues
	$in/merger.awk $out/*.mergeR_qvalueAll.txt > $out/tmp.totals_qvalueAll.txt
echo "-----------------------------------------------------------"


echo "Writting final output files ..."
	#formatting sig zscores files
	awk -F"\t" '{print $1}' $out/tmp.totals.txt > $out/tmp.id.txt
	awk -F"\t" '{print $1}' $out/tmp.totals.txt | awk -F"[/_]" '{print $4"\t"$5"\t"$6}' > $out/tmp.targets.txt
	paste $out/tmp.targets.txt $out/tmp.id.txt | sort -dsk1,1 -k2n,2 -k3nr,3 | grep -v 'ID' > $out/tmp.togeneanno.txt
	$in/merger.awk $out/*.mergeR.txt > $out/tmp2.totals.txt
	awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1'  $out/tmp2.totals.txt >  $out/$name.zscores.txt
	echo "> $out/$name.zscores.txt"
	#formatting ALL zscores files
	awk -F"\t" '{print $1}' $out/tmp.totals_zscoreall.txt > $out/tmp.id_zscoresall.txt
	awk -F"\t" '{print $1}' $out/tmp.totals_zscoreall.txt | awk -F"[/_]" '{print $4"\t"$5"\t"$6}' > $out/tmp.targets_zscoreall.txt
	paste $out/tmp.targets_zscoreall.txt $out/tmp.id_zscoresall.txt | sort -dsk1,1 -k2n,2 -k3nr,3 | grep -v 'ID' > $out/tmp.togeneanno_zscoreall.txt
	$in/merger.awk $out/*.mergeR_zscoreall.txt > $out/tmp2_zscoreall.totals_zscoreall.txt
	awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1'  $out/tmp2_zscoreall.totals_zscoreall.txt >  $out/$name.zscoresall.txt
	echo "> $out/$name.zscoresall.txt"
	
	#formatting sig pvalue files
	awk -F"\t" '{print $1}' $out/tmp.totals_pvalue.txt > $out/tmp.id_pvalue.txt
	awk -F"\t" '{print $1}' $out/tmp.totals_pvalue.txt | awk -F"[/_]" '{print $4"\t"$5"\t"$6}' > $out/tmp.targets_pvalue.txt
	paste $out/tmp.targets_pvalue.txt $out/tmp.id_pvalue.txt | sort -dsk1,1 -k2n,2 -k3nr,3 | grep -v 'ID' > $out/tmp.togeneanno_pvalue.txt
	$in/merger.awk $out/*.mergeR_pvalue.txt > $out/tmp2_pvalue.totals_pvalue.txt
	awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1'  $out/tmp2_pvalue.totals_pvalue.txt >  $out/$name.pvalue.txt
	echo "> $out/$name.pvalue.txt"
	#formatting ALL pvalue files
	awk -F"\t" '{print $1}' $out/tmp.totals_pvalueAll.txt > $out/tmp.id_pvalueAll.txt
	awk -F"\t" '{print $1}' $out/tmp.totals_pvalueAll.txt | awk -F"[/_]" '{print $4"\t"$5"\t"$6}' > $out/tmp.targets_pvalueAll.txt
	paste $out/tmp.targets_pvalueAll.txt $out/tmp.id_pvalueAll.txt | sort -dsk1,1 -k2n,2 -k3nr,3 | grep -v 'ID' > $out/tmp.togeneanno_pvalueAll.txt
	$in/merger.awk $out/*.mergeR_pvalueAll.txt > $out/tmp2_pvalueAll.totals_pvalueAll.txt
	awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1'  $out/tmp2_pvalueAll.totals_pvalueAll.txt >  $out/$name.pvalueAll.txt
	echo "> $out/$name.pvalueAll.txt"
	
	#formatting sig qvalue files
	awk -F"\t" '{print $1}' $out/tmp.totals_qvalue.txt > $out/tmp.id_qvalue.txt
	awk -F"\t" '{print $1}' $out/tmp.totals_qvalue.txt | awk -F"[/_]" '{print $4"\t"$5"\t"$6}' > $out/tmp.targets_qvalue.txt
	paste $out/tmp.targets_qvalue.txt $out/tmp.id_qvalue.txt | sort -dsk1,1 -k2n,2 -k3nr,3 | grep -v 'ID' > $out/tmp.togeneanno_qvalue.txt
	$in/merger.awk $out/*.mergeR_qvalue.txt > $out/tmp2_qvalue.totals_qvalue.txt
	awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1'  $out/tmp2_qvalue.totals_qvalue.txt >  $out/$name.qvalue.txt
	echo "> $out/$name.qvalue.txt"
	#formatting ALL qvalue files
	awk -F"\t" '{print $1}' $out/tmp.totals_qvalueAll.txt > $out/tmp.id_qvalueAll.txt
	awk -F"\t" '{print $1}' $out/tmp.totals_qvalueAll.txt | awk -F"[/_]" '{print $4"\t"$5"\t"$6}' > $out/tmp.targets_qvalueAll.txt
	paste $out/tmp.targets_qvalueAll.txt $out/tmp.id_qvalueAll.txt | sort -dsk1,1 -k2n,2 -k3nr,3 | grep -v 'ID' > $out/tmp.togeneanno_qvalueAll.txt
	$in/merger.awk $out/*.mergeR_qvalueAll.txt > $out/tmp2_qvalueAll.totals_qvalueAll.txt
	awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1'  $out/tmp2_qvalueAll.totals_qvalueAll.txt >  $out/$name.qvalueAll.txt
	echo "> $out/$name.qvalueAll.txt"	
echo "-----------------------------------------------------------"


#final cleang up 
rm -r $out/*tmp*
rm $out/0
rm $out/*mergeR*


echo "Signature is COMPLETE!"
echo "-----------------------------------------------------------"

;;
*) 
;;

esac


# futher subsetting

for celltype in $(cat $tissue);do

	echo "Processing pvalues for $celltype"
	
	echo "# separate dataset from batch" #col number is same for zscore file and pvalue file which is why it's only done once
		head -1 $out/$name.zscores.txt > $out/tmp.colnum
		COLNUM="$(awk -v val=$celltype '{for (i=1; i<=NF; i++) if ($i==val) {print i} }' $out/tmp.colnum)"
		awk -v coln=$COLNUM '{print $1,$coln}' $out/$name.zscores.txt > $out/tmp.z
		awk -v coln=$COLNUM '{print $coln}' $out/$name.pvalue.txt > $out/tmp.p
	echo "# bind zscore and pvalue"
		paste $out/tmp.z $out/tmp.p > $out/tmp.zp
	echo "# filter by positive zscore"
		awk '{if ($2 > 0) print $0;}' $out/tmp.zp > $out/tmp.zp.pos
	echo "# filter by pvalue < 0.5"
		awk -v sl=$pval '{if ($3 <= sl) print $0;}' $out/tmp.zp.pos > $out/tmp.pos.sig
	echo "# subset pvalues"
		awk '{print $1, $3}' $out/tmp.pos.sig > $out/$celltype.trans1vsAll.1MB.pvalue_pos.txt
		sed -i "1i ID \t $celltype" $out/$celltype.trans1vsAll.1MB.pvalue_pos.txt
		filesize=$(ls -l $out | grep $celltype.trans1vsAll.1MB.pvalue_pos.txt | awk '{print $5}')
		if [ $filesize == 0 ] ; then echo -e "ID \t $celltype" > $out/$celltype.trans1vsAll.1MB.pvalue_pos.txt ; fi
	echo "# filter by negative zscore"
		awk '{if ($2 < 0) print $0;}' $out/tmp.zp > $out/tmp.zp.neg
	echo "# filter by pvalue < 0.5"
		awk -v sl=$pval '{if ($3 <= sl) print $0;}' $out/tmp.zp.neg > $out/tmp.neg.sig
	echo "# subset pvalues"
		awk '{print $1, $3}' $out/tmp.neg.sig > $out/$celltype.trans1vsAll.1MB.pvalue_neg.txt
		sed -i "1i ID \t $celltype" $out/$celltype.trans1vsAll.1MB.pvalue_neg.txt
		filesize=$(ls -l $out | grep $celltype.trans1vsAll.1MB.pvalue_neg.txt | awk '{print $5}')
		if [ $filesize == 0 ] ; then echo -e "ID \t $celltype" > $out/$celltype.trans1vsAll.1MB.pvalue_neg.txt; fi
	echo "# clean up"
		rm $out/*tmp*
	echo "- - -"


	echo "Processing qvalues for $celltype"
	
	echo "# separate dataset from batch" #col number is same for zscore file and qvalue file which is why it's only done once
		head -1 $out/$name.zscores.txt > $out/tmp.colnum
		COLNUM="$(awk -v val=$celltype '{for (i=1; i<=NF; i++) if ($i==val) {print i} }' $out/tmp.colnum)"
		awk -v coln=$COLNUM '{print $1,$coln}' $out/$name.zscores.txt > $out/tmp.z
		awk -v coln=$COLNUM '{print $coln}' $out/$name.qvalue.txt > $out/tmp.q
	echo "# bind zscore and qvalue"
		paste $out/tmp.z $out/tmp.q > $out/tmp.zq
	echo "# filter by positive zscore"
		awk '{if ($2 > 0) print $0;}' $out/tmp.zq > $out/tmp.zq.pos
	echo "# filter by qvalue < 0.5"
		awk -v sl=$pval '{if ($3 <= sl) print $0;}' $out/tmp.zq.pos > $out/tmp.pos.sig
	echo "# subset qvalues"
		awk '{print $1, $3}' $out/tmp.pos.sig > $out/$celltype.trans1vsAll.1MB.qvalue_pos.txt
		sed -i "1i ID \t $celltype" $out/$celltype.trans1vsAll.1MB.qvalue_pos.txt
		filesize=$(ls -l $out | grep $celltype.trans1vsAll.1MB.qvalue_pos.txt | awk '{print $5}')
		if [ $filesize == 0 ] ; then echo -e "ID \t $celltype" > $out/$celltype.trans1vsAll.1MB.qvalue_pos.txt ; fi
	echo "# filter by negative zscore"
		awk '{if ($2 < 0) print $0;}' $out/tmp.zq > $out/tmp.zq.neg
	echo "# filter by qvalue < 0.5"
		awk -v sl=$pval '{if ($3 <= sl) print $0;}' $out/tmp.zq.neg > $out/tmp.neg.sig
	echo "# subset qvalues"
		awk '{print $1, $3}' $out/tmp.neg.sig > $out/$celltype.trans1vsAll.1MB.qvalue_neg.txt
		sed -i "1i ID \t $celltype" $out/$celltype.trans1vsAll.1MB.qvalue_neg.txt
		filesize=$(ls -l $out | grep $celltype.trans1vsAll.1MB.qvalue_neg.txt | awk '{print $5}')
		if [ $filesize == 0 ] ; then echo -e "ID \t $celltype" > $out/$celltype.trans1vsAll.1MB.qvalue_neg.txt; fi
	echo "# clean up"
		rm $out/*tmp*
	echo "- - -"


done
