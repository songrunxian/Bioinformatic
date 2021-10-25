#!/bin/bash
#

genomeindex=/home/agrogene_lq/srx/ref/index
gtf=/home/agrogene_lq/srx/ref/Oryza_sativa.IRGSP-1.0.51.gtf
output_utrs=/home/agrogene_lq/srx/ref/qapatry/output_utrs.bed
output_sequences=/home/agrogene_lq/srx/ref/qapatry/output_sequences.fa
ensembl_identifiers=/home/agrogene_lq/srx/ref/qapatry/ensembl_identifiers.txt
rmats=/home/agrogene_lq/srx/tools/rmats_turbo_v4_1_1/rmats.py

for file in `grep "RX" <(ls ${pwd})`
do
	cd $file
		#进入每个SRX提取sample name
		ls *.gz > fq_tmp1_list
		awk 'BEGIN{FS="."}{print $1}' fq_tmp1_list > fq_tmp2_list 
		awk 'BEGIN{FS="_"}{print $1}' fq_tmp2_list > fq_tmp3_list
		sort fq_tmp3_list | uniq > rfq_list
		mv rfq_list fq_list	
		rm fq_tmp*
		#STAR
		p=`pwd`/
		for sample in `cat fq_list`
		do
	#		STAR --genomeDir $genomeindex --sjdbGTFfile $gtf --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMattributes MD NH --runThreadN 12 --readFilesCommand zcat --readFilesIn $sample* --outFileNamePrefix $p$sample
	#		genomeCoverageBed -bg -ibam $sample"Aligned.sortedByCoord.out.bam" -split > $sample"Aligned.sortedByCoord.out.wig"
	#		samtools view -c -F 4 $sample"Aligned.sortedByCoord.out.bam" > $sample"Aligned.sortedByCoord.out.bam.rsl"
		        echo -e "`ls $p$sample*.bam`\t`cat $p$sample*.rsl`" >> mapping_bam_location_with_depth.txt	
		done	
		
		#featureCounts
			featureCounts -T 10 -t exon -g gene_id -a $gtf -o $file.featureCounts.txt *.bam
		
		#Dapars2
			ls *.wig > wig_tmp_list
			awk '{print $1","}' wig_tmp_list > wig_tmp2_list
			tac <(tac wig_tmp2_list | sed "1s/,//g") > wig_tmp3_list
			awk '{printf $0}' wig_tmp3_list > wig_list
			rm wig_tmp*
			echo -e "Annotated_3UTR=/home/agrogene_lq/srx/DaPars_2/rice_refseq_extracted_3UTR.bed\n\nAligned_Wig_files=`cat wig_list`\n\nOutput_directory=/lineresult/$file\n\nOutput_result_file=chromosome\n\nCoverage_threshold=1\n\nNum_Threads = 20\n\nsequencing_depth_file=mapping_bam_location_with_depth.txt" > Dapars2_configure_file
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file 1
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file 2
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file 3
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file 4
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file 5
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file 6
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file 7
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file 8
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file 9
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file 10
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file 11
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file 12
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file Mt
			python2 /home/agrogene_lq/srx/tools/DaPars2/src/Dapars2_Multi_Sample.py Dapars2_configure_file Pt
		
		#QAPA
		source /home/agrogene_lq/software/anaconda3/etc/profile.d/conda.sh
		conda activate qapa	
		for sample in `cat fq_list`
		do
			grep $sample"_1.fastq.gz" <(ls)
			if [ $? -eq 0 ];then
				gunzip $sample"_1.fastq.gz"
				gunzip $sample"_2.fastq.gz"
				salmon quant -i /home/agrogene_lq/srx/ref/qapatry/utr_library -l IU -1 $sample"_1.fastq" -2 $sample"_2.fastq" -o $sample"_transcripts_quant"
				gzip $sample"_1.fastq"
				gzip $sample"_2.fastq"
				conda info --envs	
			else
				gunzip $sample".fastq.gz"
				salmon quant -i /home/agrogene_lq/srx/ref/qapatry/utr_library -l SF -r $sample".fastq" -o $sample"_transcripts_quant"
				gzip $sample".fastq"
				conda info --envs
			fi
		done
		qapa quant --db $ensembl_identifiers *_transcripts_quant/quant.sf > pau_results.txt
		conda info --envs
		conda deactivate
		conda info --envs
		
		#rmats
		mkdir result
		mkdir tmp
		ls *.bam > bam_list
		for sample in `cat bam_list`
		do
				echo $p$sample"," >> b1_tmp.txt					
		done
		touch b2.txt
		tac <(tac b1_tmp.txt | sed "1s/,//g") > b1_tmp1.txt
		awk '{printf $0}' b1_tmp1.txt > b1.txt
		rm b1_tmp*

		source /home/agrogene_lq/software/anaconda3/etc/profile.d/conda.sh
		conda activate /home/agrogene_lq/srx/tools/rmats_turbo_v4_1_1/conda_envs/rmats
			grep $sample"_1.fastq.gz" <(ls)
			if [ $? -eq 0 ];then
				python $rmats --b1 b1.txt --b2 b2.txt --gtf $gtf --tmp ./tmp/  --od ./result/ -t paired --nthread 10 --cstat 0.0001 --readLength 150
			else
				python $rmats --b1 b1.txt --b2 b2.txt --gtf $gtf --tmp ./tmp/  --od ./result/ -t single --nthread 10 --cstat 0.0001 --readLength 150
			fi
		conda deactivate
	cd ..
done
