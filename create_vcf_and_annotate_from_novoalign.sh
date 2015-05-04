#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe threaded 2




mkdir /$1_$2/
temp=$JOB_ID

mkdir $TMPDIR/$temp

cd $TMPDIR/$temp

samtools mpileup -B -l gencode_CDS_plus_10bp_flat.bed -C 60 -f human_g1k_v37.fasta  -F 0.05  $1_novoalign.bam  $2_novoalign.bam |\
java -Xmx20g -jar VarScan.jar somatic --output-indel /$1_$2/$1_$2_varscan.indel --output-snp /$1_$2/$1_$2_varscan.snp --mpileup 1 --tumor-purity 1 --min-coverage 20 --min-var-freq 0.01 --somatic-p-value 0.00001

# reform files for annovar

awk 'BEGIN{OFS=FS="\t"}{if ((($5+$6)>=20) && (($6/($5+$6))<=0.05) && (($9+$10)>=20) && (($10/($9+$10))>=0.15)) print $1,$2,$2,$3,$4,"somatic",$7,$11,$5,$6,$9,$10}' /home/nowinski/varscan/combined_calling/$1_$2/$1_$2_varscan.snp > somatic_snp.calls

grep -v Unknown /$1_$2/$1_$2_varscan.indel | awk 'BEGIN{OFS=FS="\t"}{if ((substr($4,1,1)=="-") && (($5+$6)>=20) && (($6/($5+$6))<=0.05) && (($9+$10)>=20) && ($10/($9+$10)>=0.15)) print $1,$2,$2+(length($4)-2),substr($4,2),"-","somatic",$7,$11,$5,$6,$9,$10;
else if ((substr($4,1,1)=="+") && (($5+$6)>=20) && (($6/($5+$6))<=0.05) && (($9+$10)>=20) && ($10/($9+$10)>=0.15)) print $1,$2,$2,"-",substr($4,2),"somatic",$7,$11,$5,$6,$9,$10;}' > somatic_indel.calls

awk 'BEGIN{OFS=FS="\t"}{if ((($5+$6)>=20) && (($6/($5+$6)>0.15) && ($6/($5+$6)<0.85)) && (($9+$10)>=20) && ((($10/($9+$10))<=0.05) || (($10/($9+$10))>=0.95))) print $1,$2,$2,$3,$4,"LOH",$7,$11,$5,$6,$9,$10}' /$1_$2/$1_$2_varscan.snp > LOH_snp.calls

grep -v Unknown /$1_$2/$1_$2_varscan.indel | awk 'BEGIN{OFS=FS="\t"}{if ((substr($4,1,1)=="-") && (($5+$6)>=20) && (($6/($5+$6)>0.15) && ($6/($5+$6)<0.85)) && (($9+$10)>=20) && ((($10/($9+$10))<=0.05) || (($10/($9+$10))>=0.95))) print $1,$2,$2+(length($4)-2),substr($4,2),"-","LOH",$7,$11,$5,$6,$9,$10;
else if ((substr($4,1,1)=="+") && (($5+$6)>=20) && (($6/($5+$6)>0.15) && ($6/($5+$6)<0.85)) && (($9+$10)>=20) && ((($10/($9+$10))<=0.05) || (($10/($9+$10))>=0.95))) print $1,$2,$2,"-",substr($4,2),"LOH",$7,$11,$5,$6,$9,$10;
}' > LOH_indel.calls

cat somatic_snp.calls somatic_indel.calls LOH_snp.calls LOH_indel.calls > $1_$2.annovar

# annotate with respect to genes

annotate_variation.pl --buildver hg19 --splicing_threshold 10 -geneanno $1_$2.annovar /path...to...hg19/

awk 'BEGIN {FS=OFS="\t"} {if ($1 == "exonic;splicing") print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,substr($1,1,index($1,";")-1),substr($2,1,index($2,";")-1) ;                        
                     else if (($1 == "splicing" || $1 == "ncRNA_splicing") && index($2,"(") > 0 ) print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1,substr($2,1,index($2,"(")-1) ;    
                     else print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$1,$2 }' $1_$2.annovar.variant_function > $1_$2_1.annovar      

annotate_variation.pl --buildver hg19 --splicing_threshold 10 -geneanno $1_$2_1.annovar /path...to...hg19/
awk 'BEGIN {FS=OFS="\t"} {if ((substr($1,0,6) != "exonic") && ($1 != "splicing")) print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,".";    
                     else if ($1 == "splicing") print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$16":"substr($2,index($2,"(")+1,index($2,")")-index($2,"(")-1)}' $1_$2_1.annovar.variant_function > $1_$2_2.annovar
awk 'BEGIN {FS=OFS="\t"} {print $4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$2,$3}' $1_$2_1.annovar.exonic_variant_function >> $1_$2_2.annovar  

# cross reference with dbSNP141
/home/exome/bin/annotate_variation.pl --buildver hg19 --filter  -dbtype generic -genericdbfile dbSNP134a.annovar $1_$2_2.annovar /path...to...dbSNP141/
awk 'BEGIN {FS=OFS="\t"} {print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$2}' $1_$2_2.annovar.hg19_generic_dropped > $1_$2_3.annovar
awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,"."}' $1_$2_2.annovar.hg19_generic_filtered >> $1_$2_3.annovar

# cross reference with 1000g
/home/exome/bin/annotate_variation.pl --buildver hg19 -filter -dbtype generic -genericdbfile 1KG_2011_10.annovar $1_$2_3.annovar /path...to...1000g/

# rejig format back to input format
awk 'BEGIN {FS=OFS="\t"} {print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$2}' $1_$2_3.annovar.hg19_generic_dropped > $1_$2_4.annovar
awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,"."}' $1_$2_3.annovar.hg19_generic_filtered >> $1_$2_4.annovar