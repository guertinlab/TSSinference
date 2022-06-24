cat U656SK2407_01_L00*/*_R1_*fastq.gz > HEK_ZNF143FKBP_DMSO_rep1_PE1.fastq.gz
cat U656SK2407_02_L00*/*_R1_*fastq.gz > HEK_ZNF143FKBP_DMSO_rep2_PE1.fastq.gz
cat U656SK2407_03_L00*/*_R1_*fastq.gz > HEK_ZNF143FKBP_DMSO_rep3_PE1.fastq.gz
cat U656SK2407_04_L00*/*_R1_*fastq.gz > HEK_ZNF143FKBP_DMSO_rep4_PE1.fastq.gz
cat U656SK2407_05_L00*/*_R1_*fastq.gz > HEK_DMSO_rep1_PE1.fastq.gz
cat U656SK2407_06_L00*/*_R1_*fastq.gz > HEK_DMSO_rep2_PE1.fastq.gz
cat U656SK2407_07_L00*/*_R1_*fastq.gz > HEK_DMSO_rep3_PE1.fastq.gz
cat U656SK2407_08_L00*/*_R1_*fastq.gz > HEK_DMSO_rep4_PE1.fastq.gz
cat U656SK2407_09_L00*/*_R1_*fastq.gz > HEK_ZNF143FKBP_dTAG13_rep1_PE1.fastq.gz
cat U656SK2407_10_L00*/*_R1_*fastq.gz > HEK_ZNF143FKBP_dTAG13_rep2_PE1.fastq.gz
cat U656SK2407_11_L00*/*_R1_*fastq.gz > HEK_ZNF143FKBP_dTAG13_rep3_PE1.fastq.gz
cat U656SK2407_12_L00*/*_R1_*fastq.gz > HEK_ZNF143FKBP_dTAG13_rep4_PE1.fastq.gz
cat U656SK2407_13_L00*/*_R1_*fastq.gz > HEK_dTAG13_rep1_PE1.fastq.gz
cat U656SK2407_14_L00*/*_R1_*fastq.gz > HEK_dTAG13_rep2_PE1.fastq.gz
cat U656SK2407_15_L00*/*_R1_*fastq.gz > HEK_dTAG13_rep3_PE1.fastq.gz
cat U656SK2407_16_L00*/*_R1_*fastq.gz > HEK_dTAG13_rep4_PE1.fastq.gz

cat U656SK2407_01_L00*/*_R2_*fastq.gz > HEK_ZNF143FKBP_DMSO_rep1_PE2.fastq.gz
cat U656SK2407_02_L00*/*_R2_*fastq.gz > HEK_ZNF143FKBP_DMSO_rep2_PE2.fastq.gz
cat U656SK2407_03_L00*/*_R2_*fastq.gz > HEK_ZNF143FKBP_DMSO_rep3_PE2.fastq.gz
cat U656SK2407_04_L00*/*_R2_*fastq.gz > HEK_ZNF143FKBP_DMSO_rep4_PE2.fastq.gz
cat U656SK2407_05_L00*/*_R2_*fastq.gz > HEK_DMSO_rep1_PE2.fastq.gz
cat U656SK2407_06_L00*/*_R2_*fastq.gz > HEK_DMSO_rep2_PE2.fastq.gz
cat U656SK2407_07_L00*/*_R2_*fastq.gz > HEK_DMSO_rep3_PE2.fastq.gz
cat U656SK2407_08_L00*/*_R2_*fastq.gz > HEK_DMSO_rep4_PE2.fastq.gz
cat U656SK2407_09_L00*/*_R2_*fastq.gz > HEK_ZNF143FKBP_dTAG13_rep1_PE2.fastq.gz
cat U656SK2407_10_L00*/*_R2_*fastq.gz > HEK_ZNF143FKBP_dTAG13_rep2_PE2.fastq.gz
cat U656SK2407_11_L00*/*_R2_*fastq.gz > HEK_ZNF143FKBP_dTAG13_rep3_PE2.fastq.gz
cat U656SK2407_12_L00*/*_R2_*fastq.gz > HEK_ZNF143FKBP_dTAG13_rep4_PE2.fastq.gz
cat U656SK2407_13_L00*/*_R2_*fastq.gz > HEK_dTAG13_rep1_PE2.fastq.gz
cat U656SK2407_14_L00*/*_R2_*fastq.gz > HEK_dTAG13_rep2_PE2.fastq.gz
cat U656SK2407_15_L00*/*_R2_*fastq.gz > HEK_dTAG13_rep3_PE2.fastq.gz
cat U656SK2407_16_L00*/*_R2_*fastq.gz > HEK_dTAG13_rep4_PE2.fastq.gz


for i in *PE1.fastq.gz
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F"_PE1.fastq.gz" '{print $1}')
    echo $name
    cutadapt -m 26 -a TGGAATTCTCGGGTGCCAAGG ${name}_PE1.fastq.gz | \
       fqdedup -i - -o - | \
       fastx_trimmer -f 9 -l 38 | \
       fastx_reverse_complement -z -o ${name}_PE1.processed.fastq.gz
    cutadapt -m 26 -a GATCGTCGGACTGTAGAACTCTGAAC ${name}_PE2.fastq.gz | \
       fastx_trimmer -t 8 | \
       fastx_reverse_complement -z -o ${name}_PE2.processed.fastq.gz
    gunzip ${name}_PE1.processed.fastq.gz
    gunzip ${name}_PE2.processed.fastq.gz
    #deduplicated based on ID, UMI deduplicated effectively PE1  
    fastq_pair -t 500000 ${name}_PE1.processed.fastq ${name}_PE2.processed.fastq
    mv ${name}_PE1.processed.fastq.paired.fq ${name}_PE1.processed.paired.fastq
    mv ${name}_PE2.processed.fastq.paired.fq ${name}_PE2.processed.paired.fastq
    gzip ${name}_PE1.processed.paired.fastq
    gzip ${name}_PE2.processed.paired.fastq
    rm ${name}_PE1.processed.fastq
    rm ${name}_PE2.processed.fastq
    #testing
#    note the --rf flag below
    bowtie2 -p 3 --rf -x hg38 -1 ${name}_PE1.processed.paired.fastq.gz -2 ${name}_PE2.processed.paired.fastq.gz  | \
       samtools view -b - | \
       samtools sort - -o ${name}_PE.sorted.bam
#this also gets rid of unmated PE2--which effectively uses PE1 UMI to deduplicate.    
    samtools view -h -f 2 ${name}_PE.sorted.bam | \
	awk '!($3!=$7 && $7!="=") || $1 ~ /^@/' | \
    	samtools view -Sb - > ${name}_PE.filtered.bam
    samtools view -b -f 0x0040 ${name}_PE.filtered.bam > ${name}_PE1.bam
    samtools view -b -f 0x0080 ${name}_PE.filtered.bam > ${name}_PE2.bam
    samtools view -bh -F 20 ${name}_PE1.bam > ${name}_PE1_pro_plus.bam
    samtools view -bh -f 0x10 ${name}_PE1.bam > ${name}_PE1_pro_minus.bam
    samtools view -bh -F 20 ${name}_PE2.bam | bamToFastq -i - -fq /dev/stdout | fastx_reverse_complement -z -o ${name}_PE2_pro_plus.fastq.gz
    samtools view -bh -f 0x10 ${name}_PE2.bam | bamToFastq -i - -fq /dev/stdout | fastx_reverse_complement -z -o ${name}_PE2_pro_minus.fastq.gz
    bowtie2 -p 3 -x hg38 -U ${name}_PE2_pro_minus.fastq.gz  | \
       samtools view -b - | \
       samtools sort - -o ${name}_PE2_pro_minus.bam
    bowtie2 -p 3 -x hg38 -U ${name}_PE2_pro_plus.fastq.gz  | \
       samtools view -b - | \
       samtools sort - -o ${name}_PE2_pro_plus.bam
#the plus and minus are actually inverted here, NOTE the name change!
    seqOutBias hg38.fa ${name}_PE2_pro_minus.bam --no-scale --bw=${name}_PE2_pro_plus.bigWig --read-size=30
#the plus and minus are actually inverted here, NOTE the name change!
    seqOutBias hg38.fa ${name}_PE2_pro_plus.bam --no-scale --bw=${name}_PE2_pro_minus.bigWig --read-size=30
#named properly    
    seqOutBias hg38.fa ${name}_PE1_pro_plus.bam --no-scale --bw=${name}_PE1_pro_plus.bigWig --tail-edge --read-size=30
    seqOutBias hg38.fa ${name}_PE1_pro_minus.bam  --no-scale --bw=${name}_PE1_pro_minus.bigWig --tail-edge --read-size=30
done

for i in *PE1*.bigWig
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F".bigWig" '{print $1}')
    strand=$(echo $i | awk -F"pro_" '{print $NF}' | awk -F".bigWig" '{print $1}')
    echo $name
    echo $strand
    touch temp.txt
    bigWigToBedGraph $i ${name}_combined.bg
    if [ "$strand" == "plus" ]
    then
        echo "track type=bedGraph name=${name} color=255,0,0 alwaysZero=on visibility=full" >> temp.txt
    fi
    if [ "$strand" == "minus" ]
    then
        echo "track type=bedGraph name=${name} color=0,0,255 alwaysZero=on visibility=full" >> temp.txt
    fi
    head -49999999 ${name}_combined.bg > ${name}_2combined.bg
    cat temp.txt ${name}_2combined.bg > ${name}_combined.bedGraph
    gzip ${name}_combined.bedGraph
    rm temp.txt
done

for i in *PE2*.bigWig
do
    name=$(echo $i | awk -F"/" '{print $NF}' | awk -F".bigWig" '{print $1}')
    strand=$(echo $i | awk -F"pro_" '{print $NF}' | awk -F".bigWig" '{print $1}')
    echo $name
    echo $strand
    touch temp.txt
    bigWigToBedGraph $i ${name}_combined.bg
    if [ "$strand" == "plus" ]
    then
        echo "track type=bedGraph name=${name}_TSS color=255,0,0 alwaysZero=on visibility=full" >> temp.txt
    fi
    if [ "$strand" == "minus" ]
    then
        echo "track type=bedGraph name=${name}_TSS color=0,0,255 alwaysZero=on visibility=full" >> temp.txt
    fi
    head -49999999 ${name}_combined.bg > ${name}_2combined.bg
    cat temp.txt ${name}_2combined.bg > ${name}_combined.bedGraph
    gzip ${name}_combined.bedGraph
    rm temp.txt
done


#load in browser
#to get the merged BAM
samtools merge - *_PE.filtered.bam | \
       samtools sort - -o H9_combined.sorted.bam

#to get the merged bigWigs
cellline=XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
pluspe1files=$(ls *_PE1_pro_plus.bam)
minuspe1files=$(ls *_PE1_pro_minus.bam)
pluspe2files=$(ls *_PE2_pro_plus.bam)
minuspe2files=$(ls *_PE2_pro_minus.bam)
#notice againt the 
seqOutBias hg38.fa ${minuspe2files} --no-scale --bw=${cellline}_PE2_combined_pro_plus.bigWig --read-size=30
seqOutBias hg38.fa ${pluspe2files} --no-scale --bw=${cellline}_PE2_combined_pro_minus.bigWig --read-size=30





#mappability

#this file is generated during a seqOutBias run: hg38.tal_30.gtTxt.gz
#code from https://github.com/andrelmartins/bigWig/tree/master/calc_mappability

SEQNAMES=`gzcat hg38.fa.gz | grep ">" | sed "s/^>//g"`
gzcat hg38.tal_30.gtTxt.gz | perl tallymer2bed.pl $SEQNAMES | bedops -m  - | gzip > hg38_30mers.unmap.bed.gz

wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit

function makeBigWig {
 twoBitInfo $TWOBIT chromInfo
 gzcat $BEDPREFIX.unmap.bed.gz | grep "^chr" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,1}' | sort -k1,1 -k2,2n > $BEDPREFIX.bedGraph
 bedGraphToBigWig $BEDPREFIX.bedGraph chromInfo $BEDPREFIX.bigWig
 rm $BEDPREFIX.bedGraph chromInfo
}

TWOBIT=hg38.2bit
BEDPREFIX=hg38_30mers
makeBigWig

