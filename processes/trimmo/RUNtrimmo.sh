sampleID=$1
infolder=$2


# Find left and right reads in input folder
R1=$(find $infolder -name *R1_001.fastq.gz)
R2=$(find $infolder -name *R2_001.fastq.gz)


java -jar /local/genome/packages/trimmomatic/0.32/trimmomatic-0.32.jar PE -threads 2 $R1 $R2 $sampleID.R1.fq $sampleID.R1.unpaired.fq $sampleID.R2.fq $sampleID.R2.unpaired.fq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 CROP:120 SLIDINGWINDOW:40:15 MINLEN:36 

fastqc -t 2 -o fastqc_output $outfolder/$sampleID.R1.fq $outfolder/$sampleID.R2.fq
