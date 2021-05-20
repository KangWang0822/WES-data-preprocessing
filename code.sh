######################################
#######germline variants calling######
######################################
module load bioinfo-tools
module load FastQC/0.11.8
module load bwa/0.7.17
module load samtools/1.10
module load GATK/4.1.4.1
##############SWE################
for sample in *_1.fq;
do
  r1=$sample;
  r2=${sample/_1.fq/}_2.fq;
  bwa mem -R "@RG\\tID:readgroupx\\tPU:lanex_flowcellx\\tSM:HG00097\\tLB:libraryx\\tPL:illumina" -t 1 human_g1k_v37_chr2.fasta $r1 $r2 | samtools sort > $(basename $sample _1.fq).bam
done
##########Ordering bam##########
for sample in *.bam;
do
samtools index $sample
done
###########bam view#############
samtools view -H HG00097.bam 
#########Mark Duplicates########
for sample in *.bam;
do
gatk --java-options -Xmx7g MarkDuplicates \
      -I $sample \
      -O $(basename $sample .bam)_duplicates.bam \
      -M $(basename $sample .bam)_duplicates.txt
done  
#######Recalibrate Base Quality Scores######
for sample in *_duplicates.bam;
do
gatk --java-options -Xmx7g BaseRecalibrator -R human_g1k_v37_chr2.fasta -I $sample --known-sites /sw/courses/ngsintro/reseq/data/ref/1000G_phase1.snps.high_confidence.b37.chr2.vcf -O $(basename $sample _duplicates.bam).recal.table
gatk --java-options -Xmx7g ApplyBQSR -R human_g1k_v37_chr2.fasta -I $sample --bqsr-recal-file $(basename $sample _duplicates.bam).recal.table -O $(basename $sample _duplicates.bam).recal.bam     
done    


cd /home/kangwang/WES/reference
ln -s /sw/GATK/hg38

cd /home/kangwang/WES/data
for i in /proj/sens2019581/data/WES_promix/*
do ln -s $i   #创建数据链接 ln
done

cat >BWA.sh

#!/bin/bash -l
#SBATCH -A sens2019581
#SBATCH -p node
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH -J bwa
bwa index -a bwtsw /home/kangwang/WES/reference/genome.fa
