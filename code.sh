######################################
#######germline variants calling######
######################################
cd /proj/g2021013/nobackup/kangwang/ngsworkflow
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
#################generate VCF###############
for sample in *.recal.bam;
do
gatk --java-options -Xmx7g HaplotypeCaller -R human_g1k_v37_chr2.fasta -ERC GVCF -I $sample -O $(basename $sample .recal.bam).g.vcf 
done 
##################combine VCF###############
gatk --java-options -Xmx7g CombineGVCFs -R human_g1k_v37_chr2.fasta -V HG00097.vcf -V HG00100.vcf -V HG00101.vcf -O cohort.g.vcf
gatk --java-options -Xmx7g GenotypeGVCFs -R human_g1k_v37_chr2.fasta -V cohort.g.vcf -O cohort.vcf
#Variant filtration SNPs:
gatk --java-options -Xmx7g SelectVariants -R human_g1k_v37_chr2.fasta -V cohort.vcf --select-type-to-include SNP -O cohort.snvs.vcf
gatk --java-options -Xmx7g VariantFiltration -R human_g1k_v37_chr2.fasta -V cohort.snvs.vcf -O cohort.snvs.filtered.vcf --filter-name QDfilter --filter-expression "QD < 2.0" --filter-name MQfilter --filter-expression "MQ < 40.0"  --filter-name FSfilter --filter-expression "FS > 60.0"
#Variant filtration indels
gatk --java-options -Xmx7g SelectVariants -R human_g1k_v37_chr2.fasta -V cohort.vcf --select-type-to-include INDEL -O cohort.indels.vcf
gatk --java-options -Xmx7g VariantFiltration -R human_g1k_v37_chr2.fasta -V cohort.indels.vcf -O cohort.indels.filtered.vcf --filter-name QDfilter --filter-expression "QD < 2.0" --filter-name FSfilter --filter-expression "FS > 200.0"
#Merge filtered SNPs and indels
gatk --java-options -Xmx7g MergeVcfs -I cohort.snvs.filtered.vcf -I cohort.indels.filtered.vcf -O cohort.filtered.vcf


#########################
#SBATCH -A g2021013
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J BestPractise
## load modules
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.10
module load GATK/4.1.4.1
#Define path to reference genome
ref="/sw/courses/ngsintro/reseq/data/ref"
#Make symbolic links
ln -s /sw/courses/ngsintro/reseq/data/ref/human_g1k_v37_chr2.fasta
ln -s /sw/courses/ngsintro/reseq/data/fastq/HG00097_1.fq
ln -s /sw/courses/ngsintro/reseq/data/fastq/HG00097_2.fq
ln -s /sw/courses/ngsintro/reseq/data/fastq/HG00100_1.fq
ln -s /sw/courses/ngsintro/reseq/data/fastq/HG00100_2.fq
ln -s /sw/courses/ngsintro/reseq/data/fastq/HG00101_1.fq
ln -s /sw/courses/ngsintro/reseq/data/fastq/HG00101_2.fq
#Index reference genome
bwa index -a bwtsw human_g1k_v37_chr2.fasta
samtools faidx human_g1k_v37_chr2.fasta
gatk --java-options -Xmx7g CreateSequenceDictionary -R human_g1k_v37_chr2.fasta -O human_g1k_v37_chr2.dict
## loop through the samples:
for sample in HG00097 HG00100 HG00101;
do
  echo "Now analyzing: "$sample
  #Map the reads:
  bwa mem -R "@RG\\tID:readgroupx\\tPU:flowcellx_lanex\\tSM:"$sample"\\tLB:libraryx\\tPL:illumina" -t 1 human_g1k_v37_chr2.fasta $sample"_1.fq" $sample"_2.fq" | samtools sort > $sample".bam"
  samtools index $sample".bam"
  #Mark duplicates:
  gatk --java-options -Xmx7g MarkDuplicates -I $sample".bam" -O $sample".md.bam" -M $sample"_mdmetrics.txt"
  #Base quality score recalibration:
  gatk --java-options -Xmx7g BaseRecalibrator -R human_g1k_v37_chr2.fasta -I $sample".md.bam" --known-sites $ref"/1000G_phase1.snps.high_confidence.b37.chr2.vcf" -O $sample".recal.table"
  gatk --java-options -Xmx7g ApplyBQSR -R human_g1k_v37_chr2.fasta -I $sample".md.bam" --bqsr-recal-file $sample".recal.table" -O $sample".recal.bam"
  #HaplotypeCaller in -ERC mode:
  gatk --java-options -Xmx7g HaplotypeCaller -R human_g1k_v37_chr2.fasta -ERC GVCF -I $sample".bam" -O $sample".g.vcf" 
done
#Joint genotyping:
gatk --java-options -Xmx7g CombineGVCFs -R human_g1k_v37_chr2.fasta -V HG00097.g.vcf -V HG00100.g.vcf -V HG00101.g.vcf -O cohort.g.vcf
gatk --java-options -Xmx7g GenotypeGVCFs -R human_g1k_v37_chr2.fasta -V cohort.g.vcf -O cohort.vcf
#Variant filtration SNPs:
gatk --java-options -Xmx7g SelectVariants -R human_g1k_v37_chr2.fasta -V cohort.vcf --select-type-to-include SNP -O cohort.snvs.vcf
gatk --java-options -Xmx7g VariantFiltration -R human_g1k_v37_chr2.fasta -V cohort.snvs.vcf -O cohort.snvs.filtered.vcf --filter-name QDfilter --filter-expression "QD < 2.0" --filter-name MQfilter --filter-expression "MQ < 40.0"  --filter-name FSfilter --filter-expression "FS > 60.0"
#Variant filtration indels
gatk --java-options -Xmx7g SelectVariants -R human_g1k_v37_chr2.fasta -V cohort.vcf --select-type-to-include INDEL -O cohort.indels.vcf
gatk --java-options -Xmx7g VariantFiltration -R human_g1k_v37_chr2.fasta -V cohort.indels.vcf -O cohort.indels.filtered.vcf --filter-name QDfilter --filter-expression "QD < 2.0" --filter-name FSfilter --filter-expression "FS > 200.0"
#Merge filtered SNPs and indels
gatk --java-options -Xmx7g MergeVcfs -I cohort.snvs.filtered.vcf -I cohort.indels.filtered.vcf -O cohort.filtered.vcf
