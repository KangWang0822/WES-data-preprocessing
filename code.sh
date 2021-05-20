module load bioinfo-tools
module load FastQC/0.11.8
module load bwa/0.7.17
module load samtools/1.10
module load GATK/4.1.4.1

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
