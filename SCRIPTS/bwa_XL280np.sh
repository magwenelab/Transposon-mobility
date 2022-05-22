bwa mem -a -M ../REFERENCE/XL280_NP_NUCLEAR_FINAL_July2020.fasta ../FASTQ/JW-S1_S1_L002_R1_001.fastq.gz ../FASTQ/JW-S1_S1_L002_R2_001.fastq.gz | samtools view -b | samtools sort -o ../BAMS/JW-S1_S1_L002_.bam
samtools index ../BAMS/JW-S1_S1_L002_.bam

bwa mem -a -M ../REFERENCE/XL280_NP_NUCLEAR_FINAL_July2020.fasta ../FASTQ/JW-S1_S1_L003_R1_001.fastq.gz ../FASTQ/JW-S1_S1_L003_R2_001.fastq.gz | samtools view -b | samtools sort -o ../BAMS/JW-S1_S1_L003_.bam
samtools index ../BAMS/JW-S1_S1_L003_.bam

bwa mem -a -M ../REFERENCE/XL280_NP_NUCLEAR_FINAL_July2020.fasta ../FASTQ/JW-S1_S1_L004_R1_001.fastq.gz ../FASTQ/JW-S1_S1_L004_R2_001.fastq.gz | samtools view -b | samtools sort -o ../BAMS/JW-S1_S1_L004_.bam
samtools index ../BAMS/JW-S1_S1_L004_.bam

