# Reads selection and bowtie2 alignment
for sample in MC_4C_rep1 MC_4C_rep2 FL_4C_rep1 FL_4C_rep2
do
gunzip ${sample}.fastq.gz
grep -B 1 -A 2 -i 'TGGATGTTGCTGAAAATTGCAGTGTAATTTTCATG' ${sample}.fastq > ${sample}_bait.fastq
seqtk seq -A ${sample}_bait.fastq > ${sample}_bait.fasta

bowtie2 -p 12 -N 0 -5 35 \
        --un ${sample}_unaligned.sam \
        -x /data05/genomes/mm10/mm10 \
        -f -U ${sample}_bait.fasta \
        -S ${sample}_aligned.sam 2> ${sample}.bowtie2.log
samtools view -@ 8 -b ${sample}_aligned.sam > ${sample}_aligned.bam
samtools sort -@ 8 ${sample}_aligned.bam -o ${sample}_aligned.sorted.bam
samtools index ${sample}_aligned.sorted.bam
done
# Combine replicates for visualization
samtools merge MC_4C_combined_aligned.bam MC_4C_rep1_aligned.sorted.bam MC_4C_rep2_aligned.sorted.bam -@ 12
samtools merge FL_4C_combined_aligned.bam FL_4C_rep1_aligned.sorted.bam FL_4C_rep2_aligned.sorted.bam -@ 12
samtools sort MC_4C_combined_aligned.bam -o MC_4C_combined_aligned.sorted.bam -@ 12 
samtools sort FL_4C_combined_aligned.bam -o FL_4C_combined_aligned.sorted.bam -@ 12 
samtools index MC_4C_combined_aligned.sorted.bam
samtools index FL_4C_combined_aligned.sorted.bam
bamCoverage -b MC_4C_combined_aligned.sorted.bam -o MC_4C_combined_aligned.sorted.bw
bamCoverage -b FL_4C_combined_aligned.sorted.bam -o FL_4C_combined_aligned.sorted.bw
./bigWigToBedGraph ./MC_4C_combined_aligned.sorted.bw ./MC_4C_combined_aligned_chr11.sorted.bedGraph -chrom=chr11 
./bigWigToBedGraph ./FL_4C_combined_aligned.sorted.bw ./FL_4C_combined_aligned_chr11.sorted.bedGraph -chrom=chr11 
# Normalization 
## normalize_bedgraph.py is a custum script provided in https://github.com/porchard/normalize_bedgraph/blob/master/src/normalize_bedgraph.py
python normalize_bedgraph.py --to-number-reads 10000000 MC_4C_combined_aligned_chr11.sorted.bedGraph  > MC_4C_combined_aligned_chr11_normalized_TPM.bedGraph
python normalize_bedgraph.py --to-number-reads 10000000 FL_4C_combined_aligned_chr11.sorted.bedGraph  > FL_4C_combined_aligned_chr11_normalized_TPM.bedGraph
Rscript SmoothData.R  MC_4C_combined_aligned_chr11_normalized_TPM.bedGraph 5 MC_4C_combined_aligned_chr11_normalized_TPM_wd3  MC_4C_combined_aligned_chr11_normalized_TPM_wd5.bedGraph False
Rscript SmoothData.R  FL_4C_combined_aligned_chr11_normalized_TPM.bedGraph 5 FL_4C_combined_aligned_chr11_normalized_TPM_wd3  FL_4C_combined_aligned_chr11_normalized_TPM_wd5.bedGraph False