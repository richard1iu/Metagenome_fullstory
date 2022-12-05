
# Workflow
subread→featureCounts

## Install
```bash
# in shell
conda install -c conda-forge hisat2 stringtie
```

```R
# in R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
```

## Refence Genome files

1. Reference genome：.fasta, .fa
2. Annotation file: .gff, .gtf, .gff3

> NOTE: .fa.gz or .gff.gz is the compressed file of .fa or .gff, you can use it directly in most software or you can use `gunzip` to uncompressed it if some software  does not work.

## Quality control

```
 fastqc -t 2 sample1.fastq sample2.fastq 
```


## 
```bash
# build index
subread-buildindex -o subread_index  reference_genome.fa

# map fastq to reference genome
for R1 in $(ls *_R1.fq);
do
	sample=$(basename ${R1} _R1.fq)
	echo "Processin sample ${sampe}" 
	subjunc -i $index \
		-r $R1 \
		-R ${sample}_R2.fq \
		-T 4 -o align/${sample}_subread.bam
done
```

## quantification
```bash
#
featureCounts   -T 4 \
                -p \
                -t exon \ # only quantify exon gene
                -g gene_name \ # rowname re gene_name
                -a reference_genome.gtf \
                -o  ./exon_counts.txt  \
                 *.bam

featureCounts -T 4 -p -t exon -g gene_id -a reference_genome.gtf -o  ./gene_counts.txt   *.bam
```