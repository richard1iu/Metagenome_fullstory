
# Workflow
HISAT2→Stringtie→DESeq2 (Fast and simple)
1. Hisat2: genome alignment
2. Stringtie: transcripts alignment and qualification
3. DEseq2: differential gene expression

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

## Genome alignment
```bash
# index
hisat2-build -p 4 ./Genome/Gmax_275_v2.0.fa ./Genome/

# single end
hisat2 -p 4 -t -x ./Genome/ -U sample1.fastq -S sample1.sam --dta-cufflinks --no-unal --un-conc

# paired end
hisat2 -p 4 -t -x ./Genome/ -1 sample1.R1.fastq -2 sample1.R2.fastq -S sample1.sam

# sort
samtools view -Sb sample1.sam | samtools sort -@ 4 -o sample1.sorted.bam -
```

## Transcripts alignment and qualification
```bash
# transform gff to gtf if needed
gffread ./Genome/Gmax_275_Wm82.a2.v1.gene_exons.gff3 -T -o ./Genome/genome.gtf

# assemble and align
stringtie -p 4 -l sample -G ./Genome/genome.gtf -o sample1.gtf sample1.sorted.bam

# merge gtf of all samples
ls -1 ./*.gtf >> sample_gtf_list.txt
stringtie --merge -p 4 -o stringtie_merged.gtf -G Reference.gtf sample_gtf_list.txt

# compare sample.gtf with reference gtf
gffcompare -r ./Genome/genome.gtf -o gffcompare sample1.gtf

# quantification
stringtie -e -B -p 4 -G genome.gtf -o ballgown/sample.gtf sample.gtf
```
# convert ballgown to deseq2
Download the script of [prepDe.py](http://ccb.jhu.edu/software/stringtie/dl/prepDE.py)

```bash
# make a gtf_list with two column
ls *.gtf > path.txt
ls *.gtf | xargs -i basename {} .gtf > sampleid.txt
paste sampleid.txt path.txt > gtf_list.txt
#
conda activate py2
python ./prepDE.py -i gtf_list.txt
```

## Preprocess of gene id
If a certain gene correspond to multiple transcripts, the **labels (rownames)** of gene_id will be converted to the defualt label of stringtie **(STRG or MSTRG)**，and their counts are the sum of all transcripts.

We could use the **.annotated.gtf** file  or **.tracking** file  to reannotate  the labels of gene_count_matrix or trancript_count_matrix

```R
library(dplyr)
library(rtracklayer)

# input gene_count matrix
expr.gene <- read.csv("stringtie_gene_count_matrix.csv")

# input annotated.gtf
ensembl_anno <- rtracklayer::import('stringtie_merge.annotated.gtf') %>% as.data.frame()

# combine two file
anno_result <- dplyr::left_join(expr.gene,ensembl_anno[,(11:14)],by ="gene_id")

# remove duplicate gene id
anno_result <- anno_result[!duplicated(anno_result$gene_id),]
```

## Deseq2
```R
library(DESeq2)
database <- as.matrix(read.csv("transcript_count_matrix.csv", row.names="transcript_id"))
condition <- factor(c("control","control","KD","KD"))
coldata <- data.frame(row.names = colnames(database), condition)
countData <- countData[, rownames(colData)]
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resordered <- res[order(res$padj),]
summary(res)
write.csv(as.data.frame(resordered),file="results.csv")
```

