
## Reference materials
1. A survey of best practices for RNA-seq data analysis
2. CBC-UCONN

![RNA-seq Map](https://github.com/CBC-UCONN/RNA-Seq-Model-Organism-Arabidopsis-thaliana/blob/master/Hisat_stringtie_ballgown_WF.png)
>Figure source: CBC-UCONN

## 0.experimental design
## 1.quality control of raw reads
### 1.1 quality check
- sequence quality
- GC content
- the presence of adaptors
- overrepresented k-mers
- duplicated reads
#### software:
- FastQC on any platform
- NGSQC on Roche 454 and Illumina reads

### 1.2 quality control
- discard low-quality reads
- eliminate poor-quality bases
- trim adaptor and/or barcode sequences
#### software
- FASTX-Toolkit
- Trimmomatic
- sickle
- fastp

## 2. read alignment
### 2.1.1.genome mapping (gapped mapper)
- `HISAT2`: feature counting; lower memory (~4.3GB for human genome) than STAR (~27GB)
- `TopHat`: (1) **unspliced reads** are first mapped to locate exons, (2) **unmapped reads** are split and aligned independently to identify exon junctions. 
- `GSNAP`, `PALMapper`, `MapSplice`: optimized to identify SNPs or indels
- `STAR`, `MapSplice`: detect non-canonical splice junctions
- `GEM`: ultra-fast mapping

### 2.1.2.transcriptome mapping (ungapped mapper)
- `Salmon`: pseudo alignment 
- `Bowtie`
### 2.1.3.reference-free assembly (ungapped mapper)
1. `Trinity`: assembly
2. `Bowtie`: assembly into transcripts

The alignment step is followed by a quality evaluation and feature counting 
### 2.2.1. Quality evaluation in read alignment
1. Percentage of mapped reads,
It is the overall sequencing accuracy and the presence of contaminating DNA
2. Uniformity of read coverage on exons and the mapped strand. 
    - If reads primarily accumulate at the 3’ end of transcripts in poly(A)-selected samples, this might indicate low RNA quality in the starting material. 
    - The GC content of mapped reads may reveal PCR biases.
### software
- Picard 
- RSeQC
- Qualimap2
- rnaQUAST

### 2.2.2. quantification of expression levels of gene and transcript 

- `HTSeq`: estimate reference-genome expression levels 
- `TopHat`: Using an expectation-maximization approach that estimates transcript abundances.
- `Cufflinks`: PE reads, and may use GTF information to identify expressed transcripts, or can infer transcripts de novo from the mapping data alone.
- `RSEM`, `eXpress`, `Sailfish`, `kallisto`: Quantify expression from transcriptome mapping. Allocate multi-mapping reads among transcript and output within-sample normalized values corrected for sequencing biases.
- `NURD`: Provides an efficient way of estimating transcript expression from SE reads with a low memory and computing cost.
- EDGE-pro: prokaryotic genome 
#### quality control of count data
evaluate GC content and gene length biases
- R package: {NOISeq} 
- R package: {EDASeq}
- featureCounts 
- htseq-count

## 3.differential gene expression analysis
1. `limma` perform well and run the fastest in most circumstances. 
2. `DESeq and edgeR` perform similarly in ranking genes but are often relatively conservative or too liberal, respectively, in controlling FDR. 
3. `SAMseq` performs well in terms of FDR but presents an acceptable sensitivity when the number of replicates is relatively high, at least 10. 
4. `NOISeq` and `NOISeqBIO` (the adaptation of NOISeq for biological replication) are more efficient in avoiding false positive calls at the cost of some sensitivity but perform well with different numbers of replicates.
5. `Cuffdiff` and `Cuffdiff2` performe surprisingly poorly in the comparisons. 
6. `edgeR`, `limma-voom`, `DESeq`, `DESeq2`, and `maSigPro` can perform multiple group comparisons, include different covariates or analyze time-series data.

## 4.alternative splicing
Alternative splicing (AS) is a post-transcriptional process which **generates different transcripts from the same gene** and is vital in response to environmental stimuli by producing diverse protein products. 
- `TopHat` and its downstream tool, `FineSplice`: are the fastest tools
- `Alt Event Finder`: can detect the highest number of junctions
- `RSR`: detects the lowest number of junctions. 
- `rMATS` is faster than `rSeqDiff` but detects less differentially spliced isoforms than rSeqDiff.

## 5.Functional analysis
Characterize molecular functions or pathways in which differentially expressed genes (DEGs) are involved
### Traditional methods for microarray technology
1. Comparing a list of DEGs against the rest of the genome for overrepresented functions
2. Gene set enrichment analysis (GSEA)
RNA-seq biases such as gene length complicate the direct applications of these methods for **count data**
### RNA-seq-specific tools
1. GOseq estimates a bias effect (such as gene length) on differential expression results and adapts the traditional hypergeometric statistic used in the functional enrichment test to account for this bias. 
2. Gene Set Variation Analysis (GSVA)  or SeqGSEA packages also combine splicing and implement enrichment analyses similar to GSEA.

### Functional annotation for model species
- Gene Ontology
- Bioconductor
- DAVID 
- Babelomics 
### Functional profiling of de novo
1. annotate Protein-coding transcripts using orthology by searching for similar sequences in protein databases
    - SwissProt
    - Pfam and InterPro: contain **conserved protein domains**
2. The use of standard vocabularies such as the Gene Ontology (GO) allows for some exchangeability of functional information across orthologs.
    - Blast2GO  allow massive annotation of complete transcriptome datasets against a variety of databases and controlled vocabularies

### functional annotation non-coding RNAs functional annotation
The functional annotation of these long non-coding RNAs is more challenging as **their conservation is often less pronounced** than that of protein-coding genes.
- `Rfam` database contains most well-characterized RNA families, such as ribosomal or transfer RNAs
- `mirBase` or `Miranda` are specialized in miRNAs. 
These resources can be used for **similarity-based annotation of short non-coding RNAs**, but no standard functional annotation procedures are available yet for other RNA types such as the **long non-coding RNAs**.

## 6.gene fusion detection and eQTL mapping
### 6.1.integration with other technologies
#### 6.1.1 small and other non-coding RNAs
#### 6.1.2 gene fusion discovery
#### 6.1.3 long-read
#### 6.1.4 single-cell analysis

### 6.2.other RNA-seq
#### 6.2.1 eQTL/sQTL
#### 6.2.2 Chromatin
#### TF binding
#### Proteomics/metabolomics

# Outlook
Two major highlights in the current application of RNA-seq
1. the construction of transcriptomes from `small amounts of starting materials`
2. better transcript identification from `longer reads`.

## Single-cell RNA-seq

## Long-read sequencing