
## Reference materials
1. A survey of best practices for RNA-seq data analysis
2. CBC-UCONN
3. Systematic comparison and assessment of RNA-seq procedures for gene expression quantitative analysis

![RNA-seq Map](https://github.com/CBC-UCONN/RNA-Seq-Model-Organism-Arabidopsis-thaliana/blob/master/Hisat_stringtie_ballgown_WF.png)
>Figure source: CBC-UCONN

## 0.experimental design
## 1.Quality control of raw reads
### 1.1 Quality check
1. Sequence quality
2. GC content
3. The presence of adaptors
4. Overrepresented k-mers
5. Duplicated reads
#### Software:
1. FastQC & multiqc on any platform
2. NGSQC on Roche 454 and Illumina reads

### 1.2 quality filter
- Discard low-quality reads
- Eliminate poor-quality bases
- Trim adaptor and/or barcode sequences
#### Software
1. FASTX-Toolkit
2. Trimmomatic
3. Sickle
4. Fastp

## 2. Read Alignment
#### 2.1 Alignment
#### 2.1.1 Genome mapping (gapped mapper)
1. `HISAT2`: feature counting; lower memory (~4.3GB for human genome) than STAR (~27GB)
2. `TopHat`: (1) **unspliced reads** are first mapped to locate exons, (2) **unmapped reads** are split and aligned independently to identify exon junctions. 
3. `GSNAP`, `PALMapper`, `MapSplice`: optimized to identify SNPs or indels
4. `STAR`, `MapSplice`: detect non-canonical splice junctions
5. `GEM`: ultra-fast mapping

#### 2.1.2 Transcriptome mapping (ungapped mapper)
1. `Salmon`: pseudo alignment 
2. `Bowtie2`:
3. `RSEM`:
#### 2.1.3 Reference-free assembly (ungapped mapper)
1. `Trinity`: assembly
2. `Bowtie2`: assembly into transcripts

The alignment step is followed by a quality evaluation and gene counting 
### 2.2 Quality evaluation in read alignment
1. Percentage of mapped reads,
It is the overall sequencing accuracy and the presence of contaminating DNA
2. Uniformity of read coverage on exons and the mapped strand. 
    - If reads primarily accumulate at the 3’ end of transcripts in poly(A)-selected samples, this might indicate low RNA quality in the starting material. 
    - The GC content of mapped reads may reveal PCR biases.
#### Software
- Picard 
- RSeQC
- Qualimap2
- rnaQUAST

## 3. Quantification of expression levels of gene and transcript 

### 3.1 quantify genome
Input files: .sam/.bam and .gtf
1. `htseq-count` 
2. `featureCounts`: faster than htseq-count 
3. `EDGE-pro`: prokaryotic genome
5. `HTSeq`: estimate reference-genome expression levels 
5. `TopHat`: Using an expectation-maximization approach that estimates transcript abundances.
6. `Cufflinks`: PE reads, and may use GTF information to identify expressed transcripts, or can infer transcripts de novo from the mapping data alone.

### 3.2 quantify transcriptome
1. `RSEM`, `eXpress`, `Sailfish`, `kallisto`: Quantify expression from transcriptome mapping. Allocate multi-mapping reads among transcript and output within-sample normalized values corrected for sequencing biases.
2. `NURD`: Provides an efficient way of estimating transcript expression from SE reads with a low memory and computing cost.

### 3.3 quantify genome and/or transcriptome
- `stringtie`:

### 3.4 Normalization
1. FPKM (Fragments per Kilobase of Mapped reads)
2. RPKM (Reads per Kilobase of Mapped reads)
3. TMM (Trimmed Mean of M values )
4. TPM (Transcripts per Million )
5. RLE (Relative Log Expression)
6. UQ (Upper quartile )
7. Cov (Coverage)
8. Est_Counts (Estimated counts)
8. Eff_Counts (Effective counts)

## 4.Differential gene expression analysis

1. `limma` perform well and run the fastest in most circumstances. 
2. `DESeq and edgeR` perform similarly in ranking genes but are often relatively conservative or too liberal, respectively, in controlling FDR. 
3. `SAMseq` performs well in terms of FDR but presents an acceptable sensitivity when the number of replicates is relatively high, at least 10. 
4. `NOISeq` and `NOISeqBIO` (the adaptation of NOISeq for biological replication) are more efficient in avoiding false positive calls at the cost of some sensitivity but perform well with different numbers of replicates.
5. `Cuffdiff` and `Cuffdiff2` performe surprisingly poorly in the comparisons. 
6. `edgeR`, `limma-voom`, `DESeq`, `DESeq2`, and `maSigPro` can perform multiple group comparisons, include different covariates or analyze time-series data.

![Comparison of DGE](https://www.nature.com/articles/s41598-020-76881-x/figures/7)
>Source: Systematic comparison and assessment of RNA-seq procedures for gene expression quantitative analysis

## 5.Alternative splicing
Alternative splicing (AS) is a post-transcriptional process which **generates different transcripts from the same gene** and is vital in response to environmental stimuli by producing diverse protein products. 
- `TopHat` and its downstream tool, `FineSplice`: are the fastest tools
- `Alt Event Finder`: can detect the highest number of junctions
- `RSR`: detects the lowest number of junctions. 
- `rMATS` is faster than `rSeqDiff` but detects less differentially spliced isoforms than rSeqDiff.

## 6.Functional analysis
Characterize molecular functions or pathways in which differentially expressed genes (DEGs) are involved
### 6.1 Traditional methods for microarray technology
1. Comparing a list of DEGs against the rest of the genome for overrepresented functions
2. Gene set enrichment analysis (GSEA)
RNA-seq biases such as gene length complicate the direct applications of these methods for **count data**
### 6.2 RNA-seq-specific tools
1. `GOseq` estimates a bias effect (such as gene length) on differential expression results and adapts the traditional hypergeometric statistic used in the functional enrichment test to account for this bias. 
2. `GSVA` (Gene Set Variation Analysis)  or `SeqGSEA` packages also combine splicing and implement enrichment analyses similar to GSEA.

### 6.3 Functional annotation for model species
- Gene Ontology
- Bioconductor
- DAVID 
- Babelomics 
### 6.4 Functional profiling of de novo
1. annotate Protein-coding transcripts using orthology by searching for similar sequences in protein databases
    - SwissProt
    - Pfam and InterPro: contain **conserved protein domains**
2. The use of standard vocabularies such as the Gene Ontology (GO) allows for some exchangeability of functional information across orthologs.
    - Blast2GO: allow massive annotation of complete transcriptome datasets against a variety of databases and controlled vocabularies

### 6.5 Non-coding RNAs functional annotation
The functional annotation of these long non-coding RNAs is more challenging as **their conservation is often less pronounced** than that of protein-coding genes.
- `Rfam` database contains most well-characterized RNA families, such as ribosomal or transfer RNAs
- `mirBase` or `Miranda` are specialized in miRNAs. 
These resources can be used for **similarity-based annotation of short non-coding RNAs**, but no standard functional annotation procedures are available yet for other RNA types such as the **long non-coding RNAs**.

## 7. Gene fusion detection and eQTL mapping
### 7.1.integration with other technologies
#### 7.1.1 small and other non-coding RNAs
#### 7.1.2 gene fusion discovery
#### 7.1.3 long-read
#### 7.1.4 single-cell analysis

### 7.2.other RNA-seq
#### 7.2.1 eQTL/sQTL
#### 7.2.2 Chromatin
#### 7.2.3 TF binding
#### 7.2.4 Proteomics/metabolomics

## 8. Outlook
Two major highlights in the current application of RNA-seq
1. the construction of transcriptomes from `small amounts of starting materials`
2. better transcript identification from `longer reads`.

### 8.1 Single-cell RNA-seq

### 8.2 Long-read sequencing