# Step.1 Trimming
---------------------
## Scythe: (very slow)

uses a **Naive Bayesian** approach to classify contaminant substrings

```bash
scythe -a adapters.fasta -o SRR957824_adapt_R1.fastq SRR957824_500K_R1.fastq.gz
scythe -a adapters.fasta -o SRR957824_adapt_R2.fastq SRR957824_500K_R2.fastq.gz
```

## Sickle
uses **sliding windows** along with **quality and length thresholds** to determine when quality is sufficiently low to trim the 3'-end of reads and also determines when the quality is sufficiently high enough to trim the 5'-end of reads.

```bash
sickle pe -f SRR957824_adapt_R1.fastq -r SRR957824_adapt_R2.fastq \
    -t sanger -o SRR957824_trimmed_R1.fastq -p SRR957824_trimmed_R2.fastq \
    -s /dev/null -q 25
```

## Trimmomatic
remova adapter and low-quality reads

```bash
for gz in *_R1.fastq.gz
do

# 提取双端公共文件名，并输出检验
base=$(basename $gz _R1.fastq.gz) # base=${file/_R1.fastq.gz/}
echo $base

# 运行去接头程序
TrimmomaticPE -threads 9 \
     ${base}_R1.fastq.gz \
     ${base}_R2.fastq.gz \
     ${base}_R1.qc.fq.gz ${base}_s1_se \
     ${base}_R2.qc.fq.gz ${base}_s2_se \
     ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 \
     LEADING:2  \
     TRAILING:2 \
     SLIDINGWINDOW:4:2 \
     MINLEN:25 
done
```
## khmer

An alternative to **trimming based on the quality scores** is to **trim based on k-mer abundance** - this is known as k-mer spectral error trimming. K-mer spectral error trimming always beats quality score trimming in terms of eliminating errors

if there are low abundance k-mers in a high coverage dataset, those k-mers are almost certainly errors result.
but sometimes they may represent truly low abundance (but useful) data.

```bash
# set as env variables
export PATH="/usr/lib/khmer/bin:$PATH"

# estimate the unique (abundance=1) k-mers abundance before and after trim
abundance-dist-single.py -M 1e9 -k 21 SRR1976948_1.fastq.gz SRR1976948_1.fastq.gz.dist
abundance-dist-single.py -M 1e9 -k 21 SRR1976948_1.trimmed.fq.gz SRR1976948_1.qc.fq.gz.dist

# trim off low abundance k-mers in high coverage dataset, but keep the k-mers in low coverage dataset
interleave-reads.py SRR1976948_1.qc.fq.gz SRR1976948_2.qc.fq.gz | trim-low-abund.py -V -M 8e9 -C 3 -Z 10 - -o SRR1976948.trim.fq

# see how many k-mers are removed
unique-kmers.py SRR1976948_1.qc.fq.gz SRR1976948_2.qc.fq.gz
unique-kmers.py SRR1976948.trim.fq
```

# Step.2 Mapping/alignment
---------------------
## Bowtie2
***algorithm:** Burrows–Wheeler transform*

Before aligning the reads against a reference, it is necessary to build an index of that reference, which makes searching for patterns **much much faster**.
reference genome can be the results of assembly (e.g. the output fa of megahit)

```bash
# 1.indexing the reference genome 
bowtie2-build p0157_Sakai.fasta.gz p0157_Sakai

# 2.Aligning reads
bowtie2 -x p0157_Sakai -1 SRR957824_trimmed_R1.fastq.gz \
-2 SRR957824_trimmed_R2.fastq.gz -S SRR957824.sam
```

## bwa
***algorithm:** Burrows–Wheeler transform*

```bash
bwa index p0157_Sakai.fasta.gz -p p0157_Sakai
bwa mem -t 4 SRR957824_trimmed_R1.fastq.gz SRR957824_trimmed_R2.fastq.gz > fastq.mem.sam
```

# Option.1 Visualising the sam
---------------------
## samtools

```bash
# 1.convert SAM file into BAM, a compressed version of sam, for sort and index
samtools view -hSb -o SRR957824.bam SRR957824.sam

# 2.Sort the bam file per position in the genome and index it
samtools sort SRR957824.bam SRR2584857.sorted.bam
samtools index SRR2584857.sorted.bam

# 3.visualise
samtools tview SRR2584857.sorted.bam pO157_Sakai.fasta.gz
```

# Option.2 Variant Calling
---------------------
## samtools
finding positions where the reads are systematically different from the reference genome. 
Single nucleotide polymorphism (SNP)-based typing is particularly popular and used for a broad range of applications. 

```bash
samtools mpileup -uD -f pO157_Sakai.fasta.gz SRR2584857.sorted.bam | bcftools view - > variants.vcf
```

# Step.3 Genome assembly
---------------------
## MEGAHIT（De-novo）
if run error, try transfer the format of input file from gz to fastq and then run it again.

```bash
# megahit is more rapid and resource saving
megahit -1 ERR486840_1.fastq.gz -2 ERR486840_2.fastq.gz -o m_genitalium
```
## SPAdes
```
# SPAdes/metaSPAdes for fastq <100Mb, it is most accurate but time-cost
spades.py --pe1-1 ERR486840_1.fastq.gz --pe1-2 ERR486840_2.fastq.gz -o assembly/spades/
```

# Step.4 Quality evaluation of the Assembly
---------------------
## QUAST
Evaluate quality of assembly with quality metrics, it can be used in the case with or without reference genome
```bash
# 1.without reference
quast.py m_genitalium.fasta -o m_genitalium_report

# 2. with reference

```

## BUSCO (Find marker genes)

Marker genes are conserved across a range of species and finding intact conserved genes in our assembly would be a good indication of its quality.

1.	Firstly, constructed the minimum gene set of different species
2.	and then use HMMER, BLAST, Augustus and other tools to analyze the homologous genes in the assembly results to quantitatively evaluate whether the assembly is complete or not.

```bash
BUSCO.py -i m_genitalium.fasta -l bacteria_odb9 -o busco_genitalium -m genome
BUSCO.py -i m_genitalium.fasta -l bacteria_odb9 -o busco_genitalium -m genome -f
 ```

**Results:**
**C:** the integrity compared with the BUSCO set
**D:** the number of duplication
**M:** the number of genes that may be missing

# Option.1 Fixing misassemblies
---------------------
## Pilon

Pilon attempts to improve the input genome, including:
- Single base differences
- Small Indels
- Larger Indels or block substitution events
- Gap filling
- Identification of local misassemblies, including optional opening of new gaps

Pilon then outputs a FASTA file containing: 
- an improved representation of the genome from the read data 
- and an optional VCF file detailing variation seen between the read data and the input genome.

```bash
# 1.map our reads against the assembly
bowtie2-build assembled_contigs.fa final_contigs
bowtie2 -x final_contigs -1 ERR486840_1.fastq.gz -2 ERR486840_2.fastq.gz | samtools view -bS -o assembled_contigs.bam
samtools sort assembled_contigs.bam -o assembled_contigs.sorted.bam
samtools index assembled_contigs.sorted.bam

# 2.fix
pilon --genome assembled_contigs.fa --frags assembled_contigs.sorted.bam --output assembled_contigs_improved
```
# Step.5 Binning
----------------------
Firstly,  map the reads back against the assembly to get coverage information

```bash
# 1.build an index of that reference, the input file is the assemblyed contigs fa
bowtie2-build assembled.contigs.fa final.contigs
# 2.mapping
bowtie2 -x final.contigs -1 tara_reads_R1.fastq.gz -2 tara_reads_R2.fastq.gz | samtools view -bS -o tara_unsorted.bam
# 3*.sort the bam (NECESSARY!!!)
samtools sort tara_unsorted.bam -o tara_sorted.bam
# 4.index the bam
samtools index tara_sorted.bam
```

## concoct 1.1.0 
Unsupervised binning of metagenomic contigs by using **nucleotide composition**, **kmer frequencies**,**coverage data** in multiple samples and **linkage data** from paired end reads

```bash
# set the path of Database (NECESSARY!!!)
export CHECKM_DATA_PATH=/mnt/f/Database/checkm_data_2015_01_16

# Cut contigs into smaller parts
cut_up_fasta.py assembled_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa

# coverage depth information (input files are sorted and indexed bam)
concoct_coverage_table.py contigs_10K.bed mapping/Sample*.sorted.bam > coverage_table.tsv

# binning
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/

# Merge subcontig clustering into original contig clustering
merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv

# Extract bins as individual FASTA
mkdir concoct_output/fasta_bins
extract_fasta_bins.py assembled_contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins

```

# Step.6 Checking the quality of the bins

## checkm

If checkm **fails at the phylogeny step**, it is likely that your pc **doesn't have enough RAM**.  
pplacer requires about **35G of RAM** to place the bins in the tree of life.

**workflow:**
1. (M) > checkm tree <bin folder> <output folder>
2. (M) > checkm lineage_set <output folder> <marker file>
3. (M) > checkm analyze <marker file> <bin folder> <output folder>
4. (M) > checkm qa <marker file> <output folder>

For convenience, the 4 mandatory steps can be executed using one-step command:
> checkm lineage_wf <bin folder> <output folder>

```bash
# check the quality of bins, input file are the dir contain the fa file from binning
checkm lineage_wf -t 8 -x fa --tab_table -f bins_qa.txt concoct_output/fasta_bins/ ../checkm
```


# Step.7 Genome Annotation
---------------------
## Prokka
Prokka is a "wrapper" collecting several pieces of software from various authors

Prokka finds and annotates features (both protein coding regions and RNA genes, i.e. tRNA, rRNA) present on on a sequence.

Prokka can be used to annotate bacterial, archaeal and viral genomes quickly, generating standard output files in GenBank, EMBL and gff formats.

```bash
# setup database
prokka --setupdb

# check if there are kingdom database. if not, download the db file from github and substitute /var/lib/prokka/db
prokka --listdb

# with both .faa and assembled contigs .fa as input file 
pokkar --outdir prokka_annotation --kingdom Bacteria --proteins uniprot_mycoplasma_reviewed.faa assembled_contigs.fa

# just the assembled contigs .fa as input file
prokka --outdir prokka_annotation --kingdom Bacteria --prefix metag --metagenome  assembled_contigs.fa
```

**output:**
- The GFF and GBK files contain all the information about the features annotated (in different formats.)
- The .txt file contains a summary of the number of features annotated.
- The .faa file contains the protein sequences of the genes annotated.
- The .ffn file contains the nucleotide sequences of the genes annotated.

## prodigal (gene calling)

Protein-coding gene prediction for **prokaryotic** genomes, it is **unsupervised and rapid**

```bash
#
prodigal -i my.genome.fna -o my.genes -a my.proteins.faa
#
prodigal -i my.metagenome.fna -o my.genes -a my.proteins.faa -p meta
```

# Gene Abundance Estimation
## salmon
Salmon is a very fast RNAseq counting packages, which counts fragments without doing up-front read mapping. 
Salmon can be used with edgeR and others to do differential expression analysis (if you are quantifying RNAseq data).

```bash


# if the fastq are paried-merged, first we need split it
for file in *_merged.fq.gz
do   
tail=.fq.gz
base=${file/$tail/}
split-paired-reads.py $file -1 ${base}.1.fq -2 ${base}.2.fq;
 done  

# create index with .ffn (the nucleotide predicted protein regions) from Prokka
salmon index -t metagG.ffn -i transcript_index --type quasi -k 31

# quantify our reads against this reference
for file in *.1.fq
do
tail1=.1.fq
tail2=.2.fq
base=${file/$tail1/}
salmon quant -i transcript_index --libType IU \ # --libType must come before the read files!
      -1 ${base}_1.fq.gz -2 $base$tail2 -o $base.quant;
 done
```

