
The fullstory of analyzing RNA sequencing data when a reference genome is available

## Reference material
[CBC-UCONN](https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/)

## Guides of all steps
1. quality control of the reads
2. alignment of the reads to the reference genome
3. conversion of the files to raw counts
4. analysis of the counts with DeSeq2
5. annotation of the reads using Biomart

## All software
- Sickle
- STAR
- HTSeq
- DESeq2
- biomaRt
- fastqc

## Script
### 1.Download SRA (optional)
1. using SRA number dirctly
    - prefetch
    - fasterq-dump
2. using website
    - aria2c
    - wget
    - other tools in windows

```bash
# 1.download using sratools
prefetch SRR391536
fasterq-dump SRR391536

# 2.download via website
srapath SRR391536 > web.txt
wget -c -t -i web.txt
aria2c web.txt
```
Create sample list for batch analysis
```bash
# create sample list
ls *.fastq | while read id;do basename $id .fastq >> seq_list.txt; done
```

### 2.SRA to fastq
```bash
# single-end using fasterq-dump
fasterq-dump  SRR391536.man (download using windows)
fasterq-dump  SRR391536.SRA (download using linux)

# paired-end using fasterq-dump
fasterq-dump --split-3 SRR391536.man (download using windows)

# fastq-dump
alias fd='fastq-dump --split-3 --defline-qual '+' --defline-seq '@\\\$ac-\\\$si/\\\$ri' '
fasterq-dump SRR391536.man
```

### 3.Quality check and filter
#### 3.1Quality check

```bash
# fastqc
cat ./seq_list.txt | parallel -j 4 \
    fastqc --outdir fastqc_result raw_seq/{}_{1..2}.fastq

# multiqc
multiqc -f -o fastqc_multiqc fastqc_result
```
#### 3.2 Trim reads
1.Sickle
1. se: unpaired reads
2. -f: input paired-end forward fastq (or single)
3. -r: input paired-end reverse fastq
4. -t: quality type: solexa, illumina, sanger
5. -o: output trimmed forward fastq (or single)
6. -p: output trimmed reverse fastq
7. -q: minimum quality score
8. -l：minimum read length
9. -s: output trimmed singles fastq file
```bash
# single file
## single end
sickle se -f SRR391535.fastq -t sanger -o trimmed_SRR391535.fastq -q 35 -l 45

## paired end
sickle se -t sanger -f SRR391535_R1.fastq -r SRR391535_R2.fastq  -o trimmed_SRR391535_R1.fastq -p trimmed_SRR391535_R2.fastq -s singles_SRR391535_R1.fastq -q 35 -l 45

# multiple file

```

2.Fastp
```
cat raw_seq/seq_list.txt | parallel -j 4 \
fastp \
	--in1 raw_seq/{}_1.fastq.gz \
	--in2 raw_seq/{}_2.fastq.gz \
	--out1 trim_seq/{}_trim_1.fastq.gz \
	--out2 trim_seq/{}_trim_2.fastq.gz \
	--json fastp_result/{}_fastp.json \
	--html fastp_result/{}_fastp.html
```

3.Trimmomatic
```bash
# paired
trimmomatic PE sub1_R1.fq.gz sub1_R2.fq.gz \
               sub1_R1.trim.fq sub1_R1.unpaired.trim.fq \
               sub1_R2.trim.fq sub1_R2.unpaired.trim.fq \
               SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
```
### 4.Alignment
#### 4.1 Map to reference genome
1.STAR
>Fast but need very large memory
```bash
STAR --runMode genomeGenerate \
    --genomeDir ./genome \ # index file directory
    --genomeFastaFiles ./Gmax_275_v2 \  # input .fasta file 
    --sjdbGTFfile ./Gmax_275_Wm82.a2.v1.gene_exons \  # annotation .gtf file
    --sjdbGTFtagExonParentTranscript Parent 
    --sjdbOverhang 100 \ # defualt 100, or read length minus 1
    --runThreadN 8 \
    --readFilesIn ./sample1.fasta \ # input sample
    --outFileNamePrefix sample1 \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts # quantify gene count and output files for transripts quantification using RSEM
```

2.Hisat2
> Fast and low memory
```bash
# 
GENOME=./genome/Gmax
INDEX=./hisat2_index/indexed
INDIR=./trim_seq
OUTDIR=./alignment
SEQLIST=./raw_seq/seq_list.txt
SAMPLE=$(sed -n ${NUM}p $SEQLIST)

# 1.index
hisat2-build -p 16 $GENOME $INDEX

# 2.align
## single
hisat2 -p 8 \ # processor number
  --dta \     # report alignments results
  -x $INDEX \ # the path of index dir
  -1 ./trimmed_reads/trimmed_SRR391537_R1.fastq -2 ./trimmed_reads/trimmed_SRR391537_R2.fastq  -S ./mapping/SRR391537_R1.sam # output .sam file

## batch
hisat2 \
	-p 8 \ 
	-x $INDEX \ 
	-1 $INDIR/${SAMPLE}_trim_1.fastq.gz \
	-2 $INDIR/${SAMPLE}_trim_2.fastq.gz | \
samtools view -@ 8 -Sh -u - | \
samtools sort -@ 8 -T $SAMPLE - >$OUTDIR/$SAMPLE.bam
```
#### 4.2 De novo assembly 
##### 4.2.1 Assembly of transcriptomes
To identify which reads were derived from which genes when there are no reference genome or transcriptome

```bash
# assembly
SAM=K21
Trinity --seqType fq \
  --left ./trim_${SAM}_R1.fastq.gz \
  --right ./trim_${SAM}_R2.fastq.gz \
  --min_contig_length 300 \
  --CPU 8 \
  --max_memory 16G \
  --output trinity_${SAM} \
  --full_cleanup 

# add a sample ID prefix to each sequence name
SAM=K21
sed "s/>/>${SAM}_/g" ./trinity_${SAM}.Trinity.fasta > ./trinity_prefix_${SAM}.Trinity.fasta

# concatenate the assemblies into a single file
cat ./trinity_prefix_* > ./trinity_combine.fasta
```
##### 4.2.2 Identifying coding regions
1.Identify long open reading frames (ORFs). Many transcripts will have multiple ORFs that could correspond to a true coding sequence. 
```bash
# 1.identify long open reading frames
TransDecoder.LongOrfs -m 100 \ # at least 100 amino acids
    -t ./trinity_combine.fasta 
```
2.Identify ORFs with homology to known proteins using `hmmer`. 
Many transcripts may have multiple ORFs, which bring in evidence of homology to known proteins for the final predictions. 
```bash
# 2.identify ORFs with homology to known proteins
hmmscan --cpu 16 \
        --domtblout pfam.domtblout \
        Pfam/Pfam-A.hmm \  # very large
        trinity_combine.fasta.transdecoder_dir/longest_orfs.pep
```

3.Final predictions of real ORFs in these transcripts are 

```bash
TransDecoder.Predict -t ../Assembly/trinity_combine.fasta \
        --retain_pfam_hits pfam.domtblout \
        --cpu 16
```

##### 4.2.3 Determining and Removing Redundant Transcripts
cluster the transcripts with similar sequences and then choose one representative transcript (the centroid)
```bash
vsearch --threads 8 --log LOGFile \
        --cluster_fast ../Coding_Regions/trinity_combine.fasta.transdecoder.cds \
        --id 0.90 \
        --centroids centroids.fasta \
        --uc clusters.uc
```
#### 4.3 Reference guided transcript assembly
Stringtie is a fast and highly efficient assembler of RNA-Seq alignments into potential **transcripts**.
It can be executed in 3 different modes:
1. Exclusively reference guided : quantify the expression of known transcripts only.
2. Reference guided transcript discovery mode : quantify known transcripts and detect novel ones.
3. De-novo mode : Detect and assemble transcripts.

```bash
# generate the gtf file has information on expression levels of transcripts, exons and other features along with any novel transcripts.
stringtie -p 4 \ # threads number
          -l label \ # label used in gtf file
          -G Reference.gtf \ # Reference GTF
          -o sample.gtf \  # output gtf of sample
          sample.bam # result file of sample genome alignment, must be sorted

# Redundant transcripts across the samples should be represented once Known transcripts should hold their stable gene ID's (assigned in Ensembl)
ls -1 ./*.gtf >> sample_gtf_list.txt
stringtie --merge -p 4 -o stringtie_merged.gtf -G Reference.gtf sample_gtf_list.txt

# Compare merged GTF with Reference GTF
gffcompare -r Reference.gtf -o gffcompare stringtie_merged.gtf

# Transcript quantification
stringtie -e \ # only estimate the abundance of given reference transcripts
          -B \ # returns a Ballgown input table file
          -p 4 \ # threads
          sample.sorted.bam \ # result file of sample genome alignment
          -G stringtie_merged.gtf \ # merged_sample GTF file
          -o output.count \ # output path/file name
          -A gene_abundance.out # gene abundance estimation output file
```
### 5 Evaluating the Assembly
#### 5.1 samtools stat
```bash
## stat
cat $SEQLIST | parallel -j 12 \
	"samtools stats $INDIR/{}.bam >$OUTDIR/{}.stats"

## glimpse
grep "^SN" samtools_stats/SRR12475447.stats

```
#### 5.2 Qualimap
We need to provide the GTF annotation file so qualimap can count up how many reads map to annotated features.
```bash
# GTF annotation file
GTF=

#
cat $SEQLIST | \
parallel -j 5 \
    qualimap \
        rnaseq \
        -bam $INDIR/{}.bam \
        -gtf $GTF \
        -outdir $OUTDIR/{} \
        --java-mem-size=2G  
```

#### 5.3 rnaQUAST
```
rnaQUAST.py --transcripts ../05_Clustering/centroids.fasta \ # reference transcriptome
	--gene_mark \
  --threads 8 \
  --output_dir Genemark
```


### 6.Quantify gene expression
1.htseq  
1. -s: whether strand specific counts
2. -r: indicates the order was by alignment position
3. -t: annotation file feature
4. -i: attribute from the annotation file.
5. -f: input file types 

```bash
# single file
htseq-count \
  -s no \  # unstranded RNA-seq library.
  -r pos \ # BAM file is coordinate sorted
  —t exon \
  -i pacid \
  -f bam sample_sorted.bam \
  ./alignment/Gmax_275_Wm82.a2.v1.gene_exons > sample.counts

# multiple files
cat $SEQLIST | \
parallel -j 5 \
    "htseq-count \
        -s no \
        -r pos \
        -f bam $INDIR/{}.bam \
        $GTF \
        > $OUTDIR/{}.counts"
```

2.FeatureCount
> The fastest
```bash
#
featureCounts -p \ # paired
              -a ./genome.gtf\ # annotation .gtf file
              -o gene_counts.txt \ # output count results
              -T 4 \ # threads number
              -t exon \ # rna type
              -g gene_id \ # rowname
              sample*.bam # all input files

#
```

3.RSEM
Alignment-based **transcript** quantification
```bash
# build index using bowtie2
rsem-prepare-reference \
      -gtf annotation.gtf \
   		--bowtie2 \
      reference_genome.fa ./RSEM_index/

# align and quantify transcripts
rsem-calculate-expression \
	  -p 6 \
	  --bowtie2 \
	  --append-names \ # append name of gene and transcripts
	  --output-genome-bam \
	  --paired-end sample_R1.fastq sample_R2.fastq \
	  ./RSEM_index/ \
	  RSEM/sample_count

# only align
bowtie2 -q \
	  --phred33 \
	  --sensitive \
	  --dpad 0 \
	  --gbar 99999999 \	# omit gap
	  --mp 1,1 \
	  --np 1 \
	  --score-min L,0,-0.1 \
	  -I 1 -X 1000 --no-mixed \
	  --no-discordant \	# omit paired reads discordant
	  -p 6 \
	  -k 200 \	# output the top 200 match results
	  -x RSEM_index/ \
	  -1 sample_R1.fastq \
	  -2 sample_R2.fastq | samtools view   -S -b -o RSEM_align/sample.bam

# only quantify transcripts
rsem-calculate-expression --alignments
```

4.kallisto  
kallisto is fast, and referred to as a **"pseudo-mapper"** because it have **no base-level alignment** of reads in favor of finding reads' approximate position in the reference transcriptome.

```bash
# index reference transcriptome
kallisto index -i ./centroids.fasta.index ./centroids.fasta

# Counting reads mapping to transcripts
kallisto quant \
  -i ./centroids.fasta.index \
  -o sample1 \
  -t 8 \
  ./trim_sample1_R1.fastq.gz ./trim_sample1_R2.fastq.gz
```

5.Edge-pro  
Edge-pro can estimate gene expression levels in prokaryotic genomes from RNA-seq data. It intergrate functions of aligner and counter.
1. Index the reference genome and align the reads to this indexed reference
2. calculate RPKM values based on these alignments.  
3. Edge-Pro using Bowtie as its aligner and it is optimizing the alignments for genes that do not contain introns.

```bash
edge.pl -g ../reference_genome/NC_003210.fna \ # prokaryotic reference genome fasta file
        -p ../reference_genome/NC_003210.ptt \ # ptt file with coordinates of coding genes, in Genbank format
        -r ../reference_genome/NC_003210.rnt \ # rnt file with coordinates of rRNAs and tRNAs, in Genbank format
        -u ../sickle_quality_control/ SRR034450_trimmed.fastq \ #  the trimmed fastq file from sickle or other tools
        -o SRR034450.out \ # output files prefix
        -s /EDGE_pro/1.3.1 \ # Default: working directory
        -t 8

# convert the output of EDGE-pro to DESeq
edgeToDeseq.perl SRR034450.out.rpkm_0

# remove the duplicate row (gene) with the smallest total count
python trim_epro2deseq.py deseqFile
```