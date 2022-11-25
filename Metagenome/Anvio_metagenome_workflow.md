
## Fix your deflines
--------------------
**defline** (description line) is distinguished from the sequence data by a **greater-than (">")** symbol at the beginning

**NOTICE:** The names in contigs.fa must match the names in your BAM files when you are **mapping**.

Check your deflines prior to mapping, remove anything that is **not a digit, an ASCII letter, an underscore, or a dash character.**

**Bad deflines**
>Contig-123 length:4567
>Another defline 42
>gi|478446819|gb|JN117275.2|

**Elegant deflines**
>Contig-123
>Another_defline_42
>gi_478446819_gb_JN117275_2


```bash
# re-formatting our input FASTA 
anvi-script-reformat-fasta assembled-contigs.fa -o anvio-contigs.fa --min-len 2000 --simplify-names --report name_conversions.txt

# map our reads against the assembly
bowtie2-build anvio-contigs.fa anvio-contigs
bowtie2 -x anvio-contigs -1 tara_trimmed_R1.fq -2 tara_trimmed_R2.fq | samtools view -bS -o anvio-contigs.bam
samtools sort anvio-contigs.bam -o anvio_sorted.bam
samtools index anvio_sorted.bam
```

## create contigs database
----------------------------

#### This database contain:
- positions of open reading frames
- k-mer frequencies for each contigs
- where splits start and end
- functional and taxonomic annotation of genes, etc 

#### what they do
1. **split long contigs** (>20k): Anvi’o shows the generalized statistics for each contig (GC content etc.). For long contigs these stats are calculated across split contigs (which remain grouped)
2. Identify and locate **open reading frames** in your contigs (using Prodigal)
3. Estimate **Single Copy Gene content** (using hmmer against defined gene sets for bacteria and archaea)
4. **Calculate k-mer frequencies** for the contigs in our assemblies

#### Decorate contigs database with HMM
It will utilize multiple default **bacterial single-copy** core gene collections and identify hits among your genes to those collections using HMMER.

#### stats of contigs database
```bash
# generate database
anvi-gen-contigs-database -f anvio-contigs.fa -o anvio-contigs.db

# with default bacterial and archaea database
anvi-run-hmms -c anvio-contigs.db --num-threads 8

# --installed-hmm-profile: run a specific default HMM profile 

# --hmm-profile-dir: declare where your own .hmm are

# assess your assembly output and the number of bacterial and archaeal genomes to recover
anvi-display-contigs-stats anvio-contigs.db --server-only -P 8080
# type the address http://localhost:8080 in your own brower rather than in remote server
``` 

## Profile the bam
The profiling step explain each BAM file separately by using the information stored in the contigs database. It is one of the most critical (and also most complex and computationally demanding) steps of the metagenomic workflow

```bash
# sorted and indexed .bam (skip this step if you have done it with samtools)
anvi-init-bam tara_unsorted.bam -o tara_anvio_sorted.bam 
for bam in *.bam; do anvi-init-bam $bam -o ${bam/.bam/}_sorted.bam; done 

# see the coverge information
for file in *_sorted.bam; do anvi-profile -i $file -c anvio-contigs.db -T 8; done

# pull all information together, create a merged anvi’o profile and calculates the hierarchical relationship betwewen contigs
anvi-merge *ANVIO_PROFILE/PROFILE.db -o MERGED-SAMPLES -c anvio-contigs.db --enforce-hierarchical-clustering

# visualize our data
anvi-interactive -p MERGED-SAMPLES/PROFILE.db -c anvio-contigs.db
```

## Identifying and refining genome bins
---------------------------------------
adjust your bins that may be contaminated
```bash
# summarize the statistic information of bins
anvi-summarize -p MERGED-SAMPLES/PROFILE.db -c anvio-contigs.db -o SAMPLES-SUMMARY -C CONCOCT

# visualize it
anvi-interactive -p MERGED-SAMPLES/PROFILE.db -c anvio-contigs.db -C CONCOCT --server-only -P 8080

# use human intuition and pattern recognition to better identify contigs that should co-occur.
anvi-refine -p MERGED-SAMPLES/PROFILE.db -c anvio-contigs.db -b Bin_4 -C CONCOCT
```

## Function annotation

## Taxonomic annotation

#### workflow
1. generate your contigs database from your FASTA file
2. export your gene sequences
3. annotate them with taxonomy
4. import results back into your contigs database 

```bash
# 
anvi-get-sequences-for-gene-calls -c anvio-contigs.db -o gene_calls.fa
```