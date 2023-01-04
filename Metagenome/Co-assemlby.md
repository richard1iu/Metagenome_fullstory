## What is a co-assembly?
**Co-assembly** refers to performing an assembly for multiple samples simultaneously. This is in contrast to doing an independent assembly for each sample. 

### benefits of co-assembly: 
1) higher read depth; 
2) comparison across samples by using one reference assembly for all samples; 
3) differential coverage.

###
compare them with QUAST (for individual genome) or MetaQUAST (for metagenome assemblies). 

## Co-assembly for multiple samples
```bash
# integrate samples
R1s=`ls *_1.* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
R2s=`ls *_2.* | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`

# co-assembly for multiple samples
megahit -1 $R1s -2 $R2s --min-contig-len 1000 -t 4 -m 0.85 -o megahit_merge/ 

# glimpse the results
tail megahit_merge/log
```

## simplified headers and filter out contigs
```bash
# activate anvio env
mamba activate anvio
# format
anvi-script-reformat-fasta megahit_merge/final.contigs.fa -o contigs.fa -l 1000 --simplify-names --report name_conversions.txt
```

## Mapping
```bash
# 1.create an index 
bowtie2-build contigs.fa assembly

# 2.mapping
for sample in $(cat samples.txt)
do
## mapping samples to index
bowtie2 -x assembly -q -1 "$sample"_1.fastq.gz -2 "$sample"_2.fastq.gz --no-unal -p 4 -S "$sample".sam
samtools view -b -o "$sample".bam "$sample".sam
## sort and indexed bam
samtools sort -o "$sample".sorted.bam "$sample".bam
samtools index "$sample".sorted.bam
done
```

## Build up our contigs database
Generate an anvi’o contigs database from our co-assembly fasta file:
1. calculating tetranucleotide frequencies for each contig (uses 4-mers by default); 
2. identifies open-reading frames (“genes”) with prodigal;
3. splits long contigs into segments of roughly 20,000 bps

```bash
anvi-gen-contigs-database -f contigs.fa -o contigs.db -n "merged_metagenome"
```

## Search bacterial single-copy genes
Using the program HMMER to scan for a commonly used set of bacterial single-copy genes.
This will help us estimate genome completeness/redundancy in real-time.
```bash
# glimpse the available HMM profile
anvi-run-hmms --help

# hmmer scan
anvi-run-hmms -c contigs.db -I Campbell_et_al -T 4
```

## Functional annotation
Use NCBI COGs for functional annotation of the open-reading frames prodigal predicted. 
either BLAST or DIAMOND – DIAMOND is like a less sensitive, but faster (default is DIAMOND). 

```bash
# set up cogs database
anvi-setup-ncbi-cogs -T 4 --just-do-it

# cogs functional annotation with diamond
anvi-run-ncbi-cogs -c contigs.db  -T 4
```

## Taxonomy annotation
Assign taxonomy with a tool called Centrifuge to the open-reading frames prodigal predicted. 

```bash
#
anvi-get-sequences-for-gene-calls -c contigs.db -o gene_calls.fa

#
centrifuge -f -x centrifuge_db/nt/nt gene_calls.fa -S centrifuge_taxo.tsv -p 4

#
anvi-import-taxonomy-for-genes -c contigs.db -i centrifuge_report.tsv centrifuge_taxo.tsv -p centrifuge
```

## Profiling our samples
provide information about each sample to anvi’o so it can then integrate everything together. 
Each sample will have what’s known as a “profile database” that will keep information about that sample 

```bash
# profile each samples to contigs.db
for i in $(cat samples.txt)
do 
anvi-profile -i "$i".bam -c contigs.db -T 4
done

# merge all
anvi-merge */PROFILE.db -o merged_profile -c contigs.db
```
