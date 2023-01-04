
## reference
1. Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3
2. eLife 2021;10:e65088

nucleotide mapping and translated search to provide organism-specific gene and pathway abundance profiles from a single metagenome.
gene families are annotated using **UniRef90** identifiers and pathways using **MetaCyc** IDs

## install
```bash
# install
conda install humann -c biobakery
# test it if work well
humann_test
```

## upgrade database
```bash
# pangenome database
humann_databases --download chocophlan full /path/to/databases --update-config yes
# protein database
humann_databases --download uniref uniref90_diamond /path/to/databases --update-config yes
# annotations database
humann_databases --download utility_mapping full /path/to/databases --update-config yes
```

## qualtiy check
```bash
# a wrapper of Trimmomatic and Bowtie2
kneaddata -i raw_data_example/p144C_R1.fastq -i raw_data_example/p144C_R2.fastq \
          -o kneaddata_out \
          -db /bowtiedb/GRCh38_PhiX \ # path to host sequence of bowtie2 database
          --trimmomatic /Trimmomatic-0.36/ \ # path to Trimmomatic folder containing jarfile
          --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
          --bowtie2-options "--very-sensitive --dovetail" \
          --remove-intermediate-output \
          -t 4 
```

## humann function annotation 
```bash
# single sample
humann --input demo.fastq --output demo_fastq --threads 8 --nucleotide-database $DIR

# glimpse the result table
column -t -s $'\t' demo_fastq/demo_genefamilies.tsv | less -S

# multiple samples 
for f in *.fasta; do humann -i $f -o hmp_subset; done
## merge all results to a table
humann_join_tables -i hmp_subset -o hmp_subset_genefamilies.tsv --file_name genefamilies
## normalize result to adjust for differences in sequencing depth across the samples.
humann_renorm_table -i hmp_subset_genefamilies.tsv -o hmp_subset_genefamilies-cpm.tsv --units cpm
```

## default output
1. genefamilies.tsv
2. pathabundance.tsv
3. pathcoverage.tsv

### 1.genefamilies.tsv
This file lists the abundances of each gene family in the community in **RPK (Reads Per Kilobase)** units
This file lists the abundances of each pathway in the community, also in RPK units 