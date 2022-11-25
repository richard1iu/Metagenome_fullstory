#https://codeantenna.com/a/FVqG01T9qL
#https://www.jianshu.com/p/c72bb359f050

# 下载RDP训练集
https://www.mothur.org/wiki/RDP_reference_files#Version_14

# 训练集小写转成大写
cat trainset16_022016.pds.fasta | tr a-z A-Z > trainset16.fasta

# install vsearch, usearch and mothur

## 1.unzip files
`gzip -d *.gz`

## rename fq
```bash
rename 's/fq/fastq/' *.fq
rename 's/\.1\.fq/\_R1\.fastq/'  *.fq
rename 's/\.2\.fq/\_R2\.fastq/' *.fq
rename 's/raw\.split\.//' *.fastq
```

## 2.merge paried reads (usearch)
```bash
usearch -fastq_mergepairs *R1*.fastq -relabel @ -fastq_maxdiffs 10 \
  -fastq_pctid 80 -fastqout merged.fq
```

## 3.remove primers and barcode reads (optional)
```bash
# 1.序列随机抽样 (optional)
vsearch -fastx_subsample merged.fq -sample_size 5000 -fastqout sub5000.fq

# 2.检索引物位置 (optional)
usearch -search_oligodb sub5000.fq \
-db primers.fa \
-strand both \
-userout primer_hits.txt \
-userfields query+qlo+qhi+qstrand

# 3.查看引物位置 (optional)
head primer_hits.txt
```

## 4.quality control (using vsearch if the size of your merged.fq > 4GB)
### (1) 16s & 18s
1. stat the reads length and quality
```bash
# glimpse the length of reads (optional)
vsearch -fastq_eestats2 merged.fq -output eestats2.txt -length_cutoffs 50,500,10 -fastq_qmax 44

# glimpse the quality distributions of reads (optional)
vsearch -fastq_eestats2 merged.fq -output stats -fastq_qmax 44
```

2. filter reads based on error rates and the position of primers
```bash
vsearch -fastq_filter merged.fq \ # input merged.fa
--fastq_stripleft 19 --fastq_stripright 20 \ # cut primers
--fastq_maxlen 280 --fastq_minlen 220 \ # the max and min length of reads
--fastq_qmax 60 \ # the max q value
-fastq_maxee_rate 0.01 -fastq_maxns 0 -fastq_maxee 1 \ # the max error rate
--fastaout filtered.fa # the output
```

1. filter reads based on error rates (without primers or low quality base)
```bash
vsearch -fastq_filter merged.fq \
--fastq_stripleft 0 --fastq_stripright 0 \
--fastq_maxlen 280 --fastq_minlen 220 \
--fastq_qmax 44 \
-fastq_maxee_rate 0.01 -fastq_maxns 0 -fastq_maxee 1 \
--fastaout filtered.fa
```

### (2) ITS filter based on error rate
```bash
vsearch -fastq_filter merged.fq \
-fastq_qmax 44 \
-fastq_maxee 1 -fastq_maxns 0 -fastq_maxee_rate 0.01 \
-fastaout filtered.fa
```

## 5.dereplicate randunry reads 
```bash
# 1.vsearch
vsearch -derep_fulllength filtered.fa \
-minuniquesize 8 \
-sizein -sizeout \
-fasta_width 0 \
-uc derep.uc \
-output derep.fa

# 2.usearch (optional)
usearch -fastx_uniques filtered.fa \
-sizeout \
-relabel Uniq \
-fastaout uniques.fa
```

## 6.filter chimers (optional)
The step `unoise` (7.2) does chimers filtering automatically
### 6.1 Abundance sorting and singletons checking 
```bash
usearch -minsize 2 -sortbysize derep.fa -fastaout sorted.fa 
```

### 6.2 uchime_denove (without reference)
```bash
# vsearch
vsearch -uchime_denovo derep.fasta \
-sizein -sizeout \
-fasta_width 0 \
-nonchimeras otu_vsearch.fa

# usearch
# The input to uchime3_denovo must be denoised amplicons
usearch -uchime3_denovo sorted.fa \
-chimeras chimeras.fa \
-nonchimeras nonchimeras.fa
```

### 6.3 vsearch -uchime_ref  (with reference)
```bash
# prokaryotes
vsearch -uchime_ref filtered.fa \
-nonchimeras otu_vsearch.fa \
-db /home/kxf/Documents/database/SILVA_132_SSURef_Nr99_tax_silva.fasta 

# fungi
vsearch -uchime_ref filtered.fa \
-nonchimeras out_vsearch.fa \
-db /home/kxf/Documents/database/UNITEv6_sh_97_s.fasta 

# usearch
usearch -uchime_ref filtered.fa \
-db 16s_ref.udb \ # FASTA or UDB; SILVA for 16S or UNITE for ITS
-strand plus \ # only plus, not support paired
-mode sensitive \ # specific; balanced; high_confidence; sensitive
-uchimeout otu.txt \ # summary table
-chimeras ch.fa \ # predicted chimeras
-notmatched not.fa \ # sequences not matched to the database
-uchimealnout aln.txt \ # alignments sequences

# stat the number of reads without chimers
grep -c "^>" otu_vsearch.fa
```

## 7.cluster
### 7.1 cluter to OTU,  97.5%
```bash
usearch -cluster_otus otu_vsearch.fa -otus otus.fa -relabel OTU
```
### 7.2 cluster to ZOTU(ASV) with unoise3, 100% 
```bash
usearch -unoise3 derep.fa -zotus ZOTU.fa -minsize 9 
```

## 8.generate otu table
```bash
# 1.usearch
usearch -otutab merged.fq -otus otus.fa -otutabout otutab_raw.txt #OTU
usearch -otutab merged.fq -otus ZOTU.fa -otutabout otutab_raw.txt #OTU

# 2.vsearch
vsearch -usearch_global filtered.fa --db ZOTU.fa --id 0.99 --otutabout otutable_raw.txt

# 2.mothur
usearch -otutab filtered.fa -otus otus.fa -mothur_shared_out otutab_raw.txt -mapout map.txt -threads 12
```

## 9.construct tree
```bash
# 1.generate distance matrix
./usearch11 -calc_distmx otus.fa \
-tabbedout matrix.txt \
-maxdist 0.2 -termdist 0.3

# 2.construct tree using distance matrix
usearch -cluster_aggd matrix.txt \
-treeout clusters.tree \
-clusterout clusters.txt \
-id 0.80 -linkage min
```

## 10.annotation rep sequences
```bash
# 1.usearch
## 1.1 taxonomy predictions with bootstrap confidence
usearch -sintax otus.fa \
-db ../data/rdp_16s_v2.fa \
-strand both \
-tabbedout sintax.txt \
-sintax_cutoff 0.8

## 1.2 summary
usearch -sintax_summary sintax.txt -otutabin otutab.txt -rank p -output phylum_summary.txt

# 2.Mothur
Mothur >classify.seqs(fasta=otus.fa, 
reference=/home/kxf/Documents/database/silva.nr_v132.align, taxonomy=/home/kxf/Documents/database/silva.nr_v132.tax,
cutoff=60, processors=88) #细菌，原虫

mothur >classify.seqs(fasta=otus.fa,
reference=/home/kxf/Documents/database/RIM_DB_14_07.fasta, taxonomy=/home/kxf/Documents/database/RIM_DB_14_07_c.txt, processors=88)  #产甲烷菌

mothur >classify.seqs(fasta=otus.fa, 
reference=/home/kxf/Documents/database/UNITEv6_sh_97_s.fasta, taxonomy=/home/kxf/Documents/database/UNITEv6_sh_97_s.tax, processors=88)  #真菌ITS
```

## 11.normalize
```bash
# (1) normalize to the least number of reads
normalize.shared(shared=otutab_raw.txt)

# (2) normalize to specific number of reads / sample
usearch -otutab_norm otutab_raw.txt -sample_size 5000 -output otutab.txt
```

## 12.rarefaction (R programming)
```bash
rarefaction.single(shared=otutab.txt)
usearch -alpha_div_rare otutab.txt -output rare.txt
```

## 13.calculate α-diversity
```
usearch -alpha_div otutab.txt -output alpha.txt #32位不够
summary.single(shared=current)
```

## 14.calculate β-diversity (R programming)
```R
#distance matrix
dist.shared(shared=current, calc=thetayc-jclass-braycurtis) 
count.seqs(shared=current)

unifrac.unweighted(tree=clusters.tree, count=current, distance=lt, processors=88, random=F)

unifrac.weighted(tree=clusters.tree, count=current, distance=lt, processors=88, random=F)

# pcoa
pcoa(phylip=clusters.tree1.unweighted.phylip.dist)
pcoa(phylip=clusters.tree1.weighted.phylip.dist)
pcoa(phylip=otutab.braycurtis.$usearch.lt.dist)
pcoa(phylip=otutab.jclass.$usearch.lt.dist)
pcoa(phylip=otutab.thetayc.$usearch.lt.dist)
system(mkdir -p result/pcoa/)
system(rm -f *.rabund)
system(mv *.dist *.loadings *.axes -t result/pcoa/)
system(mv *.tree otus.fa *.rarefaction otutab.txt otutab.$usearch.count_table otutab.groups.summary *.taxonomy  -t result/)
```