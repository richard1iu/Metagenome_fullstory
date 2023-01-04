

## Script to run VSEARCH UNOISE3 workflow for 16S V4 reads on IMAC (desktop)
- Denoising: VSEARCH (version 2.15.2) UNOISE3
- Taxonomic classification: GG13_8 Naive Bayesian classifier
- Phylogenetic tree: SEPP

## Program Assumes:
 - demultiplexed F and R fastq files are placed in "seqs" folder folder on desktop
 - files are compressed (.gz)
 - primers and all non-biologic sequences have been removed
 - files end in "R1_001.fastq.gz" format
 - use following if need to modify: for file in *; do mv "$file" ${file//_R1/_R1_001}; done



## vsearch workflow
```bash
# 1.Merge reads
for R1 in *_R1_*.fastq.gz; do
    vsearch \
        --fastq_mergepairs ${R1} \
        --reverse ${R1/_R1_/_R2_} \
        --fastqout ${R1/_R1_*/merged.fq} \
        --fastaout ${R1/_R1_*/merged.fa} \
        --fastq_maxdiffs 10 --fastq_minmergelen 240 --fastq_maxmergelen 270 \
        --relabel ${R1} \
        --threads 14
done

cat *fa > my.fasta
cat *fq > my.fastq
rm *.fa
rm *.fq

sed 's/fastq.gz//g' < my.fasta > merged.fasta
sed 's/fastq.gz//g' < my.fastq > merged.fastq
rm my.fasta
rm my.fastq


# 2.Filter low quality reads
vsearch --fastq_filter merged.fastq --fastq_maxee 1.0 --fastaout filtered.fa --relabel filt --log filter_log.txt

# 3.Get unique reads
vsearch  --derep_fulllength filtered.fa --sizein --fasta_width 0 --sizeout --output uniques.fa --minuniquesize 1 --relabel unique --log uniques_log.txt

# 4.Cluster zOTUs
vsearch --cluster_unoise uniques.fa --centroids zotus.fa --uc uc_zotus.uc --log unoise_log.txt --threads 14

# 5.Remove chimeras
vsearch  --uchime3_denovo zotus.fa --nonchimeras zotus_nochime.fa --fasta_width 0 --relabel zotu --xsize --log unchime3_log.txt

# 6.Removing any masking (convert any lower to upper case)
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' zotus_nochime.fa > zotus_nochime_nomask.fa

# 7.Map reads back to centriods
vsearch --usearch_global merged.fasta --db zotus_nochime_nomask.fa --id 0.97 --otutabout unoise3_zotu_table.txt --biomout unoise3_zotu_table.biom  --threads 14 --log otutab_log.txt

# 8.Move log files
mkdir logs
mv *log.txt logs
```

## Qiime2 workflow
```bash
# 1.Load QIIME2
source activate qiime2-2022.8

# 2.Convert zOTU table and rep-seqs to q2 artifacts
mkdir q2
mkdir q2/reports

qiime tools import \
--input-path zotus_nochime_nomask.fa \
--output-path q2/rep-seqs.qza \
--type "FeatureData[Sequence]"

qiime tools import \
  --input-path unoise3_zotu_table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path q2/table.qza

# 3.Assign taxonomy using pre-trained (515F-806R) GG ref db
cd q2

qiime feature-classifier classify-sklearn \
  --i-classifier /Users/olljt2/Documents/q2_db/gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza \
  --p-confidence .8 \
  --p-n-jobs -14

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization reports/taxonomy.qzv

# 4.Fragment insertion via SEPP for tree
qiime fragment-insertion sepp \
  --i-representative-sequences rep-seqs.qza \
  --i-reference-database /Users/olljt2/Documents/q2_db/sepp-refs-gg-13-8.qza \
  --o-tree tree.qza \
  --o-placements insertion-placements.qza \
  --p-threads 14
```
 