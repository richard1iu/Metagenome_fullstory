## 1.import
### 1.1 paired end
```bash
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path all_sample.txt \
--output-path paired-demux.qza \
--input-format PairedEndFastqManifestPhred33
```

### 1.2 single end
```bash
qiime tools import --type 'SampleData[SequencesWithQuality]' \
--input-path ./xinjiang120 \
--output-path paired-demux.qza \
--input-format SingleEndFastqManifestPhred33
```

## 2.view quality
```bash
qiime demux summarize --i-data paired-demux.qza --o-visualization paired-demux.qzv
```

## 3.dada2
### chimera filtering, generating ASV,rep-seqs
```bash
qiime dada2 denoise-paired --i-demultiplexed-seqs paired-end-demux.qza \
--p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 270 --p-trunc-len-r 250 \
--o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza \
--p-n-threads 0 # all threads
```

### stat (optianal)
```bash
# ASV stat
qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

# Feature stat
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file sample-metadata.txt 

# rep-seq stat
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
```

## 4.tree
```bash
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza 
```

## 5.classification

### 5.1 use other's classfier  
qiime feature-classifier classify-sklearn \
--i-classifier Greengenes_13_8_99%_OTUs_341F-805R_classifier.qza \
--i-reads rep-seqs.qza \
--o-classification taxonomy.qza 

### 5.2 train classfier with own data
```bash
# 1.download dataset (greengenes or silva) 
wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz 
tar -zxvf gg_13_8_otus.tar.gz 

# 2.load taxonomy in dataset
qiime tools import \ 
  --type 'FeatureData[Taxonomy]' \ 
  --input-format HeaderlessTSVTaxonomyFormat \ 
  --input-path 99_otu_taxonomy.txt \ 
  --output-path ref-taxonomy.qza 
 
# 3.Extract reference reads
qiime feature-classifier extract-reads \
--i-sequences 99_otus.qza \
--p-f-primer CCTACGGGNGGCWGCAG \ # 341F primer 
--p-r-primer GACTACHVGGGTATCTAATCC \ # 805R primer 
--o-reads ref-seqs.qza 
 
# 4.Train the classifie
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy ref-taxonomy.qza \
--o-classifier Greengenes_13_8_99%_OTUs_341F-805R_classifier.qza 

# 5.classification
qiime feature-classifier classify-sklearn \
--i-classifier Greengenes_13_8_99%_OTUs_341F-805R_classifier.qza \
--i-reads rep-seqs.qza \
--o-classification taxonomy.qza 
```

