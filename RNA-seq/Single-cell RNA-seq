# Single-cell RNA sequencing (Cell Ranger)
The **cellranger** `mkfastq` pipeline is a 10x-enhanced wrapper around Illumina bcl2fastq, which demultiplexes BCL files from a sequencer into FASTQs for analysis.

## Download data (optional)
```bash
# sequences
wget http://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
tar -xvzf cellranger-tiny-bcl-1.2.0.tar.gz
# sample table
wget http://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-samplesheet-1.2.0.csv
```
## 1.Generate fastqs
```bash
cellranger mkfastq --run=./tiny_bcl \ # input bcl file
                   --samplesheet=tiny-bcl-samplesheet.csv
```

## 2.Count
cellranger count takes FASTQ files from cellranger mkfastq and performs alignment, filtering, and UMI counting.
```bash
cellranger count --cells=100 \
    --id=hgmm \ # a unique run ID string
    --transcriptome=./ refdata-cellranger-hg19-and-mm10-1.2.0 \ # path to the Cell Ranger compatible transcriptome reference
    --fastqs=./hgmm100/fastqs \ # path of the fastqs directory
    --sample=read # sample name as specified in the sample sheet
```