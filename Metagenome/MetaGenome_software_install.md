# Quality control
- Fastq quality checking: Fastqc
- Checking result integrating: MultiQC
- Filter adapter and low-quality reads: Trimmoatic
- Filtering low abundance kmer: khmer
- 
```bash
sudo apt install khmer
export PATH=/usr/lib/khmer/bin:$PATH
```

# Genome assembly
- Megahit
- Metaspades 
- Minia

```bash
conda install megahit
conda install spades
```
# Quality of assembly evaluation
- Sourmash
- QUAST
```
conda install sourmash
```


# Genome annotation
- Prokka
- Prodigal

```
conda install prokka=1.13
conda install prodigal=2.6.3
```
# Gene abundance estimation 
## without alignment
salmon
conda install salmon
## with alignment
bowtie2, samtools, bedtools
conda install bowtie2
conda install bedtools
conda install samtools

# Binning
Maxbin、MetaBAT、MetaWatt、CONCOCT、MyCC

```bash
Maxbin
curl  https://downloads.jbei.org/data/microbial_communities/MaxBin/getfile.php?MaxBin-2.2.2.tar.gz > MaxBin-2.2.2.tar.gz
tar xzvf MaxBin-2.2.2.tar.gz
cd MaxBin-2.2.2/src
make

# MetaBAT
Curl -L https://bitbucket.org/berkeleylab/metabat/downloads/metabat-static-binary-linux-x64_v0.32.4.tar.gz > metabatv0.32.4.tar.gz
tar xvf metabatv0.32.4.tar.gz
```

# Binning evaluation
checkm
conda create -n checkm checkm

# Taxa abundance estimation
Metaphlan、Kraken
```bash
# Metaphlan
wget https://bitbucket.org/biobakery/metaphlan2/get/default.zip
tar xzvf biobakery-metaphlan2-<versioned>.tar.gz
cd biobakery-metaphlan2-<versioned>/

# Kraken
conda create -n kraken=1.0
kraken db 
wget -c https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz
```

# ENA data download
enaBrowerTools: enaDataGet; enaGroupGet
https://github.com/enasequence/enaBrowserTools/releases/tag/v1.1.0
```
alias engDataGet=/home/bioinfo/ enaBrowserTools-1.1.0/python3/enaDataGet 

alias engGroupGet=/home/bioinfo/ enaBrowserTools-1.1.0/python3 enaGroupGet                                                                    
```