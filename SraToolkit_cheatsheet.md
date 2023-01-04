

## install
### 1.sra toolkit
```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz
tar zxvf sratoolkit.2.9.6-ubuntu64.tar.gz
cd sratoolkit.2.9.6-ubuntu64
#加入环境路径
echo 'export PATH=$PATH:your saved path/sratoolkit.2.9.6-ubuntu64/bin' >> ~/.bashrc
source ~/.bashrc
```
### 2. aspera high speed download (paid)

### 3. wget/aria2
```bash
# ubuntu
sudo apt install aria2
```

## download SRA
### prefetch (file size < 20GB)
```bash
# single file
prefetch SRR8956146

# multiple file
ls * > sra_list.txt
prefetch --option-file sra_list.txt

# file path
sudo updatedb
locate SRR8956146.sra
```

### fastq-dump
#### options
1. --split-3: split paired reads 
2. --fasta: output fasta
3. --gzip, --bzip2: output zip file 
4. -O|--outdir: specific the output dir
5. --defline-seq: define the format of readsID
6. --defline-qual: define the format of reads quality

```bash
alias fd='fastq-dump --split-3 --defline-qual '+' --defline-seq '@\\\$ac-\\\$si/\\\$ri' '
fd SRA_ID
```

### fasterq-dump
#### options
1. -e|threads: the number of threads
2. -p: current progress
3. -O: specific the output dir
```
fasterq-dump --split-3 --mem 16 --threads 8 ./SRR5318040 
for i in s{10..20};do fasterq-dump $i --split-3 --mem 32G -e 8 --qual-defline '+' --seq-defline '@$ac-$si/$ri'; done
```
