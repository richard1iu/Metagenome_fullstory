## Reference
[seqtk](https://github.com/lh3/seqtk)


## Install
```bash
conda install -c bioconda seqtk
```

## Convert .fq to .fa
```bash
seqtk seq -a reads.fq.gz > reads.fa
```

##  Randomly Subsample FASTQ
```bash
# Randomly Subsample single FASTQ
seqtk sample read_R1.fq 10000 > sub_R1.fq

# Randomly Subsample Paired FASTQ
# must set the same random seed to keep pairing
seqtk sample -s 123 read_R1.fq 10000 > sub_R1.fq
seqtk sample -s 123 read_R2.fq 10000 > sub_R2.fq
```
## Reverse complement FASTA/Q:
```bash
seqtk seq -r reads.fq > out.fq
```

## Extract sequences with names in file name.lst, one sequence name per line:
```bash
seqtk subseq reads.fq name.lst > out.fq
```