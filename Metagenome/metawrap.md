
## Quality control
```bash
mkdir read_qc
for fq in *_1.fastq
do
ID = $(basename $fq _1.fastq)
metawrap read_qc -1 ${ID}_1.fastq -2 ${ID}_2.fastq -t 24 -o read_qc/$ID
done

cat *_1.fastq > merged_1.fastq
cat *_2.fastq > merged_2.fastq
```

## Assembly
```bash
metawrap assembly -1 read_qc/merged_1.fastq -2 read_qc/merged_2.fastq -m 200 -t 96  -o assembly_contigs
grep ">" assembly_contigs/final_assembly.fasta | head -n3 
```

## Taxonomy annotation
```bash
metawrap kraken -o Kraken_contigs -t 96 -s 1000000 read_qc/s*.fastq assembly_contigs/final_assembly.fasta
```

## Binning based on metabat2, maxbin2, concoct
```
metawrap binning -o Binning_contigs -t 8 -a assembly_contigs/final_assembly.fasta --metabat2 --maxbin2 --concoct read_qc/s*.fastq
```

## Aggregate results of three binning
```bash
#
metawrap Bin_refinement -o bin_refine -t 8 -A Binning_contigs/metabat2_bins/ -B Binning_contigs/maxbin2_bins/ -C Binning_contigs/concoct_bins/ -c 50 -x 10
#
cat bin_refinement/metaWRAP_bins.stats
metawrap blobology -a assembly_contigs/final_assembly.fasta -t 8 -o biology_plot --bins bin_refinement/metaWRAP_bins read_qc/s*.fastq
```

## Quantification
```bash
#
metawrap quant_bins -b bin_refinement/metaWRAP_bins -t 8 -o Bin_quantification -a assembly_contigs/final_assembly.fasta read_qc/s*fastq
```

## Reassembly

```
metawrap reassemble_bins -o Bin_reassembly -1 read_qc/merged_1.fastq -2 read_qc/merged_2.fastq -t 8 -m 800 -c 50 -x 10 -b bin_refinement/metaWRAP_bins
```

## Bin Taxonomy annotation
```bash
#
metawrap classify_bins -b Bin_reassembly/reassembled_bins -o Bin_taxonomy -t 8
```

## Bin Gene annotation
```
metaWRAP annotate_bins -o Bin_function -t 8 -b Bin_reassembly/reassembled_bins/
```