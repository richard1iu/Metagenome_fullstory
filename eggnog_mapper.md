
## options
Input Data Options:
  -i FASTA_FILE         Input FASTA file containing query sequences (proteins by default;). Required unless `-m no_search`
  --itype               Type of data in the input (-i) file. Default: proteins, including: {CDS,proteins,genome,metagenome}
                        
  --translate           When --itype CDS, translate CDS to proteins before search. 
                        When --itype genome/metagenome and --genepred search, translate predicted CDS from blastx hits to proteins. (default: False)    
  --annotate_hits_table SEED_ORTHOLOGS_FILE. a .seed_orthologs file from a previous emapper.py run. Requires `-m no_search.` (default: None)
                        Annotate TSV formatted table with 4 fields: query, hit, evalue, score.
  --data_dir DIR        Path to eggnog-mapper databases. By default, "data/" or the path specified in the environment variable EGGNOG_DATA_DIR. (default: None)
optional arguments:
  --temp_dir DIR        temporary files directory.(default: current directory)
  --no_file_comments    No header lines nor stats are included in the output files (default: False)
  --resume              continue corrupt jobs
  --override            rerun corrupt jobs

## workflow of eggnog-mapper : 
1) Search orthologous sequences; CPU calculate 
2) Function annotation; Disk read and write

```bash
# For large fasta file，split it for multiprocess run to speed up
split -l 2000000 -a 3 -d input_file.faa input_file.chunk_

# 1.Search orthologous sequences
## generate run codes for cluter; --no_annot: don't annotation
for f in *.chunk_*; do echo emapper.py -m diamond --no_annot --no_file_comments --cpu 16 -i $f -o $f; done 

# 2. Function annotation
## combine all seed_orthologs files
cat *.chunk_*.emapper.seed_orthologs > input_file.emapper.seed_orthologs
## annotation
emapper.py --annotate_hits_table input.emapper.seed_orthologs --no_file_comments -o output_file --cpu 10
```


## one-step mapping
```
EGGNOG_DATA_DIR=/mnt/f/Database/eggnog_db/
emapper.py \
-i data/paried/alfafa_zm4_reference_genome/zm-4.protein.fasta  \
--temp_dir data/ \
--data_dir $EGGNOG_DATA_DIR \
-m diamond \
--dbmem \
--tax_scope Viridiplantae \
-o alfafa \
--output_dir results/annotation
```

# modify the annotation for R
```
id=alfafa # -o
sed '/^##/d' ${id}.emapper.annotations| sed 's/#//g'| awk -vFS="\t" -vOFS="\t" '{print $1,$9,$10,$12}' > ${id}.annotations
```
