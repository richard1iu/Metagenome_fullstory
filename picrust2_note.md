

## 1.empirical data output
```bash
picrust2_pipeline.py -s ZOTU_repseq.fa \
-i otutab_raw.txt \
-o picrust2_result \
-p 2

# add description columns of each functional category
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                    -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz
```


## 2.shuffle data output
It can be helpful to compare the empirical tables with tables based on shuffling the predictions. 

### 2.1.generate the EC table
```bash
# picrust version (v2.4.0)
shuffle_predictions.py -i EC_predicted.tsv.gz \
                           -o EC_predicted_shuffled \
                           -r 5 \ # how many random replicates to make
                           -s 131 # specifies a random seed
```

### 2.2. generate the pathway table
The gene family and pathway-level prediction tables can then be generated from these shuffled tables by running the standard PICRUSt2 commands. 

Below is an example of how to quickly run metagenome_pipeline.py and pathway_pipeline.py on all shuffled tables with a bash loop.

```bash

# Make folders for shuffled output
mkdir EC_metagenome_out_shuffled
mkdir pathways_out_shuffled

# script 
for i in {1..5}; do
    
    # Define in and out file paths.
    EC_SHUFFLED="EC_predicted_shuffled/EC_predicted_shuf"$i".tsv.gz"
    OUT_META="EC_metagenome_out_shuffled/rep"$i
    OUT_PATHWAYS="pathways_out_shuffled/rep"$i
    
    # PICRUSt2 scripts to get prediction abundance tables for gene and pathway levels, respectively.
    metagenome_pipeline.py -i ../table.biom -m marker_predicted_and_nsti.tsv.gz -f $EC_SHUFFLED \
                       -o $OUT_META \
                       --strat_out
    
     pathway_pipeline.py -i $OUT_META/pred_metagenome_unstrat.tsv.gz \
                         -o $OUT_PATHWAYS \
                         -p 1
done
```

## 3.compare the empirical data with shuffle data
These shuffled tables are especially helpful to get a baseline for how the predicted functional data differentiates samples (e.g. based on ordination or differential abundance testing) when the predicted ASV genomes are assigned randomly.