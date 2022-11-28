# Reference
[Harvard Chan Bioinformatics Core (HBC)](https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/experimental_planning_considerations.html)

Three special considerations can greatly affect the experiment results
1. Number and type of replicates
2. Avoiding confounding
3. Addressing batch effects


## 1.Replicate

### 1.1 biological and technical replicates
Experimental replicates can be performed as 
1. Technical replicates: use **the same biological sample** to repeat the technical or experimental steps to accurately measure **technical variation** and remove it during analysis.
2. Biological replicates use **different biological samples** of the same condition to measure the **biological variation** between samples.

![replicates tyeps](https://hbctraining.github.io/DGE_workshop_salmon_online/img/replicates.png)

## 1.2 replicates and sequencing depth
an increase in the number of replicates tends to return more differential expression genes than increasing the sequencing depth.

some general guidelines for replicates and sequencing depth:
1. General gene-level differential expression:
    - ENCODE guidelines suggest 30 million SE reads per sample (stranded).
    - 15 million reads per sample is often sufficient, if there are a good number of replicates (>3).
    - Spend money on more biological replicates, if possible.
    - Generally recommended to have read length >= 50 bp
2. Gene-level differential expression with detection of lowly-expressed genes:
    - Similarly benefits from replicates more than sequencing depth.
    - Sequence deeper with at least 30-60 million reads depending on level of expression (start with 30 million with a good number of replicates).
    - Generally recommended to have read length >= 50 bp
3. Isoform-level differential expression:
    - Of known isoforms, suggested to have a depth of at least 30 million reads per sample and paired-end reads.
    - Of novel isoforms should have more depth (> 60 million reads per sample).
    - Choose biological replicates over paired/deeper sequencing.
    - Generally recommended to have read length >= 50 bp, but longer is better as the reads will be more likely to cross exon junctions
4. Other types of RNA analyses (intron retention, small RNA-Seq, etc.):
    - Different recommendations depending on the analysis.
    - Almost always more biological replicates are better!

## 2.Confounding
A confounded experiment is one where you cannot distinguish the separate effects of two different sources of variation in the data.

For example, if all of our control mice were female and all of the treatment mice were male, then our treatment effect would be confounded by sex. We could not differentiate the effect of treatment from the effect of sex.

## 3.Batch effects
Sometimes, there are more significant differences between batchs than between treats.
![batch effects](https://hbctraining.github.io/DGE_workshop_salmon_online/img/batch_effect_pca.png)

### How to know whether you have batches?
1. Were all RNA isolations performed **on the same day**?
2. Were all library preparations performed on **the same day**?
3. Did the **same person** perform the RNA isolation/library preparation for all samples?
4. Did you use **the same reagents** for all samples?
5. Did you perform the RNA isolation/library preparation **in the same location**?

If any of the answers is **‘No’**, then you have batches.