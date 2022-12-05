```bash
# install
conda install -c conda-forge -c bioconda snakemake -y
```

## Raw data single end
NCBI Accession: SRP009826
species: Glycine max
## make config.yaml 

## make config.yaml
### make sample list
```bash
#
ls *.fq> fq_list.txt

#
echo "DATA:" > config.yaml
cat fq_list.txt|while read line;do id=$(basename $line .fq);echo $id >> config.yaml;done

#
sed '2,$ s/^/- /' config.yaml -i
```

### 