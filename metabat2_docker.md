

```bash
# pull metabat2 image
docker image pull metabat/metabat:v2.15-12-gc318b89

# run metabat2 container
docker run -it -d metabat/metabat:v2.15-12-gc318b89 bash

# get the [CONTAINER_ID]
docker ps

# copy input files into container
docker cp ./tara_sorted.bam [CONTAINER_ID]:/home
docker cp ./assembled_contigs.fa [CONTAINER_ID]:/home
# exec
docker exec -it [CONTAINER_ID] bash 

```

```bash
# generate both bins and depth.txt
runMetaBat.sh assembled_contigs.fa tara_sorted.bam

# generate depth file
jgi_summarize_bam_contig_depths --outputDepth assembled_contigs.fa.depth.txt --pairedContigs paired.txt tara_sorted.bam

# generate bins
metabat2 -i assembled_contigs.fa -a assembled_contigs.fa.depth.txt -o ./bin -v 

```







