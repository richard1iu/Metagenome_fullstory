https://bioinf.shenwei.me/taxonkit/tutorial/
https://www.jianshu.com/p/1d6edfcb4110
https://blog.csdn.net/qq_42491125/article/details/92790567

# 1.Install
```bash
# 1.download
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
tar -zxvf taxdump.tar.gz
# 2.transfer all dmp to ./taxonkit
mkdir -p $HOME/.taxonkit
cp *.dmp $HOME/.taxonkit 
# 3.check it
ll -h  ~/.taxonkit/
```

# 2*. query lineage of given taxids

```bash
# 1.query single taxids
taxonkit list --ids 9606,10090 --show-name  --show-rank
#33090: plant; 2: Bacteria; 2157: Archaea; 10239:Viruses; Eukaryota: #2759；Fungi: 4751; 9606: human being
# 2. list all taxids of a certain domain
taxonkit list --ids 2 --indent "" > bacteria.taxid.txt
# 3. annotation of all taxids
less bacteria.taxid.txt|taxonkit lineage | taxonkit reformat -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" -F | cut -f 1,3- | sed '1i\Taxid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenu\tSpecies' > bacteria_taxid_ano.txt
```

# 3.Construct bacteria sublibrary

```bash
# Find the taxonid’s corresponding accession versions of proteins 
zcat prot.accession2taxid.gz | csvtk -t -f taxid grep -P bac.taxid.txt | csvtk -t cut -f accession.version >bac.taxid.acc.txt

# NOTICE!: old-version blast has no seqidlist parameter
blastdb_aliastool -seqidlist bacteria.taxid.acc.txt -db /home/software/nr-2019-12-18/nr -out nr_bac -title nr_bac

blastp -query NR100pro.fasta -db /home/pub_guest/db/nr/nr_bac -out D1_nr.out -outfmt "6 qseqid qgi qacc qaccver qlen sseqid qseq sseq evalue  score length pident  staxids sscinames salltitles " -num_threads 16 -evalue 1e-5 -num_alignments 5
```