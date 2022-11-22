给定一段fasta格式序列，如何预测其中的基因呢？

首先需要判断该片段来自:
- 原核生物:prodigal/glimmer3
- 真核生物:augustus/snap
- 病毒序列


## 1.原核生物

能预测**蛋白编码基因**，不能预测**RNA编码基因** (barrnap:rRNA, tRNAscan-SE:tRNA)
不能预测**包含内含子**的基因
不能为预测的基因提供**功能注释**
不能处理**移码突变（indel）**的情况


```bash
prodigal -i ref.fna -a ref.pep -d ref.cds -f gff -g 11 -o ref.gff -s ref.stat  >prodigal.log

```

-i：输入文件，fasta格式
-a：基因的氨基酸序列
-d：基因的核酸序列
-f：输出文件类型gbk, gff, or sco
-g：密码子表，细菌为第11
-p：模式，单菌还是宏基因组
-o：输出结果文件，有多种格式可选
-s：统计信息

## 2.真核生物

```bash
augustus --strand=both --genemodel=partial --singlestrand=false --protein=on --introns=on --start=on --stop=on --cds=on --codingseq=on --alternatives-from-evidence=true --gff3=on --UTR=on --outfile=out.gff --species=human HS04636.fa
```

## 3.基因功能注释
给定一个fasta格式的氨基酸序列，如何得到基因的功能信息？可以使用eggnog-mapper进行分析。

```
emapper.py -i gene.fasta --output polb_bact -d bact --data_dir eggnog-mapper-1.0.3/data/
```

-i：输入文件，基因的氨基酸序列
-m：选择运行模式hmmer或者diamond
-h：输出帮助文档
–output：输出结果前缀
–output_dir：输出结果目录
–data_dir：数据库目录
–database：单独指定数据库
–dmnd_db：单独指定diamond数据库路径

## 4.rRNA预测
给定一段序列，如何找到rRNA，包括原核生物的5S，16S，23S，真核生物的5.8S，18S，28SRNA等.
由于核糖体RNA具高的保守性，因此预测准确性较高。
使用rnammer或者barrnap 软件or使用Infernal基于数据库rfam软件，直接输入fasta序列即可。

```bash
rnammer -S bac -m tsu,lsu,ssu -gff ref.gff -f ref.frn ref.fna
```
-S：物种类型，古细菌，细菌或者真菌
-m：需要rRNA类型，如果真要16S，则单独选择lsu
-gff：输出gff格式结果
-f：输出fasta格式序列

## 5.tRNA预测
给定一段序列，如何找到tRNA，可以使用tRNAscan工具

```bash
tRNAscan-SE -B -o tRNAScan.out -f tRNAScan.out.structure -m stat.list ref.fna
```
-B ：物种为细菌
-A ：物种为古细菌
-O ：输入序列为细胞器
-G ：包括全部类型
-o：输出结果
-f：tRNA二级结构
-m：统计结果