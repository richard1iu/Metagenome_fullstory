# GNU parallel
<<<<<<< HEAD

## install
>Notice: DON'T install prallel using apt 
```bash
# download and install
(wget -O - pi.dk/3 || curl pi.dk/3/) | bash
# export to env
echo "export PATH=$PATH:/home/<yourusername>/bin"
```

## |, :::, -a
```bash
seq 5 > ls.txt
# get the same result
cat ls.txt | parallel echo id_{}
parallel echo id_{} ::: $(cat ls.txt)
parallel -a ls.txt echo id_{}
```

## parallel echo
```bash
# compare the difference
echo sample.{a..z}.fq.gz
parallel echo ::: sample.{a..z}.fq.gz
for i in {a..z};do echo sample.$i.fq.gz; done

# link all paired
parallel echo ::: a b c ::: 1 2 3
#
```

##
```bash
# specific separator with -d, defualt separtor is \n
parallel echo :::: /etc/passwd
parallel -d : echo :::: /etc/passwd 

```

## some option
- --eta: Shows the remaining time to finish jobs.
- -j | --jobs: the number of jobs (integer) simultaneous running, less than the core number of your PC
- --load: The maximum CPU load (percent. e.g. 80%), if CPU load reach it, then new job will not be started
- --noswap: New jobs won't be started if the most memory of your PC are occupied  
=======
>>>>>>> 553a374 (parallel)
