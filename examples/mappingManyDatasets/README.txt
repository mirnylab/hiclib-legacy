This is an actual script that was used to re-map data from (Dixon, Nature 2012). 
Read lengths were determined by the following bash one-liner: 
for i in *.fastq; do echo -n $i; echo -n " "; head -n 2 $i | tail -n 1 | wc -m ; done > datasets.txt

