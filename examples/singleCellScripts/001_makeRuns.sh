rm runs.tsv
for i in `ls mapped-mm9 | grep -v combined | grep -v lock | sort`; do echo mapped-mm9/$i $(echo $i | head -c 12) R1 mm9 DpnII >> runs.tsv; done
for i in `ls mapped-hg19 | grep -v combined | grep -v lock | sort`; do echo mapped-hg19/$i $(echo $i | head -c 12) R1 hg19 DpnII >> runs.tsv; done
