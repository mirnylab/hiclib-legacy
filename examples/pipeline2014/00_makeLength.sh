#!/bin/bash
mkdir lengths
for i in `ls fastq`; do ./fastq-dump -Z fastq/$i | head -n 2 | tail -n 1 | wc -c > lengths/$i; done 
