#!/bin/bash

#This script makes folder which keeps lengths of all SRA files 
mkdir lengths
for i in `ls fastq`; do ./fastq-dump -Z fastq/$i | head -n 2 | tail -n 1 | wc -c > lengths/$i; done 
