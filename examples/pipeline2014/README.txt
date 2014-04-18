This is a current version of the pipeline which we use for processing many Hi-C datasets. This version is for .sra files; make necessary adjustments for fastq files.  

Usage: 

1. place .sra files in the fastq folder 
2. run 00_makeLengths.sh
3. edit 01_iterative_mapping.py. Specify genome and paths to bowtie. 
4. Run 01_iterative_mapping
5. Provide datasets.tsv file. Sample file for (Lieberman 2009) is presented. 
6. run 02_mergeDatasets.py 

