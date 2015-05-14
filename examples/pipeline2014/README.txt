This is a current version of the pipeline which we use for processing many Hi-C datasets. This version is for .sra files; make necessary adjustments for fastq files.  

Usage: 

1. place .sra files in the fastq folder 
2. run 00_makeLengths.sh  -  it will make a folder "lengths" which has lengths of reads
3. edit 01_iterative_mapping.py. Specify genome and paths to bowtie. 
4. Run 01_iterative_mapping   -  it will map reads and output .hdf5 files with mapped read data (raw, unfiltered)
5. Provide datasets.tsv file. Sample file for (Lieberman 2009) is presented. 
6. run 02_mergeDatasets.py   -  it will do several things.
    -- merge (if necessary) different .hdf5 files corresponding to one replica of one datasets  
    -- do fragment-level filtering of merged files    
    -- save heatmaps at different resolutions
    -- merge files from different replicas of the same experiment
    -- save heatmaps 

