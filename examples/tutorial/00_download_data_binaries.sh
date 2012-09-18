echo "Download Bowtie2 binaries"
mkdir bin
cd bin

wget http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.0.0-beta7/bowtie2-2.0.0-beta7-linux-x86_64.zip -O bowtie2.zip
unzip bowtie2.zip

echo "Download hg19 index for bowtie2..."
cd bowtie2-2.0.0-beta7
mkdir index
cd index
wget ftp://ftp.cbcb.umd.edu/pub/data/bowtie2_indexes/incl/hg19.1.zip
wget ftp://ftp.cbcb.umd.edu/pub/data/bowtie2_indexes/incl/hg19.2.zip
wget ftp://ftp.cbcb.umd.edu/pub/data/bowtie2_indexes/incl/hg19.3.zip
unzip hg19.1.zip
unzip hg19.2.zip
unzip hg19.3.zip
cd ../../..

echo "Download SRA binaries"
cd bin 
wget http://ftp-private.ncbi.nlm.nih.gov/sra/sdk/2.1.16/sratoolkit.2.1.16-ubuntu32.tar.gz -O sra.tar.gz
tar -xvzf sra.tar.gz
cd ..

echo "Download the data..."
mkdir data
cd data
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP001%2FSRP001261/SRR027956/SRR027956.sra -O SRR027956.sra
mkdir tmp
cd ../

echo "Download the hg19 genome FASTA database"
mkdir data/hg19
cd data/hg19
for i in $(seq 1 1 22); do wget "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr$i.fa.gz"; done
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrX.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrY.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrM.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz
gunzip *.gz
cd ../../
