echo "Download Bowtie2 binaries"
mkdir -p bin/bowtie2
cd bin/bowtie2

wget http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.0.0-beta7/bowtie2-2.0.0-beta7-linux-x86_64.zip -O bowtie2.zip
unzip bowtie2.zip
mv bowtie2-2.0.0-beta7/* ./
rm -r bowtie2-2.0.0-beta7

