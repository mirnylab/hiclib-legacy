echo "Download SRA binaries"
mkdir -p bin/sra
cd bin/sra

wget http://ftp-private.ncbi.nlm.nih.gov/sra/sdk/2.1.16/sratoolkit.2.1.16-ubuntu32.tar.gz -O sra.tar.gz
tar -xvzf sra.tar.gz

mv sratoolkit.2.1.16-ubuntu32/* ./
rm -r sratoolkit.2.1.16-ubuntu32

