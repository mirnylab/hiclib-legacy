echo "Download SRA binaries"
mkdir -p bin/sra
cd bin/sra

wget http://ftp-private.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

tar -xvzf sratoolkit.current-ubuntu64.tar.gz

mv sratoolkit.*-ubuntu64/* ./
rm -r sratoolkit.current-ubuntu64.tar.gz

