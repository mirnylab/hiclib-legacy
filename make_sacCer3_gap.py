import subprocess
import gzip
import csv
import os
import sys
import mirnylib.genome 

if len(sys.argv) != 2:
    print 'Please, supply the path to the folder with the sacCer3 genome.'
    sys.exit(1)
print sys.argv

print 'Download and parse sgdOther track...'
subprocess.call(
    'wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/database/sgdOther.txt.gz',
    shell=True)

centromere_starts = {}
centromere_ends = {}
for line in csv.reader(gzip.open('sgdOther.txt.gz'), dialect='excel-tab'):
    if line[4].startswith('CEN'):
        chr_num = int(line[4][3:])
        centromere_starts[chr_num] = min(int(line[2]), 
                                         centromere_starts.get(chr_num, 1e9))
        centromere_ends[chr_num] = max(int(line[3]), 
                                       centromere_ends.get(chr_num, -1))
os.remove('sgdOther.txt.gz')

print 'Save the centromere positions into a .gap file'
centromere_positions = {}
for i in centromere_starts:
    centromere_positions['chr' + str(i)] = (
        centromere_starts[i], centromere_ends[i])
genome_db = mirnylib.genome.Genome(sys.argv[1])
genome_db.createGapFile(centromere_positions)

