import os

import numpy as np
import Bio.Seq

import mirnylib.genome

NUM_READS = 10000
READ_LENGTH = 75
MOLECULE_LENGTH = 400
NOISE = 0.02
RANDOM_JOINT = True
STD_HEADER = 'SRR000000.{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'

genome_db = mirnylib.genome.Genome('../../fasta/sacCer3', readChrms=['#'])

if not os.path.exists('tmp'):
    os.mkdir('tmp')

def noise_seq(seq, noise):
    for i in range(len(seq)):
        if np.random.random() < noise:
            seq = seq[:i] + 'ATGC'[np.random.randint(4)] + seq[i+1:]
    return seq

fastq1 = open('./tmp/insilico_1.fastq', 'w')
fastq2 = open('./tmp/insilico_2.fastq', 'w')
for i in range(NUM_READS):
    chr1 = np.random.randint(16)
    pos1 = np.random.randint(MOLECULE_LENGTH,
                             genome_db.chrmLens[chr1] - MOLECULE_LENGTH)
    if RANDOM_JOINT:
        len1 = np.random.randint(MOLECULE_LENGTH)
    else:
        len1 = MOLECULE_LENGTH / 2
    seq1 = str(genome_db.seqs[chr1][pos1:pos1+len1].seq)
    seq1 = noise_seq(seq1, NOISE)

    chr2 = np.random.randint(16)
    pos2 = np.random.randint(MOLECULE_LENGTH, 
                             genome_db.chrmLens[chr2] - MOLECULE_LENGTH)
    len2 = MOLECULE_LENGTH - len1
    seq2 = str(genome_db.seqs[chr2][pos2:pos2+len2].seq)
    seq2 = noise_seq(seq2, NOISE)

    concatenated_seq = seq1+seq2
    
    read_id = STD_HEADER.format(i, chr1, pos1, len1, 
                                chr2, pos2+len2, len2, MOLECULE_LENGTH)

    fastq1.writelines(
            ['@' + read_id  + '\n',
             concatenated_seq[:READ_LENGTH] + '\n',
             '+' + read_id + '\n',
             'g' * READ_LENGTH + '\n'])
    fastq2.writelines(
            ['@' + read_id + '\n',
             Bio.Seq.reverse_complement(concatenated_seq[-READ_LENGTH:]) + '\n',
             '+' + read_id + '\n',
             'g' * READ_LENGTH + '\n'])

    if (i+1) % 10000 == 0:
        print '{0} random sequences generated...'.format((i+1))

