'''
mapping - map raw Hi-C reads to a genome
========================================
'''
import os,sys
import glob
import gzip
import re
import subprocess

import numpy as np

import Bio, Bio.SeqIO, Bio.Seq, Bio.Restriction
import pysam

"""
The dictionary of chromosome # to the char ID correspondence.
"""
CHRM_ID = {i:str(i) for i in range(1,23)}
CHRM_ID.update({23:'X', 24:'Y', 25:'M'})

"""
The dictionary of the chromosome char ID to # correspondence.
"""
CHRM_IDX = {str(i):i for i in range(1,23)}
CHRM_IDX.update({'X':23, 'Y':24, 'M':25})


##TODO: write some autodetection of chromosome lengthes base on genome folder
##TODO: throw an exception if no chromosomes found in chromosome folder

"""
The dictionary of chromosome # to the char ID correspondence.
"""
CHRM_ID = {i:str(i) for i in range(1,20)}
CHRM_ID.update({20:'X', 21:'Y', 22:'M'})

"""
The dictionary of the chromosome char ID to # correspondence.
"""
CHRM_IDX = {str(i):i for i in range(1,20)}
CHRM_IDX.update({'X':20, 'Y':21, 'M':22})

_CACHED_RSITES_PATH = '%s_rsites.npz'

def _detect_quality_coding_scheme(in_fastq, num_entries = 10000):
    in_file = open(in_fastq)
    max_ord = 0
    min_ord = 256
    i = 0
    while True:
        line = in_file.readline()
        if not line or i > num_entries:
            break

        if not line.startswith('@'):
            raise Exception('%s does not comply with the FASTQ standards.' % line)

        fastq_entry = [line, in_file.readline(), 
                       in_file.readline(), in_file.readline()]
        min_ord = min(min_ord, min(ord(j) for j in fastq_entry[3].strip()))
        max_ord = max(max_ord, max(ord(j) for j in fastq_entry[3].strip()))

        i += 1

    return min_ord, max_ord

def filter_fastq(ids, in_fastq, out_fastq):
    '''Filter FASTQ sequences by their IDs.

    Read entries from **in_fastq** and store in **out_fastq** only those
    the whose ID are in **ids**.
    '''
    out_file = open(out_fastq, 'w')
    in_file = open(in_fastq)
    while True:
        line = in_file.readline()
        if not line:
            break

        if not line.startswith('@'):
            raise Exception('%s does not comply with the FASTQ standards.' % line)

        fastq_entry = [line, in_file.readline(), 
                       in_file.readline(), in_file.readline()]
        read_id = line.split("/")[0][1:]
        if read_id in ids:
            out_file.writelines(fastq_entry)

def split_uniquely_mapped(in_fastq, in_sam, nonunique_fastq, unique_sam):
    '''Split unique, non-unique and unmapped reads.

    1. Read raw sequences from **in_fastq** FASTQ file and alignments from 
       **in_sam** SAM file.
    2. Save the unique aligments in **unique_sam**.
    3. Save the non-mapped and non-uniquely aligned raw sequences in 
       **nonunique_fastq** FASTQ file.
    '''
    samfile = pysam.Samfile(in_sam)

    nonunique_ids = set()
    for read in samfile:
        tags_dict = dict(read.tags)
        read_id = read.qname
        # If exists, the option 'XS' contains the score of the second 
        # best alignment. Therefore, its presence means non-unique alignment.
        if 'XS' in tags_dict or read.is_unmapped:
            nonunique_ids.add(read_id)
        else:
            unique_sam.write(read)

    filter_fastq(nonunique_ids, in_fastq, nonunique_fastq)

def iterative_mapping(bowtie_path, genome_path, fastq_path, out_sam_path,
                      min_seq_len, len_step, seq_start=0, seq_end=None,
                      nthreads = 4, bowtie_flags = ''):
    '''Map raw HiC reads iteratively with bowtie2.
    http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

    1. Truncate the sequences to the first N = **min_seq_len** base pairs,
       starting at the **seq_start** position.
    2. Map the sequences using bowtie2.
    3. Store the uniquely mapped sequences in a SAM file at **out_sam_path**.
    4. Go to the step 1, increase the truncation length N by **len_step** base 
       pairs, and map the non-mapped and non-uniquely mapped sequences,
       ...
       Stop when the 3' end of the truncated sequence reaches the **seq_end**
       position.
       
    Parameters
    ----------

    bowtie_path : str
        The path to the bowtie2 executable.

    genome_path : str
        The path to the bowtie2 genome index. Since the index consists of 
        several files with the different suffices (e.g., hg18.1.bt2, 
        hg18.2.bt.2), provide only the common part (hg18).
        
    fastq_path : str
        The path to the input FASTQ file.

    out_sam_path : str
        The path to the output SAM file.

    min_seq_len : int
        The truncation length at the first iteration of mapping.

    len_step : int
        The increase in truncation length at each iteration.

    seq_start : int, optional
        The start position of the mapped sequence in a FASTQ sequence.
        Default is 0.

    seq_end : int, optional
        The last position of the analyzed sequence in a FASTQ sequence.
        Default is None, i.e. the last position of the FASTQ sequence.

    nthreads : int, optional
        The number of Bowtie2 threads. Default is 8

    bowtie_flags : str, optional
        Extra command-line flags for Bowtie2
    '''

    raw_seq_len = len(Bio.SeqIO.parse(open(fastq_path), 'fastq').next().seq)
    if (seq_start < 0 
        or seq_start > raw_seq_len 
        or (seq_end and seq_end > raw_seq_len)):
        raise Exception('An incorrect trimming region is supplied: [%d, %d), '
                        'the raw sequence length is %d' % (
                            seq_start, seq_end, raw_seq_len))
    max_seq_len = seq_end - seq_start

    out_sam_file = None
    prev_fastq_path = fastq_path
    local_fastq_path = fastq_path
    for seq_len in range(min_seq_len, max_seq_len + 1, len_step):
        local_seq_start = seq_start
        local_seq_end = seq_end - (max_seq_len - seq_len)
        local_sam_path = '%s_%s_%s_%s.sam' % (
            os.path.splitext(fastq_path)[0],
            os.path.split(genome_path)[-1],
            str(local_seq_start),
            str(local_seq_end))

        local_5_trim = seq_start
        local_3_trim = raw_seq_len - local_seq_end
        bowtie_command = (
            ('time %s -x %s --fast --score-min L,-0.6,-0.2 '
             '-q %s -5 %s -3 %s -p %s %s > %s') % (
                bowtie_path, 
                genome_path, 
                local_fastq_path, 
                str(local_5_trim),
                str(local_3_trim),
                str(nthreads),
                bowtie_flags,
                local_sam_path))

        print 'Map reads:', bowtie_command
        subprocess.call(bowtie_command, shell=True)

        print ('Save unique aligments and send the '
               'non-unique ones to the next iteration')
        if out_sam_file is None:
            if os.path.splitext(out_sam_path)[1] in ['.bam', '.BAM']:
                write_mode = 'wb'
            else:
                write_mode = 'wh'
            out_sam_file = pysam.Samfile(
                out_sam_path, write_mode, 
                template=pysam.Samfile(local_sam_path))

        local_fastq_path = '%s_%s_%s_%s.fastq' % (
            os.path.splitext(fastq_path)[0],
            os.path.split(genome_path)[-1],
            str(local_seq_start),
            str(local_seq_end))
        split_uniquely_mapped(prev_fastq_path, local_sam_path,
                              local_fastq_path, out_sam_file)

        prev_fastq_path = local_fastq_path
     
def _find_rfrags(cuts, chrms, dirs, rsites_map, min_frag_size):
    '''Private: assign mapped reads to restriction fragments by 
    their 5' end position.
    '''
    rfrags = np.zeros(len(cuts), dtype=np.int64)
    rsites = np.zeros(len(cuts), dtype=np.int64)
    uprsites = np.zeros(len(cuts), dtype=np.int64)
    downrsites = np.zeros(len(cuts), dtype=np.int64)

    # If the fragment was not mapped.
    rfrags[chrms == 0] = -1
    rsites[chrms == 0] = -1
    uprsites[chrms == 0] = -1
    downrsites[chrms == 0] = -1

    for chrm_idx in set.difference(set(chrms), set([0])):
        chrm_id = CHRM_ID[chrm_idx]
        all_rsites = rsites_map[CHRM_ID[chrm_idx]]
        idxs = (chrms == chrm_idx)

        # Find the indexes of the restriction fragment...
        rfrags[idxs] = np.searchsorted(all_rsites, cuts[idxs]) - 1
        uprsites[idxs] = all_rsites[rfrags[idxs]]
        downrsites[idxs] = all_rsites[rfrags[idxs] + 1] 
        rsites[idxs] = np.where(dirs[idxs], downrsites[idxs], uprsites[idxs])

        too_close = (np.abs(rsites[idxs] - cuts[idxs]) <= min_frag_size)
        too_close_idxs = np.where(idxs)[0][too_close]
        rfrags[too_close_idxs] += dirs[too_close_idxs] * 2 - 1
        uprsites[too_close_idxs] = all_rsites[rfrags[too_close_idxs]]
        downrsites[too_close_idxs] = all_rsites[rfrags[too_close_idxs] + 1]
        rsites[too_close_idxs] = np.where(
            dirs[too_close_idxs],
            downrsites[too_close_idxs], 
            uprsites[too_close_idxs])
        
    return rfrags, rsites, uprsites, downrsites
    
def _cache_rsites(db_dir_path, enzyme_name, chr_file_template = 'chr%s.fa'):
    '''Private: cache the restriction maps of a genome.
    '''
    enzyme_search_func = eval('Bio.Restriction.%s.search' % enzyme_name)
    rsites_map = {}
    for fasta_path in glob.glob(
            os.path.join(db_dir_path, chr_file_template % ('*',))):
        genome_db = Bio.SeqIO.read(open(fasta_path), 'fasta')
        chrm = re.search(chr_file_template % ('(.*)',), fasta_path).group(1)
        rsites_map[chrm] = np.array(enzyme_search_func(genome_db.seq)) + 1

        # Insert the ends of a chromosome.
        rsites_map[chrm] = np.insert(rsites_map[chrm], 0, 0)
        rsites_map[chrm] = np.append(rsites_map[chrm], len(genome_db.seq))

    np.savez(os.path.join(db_dir_path, _CACHED_RSITES_PATH % (enzyme_name)),
             **rsites_map)

def fill_rsites(lib, db_dir_path, enzyme_name, 
                min_frag_size = None, chr_file_template = 'chr%s.fa'):
    '''Private: assign mapped reads to restriction fragments by 
    their 5' end position.

    Parameters
    ----------

    lib : dict
        A library of mapped Hi-C molecules. Modified by the function.

    db_dir_path: path to the genome location
        
    '''
    if enzyme_name not in Bio.Restriction.AllEnzymes:
        raise Exception('Enzyme is not found in the library: %s' % (enzyme_name,))

    # Ugly workaround for enzyme name lookup.
    enzyme_search_func = eval('Bio.Restriction.%s.search' % enzyme_name)

    # Build the restriction maps of the chromosomes.
    rsites_map = {}
    target_chrms = set(lib['chrms1']).union(set(lib['chrms2']))
    target_chrms.discard(0)

    cached_rsites_path = os.path.join(
        db_dir_path, _CACHED_RSITES_PATH % (enzyme_name))

    if not os.path.isfile(cached_rsites_path):
        _cache_rsites(db_dir_path, enzyme_name, chr_file_template)

    rsites_map = dict(np.load(cached_rsites_path))

    rsite_size = eval('len(Bio.Restriction.%s.site)' % enzyme_name)
    if min_frag_size is None:
        _min_frag_size = rsite_size / 2.0
    else:
        _min_frag_size = min_frag_size
        
    lib['rfrags1'], lib['rsites1'], lib['uprsites1'], lib['downrsites1'] = \
        _find_rfrags(lib['cuts1'], lib['chrms1'], lib['dirs1'], 
                     rsites_map, _min_frag_size)

    lib['rfrags2'], lib['rsites2'], lib['uprsites2'], lib['downrsites2'] = \
        _find_rfrags(lib['cuts2'], lib['chrms2'], lib['dirs2'], 
                     rsites_map, _min_frag_size)

    return lib

def _parse_single_sam(sam_path, target_chrms, reverse_complement=False):
    sam_file = pysam.Samfile(sam_path)
    chrms, cuts, dirs, seqs, ids = [], [], [], [], []
    for read in sam_file:
        if (read.tid + 1) in target_chrms:
            chrms.append(read.tid + 1)
            dirs.append(not read.is_reverse)
            ids.append(read.qname)
            if read.is_reverse:
                if reverse_complement:
                    seqs.append(Bio.Seq.reverse_complement(read.seq))
                else:
                    seqs.append(read.seq)
                cuts.append(read.pos + len(read.seq))
            else:
                seqs.append(read.seq)
                cuts.append(read.pos)

    return (np.array(chrms, dtype=np.int8), 
            np.array(cuts, dtype=np.int32),
            np.array(dirs, dtype=np.bool),
            np.array(seqs),
            np.array(ids))

def parse_sam(sam_path1, sam_path2, target_chrms = range(0, 25),
              reverse_complement=False):
    # TODO: Rewrite the parser using lazy reading.
    lib = {}
    chrms1, cuts1, dirs1, seqs1, ids1 = _parse_single_sam(
        sam_path1, target_chrms, reverse_complement)
    chrms2, cuts2, dirs2, seqs2, ids2 = _parse_single_sam(
        sam_path2, target_chrms, reverse_complement)
    
    all_ids = set.union(set(ids1), set(ids2))
    all_ids = list(all_ids)
    all_ids.sort()
    tot_num_reads = len(all_ids)
    ids_to_idx = {all_ids[i]:i for i in xrange(tot_num_reads)}

    sorting_1 = -1 * np.ones(tot_num_reads, dtype=np.int)
    for i in xrange(len(ids1)):
        sorting_1[ids_to_idx[ids1[i]]] = i

    sorting_2 = -1 * np.ones(tot_num_reads, dtype=np.int)
    for i in xrange(len(ids2)):
        sorting_2[ids_to_idx[ids2[i]]] = i

    num_reads1 = len(chrms1)
    lib['chrms1'] = np.append(chrms1,    [0] * (tot_num_reads-num_reads1))[sorting_1]
    lib['cuts1']  = np.append(cuts1,    [-1] * (tot_num_reads-num_reads1))[sorting_1]
    lib['dirs1']  = np.append(dirs1, [False] * (tot_num_reads-num_reads1))[sorting_1]
    lib['seqs1']  = np.append(seqs1,    [''] * (tot_num_reads-num_reads1))[sorting_1]
    lib['ids1']   = np.append(ids1,     [''] * (tot_num_reads-num_reads1))[sorting_1]

    num_reads2 = len(chrms2)
    lib['chrms2'] = np.append(chrms2,    [0] * (tot_num_reads-num_reads2))[sorting_2]
    lib['cuts2']  = np.append(cuts2,    [-1] * (tot_num_reads-num_reads2))[sorting_2]
    lib['dirs2']  = np.append(dirs2, [False] * (tot_num_reads-num_reads2))[sorting_2]
    lib['seqs2']  = np.append(seqs2,    [''] * (tot_num_reads-num_reads2))[sorting_2]
    lib['ids2']   = np.append(ids2,     [''] * (tot_num_reads-num_reads2))[sorting_2]

    return lib

