# (c) 2012 Massachusetts Institute of Technology. All Rights Reserved
# Code written by: Anton Goloborodko (golobor@mit.edu)

'''
This module contains the functions that map raw DNA sequences obtained in
the Hi-C experiment to a supplied genome.

The three main methods of this module are iterative_mapping, parse_sam and
fill_rsite.

The first, iterative_mapping() applies the bowtie2 read alignment software to
the raw reads from the sequencer. The second method, parse_sam() parses
the bowtie output, combines individual reads into pairs and converts the data
into the internal format that may be fed to the downstream functions. Finally,
fill_rsite() maps the sequences onto the restriction fragments.

-------------------------------------------------------------------------------

API Documentation
-----------------
'''

import os
import sys
import re
import glob
import gzip
import subprocess
import tempfile
import atexit
import logging
import warnings
import numpy as np
import Bio.SeqIO
import Bio.Seq
import Bio.Restriction
import pysam
import time
import gc

import mirnylib.h5dict
import mirnylib.genome

# #TODO: write some autodetection of chromosome lengthes base on genome folder
# #TODO: throw an exception if no chromosomes found in chromosome folder
# #TODO: fix #-to-ID correspondence for other species.

log = logging.getLogger(__name__)


def sleep():
    """sleep for a second, run garbage collector, sleep again.
    Sleep is split in small pieces to allow some callbacks to
    possibly terminate in between (I don't know if it makes sense, but
    it definitely does not hurt)"""
    for _ in range(10):
        time.sleep(0.1)
    gc.collect()
    for _ in range(10):
        time.sleep(0.1)



def commandExists(command):
    "checks if the bash command exists"
    command = command.split()[0]
    if subprocess.call(['which', command]) != 0:
        return False
    return True


def gzipWriter(filename):
    """
    creates a writing process with gzip or parallel gzip (pigz) attached to it
    """
    filename = os.path.abspath(filename)
    if commandExists("pigz"):
        writer = "pigz -c -4"
    else:
        writer = "gzip -c -1"
        warnings.warn("Please install 'pigz' parallel gzip for faster speed")
    bashLine = "{writer} > {filename}".format(writer=writer, filename=filename)
    log.info("""Writer created with command "{0}" """.format(bashLine))
    pwrite = subprocess.Popen([bashLine], stdin=subprocess.PIPE,
                               shell=True, bufsize=-1)
    return pwrite


def _detect_quality_coding_scheme(in_fastq, num_entries=10000):
    in_file = open(in_fastq)
    max_ord = 0
    min_ord = 256
    i = 0
    while True:
        line = in_file.readline()
        if not line or i > num_entries:
            break

        if not line.startswith('@'):
            raise Exception('%s does not comply with the FASTQ standards.')

        fastq_entry = [line, in_file.readline(),
                       in_file.readline(), in_file.readline()]
        min_ord = min(min_ord, min(ord(j) for j in fastq_entry[3].strip()))
        max_ord = max(max_ord, max(ord(j) for j in fastq_entry[3].strip()))

        i += 1

    return min_ord, max_ord


def _line_count(path):
    '''Count the number of lines in a file. The function was posted by
    Mikola Kharechko on Stackoverflow.
    '''

    f = open(path)
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read  # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    return lines


def _chunk_file(in_path, out_basename, max_num_lines):
    '''Slice lines from a large file.
    The line numbering is as in Python slicing notation.
    '''
    num_lines = _line_count(in_path)
    if num_lines <= max_num_lines:
        return [in_path, ]

    out_paths = []

    for i, line in enumerate(open(in_path)):
        if i % max_num_lines == 0:
            out_path = out_basename + '.%d' % (i // max_num_lines + 1)
            out_paths.append(out_path)
            out_file = file(out_path, 'w')
        out_file.write(line)

    return out_paths


def _filter_fastq(ids, inStream, out_fastq, in_filename="none"):
    '''Filter FASTQ sequences by their IDs.

    Read entries from **in_fastq** and store in **out_fastq** only those
    the whose ID are in **ids**.
    '''
    writingProcess = gzipWriter(out_fastq)

    num_filtered = 0
    num_total = 0
    while True:
        line = inStream.readline()

        if not line:
            break

        if not line.startswith('@'):
            raise Exception(
                '{0} does not comply with the FASTQ standards.'.format(in_filename))

        fastq_entry = [line, inStream.readline(),
                       inStream.readline(), inStream.readline()]
        read_id = line.split()[0][1:]
        if read_id in ids:
            writingProcess.stdin.writelines(fastq_entry)
            num_filtered += 1
        num_total += 1
    writingProcess.stdin.flush()
    writingProcess.stdin.close()
    writingProcess.wait()
    sleep()
    if writingProcess.returncode != 0:
        raise RuntimeError("Writing process return code {0}".format(writingProcess.returncode))
    return num_total, num_filtered


def _filter_unmapped_fastq(in_stream, in_sam, nonunique_fastq, in_filename="none"):
    '''Read raw sequences from **in_fastq** and alignments from
    **in_sam** and save the non-uniquely aligned and unmapped sequences
    to **unique_sam**.
    '''
    samfile = pysam.Samfile(in_sam)

    nonunique_ids = set()
    for read in samfile:
        tags_dict = dict(read.tags)
        read_id = read.qname

        # If exists, the option 'XS' contains the score of the second
        # best alignment. Therefore, its presence means a non-unique alignment.
        if 'XS' in tags_dict or read.is_unmapped:
            nonunique_ids.add(read_id)

    num_total, num_filtered = _filter_fastq(
        nonunique_ids, in_stream, nonunique_fastq, in_filename=in_filename)
    sleep()

    return num_total, num_filtered


def iterative_mapping(bowtie_path, bowtie_index_path, fastq_path, out_sam_path,
                      min_seq_len, len_step, **kwargs):
    '''Map raw HiC reads iteratively with bowtie2.
    http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    Iterative mapping accounts for the modification of fragments' sequences
    due to ligation.

    The algorithm of iterative correction:

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

    bowtie_index_path : str
        The path to the bowtie2 genome index. Since the index consists of
        several files with the different suffices (e.g., hg18.1.bt2,
        hg18.2.bt.2), provide only the common part (hg18).

    fastq_path : str
        The path to the input FASTQ or gzipped FASTQ file.

    out_sam_path : str
        The path to the output SAM file. If ends with .bam, then the output
        is converted to BAM.

    min_seq_len : int
        The truncation length at the first iteration of mapping.

    len_step : int
        The increase in truncation length at each iteration.

    seq_start, seq_end : int, optional
        Slice the FASTQ sequences at [seq_start:seq_end]. Default is [O:None].

    nthreads : int, optional
        The number of Bowtie2 threads. Default is 8

    bowtie_flags : str, optional
        Extra command-line flags for Bowtie2. Default is ''.

    max_reads_per_chunk : int, optional
        If positive then split input into several chunks with
        `max_reads_per_chunk` each and map them separately. Use for large
        datasets and low-memory machines.

    temp_dir : str, optional
        The path to the temporary folder. If not specified, this path is
        supplied by the OS.

    bash_reader : str, optional
        A bash application to convert the input to the FASTQ format. The
        application have to save the output into stdout.
        The default value is None, that is the app is autodetected by the
        extension (i.e. cat for .fastq, gunzip for .gz).

    '''
    bowtie_path = os.path.abspath(os.path.expanduser(bowtie_path))
    if not os.path.isfile(bowtie_path):
        raise Exception(
            'The bowtie binary is not found '
            'at the specified path: {0}.'.format(bowtie_path))
    bowtie_index_path = os.path.abspath(os.path.expanduser(bowtie_index_path))

    fastq_path = os.path.abspath(os.path.expanduser(fastq_path))
    if not os.path.isfile(fastq_path):
        raise Exception(
            'The fastq file is not found '
            'at the specified path: {0}.'.format(fastq_path))

    already_mapped = kwargs.get('already_mapped', [])
    out_sam_path = os.path.abspath(os.path.expanduser(out_sam_path))

    seq_start = kwargs.get('seq_start', 0)
    seq_end = kwargs.get('seq_end', None)
    nthreads = kwargs.get('nthreads', 4)
    max_reads_per_chunk = kwargs.get('max_reads_per_chunk', -1)
    bowtie_flags = kwargs.get('bowtie_flags', '')

    if subprocess.call(['which', 'samtools']) != 0:
        raise Exception('samtools are not installed!')

    # Check for a typo from the publication.
    assert re.search('--score-min(\s+)-L', bowtie_flags) is None, (
        'The flag --score-min -L 0.6,0.2 in the original publication was a typo. '
        'The correct notation is --score-min L ... . Please fix the supplied '
        'flags.')

    temp_dir = os.path.abspath(os.path.expanduser(
        kwargs.get('temp_dir', tempfile.gettempdir())))
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)

    bash_reader = kwargs.get('bash_reader', None)
    if bash_reader is None:
        extension = fastq_path.split('.')[-1].lower()
        if extension == 'gz':
            if commandExists("pigz"):
                bash_reader = "pigz -dc"
            else:
                bash_reader = 'gunzip -c'
        else:
            bash_reader = 'cat'
    else:
        if not commandExists(bash_reader):
            bash_reader = os.path.abspath(os.path.expanduser(bash_reader))
            if not os.path.isfile(bash_reader.split()[0]):
                raise Exception(
                    'The bash reader is not found '
                    'at the specified location {0}.'.format(bash_reader))

    reading_command = bash_reader.split() + [fastq_path, ]

    # If bash reader is not 'cat', convert file to FASTQ first and
    # run iterative_mapping recursively on the converted file.

    output_is_bam = (out_sam_path.split('.')[-1].lower() == 'bam')
    bamming_command = ['samtools', 'view', '-bS', '-'] if output_is_bam else []

    # Split input files if required and apply iterative mapping to each
    # segment separately.
    if max_reads_per_chunk > 0:
        kwargs['max_reads_per_chunk'] = -1
        log.info('Split input file %s into chunks', fastq_path)
        chunked_files = _chunk_file(
            fastq_path,
            os.path.join(temp_dir, os.path.split(fastq_path)[1]),
            max_reads_per_chunk * 4)
        log.debug('%d chunks obtained', len(chunked_files))
        for i, fastq_chunk_path in enumerate(chunked_files):
            log.info('Run iterative_mapping recursively on %s', fastq_chunk_path)
            iterative_mapping(
                bowtie_path, bowtie_index_path, fastq_chunk_path,
                out_sam_path + '.%d' % (i + 1), min_seq_len, len_step,
                **kwargs)

            # Delete chunks only if the file was really chunked.
            if len(chunked_files) > 1:
                log.info('Remove the chunks: %s', ' '.join(chunked_files))
                os.remove(fastq_chunk_path)
        return

    # Convert input relative arguments to the absolute length scale.
    reading_process = subprocess.Popen(reading_command,
                                       stdout=subprocess.PIPE)
    reading_process.stdout.readline()
    raw_seq_len = len(reading_process.stdout.readline().strip())
    log.info('The length of whole sequences in the file: %d', raw_seq_len)
    reading_process.terminate()
    sleep()

    if kwargs.get('first_iteration', True):
        has_old_files = False
        for path in sorted(glob.glob(out_sam_path + '.*')):
            try:
                mapped_len = int(path[len(out_sam_path) + 1:])
                if ((mapped_len - min_seq_len) % len_step != 0) and (mapped_len != raw_seq_len):
                    has_old_files = True
            except:
                pass

        if has_old_files:
            raise Exception(
                'The output folder contains a SAM file mapped '
                'to a different length range. '
                'Most likely, this is an artifact of previous mappings.')

    if (seq_start < 0
        or seq_start > raw_seq_len
        or (seq_end and seq_end > raw_seq_len)):
        raise Exception('An incorrect trimming region is supplied: [%d, %d), '
                        'the raw sequence length is %d' % (
                            seq_start, seq_end, raw_seq_len))
    local_seq_end = min(raw_seq_len, seq_end) if seq_end else raw_seq_len

    if min_seq_len <= local_seq_end - seq_start:
        trim_5 = seq_start
        trim_3 = raw_seq_len - seq_start - min_seq_len
        local_out_sam = out_sam_path + '.' + str(min_seq_len)
        mapping_command = [
            bowtie_path, '-x', bowtie_index_path, '-q', '-',
            '-5', str(trim_5), '-3', str(trim_3), '-p', str(nthreads)
            ] + bowtie_flags.split()

        pipeline = []
        try:
            log.info('Reading command: %s', ' '.join(reading_command))
            pipeline.append(
                subprocess.Popen(reading_command, stdout=subprocess.PIPE, bufsize=-1))
            if bamming_command:
                log.info('Mapping command: %s', ' '.join(mapping_command))
                pipeline.append(
                    subprocess.Popen(mapping_command,
                                     stdin=pipeline[-1].stdout,
                                     stdout=subprocess.PIPE, bufsize=-1))

                log.info('Output formatting command: %s', ' '.join(bamming_command))
                pipeline.append(
                    subprocess.Popen(bamming_command,
                                     stdin=pipeline[-1].stdout,
                                     stdout=open(local_out_sam, 'w'), bufsize=-1))
            else:
                log.info('Mapping command: %s', ' '.join(mapping_command))
                pipeline.append(
                    subprocess.Popen(mapping_command,
                                     stdin=pipeline[-1].stdout,
                                     stdout=open(local_out_sam, 'w'), bufsize=-1))
            pipeline[-1].wait()
        finally:
            sleep()
            for process in pipeline:
                if process.poll() is None:
                    process.terminate()

        # Check if the next iteration is required.
        if len_step <= 0:
            return
        if min_seq_len + len_step > local_seq_end - seq_start:
            return

        # Recursively go to the next iteration.
        log.info('Save the unique aligments and send the '
                     'non-unique ones to the next iteration')
        reading_process = subprocess.Popen(reading_command,
                                       stdout=subprocess.PIPE,
                                       bufsize=-1)

        unmapped_fastq_path = os.path.join(
            temp_dir, os.path.split(fastq_path)[1] + '.%d' % min_seq_len + ".fastq.gz")

        num_total, num_filtered = _filter_unmapped_fastq(
            reading_process.stdout, local_out_sam, unmapped_fastq_path, in_filename=fastq_path)

        reading_process.stdout.flush()
        reading_process.stdout.close()
        reading_process.terminate()
        sleep()


        log.info(('{0} non-unique reads out of '
                  '{1} are sent the next iteration.').format(num_filtered,
                                                             num_total))
        kwargs['first_iteration'] = False
        if commandExists("pigz"):
            kwargs["bash_reader"] = "pigz -dc"
        else:
            kwargs["bash_reader"] = 'gunzip -c'


        iterative_mapping(bowtie_path, bowtie_index_path, unmapped_fastq_path,
                          out_sam_path,
                          min_seq_len=min_seq_len + len_step,
                          len_step=len_step, **kwargs)

        os.remove(unmapped_fastq_path)



def _find_rfrags_inplace(lib, genome, min_frag_size, side):
    '''Private: assign mapped reads to restriction fragments by
    their 5' end position.
    '''
    assert isinstance(genome, mirnylib.genome.Genome)  # make Pydev happy
    side = str(side)

    chrms = lib['chrms' + side]
        # setting to zero chromosomes that are over the limit of the genome
    removeMask = chrms >= genome.chrmCount
    chrms[removeMask] = -1
    lib['chrms' + side] = chrms
    cuts = lib['cuts' + side]
    cuts[removeMask] = -1
    lib['cuts' + side] = cuts

    rfragIdxs = np.zeros(len(chrms), dtype=np.int64)
    uprsites = np.zeros(len(chrms), dtype=np.int64)
    rsites = np.zeros(len(chrms), dtype=np.int64)
    downrsites = np.zeros(len(chrms), dtype=np.int64)

    # If the fragment was not mapped.
    rfragIdxs[chrms == -1] = -1
    rsites[chrms == -1] = -1
    uprsites[chrms == -1] = -1
    downrsites[chrms == -1] = -1

    badCuts = np.nonzero(cuts >= genome.chrmLens[chrms])[0]
    if len(badCuts) > 0:
        maxDev = np.max(cuts[badCuts] - genome.chrmLens[chrms[badCuts]])
        warnings.warn(
            ('\nDetermined many ({0}) reads that map after the end of chromosome!'
             '\n Maximum deviation is {1} bp ').format(len(badCuts), maxDev))
        if maxDev > 50:
            raise StandardError("Deviation is too large. Probably, genome mismatch.")
        cuts[badCuts] = np.array(genome.chrmLens[np.array(chrms[badCuts], dtype=int)] - 1, dtype=cuts.dtype)
    if len(badCuts) > 10000:
        raise StandardError("Determined too many (%s) reads that map after "
                            "the end of chromosome!" % len(badCuts))

    strands = lib['strands' + side]
    for chrm_idx in xrange(genome.chrmCount):
        all_rsites = np.r_[0, genome.rsites[chrm_idx]]
        idxs = (chrms == chrm_idx)

        # Find the indexes of the restriction fragment...
        rfragIdxs[idxs] = np.searchsorted(all_rsites, cuts[idxs]) - 1
        uprsites[idxs] = all_rsites[rfragIdxs[idxs]]
        downrsites[idxs] = all_rsites[rfragIdxs[idxs] + 1]
        rsites[idxs] = np.where(
            strands[idxs], downrsites[idxs], uprsites[idxs])

        too_close = (np.abs(rsites[idxs] - cuts[idxs]) <= min_frag_size)
        too_close_idxs = np.where(idxs)[0][too_close]
        rfragIdxs[too_close_idxs] += strands[too_close_idxs] * 2 - 1
        uprsites[too_close_idxs] = all_rsites[rfragIdxs[too_close_idxs]]
        downrsites[too_close_idxs] = all_rsites[rfragIdxs[too_close_idxs] + 1]
        rsites[too_close_idxs] = np.where(
            strands[too_close_idxs],
            downrsites[too_close_idxs],
            uprsites[too_close_idxs])





    lib['rfragIdxs' + side] = rfragIdxs
    lib['uprsites' + side] = uprsites
    lib['downrsites' + side] = downrsites
    lib['rsites' + side] = rsites


def _find_rfrags_inplace_new(lib, genome, min_frag_size, side):
    '''Private: assign mapped reads to restriction fragments by
    their 5' end position.
    '''
    assert isinstance(genome, mirnylib.genome.Genome)  # make Pydev happy
    side = str(side)

    chrms = lib['chrms' + side]
        # setting to zero chromosomes that are over the limit of the genome
    removeMask = chrms >= genome.chrmCount
    chrms[removeMask] = -1
    del removeMask
    lib['chrms' + side] = chrms
    cuts = lib['cuts' + side]
    cuts[removeMask] = -1
    lib['cuts' + side] = cuts

    rfragIdxs = np.zeros(len(chrms), dtype=np.int64)
    rsites = np.zeros(len(chrms), dtype=np.int64)

    # If the fragment was not mapped.
    rfragIdxs[chrms == -1] = -1
    rsites[chrms == -1] = -1


    badCuts = np.nonzero(cuts >= genome.chrmLens[chrms])[0]
    if len(badCuts) > 0:
        maxDev = np.max(cuts[badCuts] - genome.chrmLens[chrms[badCuts]])
        warnings.warn(
            ('\nDetermined many ({0}) reads that map after the end of chromosome!'
             '\n Maximum deviation is {1} bp ').format(len(badCuts), maxDev))
        if maxDev > 50:
            raise StandardError("Deviation is too large. Probably, genome mismatch.")
        cuts[badCuts] = np.array(genome.chrmLens[np.array(chrms[badCuts], dtype=int)] - 1, dtype=cuts.dtype)
    if len(badCuts) > 10000:
        raise StandardError("Determined too many (%s) reads that map after "
                            "the end of chromosome!" % len(badCuts))


    strands = lib['strands' + side]
    for chrm_idx in xrange(genome.chrmCount):
        all_rsites = np.r_[0, genome.rsites[chrm_idx]]
        idxs = (chrms == chrm_idx)

        # Find the indexes of the restriction fragment...
        rfragIdxs[idxs] = np.searchsorted(all_rsites, cuts[idxs]) - 1
        uprsites[idxs] = all_rsites[rfragIdxs[idxs]]
        downrsites[idxs] = all_rsites[rfragIdxs[idxs] + 1]
        rsites[idxs] = np.where(
            strands[idxs], downrsites[idxs], uprsites[idxs])

        too_close = (np.abs(rsites[idxs] - cuts[idxs]) <= min_frag_size)
        too_close_idxs = np.where(idxs)[0][too_close]
        rfragIdxs[too_close_idxs] += strands[too_close_idxs] * 2 - 1
        uprsites[too_close_idxs] = all_rsites[rfragIdxs[too_close_idxs]]
        downrsites[too_close_idxs] = all_rsites[rfragIdxs[too_close_idxs] + 1]
        rsites[too_close_idxs] = np.where(
            strands[too_close_idxs],
            downrsites[too_close_idxs],
            uprsites[too_close_idxs])





    lib['rfragIdxs' + side] = rfragIdxs
    lib['uprsites' + side] = uprsites
    lib['downrsites' + side] = downrsites
    lib['rsites' + side] = rsites


def _parse_ss_sams(sam_basename, out_dict, genome_db,
                   max_seq_len=-1, reverse_complement=False, save_seqs=False):
    """Parse SAM files with single-sided reads.
    """
    def _for_each_unique_read(sam_basename, genome_db, action):
        sam_paths = glob.glob(sam_basename + '.*')
        if not sam_paths:
            raise Exception('No SAM/BAM files with \'%s\' basename are found.' % sam_basename)


        for sam_path in sam_paths:

            samfile = pysam.Samfile(sam_path)

            # Make Bowtie's chromosome tids -> genome_db indices dictionary.
            tid2idx = {}
            for i in xrange(len(samfile.lengths)):
                chrm_rname = samfile.getrname(i)
                chrm_label = genome_db._extractChrmLabel(chrm_rname)
                if chrm_label in genome_db.label2idx:
                    tid2idx[i] = genome_db.label2idx[chrm_label]

            for read in samfile:
                # Skip non-mapped reads.
                if read.is_unmapped:
                    continue

                # Skip non-uniquely aligned.
                for tag in read.tags:
                    if tag[0] == 'XS':
                        break
                else:
                    # Convert Bowtie's chromosome tids to genome_db indices.
                    # Skip chromosomes that are not in the genome.
                    if read.tid in tid2idx:
                        read.tid = tid2idx[read.tid]
                        action(read)

    # Calculate reads statistics.
    def _count_stats(read):
        # In Python, function is an object and can have an attribute.
        # We are using the .cache attribute to store the stats.
        _count_stats.id_len = max(_count_stats.id_len,
                                  len(read.qname))
        _count_stats.seq_len = max(_count_stats.seq_len,
                                   len(read.seq))
        _count_stats.num_reads += 1
    _count_stats.id_len = 0
    _count_stats.seq_len = 0
    _count_stats.num_reads = 0
    _for_each_unique_read(sam_basename, genome_db, _count_stats)
    sam_stats = {'id_len': _count_stats.id_len,
                 'seq_len': _count_stats.seq_len,
                 'num_reads': _count_stats.num_reads}
    log.info(
        'Parsing SAM files with basename {0}, # of reads: {1}'.format(
            sam_basename, sam_stats['num_reads']))
    if max_seq_len > 0:
        sam_stats['seq_len'] = min(max_seq_len, sam_stats['seq_len'])

    if sam_stats['num_reads'] == 0:
        out_dict.update(
            {'chrms': [], 'strands': [], 'cuts': [], 'seqs': [], 'ids': []})
        return out_dict

    # Read and save each type of data separately.
    def _write_to_array(read, array, value):
        array[_write_to_array.i] = value

    def inc(function):
        function.i += 1

    # ...chromosome ids
    chrmBuf = np.zeros((sam_stats['num_reads'],), dtype=np.int8)
    strandBuf = np.zeros((sam_stats['num_reads'],), dtype=np.bool)
    cutBuf = np.zeros((sam_stats['num_reads'],), dtype=np.int64)
    idBuf = np.zeros(
        (sam_stats['num_reads'],), dtype='|S%d' % sam_stats['id_len'])

    _write_to_array.i = 0

    if save_seqs:
        seqBuf = np.zeros(
            (sam_stats['num_reads'],), dtype='|S%d' % sam_stats['seq_len'])

        _for_each_unique_read(sam_basename, genome_db,
            action=lambda read: (_write_to_array(read, chrmBuf, read.tid),
                                 _write_to_array(read, strandBuf, not read.is_reverse),
                                 _write_to_array(read, cutBuf, read.pos + (len(read.seq) if read.is_reverse else 0)),
                                 _write_to_array(read, idBuf, read.qname[:-2] if read.qname.endswith('/1') or read.qname.endswith('/2') else read.qname),
                                 _write_to_array(read, seqBuf, Bio.Seq.reverse_complement(read.seq) if read.is_reverse and reverse_complement else read.seq),
                                 inc(_write_to_array)))

        out_dict['chrms'] = chrmBuf
        out_dict["strands"] = strandBuf
        out_dict["cuts"] = cutBuf
        out_dict["ids"] = idBuf
        out_dict['seqs'] = seqBuf
    else:
        print "In a recent update by default we're not saving sequences!!!"
        print "use parse_sams(save_seqs=True) to save sequences"
        warnings.warn(RuntimeWarning("Since 14-01-20 we're not saving sequences by default"))
        _for_each_unique_read(sam_basename, genome_db,
            action=lambda read: (_write_to_array(read, chrmBuf, read.tid),
                                 _write_to_array(read, strandBuf, not read.is_reverse),
                                 _write_to_array(read, cutBuf, read.pos + (len(read.seq) if read.is_reverse else 0)),
                                 _write_to_array(read, idBuf, read.qname[:-2] if read.qname.endswith('/1') or read.qname.endswith('/2') else read.qname),
                                 inc(_write_to_array)))

        out_dict['chrms'] = chrmBuf
        out_dict["strands"] = strandBuf
        out_dict["cuts"] = cutBuf
        out_dict["ids"] = idBuf


    return out_dict


def parse_sam(sam_basename1, sam_basename2, out_dict, genome_db, save_seqs=False, **kwargs):
    '''Parse SAM/BAM files with HiC reads.

    Parameters
    ----------

    sam_basename1 : str
        A basename of SAM files with the mapped sequences of the first
        side of Hi-C molecules.

    sam_basename2 : str
        A basename of SAM files with the mapped sequences of the second
        side of Hi-C molecules.

    out_dict : dict-like
        A dict-like structure to store the library of matched HiC reads.

    genome_db : str or mirnylib.genome.genome
        A path to a folder with FASTA files or a genome object. It is used
        to convert Bowtie chromosome indices to internal indices.

    max_seq_len : int, optional
        The length the sequences are truncated to before saving
        into the library. The default value is -1, i.e. the sequences are
        not truncated.

    reverse_complement : bool, optional
        If True then the sequences of reads on the reversed strand will be
        reverse complemented. False by default.

    keep_ids : bool, optional
        If True then the IDs of reads are stored. False by default.

    enzyme_name : str, optional
        If specified, assign the reads to the restriction fragments with
        the fill_rsites() function.

        The name of the restriction enzyme. The full list of possible names
        can be found in Bio.Restriction.AllEnzymes. If 'auto' and genome_db
        has an enzyme set then use this enzyme.

    min_frag_size : int, optional
        The minimal distance between a cut site and a restriction site.
        Used only if enzyme_name is specified.

        If the actual distance is less or equal than minimal then the ultra-sonic
        fragment is assigned to the next restriction fragment in the direction
        of the read. Default is None, which means it is set to a half
        of the length of the restriction motif.
    '''

    max_seq_len = kwargs.get('max_seq_len', -1)
    reverse_complement = kwargs.get('reverse_complement', False)
    keep_ids = kwargs.get('keep_ids', False)
    enzyme_name = kwargs.get('enzyme_name', None)
    min_frag_size = kwargs.get('min_frag_size', None)

    if isinstance(genome_db, str):
        genome_db = mirnylib.genome.Genome(genome_db)
    assert isinstance(genome_db, mirnylib.genome.Genome)

    # Parse the single-sided  reads.
    ss_lib = {}
    ss_lib[1] = mirnylib.h5dict.h5dict()
    ss_lib[2] = mirnylib.h5dict.h5dict()

    log.info('Parse the first side of the reads from %s' % sam_basename1)
    _parse_ss_sams(sam_basename1, ss_lib[1], genome_db,
                   1 if not max_seq_len else max_seq_len, reverse_complement, save_seqs=save_seqs)

    log.info('Parse the second side of the reads from %s' % sam_basename2)
    _parse_ss_sams(sam_basename2, ss_lib[2], genome_db,
                   1 if not max_seq_len else max_seq_len, reverse_complement, save_seqs=save_seqs)

    # Determine the number of double-sided reads.
    all_ids = np.unique(np.concatenate((ss_lib[1]['ids'], ss_lib[2]['ids'])))
    tot_num_reads = all_ids.shape[0]
    if tot_num_reads == 0:
        log.warning(
            'The SAM files %s and %s do not contain unique double sided reads' %
            (sam_basename1, sam_basename2))

    # Pair single-sided reads and write into the output.
    for i in [1, 2]:
        sorting = np.searchsorted(all_ids, ss_lib[i]['ids'])
        for key in ss_lib[i].keys():
            # Create empty arrays if input is empty.
            if tot_num_reads == 0:
                out_dict[key + str(i)] = []
                continue

            # Don't save ids and seqs if not requested.
            if key == 'ids' and not keep_ids:
                continue
            if key == 'seq' and not max_seq_len:
                continue

            # The default value is -1 for an undefined cut site and chromosome
            # and 0 for other data.
            if key == 'cuts' or key == 'chrms':
                buf = -1 * np.ones(shape=tot_num_reads,
                                   dtype=ss_lib[i].value_dtype(key))
            else:
                buf = np.zeros(shape=tot_num_reads,
                               dtype=ss_lib[i].value_dtype(key))

            buf[sorting] = ss_lib[i][key]
            out_dict[key + str(i)] = buf
            del buf

    misc_dict = {}
    misc_dict['genome'] = {}
    misc_dict['genome']['idx2label'] = dict(genome_db.idx2label)
    misc_dict['genome']['label2idx'] = dict(genome_db.label2idx)
    out_dict['misc'] = misc_dict

    if not (enzyme_name is None):
        fill_rsites(out_dict, genome_db, enzyme_name, min_frag_size)

    return out_dict


def fill_rsites(lib, genome_db, enzyme_name='auto', min_frag_size=None):
    '''Assign the mapped reads to the restriction fragments.

    Parameters
    ----------

    lib : dict
        A library of mapped Hi-C molecules. Gets modified by the function.

    genome_db : str or mirnylib.genome.genome
        A path to the folder with genome sequences in FASTA format or
        a mirnylib.genome.genome object.

    enzyme_name : str
        A name of the restriction enzyme. The full list of possible names
        can be found in Bio.Restriction.AllEnzymes.

    min_frag_size : int
        The minimal distance between a cut site and a restriction site.
        If the actual distance is less than minimal then the ultra-sonic
        fragment is assigned to the next restriction fragment in the direction
        of the read.
    '''

    if isinstance(genome_db, str):
        genome_db = mirnylib.genome.Genome(genome_db)
    assert isinstance(genome_db, mirnylib.genome.Genome)

    if len(lib['chrms1']) == 0:
        return lib

    if enzyme_name == 'auto':
        if not genome_db.hasEnzyme():
            raise Exception('Set a restriction enzyme in the genome object or '
                            'supply its name')
    else:
        if enzyme_name not in Bio.Restriction.AllEnzymes:
            raise Exception('Enzyme is not found in the library: %s' %
                (enzyme_name,))
        genome_db.setEnzyme(enzyme_name)

    rsite_size = eval('len(Bio.Restriction.%s.site)' % genome_db.enzymeName)
    if min_frag_size is None:
        _min_frag_size = rsite_size / 2.0
    else:
        _min_frag_size = min_frag_size

    _find_rfrags_inplace(lib, genome_db, _min_frag_size, 1)
    _find_rfrags_inplace(lib, genome_db, _min_frag_size, 2)

    return lib
