import os
import logging

import numpy as np

import mirnylib.genome
import mirnylib.h5dict

lib = mirnylib.h5dict.h5dict('./tmp/insilico_mapped_reads.hdf5')

# Read the true positions of reads from the IDs.
chrms1 = np.zeros(shape=lib['chrms1'].shape, dtype=int)
chrms2 = np.zeros(shape=lib['chrms1'].shape, dtype=int)
cuts1 = np.zeros(shape=lib['chrms1'].shape, dtype=int)
cuts2 = np.zeros(shape=lib['chrms1'].shape, dtype=int)
lens1 = np.zeros(shape=lib['chrms1'].shape, dtype=int)
lens2 = np.zeros(shape=lib['chrms1'].shape, dtype=int)
ids = np.where(lib['ids1'] != '', lib['ids1'], lib['ids2'])

for i in range(len(lib['chrms1'])):
    _, chrm1, cut1, len1, chrm2, cut2, len2, molecule_length = ids[i].split('_')
    chrms1[i] = int(chrm1)
    chrms2[i] = int(chrm2)
    cuts1[i] = int(cut1)
    cuts2[i] = int(cut2)
    lens1[i] = int(len1)
    lens2[i] = int(len2)
molecule_length = int(molecule_length)

first_mapped = (lib['chrms1'] >= 0)
first_mapped_corr = (chrms1==lib['chrms1']) * (cuts1==lib['cuts1'])
first_mismapped = first_mapped * (first_mapped_corr == False) 
first_same_frag = (first_mismapped * (lib['chrms1'] == chrms2)
                   * (np.abs(lib['cuts1'] - cuts2) < 1.10 * molecule_length))
second_mapped = (lib['chrms2'] >= 0)
second_mapped_corr = (chrms2==lib['chrms2']) * (cuts2==lib['cuts2'])
second_mismapped = second_mapped * (second_mapped_corr == False) 
second_same_frag = (second_mismapped * (lib['chrms2'] == chrms1)
                   * (np.abs(lib['cuts2'] - cuts1) < 1.10 * molecule_length))

print 'The first side of Hi-C molecules:'
print 'Out of {0} mapped reads, {1} are mapped correctly.'.format(
      first_mapped.sum(), first_mapped_corr.sum())
print ('Out of {0} mismapped reads, {1} were too short and were '
       'mapped to the same restriction fragment as the second side. '
       'These reads will be filtered at a later stage of data analysis.').format(
       first_mismapped.sum(), first_same_frag.sum())
print

print 'The second side of Hi-C molecules:'
print 'Out of {0} mapped reads, {1} are mapped correctly.'.format(
      second_mapped.sum(), second_mapped_corr.sum())
print ('Out of {0} mismapped reads, {1} were too short and were '
       'mapped to the same restriction fragment at the first side. '
       'These reads will be filtered at a later stage of data analysis.').format(
       second_mismapped.sum(), second_same_frag.sum())
print

ds_mapped = (first_mapped * second_mapped)
ds_mapped_corr = (first_mapped_corr * second_mapped_corr)
ds_mapped_same_frag = (first_same_frag + second_same_frag)
ds_mismapped = (first_mapped * second_mapped  
    * ((first_mapped_corr == False) * (first_same_frag == False)
        + (second_mapped_corr == False) * (second_same_frag == False)))
ss = (first_mapped == False) + (second_mapped == False)

print 'Double-sided reads:'
print '{0} double-sided reads are mapped correctly,'.format(ds_mapped_corr.sum())
print '{0} have both sides mapped to the same restriction fragment and will be filtered out later,'.format(ds_mapped_same_frag.sum())
print '{0} are mismapped for other reason,'.format(ds_mismapped.sum())
print '{0} have only one side mapped.'.format(ss.sum())

