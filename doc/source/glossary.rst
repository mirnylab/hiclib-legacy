Glossary
========

General Hi-C terms
------------------

restriction fragment - 
    a fragment between two sites of restriction
    There are around 800k fragments in the Human or Mouse genome. 

ligation junction - 
    a site of enzymatic ligation between two restriction 
    fragments, forming an interaction pair.
    For strictly enzymatic cuts, it has a canonical sequence,
    determined by the restrictase. In reality, junction sites 
    may have different sequences due to random breaks. 

ultra-sonic cut site - 
    position of ultra-sonic shear site of a Hi-C molecule. 
    Also, the first position of the mapped sequence for each side of the read. 

ultra-sonic fragment (a.k.a. Hi-C molecule) -
    a DNA fragment sequenced in Hi-C experiment. The Hi-C experiment is designed
    to sequence only fragments that are formed by two ultra-sonic cuts and
    contain a ligation junction. However, the specificity is not 100%,
    so the collection of sequenced molecules may contain molecular by-products such as 
    dangling ends or self-circles.

Hi-C library - a collection of sequenced Hi-C molecules

read - the sequences of both ends (sides) of a Hi-C molecule

sides of a read - 
    a sequence of one end of a Hi-C molecule. Usually, the sides
    are numbered, e.g. "the first side" and "the second side". 
    Usually, the second side is slightly lower-quality due to sequencing 
    protocol. 

strand - 
    the genomic strand to which the side of a read was mapped. 
    Strand = True correspond to increase of genomic position in the direction of the read, 
    strand = False - decreasing. 

upstream, downstream - here, indicates genomic position with higher/lower coordinate along the chromosome. 
        

single-sided read - 
    read with one uniquely-mapped side. We don't distinguish if the other side was
    not-mappable, non-uniquely-mappable or error (e.g. low-quality). 

double-sided read - read with both sides mapped

dangling end -  Double-sided read, most likely corresponding to an un-ligated molecular product. 

self-circle - a self-circulized ligation product, i.e. a double-sided read, 
connecting two ends of a restriction fragment. 
        

Programming terms
-----------------

chromosome label - basename of a chromosomal fastq file, i.e. chr1, chr2, ... , chrX, chrM, chr4_random

chromosome ID - 
    zero-based number of chromosome in the list of chromosome labels. 
    E.g. for human genome: chr1 - 0; chr2 - 1;...  chr22 - 21; chrX - 22; chrY - 23; chrM - 24. 
    If any chromosomes are not used for the analysis, above numbers will change.
        
restriction fragment ID - 
    unique ID of a restriction fragment in the genome.
    Is calculated using the following formula: 
    rfragID = (fragStart + fragEnd)/2  + fragChrom * (max(chromLens) + 1000)
    Note that we use integer division.
    Extra 1000 has no particular meaning.
        

Binned data terms
-----------------

Genomic position - zero-based position along the chromosome

Bin - 
    region of a genome at a given resolution. Bins are zero-based, and start/end at exact multiples of resolution. 
    At 1MB resolutions bins for each would be [0,999999], [1000000,1999999],[x000000,chromLength-1]


        
