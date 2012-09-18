from mirnylib import genome
from hiclib import fragmentHiC 

genome_db = genome.Genome('data/hg19', readChrms=['#', 'X'])

fragments = fragmentHiC.HiCdataset(
    filename='data/fragment_dataset.hdf5',
    genome=genome_db,
    maximumMoleculeLength=500,
    override=True)

fragments.parseInputData(
    dictLike='data/mapped_reads.hdf5')

fragments.filterRsiteStart(offset=5)
fragments.filterDuplicates()
fragments.filterLarge()
fragments.filterExtreme(cutH=0.005, cutL=0)

fragments.saveHeatmap('heatmap-res-500k.hdf5', resolution=500000)
fragments.saveHeatmap('heatmap-res-1M.hdf5', resolution=1000000)
fragments.saveHeatmap('heatmap-res-2M.hdf5', resolution=2000000)

# Requires at least 8Gb of RAM:
fragments.saveHeatmap('heatmap-res-200k.hdf5', resolution=200000) 
