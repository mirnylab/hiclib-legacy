from setuptools import setup

setup(
    name='hiclib',
    url='http://mirnylab.bitbucket.org/hiclib/index.html',
    description=('Hi-C data analysis library.'),
    package_dir={'': 'src'},
    packages=['hiclib'],
    install_requires=[
        'biopython',
        'mirnylib',
        'pysam',
    ],
    dependency_links=[
        'https://bitbucket.org/mirnylab/mirnylib/get/tip.tar.bz2'
    ],
)