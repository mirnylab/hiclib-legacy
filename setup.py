from setuptools import setup

setup(
    name='hiclib',
    url='http://mirnylab.bitbucket.org/hiclib/index.html',
    description=('Hi-C data analysis library.'),
    packages=['hiclib'],
    install_requires=[
        'numpy>=1.6',
        'scipy',
        'matplotlib',
        'biopython',
        'mirnylib',
        'pysam',
    ],
    dependency_links=[
        'https://bitbucket.org/mirnylab/mirnylib/get/tip.tar.bz2'
    ],
)