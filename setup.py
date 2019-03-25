from setuptools import setup, Extension
import numpy 
from Cython.Build import cythonize

ext_modules = []
from Cython.Distutils import build_ext



setup(
    name='hiclib',
    url='http://mirnylab.bitbucket.org/hiclib/index.html',
    description=('Hi-C data analysis library.'),
    ext_modules=cythonize(ext_modules),    
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
