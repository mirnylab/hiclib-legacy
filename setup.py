from setuptools import setup, Extension
import numpy 
from Cython.Build import cythonize

ext_modules = []
from Cython.Distutils import build_ext

ext_modules += [Extension("hiclib.fastBinSearch", [ "binarySearch/fastBinSearch.pyx", "binarySearch/mycode.cpp"],
        language = "c++",
         include_dirs=[numpy.get_include()], 
         extra_compile_args = [ "-Ofast", "-fopenmp", "-std=c++11"],
         extra_link_args = [ "-Ofast", "-lgomp", "-std=c++11"]),   
         ]


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
