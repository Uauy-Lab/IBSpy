import setuptools
from distutils.core import setup, Extension
from Cython.Build import cythonize
import unittest

def unit_tests():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='test_*.py')
    return test_suite

# See https://stackoverflow.com/a/51272967/5188860
module1 = Extension('kmerGWAS',
                   #sources = ['IBSpy/kmerGWAS/kmer_gwas.pyx']
                     sources = ['IBSpy/kmerGWAS/nucleotide.c',
                                'IBSpy/kmerGWAS/kmer_gwas.pyx',
                                'IBSpy/kmerGWAS/kmer_general.c',
                                'IBSpy/kmerGWAS/kmer_db.c']
                               )

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="IBSpy", # Replace with your own username
    version="0.0.6",
    author="Ricardo H. Ramirez-Gonzalez",
    author_email="ricardo.ramirez-gonzalez@jic.ac.uk",
    description="A package to detect IBS regions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Uauy-Lab/IBSpy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    test_suite='setup.unit_tests',
    #ext_modules=[module1]
    ext_modules = cythonize(
        [module1],
        compiler_directives={'language_level': "3"}),
    entry_points={  # Optional
        'console_scripts': [
            'IBSpy_wcount=IBSpy.bin.IBSpy_window_count:main',
            'IBSpy=IBSpy:main'
        ],
    }
)
