import setuptools
from distutils.core import setup, Extension
from Cython.Build import cythonize

# See https://stackoverflow.com/a/51272967/5188860
module1 = Extension('kmerGWAS',
                    sources = ['IBSpy/kmerGWAS/nucleotide.c', 'IBSpy/kmerGWAS/kmer_general.pyx',
                               'IBSpy/kmerGWAS/kmer_general_c.c'])

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="IBSpy", # Replace with your own username
    version="0.0.1",
    author="Ricardo H. Ramirez-Gonzalez",
    author_email="ricardo.ramirez-gonzalez@jic.ac.ul",
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
    #ext_modules=[module1]
    ext_modules = cythonize(
        [module1],
        compiler_directives={'language_level': "3"})
)
