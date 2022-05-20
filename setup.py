import setuptools
# from distutils.core import setup, Extension
from Cython.Build import cythonize
import unittest
import distutils.cmd
from setuptools import Extension
import subprocess

class PyreverseCommand(distutils.cmd.Command):
    description = 'run pyreverse to update the UML diagram'
    user_options = []

    def initialize_options(self):
        """Set default values for options."""
        # Each user option must be listed here with their default value.
        #self.pylint_rcfile = ''
        pass

    def finalize_options(self):
        """Post-process options."""
        pass 

    def run(self):
        """Run command."""
        command = ['pyreverse']
        command.append('-ddocs')
        command.append('-opng')
        command.append('IBSpy')
        self.announce(
            'Running command: %s' % str(command),
            level=distutils.log.INFO)
        subprocess.check_call(command)

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
    name="IBSpy",
    version="0.4.10",
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
    cmdclass={
        'uml': PyreverseCommand,
    },
    # test_suite='tests',
    #ext_modules=[module1]
    ext_modules = cythonize(
        [module1],
        compiler_directives={'language_level': "3"}),
    entry_points={  # Optional
        'console_scripts': [
            'IBSpy_wcount=IBSpy.IBSpy_window_count:main',
            'IBSpy=IBSpy:main',
            'IBSplot=IBSpy.IBSpy_plot_window_results:main',
            'IBSpy_matrix=IBSpy.build_score_matrix:main',
            'IBSpy_haplotypes=IBSpy.call_haplotypes:main',
            'IBSpy_extract=IBSpy.extract_window:main'
        ],
    }
)
