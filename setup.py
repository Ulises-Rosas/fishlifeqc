#!/usr/bin/env python

import sys
import platform
import setuptools
from distutils.core import setup, Extension


if platform.architecture()[0] != '64bit':
    sys.stderr.write('Architecture requires 64bit')
    sys.stderr.flush()
    exit()

myos = sys.platform

if myos == 'darwin':
    bins = [
        './ext_bin/blast/blastn_Darwin_64bit',
        './ext_bin/blast/makeblastdb_Darwin_64bit',
        './ext_bin/raxml/raxmlHPC-PTHREADS-SSE3_Darwin_64bit',
        './ext_bin/consel/darwin/seqmt',
        './ext_bin/consel/darwin/makermt',
        './ext_bin/consel/darwin/consel',
        './ext_bin/consel/darwin/catpv'
            ]

elif myos == 'linux' or myos == "linux2":
    bins = [
        './ext_bin/blast/blastn_Linux_64bit',
        './ext_bin/blast/makeblastdb_Linux_64bit',
        './ext_bin/raxml/raxmlHPC-PTHREADS-SSE3_Linux_64bit',
        './ext_bin/consel/linux/seqmt',
        './ext_bin/consel/linux/makermt',
        './ext_bin/consel/linux/consel',
        './ext_bin/consel/linux/catpv'
            ]

elif myos == 'win32':
    bins = [
        './ext_bin/blast/blastn.exe',
        './ext_bin/blast/makeblastdb.exe',
        './ext_bin/raxml/raxmlHPC-PTHREADS-SSE3.exe'
            ]
else:
    sys.stderr.write('Package does not work with %s operative system'  % myos)
    sys.stderr.flush()
    exit()

ext_modules = [
              Extension('runshell', sources = ["./fishlifeqc/runshell.c"], extra_compile_args=['-std=c99']),
              Extension('fishlifeseq', sources= ["./fishlifeqc/sequtils.c"], extra_compile_args=['-std=c99'])
              ]

dependencies = [
                "boldminer", # own package
                "fishlifetraits>=0.4.2", # own package
                'dendropy==4.4.0'
                ]

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(name="fishlifeqc",
      version='1.2.3',
      author='Ulises Rosas',
    #   long_description = readme,
    #   long_description_content_type='text/markdown',
      url='https://github.com/Ulises-Rosas/fishlifeqc',
      packages=['fishlifeqc', 'qcutil'],
      ext_modules=ext_modules,
      data_files = [ ('bin', bins) ],
      install_requires = dependencies,
      zip_safe = False,
      entry_points={
        'console_scripts': [
            'fishlifeqc   = fishlifeqc.core_fishlife:main',
            'qcutil = qcutil.cli:main'
            ]
        },
      scripts=[
          './scripts/splitexonfiles.py',
          './scripts/merge.py',
          './fishlifeqc/concatenate.py',
          './scripts/codon_partition.py',
          './scripts/cons_trees.py',
          './scripts/collapse.py'
          ],
      classifiers=[
          'Programming Language :: Python :: 3'
      ]
    )
