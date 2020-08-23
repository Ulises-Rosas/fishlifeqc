
#!/usr/bin/env python3

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
    bins = ['./ext_bin/raxml/raxmlHPC-SSE3_Darwin_64bit',
            './ext_bin/raxml/raxmlHPC-PTHREADS-SSE3_Darwin_64bit',
            './ext_bin/blast/blastn_Darwin_64bit',
            './ext_bin/blast/makeblastdb_Darwin_64bit']

elif myos == 'linux' or myos == "linux2":
    bins = ['./ext_bin/raxml/raxmlHPC-SSE3_Linux_64bit',
            './ext_bin/raxml/raxmlHPC-PTHREADS-SSE3_Linux_64bit',
            './ext_bin/blast/blastn_Linux_64bit',
            './ext_bin/blast/makeblastdb_Linux_64bit']

elif myos == 'win32':
    bins = ['./ext_bin/raxml/raxmlHPC-SSE3.exe',
            './ext_bin/raxml/raxmlHPC-PTHREADS-SSE3.exe',
            './ext_bin/blast/blastn.exe',
            './ext_bin/blast/makeblastdb.exe']
else:
    sys.stderr.write('Package does not work with %s operative system'  % myos)
    sys.stderr.flush()
    exit()

ext_modules = [
              Extension('runshell', sources = ["./fishlifeqc/runshell.c"]),
              Extension('fishlifeseq', sources= ["./fishlifeqc/sequtils.c"])
              ]

dependencies = [
                "boldminer",
                ]

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(name="fishlifeqc",
      version='0.2',
      long_description = readme,
      long_description_content_type='text/markdown',
      url='https://github.com/Ulises-Rosas/fishlifeqc',
      packages=['fishlifeqc'],
      ext_modules=ext_modules,
      data_files = [ ('bin', bins) ],
      install_requires = dependencies,
      entry_points={
        'console_scripts': [
            # 'pairedblast  = fishlifeqc.pairedblast:main',
            'fishlifeqc   = fishlifeqc.core_fishlife:main',
            ]
        },
      scripts=['./scripts/splitexonfiles.py'],
      classifiers=[
          'Programming Language :: Python :: 3'
      ]
    )
