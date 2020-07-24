
#!/usr/bin/env python3

import sys
import platform
import setuptools
from distutils.core import setup, Extension

module = Extension('runshell', sources= ["./fishlifeqc/runshell.c"])

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

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(name="fishlifeqc",
      version='0.2',
      long_description = readme,
      long_description_content_type='text/markdown',
      url='https://github.com/Ulises-Rosas/fishlifeqc',
      packages=['fishlifeqc'],
      ext_modules=[module],
      data_files = [ ('bin', bins) ],
      classifiers=[
          'Programming Language :: Python :: 3'
      ]
    )
