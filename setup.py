
#!/usr/bin/env python3

import setuptools
from distutils.core import setup

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(name="fishlifeqc",
      version='0.1',
      long_description = readme,
      long_description_content_type='text/markdown',
      url='https://github.com/Ulises-Rosas/mfeprimer-py3',
      packages=['fishlifeqc'],
      data_files=[
          ('bin', [ 'ext_bin/raxmlHPC-SSE3_Darwin_64bit',
                    'ext_bin/raxmlHPC-PTHREADS-SSE3_Darwin_64bit',
                    'ext_bin/raxmlHPC-SSE3_Linux_64bit',
                    'ext_bin/raxmlHPC-PTHREADS-SSE3_Linux_64bit',
                    ])
      ],
      classifiers=[
          'Programming Language :: Python :: 3'
      ]
    )
