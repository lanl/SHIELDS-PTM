#!/usr/bin/env python

from setuptools import setup, find_packages

long_description = '# ptm_python\n'+\
                   '## SHIELDS PTM (Particle Tracing Model) preprocessing and post-processing'+\
                   'tools\n'+\
                   '\n'+\
                   '### Preprocessing\n'+\
                   '- Conversion of SWMF 3d fields to PTM format\n'+\
                   '\n'+\
                   '### Postprocessing\n'+\
                   '- Flux mapping tools\n'

setup(name='ptm_python',
      version='0.99.0',
      description='Pre- and post-processing library for PTM Particle Tracing Model',
      author='Steve Morley <smorley@lanl.gov>, Jesse Woodroffe',
      install_requires=['numpy', 'scipy', 'matplotlib', 'spacepy'],
      long_description=long_description,
      long_description_content_type="text/markdown",
      packages=find_packages(exclude=['tests']),
      classifiers=['Development Status :: 4 - Beta',
                   'License :: OSI Approved :: BSD License',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Atmospheric Science',
                   'Topic :: Scientific/Engineering :: Information Analysis',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   'Programming Language :: Python :: 3.8',
                   ],
     test_suite='test_ptm.py'
     )
