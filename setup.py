#!/usr/bin/env python

import os

try:
    from setuptools import setup
except ImportError as e:
    from distutils.core import setup

requirements = [
    'numpy>=1.14',
    'scipy',
    'astropy',
    'future-fstrings',
]

PACKAGE_NAME = 'pk_common_functions'
__version__ = '0.0.1'


with open("README.md", "r") as fh:
    long_description = fh.read()


setup(name=PACKAGE_NAME,
      version=__version__,
      description="Development Status :: 1 - Beta",
      long_description=long_description,
      long_description_content_type="text/markdown",
      author="P. Kamphuis",
      author_email="peterkamphuisastronomy@gmail.com",
      url="https://github.com/PeterKamphuis/Extract_Scaleheights",
      packages=[PACKAGE_NAME],
      python_requires='>=3.6',
      install_requires=requirements,
      include_package_data=False,
      # package_data - any binary or meta data files should go into MANIFEST.in
    
      license="GNU GPL v3",
      classifiers=[
          "Development Status :: 4 - Beta",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: POSIX :: Linux",
          "Programming Language :: Python :: 3",
          "Topic :: Scientific/Engineering :: Astronomy"
      ]
      )
