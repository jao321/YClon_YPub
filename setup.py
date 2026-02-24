#!/usr/bin/env python3
import os
import sys

if sys.version_info < (3,8,0):
    sys.exit('Python 3.8.0 or newer is required.\n')

try:
    from setuptools import setup, find_packages
except ImportError:
    sys.exit('You need to install setuptools so you can install YClon.\n')

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
setup(
    name='YClon',
    version='2.0.2',
    description='Set scripts to cluster AIRR repertoire into clonotypes',
    maintainer=['Jay Gervasio'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    maintainer_email=['joaodgervasio@@ufmg.br'],
    include_package_data=True,
    packages=find_packages(),
    install_requires=[
        'alive_progress',
        'joblib',
        'scipy>=1.11',
        'pandas>=2.1',
        'numpy>=1.26',
        'scikit-learn>=1.3',
        'rich'
    ],
    entry_points={
      'console_scripts': [
         'YClon=YClon.YClon:main',
         'YPub=YClon.YPub:main',
      ],
    }
)
