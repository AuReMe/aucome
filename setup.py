#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

setup(name='aucome',
      url='https://github.com/aureme/aucome',
      license='GPLv3+',
      description='Automatic Comparison of Metabolism',
      long_description=
      'AuCoMe is a Python3 workflow that aims at reconstructing homogeneous metabolic networks from genomes with heterogeneous levels of annotations. More information on usage and troubleshooting on Github: https://github.com/aureme/aucome',
      author='AuReMe',
      author_email='gem-aureme@inria.fr',
      python_requires='>=3.6',
      classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
          'Development Status :: 3 - Alpha',       
        # Audience
          'Intended Audience :: Science/Research',
          'Intended Audience :: Developers',
          'License :: OSI Approved',
          'Natural Language :: English',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          
        # Environnement, OS, languages
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7'
      ],
      packages=['aucome'],
      install_requires=['matplotlib', 'mpwt', 'padmet', 'rpy2==3.0.5',
                        'seaborn', 'supervenn', 'tzlocal'],      
      entry_points={
          'console_scripts': [
              'aucome = aucome.__main__:main'
          ]
      },
)
