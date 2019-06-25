import os

from io import open
from setuptools import setup

setup_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(setup_directory, 'README.rst'), encoding='utf-8') as readme_file:
    long_description = readme_file.read()

setup(name='aucome',
      description='Automatic Comparison of Metabolism',
      long_description=long_description,
      version='0.0.4',
      url='https://github.com/AuReMe/aucome',
      author='A. Belcour',
      author_email='arnaud.belcour@gmail.com',
      classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Audience
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',

        # Environnement, OS, languages
        'Programming Language :: Python :: 3'
      ],
      packages=['aucome'],
      install_requires=[
            'docopt',
            'eventlet',
            'mpwt',
            'padmet',
            'requests'

      ],
      entry_points={
          'console_scripts': [
              'aucome = aucome.__main__:main'
          ]
      },
)