import os

from io import open
from distutils.util import convert_path
from setuptools import setup, find_packages

setup_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(setup_directory, 'README.rst'), encoding='utf-8') as readme_file:
    long_description = readme_file.read()

init_data = {}
init_pathname = convert_path('aucome/__init__.py')
with open(init_pathname) as init_file:
    exec(init_file.read(), init_data)

setup(name='aucome',
      description='Automatic Comparison of Metabolism',
      long_description=long_description,
      version=init_data['__version__'],
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
            'mpwt==0.6.0',
            'padmet==4.0',
            'eventlet==0.31.0',
            'requests==2.22.0',
            'rpy2==3.0.5',
            'supervenn==0.1.13',
      ],
      entry_points={
          'console_scripts': [
              'aucome = aucome.__main__:main'
          ]
      },
)
