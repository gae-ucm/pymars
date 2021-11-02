from setuptools import setup, find_packages
from os import path
from pymars.version import *

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='pymars',
      version=get_version_pypi(),
      author="gae-ucm",
      author_email="tmiener@ucm.es",
      description='Python-based package for IACT event analysis of the MAGIC telescopes.',
      long_description=long_description,
      long_description_content_type='text/x-rst',
      url='https://github.com/gae-ucm/pymars',
      license='GPL-3.0',
      packages=['pymars'],
      install_requires=[
          'matplotlib',
          'numpy>=1.15.0',
          'scipy',
          'jupyter',
          'pandas',
          'pint-pulsar',
          'ctaplot',
          'uproot',
          'tables',
          ],
      entry_points = {
        'console_scripts': ['pymars-podie=pymars.podie:main',
                            'pymars-pulsar=pymars.pulsar:main',
                            'pymars-merger=pymars.merger:main',
                            'pymars-plotRES=pymars.plotRES:main',
                            'pymars-plotROC=pymars.plotROC:main',
                            'pymars-plotCRABRES=pymars.plotCRABRES:main']
      },
      include_package_data=True,
      dependencies=[],
      dependency_links=[],
      zip_safe=False)
