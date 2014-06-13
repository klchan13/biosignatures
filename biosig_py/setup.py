"""

biosignatures: elemental signature analysis of cyanobacteria
------------------------------------------------------------

"""
from distutils.core import setup

setup(name='biosignatures',
      packages = ['biosignatures'],
      package_dir = {'biosignatures':'.'},
      package_data={'biosignatures': ['./*.py']})
