#from distutils.core import setup
#from distutils.extension import Extension
#
from setuptools import setup
from setuptools.extension import Extension

import subprocess

import os

#update version
args = 'git describe --tags'
p = subprocess.Popen(args.split(), stdout=subprocess.PIPE)
long_version = p.communicate()[0].decode("utf-8").strip()
spl = long_version.split('-')

if len(spl) == 3:
    main_version = spl[0]
    commit_number = spl[1]
    version_hash = spl[2]
    version = f'{main_version}.dev{commit_number}'
else:
    version_hash = '---'
    version = long_version
    
# version = "1.0" # pipeline with database

# Set this to true to add install_requires to setup
# Turned off for incremental builds as it kills "reload(mastquery.query)" 
if 0:
    install_requires=[
         'astropy>=2.0.0',
         'scipy',
         'numpy>=1.10.2',
         'matplotlib>=2.0.2',
         'mpdaf>=1.0']
else:
    install_requires = []    
    
#lines = open('grizli/version.py').readlines()
version_str =f"""# git describe --tags
__version__ = "{version}"
__long_version__ = "{long_version}"
__version_hash__ = "{version_hash}"
"""

fp = open('mospipe/version.py','w')
fp.write(version_str)
fp.close()
print('Git version: {0}'.format(version))

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "mospipe",
    version = version,
    author = "Gabriel Brammer",
    author_email = "gbrammer@gmail.com",
    description = "MOSFIRE pipeline",
    license = "MIT",
    url = "https://github.com/gbrammer/mospipe",
    download_url = "https://github.com/gbrammer/mospipe/tarball/{0}".format(version),
    packages=['mospipe'],
    classifiers=[
        "Development Status :: 1 - Planning",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    install_requires=install_requires,
    package_data={'mospipe': ['data/*']},
)
