

# relevant docs for making python packages:
# https://packaging.python.org/tutorials/distributing-packages/
# https://pip.pypa.io/en/stable/user_guide/#installing-packages


# to run the osprey module in 'develop' mode (ie, where python points directly to the source files),
# run this script directly:
# $ python setup.py develop
# to uninstall the development version, run:
# $ python setup.py develop --uninstall

rootDir = '../'

# the jerkar build script also runs this setup.py to build the python distribution files
# but it will change the rootDir variable so paths still resolve correctly
# the jerkar build script copies files as needed to create the following dir structure for this script:
# build/output/python/
#    osprey/
#       *.py
#       osprey-version.jar
#       natives/
#    setup.py
#    README.rst
#    LICENSE.txt


import setuptools
import os


# read the osprey version
with open(os.path.join(rootDir, 'resources/config/version'), 'r') as file:
	version = file.read()


setuptools.setup(
	name='osprey',
	version=version,
	description='Open-Source Protein Redesign for You',
	url='https://github.com/donaldlab/OSPREY_refactor',
	packages=['osprey'],
	install_requires=['JPype1'],
	package_data={
		'osprey': ['*.jar', 'natives/*', '*.rst']
	},
	classifiers=[
		'Programming Language :: Python :: 2.7',
		'Programming Language :: Python :: 3',
		'Operating System :: OS Independent'
	]
)

# TODO: get a pypi account and upload our python package when it's ready for release
# then installs can be as easy as `pip install osprey`

