## This file is part of OSPREY.
## 
## OSPREY is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
## 
## OSPREY is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.

import setuptools
import os


# relevant docs for making python packages:
# https://packaging.python.org/tutorials/distributing-packages/
# https://pip.pypa.io/en/stable/user_guide/#installing-packages


# this script gets run in multiple folders by different gradle tasks,
# so we need a configurable root
rootDir = '../'

# when run in the current folder by the gradle task `pythonDevelop`, we keep this rootDir
# the gradle task `pythonBdist` will re-write rootDir when it copies this script to the build dir

# the build dir should look like this when setup.py is called:
# build/python/bdist/
#    osprey/
#       *.py
#       lib/
#    setup.py
#    README.rst
#    LICENSE.txt


# read the osprey version
with open(os.path.join(rootDir, 'resources/config/version'), 'r') as file:
	version = file.read()


setuptools.setup(
	name='osprey',
	version=version,
	description='Open-Source Protein Redesign for You',
	url='https://github.com/donaldlab/OSPREY_refactor',
	packages=['osprey'],
	python_requires='>=2.7,<3',
	install_requires=['JPype-py2>=0.5.8'],
	package_data={
		'osprey': ['lib/*.jar', 'README.rst', 'LICENSE.txt']
	},
	classifiers=[
		'Programming Language :: Python :: 2.7',
		'Operating System :: OS Independent'
	]
)
