## This file is part of OSPREY 3.0
## 
## OSPREY Protein Redesign Software Version 3.0
## Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
## 
## OSPREY is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License version 2
## as published by the Free Software Foundation.
## 
## You should have received a copy of the GNU General Public License
## along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
## 
## OSPREY relies on grants for its development, and since visibility
## in the scientific literature is essential for our success, we
## ask that users of OSPREY cite our papers. See the CITING_OSPREY
## document in this distribution for more information.
## 
## Contact Info:
##    Bruce Donald
##    Duke University
##    Department of Computer Science
##    Levine Science Research Center (LSRC)
##    Durham
##    NC 27708-0129
##    USA
##    e-mail: www.cs.duke.edu/brd/
## 
## <signature of Bruce Donald>, Mar 1, 2018
## Bruce Donald, Professor of Computer Science

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
