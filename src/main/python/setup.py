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
import sys


# relevant docs for making python packages:
# https://packaging.python.org/tutorials/distributing-packages/
# https://pip.pypa.io/en/stable/user_guide/#installing-packages
# https://pypi.org/classifiers/


# this script gets run in multiple folders by different gradle tasks,
# so we need a configurable root
rootDir = '../../../'

# when run in the current folder by the gradle task `pythonDevelop`, we keep this rootDir,
# but the gradle task `pythonWheel` will re-write rootDir when it copies this script to the build dir


# read the osprey version
with open(os.path.join(rootDir, 'build/osprey-version'), 'r') as file:
	version = file.read()


# collect the python package classifiers
classifiers = [
	'Programming Language :: Python :: 3',
	'License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)'
]

# get OS details
if sys.platform in ('win32', 'cygwin'):
	classifiers += ['Operating System :: Microsoft :: Windows']
elif sys.platform == 'darwin':
	classifiers += ['Operating System :: MacOS']
else:
	classifiers += ['Operating System :: POSIX :: Linux']


setuptools.setup(
	name='osprey',
	version=version,
	description='Open-Source Protein Redesign for You',
	url='https://github.com/donaldlab/OSPREY_refactor',
	packages=['osprey'],
	python_requires='>=3',
	install_requires=['JPype1==1.5.0'],
	package_data={
		'osprey': [
			'lib/*.jar',
			'*.md',
			'*.txt',
			# the runtime's support for globs leaves much to be desired...
			#'jre/**/*'
			'jre/*', 'jre/*/*', 'jre/*/*/*', 'jre/*/*/*/*', 'jre/*/*/*/*/*', 'jre/*/*/*/*/*/*',
			#'progs/**/*'
			'progs/*', 'progs/*/*', 'progs/*/*/*', 'progs/*/*/*/*', 'progs/*/*/*/*/*', 'progs/*/*/*/*/*/*'
		]
	},
	classifiers=classifiers
)
