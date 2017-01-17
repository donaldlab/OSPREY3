
# TODO: make python installable package with dependencies
# depends on JPype1 (e.g., pip install JPype1)
# requires python3
# see: https://packaging.python.org/distributing/

import jvm, jpype

# NOTE: this var gets set by the build system during packaging
# so the release version of this script will point to the final jar file for osprey
# instead of the development classes folder
_ospreyPath = '../../bin'

c = None


def start(heapSizeMB=1024, enableAssertions=False):
	
	# start the jvm
	jvm.addClasspath(_ospreyPath)
	jvm.start(heapSizeMB, enableAssertions)

	# set up class factories
	global c
	c = jvm.Packages()
	c = jpype.JPackage('edu.duke.cs.osprey')

	# print the preamble
	print("OSPREY %s" % c.control.Main.Version)

