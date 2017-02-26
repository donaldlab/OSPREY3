
:orphan:

.. _constructors:

Accessing Java Classes from Python
==================================

Almost all of Osprey's code is written in Java because of Java's runtime performance
and its rich development tools. This Python API was created to make Osprey's sophisticated
computational libraries available with the ease and flexibilty of Python scripting.

Most users of Osprey will not have to call constructors on Java classes directly.
The ``osprey`` module (see the :ref:`api_reference`) provides many convenient
functions for creating Osprey objects (with appropriately Pythonic conventions)
and are sufficient for completing most designs.

However, more advanced users may wish to invoke some of Osprey's more 'behind-the-scenes' code
from a Python script, even though that code has not been explicitly designed for Python access.
To create an instance of an Osprey Java class in a Python script, use the
class factory in the osprey module, which is named simply ``c``::

	import osprey
	instance = osprey.c.SomeJavaClass()

Osprey uses a library called `JPype <http://jpype.readthedocs.io/en/latest/>`_ to allow Java classes
to be called from Python. The ``osprey.c`` object is an instance of JPype's JPackage class that points
to the Java package ``edu.duke.cs.osprey``. For example, to instantiate a hypothetical Osprey class
``edu.duke.cs.osprey.package.Class``::

	instance = osprey.c.package.Class()

For lots of instantiations, try assigning the class constructor to a variable::

	Class = osprey.c.package.Class
	instance1 = Class()
	instance2 = Class()
	instance3 = Class()
	...

Java class constructors can be called with arguments too::

	ParameterizedClass = osprey.c.package.ParameterizedClass
	instance1 = ParameterizedClass('argument')

.. warning:: The Java language supports method overloading, but the Python language does not.
	This means it's tricky to figure out which Java method overload should be called by the
	Python script. JPype does a lot of work to figure this out based on the Python arguments
	passed to the method, but JPype might not always find a match (and raise a
	``RuntimeError: No matching overloads found``), or it might not pick the method you
	intended. For more information on how to get the right method overload, read
	`how JPype does automatic type conversion
	<http://jpype.readthedocs.io/en/latest/userguide.html#type-conversion>`_.


Accessing Java Inner Classes
----------------------------

The convention for referencing an inner class in Java is::

	package.name.OuterClass$InnerClass

but this is not valid Python syntax, so it can't be used with the JPackage classes.
To access inner classes, Osprey provides a helper method.

.. autofunction:: osprey.jvm.getInnerClass

To call a constructor on the hypothetical Osprey class ``edu.duke.cs.osprey.package.Outer$Inner``::

	import osprey
	instance = osprey.jvm.getInnerClass(osprey.c.package.Outer, 'Inner')

