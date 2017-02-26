
OSPREY
======

*Open-Source Protein REdesign for You*

OSPREY is developed and maintained by the `Donald Lab`_
in the `Department of Computer Science`_
at `Duke University`_.

.. _Donald Lab: http://www.cs.duke.edu/donaldlab/home.php
.. _Department of Computer Science: http://www.cs.duke.edu
.. _Duke University: https://www.duke.edu/

Osprey is released under a permissive GPL-style license. See the
`full Osprey license`_ for more information. 

.. _full Osprey license: http://www.cs.duke.edu/donaldlab/software/osprey/osprey.2.2/license.pdf


Installation
------------

Download an Osprey distribution zip, extract it to your favorite folder,
then run the ``setup.py`` install script with the ``install`` command::

    $ python setup.py install

.. note:: You may need super-user privileges to install Python packages. In linux, try the ``sudo`` command.


Running Osprey
--------------

using Python scripts
~~~~~~~~~~~~~~~~~~~~

Python scripting is the preferred way of using Osprey due to its simplicity and flexibilty.
To run Osprey from a Python script:

.. code:: python

	import osprey
	osprey.start()
	
	# run osprey commands, e.g.
	osprey.printGpuInfo()
	
For more information about Python scripting with Osprey, see the `Python documentation`_.

.. _Python documentation: TODO

.. note:: TODO: make the docs link work


using the command-line interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Python interface to Osprey represents a significant improvement in the user interface over the
older command-line interface, and new Osprey projects should consider using the Python interface
rather than the command-line interface.

However, for backwards compatibility, the command-line interface is still provided, although
it may not receive feature updates in the future.

Run Osprey at the command line by navigating to the folder where you extracted the Osprey distribution.
Then navigate to the ``osprey`` folder. Then entering the following command into a shell::

    $ java -jar /path/to/osprey.jar [commands]
    
.. note:: To use GPU acceleration, you'll need to need to tell Java where to find the operating
	system-specific GPU libraries. The Python interface does this automatically, but on the command line,
	you'll need to supply the additional JVM argument::
	
		-Djava.library.path=natives
		
	For example, to run the ``GpuInfo`` command::
	
		$ java -Djava.library.path=natives -jar osprey.jar GpuInfo

Try the ``help`` command to show a list of commands recognized by Osprey.

To show the version of your Osprey installation, try::

    $ java -jar osprey.jar version

To run a GMEC-based protein design, try::

    $ java -jar osprey.jar FindGMEC /path/to/config1 /path/to/config2 ...
    
For more information about the Osprey command-line interface and config file format, see
the `command-line interface documentation`_.

.. _command-line interface documentation: TODO

.. note:: TODO: make this doc link work too


Contributing
------------

Osprey is open-source software and contributions are welcome.

Download project from Github
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, clone the git repository to a folder of your choosing::

	$ mkdir osprey
	$ cd osprey
	$ git clone git@github.com:donaldlab/OSPREY_refactor.git .


Setup Eclipse build path
~~~~~~~~~~~~~~~~~~~~~~~~

Osprey uses the `Jerkar`_ build system. First `install Jerkar if needed`_,
then run this in a shell at the project root::

    $ jerkar eclipse#generateFiles

This will download the Java libraries needed by Osprey and generate the ``.classpath``
file for Eclipse. Then Eclipse can build the project and you’ll be able to launch Osprey from
within Eclipse for debugging.

Jerkar also provides an `Eclipse plugin`_ to make running Jerkar easier.

.. _Jerkar: http://project.jerkar.org
.. _install Jerkar if needed: http://project.jerkar.org/documentation/latest/getting_started.html
.. _Eclipse plugin: https://github.com/jerkar/eclipsePlugin4Jerkar


Build Osprey distribution
~~~~~~~~~~~~~~~~~~~~~~~~~

To build an Osprey distribution zip file, simply run in the Osprey project folder::

	$ jerkar doDist

The distribution zip file will be saved to ``build/output``.


Documentation
~~~~~~~~~~~~~

Osprey is mostly implemented in Java, but exposes a Python API as its user interface.
For the most part, this Python API is a very thin wrapper around the Java classes and is
dynamically generated using `JPype`_. The goal of this arrangement is use the Java language
and the JVM to enable rapid development of high-performance code, but provide end users with
the ease and flexibility of Python scripting. As much as possible, Osprey should look and
feel like a native Python application, even though it's really not.

.. _JPype: http://jpype.readthedocs.io/en/latest/

Most of the documentation for Osprey exists in the form of Javadoc comments, to make it easier
to keep the documentation up-to-date with code changes. However, the documentation the users see
is created with the typical Python documentation toolchain, `Sphinx`_. Following Python conventions,
Osprey's documentation outside of javadoc comments is written in the `ReStructured Text (RST)`_
format and can be found in the ``python/doc`` folder. For the javadoc comments, Osprey contains
a custom Sphinx extension (at ``python/doc/javadoc.py``) that converts javadoc comments into RST
documentation much like the `autodoc extension to Sphinx`_.

.. _Sphinx: http://www.sphinx-doc.org/en/stable/
.. _ReStructured Text (RST): https://en.wikipedia.org/wiki/ReStructuredText
.. _autodoc extension to Sphinx: http://www.sphinx-doc.org/en/stable/ext/autodoc.html

To build the documentation for Osprey, run the Sphinx tool from the ``doc`` folder::

	$ cd python/doc
	$ make html

.. note:: For quick edit-compile-test cycles when editing documentation, it's helpful
	to run ``make clean`` before ``make html`` which makes sure all documentation is refreshed
	regardless of which RST documents have been recently edited. e.g.::
	
		$ make clean && make html
		
.. warning:: Sphinx can detect problems with the documentation during building. When this happens,
	these problems will be reported to the console, usually in red text.
	These warning messages usually indicate something is missing or incorrect
	in the documentation, and that the underlying problems should be fixed before
	the documentation is released.

Then open the ``python/doc/_build/html/index.html`` file in your browser to view the documentation.

Osprey's javadoc extension to Sphinx provides a few directives and roles to allow referring to
Java classes, fields, methods, and javadoc comments from the RST documentation:


Sphinx Directives
~~~~~~~~~~~~~~~~~

**.. javaclass:: java_class_reference**

	where ``java_class_reference`` is the fully-qualified name of a Java class, e.g.::
	
		package.Class
		package.OuterClass$InnerClass
		
	This directive will automatically scan the source code for the specified class and show
	all the public constructors, methods, and fields for the class. Javadoc comments will be
	shown with the constructors, methods, arguments, fields, etc, and Java type information
	will be shown in the documentation where possible.
	
	.. note:: When the java reference is prefixed with a ``.``, the package ``edu.duke.cs.osprey``
		is automatically inferred. Therefore, references to Osprey java classes can be shortened
		from, .e.g.::
			
			edu.duke.cs.osprey.subpackage.Class
			
		to::
		
			.subpackage.Class
	

Sphinx Roles
~~~~~~~~~~~~

**:java:ref:`java_reference`**

	where ``java_reference`` is the fully-qualified name to a Java class, method, or field, e.g.::
	
		package.Class
		package.OuterClass$InnerClass
		package.Class#method
		package.Class#field
		package.OuterClass$InnerClass$ReallyInnerClass#field
		
	This role will create a clickable link to the RST documentation for the referenced Java class,
	method, field, etc.
	
	.. note:: the ````` characters are not single quotes ``'``, but rather grave characters, or backticks.
		
		
**:java:classdoc:`** ``java_class_reference`` **`**

**:java:methoddoc:`** ``java_method_reference`` **`**

**:java:fielddoc:`** ``java_field_reference`` **`**

	where ``java_class_reference`` is any reference allowed by **.. javaclass::**, and
	``java_method_reference`` and ``java_field_reference`` refer to a Java class method or
	field respectively using the ``#`` notation described by **:java:ref:``**

	This role will copy the javadoc comment for the referenced class, method, or field
	into the RST documentation.
	
	
Python Docstring field extensions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`Builder`_ classes in Java are a graceful way to handle class constructors that have
many optional arguments, even though the Java language offers no explicit support for
optional method arguments. Since the Python language *does* support explicit optional
function arguments, Osprey's Python module provides custom builder functions that wrap
these Java builder classes and make Osprey's Python API seem more 'Pythonic'.

.. _Builder: https://en.wikipedia.org/wiki/Builder_pattern#Java

Osprey adds new docstring fields to help translate the javadoc comments for these builder
classes into the Python builder functions documentation.

**:default** ``argname`` **:** ``value``

	This docstring field causes the documentation to display ``value`` as the default value for
	the function or method argument named ``argname``, instead of the default value in the
	Python code itself.
	
	This extension is used internally by the **:builder_option:** docstring field,
	but is also useful on its own.
	
**:builder_option** ``argname`` ``java_field_ref`` **:**

	This extension generates documentation for the builder function argument named ``argname``
	that represents the field referred to by ``java_field_ref`` in a Java builder class.
	The documentation will show the javadoc comment for the field (if any exists) and the type
	of the field. If a value is assigned in the field initializer, then the default value
	will be shown in the Python documentation as well.


**:builder_return** ``java_class_ref`` **:**

	This extension automatically creates an **:rtype:** docstring field based on
	the ``build()`` method of the Java Builder class referenced by ``java_class_ref``.
	

Javadoc extensions
~~~~~~~~~~~~~~~~~~

Since Osprey's documentation toolchain renders javadoc comments into RST, we can easily
define a few new javadoc tags that invoke RST features that wouldn't otherwise be present
in javadoc-based documentation.

**@note** ``message``

	This javadoc tag causes ``message`` to appear inside an RST ``note`` directive, like so:
	
	.. note:: ``message``
	
**@warn** ``message``

	This javadoc tag causes ``message`` to appear inside an RST ``warning`` directive, like so:
	
	.. warning:: ``message``
	
**@cite** ``KEY`` ``citation``

	This javadoc tag renders a citation using ``KEY`` as a unique key, like so: [KEY]_
	
	.. [KEY] ``citation``
	