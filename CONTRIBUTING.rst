
Osprey guide for contributors
=============================

Prerequisites
~~~~~~~~~~~~~

To devleop Osprey, you'll need to have the following software packages installed:

 * `Python 2.7`_ (not Python 3+)
 * `pip`_
 * Sphinx (for building documentation. see `Documentation`_ section below)

.. _Python 2.7: https://www.python.org/download/releases/2.7/
.. _pip: https://pip.pypa.io/en/stable/


Clone project from Github
~~~~~~~~~~~~~~~~~~~~~~~~~

The Osprey code repository is hosted at `Github`_.
To develop Osprey, first clone the git repository to a folder of your choosing::

	$ mkdir osprey
	$ cd osprey
	$ git clone git@github.com:donaldlab/OSPREY_refactor.git .

.. _Github: https://github.com/donaldlab/OSPREY_refactor


Setup development environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Osprey uses the `Gradle`_ build system with the `Kotlin DSL`_ to automate certain development tasks like:

* Download Java libraries (e.g. jar files)
* Maintain the Java classpath needed to compile or run Osprey
* Compile code and build the Osprey jar file
* Compile Osprey documentation
* Build the Osprey distribution zip files
* Configure Java IDEs
* Configure python development environment

.. _Gradle: https://gradle.org/
.. _Kotlin DSL: https://blog.gradle.org/kotlin-meets-gradle

You don't need to install Gradle, since it's already included in the Osprey repo.
To begin developing Osprey, just use your favorite Java IDE to import the Gradle project.
Gradle has built-in or plugin support in many popular IDEs, including `Eclipse`_, `IntelliJ IDEA`_,
and `Netbeans`_. The import process in your IDE should walk you through the steps needed
to setup Osprey for compilation and testing.

.. _Eclipse: https://www.eclipse.org/
.. _IntelliJ IDEA: https://www.jetbrains.com/idea/
.. _Netbeans: https://netbeans.org/

Osprey's build scripts use a bunch of Python build tools, so you'll need to install those::

    $ pip install setuptools wheel sphinx guzzle_sphinx_theme javalang


Once your build environment is set up, setup the Osprey Python development environment with the
Gradle task ``pythonDevelop``::

	$ ./gradlew pythonDevelop

This task 'installs' the Osprey Python package in editable mode, so your changes to the Python source files will
be immediately visible to Python scripts. This task will also install the `JPype`_ dependency
needed by Osprey.

To undo development mode (ie, uninstall Osprey's python package), use the Gradle task
``pythonUndevelop``::

	$ ./gradlew pythonUndevelop

This will remove the Osprey package from python. Removing the development verison of Osprey
is necessary if you wish to install the Python package included in Osprey distribution.


Build Osprey distribution
~~~~~~~~~~~~~~~~~~~~~~~~~

To build the Osprey distribution files, simply run the Gradle task ``assemble`` from your
IDE's Gradle UI.

Or, if you're not using an IDE, you can run Gradle from the command line::

	$ ./gradlew assemble

**WARNING:** Using the usual ``build`` task in Gradle will indeed build the distribution,
but it will also run the entire suite of Osprey unit tests. The testing suite is rather
extensive and can take up to an hour to run depending on available hardware. To build Osprey
without running the test suite, use the ``assemble`` task instead.

The distribution files will be saved to ``build/distributions``.


Documentation
~~~~~~~~~~~~~

Osprey is mostly implemented in Java, but exposes a Python API as its user interface.
For the most part, this Python API is a very thin wrapper around the Java classes and is
dynamically generated using `JPype`_. The goal of this arrangement is use the Java language
and the JVM to enable rapid development of high-performance code, but provide end users with
the ease and flexibility of Python scripting. As much as possible, Osprey should look and
feel like a real Python application, even though it's really not.

.. _JPype: http://jpype.readthedocs.io/en/latest/

Most of the documentation for Osprey exists in the form of Javadoc comments, to make it easier
to keep the documentation up-to-date with code changes. However, the documentation the users see
is created with the typical Python documentation toolchain, `Sphinx`_. Following Python conventions,
Osprey's documentation outside of javadoc comments is written in the `ReStructured Text (RST)`_
format and can be found in the ``python/doc`` folder. For the javadoc comments, Osprey contains
a custom Sphinx extension (at ``python/doc/javadoc.py``) that converts javadoc comments into RST
documentation much like the `autodoc extension to Sphinx`_.

.. _Sphinx: http://www.sphinx-doc.org
.. _ReStructured Text (RST): https://en.wikipedia.org/wiki/ReStructuredText
.. _autodoc extension to Sphinx: http://www.sphinx-doc.org/en/stable/ext/autodoc.html

To build the documentation for Osprey, run the Sphinx tool using the Gradle task ``makeDoc``::

	$ ./gradlew makeDoc

For quick edit-compile-test cycles when editing documentation, it's helpful
to run the ``cleanDoc`` task before ``makeDoc`` task. This makes sure all documentation
is refreshed regardless of which RST documents have been recently edited. Sphinx won't know if
you've updated Javadoc comments, for instance. The Gradle task ``remakeDoc`` chains the two
commands automatically::

    $ ./gradlew remakeDoc

**NOTE:** Sphinx can detect problems with the documentation during building.
When this happens, these problems will be reported to the console, usually in red text.
These warning messages usually indicate something is missing or incorrect
in the documentation, and that the underlying problems should be fixed before
the documentation is released.

Documentation is built to the ``build/python/doc`` folder. Open in the ``index.html``
file there in your Browser to view the documentation.


Extensions to Sphinx to read Javadoc comments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Osprey's javadoc extension to Sphinx provides a few directives and roles to allow referring to
Java classes, fields, methods, and javadoc comments from the RST documentation:


Sphinx Directives
-----------------


**.. javaclass:: java_class_reference**
    
    where ``java_class_reference`` is the fully-qualified name of a Java class, e.g.::
    
    	package.Class
    	package.OuterClass$InnerClass
    
    This directive will automatically scan the source code for the specified class
    and show all the public constructors, methods, and fields for the class. Javadoc
    comments will be shown with the constructors, methods, arguments, fields, etc,
    and Java type information will be shown in the documentation where possible.
    
        **NOTE:** When the java reference is prefixed with a ``.``, the package
        ``edu.duke.cs.osprey`` is automatically inferred. Therefore, references
        to Osprey java classes can be shortened from, .e.g.::
        	
        	edu.duke.cs.osprey.subpackage.Class
        	
        to::
        
        	.subpackage.Class


Sphinx Roles
------------

**:java:ref:`java_reference`**
    
    where ``java_reference`` is the fully-qualified name to a Java class, method, or field, e.g.::
    
    	package.Class
    	package.OuterClass$InnerClass
    	package.Class#method
    	package.Class#field
    	package.OuterClass$InnerClass$ReallyInnerClass#field
    
    This role will create a clickable link to the RST documentation for the referenced Java class,
    method, field, etc.
    
        **NOTE:** the ````` characters are not single quotes ``'``, but rather grave
        characters, or backticks.
    

**:java:classdoc:`** ``java_class_reference`` **`**

**:java:methoddoc:`** ``java_method_reference`` **`**

**:java:fielddoc:`** ``java_field_reference`` **`**
    
    where ``java_class_reference`` is any reference allowed by **.. javaclass::**, and
    ``java_method_reference`` and ``java_field_reference`` refer to a Java class method or
    field respectively using the ``#`` notation described by **:java:ref:``**
    
    This role will copy the javadoc comment for the referenced class, method, or field
    into the RST documentation.
	
	
Python Docstring field extensions
---------------------------------

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
	

Extensions to Javadoc enabled by Sphinx
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since Osprey's documentation toolchain renders javadoc comments into RST, we can easily
define a few new javadoc tags that invoke RST features that wouldn't otherwise be present
in javadoc-based documentation.

**@note** ``message``
    
    This javadoc tag causes ``message`` to appear inside an RST ``note`` directive, like so:
    
        **NOTE:** ``message``
	
**@warn** ``message``
    
    This javadoc tag causes ``message`` to appear inside an RST ``warning`` directive, like so:
    
        **WARNING:** ``message``
	
**@cite** ``KEY`` ``citation``

	This javadoc tag renders a citation using ``KEY`` as a unique key, like so: [KEY]_
	
	.. [KEY] ``citation``
	
