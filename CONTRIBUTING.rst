
Osprey guide for contributors
=============================


Clone project from Github
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Osprey code repository is hosted at `Github`_.
To develop Osprey, first clone the git repository to a folder of your choosing::

	$ mkdir osprey
	$ cd osprey
	$ git clone git@github.com:donaldlab/OSPREY_refactor.git .

.. _Github: https://github.com/donaldlab/OSPREY_refactor


Setup Eclipse build path
~~~~~~~~~~~~~~~~~~~~~~~~

Osprey uses the `Jerkar`_ build system. You don't need to install Jerkar, since it's
included in the Osprey repo. Just run this in a shell at the project root::

    $ ./jerkar eclipse#generateFiles

This will download the Java libraries needed by Osprey and generate the ``.classpath``
file for Eclipse. Then Eclipse can build the project and you’ll be able to launch
Osprey from within Eclipse for debugging.

.. _Jerkar: http://project.jerkar.org

\
    **TODO:** Jerkar supports other IDEs too (eg, Netbeans, IntelliJ), but I don't know
    how to use those IDEs. Anyone want to help write docs for those IDEs?


Build Osprey distribution
~~~~~~~~~~~~~~~~~~~~~~~~~

To build an Osprey distribution zip file, simply run in the Osprey project folder::

	$ ./jerkar doDist

The distribution zip file will be saved to ``build/output``.

To build the documentation, you'll need to install `Sphinx`_.
If you want to skip building the documentation, try this command instead::

    $ ./jerkar doDist -makeDocs=false

.. _Sphinx: http://www.sphinx-doc.org


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

.. _ReStructured Text (RST): https://en.wikipedia.org/wiki/ReStructuredText
.. _autodoc extension to Sphinx: http://www.sphinx-doc.org/en/stable/ext/autodoc.html

To build the documentation for Osprey, run the Sphinx tool from the ``doc`` folder::

	$ cd python/doc
	$ make html

For quick edit-compile-test cycles when editing documentation, it's helpful
to run ``make clean`` before ``make html`` which makes sure all documentation
is refreshed regardless of which RST documents have been recently edited. e.g.::

    $ make clean && make html

\
    **WARNING:** Sphinx can detect problems with the documentation during building.
    When this happens,
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
    
        **NOTE:** ``message``
	
**@warn** ``message``
    
    This javadoc tag causes ``message`` to appear inside an RST ``warning`` directive, like so:
    
        **WARNING:** ``message``
	
**@cite** ``KEY`` ``citation``

	This javadoc tag renders a citation using ``KEY`` as a unique key, like so: [KEY]_
	
	.. [KEY] ``citation``
	
