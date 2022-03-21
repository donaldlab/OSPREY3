+++
title = "Documentation"
weight = 2
+++


TODO: explain how to edit the documentation



### Prerequitites for building the documentation

[hugo](https://gohugo.io)
* [install release](https://github.com/gohugoio/hugo/releases)
  (get the "extended" version, needed by pydoc-markdown)

[pydoc-markdown](https://github.com/NiklasRosenstein/pydoc-markdown):
* [installation](https://pydoc-markdown.readthedocs.io/en/latest/docs/getting-started/)



TODO: update this for the new doc system

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
