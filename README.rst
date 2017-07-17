
OSPREY
======

*Open-Source Protein REdesign for You*

OSPREY is developed and maintained by the `Donald Lab`_
in the `Department of Computer Science`_
at `Duke University`_.

.. _Donald Lab: http://www.cs.duke.edu/donaldlab/home.php
.. _Department of Computer Science: http://www.cs.duke.edu
.. _Duke University: https://www.duke.edu/


Development Version
-------------------

This is the development version of Osprey. For stable versions, please visit:
http://www.cs.duke.edu/donaldlab/osprey.php


License
-------

Osprey is released under a permissive GPL-style license. See the
`full Osprey license`_ for more information. 

.. _full Osprey license: LICENSE.txt


Prerequisites
-------------

 * `Python 2.7`_ (not Python 3+)
 * `pip`_

.. _Python 2.7: https://www.python.org/download/releases/2.7/
.. _pip: https://pip.pypa.io/en/stable/


Installation
------------

Download the Osprey Python distribution zip, extract it to your favorite folder,
then use the provided shell script to install::

    $ ./install.sh

if using Linux/Mac, or::

	$ install.bat

if using Windows

The install scripts use ``pip`` internally to install the Python package. If you want to customize
the installation of the python package, you can ingore the install scripts and call ``pip`` directly::

	$ pip2 install osprey --no-index --use-wheel --find-link=wheelhouse


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
	
For more information about Python scripting with Osprey, see the tutorial at ``doc/tutorial.html``
in the distribution zip or the Python documentation at ``doc/api.osprey.html``.

\
    **TODO:** add links to online docs (eg tutorials, references)


using the command-line interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Python interface to Osprey represents a significant improvement in the user interface over the
older command-line interface, and new Osprey projects should consider using the Python interface
rather than the command-line interface.

However, for backwards compatibility, the command-line interface is still provided, although
it may not receive feature updates in the future. Or it may be removed entirely.

To access the command-line interface, download the Osprey Java distribution.
Extract it to your favorite folder, then enter the following command into a shell::

    $ cd bin
    $ ./osprey [commands]

where ``[commands]`` are the Osprey commands you want to run. You can run Osprey without
``[commands]`` and Osprey will print a list of the available commands.

To show the version of your Osprey installation, try::

    $ ./osprey version

To run a GMEC-based protein design, try::

    $ ./osprey FindGMEC /path/to/config1 /path/to/config2 ...

To show GPU informatino, try::

    $ ./osprey GPUInfo

The ``GpuInfo`` command prints info about available GPUs in the system, and which
ones Osprey can use.


Contributing
------------

Osprey is open-source software and contributions are welcome.

See the `guide for contributors`_ to see how to compile and package Osprey.

.. _guide for contributors: CONTRIBUTING.rst
