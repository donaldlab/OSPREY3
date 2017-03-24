
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
	
For more information about Python scripting with Osprey, see the tutorial at ``doc/tutorial.html``
in the distribution zip or the Python documentation at ``doc/api.osprey.html``.

.. note:: TODO: add links to online docs (eg tutorials, references)


using the command-line interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Python interface to Osprey represents a significant improvement in the user interface over the
older command-line interface, and new Osprey projects should consider using the Python interface
rather than the command-line interface.

However, for backwards compatibility, the command-line interface is still provided, although
it may not receive feature updates in the future.

Run Osprey at the command line by navigating to the folder where you extracted the Osprey distribution. Then enter the following command into a shell::

    $ java -jar osprey/osprey-<version>.jar [commands]
    
where ``<version>`` is the version of Osprey and ``[commands]`` are the Osprey commands
you want to run. You can run Osprey without ``[commands]`` and Osprey will print a list
of the available commands.
    
To show the version of your Osprey installation, try::

    $ java -jar osprey.jar version

To run a GMEC-based protein design, try::

    $ java -jar osprey.jar FindGMEC /path/to/config1 /path/to/config2 ...
    

.. note:: To use GPU acceleration, you'll need to need to tell Java where to find the operating
	system-specific GPU libraries. The Python interface does this automatically, but on the command line,
	you'll need to supply the additional JVM argument::
	
		-Djava.library.path=osprey/natives
		
	For example, to run the ``GpuInfo`` command::
	
		$ java -Djava.library.path=osprey/natives -jar osprey/osprey-3.0.jar GpuInfo
		
	The ``GpuInfo`` command prints info about available GPUs in the system, and which
	ones Osprey can use.
	
	If the `java.library.path` is not correctly set, Osprey will think the
	native GPU libraries (e.g. CUDA, OpenCL) are not installed.


Contributing
------------

Osprey is open-source software and contributions are welcome.

See the `guide for contributors`_ to see how to compile and package Osprey.

.. _guide for contributors: CONTRIBUTING.rst
