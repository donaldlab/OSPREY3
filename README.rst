
OSPREY
======

*Open-Source Protein REdesign for You*

OSPREY is developed and maintained by the `Donald Lab`_
in the `Department of Computer Science`_
at `Duke University`_.

.. _Donald Lab: http://www.cs.duke.edu/donaldlab/home.php
.. _Department of Computer Science: http://www.cs.duke.edu
.. _Duke University: https://www.duke.edu/

For an introduction to OSPREY 3.0 and its new features please read this paper: 

*Journal of Computational Chemistry* 2018; 39(30): 2494-2507 `Cover article`_.

.. _Cover article: http://www.cs.duke.edu/brd/papers/jcc18-osprey3point0/cover-jcc.25043.pdf

Available here:

`Journal of Computational Chemistry`_.

`Cover Image \(Osprey)`_ 

`PDF of paper`_

.. _Journal of Computational Chemistry: https://onlinelibrary.wiley.com/doi/10.1002/jcc.25522
.. _Cover Image (Osprey): http://www.cs.duke.edu/brd/papers/jcc18-osprey3point0/cover-jcc.25043.pdf
.. _PDF of paper: http://www.cs.duke.edu/brd/papers/jcc18-osprey3point0/jcc18-osprey-donald.pdf



Citation requirements
~~~~~~~~~~~~~~~~~~~~~
We require everyone who publishes or presents results from OSPREY to please mention the name "OSPREY," and to cite our papers as described in CITING_OSPREY.txt (especially our new paper introducing OSPREY 3.0). 


License
~~~~~~~

`GPLv2`_

Copyright (C) 2017 Duke University

This program is free software; you can redistribute it and/or
modify it under the terms of the `GNU General Public License version 2`_
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The full text of the GPLv2 is included in the accompanying `LICENSE.txt`_

.. _GPLv2: https://www.gnu.org/licenses/gpl-2.0.html
.. _GNU General Public License version 2: https://www.gnu.org/licenses/gpl-2.0.html
.. _LICENSE.txt: LICENSE.txt


Installation
~~~~~~~~~~~~

on OSX:
-------

 #. Install `Python`_.
     * Both v2 and v3 are supported.
     * Install the `64-bit` verison, not the `32-bit` version.
 #. Download the `newest Osprey Python release`_.
     * Download either `osprey-osx-python2` or `osprey-osx-python3` to match your Python version.
     * Extract the acrhive to your favorite folder.
 #. Run the ``install.sh`` script to install Osprey.


on Windows:
-----------

Make sure you're running 64-bit Windows. Osprey is not supported on 32-bit Windows.

 #. Install `Python`_.
     * Both v2 and v3 are supported.
     * Install the `x86-64` verison, not the `x86` version.
     * During Python installation, enable the option to ``Add python.exe to Path``.
 #. Download the `newest Osprey Python release`_.
     * Download either `osprey-win-python2` or `osprey-win-python3` to match your Python version.
     * Extract the acrhive to your favorite folder.
 #. Run the ``install`` batch script to install Osprey.


on Linux:
---------------------

*Including distributions like Ubuntu and Mint*

Make sure you're running 64-bit Linux. Osprey is not supported on 32-bit Linux.

 #. Install `Python`_.
     * For Debian-like distributions, including Ubuntu and Mint, run::

       $ sudo apt install python3 python3-pip

       or::

       $ sudo apt install python2.7 python-pip

       depending on which version of Python you want to use.

 #. Download the `newest Osprey Python release`_.
     * Download either `osprey-linux-python2` or `osprey-linux-python3` to match your Python version.
     * Extract the acrhive to your favorite folder.
 #. Run the ``install`` script to install Osprey.

 	$ ./install.sh


.. _Python: https://www.python.org/downloads/
.. _newest Osprey Python release: https://github.com/donaldlab/OSPREY_refactor/releases


To check your installation:
---------------------------

From the command line, run:

.. code:: pycon

	python3
	>>> import osprey
	>>> osprey.start()

If successful, should should be greeted with a message something like the following::

	OSPREY 3.2-beta1-dev, Python 3.6.9, Java 11.0.6, Linux-4.15.0-91-generic-x86_64-with-LinuxMint-19-tara
	Using up to 1024 MiB heap memory: 128 MiB for garbage, 896 MiB for storage


Upgrading from an older version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you're upgrading from an older installation of Osprey, just run the install script. There's no need
to explicitly uninstall the older version.


Uninstallation
~~~~~~~~~~~~~~

To uninstall Osprey, use the provided shell script.

on Windows::

	> uninstall.bat

on Linux or Mac::

	$ ./uninstall.sh


Running Osprey
~~~~~~~~~~~~~~

using Python scripts
--------------------

Python scripting is the preferred way of using Osprey due to its simplicity and flexibilty.
To run Osprey from a Python script:

.. code:: python

	import osprey
	osprey.start()
	
	# run osprey commands, e.g.
	osprey.printGpuInfo()
	
For more information about Python scripting with Osprey, see the tutorial at ``doc/tutorial.html``
(in the downloaded zip file) or the Python documentation at ``doc/api.osprey.html``.

Many Osprey features are explained in example scripts
which can be found in the downloaded zip file at ``examples/python.*/*.py``.

A comprehensive manual for Osprey has yet to be written,
but these example scripts can help you get started with common design tasks.


using the command-line interface
--------------------------------

The Python interface to Osprey represents a significant improvement in the user interface over the
older command-line interface, and new Osprey projects should consider using the Python interface
rather than the command-line interface.

However, for backwards compatibility, the command-line interface is still provided, although
it may not receive feature updates in the future. It may eventually be removed from Osprey.

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
~~~~~~~~~~~~

Osprey is open-source software and contributions are welcome.

See the `guide for contributors`_ to see how to compile and package Osprey.

.. _guide for contributors: CONTRIBUTING.rst
