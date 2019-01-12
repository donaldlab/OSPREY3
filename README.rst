
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
~~~~~~~
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

on Windows:
-----------

Make sure you're running 64-bit Windows. Osprey is not supported on 32-bit Windows.

 1. Install `Python 2.7 x86-64`_, choose the ``Windows x86-64 MSI Installer`` option, not Python 3+ or 32-bit version.
 2. During Python installation, enable the option to ``Add python.exe to Path``.
 3. Install `Java 8 64-bit`_, choose the ``Windows Offline (64-bit)`` option, not 32-bit version.
 4. After installing Java, add the ``C:\Program Files\Java\jre1.8.0_151\bin`` folder to your ``PATH`` environment variable.
    (`See how to set the PATH Environment Variable`_) Be sure to replace the ``jre1.8.0_151`` part with the actual Java
    installation folder on your computer. Tragically, the Java installer does not do this for you.
 5. Download the `newest Osprey Python release`_ (not the source files) and extract it to your favorite folder.
 6. Run the ``install`` batch script to install Osprey.

.. _Python 2.7 x86-64: https://www.python.org/downloads/release/python-2714/
.. _pip: https://pip.pypa.io/en/stable/
.. _Java 8 64-bit: https://www.java.com/en/download/manual.jsp
.. _See how to set the PATH Environment Variable: https://www.java.com/EN/DOWNLOAD/HELP/PATH.XML
.. _newest Osprey Python release: https://github.com/donaldlab/OSPREY_refactor/releases

on Debian-like Linux:
---------------------

*Including distributions like Ubuntu and Mint*

Make sure you're running 64-bit Linux. Osprey is not supported on 32-bit Linux.

 1. Install prerequisites::

	$ sudo apt-get install python2.7 python-pip openjdk-8-jre

 2. Download the `newest Osprey Python release`_ (not the source files) and extract it to your favorite folder.
 3. Run the install shell script to install Osprey::

 	$ ./install.sh


manually using ``pip``:
-----------------------

The install scripts use ``pip`` internally to install the Python package. If you want to customize
the installation of the python package, you can ingore the install scripts and call ``pip`` directly.
First download the `newest Osprey Python release`_ and extract it to your favorite folder. Then call ``pip``::

	$ pip2 install osprey --no-index --use-wheel --find-link=wheelhouse


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
in the distribution zip or the Python documentation at ``doc/api.osprey.html``.

\
    **TODO:** add links to online docs (eg tutorials, references)


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
