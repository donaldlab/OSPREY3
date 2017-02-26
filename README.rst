
OSPREY
======

*Open-Source Protein REdesign for You*

OSPREY is developed and maintained by the `Donald Lab <http://www.cs.duke.edu/donaldlab/home.php>`_
in the `Department of Computer Science <http://www.cs.duke.edu/>`_
at `Duke University <https://www.duke.edu/>`_.

Osprey is released under a permissive GPL-style license. See the
`full Osprey license <http://www.cs.duke.edu/donaldlab/software/osprey/osprey.2.2/license.pdf>`_
for more information. 


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
	
For more information about Python scripting with Osprey, see the `documentation <docs>`_.

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
the `command-line interface documentation <docs>`_.

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

Osprey uses the `Jerkar <http://project.jerkar.org/>`_ build system.
First `install Jerkar if needed <http://project.jerkar.org/documentation/latest/getting_started.html>`_,
then run this in a shell at the project root::

    $ jerkar eclipse#generateFiles

This will download the Java libraries needed by Osprey and generate the ``.classpath``
file for Eclipse. Then Eclipse can build the project and you’ll be able to launch Osprey from
within Eclipse for debugging.

Jerkar also provides an `Eclipse plugin <https://github.com/jerkar/eclipsePlugin4Jerkar>`_ to make running Jerkar easier.


Build Osprey distribution
~~~~~~~~~~~~~~~~~~~~~~~~~

To build an Osprey distribution zip file, simply run in the Osprey project folder::

	$ jerkar doDist

The distribution zip file will be saved to ``build/output``.
