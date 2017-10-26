
Tutorial
========

Hello there!

If you're wondering how to write your first Python script to do protein redesign
using Osprey, you've come to the right place. If you haven't installed Osprey yet,
you'll have to do that first. Then come back here.

So, without further ado, here's our first Python script.

.. code-block:: python
	:linenos:

	import osprey

	osprey.start()

	# define a strand
	strand = osprey.Strand('1CC8.ss.pdb')
	strand.flexibility['A2'].setLibraryRotamers('ALA', 'GLY')
	strand.flexibility['A3'].setLibraryRotamers(osprey.WILD_TYPE, 'VAL')
	strand.flexibility['A4'].setLibraryRotamers(osprey.WILD_TYPE)

	# make the conf space
	confSpace = osprey.ConfSpace(strand)

	# choose a forcefield
	ffparams = osprey.ForcefieldParams()

	# how should we compute energies of molecules?
	ecalc = osprey.EnergyCalculator(confSpace, ffparams)

	# how should we define energies of conformations?
	confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)

	# how should confs be ordered and searched?
	emat = osprey.EnergyMatrix(confEcalc)
	astar = osprey.AStarMPLP(emat, confSpace)

	# find the best sequence and rotamers
	gmec = osprey.GMECFinder(astar, confEcalc).find()

This is an example script in the Osprey distribution.

You can find it at ``examples/python.GMEC/findGMEC.py``
It's a bit much to take in all at once though, so let's break it down line-by-line.

.. code-block:: python
	:lineno-start: 1

	import osprey


This line imports the ``osprey`` module into the Pyton script. It's pretty standard for python scripts.

.. code-block:: python
	:lineno-start: 3

	osprey.start()

This line is special to Osprey though. Under-the-hood, Osprey is doing most of the
computational work using Java, rather than Python. The :py:func:`osprey.start`
starts the Java Virtual Machine (JVM) and gets it ready to accept commands
from the Python script.

.. code-block:: python
	:lineno-start: 5

	# define a strand
	strand = osprey.Strand('1CC8.ss.pdb')
	strand.flexibility['A2'].setLibraryRotamers('ALA', 'GLY')
	strand.flexibility['A3'].setLibraryRotamers(osprey.WILD_TYPE, 'VAL')
	strand.flexibility['A4'].setLibraryRotamers(osprey.WILD_TYPE)

These lines are where much of the configuration for your design will take place, since
this is where we define what flexibility is allowed in the design.

An :py:func:`osprey.Strand` is molecule fragment that we have annotated with information
about sequence mutations and conformational flexibility. In this case, we've defined
a strand based on the molecule in PDB file ``1CC8.ss.pdb``.

For residues 2, 3, and 4 in chain A, we've called the :java:ref:`.confspace.Strand$ResidueFlex#setLibraryRotamers`
function to define flexibility information. The *Template Library* in Osprey contains rotamers for
each amino acid, so calling ``setLibraryRotamers`` sets the flexibility for that residue to
the rotamers from the specified amino acids.

For residue ``A2`` specifically, we're forcing the sequence to mutate to either Alanine or Glycine.

For residue ``A3``, we're allowing a mutation to Valine, but the sequence can also stay at whatever
amino acid is present in the PDB file by using the magic constant, :py:const:`osprey.WILD_TYPE`.

For residue ``A4``, we're not allowing any mutations this time, but, by using :py:const:`osprey.WILD_TYPE`
for the amino acid, we're still allowing flexibilty at this design position. Only the rotamers
defined for the wild-type amino acid will be used at this position though.

The rest of the residues have no flexibilty specified, so they will remain completely fixed
in the conformation specified by the PDB file throughout all of Osprey's analyses.

.. code-block:: python
	:lineno-start: 11

	# make the conf space
	confSpace = osprey.ConfSpace(strand)

Here is where we make our *Conformation Space* for the design. The :py:func:`osprey.ConfSpace` is the
object where we collect all of the information about design flexibility. In more complicated designs,
it can hold information about multiple strands, but for now, we have just the one strand.
This information will remain essentially constant for the rest of the script. Other parts of
Osprey will refer to the conf space to see what flexibility is defined for this design.

.. code-block:: python
	:lineno-start: 14

	# choose a forcefield
	ffparams = osprey.ForcefieldParams()

The *Forcefield* and its parameters tells Osprey how to calculate the energy for a molecule in a
specific conformation. For now, we'll just use the default forcefield by calling
:py:func:`osprey.ForcefieldParams` without any arguments.


.. code-block:: python
	:lineno-start: 17
	
	# how should we compute energies of molecules?
	ecalc = osprey.EnergyCalculator(confSpace, ffparams)

This line sets up the energy calculator for Osprey. The job of energy calculator is to calculate the
physical energy of a molecule based on the forcefield parameters, and minimize over the continuous degrees
of freedom of the molecule as needed. Much of Ospey's run time is spent inside the energy calculator, so
lots of work has been done to keep it running as fast as possible. Osprey's energy calculator supports
many different kinds of parallel hardware (including GPUs), but we'll just use the default single-threaded
CPU calculator for this toy design. For larger designs, you'll definitely want to explore the options in
the :py:func:`osprey.EnergyCalculator` to reduce the runtimes of your designs.

.. code-block:: python
	:lineno-start: 20
	
	# how should we define energies of conformations?
	confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)

This line combines the energy calculator with information about your design's conformation space.
The :py:func:`osprey.ConfEnergyCalculator` allows different options for how conformation energies
are created from interactions between residues, but again, we'll just stick with the defaults this time.

.. code-block:: python
	:lineno-start: 23
	
	# how should confs be ordered and searched?
	emat = osprey.EnergyMatrix(confEcalc)
	astar = osprey.AStarMPLP(emat, confSpace)

In the most abstract and simplest sense, Osprey computes designs by executing two steps:

	1.	Define a sort order for all conformations in the conformation space,
		and then start enumerating conformations in that order.

	2.	Stop enumerating conformations when we can prove that we've already
		enumerated the one with the lowest energy.

Many of the sophisticated algorithms implemented in Osprey are variations on this simple theme.

The two lines of code above define how Osprey should do step 1. Since there can be an astronomical
number of conformations in the conformation space, explicitly sorting a list of all them
is generally far too expensive to do in practice. Instead, Osprey uses `A* search`_ to find
the first few conformations in the sort order very quickly.

Osprey's A* implementation uses a matrix of energies between pairs of residue conformations
to define the sort order for the conformation space. The :py:func:`osprey.EnergyMatrix` function
is used to compute the energy matrix using the specified forcefield parameters. Then we call
the :py:func:`osprey.AStarMPLP` function to create the object that performs the A* search
on the conformation space.

.. _A* search: https://en.wikipedia.org/wiki/A*_search_algorithm

.. tip:: add the ``cacheFile='path'`` argument to :py:func:`osprey.EnergyMatrix` to reuse
	the energy matrix between runs of your design. If the energy matrix takes a long time to
	compute, this can save you a lot of time.
	
To perform step 2, Osprey uses the conformation energy calculator on each conformation returned
by the A* search. Different algorithms in Osprey use the resulting energies with various criteria
to choose when to stop enumerating conformations, but this simple discrete design will just stop
after the first conformation is reached.

.. code-block:: python
	:lineno-start: 27

	# find the best sequence and rotamers
	gmec = osprey.GMECFinder(confSpace, astar, ecalc).find()

Finally, this line of code runs the design and computes the best sequence in the conformation
space: the *Global Minimum Energy Conformation*, or GMEC. This is the part of Osprey that
actually computes steps 1 and 2, whereas before we just defined how they should be computed.
This function can take a long time to run, depending on how large the design is, so check the
script's output to get progress information.

To configure parameters for the GMEC search, try the API function :py:func:`osprey.GMECFinder`.

For our simple example though, the whole script should take only a few seconds to complete.
When it's done, you should be greeted something very similar to the following output::

	OSPREY 3.0
	read PDB file from file: 1CC8.ss.pdb
	Calculating energy matrix with 133 entries...
	Searching for min score conformation...
		(among 154.0 possibilities)
	Found min score conformation in 8.9 ms
	Computing energy of min score conf...
	Found GMEC!
		Residue Conf Ids       1   3   4
		Residue types        GLY GLU ILE
		Rotamer numbers        L  L3  L4
		Energy               -30.705504
		Score                -30.705504 (gap: 0.000000)

This simple example only shows off a few of Osprey's features. If you'd like to learn more,
browse the :ref:`api_reference`.

