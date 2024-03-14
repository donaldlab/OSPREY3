+++
title = "Conformations and Conformation Spaces"
menuTitle = "Conformation Spaces"
weight = 1
+++


A *conformation space* is the set of all mutations we can make to a molecule,
and all the ways we allow the atom coordinates to change within a mutation.

A conformation space describes its mutations to a molecule by dividing up the atoms
into non-overlapping groups called *design positions*. A design position is a small
subset of atoms in the molecule that we can replace with other sets of atoms.

In a protein, a design position might include all the atoms at, say, residue 42.
That design position might also allow residue 42 to mutate from its wild-type Alanine
amino acid to a Valine amino acid. The design position might also include three
different sets of atom coordinates for the Valine amino acid, representing three
different conformational states for the residue.

Here's a description of a conformation space that includes the design
position we described above that contains the mutation from Alanine to Valine,
and also another design position at residue 17 that only has discrete flexibility
for the wild-type amino acid, Tyrosine:

* Residue 42
    * Alanine (wild-type)
        * conformation 1
    * Valine (mutation)
        * conformation 1
        * conformation 2
        * conformation 3
* Residue 17
    * Tyrosine (wild-type)
        * conformation 1
        * conformation 2
        * ...
        * conformation 8

{{% notice note %}}
The term *conformation* is slightly overloaded here. OSPREY uses several different
senses of the word conformation, depending on the context. A conformation can be:

1: An entire molecule where sets of atoms have been swapped out with others or moved from their original positions.

2: The mathematical representation of such a molecule in a conformation space.

3: A molecular fragment that can be exchanged with the existing atoms at a design position in a conformation space.
In the specific case of protein molecules, these are sometimes called *residue conformations*.
{{% /notice %}}

In this example conformation space, there are a total of 4 * 8 = 32 conformations,
since Residue 42 can be one of 4 different options, Reisdue 17 can be one of 8 different
options, and both of these choices are completely independent.

Also in this example conformation, there are a total of 2 * 1 = 2 sequences. The sequences
are (in no particular order):
* Alanine @ Residue 42, Tyrosine @ Residue 17
    * (this sequence has 1*8 = 8 conformations)
* Valine @ Residue 42, Tyrosine @ Residue 17
    * (this sequence has 3*8 = 24 conformations)

To describe these conformations compactly, we index everything in the conformation space
with integers. So the example conformation space becomes:

* `posi=0` Residue 42
    * Alanine (wild-type)
        * `confi=0` conformation 1
    * Valine (mutation)
        * `confi=1` conformation 1
        * `confi=2` conformation 2
        * `confi=3` conformation 3
* `posi=1` Residue 17
    * Tyrosine (wild-type)
        * `confi=0` conformation 1
        * `confi=1` conformation 2
        * ...
        * `confi=7` conformation 8

Now we can describe an entire conformation with a simple list of integers. For example,
the conformation `[0,0]` represents the molecule where Residue 42 is Alanine,
conformation 1 and Residue 17 is Tyrosine, conformation 1.
Another conformation is `[3,7]`.

Sometimes a conformation will only be partially assigned. In that case, the conformation
index at that design position will be `-1`. So the conformation `[-1,-1]` describes only
the atoms in the molecule that lie outside any design positions. The state of the atoms
inside the design positions is undefined. Or rather, unassigned.


## How Conformation Spaces are created

In OSPREY, there are four overall steps to doing a design:

1. Import your molecule(s) of interest
2. Configure a conformation space with your design parameters
3. Compile the conformation space into an optimized format
4. Run the design using the compiled conformation space

The next sections will explain the purpose of each of these steps in a little more detail
and explain which parts of the code are responsible for them.


### 1. Import molecules

Most designers using OSPREY start with `PDB` files containing their molecules of interest
-- usually proteins, but also sometimes small molecules, lipids, carbohydrates, or even nucleic acids.
Tragically, the `PDB` file doesn't natively contain all the information about the molecules OSPREY
needs to do a design. Basic information like the presense of covalent bonds or hydrogen atoms must
be inferred to complete the molecule description for the purposes of design.

OSPREY has tools to make these inferences, but they tend to be very error-prone when fully automated.
To help correct the errors made by automated tools, OSPREY also has manual tools to perform the same tasks.
Typically, the design prep will start with running the automated tools, inspecting the results,
and then using the manual tools to correct any errors.

On the OSPREY Desktop side, the molecule import tools are all collected in the Kotlin `.gui.prep.MoleculePrep` class
under the `Prepare` menu section of the code in the Kotlin GUI code.

On the OSPREY Server side, we have a Python API to handle molecule preparation that does all the same things
as the Desktop software at `/src/main/python/osprey/prep.py`.

Once configured, OSPREY molecules are saved in a custom file format and usually given an `.omol` extension.
The `omol` files save all the information we added that was missing from the `PDB` files. On the Kotlin side,
the code for managing `omol` files is in the Kotlin file `gui/io/OMOL.kt`. On the Java side, the `omol` files
are managed by the `.structure.OMOLIO` class.


### 2. Configure conformation spaces

Once the molecules are saved in the `omol` format, then the designer will build the conformation space(s)
for the design, including configing any mutations and flexibility.

On the OSPREY Desktop side, the main Kotlin class for managing conformation space prep is
`.gui.prep.ConfSpacePrep`. On the OSPREY Server side, the Python API that handles conformation space
preparation is also at `/src/main/python/osprey/prep.py`.

It's worth noting here that, unlike older versions of OSPREY, any kind molecule can be designed, not
just proteins. Which means much of the conceptual parts of a design have been generalized from their
older protein-specific counterparts. For example, conformation space design positions need no longer correspond
exactly to protein residues, since, e.g. small molecules have no residues. The conformation space system
is general enough to design all kinds of molecules now, so defining a design position takes quite a bit more
work than it used to. See the Kotlin `.gui.prep.DesignPosition` class for details on what it takes now to
define a design position. OSPREY now uses a general anchor system to define how to remove, replace, and
align molecular fragments for designs. The Kotlin `.gui.prep.ConfSwitcher` class is responsible for actually
performing the molecular fragment replacements during the prep phase of designing.
In the specific case of proteins, OSPREY has convenient shortcuts using the familiar residue-centric paradigm
in the Kotlin `.gui.prep.Proteins` object.

OSPREY also uses a molecular fragment library to provide a collection of commonly-used fragments in designs,
like amino acids for proteins, along with their rotameric conformations. The default conformation libraries
in OSPREY are saved at `/src/main/resources/edu/duke/cs/osprey/gui/conflib`. In particular, OSPREY's
traditional amino acid rotamer library, based on the work of Lovell at al, has been converted to the new
conformation library format and saved in the `lovell.conflib` file. Conformation libraries are handled by
the Kotlin `.gui.io.ConfLib` class, including reading/writing them from/to files.


### 3. Compile the conformation spaces

The conformation space preparation system in OSPREY is very flexible and very powerful, but it's
not necessarily performant. Since the success of OSPREY designs often depends on how many conformations
can be effectively searched, great care has been put into making sure OSPREY can search conformations
as efficiently as possible. One part of that work is the conformation compilation system.

The main purpose of the conformation space compilation system is two-fold:

1. Do as much pre-computation as possible in the design preparation phase,
   so the design phase itself can be simpler and run more efficently.
2. Provide immediate feedback to users about design preparation errors before
   investing time and other resources running the design.

The conformation space compiler itself lives in the Kotlin `.gui.compiler.ConfSpaceCompiler` class.
The compiler does quite a bit of work to make the conformation spaces ready for design. Perhaps the most
important step is collecting the forcefield parameters for the design molecules. OSPREY currently provides
the EEF1 implicit solvation forcefield and also various flavors of Amber forcefields. The forcefield system
in OSPREY is general enough to allow integrating with (hopefully) any forcefield system, so it should be possible
to integrate with other forcefiled systems in the future, should such a thing ever be desired.

All of the preparation-side forcefield code is in the Kotlin package `.gui.forcefield`.
All of the available OSPREY forcefield implementations are collected in the Kotlin `.gui.forcefiled.Forcefield`
class. Specifically, the Kotlin `.gui.forcefield.ForcefieldParams` interface is the main entry point for OSPREY
to integrate with forcefield systems. You'll find subclasses of that interface for both the Amber and EEF1
forcefield implementations. This interface is primarily used by the conformation space compiler itself
to collect the forcefield parameters for a design, including forcefield parameters for all of the possible
mutations for the molecules being designed.

Once all the forcefield parameters are assembled, the compiler writes out the conformation spaces's conformations and
forcefield parameters into an optimized binary file format and then into a file usually given a `.ccs` extension
(for compiled conformation space). These compiled files can be very large, so they're usually compressed
and written to a `.ccsx` file (the `x` represents the compression). These `.ccsx` files are then given to
OSPREY Server as the definitions of the conformation spaces for a design. While still in memory, the compiled
conformation space is represented by the Kotlin `.gui.compiler.CompiledConfSpace` class. Then the Kotlin code at
`gui/io/CompiledConfSpace.kt` renders the compiled conformation space into an efficient binary encoding.


### 4. Run the design

All of the OSPREY Server design algorithms require conformation spaces as input. For compiled conformation spaces,
read the files using the Java `.confspace.compiled.ConfSpace` class, specifically the `fromBytes` method (of
course, read the file into a byte array first).

From there, the various `assign` methods in the `ConfSpace` class will create conformations, and the
`makeCoords` method will build the atomic coordinates for a conformation using an `AssignedCoords` class instance,
along with any continuous motions that were configured in the conformation space.

Then, just run your favorite design algorithms using the loaded conformation spaces. The details for that
will be specific to the particular design algorithm you're using.
