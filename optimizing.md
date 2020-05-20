
# Optimizing Osprey

First of all, thanks for taking an interest in improving Osprey!
We're so delighted you want to help. =D

This document will briefly explain how Osprey works at a high level,
and where the bottlenecks are, and how to navigate and build the code
so you can get started optimizing quickly.


## Part 1. A high-level overview of how Osprey designs molecules

Osprey's purpose in life is to find protein sequences that have properties we're
interested in, like stability or affinity with other molecules. Osprey does this
by doing two basic tasks over and over:

1. Find a conformation
2. Compute the minimized energy of that conformation

In this case, a *conformation* is a molecule, but we've replaced some of the atoms
in it with a different set of atoms. This allows us to describe mutations to a
molecule and discrete flexibility in the atom coordinates in an efficent, compact way.
Osprey finds these conformations using various graph algorithms operating on a
*conformation space*. A conformation space is the set of all mutations we can make
to a molecule, and all the ways we allow the atom coordinates to change within a mutation.

Once Osprey has a conformation, it computes the energy of that conformation according to
a physical forcefield to update its scoring model about that conformation, and the
sequence of which the conformation is a part. To model the possibility that the
conformation we actually want isn't perfectly represented by the atom coordinates
we started with, Osprey allows the atoms to move slightly from their original positions.
Osprey uses a *continuous motion* to move the atoms, for example, by a dihedral angle
rotation around a rotatable bond. Using the forcefield and the continuous motions as
an objective function, Osprey minimizes the conformation to find the energy of the best
nearby atom coordinates.

Once Osprey has enough information about different conformations and their minimized
energies, it can claim that a certain sequence has or doesn't have the properties
that the design project seeks.

For now, 2 is the main performance bottleneck in the Osprey code.
This document will describe enough background to understand how 2
works, and give an overview of the code that implements it.


### Conformations and conformation spaces in more detail

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
and also another design position at residue 17 that only has discrete flexibility:

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

Note that the term *conformation* is slightly overloaded here. The term *conformation*
generally means an entire molecule where sets of atoms have been swapped.
However, In the context of a design position, the term *conformation* means a set of atoms
that can be swapped in.

In this example conformation space, there are a total of 4 * 8 = 32 conformations,
since Residue 42 can be one of 4 different options, Reisdue 17 can be one of 8 different
options, and both of these choices are completely independent.

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


### Conformations and forcefields in more detail

Once we have a conformation, we want to score it to see how biophysically likely it
is to occur in some assumed environment, like the aqueous interior of a cell.

Osprey uses two different forcefields to do this:
 * [AMBER](https://en.wikipedia.org/wiki/AMBER)
 * [EEF1 or Effective Energy Function 1](http://charmm.sunhwanj.com/html/eef1.html)
 
For the AMBER forecfield, we only use the non-bonded terms modeling van der Waals
and electrostatic interactions. EEF1 models the interactions between the molecule
and an implicit bulk solvent (eg water).

The forcefields are essentially functions that map a set of atoms into a real
scalar value called *energy*. Generally speaking, the lower the energy of a conformation,
the more likey we are to observe that conformation occurring in nature.

These particular forcefields work by treating the conformation (which could be a single
molecule by itself, or multiple molecules interacting with each other) as a list of pairs
of atoms, where each pair of atoms is associated with *forcefield parameters*.

For example, in AMBER, every pair of atoms with at least two bonds between them
(or no bonds between them at all) is associated with three values:
`vdwA`, `vdwB`, and `esQ`.

The energy of that atom pair is represented by an equation similar to:
```
E = eqQ/(r^2) + vdwA/(r^12) - vdwB/(r^6)
```
where `r` is the radius (aka distance) between the two atom positions.
The total energy for the conformation is the sum of all the energies of all
the atom pairs.

Therefore, computing the energy of a conformation consists of enumerating all
the atom pairs, computing the atom pair energies, and then computing the sum of all
the atom pair energies.

The conformation space also contains the lists of atom pairs for each interaction
between design positions, and the forcefield parameters for each atom pair.


### Energetic interactions between design positions

Since atoms in the molecule change depending on which conformation we're looking at,
the list of atom pairs and forcefield parameters must change too. Groups of atom pairs
in the conformation space are organized by pairs of design positions, and the conformations
those design positions can have.

For example, `posi1=0,posi2=1` represents all the atom pairs between the atoms in
position `0` and position `1` in a conformation. If the conformation in this case is `[0,0]`,
then we'll lookup the atom pairs for this position interaction in the conformation space
at position `0` conformation `0`, and position `1` conformation `0`.

Positions can also be paired with themselves. So `posi1=0,posi2=0` refers to all the atom
pairs *within* position `0`.

Positions can also be paired with the *static* position, which is all of the atoms outside
of any design position. eg `posi1=0,posi2=static`.

Finally, there's the interaction within the static position,
eg `posi1=static,posi2=static`. This represents the internal energy of all the atoms
not in any design positions.

The `static` position was so named because its atoms aren't changed by any design position.
However, other mechanisms can move the `static` atoms, such as dihedral angles defined
for a whole molecule, or translation and rotation motions for a whole molecule.
So it's not correct to assume that the atom coordinates of the `static` position are
constant.

When we're evaluating a conformation, we'll want the complete set of all position
interactions, including interactions within design positions, interactions between all
pairs of design positions, interactions between all positions and the static position,
and interactions within the static position.

Osprey's design algorithm that try to find conformations to evaluate sometimes only
care about interactions between some of the positions. For example, a design algorithm
in Osprey may want to compute the pairwise energy between positions `0` and `1`
independently of the other interactions. So our energy calculators need to support
compute energies on aritrary subsets of the position interactions.


### Computing minimized energies for conformations 

Computing the energy from the atom coordinates exactly as described by the conformation
space is a crude model of molecular conformational flexibility. Often, small overlaps
between atoms that have just been mutated would lead to extremely high forces between
these atoms (recall the `r^12` terms in the van der Waals equations!) and hence extremely
poor energies. Osprey improves its biophysical modeling by allowing the atom coordinates
to "relax" before recording the conformation energy. That is, the atoms are allowed to
move a little closer to an equilibrium state where the atomic forces are more in balance.
As long as the atoms don't move too far away from their original positions, it's sound
to say the "relaxed" conformation is a better physical model of the origial conformation
before relaxation.

Osprey performs this relaxtion by minimizing the energy function over continuous motions
of the atoms, like dihedral rotations around rotatable bonds. The parameters for these
dihedral rotations are also defined in the conformation space. In general, there are
many different kinds of continuous motions we need to model for molecular flexiblity,
and dihedral angles are one of these. More continuous motions are planned for the future,
but for now, dihedral angle motions are the only motions implemented in the newest
minimizer code because they're the ones we use most often by far.

Our minimization algorithm of choice in Osprey is
[Cyclic Coordinate Descent (CCD)](https://en.wikipedia.org/wiki/Coordinate_descent) due
to its simplicity and speed. All the continuous motions are translated into scalar
[*degrees of freedom*](https://en.wikipedia.org/wiki/Degrees_of_freedom_(statistics)).
For dihedral angle rotations, there's only one degree of freedom:
the rotation angle. But in general, a continuous motion can map to many degrees of freedom.
For example, translation and rotation motions map to 6 degees of freedom.

Each degree of freedom is bounded to lie within an interval. These bounding intervals
are also described in the conformation space.

Since the forcefield can be broken down by position interaction, the line search
subroutine inside CCD can focus only on the atom pairs affected by the degree of freedom
being searched. This single optimization causes Osprey's CCD implementation to consistently
out-perform other gradient-based minimization algorithms we've tried in the past, like
[Newtonian methods](https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization).


## Part 2. The code

Now that we've covered some of the high-level concepts of Osprey, we can start to
explore how these are implemented in the code.


### Overall architecture

Osprey is almost entirely written in Java (see `/src`).
There's a thin layer of Python on top (see `/python/osprey`) acting as a scripting API
for designers and data scientists. Chances are, you won't need to use the Python API
at all.
A small part of performance-critical code has been copied into to C++ (see `/native`)
and optimized for various hardware platforms, including CPUs and GPUs.

Osprey is also a very old project full of lots of legacy code written mostly by students who
are not professional software engineers. This is probably the main reason why Osprey is in
Java: to remain accessible to students who are trying to conduct research projects in
computational biology and bioinformatics, or structural biology, and not necessarily
in software engineering.

Because of that, Osprey has naturally adopted a copy-on-write approach to improving
systems. Rather than do painful refactoring that may distrupt a student's progress
towards a research project, it's generally far easier to write new independent code
to perform a new feature, and convince students to upgrade to the new version after a
research project has finished, rather than during the project.

This copy-on-write philosophy has led to the proliferation of many different implemenations
of similar features in the codebase over time. And due to the needs of academic
reproducibility, it's essentially impossible to delete older implementations that are
no longer actively in use. So Osprey's code is littered with a complete history of
various implementations for concepts like conformation spaces, energy functions, and
minimizers. This can all be extremely confusing to newcomers to the code, so one
of the main purposes of this document is to point you quickly in the direction of
the modern and actively in-use implementaions of Osprey's basic features.


### Where's the active code?

Here's a quick roadmap of the Java classes that implement the most modern versions of
Osprey's basic features. The class names will be given by Java's verbose Fully-Qualified
class name, but the common package `edu.duke.cs.osprey` has been omitted.

These classes are all in the `jdk14` branch of the git repo for reasons that will be
explained later.

 * **Conformation Space** `confspace.compiled.ConfSpace`
   * This class reads the conformation space description from a file that's
     prepared by an entirely separate [GUI application](https://github.com/donaldlab/osprey-gui).
     If you wanted to prepare a new conformation space from scratch, you'll need to run
     the GUI program. But for the purposes of optimization, it's sufficient to use the
     previously-prepared conformation spaces that are built into Osprey's test suite
     and benchmarks.
     
 * **Conformation** `confspace.compiled.AssignedCoords`
   * The atom coordinates and degrees of freedom that embody a conformation.
     
 * **Conformation energy calculator** `energy.compiled.CPUConfEnergyCalculator`
   * This is the Java implementation of the newest energy calculator and minimizer.
     It's the slowest option available, but it's accessible to students,
     and makes a good baseline for benchmarks.
     
 * **"native" energy calculator** `energy.compiled.NativeConfEnergyCalculator`
   * This is the reference implementation of the C++ version of the energy calculator
     and minimizer. It's basically a straight port of the Java code into C++ and it's
     totally unoptimized. It also makes a good baseline for benchmarking.
   * It was created originally as a stepping stone for the CUDA implementation,
     so it's written in a very C-style subset of C++, and uses exactly zero libraries.
   * Also, C++ isn't my main language these days, so I'm a little rusty.
     If you have feedback on how to write better C++ code, please share.
     I'd love to learn more!
     
 * **CUDA energy calculator** `energy.compiled.CudaConfEnergyCalculator`
   * This is the CUDA version of the energy calculator. It's been optimized quite a bit
     and is currently the fastest option available by far.
   * The approach to parallelism here is to run one minimization per SM (symmetric
     multiprocessor) on the GPU and rely heavily on the `syncthreads()` mechanism to
     allow us to mix parallel and serial code and handle temporal dependencies.
   * Newer CUDA versions have introduced cooperative kernels that might give us
     better options for parallelism, but I haven't had time to look into that yet.
     It's not clear yet if the multi-block synchronization mechanisms are fast enough
     for our purposes. We're never huring for more conformations to minimize concurrently,
     so one SM per minimization has worked really well so far.
   * We saturate the compute capability of the GPUs by submitting many conformation
     minimizations in parallel, so all the SMs are kept busy.
   * NVidia's SIMT approach works well on the energy function, since computing atom pair
     energies is a bunch of sequential reads and math. But SIMT doesn't work well on the
     CCD part, since CCD is essentially a serial algorithm.
     
 * **Intel energy calculator** `energy.compiled.IntelConfEnergyCalculator`
   * For now, this is nearly a verbatim copy of the "native" conf energy calculator code,
     but I've updated cmake to use `icc` instead of the default `gcc`.
   * I've just started learning how to use `icc`, so I haven't made much progress
     vectorizing the code yet.
   * **If you're looking for a good place to start optimizing, this is it!**

 * **Correctness/Regession tests** `energy.compiled.TestNativeConfEnergyCalculator` (in `/test` not `/src`)
   * Run these tests to make sure the computed energies are correct.
   * Your IDE should give you convenient options to run the tests
     (I recommend [IntelliJ IDEA](https://www.jetbrains.com/idea/)).
     
 * **Benchmarks** `BenchmarkEnergies` (in `/test` not `/src`)
   * Run the `main()` method
   * You'll want to comment/uncomment various sections to disable/enable certain code paths.
   * Use the `benchmarkMinimize()` code path to run the benchmarks for the minimizer.
   * If you want, use the `nativeLab()` code path to isolate certain functions for testing.
   * The benchmarks refer to a `classic` implemenation of Osprey's basic features.
     This is an older way of perparing conformation spaces that we're going to start
     (hopefully) phasing out when the new conformation space tools are more
     production-ready.
   
 * **Thread Pool** `parallelism.ThreadPoolTaskExecutor`
   * Osprey's main tool for task-level parallelism and distributing work across
     parallel hardware.
   * In a nutshell, a task is shipped off to a worker thread where it gets processed.
     Then the task's result is put on a queue where a single listener thread drains the queue
     and sends the result to a callback function specified by the submitter of the task.
   * The benchmarks don't actually use this mechanism now that I think about it.
     They use something simpler that makes repeating and timing the workload easier.

 * **Parallelism** `parallelism.Parallelism`
   * A simpler helper class that describes options for parallelism, like
     the desired size of the thread pool.
   * I don't think the minimizer tests or the benchmarks actually use this,
     but all the design algorithms themselves do.
   * The `streamsPerGpu` option isn't used anymore in the latest CUDA code.
     Stream management is handled internally by the CUDA code now based on
     hardware specs queried at runtime.


### How does Osprey handle high-level concurrency and parallelism?

At a high level, Osprey uses task parallelism to distribute work among CPU cores
and GPU SMs. The `ThreadPoolTaskExecutor` is the main class that handles distributing
parallel workloads using a thread pool.

Most inputs to these parallel tasks are constant objects like the `ConfSpace` instances.
These objects are designed to be created once and then remain constant
throughout the rest of the program, serving as lookup tables for later computations.
Parallel tasks read these objects without sychronizing since there's no possiblity of
a race.

For varying objects like conformations, or lists of position interations, these
objects are created as part of a task definition and moved across threads with the task,
so no two threads can access it at once (outside of the `TaskExecutor` system).
Other varying objects are created wholly whithin the task as it executes, and are
essentially thread-local.


### How does the Java code interact with the C++ code?

Tragically, the object-oriented programming paradigm fails miserably when you really want
your code to go fast. So some of Osprey's performance-critical functions are implemented
in C++ in more of a data-oriented programming style.

Osprey has a few Java classes that translate the data needed for the C++ implementations
into contiguous blocks of memory suitable for C++ code to consume. Once the data buffers
are ready, the Java code sends the pointers to the C++ functions (defined via `extern "C"`)
over one of Java's FFI (foreign function interface) mechanisms called
[JNA (Java native access)](https://github.com/java-native-access/jna).
JNA loads the `.so` libraries at runtime into the JVM (Java virutal machine) process and
binds the native functions to Java function stubs. It also handles type marshalling
for primitive types and pointers across the Java/C++ boundary.

The `NativeConfEnergyCalculator`, `CUDAConfEnergyCalculator`, and `IntelConfEnergyCalculator`
classes are all the newest versions of the Java/C++ bridge class idea.
They're based on a [Foreign Memory Access API](https://openjdk.java.net/jeps/370)
that's a new feature of JDK14 (Java development kit, ie the compiler and runtime for Java).
These classes handle all the memory layouts and alignments
for the constant data that's needed by the C++ code. They use a crude struct definition
API I wrote on the Java side to help read/write the data from/into the raw memory buffers.


## Part 3. Ok great, how do I compile?

First clone the git repo. Then checkout the `jdk14` branch.
The newest code isn't in the `master` branch yet, because not all
the Osprey devs are using JDK14 (Java development kit, ie the compiler
and runtime for Java) yet.

### Java side

The java build process is complicated. It's all managed by a build system called
[Gradle](https://gradle.org/) though. If you're using the
[IntelliJ IDE](https://www.jetbrains.com/idea/) (which I highly recommend),
IntelliJ has really good integration with Gradle.
Just open the root folder of the osprey project with IntelliJ to import everyting.
Then tell IntelliJ to run the `main()` methods you want (little green arrow in the gutter),
make sure it delegates the build process to Gradle in the IntelliJ settings, and
IntelliJ will handle all the details for you.

Like I mentioned above, the newest code requires JDK14, which isn't the default
Java for many systems yet. The easiest place to get it is
[AdoptOpenJDK](https://adoptopenjdk.net/). Make sure IntelliJ is configured to use it
for the Osprey project.

There are ways to build and run on the command line too, if you want to, for remote
profiling/debugging tools etc. It's hard to get gradle to run a specific `main()`
method directly on the command line though. IntelliJ apparently implemented some special
Gradle integration to get that working. But once you've run a `main()` method in the
IDE, I've configured Gradle to spit out the entire command line invocation to stdout.
Just copy that into your shell or analysis tools.

Fair warning though: The Java class path needed for Osprey is ridiculously long,
so the command line invocation is also ridiculously long.

For more information about the build system, consult `/contributing.rst`.


### C++ side

The C++ side is much easier (so far). Each `ConfEnergyCalculator` implementation
is a separate subfolder inside the `/native` folder. Just open the project folder
with your favorite C++ development environment and run `cmake`.

The main `cmake` task (named something like `XXXConfEcalc`) will build the shared
library for the project. There's a second task (named something like
`XXXConfEcalc_CopyLibs`) that will build the shared library and also copy it to
the place where the Java side expects to find it. Use that task if you want to
run the Java code using your newly-built library.


## Part 4. Getting more info

If you have any questions, want some clarification, or need any more information,
don't hesitate to send a message to `osprey@cs.duke.edu`.
It's a mailing list monitored by several people who will be able
to answer your questions.
