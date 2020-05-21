
# Overview of Task 2: Evaluate the conformation

*This document focuses on just the minimization bottlenecks in Osprey.
You'll want to read [the general overview](optimizing.md) first.*


## Part 1. Theory

### Forcefields in more detail

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

The class names will be given by Java's verbose Fully-Qualified
class name, but the common package `edu.duke.cs.osprey` has been omitted.
All of these classes are in the `/src` folder.

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


### Testing and Benchmarking

All of these classes are in the `/test` folder.

 * **Correctness/Regession tests** `energy.compiled.TestNativeConfEnergyCalculator`
   * Run these tests to make sure the computed energies are correct.
   * Your IDE should give you convenient options to run the tests
     (I recommend [IntelliJ IDEA](https://www.jetbrains.com/idea/)).
     
 * **Benchmarks** `BenchmarkEnergies`
   * Run the `main()` method
   * You'll want to comment/uncomment various sections to disable/enable certain code paths.
   * Use the `benchmarkMinimize()` code path to run the benchmarks for the minimizer.
   * If you want, use the `nativeLab()` code path to isolate certain functions for testing.
   * The benchmarks refer to a `classic` implemenation of Osprey's basic features.
     This is an older way of perparing conformation spaces that we're going to start
     (hopefully) phasing out when the new conformation space tools are more
     production-ready.


### Interaction between Java code and C++ code

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


### Compiling the C++ code

The C++ side is much easier to compile than the Java side (so far).
Each `ConfEnergyCalculator` implementation
is a separate subfolder inside the `/native` folder. Just open the project folder
with your favorite C++ development environment and run `cmake`.

The main `cmake` task (named something like `XXXConfEcalc`) will build the shared
library for the project. There's a second task (named something like
`XXXConfEcalc_CopyLibs`) that will build the shared library and also copy it to
the place where the Java side expects to find it. Use that task if you want to
run the Java code using your newly-built library.
