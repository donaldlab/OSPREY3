
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

**Task 1:** [Find a conformation](optimizing-graphs.md)

**Task 2:** [Evaluate the conformation](optimizing-minimizers.md)

In this case, a *conformation* is a molecule, but we've replaced some of the atoms
in it with a different set of atoms. This allows us to describe mutations to a
molecule and discrete flexibility in the atom coordinates in an efficent, compact way.
Osprey finds these conformations using various graph algorithms operating on a
*conformation space*. A conformation space is the set of all mutations we can make
to a molecule, and all the ways we allow the atom coordinates to change within a mutation.

Once Osprey has a conformation, it evaluates that conformation according to
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

For now, task 2 is the main performance bottleneck in the Osprey code.
Making Osprey faster at task 2 usually translates well into overall reductions
in running time.

Task 1 can be very memory hungry though, since the space of all possible conformations
is exponentially large. So as designs include more and more choices for mutations,
we can easily exaust all the RAM available in even large server machines trying
to find conformations before even getting to task 2.

This document will describe enough background to understand how both tasks
work, and give an overview of the code that implements them.

But first, some preliminaries that are common to both tasks:


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

Note that the term *conformation* is slightly overloaded here. The term *conformation*
generally means an entire molecule where sets of atoms have been swapped.
However, In the context of a design position, the term *conformation* means a set of atoms
that can be swapped in.

In this example conformation space, there are a total of 4 * 8 = 32 conformations,
since Residue 42 can be one of 4 different options, Reisdue 17 can be one of 8 different
options, and both of these choices are completely independent.

Also in this example conformation, there are a total of 2*1 = 1 sequences. The sequences
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


## Specializing into graphs or minimizations

That's it for the preliminaries.
From here, head over to either task to dive deeper.

**Task 1:** [Find a conformation](optimizing-graphs.md)

**Task 2:** [Evaluate the conformation](optimizing-minimizers.md)


## Part 2. The code

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
explained in the later section on compilation.


### Code shared between both tasks

 * **Conformation Space** `confspace.compiled.ConfSpace`
   * This class reads the conformation space description from a file that's
     prepared by an entirely separate [GUI application](https://github.com/donaldlab/osprey-gui).
     If you wanted to prepare a new conformation space from scratch, you'll need to run
     the GUI program. But for the purposes of optimization, it's sufficient to use the
     previously-prepared conformation spaces that are built into Osprey's test suite
     and benchmarks.
     
 * **Conformation** `confspace.compiled.AssignedCoords`
   * The atom coordinates and degrees of freedom that embody a conformation.

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


## Part 3. Ok great, how do I compile?

First clone the git repo. Then checkout the `jdk14` branch.
The newest code isn't in the `master` branch yet, because not all
the Osprey devs are using JDK14 (Java development kit, ie the compiler
and runtime for Java) yet.

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


## Part 4. Getting more info

If you have any questions, want some clarification, or need any more information,
don't hesitate to send a message to `osprey@cs.duke.edu`.
It's a mailing list monitored by several people who will be able
to answer your questions.
