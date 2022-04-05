+++
title = "The Architecture of OSPREY"
menuTitle = "Architecture"
weight = 7
+++

This section describes how OSPREY itself is organized internally.
It might serve as a nice introduction to the OSPREY code for a developer who is new to the project.


## Overview

![High-level overview of OSPREY components](overview.svg?height=800px)

OSPREY's purpose in life is to find protein sequences that have properties we're
interested in, like stability or affinity with other molecules. OSPREY does this
by doing two basic tasks over and over:

**Task 1:** [Find a conformation using graph theory](graph-search)

**Task 2:** [Evaluate a conformation using physics](conf-energy)

In this case, a [*conformation*](conformation-spaces) is a molecule, but we've replaced some of the atoms
in it with a different set of atoms. This allows us to describe mutations to a
molecule and discrete flexibility in the atom coordinates in an efficent, compact way.
OSPREY finds these conformations using various graph algorithms operating on a
[*conformation space*](conformation-spaces). A conformation space is the set of all mutations we can make
to a molecule, and all the ways we allow the atom coordinates to change within a mutation.

Once OSPREY has a conformation, it evaluates that conformation according to
a physical forcefield to update its scoring model about that conformation, and the
sequence of which the conformation is a part. To model the possibility that the
conformation we actually want isn't perfectly represented by the atom coordinates
we started with, OSPREY allows the atoms to move slightly from their original positions.
OSPREY uses a *continuous motion* to move the atoms, for example, by a dihedral angle
rotation around a rotatable bond. Using the forcefield and the continuous motions as
an objective function, OSPREY minimizes the conformation to find the energy of the best
nearby atom coordinates.

Once OSPREY has enough information about different conformations and their minimized
energies, it can claim that a certain sequence has or doesn't have the properties
that the design project seeks.


## The code

### Overall architecture

OSPREY is almost entirely written in [Java][java] (see `/src/main/java`). Newer parts of OSPREY, like the Desktop software,
are written in [Kotlin][kotlin] though, which in many ways is like a successor to the Java language.
There's a thin layer of [Python][python] on top (see `/src/main/python`) acting as a scripting API
for designers and data scientists.
A small part of performance-critical code has been copied into to C++ (see `/src/main/cc` and `/src/main/cu`)
and optimized for various hardware platforms, including CPUs and GPUs.

[java]: https://www.oracle.com/java/technologies/
[kotlin]: https://kotlinlang.org/
[python]: https://www.python.org/

Since OSPREY is largely developed by academic researchers working on projects in biology and chemistry,
OSPREY has naturally adopted a copy-on-write approach to improving systems. Rather than upgrade
an existing system in-place with a painful refactoring that may distrupt a researcher's progress towards
a project (like a paper, or a PhD), it's generally far easier to implement new features in independent code,
and then convince other researchers to upgrade to the new version between research projects, rather than during them.

This copy-on-write philosophy has led to the proliferation of many different implemenations
of similar features in OSPREY's code over time. And due to the needs of academic
reproducibility, it's essentially impossible to delete older implementations that are
no longer actively in use. So OSPREY's code is littered with a complete history of
various implementations for concepts like conformation spaces, energy functions, and
minimizers. This can all be extremely confusing to new OSPREY developers, so one
of the purposes of this document is to point you quickly in the direction of
the modern and actively in-use implementaions of OSPREY's basic features.


### Where's the active code?

Here's a quick roadmap of the source code that implements the most modern versions of
OSPREY's basic features. The sources will be given by the source root folder and the
verbose Fully-Qualified name of the relevant class, but the common package `edu.duke.cs.osprey` has been omitted.


#### Conformations and Conformation Spaces

* **Conformation Space:** `/src/main/java` / `.confspace.compiled.ConfSpace`
    * This class reads the conformation space description from a file that has been
      prepared by the conformation space preparation tools.

* **Conformation:** `/src/main/java` / `.confspace.compiled.AssignedCoords`
    * The atom coordinates and degrees of freedom that embody a conformation.


#### Conformation Search

* **Energy Matrix:** `/src/main/java` / `ematrix.EnergyMatrix`
    * Implementation of the lookup table for the A* scoring function

* **Energy matrix calculator:** `/src/main/java` / `ematrix.compiled.EmatCalculator`
    * Functions to compute the energy matrix from a conformation space
      and a conformation energy calculator instance.

* **A\* search:** `/src/main/java` / `astar.conf.ConfAStarTree`
    * Osprey's implemenation of A* search
    * There's also an implementation of a memory-bounded version of A*
      in there too, that trades time for space.
    * There are options to choose different heuristics for A* too.


#### Energy Calculation and Minimization

* **Conformation energy calculator:** `/src/main/java` / `energy.compiled.CPUConfEnergyCalculator`
    * This is the Java implementation of the newest energy calculator and minimizer.
      It's the slowest option available, but the design is simple, it's pure Java,
      and it makes a good baseline for benchmarks.

* **"native" energy calculator:** `/src/main/java` / `energy.compiled.NativeConfEnergyCalculator`
    * This is the reference implementation of the C++ version of the energy calculator
      and minimizer. It's basically a straight port of the Java code into C++ and it's
      totally unoptimized. It also makes a good baseline for benchmarking.
    * It was created originally as a stepping stone for the CUDA implementation,
      so it's written in a very C-style subset of C++, and uses exactly zero libraries.
    * The native side of the code is at `/src/main/cc/ConfEcalc`.

* **CUDA energy calculator:** `/src/main/java` / `energy.compiled.CudaConfEnergyCalculator`
    * This is the CUDA version of the energy calculator. It's been optimized quite a bit
      and is currently the fastest option available by far.
    * Of course, it requires CUDA-capable GPU hardware to run.
    * The native side of the code is at `/src/main/cu/CudaConfEcalc`.


#### TODO: other categories? design prep? design algorithms?


#### Parallelism and high-performance code

* **Thread Pool:** `/src/main/java` / `.parallelism.ThreadPoolTaskExecutor`
    * OSPREY's main tool for task-level parallelism and distributing work across
      multi-core CPUs.
    * In a nutshell, a task is shipped off to a worker thread where it gets processed.
      Then the task's result is put on a queue where a single listener thread drains the queue
      and sends the result to a callback function specified by the submitter of the task.

* **Parallelism:** `/src/main/java` / `.parallelism.Parallelism`
    * A simpler helper class that describes options for parallelism, like
      the desired size of the thread pool.
    * The `streamsPerGpu` option isn't used anymore in the latest CUDA code.
      Stream management is handled internally by the CUDA code now based on
      hardware specs queried at runtime.

* **Cluster Communication:** `/src/main/java` / `.parallelism.Cluster`
    * OSPREY's main tool for using Hazelcast to coordinate communication between machines
      in multi-machine jobs. Supports integration with SLURM as well.
