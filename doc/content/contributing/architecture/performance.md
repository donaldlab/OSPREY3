+++
title = "High-performance code in OSPREY"
menuTitle = "Performance"
weight = 9
+++


OSPREY relies on searching exponentially-sized combinatorial spaces to find sequences
and conformations of interest for molecular designs. Since these search spaces are so huge,
the runtime performance of OSPREY has a big impact on the quality of the results that can be found.
If OSPREY runs faster, the designer can design larger conformation spaces with more interesting options
to explore, and can hopefully get more interesting and biologically relevant results.
So there is a large incentive to invest in the engineering of OSPREY to improve runtime
performance as much as possible.


## What's the current bottleneck?

If we take the [two-task view of OSPREY]({{< ref "/#overview" >}}), for now,
task 2 is the main performance bottleneck in terms of time.
Making Osprey faster at task 2 usually translates well into overall reductions
in running time.

Task 1 can be very memory hungry though, since the space of all possible conformations
is exponentially large. So as designs include more and more choices for mutations,
we can easily exhaust all the RAM available in even large server machines trying
to find conformations before even getting to task 2.


## How does OSPREY handle high-level concurrency and parallelism?

### Within a single machine

At a high level, OSPREY uses task parallelism to distribute work among CPU cores
and GPU SMs. The `ThreadPoolTaskExecutor` is the main class that handles distributing
parallel workloads using a thread pool.

Most inputs to these parallel tasks are constant objects like the `ConfSpace` instances.
These objects are designed to be created once and then remain constant
throughout the rest of the program, serving as lookup tables for later computations.
Parallel tasks read these objects without sychronizing since there's no possiblity of
a race.

For varying objects like conformations, or lists of position interactions, these
objects are created as part of a task definition and moved across threads with the task,
so no two threads can access it at once (outside of the `TaskExecutor` system).
Other varying objects are created wholly whithin the task as it executes, and are
essentially thread-local.


### Across multiple machines in a cluster

OSPREY has also recently (re-)gained the ability to coordinate computations across a cluster
of networked machines, rather that just within a singe machine. OSPREY uses the [Hazelcast][hazelcast]
library in much of the newer code to implement low-level message-passing between cluster members.
For this application, Hazelcast is a more modern and easy-to-use alternative to the classic [MPI][mpi] systems
of yesteryear for the Java ecosystem.
Hazelcast also provides high-level distributed data structures (like lists and maps) that can be used
in distributed design algorithms.

[hazelcast]: https://hazelcast.com/open-source-projects/imdg/

OSPREY also has some integration with the [SLURM][slurm] job scheduler. OSPREY can use information from SLURM
to setup the inter-machine communication, even when launching multiple multi-machine SLURM jobs at the same time.

[slurm]: https://slurm.schedmd.com/overview.html


## OSPREY uses "native" code for performance-critical functions

In the Java ecosystem, any code that runs directly on the operating system and not the
Java Virtual Machine (JVM) is called *native* code.

Tragically, the object-oriented programming paradigm fails miserably when you really want
your code to go fast. So some of OSPREY's performance-critical functions are implemented
in C++ in more of a data-oriented programming style.

OSPREY has a few Java classes that translate the data needed for the C++ implementations
into contiguous blocks of memory suitable for C++ code to consume. Once the data buffers
are ready, the Java code sends the pointers to the C++ functions (defined via `extern "C"`)
over one of Java's FFI (foreign function interface) mechanisms called
[JNA (Java native access)](https://github.com/java-native-access/jna).
JNA loads the `.so`/`.dll`/`.dylib` libraries at runtime into the JVM (Java virutal machine) process and
binds the native functions to Java function stubs. It also handles type marshalling
for primitive types and pointers across the Java/C++ boundary.

The `NativeConfEnergyCalculator` and `CUDAConfEnergyCalculator`
classes are all the newest versions of the Java/C++ bridge class idea.
They're based on a [Foreign Memory Access API](https://openjdk.java.net/jeps/370)
that's a new feature of JDK14 (Java development kit, i.e., the compiler and runtime for Java).
These classes handle all the memory layouts and alignments
for the constant data that's needed by the C++ code. They use a custom struct definition
API on the Java side to help read/write the data from/into the raw memory buffers.


### Including on GPUs

OSPREY implements some of its energy calculations using very high-performance computing hardware,
like [CUDA][cuda]-capable [Graphics Processing Units (GPUs)][nvidia-gpus].

[cuda]: https://en.wikipedia.org/wiki/CUDA
[nvidia-gpus]: https://en.wikipedia.org/wiki/List_of_Nvidia_graphics_processing_units

Energy calculations are a bit tricky to parallelize on GPU hardware.
OSPREY's minimized energy calculations use [Cyclic Coordinate Descent (CCD)][ccd] which is a notoriously
linear algorithm. However, to make efficient use of GPU hardware, it's best to have thousands of calculations
running in parallel. Luckily, the forcefield calculations themselves are [embarassingly parallel][embarassingly-parallel]
with thousands (and more) independent elements to compute. So getting good performance for OSPREY on GPUs comes down
to using the parallelism inherent in the forcefield calculations and hiding the linearity of
the CCD algorithm as much as possible.

[ccd]: https://en.wikipedia.org/wiki/Coordinate_descent
[embarassingly-parallel]: https://en.wikipedia.org/wiki/Embarrassingly_parallel

The parallelization approach to energy calculation that has worked best for OSPREY so far
is to run one minimization per streaming multiprocessor (SM) on the GPU and rely heavily
on the `syncthreads()` mechanism to allow us to mix parallel and serial code and handle
temporal dependencies. Modern GPUs typically have tens of SMs on-chip, so not only
can OSPREY use the parallelism within a single minimization to saturate the SM hardware,
OSPREY can minimize multiple conformations in parallel to saturate all the SMs on a GPU
-- or even multiple GPUs. Since OSPREY is usually trying to minimize astronomical numbers
of conformations, this leaves us with no shortage of parallelism to fill up available GPU hardware.

Newer CUDA versions have introduced cooperative kernels that might give us
better options for parallelism, this hasn't yet been explored in OSPREY.
It's not clear yet if the multi-block synchronization mechanisms are fast enough
for our purposes. We're never hurting for more conformations to minimize concurrently,
so one SM per minimization has worked really well so far.
