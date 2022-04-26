+++
title = "COFFEE"
weight = 1
+++

## Cluster Operations For Free Energy Estimation

COFFEE is a design tool that calculates the free energy of a conformation space (which is,
equivalently, a partition function calculator) that can efficiently use the resources of a
whole compute cluster, rather than just one compute server.

The name "COFFEE" itself is just a temporary code name for the project, and will likely be changed
as the tool becomes more mature.

COFFEE uses ideas from previous design tools in OPSREY:
 * MARK\*: recursive K*
 * SOFEA: computing in bounded memory, early stopping with partial results
 * COMETS: multi-state conformation spaces

There are also a few key ideas that distinguish COFFEE from previous design tools in OSPREY:
 * COFFEE uses distributed computing to use compute cluster resources efficiently.
 * COFFEE trades off the epsilion approximation of the partition function to achieve the bounded-memory goal.
 * COFFEE relaxes the best-first approach of A\* for conformational search into a good-enough-soon approach.

The main goal of COFFEE is to be a very practical design tool. Theoretical and mathematical
guarantees are sacrificed to provide conveniences for the designer. For example, many previous design tools
in OSPREY attempt to compute rigorous approximations to partition function values. To compute the
approximations in a reasonable amount of time, the designer has to hope that sufficient computational
resources are available, and that the design calculation completes without exhausting memory resources.
The designer usually has no a priori knowledge of how much time or memory is needed to complete design calculations,
and stopping the computation early usually gives no useful results. So running a design becomes a gamble of
whether or not the design can finish in a reasonable amount of time with the compute resources available,
but the designer can't know if the computation will finish until after spending the resources, both time and compute.
The best strategy we can offer to OSPREY designers is to start with the smallest possible version of the design
and try to get the computation to finish. Then increase the size of the design incrementally until the design no
longer finishes in a reasonable amount of time or runs out of memory.

In practice, regardless of the design problem, a designer generally only has access to a fixed amount of computing
resources, and is only patient enough to wait a week or two for results. COFFEE attempts to accommodate these
realities in a couple ways. First, the design of COFFEE handles the most memory-hungry operations (the conformation search)
with a bounded-memory system called the NodeDB so that the computation can continue even when the memory
allocated to the NodeDB has been exhausted. The trade-off to using bounded memory is that after NodeDB runs out
of memory, the free energy estimation loses precision. So the epsilon approximation is sacrificed for the convenience of
having a computation that provides a result even after running out of memory.

COFFEE also maintains partial results during the computation in a database called SeqDB, so stopping the
computation early will still provide useful results to the designer. Even with bounded-memory and early-stopping
properties, COFFEE can't give the designer any useful predictions for how long or how many resources
a design computation will need to reach a desired precision, but those properties take much of the gambling out
of designs because the computation will (barring any bugs in the code) end with a concrete result that might be
still precise enough to be useful too. So designers can start designs using COFFEE knowing the effort won't be
totally wasted, and if the results aren't as precise as needed, they can try again with a smaller design. 

NodeDB and SeqDB are explained in more detail in the following sections, along with references to the code that
implements them, and other components of COFFEE too.


## NodeDB: The Node Database

NodeDB is a custom distributed database that holds the state of the conformation search for a free
energy computation. On each member of the compute cluster, NodeDB uses a fixed amount of memory and a system
to manage that memory, to handle the case when free space is exhausted. NodeDB is implemented in the Java
`.coffee.nodedb.NodeDB` class.


### Memory Management

Conformation search in COFFEE using NodeDB is very similar to the recursive K\* approach used in MARK\*.
The conformation search is represented as a tree where nodes represent conformation assignments to design
positions, and nodes are scored using bounds on contributions to the partition function value.

Traditionally, the conformation search tree is searched using an A\*-like approach using a priority queue
to efficiently remove the best-scoring nodes first. But in COFFEE, to implement the bounded memory system, NodeDB
uses a double-ended priority queue that can efficiently remove both the best-scoring and worst-scoring nodes.
When NodeDB needs to store new nodes but has no more available space, NodeDB will just discard some of the
worst-scoring nodes to free up extra space. The hope is that the worst-scoring nodes represent negligible
contributions to the partition function value, so their absence will not be missed. However, by discarding
nodes from the conformation search, COFFEE is unable to guarantee any epsilon approximation to the partition
function value. In the trivial, but extreme, case where epsilon is set to exactly zero, all nodes must be
considered, and none can be discarded entirely.

The double-ended priority queue and memory management is implemented in the Java
`.coffee.nodedb.PriorityDequeFixedIndex` class.


### Partition function calculation is easier than GMEC search in some ways

For older GMEC-style conformational searches, having a strict iteration order of conformations in increasing
order of conformation energy (or bounds on energy) was vitally important to the correctness of the search.
However, for partition function calculations, conformation search order is important only for the performance
of the approximation, but not for correctness. In the trivial, but extreme, case when epsilon is exactly zero,
all conformations must be considered to correctly compute a partition function value, but since partition function
calculation is inherently a sum, the order in which those conformations are computed has no effect on the
resulting partition function value.  However, when epsilon is non-zero, computing conformations with low energies
first gives an accurate approximation to the sum faster, since conformations with high energies make negligible
contributions to the sum.

Therefore, an efficient partition function calculator need only to avoid computing conformations with high
energies, while still computing all the conformations with low energies, but strict ordering of conformations
by energy is unnecessary, especially when computing that strict ordering imposes high computational costs.
If A\* search is a best-first search, then we'll call this more relaxed philosophy a good-enough-soon search.
Computing the best (i.e. absolute lowest-energy) conformation first isn't as important for partition function
calculation as computing the good-enough (i.e. low-energy) conformations early in the search (i.e. soon), since
the summation nature of partition function calculation is insensitive to the order in which the low-energy
conformations are computed.


### Distributed Database

To make efficient use of cluster resources, NodeDB can be distributed across a cluster of multiple machines.
Each member of the cluster allocates a separate pool of memory for the NodeDB, and each node is stored in just
one NodeDB member. At the start of the search, all NodeDB members are empty since they're not storing
any nodes. As the partition function computation proceeds, NodeDB performs two primary tasks:

 * Store nodes that are created by the search
 * Retreive nodes corresponding to low-energy conformations

When storing new nodes, NodeDB prefers to store the node on the cluster member where the node was created, if
space is available. If the local NodeDB member doesn't have enough free space to store the node, then the
node is sent over the network to another NodeDB member with available free space. The particular choice of NodeDB
member isn't super important, but heuristically, the NodeDB member with the most available free space is chosen.
These heuristics require that each NodeDB member keep statistics on the state of all other NodeDB members.
The statistics tracking is implemented in the Java `.coffee.nodedb.Neighbors` class.
Finally, if all NodeDB members have exhausted their storage space, the local NodeDB member will remove high-energy
nodes from its own local storage to make extra space for the new nodes.

When a cluster member wants to compute new conformations, it will query NodeDB to get a batch of nodes that
correspond to low-energy conformations. Since the strictly lowest-energy nodes may reside in the storage of
multiple different NodeDB members, it's not computationally efficient to retrieve them all exactly. Instead,
the cluster member simply requests a batch of lowest-energy nodes from the single NodeDB member whose current
statistics advertise it as having the lowest-energy nodes in the cluster. Then the cluster member processes those
nodes and hopes for the best. Once all the nodes in the batch are processed, and any resulting new nodes are
stored in NodeDB, the process repeats by requesting more low-energy nodes until the conformation search is
completed. Depending on how often NodeDB member statistics are refreshed, the retrieved block of nodes may or may
not actually contain the lowest-energy nodes currently known to the conformation search. But probably the retrieved
nodes will be close to the truly lowest-energy node at any given moment. In this distributed environment where the
strictly lowest-energy nodes may be difficult to retrieve exactly, the strict node ordering required by a traditional
best-first A\* approach may have significant performance penalties. But the more relaxed good-enough-soon
approach is naturally a good fit here, since batches of low-energy nodes can be retrieved efficiently from single
NodeDB members.


## SeqDB: The Sequence Database

Another major component of COFFEE is SeqDB. The sequence database records the results of the partition function
calculation. Once a cluster member has processed a batch of nodes from NodeDB, the resulting bounds on
the partition function value are recorded in SeqDB. SeqDB is not distributed like NodeDB is. SeqDB is hosted
on just one member of the cluster. In terms of data structures, SeqDB is essentially a map from sequences
(which may be partially assigned) to bounds on partition function values.
SeqDB is implemented in the Java `.coffee.seqdb.SeqDB` class.

At the start of the computation, SeqDB is empty. The first step of the computation computes partition function
value bounds for the totally unassigned sequence corresponding to the root node of the conformation tree.
After they are computed, the root node bounds are saved to SeqDB at the totally unassigned sequence.
As the computation progresses and the interior of the conformation tree gets explored, bounds get computed
and refined for the partially-assigned sequences corresponding to the interior nodes of the conformation tree.
The bounds stored at any partially-assigned sequence represent the bounds of any unexplored conformation subtree
that matches the partial sequence. As the conformation subtree gets explored, the uncertainty of the
partition function value is reduced. When a node is processed, its bounds get subtracted from the corresponding
sequence assignment in SeqDB. Then, when the new children nodes are created, their bounds are added to their
corresponding sequence assignment in SeqDB, which will have one more assignment than the parent sequence.
Typically the sum of the bounds of all the child nodes are smaller than the bound of the parent node, so node
processing reduces the overall uncertainty of the partition function value. Finally, as the leaf nodes in the
conformation tree are searched, energy minimization of the corresponding conformation results in exactly zero
uncertainty of the partition function contribution of that conformation.

To get the total bound on the partition function value for a sequence, you need to not only look at the bound
in SeqDB for the fully-assigned sequence, but also combine it with bounds at all the preceding partially assigned
sequences (all the way to the totally-unassigned sequence) to account for the uncertainty of the unexplored subtrees
of the conformation tree for that sequence.

The partition function computation is complete when the bounds on the sequence are sufficiently precise.
SeqDB itself is saved as an output of the COFFEE computation, so the partition function value bounds persist
even after the computation finishes, or is aborted. Therefore, if the computation is stopped early, SeqDB can
still be queried later to get the best-known partition function value bounds on sequences at the time the
computation was stopped.


### Multi-state and multi-sequence

Not only can SeqDB hold the partition function value bounds for one sequence, it can also hold bounds
for multiple sequences and multiple design states simultaneously. For example, the typical affinity design
has three design states: the designed molecule, the target molecule, and the molecular complex.
A multi-state design may have even more states for positive and negative design.

The conformation tree used in the calculation need not be limited to a single sequence either. If a conformation
spaces includes mutations, the full conformation tree for the conformation space will also describe multiple
sequences. COFFEE can compute partition functions for all of the sequences simultaneously and store all the
bounds in SeqDB. In a multi-sequence conformation tree, COFFEE will not search the sequences in any explicit
order. Instead, bounds on sequence partition function values will be refined when the conformation nodes
for those sequences get retrieved as the low-energy conformations from NodeDB. So in effect, COFFEE will tend
to search sequences that have higher partition function values (or lower free energies) first.

COFFEE can't search different design states simultaneously though. Partition function values for different
design states must be computed serially. The next component of COFFEE we'll talk about is the Director, which,
among other things, is responsible for telling COFFEE which design state to search at any given time.


## The Director

Another major component of COFFEE is the Director. The director coordinates the node processing work among the
cluster members in the COFFEE computation. The metaphor here is the director of a band (think big band, not
rock band). Communicating with each member of the band individually is inefficient, so the director sends the same
signal to all members of the band at once. Each member of the band receives the shared signal and decides what to
do independently of the other members of the band. The director in COFFEE works much the same way. It sends one
signal to all the members of the cluster, and the cluster members use that signal to decide how to participate
in the COFFEE computation. The director is defined by the Java `.coffee.Coffee.Director` interface, and
director implementations live in the Java `.coffee.directors` package.

Computational work in COFFEE is divided between two different systems. The director send signals (aka directions)
to the whole cluster. Those directions are received by the node processors. One node processor runs on each member
of the cluster, and one director runs on just one member in the cluster. The node processors are responsible for
executing the directions given by the director, including retrieving nodes from NodeDB, computing bounds on
partition function values, and minimizing conformation energies. The node processor is implemented in the
Java `.coffee.NodeProcessor` class.

COFFEE was designed to support different implementations of the director, so different kinds of design computations
could be performed using the same NodeDB, SeqDB, and node processors. For example, the Java
`.coffee.directors.PfuncDirector` class directs a COFFEE cluster to compute just a single partition function
from a conformation space. The Java `.coffee.directors.PfuncsDirector` class directs the cluster to compute
a list of partition functions from the same conformation space. The Java `.coffee.directors.KStar` is a complete
re-implementation of the K\* design algorithm using COFFEE. The Java `.coffee.directors.AffinityDirector` is
a very naive attempt to direct COFFEE to do a sequence search for an affinity design using COFFEE's multi-sequence
partition function calculation capabilities.
