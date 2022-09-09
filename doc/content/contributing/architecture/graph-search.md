+++
title = "Graph Search"
weight = 2
+++


## Task 1: Find a conformation using graph theory

Osprey uses several graph algorithms to find conformations to evaluate in the conformation
space, including primarily [A* search](https://en.wikipedia.org/wiki/A*_search_algorithm).

Once we know enough of the lowest-energy conformations in the conformation space for a
given sequence, that gives us enough information to characterize that sequence.
So the main goal of task 1 is to quickly find the lowest energy conformations in the
conformation space. Either in the whole space, or in the subspace that describes a
particular sequence.

### How Osprey uses A* search

If we totally order the design positions in the conformation space (an arbitrary decision,
there's no natural order to the design positions), the conformations at each design position
become the nodes at one level of the tree. Our example conformation space can be
represented by the following tree.

 * root
   * Alanine, conf 1 `posi=0,confi=0`
     * Tyrosine, conf 1 `posi=1,confi=0`
     * Tyrosine, conf 2 `posi=1,confi=1`
     * ...
     * Tyrosine, conf 8 `posi=1,confi=7`
   * Valine, conf 1 `posi=0,confi=1`
     * Tyrosine, conf 1 `posi=1,confi=0`
     * Tyrosine, conf 2 `posi=1,confi=1`
     * ...
     * Tyrosine, conf 8 `posi=1,confi=7`
   * Valine, conf 2 `posi=0,confi=2`
     * Tyrosine, conf 1 `posi=1,confi=0`
     * Tyrosine, conf 2 `posi=1,confi=1`
     * ...
     * Tyrosine, conf 8 `posi=1,confi=7`
   * Valine, conf 3 `posi=0,confi=3`
     * Tyrosine, conf 1 `posi=1,confi=0`
     * Tyrosine, conf 2 `posi=1,confi=1`
     * ...
     * Tyrosine, conf 8 `posi=1,confi=7`

If we had chosen the opposite order for the design positions, we would have gotten
a different tree, but the set of conformations represented (and crucially their energies)
are still the same. So the order of the tree does not affect the solution to the problem,
so therefore we can optimize the position ordering in any way we like using heuristics.

The path from the root node to any leaf node in the tree describes a unique conformation.
For example, The path `root` -> `posi=0,confi=0` -> `posi=1,confi=7`
corresponds to the conformation `[0,7]` where the implicit order of design positions
is `[0,1]`.

Osprey uses A* search to essentially sort the leaf nodes (corresponding to conformations)
according to a scoring function that is a *lower bound* on the minimized energy
of the conformation. If you run A* the first time, you'll get the optimal leaf node in
the tree (assuming the heuristic function is sound). But if you keep iterating A* search
beyond the optimal solution, you'll get the *next* lowest scoring node, and the next
lowest scoring node, and so on in order.

The hope is that by sorting the conformations by lower bounds on their energies,
the lower bounds will be close to the real energies, and we can find the low-energy
conformations quickly.

However, in practice, it's often the case the the lower bounds are very loose.
We may do all the work of sorting the conformations by lower bound and enumerating
them in order only to find their minimized energies are much higher than the bound.
We may have to enumerate and minimzie many many conformations before we find the
conformations that are truly low-energy.


### Graph structure of the A* scoring function

The scoring function we use to sort leaf nodes is based on a different graph.
Each node in the graph corresponds to a design position in the conformation space,
plus we'll add one more node for a "static" design position, which represents
all the atoms that aren't in any of the design positions. The graph for our
example conformation space looks like this:
```
   S
  /  \
 0 -- 1
```
Where `S` is the "static" design position, and `0` and `1` are the design positions
corresponding to residues 42 and 17 in the protein respectively.

The edges in the graph represent interactions between the design positions.
The edges between the numbered design positions form a complete graph in general.
And each numbered design position has an edge to the static design position.
There are also self edges for each node, but I haven't shown them here
because ASCII art is hard.

In this graph, a conformation is represented by making an assignment to each design
position. For example, the conformation `[3,7]` induces this graph:
```
     S
  /     \
 0=3 -- 1=7
```
Therefore, each edge in the graph is represented not only by its incident nodes,
but also the assignments at those nodes.

The scoring function works by computing a scalar weight for each possible edge in the
graph. For example:
 * the `(S,0=3)` edge might have a weight of `5.0`
 * the `(S,0=5)` edge might have a weight of `4.2`
 * the `(0=3,1=4)` edge might have a weight of `3.5`
 * the `(0=3,0=3)` edge might have a weight of `8.3`
 
and so on. However, there are no self edges between different assignments of the same
design position. For example, there is no `(0=3,0=4)` edge.

The score of a conformation is simply the sum of all edge weights after all
design positions have been assigned. Due to the way we compute the edge
weights for this graph, the sum of weights is guaranteed to be a lower bound
on the minimized energy of the entire conformation.


### The computed edge weights form an "energy matrix"

Edge weights are computed by miniming only part of a conformation, rather than
minimizing an entire conformation. For example, for the edge `(S,0=3)`, Osprey
creates a partial molecule containing the atoms for position `S` and position `0`.
Then osprey computes the minimized energy for this partial molecule and stores the
result in a lookup table. This lookup table is called an "energy matrix".

For more details on the energy minimization process, see the later section on Task 2.

The energy matrix is a simplfication of the graph where the "static" position/node
has been removed. Therefore, we have to distribute the edge weights from the edges
incident to the "static" position to other edges in the graph.

The simplest way Osprey does this is by moving all the static-incident edge weights
onto the self edges. So the energy matrix is a collection of two different types of
edge weights:
 * **singles:** (or "one body" energies)
   * for example `(S,0=3) + (0=3,0=3)`
   * indexed in two dimensions: `posi1,confi1`
 * **pairs:** (or "pairwise" energies)
   * for example `(0=3,1=2)`
   * indexed in four dimensions: `posi1,confi1,posi2,confi2`
   
Osprey also has settings to use different distributions of edges onto the energy
matrix that can result in much tighter lower bounds for the A* search. Older parts of the code
call this process "energy partitioning" (embodied in the `EnergyPartition` class), and
the newest parts of the code calls this process the "position interation distribution"
(embodied in the `PosInterDist` class). The idea there is that when the energy function is
described at a high level using the interations between design positions in the conformation
space, those interations can be mapped to different "distributions" of low-level interations
between individual atom pairs. In the newest energy calculation code, these distribution
settings can be configured by any class that accepts a position interaction generator
(embodied in the `PosInterGen` class).

Once the energy matrix is computed, Osprey has enough information to run the A* algorithm
and sort the conformations by their A* scores. All the edge weights are readily available
in a lookup table.


### The energy matrix is the basis for the A* scoring function

The A* scoring function is broken into two parts:
 * "G" score
   * The G score is simply the summed energy matrix entries
     corresponding to the assigned design positions.
 * "H" score, or "heuristic" score
   * The H score is a lower bound on the summed energy matrix entries
     where at least one design position is unassigned.
   * This is the "heuristic" function of A* search.

The G score simply reads the values out of the energy matrix
corresponding to the node assignments. The G scoring function implementation
is linear in the number of design positions and is very fast to compute.

Opsrey's implementation of the H score is polynomial in the number of design positions
and is by far the bottleneck of the A* computation. The H function must optimize over
the possible choices avaiable to each design position to compute the lower bound on the
scoring function for the whole conformation.

The heuristic function used to be the main bottleneck in Osprey overall,
but much work has gone into optimizing the Java code, so it's no longer the main
bottleneck in modern versions of Osprey.


### A* is a total memory hog

In a design with many design positions, A* can eat an insanely huge amount of
memory very quickly, so we've tried lots of different tricks to deal with that.

We've optimized our A* node implementation to save as much memory as possible
(at least as well as you can do in Java, which is admittedly not the best
language for memory efficiency).
We've integrated an I/O-efficient implementation
of a priority queue into Osprey to spill A* memory to local storage (eg SSD, NVME).

Newer attempts at conformation search in OSPREY deviate slightly from the traditional A*-based approaches.
One experimental approach, codenamed "COFFEE" for now, tries to relax the requirements
of "best-first" search in the hopes that all we really need is "good enough-soon" search.
Once strict sorting of possible conformations in the search is no longer required,
that gives us more options to store the possible conformatations efficiently,
and distribute the computation across parallel hardware -- even across an entire compute cluster.

TODO: add a link here to the upcoming COFFEE section.
