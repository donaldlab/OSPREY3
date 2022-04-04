+++
title = "Optimizing Osprey"
weight = 9
+++

TODO: update this to make more sense in the new context

For now, task 2 is the main performance bottleneck in the Osprey code.
Making Osprey faster at task 2 usually translates well into overall reductions
in running time.

Task 1 can be very memory hungry though, since the space of all possible conformations
is exponentially large. So as designs include more and more choices for mutations,
we can easily exaust all the RAM available in even large server machines trying
to find conformations before even getting to task 2.



## Specializing into graphs or minimizations

That's it for the preliminaries.
From here, head over to either task to dive deeper.

**Task 1:** [Find a conformation using graph theory](optimizing-graphs.md)

**Task 2:** [Evaluate a conformation using physics](optimizing-minimizers.md)
