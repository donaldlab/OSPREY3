+++
title = "High-performance code in OSPREY"
menuTitle = "Performance"
weight = 9
+++


TODO: introduce this section?


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

For varying objects like conformations, or lists of position interations, these
objects are created as part of a task definition and moved across threads with the task,
so no two threads can access it at once (outside of the `TaskExecutor` system).
Other varying objects are created wholly whithin the task as it executes, and are
essentially thread-local.


### Across multiple machines in a cluster

TODO: describe hazelcast and OSPREY's cluster tools
