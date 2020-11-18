# PsoPath: Particle Swarm Optimisation algorithm for the shortest path problem

```python
Copyright (C) 2020 by Constantine Kyriakopoulos
Version: 1.0.2
License: GNU GPL Version 3
```


## About the project

The shortest path problem is solved by many methods. Heuristics offer lower complexity 
in expense of accuracy. There are many use-cases where the lower accuracy is acceptable 
in return of lower consumption of computing resources.

Particle Swarm Optimisation (PSO) emulates the social behaviour of, e.g., a flock of 
birds, as a stochastic optimisation method [1]. Specifically, a particle is an entity 
representing a solution in the search space. Several particles cooperate inside an 
algorithmic flow to occupy positions close to the best solution. So, when a number of 
predefined iterations for evaluation ends (depending on topology resources like 
the node-number), the particle which provides the best position is the one to provide 
the solution to the problem. To achieve this goal, a fitness function evaluates every 
particle’s position at any time. Its velocity depends on neighbours’ positions, and 
the current and previous best position the particle had at the moment of evaluation.

This is a heuristic method, i.e., optimal results are not always feasible. According to 
topology's resources like the node and edge numbers, the proper numbers of iterations 
and particles must be used. Large numbers lead to paths with higher probability of 
being optimal but more computational resources are consumed.


## Prerequisites to build

There are only two requirements, i.e., the Boost Library and the availability of C++20 
or C++17 standard. Boost is utilised for parsing the JSON representation of the topology. 
Tested with clang 11 and libc++ from the LLVM project.


## Usage

Create a new instance of <em>PsoSystem</em> in your code passing as parameters the JSON 
topology file and the numbers of iterations and particles (default values are also provided 
but the shortest paths aren't returned under every topology size). Next, execute the method 
<em>path(src, dest)</em> where <em>src</em> is the source node and <em>dest</em> the 
destination to reach. This returns the valid path that particles converge to.


## Related work

```python
[1] Zhang, Y., Wang, S. and Ji, G. A comprehensive survey on particle swarm optimization 
algorithm and its applications. Mathematical Problems in Engineering, 2015.
```

## Changelog

<pre>
1.0.2  2020-1X-XX    Extended the interface for individual edge insertion
1.0.1  2020-11-17    Execution time improvement
1.0    2020-11-13    Initial public release
</pre>
