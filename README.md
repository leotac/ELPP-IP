Solving Elementary Shortest Path problems with Integer Programming
====

Solving Elementary Shortest Path, a.k.a. Longest Path problems, via mixed-integer programming
over directed or undirected graphs.

The Longest Path problem is a classic NP-hard problem, also known as any of the following:
- Longest [Simple|Elementary] Path (with positive cycles)
- Shortest [Simple|Elementary] Path (with negative cycles)

*If you want to nitpick, some people may say that simple != elementary.


The provided C++ class ElppSolver has a constructor that solves a Mixed-Integer Programming model
to compute a longest path between two nodes for a given directed graph.

Before solving a problem, update_problem() must be called once to set the objective function
(arc costs) and, optionally, bounds on the arc variables.
This allows the solver to be called repeatedly, more than once, with modifications
only in the variable bounds and objective coefficients, without the need to rebuild
the whole model. This also means that the solver can use MIP warm start information.
This is especially handy in the context, e.g., of a branch-and-price, where
the pricing solver has to be called several times.


*External dependencies:* [IBM Cplex]() and [Lemon graph library]()

Requires compiler with C++11 support.
Tested on Linux machines with GCC 4.6 and above.
