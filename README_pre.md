Solving Longest Path problems with Integer Programming
====

This repository contains code that uses several Integer Programming formulations to solve
Longest Path problems over directed or undirected graphs.

The problem
---

Given a graph *G=(V,A)*, an origin node *s* and a destination *t*, the Longest Path problem 
requires to find the longest path connecting *s* and *t*.
Unlike its minimization counterpart, for which well-known polynomial algorithms like Dijkstra exist, 
the Longest Path problem is NP-hard, and related to problems such as the Traveling Salesman Problem. 

Several variants exist: the graph can be weighted or unweighted, directed or undirected.
The problem is also sometimes known as the Elementary Shortest Path Problem,
where an elementary path is a path where no node is visited more than once. 
Specifying that the path must be elementary hints at the fact that we are solving a problem where
cycles of negative cost may arise, and we are focusing on paths without cycles.

From now on, I'll stick to the maximization version, and refer to the problem as ELPP (Elementary Longest Path Problem).

NP-hardness is no excuse to avoid solving a problem to optimality (or, at least, attempting to!).
One way to solve ELPP with an exact algorithm, i.e., to **proven global optimality**, 
is by using the awesomeness of Integer Programming (IP) and modern branch-and-cut algorithms.

The core of a simple IP formulation is as simple as this:

\begin{align}
   &\max \sum_{(i,j)\in A}c_{ij}x_{ij}\label{espp:objective}\\
   & \sum_{\makebox[20pt]{$\substack{(i,j) \in \delta^+(i)}$}}x_{ij} - \sum_{\makebox[20pt]{$\substack{(j,i) \in \delta^-(i)}$}}x_{ji} = 
\begin{cases}
   1 & \text{if } i = s \\
   -1 & \text{if } i = t \\
   0   & \text{else} \end{cases} & \forall i \in V \label{espp:balance}\\
&\sum_{(i,j) \in \delta^+(i)}x_{ij} \leq 1 & \forall i \in V \label{espp:degree}\\
&  x_{ij}\in\{0,1\} & \forall (i,j) \in A, \label{espp:integrality}
\end{align}
where $c_{ij}\in \mathbb{R}$ are the arc costs, and
$x_{ij}$ are binary arc variables that take value 1 if the arc $(i,j)$ belongs to the path. 

This formulation, however, is not sufficient to guarantee that the optimal path will be *elementary*.

There are mainly two ways to ensure this is the case:
- polynomial-size formulations, that have a size (number of variables and constraints) which is polynomial w.r.t. the size of the graph *G*
- exponential-size formulations, that have an exponential number of **subtour elimination constraints**.

In the first case, the whole Integer Programming model can be built statically, for reasonably large graphs. 
However, these formulations turn out to be pretty challenging to solve.

In the second case, the whole model cannot be built (except for tiny instances). Nevertheless, it is possible to generate
the subtour elimination constraints dynamically, in a phase called separation. 
In practice, this means starting with no subtour elimination constraints. 
Then, during the branch-and-bound search, every time we find a solution x, we check whether it violates
any of the constraints we left out of the model. If this is the case, we add those constraints, and continue.

Finding a violated constraint is done in polynomial time solving a MinCut problem, or finding Strongly Connected components
in the support graph defined by the solution.

This kind of formulation is the way to go for similar problems (TSP in primis), and it works quite well also
for Longest Path problems.

Installation and usage
---

The class `ElppSolver` provides methods to build and solve a Longest Path problem with one of several formulations.

Before solving a problem, `update_problem()` must be called once to set the objective function
(arc costs) and, optionally, bounds on the arc variables.
This allows the solver to be called repeatedly, with modifications
only in the variable bounds and/or objective coefficients, without the need to rebuild
the whole model. It also means that the solver can use MIP warm start information,
that might be useful is the cost coefficients haven't changed much.
The possibility to modify the problem and calling it again might be handy in the context, 
e.g., of a branch-and-price, where the pricing solver has to be called several times.

In addition to the basic Elementary Longest/Shortest Path Problem, it is easy to add additional side constraints
via the method `add_constraint()`. 
Although in presence of additional constraints, you might want to consider a specialized algorithm,
a simple, out-of-the-box Integer Programming approach can give you pretty good results, and it's extremely flexible.

External dependencies
---

- [IBM Cplex Optimizer](https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) 
- [Boost library](http://www.boost.org)
- [Lemon graph library](http://lemon.cs.elte.hu/)

For Boost and Lemon, please follow the instruction provided on their website.
For Cplex, if you are a student or faculty (an academic email address should be sufficient), you can get access to the software via their Academic Initiative:
for students [here](https://ibm.onthehub.com/WebStore/OfferingDetails.aspx?o=9b4eadea-9776-e611-9421-b8ca3a5db7a1)
and for professors or researchers [here](https://ibm.onthehub.com/WebStore/OfferingDetails.aspx?o=6fcc1096-7169-e611-9420-b8ca3a5db7a1).

Requires a compiler with C++11 support. Tested on Linux machines with GCC 4.6 and above.

Using the code/data and citing
--

Feel free to use/modify this code (that comes with no guarantees :-). 
If you use it for research, consider citing the 
[article](http://dx.doi.org/10.1016/j.ejor.2016.01.003) where I used it to compare several Integer Programming formulations.
This repository also includes the datasets used in the article.

``` tex
@article{taccari2016espp,
  title={Integer programming formulations for the elementary shortest path problem},
  author={Taccari, Leonardo},
  journal={European Journal of Operational Research},
  volume={252},
  number={1},
  pages={122--130},
  year={2016},
  publisher={Elsevier}
}
```

A preprint version can be found on Optimization Online [here](http://www.optimization-online.org/DB_FILE/2014/09/4560.pdf).
