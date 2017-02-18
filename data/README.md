Data
===

**Note:** all instances assume that an *(elementary) longest path* is sought, so the problem must be formulated as a maximization or one might need to flip the sign of the optimization or of the weight.

Each .dat file contains a network topology, with, respectively:
+ number of nodes and arcs
+ nodes
+ arcs and weights
Each .st file contains one or more (origin,destination) pairs.

**drexl:** pool containing four subsets: *prc-f*, *prc-p*, *prc-l*, *rnd-d* from Drexl.   
For each subset, there are 3 graph sizes: 25, 50, 100.  
Note that the actual number of nodes in the pricing instances is N + 2. Not sure why Drexl didn't count those. Also, sometimes the graphs are complete, sometimes a few arcs (incoming or outgoing from source or destination) are missing. 
For each subset and size, there are 30 instances with different arc weights (both positive and negative).  
Each instance corresponds to a problem where one wants to determine the ELPP between *s* = 1 and  *t* = largest_node_index.  

**rand**: pool containing randomly generated network topologies with arc density 1/15 and size (nodes): 50, 100, 200, 500, 1000.  
Since the graphs are randomly generated, there are no prescribed source and target nodes, and source and target can be picked randomly.
The (s,t) pairs used in my paper are in the .st file.

**snd_pricing:** pool containing pricing subproblems on graphs with both positive and negative weights.  
There are 5 subsets, corresponding to topologies from the SND library: atlanta, france, geant, germany, nobel-us.  
For each of them, 30 instances are collected. 
Most of the graphs are very sparse and should be easy to solve.

**sp:** smallest instance ("rome") from the shortest path DIMACS challenge. All weights are positive.
