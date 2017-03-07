/**@file   main.cpp
 * @brief  Main file for ELPP tests
 * @author Leonardo Taccari
  */

/* standard library includes */
#include <iostream>
#include <fstream>
#include <string>

/* user defined includes */
#include "type.h"
#include "elpp.h"
#include "graph.h"
#include "readers.h"
#include <algorithm>
#include <unordered_map>
#include <chrono>

/* namespace usage */
using namespace std;


int main(int argc, char** argv)
{
   srand(2);
   string filename = "";
   bool help = false;
   bool relax = false;
   bool pruning = false;
   bool has_origin = false, has_destination = false;
   NODE s, t;
   ElppForm formulation;
   int timelimit = 1200; 
   int max_cuts = 1;
   double epsilon = 0.001;
   
   int o = 0;
   int tmp = -1;
   while ((o = getopt(argc, argv, "f:p::c:k:m:e:n:s:t:T:a::b::r::")) != -1)
      switch (o)
      {
         case 'f':
            filename = optarg;
            break;
         case 'p':
            pruning = true;
            break;
         case 'c':
            tmp = atoi(optarg);
            if(tmp >= 0 && tmp < ElppFormulationName.size())
                formulation = ElppForm(tmp);
            break;
         case 'T':
            timelimit = atoi(optarg);
            break;
         case 'm':
            max_cuts = atoi(optarg);
            break;
         case 'e':
            epsilon = atof(optarg);
            break;
         case 's':
            s = atoi(optarg);
            has_origin = true;
            break;
         case 't':
            t = atoi(optarg);
            has_destination = true;
            break;
         case 'r':
            relax = true; 
            break;
         default:
            help = true;
            break;
      }

   if ( filename == "" || help == true )
   {
      cout << "Usage: " << argv[0] << " -f <file_name> [other options]" << endl;
      cout << "Note: if no origin(s) and/or destination(s) are specified, all combinations are solved." << endl;
      cout << "Options:" << endl;
      cout << "-f\t(mandatory) name of data file." << endl;
      cout << "-s\torigin" << endl;
      cout << "-t\tdestination" << endl;
      cout << "-T\ttime limit (default: 1200 seconds)." << endl;
      cout << "-c\tsubtour elimination constraints:\n\t   0: None\n\t   1: MCF\n\t   2: Separation via Strong Components [default]\n";
      cout << "\t   3: Separation via Min-Cut\n\t   4: SF\n\t   5: RLT\n\t   6: MTZ\n\t   7: DL\n\t   8: DFJ\n\t   9: MCFsep" << endl;
      cout << "-e\tcut violation tolerance" << endl;
      cout << "-m\tmax number of cuts to be added in a callback. set -1 to add all violated inequalities (default: 1)" << endl;
      cout << "-p\tprune nodes that cannot be on an (s,t) path" << endl;
      cout << "-r\tsolve LP relaxation" << endl;
      cout << endl << "Examples:" << endl;
      cout << "\t" << argv[0] << " -f data/toy.dat -s 1 -t 5 -c 2" << endl;
      cout << "\t\tcompute the Longest Path over the graph in toy.dat" << endl;
      cout << "\t\tfrom node 1 to node 5 with formulation 2 (SC)" << endl;
      exit(0);
   }

   switch(formulation)
   {
       case NONE:
           cout << formulation << "- No subtour elimination" << endl;
           break;
       case MCF:
           cout << formulation << "- Using MCF formulation" << endl;
           break;
       case SC:
           cout << formulation << "- Using dynamic SEC separation via Strong Components" << endl;
           break;
       case MinCut:
           cout << formulation << "- Using dynamic SEC separation via Min-Cut" << endl;
           break;
       case SF:
           cout << formulation << "- Using SF formulation" << endl;
           break;
       case RLT:
           cout << formulation << "- Using RLT formulation" << endl;
           break;
       case MTZ:
           cout << formulation << "- Using MTZ formulation" << endl;
           break;
       case DL:
           cout << formulation << "- Using DL formulation" << endl;
           break;
       case DFJ:
           cout << formulation << "- Using DFJ formulation" << endl;
           break;
       case MCFsep:
           cout << formulation << "- Using MCF formulation with row-generation" << endl;
           break;
       default:
           formulation = ElppForm::SC; 
           cout << "No valid formulation selected. Using default, Strong Component (SC)." << endl;
   }

   unordered_map<NODE_PAIR, double> cost;
   std::shared_ptr<Graph> G = std::make_shared<Graph>();
   if(! Readers::read_graph(filename, G, cost) ) {
       cout << "Something wrong reading from file." << endl;
       exit(1);
   }
   if(!has_origin || !has_destination) {
       cout << "Both s and t must be specified." << endl;
       exit(1);
   }
   if(!G->has_node(s)) {
       cout << s << " not in graph!" << endl;
       exit(1);
   }
   if(!G->has_node(t)) {
       cout << t << " not in graph!" << endl;
       exit(1);
   }
   
   cout << "Compute longest path from node " << s << " to node " << t << "." << endl;
   IloEnv env;
   if(pruning) {
       cout << "Before pruning: " << G->num_nodes() << " nodes and " << G->num_arcs() << " arcs." << endl;
       G->prune(s,t);
       cout << "After pruning: " << G->num_nodes() << " nodes and " << G->num_arcs() << " arcs." << endl;
   }
   ElppSolver elpp_solver = ElppSolver(env, NODE_PAIR(s,t), G, formulation, relax, timelimit, epsilon, max_cuts);
   elpp_solver.update_problem(cost);
   auto timer_start = std::chrono::high_resolution_clock::now();
   if(relax)
       elpp_solver.solveLP();
   else
       elpp_solver.solve();
   auto timer_end = std::chrono::high_resolution_clock::now();
   double duration = std::chrono::duration<double>(timer_end - timer_start).count();

   cout << ElppFormulationName[formulation] << "\t: ";
   if(elpp_solver.getStatus() == IloAlgorithm::Optimal)
       cout << elpp_solver.getObjValue() << " " << duration << "s" << endl;
   else
       cout << "Status: " << elpp_solver.getStatus() << ", not solved to optimality! " << duration << endl;
   
   env.end();
   return 0;
}
