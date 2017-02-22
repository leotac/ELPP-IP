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
#include <algorithm>
#include <random>
#include <map>
#include <unordered_map>
#include <lemon/time_measure.h>

/* namespace usage */
using namespace std;

void read_graph(std::string filename, std::shared_ptr<Graph> G, unordered_map<NODE_PAIR, double>& cost)
{
   /* Read data file */
   ifstream infile(filename);
   if (!infile.good())
   {
      cout << "Can't read file " << filename <<endl;
      exit(0);
   }

   int n, m;
   infile >> n >> m;
   cout << n << " " << m << endl;

   for(int k=0;k<n;k++)
   {
      NODE i;
      infile >> i;
      G->add_node(i);
   }

   for(int k=0;k<m;k++)
   {
      NODE i,j;
      double c;
      infile >> i >> j >> c;
      G->add_arc(i,j);
      cost[NODE_PAIR(i,j)] = c;
   }
}


int main(int argc, char** argv)
{
   srand(2);
   string filename = "";
   string st_filename ="";
   bool help = false;
   bool relax = false;
   vector<NODE> origins;
   vector<NODE> destinations;
   vector<pair<NODE,NODE>> st_pairs;
   vector<ElppForm> formulations;
   int pairs = 0;
   int K = 0;
   int timelimit = 1200; 
   int max_cuts = 1;
   double epsilon = 0.001;
   
   int o = 0;
   int tmp = -1;
   while ((o = getopt(argc, argv, "f:p:c:k:m:e:n:s:t:T:a::b::r::")) != -1)
      switch (o)
      {
         case 'f':
            filename = optarg;
            break;
         case 'p':
            st_filename = optarg;
            break;
         case 'c':
            tmp = atoi(optarg);
            if(tmp >= 0 && tmp < ElppFormulationName.size())
                formulations.push_back(ElppForm(tmp));
            break;
         case 'k':
            K = atoi(optarg);
            break;
         case 'T':
            timelimit = atoi(optarg);
            break;
         case 'n':
            pairs = atoi(optarg);
            break;
         case 'm':
            max_cuts = atoi(optarg);
            break;
         case 'e':
            epsilon = atof(optarg);
            break;
         case 's':
            origins.push_back(atoi(optarg));
            break;
         case 't':
            destinations.push_back(atoi(optarg));
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
      cout << "-p\tname of file with a list of (s,t) pairs." << endl;
      cout << "-s\torigin(s)\t(can specify more than one)" << endl;
      cout << "-t\tdestination(s)\t(can specify more than one)" << endl;
      cout << "-n\tuse a set of random origins and/or destinations of cardinality n and solve for all combinations [O(n^2) problems]" << endl;
      cout << "-T\ttime limit (default: 1200 seconds)." << endl;
      cout << "-c\tsubtour elimination constraints (can specify more than one):\n\t   0: None\n\t   1: MCF\n\t   2: Separation via Strong Components [default]\n";
      cout << "\t   3: Separation via Min-Cut\n\t   4: SF\n\t   5: RLT\n\t   6: MTZ\n\t   7: DL\n\t   8: DFJ(?)\n\t   9: MCFsep" << endl;
      cout << "-k\tstop after k problems have been solved" << endl;
      cout << "-e\tcut violation tolerance" << endl;
      cout << "-m\tmax number of cuts to be added in a callback. set -1 to add all violated inequalities (default: 1)" << endl;
      cout << "-r\tsolve LP relaxation" << endl;
      cout << endl << "Examples:" << endl;
      cout << "\t" << argv[0] << " -f data/toy.dat" << endl;
      cout << "\t\tcompute the Longest Path over the graph in toy.dat" << endl;
      cout << "\t\tfor all possible (s,t) pairs with the default formulation (SC)" << endl;
      cout << "\t" << argv[0] << " -f data/toy.dat -s 1 -t 5 -c 2" << endl;
      cout << "\t\tcompute the Longest Path over the graph in toy.dat" << endl;
      cout << "\t\tfrom node 1 to node 5 with formulation 2 (SC)" << endl;
      cout << "\t" << argv[0] << " -f data/toy.dat -s 1 -n 4 -c 1 -c 3" << endl;
      cout << "\t\tcompute the Longest Path over the graph in toy.dat" << endl;
      cout << "\t\tfrom node 1 to 4 other randomly selected nodes" << endl;
      cout << "\t\twith formulation 1 (MCF) and formulation 3 (Min-Cut)" << endl;
      cout << "\t" << argv[0] << " -f data/toy.dat -p data/toy.st" << endl;
      cout << "\t\tcompute the Longest Path over the graph in toy.dat" << endl;
      cout << "\t\tfor all the (s,t) pairs loaded from toy.st" << endl;
      cout << "\t\twith the default formulation (SC)" << endl;
      cout << "\t" << argv[0] << " -f data/toy.dat -n 4" << endl;
      cout << "\t\tcompute the Longest Path over the graph in toy.dat" << endl;
      cout << "\t\tfor all the (s,t) pairs connecting 4 randomly selected" << endl;
      cout << "\t\tnodes to 4 other randomly selected nodes" << endl;
      cout << "\t\twith formulation 1 (MCF) and formulation 3 (Min-Cut)" << endl;

      exit(0);
   }


   /* instance name */
   size_t lastdot = filename.find_first_of(".");
   size_t lastslash = filename.find_last_of("/");

   string stripped = filename.substr(0,lastdot);
   string instance_name;
   if (lastslash != string::npos)
      instance_name = stripped.substr(lastslash + 1, string::npos);
   else
      instance_name = stripped; //bah

   cout << "Selected formulations:"<<endl;
   for(ElppForm form : formulations)
      switch(form)
      {
         case NONE:
            cout << form << "- No subtour elimination" << endl;
            break;
         case MCF:
            cout << form << "- Using MCF formulation" << endl;
            break;
         case SC:
            cout << form << "- Using dynamic SEC separation via Strong Components" << endl;
            break;
         case MinCut:
            cout << form << "- Using dynamic SEC separation via Min-Cut" << endl;
            break;
         case SF:
            cout << form << "- Using SF formulation" << endl;
            break;
         case RLT:
            cout << form << "- Using RLT formulation" << endl;
            break;
         case MTZ:
            cout << form << "- Using MTZ formulation" << endl;
            break;
         case DL:
            cout << form << "- Using DL formulation" << endl;
            break;
         case DFJ:
            cout << form << "- Using DFJ formulation" << endl;
            break;
         case MCFsep:
            cout << form << "- Using MCF formulation with row-generation" << endl;
            break;
         default:
            cout << form << "- Unknown. Should not happen." << endl;
      }
   if(formulations.size() == 0)
   {
      cout << "No formulation selected. Using Strong Component (SC)." << endl;
      formulations.push_back(ElppForm::SC);
   }

   unordered_map<NODE_PAIR, double> cost;
   std::shared_ptr<Graph> G = std::make_shared<Graph>();
   read_graph(filename, G, cost);

   /* Add (s,t) pairs from .st file */
   if(st_pairs.size() == 0 && destinations.size() == 0)
   {
      ifstream st_file;
      st_file.open(st_filename);
      NODE s, t;
      while(st_file >> s >> t)
         if(find(G->nodes().begin(), G->nodes().end(), s) != G->nodes().end() && find(G->nodes().begin(), G->nodes().end(), t) != G->nodes().end())
            st_pairs.push_back(pair<NODE,NODE>(s,t));
      st_file.close();
   }


   mt19937 rnd_engine(12);
   
   /* Build set of origins*/
   bool ok = true;
   for(NODE s : origins)
      if (find(G->nodes().begin(), G->nodes().end(), s) == G->nodes().end())
      {
         cout << s << " not in graph!" << endl;
         ok = false;
         break;
      }
   if((st_pairs.size() == 0 && origins.size() == 0) || !ok)
   {
      origins = G->nodes();
      if(pairs>0 && pairs<int(origins.size()))
      {
         shuffle(origins.begin(), origins.end(), rnd_engine);
         origins.resize(pairs);
      }
   }

   /* Build set of destinations*/
   ok = true;
   for(NODE t : destinations)
      if (std::find(G->nodes().begin(), G->nodes().end(), t) == G->nodes().end())
      {
         cout << t << " not in graph!" << endl;
         ok = false;
         break;
      }
   if((st_pairs.size() == 0 && destinations.size() == 0) || !ok)
   {
      destinations = G->nodes();
      if(pairs>0 && pairs<int(destinations.size()))
      {
         shuffle(destinations.begin(), destinations.end(), rnd_engine);
         destinations.resize(pairs);
      }
   }
 

   /* Add all combinations to set of (s,t) pairs */ 
   for(NODE s : origins)
      for(NODE t : destinations)
         st_pairs.push_back(pair<NODE,NODE>(s,t)); 

   cout << "Set of " << st_pairs.size() << " s-t pairs." << endl;
   if(K == 0) 
      K = int(st_pairs.size());
   else
      cout << "[Solving at most " << K << ".]" << endl;

   int k=0;
   bool stop=false;
   /** Call solver for given origin-destination pairs, with the selected formulations */
   for(auto st : st_pairs)
   {
      NODE s = st.first;
      NODE t = st.second;
      if(s != t && !stop)
      {
         cout << "----" << s << " " << t << "----" << k << endl;
         for(ElppForm form : formulations)
         {
            IloEnv env;
            lemon::Timer time;
            if((/*MCF*/ form == ElppForm::MCF && G->num_nodes() >= 100) 
                  || (/*MCFsep*/ relax && form == ElppForm::MCFsep && G->num_nodes() >= 500))
            {
               cout << ElppFormulationName[form] << "\t: - -" << endl;
            }
            else{
               ElppSolver elpp_solver = ElppSolver(env, NODE_PAIR(s,t), G, form, relax, timelimit, epsilon, max_cuts);
               elpp_solver.update_problem(cost);
               time.start();
               if(relax)
                  elpp_solver.solveLP();
               else
                  elpp_solver.solve();
               time.stop();
               cout << "###" << ElppFormulationName[form] << "\t: ";
               if(elpp_solver.getStatus() == IloAlgorithm::Optimal)
                  cout << elpp_solver.getObjValue() << " " << time.userTime() << endl;
               else
               {
                  cout << "Status: " << elpp_solver.getStatus() << ", not solved to optimality! " << time.userTime() << endl;
               }
            }
            env.end();
         }
         if(++k == K)
            stop=true;
      }   
   }
   return 0;
}
