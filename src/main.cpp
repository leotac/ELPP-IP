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
#include <algorithm>
#include <random>
#include <map>
#include <unordered_map>
#include <lemon/time_measure.h>

/* namespace usage */
using namespace std;


int main(int argc, char** argv)
{
   srand(2);
   string filename = "";
   string infodir = "";
   string st_filename ="";
   bool bounds = false;
   //int formulation = 1;
   bool help = false;
   bool relax = false;
   bool all = false;
   vector<NODE> origins;
   vector<NODE> destinations;
   vector<pair<NODE,NODE>> st_pairs;
   vector<ElppForm> formulations;
   int pairs = 0;
   int K = 0;
   int timelimit = 1200; 
   int max_cuts = 1;
   double epsilon = 0.001;
   
   bool append = false;

   int o = 0;
   while ((o = getopt(argc, argv, "f:p:i:c:k:m:e:n:s:t:T:a::b::r::x::")) != -1)
      switch (o)
      {
         case 'f':
            filename = optarg;
            break;
         case 'p':
            st_filename = optarg;
            break;
         case 'i':
            infodir = optarg;
            break;
         case 'a':
            all = true;
            break;
         case 'c':
            formulations.push_back((ElppForm) atoi(optarg));
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
         case 'b':
            bounds = true; 
            break;
         case 'r':
            relax = true; 
            break;
         case 'x':
            append = true; 
            break;
         default:
            help = true;
            break;
      }

   if ( filename == "" || help == true )
   {
      cout << "Usage: " << argv[0] << " -f <file_name> [other options]" << endl;
      cout << "If no origin(s) and/or destination(s) are specified, all combinations are solved." << endl;
      cout << "Options:" << endl;
      cout << "-f\t(mandatory) name of data file." << endl;
      cout << "-p\t(optional) name of file with a list of (s,t) pairs. If not provided, try to infer from the data file." << endl;
      cout << "-i\tdirectory of the info file (default: none)." << endl;
      cout << "-T\ttime limit (default: 1200 seconds)." << endl;
      cout << "-c\tsubtour elimination constraints:\n\t   0: None\n\t   1: MCF [default]\n\t   2: Separation via Strong Components\n";
      cout << "\t   3: Separation via Min-Cut\n\t   4: SF\n\t   5: RLT\n\t   6: MTZ\n\t   7: DL\n\t   8: DFJ(?)\n\t   9: MCFsep" << endl;
      cout << "-k\tstop after k problems have been solved" << endl;
      cout << "-e\tviolation tolerance" << endl;
      cout << "-m\tmax number of cuts to be added in a callback. set -1 to add all violated inequalities (default: 1)" << endl;
      cout << "-b\tread arc variables bounds from data file" << endl;
      cout << "-s\torigin(s)\t(can specify more than one)" << endl;
      cout << "-t\tdestination(s)\t(can specify more than one)" << endl;
      cout << "-a\ttry all-to-all (s,t) pairs" << endl;
      cout << "-n\tuse a set of random origins and/or destinations of cardinality n and solve for all combinations [O(n^2) problems]" << endl;
      cout << "-r\tsolve LP relaxation" << endl;
      cout << "-x\tappend to existing info file" << endl;
      exit(0);
   }

   unordered_map<NODE_PAIR, double> cost;
   map<NODE, vector<NODE>> out_adj_list;
   map<NODE, vector<NODE>> in_adj_list;

   /* instance name */
   size_t lastdot = filename.find_first_of(".");
   size_t lastslash = filename.find_last_of("/");

   string stripped = filename.substr(0,lastdot);
   string instance_name;
   if (lastslash != string::npos)
      instance_name = stripped.substr(lastslash + 1, string::npos);
   else
      instance_name = stripped; //bah

   /* try to infer st filename */
   if(st_filename == "")
      st_filename = stripped + ".st";


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
            cout << form << "- Unknown." << endl;
      }
   if(formulations.size() == 0)
   {
      cout << "None selected. Using MCF." << endl;
      formulations.push_back(ElppForm::MCF);
   }

   /* Read data file */
   ifstream infile(filename);
   if (!infile.good())
   {
      cout << "Can't read file " << filename <<endl;
      exit(0);
   }

   int n, m;
   infile >> n >> m;
   //cout << n << " " << m << endl;
   
   vector<NODE> nodes;
   nodes.reserve(n);
   vector<NODE_PAIR> arcs;
   arcs.reserve(m);
   map<NODE_PAIR, double> ubs, lbs;
   
   for(int k=0;k<n;k++)
   {
      NODE i;
      infile >> i;
      //cout << i << endl;
      nodes.push_back(i);
      out_adj_list[i] = vector<NODE>();
      in_adj_list[i] = vector<NODE>();
   }

   for(int k=0;k<m;k++)
   {
      NODE i,j;
      double c;
      infile >> i >> j >> c;
      //cout << i << " "<< j << " " << c <<endl;
      out_adj_list[i].push_back(j);
      in_adj_list[j].push_back(i);
      cost[NODE_PAIR(i,j)] = c;
      arcs.push_back(NODE_PAIR(i,j));
   }

   if(bounds) 
      for(int k=0;k<m;k++)
      {
         NODE i,j;
         double lb = 0;
         double ub = 1;
         infile >> i >> j >> lb >> ub;
         //cout << i << " "<< j << " " << lb << " " << ub <<endl;
         if(cost.count(NODE_PAIR(i,j))==0) 
         {
            cout << "Unexpected end of file or non-existing arc: ( " << i << ", " << j << " )" <<endl;
            exit(-1);
         }
         lbs[NODE_PAIR(i,j)] = lb;
         ubs[NODE_PAIR(i,j)] = ub;
      }

   /* Add (s,t) pairs from .st file */
   if(!all && st_pairs.size() == 0 && destinations.size() == 0)
   {
      ifstream st_file;
      st_file.open(st_filename);
      NODE s, t;
      while(st_file >> s >> t)
         if(find(nodes.begin(), nodes.end(), s) != nodes.end() && find(nodes.begin(), nodes.end(), t) != nodes.end())
            st_pairs.push_back(pair<NODE,NODE>(s,t));
      st_file.close();
   }


   mt19937 rnd_engine(12);
   
   /* Build set of origins*/
   bool ok = true;
   for(NODE s : origins)
      if (find(nodes.begin(), nodes.end(), s) == nodes.end())
      {
         cout << s << " not in graph!" << endl;
         ok = false;
         break;
      }
   if((st_pairs.size() == 0 && origins.size() == 0) || !ok)
   {
      origins = nodes;
      if(pairs>0 && pairs<(int)origins.size())
      {
         shuffle(origins.begin(), origins.end(), rnd_engine);
         origins.resize(pairs);
      }
   }

   /* Build set of destinations*/
   ok = true;
   for(NODE t : destinations)
      if (std::find(nodes.begin(), nodes.end(), t) == nodes.end())
      {
         cout << t << " not in graph!" << endl;
         ok = false;
         break;
      }
   if((st_pairs.size() == 0 && destinations.size() == 0) || !ok)
   {
      destinations = nodes;
      if(pairs>0 && pairs<(int)destinations.size())
      {
         shuffle(destinations.begin(), destinations.end(), rnd_engine);
         destinations.resize(pairs);
      }
   }
 
   
   /* Add all combinations to set of (s,t) pairs */ 
   for(NODE s : origins)
      for(NODE t : destinations)
         st_pairs.push_back(pair<NODE,NODE>(s,t)); 



   /* */
   string info_filename = "";
   if(infodir != "")
   {
      info_filename = infodir + "/" + instance_name + ".nfo";
      if(relax)
         info_filename += ".relax";

      ofstream infofile;
      if(!ifstream(info_filename) || !append)
      {
         infofile.open(info_filename);
         infofile << "inst\t";
         if(relax)
            infofile << "s\tt\tform\topt\ttime\tticks\tsolveT\tcutT\tnodes\tcuts\tint" << endl;
         else
            infofile << "s\tt\tform\topt\ttime\tticks\tsolveT\tcutT\tnodes\tcuts\ttol\tmaxc" << endl;
         infofile.close();
      }
   }
   
   cout << "Set of " << st_pairs.size() << " s-t pairs." << endl;
   if(K == 0) 
      K = (int) st_pairs.size();
   else
      cout << "[Solving at most " << K << ".]" << endl;

   int k=0;
   bool stop=false;
   /** Call solver for given origin-destination pairs */
   for(auto st : st_pairs)
   {
      NODE s = st.first;
      NODE t = st.second;
      if(s != t && !stop)
      {
         cout << "----" << s << " " << t << "----" << k << endl;
         ++k;
         for(ElppForm form : formulations)
         {
            IloEnv env;
            lemon::Timer time;
            if((/*MCF*/ form == ElppForm::MCF && n >= 100) 
                  || (/*MCFsep*/ relax && form == ElppForm::MCFsep && n >= 500) 
                  //|| (/*MCFsep*/ !relax && form == ElppForm::MCFsep && n >= 500) 
                  
                  )
            {
               cout << which(form) << "\t: - -" << endl;
               if(info_filename != "")
               {
                  ofstream infofile(info_filename, std::ios_base::app);
                  infofile << instance_name << "\t" << s << "\t" << t << "\t" << which(form) << "\t-\t-\t-\t-\t-\t-\t-";
                  if(relax)
                     infofile << "\t0" << endl;
                  else
                     infofile << "\t-\t-" << endl;
                  infofile.close();
               }
            }
            else{

               ElppSolverCplex elpp_solver = ElppSolverCplex(env, NODE_PAIR(s,t), nodes, arcs, out_adj_list, in_adj_list, form, relax, timelimit, epsilon, max_cuts);
               if(bounds)
                  elpp_solver.update_problem(cost, lbs, ubs);
               else
                  elpp_solver.update_problem(cost);
               time.start();
               if(relax)
                  elpp_solver.solveLP();
               else
                  elpp_solver.solve();
               time.stop();
               if(info_filename != "")
                  elpp_solver.append_info(info_filename, instance_name);
               cout << "###" << which(form) << "\t: ";
               if(elpp_solver.getStatus() == IloAlgorithm::Optimal)
                  cout << elpp_solver.getObjValue() << " " << time.userTime() << endl;
               else
               {
                  cout << "Status: " << elpp_solver.getStatus() << ", not solved to optimality! " << time.userTime() << endl;
               }
            }
            env.end();
         }
         if(k == K)
            stop=true;


      }   
   }
   return 0;
}
