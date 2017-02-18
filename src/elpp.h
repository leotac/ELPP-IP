/**@file elpp.h
 * @brief Elementary Longest Path Problem solver
 * @author Leonardo Taccari 
 */

#ifndef __ELPP_H__
#define __ELPP_H__

#include "type.h"
#include <map>
#include <unordered_map>
#include <string>

using namespace std;
typedef unordered_map<NODE_PAIR,IloNumVar> ARC_VARS;
typedef map<NODE,vector<NODE>> ADJ_LIST;
typedef unordered_map<NODE_PAIR,int> INDEX_MAP;

class ElppSolverCplex
{
   public:
      ElppSolverCplex(){};
      ElppSolverCplex(IloEnv env,
            NODE_PAIR st_,
            const vector<NODE>& nodes_, 
            const vector<NODE_PAIR>&  arcs_,
            const map<NODE, vector<NODE>>& out_adj_list_,
            const map<NODE, vector<NODE>>& in_adj_list_,
            ElppForm formulation,
            bool relax,
            int timelimit,
            double epsilon,
            int max_cuts);

      void update_problem(
            const unordered_map<NODE_PAIR, IloNum>& obj_coeff 
            );
      void update_problem(
            const unordered_map<NODE_PAIR, IloNum>& obj_coeff, 
            const map<NODE_PAIR, IloNum>& lbs,
            const map<NODE_PAIR, IloNum>& ubs
            );

      void update_problem(
            const unordered_map<NODE_PAIR, IloNum>& obj_coeff, 
            const map<NODE_PAIR, IloNum>& lbs,
            const map<NODE_PAIR, IloNum>& ubs,
            const unordered_map<NODE_PAIR, IloNum>& lhs, 
            IloNum rhs,
            bool use_extra_con
            );

      void solve();
      void solveRoot();
      void solveLP();

      void clear();

      void append_info(string filename, string instance_name);

      void writeLP(string filename);

      void printInstance(
            string filename,
            const unordered_map<NODE_PAIR, IloNum>& obj_coeff, 
            const map<NODE_PAIR, IloNum>& lbs,
            const map<NODE_PAIR, IloNum>& ubs
            );

      IloAlgorithm::Status getStatus();
      IloNum getObjValue();
      IloNum getBestObjValue();
      IloNum getValue(NODE_PAIR arc);
      IloNum getValue(NODE node);
      bool isInteger();
      int pathLength();

   private:
      void build_problem_sec();
      void build_problem_mf(bool sep);
      void build_problem_sf();
      void build_problem_rlt();
      void build_problem_mtz();
      void build_problem_dl();

      unordered_map<NODE_PAIR, IloNumVar>          sigma_vars;
      unordered_map<NODE_PAIR,IloNumVar>           u_var;
      unordered_map<TRIPLET,IloNumVar>           qq_var;
      unordered_map<NODE,IloNumVar>                zz_var;
      unordered_map<NODE,IloNumVar>                p_var;

      IloNumVarArray x_vararray;
      unordered_map<NODE_PAIR, int> index;

      IloModel                model;
      IloCplex                cplex;
      IloRange               extra_con;
      IloObjective            objective;
   

      int timelimit;
      ElppForm formulation;
      bool relax;
      double elapsed_time;
      double elapsed_ticks;
      int ncuts;
      double tol;
      int max_cuts;

      //Graph structure
      NODE_PAIR st;
      vector<NODE>                       nodes;         /* array of nodes */
      vector<NODE_PAIR>                  arcs;          /* arcs */
      map<NODE, vector<NODE>>            out_adj_list;      /* adjacency list */
      map<NODE, vector<NODE>>            in_adj_list;      /* adjacency list */

};


#endif
