/**@file elpp.cpp
 * @brief Elementary Longest Path problem solver
 * @author Leonardo Taccari 
 */

#include "elpp.h"
#include "unionfind.h"
#include "separation.h"
#include "callbacks.h"
#include <iostream>
#include <fstream>
#include <algorithm> 
#include <lemon/time_measure.h>
#define LOG if(false) cerr
#define TOL 0.001

using namespace std;
using namespace lemon;

/*******************************************************/
/************* Build the model *************************/
/*******************************************************/

/** Constructor */
ElppSolver::ElppSolver(
      IloEnv env,
      NODE_PAIR st_,
      const vector<NODE>& nodes_, 
      const vector<NODE_PAIR>&  arcs_,
      const map<NODE, vector<NODE>>& out_adj_list_,
      const map<NODE, vector<NODE>>& in_adj_list_,
      ElppForm formulation_,
      bool relax_,
      int timelimit_,
      double epsilon,
      int max_cuts_)
{
   // Initialize graph structures
   st = st_;
   nodes = nodes_;
   arcs = arcs_;
   out_adj_list = out_adj_list_; 
   in_adj_list = in_adj_list_; 
   // Initialize Cplex structures
   model = IloModel(env);
   extra_con = IloRange();
   objective = IloObjective();

   // Settings
   formulation = formulation_;
   relax = relax_;
   timelimit = timelimit_;
   tol = epsilon;
   max_cuts = max_cuts_;

   /* Add x variables */
   // sigma[i,j] 
   // Whether arc (i,j) is in the longest path 
   char var_name[255];
   x_vararray = IloNumVarArray(env);
   int idx = 0;
   for (auto arc : arcs)
   {   
      IloNumVar var;
      snprintf(var_name, 255, "sigma_%d_%d", arc.first,arc.second);
      if(relax)
         var = IloNumVar(env, 0,1, IloNumVar::Float, var_name);
      else
         var = IloNumVar(env, 0,1, IloNumVar::Bool, var_name);
      sigma_vars[arc] = var;
      x_vararray.add(var);
      model.add(var);
      index[arc] = idx;
      ++idx;
   }

   LOG << "Building problem...";
   //LOG << 4*sizeof(IloNumVar::ImplClass)*((nodes.size()-1)*arcs.size() + arcs.size() + nodes.size())<< " Bytes" << endl;
   switch(formulation)
   {
      case MCF:
         build_problem_mf(false);
         break;
      case MCFsep:
         build_problem_mf(true);
         break;
      case NONE:
      case SC:
      case MinCut:
      case DFJ:
         build_problem_sec();
         break;
      case SF:
         build_problem_sf();
         break;
      case RLT:
         build_problem_rlt();
         break;
      case MTZ:
         build_problem_mtz();
         break;
      case DL:
         build_problem_dl();
         break;
      default:
         build_problem_mf(false);
   }
   LOG << "done." << endl;

   // Cplex settings
   // Need to load the model just once: it is automatically updated
   cplex = IloCplex(model);
   //cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);

   //cplex.setOut(env.getNullStream());
   cplex.setParam(IloCplex::MIPDisplay, 2);
   //cplex.setParam(IloCplex::SimDisplay, 2);

   // Disable MIPstarts in presence of lazy constraints. 
   // Needed when solving the same problem more than once (otherwise, problems with lazy cons).
   // Currently disabled for all formulations, in order not to change computational results when re-optimizing
   if(formulation > 0)
      cplex.setParam(IloCplex::AdvInd, 0);

   // Disable presolve? (works well in pricing)
   //cplex.setParam(IloCplex::PreInd, false);

   cplex.setParam(IloCplex::EpGap, 1e-09); //default:1e-04, 0.1%
   cplex.setParam(IloCplex::Threads, 1);
   cplex.setParam(IloCplex::TreLim, 2048);
   cplex.setParam(IloCplex::WorkMem, 4096);
   cplex.setParam(IloCplex::TiLim, timelimit);
   cplex.setParam(IloCplex::EpInt, 1e-06);

   //double tol = TOL;
   //int max_cuts = 1;
   /* User cuts / lazy callbacks*/
   switch(formulation)
   {
      case NONE: /* No subtour elimination constraints */
      case MCF: /* MCF */
      case SF:
      case RLT:
      case MTZ:
      case DL:
         break;
      case SC: /* Separation with Strong Component*/
      case MinCut: /* Separation with Min-Cut*/
      case MCFsep: /* Separation with Min-Cut for MCF*/
      case DFJ: /* Separation with SC for DFJ*/
         cplex.use(StrongComponentLazyCallback(env,nodes, arcs, st, sigma_vars, qq_var, zz_var, out_adj_list, in_adj_list, x_vararray, index, tol, max_cuts, formulation));
         cplex.use(ElppCutCallback(env,nodes, arcs, st, sigma_vars, qq_var, zz_var, out_adj_list, in_adj_list, x_vararray, index, tol, max_cuts, formulation));
         break;
      default:
         break;
   }

   // Set a solution limit?
   //cplex.setParam(IloCplex::IntSolLim, 1);
}


void ElppSolver::build_problem_sec()
{
   /***************/ 
   /* Build model */
   /***************/ 
   NODE s = st.first;
   NODE t = st.second;

   IloEnv env = model.getEnv();

   /*******************/
   /* Add constraints */
   /*******************/

   // Flow conservation constraints (on sigma)
   for (NODE node : nodes)
   {
      IloRange con;

      IloNum coeff;
      if( node == s)
         coeff = 1.0;
      else if ( node == t)
         coeff = -1.0;
      else
         coeff = 0.0;

      con = IloRange(env, coeff, coeff);

      for(auto j : out_adj_list.at(node))
         con.setLinearCoef(sigma_vars[NODE_PAIR(node,j)],1.0);
      for(auto j : in_adj_list.at(node))
         con.setLinearCoef(sigma_vars[NODE_PAIR(j,node)],-1.0);
      model.add(con);
   }

   // Incoming s = 0 
   IloRange con = IloRange(env, 0, 0);
   for(auto j : in_adj_list.at(s)) 
      con.setLinearCoef(sigma_vars[NODE_PAIR(j,s)],1.0);
   model.add(con);

   // Outgoing t = 0
   con = IloRange(env, 0, 0);
   for(auto j : out_adj_list.at(t)) 
      con.setLinearCoef(sigma_vars[NODE_PAIR(t,j)],1.0);
   model.add(con);

   // Outgoing degree constraints (on sigma)
   for (NODE node : nodes)
     {
     IloRange con;
     con = IloRange(env, -IloInfinity, 1.0);

     for(auto j : out_adj_list.at(node)) 
        con.setLinearCoef(sigma_vars[NODE_PAIR(node,j)],1.0);
     model.add(con);
     }

}

void ElppSolver::build_problem_mf(bool sep)
{
   /***************/ 
   /* Build model */
   /***************/ 
   NODE s = st.first;
   NODE t = st.second;


   LOG << arcs.size() << endl;
   LOG << nodes.size() << endl;
   LOG << "Number of variables to build: " << arcs.size() + arcs.size()*(nodes.size()-1) +  (nodes.size()-1)<< endl;
   LOG << "Number of constraints to add each time: ~" << arcs.size() + nodes.size() << endl;

   IloEnv env = model.getEnv();
   char var_name[255];

   // q[i,j,h]
   // (Imaginary) flow assigned to node h on arc (i,j)
   //unordered_map<TRIPLET,IloNumVar> qq_var;
   for (auto arc : arcs)
   {   
      for (NODE h : nodes)
         if ( h != s )
         {
            snprintf(var_name, 255, "q_%d_%d_%d", arc.first,arc.second,h);
            IloNumVar var;
            var = IloNumVar(env,var_name);
            var.setBounds(0, 1);
            qq_var[TRIPLET(h,arc)] = var;
            model.add(var); //Add all variables even when they appear in no constraints yet (MCFsep)
         }
   }

   // z[h]
   // Whether node h is used
   //unordered_map<NODE,IloNumVar> zz_var;
   for (NODE h : nodes)
      if ( h != s) // Yes!
      {
         snprintf(var_name, 255, "z_%d", h);
         IloNumVar var;
         var = IloNumVar(env, var_name);
         var.setBounds(0, 1);
         zz_var[h] = var;
         model.add(var); //Add all variables even when they appear in no constraints yet (MCFsep)
      }

   /*******************/
   /* Add constraints */
   /*******************/

   // Flow conservation constraints (on sigma)
   for (NODE node : nodes)
   {
      IloRange con;

      IloNum coeff;
      if( node == s)
         coeff = 1.0;
      else if ( node == t)
         coeff = -1.0;
      else
         coeff = 0.0;

      con = IloRange(env, coeff, coeff);

      for(auto j : out_adj_list.at(node))
         con.setLinearCoef(sigma_vars[NODE_PAIR(node,j)],1.0);
      for(auto j : in_adj_list.at(node))
         con.setLinearCoef(sigma_vars[NODE_PAIR(j,node)],-1.0);
      model.add(con);
   }

   // Outgoing degree constraints (on sigma)
   for (NODE node : nodes)
   {
      IloRange con;
      con = IloRange(env, -IloInfinity, 1.0);

      for(auto j : out_adj_list.at(node)) 
         con.setLinearCoef(sigma_vars[NODE_PAIR(node,j)],1.0);
      model.add(con);
   }

   // Incoming s = 0 
   IloRange con = IloRange(env, 0, 0);
   for(auto j : in_adj_list.at(s)) 
      con.setLinearCoef(sigma_vars[NODE_PAIR(j,s)],1.0);
   model.add(con);

   // Outgoing t = 0
   con = IloRange(env, 0, 0);
   for(auto j : out_adj_list.at(t)) 
      con.setLinearCoef(sigma_vars[NODE_PAIR(t,j)],1.0);
   model.add(con);

   if(!sep)
   {
      // Ham1
      // qq[h,i,j] <= sigma[i,j] (s,t)
      for (auto arc : arcs)
      {   
         for (NODE h : nodes)
         {
            if(h != s)
            { 
               IloRange con;
               con = IloRange(env, -IloInfinity, 0.0);

               con.setLinearCoef(qq_var[TRIPLET(h,arc)],1.0);
               con.setLinearCoef(sigma_vars[arc],-1.0);
               model.add(con);
            }
         }
      }

      // Ham2
      // sum{(s,j) in A}qq[h,s,j] = z[h] (s,t)
      for (NODE h : nodes)
      {
         if ( h != s)
         {
            IloRange con;
            con = IloRange(env, 0.0, 0.0);

            con.setLinearCoef(zz_var[h],-1.0);

            for(auto j : out_adj_list.at(s)) 
               con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(s,j))],1.0);

            model.add(con);
         }
      }

      // Ham3
      // sum{(i,s) in A}qq[h,i,s] = 0
      for (NODE h : nodes)
      {
         if ( h != s)
         {
            IloRange con;
            con = IloRange(env, 0.0, 0.0);

            for(auto j : in_adj_list.at(s)) 
               con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(j,s))],1.0);

            model.add(con);
         }
      }

      // Ham4
      // sum{(i,h) in A}q[h,i,h] = z[h]
      for (NODE h : nodes)
      {
         if ( h != s)
         {
            IloRange con;
            con = IloRange(env, 0.0, 0.0);

            con.setLinearCoef(zz_var[h],-1.0);

            for(auto j : in_adj_list.at(h)) 
               con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(j,h))],1.0);

            model.add(con);
         }
      }  

      // Ham5
      // sum{(h,j) in A}q[h,h,j] = 0 
      for (NODE h : nodes)
      {
         if ( h != s)
         {
            IloRange con;
            con = IloRange(env, 0.0, 0.0);

            for(auto j : out_adj_list.at(h)) 
               con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(h,j))],1.0);

            model.add(con);
         }
      }

      // Ham6
      // sum{(i,j) in A}qq[h,i,j] - sum{(j,i) in A}qq[h,j,i] = 0
      for (NODE node : nodes)
      {
         for (NODE h : nodes)
         {
            if ( h != s && node != h && node != s)
            {
               IloRange con;
               con = IloRange(env, 0.0, 0.0);
               for(auto j : out_adj_list.at(node))
                  con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(node,j))],-1.0);
               for(auto j : in_adj_list.at(node))
                  con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(j,node))],1.0);
               model.add(con);
            }
         }  
      }

      // Ham7
      // sum{(i,h) in A}sigma[i,h] = z[h]
      // if node h is reached, z[h] = 1
      for (NODE h : nodes)
      {
         if ( h != s)
         {
            IloRange con;
            con = IloRange(env, 0.0, 0.0);

            con.setLinearCoef(zz_var[h],-1.0);

            for(auto j : in_adj_list.at(h)) 
               con.setLinearCoef(sigma_vars[NODE_PAIR(j,h)],1.0);
            model.add(con);
         }
      }
   }
}

/** Single flow formulation */
void ElppSolver::build_problem_sf()
{
   /***************/ 
   /* Build model */
   /***************/ 
   NODE s = st.first;
   NODE t = st.second;

   IloEnv env = model.getEnv();
   /* Add variables */
   char var_name[255];

   // u[i,j] 
   // (Imaginary) flow assigned on arc (i,j)
   //unordered_map<NODE_PAIR,IloNumVar> u_var;
   for (auto arc : arcs)
   {   
      snprintf(var_name, 255, "q_%d_%d", arc.first,arc.second);
      IloNumVar var;
      var = IloNumVar(env,var_name);
      var.setBounds(0, (int)nodes.size()-1);
      u_var[arc] = var;
   }

   // z[h]
   // Whether node h is used
   //unordered_map<NODE,IloNumVar> zz_var;
   for (NODE h : nodes)
      if ( h != s )
      {
         snprintf(var_name, 255, "z_%d", h);
         IloNumVar var;
         var = IloNumVar(env, var_name);
         var.setBounds(0, 1);
         zz_var[h] = var;
      }

   /*******************/
   /* Add constraints */
   /*******************/
   // Flow conservation constraints (on sigma)
   for (NODE node : nodes)
   {
      IloRange con;

      IloNum coeff;
      if( node == s)
         coeff = 1.0;
      else if ( node == t)
         coeff = -1.0;
      else
         coeff = 0.0;

      con = IloRange(env, coeff, coeff);

      for(auto j : out_adj_list.at(node))
         con.setLinearCoef(sigma_vars[NODE_PAIR(node,j)],1.0);
      for(auto j : in_adj_list.at(node))
         con.setLinearCoef(sigma_vars[NODE_PAIR(j,node)],-1.0);
      model.add(con);
   }

   // Outgoing degree constraints (on sigma)
   for (NODE node : nodes)
   {
      IloRange con;
      con = IloRange(env, -IloInfinity, 1.0);

      for(auto j : out_adj_list.at(node)) 
         con.setLinearCoef(sigma_vars[NODE_PAIR(node,j)],1.0);
      model.add(con);
   }

   // Incoming s = 0 
   IloRange con = IloRange(env, 0, 0);
   for(auto j : in_adj_list.at(s)) 
      con.setLinearCoef(sigma_vars[NODE_PAIR(j,s)],1.0);
   model.add(con);

   // Outgoing t = 0
   con = IloRange(env, 0, 0);
   for(auto j : out_adj_list.at(t)) 
      con.setLinearCoef(sigma_vars[NODE_PAIR(t,j)],1.0);
   model.add(con);

   // SF1
   // q[i,j] <= (n-1)*sigma[i,j] (s,t)
   for (auto arc : arcs)
   {   
      IloRange con;
      con = IloRange(env, -IloInfinity, 0.0);

      con.setLinearCoef(u_var[arc],1.0);
      con.setLinearCoef(sigma_vars[arc],-((int)nodes.size()-1));
      model.add(con);
   }

   // SF2
   // sum{(s,j) in A}qq[s,j] = sum{h in V: h!=s}z[h] (s,t)
   con = IloRange(env, 0.0, 0.0);

   for (NODE h : nodes)
      if ( h != s)
         con.setLinearCoef(zz_var[h],-1.0);

   for(auto j : out_adj_list.at(s)) 
      con.setLinearCoef(u_var[NODE_PAIR(s,j)],1.0);

   model.add(con);

   // SF2-bis
   // sum{(i,s) in A}qq[i,s] = 0
   con = IloRange(env, 0.0, 0.0);
   for(auto j : in_adj_list.at(s)) 
      con.setLinearCoef(u_var[NODE_PAIR(j,s)],1.0);
   model.add(con);

   // SF3
   // sum{(i,h) in A}q[i,h] - sum{(h,j) in A}q[h,j] = z[h]
   for (NODE h : nodes)
   {
      if ( h != s)
      {
         IloRange con;
         con = IloRange(env, 0.0, 0.0);

         con.setLinearCoef(zz_var[h],-1.0);

         for(auto j : in_adj_list.at(h))
            con.setLinearCoef(u_var[NODE_PAIR(j,h)],1.0);  // incoming
         for(auto j : out_adj_list.at(h))
            con.setLinearCoef(u_var[NODE_PAIR(h,j)],-1.0); // outgoing

         model.add(con);
      }
   }  

   // SF4
   // sum{(i,h) in A}sigma[i,h] = z[h]
   // if node h is reached, z[h] = 1
   for (NODE h : nodes)
   {
      if ( h != s)
      {
         IloRange con;
         con = IloRange(env, 0.0, 0.0);

         con.setLinearCoef(zz_var[h],-1.0);

         for(auto j : in_adj_list.at(h)) 
            con.setLinearCoef(sigma_vars[NODE_PAIR(j,h)],1.0);
         model.add(con);
      }
   }
}


/** RLT formulation */
void ElppSolver::build_problem_rlt()
{
   /***************/ 
   /* Build model */
   /***************/ 
   NODE s = st.first;
   NODE t = st.second;

   IloEnv env = model.getEnv();

   /* Add variables */
   char var_name[255];

   // alpha[i,j] 
   unordered_map<NODE_PAIR,IloNumVar> alpha_var;
   for (auto arc : arcs)
   {   
      snprintf(var_name, 255, "alpha_%d_%d", arc.first,arc.second);
      IloNumVar var;
      var = IloNumVar(env,var_name);
      var.setBounds(0, IloInfinity);
      alpha_var[arc] = var;
   }

   // beta[i,j] 
   unordered_map<NODE_PAIR,IloNumVar> beta_var;
   for (auto arc : arcs)
   {   
      snprintf(var_name, 255, "beta_%d_%d", arc.first,arc.second);
      IloNumVar var;
      var = IloNumVar(env,var_name);
      var.setBounds(0, IloInfinity);
      beta_var[arc] = var;
   }


   /*******************/
   /* Add constraints */
   /*******************/
   // Flow conservation constraints (on sigma)
   for (NODE node : nodes)
   {
      IloRange con;

      IloNum coeff;
      if( node == s)
         coeff = 1.0;
      else if ( node == t)
         coeff = -1.0;
      else
         coeff = 0.0;

      con = IloRange(env, coeff, coeff);

      for(auto j : out_adj_list.at(node))
         con.setLinearCoef(sigma_vars[NODE_PAIR(node,j)],1.0);
      for(auto j : in_adj_list.at(node))
         con.setLinearCoef(sigma_vars[NODE_PAIR(j,node)],-1.0);
      model.add(con);
   }

   // Outgoing degree constraints (on sigma)
   for (NODE node : nodes)
   {
      IloRange con;
      con = IloRange(env, -IloInfinity, 1.0);

      for(auto j : out_adj_list.at(node)) 
         con.setLinearCoef(sigma_vars[NODE_PAIR(node,j)],1.0);
      model.add(con);
   }

   // Incoming s = 0 
   IloRange con = IloRange(env, 0, 0);
   for(auto j : in_adj_list.at(s)) 
      con.setLinearCoef(sigma_vars[NODE_PAIR(j,s)],1.0);
   model.add(con);

   // Outgoing t = 0
   con = IloRange(env, 0, 0);
   for(auto j : out_adj_list.at(t)) 
      con.setLinearCoef(sigma_vars[NODE_PAIR(t,j)],1.0);
   model.add(con);


   // RLT1
   // alpha[i,j] = beta[i,j] + x[i,j]
   for (auto arc : arcs)
   {   
      IloRange con;
      con = IloRange(env, 0.0, 0.0);

      con.setLinearCoef(alpha_var[arc], 1);
      con.setLinearCoef(beta_var[arc], -1);
      con.setLinearCoef(sigma_vars[arc], -1);
      model.add(con);
   }

   // RLT2
   // x[s,j] + sum{(i,j) in A: i<>s}alpha[i,j] - sum{(j,i) in A}beta[j,i] = 0;
   for(auto j : out_adj_list.at(s))
      if(j!=t)
      {
         IloRange con;
         con = IloRange(env, 0.0, 0.0);
         con.setLinearCoef(sigma_vars[NODE_PAIR(s,j)], 1);

         for(auto i : in_adj_list.at(j))
            if(i!=s)
               con.setLinearCoef(alpha_var[NODE_PAIR(i,j)], 1);

         for(auto i : out_adj_list.at(j))
            con.setLinearCoef(beta_var[NODE_PAIR(j,i)], -1);

         model.add(con);
      }

   // RLT3
   // sum{(i,j) in A}alpha[i,j] - sum{(i,j) in A}beta[i,j]=0
   for(auto j : nodes)
      if( j!=s && j!=t && (find(out_adj_list.at(s).begin(),out_adj_list.at(s).end(),j) == out_adj_list.at(s).end()))
      {
         IloRange con;
         con = IloRange(env, 0.0, 0.0);
         for(auto i : in_adj_list.at(j))
            con.setLinearCoef(alpha_var[NODE_PAIR(i,j)], 1);

         for(auto i : out_adj_list.at(j))
            con.setLinearCoef(beta_var[NODE_PAIR(j,i)], -1);

         model.add(con);
      }

//   // RLT4
//   // x[i,j] <= alpha[i,j]
//   for (auto arc : arcs)
//      if(arc.first!=s && arc.second!=s)
//      {   
//         IloRange con;
//         con = IloRange(env, -IloInfinity, 0.0);
//         con.setLinearCoef(sigma_vars[arc], 1);
//         con.setLinearCoef(alpha_var[arc], -1);
//         model.add(con);
//      }

   // RLT5
   // alpha[i,j] <= (n-1)*x[i,j]
   for (auto arc : arcs)
      if(arc.first!=s)
      {   
         IloRange con;
         con = IloRange(env, -IloInfinity, 0.0);
         con.setLinearCoef(alpha_var[arc], 1);
         con.setLinearCoef(sigma_vars[arc],1.0-(int)nodes.size());
         model.add(con);
      }

   // RLT6
   // x[i,j] <= beta[i,j]
   for (auto arc : arcs)
      if(arc.first!=s && arc.second!=s)
      {   
         IloRange con;
         con = IloRange(env, -IloInfinity, 0.0);
         con.setLinearCoef(sigma_vars[arc], 1);
         con.setLinearCoef(beta_var[arc], -1);
         model.add(con);
      }

//   // RLT7
//   // beta[i,j] <= (n-1)*x[i,j]
//   for (auto arc : arcs)
//      if(arc.first!=s)
//      {   
//         IloRange con;
//         con = IloRange(env, -IloInfinity, 0.0);
//         con.setLinearCoef(beta_var[arc], 1);
//         con.setLinearCoef(sigma_vars[arc],1.0-(int)nodes.size());
//         model.add(con);
//      }
//
//   // RLT8
//   // alpha[i,j] >= beta[i,j]
//   for (auto arc : arcs)
//   {   
//      IloRange con;
//      con = IloRange(env, -IloInfinity, 0.0);
//      con.setLinearCoef(beta_var[arc], 1);
//      con.setLinearCoef(alpha_var[arc],-1);
//      model.add(con);
//   }

}  /* RLT */

/** MTZ formulation */
void ElppSolver::build_problem_mtz()
{
   /***************/ 
   /* Build model */
   /***************/ 
   NODE s = st.first;
   NODE t = st.second;

   IloEnv env = model.getEnv();
   /* Add variables */
   char var_name[255];

   // p[i] 
   // Label of node i 
   for (auto i : nodes)
   {   
      snprintf(var_name, 255, "p_%d", i);
      IloNumVar var;
      var = IloNumVar(env,var_name);
      var.setBounds(0, (int)nodes.size());
      p_var[i] = var;
   }
   /*******************/
   /* Add constraints */
   /*******************/
   // Flow conservation constraints (on sigma)
   for (NODE node : nodes)
   {
      IloRange con;

      IloNum coeff;
      if( node == s)
         coeff = 1.0;
      else if ( node == t)
         coeff = -1.0;
      else
         coeff = 0.0;

      con = IloRange(env, coeff, coeff);

      for(auto j : out_adj_list.at(node))
         con.setLinearCoef(sigma_vars[NODE_PAIR(node,j)],1.0);
      for(auto j : in_adj_list.at(node))
         con.setLinearCoef(sigma_vars[NODE_PAIR(j,node)],-1.0);
      model.add(con);
   }

   // Outgoing degree constraints (on sigma)
   for (NODE node : nodes)
   {
      IloRange con;
      con = IloRange(env, -IloInfinity, 1.0);

      for(auto j : out_adj_list.at(node)) 
         con.setLinearCoef(sigma_vars[NODE_PAIR(node,j)],1.0);
      model.add(con);
   }

   // Incoming s = 0 
   IloRange con = IloRange(env, 0, 0);
   for(auto j : in_adj_list.at(s)) 
      con.setLinearCoef(sigma_vars[NODE_PAIR(j,s)],1.0);
   model.add(con);

   // Outgoing t = 0
   con = IloRange(env, 0, 0);
   for(auto j : out_adj_list.at(t)) 
      con.setLinearCoef(sigma_vars[NODE_PAIR(t,j)],1.0);
   model.add(con);

   // MTZ
   // p[i] - p[j] + (n-1)*x[i,j] >= (n-2);
   for (NODE i : nodes)
   { 
      if(i!=s)
         for(auto j : out_adj_list.at(i))
            if(j!=t)
            {
               IloRange con;
               con = IloRange(env, -IloInfinity, (int)nodes.size()-2);

               con.setLinearCoef(p_var[i],1.0);
               con.setLinearCoef(p_var[j],-1.0);
               con.setLinearCoef(sigma_vars[NODE_PAIR(i,j)], ((int)nodes.size()) - 1.0);
               model.add(con);
            }
   }


} /* MTZ */

/** DL formulation */
void ElppSolver::build_problem_dl()
{
   /***************/ 
   /* Build model */
   /***************/ 
   NODE s = st.first;
   NODE t = st.second;

   IloEnv env = model.getEnv();
   /* Add variables */
   char var_name[255];

   // p[i] 
   // Label of node i 
   for (auto i : nodes)
   {   
      snprintf(var_name, 255, "p_%d", i);
      IloNumVar var;
      var = IloNumVar(env,var_name);
      var.setBounds(0, (int)nodes.size());
      p_var[i] = var;
   }
   /*******************/
   /* Add constraints */
   /*******************/
   // Flow conservation constraints (on sigma)
   for (NODE node : nodes)
   {
      IloRange con;

      IloNum coeff;
      if( node == s)
         coeff = 1.0;
      else if ( node == t)
         coeff = -1.0;
      else
         coeff = 0.0;

      con = IloRange(env, coeff, coeff);

      for(auto j : out_adj_list.at(node))
         con.setLinearCoef(sigma_vars[NODE_PAIR(node,j)],1.0);
      for(auto j : in_adj_list.at(node))
         con.setLinearCoef(sigma_vars[NODE_PAIR(j,node)],-1.0);
      model.add(con);
   }

   // Outgoing degree constraints (on sigma)
   for (NODE node : nodes)
   {
      IloRange con;
      con = IloRange(env, -IloInfinity, 1.0);

      for(auto j : out_adj_list.at(node)) 
         con.setLinearCoef(sigma_vars[NODE_PAIR(node,j)],1.0);
      model.add(con);
   }

   // Incoming s = 0 
   IloRange con = IloRange(env, 0, 0);
   for(auto j : in_adj_list.at(s)) 
      con.setLinearCoef(sigma_vars[NODE_PAIR(j,s)],1.0);
   model.add(con);

   // Outgoing t = 0
   con = IloRange(env, 0, 0);
   for(auto j : out_adj_list.at(t)) 
      con.setLinearCoef(sigma_vars[NODE_PAIR(t,j)],1.0);
   model.add(con);

   // DL
   // p[i] - p[j] + (n-1)*x[i,j] + (n-3)*x[j,i] >= (n-2) [important: s!=i, t!=j]
   for (NODE i : nodes)
   { 
      if(i!=s)
         for(auto j : out_adj_list.at(i))
            if(t!=j)
            {
               IloRange con;
               con = IloRange(env, -IloInfinity, (int)nodes.size()-2);

               con.setLinearCoef(p_var[i],1.0);
               con.setLinearCoef(p_var[j],-1.0);
               con.setLinearCoef(sigma_vars[NODE_PAIR(i,j)], ((int)nodes.size()) - 1.0);
               // If (j,i) in the graph
               if(!(find(in_adj_list.at(i).begin(),in_adj_list.at(i).end(),j) == in_adj_list.at(i).end()))
                  con.setLinearCoef(sigma_vars[NODE_PAIR(j,i)], ((int)nodes.size()) - 3.0 );
               model.add(con);
            }
   }

} /* DL */


/** Update the objective fuction */
void ElppSolver::update_problem(
      const unordered_map<NODE_PAIR, IloNum>& obj_coeff 
      )
{
   /***************/ 
   /* Build model */
   /***************/ 
   NODE s = st.first;
   NODE t = st.second;

   /* Build objective function */
   IloEnv env = model.getEnv();
   IloExpr totalCost(env);
   IloNum objcoeff = 0.0;
   IloNumVar var;
   for (auto arc : arcs)
   {   
      var = sigma_vars[arc];
      objcoeff = obj_coeff.at(arc);
      //objective.setLinearCoef(var, objcoeff);
      totalCost += var*objcoeff;
   }

   objective = IloObjective(env, totalCost, IloObjective::Maximize);
   model.add(objective);

}

/** Update the problem with bounds (e.g., branching info) */
void ElppSolver::update_problem(
      const unordered_map<NODE_PAIR, IloNum>& obj_coeff, 
      const map<NODE_PAIR, IloNum>& lbs,
      const map<NODE_PAIR, IloNum>& ubs
      )
{
   /***************/ 
   /* Build model */
   /***************/ 
   NODE s = st.first;
   NODE t = st.second;

   IloEnv env = model.getEnv();

   if(objective.getImpl() == NULL)
   {
      objective = IloObjective(env, 0.0, IloObjective::Maximize);
      model.add(objective);
   }

   /* Update objective function */
   for (auto arc : arcs)
   {   
      IloNumVar var;
      var = sigma_vars[arc];

      IloNum objcoeff = 0.0;
      objcoeff = obj_coeff.at(arc);
      objective.setLinearCoef(var, objcoeff);

      // Update bounds
      IloNum upper = ubs.at(arc);
      IloNum lower = lbs.at(arc);
      var.setBounds(lower, upper);
   }
}


/** Update the problem */
void ElppSolver::update_problem(
      const unordered_map<NODE_PAIR, IloNum>& obj_coeff, 
      const map<NODE_PAIR, IloNum>& lbs,
      const map<NODE_PAIR, IloNum>& ubs,
      const unordered_map<NODE_PAIR, IloNum>& lhs, 
      IloNum rhs,
      bool use_extra_con
      )
{
   /***************/ 
   /* Build model */
   /***************/ 
   NODE s = st.first;
   NODE t = st.second;

   IloEnv env = model.getEnv();

   if(objective.getImpl() == NULL)
   {
      objective = IloObjective(env, 0.0, IloObjective::Maximize);
      model.add(objective);
   }

   /* Update objective function */
   for (auto arc : arcs)
   {   
      IloNumVar var;
      var = sigma_vars[arc];

      IloNum objcoeff = 0.0;
      objcoeff = obj_coeff.at(arc);
      objective.setLinearCoef(var, objcoeff);

      // Update branching information
      IloNum upper = ubs.at(arc);
      IloNum lower = lbs.at(arc);
      var.setBounds(lower, upper);
   }

   /* Lower bound on *real* objective function value */
   if ( use_extra_con )
   {
      if(extra_con.getImpl() == NULL)
         model.remove(extra_con);
      /** Reduced cost (lhs) must be greater than w (rhs) to be an improving column */
      extra_con.end();
      extra_con = IloRange(env, rhs + 1e-05, IloInfinity);

      for (auto arc : arcs)
      {   
         extra_con.setLinearCoef(sigma_vars[arc], lhs.at(arc));
      }
      model.add(extra_con);
   }

}

/** Solve problem */
void ElppSolver::solve()
{
   LOG << "------------------ SOLVING ELPP FOR " << st.first << " - " << st.second << "-------------" << endl;

   double start_time = cplex.getCplexTime();
   double start_ticks = cplex.getDetTime();

   cplex.solve();

   //   for(NODE_PAIR arc : arcs)
   //      if(cplex.getValue(sigma_vars[arc])>0)
   //      {
   //         cout << arc << " " << cplex.getValue(sigma_vars[arc]);
   //         if(u_var.count(arc) == 1)
   //            cout << " (" << cplex.getValue(u_var[arc]) << ")";
   //         cout << endl;
   //      }

   elapsed_time = cplex.getCplexTime() - start_time;
   elapsed_ticks = cplex.getDetTime() - start_ticks;

}

/** Solve root node (LP relaxation) */
void ElppSolver::solveRoot()
{
   LOG << "------------------ Solving ELPP root node for " << st.first << " - " << st.second << "-------------" << endl;

   // Disable all Cplex cuts
   cplex.setParam(IloCplex::ImplBd, -1);
   cplex.setParam(IloCplex::Cliques, -1);
   cplex.setParam(IloCplex::Covers, -1);
   cplex.setParam(IloCplex::DisjCuts, -1);
   cplex.setParam(IloCplex::FlowPaths, -1);
   cplex.setParam(IloCplex::FlowCovers, -1);
   cplex.setParam(IloCplex::FracCuts, -1);
   cplex.setParam(IloCplex::GUBCovers, -1);
   //cplex.setParam(IloCplex::LiftProjCuts, -1); //enum not found? //SplitCuts?
   cplex.setParam(IloCplex::MCFCuts, -1);
   cplex.setParam(IloCplex::MIRCuts, -1);
   cplex.setParam(IloCplex::ZeroHalfCuts, -1);

   // Disable presolve
   cplex.setParam(IloCplex::PreInd, false);

   // Limit to 0-th node
   cplex.setParam(IloCplex::NodeLim, 0);

   double start_time = cplex.getCplexTime();
   double start_ticks = cplex.getDetTime();

   if((formulation!=1) || nodes.size() < 100)
   {
      cplex.solve();
      elapsed_time = cplex.getCplexTime() - start_time;
      elapsed_ticks = cplex.getDetTime() - start_ticks;
   }

}

/** Solve LP relaxation manually */
void ElppSolver::solveLP()
{
   LOG << "------------------ Solving LP for " << st.first << " - " << st.second << "-------------" << endl;

   double start_time = cplex.getCplexTime();
   double start_ticks = cplex.getDetTime();

   bool separated = true;
   ncuts = 0;
   while(separated)
   {
      separated = false;
      LOG << "Time elapsed: " << cplex.getCplexTime() - start_time << endl;
      LOG << "Set new time limit of: " << timelimit - (cplex.getCplexTime() - start_time) << endl;
      cplex.setParam(IloCplex::TiLim, max(0.0,timelimit - (cplex.getCplexTime() - start_time))); //TILIM - elapsed time
      LOG << "Solve updated LP" << endl;
      cplex.solve();
      if(cplex.getStatus() != IloAlgorithm::Optimal)
         break;
      LOG << "Current obj value: " << cplex.getObjValue()<< endl;

      IloEnv env = cplex.getEnv();
      IloNumArray val = IloNumArray(env, sigma_vars.size());
      cplex.getValues(val, x_vararray);

      unordered_map<NODE_PAIR, IloNum> xSol;
      for(NODE_PAIR arc : arcs)
      {
         xSol[arc] = val[index[arc]];
         //LOG << arc.first << " " << arc.second << ": " <<  xSol[arc] << endl;
      }

      vector<IloExpr> cutLhs, cutRhs;
      vector<IloRange> cons;
      vector<double> violation;

      switch(formulation)
      {
         case SC: /* only SC */
            if(!separate_weak(env, xSol, nodes, arcs, st, sigma_vars, out_adj_list, in_adj_list, cutLhs, cutRhs, violation))
               separate_sc(env, xSol, nodes, arcs, st, sigma_vars, out_adj_list, in_adj_list, cutLhs, cutRhs, violation);
            break;
         case MinCut:
            if(!separate_sc(env, xSol, nodes, arcs, st, sigma_vars, out_adj_list, in_adj_list, cutLhs, cutRhs, violation))
               separate_min_cut(env, xSol, nodes, arcs, st, sigma_vars, out_adj_list, in_adj_list, cutLhs, cutRhs, violation);
            break;
         case DFJ:
            if(!separate_sc_dfj(env, xSol, nodes, arcs, st, sigma_vars, out_adj_list, in_adj_list, cutLhs, cutRhs, violation))
               separate_min_cut_dfj(env, xSol, nodes, arcs, st, sigma_vars, out_adj_list, in_adj_list, cutLhs, cutRhs, violation);
            break;
         case MCFsep:
            if(!separate_sc_mf(env, xSol, nodes, arcs, st, sigma_vars, qq_var, zz_var, out_adj_list, in_adj_list, cons))
               separate_min_cut_mf(env, xSol, nodes, arcs, st, sigma_vars, qq_var, zz_var, out_adj_list, in_adj_list, cons);
            break;
      }

      switch(formulation)
      {
         case SC:
         case MinCut: 
         case DFJ:
            {  // Only need to get the max_cuts maximally-violated inequalities
               vector<int> p(violation.size()); /* vector with indices */
               iota(p.begin(), p.end(), 0);  /* increasing */
               bool sorted = false;

               int attempts;
               if(max_cuts < 0)
                  attempts = violation.size();
               else
               {
                  attempts = min(max_cuts, (int) violation.size());
                  partial_sort(p.begin(), p.begin() + attempts, p.end(), [&](int i, int j){ return violation[i] > violation[j]; }); /* sort indices according to violation */
                  sorted = true;
               }

               for(unsigned int i=0; i<attempts; ++i)
               {
                  if(violation[p[i]] >= TOL)
                  { 
                     LOG << "Adding user cut for the " << i+1 << "-th maximally violated constraint. Violation: " << violation[p[i]] << endl;
                     try{
                        LOG << (cutLhs[p[i]] >= cutRhs[p[i]]) << endl;
                        model.add(cutLhs[p[i]] >= cutRhs[p[i]]);
                        separated = true;
                        ++ncuts;
                     }
                     catch(IloException e)
                     {
                        cerr << "Cannot add cut" << endl;
                     }
                  }
                  else /* sorted, so no further violated ineq exist */
                     if(sorted)
                        break;
               }
               for(unsigned int i=0; i<cutLhs.size(); ++i)
               {
                  cutLhs[i].end();
                  cutRhs[i].end();
               }
            }
            break;

         case MCFsep: /*MCFsep -- only add 1 set of MCF constraints each time..*/
            for(unsigned int i=0; i<cons.size();++i)
            {
               LOG << cons[i] << endl;
               model.add(cons[i]);
               separated = true;
               ++ncuts;
            }
            break;

         default:
            break;
      }

      val.end();
   } /* while(separated) */
   
   //   cplex.solve();
   //   cerr << "Final obj value: " << cplex.getObjValue()<< endl;
   //   IloEnv env = cplex.getEnv();
   //   IloNumArray val = IloNumArray(env, sigma_vars.size());
   //   cplex.getValues(val, x_vararray);
   //
   //   unordered_map<NODE_PAIR, IloNum> xSol;
   //   for(NODE_PAIR arc : arcs)
   //   {
   //      xSol[arc] = val[index[arc]];
   //      LOG << arc.first << " " << arc.second << ": " <<  xSol[arc] << endl;
   //   }

   elapsed_time = cplex.getCplexTime() - start_time;
   elapsed_ticks = cplex.getDetTime() - start_ticks;

}


void ElppSolver::clear()
{
   LOG << "Clearing cuts and lazy constraints" << endl;
   cplex.clearUserCuts();
   cplex.clearLazyConstraints();
}

IloAlgorithm::Status ElppSolver::getStatus()
{
   return cplex.getStatus();
}

IloNum ElppSolver::getBestObjValue()
{
   return cplex.getBestObjValue();
}

IloNum ElppSolver::getObjValue()
{
   return cplex.getObjValue();
}

IloNum ElppSolver::getValue(NODE_PAIR arc)
{
   return cplex.getValue(sigma_vars[arc]);
}

IloNum ElppSolver::getValue(NODE node)
{
   return cplex.getValue(p_var[node]);
}

/* check if the current solution is feasible for the original problem */
bool ElppSolver::isInteger()
{
   if(cplex.getStatus() != IloAlgorithm::Optimal)
      return false;

   IloEnv env = cplex.getEnv();
   IloNumArray val = IloNumArray(env, sigma_vars.size());
   cplex.getValues(val, x_vararray);
   for (auto arc : arcs)
      if(val[index[arc]] <= 1.0 - TOL && val[index[arc]] >= 0.0 + TOL)
         return false;
   return true;
}

int ElppSolver::pathLength()
{
   if(cplex.getStatus() != IloAlgorithm::Optimal)
      return false;

   IloEnv env = cplex.getEnv();
   IloNumArray val = IloNumArray(env, sigma_vars.size());
   cplex.getValues(val, x_vararray);
   int counter = 0;
   for (auto arc : arcs)
      if(val[index[arc]] >= 1.0 - TOL)
         ++counter;
   return counter;
}
void ElppSolver::printInstance(
      string filename,
      const unordered_map<NODE_PAIR, IloNum>& obj_coeff, 
      const map<NODE_PAIR, IloNum>& lbs,
      const map<NODE_PAIR, IloNum>& ubs)
{
   ofstream outfile(filename);
   outfile << nodes.size() << " " << arcs.size() << endl;
   outfile << st.first << " " << st.second << endl;
   for(NODE node : nodes)
      outfile << node << endl;
   for(NODE_PAIR arc : arcs)
      outfile << arc.first << " " << arc.second << " " << obj_coeff.at(arc) << endl;
   for(NODE_PAIR arc : arcs)
      outfile << arc.first << " " << arc.second << " " << lbs.at(arc) << " " << ubs.at(arc) << endl;
}

void ElppSolver::writeLP(string filename)
{
   string ciccio = filename+".lp";
   cplex.exportModel(ciccio.c_str());
}

/* Print info file */
void ElppSolver::append_info(string info_filename, string instance_name)
{
   ofstream infofile(info_filename, std::ios_base::app);
   infofile << instance_name << "\t";
   infofile << st.first << "\t";
   infofile << st.second << "\t";
  
   infofile << which(formulation)  << "\t";

   if(cplex.getStatus() == IloAlgorithm::Optimal)
      infofile << cplex.getObjValue() << "\t";
   else 
      infofile << "-" << "\t";
   infofile << elapsed_time << "\t";
   infofile << elapsed_ticks << "\t";
   infofile << cplex.getNnodes() << "\t";
   if(relax)
      infofile << ncuts << "\t" << isInteger() << endl;
   else
   {
      infofile << cplex.getNcuts(IloCplex::CutUser) << "\t";
      infofile << tol << "\t";
      infofile << max_cuts << endl;
      //infofile << pathLength() << endl;
   }

}

