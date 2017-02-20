/**@file elpp.cpp
 * @brief Elementary Longest Path problem solver
 * @author Leonardo Taccari 
 */

#include "elpp.h"
#include "separation.h"
#include <iostream>
#include <algorithm> 
#define LOG if(false) cerr

using namespace std;
using namespace lemon;

// Lazy constraint: called only when integer feasible incumbent is found
class StrongComponentLazyCallbackI : public IloCplex::LazyConstraintCallbackI {
   vector<NODE> nodes;
   vector<NODE_PAIR> arcs;
   NODE_PAIR st;
   ARC_VARS sigma_vars;
   unordered_map<TRIPLET,IloNumVar> qq_var;
   unordered_map<NODE,IloNumVar> zz_var;
   ADJ_LIST out_adj_list;
   ADJ_LIST in_adj_list;
   IloNumVarArray x_vararray;
   INDEX_MAP index;
   const double tol;
   const int max_cuts;
   const ElppForm form;

   public:
      ILOCOMMONCALLBACKSTUFF(StrongComponentLazyCallback)
      StrongComponentLazyCallbackI(IloEnv env, vector<NODE> nodes_, vector<NODE_PAIR> arcs_, NODE_PAIR st_, ARC_VARS sigma_vars_, unordered_map<TRIPLET,IloNumVar> qq_var_, unordered_map<NODE,IloNumVar> zz_var_, ADJ_LIST out_adj_list_, ADJ_LIST in_adj_list_, IloNumVarArray x_vararray_, INDEX_MAP index_, double tol_, int max_cuts_, ElppForm form_)
      :  IloCplex::LazyConstraintCallbackI(env), nodes(nodes_), arcs(arcs_), st(st_), sigma_vars(sigma_vars_), qq_var(qq_var_), zz_var(zz_var_), out_adj_list(out_adj_list_), 
      in_adj_list(in_adj_list_), x_vararray(x_vararray_), index(index_), tol(tol_), max_cuts(max_cuts_), form(form_) {}

   void main();
};

IloCplex::Callback StrongComponentLazyCallback(IloEnv env, vector<NODE> nodes, vector<NODE_PAIR> arcs, NODE_PAIR st, ARC_VARS sigma_vars, unordered_map<TRIPLET,IloNumVar> qq_var, unordered_map<NODE,IloNumVar> zz_var, ADJ_LIST out_adj_list, ADJ_LIST in_adj_list, IloNumVarArray x_vararray, INDEX_MAP index, double tol, int max_cuts, ElppForm form) {
   return (IloCplex::Callback(new (env) StrongComponentLazyCallbackI(env, nodes, arcs, st, sigma_vars, qq_var, zz_var, out_adj_list, in_adj_list, x_vararray, index, tol, max_cuts, form)));}

void StrongComponentLazyCallbackI::main()
{
   global_timer.start();
   IloEnv masterEnv = getEnv();

   LOG << "--STRONGCOMPONENT-LAZYCONSTR--" << endl;

   IloNumArray val = IloNumArray(masterEnv, sigma_vars.size());
   getValues(val, x_vararray);

   unordered_map<NODE_PAIR, IloNum> xSol;
   for(NODE_PAIR arc : arcs)
   {
      //LOG << arc.first << " " << arc.second << endl;
      xSol[arc] = val[index[arc]];
      LOG << arc.first << " " << arc.second << ": " <<  xSol[arc] << endl;
   }

   vector<IloExpr> cutLhs, cutRhs;
   vector<double> violation;
   vector<IloRange> cons;
   switch(form)
   {
      case SC:
      case MinCut: /* if MinCut, try SC separations, then MinCut */
         {
            separate_sc(masterEnv, xSol, nodes, arcs, st, sigma_vars, out_adj_list, in_adj_list, cutLhs, cutRhs, violation);

            // Only need to get the max_cuts maximally-violated inequalities
            vector<int> p(violation.size()); /* vector with indices */
            iota(p.begin(), p.end(), 0);  /* increasing */
            bool sorted = false;

            int attempts = 0;
            if(max_cuts < 0)
               attempts = (int) violation.size();
            else
            {
               attempts = min(max_cuts, (int) violation.size());
               partial_sort(p.begin(), p.begin() + attempts, p.end(), [&](int i, int j){ return violation[i] > violation[j]; }); /* sort indices according to violation */
               sorted = true;
            }

            for(unsigned int i=0; i<attempts; ++i)
            {
               LOG << violation[p[i]] << endl;
               if(violation[p[i]] >= tol)
               { 
                  LOG << "Adding lazy constraint for the " << i+1 << "-th maximally violated constraint." << endl;
                  try{
                     LOG << (cutLhs[p[i]] >= cutRhs[p[i]]) << endl;
                     add(cutLhs[p[i]] >= cutRhs[p[i]]).end();
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
            for(unsigned int i=0;i<cutLhs.size();++i)
            {
               cutLhs[i].end();
               cutRhs[i].end();
            }
         }
         break;
      
      case MCFsep:
         separate_sc_mf(masterEnv, xSol, nodes, arcs, st, sigma_vars, qq_var, zz_var, out_adj_list, in_adj_list, cons);
         for(unsigned int i=0; i<cons.size();++i)
         {
            LOG << i << ": " << cons[i] << endl;
            add(cons[i]).end();
         }
         break;
   
      default: 
         break;
   }
   
   LOG << "--END ILOLAZYCONSTRAINTCALLBACK--" << global_timer.userTime() << endl;
   global_timer.stop();
   return;
}

// User cut: called also on fractional solutions 
class ElppCutCallbackI : public IloCplex::UserCutCallbackI {
   vector<NODE> nodes;
   vector<NODE_PAIR> arcs;
   NODE_PAIR st;
   ARC_VARS sigma_vars;
   unordered_map<TRIPLET,IloNumVar> qq_var;
   unordered_map<NODE,IloNumVar> zz_var;
   ADJ_LIST out_adj_list;
   ADJ_LIST in_adj_list;
   IloNumVarArray x_vararray;
   INDEX_MAP index;
   const double tol;
   const int max_cuts;
   const ElppForm form;

   public:
      ILOCOMMONCALLBACKSTUFF(ElppCutCallback)
      ElppCutCallbackI(IloEnv env, vector<NODE> nodes_, vector<NODE_PAIR> arcs_, NODE_PAIR st_, ARC_VARS sigma_vars_, unordered_map<TRIPLET,IloNumVar> qq_var_, unordered_map<NODE,IloNumVar> zz_var_, ADJ_LIST out_adj_list_, ADJ_LIST in_adj_list_, IloNumVarArray x_vararray_, INDEX_MAP index_, double tol_, int max_cuts_, ElppForm form_)
      :  IloCplex::UserCutCallbackI(env), nodes(nodes_), arcs(arcs_), st(st_), sigma_vars(sigma_vars_), qq_var(qq_var_), zz_var(zz_var_), out_adj_list(out_adj_list_), 
      in_adj_list(in_adj_list_), x_vararray(x_vararray_), index(index_), tol(tol_), max_cuts(max_cuts_), form(form_) {}

   void main();
};

IloCplex::Callback ElppCutCallback(IloEnv env, vector<NODE> nodes, vector<NODE_PAIR> arcs, NODE_PAIR st, ARC_VARS sigma_vars, unordered_map<TRIPLET,IloNumVar> qq_var, unordered_map<NODE,IloNumVar> zz_var, ADJ_LIST out_adj_list, ADJ_LIST in_adj_list, IloNumVarArray x_vararray, INDEX_MAP index, double tol, int max_cuts, ElppForm form) {
   return (IloCplex::Callback(new (env) ElppCutCallbackI(env, nodes, arcs, st, sigma_vars, qq_var, zz_var, out_adj_list, in_adj_list, x_vararray, index, tol, max_cuts, form)));}

void ElppCutCallbackI::main()
{
   // Skip the separation if not at the end of the cut loop
   if ( !isAfterCutLoop() )
      return;

   global_timer.start();
   LOG << "--ELPP USERCUT--" << endl;

   IloEnv masterEnv = getEnv();
   IloNumArray val = IloNumArray(masterEnv, sigma_vars.size());
   getValues(val, x_vararray);

   unordered_map<NODE_PAIR, IloNum> xSol;
   for(NODE_PAIR arc : arcs)
   {
      //LOG << arc.first << " " << arc.second << endl;
      xSol[arc] = val[index[arc]];
      //LOG << arc.first << " " << arc.second << ": " <<  xSol[arc] << endl;
   }

   vector<IloExpr> cutLhs, cutRhs;
   vector<double> violation;
   vector<IloRange> cons;

   switch(form)
   {
      case SC:
      case MinCut: /* if MinCut, try SC separations, then MinCut */
         {
            if(!separate_sc(masterEnv, xSol, nodes, arcs, st, sigma_vars, out_adj_list, in_adj_list, cutLhs, cutRhs, violation) && form == MinCut) 
               separate_min_cut(masterEnv, xSol, nodes, arcs, st, sigma_vars, out_adj_list, in_adj_list, cutLhs, cutRhs, violation);

            // Only need to get the max_cuts maximally-violated inequalities
            vector<int> p(violation.size()); /* vector with indices */
            iota(p.begin(), p.end(), 0);  /* increasing */
            bool sorted = false;

            int attempts = 0;
            if(max_cuts < 0)
               attempts = (int) violation.size();
            else
            {
               attempts = min(max_cuts, (int) violation.size());
               partial_sort(p.begin(), p.begin() + attempts, p.end(), [&](int i, int j){ return violation[i] > violation[j]; }); /* sort indices according to violation */
               sorted = true;
            }

            for(unsigned int i=0; i<attempts; ++i)
            {
               if(violation[p[i]] >= tol)
               {
                  LOG << "Adding user cut for the " << i+1 << "-th maximally violated constraint. Violation (" << p[i] << "): " << violation[p[i]] << endl;
                  try{
                     LOG << (cutLhs[p[i]] >= cutRhs[p[i]]) << endl;
                     add(cutLhs[p[i]] >= cutRhs[p[i]]).end();
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
         if(!separate_sc_mf(masterEnv, xSol, nodes, arcs, st, sigma_vars, qq_var, zz_var, out_adj_list, in_adj_list, cons))
            separate_min_cut_mf(masterEnv, xSol, nodes, arcs, st, sigma_vars, qq_var, zz_var, out_adj_list, in_adj_list, cons);

         for(unsigned int i=0; i<cons.size();++i)
         {
            LOG << cons[i] << endl;
            add(cons[i]);
         }
         break;

      default:
         break;
   }

   LOG << "--END ILOUSERCUTCALLBACK--" << global_timer.userTime() << endl;
   global_timer.stop();
   return;
}

