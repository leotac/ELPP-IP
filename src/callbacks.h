/**@file callbacks.h
 * @brief Elementary Longest Path problem solver
 * @author Leonardo Taccari 
 */

#include "elpp.h"
#include "graph.h"
#include "separation.h"
#include <iostream>
#include <algorithm> 
#define LOG if(false) cerr

using namespace std;
using namespace lemon;

// Lazy constraint: called only when integer feasible incumbent is found
class StrongComponentLazyCallbackI : public IloCplex::LazyConstraintCallbackI {
   NODE_PAIR st;
   std::shared_ptr<Graph> G;
   ARC_VARS sigma_vars;
   unordered_map<TRIPLET,IloNumVar> qq_var;
   unordered_map<NODE,IloNumVar> zz_var;
   IloNumVarArray x_vararray;
   INDEX_MAP index;
   const double tol;
   const int max_cuts;
   const ElppForm form;

   public:
      ILOCOMMONCALLBACKSTUFF(StrongComponentLazyCallback)
      StrongComponentLazyCallbackI(IloEnv env, std::shared_ptr<Graph> graph, NODE_PAIR st_, ARC_VARS sigma_vars_, unordered_map<TRIPLET,IloNumVar> qq_var_, unordered_map<NODE,IloNumVar> zz_var_, IloNumVarArray x_vararray_, INDEX_MAP index_, double tol_, int max_cuts_, ElppForm form_)
      :  IloCplex::LazyConstraintCallbackI(env), G(graph), st(st_), sigma_vars(sigma_vars_), qq_var(qq_var_), zz_var(zz_var_), 
      x_vararray(x_vararray_), index(index_), tol(tol_), max_cuts(max_cuts_), form(form_) {}

   void main();
};

IloCplex::Callback StrongComponentLazyCallback(IloEnv env, std::shared_ptr<Graph> graph, NODE_PAIR st, ARC_VARS sigma_vars, unordered_map<TRIPLET,IloNumVar> qq_var, unordered_map<NODE,IloNumVar> zz_var, IloNumVarArray x_vararray, INDEX_MAP index, double tol, int max_cuts, ElppForm form) {
   return (IloCplex::Callback(new (env) StrongComponentLazyCallbackI(env, graph, st, sigma_vars, qq_var, zz_var, x_vararray, index, tol, max_cuts, form)));}

void StrongComponentLazyCallbackI::main()
{
   IloEnv masterEnv = getEnv();

   LOG << "--STRONGCOMPONENT-LAZYCONSTR--" << endl;

   IloNumArray val = IloNumArray(masterEnv, sigma_vars.size());
   getValues(val, x_vararray);

   unordered_map<NODE_PAIR, double> xSol;
   for(NODE_PAIR arc : G->arcs())
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
            separate_sc(masterEnv, xSol, G, st, sigma_vars, cutLhs, cutRhs, violation);

            // Only need to get the max_cuts maximally-violated inequalities
            vector<int> p(violation.size()); /* vector with indices */
            iota(p.begin(), p.end(), 0);  /* increasing */
            bool sorted = false;

            int attempts = 0;
            if(max_cuts < 0)
               attempts = int(violation.size());
            else
            {
               attempts = min(max_cuts, int(violation.size()));
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
                     LOG << "Cannot add cut" << endl;
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
         separate_sc_mf(masterEnv, xSol, G, st, sigma_vars, qq_var, zz_var, cons);
         for(unsigned int i=0; i<cons.size();++i)
         {
            LOG << i << ": " << cons[i] << endl;
            add(cons[i]).end();
         }
         break;
   
      default: 
         break;
   }
   
   LOG << "--END ILOLAZYCONSTRAINTCALLBACK--" << endl;
   return;
}

// User cut: called also on fractional solutions 
class ElppCutCallbackI : public IloCplex::UserCutCallbackI {
   std::shared_ptr<Graph> G;
   NODE_PAIR st;
   ARC_VARS sigma_vars;
   unordered_map<TRIPLET,IloNumVar> qq_var;
   unordered_map<NODE,IloNumVar> zz_var;
   IloNumVarArray x_vararray;
   INDEX_MAP index;
   const double tol;
   const int max_cuts;
   const ElppForm form;

   public:
      ILOCOMMONCALLBACKSTUFF(ElppCutCallback)
      ElppCutCallbackI(IloEnv env, std::shared_ptr<Graph> graph, NODE_PAIR st_, ARC_VARS sigma_vars_, unordered_map<TRIPLET,IloNumVar> qq_var_, unordered_map<NODE,IloNumVar> zz_var_, IloNumVarArray x_vararray_, INDEX_MAP index_, double tol_, int max_cuts_, ElppForm form_)
      :  IloCplex::UserCutCallbackI(env), G(graph), st(st_), sigma_vars(sigma_vars_), qq_var(qq_var_), zz_var(zz_var_), 
      x_vararray(x_vararray_), index(index_), tol(tol_), max_cuts(max_cuts_), form(form_) {}

   void main();
};

IloCplex::Callback ElppCutCallback(IloEnv env,std::shared_ptr<Graph> graph, NODE_PAIR st, ARC_VARS sigma_vars, unordered_map<TRIPLET,IloNumVar> qq_var, unordered_map<NODE,IloNumVar> zz_var, IloNumVarArray x_vararray, INDEX_MAP index, double tol, int max_cuts, ElppForm form) {
   return (IloCplex::Callback(new (env) ElppCutCallbackI(env, graph, st, sigma_vars, qq_var, zz_var, x_vararray, index, tol, max_cuts, form)));}

void ElppCutCallbackI::main()
{
   // Skip the separation if not at the end of the cut loop
   if ( !isAfterCutLoop() )
      return;

   LOG << "--ELPP USERCUT--" << endl;

   IloEnv masterEnv = getEnv();
   IloNumArray val = IloNumArray(masterEnv, sigma_vars.size());
   getValues(val, x_vararray);

   unordered_map<NODE_PAIR, double> xSol;
   for(NODE_PAIR arc : G->arcs())
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
            if(!separate_sc(masterEnv, xSol, G, st, sigma_vars, cutLhs, cutRhs, violation) && form == MinCut) 
               separate_min_cut(masterEnv, xSol, G, st, sigma_vars, cutLhs, cutRhs, violation);

            // Only need to get the max_cuts maximally-violated inequalities
            vector<int> p(violation.size()); /* vector with indices */
            iota(p.begin(), p.end(), 0);  /* increasing */
            bool sorted = false;

            int attempts = 0;
            if(max_cuts < 0)
               attempts = int(violation.size());
            else
            {
               attempts = min(max_cuts, int(violation.size()));
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
                     LOG << "Cannot add cut" << endl;
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
         if(!separate_sc_mf(masterEnv, xSol, G, st, sigma_vars, qq_var, zz_var, cons))
            separate_min_cut_mf(masterEnv, xSol, G, st, sigma_vars, qq_var, zz_var, cons);

         for(unsigned int i=0; i<cons.size();++i)
         {
            LOG << cons[i] << endl;
            add(cons[i]);
         }
         break;

      default:
         break;
   }

   LOG << "--END ILOUSERCUTCALLBACK--" << endl;
   return;
}

