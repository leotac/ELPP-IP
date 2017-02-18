#ifndef __separation_h__
#define __separation_h__

#include "elpp.h"
#include "unionfind.h"
#include <iostream>
#include <fstream>
#include <algorithm> 
#include <lemon/smart_graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/connectivity.h>
#include <lemon/preflow.h>
#include <lemon/time_measure.h>

using namespace std;
using namespace lemon;

typedef SmartDigraph::Node LemonNode;
typedef SmartDigraph::Arc LemonArc;

void build_support_graph(SmartDigraph& support_graph,  unordered_map<NODE,LemonNode>& v_nodes,  map<LemonNode,NODE>& rev_nodes, 
      const unordered_map<NODE_PAIR, IloNum>& xSol, const vector<NODE>& nodes, const ADJ_LIST& out_adj_list, NODE s, NODE t);

void build_cap_graph(SmartDigraph& cap_graph, SmartDigraph::ArcMap<double>& x_capacities, unordered_map<NODE,LemonNode>& v_nodes,  map<LemonNode,NODE>& rev_nodes, 
      const unordered_map<NODE_PAIR, IloNum>& xSol, const vector<NODE>& nodes, const ADJ_LIST& out_adj_list, NODE s, NODE t);

bool separate_sc(IloEnv masterEnv, const unordered_map<NODE_PAIR, IloNum>& xSol, const vector<NODE>& nodes, const vector<NODE_PAIR>& arcs, 
      NODE_PAIR st, const ARC_VARS& sigma_vars, const ADJ_LIST& out_adj_list, const ADJ_LIST& in_adj_list,
      vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation);

bool separate_sc_dfj(IloEnv masterEnv, const unordered_map<NODE_PAIR, IloNum>& xSol, const vector<NODE>& nodes, const vector<NODE_PAIR>& arcs, 
      NODE_PAIR st, const ARC_VARS& sigma_vars, const ADJ_LIST& out_adj_list, const ADJ_LIST& in_adj_list,
      vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation);

bool separate_sc_mf(IloEnv masterEnv, const unordered_map<NODE_PAIR, IloNum>& xSol, const vector<NODE>& nodes, const vector<NODE_PAIR>& arcs, 
      NODE_PAIR st, ARC_VARS& sigma_var, unordered_map<TRIPLET,IloNumVar>& qq_var, unordered_map<NODE,IloNumVar>& zz_var, const ADJ_LIST& out_adj_list, 
      const ADJ_LIST& in_adj_list, vector<IloRange>& cons);

bool separate_weak(IloEnv masterEnv, const unordered_map<NODE_PAIR, IloNum>& xSol, const vector<NODE>& nodes, const vector<NODE_PAIR>& arcs, 
      NODE_PAIR st, const ARC_VARS& sigma_vars, const ADJ_LIST& out_adj_list, const ADJ_LIST& in_adj_list,
      vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation);

bool separate_min_cut(IloEnv masterEnv, const unordered_map<NODE_PAIR, IloNum>& xSol, const vector<NODE>& nodes, const vector<NODE_PAIR>& arcs, 
      NODE_PAIR st, const ARC_VARS& sigma_vars, const ADJ_LIST& out_adj_list, const ADJ_LIST& in_adj_list,
      vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation);

bool separate_min_cut_dfj(IloEnv masterEnv, const unordered_map<NODE_PAIR, IloNum>& xSol, const vector<NODE>& nodes, const vector<NODE_PAIR>& arcs, 
      NODE_PAIR st, const ARC_VARS& sigma_vars, const ADJ_LIST& out_adj_list, const ADJ_LIST& in_adj_list,
      vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation);

bool separate_min_cut_mf(IloEnv masterEnv, const unordered_map<NODE_PAIR, IloNum>& xSol, const vector<NODE>& nodes, const vector<NODE_PAIR>& arcs, 
      NODE_PAIR st, ARC_VARS& sigma_var, unordered_map<TRIPLET,IloNumVar>& qq_var, unordered_map<NODE,IloNumVar>& zz_var, const ADJ_LIST& out_adj_list, 
      const ADJ_LIST& in_adj_list, vector<IloRange>& cons);

#endif
