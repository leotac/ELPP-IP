/**@file elpp.cpp
 * @brief Elementary Longest Path problem solver
 * @author Leonardo Taccari 
 */

#include "elpp.h"
#include "graph.h"
#include "unionfind.h"
#include "separation.h"
#include <iostream>
#include <fstream>
#include <algorithm> 
#include <lemon/smart_graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/connectivity.h>
#include <lemon/preflow.h>
#include <lemon/time_measure.h>
#define LOG if(false) cerr
#define TOL 0.001

using namespace std;
using namespace lemon;


void build_support_graph(SmartDigraph& support_graph,  unordered_map<NODE,LemonNode>& v_nodes,  map<LemonNode,NODE>& rev_nodes, 
      const unordered_map<NODE_PAIR, double>& xSol, std::shared_ptr<Graph> G, NODE s, NODE t)
{
   LemonNode a, b;
   for(NODE i : G->nodes())
   {
      if(i!=s && i!=t)
         for(NODE j : G->outgoing_from(i))
         {
            if(j!=s && j!=t && xSol.at(NODE_PAIR(i,j)) > TOL)
            {
               if(v_nodes.count(i) == 0)
               {
                  a = support_graph.addNode();
                  v_nodes[i] = a;
                  rev_nodes[a] = i;
               }
               if(v_nodes.count(j) == 0)
               {
                  b = support_graph.addNode();
                  v_nodes[j] = b;
                  rev_nodes[b] = j;
               }
               support_graph.addArc(v_nodes[i],v_nodes[j]);
               LOG << "added arc: " << i << " " << j << endl;
            }
         }
   }
}

void build_cap_graph(SmartDigraph& cap_graph, SmartDigraph::ArcMap<double>& x_capacities, unordered_map<NODE,LemonNode>& v_nodes,  map<LemonNode,NODE>& rev_nodes, 
      const unordered_map<NODE_PAIR, double>& xSol, std::shared_ptr<Graph> G, NODE s, NODE t)
{
   LemonNode a, b;
   LemonArc arc;
   for(NODE i : G->nodes())
   {
      if(i!=s && i!=t) //no need for outgoing from t
         for(NODE j : G->outgoing_from(i))
         {
            if(j!=s /*&& j!=st.second*/) // && xSol.at(NODE_PAIR(i,j)) > 0.001)
            {
               if(v_nodes.count(i) == 0)
               {
                  a = cap_graph.addNode();
                  v_nodes[i] = a;
                  rev_nodes[a] = i;
               }
               if(v_nodes.count(j) == 0)
               {
                  b = cap_graph.addNode();
                  v_nodes[j] = b;
                  rev_nodes[b] = j;
               }
               arc = cap_graph.addArc(v_nodes[i],v_nodes[j]);
               x_capacities[arc] = xSol.at(NODE_PAIR(i,j));
               LOG << "added arc: " << i << " " << j;
               LOG << " with capacity: " << x_capacities[arc] << endl;
            }
         }
   }
}


/***************************************************
            SEPARATION PROCEDURES
 ****************************************************/

// Weak Component separation (only find disconnected components)
bool separate_weak(IloEnv masterEnv, const unordered_map<NODE_PAIR, double>& xSol, std::shared_ptr<Graph> G, 
      NODE_PAIR st, const ARC_VARS& sigma_vars,
      vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation)
{
   bool ret = false;

   /* Find (weakly) connected components - UnionFind*/
   cutLhs = vector<IloExpr>();
   cutRhs = vector<IloExpr>();
   violation = vector<double>();

   UnionFind<NODE> forest(G->nodes());
   NODE u,v;
   unordered_map<NODE,bool> reached;   //default: false
   for(auto& arc : G->arcs()){
      u = arc.first;
      v = arc.second;
      if(xSol.at(arc) > TOL)
      {
         reached[u] = true;
         reached[v] = true;
         if(forest.find_set(u) != forest.find_set(v))
            forest.join(u,v);
      }
   }
   unordered_map<NODE, IloExpr> cutLhsMap;
   unordered_map<NODE,IloExpr> cutRhsMap;
   unordered_map<NODE, int> cardinality;
   unordered_map<NODE, double> out_degree;
   unordered_map<NODE, double> max_node_degree;
   unordered_map<NODE, double> node_out_degree;
   for(NODE i : G->nodes())
   {
      // Not connected to s!
      if(reached[i] == true && forest.find_set(i) != forest.find_set(st.first))
      {

         LOG << "--------------- node " << i << endl;
         NODE comp = forest.find_set(i);
         if(0 == cardinality[comp]++)
         {  LOG << "Initialized lhs" << endl;
            cutLhsMap[comp] = IloExpr(masterEnv);
         }
         for (NODE j : G->outgoing_from(i))
         {
            node_out_degree[i] += xSol.at(NODE_PAIR(i,j));
            if( comp != forest.find_set(j))
            {
               out_degree[comp] += xSol.at(NODE_PAIR(i,j));
               cutLhsMap[comp] += (sigma_vars.at(NODE_PAIR(i,j)));
            }
         }

         LOG << "degree: " << node_out_degree[i] << endl;
         if(node_out_degree[i] > max_node_degree[comp])
         {
            LOG << "add rhs " << endl;
            max_node_degree[comp] = node_out_degree[i];
            cutRhsMap[comp] = IloExpr(masterEnv);
            for (NODE j : G->outgoing_from(i))
               cutRhsMap[comp] += (sigma_vars.at(NODE_PAIR(i,j)));
         }
      }
   }

   for(auto element : cutLhsMap)
   {
      NODE comp = element.first;
      IloExpr lhs = element.second;
      IloExpr rhs = cutRhsMap[comp];
      double v = max_node_degree[comp] - out_degree[comp];
      if(v >= TOL)
      {
         cutLhs.push_back(lhs);
         cutRhs.push_back(rhs);
         violation.push_back(v);
         ret = true;
      }
   }

   return ret;
}


// Strong Component separation
bool separate_sc(IloEnv masterEnv, const unordered_map<NODE_PAIR, double>& xSol, std::shared_ptr<Graph> G, 
      NODE_PAIR st, const ARC_VARS& sigma_vars,
      vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation)
{
   /** Strong component separation */
   bool ret = false;
   
   // Build support graph
   SmartDigraph support_graph;
   unordered_map<NODE,LemonNode> v_nodes;
   map<LemonNode,NODE> rev_nodes;
   NODE s = st.first;
   NODE t = st.second;
   build_support_graph(support_graph, v_nodes, rev_nodes, xSol, G, s, t);

   // Search for strong components
   SmartDigraph::NodeMap<int> nodemap(support_graph);
   int components = stronglyConnectedComponents(support_graph, nodemap);
   //LOG << components << " strong components" << endl;
   //NODE N;
   //LemonNode n;
   //for(auto x : v_nodes)
   //{
   //   N = x.first;
   //   n = x.second;
   //   LOG << N << ": " << nodemap[n] << endl;
   //}

   cutLhs = vector<IloExpr>(components);
   cutRhs = vector<IloExpr>(components);
   violation = vector<double>(components);
   vector<double> cardinality(components,0);
   vector<double> out_degree(components,0);
   vector<double> max_node_degree(components,0);
   SmartDigraph::NodeMap<double> node_out_degree(support_graph,0);
   for(SmartDigraph::NodeIt i(support_graph); i!=INVALID; ++i)
   {
      LOG << "--------------- node " << rev_nodes[i] << endl;
      int comp = nodemap[i];
      if(0 == cardinality[comp]++)
      {  LOG << "Initialized lhs" << endl;
         cutLhs[comp] = IloExpr(masterEnv);
      }
      for (NODE j : G->outgoing_from(rev_nodes[i]))
      {
         node_out_degree[i] += xSol.at(NODE_PAIR(rev_nodes[i],j));
         if( v_nodes.count(j) == 0 || comp != nodemap[v_nodes[j]])
         {
            out_degree[comp] += xSol.at(NODE_PAIR(rev_nodes[i],j));
            cutLhs[comp] += (sigma_vars.at(NODE_PAIR(rev_nodes[i],j)));
         }
      }

      LOG << "degree: " << node_out_degree[i] << endl;
      if(node_out_degree[i] >= max_node_degree[comp] + TOL)
      {
         LOG << "add rhs " << endl;
         max_node_degree[comp] = node_out_degree[i];
         cutRhs[comp] = IloExpr(masterEnv);
         for (NODE j : G->outgoing_from(rev_nodes[i]))
            //for (SmartDigraph::OutArcIt arc(support_graph, i); arc!=INVALID; ++arc)
            //   if(comp != nodemap[support_graph.target(arc)])
            cutRhs[comp] += (sigma_vars.at(NODE_PAIR(rev_nodes[i],j)));
      }
   }

   for(int i=0;i<components;++i)
   {
      LOG << "component " << i << ": " << cardinality[i] << endl;
      LOG << "delta(S) " << out_degree[i] << endl;
      LOG << "max(delta(i)) " << max_node_degree[i] << endl;
      LOG << "lhs: " << cutLhs[i] << endl;
      LOG << "rhs: " << cutRhs[i] << endl;
      violation[i] = max_node_degree[i] - out_degree[i];
      if(violation[i] >= TOL)
         ret = true;        
   }

   //   for(SmartDigraph::NodeIt i(support_graph); i!=INVALID; ++i)
   //      LOG << "node " << rev_nodes[i] << ": " << "delta=" << node_out_degree[i] << endl;

   return ret;
}


// Min-cut separation
bool separate_min_cut(IloEnv masterEnv, const unordered_map<NODE_PAIR, double>& xSol, std::shared_ptr<Graph> G,
      NODE_PAIR st, const ARC_VARS& sigma_var, 
      vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation)
{
   // Build graph with x values as capacities
   SmartDigraph cap_graph;
   SmartDigraph::ArcMap<double> x_capacities(cap_graph);
   unordered_map<NODE,LemonNode> v_nodes;
   map<LemonNode,NODE> rev_nodes;
   NODE s = st.first;
   NODE t = st.second;
   build_cap_graph(cap_graph, x_capacities, v_nodes, rev_nodes, xSol, G, s, t);

   LOG << "Built graph" << endl;
   cutLhs = vector<IloExpr>();
   cutRhs = vector<IloExpr>();
   violation = vector<double>();
   IloExpr newCutLhs;
   IloExpr newCutRhs;
   double newViolation;
   double min_cut_value;
   for(NODE v : G->nodes())
      if(v!=s && v!=t)
      {
         Preflow<SmartDigraph, SmartDigraph::ArcMap<double>> min_cut(cap_graph, x_capacities, v_nodes[v], v_nodes[t]);
         min_cut.runMinCut();
         min_cut_value = min_cut.flowValue();

         LOG << "Ran min-cut" << endl;
         // Compute out degree of v
         double node_out_degree = 0;
         for (NODE j : G->outgoing_from(v))
            node_out_degree += xSol.at(NODE_PAIR(v,j));

         LOG << v << endl;
         LOG << "Min-cut " << min_cut_value << endl;
         LOG << "Node degree " << node_out_degree << endl;

         //if mincut < degree(v), add delta+(S) >= delta+(v)
         if ( node_out_degree > min_cut_value )
         {
            newCutLhs = IloExpr(masterEnv);
            newCutRhs = IloExpr(masterEnv);
            newViolation = node_out_degree - min_cut_value;

            // for all nodes i in S
            for(SmartDigraph::NodeIt i(cap_graph); i!=INVALID; ++i)
            {
               //if in S
               if(min_cut.minCut(i))
                  // for all nodes j adjacent to i, not in S
                  for (NODE j : G->outgoing_from(rev_nodes[i]))
                  {
                     if( v_nodes.count(j) == 0 || !min_cut.minCut(v_nodes[j]))
                     {
                        newCutLhs += (sigma_var.at(NODE_PAIR(rev_nodes[i],j)));
                     }
                  }
            }

            for (NODE j : G->outgoing_from(v))
               newCutRhs += (sigma_var.at(NODE_PAIR(v,j)));

            cutLhs.push_back(newCutLhs);
            cutRhs.push_back(newCutRhs);
            violation.push_back(newViolation);

            LOG << "node " << v << endl;
            LOG << "cut " << cutLhs.size() << endl;
            LOG << "delta(S) " << min_cut_value << endl;
            LOG << "delta(v) " << node_out_degree << endl;
            LOG << "lhs: " << newCutLhs << endl;
            LOG << "rhs: " << newCutRhs << endl;
         }
      }   

   if(cutLhs.size() > 0)
      return true;
   else 
      return false;
}

// Min-cut separation for DFJ
bool separate_min_cut_dfj(IloEnv masterEnv, const unordered_map<NODE_PAIR, double>& xSol, std::shared_ptr<Graph> G, 
      NODE_PAIR st, const ARC_VARS& sigma_var,
      vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation)
{
   // Build graph with x values as capacities
   SmartDigraph cap_graph;
   SmartDigraph::ArcMap<double> x_capacities(cap_graph);
   unordered_map<NODE,LemonNode> v_nodes;
   map<LemonNode,NODE> rev_nodes;
   NODE s = st.first;
   NODE t = st.second;
   build_cap_graph(cap_graph, x_capacities, v_nodes, rev_nodes, xSol, G, s, t);
   
   LOG << "Built graph" << endl;
   cutLhs = vector<IloExpr>();
   cutRhs = vector<IloExpr>();
   violation = vector<double>();
   IloExpr newCutLhs;
   IloExpr newCutRhs;
   double newViolation;
   double min_cut_value;
   for(NODE v : G->nodes())
      if(v!=s && v!=t)
      {
         Preflow<SmartDigraph, SmartDigraph::ArcMap<double>> min_cut(cap_graph, x_capacities, v_nodes[v], v_nodes[t]);
         min_cut.runMinCut();
         min_cut_value = min_cut.flowValue();

         LOG << "Ran min-cut" << endl;

         LOG << v << endl;
         LOG << "Min-cut " << min_cut_value << endl;

         //if mincut < degree(v), add A(S) <= card(S)-1)
         if(min_cut_value < 1)
         {
            newCutLhs = IloExpr(masterEnv);
            newCutRhs = IloExpr(masterEnv);
            double lhsValue = 0.;

            // for all nodes i in S
            int cardS = 0;
            for(SmartDigraph::NodeIt i(cap_graph); i!=INVALID; ++i)
               if(min_cut.minCut(i))
                  ++cardS;

            for(SmartDigraph::ArcIt e(cap_graph); e!=INVALID; ++e)
            {
               LemonNode i = cap_graph.source(e);
               LemonNode j = cap_graph.target(e);
               if(min_cut.minCut(i) && min_cut.minCut(j))
               {
                  newCutRhs += sigma_var.at(NODE_PAIR(rev_nodes[i],rev_nodes[j]));
                  lhsValue += xSol.at(NODE_PAIR(rev_nodes[i],rev_nodes[j]));
               }
            }

            newViolation = lhsValue - (cardS-1);
            newCutLhs += cardS - 1;

            if(newViolation >= TOL)
            {
               cutLhs.push_back(newCutLhs);
               cutRhs.push_back(newCutRhs);
               violation.push_back(newViolation);
            }
   
            LOG << "node " << v << endl;
            LOG << "violation " << newViolation << endl;
            LOG << "cut " << cutRhs.size() << endl;
            LOG << "lhs: " << newCutLhs << endl;
            LOG << "rhs: " << newCutRhs << endl;
         }
      }   

   if(cutLhs.size() > 0)
      return true;
   else 
      return false;
}


// Strong Component separation for DFJ
bool separate_sc_dfj(IloEnv masterEnv, const unordered_map<NODE_PAIR, double>& xSol, std::shared_ptr<Graph> G, 
      NODE_PAIR st, const ARC_VARS& sigma_vars, 
      vector<IloExpr>& cutLhs, vector<IloExpr>& cutRhs, vector<double>& violation)
{
   /** Strong component separation */
   bool ret = false;

   // Build support graph
   SmartDigraph support_graph;
   unordered_map<NODE,LemonNode> v_nodes;
   map<LemonNode,NODE> rev_nodes;
   NODE s = st.first;
   NODE t = st.second;
   build_support_graph(support_graph, v_nodes, rev_nodes, xSol, G, s, t);

   // search for strong components
   SmartDigraph::NodeMap<int> nodemap(support_graph);
   int components = stronglyConnectedComponents(support_graph, nodemap);

   cutLhs = vector<IloExpr>(components);
   cutRhs = vector<IloExpr>(components);
   violation = vector<double>(components);
   vector<double> out_degree(components,0);
   vector<double> lhsValue(components,0);

   vector<int> cardS(components, 0);
   for(SmartDigraph::NodeIt i(support_graph); i!=INVALID; ++i)
   {
      int comp = nodemap[i];
      if(0 == cardS[comp]++) //initialize lhs and increase cardinality
      {  LOG << "Initialized lhs and rhs" << endl;
         cutLhs[comp] = IloExpr(masterEnv);
         cutRhs[comp] = IloExpr(masterEnv);
      }
   }

   for(SmartDigraph::ArcIt e(support_graph); e!=INVALID; ++e)
   {
      LemonNode i = support_graph.source(e);
      LemonNode j = support_graph.target(e);
      if(nodemap[i] == nodemap[j]) /* same component */
      {
         int comp = nodemap[i];
         cutRhs[comp] += sigma_vars.at(NODE_PAIR(rev_nodes[i],rev_nodes[j]));
         lhsValue[comp] += xSol.at(NODE_PAIR(rev_nodes[i],rev_nodes[j]));
      }
   }

   for(int i=0;i<components;++i)
   {
      cutLhs[i] += cardS[i] - 1;
      violation[i] = lhsValue[i] - (cardS[i] - 1);
      LOG << "-- Strong Component " << i << ": " << cardS[i] << endl;
      LOG << "lhs: " << cutLhs[i] << endl;
      LOG << "rhs: " << cutRhs[i] << endl;
      LOG << "violation: " << violation[i] << endl;
      if(violation[i] >= TOL)
         ret = true;

   }

   return ret;
}


// Strong Component separation for MCF
bool separate_sc_mf(IloEnv masterEnv, const unordered_map<NODE_PAIR, double>& xSol, const std::shared_ptr<Graph> G,
      NODE_PAIR st,  ARC_VARS& sigma_vars, unordered_map<TRIPLET,IloNumVar>& qq_var, unordered_map<NODE,IloNumVar>& zz_var,
      vector<IloRange>& cons)
{
   /** Strong component separation */
   bool ret = false;
   
   // Build support graph
   SmartDigraph support_graph;
   unordered_map<NODE,LemonNode> v_nodes;
   map<LemonNode,NODE> rev_nodes;
   NODE s = st.first;
   NODE t = st.second;
   build_support_graph(support_graph, v_nodes, rev_nodes, xSol, G, s, t);

   // Search for strong components
   SmartDigraph::NodeMap<int> nodemap(support_graph);
   int components = stronglyConnectedComponents(support_graph, nodemap);

   cons = vector<IloRange>();
   double violation;
   vector<int> cardinality(components,0);
   vector<double> out_degree(components,0);
   vector<double> max_node_degree(components,0);
   vector<NODE> max_node(components,0);
   SmartDigraph::NodeMap<double> node_out_degree(support_graph,0);
   for(SmartDigraph::NodeIt i(support_graph); i!=INVALID; ++i)
   {
      LOG << "--------------- node " << rev_nodes[i] << endl;
      int comp = nodemap[i];
      ++cardinality[comp];
      
      /* compute the out degree of each node, while also updating delta(S) of its component */
      for (NODE j : G->outgoing_from(rev_nodes[i]))
      {
         node_out_degree[i] += xSol.at(NODE_PAIR(rev_nodes[i],j));
         if( v_nodes.count(j) == 0 || comp != nodemap[v_nodes[j]])
            out_degree[comp] += xSol.at(NODE_PAIR(rev_nodes[i],j));
      }

      LOG << "degree: " << node_out_degree[i] << endl;
      /* identify the maximally violated node in each component */
      if(node_out_degree[i] > max_node_degree[comp])
      {
         max_node_degree[comp] = node_out_degree[i];
         max_node[comp] = rev_nodes[i];
      }
   }

   for(int i=0;i<components;++i)
   {
      LOG << "component " << i << ": " << cardinality[i] << endl;
      LOG << "delta(S) " << out_degree[i] << endl;
      LOG << "max(delta(i)) " << max_node_degree[i] << endl;
      violation = max_node_degree[i] - out_degree[i];
      if(violation >= TOL)
      {
         LOG << "violation " << violation << endl;
         NODE h = max_node[i];
         IloRange con;
         // Ham1
         // qq[h,i,j] <= sigma[i,j] (s,t)
         for (auto& arc : G->arcs())
         {   
            con = IloRange(masterEnv, -IloInfinity, 0.0);
            con.setLinearCoef(qq_var[TRIPLET(h,arc)],1.0);
            con.setLinearCoef(sigma_vars[arc],-1.0);
            cons.push_back(con);
         }

         // Ham2
         // sum{(s,j) in A}qq[h,s,j] = z[h] (s,t)
         con = IloRange(masterEnv, 0.0, 0.0);
         con.setLinearCoef(zz_var[h],-1.0);
         for(auto j : G->outgoing_from(s)) 
            con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(s,j))],1.0);
         cons.push_back(con);

         // Ham3
         // sum{(i,s) in A}qq[h,i,s] = 0
         con = IloRange(masterEnv, 0.0, 0.0);
         for(auto j : G->incoming_to(s)) 
            con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(j,s))],1.0);
         cons.push_back(con);

         // Ham4
         // sum{(i,h) in A}q[h,i,h] = z[h]
         con = IloRange(masterEnv, 0.0, 0.0);
         con.setLinearCoef(zz_var[h],-1.0);
         for(auto j : G->incoming_to(h)) 
            con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(j,h))],1.0);
         cons.push_back(con);

         // Ham5
         // sum{(h,j) in A}q[h,h,j] = 0 
         con = IloRange(masterEnv, 0.0, 0.0);
         for(auto j : G->outgoing_from(h)) 
            con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(h,j))],1.0);
         cons.push_back(con);

         // Ham6
         // sum{(i,j) in A}qq[h,i,j] - sum{(j,i) in A}qq[h,j,i] = 0
         for (NODE node : G->nodes())
         {
            if (node != h && node != s)
            {
               con = IloRange(masterEnv, 0.0, 0.0);
               for(auto j : G->outgoing_from(node))
                  con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(node,j))],-1.0);
               for(auto j : G->incoming_to(node))
                  con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(j,node))],1.0);
               cons.push_back(con);
            }
         }

//         // Ham7
//         // sum{(t,j) in A}q[h,t,j] = 0  //don't see why needed TODO
//         con = IloRange(masterEnv, 0.0, 0.0);
//         for(auto j : G->outgoing_from(t))
//            con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(t,j))],1.0);
//         cons.push_back(con);

         // Ham8
         // sum{(i,h) in A}sigma[i,h] = z[h]
         // if node h is reached, z[h] = 1
         con = IloRange(masterEnv, 0.0, 0.0);
         con.setLinearCoef(zz_var[h],-1.0);
         for(auto j : G->incoming_to(h)) 
            con.setLinearCoef(sigma_vars[NODE_PAIR(j,h)],1.0);
         cons.push_back(con);

         /* if you want to separate only one violation */
         //return true;
         
         ret = true;
      }

   }

   return ret;
}

// Min-cut separation for MCF
bool separate_min_cut_mf(IloEnv masterEnv, const unordered_map<NODE_PAIR, double>& xSol, std::shared_ptr<Graph> G, 
      NODE_PAIR st, ARC_VARS& sigma_vars, unordered_map<TRIPLET,IloNumVar>& qq_var, unordered_map<NODE,IloNumVar>& zz_var, 
      vector<IloRange>& cons)
{
   // Build graph with x values as capacities
   SmartDigraph cap_graph;
   SmartDigraph::ArcMap<double> x_capacities(cap_graph);
   unordered_map<NODE,LemonNode> v_nodes;
   map<LemonNode,NODE> rev_nodes;
   NODE s = st.first;
   NODE t = st.second;
   build_cap_graph(cap_graph, x_capacities, v_nodes, rev_nodes, xSol, G, s, t);

   LOG << "Built graph" << endl;
   cons = vector<IloRange>();
   double min_cut_value;
   for(NODE v : G->nodes())
      if(v!=s && v!=t)
      {
         Preflow<SmartDigraph, SmartDigraph::ArcMap<double>> min_cut(cap_graph, x_capacities, v_nodes[v], v_nodes[t]);
         min_cut.runMinCut();
         min_cut_value = min_cut.flowValue();

         LOG << "Ran min-cut" << endl;
         // Compute out degree of v
         double node_out_degree = 0;
         for (NODE j : G->outgoing_from(v))
            node_out_degree += xSol.at(NODE_PAIR(v,j));

         LOG << v << endl;
         LOG << "Min-cut " << min_cut_value << endl;
         LOG << "Node degree " << node_out_degree << endl;

         //if mincut < degree(v), add delta+(S) >= delta+(v)
         if ( node_out_degree >= min_cut_value + TOL)
         {
            NODE h = v;
            IloRange con;
            // Ham1
            // qq[h,i,j] <= sigma[i,j] (s,t)
            for (auto& arc : G->arcs())
            {   
               con = IloRange(masterEnv, -IloInfinity, 0.0);
               con.setLinearCoef(qq_var[TRIPLET(h,arc)],1.0);
               con.setLinearCoef(sigma_vars[arc],-1.0);
               cons.push_back(con);
            }

            // Ham2
            // sum{(s,j) in A}qq[h,s,j] = z[h] (s,t)
            con = IloRange(masterEnv, 0.0, 0.0);
            con.setLinearCoef(zz_var[h],-1.0);
            for(auto j : G->outgoing_from(s)) 
               con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(s,j))],1.0);
            cons.push_back(con);

            // Ham3
            // sum{(i,s) in A}qq[h,i,s] = 0
            con = IloRange(masterEnv, 0.0, 0.0);
            for(auto j : G->incoming_to(s)) 
               con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(j,s))],1.0);
            cons.push_back(con);

            // Ham4
            // sum{(i,h) in A}q[h,i,h] = z[h]
            con = IloRange(masterEnv, 0.0, 0.0);
            con.setLinearCoef(zz_var[h],-1.0);
            for(auto j : G->incoming_to(h)) 
               con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(j,h))],1.0);
            cons.push_back(con);

            // Ham5
            // sum{(h,j) in A}q[h,h,j] = 0 
            con = IloRange(masterEnv, 0.0, 0.0);
            for(auto j : G->outgoing_from(h)) 
               con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(h,j))],1.0);
            cons.push_back(con);

            // Ham6
            // sum{(i,j) in A}qq[h,i,j] - sum{(j,i) in A}qq[h,j,i] = 0
            for (NODE node : G->nodes())
            {
               if (node != h && node != s)
               {
                  con = IloRange(masterEnv, 0.0, 0.0);
                  for(auto j : G->outgoing_from(node))
                     con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(node,j))],-1.0);
                  for(auto j : G->incoming_to(node))
                     con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(j,node))],1.0);
                  cons.push_back(con);
               }
            }

//            // Ham7
//            // sum{(t,j) in A}q[h,t,j] = 0  //don't see why needed TODO
//            con = IloRange(masterEnv, 0.0, 0.0);
//            for(auto j : G->outgoing_from(t))
//               con.setLinearCoef(qq_var[TRIPLET(h,NODE_PAIR(t,j))],1.0);
//            cons.push_back(con);

            // Ham8
            // sum{(i,h) in A}sigma[i,h] = z[h]
            // if node h is reached, z[h] = 1
            con = IloRange(masterEnv, 0.0, 0.0);
            con.setLinearCoef(zz_var[h],-1.0);
            for(auto j : G->incoming_to(h)) 
               con.setLinearCoef(sigma_vars[NODE_PAIR(j,h)],1.0);
            cons.push_back(con);

            /* if you want to separate only one violation */
            //return true;
         }
      }   

   if(cons.size() > 0)
      return true;
   else 
      return false;
}

