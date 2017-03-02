/**@file graph.h
 * @brief Graph structure 
 * @author Leonardo Taccari 
 */

#ifndef __GRAPH_H__
#define __GRAPH_H__

#include "type.h"
#include <map>
#include <unordered_map>
#include <vector>
#include <numeric>
#include <algorithm>

class Graph
{
   public:
      Graph() {}

      void add_node(NODE n) { 
         if(!has_node(n))
            nodes_m.push_back(n);
         if(out_adj_list_m.count(n) == 0)
            out_adj_list_m[n] = vector<NODE>();
         if(in_adj_list_m.count(n) == 0)
            in_adj_list_m[n] = vector<NODE>();
      }

      void add_arc(NODE i, NODE j) {
         out_adj_list_m[i].push_back(j);
         in_adj_list_m[j].push_back(i);
         arcs_m.push_back(NODE_PAIR(i,j));
      }
      
      bool has_node(NODE n) {
          return std::find(std::begin(nodes_m), std::end(nodes_m), n) != std::end(nodes_m);
      }

      const vector<NODE>& nodes() const { return nodes_m; }
      const vector<NODE_PAIR>& arcs() const { return arcs_m; }
      
      const map<NODE, vector<NODE>>& out_adj_list() const { return out_adj_list_m; }
      const vector<NODE>& outgoing_from(NODE i) const { return out_adj_list_m.at(i); }
      
      const map<NODE, vector<NODE>>& in_adj_list() const { return in_adj_list_m; }
      const vector<NODE>& incoming_to(NODE i) const { return in_adj_list_m.at(i); }

      unsigned num_nodes() { return nodes_m.size(); }
      unsigned num_arcs() { return arcs_m.size(); }
  
      bool check() { 
         // check number of nodes in adjacency lists is ok 
         if(out_adj_list_m.size() != nodes_m.size() || in_adj_list_m.size() != nodes_m.size()) {
            cout << "Check failed: number of nodes." << endl;
            return false;
         }
         
         // check number of arcs computed from the adjacency lists is ok
         auto count_arcs = [](unsigned res, const pair<NODE,vector<NODE>>& key_value_pair){ 
                            return res + key_value_pair.second.size(); };
         unsigned num_out_arcs = std::accumulate(out_adj_list_m.begin(), out_adj_list_m.end(), unsigned(0), count_arcs);
         unsigned num_in_arcs = std::accumulate(in_adj_list_m.begin(), in_adj_list_m.end(), unsigned(0), count_arcs);
         if( num_in_arcs != arcs_m.size() ||  num_out_arcs != arcs_m.size() ) {
            cout << "Check failed: number of arcs." << arcs_m.size() << num_out_arcs << num_in_arcs << endl;
            return false;
         }
         
         return true; 
      }

      void prune(NODE s, NODE t); //Prune nodes that canno be in a (s,t) path

      bool are_connected(NODE s, NODE t);
      void dijkstra(NODE root,
         std::unordered_map<NODE, double>& distance,
         std::unordered_map<NODE, NODE>& parent,
         bool outgoing = true // true whether you want to compute all the shortest paths *from* the root node.
                              // otherwise, compute shortest paths *to* the root node.
      );

   private:
      vector<NODE>                       nodes_m;         /* array of nodes */
      vector<NODE_PAIR>                  arcs_m;          /* arcs */
      map<NODE, vector<NODE>>            out_adj_list_m;      /* adjacency list */
      map<NODE, vector<NODE>>            in_adj_list_m;      /* adjacency list */
};

#endif
