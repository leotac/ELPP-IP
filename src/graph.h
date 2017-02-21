/**@file graph.h
 * @brief Graph structure 
 * @author Leonardo Taccari 
 */

#ifndef __GRAPH_H__
#define __GRAPH_H__

#include "type.h"
#include <map>
#include <vector>

class Graph
{
   public:
      Graph() {}

      void add_node(NODE n) { 
         nodes_m.push_back(n);
         out_adj_list_m[n] = vector<NODE>();
         in_adj_list_m[n] = vector<NODE>();
      }

      void add_arc(NODE i, NODE j) {
         out_adj_list_m[i].push_back(j);
         in_adj_list_m[j].push_back(i);
         arcs_m.push_back(NODE_PAIR(i,j));
      }

      const vector<NODE>& nodes() const { return nodes_m; }
      const vector<NODE_PAIR>& arcs() const { return arcs_m; }
      
      const map<NODE, vector<NODE>>& out_adj_list() const { return out_adj_list_m; }
      const vector<NODE>& outgoing_from(NODE i) const { return out_adj_list_m.at(i); }
      
      const map<NODE, vector<NODE>>& in_adj_list() const { return in_adj_list_m; }
      const vector<NODE>& incoming_to(NODE i) const { return in_adj_list_m.at(i); }

      unsigned num_nodes() { return nodes_m.size(); }
      unsigned num_arcs() { return arcs_m.size(); }
   
   private:
      vector<NODE>                       nodes_m;         /* array of nodes */
      vector<NODE_PAIR>                  arcs_m;          /* arcs */
      map<NODE, vector<NODE>>            out_adj_list_m;      /* adjacency list */
      map<NODE, vector<NODE>>            in_adj_list_m;      /* adjacency list */
};

#endif
