#include <unordered_map>
#include <queue>
#include "graph.h"
#include "boost/heap/fibonacci_heap.hpp"

typedef boost::heap::fibonacci_heap<std::pair<double,NODE>, boost::heap::compare<greater<pair<double,NODE>>>> pq_type;

template<typename T>
void erase_by_value(vector<T>& v,  T val) {
   v.erase(std::remove(v.begin(), v.end(), val), v.end());
}

void Graph::prune(NODE s, NODE t) {

   std::unordered_map<NODE, double> distance_from_s;
   std::unordered_map<NODE, NODE> parent_from_s;
   dijkstra(s, distance_from_s, parent_from_s, true);

   std::unordered_map<NODE, double> distance_to_t;
   std::unordered_map<NODE, NODE> parent_to_t;
   dijkstra(t, distance_to_t, parent_to_t, false);

   for(NODE n : nodes_m) {
      if(distance_from_s[n] < std::numeric_limits<double>::infinity() && distance_to_t[n] < std::numeric_limits<double>::infinity())
         continue;
      //For each v reachable from n, remove n from the incoming star
      for(NODE v : out_adj_list_m[n]) {
         erase_by_value(in_adj_list_m[v], n);
         erase_by_value(arcs_m, NODE_PAIR(n,v));
      }
      out_adj_list_m.erase(n);

      //For each v that n is reachable from, remove n from the outgoing star
      for(NODE v : in_adj_list_m[n]) {
         erase_by_value(out_adj_list_m[v], n);
         erase_by_value(arcs_m, NODE_PAIR(v,n));
      }
      in_adj_list_m.erase(n);
   }
  
   nodes_m.erase(std::remove_if(std::begin(nodes_m), std::end(nodes_m), [&distance_from_s, &distance_to_t](NODE n){ 
            if(distance_from_s[n] < std::numeric_limits<double>::infinity() 
               && distance_to_t[n] < std::numeric_limits<double>::infinity())
               return false;
            return true;
            }), nodes_m.end());
}

bool Graph::are_connected(NODE s, NODE t)
{
    std::unordered_map<NODE, bool> reached;
    for (NODE node : nodes_m)
        reached[node] = false;
    std::priority_queue<NODE> queue;

    NODE cur_node;
    queue.push(s);
    reached[s] = true;
    while ( !queue.empty() ) {
        cur_node = queue.top();
        queue.pop();
        for (NODE v : outgoing_from(cur_node))
            if ( reached[v] == false ) {
                queue.push(v);
                reached[v] = true;
            }
    }
    return reached[t];
}

/** Simple implementation of Dijkstra algorithm. */ 
void Graph::dijkstra(NODE origin,
      std::unordered_map<NODE, double>& distance,
      std::unordered_map<NODE, NODE>& parent,
      bool outgoing // true whether you want to compute shortest paths from s. otherwise, compute shortest paths to s.
      )
{
    std::unordered_map<NODE, pq_type::handle_type> pq_handle;
    std::unordered_map<NODE, bool> visited;

    pq_type pq;
    std::pair<double, NODE> item;
    NODE cur_node = origin;

    // Initialize structures for Dijkstra
    for (NODE node : nodes_m) {
        if(node == origin)
            distance[node] = 0;
        else
           distance[node] = std::numeric_limits<double>::infinity();
        parent[node] = node;
        visited[node] = false;
        pq_handle[node] = pq.push(std::make_pair(distance[node], node));
    }

    vector<NODE> adjacent_nodes;
    while ( !pq.empty() ) {
        item = pq.top();
        cur_node = item.second;
        visited[cur_node] = true;
        pq.pop();
        if(outgoing) 
            adjacent_nodes = outgoing_from(cur_node);
        else
            adjacent_nodes = incoming_to(cur_node);
        
        for (NODE v : adjacent_nodes)
        {  
            double arc_weight = 1;
            if ( visited[v] == false
                    && (distance[v] > distance[cur_node] + arc_weight))
            {
                distance[v] = distance[cur_node] + arc_weight;

                // Update distance with increase(), never with decrease()
                // (in shortest paths, we actually decrease it - but stdlib uses max-heap by default,
                // while we are using a min-heap so that we can use it)
                // TODO: make sure the update is in the correct direction
                // wrt the priority queue. Does it actually change the performance?
                pq.increase(pq_handle[v], std::make_pair(distance[v],v));
                parent[v] = cur_node;
            }
        }
    }
}   
