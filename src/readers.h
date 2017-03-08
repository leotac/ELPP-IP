/**@file   readers.h
 * @brief  File readers 
 * @author Leonardo Taccari
  */

/* standard library includes */
#include <iostream>
#include <fstream>
#include <string>
#include <memory>

/* user defined includes */
#include "type.h"
#include "graph.h"

/* namespace usage */
using namespace std;

class Readers
{
    public:
        static bool read_graph(std::string filename, std::shared_ptr<Graph> G, unordered_map<NODE_PAIR, double>& cost)
        {
            /* Read data file */
            ifstream infile(filename);
            if (!infile.good())
            {
                cout << "Can't read file " << filename <<endl;
                return false;
            }

            /* suffix */
            size_t lastdot = filename.find_last_of(".");
            string suffix = filename.substr(lastdot + 1, string::npos);

            if (suffix == "dat")
                return read_datfile(infile, G, cost);
            if (suffix == "dimacs")
                return read_dimacsfile(infile, G, cost);
            if (suffix == "graph")
                return read_graphfile(infile, G, cost);

            return false;
        };

    private:
        static bool read_datfile(std::ifstream& infile, std::shared_ptr<Graph> G, unordered_map<NODE_PAIR, double>& cost)
        {
            int n, m;
            infile >> n >> m;

            for(int k=0;k<n;k++)
            {
                NODE i;
                infile >> i;
                if(G->has_node(i)) {
                    cout << "Node " << i << " added twice." << endl;
                    return false;
                }
                G->add_node(i);
            }

            for(int k=0;k<m;k++)
            {
                NODE i,j;
                double c;
                infile >> i >> j >> c;

                if(!G->has_node(i) || !G->has_node(j)) {
                    cout << "Arc between non-existing nodes (either " << i << " or " << j << ")." << endl;
                    return false;
                }

                if(i == j) {
                    cout << "Self-loop detected." << endl;
                    return false;
                }

                G->add_arc(i,j);

                if(cost.count(NODE_PAIR(i,j)) > 0) {
                    cout << "Arc (" << i << "," << j << ") added twice." << endl;
                    return false;
                }
                cost[NODE_PAIR(i,j)] = c;
            }

            return G->check();

        };

        static bool read_graphfile(std::ifstream& infile, std::shared_ptr<Graph> G, unordered_map<NODE_PAIR, double>& cost)
        {
           // Format used in the 10th DIMACS Implementation Challenge
           string line;
           int n, m;
           int weight_type = 0;
           bool has_node_weights = false, has_arc_weights = false; 
           std::getline(infile, line);
           std::stringstream ss(line);
           ss >> n >> m >> weight_type;

           if(weight_type == 1)
              has_arc_weights = true;
           if(weight_type == 10)
              has_node_weights = true;
           if(weight_type == 11) {
              has_node_weights = true;
              has_arc_weights = true;
           }

           int j;
           double c;
           for(int i=1; i<=n; ++i)
           {
              G->add_node(i);
              std::getline(infile, line);
              std::stringstream ss(line);
              
              if(has_node_weights) //skip node weight, if there
                 ss >> j;

              while(ss >> j)
              {
                 cout << j << "-";
                 cout << j << endl;

                 if( j > n || j < 1) {
                    cout << "Node " << j << " out of bounds." << endl;
                    return false;
                 }

                 if(i == j) {
                    cout << "Self-loop " << i << " detected." << endl;
                    return false;
                 }

                 G->add_arc(i,j);

                 if(cost.count(NODE_PAIR(i,j)) > 0) {
                    cout << "Arc (" << i << "," << j << ") added twice." << endl;
                    return false;
                 }

                 if(has_arc_weights) {
                    ss >> c; 
                    cost[NODE_PAIR(i,j)] = c;
                 }
                 else
                    cost[NODE_PAIR(i,j)] = 1;
              }
           }
           return G->check();
        };
        
        static bool read_dimacsfile(std::ifstream& infile, std::shared_ptr<Graph> G, unordered_map<NODE_PAIR, double>& cost){
           // Format used in the 9th DIMACS Implementation Challenge
           string line;
           int n, m;
           char line_type;
           string dummy;
           bool read_problem_line = false;
           while(std::getline(infile, line)) {
              std::stringstream ss(line);
               ss >> line_type;
               if(line_type == 'c')
                  continue;
               if(line_type == 'p') {
                  ss >> dummy >> n >> m;
                  read_problem_line = true;
               }
               if(!read_problem_line && line_type == 'a') {
                  cout << "Cannot have an arc line before the problem line has been read." << endl;
                  return false;
               }
               if(read_problem_line && line_type == 'a') {
                  int i,j;
                  double c;
                  ss >> i >> j >> c;
                  G->add_node(i);
                  G->add_node(j);
                  G->add_arc(i,j);
                  cost[NODE_PAIR(i,j)] = c;
               }
           }
           if(G->num_nodes() != n || G->num_arcs() != m) {
              //cout << G->num_nodes() << endl;
              //cout << n << endl;
              //cout << G->num_arcs() << endl;
              //cout << m << endl;
              cout << "Warning: the number of read nodes or arcs does not match the problem line." << endl;
           }
           return G->check();
        };

};
