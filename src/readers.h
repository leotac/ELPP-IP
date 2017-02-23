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
            size_t lastdot = filename.find_first_of(".");
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
            cout << n << " " << m << endl;

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

        static bool read_graphfile(std::ifstream& infile, std::shared_ptr<Graph> G, unordered_map<NODE_PAIR, double>& cost){};
        static bool read_dimacsfile(std::ifstream& infile, std::shared_ptr<Graph> G, unordered_map<NODE_PAIR, double>& cost){};

};
