#ifndef __TYPE_H__
#define __TYPE_H__

#include <string>
#include <vector>
#include <array>

#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wredundant-decls"
#include <ilcplex/ilocplex.h>
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop
#pragma GCC diagnostic pop


using namespace std;

typedef int NODE;
typedef pair<NODE, NODE> NODE_PAIR;
typedef pair<NODE, NODE_PAIR> TRIPLET;
typedef pair<NODE_PAIR, NODE_PAIR> QUADRUPLET;
typedef pair<QUADRUPLET, NODE> QUINTUPLET;

enum ElppForm{
      NONE     = 0,
      MCF      = 1,
      SC       = 2,
      MinCut   = 3,
      SF       = 4,
      RLT      = 5,
      MTZ      = 6,
      DL       = 7,
      DFJ      = 8,
      MCFsep   = 9
};

static std::vector<std::string> ElppFormulationName = { "None", "MCF", "SC", "MinCut", "SF", "RLT", "MTZ", "DL", "DFJ", "MCFsep" };

template <class T1, class T2>
std::ostream& operator<<(std::ostream& out, const std::pair<T1, T2>& value)
{
   out << value.first << ' ' << value.second;
   return out;
}

/** Hash function for type pair< , > */
namespace std
{
template<typename a, typename b>
struct hash< pair<a, b> > {
private:
   const hash<a> ah;
   const hash<b> bh;
public:
   hash() : ah(), bh() {}
   size_t operator()(const pair<a, b> &p) const {
      size_t seed = ah(p.first);
      return bh(p.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);
   }
};
}


#endif
