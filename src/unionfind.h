#ifndef __UNIONFIND_H__
#define __UNIONFIND_H__
#include <vector>
#include <unordered_map>

using namespace std;

template <class T>
struct UnionFind{
   unordered_map<T, T> parent; //parent for each node
   unordered_map<T, int> rank;  //rank for each node

   UnionFind(const vector<T> elements)
   {
      for(T i : elements)
         make_set(i);
   }

   void make_set(T i)
   {
      parent[i]=i;
      rank[i]=0;
   }

   T find_set(T i) 
   {
      if(i!=parent[i])
         parent[i] = find_set(parent[i]); //path compression
      return parent[i];
   }

   void join(T i, T j)	
   {
      i=find_set(i);
      j=find_set(j);
      if (i==j)
         return;
      if (rank[i]>rank[j])  //i longer. attach j to i
         parent[j]=i;
      else{
         parent[i]=j; //i not longer. attach i to j
         if(rank[i]==rank[j]) //rank grows if equal length
            rank[j]+=1;
      }
   }

};

#endif
