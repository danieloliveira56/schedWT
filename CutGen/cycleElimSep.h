#include "cutInterface.h"

#include <vector>

#define ECEC_VIOLATED_EPS      0.003

class ExtCycleEliminator
{
   static const double Infinity;

   static const double Tolerance;

   static const double MinViolation;

   struct Vertex
   {
      // index of the capacitated arc in the LP solution
      int index;

      // a mark number for BFS
      int mark;

      // index of the previous vertex in the BFS
      int previous;

      // maximum flow passed to ths vertex by the BFS
      double maxFlow;

      // vector of adjacent output vertices
      std::vector<int> outVertices;

      // vector of indices of returning arcs in the output vertices
      std::vector<int> outReturns;

      // vector of capacities for output arcs
      std::vector<double> outCaps;
   };

public:
   ExtCycleEliminator( LpSolution* solution, int* demands, int capacity,
         int* setList, int setSize );

   // return true if the generated cut is violated
   // the vectors "heads", "tails" and "demands" describe the non-saturated
   // entering arcs and all compatible. The generated cut is the sum of
   // non-saturated arcs less than or equal to the sum of compatible leaving
   // arcs.
   bool generateCut( std::vector<int>& heads,  std::vector<int>& tails,
         std::vector<int>& demands, std::vector<double>& coeffs );

private:
   // input data for problem description
   LpSolution* solution_;

   // set of vertex for the cut generation
   int* setList_;
   int setSize_;

   // number of vertices
   int n;

   // vector of graph vertices to solve the minimum cut problem
   std::vector<Vertex> vertices;

   // index of the first leaving arc in the network
   int firstLeaving;

   // mark number used in the BFS that finds augmenting paths
   int currentMark;

   // minimum cut allowaed for the network (otherwise an ECE cut is violated)
   double minCut;

   // Vector of booleans that indicate the possible differences of
   // capacities of entering arcs and capacities of leaving arcs
   // diffCap[i][j][d] = true means that the capacity of an arc entering in
   // setList[i] may differ from that of an arc leaving setList[j] by d.
   // Since diffCap[i][j][x] = diffCap[j][i][x], we only initialize diffCap
   // for i <= j.
   std::vector< std::vector< std::vector<bool> > > diffCap;

   // send flow through an augmenting path from s to t
   // return the amount of flow sent
   double sendFlow();
};
