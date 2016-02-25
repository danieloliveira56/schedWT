/** Header File:
 ** Branch-and-bound node data sructures
 **/

#ifndef _B_B_NODE_H_
#define _B_B_NODE_H_

#include "cutInterface.h"
#include "dwMaster.hpp"

class DWMaster;
class CPUTimer;
class PricingSolver;
class Instance;
struct CutList;

class BBNode
{
public:
   // constructor/destructor
   BBNode( Instance *inst, const double ub, PricingSolver* pricing,
         const double *_duals );
   BBNode( BBNode& other );
   ~BBNode();

   // solve the current node branching and recursively solving the child nodes
   // if necessary
   void solve( CPUTimer& t, double& secsLP, double& secsNode, int& cutRounds,
         bool doBranch = true );

   // getters/setters
   inline int getNodeCount() { return nodeCount; }
   inline double getUB() { return ub; }
   inline double getFirstLB() { return firstLB; }
   inline double getLB() { return dwm->getSolver()->getObjValue(); }
   inline int getFirstNumIter() { return firstNumIter; }
   inline int getFirstArcs() { return firstArcs; }
   inline int getRemainArcs() { return dwm->getArcCount(); }
   inline int getMissPricings() { return missPricings; }
   inline int getStabChanges() { return stabChanges; }
   void write_sol(char fname[]);

private:
   // Dantzig-Wolf Master LP for this node
   DWMaster *dwm;

   // pointer to the instance data
   Instance *inst;

   // current DWM solution expanded to arc and capacitated arc variables
   LpSolution sol;

   // indicates that "sol" is an integer solution
   bool isIntegerSol;

   // total number of nodes used in the subtree rooted at this node
   int nodeCount;

   // number of iterations of the DWM LP
   int iteration;

   // number of this node in the whole tree
   int nodeNumber;

   // current level in the branch-and-bound tree
   int level;

   // current upper bound
   double ub;

   // lower bound obtained by the first DWM LP, with no cuts
   double firstLB;

   // number of iterations to solve the first DWM LP, with no cuts
   int firstNumIter;

   // other statistics for the first DWM LP, with no cuts
   int missPricings;
   int stabChanges;

   // remaining arcs after solving the first DWM LP, with no cuts
   double firstArcs;

   // call the separation rotines to generate cuts in the DWM LP
   int GenerateCuts(int cutRound,  CPUTimer& t, double currLB);

   // expand the q-route LP solution into capacitated arc variable value
   void ExpandSol();

   // translate a list o cuts into DWM constraints
   void ProcessCutList( CutList* cutList );

   // find a job and a time for branching
   void setBranchingJobTime( int& job, int& jtime, bool& fathomLeft,
         bool& fathomRight, std::vector<bool>& jobSet );

   // find a set of arcs for branching
   void setBranchingArcs( std::vector< std::vector<bool> >& arcs,
         bool& fathomLeft, bool& fathomRight );
};

#endif // _B_B_NODE_H_
