#ifndef __CUT_GENERATOR_H__
#define __CUT_GENERATOR_H__

#include "cutInterface.h"


#include <stdio.h>
#include <string.h>

#include <vector>

class ProbCharact
{
public:
   ProbCharact();

   int has_root;
   int is_fcnf;

   int min_root_indeg;
   int max_root_indeg;

   int min_root_outdeg;
   int max_root_outdeg;

   int min_nonroot_indeg;
   int max_nonroot_indeg;

   int min_nonroot_outdeg;
   int max_nonroot_outdeg;

   int arcsin_equal_arcsout;

};


// Object for generating cuts
class CutGenerator
{
public:
   // Constructor
   CutGenerator( InstanceInfo* instance );

   // Destructor
   ~CutGenerator();

   // set the tolerance for the variable values
   inline void setIntEps( double value ) { intEps = value; }

   // set the LP solution pointer
   inline void setLpSolution( LpSolution* sol ) { solution = sol; }

   // set problem type
   inline void setProbType( ProblemType pType ) { probType = pType; }

   // set problem class
   int setProbClass( ProblemClass pClass );

   // set the maximum number of cuts generated in each batch
   inline void setCutBatch( int value ) { cutBatch = value; }

   // set the maximum depth to be reached in the subgraph enumeration
   inline void setMaxDepth( int value ) { maxDepth = value; }

   // set the maximum number of recursive calls during the subgraph enumeration
   inline void setMaxCalls( int value ) { maxCalls = value; }

   // set the minimum set size for a generated cut
   inline void setMinSetSize( int value ) { minSetSize = value; }

   /** EXTENDED CAPACITY CUT ROUTINES **/

   // generate Extended Capacity cuts (ECCs) by heuristic
   void extCapCutGenByHeur( CutList* cuts );

   // generate a single ECC from a set of vertices and r = a/b multiplier
   bool genSingleECCCut(CutList* cuts, int num, int den, int subtr,
         int *SetList, int SetListSize, double minViolation);

   // generate a single Triclique Cut from a set of 3 vertices
   bool genSingleTriCliqueCut( CutList* cuts, int S_i, int S_j, int S_k, double minViolation );
   
   // generate a single Generic Clique Cut from the set of x i,j,d indexs
   bool genGenericCliqueCut( CutList* cuts, std::vector<int> x_i, std::vector<int> x_j, std::vector<int> x_d, double minViolation );

   // generate Overload Elimination Cuts (OECs) by heuristic
   void extOECGenByHeur( CutList* cuts );   

   // generate a single ECC from a set of vertices and time t
   bool genSingleOECCut( CutList* cuts, int t, int m, 								 
      int *SetList, int SetListSize, double minViolation );

   bool genSingleECCCut_s( CutList* cuts, int num, int den, int s_num,
      int subtr, int *SetList, int SetListSize, double minViolation);

   bool isExtCapCutValid();

private:

   // auxiliary arithmetic functions
   double eccfrac(long long num, int den );
   int eccceil(long long num, int den );

   double eccfrac_s(long long num, int den, double s );

   /** EXTENDED CYCLE ELIMINATION CUT ROUTINES **/

public:

   // generate Extended Cycle Elimination Cuts by heuristic
   void extCyElimCutGenByHeur( CutList* cuts );

   bool isExtCyElimValid();

   /** OTHER AUXILIARY ROUTINES **/

private:
   void copyInstance( InstanceInfo* instance );

   // Convert the vertex index from the standard convention to the
   // convention used in the SCCs and ECCs
   inline int vertexToGsec( int vertex )
   { return ((vertex == 0)? instance_.numNodes: vertex); }

   // Convert the vertex index from the convention used in the SCCs and ECCs
   // to the standard convention
   inline int gsecToVertex( int vertex )
   { return ((vertex == instance_.numNodes)? 0: vertex); }

   // Check if a given generic cut is violated
   bool isViolated( Cut* cut, double tolerance );

   // Compute the lhs of the cut (activity)
   double computeActivity( Cut* cut );

   // Tolerance for the variable values
   double intEps;

   // Instance data
   InstanceInfo instance_;

   // current LP solution pointer
   LpSolution* solution;

   // problem type assumed when generating the cuts
   ProblemType probType;

   // Problem characteristics
   ProbCharact probProperties;

   // Problem class
   ProblemClass probClass;

   // Maximum number of cuts generated in each batch
   int cutBatch;

   // Maximum depth to be reached in the subgraph enumeration
   int maxDepth;

   // Maximum number of recursive calls during the subgraph enumeration
   int maxCalls;

   // Minimum set size for a generated cut
   int minSetSize;

   // Next vertex to start the set construction for ECCs
   int nextEccVertex;
};



double ExtCC_GetArcCoeff( void* data, int i, int j );

double ExtCC_GetArcCapCoeff( void* data, int i, int j, int d );

double InDegree_GetArcCoeff( void* data, int i, int j );

double InDegree_GetArcCapCoeff( void* data, int i, int j, int d );

double ArcFix_GetArcCoeff( void* data, int i, int j );

double ArcFix_GetArcCapCoeff( void* data, int i, int j, int d );

#endif // __CUT_GENERATOR_H__
