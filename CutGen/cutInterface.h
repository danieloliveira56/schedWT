#ifndef __CUT_INTERFACE_H__
#define __CUT_INTERFACE_H__

#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
#include <hash_map>
#else
#include <ext/hash_map>
#endif
#include <map>

#ifdef _WIN32
using std::hash;
using std::hash_map;
#else
using __gnu_cxx::hash;
using __gnu_cxx::hash_map;
#endif


typedef enum
{
      PROB_FCNF, PROB_FCNF_SS, PROB_TREE, PROB_PATH, PROB_TSP
} ProblemType;


// Problem class will be the class of problems we will be solving.
// This is to be very specific (in contrast with ProblemType)
// For instance for each problem we solve, we will have a new ProblemClass
// CVRP, TSP, TDTSP, etc.
// The idea is that we will have to write one routine for each problemClass to set their
// characteristics
// Even though this may repeat some code, the effort should not be big and with that if some new characteristic
// is to be included in the future, one would only need to add the implementation for each problem Class
// and the cut generation routines would take care of it immediately.
typedef enum
{
    PCLASS_UNKNOWN, PCLASS_CVRP, PCLASS_TDTSP, PCLASS_SDVRP
} ProblemClass;


typedef enum
{
   CUT_UNKNOWN, CUT_SCC, CUT_PATH_EQ, CUT_EXTCC, CUT_GENERIC, CUT_INDEGREE, CUT_ARCFIX
} CutType;

// when the function "GetArcCapCoeff" is called with the last argument equal
// to this constant, it returns non-zero if there is a non-zero value for the
// given arc with any demand.
const int AnyDemand = -1;

typedef struct
   {
      int capacity;
      int numNodes;
         // root node has index zero and other nodes range from 1 to numNodes-1
      int* demand;   // node demands

      // For problems like CVRP, specifies how many vehicles come out of the root
      int nrootbranches;

   } InstanceInfo;

typedef struct
   {
      int i;
      int j;
      double value;
   } CUTS_ArcVariable;

typedef struct
   {
      int i;
      int j;
      int d;
      double value;
   } CUTS_ArcCapVariable;

struct LpSolution
{
   CUTS_ArcVariable* arcs;   // all existing arc variables
   int numArcs;
   CUTS_ArcCapVariable* arcsCap;   // all existing arc-capacity variables
   int numArcsCap;
};


/*****************************************************************************************/
/*     Data structures used for GenericCutData                                           */
/*****************************************************************************************/

struct ArcHashKey
{
   int i;
   int j;
#ifdef _WIN32
   operator size_t() const
   {
      return (i + j*300);
   }
#endif
   bool operator< ( const ArcHashKey &arc ) const
   {
      if (i != arc.i) return (i < arc.i);
      return (j < arc.j);
   }
};

#ifndef _WIN32
struct eqArcKey
{
   bool operator() ( const ArcHashKey &arc1, const ArcHashKey &arc2 ) const
   {
      return( (arc1.i == arc2.i) && (arc1.j == arc2.j) );
   }
};

struct hashArcKey
{
   size_t operator() ( const ArcHashKey &arc ) const
   {
      hash<int> H;
      return( H(arc.i*1000000 + arc.j ) );
   }
};
#endif

struct ArcCapHashKey
{
   int i;
   int j;
   int d;
#ifdef _WIN32
   operator size_t() const
   {
      return (i + j*300 + d*90000);
   }
#endif
   bool operator< ( const ArcCapHashKey &arc ) const
   {
      if (i != arc.i) return (i < arc.i);
      if (j != arc.j) return (j < arc.j);
      return (d < arc.d);
   }
};

#ifndef _WIN32
struct eqArcCapKey
{
   bool operator() ( const ArcCapHashKey &arc1, const ArcCapHashKey &arc2 ) const
   {
      return( (arc1.i == arc2.i) && (arc1.j == arc2.j) && (arc1.d == arc2.d) );
   }
};

struct hashArcCapKey
{
   size_t operator() ( const ArcCapHashKey &arc ) const
   {
      hash<int> H;
      return( H(arc.i + arc.j*1000 + arc.d*1000000) );
   }
};
#endif

#ifdef _WIN32
typedef hash_map< ArcCapHashKey, double > ArcCapCoeffHash;
#else
typedef hash_map< ArcCapHashKey, double, hashArcCapKey, eqArcCapKey > ArcCapCoeffHash;
#endif
// typedef std::map< ArcCapHashKey, double > ArcCapCoeffHash;


/*****************************************************************************************/
/*****************************************************************************************/

struct StrenCCData
{
   std::vector<bool> inSet;
   int setDem;
   int capacity;
   int minDem;
};

struct ExtCCData
{
   int weightS;
   std::vector<bool> S;

   int capacity;
   std::vector<double> incapcoeff;
   std::vector<double> outcapcoeff;
};

struct GenericCutData
{
   ArcCapCoeffHash coeffs;
};


struct ArcFixCutData
{
   int i;
   int j;
};

struct InDegreeCutData
{
   int vertex;
};

struct Cut
{
      void* data;
      double rhs;
      char sense;    // '<', '=' or '>'
      double violation;
      double (*GetArcCoeff)( void* data, int i, int j );
      double (*GetArcCapCoeff)( void* data, int i, int j, int d );
      bool (*IsEqual)( void* data1, void* data2 );
      void* (*CopyData)( void* data );
      void (*DestroyData)( void* data );
      char templabel[100];
};


struct CutList
{
   Cut* cuts;
   int numCuts;
   int maxNumCuts;
   double diversity;
};

// Return a pointer to a new cut generator object
void* InitCutGenerator( InstanceInfo* instance );

// destroy a cut generator object
void DestroyCutGenerator( void* cutGenPointer );

//----- All other functions return an error code or 0 if success -----

int SeparateECCbyHeuristic( void* cutGen, LpSolution* sol, CutList* cuts,
      int maxCuts, ProblemType prob, int minSetSize );

int SeparateExtCyElimByHeur( void* cutGen, LpSolution* sol, CutList* cuts,
      int maxCuts, ProblemType prob );

int SeparateOECbyHeuristic( void* cutGen, LpSolution* sol, CutList* cuts,
      int maxCuts, ProblemType prob, int minSetSize );


void ClearCutList( CutList* list );

void ResetCut( Cut* cut );

#endif // __CUT_INTERFACE_H__
