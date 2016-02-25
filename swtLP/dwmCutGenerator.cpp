#include "dwmCutGenerator.hpp"
#include "../Instance.hpp"
#include "dwMaster.hpp"
#ifdef _WIN32
#include <OsiSolverInterface.hpp>
#else
#include <coin/OsiSolverInterface.hpp>
#endif

using namespace std;

const double EPS = 1e-5;

DWMCutGenerator::DWMCutGenerator( Instance *_inst, DWMaster *_dwm )
   :
   cutGenerator(0),
   inst(_inst),
   dwm(_dwm),
   solCapArcs(0),
   solCapCapArcs(0)
{
   cutList.cuts = 0;
   cutList.numCuts = 0;

   sol.arcs       = 0;
   sol.arcsCap    = 0;
   sol.numArcs    = 0;
   sol.numArcsCap = 0;
}

int DWMCutGenerator::generateCuts()
{
   int result = 0;

   allocCutGenerator();

   ProblemType pt = PROB_PATH;

   printf("Starting cut generation...\n");
   printf("StrenCCsH = ");
   //nprevcuts = cutList.numCuts;
   //status = SeparateStrenCCbyHeuristic( cutGenPointer,&sol,&cutList,10, problemType,8);
   //numberSCCH = cutList.numCuts - nprevcuts;
   //printf("%d\n",numberSCCH);


   return result;
}

void DWMCutGenerator::allocCutGenerator()
{
   // cut generator initialization
   if (!cutGenerator)
   {
      int *demand = new int[ inst->jobs() + 2 ];
      demand[0] = 0.0;
      for ( int j=1 ; (j<=inst->jobs()) ; j++ )
         demand[j] = inst->ptime()[j];
      demand[inst->jobs()+1] = 0.0;

      InstanceInfo ii;
      ii.capacity      = inst->T();
      ii.numNodes      = inst->jobs();
      ii.nrootbranches = -1;
      ii.demand        = demand;

      cutGenerator = InitCutGenerator( &ii );

      delete[] demand;
   }

}

CUTS_ArcVariable *DWMCutGenerator::getSolArcVar( const int index )
{
   if ( index >= solCapArcs )
   {
      solCapArcs = max( 5000, ((int)(((double)solCapArcs)*1.5)) );
      size_t space = sizeof(CUTS_ArcVariable)*solCapArcs;

      if ( sol.arcs )
      {
         CUTS_ArcVariable *tmp = (CUTS_ArcVariable*) realloc( sol.arcs, space );
         if (!tmp)
         {
            fprintf(stderr, "No memory to store fractional solution.");
            exit(EXIT_FAILURE);
         }
         sol.arcs = tmp;
      }
      else
      {
         sol.arcs = (CUTS_ArcVariable*) malloc( space );
         if (!sol.arcs)
         {
            fprintf(stderr, "No memory to store fractional solution.");
            exit(EXIT_FAILURE);
         }
      }
   }

   return &(sol.arcs[index]);
}

CUTS_ArcCapVariable *DWMCutGenerator::getSolCapacitatedArcVar( const int index )
{
   if ( index >= solCapCapArcs )
   {
      solCapCapArcs = max( 5000, ((int)(((double)solCapCapArcs)*1.5)) );
      size_t space = sizeof(CUTS_ArcCapVariable)*solCapCapArcs;

      if ( sol.arcsCap )
      {
         CUTS_ArcCapVariable *tmp = (CUTS_ArcCapVariable*) realloc( sol.arcsCap, space );
         if (!tmp)
         {
            fprintf(stderr, "No memory to store fractional solution.");
            exit(EXIT_FAILURE);
         }
         sol.arcsCap = tmp;
      }
      else
      {
         sol.arcsCap = (CUTS_ArcCapVariable*) malloc( space );
         if (!sol.arcsCap)
         {
            fprintf(stderr, "No memory to store fractional solution.");
            exit(EXIT_FAILURE);
         }
      }
   }

   return &(sol.arcsCap[index]);
}

void DWMCutGenerator::getSolFromCurrentLP()
{
   sol.numArcs = 0;
   sol.numArcsCap = 0;

   OsiSolverInterface *solver = dwm->getSolver();
   const double *x = solver->getColSolution();
   const int cols  = solver->getNumCols();

   for ( int j=0 ; (j<cols) ; j++ )
   {
      if ( x[j] <= EPS )
         continue;

      const string cname = solver->getColName( j );
      if ( cname.find(COL_PREFIX_LBDA + string("_")) != cname.npos )
         continue;
      size_t sizePrefix = string(COL_PREFIX_LBDA + string("_")).size();
      size_t lenNum = cname.size() - sizePrefix + 1;
      string strIdx = cname.substr( sizePrefix, lenNum );
      int colIndex = atoi( strIdx.c_str() );
   }
}

DWMCutGenerator::~DWMCutGenerator()
{
   if (cutGenerator)
      DestroyCutGenerator( cutGenerator );
}
