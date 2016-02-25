
#ifndef DWMASTER_H
#define DWMASTER_H

/*
 * dwMaster.hpp
 */

#ifdef _WIN32
#include <OsiSolverInterface.hpp>
#else
#include <coin/OsiSolverInterface.hpp>
#endif
#include "PricingSolver.hpp"
#include "roundEccSep.h"
#include <string>

// tolerances
const double RedCostEps = 1e-10;
const double DwmValueEps = IntegerTolerance / 2.0;
const double DwmStabEps = 0.02;
const double DwmCoeffEps = CoeffTolerance;
const double CutCleanDualEps = 1e-5;
const double CutCleanSlackEps = 0.01;
const double ColCleanXEps = 1e-10;
const double ColExpandEps = COL_EXPAND_EPS;
const double SolCheckEps = 1e-8;
const double ObjCheckFactor = 1e-5;

class Instance;
class OsiSolverInterface;
class Model;
#ifdef CPX
class OsiCpxSolverInterface;
#endif
#ifdef CLP
class OsiClpSolverInterface;
#endif
#ifdef GLP
class OsiGlpkSolverInterface;
#endif

const std::string COL_PREFIX_LBDA = "lbda";

class DWMaster
{
public:
   // constructors
   DWMaster( Instance *inst, const double ub, PricingSolver* pricing,
         const double *_duals );
   DWMaster( DWMaster& other );


   // solve the LP relaxation
   void solveRelaxation( bool showLog );

   /**
    * Does the pricing, computing the best reduced cost (brcc)
    * and the respective route, if brc is negative, this route is
    * inserted into the LP and it returns true, otherwise returns zero.
    * If update is true, the constraints changed since the last pricing.
    **/
   bool doPricing( bool showLog, bool update, const double *_duals = 0 );

   /**
    * execution status
    **/
   bool isAbandoned() const;
   bool isProvenOptimal() const;
   bool isInfeasible() const;

   /**
    * getters/setters
    **/
   inline OsiSolverInterface *getSolver() { return solver; }
   inline int getNumCols() { return solver->getNumCols(); }
   inline Model* getModel() { return model; }
   inline int getArcCount() { return pricingSolver->getArcCount(); }
   inline double getUB() { return _ub; }
   inline void setUB(double ub) { _ub = ub; }
   inline int getMissPricings() { return missPricings; }
   inline int getStabChanges() { return stabChanges; }
   inline bool getIsArcTimeIndexed() { return pricingSolver->isArcTimeIndexed(); }

   /**
    * querying constraint and variable
    * indexes
    **/
   inline int getIdxDegreeConstraints() { return idxConsSelRouteJob; }
   inline int getIdxConstraintMachines() { return idxConsMachines; }
   inline int getIdxCQRouteVars() { return idxCQRouteVars; }

   /**
    * generate an LP file equivalent to the current problem
    **/
   inline void generateLpFile( int nodeNumber )
   {
      pricingSolver->updateConstraints( model );
      pricingSolver->generateLpFile( nodeNumber );
   }

   inline void generateTimeIndexedLpFile( int nodeNumber )
   {
      pricingSolver->updateConstraints( model );
      pricingSolver->generateTimeIndexedLpFile( nodeNumber );
   }

   inline void generateArcTimeFixDataFile( int nodeNumber )
   {
      pricingSolver->updateConstraints( model );
      pricingSolver->generateArcTimeFixDataFile( nodeNumber );
   }

   
   

   /**
    * insert rows in the DWM LP, assuming that the constraints are already
    * inserted in the model
    **/
   void insertRowsInLP();

   /**
    * reduced cost of columns in the LP
    */
   const double *getRCs() const;

   /**
    * returns the contents
    * of the idx-th qroute
    **/
   CapArcKey *getCQRoute( const int idx, int& length );

   /**
    * return dual costs from constraints in the model
    **/
   const double *getRowPrice() const;

   /**
    * return the primal solution from variables in the model
    **/
   const double *getColSolution() const;

   /**
    * remove unused rows (whose dual variables are zero)
    **/
   void removeUnusedRows();

   /**
    * remove unused columns (only those whose reduced costs are greater than
    * zero, NOT those whose primal values are zero)
    **/
   void removeUnusedCols();

   /**
    * destroy the current pricing object and create a new one using the current
    * dual variables to fix variables by reduced cost
    **/
   void changePricing();

   // Destructor
   virtual ~DWMaster();

private:
   void computeInfoQRoute( const double *_duals );
   // pricing solution information
   double* lhsQR;    // left-hand side of each constraint in qRoute
   std::vector<CapArcKey> cqRoute; // contents of the last qroute discovered

   double prsRC;
   double prsCost;
   std::vector<int>    ncIdx;   // indexes and
   std::vector<double> ncCoef;  // coefs for new column

   // best lower bound and corresponding dual vector found so far
   double bestLB;
   std::vector<double> bestDuals;

   /**
    * inserts into LP the newly discovered qRoute and
    * returns column name
    ***/
   int insertQRoute();

   // contents of qroutes
   std::vector< std::vector<CapArcKey> > cqRoutes;

   Instance *_inst;
   OsiSolverInterface *solver;
   void createConsMachines();
   void createConsSelQRouteJob();

   // constraint inidices
   int idxConsMachines;
   int idxConsSelRouteJob;
   int idxConsExtCapCuts;

   // costs and indices of artificial variables
   // related to jobs and cuts
   double *aCostJobsAndCuts;
   void setACostJobsAndCuts(double ub);
   int idxAVarsJobs;
   int idxAVarsCuts;
   int idxNextAVarCut;

   // artificial variable related to machines constraints
   double aCostMachines;
   int idxACostMachines;

   // artificial vars for jobs and cuts
   void createAVarsJobsAndCuts();
   int idxCQRouteVars;

   // remove invalid columns (whose arcs hav been fixed) from LP
   void removeInvalidCols();

   double _ub;

   Model *model;
   PricingSolver *pricingSolver;
#ifdef CLP
   OsiClpSolverInterface *clpSolver;
#endif
#ifdef CPX
   OsiCpxSolverInterface *cpxSolver;
#endif
#ifdef GLP
   OsiGlpkSolverInterface *glpSolver;
#endif

   // index of the next lambda variable
   int nextLbdaIdx;

   // index of the next cut
   int nextCutIdx;

   // pricing iteration counters
   int pricingIt;
   int previousFix;
   int missPricings;
   int stabChanges;

   // current stabilization factor
   double currStabFactor;

   // LB in the next LP clean-up
   double cleanUpLB;
};

#endif /* ifndef DWMASTER_H */

