/** Header File:
 ** Object the solves the master problem
 **/

#ifndef _MASTER_SOLVER_H_
#define _MASTER_SOLVER_H_

class Instance;
class Model;
class PricingSolver;

class MasterSolver
{
public:
   // Constructor
   // n is the number of nodes including the root
   // dists is the matrix of distances
   MasterSolver( Instance *inst, int ub );

   // Destructor
   ~MasterSolver();

   // change the pricing solver to the graph based one
   // - delete the old one but keep the dual variables
   void changePricing( int ub );

   // solve the problem exactly (using the volume algorithm
   void solve( int& ub, double minGap = 0.0 );

   // check if all arcs of a given solution are non-fixed
   void checkSolution( int* sol, int size );

   // get the array of dual variable values
   inline const double *getDuals() const { return duals; }

   // get the last lower bound
   inline double getLowerBound() { return lowerBound; }

   // set the array of dual variables
   inline void setDuals( const double* newDuals )
   {
      for (int i = 0; i < _inst->jobs(); i++)
         duals[i] = newDuals[i];
   }

   // get the pointer to the pricing solver and zeroize this pointer
   inline PricingSolver* getPricingSolver()
   {
      PricingSolver* ret = pricing;
      pricing = 0;
      return ret;
   }

private:
   Instance *_inst;

   // pointer to the model to be solved
   Model* model;

   // pointer to the pricing solver
   PricingSolver* pricing;

   // array with the current values of the dual variables
   double* duals;

   // last calculated lower bound when solving the relaxation
   double lowerBound;

   // fraction of initial GAP for the step size
   double stepSizeGapFrac;

   // change the current dual values following the subgradient direction
   bool subgradientStep( double* lhs, double* currLhs, double initialGap );
};

#endif // _MASTER_SOLVER_H_
