/** Header File:
 ** Object the solves the pricing subproblem
 **/

#ifndef _PRICING_SOLVER_H_
#define _PRICING_SOLVER_H_

#include <vector>

const double DualTolerance = 1E-9;
const double CoeffTolerance = 1E-11;
const double IntegerTolerance = 1E-3;

const double DoubleInfinity = 1E+9;

struct CapArcKey
{
   int i;
   int j;
   int d;
};

struct CapArcKeyComp {
   bool operator() (const CapArcKey& lhs, const CapArcKey& rhs) const
   {
      if (lhs.i != rhs.i) return (lhs.i < rhs.i);
      if (lhs.j != rhs.j) return (lhs.j < rhs.j);
      return (lhs.d < rhs.d);
   }
};

class Instance;
class Model;
class ArcNode;
class CapArcNode;
class ModCapArcNode;
class ConstrShadow;
class ModCapArcShadow;

class PricingSolver
{
public:
   // Constructor/Destructor
   PricingSolver() {}
   virtual ~PricingSolver() {}

   // Create a copy of this
   virtual PricingSolver* copy() = 0;

   // update the data structures for speeding up the reduced cost calculations
   virtual void updateConstraints(Model* model) = 0;

   // update the reduced costs of the original variables
   virtual void updateReducCosts(const double* duals) = 0;

   /**
    * solve the subproblem: return the obtained lower bound
    **/
   virtual void solve( double &rc, double &lowerBound ) = 0;

   /**
    * try to fix all modulus arc variables to zero
    * remRC: reduced cost for ONE machine
    * without any fixation
    **/
   virtual void fixArcs( const double ub, const double remRC ) = 0;

   // get the optimal subproblem solution as an array of capacitated arc keys
   virtual CapArcKey* getSolution(int& size) = 0;

   // check if all arcs of a given solution are non-fixed
   virtual bool checkSolution(int* sol, int size, bool mustCheck) = 0;

   // search for a better feasible solution using the non-fixed arcs
   // -> the number of arcs in this solution is always numJobs + numMachines
   virtual double searchFeasible( const double ub, CapArcKey* sol ) = 0;

   // generate an LP file equivalent to the current problem
   virtual void generateLpFile( int nodeNumber ) = 0;

   // generate an LP file equivalent to the current problem
   virtual void generateTimeIndexedLpFile( int nodeNumber ) = 0;

   // pure virtual method implementation
   virtual void generateArcTimeFixDataFile( int nodeNumber ) = 0;

   // whether the formulation is arc-time indexed
   virtual bool isArcTimeIndexed() = 0;

   // Getters/Setters
   inline int getArcCount() { return arcCount; }

protected:

   // counter for non-fixed dynamic programming arcs
   int arcCount;
};

struct SolArrayDescr
{
    double cost[2];
    int job[2];
    int prev[2];
};

// Main pricing solver where dynamic programming transitions are stored in a
// graph.
class PricingSolverGraph : public PricingSolver
{
public:
   /**
    * Constructor
    * duals can ONLY be used when solving one machine problem
    **/
   PricingSolverGraph( Instance *inst, const double* duals, double ub );

   // Alternative contructor for copied objects
   PricingSolverGraph( Instance *inst, ArcNode* root_,
         CapArcNode* capStart_, CapArcNode* capEnd_,
         ModCapArcNode* modCapStart_, ModCapArcNode* modCapEnd_ );

   // Destructor
   ~PricingSolverGraph();

   // pure virtual method implementation
   PricingSolver* copy();

   // pure virtual method implementation
   void updateConstraints(Model* model);

   // pure virtual method implementation
   void updateReducCosts(const double* duals);

   // pure virtual method implementation
   virtual void solve( double &rc, double &lowerBound );

   // pure virtual method implementation
   void fixArcs( const double ub, const double remRC );

   // pure virtual method implementation
   CapArcKey* getSolution(int& size);

   // pure virtual method implementation
   bool checkSolution(int* sol, int size, bool mustCheck);

   // pure virtual method implementation
   double searchFeasible( const double ub, CapArcKey* sol );

   // pure virtual method implementation
   void generateLpFile( int nodeNumber );

   // pure virtual method implementation
   void generateTimeIndexedLpFile( int nodeNumber );

   
   // pure virtual method implementation
   void generateArcTimeFixDataFile( int nodeNumber );

   

   // pure virtual method implementation
   bool isArcTimeIndexed() { return true; }

private:
   Instance *_inst;

   // entry/exit points for the graph of original variables
   ArcNode* root;
   CapArcNode* capStart;
   CapArcNode* capEnd;
   ModCapArcNode* modCapStart;
   ModCapArcNode* modCapEnd;

   // array of constraint descriptors for speeding up the reduced cost
   // calculations
   ConstrShadow* constraints;
   int numConstr;

   // constant term of the lower bound
   double constLB;

   // counter for searches in the original variable graph
   int searchCount;

   // array of arcs to be deleted
   std::vector<ModCapArcShadow*> toBeDel;

   // check if and arc (i,j,t) exists, considering the preprocessing rule with
   // tie-breaking
   inline bool checkArc(int i, int j, int t);

   // fill the forward and backward lower bound arrays
   void updateLowerBounds(double** fwd, double** bwd,
         const double* duals);

   // solve the subproblem for the modulus node passed as an argument
   long double recursiveSolve(ModCapArcNode* start);

   // try to fix arcs until the modulus node passed as an argument
   long double recursiveFixArcs( ModCapArcNode* end, const double ub, const double remRC );
};

// Auxiliary pricing solver implementing a simpler dynamic programming
// algorithm that uses an array of states it assumes that:
// - the model contains only the standard degree constraints (all vertices
//   but the root have degree constraints equal to one)
// - there is only one machine
class PricingSolverArray : public PricingSolver
{
public:
   // Constructor/Destructor
   PricingSolverArray( Instance *inst );
   ~PricingSolverArray();

   // pure virtual method implementation
   PricingSolver* copy();

   // pure virtual method implementation
   void updateConstraints(Model* model) {}

   // pure virtual method implementation
   void updateReducCosts( const double* duals );

   // pure virtual method implementation
   virtual void solve( double &rc, double &lowerBound );

   // pure virtual method implementation
   void fixArcs( const double ub, const double remRC );

   // pure virtual method implementation
   CapArcKey* getSolution(int& size);

   // pure virtual method implementation
   bool checkSolution(int* sol, int size, bool mustCheck);

   // pure virtual method implementation
   double searchFeasible( const double ub, CapArcKey* sol );

   // pure virtual method implementation
   void generateLpFile( int nodeNumber );
   
   // pure virtual method implementation
   void generateTimeIndexedLpFile( int nodeNumber );

   // pure virtual method implementation
   void generateArcTimeFixDataFile( int nodeNumber );

   // pure virtual method implementation
   bool isArcTimeIndexed() { return false; }

private:
   Instance* _inst;

   // array of dual variables (assuming one per vertex)
   double* _duals;

   // array of solutions (2 per capacity)
   SolArrayDescr* solution;

   // constant term of the lower bound
   double constLB;
};

#endif // _PRICING_SOLVER_H_
