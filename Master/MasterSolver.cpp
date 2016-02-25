/** Source Code File:
 ** Object the solves the master problem
 **/

#include "../Instance.hpp"
#include "MasterSolver.hpp"
#include "PricingSolver.hpp"
#include "Model.hpp"

#include <cstdio>
#include <cmath>

// fraction for updating the current solution
static const double maxCurrSolUpdFrac = 0.01;

// initial fraction of GAP for the step size
static const double initStepSizeGapFrac = 10.0;

// number of non-improments before reducing the step size GAP fraction
static const int nonImproveReduc = 20;

// factor of reduction for the step size GAP fraction
static const double gapFracReduc = 0.66;

// factor of increase for the step size GAP fraction
static const double gapFracInc = 1.1;

// minimum step size GAP fraction to continue
static const double minStepSize = 1E-6;

MasterSolver::MasterSolver( Instance *inst, int ub ) :
   _inst( inst ),
   lowerBound( -DoubleInfinity )
{
   // create the model and the pricing subproblem solver
   model = new Model( inst->jobs()+1 );
   if (inst->machines() > 1)
      pricing = new PricingSolverGraph( _inst, (const double*) NULL, ub );
   else
      pricing = new PricingSolverArray( inst );

   // create and zeroize the array of duals
   duals = new double[ inst->jobs() ];
   for ( int i = 0; i < inst->jobs() ; i++ )
      duals[i] = double(((i*104729)%inst->jobs())*ub)/double(inst->jobs()*inst->T());
}

MasterSolver::~MasterSolver()
{
   // release the model and the pricing subproblem solver
   delete model;
   if (pricing != 0) delete pricing;

   // release the array of duals
   delete [] duals;
}

void MasterSolver::changePricing( int ub )
{
   delete pricing;
   pricing = new PricingSolverGraph( _inst, duals, ub );
}

void MasterSolver::solve( int& ub, double minGap )
{
   // initializations
   int m = model->getNumConstraints();
   double* currLhs = new double[m];
   for ( int i = 0; i < m; i++ ) currLhs[i] = 0.0;
   double upperBound = ub;
   lowerBound = -DoubleInfinity;
   pricing->updateConstraints( model );

   // make a copy of the array of duals
   double* bestDuals = new double[m];
   for ( int i = 0; i < m; i++ ) bestDuals[i] = duals[i];

   // the Volume Algorithm main loop
   stepSizeGapFrac = initStepSizeGapFrac;
   int count = 0;
   int totalIter = 0;
   double initialGap = 0.0;
   bool good = false;
   double* lhs = new double[m];
   int* qroute = new int[ _inst->T()+1 ];
   do
   {
      // optimize the subproblem
      for ( int i = 0; i < m; i++ ) lhs[i] = 0.0;
      pricing->updateReducCosts( duals );

      double rc, cost_;
      pricing->solve( rc, cost_ );

      int aux;
      CapArcKey* solArcs = pricing->getSolution( aux );
      for ( int i = 0; i < aux ; i++ ) qroute[i] = solArcs[i].j;

      // fix arcs if its time
      if ( ((totalIter + count) % 5) == 0 )
      {
         int prevArcCount = pricing->getArcCount();
         pricing->fixArcs( upperBound, rc );
         if (prevArcCount > pricing->getArcCount())
            fprintf( stderr, "Fixed %d arcs out of %d.\n",
                  prevArcCount - pricing->getArcCount(), prevArcCount );
      }
      //printf( "Relaxed solution:" );
      //for (int i = 0; i < aux; i++) printf( " %d", qroute[i] );
      //printf( "\n" );

      // calculate the LHS's
      for ( int j = 0; j < m ; j++ ) lhs[j] = 0.0;
      for ( int i = 0; i < aux ; i++ )
         for ( int j = 0; j < m; j++ )
         {
            lhs[j] += ((double)_inst->machines()) * model->getConstraint( j )->getVarCoeff(
                         solArcs[i].i, solArcs[i].j );
            lhs[j] += ((double)_inst->machines()) * model->getConstraint( j )->getCapVarCoeff(
                         solArcs[i].i, solArcs[i].j, solArcs[i].d );
         }

      // check if the relaxed solution is feasible
      bool optimal_ = true;
      for ( int j = 0; j < m; j++ )
      {
         if ( fabs( lhs[j] - model->getConstraint( j )->getRhs() )
               > CoeffTolerance )
            optimal_ = false;
      }
      if ( optimal_ )
      {
         fprintf( stderr, "Optimal solution of %g found!\n", cost_ );
         for ( int i = 0; i < aux ; i++ ) fprintf( stderr, "%d ", qroute[i] );
         fprintf( stderr, "\n" );
         lowerBound = upperBound = ub = cost_;
         count++;
         delete [] solArcs;
         break;
      }

      // check if the best lower bound has been improved
      if ( ( ( cost_ -  lowerBound ) > ( upperBound - lowerBound ) * 0.00001 )
            || ( lowerBound == -DoubleInfinity ) )
      {
         // make a copy of the array of duals
         for ( int i = 0; i < m; i++ ) bestDuals[i] = duals[i];

         // do other updates
         totalIter += count + 1;
         count = 0;
         fprintf( stderr, "LB = %g <- %g (%c), %d iterations\n",
                 cost_, lowerBound, good? 'G': 'Y', totalIter );
         if ( good ) stepSizeGapFrac *= gapFracInc;
         lowerBound = cost_;
         if ( initialGap == 0.0 ) initialGap = upperBound - lowerBound;
      }
      else
      {
         // restore the array of duals
         for ( int i = 0; i < m; i++ ) duals[i] = bestDuals[i];

         // check if update the step size
         count++;
         if ( ( count % nonImproveReduc ) == 0 )
         {
            stepSizeGapFrac *= gapFracReduc;
            fprintf( stderr, "Step reduced to %g\n", stepSizeGapFrac );
         }
      }

      // update the lagrangian multipliers
      good = subgradientStep( lhs, currLhs, initialGap );

      delete [] solArcs;
   }
   while ( ( stepSizeGapFrac > minStepSize ) &&
           ( lowerBound + minGap <= upperBound - 1.0 + IntegerTolerance ) &&
           ( count < 200 * _inst->jobs() ) );

   // restore the array of duals
   for ( int i = 0; i < m; i++ ) duals[i] = bestDuals[i];

   // release the auxiliary arrays
   delete [] currLhs;
   delete [] bestDuals;
   delete [] lhs;
   delete [] qroute;
}

void MasterSolver::checkSolution( int* sol, int size )
{
   pricing->checkSolution( sol,size,true );
}

/*** PRIVATE METHODS ***/

bool MasterSolver::subgradientStep( double* lhs, double* currLhs,
                                    double initialGap )
{
   // initialize the current LHSs if necessary
   int m = model->getNumConstraints();
   bool currZero = true;
   for ( int j = 0; j < m; j++ )
      if ( currLhs[j] != 0.0 ) currZero = false;
   if ( currZero )
      for ( int j = 0; j < m; j++ ) currLhs[j] = lhs[j];

   // calculate the subgradient
   double improve = 0.0;
   double dividend = 0.0;
   double divisor = 0.0;
   double product = 0.0;
   for ( int j = 0; j < m; j++ )
   {
      double rhs = model->getConstraint( j )->getRhs();
      double viol = currLhs[j] - rhs;
      double viol2 = lhs[j] - rhs;
      dividend += viol * ( viol - viol2 );
      divisor += ( viol - viol2 ) * ( viol - viol2 );
      product += viol * viol2;
      improve += viol * viol;
   }

   // calculate the step size
   double stepSize = stepSizeGapFrac * initialGap / improve;

   // update the current objective and the previous step
   for ( int j = 0; j < m; j++ )
   {
      double rhs = model->getConstraint( j )->getRhs();
      duals[j] += ( rhs - currLhs[j] ) * stepSize;
   }

   // combine with the current solution (for the Volume algorithm)
   double currSolUpdateFrac = maxCurrSolUpdFrac;
   for ( int j = 0; j < m; j++ )
   {
      currLhs[j] *= 1.0 - currSolUpdateFrac;
      currLhs[j] += lhs[j] * currSolUpdateFrac;
   }
   return ( product >= 0 );
}
