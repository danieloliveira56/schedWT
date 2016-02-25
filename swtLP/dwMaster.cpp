/*
 * dwMaster.cpp
 */

// External includes
#include <algorithm>
#include <vector>
#include <string>
#include <cstring>
#include <cassert>
#include <cmath>
#include <map>
#ifdef _WIN32
#include <OsiSolverInterface.hpp>
#else
#include <coin/OsiSolverInterface.hpp>
#endif
#ifdef CLP
#ifdef _WIN32
#include <OsiClpSolverInterface.hpp>
#else
#include <coin/OsiClpSolverInterface.hpp>
#endif
#endif
#ifdef CPX
#ifdef _WIN32
#include <OsiCpxSolverInterface.hpp>
#else
#include <coin/OsiCpxSolverInterface.hpp>
#endif
#include <cplex.h>
#endif
#ifdef GLP
#ifdef _WIN32
#include <OsiGlpkSolverInterface.hpp>
#else
#include <coin/OsiGlpkSolverInterface.hpp>
#endif
#endif
// Local includes
#include "dwMaster.hpp"
#include "../Instance.hpp"
#include "../Pricing/PricingSolver.hpp"
#include "../Pricing/ConstrShadow.hpp"
#include "../Pricing/OrigVariable.hpp"
#include "../Master/Model.hpp"

//#define NO_STABILIZATION
//#define NO_FIXING

// minimum, maximum and step for the stabilization factor maintanance
const double MaxStabFactor = 0.1;
const double MinStabFactor = 0.1;
const double StepStabFactor = 1.33;

// maximum total number of cuts added
const int MaxTotalCuts = 10000;

// maximum ratio between the number of columns and the number of rows
const int MaxColsRowsRatio = 6;

// maximum numbe of arcs to run the heuristic
const int MaxArcHeuristic = 100000;

using namespace std;

DWMaster::DWMaster( Instance *inst, const double ub, PricingSolver* pricing,
      const double *_duals )
   :
     lhsQR( new double[inst->jobs()] ),
     _inst(inst),
     _ub(ub),
     model( new Model(inst->jobs()+1) ),
     nextLbdaIdx(0),
     nextCutIdx(0),
     pricingIt(0),
     missPricings(0),
     stabChanges(0)
{

#ifdef CLP
   clpSolver = new OsiClpSolverInterface();
   clpSolver->messageHandler()->setLogLevel(0);
   clpSolver->getModelPtr()->setLogLevel(0);
   clpSolver->getModelPtr()->scaling( 2/* geometric */);
   solver = dynamic_cast< OsiSolverInterface* >( clpSolver );
#endif
#ifdef GLP
   glpSolver = new OsiGlpkSolverInterface();
   solver = dynamic_cast< OsiSolverInterface* >( cpxSolver );
#endif
#ifdef CPX
   cpxSolver = new OsiCpxSolverInterface();
   solver = dynamic_cast< OsiSolverInterface* >( cpxSolver );
   CPXENVptr cpxEnv = cpxSolver->getEnvironmentPtr();
   CPXsetdblparam(cpxEnv, CPX_PARAM_EPRHS, 1e-8);
#endif
   solver->messageHandler()->setLogLevel(0);
   solver->messageHandler()->setPrefix( false );
   solver->setHintParam(OsiDoReducePrint,true,OsiHintTry);

   // to store column names in the LP
   solver->setIntParam(OsiNameDiscipline, 1);

   aCostJobsAndCuts = new double[ inst->jobs() + MaxTotalCuts ];
   setACostJobsAndCuts(ub);

   aCostMachines = (ub + 1.0) * 100.0;
   if (_inst->machines() > 1)
      createConsMachines();

   createConsSelQRouteJob();
   createAVarsJobsAndCuts();

   //solver->setHintParam(OsiDoScale);
   solver->initialSolve();
#ifdef DEBUG
   solver->writeLp( "dwmEmpty" );
#endif
   ncIdx.reserve( solver->getNumRows() );
   ncCoef.reserve( solver->getNumRows() );

   currStabFactor = MaxStabFactor;

   bestLB = 0;
   cleanUpLB = 0;
   if (pricing != 0)
   {
      pricingSolver = pricing;
      pricingIt = previousFix = 0;
      doPricing(false, true, _duals);
      solveRelaxation( false );
   }
   else
   {
      if (_duals == 0)
         pricingSolver = new PricingSolverArray( inst );
      else
         pricingSolver = new PricingSolverGraph( inst, _duals , ub );
      //if (_inst->machines() > 1)
      doPricing(false, true);
#ifdef NO_FIXING
      solveRelaxation( false );
#endif
   }
}

DWMaster::DWMaster( DWMaster& other )
      :
      _inst(other._inst),
      idxConsMachines(other.idxConsMachines),
      idxConsSelRouteJob(other.idxConsSelRouteJob),
      idxConsExtCapCuts(other.idxConsExtCapCuts),
      idxAVarsJobs(other.idxAVarsJobs),
      idxAVarsCuts(other.idxAVarsCuts),
      idxNextAVarCut(other.idxNextAVarCut),
      idxACostMachines(other.idxACostMachines),
      idxCQRouteVars(other.idxCQRouteVars),
      nextLbdaIdx(other.nextLbdaIdx),
      nextCutIdx(other.nextCutIdx),
      pricingIt(0),
      previousFix(0),
      cleanUpLB(other.cleanUpLB)
{
   bestLB = other.bestLB;

   _ub = other._ub;

   currStabFactor = other.currStabFactor;

   solver = other.solver->clone();
#ifdef CLP
   clpSolver = dynamic_cast< OsiClpSolverInterface* >( solver );
#endif
#ifdef GLP
   glpSolver = dynamic_cast< OsiGlpSolverInterface* >( solver );
#endif
#ifdef CPX
   cpxSolver = dynamic_cast< OsiCpxSolverInterface* >( solver );
#endif

   lhsQR = new double[other.model->getNumConstraints()];

   model = new Model( *(other.model) );

   aCostJobsAndCuts = new double[ _inst->jobs() + MaxTotalCuts ];
   setACostJobsAndCuts(_ub);

   aCostMachines = (_ub + 1.0) * 100.0;

   bestDuals = other.bestDuals;

   cqRoutes = other.cqRoutes;

   pricingSolver = other.pricingSolver->copy();
}

void DWMaster::createConsSelQRouteJob()
{
   idxConsSelRouteJob = solver->getNumRows();

   const int rows = _inst->jobs();

   vector< int > starts( rows+1, 0 );
   vector< double > lbANDub( rows, 1.0 );

   solver->addRows( rows, &(starts[0]), NULL, NULL, &(lbANDub[0]), &(lbANDub[0]) );
   idxConsExtCapCuts = solver->getNumRows();

   vector< string > names;
   char name[ 256 ];
   for ( int j=1 ; (j<=_inst->jobs()) ; j++ )
   {
      sprintf( name, "cj%05d", j );
      names.push_back( string(name) );
   }

   solver->setRowNames( names, 0, names.size(), idxConsSelRouteJob );
}

void DWMaster::createConsMachines()
{
   idxConsMachines = solver->getNumRows();
   double lbANDub = _inst->machines();
   int starts[2] = { 0, 0 };

   solver->addRows( 1, &(starts[0]) , NULL, NULL, &lbANDub, &lbANDub );
   solver->setRowName( solver->getNumRows()-1, "cMachines");

   // adding artificial var for this row
   double cLB = 0.0, cUB = solver->getInfinity(), cOBJ = aCostMachines;
   idxACostMachines = solver->getNumCols();
   solver->addCol( 1, &idxConsMachines, &lbANDub, cLB, cUB, cOBJ );
   solver->setColName( idxACostMachines, "aMachines" );
}

void DWMaster::setACostJobsAndCuts(double ub)
{
   const int T = _inst->T() + 1;
   for (int j = 1; j <= _inst->jobs(); j++)
      aCostJobsAndCuts[j-1] = (_inst->getCost(j, T) + 1.0) * 100.0;
   for (int c = 0; c < MaxTotalCuts; c++)
      aCostJobsAndCuts[_inst->jobs() + c] = (ub + 1.0) * 100.0;
}

void DWMaster::createAVarsJobsAndCuts()
{
   // set the artificial variabel indices
   idxAVarsJobs = solver->getNumCols();
   idxAVarsCuts = idxAVarsJobs + _inst->jobs();
   idxNextAVarCut = idxAVarsCuts;
   int numAVars = _inst->jobs() + MaxTotalCuts;
   idxCQRouteVars = idxAVarsJobs + numAVars;

   vector< double > lb( numAVars, 0.0 );
   vector< double > ub( numAVars, solver->getInfinity() );

   vector< int > starts( numAVars+1 );
   for (int j = 0 ; j < _inst->jobs(); j++)
      starts[j] = j;
   for (int c = 0; c <= MaxTotalCuts; c++)
      starts[_inst->jobs()+c] = _inst->jobs();

   vector< int > rows( _inst->jobs() );
   for (int j = 0; j < _inst->jobs(); j++)
      rows[j] = idxConsSelRouteJob+j;

   vector< double > elements( _inst->jobs(), 1.0 );

   solver->addCols( numAVars, &(starts[0]), &(rows[0]), &(elements[0]), &(lb[0]), &(ub[0]),
         aCostJobsAndCuts );

   fprintf( stderr, "After add cols: %d rows: %d\n", solver->getNumCols(), solver->getNumRows() );

   vector< string > names;
   char name[256];
   for (int j=0; j < _inst->jobs(); j++)
   {
      sprintf( name, "aj%05d", j );
      names.push_back( name );
   }
   for (int c = 0; c < MaxTotalCuts; c++)
   {
      sprintf( name, "ac%05d", c );
      names.push_back( name );
   }

   solver->setColNames( names, 0, names.size(), idxAVarsJobs );
}

bool DWMaster::isAbandoned() const
{
   return solver->isAbandoned();
}

bool DWMaster::isProvenOptimal() const
{
   return solver->isProvenOptimal();
}

bool DWMaster::isInfeasible() const
{
   return ( (solver->isProvenDualInfeasible()) || (solver->isProvenPrimalInfeasible()) );
}

void DWMaster::solveRelaxation( bool showLog )
{
   // solver->writeLp( "DWM" );
   solver->resolve();

#ifdef CLP
   clpSolver->getModelPtr()->factorize();
   clpSolver->getModelPtr()->primal();
   clpSolver->getModelPtr()->dual();
#elif !defined(CPX)
   #error *** Solver still not supported ***
#endif

   if ( solver->isAbandoned() )
   {
      fprintf( stderr, "Linear Program could not be resolved. Aborting." );
      exit( EXIT_FAILURE );
   }

   if ( (solver->isProvenPrimalInfeasible()) || (solver->isProvenDualInfeasible()) )
   {
      fprintf( stderr, "Linear Program infeasible. Aborting." );
      exit( EXIT_FAILURE );
   }

   if ( !(solver->isProvenOptimal()) )
   {
      fprintf( stderr, "Optimal solution to linear program not reached. Aborting." );
      exit( EXIT_FAILURE );
   }

   if (showLog) fprintf( stderr,", LP: %9.2f\n", solver->getObjValue());

   // check the LP solution for infeasibilities and inconsistencies
   /* const double* x = &(solver->getColSolution()[0]);
   const double* objx = &(solver->getObjCoefficients()[0]);
   double obj = 0.0;
   bool ok = true;
   for (int c = 0; c < idxNextAVarCut; c++)
   {
      if (x[c] < -SolCheckEps)
      {
         //fprintf( stderr, "ERROR: x_%d = %lg\n", c, x[c] );
         ok = false;
         break;
      }
      obj += objx[c] * x[c];
   }
   for (int c = idxCQRouteVars; c < solver->getNumCols(); c++)
   {
      double ub = 1.0;
      if ((int)cqRoutes[c-idxCQRouteVars].size() == 0) ub = 0.0;
      if ((x[c] < -SolCheckEps) || (x[c] > ub + SolCheckEps))
      {
         //fprintf( stderr, "ERROR: x_%d = %lg\n", c, x[c] );
         ok = false;
         break;
      }
      obj += objx[c] * x[c];
   }
#ifdef CPX
   if (!ok)
   {
      CPXENVptr cpxEnv = cpxSolver->getEnvironmentPtr();
      CPXLPptr cpxLp = cpxSolver->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_PROBLEM);
      CPXbaropt( cpxEnv, cpxLp );
      fprintf( stderr, "#" );
   }
#endif */

   // check for lower bound error
   if ((solver->getObjValue() < bestLB - bestLB * 1e-6) &&
         (solver->getObjValue() < bestLB - DwmValueEps) &&
         (bestLB < _ub - 1.0 + IntegerTolerance))
   {
#ifdef CLP
      clpSolver->getModelPtr()->dual();

      fprintf( stderr, "*" );
      if ((solver->getObjValue() < bestLB - bestLB * 1e-6) &&
            (solver->getObjValue() < bestLB - DwmValueEps))
#elif defined(CPX)
      CPXENVptr cpxEnv = cpxSolver->getEnvironmentPtr();
      CPXLPptr cpxLp = cpxSolver->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_PROBLEM);
      CPXbaropt( cpxEnv, cpxLp );
      fprintf( stderr, "*" );

      if ((solver->getObjValue() < bestLB - bestLB * 1e-6) &&
            (solver->getObjValue() < bestLB - DwmValueEps))
#endif
      {
         fprintf( stderr, "Error: LB (%11.4lf) > LP (%11.4lf)\n",
               bestLB, solver->getObjValue() );
         /* const double* pi = solver->getRowPrice();
         fprintf( stderr, "LP dual vars:\n" );
         for (int j = 0; j < solver->getNumRows(); j++)
         {
            fprintf(stderr, "pi_%s=%10.3f\n", j, pi[j]);
         } */
         fprintf(stderr, "Best duals:\n" );
         for (int j = 0; j < (int)bestDuals.size(); j++)
         {
            fprintf(stderr, "pi_%d=%10.3f\n", j, bestDuals[j]);
         }
         pricingSolver->updateConstraints( model );
         pricingSolver->updateReducCosts( &bestDuals[0] );
         double rc, lb;
         pricingSolver->solve( rc, lb );
         fprintf(stderr, "recalculated LB=%10.3f\n", lb);
         /* const double* xsol = solver->getColSolution();
         for (int i = idxCQRouteVars; i < solver->getNumCols(); i++)
         {
            if (fabs(xsol[i]) <= 1E-5) continue;
            CapArcKey* qroute = &(cqRoutes[i-idxCQRouteVars][0]);
            int length = (int)cqRoutes[i-idxCQRouteVars].size();

            // check the lower bound
            lb = solver->getObjCoefficients()[i];
            for (int j = 0; j < (int)bestDuals.size(); j++)
            {
               double coeff = solver->getMatrixByCol()->getCoefficient(j,i);
               lb += model->getConstraint(j)->getRhs() * bestDuals[j];
               lb -= coeff * bestDuals[j];
               double check = 0.0;
               for (int k = 0; k < length; k++)
               {
                  check += model->getConstraint(j)->getVarCoeff(qroute[k].i,qroute[k].j) +
                        model->getConstraint(j)->getCapVarCoeff(qroute[k].i,qroute[k].j,qroute[k].d);
               }
               if (fabs(check - coeff) > 1E-15)
               {
                  fprintf(stderr, "Coeff of col %s in row %s is %g insted of %g\n",
                        solver->getColNames()[i].c_str(),
                        solver->getRowNames()[j].c_str(),
                        coeff, check );
               }
            }
            fprintf(stderr, "LB for column %s=%10.3f\n", solver->getColNames()[i].c_str(), lb);
         } */
         solver->writeLp( "DWM" );
#ifdef CPX
         CPXENVptr cpxEnv = cpxSolver->getEnvironmentPtr();
         CPXLPptr cpxLp = cpxSolver->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_PROBLEM);
         CPXwriteprob(cpxEnv, cpxLp, "DWM_CPX.lp", "lp");
         double val;
         CPXgetobjval(cpxEnv, cpxLp, &val);
         fprintf(stderr, "        CPLEX obj val = %11.4lf\n", val);
#endif
         exit( EXIT_FAILURE );
      }
   }
/*
   char lpname[10];
   sprintf(lpname, "DWM_%d", pricingIt%2);
   solver->writeLp( lpname );
*/

   //if (solver->getObjValue() < bestLB + 0.01)
   //solver->writeLp( "DWM" );
}

bool DWMaster::doPricing( bool showLog, bool update, const double *_duals )
{
   // if new constraints have been added, update them
   bool cleaned = false;
   if (update)
   {
      //printf( "duals:\n" );
      //for (int r = 0; r < bestDuals.size(); r++)
      //   printf( "%10.3f\n", bestDuals[r] );
      pricingSolver->updateConstraints( model );
      bestDuals.resize(model->getNumConstraints(), 0);
   }
   else
   {
      // clean the LP if too many columns have been inserted
      if ((solver->getNumCols() - idxCQRouteVars) > MaxColsRowsRatio * solver->getNumRows())
      {
         if (cleanUpLB < bestLB - DwmValueEps)
         {
            removeUnusedCols();
            cleanUpLB = bestLB;
            cleaned = true;
         }
      }
   }

   // get the array of dual variables from DWM or from the argument
   const double *selDuals = _duals;
   if (selDuals == 0)
      selDuals = solver->getRowPrice()+idxConsSelRouteJob;

   // calculate the stabilized duals for pricing
   std::vector<double> pricingDuals( bestDuals.size() );
   if ((_duals != 0) || (bestLB >= solver->getObjValue() - DwmStabEps))
   {
      for (int r = 0; r < model->getNumConstraints(); r++)
         pricingDuals[r] = selDuals[r];
   }
   else
   {
#ifdef NO_STABILIZATION
      for (int r = 0; r < model->getNumConstraints(); r++)
         pricingDuals[r] = selDuals[r];
#else
      double factor = currStabFactor;
      if (MaxStabFactor * (_ub - bestLB) / (solver->getObjValue() - bestLB) < factor)
         factor = MaxStabFactor * (_ub - bestLB) / (solver->getObjValue() - bestLB);
      for (int r = 0; r < model->getNumConstraints(); r++)
         pricingDuals[r] = bestDuals[r] * (1.0 - factor) + selDuals[r] * factor;
#endif
   }

   // solve the pricing subproblem
   pricingSolver->updateReducCosts( &pricingDuals[0] );
   double rc, lb;
   pricingSolver->solve( rc, lb );

   // calculate the cost and true reduced cost of the generated column
   computeInfoQRoute( solver->getRowPrice() );
   if (showLog)
      fprintf( stderr,"S: % 7d, RC: %9.2f, LB: %9.2f (%9.2f)",
            pricingSolver->getArcCount(), prsRC, lb, bestLB );

/*
   fprintf(stderr, "idxMach=%d, idxJobs=%d\n", idxConsMachines, idxConsSelRouteJob);
   fprintf(stderr, "Duals:\nDWM_000=%10.3f\n", solver->getRowPrice()[0] );
   for (int r = 0; r < model->getNumConstraints(); r++)
      fprintf(stderr, "Duals:\nDWM_%03d=%10.3f, price_%03d=%10.3f\n",
            r+1, solver->getRowPrice()[r+1], r+1, pricingDuals[r] );
*/
   // check if should continue solving the DWM LP
   bool continueLP = ((bestLB < solver->getObjValue() - DwmValueEps) || update);

   // store the array of best duals if necessary
   bool improveLB = (lb > bestLB);
   if (lb > bestLB)
   {
      for (int r = 0; r < model->getNumConstraints(); r++)
         bestDuals[r] = pricingDuals[r];
      bestLB = lb;
      stabChanges++;

      // increase the stabilization factor (look farer from the best dual)
      currStabFactor *= StepStabFactor;
      if (currStabFactor > MaxStabFactor) currStabFactor = MaxStabFactor;

      // fix arcs after 20 iterations if the bound has improved
#ifdef NO_FIXING
      if (!continueLP)
#else
      if ((pricingIt == 0) || (previousFix <= pricingIt - 20) || !continueLP)
#endif
      {
#ifndef NO_FIXING
         pricingSolver->fixArcs( _ub, rc );
#endif
         removeInvalidCols();
         continueLP = ((bestLB < solver->getObjValue() - DwmValueEps) || update);
         previousFix = pricingIt;
      }
   }
   else
   {
      // reduce the stabilization factor (look closer to the best dual)
      currStabFactor /= StepStabFactor;
      if (currStabFactor < MinStabFactor) currStabFactor = MinStabFactor;

      if (!continueLP)
      {
#ifndef NO_FIXING
         pricingSolver->fixArcs( _ub, rc );
#endif
         removeInvalidCols();
         continueLP = ((bestLB < solver->getObjValue() - DwmValueEps) || update);
         previousFix = pricingIt;
      }
   }
   if ((!continueLP) && (pricingSolver->getArcCount() <= MaxArcHeuristic))
   {
      std::vector<CapArcKey> sol( _inst->jobs() + _inst->machines() );
      _ub = pricingSolver->searchFeasible( _ub, &sol[0] );
   }
   pricingIt++;

   // insert the column in the DWM if it has a negative reduced cost (using
   // the duals from the DWM)
   if ( (prsRC < -RedCostEps) && continueLP  &&
         (lb < _ub - 1.0 + IntegerTolerance) )
      insertQRoute();
   else if ((bestLB < solver->getObjValue() - DwmValueEps) &&
         (lb < _ub - 1.0 + IntegerTolerance) &&
         (pricingIt > 2) && !cleaned)
   {
      missPricings++;
      //fprintf( stderr, "Miss pricing: it = %d, fix = %d, lb = %g, gap = %g\n",
      //      pricingIt, previousFix, bestLB, solver->getObjValue() - bestLB );
   }

   if ( (prsRC >= -RedCostEps)  && (pricingIt > previousFix + 1) &&
         (lb < _ub - 1.0 + IntegerTolerance) && continueLP && !update &&
         !improveLB )
   {
      fprintf( stderr, "ERROR: rc=%lf, LB=%lf(%lf) (not improved), DWM=%lf.\n",
            prsRC, lb, bestLB, solver->getObjValue() );
      const double* duals = solver->getRowPrice();
      fprintf( stderr, "it = %d, fix = %d, duals:\n", pricingIt, previousFix );
      for (int r = 0; r < model->getNumConstraints(); r++)
      {
         double dual = duals[ idxConsSelRouteJob+r ];
         double pdual = pricingDuals[r];
         double bdual = bestDuals[r];
         fprintf( stderr, "pi_%d=%lf<-(%lf)-%lf\n", r, dual, pdual, bdual );
      }
      pricingSolver->updateReducCosts( &bestDuals[0] );
      pricingSolver->solve( rc, lb );
      fprintf( stderr, "recalculated best LB = %lf\n", lb );
      pricingSolver->updateReducCosts( &pricingDuals[0] );
      pricingSolver->solve( rc, lb );
      fprintf( stderr, "recalculated pricing LB = %lf\n", lb );
      solver->writeLp( "DWM" );
      throw(-908934);
   }

   // stop only if the LB reaches the DWM value
   return (continueLP && (lb < _ub - 1.0 + IntegerTolerance));
}

void DWMaster::computeInfoQRoute( const double *_duals )
{
   std::fill( lhsQR, lhsQR+model->getNumConstraints(), 0 );

   int length;
   CapArcKey *qRoute = pricingSolver->getSolution( length );

#ifdef DEBUG
   assert( (length>0) );
   printf("best route from pricing:\n");
   for ( int i=0 ; (i<length) ; i++ )
      if ( ( qRoute[i].j > 0 ) && ( qRoute[i].j <= this->_inst->jobs() ) )
         printf("%d ", qRoute[i].j );
   printf("\n");
#endif

   // calculate the original cost and make a copy of the current q-route
   cqRoute.resize(length);
   prsCost = 0.0;
   for (int k = 0; k < length; k++)
   {
      cqRoute[k] = qRoute[k];
      const double jobCost = _inst->getCost( qRoute[k].j, qRoute[k].d );
      prsCost += jobCost;
   }

   // calculate the reduced cost of the current q-route
   prsRC   = prsCost;
   for (int r = 0; r < model->getNumConstraints(); r++)
   {
      double dual = _duals[ idxConsSelRouteJob+r ];
      Constraint* constr = model->getConstraint(r);
      for (int k = 0; k < length; k++)
      {
         // consider the coefficient of the arc variable
         double coeff = constr->getVarCoeff(qRoute[k].i, qRoute[k].j);
         if (fabs(coeff) > DwmCoeffEps)
         {
            lhsQR[r] += coeff;
            prsRC -= coeff * dual;
         }

         // consider the coefficient of the extended arc variable
         coeff = constr->getCapVarCoeff(qRoute[k].i, qRoute[k].j, qRoute[k].d);
         if (fabs(coeff) > DwmCoeffEps)
         {
            lhsQR[r] += coeff;
            prsRC -= coeff * dual;
         }
      }
   }

   if (_inst->machines() > 1)
      prsRC -= (_duals[ idxConsMachines ]);
   delete [] qRoute;
}

void DWMaster::insertRowsInLP()
{
   int firstRow = solver->getNumRows();
   int numRows = idxConsSelRouteJob + model->getNumConstraints() - firstRow;

   assert( numRows <= (idxCQRouteVars - idxNextAVarCut) );

   // resize the LHS array
   delete[] lhsQR;
   lhsQR = new double[model->getNumConstraints()];

   // initialize the arrays used to insert the rows
   vector< int > starts( numRows+1, 0 );
   vector< double > lbs( numRows );
   vector< double > ubs( numRows );
   vector< int > columns;
   vector< double > coeffs;

   // calculate the row coeficients
   int pos = 0;
   for (int r = 0; r < numRows; r++)
   {
      // set the RHS and the coefficient for the slack variable
      Constraint* constr = model->getConstraint(firstRow+r-idxConsSelRouteJob);
      if (constr->getType() == '>')
      {
         lbs[r] = constr->getRhs();
         ubs[r] = solver->getInfinity();
      }
      else
      {
         lbs[r] = -solver->getInfinity();
         ubs[r] = constr->getRhs();
      }
      starts[r] = pos;
      if (constr->getRhs() != 0.0)
      {
         pos++;
         columns.push_back(idxNextAVarCut++);
         if (constr->getType() == '>')
            coeffs.push_back(1.0);
         else
            coeffs.push_back(-1.0);
      }

      // set the coefficients for the qroutes
      for (int c = idxCQRouteVars; c < solver->getNumCols(); c++)
      {
         // skip the fixed qroutes
         int length = cqRoutes[c-idxCQRouteVars].size();
         if (length == 0) continue;
         CapArcKey* qRoute = &(cqRoutes[c-idxCQRouteVars][0]);

         // calculate the coefficient
         double coeff = 0.0;
         for (int k = 0; k < length; k++)
         {
            coeff += constr->getVarCoeff(qRoute[k].i, qRoute[k].j);
            coeff += constr->getCapVarCoeff(qRoute[k].i, qRoute[k].j, qRoute[k].d);
         }

         // add to the DWM LP if necessary
         if (fabs(coeff) > DwmCoeffEps)
         {
            columns.push_back(c);
            coeffs.push_back(coeff);
            pos++;
         }
      }
   }
   starts[numRows] = pos;

   // insert the rows
   solver->addRows( numRows, &(starts[0]), &(columns[0]), &(coeffs[0]), &(lbs[0]), &(ubs[0]) );
   fprintf( stderr, "Added %d rows, having %d rows after it\n", numRows, solver->getNumRows() );

   vector< string > names;
   char name[ 256 ];
   for (int j = 0; j < numRows; j++)
   {
      sprintf( name, "cut%05d", nextCutIdx++ );
      names.push_back( string(name) );
   }

   solver->setRowNames( names, 0, names.size(), firstRow );
   //solver->writeLp( "DWM_w_cuts" );
}

int DWMaster::insertQRoute()
{
   ncIdx.clear(); ncCoef.clear();
   for (int j = 0; j < model->getNumConstraints(); j++)
   {
      if (fabs(lhsQR[j]) > DwmCoeffEps)
      {
         ncIdx.push_back( idxConsSelRouteJob+j );
         ncCoef.push_back( lhsQR[j] );
      }
   }

   if (_inst->machines() > 1)
   {
      ncIdx.push_back( idxConsMachines );
      ncCoef.push_back( 1.0 );
   }

   int cIdx = solver->getNumCols();
   char name[80];
   sprintf( name, "qroute_%d", nextLbdaIdx++ );
   //sprintf( name, "%s_%d", COL_LBDA_PREFIX.c_str(), nextLbdaIdx++ );
#ifdef DEBUG
   printf("inserted column %s\n", name); fflush(stdout);
#endif
   solver->addCol( ncIdx.size(), &(ncIdx[0]), &(ncCoef[0]), 0.0, 1.0, prsCost, string(name) );
   //solver->setInteger( cIdx );

   cqRoutes.push_back(cqRoute);

   return cIdx;
}

void DWMaster::removeInvalidCols()
{
   int count = 0;
   for (int c = idxCQRouteVars; c < solver->getNumCols(); c++)
   {
      // skip the fixed qroutes
      int length = cqRoutes[c-idxCQRouteVars].size();
      if (length == 0) continue;
      CapArcKey* qRoute = &(cqRoutes[c-idxCQRouteVars][0]);

      // if the current route is invalid...
      std::vector<int> sol;
      for (int k = 0; k < length-1; k++)
         sol.push_back( qRoute[k].j );
      if ( !pricingSolver->checkSolution(&sol[0], sol.size(), false) )
      {
         // remove from the LP
         solver->setColUpper(c, 0.0);
         cqRoutes[c-idxCQRouteVars].clear();
         count++;
      }
   }
   solveRelaxation( false );
   //fprintf( stderr, "%d columns deleted\n", count );
}

void DWMaster::removeUnusedCols()
{
   int count = 0;
   for (int c = idxCQRouteVars; c < solver->getNumCols(); c++)
   {
      // skip the fixed qroutes
      int length = cqRoutes[c-idxCQRouteVars].size();
      if (length == 0) continue;

      // if the column is unused...
      double rc = solver->getReducedCost()[c];
      //double x = solver->getColSolution()[c];
      if ( rc > RedCostEps )
      //if (x < ColCleanXEps)
      {
         // remove from the LP
         cqRoutes[c-idxCQRouteVars].clear();
         count++;
      }
   }
   fprintf( stderr, "%d columns deleted\n", count );

   // really clean the LP
   int cc = idxCQRouteVars;
   int numCols = solver->getNumCols();
   for (int c = idxCQRouteVars; c < numCols; c++)
   {
      // remove and skip the fixed qroutes
      int length = cqRoutes[c-idxCQRouteVars].size();
      if (length == 0)
      {
         // really remove from the LP
         solver->deleteCols(1, &cc);
         continue;
      }

      // if the column is shifted...
      if ( c > cc )
         cqRoutes[cc-idxCQRouteVars] = cqRoutes[c-idxCQRouteVars];
      cc++;
   }
   cqRoutes.resize( solver->getNumCols() - idxCQRouteVars );
   fprintf( stderr, "Clean LP has %d columns\n", solver->getNumCols() );
}

void DWMaster::changePricing()
{
   fprintf( stderr,"\nChanging the pricing solver to arc-time indexed...\n" );
   delete pricingSolver;
   pricingSolver = new PricingSolverGraph( _inst, &bestDuals[0] , _ub );

   // wt40-4m-36: cost 6420
   //int r1[] = { 14, 7, 19, 22, 2, 40, 15, 1, 29 };
   //int r2[] = { 37, 28, 24, 13, 3, 4, 18, 20, 34, 27, 23, 5 };
   //int r3[] = { 38, 26, 9, 6, 33, 39, 8, 25, 11, 35, 31 };
   //int r4[] = { 12, 10, 21, 36, 17, 32, 30, 16 };

   // wt40-2m-1: cost 606
   //int r1[] = { 39, 8, 18, 32, 21, 4, 9, 15, 16, 31, 10, 33, 27, 1, 17, 7, 35, 34, 25, 22, 26, 6, 38 };
   //int r2[] = { 13, 40, 29, 24, 3, 5, 28, 14, 30, 2, 11, 20, 12, 23, 36, 19, 37 };

   /* int cost = 0;
   int t = 0;
   for (int k = int(sizeof(r1)/sizeof(int)) - 1; k>= 0; k--)
   {
      t += _inst->ptime()[r1[k]];
      cost += _inst->getCost(r1[k], t);
   }
   t = 0;
   for (int k = int(sizeof(r2)/sizeof(int)) - 1; k >= 0; k--)
   {
      t += _inst->ptime()[r2[k]];
      cost += _inst->getCost(r2[k], t);
   }
   t = 0;
   for (int k = int(sizeof(r3)/sizeof(int)) - 1; k >= 0; k--)
   {
      t += _inst->ptime()[r3[k]];
      cost += _inst->getCost(r3[k], t);
   }
   t = 0;
   for (int k = int(sizeof(r4)/sizeof(int)) - 1; k >= 0; k--)
   {
      t += _inst->ptime()[r4[k]];
      cost += _inst->getCost(r4[k], t);
   }

   fprintf( stderr, "Checking solution with cost %d...\n", cost );
   pricingSolver->checkSolution( r1, sizeof(r1)/sizeof(int), true );
   pricingSolver->checkSolution( r2, sizeof(r2)/sizeof(int), true );
   pricingSolver->checkSolution( r3, sizeof(r3)/sizeof(int), true );
   pricingSolver->checkSolution( r4, sizeof(r4)/sizeof(int), true ); */

   removeInvalidCols();
   fprintf( stderr,"\n" );
   doPricing(false, true, &bestDuals[0]);
}

void DWMaster::removeUnusedRows()
{
   int numRows = solver->getNumRows();
   const double* osilhs = solver->getRowActivity();
   vector<double> lhs(numRows);
   for (int r = 0; r < numRows; r++) lhs[r] = osilhs[r];
   int rr = idxConsExtCapCuts;
   for (int r = idxConsExtCapCuts; r < numRows; r++)
   {
      // if the dual variable is zero, remove and skip the row
      //if ( (fabs(bestDuals[r-idxConsSelRouteJob]) < CutCleanDualEps) &&
      //      model->getConstraint(r-idxConsSelRouteJob)->canBeDeleted() )

      // if the constraint is not tight, remove and skip the row
      double slack = lhs[r] -
            model->getConstraint(r-idxConsSelRouteJob)->getRhs();
      if (model->getConstraint(r-idxConsSelRouteJob)->getType() == '<')
         slack = -slack;
      if ( (slack > CutCleanSlackEps) &&
            (fabs(bestDuals[r-idxConsSelRouteJob]) < CutCleanDualEps) &&
            model->getConstraint(r-idxConsSelRouteJob)->canBeDeleted() )
      {
         solver->deleteRows(1, &rr);
         delete model->getConstraint(r-idxConsSelRouteJob);
         continue;
      }

      // move the row information to the correct place if necessary
      if (rr < r)
      {
         bestDuals[rr-idxConsSelRouteJob] = bestDuals[r-idxConsSelRouteJob];
         model->setConstraint(rr-idxConsSelRouteJob,
               model->getConstraint(r-idxConsSelRouteJob));
      }
      rr++;
   }
   model->resizeConstraintArray(rr-idxConsSelRouteJob);
   bestDuals.resize(rr-idxConsSelRouteJob);
   fprintf( stderr, "Clean LP has %d rows\n", solver->getNumRows() );
}

const double *DWMaster::getRCs() const
{
   return solver->getReducedCost();
}

const double *DWMaster::getRowPrice() const
{
   return solver->getRowPrice();
}

const double *DWMaster::getColSolution() const
{
   return &(solver->getColSolution()[0]);
}

CapArcKey *DWMaster::getCQRoute( const int idx, int& length )
{
   assert( idx < (int)cqRoutes.size()+idxCQRouteVars );

   if (idx < idxCQRouteVars)
   {
      length = 0;
      return 0;
   }
   else
   {
      length = cqRoutes[ idx-idxCQRouteVars ].size();
      return &(cqRoutes[ idx-idxCQRouteVars ][0]);
   }
}

DWMaster::~DWMaster()
{
   delete[] lhsQR;
   delete pricingSolver;
   delete model;

#ifdef CLP
   delete clpSolver;
#endif
#ifdef CPX
   delete cpxSolver;
#endif
#ifdef GLP
   delete glpSolver;
#endif

   delete[] aCostJobsAndCuts;
}
