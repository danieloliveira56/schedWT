/** Source Code File:
 ** Object the solves the pricing subproblem
 **/

#include "../Instance.hpp"
#include "PricingSolver.hpp"
#include "Model.hpp"
#include "ConstrShadow.hpp"
#include "OrigVariable.hpp"

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <climits>
#include <set>
#include <list>

#define FIX_ROOT_ARCS

// number of tries when searching for good feasible solutions
const int numTriesFeasible = 1000;

// tolerance for coefficients of the generated LP
const double LpCoeffEps = 1E-12;

/*****************************************************************************
 *                             PricingSolverGraph                            *
 *****************************************************************************/

PricingSolverGraph::PricingSolverGraph( Instance *inst, const double* duals,
      double ub ) :
      _inst( inst )
{
   // all jobs, including dummy
   const int jobs    = inst->jobs() + 1;
   const int lastJob = jobs-1;
   // all times, including zero
   const int T       = inst->T() + 1;
   const int lastT   = T-1;
   const int idx     = -1;

   // if there is a vector of dual variables
   double** fwd = 0;
   double** bwd = 0;
   double bestRc = 0.0;
   if (duals != 0)
   {
      // fill arrays of lower bounds by time forward and backward
      fwd = new double*[T];
      bwd = new double*[T];
      for (int t = 0; t < T; t++)
      {
         fwd[t] = new double[jobs];
         bwd[t] = new double[jobs];
      }
      updateLowerBounds(fwd, bwd, duals);
      constLB = 0;
      for (int j = 0; j < lastJob; j++)
         constLB += duals[j];
      if (_inst->machines() > 1)
         for (int t = 1; t < T; t++)
            if (bestRc > fwd[t][0])
               bestRc = fwd[t][0];
      //for (int j = 1; j <= lastJob; j++)
      //   printf( "fwd_%d = %g, bwd_%d = %g\n",
      //         j, fwd[lastT][j] + constLB, j, bwd[0][j] + constLB );
   }

   /*** create the graph of original variables ***/
   // create each node
   ArcNode** nodes = new ArcNode*[ jobs ];
   for ( int i = 0; ( i < jobs ); i++ )
      nodes[i] = new ArcNode( i, lastJob, lastJob );
   root = nodes[0];

   // create each capacitated node
   CapArcNode*** capNodes = new CapArcNode**[ T ];
   for ( int t = 0; (t<T) ; t++ )
      capNodes[t] = new CapArcNode*[ jobs ];
   capNodes[ 0 ][ 0 ] = new CapArcNode( 0, 0, lastJob, 1 );
   capNodes[ lastT ][ 0 ] = new CapArcNode( lastT, 0, 1, T*lastJob );
   for ( int t = 1; ( t < T ) ; t++ )
   {
      for ( int i = 1 ; ( i< jobs ) ; i++ )
      {
         int maxIn, maxOut;
         maxIn = maxOut = jobs;

         // special cases
         if ( t <= _inst->ptime()[i] )
            maxOut = 1;
         else if (t == lastT)
            maxIn = 1;
         else if (duals != 0)
         {
            maxIn = 0;
            for ( int j = 0 ; ( j< jobs ) ; j++ )
            {
               if (i == j) continue;
               if ((j == 0) && (_inst->machines() == 1)) continue;

               int currJobJ = i;
               int currJobI = j;
               int endp = t + _inst->ptime()[currJobI];
               if (endp > lastT) continue;
               int startp = t - _inst->ptime()[currJobJ];
               if (startp < 0) continue;
               if (duals != 0)
               {
                  // calculate the lower bound on the cost including this arc
                  double LB = fwd[t][currJobJ] + _inst->getCost(currJobJ, t) - 
                        duals[idx + currJobJ] + bwd[t][currJobI] +
                        (_inst->machines() - 1) * bestRc + constLB;

                  // skip if it is less than one unit bellow the upper bound
                  if (LB > ub - 1.0 + IntegerTolerance) continue;
               }

               if ((j != 0) && !checkArc(currJobI, currJobJ, t) ) continue;
               maxIn++;
            }
            maxOut = 0;
            for ( int j = 1 ; ( j< jobs ) ; j++ )
            {
               if (i == j) continue;

               int currJobJ = j;
               int currJobI = i;
               int tt = t - _inst->ptime()[currJobI];
               int startp = tt - _inst->ptime()[currJobJ];
               if (startp < 0) continue;
               if (duals != 0)
               {
                  // calculate the lower bound on the cost including this arc
                  double LB = fwd[tt][currJobJ] + _inst->getCost(currJobJ, tt) - 
                        duals[idx + currJobJ] + bwd[tt][currJobI] +
                        (_inst->machines() - 1) * bestRc + constLB;

                  // skip if it is less than one unit bellow the upper bound
                  if (LB > ub - 1.0 + IntegerTolerance) continue;
               }

               if ( !checkArc(currJobI, currJobJ, tt) ) continue;
               maxOut++;
            }
         }

         if ((maxIn > 0) && (maxOut > 0))
            capNodes[ t ][ i ] = new CapArcNode( t, i, maxIn, maxOut );
         else
            capNodes[ t ][ i ] = 0;
      } // all jobs
   } // all capacities

   capStart = capNodes[ lastT ][ 0 ];
   capEnd   = capNodes[ 0 ][ 0 ];

   // create each modulus node
   ModCapArcNode*** modNodes = new ModCapArcNode**[ T ];
   for ( int t = 0; (t < T) ; t++ )
      modNodes[t] = new ModCapArcNode*[ jobs ];
   modNodes[0][0] = new ModCapArcNode( 0, 0, 0, lastJob, 1 );
   modNodes[lastT][0] = new ModCapArcNode( lastT, 0, 0, 1, T*lastJob );
   for ( int t = 1; ( t < T ) ; t++ )
      for ( int i = 1; ( i < jobs ) ; i++ )
      {
         int maxIn, maxOut;
         maxIn = maxOut = jobs;

         // special cases
         if ( t <= _inst->ptime()[i] )
            maxOut = 1;
         else if (t == lastT)
            maxIn = 1;
         else if (duals != 0)
         {
            maxIn = 0;
            for ( int j = 0 ; ( j< jobs ) ; j++ )
            {
               if (i == j) continue;
               if ((j == 0) && (_inst->machines() == 1)) continue;

               int currJobJ = i;
               int currJobI = j;
               int endp = t + _inst->ptime()[currJobI];
               if (endp > lastT) continue;
               int startp = t - _inst->ptime()[currJobJ];
               if (startp < 0) continue;
               if (duals != 0)
               {
                  // calculate the lower bound on the cost including this arc
                  double LB = fwd[t][currJobJ] + _inst->getCost(currJobJ, t) - 
                        duals[idx + currJobJ] + bwd[t][currJobI] +
                        (_inst->machines() - 1) * bestRc + constLB;

                  // skip if it is less than one unit bellow the upper bound
                  if (LB > ub - 1.0 + IntegerTolerance) continue;
               }

               if ((j != 0) && !checkArc(currJobI, currJobJ, t) ) continue;
               maxIn++;
            }
            maxOut = 0;
            for ( int j = 1 ; ( j< jobs ) ; j++ )
            {
               if (i == j) continue;

               int currJobJ = j;
               int currJobI = i;
               int tt = t - _inst->ptime()[currJobI];
               int startp = tt - _inst->ptime()[currJobJ];
               if (startp < 0) continue;
               if (duals != 0)
               {
                  // calculate the lower bound on the cost including this arc
                  double LB = fwd[tt][currJobJ] + _inst->getCost(currJobJ, tt) - 
                        duals[idx + currJobJ] + bwd[tt][currJobI] +
                        (_inst->machines() - 1) * bestRc + constLB;

                  // skip if it is less than one unit bellow the upper bound
                  if (LB > ub - 1.0 + IntegerTolerance) continue;
               }

               if ( !checkArc(currJobI, currJobJ, tt) ) continue;
               maxOut++;
            }
         }

         if ((maxIn > 0) && (maxOut > 0))
            modNodes[t][i] = new ModCapArcNode( t, 0, i, maxIn, maxOut );
         else
            modNodes[t][i] = 0;
      }
   modCapStart = modNodes[lastT][0];
   modCapEnd = modNodes[0][0];

   // create each arc, capacitated arc and modulus arc
   arcCount = 0;
   ArcVariable* arc;
   CapArcVariable* capArc;
   for ( int i = 0; ( i < jobs ) ; i++ )
      for ( int j = 0; ( j < jobs ) ; j++ )
      {
         // create the arc
         if ( i == j ) continue;

         arc = new ArcVariable( nodes[j], nodes[i], T );

         // create the capacitated arcs and modulus arcs
         if ( i == 0 )
         {
            for ( int t=inst->Tmin() ; ( t<T ) ; t++ )
            {
               //printf(" i: %d j: %d t: %d cost: %d\n", i, j, t, inst->getCost( j, t ) );
               //fflush(stdout);
               if (!capNodes[ t ][ j ])
                  continue;
               capArc = new CapArcVariable( t, capNodes[ t ][ j ],
                                            capStart, arc, inst->getCost( j, t ), 1 );

               new ModCapArcVariable( t, 0, modNodes[ t ][ j ],
                                      modCapStart, capArc );
               arcCount++;
            }
         }
         else if ( j == 0 )
         {
            capArc = new CapArcVariable( 0, capNodes[0][j], capNodes[ inst->ptime()[i] ][i],
                                         arc, 0, 1 );
            new ModCapArcVariable( 0, 0, modNodes[0][j], modNodes[ inst->ptime()[i] ][i],
                                   capArc );
            arcCount++;
         }
         else
         {
            const int currJobJ = j;
            const int currJobI = i;

            for ( int t = 1; ( t < lastT ) ; t++ )
            {
               const int endp = t + _inst->ptime()[currJobI];
               if (endp > lastT) continue;
               const int startp = t - _inst->ptime()[currJobJ];
               if (startp < 0) continue;

               if (duals != 0)
               {
                  // calculate the lower bound on the cost including this arc
                  double LB = fwd[t][j] + _inst->getCost(j, t) - 
                        duals[idx + j] + bwd[t][i] +
                        (_inst->machines() - 1) * bestRc + constLB;

                  // skip if it is less than one unit bellow the upper bound
                  if (LB > ub - 1.0 + IntegerTolerance) continue;
               }

               if ( !checkArc(currJobI, currJobJ, t) ) continue;

               //printf(" i: %d j: %d t: %d endp: %d cost: %d\n", i, j, t, endp, costJinHead );

               capArc = new CapArcVariable( t, capNodes[t][j],
                     capNodes[endp][i], arc, _inst->getCost(currJobJ, t), 1 );
               new ModCapArcVariable( t, 0, modNodes[t][j],
                     modNodes[endp][i], capArc );
               arcCount++;
            }
         }

         // remove the arcs with no capacitated arc associated
         if (arc->getNumCaps() == 0)
            delete arc;
      }

   std::set< ModCapArcNode* > toBeDeleted;
   // checking for non-reachable nodes
   for ( int t = 1; ( t < T ) ; t++ )
      for ( int i = 1; ( i < jobs ) ; i++ )
      {
         if (modNodes[t][i] == 0) continue;
         ModCapArcNode *mcArcNode = modNodes[t][i];
         if ( mcArcNode->getOutDegree() == 0 )
         {
            if ( mcArcNode->getInDegree() == 0 )
            {
               delete mcArcNode;
               delete capNodes[t][i];
            }
            else
               toBeDeleted.insert( mcArcNode );
         }
      }

   for ( std::set< ModCapArcNode* >::iterator sIt=toBeDeleted.begin() ; (sIt!=toBeDeleted.end()) ; sIt++ )
      delete (*sIt);

   // release the auxiliary data
   delete [] nodes;
   for ( int d = 0; ( d < T ) ; d++ )
      delete [] capNodes[d];
   delete [] capNodes;
   for ( int d = 0; ( d < T ) ; d++ )
      delete [] modNodes[d];
   delete [] modNodes;
   if (fwd != 0)
   {
      for (int t = 0; t < T; t++)
         delete [] fwd[t];
      delete [] fwd;
   }
   if (bwd != 0)
   {
      for (int t = 0; t < T; t++)
         delete [] bwd[t];
      delete [] bwd;
   }
   fprintf( stderr, "Created %d time-indexed arcs.\n", arcCount );

   /*** initialize other data ***/
   numConstr = 0;
   constraints = 0;
   constLB = 0.0;
   searchCount = 0;
   srand( 123456 );
}

PricingSolverGraph::PricingSolverGraph( Instance *inst, ArcNode* root_,
      CapArcNode* capStart_, CapArcNode* capEnd_,
      ModCapArcNode* modCapStart_, ModCapArcNode* modCapEnd_ ) :
      _inst( inst ),
      root( root_ ),
      capStart( capStart_ ),
      capEnd( capEnd_ ),
      modCapStart( modCapStart_ ),
      modCapEnd( modCapEnd_ )
{
   // all jobs, including dummy
   const int jobs    = inst->jobs() + 1;
   const int lastJob = jobs-1;
   // all times, including zero
   const int T       = inst->T() + 1;
   const int lastT   = T-1;

   /*** find the references to existing nodes ***/
   // create each node
   ArcNode** nodes_ = new ArcNode*[ jobs ];
   for ( int i = 0; ( i < jobs ); i++ )
      nodes_[i] = root_;

   // create each capacitated node
   CapArcNode*** capNodes_ = new CapArcNode**[ T ];
   for ( int t = 0; (t<T) ; t++ )
      capNodes_[t] = new CapArcNode*[ jobs ];
   capNodes_[ 0 ][ 0 ] = capEnd_;
   capNodes_[ lastT ][ 0 ] = capStart_;
   for ( int t = 1; ( t < T ) ; t++ )
      for ( int i = 1 ; ( i< jobs ) ; i++ )
         capNodes_[ t ][ i ] = 0;

   // create each modulus node
   ModCapArcNode*** modNodes_ = new ModCapArcNode**[ T ];
   for ( int t = 0; (t < T) ; t++ )
      modNodes_[t] = new ModCapArcNode*[ jobs ];
   modNodes_[0][0] = modCapEnd_;
   modNodes_[lastT][0] = modCapStart_;
   for ( int t = 1; ( t < T ) ; t++ )
      for ( int i = 1; ( i < jobs ) ; i++ )
         modNodes_[t][i] = 0;

   // traverse the graph saving the references to existing nodes
   std::list<ModCapArcNode*> nodeQueue;
   nodeQueue.push_back( modCapStart_ );
   while ( !nodeQueue.empty() )
   {
      ModCapArcNode* currMod = nodeQueue.front();
      nodeQueue.pop_front();

      // take the adjacent arcs
      for (int k = 0; k < currMod->getOutDegree(); k++)
      {
         ModCapArcVariable* modArc = currMod->getOut(k);
         ModCapArcNode* nextMod = modArc->getHead();
         int t = nextMod->getCap();
         int i = nextMod->getIndex();
         if (modNodes_[t][i] == 0)
         {
            // save the references to the adjacent node
            modNodes_[t][i] = nextMod;
            capNodes_[t][i] = modArc->getCapArc()->getHead();
            nodes_[i] = modArc->getCapArc()->getArc()->getHead();

            // store the adjacent node for future inspection
            nodeQueue.push_back( nextMod );
         }
      }
   }

   /*** create the copy of the graph of original variables ***/
   // create each node
   ArcNode** nodes = new ArcNode*[ jobs ];
   for ( int i = 0; ( i < jobs ); i++ )
   {
      if (nodes_[i] != 0)
         nodes[i] = new ArcNode( i, lastJob, lastJob );
      else
         nodes[i] = 0;
   }
   root = nodes[0];

   // create each capacitated node
   CapArcNode*** capNodes = new CapArcNode**[ T ];
   for ( int t = 0; (t<T) ; t++ )
      capNodes[t] = new CapArcNode*[ jobs ];
   capNodes[ 0 ][ 0 ] = new CapArcNode( 0, 0, lastJob, 1 );
   capNodes[ lastT ][ 0 ] = new CapArcNode( lastT, 0, 1, T*lastJob );
   for ( int t = 1; ( t < T ) ; t++ )
   {
      for ( int i = 1 ; ( i< jobs ) ; i++ )
      {
         CapArcNode* capNode = capNodes_[ t ][ i ];
         if (capNode != 0)
         {
            capNodes[ t ][ i ] = new CapArcNode( t, i,
                  capNode->getInDegree(), capNode->getOutDegree() );
         }
         else
            capNodes[ t ][ i ] = 0;
      } // all jobs
   } // all capacities

   capStart = capNodes[ lastT ][ 0 ];
   capEnd   = capNodes[ 0 ][ 0 ];

   // create each modulus node
   ModCapArcNode*** modNodes = new ModCapArcNode**[ T ];
   for ( int t = 0; (t < T) ; t++ )
      modNodes[t] = new ModCapArcNode*[ jobs ];
   modNodes[0][0] = new ModCapArcNode( 0, 0, 0, lastJob, 1 );
   modNodes[lastT][0] = new ModCapArcNode( lastT, 0, 0, 1, T*lastJob );
   for ( int t = 1; ( t < T ) ; t++ )
      for ( int i = 1; ( i < jobs ) ; i++ )
      {
         ModCapArcNode* modNode = modNodes_[t][i];
         if (modNode != 0)
            modNodes[t][i] = new ModCapArcNode( t, 0, i,
                  modNode->getInDegree(), modNode->getOutDegree() );
         else
            modNodes[t][i] = 0;
      }
   modCapStart = modNodes[lastT][0];
   modCapEnd = modNodes[0][0];

   // create each arc, capacitated arc and modulus arc
   arcCount = 0;
   ArcVariable* arc;
   ArcVariable* arc_;
   CapArcVariable* capArc_;
   CapArcVariable* capArc;
   for ( int i = 0; ( i < jobs ) ; i++ )
   {
      for (int k = 0; k < nodes_[i]->getOutDegree(); k++)
      {
         // create the arc
         arc_ = nodes_[i]->getOut( k );
         int j = arc_->getHead()->getIndex();
         arc = new ArcVariable( nodes[j], nodes[i], arc_->getNumCaps() );

         // create the capacitated and modulus arcs
         for (int c = 0; c < arc_->getNumCaps(); c++)
         {
            capArc_ = arc_->getCap( c );
            int t = capArc_->getCap();
            int ti = capArc_->getTail()->getCap();
            int tj = capArc_->getHead()->getCap();
            capArc_->resetReducCost();
            capArc = new CapArcVariable( t, capNodes[tj][j], capNodes[ti][i],
                  arc, capArc_->getReducCost(), 1 );
            new ModCapArcVariable( t, 0, modNodes[tj][j], modNodes[ti][i],
                  capArc );
            arcCount++;
         }
      }
   }

   // release the auxiliary data
   delete [] nodes;
   delete [] nodes_;
   for ( int d = 0; ( d < T ) ; d++ )
   {
      delete [] capNodes[d];
      delete [] capNodes_[d];
   }
   delete [] capNodes;
   delete [] capNodes_;
   for ( int d = 0; ( d < T ) ; d++ )
   {
      delete [] modNodes[d];
      delete [] modNodes_[d];
   }
   delete [] modNodes;
   delete [] modNodes_;
   fprintf( stderr, "Created %d time-indexed arcs.\n", arcCount );

   /*** initialize other data ***/
   numConstr = 0;
   constraints = 0;
   constLB = 0.0;
   searchCount = 0;
}

PricingSolverGraph::~PricingSolverGraph()
{
   // release the original variable graph
   delete modCapStart; // almost all graph is released with this deletion
   delete modCapEnd;
   delete capStart;
   delete capEnd;
   delete root;

   // release the array of constraint shadows and other arrays
   if ( constraints != 0 ) delete [] constraints;
}

PricingSolver* PricingSolverGraph::copy()
{
   return new PricingSolverGraph( _inst, root, capStart, capEnd,
         modCapStart, modCapEnd );
}

void PricingSolverGraph::updateConstraints( Model* model )
{
   // release the array of constraint shadows
   if ( constraints != 0 ) delete [] constraints;

   // create a new one and associate with the true constraints
   numConstr = model->getNumConstraints();
   constraints = new ConstrShadow [numConstr];
   for ( int c = 0; c < model->getNumConstraints(); c++ )
      constraints[c].setConstr( model->getConstraint( c ) );

   // search the graph at the arc level
   searchCount++;
   ArcNode* queue = root;
   root->setAuxNext( root );
   root->setAuxCount( searchCount );
   for ( ;; )
   {
      // get the first node in the queue
      ArcNode* first = queue->getAuxNext();
      first->setAuxCount( searchCount );

      // check all adjacent arcs
      for ( int a = 0; a < first->getOutDegree(); a++ )
      {
         ArcVariable* arc = first->getOut( a );
         ArcNode* head = arc->getHead();
         int i = first->getIndex();
         int j = head->getIndex();

         // check all constraints
         for ( int c = 0; c < numConstr; c++ )
         {
            // check the arc without capacity
            Constraint* constr = model->getConstraint( c );
            double coeff = constr->getVarCoeff( i, j );
            CoeffShadow* aux;
            if ( fabs( coeff ) > CoeffTolerance )
            {
               aux = new CoeffShadow( 0, arc->getShadow(), coeff );
               constraints[c].addCoeff( aux );
            }

            // check if a capacitated arc exists
            if ( constr->checkCapVarCoeff( i, j ) )
            {
               for ( int ca = 0; ca < arc->getNumCaps(); ca++ )
               {
                  CapArcVariable* capArc = arc->getCap( ca );
                  int d = capArc->getCap();

                  // check the capacitated arc
                  coeff = constr->getCapVarCoeff( i, j, d );
                  if ( fabs( coeff ) > CoeffTolerance )
                  {
                     aux = new CoeffShadow( capArc->getShadow(), 0,
                                            coeff );
                     constraints[c].addCoeff( aux );
                  }
               }
            }
         }

         // mark and save the adjacent nodes that are still not explored
         if ( head->getAuxCount() != searchCount )
         {
            head->setAuxCount( searchCount );
            head->setAuxNext( first );
            queue->setAuxNext( head );
            queue = head;
         }
      }

      // remove the current node from the queue
      if ( first == queue ) break;
      queue->setAuxNext( first->getAuxNext() );
   }
}

void PricingSolverGraph::updateReducCosts( const double* duals )
{
   // search the graph at the arc level to reset the reduced costs
   searchCount++;
   ArcNode* queue = root;
   root->setAuxNext( root );
   root->setAuxCount( searchCount );
   for ( ;; )
   {
      // get the first node in the queue
      ArcNode* first = queue->getAuxNext();
      first->setAuxCount( searchCount );

      // check all adjacent arcs
      for ( int a = 0; a < first->getOutDegree(); a++ )
      {
         // reset the reduced cost of the arc without capacity
         ArcVariable* arc = first->getOut( a );
         arc->resetReducCost();

         for ( int ca = 0; ca < arc->getNumCaps(); ca++ )
         {
            // reset the reduced cost of the capacitated arc
            CapArcVariable* capArc = arc->getCap( ca );
            capArc->resetReducCost();
         }

         // mark and save the adjacent nodes that are still not explored
         ArcNode* head = arc->getHead();
         if ( head->getAuxCount() != searchCount )
         {
            head->setAuxCount( searchCount );
            head->setAuxNext( first );
            queue->setAuxNext( head );
            queue = head;
         }
      }

      // remove the current node from the queue
      if ( first == queue ) break;
      queue->setAuxNext( first->getAuxNext() );
   }

   // traverse all constraints with non-zero dual values
   //printf( "Duals:\n" );
   constLB = 0.0;
   for ( int c = 0; c < numConstr; c++ )
   {
      if ( fabs( duals[c] ) <= DualTolerance ) continue;
      //printf( "dual_%d = %g\n", c, duals[c] );

      constLB += duals[c] * constraints[c].getConstr()->getRhs();

      // traverse the constraint coefficients
      CoeffShadow* prev = 0;
      CoeffShadow* coeff = constraints[c].getCoeffList();
      while ( coeff != 0 )
      {
         // remove the coefficient if the original variable has been deleted
         if ( ( coeff->getArc() == 0 ) && ( coeff->getCapArc() == 0 ) )
         {
            constraints[c].removeCoeff( prev );
            delete coeff;
            if ( prev == 0 )
               coeff = constraints[c].getCoeffList();
            else
               coeff = prev->getNext();
         }

         // update the reduced cost otherwise
         else
         {
            if ( coeff->getArc() != 0 )
               coeff->getArc()->addToReducCost(
                  -coeff->getCoeff() * duals[c] );
            if ( coeff->getCapArc() != 0 )
               coeff->getCapArc()->addToReducCost(
                  -coeff->getCoeff() * duals[c] );
            prev = coeff;
            coeff = coeff->getNext();
         }
      }
   }
}

void PricingSolverGraph::solve(  double &rc, double &lowerBound )
{
   // call the recursive solving procedure
   arcCount = 0;
   searchCount++;
   rc = recursiveSolve( modCapStart );
   lowerBound = rc * ((double)_inst->machines()) + constLB;
}

void PricingSolverGraph::fixArcs( const double ub, const double remRC )
{
   //printf(" ub: %g\n ", ub );
   // allocate the array of arcs to be deleted
   toBeDel.reserve( arcCount );
   toBeDel.clear();

   // call the recursive fixing procedure
   int prevArcCount = arcCount;
   arcCount = 0;
   searchCount++;
   recursiveFixArcs( modCapEnd, ub, remRC );

   // delete the fixed arcs
   for ( int a = 0; a < arcCount; a++ )
   {
      ModCapArcVariable* arc = toBeDel[a]->getModArc();
      if ( arc != 0 ) delete arc;
      toBeDel[a]->deleteRef();
   }
   arcCount = prevArcCount - arcCount;
}

CapArcKey* PricingSolverGraph::getSolution( int& size )
{
   // allocate memory for the solution
   size = capStart->getCap() + 1;
   CapArcKey* sol = new CapArcKey[size];

   // extract the solution from the graph
   ModCapArcNode* node = modCapStart;
   int pos = 0;
   int k;
   for ( k = 0; k < size; k++ )
   {
      ModCapArcVariable* arc = node->getSolOut( pos ).lastArc;
      if (!arc)
         break;
      assert(arc->getTail());
      assert(arc->getHead());

      pos = node->getSolOut( pos ).nextPos;
      sol[k].i = node->getIndex();
      sol[k].j = arc->getHead()->getIndex();
      sol[k].d = arc->getCap();
      node = arc->getHead();
   }
   size = k;

   // return the solution
   return sol;
}

bool PricingSolverGraph::checkSolution( int* sol, int size, bool mustCheck )
{
   // find the solution in the graph
   ModCapArcNode* node = modCapEnd;
   for ( int k = size-1; k >= -1; k-- )
   {
      // get the next solution node
      int j = 0;
      if ( k >= 0 ) j = sol[k];

      // find the arc that leads to the next solution node
      int pos;
      ModCapArcVariable* arc = 0;
      for ( pos = 0; pos < node->getInDegree(); pos++ )
      {
         arc = node->getIn( pos );
         if ( arc->getTail()->getIndex() == j ) break;
      }

      // throw an exception if did not find
      if ( pos == node->getInDegree() )
      {
         if (mustCheck)
         {
            fprintf( stderr, "Checking Error: arc (%d,%d,%d) of solution fixed to "
                    "zero.\n", j, node->getIndex(), node->getCap() );
            if (!checkArc(j, node->getIndex(), node->getCap()))
               fprintf( stderr, "Warning: the arc is preprocessed!\n" );
            exit( -1 );
         }
         else
            return false;
      }

      // go to the next node
      node = arc->getTail();
   }
   return true;
}

double PricingSolverGraph::searchFeasible( const double ub, CapArcKey* sol )
{
   double newUB = ub;
   bool found = false;
   int size = _inst->jobs() + _inst->machines();
   std::vector<bool> visited( _inst->jobs()+1 );
   std::vector<bool> blocked( _inst->jobs()+1 );
   std::vector<double> blockCount( _inst->jobs()+1, 1.0 );
   std::vector<ModCapArcNode*> currNodes( _inst->machines() );
   int bestCap = _inst->T();

   // repeat a number of tries
   for (int r = 0; r < numTriesFeasible; r++)
   {
      // clean the array of visited jobs
      for (int j = 0; j <= _inst->jobs(); j++)
         visited[j] = false;

      // start all machines at the root node
      for (int i = 0; i < _inst->machines(); i++)
         currNodes[i] = modCapStart;

      // initialized the current costs
      double origCost = 0;
      double lagCost = constLB;

      // build the solution arc by arc
      found = true;
      for (int a = 0; a < size; a++)
      {
         // find the best node to extend and calculate the current lower bound
         double lb = lagCost + currNodes[0]->getSolOut(0).cost;
         double origLB = origCost + currNodes[0]->getSolOut(0).origCost;
         int ii = 0;
         for (int i = 1; i < _inst->machines(); i++)
         {
            lb += currNodes[i]->getSolOut(0).cost;
            origLB += currNodes[i]->getSolOut(0).origCost;
            if (currNodes[i]->getCap() > currNodes[ii]->getCap())
               ii = i;
         }
         //fprintf( stderr, "Extending node %d^%d, LB = %g, UB = %g\n",
         //      currNodes[ii]->getIndex(), currNodes[ii]->getCap(), lb, ub );
         lb -= currNodes[ii]->getSolOut(0).cost;
         origLB -= currNodes[ii]->getSolOut(0).origCost;

         // clean the array of blocked jobs
         for (int j = 0; j <= _inst->jobs(); j++)
            blocked[j] = false;

         // find the best arc to use in the extension
         ModCapArcVariable* bestArc = 0;
         double bestCost = 1000000.0;
         for (int k = 0; k < currNodes[ii]->getOutDegree(); k++)
         {
            ModCapArcVariable* arc = currNodes[ii]->getOut(k);
            double cost;
            cost = origLB + arc->getCapArc()->getOrigCost() +
                  arc->getHead()->getSolOut(0).origCost;
            if (cost > ub - 1.0 + IntegerTolerance ) continue;
            cost = lb + arc->getCapArc()->getReducCost() +
                  arc->getCapArc()->getArc()->getReducCost() +
                  arc->getHead()->getSolOut(0).cost;
            if (cost > ub - 1.0 + IntegerTolerance ) continue;
            if (visited[arc->getHead()->getIndex()])
            {
               blocked[arc->getHead()->getIndex()] = true;
               continue;
            }
            //cost += double(rand() % 1000) * (ub - cost) / 1000.0;
            cost += double(rand() % 1000) * (ub - cost) * (0.5*blockCount[arc->getHead()->getIndex()]+0.5) / 1000.0;
            //cost = double(rand() % 1000) * blockCount[arc->getHead()->getIndex()];
            if (cost < bestCost)
            {
               bestCost = cost;
               bestArc = arc;
            }
         }
         if (bestArc == 0)
         {
            found = false;
            //fprintf( stderr, "Stopped!\n" );

            // update the blocking counters and stop
            for (int j = 0; j <= _inst->jobs(); j++)
            {
               blockCount[j] *= 0.9;
               if (blocked[j]) blockCount[j] += 0.1;
            }
            break;
         }
         //fprintf( stderr, "Found arc to %d^%d, lb = %g\n",
         //      bestArc->getHead()->getIndex(), bestArc->getHead()->getCap(),
         //      bestCost );

         // update the current state and costs
         lagCost += bestArc->getCapArc()->getReducCost() +
            bestArc->getCapArc()->getArc()->getReducCost();
         origCost += bestArc->getCapArc()->getOrigCost();
         currNodes[ii] = bestArc->getHead();
         if (bestArc->getHead()->getIndex() != 0)
            visited[bestArc->getHead()->getIndex()] = true;
         sol[a].i = bestArc->getTail()->getIndex();
         sol[a].j = bestArc->getHead()->getIndex();
         sol[a].d = bestArc->getCap();
         if (sol[a].d < bestCap) bestCap = sol[a].d;
      }
      if ((found) && (origCost < ub - 1.0 + IntegerTolerance))
      {
         newUB = origCost;

         // print the solution
         fprintf( stderr, "NEW BEST FEASIBLE SOLUTION OF %g FOUND!\n", origCost );
         int lastStart = -1;
         for (int i = 0; i < _inst->machines(); i++)
         {
            // find the next start
            for (int a = lastStart+1; a < size; a++)
               if (sol[a].i == 0)
               {
                  lastStart = a;
                  break;
               }

            // print the solution for the current machine
            int aa = lastStart;
            while (sol[aa].j != 0)
            {
               fprintf( stderr, " %d", sol[aa].j );

               // find the next job executed
               for (int a = aa+1; a < size; a++)
                  if (sol[a].i == sol[aa].j)
                  {
                     aa = a;
                     break;
                  }
            }
            fprintf( stderr, ".\n" );
         }
         break;
      }
   }

   // return the new upper bound if found
   if (found)
      return newUB;
   else
   {
      fprintf( stderr, "\nHeuristic got infeasible until t=%d (from T=%d)"
            " in %d tries.\n", bestCap, _inst->T(), numTriesFeasible );
      return ub;
   }
}

//#define CONV(x)   (((x)==0)? _inst->jobs(): ((x)-1))
//#define CONV(x)   (((x)==0)? _inst->jobs()+1: (x))
#define CONV(x)   (x)

void PricingSolverGraph::generateLpFile( int nodeNumber )
{
   // add the node number to the LP filename
   char filename[120];
   sprintf( filename, "%s_node-%d.lp", _inst->getName(), nodeNumber );

   // open the LP file
   FILE* f = fopen( filename, "wt" );
   if (f == NULL)
   {
      fprintf( stderr, "Cannot open LP file for writting!\n" );
      throw( -3 );
   }

   // initialize the LP and auxiliary data structures
   fprintf( f, "Minimize\n  Obj: " );
   std::vector<CapArcNode*> capNodes;
   std::vector<ArcVariable*> arcVars;

   // write the objective function and store the capacitated nodes and arcs
   bool isFirst = true;
   searchCount++;
   int lineCount = 0;
   ArcNode* queue = root;
   root->setAuxNext( root );
   root->setAuxCount( searchCount );
   for ( ;; )
   {
      // get the first node in the queue
      ArcNode* first = queue->getAuxNext();
      first->setAuxCount( searchCount );

      // check all adjacent arcs
      for ( int a = 0; a < first->getOutDegree(); a++ )
      {
         // get all capacitated arcs associated to the arc without capacity
         ArcVariable* arc = first->getOut( a );
         for ( int ca = 0; ca < arc->getNumCaps(); ca++ )
         {
            // write the cost of the current arc with capacity if non-zero
            CapArcVariable* capArc = arc->getCap( ca );
            if (fabs(capArc->getOrigCost()) > LpCoeffEps)
            {
               if (!isFirst) fprintf( f, " + " );
               fprintf( f, "%g w%d_%d_%d", capArc->getOrigCost(),
                     CONV(capArc->getTail()->getIndex()),
                     CONV(capArc->getHead()->getIndex()), capArc->getCap() );
               isFirst = false;
               lineCount++;
               if (lineCount == 10)
               {
                  fprintf( f, "\n  " );
                  lineCount = 0;
               }
            }

            // store the adjacent CapArcNodes
            if (capArc->getMod(0)->getHead()->getAuxCount() != searchCount)
            {
               capNodes.push_back( capArc->getHead() );
               capArc->getMod(0)->getHead()->setAuxCount( searchCount );
            }
            if (capArc->getMod(0)->getTail()->getAuxCount() != searchCount)
            {
               capNodes.push_back( capArc->getTail() );
               capArc->getMod(0)->getTail()->setAuxCount( searchCount );
            }
         }

         // store the ArcVariables
         arcVars.push_back( arc );

         // mark and save the adjacent nodes that are still not explored
         ArcNode* head = arc->getHead();
         if ( head->getAuxCount() != searchCount )
         {
            head->setAuxCount( searchCount );
            head->setAuxNext( first );
            queue->setAuxNext( head );
            queue = head;
         }
      }

      // remove the current node from the queue
      if ( first == queue ) break;
      queue->setAuxNext( first->getAuxNext() );
   }

   // write the flow conservation constraints
   fprintf( f, "\n\nSubject To\n\n" );
   for (int cn = 0; cn < (int)capNodes.size(); cn++)
   {
      // skip the root node
      if (capNodes[cn]->getIndex() == 0) continue;

      // initialize the contraint
      isFirst = true;
      fprintf( f, "flow%d_%d: ", CONV(capNodes[cn]->getIndex()),
            capNodes[cn]->getCap() );
      lineCount = 0;

      // write the out going flow terms
      for (int out = 0; out < capNodes[cn]->getOutDegree(); out++)
      {
         CapArcVariable* capArc = capNodes[cn]->getOut(out);
         if (!isFirst) fprintf( f, " + " );
         isFirst = false;
         fprintf( f, "w%d_%d_%d", CONV(capArc->getTail()->getIndex()),
               CONV(capArc->getHead()->getIndex()), capArc->getCap() );
         lineCount++;
         if (lineCount == 10)
         {
            fprintf( f, "\n  " );
            lineCount = 0;
         }
      }

      // write the incomming flow terms
      for (int in = 0; in < capNodes[cn]->getInDegree(); in++)
      {
         CapArcVariable* capArc = capNodes[cn]->getIn(in);
         isFirst = false;
         fprintf( f, " - w%d_%d_%d", CONV(capArc->getTail()->getIndex()),
               CONV(capArc->getHead()->getIndex()), capArc->getCap() );
         lineCount++;
         if (lineCount == 10)
         {
            fprintf( f, "\n  " );
            lineCount = 0;
         }
      }

      // finish the constraint
      fprintf( f, " = 0\n" );
   }
   fprintf( f, "\n" );

   // write the arc constaints
   for (int a = 0; a < (int)arcVars.size(); a++)
   {
      // initialize the contraint
      int i = CONV(arcVars[a]->getTail()->getIndex());
      int j = CONV(arcVars[a]->getHead()->getIndex());
      fprintf( f, "arc%d_%d: x%d_%d", i, j, i, j );
      lineCount = 0;

      // get all capacitated arcs associated to the arc without capacity
      for ( int ca = 0; ca < arcVars[a]->getNumCaps(); ca++ )
      {
         // write the term of the current capacitated arc
         CapArcVariable* capArc = arcVars[a]->getCap( ca );
         fprintf( f, " - w%d_%d_%d", i, j, capArc->getCap() );
         lineCount++;
         if (lineCount == 10)
         {
            fprintf( f, "\n  " );
            lineCount = 0;
         }
      }
      fprintf( f, " = 0\n" );
   }
   fprintf( f, "\n" );

   // write the constraint on the number of machines if necessary
   if (_inst->machines() > 1)
   {
      // initialize the constraint
      fprintf( f, "nmachines: " );
      isFirst = true;
      lineCount = 0;

      // check all adjacent arcs
      for ( int a = 0; a < root->getOutDegree(); a++ )
      {
         // get all capacitated arcs associated to the arc without capacity
         if (!isFirst) fprintf( f, " + " );
         isFirst = false;
         int i = CONV(root->getOut(a)->getTail()->getIndex());
         int j = CONV(root->getOut(a)->getHead()->getIndex());
         fprintf( f, "x%d_%d", i, j );
         lineCount++;
         if (lineCount == 10)
         {
            fprintf( f, "\n  " );
            lineCount = 0;
         }
      }
      fprintf( f, " = %d\n\n", _inst->machines() );
   }

   // write the model constraints:
   // - constraint on the number of machines
   // - degree constraints
   // - constraints associated to cutting planes
   // - constraints used for branching
   for ( int c = 0; c < numConstr; c++ )
   {
      // initialize the contraint
      isFirst = true;
      fprintf( f, "model%d: ", c );
      lineCount = 0;

      // traverse the constraint coefficients
      CoeffShadow* prev = 0;
      CoeffShadow* coeff = constraints[c].getCoeffList();
      while ( coeff != 0 )
      {
         // remove the coefficient if the original variable has been deleted
         if ( ( coeff->getArc() == 0 ) && ( coeff->getCapArc() == 0 ) )
         {
            constraints[c].removeCoeff( prev );
            delete coeff;
            if ( prev == 0 )
               coeff = constraints[c].getCoeffList();
            else
               coeff = prev->getNext();
         }

         // write the coefficient otherwise
         else
         {
            // write the coefficient sign
            if (!isFirst)
            {
               if (coeff->getCoeff() < 0.0)
                  fprintf( f, " - " );
               else
                  fprintf( f, " + " );
            }
            else if (coeff->getCoeff() < 0.0)
               fprintf( f, " -" );
            isFirst = false;

            // write the coefficient value and the variable
            if ( coeff->getArc() != 0 )
            {
               fprintf( f, "%.10lg x%d_%d", fabs(coeff->getCoeff()),
                     CONV(coeff->getArc()->getTail()->getIndex()),
                     CONV(coeff->getArc()->getHead()->getIndex()) );
            }
            if ( coeff->getCapArc() != 0 )
            {
               fprintf( f, "%.10lg w%d_%d_%d", fabs(coeff->getCoeff()),
                     CONV(coeff->getCapArc()->getTail()->getIndex()),
                     CONV(coeff->getCapArc()->getHead()->getIndex()),
                     coeff->getCapArc()->getCap() );
            }
            lineCount++;
            if (lineCount == 10)
            {
               fprintf( f, "\n  " );
               lineCount = 0;
            }

            // move to the next coefficient
            prev = coeff;
            coeff = coeff->getNext();
         }
      }

      // finish the constraint
      char ctrType = constraints[c].getConstr()->getType();
      switch( ctrType )
      {
      case '<': fprintf( f, " <=" ); break;
      case '>': fprintf( f, " >=" ); break;
      case '=': fprintf( f, " =" ); break;
      default:
         fprintf(stderr, "Unknown constraint type %c.\n", ctrType);
         throw(-4);
      }
      fprintf( f, " %.10lg\n", constraints[c].getConstr()->getRhs() );
   }
   fprintf( f, "\n" );

   // write the demands
   for (int i = 1; i <= _inst->jobs(); i++)
   {
      fprintf( f, "demand%d: d%d_%d = 1\n", CONV(i), CONV(i),
            _inst->ptime()[i] );
   }

   // write the variable types
   fprintf( f, "\nBinaries\n" );
   for (int a = 0; a < (int)arcVars.size(); a++)
   {
      // write the arc without capacity
      int i = CONV(arcVars[a]->getTail()->getIndex());
      int j = CONV(arcVars[a]->getHead()->getIndex());
      fprintf( f, "x%d_%d\n", i, j );

      // get all capacitated arcs associated to the arc without capacity
      for ( int ca = 0; ca < arcVars[a]->getNumCaps(); ca++ )
      {
         // write the capacitated arc
         CapArcVariable* capArc = arcVars[a]->getCap( ca );
         fprintf( f, "w%d_%d_%d\n", i, j, capArc->getCap() );
      }
   }
   fprintf( f, "\nEnd\n" );

   // close the LP file
   fclose( f );
   fprintf(stderr, "GENERATED ARC-TIME-INDEXED LP FILE %s WITH %d ARCS.\n", filename, arcCount );
}

/*** PRIVATE METHODS ***/

inline bool PricingSolverGraph::checkArc(int i, int j, int t)
{
   const int deltai   = _inst->ptime()[i] - _inst->ptime()[j];
   const int endp = t + _inst->ptime()[i];
   const int costJinHead = _inst->getCost( j, t );
   const int costIinHead = _inst->getCost( i, t + deltai );

   const int costJinTail = _inst->getCost( j, endp );
   const int costIinTail = _inst->getCost( i, endp );

   if ( (costJinHead+costIinTail) > (costIinHead+costJinTail) )
      return false;
   else
   {
      // breaking ties if necessary
      if ( (costJinHead+costIinTail) == (costIinHead+costJinTail) )
      {
         if ( _inst->duedate()[j] > _inst->duedate()[i] )
            return false;
         else
            if ( _inst->duedate()[j] == _inst->duedate()[i] )
               if ( j > i )
                  return false;
      }
   }
   return true;
}

void PricingSolverGraph::updateLowerBounds(double** fwd, double** bwd,
      const double* duals)
{
   // all jobs, including dummy
   const int jobs    = _inst->jobs() + 1;
   const int lastJob = jobs-1;
   // all times, including zero
   const int T       = _inst->T() + 1;
   const int lastT   = T-1;

   // fill all states forward
   fwd[0][0] = 0;
   for (int t = lastT; t > 0; t--)
      fwd[t][0] = DoubleInfinity;
   for (int t = 1; t < T; t++)
   {
      for (int i = 1; i <= lastJob; i++)
      {
         fwd[t][i] = DoubleInfinity;
         if (t < _inst->ptime()[i]) continue;
         for (int j = 0; j <= lastJob; j++)
         {
            double cost = 0.0;
            if ((j == 0) && (t > _inst->ptime()[i])) continue;
            if ((j > 0) && (t == _inst->ptime()[i])) continue;
            if (j > 0) cost = _inst->getCost(j, t - _inst->ptime()[i]) - duals[j-1];
            if (t == _inst->ptime()[i])
            {
               fwd[t][i] = 0.0;
               continue;
            }
            if (i == j) continue;
            if ( !checkArc(i, j, t - _inst->ptime()[i]) ) continue;
            double cost2 = cost + fwd[t - _inst->ptime()[i]][j];
            if (cost2 < fwd[t][i])
               fwd[t][i] = cost2;
         }
      }
   }

   // fill all states backward
   for (int t = lastT; t >= 0; t--)
      bwd[t][0] = 0;
   for (int t = lastT-1; t > 0; t--)
   {
      for (int j = 1; j <= lastJob; j++)
      {
         bwd[t][j] = DoubleInfinity;
         if ((lastT - t) < _inst->ptime()[j]) continue;
         double cost = _inst->getCost(j, t + _inst->ptime()[j]) - duals[j-1];
         if ((_inst->machines() > 1) || ((lastT - t) == _inst->ptime()[j]))
         {
            bwd[t][j] = cost;

            // fwd[t][0] should store the best pseudo-schedule that finishes at time t
            if (fwd[t + _inst->ptime()[j]][0] > cost + fwd[t + _inst->ptime()[j]][j])
               fwd[t + _inst->ptime()[j]][0] = cost + fwd[t + _inst->ptime()[j]][j];

            if ((lastT - t) == _inst->ptime()[j]) continue;
         }
         for (int i = 1; i <= lastJob; i++)
         {
            if (i == j) continue;
            if ( !checkArc(i, j, t + _inst->ptime()[j]) ) continue;
            double cost2 = cost + bwd[t + _inst->ptime()[j]][i];
            if (cost2 < bwd[t][j])
               bwd[t][j] = cost2;
         }
      }
   }
}

long double PricingSolverGraph::recursiveSolve( ModCapArcNode* start )
{
   // process the node if still not processed
   if ( start->getAuxCount() != searchCount )
   {
      // initialize de current cost as empty
      SolutionDescr sol;
      sol.cost = DoubleInfinity;
      sol.origCost = DoubleInfinity;
      sol.lastArc = 0;
      sol.nextPos = 0;

      // handle the base case
      if ( start == modCapEnd )
      {
         sol.cost = 0;
         sol.origCost = 0;
      }
      else
      {
         // check all outgoing arcs
         for ( int pos = 0; pos < start->getOutDegree(); pos++ )
         {
            ModCapArcVariable* arc = start->getOut( pos );
            ModCapArcNode* node = arc->getHead();
            arcCount++;

            // process the adjacent node
            recursiveSolve( node );

            // check if the current solution is better
            long double cost = node->getSolOut( 0 ).cost +
                          arc->getCapArc()->getReducCost() +
                          arc->getCapArc()->getArc()->getReducCost();
            if ( cost < sol.cost )
            {
               sol.cost = cost;
               sol.lastArc = arc;
            }

            // check if the current solution is better with original costs
            double origCost = node->getSolOut( 0 ).origCost +
                          arc->getCapArc()->getOrigCost();
            if ( origCost < sol.origCost )
               sol.origCost = origCost;
         }
      }

      // set the solution and mark the node as processed
      start->setSolOut( 0, sol );
      start->setAuxCount( searchCount );
   }

   // return the solution cost
   return start->getSolOut( 0 ).cost;
}

long double PricingSolverGraph::recursiveFixArcs( ModCapArcNode* end,
      const double ub, const double remRC )
{
   // process the node if still not processed
   if ( end->getAuxCount() != searchCount )
   {
      // initialize de current cost as empty
      SolutionDescr sol;
      sol.cost = DoubleInfinity;
      sol.origCost = DoubleInfinity;
      sol.lastArc = 0;
      sol.nextPos = 0;
      double remOrigCost =
            modCapStart->getSolOut(0).origCost * double(_inst->machines()-1);

      // handle the base case
      if ( end == modCapStart )
      {
         sol.cost = constLB + remRC*((double)(_inst->machines()-1));
         sol.origCost = remOrigCost;
      }
      else
      {
         // check all incoming arcs
         for ( int pos = 0; pos < end->getInDegree(); pos++ )
         {
            ModCapArcVariable* arc = end->getIn( pos );
            ModCapArcNode* node = arc->getTail();

            // process the adjacent node
            long double cost = recursiveFixArcs( node, ub, remRC );

            // check if the arc variable can be fixed to zero
            cost += end->getSolOut(0).cost +
                  arc->getCapArc()->getReducCost() +
                  arc->getCapArc()->getArc()->getReducCost();
            double origCost = node->getSolIn(0).origCost +
                  end->getSolOut(0).origCost +
                  arc->getCapArc()->getOrigCost();
            //if (cost > ub - 1.0 + IntegerTolerance)
            if ( ((cost > ub - 1.0 + IntegerTolerance) ||
                  (origCost > ub - 1.0 + IntegerTolerance))
#ifndef FIX_ROOT_ARCS
                  && (node->getIndex() != 0)
#endif
                  )
            {
               toBeDel.push_back(arc->getShadow());
               arc->getShadow()->addRef();
               arcCount++;
            }
            else
            {
               // check if the current solution is better
               cost -= end->getSolOut( 0 ).cost;
               if ( cost < sol.cost )
               {
                  sol.cost = cost;
                  sol.lastArc = arc;
               }

               // check if the current solution is better with original cost
               origCost -= end->getSolOut(0).origCost;
               if ( origCost < sol.origCost )
                  sol.origCost = origCost;
            }
         }
      }

      // set the solution and mark the node as processed
      end->setSolIn( 0, sol );
      end->setAuxCount( searchCount );
      assert( end );
      //printf(" j %d t %d %g \n", end->getIndex(), end->getCap(), end->getSolIn( 0 ).cost);
   }

   // return the solution cost
   return end->getSolIn( 0 ).cost;
}

/*****************************************************************************
 *                             PricingSolverArray                            *
 *****************************************************************************/

PricingSolverArray::PricingSolverArray( Instance *inst ) :
      _inst( inst )
{
   // all jobs, including dummy
   const int jobs    = inst->jobs() + 1;
   const int lastJob = jobs-1;
   // all times, including zero
   const int T       = inst->T() + 1;

   // allocate all arrays
   _duals = new double[lastJob];
   solution = new SolArrayDescr[T];

   // other initializations
   constLB = 0;
   arcCount = lastJob * (T-1);
}

PricingSolverArray::~PricingSolverArray()
{
   // release all arrays
   delete [] _duals;
   delete [] solution;
}

PricingSolver* PricingSolverArray::copy()
{
   return new PricingSolverArray(_inst);
}

void PricingSolverArray::updateReducCosts( const double* duals)
{
   // all jobs, including dummy
   const int jobs    = _inst->jobs() + 1;
   const int lastJob = jobs-1;
   constLB = 0;
   for (int j = 0; j < lastJob; j++)
   {
      _duals[j] = duals[j];
      constLB += _duals[j];
   }
}

void PricingSolverArray::solve( double &rc, double &lowerBound )
{
   // all jobs, including dummy
   const int jobs    = _inst->jobs() + 1;
   const int lastJob = jobs-1;
   // all times, including zero
   const int T       = _inst->T() + 1;
   const int lastT   = T-1;

   // fill all states
   solution[0].cost[0] = 0;
   solution[0].cost[1] = 0;
   solution[0].job[0] = 0;
   solution[0].job[1] = 0;
   solution[0].prev[0] = 0;
   solution[0].prev[1] = 0;
   for (int t = 1; t < T; t++)
   {
      solution[t].cost[0] = DoubleInfinity;
      solution[t].cost[1] = DoubleInfinity;
      solution[t].job[0] = 0;
      solution[t].job[1] = 0;
      solution[t].prev[0] = 0;
      solution[t].prev[1] = 0;
      for (int j = 1; j <= lastJob; j++)
      {
         if (t < _inst->ptime()[j]) continue;
         double cost = _inst->getCost(j, t) - _duals[j-1];
         int prev;
         if (solution[t - _inst->ptime()[j]].job[0] == j)
            prev = 1;
         else
            prev = 0;
         cost += solution[t - _inst->ptime()[j]].cost[prev];
         if (cost < solution[t].cost[0])
         {
            solution[t].cost[1] = solution[t].cost[0];
            solution[t].job[1] = solution[t].job[0];
            solution[t].prev[1] = solution[t].prev[0];
            solution[t].cost[0] = cost;
            solution[t].job[0] = j;
            solution[t].prev[0] = prev;
         }
         else if (cost < solution[t].cost[1])
         {
            solution[t].cost[1] = cost;
            solution[t].job[1] = j;
            solution[t].prev[1] = prev;
         }
      }
   }

   // Find the last completion time of the best solution
   int bestT = lastT;
   if (_inst->machines() > 1)
      for (int t = 0; t < lastT; t++)
      {
         if (solution[t].cost[0] < solution[bestT].cost[0])
            bestT = t;
      }
   rc = solution[bestT].cost[0];
   lowerBound = rc * ((double)_inst->machines()) + constLB;
}

void PricingSolverArray::fixArcs( const double ub, const double remRC )
{
   // Do nothing!
}

CapArcKey* PricingSolverArray::getSolution(int& size)
{
   // Find the last completion time of the best solution
   int bestT = _inst->T();
   if (_inst->machines() > 1)
      for (int t = 0; t < _inst->T(); t++)
      {
         if (solution[t].cost[0] < solution[bestT].cost[0])
            bestT = t;
      }

   // count the number of jobs in the subproblem solution
   size = 1;
   int t = bestT;
   int ii = 0;
   int prev;
   while (t > 0)
   {
      size++;
      prev = solution[t].prev[ii];
      t -= _inst->ptime()[solution[t].job[ii]];
      ii = prev;
   }

   // recover the solution and return
   CapArcKey* sol = new CapArcKey[size];
   int k = 0;
   t = bestT;
   ii = 0;
   int i = 0;
   while (t > 0)
   {
      sol[k].d = t;
      sol[k].i = i;
      sol[k].j = solution[t].job[ii];
      k++;
      i = solution[t].job[ii];
      prev = solution[t].prev[ii];
      t -= _inst->ptime()[solution[t].job[ii]];
      ii = prev;
   }
   sol[k].d = 0;
   sol[k].i = i;
   sol[k].j = 0;
   return sol;
}

bool PricingSolverArray::checkSolution(int* sol, int size, bool mustCheck)
{
   // Do nothing!

   return true;
}

double PricingSolverArray::searchFeasible( const double ub, CapArcKey* sol )
{
   // Do nothing!

   return ub;
}

void PricingSolverArray::generateLpFile( int nodeNumber )
{
   fprintf( stderr, "LP generation not implemented for time-indexed pricing.\n" );
   throw( -2 );
}


void PricingSolverArray::generateTimeIndexedLpFile( int nodeNumber )
{
   fprintf( stderr, "LP generation not implemented for time-indexed pricing.\n" );
   throw( -2 );
}

void PricingSolverArray::generateArcTimeFixDataFile( int nodeNumber )
{
   fprintf( stderr, "LP generation not implemented for time-indexed pricing.\n" );
   throw( -2 );
}


//void PricingSolverGraph::generateTimeIndexedLpFile( int nodeNumber )
//{
//   // add the node number to the LP filename
//   char filename[120];
//   sprintf( filename, "%s_node-%d_time_indexed.lp", _inst->getName(), nodeNumber );
//
//   // open the LP file
//   FILE* f = fopen( filename, "wt" );
//   if (f == NULL)
//   {
//      fprintf( stderr, "Cannot open LP file for writting!\n" );
//      throw( -3 );
//   }
//
//   // initialize auxiliary data structures
//   std::vector<CapArcNode*> capNodes;
//   std::vector<ArcVariable*> arcVars;
//   
//   std::vector<std::vector<bool> > unfixedJobTime(_inst->jobs() + 1);
//   std::vector<std::vector<double> > unfixedJobTimeCost(_inst->jobs() + 1);
//
//   for (int i = 1; i < (int)unfixedJobTimeCost.size(); i++)
//		unfixedJobTime[i].resize(_inst->T() + 1);
//	
//   // write the objective function and store the capacitated nodes and arcs
//   searchCount++;
//   ArcNode* queue = root;
//   root->setAuxNext( root );
//   root->setAuxCount( searchCount );
//   for ( ;; )
//   {
//      // get the first node in the queue
//      ArcNode* first = queue->getAuxNext();
//      first->setAuxCount( searchCount );
//
//      // check all adjacent arcs
//      for ( int a = 0; a < first->getOutDegree(); a++ )
//      {
//         // get all capacitated arcs associated to the arc without capacity
//         ArcVariable* arc = first->getOut( a );
//         for ( int ca = 0; ca < arc->getNumCaps(); ca++ )
//         {
//            // write the cost of the current arc with capacity if non-zero
//            CapArcVariable* capArc = arc->getCap( ca );
//
//			//if (CONV(capArc->getHead()->getIndex()) != 0) //verifica indice i do arco
//			//	unfixedJobTime[CONV(capArc->getHead()->getIndex())][capArc->getCap()] = true;
//
//            // store the adjacent CapArcNodes
//            if (capArc->getMod(0)->getHead()->getAuxCount() != searchCount)
//            {
//               capNodes.push_back( capArc->getHead() );
//               capArc->getMod(0)->getHead()->setAuxCount( searchCount );
//            }
//            if (capArc->getMod(0)->getTail()->getAuxCount() != searchCount)
//            {
//               capNodes.push_back( capArc->getTail() );
//               capArc->getMod(0)->getTail()->setAuxCount( searchCount );
//            }
//         }
//
//         // store the ArcVariables
//         arcVars.push_back( arc );
//
//         // mark and save the adjacent nodes that are still not explored
//         ArcNode* head = arc->getHead();
//         if ( head->getAuxCount() != searchCount )
//         {
//            head->setAuxCount( searchCount );
//            head->setAuxNext( first );
//            queue->setAuxNext( head );
//            queue = head;
//         }
//      }
//
//      // remove the current node from the queue
//      if ( first == queue ) break;
//      queue->setAuxNext( first->getAuxNext() );
//   }
//
//   for (int a = 0; a < (int)arcVars.size(); a++)
//   {
//      int i = CONV(arcVars[a]->getHead()->getIndex());
//      int j = CONV(arcVars[a]->getTail()->getIndex());
//
//      for ( int ca = 0; ca < arcVars[a]->getNumCaps(); ca++ )
//      {
//         CapArcVariable* capArc = arcVars[a]->getCap( ca );
//		 if (i != 0)//if (j != 0) //verifica indice j do arco
//			unfixedJobTime[i][capArc->getCap()] = true;
//      }
//   }
//
//   fprintf( f, "Minimize\n  Obj: " );
//
//   bool isFirst = true;
//   int lineCount = 0;
//
//   //inserido por Daniel em 13/08/2015
//   for (int j = 1; j <= _inst->jobs(); j++)
//   {
//	   for (int t = 0; t <= _inst->T(); t++)
//	   {
//		   if ( (unfixedJobTime[j][t]) && (_inst->getCost(j,t) > LpCoeffEps) ) //+ _inst->ptime()[j]
//		   {
//				if (!isFirst)
//					fprintf( f, " + " );
//				
//				fprintf( f, "%d y%d_%d", _inst->getCost(j,t), j, t);
//		   
//				isFirst = false;
//				lineCount++;
//				if (lineCount == 10)
//				{
//					fprintf( f, "\n  " );
//					lineCount = 0;
//				} 
//		   }
//	   }
//   }
//
//   //inserido por Daniel em 13/08/2015
//   // write the flow conservation constraints
//	fprintf( f, "\n\nSubject To\n\n" );
//
//	// initialize the contraint
//	isFirst = true;
//	fprintf( f, "m: ");
//	lineCount = 0;
//	for (int j = 1; j <= _inst->jobs(); j++)
//	{
//		if (unfixedJobTime[j][_inst->ptime()[j]])
//		{
//			if (!isFirst) fprintf( f, " + " );
//			fprintf( f, "y%d_%d", j, _inst->ptime()[j] );
//				
//			isFirst = false;
//			lineCount++;
//			if (lineCount == 10)
//			{
//			fprintf( f, "\n  " );
//			lineCount = 0;
//			}
//		}
//					
//	}
//	// finish the constraint
//	fprintf( f, " = %d\n", _inst->machines() );
//
//	//write demand constraints
//	for (int j = 1; j < unfixedJobTime.size(); j++)
//	{
//		// initialize the contraint
//		isFirst = true;
//		fprintf( f, "demand%d: ", j);
//		lineCount = 0;
//		for (int t = 0; t < unfixedJobTime[j].size(); t++)
//		{     
//			if (unfixedJobTime[j][t])
//			 {     
//				if (!isFirst) fprintf( f, " + " );
//				fprintf( f, "y%d_%d", j, t );
//				
//				isFirst = false;
//				lineCount++;
//				if (lineCount == 10)
//				{
//				fprintf( f, "\n  " );
//				lineCount = 0;
//				}
//			}
//		}
//		// finish the constraint
//		fprintf( f, " = 1\n" );
//	}
//
//	//write flow constraints
//	for (int t = 1; t <= _inst->T(); t++)
//	{ 
//		// initialize the contraint
//		isFirst = true;
//		fprintf( f, "flow%d: ", t);
//		lineCount = 0;
//		for (int j = 1; j < unfixedJobTime.size(); j++)
//		{
//			if (unfixedJobTime[j][t])
//			{ 
//				if (!isFirst) fprintf( f, " + " );
//				fprintf( f, "y%d_%d", j, t );
//				
//				isFirst = false;
//				lineCount++;
//				if (lineCount == 10)
//				{
//				fprintf( f, "\n  " );
//				lineCount = 0;
//				}
//			}
//
//			if ( (t + _inst->ptime()[j] <= _inst->T() ) && (unfixedJobTime[j][t + _inst->ptime()[j]] ) )
//			{ 
//				fprintf( f, " - y%d_%d", j, t + _inst->ptime()[j] );
//				
//				isFirst = false;
//				lineCount++;
//				if (lineCount == 10)
//				{
//					fprintf( f, "\n  " );
//					lineCount = 0;
//				}
//			}
//		}
//		// finish the constraint
//		fprintf( f, " > 0\n" );
//	}
//
//   int yCount = 0;
//
//   // write the variable types
//   fprintf( f, "\nBinaries\n" );
//
//	for (int j = 1; j <= _inst->jobs(); j++)
//	{	
//		for (int t = 0; t <= _inst->T(); t++)
//		{     
//			if (unfixedJobTime[j][t])			  
//			{
//				fprintf( f, "y%d_%d\n", j, t );						
//				yCount++;
//			}	
//		}
//	}
//
//   fprintf( f, "\nEnd\n" );
//
//   // close the LP file
//   fclose( f );
//   fprintf(stderr, "GENERATED TIME INDEXED LP FILE %s WITH %d VARIABLES.\n", filename, yCount );
//}



void PricingSolverGraph::generateTimeIndexedLpFile( int nodeNumber )
{
   // add the node number to the LP filename
   char filename[120];
   sprintf( filename, "schedwtdata/%s_node-%d_time_indexed.lp", _inst->getName(), nodeNumber );

   // open the LP file
   FILE* f = fopen( filename, "wt" );
   if (f == NULL)
   {
      fprintf( stderr, "Cannot open LP file for writting!\n" );
      throw( -3 );
   }

   // initialize auxiliary data structures
   std::vector<CapArcNode*> capNodes;
   std::vector<ArcVariable*> arcVars;
   
   std::vector<std::vector<bool> > unfixedJobTime(_inst->jobs() + 1);
   //std::vector<std::vector<double> > unfixedJobTimeCost(_inst->jobs() + 1);

   for (int i = 1; i <= _inst->jobs(); i++)
		unfixedJobTime[i].resize(_inst->T() + 1);
	
   // write the objective function and store the capacitated nodes and arcs
   searchCount++;
   ArcNode* queue = root;
   root->setAuxNext( root );
   root->setAuxCount( searchCount );
   for ( ;; )
   {
      // get the first node in the queue
      ArcNode* first = queue->getAuxNext();
      first->setAuxCount( searchCount );

      // check all adjacent arcs
      for ( int a = 0; a < first->getOutDegree(); a++ )
      {
         // get all capacitated arcs associated to the arc without capacity
         ArcVariable* arc = first->getOut( a );
         for ( int ca = 0; ca < arc->getNumCaps(); ca++ )
         {
            // write the cost of the current arc with capacity if non-zero
            CapArcVariable* capArc = arc->getCap( ca );

			//if (CONV(capArc->getHead()->getIndex()) != 0) //verifica indice i do arco
			//	unfixedJobTime[CONV(capArc->getHead()->getIndex())][capArc->getCap()] = true;

            // store the adjacent CapArcNodes
            if (capArc->getMod(0)->getHead()->getAuxCount() != searchCount)
            {
               capNodes.push_back( capArc->getHead() );
               capArc->getMod(0)->getHead()->setAuxCount( searchCount );
            }
            if (capArc->getMod(0)->getTail()->getAuxCount() != searchCount)
            {
               capNodes.push_back( capArc->getTail() );
               capArc->getMod(0)->getTail()->setAuxCount( searchCount );
            }
         }

         // store the ArcVariables
         arcVars.push_back( arc );

         // mark and save the adjacent nodes that are still not explored
         ArcNode* head = arc->getHead();
         if ( head->getAuxCount() != searchCount )
         {
            head->setAuxCount( searchCount );
            head->setAuxNext( first );
            queue->setAuxNext( head );
            queue = head;
         }
      }

      // remove the current node from the queue
      if ( first == queue ) break;
      queue->setAuxNext( first->getAuxNext() );
   }

   for (int a = 0; a < (int)arcVars.size(); a++)
   {
      int i = CONV(arcVars[a]->getHead()->getIndex());
      int j = CONV(arcVars[a]->getTail()->getIndex());

      for ( int ca = 0; ca < arcVars[a]->getNumCaps(); ca++ )
      {
         CapArcVariable* capArc = arcVars[a]->getCap( ca );
		 if (i != 0)//if (j != 0) //verifica indice j do arco
			unfixedJobTime[i][capArc->getCap()] = true;
      }
   }

   //populate the dominho
	std::vector<std::vector<int> > auxUnfixedTimes(_inst->jobs() + 1);	
	for (int j = 1; j <= _inst->jobs(); j++)
	{
		auxUnfixedTimes[j].push_back(_inst->ptime()[j]-1);
		for (int t = 1; t <= _inst->T(); t++)
		{				
			if (unfixedJobTime[j][t])			
				auxUnfixedTimes[j].push_back(t);			
		}
	}

   fprintf( f, "Minimize\n  Obj: " );

   for (int j = 1; j <= _inst->jobs(); j++)
   {
		for (int i = 1; i < auxUnfixedTimes[j].size(); i++)
		{
		if (_inst->getCost(j,auxUnfixedTimes[j][i]) > LpCoeffEps)  //+ _inst->ptime()[j]		   				
			fprintf( f, " + %d z%d_%d - %d z%d_%d", _inst->getCost(j,auxUnfixedTimes[j][i]), j, auxUnfixedTimes[j][i],
				_inst->getCost(j,auxUnfixedTimes[j][i]), j, auxUnfixedTimes[j][i-1] ); 
		}
   }

   // write the flow conservation constraints
	fprintf( f, "\n\nSubject To\n\n" );

	//write the domino constraints -ok
	for (int j = 1; j <= _inst->jobs(); j++)
	{		
		for (int i = 1; i < auxUnfixedTimes[j].size(); i++)
		{		
			// initialize the contraint
			fprintf( f, "domino%d_%d: ", j, i);
			
			fprintf( f, "z%d_%d", j, auxUnfixedTimes[j][i-1]);
			fprintf( f, " - z%d_%d", j, auxUnfixedTimes[j][i]);

			// finish the constraint
			fprintf( f, " < 0\n" );
		}
	}
	// end - write the domino constraints
	 
	//write the kick flow constraint -ok
	fprintf( f, "m: ");
	for (int j = 1; j <= _inst->jobs(); j++)
	{
		if (unfixedJobTime[j][_inst->ptime()[j]])
			fprintf( f, " + z%d_%d - z%d_%d", j, _inst->ptime()[j], j, _inst->ptime()[j]-1 );					
	}
	fprintf( f, " = %d\n", _inst->machines() );
	//end - write the kick flow constraint

	//write demand constraints -ok
	for (int j = 1; j <= _inst->jobs(); j++)
	{
		// initialize the contraint
		fprintf( f, "demand%d: ", j); 			 	
		fprintf( f, " z%d_%d - z%d_%d", j, auxUnfixedTimes[j][auxUnfixedTimes[j].size()-1], 
			j, auxUnfixedTimes[j][0] );			
		
		// finish the constraint
		fprintf( f, " = 1\n" );
	}
	//end - write demand constraints

	//write flow constraints
	for (int t = 1; t <= _inst->T(); t++)
	{
		// initialize the contraint
		fprintf( f, "flow%d: ", t);
		for (int j = 1; j <= _inst->jobs(); j++)
		{
			if (unfixedJobTime[j][t])	
			{
				fprintf( f, " + z%d_%d", j, t );

				//porcaria - trocar depois - busca o t nao fixado anterior
				for (int i = 0; i < (int)auxUnfixedTimes[j].size(); i++)
				{
					if (auxUnfixedTimes[j][i] == t)
					{					
						fprintf( f, " - z%d_%d", j, auxUnfixedTimes[j][i-1] );
						break;
					}
				}
				//porcaria - trocar depois
			}
		}
		fprintf( f, " ");

		for (int j = 1; j <= _inst->jobs(); j++)
		{
			if ( (t + _inst->ptime()[j] <= _inst->T() ) && (unfixedJobTime[j][t + _inst->ptime()[j]] ) )		
			{
				fprintf( f, " - z%d_%d", j, t + _inst->ptime()[j] );	
				
				//porcaria - trocar depois - busca o t nao fixado anterior
				for (int i = 0; i < (int)auxUnfixedTimes[j].size(); i++)
				{
					if (auxUnfixedTimes[j][i] == t + _inst->ptime()[j])
					{					
						fprintf( f, " + z%d_%d", j, auxUnfixedTimes[j][i-1] );
						break;
					}
				}
				//porcaria - trocar depois
			}
		}
		// finish the constraint
		fprintf( f, " > 0\n" );
	}
	//end - write flow constraints
	       
	int varCount = 0;

   // write the variable types
   fprintf( f, "\nBinaries\n" );

	for (int j = 1; j <= _inst->jobs(); j++)
	{			
		for (int i = 0; i < (int)auxUnfixedTimes[j].size(); i++)	  
		{				
			fprintf( f, "z%d_%d\n", j, auxUnfixedTimes[j][i] );	
			varCount++;
		}
	}

   fprintf( f, "\nEnd\n" );

   // close the LP file
   fclose( f );
   fprintf(stderr, "GENERATED TIME INDEXED LP FILE %s WITH %d VARIABLES.\n", filename, varCount );
}

















void PricingSolverGraph::generateArcTimeFixDataFile( int nodeNumber )
{
   // add the node number to the LP filename
   char filename[120];
   sprintf( filename, "schedwtdata/%s_node-%d_arc_time_fix_data.txt", _inst->getName(), nodeNumber );

   // open the LP file
   FILE* f = fopen( filename, "wt" );
   if (f == NULL)
   {
      fprintf( stderr, "Cannot open LP file for writting!\n" );
      throw( -3 );
   }

   // initialize auxiliary data structures
   std::vector<CapArcNode*> capNodes;
   std::vector<ArcVariable*> arcVars;
   
   // write the objective function and store the capacitated nodes and arcs
   searchCount++;
   ArcNode* queue = root;
   root->setAuxNext( root );
   root->setAuxCount( searchCount );
   for ( ;; )
   {
      // get the first node in the queue
      ArcNode* first = queue->getAuxNext();
      first->setAuxCount( searchCount );

      // check all adjacent arcs
      for ( int a = 0; a < first->getOutDegree(); a++ )
      {
         // get all capacitated arcs associated to the arc without capacity
         ArcVariable* arc = first->getOut( a );
         for ( int ca = 0; ca < arc->getNumCaps(); ca++ )
         {
            // write the cost of the current arc with capacity if non-zero
            CapArcVariable* capArc = arc->getCap( ca );

			//if (CONV(capArc->getHead()->getIndex()) != 0) //verifica indice i do arco
			//	unfixedJobTime[CONV(capArc->getHead()->getIndex())][capArc->getCap()] = true;

            // store the adjacent CapArcNodes
            if (capArc->getMod(0)->getHead()->getAuxCount() != searchCount)
            {
               capNodes.push_back( capArc->getHead() );
               capArc->getMod(0)->getHead()->setAuxCount( searchCount );
            }
            if (capArc->getMod(0)->getTail()->getAuxCount() != searchCount)
            {
               capNodes.push_back( capArc->getTail() );
               capArc->getMod(0)->getTail()->setAuxCount( searchCount );
            }
         }

         // store the ArcVariables
         arcVars.push_back( arc );

         // mark and save the adjacent nodes that are still not explored
         ArcNode* head = arc->getHead();
         if ( head->getAuxCount() != searchCount )
         {
            head->setAuxCount( searchCount );
            head->setAuxNext( first );
            queue->setAuxNext( head );
            queue = head;
         }
      }

      // remove the current node from the queue
      if ( first == queue ) break;
      queue->setAuxNext( first->getAuxNext() );
   }

   
	int varCount = 0;
	for (int a = 0; a < (int)arcVars.size(); a++)
	{
		int i = CONV(arcVars[a]->getHead()->getIndex());
		int j = CONV(arcVars[a]->getTail()->getIndex());

		for ( int ca = 0; ca < arcVars[a]->getNumCaps(); ca++ )
		{
			CapArcVariable* capArc = arcVars[a]->getCap( ca );
			
			fprintf( f, "%d %d %d\n", i, j, capArc->getCap() );
			varCount++;
		}	
	}


   // close the LP file
   fclose( f );
   fprintf(stderr, "GENERATED ARC TIME INDEXED FIXED VARIABLES DATA %s WITH %d VARIABLES.\n", filename, varCount );
}