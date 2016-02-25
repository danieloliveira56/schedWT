#include "cycleElimSep.h"

#include <list>

const double ExtCycleEliminator::Infinity = 1000.0;

const double ExtCycleEliminator::Tolerance = 1E-5;

const double ExtCycleEliminator::MinViolation = ECEC_VIOLATED_EPS;

ExtCycleEliminator::ExtCycleEliminator( LpSolution* solution, int* demands,
      int capacity, int* setList, int setSize )
{
   solution_ = solution;
   setList_ = setList;
   setSize_ = setSize;

   /** build the graph to calculate the minimum cut **/
   // put the source and destination
   Vertex v;
   v.index = -1;
   v.mark = 0;
   vertices.push_back( v );
   v.index = -2;
   v.mark = 0;
   vertices.push_back( v );

   // find the value of n search for the second demand zero.
   n = 1;
   while (demands[n] != 0) n++;
   // printf( "n = %d, cap = %d, set = (", n, capacity );

   // build a vector of indices of nodes in the set list
   std::vector<int> inSet;
   inSet.resize( n, 0 );
   int i, j, k, d;
   int totalDem = 0;
   for (k = 1; k <= setSize; k++)
   {
      // printf( "%d(%d)", setList[k], demands[setList[k]] );
      // if (k == setSize) printf(")\n");
      // else printf(",");
      inSet[setList[k]] = k;
      totalDem += demands[setList[k]];
   }

   // build a vector of booleans that indicate the possible differences of
   // capacities of entering arcs and capacities of leaving arcs
   // diffCap[i][j][d] = true means that the capacity of an arc entering in
   // setList[i] may differ from that of an arc leaving setList[j] by d.
   // Since diffCap[i][j][x] = diffCap[j][i][x], we only initialize diffCap
   // for i <= j.
   diffCap.resize(setSize+1);
   for (i = 1; i <= setSize; i++)
   {
      diffCap[i].resize(setSize+1);
      for (j = i; j <= setSize; j++)
      {
         if (totalDem < capacity)
            diffCap[i][j].resize(totalDem+1, false);
         else
            diffCap[i][j].resize(capacity+1, false);
         if (i == j)
            diffCap[i][j][demands[setList[i]]] = true;
         else
         {
            // set all demands using dynamic programming
            d = demands[setList[i]]+demands[setList[j]];
            if ( d < (int) diffCap[i][j].size())
               diffCap[i][j][d] = true;
            for (k = 1; k <= setSize; k++)
            {
               if ((i == k) || (j == k)) continue;
               for (d = diffCap[i][j].size()-1; d >= 0; d--)
               {
                  if ( diffCap[i][j][d] && ((d+demands[setList[k]])
                        < (int) diffCap[i][j].size()) )
                     diffCap[i][j][d+demands[setList[k]]] = true;
               }
            }
         }
      }
   }

   // put the arcs that enter the set
   double x;
   minCut = 0.0;
   for (k = 0; k < solution_->numArcsCap; k++)
   {
      i = solution_->arcsCap[k].i;
      j = solution_->arcsCap[k].j;
      d = solution_->arcsCap[k].d;
      x = solution_->arcsCap[k].value;
      if ( (inSet[j] > 0) && (inSet[i] == 0) )
      {
         // add the arc from the source to the new vertex
         vertices[0].outVertices.push_back( vertices.size() );
         vertices[0].outReturns.push_back( -1 );
         vertices[0].outCaps.push_back( x );

         v.index = k;
         v.mark = 0;
         vertices.push_back( v );

         minCut += x;
         // printf( "entering arc (%d,%d)^%d with flow %g\n", i, j, d, x );
      }
   }

   // put the arcs that leave the set
   firstLeaving = vertices.size();
   int l, ii, jj, dd, pos;
   for (k = 0; k < solution_->numArcsCap; k++)
   {
      i = solution_->arcsCap[k].i;
      j = solution_->arcsCap[k].j;
      d = solution_->arcsCap[k].d;
      x = solution_->arcsCap[k].value;
      if ((inSet[i] > 0) && (inSet[j] == 0))
      {
         v.index = k;
         v.mark = 0;
         pos = vertices.size();
         vertices.push_back( v );

         // add the arc from the new vertex to the sink
         vertices[pos].outVertices.push_back( 1 );
         vertices[pos].outReturns.push_back( -1 );
         vertices[pos].outCaps.push_back( x );
         // printf( "leaving arc (%d,%d)^%d with flow %g\n", i, j, d, x );

         // add the arcs from matching entering arcs
         for (l = 2; l < firstLeaving; l++)
         {
            ii = solution_->arcsCap[vertices[l].index].i;
            jj = solution_->arcsCap[vertices[l].index].j;
            dd = solution_->arcsCap[vertices[l].index].d;

            // check if the arc can match the current leaving arc
            int iii, jjj;
            iii = inSet[jj];
            jjj = inSet[i];
            if (iii > jjj)
            {
               iii = inSet[i];
               jjj = inSet[jj];
            }
            if ( ((dd - d) > 0) && ((dd - d) < (int) diffCap[iii][jjj].size()) )
               if (diffCap[iii][jjj][dd - d] && ((ii != j) || (ii == 0)) )
               {
                  // add an arc of infinity capacity and a backward arc of
                  // null capacity
                  int aux = vertices[l].outVertices.size();
                  vertices[l].outVertices.push_back(pos);
                  vertices[l].outReturns.push_back(
                        vertices[pos].outVertices.size() );
                  vertices[l].outCaps.push_back(Infinity);
                  vertices[pos].outVertices.push_back(l);
                  vertices[pos].outReturns.push_back(aux);
                  vertices[pos].outCaps.push_back(0);
                  // if (j == 0) printf( "  (%d,%d)^%d", ii, jj, dd );
               }
         }
      }
   }
}

bool ExtCycleEliminator::generateCut( std::vector<int>& heads,
      std::vector<int>& tails, std::vector<int>& demands,
      std::vector<double>& coeffs )
{
   currentMark = 0;
   // printf( "generating...\n" );

   // send flow from the source to the sink while possible
   double totalFlow = sendFlow();
   double flow = totalFlow;
   // printf( "sent flow of %g\n", flow );
   while ( (flow >= Tolerance) && ((minCut - totalFlow) >= MinViolation) )
   {
      flow = sendFlow();
      totalFlow += flow;
      // printf( "sent flow of %g\n", flow );
   }

   // check if the cut is violated
   if ((minCut - totalFlow) < MinViolation) return false;

   // printf( "total flow of %g out of %g\n", totalFlow, minCut );

   // obtain the entering arcs
   std::vector<int> heads_;
   std::vector<int> tails_;
   std::vector<int> demands_;
   int k, d;
   int maxCap = 0;
   // printf( "num entering = %d out of %d (%d)\n", firstLeaving-2, vertices.size() );
   for (k = 2; k < firstLeaving; k++)
   {
      // printf( "In flow of (%d,%d)^%d is %g out of %g\n",
      //       solution_->arcsCap[vertices[k].index].i,
      //       solution_->arcsCap[vertices[k].index].j,
      //       solution_->arcsCap[vertices[k].index].d,
      //       vertices[0].outCaps[k-2],
      //       solution_->arcsCap[vertices[k].index].value );
      if (vertices[k].mark != currentMark) continue;
      tails_.push_back( solution_->arcsCap[vertices[k].index].i );
      heads_.push_back( solution_->arcsCap[vertices[k].index].j );
      d = solution_->arcsCap[vertices[k].index].d;
      demands_.push_back( d );
      if (d > maxCap) maxCap = d;
   }
   // for (k = firstLeaving; k < vertices.size(); k++)
   // {
      // printf( "Out flow of (%d,%d)^%d is %g out of %g\n",
      //       solution_->arcsCap[vertices[k].index].i,
      //       solution_->arcsCap[vertices[k].index].j,
      //       solution_->arcsCap[vertices[k].index].d,
      //       vertices[k].outCaps[0],
      //       solution_->arcsCap[vertices[k].index].value );
   // }

   // build a matrix of marks for leaving vertices and capacities
   std::vector< std::vector<int> > leaving;
   leaving.resize( setSize_+1 );
   for (k = 1; k <= setSize_; k++)
      leaving[k].resize( maxCap+1, -1 );

   // build a vector of indices of nodes in the set list
   std::vector<int> inSet;
   inSet.resize( n, 0 );
   for (k = 1; k <= setSize_; k++)
      inSet[setList_[k]] = k;

   // build a matrix of marks for entering vertices and capacities
   std::vector< std::vector<int> > entering;
   entering.resize( setSize_+1 );
   for (k = 1; k <= setSize_; k++)
      entering[k].resize( maxCap+1, 0 );

   // set the leaving arcs: for all entering arcs
   int i, j, l, p, w;
   heads.clear();
   tails.clear();
   demands.clear();
   for (k = 0; k < (int) heads_.size(); k++)
   {
      // update the entering arc counter
      entering[inSet[heads_[k]]][demands_[k]]++;

      // for all leaving vertices
      for (l = 1; l <= setSize_; l++)
      {
         i = inSet[heads_[k]];
         j = l;
         if (i > l)
         {
            j = i;
            i = l;
         }

         // for all compatible demands
         for (p = 0; (p < (int) diffCap[i][j].size()) && (p <= demands_[k]); p++)
         {
            if ( !diffCap[i][j][p] ) continue;
            d = demands_[k] - p;

            // all leaving arcs have already been generated
            if ( leaving[l][d] == n ) continue;

            // for all vertices outside the set that do not form cycle
            for (w = 0; w < n; w++)
            {
               // this leaving arc has already been generated
               if ((leaving[l][d] >= 0) && (leaving[l][d] != w)) continue;

               // this leaving arc forms cycle (not with the root
               if ((w == tails_[k]) && (tails_[k] != 0)) continue;

               // this leaving arc goes to the root with non-zero demand
               if ((w == 0) && (d > 0)) continue;

               // this leaving arc has demand zero and does not go to the root
               if ((w > 0) && (d == 0)) continue;

               // this leaving arc goes inside the set
               if (inSet[w] > 0) continue;

               tails.push_back( setList_[l] );
               heads.push_back( w );
               demands.push_back( d );
               coeffs.push_back( -1.0 );
               // printf( " + (%d,%d)^%d", setList_[l], w, d );
            }
            if ( ((leaving[l][d] == -1) || (leaving[l][d] == tails_[k])) &&
                  (tails_[k] != 0) )
               // generated all leaving arcs except to "tails_[k]"
               leaving[l][d] = tails_[k];
            else
               // generated the remaining leaving arc
               leaving[l][d] = n;
         }
      }
   }

   // printf( " >= " );

   // set the lifted entering arcs
   for (k = 0; k < (int) heads_.size(); k++)
   {
      // if more than one entering arc in the same vertex with the same demand
      if ( (entering[inSet[heads_[k]]][demands_[k]] > 1) ||
            ((entering[inSet[heads_[k]]][demands_[k]] > 0) &&
            (tails_[k] == 0)) )
      {
         // set all entering arcs in this vertex with this demand
         for (w = 0; w < n; w++)
         {
            // this leaving arc goes inside the set
            if (inSet[w] > 0) continue;

            heads.push_back( heads_[k] );
            tails.push_back( w );
            demands.push_back( demands_[k] );
            coeffs.push_back( 1.0 );
            // printf( " + (%d,%d)^%d",
            //       w, heads_[k], demands_[k] );
         }
         entering[inSet[heads_[k]]][demands_[k]] = 0;
      }
      else if (entering[inSet[heads_[k]]][demands_[k]] > 0)
      {
         heads.push_back( heads_[k] );
         tails.push_back( tails_[k] );
         demands.push_back( demands_[k] );
         coeffs.push_back( 1.0 );
         // printf( " + (%d,%d)^%d",
         //       tails_[k], heads_[k], demands_[k] );
      }
   }
   // printf( "\n" );
   return true;
}

/*************************************************************************
*                        AUXILIARY PRIVATE METHODS                       *
*************************************************************************/

double ExtCycleEliminator::sendFlow()
{
   int i, j, k;
   double cap;
   std::list<int> activeVertices;
   currentMark++;

   // find a path from the source to the sink
   activeVertices.push_back(0);
   vertices[0].previous = -1;
   vertices[0].maxFlow = Infinity;
   vertices[0].mark = currentMark;
   while ( !activeVertices.empty() )
   {
      // check all adjacent vertices
      i = activeVertices.front();
      activeVertices.pop_front();
      for (k = 0; k < (int) vertices[i].outVertices.size(); k++)
      {
         j = vertices[i].outVertices[k];
         cap = vertices[i].outCaps[k];

         // skip arcs with no capacity
         if (cap < Tolerance) continue;

         // skip the already visited vertices
         if (vertices[j].mark == currentMark) continue;

         // set the previous, the maxFlow and the mark
         vertices[j].previous = i;
         vertices[j].maxFlow = vertices[i].maxFlow;
         if (vertices[j].maxFlow > cap)
            vertices[j].maxFlow = cap;
         vertices[j].mark = currentMark;

         // check if the sink has been reached
         if (j == 1)
         {
            // send flow to the sink and return
            double flow = vertices[j].maxFlow;
            i = vertices[j].previous;
            while (i != -1)
            {
               // update the flows
               for (k = 0; k < (int) vertices[i].outVertices.size(); k++)
                  if (vertices[i].outVertices[k] == j) break;
               vertices[i].outCaps[k] -= flow;
               if (vertices[i].outReturns[k] >= 0)
                  vertices[j].outCaps[vertices[i].outReturns[k]] += flow;

               // got to the previous vertex
               j = i;
               i = vertices[j].previous;
            }
            return flow;
         }

         // add the current vertex to the queue
         activeVertices.push_back(j);
      }
   }
   return 0.0;
}
