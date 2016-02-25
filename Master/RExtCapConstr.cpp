/** Source Code File:
 ** Data structure for the rounded extended capacity constraint
 **/

#include "RExtCapConstr.hpp"
#include <assert.h>

bool RExtCapConstr::checkCapVarCoeff( int i, int j )
{
//   fprintf(stderr," i = %d, j = %d, S.size() = %d\n", i, j, eccData->S.size() );
   assert( eccData != NULL );
   assert( i < (int)eccData->S.size() );
   assert( j < (int)eccData->S.size() );
   assert( i >= 0 );
   assert( j >= 0 );

   if (!withIntArcs)
   {
      if ( eccData->S[j] != eccData->S[i] )
         return 1.0;
      else
         return 0.0;
   }
   else
   {
      if ( eccData->S[j] || eccData->S[i] )
         return 1.0;
      else
         return 0.0;
   }
}

double RExtCapConstr::getCapVarCoeff( int  i, int j, int d )
{
//   fprintf(stderr," i = %d, j = %d, S.size() = %d\n", i, j, eccData->S.size() );
   assert( eccData != NULL );
   assert( i < (int)eccData->S.size() );
   assert( j < (int)eccData->S.size() );
   assert( i >= 0 );
   assert( j >= 0 );

   if (d > eccData->capacity)
      return 0.0;
   if (!withIntArcs)
   {
      if ( eccData->S[j] && !eccData->S[i] )
         // return (eccData->incapcoeff[d] / trueRhs);
         return eccData->incapcoeff[d];
      else if( !eccData->S[j] && eccData->S[i] )
         // return (eccData->outcapcoeff[d] / trueRhs);
         return eccData->outcapcoeff[d];
      else
         return 0.0;
   }
   else
   {
      double coeff = 0.0;
      if ( eccData->S[j] )
         // return (eccData->incapcoeff[d] / trueRhs);
         coeff += eccData->incapcoeff[d];
      else if( eccData->S[i] )
         // return (eccData->outcapcoeff[d] / trueRhs);
         coeff += eccData->outcapcoeff[d];
      return coeff;
   }
}

double RExtCapConstr::getVarCoeff( int i, int j )
{
   return 0.0;
}

Constraint* RExtCapConstr::copy()
{
   ExtCCData* data = new ExtCCData;
   *data = *eccData;
   return new RExtCapConstr( data, trueRhs, withIntArcs );
}

RExtCapConstr::RExtCapConstr(ExtCCData* data, double rhs, bool wia)
//: Constraint('>', 1.0)
: Constraint('>', rhs), withIntArcs(wia)
{
    eccData = data;
    trueRhs = rhs;
    canBeDeleted_ = true;
}

RExtCapConstr::~RExtCapConstr()
{
   delete eccData;
}
