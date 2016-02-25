/** Source Code File:
 ** Data structure for the generic constraint
 **/

#include "GenericConstr.hpp"
#include <assert.h>

bool GenericConstr::checkCapVarCoeff( int i, int j )
{
   ArcCapCoeffHash::iterator it;
   ArcCapHashKey k;

   k.i = i;
   k.j = j;
   k.d = AnyDemand;
   it = cutData->coeffs.find( k );
   if( it == cutData->coeffs.end() )
      return false;
   else
      return true;
}

double GenericConstr::getCapVarCoeff( int  i, int j, int d )
{
   ArcCapCoeffHash::iterator it;
   ArcCapHashKey k;

   k.i = i;
   k.j = j;
   k.d = d;
   it = cutData->coeffs.find( k );
   if( it == cutData->coeffs.end() )
      return 0.0;
   else
      return (*it).second;
}

double GenericConstr::getVarCoeff( int i, int j )
{
   return 0.0;
}

Constraint* GenericConstr::copy()
{
   GenericCutData* data = new GenericCutData;
   *data = *cutData;
   return new GenericConstr( data, type_, rhs_ );
}

GenericConstr::GenericConstr(GenericCutData* data, char type, double rhs)
: Constraint(type, rhs)
{
    cutData = data;
    canBeDeleted_ = true;
}

GenericConstr::~GenericConstr()
{
   delete cutData;
}
