/** Source Code File:
 ** Data structure for the degree constraint
 **/

#include "DegreeConstr.hpp"

bool DegreeConstr::checkCapVarCoeff( int i, int j )
{
    return false;
}

double DegreeConstr::getCapVarCoeff( int  i, int j, int d )
{
    return 0.0;
}

double DegreeConstr::getVarCoeff( int i, int j )
{
   if (j == j_) return 1.0;
      else return 0.0;
}

Constraint* DegreeConstr::copy()
{
   return new DegreeConstr( j_ );
}

DegreeConstr::DegreeConstr(int j)
: Constraint('=', 1.0)
{
    j_ = j;
    canBeDeleted_ = false;
}

DegreeConstr::~DegreeConstr()
{
}
