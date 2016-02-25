/** Source Code File:
 ** Data structure for the arc constraints with binary coefficients
 **/

#include "BinaryArcConstr.hpp"
#include <assert.h>

bool BinaryArcConstr::checkCapVarCoeff( int i, int j )
{
   return false;
}

double BinaryArcConstr::getCapVarCoeff( int  i, int j, int d )
{
   return 0.0;
}

double BinaryArcConstr::getVarCoeff( int i, int j )
{
   if (arcs[i][j])
      return 1.0;
   else
      return 0.0;
}

Constraint* BinaryArcConstr::copy()
{
   return new BinaryArcConstr( arcs, type_, rhs_ );
}

BinaryArcConstr::BinaryArcConstr(std::vector< std::vector<bool> >& arcs_, char type, double rhs)
: Constraint(type, rhs)
{
    arcs = arcs_;
    canBeDeleted_ = false;
}

BinaryArcConstr::~BinaryArcConstr()
{
}
