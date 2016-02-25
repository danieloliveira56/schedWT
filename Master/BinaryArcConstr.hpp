/** Header File:
 ** Data structure for the arc constraints with binary coefficients
 **/

#ifndef _BINARY_ARC_CONSTR_H_
#define _BINARY_ARC_CONSTR_H_

#include "Model.hpp"
#include "cutInterface.h"

#include <vector>

class BinaryArcConstr : public Constraint
{
public:
    bool checkCapVarCoeff( int i, int j );

    double getCapVarCoeff( int  i, int j, int d );

    double getVarCoeff( int i, int j );

    Constraint* copy();

    // Constructor
    BinaryArcConstr(std::vector< std::vector<bool> >& arcs_, char type, double rhs);

    // Destructor
    ~BinaryArcConstr();

private:
    std::vector< std::vector<bool> > arcs;
};

#endif  // _BINARY_ARC_CONSTR_H_
