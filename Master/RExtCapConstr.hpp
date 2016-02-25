/** Header File:
 ** Data structure for the rounded extended capacity constraint
 **/

#ifndef _R_EXT_CAP_CONSTR_H_
#define _R_EXT_CAP_CONSTR_H_

#include "Model.hpp"
#include "cutInterface.h"

class RExtCapConstr : public Constraint
{
public:
    bool checkCapVarCoeff( int i, int j );

    double getCapVarCoeff( int  i, int j, int d );

    double getVarCoeff( int i, int j );

    Constraint* copy();

    // Constructor
    RExtCapConstr(ExtCCData* data, double rhs, bool wia = false);

    // Destructor
    ~RExtCapConstr();

private:
    ExtCCData* eccData;

    double trueRhs;

    bool withIntArcs;
};

#endif  // _R_EXT_CAP_CONSTR_H_
