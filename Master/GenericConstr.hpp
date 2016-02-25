/** Header File:
 ** Data structure for the generic constraint
 **/

#ifndef _GENERIC_CONSTR_H_
#define _GENERIC_CONSTR_H_

#include "Model.hpp"
#include "cutInterface.h"

class GenericConstr : public Constraint
{
public:
    bool checkCapVarCoeff( int i, int j );

    double getCapVarCoeff( int  i, int j, int d );

    double getVarCoeff( int i, int j );

    Constraint* copy();

    // Constructor
    GenericConstr(GenericCutData* data, char type, double rhs);

    // Destructor
    ~GenericConstr();

private:
    GenericCutData* cutData;
};

#endif  // _GENERIC_CONSTR_H_
