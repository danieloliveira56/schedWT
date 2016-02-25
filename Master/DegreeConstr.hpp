/** Header File:
 ** Data structure for the degree constraint
 **/

#ifndef _DEGREE_CONSTR_H_
#define _DEGREE_CONSTR_H_

#include "Model.hpp"

class DegreeConstr : public Constraint
{
public:
    bool checkCapVarCoeff( int i, int j );

    double getCapVarCoeff( int  i, int j, int d );

    double getVarCoeff( int i, int j );

    Constraint* copy();

    // Constructor
    DegreeConstr(int j);

    // Destructor
    ~DegreeConstr();

private:
    int j_;
};

#endif  // _DEGREE_CONSTR_H_
