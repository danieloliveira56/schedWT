/** Source Code File:
 ** Data structures for LP model used by the master problem
 **/

#include "Model.hpp"
#include "DegreeConstr.hpp"

Constraint::Constraint(char type, double rhs)
{
    type_ = type;
    rhs_ = rhs;
}

Constraint::~Constraint()
{
}

Model::Model(int n)
{
    // Create an LP model with all degree constraints
    constraintArray = new Constraint*[n-1];
    for (int i = 0; i < n-1; i++)
        constraintArray[i] = new DegreeConstr(i+1);
    numConstraints = n-1;
}

Model::Model( Model& other ) :
      constraintArray( new Constraint*[other.numConstraints] )
{
    for (int i = 0; i < other.numConstraints; i++)
        constraintArray[i] = other.constraintArray[i]->copy();
   numConstraints = other.numConstraints;
}

Model::~Model()
{
    // release the array of constraints
    for (int i = 0; i < numConstraints; i++)
        if (constraintArray[i] != 0)
            delete constraintArray[i];
    delete [] constraintArray;
}

void Model::resizeConstraintArray(int newSize)
{
    // only change the size if the new size is not greater
    if (newSize <= numConstraints)
    {
        numConstraints = newSize;
        return;
    }

    // allocate the new Array and copy the constraints
    Constraint** newArray = new Constraint*[newSize];
    for (int i = 0; i < numConstraints; i++)
        newArray[i] = constraintArray[i];
    for (int i = numConstraints; i < newSize; i++)
        newArray[i] = 0;

    // replace the old array
    delete [] constraintArray;
    constraintArray = newArray;
    numConstraints = newSize;
}
