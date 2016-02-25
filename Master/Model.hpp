/** Header File:
 ** Data structures for LP model used by the master problem
 **/

#ifndef _MODEL_H_
#define _MODEL_H_

class Constraint
{
public:
    // Getters/Setters
    inline char getType() { return type_; }
    inline double getRhs() { return rhs_; }
    inline bool canBeDeleted() { return canBeDeleted_; }
    inline void canBeDeleted(bool value) { canBeDeleted_ = value; }

    // Return true if "getCapVarCoeff" will return non-zero for at least one
    // demand and false otherwise
    virtual bool checkCapVarCoeff( int i, int j ) = 0;

    // Return the coefficient of the original variable associated with the arc
    // (i,j) with capacity d
    virtual double getCapVarCoeff( int  i, int j, int d ) = 0;

    // Return the coefficient of the original variable associated with the arc
    // (i,j) with no capacity
    virtual double getVarCoeff( int i, int j ) = 0;

    // Make a copy of this constraint
    virtual Constraint* copy() = 0;

    // Destructor
    virtual ~Constraint();

    // Constructor
    Constraint(char type, double rhs);

protected:
    // type: '<', '=' or '>'
    char type_;

    // right-hand side
    double rhs_;

    // indicate the this constraint can be deleted
    bool canBeDeleted_;
};

class Model
{
public:
    // Constructor
    // n is the number of nodes including the root
    Model(int n);

    // Copy constructor
    Model(Model& other);

    // Destructor
    ~Model();

    // Resize the constraint array
    void resizeConstraintArray(int newSize);

    // Getters/Setters
    inline Constraint* getConstraint(int pos) {  return constraintArray[pos]; }
    inline int getNumConstraints() { return numConstraints; }
    inline void setConstraint(int pos, Constraint* constr)
    {  constraintArray[pos] = constr; }

private:
    // Array of pointers to Constraint objects
    Constraint** constraintArray;

    // number of constraints
    int numConstraints;
};

#endif // _MODEL_H_
