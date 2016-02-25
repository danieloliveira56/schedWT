/** Header File:
 ** Data structures to speed up the reduced cost calculation
 **/

#ifndef _CONSTR_SHADOW_H_
#define _CONSTR_SHADOW_H_

class ArcVariable;

// this class refers to an ArcVarible object but is kept alive
// after the deletion of the referred object to indicate this deletion
// to other data structures that refer to it.
class ArcShadow
{
public:
    inline ArcVariable* getArc() { return arc; }
    inline void deleteArc()
    {
        arc = 0;
        if (numRefs == 0) delete this;
    }
    inline void addRef() { numRefs++; }
    inline void deleteRef()
    {
        numRefs--;
        if ((arc == 0) && (numRefs == 0)) delete this;
    }

    // Constructor
    ArcShadow(ArcVariable* a);
    ~ArcShadow();
private:
    ArcVariable* arc;
    int numRefs;
};

class CapArcVariable;

// this class refers to a CapArcVarible object but is kept alive
// after the deletion of the referred object to indicate this deletion
// to other data structures that refer to it.
class CapArcShadow
{
public:
    inline CapArcVariable* getCapArc() { return capArc; }
    inline void deleteCapArc()
    {
        capArc = 0;
        if (numRefs == 0) delete this;
    }
    inline void addRef() { numRefs++; }
    inline void deleteRef()
    {
        numRefs--;
        if ((capArc == 0) && (numRefs == 0)) delete this;
    }

    // Constructor
    CapArcShadow(CapArcVariable* ca);
    ~CapArcShadow();

private:
    CapArcVariable* capArc;
    int numRefs;
};

class ModCapArcVariable;

// this class refers to a ModCapArcVarible object but is kept alive
// after the deletion of the referred object to indicate this deletion
// to other data structures that refer to it.
class ModCapArcShadow
{
public:
    inline ModCapArcVariable* getModArc() { return modArc; }
    inline void deleteModArc()
    {
        modArc = 0;
        if (numRefs == 0) delete this;
    }
    inline void addRef() { numRefs++; }
    inline void deleteRef()
    {
        numRefs--;
        if ((modArc == 0) && (numRefs == 0)) delete this;
    }

    // Constructor
    ModCapArcShadow(ModCapArcVariable* ma);
    ~ModCapArcShadow();

private:
    ModCapArcVariable* modArc;
    int numRefs;
};

// this class is an element of a list of coefficients in a constraint shadow
class CoeffShadow
{
public:
    inline void setNext(CoeffShadow* n) { next = n; }

    inline CapArcVariable* getCapArc()
    {
        if (capArc == 0) return 0;
        else return capArc->getCapArc();
    }
    inline ArcVariable* getArc()
    {
        if (arc == 0) return 0;
        else return arc->getArc();
    }
    inline double getCoeff() { return coeff; }
    inline CoeffShadow* getNext() { return next; }

    // Constructors/Destructor
    CoeffShadow(CapArcShadow* cas, ArcShadow* as, double c);
    ~CoeffShadow();

private:
    CapArcShadow* capArc;
    ArcShadow* arc;
    double coeff;
    CoeffShadow* next;
};

class Constraint;

// this class stores additional information about a Constraint to speed up the
// reduced cost calculations
class ConstrShadow
{
public:
    inline void setConstr(Constraint* c) { constr = c; }
    inline void addCoeff(CoeffShadow* cs)
    {
        cs->setNext(coeffList);
        coeffList = cs;
    }
    inline void removeCoeff(CoeffShadow* prev)
    {
        if (prev == 0)
            coeffList = coeffList->getNext();
        else
            prev->setNext(prev->getNext()->getNext());
    }
    inline CoeffShadow* getCoeffList() { return coeffList; }
    inline Constraint* getConstr() { return constr; }

    // update the reduced costs (and remove shadows of deleted coefficients)
    void updateRedCosts(double dual);

    // Constructor/Destructor
    ConstrShadow();
    ~ConstrShadow();

private:
    CoeffShadow* coeffList;
    Constraint* constr;
};

#endif // _CONSTR_SHADOW_H_
