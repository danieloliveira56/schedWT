/** Source Code File:
 ** Data structures to speed up the reduced cost calculation
 **/

#include "Model.hpp"
#include "OrigVariable.hpp"
#include "ConstrShadow.hpp"

//extern int count;

ArcShadow::ArcShadow(ArcVariable* a)
{
    arc = a;
    numRefs = 0;
    //count++;
}

ArcShadow::~ArcShadow()
{
    //count--;
}

CapArcShadow::CapArcShadow(CapArcVariable* ca)
{
    capArc = ca;
    numRefs = 0;
    //count++;
}

CapArcShadow::~CapArcShadow()
{
    //count--;
}

ModCapArcShadow::ModCapArcShadow(ModCapArcVariable* ma)
{
    modArc = ma;
    numRefs = 0;
    //count++;
}

ModCapArcShadow::~ModCapArcShadow()
{
    //count--;
}

CoeffShadow::CoeffShadow(CapArcShadow* cas, ArcShadow* as, double c)
{
    capArc = cas;
    if (cas != 0) cas->addRef();
    arc = as;
    if (as != 0) as->addRef();
    coeff = c;
    next = 0;
    //count++;
}

CoeffShadow::~CoeffShadow()
{
    if (capArc != 0) capArc->deleteRef();
    if (arc != 0) arc->deleteRef();
    //count--;
}

void ConstrShadow::updateRedCosts(double dual)
{
    // traverse the list of constraint coefficients
    CoeffShadow* prev = 0;
    CoeffShadow* aux = coeffList;
    while (aux != 0)
    {
        // remove the coefficients of deleted variables
        if ((aux->getCapArc() == 0) && (aux->getArc() == 0))
        {
            if (prev == 0)
            {
                coeffList = aux->getNext();
                delete aux;
                aux = coeffList;
            }
            else
            {
                prev->setNext(aux->getNext());
                delete aux;
                aux = prev->getNext();
            }
            continue;
        }

        // update the reduced cost associated to the current coefficient
        if (aux->getArc() == 0)
            aux->getCapArc()->addToReducCost(aux->getCoeff()*dual);
        else
            aux->getArc()->addToReducCost(aux->getCoeff()*dual);

        // go to the next
        prev = aux;
        aux = aux->getNext();
    }
}

ConstrShadow::ConstrShadow()
{
    constr = 0;
    coeffList = 0;
    //count++;
}

ConstrShadow::~ConstrShadow()
{
    CoeffShadow* aux = coeffList;
    CoeffShadow* next;
    while (aux != 0)
    {
        next = aux->getNext();
        delete aux;
        aux = next;
    }
    //count--;
}
