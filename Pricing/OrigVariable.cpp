/** Source Code File:
 ** Data structures for the original problem variables
 **/

#include "OrigVariable.hpp"
#include <cstdio>

// number of objects which were
// created and were not freed
int ModCapArcVariable::activeObjects = 0;
int ModCapArcNode::activeObjects     = 0;
int CapArcVariable::activeObjects    = 0;
int CapArcNode::activeObjects        = 0;
int ArcVariable::activeObjects       = 0;
int ArcNode::activeObjects           = 0;

ModCapArcVariable::ModCapArcVariable( int c, int m, ModCapArcNode* h,
      ModCapArcNode* t, CapArcVariable* ca )
{
   cap = c;
   mod = m;
   head = h;
   h->addIn(this);
   tail = t;
   t->addOut(this);
   capArc = ca;
   ca->addMod(this);
   shadow = new ModCapArcShadow(this);

   ModCapArcVariable::activeObjects++;
}

ModCapArcVariable::~ModCapArcVariable()
{
   head->removeIn(posHead);
   tail->removeOut(posTail);
   capArc->removeMod(posCapArc);
   if ((head->getIndex() != 0) && (head->getInDegree() == 0)
         && !head->isDeleted())
      delete head;
   if ((tail->getIndex() != 0) && (tail->getOutDegree() == 0)
         && !tail->isDeleted())
      delete tail;
   if ((capArc->getNumMods() == 0) && !capArc->isDeleted())
      delete capArc;
   shadow->deleteModArc();

   ModCapArcVariable::activeObjects--;
}

ModCapArcNode::ModCapArcNode( int c, int m, int i, int maxIn, int maxOut )
{
   cap = c;
   mod = m;
   index = i;
   inDegree = 0;
   in = new ModCapArcVariable*[maxIn+1];
   outDegree = 0;
   out = new ModCapArcVariable*[maxOut+1];
   deleted = false;
   solIn = new SolutionDescr[1];
   numSolIn = 1;
   solOut = new SolutionDescr[1];
   solOut->cost = 1e9;
   numSolOut = 1;
   auxNext = 0;
   auxCount = 0;

   ModCapArcNode::activeObjects++;
}

ModCapArcNode::~ModCapArcNode()
{
   deleted = true;
   while (inDegree > 0)
      delete in[0];
   delete [] in;
   while (outDegree > 0)
      delete out[0];
   delete [] out;
   delete [] solIn;
   delete [] solOut;

   ModCapArcNode::activeObjects--;
}

void ModCapArcNode::setNumSolIn(int num)
{
   delete [] solIn;
   solIn = new SolutionDescr[num];
   numSolIn = num;
}

void ModCapArcNode::setNumSolOut(int num)
{
   delete [] solOut;
   solOut = new SolutionDescr[num];
   numSolOut = num;
}

CapArcVariable::CapArcVariable( int c, CapArcNode* h, CapArcNode* t,
      ArcVariable* a, double oc, int maxMods )
{
   cap = c;
   head = h;
   h->addIn(this);
   tail = t;
   t->addOut(this);
   arc = a;
   a->addCap(this);
   origCost = oc;
   numMods = 0;
   mods = new ModCapArcVariable*[maxMods];
   shadow = new CapArcShadow(this);
   deleted = false;

   CapArcVariable::activeObjects++;
}

CapArcVariable::~CapArcVariable()
{
   deleted = true;
   head->removeIn(posHead);
   tail->removeOut(posTail);
   arc->removeCap(posArc);
   if ((head->getIndex() != 0) && (head->getInDegree() == 0)
         && (head->getOutDegree() == 0) && !head->isDeleted())
      delete head;
   if ((tail->getIndex() != 0) && (tail->getOutDegree() == 0)
         && (tail->getInDegree() == 0) && !tail->isDeleted())
      delete tail;
   if ((arc->getNumCaps() == 0) && !arc->isDeleted())
      delete arc;
   while (numMods > 0)
      delete mods[0];
   delete [] mods;
   shadow->deleteCapArc();

   CapArcVariable::activeObjects--;
}


CapArcNode::CapArcNode( int c, int i, int maxIn, int maxOut )
{
   cap = c;
   index = i;
   inDegree = 0;
   in = new CapArcVariable*[maxIn+1];
   outDegree = 0;
   out = new CapArcVariable*[maxOut+1];
   deleted = false;

   CapArcNode::activeObjects++;
}

CapArcNode::~CapArcNode()
{
   deleted = true;
   while (inDegree > 0)
      delete in[0];
   delete [] in;
   while (outDegree > 0)
      delete out[0];
   delete [] out;

   CapArcNode::activeObjects--;
}


ArcVariable::ArcVariable( ArcNode* h, ArcNode* t, int maxCaps )
{
   head = h;
   h->addIn(this);
   tail = t;
   t->addOut(this);
   numCaps = 0;
   caps = new CapArcVariable*[maxCaps];
   shadow = new ArcShadow(this);
   deleted = false;

   ArcVariable::activeObjects++;
}

ArcVariable::~ArcVariable()
{
   deleted = true;
   head->removeIn(posHead);
   tail->removeOut(posTail);
   if ((head->getIndex() != 0) && (head->getInDegree() == 0)
         && (head->getOutDegree() == 0) && !head->isDeleted())
      delete head;
   if ((tail->getIndex() != 0) && (tail->getOutDegree() == 0)
         && (tail->getInDegree() == 0) && !tail->isDeleted())
      delete tail;
   while (numCaps > 0)
      delete caps[0];
   delete [] caps;
   shadow->deleteArc();

   ArcVariable::activeObjects--;
}


ArcNode::ArcNode( int i, int maxIn, int maxOut )
{
   index = i;
   inDegree = 0;
   in = new ArcVariable*[maxIn+1];
   outDegree = 0;
   out = new ArcVariable*[maxOut+1];
   deleted = false;
   auxNext = 0;
   auxCount = 0;

   ArcNode::activeObjects++;
}

ArcNode::~ArcNode()
{
   deleted = true;
   while (inDegree > 0)
      delete in[0];
   delete [] in;
   while (outDegree > 0)
      delete out[0];
   delete [] out;

   ArcNode::activeObjects--;
}

void printLeaks()
{
   fprintf( stderr,"Active Objects by class:\n");
   fprintf( stderr,"ModCapArcVariable  : %d\n", ModCapArcVariable::activeObjects );
   fprintf( stderr,"ModCapArcNode      : %d\n", ModCapArcNode::activeObjects     );
   fprintf( stderr,"CapArcVariable     : %d\n", CapArcVariable::activeObjects    );
   fprintf( stderr,"CapArcNode         : %d\n", CapArcNode::activeObjects        );
   fprintf( stderr,"ArcVariable        : %d\n", ArcVariable::activeObjects       );
   fprintf( stderr,"ArcNode            : %d\n", ArcNode::activeObjects           );
}

