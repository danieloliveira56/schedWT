/** Header File:
 ** Data structures for the original problem variables
 **/

#ifndef _ORIG_VARIABLE_H_
#define _ORIG_VARIABLE_H_

#include "ConstrShadow.hpp"

class ModCapArcVariable;

// a partial subproblem solution descriptor
struct SolutionDescr
{
    ModCapArcVariable* lastArc;
    int nextPos;
    long double cost;
    double origCost;
};

class ModCapArcNode;

class CapArcVariable;

class ModCapArcVariable
{
   public:
      inline int getCap() {
         return cap;
      }
      inline int getMod() {
         return mod;
      }
      inline ModCapArcNode* getHead() {
         return head;
      }
      inline ModCapArcNode* getTail() {
         return tail;
      }
      inline void setPosHead( int pos ) {
         posHead = pos;
      }
      inline void setPosTail( int pos ) {
         posTail = pos;
      }
      inline CapArcVariable* getCapArc() {
         return capArc;
      }
      inline void setPosCapArc( int pos ) {
         posCapArc = pos;
      }
      inline ModCapArcShadow* getShadow() {
         return shadow;
      }

      // Constructor/Destructor
      ModCapArcVariable( int c, int m, ModCapArcNode* h, ModCapArcNode* t,
            CapArcVariable* ca );
      ~ModCapArcVariable();

      // number of objects which were
      // created and were not freed
      static int activeObjects;
   private:
      int cap;
      int mod;
      ModCapArcNode* head;
      ModCapArcNode* tail;
      int posHead;
      int posTail;
      CapArcVariable* capArc;
      int posCapArc;
      ModCapArcShadow* shadow;

};

class ModCapArcNode
{
   public:
      inline int getCap() {
         return cap;
      }
      inline int getMod() {
         return mod;
      }
      inline int getIndex() {
         return index;
      }
      inline int getInDegree() {
         return inDegree;
      }
      inline int getOutDegree() {
         return outDegree;
      }
      inline ModCapArcVariable* getIn( int pos ) {
         return in[pos];
      }
      inline ModCapArcVariable* getOut( int pos ) {
         return out[pos];
      }
      inline void addIn( ModCapArcVariable* modCapArc )
      {
         in[inDegree] = modCapArc;
         modCapArc->setPosHead( inDegree );
         inDegree++;
      }
      inline void addOut( ModCapArcVariable* modCapArc )
      {
         out[outDegree] = modCapArc;
         modCapArc->setPosTail( outDegree );
         outDegree++;
      }
      inline void removeIn( int pos )
      {
         if ( pos != inDegree-1 )
         {
            in[pos] = in[inDegree-1];
            in[pos]->setPosHead( pos );
         }
         inDegree--;
      }
      inline void removeOut( int pos )
      {
         if ( pos != outDegree-1 )
         {
            out[pos] = out[outDegree-1];
            out[pos]->setPosTail( pos );
         }
         outDegree--;
      }
      inline bool isDeleted() {
         return deleted;
      }
      inline void setDeleted( bool del ) {
         deleted = del;
      }

      // Getters/Setters for the solution pointers
      inline int getNumSolIn() {
         return numSolIn;
      }
      inline SolutionDescr getSolIn( int pos ) {
         return solIn[pos];
      }
      inline void setSolIn( int pos, SolutionDescr& sol ) {
         solIn[pos] = sol;
      }
      void setNumSolIn( int num );  // also reallocate the array "solIn"
      inline int getNumSolOut() {
         return numSolOut;
      }
      inline SolutionDescr getSolOut( int pos ) {
         return solOut[pos];
      }
      inline void setSolOut( int pos, SolutionDescr& sol ) {
         solOut[pos] = sol;
      }
      void setNumSolOut( int num );  // also reallocate the array "solOut"

      // auxiliary getters and setters
      inline ModCapArcNode* getAuxNext() {
         return auxNext;
      }
      inline void setAuxNext( ModCapArcNode* a ) {
         auxNext = a;
      }
      inline int getAuxCount() {
         return auxCount;
      }
      inline void setAuxCount( int c ) {
         auxCount = c;
      }

      // Constructor/Destructor
      ModCapArcNode( int c, int m, int i, int maxIn, int maxOut );
      ~ModCapArcNode();

      static int activeObjects;
   private:
      int cap;
      int mod;
      int index;
      int inDegree;
      ModCapArcVariable** in;
      int outDegree;
      ModCapArcVariable** out;

      // when true, executing the destructor
      bool deleted;

      // current best solutions using input and output arcs
      int numSolIn;
      SolutionDescr* solIn;
      int numSolOut;
      SolutionDescr* solOut;

      // auxiliary pointer and counter used to search in the graph
      ModCapArcNode* auxNext;
      int auxCount;

};

class CapArcNode;

class ArcVariable;

class CapArcVariable
{
   public:
      inline int getCap() {
         return cap;
      }
      inline CapArcNode* getHead() {
         return head;
      }
      inline CapArcNode* getTail() {
         return tail;
      }
      inline void setPosHead( int pos ) {
         posHead = pos;
      }
      inline void setPosTail( int pos ) {
         posTail = pos;
      }
      inline ArcVariable* getArc() {
         return arc;
      }
      inline void setPosArc( int pos ) {
         posArc = pos;
      }
      inline void resetReducCost() {
         reducCost = origCost;
      }
      inline void addToReducCost( double cost ) {
         reducCost += (long double)cost;
      }
      inline long double getReducCost() {
         return reducCost;
      }
      inline double getOrigCost() {
         return origCost;
      }
      inline ModCapArcVariable* getMod( int pos ) {
         return mods[pos];
      }
      inline int getNumMods() {
         return numMods;
      }
      inline void addMod( ModCapArcVariable* modCapArc )
      {
         mods[numMods] = modCapArc;
         modCapArc->setPosCapArc( numMods );
         numMods++;
      }
      inline void removeMod( int pos )
      {
         if ( pos != numMods-1 )
         {
            mods[pos] = mods[numMods-1];
            mods[pos]->setPosCapArc( pos );
         }
         numMods--;
      }
      inline CapArcShadow* getShadow() {
         return shadow;
      }
      inline bool isDeleted() {
         return deleted;
      }

      // Constructor/Destructor
      CapArcVariable( int c, CapArcNode* h, CapArcNode* t, ArcVariable* a,
            double oc, int maxMods );
      ~CapArcVariable();

      static int activeObjects;
   private:
      int cap;
      CapArcNode* head;
      CapArcNode* tail;
      int posHead;
      int posTail;
      ArcVariable* arc;
      int posArc;
      double origCost;
      long double reducCost;
      int numMods;
      ModCapArcVariable** mods;
      CapArcShadow* shadow;

      // when true, executing the destructor
      bool deleted;

};


class CapArcNode
{
   public:
      inline int getCap() {
         return cap;
      }
      inline int getIndex() {
         return index;
      }
      inline int getInDegree() {
         return inDegree;
      }
      inline int getOutDegree() {
         return outDegree;
      }
      inline CapArcVariable* getIn( int pos ) {
         return in[pos];
      }
      inline CapArcVariable* getOut( int pos ) {
         return out[pos];
      }
      inline void addIn( CapArcVariable* capArc )
      {
         in[inDegree] = capArc;
         capArc->setPosHead( inDegree );
         inDegree++;
      }
      inline void addOut( CapArcVariable* capArc )
      {
         out[outDegree] = capArc;
         capArc->setPosTail( outDegree );
         outDegree++;
      }
      inline void removeIn( int pos )
      {
         if ( pos != inDegree-1 )
         {
            in[pos] = in[inDegree-1];
            in[pos]->setPosHead( pos );
         }
         inDegree--;
      }
      inline void removeOut( int pos )
      {
         if ( pos != outDegree-1 )
         {
            out[pos] = out[outDegree-1];
            out[pos]->setPosTail( pos );
         }
         outDegree--;
      }
      inline bool isDeleted() {
         return deleted;
      }

      // Constructor/Destructor
      CapArcNode( int c, int i, int maxIn, int maxOut );
      ~CapArcNode();

      static int activeObjects;
   private:
      int cap;
      int index;
      int inDegree;
      CapArcVariable** in;
      int outDegree;
      CapArcVariable** out;

      // when true, executing the destructor
      bool deleted;

};

class ArcNode;

class ArcVariable
{
   public:
      inline ArcNode* getHead() {
         return head;
      }
      inline ArcNode* getTail() {
         return tail;
      }
      inline void setPosHead( int pos ) {
         posHead = pos;
      }
      inline void setPosTail( int pos ) {
         posTail = pos;
      }
      inline void resetReducCost() {
         reducCost = 0.0;
      }
      inline void addToReducCost( double cost ) {
         reducCost += (long double)cost;
      }
      inline long double getReducCost() {
         return reducCost;
      }
      inline CapArcVariable* getCap( int pos ) {
         return caps[pos];
      }
      inline int getNumCaps() {
         return numCaps;
      }
      inline void addCap( CapArcVariable* capArc )
      {
         caps[numCaps] = capArc;
         capArc->setPosArc( numCaps );
         numCaps++;
      }
      inline void removeCap( int pos )
      {
         if ( pos != numCaps-1 )
         {
            caps[pos] = caps[numCaps-1];
            caps[pos]->setPosArc( pos );
         }
         numCaps--;
      }
      inline ArcShadow* getShadow() {
         return shadow;
      }
      inline bool isDeleted() {
         return deleted;
      }

      // Constructor/Destructor
      ArcVariable( ArcNode* h, ArcNode* t, int maxCaps );
      ~ArcVariable();

      static int activeObjects;
   private:
      ArcNode* head;
      ArcNode* tail;
      int posHead;
      int posTail;
      int numCaps;
      long double reducCost;
      CapArcVariable** caps;
      ArcShadow* shadow;

      // when true, executing the destructor
      bool deleted;

};


class ArcNode
{
   public:
      inline int getIndex() {
         return index;
      }
      inline int getInDegree() {
         return inDegree;
      }
      inline int getOutDegree() {
         return outDegree;
      }
      inline ArcVariable* getIn( int pos ) {
         return in[pos];
      }
      inline ArcVariable* getOut( int pos ) {
         return out[pos];
      }
      inline void addIn( ArcVariable* arc )
      {
         in[inDegree] = arc;
         arc->setPosHead( inDegree );
         inDegree++;
      }
      inline void addOut( ArcVariable* arc )
      {
         out[outDegree] = arc;
         arc->setPosTail( outDegree );
         outDegree++;
      }
      inline void removeIn( int pos )
      {
         if ( pos != inDegree-1 )
         {
            in[pos] = in[inDegree-1];
            in[pos]->setPosHead( pos );
         }
         inDegree--;
      }
      inline void removeOut( int pos )
      {
         if ( pos != outDegree-1 )
         {
            out[pos] = out[outDegree-1];
            out[pos]->setPosTail( pos );
         }
         outDegree--;
      }
      inline bool isDeleted() {
         return deleted;
      }

      // auxiliary getters and setters
      inline ArcNode* getAuxNext() {
         return auxNext;
      }
      inline void setAuxNext( ArcNode* a ) {
         auxNext = a;
      }
      inline int getAuxCount() {
         return auxCount;
      }
      inline void setAuxCount( int c ) {
         auxCount = c;
      }

      // Constructor/Destructor
      ArcNode( int i, int maxIn, int maxOut );
      ~ArcNode();

      static int activeObjects;
   private:
      int index;
      int inDegree;
      ArcVariable** in;
      int outDegree;
      ArcVariable** out;

      // when true, executing the destructor
      bool deleted;

      // auxiliary pointer and counter used to search in the graph
      ArcNode* auxNext;
      int auxCount;

};

void printLeaks();

#endif  // _ORIG_VARIABLE_H_
