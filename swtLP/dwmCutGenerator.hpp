#ifndef DWMCUTGENERATOR_HPP_INCLUDED
#define DWMCUTGENERATOR_HPP_INCLUDED

class Instance;
class DWMaster;
#include "cutInterface.h"

class DWMCutGenerator
{
public:
   DWMCutGenerator( Instance *_inst, DWMaster *_dwm );
   virtual ~DWMCutGenerator();

   int generateCuts();
private:
   void allocCutGenerator();

   void getSolFromCurrentLP();

   void *cutGenerator;
   CutList cutList;

   Instance *inst;
   DWMaster *dwm;

   LpSolution sol;
   /**
    * returns the index-th arc variable in sol
    * (allocates space if necessary)
    **/
   CUTS_ArcVariable *getSolArcVar( const int index );
   /**
    * returns the index-th capacitated arc variable in sol
    * (allocates space if necessary)
    **/
   CUTS_ArcCapVariable *getSolCapacitatedArcVar( const int index );

   int solCapArcs;     // capacity for arcs in sol
   int solCapCapArcs;  // capacity for capacitated arcs
};

#endif // DWMCUTGENERATOR_HPP_INCLUDED
