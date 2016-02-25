#include "cutInterface.h"
#include "CutGenerator.h"

void* InitCutGenerator( InstanceInfo* instance )
{
   return new CutGenerator( instance );
}

void DestroyCutGenerator( void* cutGenPointer )
{
   CutGenerator* aux = (CutGenerator*)cutGenPointer;
   delete aux;
}

int SeparateECCbyHeuristic( void* cutGen, LpSolution* sol, CutList* cuts, int maxCuts,
      ProblemType prob, int minSetSize )
{
   CutGenerator* cutGenObj = (CutGenerator*) cutGen;
   cutGenObj->setProbType( prob );
   if( ! cutGenObj->isExtCapCutValid() )
      return 1;
   cutGenObj->setLpSolution( sol );
   cutGenObj->setCutBatch( maxCuts );
   cutGenObj->setMinSetSize( minSetSize );
   cutGenObj->extCapCutGenByHeur( cuts );
   return 0;
}

int SeparateExtCyElimByHeur( void* cutGen, LpSolution* sol, CutList* cuts, int maxCuts,
      ProblemType prob )
{
   CutGenerator* cutGenObj = (CutGenerator*) cutGen;
   cutGenObj->setProbType( prob );
   if( ! cutGenObj->isExtCyElimValid() )
      return 1;
   cutGenObj->setLpSolution( sol );
   cutGenObj->setCutBatch( maxCuts );
   cutGenObj->extCyElimCutGenByHeur( cuts );
   return 0;
}

void ClearCutList( CutList* list )
{
   int i;
   for( i = 0; i < list->numCuts; i++ )
   {
      //fprintf(stderr,"Resetting cut %d / %d", i, list->numCuts );
      //fprintf(stderr,"  Data is %s\n", list->cuts[i].data ? " NOT NULL " : " NULL " );
      ResetCut( &( list->cuts[i] ) );
   }
   list->numCuts = 0;
}

void ResetCut( Cut* cut )
{
   if( cut->data )
      cut->DestroyData( cut->data );
   cut->data = 0;

   cut->rhs = 0.0;
   cut->sense = '?';
   strcpy( cut->templabel, "" );
}



