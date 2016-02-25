#include <cstdio>
#include <cmath>
#include "../Instance.hpp"
#include "dwMaster.hpp"
#include "MasterSolver.hpp"
#include "../CPUTimer.h"
#include "BBNode.hpp"

#include <sstream>
#include "latex_daniel.h"
#define DANIEL_XL_OUTPUT

// node count used when the B&B is interrupted
int runTimeNodeCount = 0;

int main( int argc, char **argv )
{
   // print the all command line
   fprintf( stderr, "Command line:" );
   for (int i = 0; i < argc; i++)
      fprintf( stderr, " %s", argv[i] );
   fprintf( stderr, "\n\n" );

   if ( argc<3 )
   {
      fprintf( stderr, "usage:\n" );
      fprintf( stderr, "\tswtLP instanceFile upperBound\n\n" );
      exit( EXIT_FAILURE );
   }

   Instance inst( argv[1] );

   int ub = atoi( argv[2] );

   CPUTimer t;
   t.start();
   double endvoltime = 0.0;

   const int nDualVars = inst.jobs();
   double *dualsVol = new double[ nDualVars ];
   bool hasDuals = false;
   PricingSolver* pricing = 0;
   double lb = 0;

   // The Volume Algorithm provides dual values and a lower bound
   // for one machine
   if (inst.machines()==1)
   {
      MasterSolver master( &inst, ub );
      hasDuals = true;

      // get the dual variable values from a file if it exists
      char fDualsName[256];
      strcpy( fDualsName, inst.getName() );
      strcat( fDualsName, "_duals.txt" );
      FILE *fDuals = fopen( fDualsName, "r" );
      if ( !fDuals )
      {
         // solve using the fast column generation
         master.solve( ub, ub * 0.001 ); // stop with 0.1% of gap
         t.stop();
         endvoltime = t.getCPUTotalSecs();
         t.start();

         // save a copy of the dual variables
         memcpy( dualsVol, master.getDuals(), sizeof(double)*inst.jobs() );
         fDuals = fopen( fDualsName, "w" );
         if (!fDuals)
         {
            fprintf( stderr, "Could not open a file !.\n" );
            exit( EXIT_FAILURE );
         }
         for ( int i=0 ; (i<inst.jobs()) ; ++i )
            fprintf( fDuals, "%g\n", dualsVol[i] );
         fprintf( fDuals, "\n");
         fclose(fDuals);

         // continue solving using the more expensive subproblem
         if (master.getLowerBound() < ub - 1.0 + IntegerTolerance)
         {
            master.changePricing( ub );
            master.solve( ub );
         }
      }
      else
      {
         // read the dual variable values
         float f;
         fprintf( stderr,"Getting dual values from cache.\n");
         for ( int i=0 ; (i<inst.jobs()) ; i++ )
         {
            if ( fscanf( fDuals, "%g", &f ) != 1 )
            {
               fprintf( stderr, "missing data on file of dual values.\n");
               exit(EXIT_FAILURE);
            }
            dualsVol[i] = f;
         }
         fclose(fDuals);
         fflush(stdout); fflush(stderr);

         // continue solving using the more expensive subproblem
         master.setDuals(dualsVol);
         master.changePricing( ub );
         master.solve( ub );
      }
      pricing = master.getPricingSolver();
      lb = master.getLowerBound() - 1e-3;

      // save a copy of the dual variables
      memcpy( dualsVol, master.getDuals(), sizeof(double)*inst.jobs() );
   }
   else
   {
      delete [] dualsVol;
      dualsVol = 0;
   }

   t.stop();
   double endvolarc = t.getCPUTotalSecs();
   t.start();

#ifdef DEBUG
   //master.checkSolution();
#endif

   // Solve by Branch-and-Bound if necessary
   BBNode* root = 0;
   int cutRounds = 0;
   int numIter = 0;
   double endfirstLP = endvolarc;
   double endcuts = endvolarc;
   double firstLB = 0.0;
   int firstArcs = 0;
   int remainArcs = 0;
   int missPricings = 0;
   int stabChanges = 0;
   bool stopTime = false;
   if (lb < ub - 1.0 + IntegerTolerance)
   {
      fprintf( stderr,"Creating initial linear program and pricing structures...\n");
      root = new BBNode( &inst, ub, pricing, dualsVol );
      fprintf( stderr,"...done.\n"); fflush(stdout);

      // solve the problem but Branch-and-Bound
      try
      {
         root->solve(t, endfirstLP, endcuts, cutRounds);
      }
      catch( char except )
      {
         if (except == 'T')
            stopTime = true;
         else
         {
            fprintf( stderr, "ERROR: Unknown exception!\n" );
            return (-734837);
         }
      }

      // get the statistical data for report
      ub = root->getUB()+1e-5;
      lb = root->getLB() - 1e-3;
      if (lb > ub) lb = ub;
      firstLB = root->getFirstLB() - 1e-3;
      if (firstLB > ub) firstLB = ub;
      numIter = root->getFirstNumIter();
      firstArcs = root->getFirstArcs();
      remainArcs = root->getRemainArcs();
      missPricings = root->getMissPricings();
      stabChanges = root->getStabChanges();
   }

   // update the solution time and get the node count
   int nodeCount = 1;
   if (root != 0)
   {
      nodeCount = root->getNodeCount();
      delete root;
   }
   t.stop();
   double endBBound = t.getCPUTotalSecs();

   if (dualsVol != 0) delete[] dualsVol;

   // adjust the arc counts
   if (int(ceil(firstLB)) >= int(ub)) firstArcs = 0;
   if (int(ceil(lb)) >= int(ub)) remainArcs = 0;

   //inserido por Daniel em 06-06-2015 - imprime a solução inteira encontrada
   //root->write_sol("int.txt");

   // print the statistics
   fprintf( stderr,"Successfully executed.\n");
   fprintf( stderr,"Time-indexed volume time     = %7.1lf secs\n",
         endvoltime );
   fprintf( stderr,"Arc-time-indexed volume time = %7.1lf secs\n",
         endvolarc - endvoltime );
   fprintf( stderr,"First LP lower bound         = %d\n", int(ceil(firstLB)) );
   fprintf( stderr,"First LP gap                 = %6.3lf%c\n",
            100.0 * double(int(ub) - int(ceil(firstLB)))/double(int(ub)), '%' );
   fprintf( stderr,"First LP iterations          = %d\n", numIter );
   fprintf( stderr,"First LP misprices           = %d\n", missPricings );
   fprintf( stderr,"First LP st. center changes  = %d\n", stabChanges );
   fprintf( stderr,"First LP time                = %7.1lf secs\n",
         endfirstLP - endvolarc );
   fprintf( stderr,"Remaining arcs after 1st LP  = %d\n", firstArcs );
   fprintf( stderr,"Root lower bound             = %d\n", int(ceil(lb)) );
   fprintf( stderr,"Root gap                     = %6.3lf%c\n",
            100.0 * double(int(ub) - int(ceil(lb)))/double(int(ub)), '%' );
   fprintf( stderr,"Root cut rounds              = %d\n", cutRounds );
   fprintf( stderr,"Root remaining time          = %7.1lf secs\n",
         endcuts - endfirstLP );
   fprintf( stderr,"Remaining arcs after root    = %d\n", remainArcs );
   if (!stopTime)
   {
      fprintf( stderr,"Branch-and-Bound nodes       = %d\n", nodeCount );
      fprintf( stderr,"Branch-and-Bound time        = %7.1lf secs\n",
            endBBound - endcuts );
      fprintf( stderr,"Overall time                 = %7.1lf secs\n", endBBound );
      fprintf( stderr,"Optimal value                = %d\n", int(ub) );
   }
   else
   {
      fprintf( stderr,"Branch-and-Bound nodes       = >%d\n", runTimeNodeCount );
      fprintf( stderr,"Branch-and-Bound time        = >%7.1lf secs\n",
            endBBound - endcuts );
      fprintf( stderr,"Overall time                 = >%7.1lf secs\n", endBBound );
      fprintf( stderr,"Best upper bound             = %d\n", int(ub) );
   }

#define WIN_PAUSE
#ifdef WIN32 && WIN_PAUSE
   system("pause");
#endif

#ifdef DANIEL_XL_OUTPUT
    LatexTable outputTable;
   
   outputTable.setCaption("Resumo - Teste Root dd/mm/aaaa");
   outputTable.setLabel("teste_dd_mm_resumo");

   outputTable.appendHeader("Inst"); //1
   outputTable.appendHeader("Heu UB"); //2
   outputTable.appendHeader("First LB"); //3
   outputTable.appendHeader("Iter"); //4
   outputTable.appendHeader("Time"); //5
   outputTable.appendHeader("R. Arcs"); //6
   outputTable.appendHeader("Root LB"); //7
   outputTable.appendHeader("CutRounds"); //8
   outputTable.appendHeader("Time"); //9
   outputTable.appendHeader("R. Arcs"); //10   
   outputTable.appendHeader("BCP Nodes"); //11   
   outputTable.appendHeader("BCP Time"); //12   
   outputTable.appendHeader("Opt"); //13   
   outputTable.appendHeader("Status"); //14  

   outputTable.newLine();

   outputTable.appendField(inst.getName()); //1
   outputTable.appendField(ub); //2
   outputTable.appendField(firstLB);  //3
   outputTable.appendField(numIter); //4
   outputTable.appendField(endfirstLP - endvolarc); //5   
   outputTable.appendField(firstArcs);  //6
   outputTable.appendField(lb);  //7
   outputTable.appendField(cutRounds); //8
   outputTable.appendField(endcuts - endfirstLP); //9
   outputTable.appendField(remainArcs); //10
   outputTable.appendField(runTimeNodeCount); //11
   outputTable.appendField(endBBound - endcuts); //12
   outputTable.appendField(ub); //13
	if (!stopTime)
	{
		outputTable.appendField("Optimal"); //14
	}
	else
	{
		outputTable.appendField("Stopped"); //14
	}


   bool printHeader = false;
   bool printFormat = false;
   
   char fname[100];

   //imprime excel
   printHeader = false;
   sprintf(fname, "xl_tbl_output_final.txt");
   FILE *f2  = fopen(fname, "rt");
   if (f2 == NULL)
      {
      printHeader = true;
   } else
   {
      fclose(f2);
   }
   
   f2 = fopen(fname, "at");
   fprintf(f2, "%s", outputTable.printExcel(printHeader).c_str() );
   fclose(f2);
#endif

#ifdef LATEX_OUTPUT
   // printf the table line for the paper
   if (!stopTime)
      fprintf( stderr,"IPCO: %s & %d & %d & %d & %d & %6.1lf & %d & %d & %d &"
            " %6.1lf & %d & %d & %6.1lf & %d\n",
            inst.getName(),
            int(ceil(firstLB)), numIter, missPricings, stabChanges,
            endfirstLP - endvolarc, firstArcs, int(ceil(lb)), cutRounds,
            endcuts - endfirstLP, remainArcs, nodeCount, endBBound, int(ub) );
   else
      fprintf( stderr,"IPCO: %s & %d & %d & %d & %d & %6.1lf & %d & %d & %d &"
            " %6.1lf & %d & %d & >%6.1lf & %d\n",
            inst.getName(),
            int(ceil(firstLB)), numIter, missPricings, stabChanges,
            endfirstLP - endvolarc, firstArcs, int(ceil(lb)), cutRounds,
            endcuts - endfirstLP, remainArcs, nodeCount, endBBound, int(ub) );
#endif

   exit( EXIT_SUCCESS );
}