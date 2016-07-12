/** Source Code File:
 ** Branch-and-bound node routines
 **/

#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <sstream>

#include "BBNode.hpp"
#include "../Instance.hpp" 
#include "dwMaster.hpp" 
#include "../CPUTimer.h"
#include "Model.hpp"
#include "RExtCapConstr.hpp"
#include "GenericConstr.hpp"
#include "BinaryArcConstr.hpp"
#include "CutGenerator.h"
#include "latex_daniel.h"

#include <map>
#include <algorithm>

using std::sort;

//#define NO_CUT
#define STRONG_BRANCHING
//#define SB_BRANCH_ON_TIME
///#define SB_BRANCH_ON_ACCUM_TIME
//#define SB_BRANCH_ON_SINGLE_TIME
//#define DANIEL_BRANCHING
#define LOG_ON_TREE
//#define ROOT_ONLY

#define ORIG_ECCSEP 
//#define GENETIC_ECC_SEP
//#define GENETIC_ECC_s_SEP

#define TRICLIQUECUT_SEP
//#define GENERIC_CLIQUECUT_SEP

#define INTEGRATED_OEC_SEP
#define OEC_SEP

//#define PRINT_TEX
#define PRINT_XL
//#define PRINT_TEX_GERAL
#define PRINT_XL_GERAL


// cpu time limit for the whole execution in seconds
#ifdef ROOT_ONLY
// time limits at node level (only stops after finishing the current node)
const double CpuTimeTryLimit = 0.0;         // stop if not enough completed
const double CpuTimeLastLimit = 0.0;        // stop anyway

#ifdef NO_CUT
// time limits at pricing level (may stop with an unfinished node)
const double CpuTimeHardLimit = 14400.0;    // stop at any time
#endif
#else
// time limits at node level (only stops after finishing the current node)
const double CpuTimeTryLimit = 14400.0;     // stop if not enough completed
const double CpuTimeLastLimit = 86400.0;    // stop anyway
#endif
double CpuTimeLimit = CpuTimeTryLimit;

// B&B level to consider enough completed after "CpuTimeTryLimit" seconds
const int LevelToContinue = 4;

// maximum number of non-fixed arcs required to generate an LP file for
// post-solving a node
const int MaxArcsLpGeneration = 200000;
//const int MaxArcsLpGeneration = 130400;

// minimum GAP reduction to continue with cuts
//Artur used 0.001
const double MinGapReducCut = 0.001;

// minimum reduction on the number of arcs to continue with cuts
//artur used 1
const double MinArcReducCut = 1;

// maximum instance size for using cuts in non-root branch-and-bound nodes
//Artur used 50, but didn't use the parameter below
const int MaxSizeNonRootCuts = 100;

// maximum instance size for using RHECC cuts in non-root branch-and-bound nodes
//defined by Daniel
const int MaxSizeNonRootRHECCCuts = 50;

// maximum number of cut round in non-root branch-and-bound nodes
//Artur used 1
const int MaxCutRoundsNonRoot = 1; 

// fraction of the upper bound per job added to all arc costs for branching
const double FracCostAddBranch = 1.0;

// number of candidates when doing strong branching
const int NumStrBrCandidates = 16;

// maximum number of column generations when evaluating the strong branching
// after the DWM LP goes bellow the upper bound
const int MaxColGenStrBranch = 20;

// maximum number of separated cuts per cut round
const int MaxSepECECsPerRound = 50;
const int MaxSepECCsPerRoundSingle = 20;
const int MaxSepECCsPerRoundMulti = 50;
int MaxSepCutsPerRound;

// maximum number of cuts inserted in the relaxation per cut round
const int MaxInsECECsPerRound = 50;
const int MaxInsECCsPerRoundSingle = 20;
const int MaxInsECCsPerRoundMulti = 50;
int MaxInsCutsPerRound;

// minimum improvement over the lower bound to allow removing rows from LP
const double MinImproveToCleanLP = 0.1;

// target variable value for the time-based branching rules
const double TargetBrTimeValue = 0.5;

// node count used when the B&B is interrupted
extern int runTimeNodeCount;

// A pointer to the cut generation object
void* cutGenPointer = 0;

//due to the randomness of the OEC separation, a tolerance of bad cut rounds may be of help
int badRoundsLimit = 2;

BBNode::BBNode( Instance *_inst, double _ub, PricingSolver* pricing,
      const double *_duals )
{
   // Create the Dantzig-Wolf Master
   dwm = new DWMaster(_inst, _ub, pricing, _duals);

   // do other initializations
   nodeCount = 1;
   nodeNumber = 0;
   level = 0;
   inst = _inst;
   ub = _ub;
   sol.arcs = 0;
   sol.arcsCap = 0;
   isIntegerSol = false;
   firstLB = 0.0;
   firstNumIter = 0;
   firstArcs = 0;
   missPricings = 0;
   stabChanges = 0;
}

BBNode::BBNode( BBNode& other )
{
   // Copy the Dantzig-Wolf Master
   dwm = new DWMaster(*other.dwm);

   // do other initializations
   nodeCount = 1;
   nodeNumber = -1;
   level = other.level + 1;
   inst = other.inst;
   ub = other.ub;
   sol.arcs = 0;
   sol.arcsCap = 0;
   isIntegerSol = false;
   firstLB = 0.0;
   firstNumIter = 0;
   firstArcs = 0;
   missPricings = 0;
   stabChanges = 0;
}

BBNode::~BBNode()
{
   // Release the Dantzig-Wolf Master
   delete dwm;

   // release the previous expanded solution is any
   if (sol.arcs != 0) delete [] sol.arcs;
   if (sol.arcsCap != 0) delete [] sol.arcsCap;
}

void BBNode::write_sol(char fname[]) {
   std::string instName = inst->getName();
   instName.replace(instName.end()-4,instName.end(),""); //removes the .txt at end of string

   // write the fractional arc-time solution to a file
   //   char fname[100];

     ////std::cout << instName << endl;

   //   sprintf(fname, "Sols/%s-int.txt", instName.c_str());

      FILE *f = fopen(fname, "wt");
      if (f == NULL)
      {
         fprintf(stderr, "ERROR: Can't open the file %s for writing\n", fname);
     } 
     else
     {
        fprintf( f, "%d\n", sol.numArcsCap );
        for (int k = 0; k < sol.numArcsCap; k++)
          fprintf( f, " %d %d %d %lg\n", sol.arcsCap[k].i, sol.arcsCap[k].j,
               sol.arcsCap[k].d, sol.arcsCap[k].value );
        fclose(f);
     }      
}

void BBNode::solve( CPUTimer& t, double& secsLP, double& secsNode, int& cutRounds,
      bool doBranch )
{
   bool hasCuts = (level > 0);
   double currLB = 0.0;
   double prevLB = 0.0;
   int prevArcCnt = 1000000000;
   
//#define DANIEL_COUNT_FULL_ATIF_VARS
#ifdef DANIEL_COUNT_FULL_ATIF_VARS
   dwm->changePricing();
   firstArcs = dwm->getArcCount();
   fprintf( stderr,"inst:\t %s\n", inst->getName());
   fprintf( stderr,"firstArcs:\t %7.1lf\n", firstArcs);
   throw -1;
#endif
   double firstArcCount = double(dwm->getArcCount());
   CPUTimer tCuts;
   
   // Iterate adding cuts until no cut is added
   do
   {
      iteration = 0;
      int sbIter = 0;

      // Iterate adding columns until no column is added
#ifndef LOG_ON_TREE
      bool showLog = (level == 0);
#else
      bool showLog = true;
#endif
      if (showLog) fprintf( stderr,"%d: ", iteration);
      bool hasCols = dwm->doPricing( showLog, hasCuts );
      while ( hasCols || !dwm->getIsArcTimeIndexed() )
      {
         // If there is no column to insert, change the pricing
         if (!hasCols) dwm->changePricing();

         hasCuts = false;
         iteration++;
         dwm->solveRelaxation( showLog );
         currLB = dwm->getSolver()->getObjValue();
#ifndef LOG_ON_TREE
         showLog = ((level == 0) && ((iteration % 5) == 0));
#else
         showLog = ((iteration % 5) == 0);
#endif
         if (showLog) fprintf( stderr,"%d: ", iteration);
         if (iteration > 1000000) break;
         if (currLB < ub - 1.0 + IntegerTolerance)
         {
            sbIter++;
            if ((sbIter > MaxColGenStrBranch) && !doBranch) break;
         }

#ifdef ROOT_ONLY
#ifdef NO_CUT
         // stop if the time limit is reached
         t.stop();
         if (t.getCPUTotalSecs() > CpuTimeHardLimit)
         {
            // update the first LP time, number of iterations, LB and number of arcs
            secsLP = t.getCPUTotalSecs();
            firstNumIter = iteration;
            firstLB = currLB;
            firstArcs = dwm->getArcCount();
            missPricings = dwm->getMissPricings();
            stabChanges = dwm->getStabChanges();
            secsNode = t.getCPUTotalSecs();
            t.start();
            fprintf( stderr, "\nLEVEL %d: STOPPED BY TIME LIMIT.\n", level );
            throw( 'T' );
         }
         t.start();
#endif
#endif
         hasCols = dwm->doPricing( showLog, hasCuts );
      }
      if (doBranch)
         fprintf( stderr, "\nLEVEL %d: DWM LB = %10.3f, %d ARCS.\n", level,
               currLB, dwm->getArcCount() );

      if (cutRounds == 0)
      {
         // force expanding the fractional solution if no cut round has been done
         ExpandSol();

         // update the first LP time, number of iterations, LB and number of arcs
         t.stop();
         secsLP = t.getCPUTotalSecs();
         t.start();
         firstNumIter = iteration;
         firstLB = currLB;
         firstArcs = dwm->getArcCount();
         missPricings = dwm->getMissPricings();
         stabChanges = dwm->getStabChanges();

		 tCuts.start();
		//fprintf( stderr,"\nStarting cut rounds timer\n");
      }
      ub = dwm->getUB();

      // remove and add cuts
      if ((currLB < ub - 1.0 + IntegerTolerance) && doBranch)
      {
         // expand the fractional solution for this cut round
         //if (currLB >= prevLB + MinImproveToCleanLP)
         if ( (((currLB - prevLB) >= 0.2) && ((ub - currLB) < (ub - prevLB)*(1.0-MinGapReducCut))) ||
               (double(dwm->getArcCount()) < prevArcCnt*(1.0-MinArcReducCut)) )
            ExpandSol();
         dwm->removeUnusedCols();
         dwm->solveRelaxation(false);
         dwm->removeUnusedRows();
         dwm->solveRelaxation(false);
		 //Artur used >= 0.2 for the minimum gap decrease
         if ( (((currLB - prevLB) >= 0.2) && ((ub - currLB) < (ub - prevLB)*(1.0-MinGapReducCut))) ||
              (double(dwm->getArcCount()) < prevArcCnt*(1.0-MinArcReducCut)) )
		 {
            hasCuts = false;

#ifdef NO_CUT
            if (0)
#else
            if ((!isIntegerSol) && ((level == 0) ||
                  ((inst->jobs() <= MaxSizeNonRootCuts) &&
                   (cutRounds < MaxCutRoundsNonRoot))) )
#endif
            {
               // write the fractional arc-time solution to a file
               if (level == 0)
               {
                  char fname[100];
                  sprintf(fname, "lpSol%d.txt", cutRounds);
                  FILE *f = fopen(fname, "wt");
                  if (f != NULL)
                  {
                     fprintf( f, "%d\n", sol.numArcsCap );
                     for (int k = 0; k < sol.numArcsCap; k++)
                        fprintf( f, " %d %d %d %lg\n", sol.arcsCap[k].i, sol.arcsCap[k].j,
                              sol.arcsCap[k].d, sol.arcsCap[k].value );
                     fclose(f);
                  }
                  else
                     printf("WARNING: Can't open the file %s for writing\n", fname);
               }
            
               hasCuts = (GenerateCuts(cutRounds, t, currLB) > 0);
               //dwm->generateLpFile( nodeNumber );
               cutRounds++;
            }
         }
      }
      prevLB = currLB;
      prevArcCnt = dwm->getArcCount();
   }
   while (hasCuts);

   tCuts.stop();
   //fprintf( stderr,"Stopping cut rounds timer");

   //output cut generation results to latex table
    std::string instName = inst->getName();
   instName.replace(instName.end()-4,instName.end(),""); //removes the .txt at end of string
   double lb = dwm->getSolver()->getObjValue();
   if (lb > ub) lb = ub;

   LatexTable outputTable;
   
   outputTable.setCaption("Resumo - Teste Root dd/mm/aaaa");
   outputTable.setLabel("teste_dd_mm_resumo");

   outputTable.appendHeader("Instance"); //1
   outputTable.appendHeader("Heu UB"); //2
   outputTable.appendHeader("First LB"); //3
   outputTable.appendHeader("R. Arcs"); //4
   outputTable.appendHeader("Root LB"); //5
   outputTable.appendHeader("CutRounds"); //6
   outputTable.appendHeader("Time"); //7
   outputTable.appendHeader("R. Arcs"); //8
      
   outputTable.newLine();
   
   outputTable.appendField(instName); //1
   outputTable.appendField(ub); //2
   outputTable.appendField(firstLB);  //3
   outputTable.appendField(firstArcCount);  //4
   outputTable.appendField(lb);  //5
   outputTable.appendField(cutRounds); //6
   outputTable.appendField(tCuts.getCPUTotalSecs()); //7
   outputTable.appendField(double(dwm->getArcCount())); //8
 
   bool printHeader = false;
   bool printFormat = false;
   
   char fname[100];
#ifdef PRINT_TEX_GERAL
   sprintf(fname, "latex_tbl_output_geral.tex");
   FILE *f  = fopen(fname, "rt");
   if (f == NULL)
      {
      printHeader = true;
      printFormat = true;
   } else
   {
      fclose(f);
   }
   
   f = fopen(fname, "at");
   fprintf(f, "%s", outputTable.print(printHeader,printFormat).c_str() );
   fclose(f);
#endif //PRINT_TEX_GERAL

#ifdef PRINT_XL_GERAL
   //imprime excel
   printHeader = false;
   sprintf(fname, "xl_tbl_output_geral.txt");
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
#endif //PRINT_XL_GERAL

   // update the cut generation time
   t.stop();
   secsNode = t.getCPUTotalSecs();
   t.start();

   // stop if it must not branch
   if (!doBranch) return;

   // finished one more node
   runTimeNodeCount++;

   // stop if the upper bound cannot be improved
   if (currLB >= ub - 1.0 + IntegerTolerance)
   {
      if (level <= LevelToContinue) CpuTimeLimit = CpuTimeLastLimit;
      return;
   }
   
   if (level == 0)
   {
		dwm->generateArcTimeFixDataFile( nodeNumber );

		dwm->generateTimeIndexedLpFile( nodeNumber );

		char filename[120];
		sprintf( filename, "schedwtdata/%s_node-%d_time_indexed.lp", inst->getName(), nodeNumber );
		fprintf(stderr, "\nCPLEX COMMANDS:\n"
			"  set mip stra var 3\n"
			"  set mip tol mip 0.000001\n"
			"  r %s\n"
			"  set mip to u %d\n"
			"  opt\n"
			"  dis sol v -\n"
			"  quit\n\n",
			filename, int(ub));
		fprintf( stderr, "\nLEVEL %d: TIME INDEXED FORMULATION LP FILE GENERATED.\n", level );

		// check if can generate the LP for post-solving
		if (dwm->getArcCount() < MaxArcsLpGeneration)
		{	 
			dwm->generateLpFile( nodeNumber );
	 
			char filename[120];
			sprintf( filename, "%s_node-%d.lp", inst->getName(), nodeNumber );
			fprintf(stderr, "\nCPLEX COMMANDS:\n"
				"  set mip stra var 3\n"
				"  set mip tol mip 0.000001\n"
				"  r %s\n"
				"  set mip to u %d\n"
				"  opt\n"
				"  dis sol v -\n"
				"  quit\n\n",
				filename, int(ub));
			fprintf( stderr, "\nLEVEL %d: ARC-TIME INDEXED FORMULATION LP FILE GENERATED.\n", level );
		}
   }

		
  
   // stop if the time limit is reached
   if (secsNode > CpuTimeLimit)
   {
      fprintf( stderr, "\nLEVEL %d: STOPPED BY TIME LIMIT.\n", level );
      throw( 'T' );
   }

   // prepare the branching constraints
   fprintf( stderr, "\nLEVEL %d: BRANCHING... UB = %g\n\n", level, ub );
   int numConstr = dwm->getModel()->getNumConstraints();
   std::vector< std::vector<bool> > branchArcs;
   bool fathomLeft, fathomRight;
#if defined(STRONG_BRANCHING) && defined(SB_BRANCH_ON_TIME)
   int brJob, brTime;
   std::vector<bool> brSet;
#endif
   for (;;)
   {
      // WARNNING: branching performs better when using the same fractional
      // solution as the last cut separation (not updated after the cut insertions).
#if defined(STRONG_BRANCHING) && defined(SB_BRANCH_ON_TIME)
      setBranchingJobTime(brJob, brTime, fathomLeft, fathomRight, brSet);
#else
      setBranchingArcs(branchArcs, fathomLeft, fathomRight);
#endif
      if (fathomLeft && fathomRight)
      {
         if (level <= LevelToContinue) CpuTimeLimit = CpuTimeLastLimit;
         fprintf( stderr, "\nLEVEL %d: FINISHED BOTH BRANCHES WITH %d NODES. UB = %g\n\n",
               level, nodeCount, ub );
         return;
      }
      if (fathomLeft || fathomRight)
      {
         int iteration = 0;

         // add the branching constraint as a cut
         if (fathomLeft)
         {
            fprintf( stderr, "\nLEVEL %d: ADDING RIGHT BRANCH AS A CUT.\n",
                  level );
            dwm->getModel()->resizeConstraintArray(numConstr+1);
#if defined(STRONG_BRANCHING) && defined(SB_BRANCH_ON_TIME)
            ExtCCData* data = new ExtCCData;
            data->capacity = inst->T()-1;
            data->incapcoeff.resize(inst->T(), 0.0);
            data->outcapcoeff.resize(inst->T(), 0.0);
            data->S = brSet;
            data->weightS = inst->ptime()[brJob];
#ifdef SB_BRANCH_ON_SINGLE_TIME
            data->incapcoeff[brTime] = 1.0;
            Constraint* ctr = new RExtCapConstr(data, 1.0, true);
#else
            for (int t = brTime; t < inst->T(); t++)
               data->incapcoeff[t] = -1.0;
            Constraint* ctr = new RExtCapConstr(data, 0.0, true);
#endif
            ctr->canBeDeleted(false);
#else
            Constraint* ctr = new BinaryArcConstr(branchArcs, '>', 1.0);
#endif
            dwm->getModel()->setConstraint(numConstr, ctr);
            dwm->insertRowsInLP();
            numConstr = dwm->getModel()->getNumConstraints();
         }
         else
         {
            fprintf( stderr, "\nLEVEL %d: ADDING LEFT BRANCH AS A CUT.\n",
                  level );
            dwm->getModel()->resizeConstraintArray(numConstr+1);
#if defined(STRONG_BRANCHING) && defined(SB_BRANCH_ON_TIME)
            ExtCCData* data = new ExtCCData;
            data->capacity = inst->T()-1;
            data->incapcoeff.resize(inst->T(), 0.0);
            data->outcapcoeff.resize(inst->T(), 0.0);
            data->S = brSet;
            data->weightS = inst->ptime()[brJob];
#ifdef SB_BRANCH_ON_ACCUM_TIME
            for (int t = brTime; t < inst->T(); t++)
               data->incapcoeff[t] = 1.0;
            Constraint* ctr = new RExtCapConstr(data, 1.0, true);
#else
#ifdef SB_BRANCH_ON_SINGLE_TIME
            data->incapcoeff[brTime] = -1.0;
#else
            for (int t = 0; t < brTime; t++)
               data->incapcoeff[t] = -1.0;
#endif
            Constraint* ctr = new RExtCapConstr(data, 0.0);
#endif
            ctr->canBeDeleted(false);
#else
            Constraint* ctr = new BinaryArcConstr(branchArcs, '<', 0.0);
#endif
            dwm->getModel()->setConstraint(numConstr, ctr);
            dwm->insertRowsInLP();
            numConstr = dwm->getModel()->getNumConstraints();
         }
         hasCuts = true;

         // Iterate adding columns until no column is added
#ifndef LOG_ON_TREE
         bool showLog = (level == 0);
#else
         bool showLog = true;
#endif
         if (showLog) fprintf( stderr,"%d: ", iteration);
         while ( dwm->doPricing( showLog, hasCuts ) )
         {
            hasCuts = false;
            iteration++;
            dwm->solveRelaxation( showLog );
            currLB = dwm->getSolver()->getObjValue();
            showLog = ((level == 0) && ((iteration % 5) == 0));
            if (showLog) fprintf( stderr,"%d: ", iteration);
            if (iteration > 1000000) break;
         }
         ExpandSol();
         fprintf( stderr, "\nLEVEL %d: DWM LB = %10.3f, %d ARCS.\n", level,
               currLB, dwm->getArcCount() );
         ub = dwm->getUB();
      }
      else
         break;

      // stop if the upper bound cannot be improved
      if (currLB >= ub - 1.0 + IntegerTolerance)
      {
         if (level <= LevelToContinue) CpuTimeLimit = CpuTimeLastLimit;
         return;
      }
   }

   // initialize auxiliary variables
   double _secsLP;
   double _secsNode;
   int _cutRounds;

   // build the left node and solve it
   {
      BBNode left(*this);
      left.nodeNumber = nodeNumber + nodeCount;
      left.dwm->getModel()->resizeConstraintArray(numConstr+1);
#if defined(STRONG_BRANCHING) && defined(SB_BRANCH_ON_TIME)
      ExtCCData* data = new ExtCCData;
      data->capacity = inst->T()-1;
      data->incapcoeff.resize(inst->T(), 0.0);
      data->outcapcoeff.resize(inst->T(), 0.0);
      data->S = brSet;
      data->weightS = inst->ptime()[brJob];
#ifdef SB_BRANCH_ON_ACCUM_TIME
      for (int t = brTime; t < inst->T(); t++)
         data->incapcoeff[t] = 1.0;
      Constraint* ctr = new RExtCapConstr(data, 1.0, true);
#else
#ifdef SB_BRANCH_ON_SINGLE_TIME
      data->incapcoeff[brTime] = -1.0;
#else
      for (int t = 0; t < brTime; t++)
         data->incapcoeff[t] = -1.0;
#endif
      Constraint* ctr = new RExtCapConstr(data, 0.0);
#endif
      ctr->canBeDeleted(false);
#else
      Constraint* ctr = new BinaryArcConstr(branchArcs, '<', 0.0);
#endif
      left.dwm->getModel()->setConstraint(numConstr, ctr);
      left.dwm->insertRowsInLP();
      _cutRounds = 0;
      left.solve(t, _secsLP, _secsNode, _cutRounds);
      ub = left.ub;
      dwm->setUB( ub );
      nodeCount += left.nodeCount;
      fprintf( stderr, "\nLEVEL %d: FINISHED LEFT BRANCH WITH %d NODES. UB = %g\n\n",
            level, left.nodeCount, ub );
   }

   // wait more if a node at a low level is finished
   if (level+1 <= LevelToContinue) CpuTimeLimit = CpuTimeLastLimit;

   // build the right node if necessary and solve it
   if (currLB < ub - 1.0 + IntegerTolerance)
   {
      BBNode right(*this);
      right.nodeNumber = nodeNumber + nodeCount;
      right.dwm->getModel()->resizeConstraintArray(numConstr+1);
#if defined(STRONG_BRANCHING) && defined(SB_BRANCH_ON_TIME)
      ExtCCData* data = new ExtCCData;
      data->capacity = inst->T()-1;
      data->incapcoeff.resize(inst->T(), 0.0);
      data->outcapcoeff.resize(inst->T(), 0.0);
      data->S = brSet;
      data->weightS = inst->ptime()[brJob];
#ifdef SB_BRANCH_ON_SINGLE_TIME
      data->incapcoeff[brTime] = 1.0;
      Constraint* ctr = new RExtCapConstr(data, 1.0, true);
#else
      for (int t = brTime; t < inst->T(); t++)
         data->incapcoeff[t] = -1.0;
      Constraint* ctr = new RExtCapConstr(data, 0.0, true);
#endif
      ctr->canBeDeleted(false);
#else
      Constraint* ctr = new BinaryArcConstr(branchArcs, '>', 1.0);
#endif
      right.dwm->getModel()->setConstraint(numConstr, ctr);
      right.dwm->insertRowsInLP();
      _cutRounds = 0;
      right.solve(t, _secsLP, _secsNode, _cutRounds);
      ub = right.ub;
      dwm->setUB( ub );
      nodeCount += right.nodeCount;
      fprintf( stderr, "\nLEVEL %d: FINISHED RIGHT BRANCH WITH %d NODES. UB = %g\n\n",
            level, right.nodeCount, ub );
   }
}

int BBNode::GenerateCuts(int cutRound,  CPUTimer& t, double currLB)
{
   InstanceInfo ii;
   int* demand = 0;
   int numberECCH = 0, //ORIG_ECCSEP
      numberECECH = 0, 
      numberECCG = 0, //GENETIC_ECC_SEP e GENETIC_ECC_s_SEP
      numberTCC = 0, //TRICLIQUECUT_SEP
      numberGCC = 0, //GENERIC_CLIQUECUT_SEP
      numberOEC = 0; //OEC_SEP
   double sepDuration = 0;
   double diversidade = 0;
   int status = 0;
   CutList cutList;

   cutList.maxNumCuts = 2000;
   cutList.numCuts = 0;
   cutList.cuts = new Cut[2000];
   demand = new int[ inst->jobs()+1 ];
   demand[0] = 0;
   for( int i = 1; i <= inst->jobs(); i++ )
      demand[i] = inst->ptime()[i];


   if( cutGenPointer == 0 )
   {
      // create a cut generator object passing the instance information
      ii.capacity = inst->T();
      ii.numNodes = inst->jobs()+1;
      ii.demand = demand;
      ii.nrootbranches = inst->machines();

      cutGenPointer = InitCutGenerator( &ii );
   }


   // separate ECCs by both enumeration and heuristic
   ProblemType problemType = PROB_PATH;

   double avgViol = 0.0;
   double bestViol = 0.0;
   
   MaxSepCutsPerRound = MaxSepECECsPerRound;
   MaxInsCutsPerRound = MaxInsECECsPerRound;

   fprintf(stderr,"ECECsH = ");
   status = SeparateExtCyElimByHeur(cutGenPointer,&sol,&cutList,MaxSepCutsPerRound, problemType);
   numberECECH = cutList.numCuts;   
   for (int c = 0; c < cutList.numCuts; c++)
   {
      avgViol += cutList.cuts[c].violation / double(cutList.numCuts);
      if (bestViol < cutList.cuts[c].violation)
         bestViol = cutList.cuts[c].violation;
   }
   fprintf(stderr,"%d, viol: avg = %lg, best = %lg\n",numberECECH, avgViol, bestViol);
   ProcessCutList(&cutList);

   if (inst->machines() == 1)
   {
      MaxSepCutsPerRound = MaxSepECCsPerRoundSingle;
      MaxInsCutsPerRound = MaxInsECCsPerRoundSingle;
   }
   else
   {
      MaxSepCutsPerRound = MaxSepECCsPerRoundMulti;
      MaxInsCutsPerRound = MaxInsECCsPerRoundMulti;
   }

   CPUTimer tEcc;

   std::string instName = inst->getName();
   instName.replace(instName.end()-4,instName.end(),""); //removes the .txt at end of string

   //code to print a summary table of the cut separation round
   LatexTable outputTable;

   outputTable.setCaption("Resultados- Teste Root dd/mm/aaaa"); //set tex caption
   outputTable.setLabel("teste_dd_mm"); //set tex label

   //set table headers
   outputTable.appendHeader("Instance");
   outputTable.appendHeader("CutRound");
   outputTable.appendHeader("Curr LB");
   outputTable.appendHeader("Cut Family");
   outputTable.appendHeader("n");
   outputTable.appendHeader("Avg");
   outputTable.appendHeader("Best");
   outputTable.appendHeader("Time");
   outputTable.appendHeader("Diversity");

#ifdef ORIG_ECCSEP   
   if ( (level == 0) || (inst->jobs() <= MaxSizeNonRootRHECCCuts) )
   {
	   // write the fractional arc-time solution to a file
	   do
	   {
		 char fname[100];

		 //std::cout << instName << endl;

		  sprintf(fname, "Sols/%s-lpSol%d_%d.txt", instName.c_str(), nodeCount - 1, cutRound);

		  FILE *f = fopen(fname, "wt");
		  if (f == NULL)
		  {
			 fprintf(stderr, "ERROR: Can't open the file %s for writing\n", fname);
			 break;
		  }
		  fprintf( f, "%d\n", sol.numArcsCap );
		  for (int k = 0; k < sol.numArcsCap; k++)
			 fprintf( f, " %d %d %d %lg\n", sol.arcsCap[k].i, sol.arcsCap[k].j,
				   sol.arcsCap[k].d, sol.arcsCap[k].value );
		  fclose(f);
	   }     
	   while (false);

	   fprintf(stderr,"ECCsH = ");
	   tEcc.start();
	   status = SeparateECCbyHeuristic(cutGenPointer,&sol,&cutList,MaxSepCutsPerRound, problemType,2);
	   tEcc.stop();

	   numberECCH = cutList.numCuts;
	   avgViol = 0.0;
	   bestViol = 0.0;
	   for (int c = 0; c < cutList.numCuts; c++)
	   {
		  avgViol += cutList.cuts[c].violation / double(cutList.numCuts);
		  if (bestViol < cutList.cuts[c].violation)
			 bestViol = cutList.cuts[c].violation;
	   }
	   sepDuration = tEcc.getCPUTotalSecs();
	   diversidade = cutList.diversity;

	   fprintf(stderr,"%d, viol: avg = %lg, best = %lg, time = %lg, diversity = %lg\n",
			 numberECCH, avgViol, bestViol, sepDuration, diversidade);

	   ProcessCutList(&cutList);
      
	   outputTable.newLine();
   
	   outputTable.appendField(instName);
	   outputTable.appendField(cutRound);
	   outputTable.appendField(currLB);   
	   outputTable.appendField("ORIG_ECCSEP");
	   outputTable.appendField(numberECCH);
	   outputTable.appendField(avgViol);
	   outputTable.appendField(bestViol);
	   outputTable.appendField(sepDuration);
	   outputTable.appendField(diversidade);
   
   }
#endif //ORIG_ECCSEP

#ifdef GENETIC_ECC_SEP
   do
   {
      // write the fractional arc-time solution to a file
      char fname[100];

     //std::cout << instName << endl;

      sprintf(fname, "Sols/%s-lpSol%d_%d.txt", instName.c_str(), nodeCount - 1, cutRound);

      FILE *f = fopen(fname, "wt");
      if (f == NULL)
      {
         fprintf(stderr, "ERROR: Can't open the file %s for writing\n", fname);
         break;
      }
      fprintf( f, "%d\n", sol.numArcsCap );
      for (int k = 0; k < sol.numArcsCap; k++)
         fprintf( f, " %d %d %d %lg\n", sol.arcsCap[k].i, sol.arcsCap[k].j,
               sol.arcsCap[k].d, sol.arcsCap[k].value );
      fclose(f);
      
      // call the genetic ECC separator
      char cmd[200];
      sprintf(cmd, "GenSep/gensep %s %s %d >> sep.txt", inst->getFileName(), fname, inst->machines());
      fprintf(stderr, "Running %s...\nECCsG = ", cmd);
      tEcc.start();
      system(cmd);
      tEcc.stop();
      
      // retrieve the cuts from the file
      sprintf(fname, "Cuts/%s-lpSol%d_%d_cut.txt", instName.c_str(), nodeCount - 1, cutRound);
      f = fopen(fname, "rt");
      if (f == NULL)
      {
         fprintf(stderr, "ERROR: Can't open the file %s for reading\n", fname);
         break;
      }
      fscanf( f, "%d", &numberECCG );
      for (int c = 0; c < numberECCG; c++)
      {
         int num, den, setSize;
         fscanf( f, "%d %d %d", &num, &den, &setSize );
         std::vector<int> setDesc(setSize);
         for (int k = 0; k < setSize; k++) fscanf( f, "%d", &setDesc[k] );
         if (((CutGenerator*)cutGenPointer)->genSingleECCCut(&cutList, num, den, 0, &setDesc[0], setSize, 0.05))
         {
            if (cutList.numCuts >= MaxInsCutsPerRound) break;
         }
      }
     fscanf( f, "%lg", &sepDuration ); //inserida por Daniel Dias em 15-04-2015
     fscanf( f, "%lg", &diversidade ); //inserida por Daniel Dias em 17-04-2015
     t.stop();
     t.increaseCPUTotalSecs( sepDuration );
     t.start();
      fclose(f);
   }
   while (false);

   numberECCG = cutList.numCuts;
   avgViol = 0.0;
   bestViol = 0.0;
   for (int c = 0; c < cutList.numCuts; c++)
   {
      avgViol += cutList.cuts[c].violation / double(cutList.numCuts);
      if (bestViol < cutList.cuts[c].violation)
         bestViol = cutList.cuts[c].violation;
   }
   /*fprintf(stderr,"%d, viol: avg = %lg, best = %lg, time = %.1lf\n",
         numberECCG, avgViol, bestViol, tEcc.getCPUTotalSecs());*/ // comentado por Daniel Dias em 15-04-2015 
   
   fprintf(stderr,"%d, viol: avg = %lg, best = %lg, time = %lg, diversity = %lg\n",
         numberECCG, avgViol, bestViol, sepDuration, diversidade);  // codigo por Daniel Dias em 17-04-2015 

   ProcessCutList(&cutList);
      
   outputTable.newLine();
   
   outputTable.appendField(instName);
   outputTable.appendField(cutRound);
   outputTable.appendField(currLB);   
   outputTable.appendField("GEN_RECC");
   outputTable.appendField(numberECCG);
   outputTable.appendField(avgViol);
   outputTable.appendField(bestViol);
   outputTable.appendField(sepDuration);
   outputTable.appendField(diversidade);
#endif //GENETIC_ECCSEP

#ifdef TRICLIQUECUT_SEP
   do
   {
      // write the fractional arc-time solution to a file
      char fname[100];

     //std::cout << instName << endl;

		#ifdef WIN32
			sprintf(fname, "Sols\\%s-lpSol%d_%d.txt", instName.c_str(), nodeCount - 1, cutRound);
		#else
			sprintf(fname, "Sols/%s-lpSol%d_%d.txt", instName.c_str(), nodeCount - 1, cutRound);
		#endif

      FILE *f = fopen(fname, "wt");
      if (f == NULL)
      {
         fprintf(stderr, "ERROR: Can't open the file %s for writing\n", fname);
         break;
      }
      fprintf( f, "%d\n", sol.numArcsCap );
      for (int k = 0; k < sol.numArcsCap; k++)
         fprintf( f, " %d %d %d %lg\n", sol.arcsCap[k].i, sol.arcsCap[k].j,
               sol.arcsCap[k].d, sol.arcsCap[k].value );
      fclose(f);
      
      // call the triangle clique cuts separator
      char cmd[200];
		#ifdef WIN32
		sprintf(cmd, "TriCliqueSep\\tricliquesep %s %s %d >> sep.txt", inst->getFileName(), fname, inst->machines());
		#else
		sprintf(cmd, "TriCliqueSep/tricliquesep %s %s %d >> sep.txt", inst->getFileName(), fname, inst->machines());
		#endif

      fprintf(stderr, "Running %s...\nTCsG = ", cmd);
      tEcc.start();
      system(cmd);
      tEcc.stop();
      
      // retrieve the cuts from the file
	  #ifdef WIN32
		sprintf(fname, "Cuts\\%s-lpSol%d_%d_tccut.txt", instName.c_str(), nodeCount - 1, cutRound);
		#else
		sprintf(fname, "Cuts/%s-lpSol%d_%d_tccut.txt", instName.c_str(), nodeCount - 1, cutRound);
		#endif

      f = fopen(fname, "rt");
      if (f == NULL)
      {
         fprintf(stderr, "ERROR: Can't open the file %s for reading\n", fname);
         break;
      }
      fscanf( f, "%d", &numberTCC );
      for (int c = 0; c < numberTCC; c++)
      {
         int S_i, S_j, S_k;
         fscanf( f, "%d %d %d", &S_i, &S_j, &S_k );
         if (((CutGenerator*)cutGenPointer)->genSingleTriCliqueCut(&cutList, S_i, S_j, S_k, 0.05))
         {
            if (cutList.numCuts >= MaxInsCutsPerRound) break;
         }
      }
     fscanf( f, "%lg", &sepDuration ); //inserida por Daniel Dias em 15-04-2015
     t.stop();
     t.increaseCPUTotalSecs( sepDuration );
     t.start();
      fclose(f);
   }
   while (false);

   numberTCC = cutList.numCuts;
   avgViol = 0.0;
   bestViol = 0.0;
   for (int c = 0; c < cutList.numCuts; c++)
   {
      avgViol += cutList.cuts[c].violation / double(cutList.numCuts);
      if (bestViol < cutList.cuts[c].violation)
         bestViol = cutList.cuts[c].violation;
   }

   fprintf(stderr,"%d, viol: avg = %lg, best = %lg, time = %lg\n",
         numberTCC, avgViol, bestViol, sepDuration);  // codigo por Daniel Dias em 17-04-2015 

   ProcessCutList(&cutList);

   outputTable.newLine();
   
   outputTable.appendField(instName);
   outputTable.appendField(cutRound);
   outputTable.appendField(currLB);   
   outputTable.appendField("TRICLIQUE");
   outputTable.appendField(numberTCC);
   outputTable.appendField(avgViol);
   outputTable.appendField(bestViol);
   outputTable.appendField(sepDuration);
   outputTable.appendField(diversidade);
#endif //TRICLIQUECUT_SEP

#ifdef GENERIC_CLIQUECUT_SEP
   do
   {
      // write the fractional arc-time solution to a file
      char fname[100];

     //std::cout << instName << endl;

      sprintf(fname, "Sols/%s-lpSol%d_%d.txt", instName.c_str(), nodeCount - 1, cutRound);

      FILE *f = fopen(fname, "wt");
      if (f == NULL)
      {
         fprintf(stderr, "ERROR: Can't open the file %s for writing\n", fname);
         break;
      }
      fprintf( f, "%d\n", sol.numArcsCap );
      for (int k = 0; k < sol.numArcsCap; k++)
         fprintf( f, " %d %d %d %lg\n", sol.arcsCap[k].i, sol.arcsCap[k].j,
               sol.arcsCap[k].d, sol.arcsCap[k].value );
      fclose(f);
      
      // call the genetic ECC separator
      char cmd[200];
      sprintf(cmd, "GenCliqueSep/gencliquesep %s %s %d >> sep.txt", inst->getFileName(), fname, inst->machines());
      fprintf(stderr, "Running %s...\nGCsG = ", cmd);
      tEcc.start();
      system(cmd);
      tEcc.stop();
 
      // retrieve the cuts from the file
      sprintf(fname, "Cuts/%s-lpSol%d_%d_gccut.txt", instName.c_str(), nodeCount - 1, cutRound);
      f = fopen(fname, "rt");
      if (f == NULL)
      {
         fprintf(stderr, "ERROR: Can't open the file %s for reading\n", fname);
         break;
      }
      fscanf( f, "%d", &numberGCC );
      for (int c = 0; c < numberGCC; c++)
      {
         std::vector<int> x_i;
         std::vector<int> x_j;
         std::vector<int> x_d;
         int num_coef;
         fscanf( f, "%d", &num_coef );
       for (int n = 0; n < num_coef; n++)
       {
          int temp_i, temp_j, temp_d;
          fscanf( f, "%d %d %d", &temp_i, &temp_j, &temp_d );
          x_i.push_back(temp_i);
          x_j.push_back(temp_j);
          x_d.push_back(temp_d);
       }

       //genGenericCliqueCut( CutList* cuts, std::vector<int> x_i, std::vector<int> x_j, std::vector<int> x_d, double minViolation )
         if (((CutGenerator*)cutGenPointer)->genGenericCliqueCut(&cutList, x_i, x_j, x_d, 0.05))
         {
            if (cutList.numCuts >= MaxInsCutsPerRound) break;
         }
      }
     fscanf( f, "%lg", &sepDuration ); //inserida por Daniel Dias em 15-04-2015
     t.stop();
     t.increaseCPUTotalSecs( sepDuration );
     t.start();
      fclose(f);
   }
   while (false);

   numberGCC = cutList.numCuts;
   avgViol = 0.0;
   bestViol = 0.0;
   for (int c = 0; c < cutList.numCuts; c++)
   {
      avgViol += cutList.cuts[c].violation / double(cutList.numCuts);
      if (bestViol < cutList.cuts[c].violation)
         bestViol = cutList.cuts[c].violation;
   }
   fprintf(stderr,"%d, viol: avg = %lg, best = %lg, time = %lg\n",
         numberGCC, avgViol, bestViol, sepDuration);  // codigo por Daniel Dias em 17-04-2015 

   ProcessCutList(&cutList);

   outputTable.newLine();
   
   outputTable.appendField(instName);
   outputTable.appendField(cutRound);
   outputTable.appendField(currLB);   
   outputTable.appendField("GENERIC_CLIQUE");
   outputTable.appendField(numberGCC);
   outputTable.appendField(avgViol);
   outputTable.appendField(bestViol);
   outputTable.appendField(sepDuration);
   outputTable.appendField(diversidade);
#endif //GENERIC_CLIQUECUT_SEP

#ifdef GENETIC_ECC_s_SEP
   do
   {
      // write the fractional arc-time solution to a file
      char fname[100];
	  
#ifdef WIN32
      sprintf(fname, "Sols\\%s-lpSol%d_%d.txt", instName.c_str(), nodeCount - 1, cutRound);
#else
      sprintf(fname, "Sols/%s-lpSol%d_%d.txt", instName.c_str(), nodeCount - 1, cutRound);
#endif

      FILE *f = fopen(fname, "wt");
      if (f == NULL)
      {
         fprintf(stderr, "ERROR: Can't open the file %s for writing\n", fname);
         break;
      }
      fprintf( f, "%d\n", sol.numArcsCap );
      for (int k = 0; k < sol.numArcsCap; k++)
         fprintf( f, " %d %d %d %lg\n", sol.arcsCap[k].i, sol.arcsCap[k].j,
               sol.arcsCap[k].d, sol.arcsCap[k].value );
      fclose(f);
      
      // call the genetic ECC separator
      char cmd[200];
	#ifdef WIN32
      sprintf(cmd, "RECC-Sep\\reccsep %s %s %d >> sep.txt", inst->getFileName(), fname, inst->machines());
	#else
      sprintf(cmd, "RECC-Sep/reccsep %s %s %d >> sep.txt", inst->getFileName(), fname, inst->machines());
	#endif
     
	  fprintf(stderr, "Running %s...\nECCsG_s = ", cmd);
      tEcc.start();
      system(cmd);
      tEcc.stop();
      
      // retrieve the cuts from the file
		#ifdef WIN32
		sprintf(fname, "Cuts\\%s-lpSol%d_%d_recc_cut.txt", instName.c_str(), nodeCount - 1, cutRound);
		#else
		sprintf(fname, "Cuts/%s-lpSol%d_%d_recc_cut.txt", instName.c_str(), nodeCount - 1, cutRound);
		#endif

      f = fopen(fname, "rt");
      if (f == NULL)
      {
         fprintf(stderr, "ERROR: Can't open the file %s for reading\n", fname);
         break;
      }
      fscanf( f, "%d", &numberECCG );
      for (int c = 0; c < numberECCG; c++)
      {
         int num, den, s_num, setSize;
         fscanf( f, "%d %d %d %d", &num, &den, &s_num, &setSize );
         std::vector<int> setDesc(setSize);
         for (int k = 0; k < setSize; k++) fscanf( f, "%d", &setDesc[k] );
         if (((CutGenerator*)cutGenPointer)->genSingleECCCut_s(&cutList, num, den, s_num, 0, &setDesc[0], setSize, 0.05))
         {
            if (cutList.numCuts >= MaxInsCutsPerRound) break;
         }
      }
     fscanf( f, "%lg", &sepDuration ); //inserida por Daniel Dias em 15-04-2015
     fscanf( f, "%lg", &diversidade ); //inserida por Daniel Dias em 17-04-2015
     t.stop();
     t.increaseCPUTotalSecs( sepDuration );
     t.start();
      fclose(f);
   }
   while (false);

   numberECCG = cutList.numCuts;
   avgViol = 0.0;
   bestViol = 0.0;
   for (int c = 0; c < cutList.numCuts; c++)
   {
      avgViol += cutList.cuts[c].violation / double(cutList.numCuts);
      if (bestViol < cutList.cuts[c].violation)
         bestViol = cutList.cuts[c].violation;
   }
   fprintf(stderr,"%d, viol: avg = %lg, best = %lg, time = %lg, diversity = %lg\n",
         numberECCG, avgViol, bestViol, sepDuration, diversidade);  // codigo por Daniel Dias em 17-04-2015 

   ProcessCutList(&cutList);
   
   outputTable.newLine();
   
   outputTable.appendField(instName);
   outputTable.appendField(cutRound);
   outputTable.appendField(currLB);   
   outputTable.appendField("GEN_RECC_s");
   outputTable.appendField(numberECCG);
   outputTable.appendField(avgViol);
   outputTable.appendField(bestViol);
   outputTable.appendField(sepDuration);
   outputTable.appendField(diversidade);
#endif //GENETIC_ECC_s_SEP
      

#ifdef OEC_SEP
   do
   {
      // write the fractional arc-time solution to a file
      char fname[100];

     //std::cout << instName << endl;

	#ifdef WIN32
	sprintf(fname, "Sols\\%s-lpSol%d_%d.txt", instName.c_str(), nodeCount - 1, cutRound);
	#else
	sprintf(fname, "Sols/%s-lpSol%d_%d.txt", instName.c_str(), nodeCount - 1, cutRound);
	#endif

      FILE *f = fopen(fname, "wt");
      if (f == NULL)
      {
         fprintf(stderr, "ERROR: Can't open the file %s for writing\n", fname);
         break;
      }
      fprintf( f, "%d\n", sol.numArcsCap );
      for (int k = 0; k < sol.numArcsCap; k++)
         fprintf( f, " %d %d %d %lg\n", sol.arcsCap[k].i, sol.arcsCap[k].j,
               sol.arcsCap[k].d, sol.arcsCap[k].value );
      fclose(f);
      
      // call the genetic ECC separator
#ifndef INTEGRATED_OEC_SEP
      char cmd[200];
	  #ifdef WIN32
      sprintf(cmd, "OEC-Sep\\oecsep %s %s %d >> sep.txt", inst->getFileName(), fname, inst->machines());
	  #else
      sprintf(cmd, "OEC-Sep/oecsep %s %s %d >> sep.txt", inst->getFileName(), fname, inst->machines());
	  #endif
      fprintf(stderr, "Running %s...\nOECsG = ", cmd);
      tEcc.start();
      system(cmd);
      tEcc.stop();
#else
      tEcc.start();      
	  status = SeparateOECbyHeuristic(cutGenPointer,&sol,&cutList,MaxSepCutsPerRound,problemType,2);
      tEcc.stop();	 
#endif
 
    // retrieve the cuts from the file
#ifndef INTEGRATED_OEC_SEP
	#ifdef WIN32
      sprintf(fname, "Cuts\\%s-lpSol%d_%d_oec_cut.txt", instName.c_str(), nodeCount - 1, cutRound);
	#else
      sprintf(fname, "Cuts/%s-lpSol%d_%d_oec_cut.txt", instName.c_str(), nodeCount - 1, cutRound);
	#endif
      f = fopen(fname, "rt");
      if (f == NULL)
      {
         fprintf(stderr, "ERROR: Can't open the file %s for reading\n", fname);
         break;
      }
      fscanf( f, "%d", &numberOEC );
      for (int c = 0; c < numberOEC; c++)
      {
         int tCorte;
         int setSize;
         fscanf( f, "%d %d", &tCorte, &setSize );
         std::vector<int> setDesc(setSize);
         for (int k = 0; k < setSize; k++) fscanf( f, "%d", &setDesc[k] );         

         //bool genSingleOECCut( CutList* cuts, int t, int *SetList, int SetListSize, double minViolation );
         if (((CutGenerator*)cutGenPointer)->genSingleOECCut(&cutList, tCorte, inst->machines(), &setDesc[0], setSize, 0.05))
         {
         if (cutList.numCuts >= MaxInsCutsPerRound) break;
         }
      }
     fscanf( f, "%lg", &sepDuration ); //inserida por Daniel Dias em 15-04-2015
     t.stop();
     t.increaseCPUTotalSecs( sepDuration );
     t.start();
     fclose(f);
#endif
   }
   while (false);

   numberOEC = cutList.numCuts;
   avgViol = 0.0;
   bestViol = 0.0;
   for (int c = 0; c < cutList.numCuts; c++)
   {
	 // fprintf(stderr,"OEC(%d) viol: %lg\n",
		//  c, cutList.cuts[c].violation);     
      avgViol += cutList.cuts[c].violation / double(cutList.numCuts);
      if (bestViol < cutList.cuts[c].violation)
         bestViol = cutList.cuts[c].violation;
   }
 
   fprintf(stderr,"%d, viol: avg = %lg, best = %lg, time = %lg\n",
         numberOEC, avgViol, bestViol, sepDuration);  // codigo por Daniel Dias em 17-04-2015 

   ProcessCutList(&cutList);
   
   outputTable.newLine();
   
   outputTable.appendField(instName);
   outputTable.appendField(cutRound);
   outputTable.appendField(currLB);   
   outputTable.appendField("OEC");
   outputTable.appendField(numberOEC);
   outputTable.appendField(avgViol);
   outputTable.appendField(bestViol);
   outputTable.appendField(sepDuration);
   outputTable.appendField(diversidade);
#endif //OEC_SEP

#ifdef PRINT_TEX
   bool printHeader = false;
   bool printFormat = false;

   char fname[100];
    sprintf(fname, "latex_tbl_output.tex");
   FILE *f  = fopen(fname, "rt");
    if (f == NULL)
    {
        printHeader = true;
        printFormat = true;
   } else
   {
      fclose(f);
   }
#endif //PRINT_TEX
   
#ifdef PRINT_XL      //imprime excel   
   bool printHeader_xl = false;
   bool printFormat_xl = false;

   char fname_xl[100];   

   printHeader_xl = false;
   sprintf(fname_xl, "xl_tbl_output.txt");
   FILE *f_xl = fopen(fname_xl, "rt");
   if (f_xl == NULL)
   {
      printHeader_xl = true;
   } 
   else
   {
      fclose(f_xl);
   }
   
   f_xl = fopen(fname_xl, "at");
   fprintf(f_xl , "%s", outputTable.printExcel(printHeader_xl).c_str() );
   fclose(f_xl);
#endif //PRINT_XL
   
   delete [] demand;
   delete [] cutList.cuts;

   // add the cuts to the DWM LP if any, and return
   if (numberECCH + numberECECH + numberECCG + numberTCC + numberGCC + numberOEC > 0) dwm->insertRowsInLP();
   return numberECCH + numberECECH + numberECCG + numberTCC + numberGCC + numberOEC;
}

void BBNode::ExpandSol()
{
   sol.numArcs = 0;
   sol.numArcsCap = 0;
   const double* xsol = dwm->getColSolution();
   std::map<CapArcKey,double,CapArcKeyComp> arcMap;
   std::map<CapArcKey,double,CapArcKeyComp> capArcMap;
   std::map<CapArcKey,double,CapArcKeyComp>::iterator it;
   fprintf( stderr, "qroutes start at col %d out of %d cols\n",
         dwm->getIdxCQRouteVars(), dwm->getNumCols() );

   // count the fractional routes
   //double routeCount = 0.0;
   //for(int c = dwm->getIdxCQRouteVars(); c < dwm->getNumCols(); c++ )
   //   if( xsol[c] > ColExpandEps )
   //      routeCount += xsol[c];
   //double factor = double(inst->machines()) / routeCount;
   double factor = 1.0;

   // for all the columns
   isIntegerSol = true;
   for(int c = dwm->getIdxCQRouteVars(); c < dwm->getNumCols(); c++ )
   {
      // take only the ones with non-zero primal values
      if( xsol[c] > ColExpandEps )
      {
         // check if the solution is not integer
         if (xsol[c] < 1.0 - ColExpandEps) isIntegerSol = false;
         //fprintf( stderr, "x_%d = %10.8f\n", c, xsol[c] );

         // for all arcsin the corresponding qroute
         int length;
         CapArcKey* qroute = dwm->getCQRoute(c, length);
         if (length == 0)
         {
            fprintf( stderr, "Error: artificial var %d with val %g\n", c, xsol[c] );
            throw -984836;
         }
         for (int a = 0; a < length; a++)
         {
            // add the arc in the map if necessary
            CapArcKey key;
            key.i = qroute[a].i;
            key.j = qroute[a].j;
            key.d = 0;
            it = arcMap.find( key );
            if (it != arcMap.end())
               it->second += xsol[c] * factor;
            else
               arcMap.insert( std::pair<CapArcKey,double>(key,xsol[c] * factor) );

            // add the capacitated arc in the map if necessary
            key.d = qroute[a].d;
            it = capArcMap.find( key );
            if (it != capArcMap.end())
               it->second += xsol[c] * factor;
            else
               capArcMap.insert( std::pair<CapArcKey,double>(key,xsol[c] * factor) );
         }
      }
   }

   bool printOpt = false;
   // print a message if the solution is integer
   if (isIntegerSol && (dwm->getSolver()->getObjValue() < ub - IntegerTolerance))
   {
      ub = dwm->getSolver()->getObjValue();
      dwm->setUB( ub );
      fprintf( stderr, "INTEGER SOLUTION OF %g FOUND!\n", ub );
      for(int c = dwm->getIdxCQRouteVars(); c < dwm->getNumCols(); c++)
      {
         // take only the ones with non-zero primal values
         if( xsol[c] > ColExpandEps )
         {
            // for all arcsin the corresponding qroute
            int length;
            CapArcKey* qroute = dwm->getCQRoute(c, length);
            for (int a = length-1; a > 0; a--)
               fprintf( stderr, " %d", qroute[a].i );
            fprintf( stderr, ".\n" );
         }
      }
	  printOpt = true;
   }

   // release the previous expanded solution is any
   if (sol.arcs != 0) delete [] sol.arcs;
   sol.arcs = 0;
   if (sol.arcsCap != 0) delete [] sol.arcsCap;
   sol.arcsCap = 0;

   // store the fractional primal solution
   //printf( "\nArcs in sol:\n" );
   sol.arcs = new CUTS_ArcVariable[arcMap.size()];
   for (it = arcMap.begin(); it != arcMap.end(); it++)
   {
      sol.arcs[sol.numArcs].value = it->second;
      sol.arcs[sol.numArcs].i = it->first.i;
      sol.arcs[sol.numArcs].j = it->first.j;
      //printf( "x(%3d,%3d)=%4.2f ", sol.arcs[sol.numArcs].i, sol.arcs[sol.numArcs].j,
      //      sol.arcs[sol.numArcs].value );
      sol.numArcs++;
   }
   //printf( "\nExt arcs in sol:\n" );
   sol.arcsCap = new CUTS_ArcCapVariable[capArcMap.size()];
   for (it = capArcMap.begin(); it != capArcMap.end(); it++)
   {
      sol.arcsCap[sol.numArcsCap].value = it->second;
      sol.arcsCap[sol.numArcsCap].i = it->first.i;
      sol.arcsCap[sol.numArcsCap].j = it->first.j;
      sol.arcsCap[sol.numArcsCap].d = it->first.d;
      //printf( "(%3d,%3d,%4d)=%4.2f ", sol.arcsCap[sol.numArcsCap].i, sol.arcsCap[sol.numArcsCap].j,
      //      sol.arcsCap[sol.numArcsCap].d, sol.arcsCap[sol.numArcsCap].value );
      sol.numArcsCap++;
   }
   //printf( "\n" ); fflush(stdout);

    // writes solution if the solution is integer
   if (printOpt)
   {
		std::string instName = inst->getName();
		instName.replace(instName.end()-4,instName.end(),""); //removes the .txt at end of string
		char fname[100];

		#ifdef WIN32
		sprintf(fname, "Sols\\%s-opt.txt", instName.c_str());
		#else
		sprintf(fname, "Sols/%s-opt.txt", instName.c_str());
		#endif

		FILE *f = fopen(fname, "wt");
		if (f == NULL)
		{
			fprintf(stderr, "ERROR: Can't open the file %s for writing\n", fname);
		}
		else
		{
			fprintf( f, "%d\n", sol.numArcsCap );
			for (int k = 0; k < sol.numArcsCap; k++)
				fprintf( f, " %d %d %d %lg\n", sol.arcsCap[k].i, sol.arcsCap[k].j,
					sol.arcsCap[k].d, sol.arcsCap[k].value );
			fclose(f);
		}
   }
}

struct CutComparison
{
   bool operator()(Cut a, Cut b)
   {
      return (a.violation > b.violation);
   }
};

void BBNode::ProcessCutList( CutList* cutList )
{
   // sort the cut list and remove the less violated ones if necessary
   if (cutList->numCuts > MaxInsCutsPerRound)
   {
      sort( &(cutList->cuts[0]), &(cutList->cuts[cutList->numCuts]), CutComparison());
      for (int i = MaxInsCutsPerRound; i < cutList->numCuts; i++)
         ResetCut( &(cutList->cuts[i]) );
      cutList->numCuts = MaxInsCutsPerRound;
   }

   // allocate memory for the new constraints
   int nrows = dwm->getModel()->getNumConstraints();
   dwm->getModel()->resizeConstraintArray(nrows + cutList->numCuts);

   // insert the new constraints in the DWM
   for (int i = 0; i < cutList->numCuts; i++)
   {
      if (strcmp(cutList->cuts[i].templabel, "HECC") == 0)
      {
         ExtCCData* eccData = (ExtCCData*)cutList->cuts[i].data;
         cutList->cuts[i].data = 0;
         Constraint* constr = new RExtCapConstr(eccData, cutList->cuts[i].rhs);
         dwm->getModel()->setConstraint(nrows + i, constr);
      }
      else
      {
         GenericCutData* cutData = (GenericCutData*)cutList->cuts[i].data;
         cutList->cuts[i].data = 0;
         Constraint* constr = new GenericConstr(cutData, cutList->cuts[i].sense,
               cutList->cuts[i].rhs);
         dwm->getModel()->setConstraint(nrows + i, constr);
      }
      // fprintf( stderr, "Inserting cut %d with violation %g\n", i, cutList->cuts[i].violation );
   }

   // clear the list of cuts
   ClearCutList( cutList );
}

struct ArcComparison
{
   ArcComparison( std::vector<double>* costs )
   {
      _costs = costs;
   }

   bool operator()(CUTS_ArcVariable a, CUTS_ArcVariable b)
   {
      return (_costs[a.i][a.j] > _costs[a.i][a.j]);
   }

   std::vector<double>* _costs;
};



#ifdef DANIEL_BRANCHING

void BBNode::setBranchingArcs( std::vector< std::vector<bool> >& arcs,
         bool& fathomLeft, bool& fathomRight )
{
   fathomLeft = fathomRight = false;

   // do strong branching if the gap is not too small
   if ((ub - dwm->getSolver()->getObjValue()) >= (1.0 + IntegerTolerance))
      fprintf( stderr, "\nDOING DANIEL BRANCHING...\n\n" );

   int menor_t = sol.arcsCap[0].d;
   double maiorValor = 0;

   fprintf( stderr, "Listando variaveis\n");
   for (int a = 0; a < sol.numArcsCap; a++)
   {
      int i = sol.arcsCap[a].i;
      int j = sol.arcsCap[a].j;
      int t = sol.arcsCap[a].d;
      double valor = sol.arcsCap[a].value;
      
      fprintf( stderr, "x%d_%d_%d = %g\n", j, i, t, valor);
   }

   for (int a = 0; a < sol.numArcsCap; a++)
   {
      int i = sol.arcsCap[a].i;
      int j = sol.arcsCap[a].j;
      int t = sol.arcsCap[a].d;
      double valor = sol.arcsCap[a].value;

     if (valor < 1)
     {
        if (t < menor_t)
        {
           menor_t = t;
           fprintf( stderr, "menor_t: %d\n", menor_t);
           fprintf( stderr, "x%d_%d_%d = %g\n", j, i, t, valor);
        }
     }
   }   

   int fixa_a = 0;
   for (int a = 0; a < sol.numArcsCap; a++)
   {
      int i = sol.arcsCap[a].i;
      int j = sol.arcsCap[a].j;
      int t = sol.arcsCap[a].d;
      double valor = sol.arcsCap[a].value;

     if (t == menor_t) 
     {
        if ((valor < 1) && (valor > maiorValor))
        {
         maiorValor = sol.arcsCap[a].value;
         fixa_a = a;
         fprintf( stderr, "fixa_a: %d\n", a);
         fprintf( stderr, "x%d_%d_%d = %g\n", j, i, t, valor);
        }
     }
   }   


   // initialize the boolean matrix of arcs
   arcs.clear();
   arcs.resize(inst->jobs()+1);
   for (int i = 0; i <= inst->jobs(); i++)
      arcs[i].resize(inst->jobs()+1, false);

   int i = sol.arcsCap[fixa_a].i;
   int j = sol.arcsCap[fixa_a].j;

   arcs[i+1][j+1] = true;

  fprintf( stderr, "Branching arco: (%d, %d) \n", j, i);
}

#else
#ifdef STRONG_BRANCHING

struct BranchJobComp
{
   std::vector<double> avgTimes;
   bool operator() (int i, int j) { return (avgTimes[i] < avgTimes[j]); }
};

void BBNode::setBranchingJobTime( int& job, int& jtime, bool& fathomLeft,
      bool& fathomRight, std::vector<bool>& jobSet )
{
   fathomLeft = fathomRight = false;

   // do strong branching if the gap is not too small
   if ((ub - dwm->getSolver()->getObjValue()) >= (1.0 + IntegerTolerance))
      fprintf( stderr, "\nDOING STRONG BRANCHING...\n\n" );

   // calculate the fraction of each job finishing at each time in the
   // relaxation solution
   std::vector< std::vector<double> > jobTimes( inst->jobs() + 1 );
   for (int i = 0; i <= inst->jobs(); i++)
      jobTimes[i].resize( inst->T(), 0.0 );
   for (int a = 0; a < sol.numArcsCap; a++)
   {
      // add the value of the current capacitated arc to the corresponding
      // job and time
      int i = sol.arcsCap[a].i;
      int j = sol.arcsCap[a].j;
      int t = sol.arcsCap[a].d;
      jobTimes[j][t] += sol.arcsCap[a].value;
   }

   // initialize the order for evaluating the branching jobs
   std::vector<int> ord( inst->jobs());
   for (int k = 0; k < inst->jobs(); k++) ord[k] = k+1;

   // initialize the middle times and scores used for branching
   std::vector<int> middleTime( inst->jobs() + 1, -1 );
   std::vector<double> brScores( inst->jobs() + 1, 0.0 );

#ifdef SB_BRANCH_ON_ACCUM_TIME
   // sort the branching jobs by the average finishing time
   BranchJobComp brJobCmp;
   brJobCmp.avgTimes.resize( inst->jobs() + 1, 0.0 );
   for (int i = 1; i <= inst->jobs(); i++)
      for (int t = 0; t < inst->T(); t++)
         brJobCmp.avgTimes[i] += double(t) * jobTimes[i][t];
   std::sort(ord.begin(), ord.end(), brJobCmp);

   // calculate the middle time and the score
   std::vector<double> accum( inst->jobs() + 1, 0.0 );
   for (int t = inst->T() - 1; t >= 0; t--)
   {
      double sum = 0.0;
      double score = 0.0;
      for (int k = 0; k < inst->jobs(); k++)
      {
         int i = ord[k];
         if (middleTime[i] != -1) break;
         sum += jobTimes[i][t];
         accum[i] += sum;
         brScores[i] += sum * double(t);
         if (accum[i] >= TargetBrTimeValue)
         {
            if (((accum[i] - sum) > IntegerTolerance) && (accum[i] > 1.0 - IntegerTolerance))
            {
               accum[i] -= sum;
               middleTime[i] = t+1;
               brScores[i] = (brScores[i] - accum[i] * double(t)) * (1.0 - accum[i]); //  * double(t+1);
            }
            else
            {
               middleTime[i] = t;
               brScores[i] = (brScores[i] - accum[i] * double(t-1)) * (1.0 - accum[i]); //  * double(t);
            }
         }
      }
   }
#else
#ifdef SB_BRANCH_ON_SINGLE_TIME
   double avgCost = double(ub) / double(inst->jobs());
#endif
   // calculate the median time and the average distance to the median for
   // each job (use as a branching score)
   for (int k = 0; k < inst->jobs(); k++)
   {
      int i = ord[k];
#ifdef SB_BRANCH_ON_SINGLE_TIME
      middleTime[i] = 0;
      for (int t = 1; t < inst->T(); t++)
      {
         if (jobTimes[i][t] > jobTimes[i][middleTime[i]])
            middleTime[i] = t;
      }
      if (jobTimes[i][middleTime[i]] < TargetBrTimeValue)
         brScores[i] = (double(inst->getCost(i,middleTime[i])) + avgCost) *
               jobTimes[i][middleTime[i]] / TargetBrTimeValue;
      else
         brScores[i] = (double(inst->getCost(i,middleTime[i])) + avgCost) *
               (1.0 - jobTimes[i][middleTime[i]]) / (1.0 - TargetBrTimeValue);
#else
      int prev = -1;
      double accum = 0.0;
      double distZero = 0.0;
      for (int t = 0; t < inst->T(); t++)
      {
         accum += jobTimes[i][t];
         if ((accum >= (1.0 - TargetBrTimeValue)) && (jobTimes[i][t] > IntegerTolerance) &&
               (prev != -1) && (middleTime[i] == -1))
         {
            middleTime[i] = (t + prev + 1) / 2;
            brScores[i] = double(middleTime[i]) * accum - distZero;
         }
         if (middleTime[i] != -1)
            brScores[i] += double(t - middleTime[i] + 1) * jobTimes[i][t];
         distZero += double(t) * jobTimes[i][t];
         if (jobTimes[i][t] > IntegerTolerance)
            prev = t;
      }
#endif
   }
#endif

   // choose the candidates for strong branching
   std::vector<double> bestScores(NumStrBrCandidates, IntegerTolerance);
   std::vector<int> bestJobs(NumStrBrCandidates, -1);
   for (int i = 1; i <= inst->jobs(); i++)
   {
      // update the array of best minimum estimated gains
      int worse = 0;
      for (int k = 1; k < NumStrBrCandidates; k++)
      {
         if (bestScores[k] < bestScores[worse])
            worse = k;
      }
      if (bestScores[worse] < brScores[i])
      {
         bestScores[worse] = brScores[i];
         bestJobs[worse] = i;
      }
   }

   // for each vertex...
   double bestMinGain = 0.0;
   int bestJob = 0;
   int bestTime = 0;
   for (int k = 0; k < NumStrBrCandidates; k++)
   {
      int i = bestJobs[k];
      if (i == -1) continue;
#ifdef SB_BRANCH_ON_ACCUM_TIME
      double leftGain = (1.0 - accum[i]);
      double rightGain = accum[i];
#else
#ifdef SB_BRANCH_ON_SINGLE_TIME
      double leftGain = (double(inst->getCost(i,middleTime[i])) + avgCost) *
            jobTimes[i][middleTime[i]] / TargetBrTimeValue;
      double rightGain = (double(inst->getCost(i,middleTime[i])) + avgCost) *
            (1.0 - jobTimes[i][middleTime[i]]) / (1.0 - TargetBrTimeValue);
#else
      double leftGain = 0.0;
      double rightGain = 0.0;
      for (int t = 0; t < inst->T(); t++)
      {
         if (t < middleTime[i])
            leftGain += jobTimes[i][t];
         else
            rightGain += jobTimes[i][t];
      }
#endif
#endif

      // use the estimated gains if the GAP is too small
      if ((ub - dwm->getSolver()->getObjValue()) < (1.0 + IntegerTolerance))
      {
         double minGain = leftGain;
         if (minGain > rightGain) minGain = rightGain;
         if (minGain > bestMinGain)
         {
            bestMinGain = minGain;
            bestJob = i;
            bestTime = middleTime[i];
         }
      }

      // solve the root LP for each created branch
      else
      {
         // initialize auxiliary variables
         CPUTimer _t;
         _t.start();
         double _secsLP;
         double _secsNode;
         int _cutRounds = 0;
         int numConstr = dwm->getModel()->getNumConstraints();

         // build the left node and solve its root LP only
         {
            BBNode left(*this);
            left.dwm->getModel()->resizeConstraintArray(numConstr+1);
            ExtCCData* data = new ExtCCData;
            data->capacity = inst->T()-1;
            data->incapcoeff.resize(inst->T(), 0.0);
            data->outcapcoeff.resize(inst->T(), 0.0);
            data->weightS = inst->ptime()[i];
            data->S.resize(inst->jobs()+1, false);
            data->S[i] = true;
#ifdef SB_BRANCH_ON_ACCUM_TIME
            for (int kk = 0; ord[kk] != i; kk++) data->S[ord[kk]] = true;
            for (int t = middleTime[i]; t < inst->T(); t++)
               data->incapcoeff[t] = 1.0;
            Constraint* ctr = new RExtCapConstr(data, 1.0, true);
#else
#ifdef SB_BRANCH_ON_SINGLE_TIME
            data->incapcoeff[middleTime[i]] = -1.0;
#else
            for (int t = 0; t < middleTime[i]; t++)
               data->incapcoeff[t] = -1.0;
#endif
            Constraint* ctr = new RExtCapConstr(data, 0.0);
#endif
            left.dwm->getModel()->setConstraint(numConstr, ctr);
            left.dwm->insertRowsInLP();
            left.solve(_t, _secsLP, _secsNode, _cutRounds, false);
            ub = left.ub;
            dwm->setUB( ub );
            double approx = leftGain;
            leftGain = left.dwm->getSolver()->getObjValue();
            fprintf( stderr, "STRONG BRANCHING LEFT PROBE: j=%d, t=%d,"
                  " DWM LB = %9.2f (%g)\n\n", i, middleTime[i], leftGain,
                  approx );
            if (leftGain >= ub - 1.0 + IntegerTolerance)
               fathomLeft = true;
         }

         // build the right node and solve its root LP only
         {
            BBNode right(*this);
            right.dwm->getModel()->resizeConstraintArray(numConstr+1);
            ExtCCData* data = new ExtCCData;
            data->capacity = inst->T()-1;
            data->incapcoeff.resize(inst->T(), 0.0);
            data->outcapcoeff.resize(inst->T(), 0.0);
            data->S.resize(inst->jobs()+1, false);
            data->S[i] = true;
#ifdef SB_BRANCH_ON_ACCUM_TIME
            for (int kk = 0; ord[kk] != i; kk++) data->S[ord[kk]] = true;
#endif
            data->weightS = inst->ptime()[i];
#ifdef SB_BRANCH_ON_SINGLE_TIME
            data->incapcoeff[middleTime[i]] = 1.0;
            Constraint* ctr = new RExtCapConstr(data, 1.0, true);
#else
            for (int t = middleTime[i]; t < inst->T(); t++)
               data->incapcoeff[t] = -1.0;
            Constraint* ctr = new RExtCapConstr(data, 0.0, true);
#endif
            right.dwm->getModel()->setConstraint(numConstr, ctr);
            right.dwm->insertRowsInLP();
            right.solve(_t, _secsLP, _secsNode, _cutRounds, false);
            ub = right.ub;
            dwm->setUB( ub );
            _t.stop();
            double approx = rightGain;
            rightGain = right.dwm->getSolver()->getObjValue();
            fprintf( stderr, "STRONG BRANCHING RIGHT PROBE: j=%d, t=%d,"
                  " DWM LB = %9.2f (%g)\n\n", i, middleTime[i], rightGain,
                  approx );
            if (rightGain >= ub - 1.0 + IntegerTolerance)
               fathomRight = true;
         }

         // stop if fathomed in any side
         if (fathomLeft || fathomRight)
         {
            job = i;
            jobSet.clear();
            jobSet.resize(inst->jobs() + 1, false);
            jobSet[i] = true;
#ifdef SB_BRANCH_ON_ACCUM_TIME
            for (int kk = 0; ord[kk] != i; kk++) jobSet[ord[kk]] = true;
#endif
            jtime = middleTime[i];
            return;
         }

         // update the branching choice
         double minGain = leftGain;
         if (minGain > rightGain) minGain = rightGain;
         if (minGain > bestMinGain)
         {
            bestMinGain = minGain;
            bestJob = i;
            bestTime = middleTime[i];
         }
      }
   }

   if (bestJob == 0)
   {
      fprintf( stderr, "ERROR: no branching found!\n" );
      for (int k = 0; k < inst->jobs(); k++)
      {
         int j = ord[k];
         fprintf( stderr, "j=%d:", j );
         for (int t = 0; t < inst->T(); t++)
            if (jobTimes[j][t] > 1e-12)
               fprintf( stderr, " (%d,%lg)", t, jobTimes[j][t] );
         fprintf( stderr, "\n" );
      }

      exit(-1);
   }

   job = bestJob;
   jobSet.clear();
   jobSet.resize(inst->jobs() + 1, false);
   jobSet[bestJob] = true;
#ifdef SB_BRANCH_ON_ACCUM_TIME
   for (int kk = 0; ord[kk] != bestJob; kk++) jobSet[ord[kk]] = true;
#endif
   jtime = bestTime;
   fprintf( stderr, "Branching choice: j=%d, t=%d\n", bestJob, bestTime );
}

void BBNode::setBranchingArcs( std::vector< std::vector<bool> >& arcs,
         bool& fathomLeft, bool& fathomRight )
{
   fathomLeft = fathomRight = false;

   // do strong branching if the gap is not too small
   if ((ub - dwm->getSolver()->getObjValue()) >= (1.0 + IntegerTolerance))
      fprintf( stderr, "\nDOING STRONG BRANCHING...\n\n" );

   // calculate the cost of each arc preferring the ones having smaller or
   // larger original costs
   std::vector< std::vector<double> > costsSmall( inst->jobs() + 1 );
   std::vector< std::vector<double> > costsLarge( inst->jobs() + 1 );
   for (int i = 0; i <= inst->jobs(); i++)
   {
      costsSmall[i].resize( inst->jobs() + 1, 1E-5 );
      costsLarge[i].resize( inst->jobs() + 1, 1E-5 );
   }
   double costJob = FracCostAddBranch * ub / double(inst->jobs());
   for (int a = 0; a < sol.numArcsCap; a++)
   {
      // add the cost of the current capacitated arc to the corresponding arc
      int i = sol.arcsCap[a].i;
      int j = sol.arcsCap[a].j;
      int t = sol.arcsCap[a].d;
      costsSmall[i][j] += (costJob / (costJob + inst->getCost(j, t))) *
            sol.arcsCap[a].value;
      costsLarge[i][j] += ((costJob + inst->getCost(j, t)) / costJob) *
            sol.arcsCap[a].value;
   }

   // choose the candidates for strong branching
   std::vector<double>* costs;
   std::vector<CUTS_ArcVariable> adjArcs(inst->jobs());
   std::vector<double> bestSmallest(NumStrBrCandidates, 0.0);
   std::vector<double> bestLargest(NumStrBrCandidates, 0.0);
   double* bestMinGains;
   for (int i = 1; i <= inst->jobs(); i++)
   {
      for (int out = 0; out < 2; out++)
      {
         for (int pref = 0; pref < 2; pref++)
         {
            // choose the criterion to be considered
            if (pref == 0)
            {
               costs = &costsSmall[0];
               bestMinGains = &bestSmallest[0];
            }
            else
            {
               costs = &costsLarge[0];
               bestMinGains = &bestLargest[0];
            }

            // build an array of adjacent arcs sorted non-increasing by costs
            int k = 0;
            for (int j = 0; j <= inst->jobs(); j++)
            {
               if (j == i) continue;
               if (out == 0)
               {
                  adjArcs[k].i = j;
                  adjArcs[k].j = i;
               }
               else
               {
                  adjArcs[k].i = i;
                  adjArcs[k].j = j;
               }
               k++;
            }
            sort(adjArcs.begin(), adjArcs.end(), ArcComparison(costs));

            // divide the costs between left and right branching
            // -> sorting helps the balance being not so bad although suboptimal
            double leftGain = 0.0;
            double rightGain = 0.0;
            for (k = 0; k < inst->jobs(); k++)
            {
               if (leftGain <= rightGain)
                  leftGain += costs[adjArcs[k].i][adjArcs[k].j];
               else
                  rightGain += costs[adjArcs[k].i][adjArcs[k].j];
            }

            // update the array of best minimum estimated gains
            double minGain = leftGain;
            if (minGain > rightGain) minGain = rightGain;
            int worse = 0;
            for (k = 1; k < NumStrBrCandidates; k++)
            {
               if (bestMinGains[k] < bestMinGains[worse])
                  worse = k;
            }
            if (bestMinGains[worse] < minGain)
               bestMinGains[worse] = minGain;
         }
      }
   }

   // calculate the worse gain that should be considered as candidate for each
   // criterion
   double worseGain[2];
   for (int pref = 0; pref < 2; pref++)
   {
      // choose the criterion to be considered
      if (pref == 0)
         bestMinGains = &bestSmallest[0];
      else
         bestMinGains = &bestLargest[0];

      // find the worse relevant gain
      worseGain[pref] = bestMinGains[0];
      for (int k = 1; k < NumStrBrCandidates; k++)
      {
         if (bestMinGains[k] < worseGain[pref])
            worseGain[pref] = bestMinGains[k];
      }
   }

   // for each vertex, input or output arcs, and each criterion...
   double bestMinGain = 0.0;
   int bestJob = 0;
   int bestPref = 0;
   bool bestIsOut = false;
   for (int i = 1; i <= inst->jobs(); i++)
   {
      for (int out = 0; out < 2; out++)
      {
         for (int pref = 0; pref < 2; pref++)
         {
            // choose the criterion to be considered
            if (pref == 0)
            {
               costs = &costsSmall[0];
               bestMinGains = &bestSmallest[0];
            }
            else
            {
               costs = &costsLarge[0];
               bestMinGains = &bestLargest[0];
            }

            // build an array of adjacent arcs sorted non-increasing by costs
            int k = 0;
            for (int j = 0; j <= inst->jobs(); j++)
            {
               if (j == i) continue;
               if (out == 0)
               {
                  adjArcs[k].i = j;
                  adjArcs[k].j = i;
               }
               else
               {
                  adjArcs[k].i = i;
                  adjArcs[k].j = j;
               }
               k++;
            }
            sort(adjArcs.begin(), adjArcs.end(), ArcComparison(costs));

            // initialize the boolean matrix of arcs
            arcs.clear();
            arcs.resize(inst->jobs()+1);
            for (int j = 0; j <= inst->jobs(); j++)
               arcs[j].resize(inst->jobs()+1, false);

            // divide the costs between left and right branching
            // -> sorting helps the balance being not so bad although suboptimal
            double leftGain = 0.0;
            double rightGain = 0.0;
            for (k = 0; k < inst->jobs(); k++)
            {
               if (leftGain <= rightGain)
               {
                  leftGain += costs[adjArcs[k].i][adjArcs[k].j];
                  arcs[adjArcs[k].i][adjArcs[k].j] = true;
               }
               else
                  rightGain += costs[adjArcs[k].i][adjArcs[k].j];
            }

            // use the estimated gains if the GAP is too small
            if ((ub - dwm->getSolver()->getObjValue()) < (1.0 + IntegerTolerance))
            {
               // update the branching choice if the first criterion is used
               if (pref == 0)
               {
                  double minGain = leftGain;
                  if (minGain > rightGain) minGain = rightGain;
                  if (minGain > bestMinGain)
                  {
                     bestMinGain = minGain;
                     bestJob = i;
                     bestIsOut = (out == 1);
                     bestPref = pref;
                  }
               }
            }

            // solve the root LP for each created branch
            else if ((leftGain >= worseGain[pref]) && (rightGain >= worseGain[pref]))
            {
               // initialize auxiliary variables
               CPUTimer _t;
               _t.start();
               double _secsLP;
               double _secsNode;
               int _cutRounds = 0;
               int numConstr = dwm->getModel()->getNumConstraints();

               // build the left node and solve its root LP only
               {
                  BBNode left(*this);
                  left.dwm->getModel()->resizeConstraintArray(numConstr+1);
                  Constraint* ctr = new BinaryArcConstr(arcs, '<', 0.0);
                  left.dwm->getModel()->setConstraint(numConstr, ctr);
                  left.dwm->insertRowsInLP();
                  left.solve(_t, _secsLP, _secsNode, _cutRounds, false);
                  ub = left.ub;
                  dwm->setUB( ub );
                  double approx = leftGain;
                  leftGain = left.dwm->getSolver()->getObjValue();
                  fprintf( stderr, "STRONG BRANCHING LEFT PROBE: j=%d(%s, %s),"
                        " DWM LB = %9.2f (%g)\n\n", i, (out == 1)? "out": "in",
                        (pref == 0)? "small": "large", leftGain, approx );
                  if (leftGain >= ub - 1.0 + IntegerTolerance)
                     fathomLeft = true;
               }

               // build the right node and solve its root LP only
               {
                  BBNode right(*this);
                  right.dwm->getModel()->resizeConstraintArray(numConstr+1);
                  Constraint* ctr = new BinaryArcConstr(arcs, '>', 1.0);
                  right.dwm->getModel()->setConstraint(numConstr, ctr);
                  right.dwm->insertRowsInLP();
                  right.solve(_t, _secsLP, _secsNode, _cutRounds, false);
                  ub = right.ub;
                  dwm->setUB( ub );
                  _t.stop();
                  double approx = rightGain;
                  rightGain = right.dwm->getSolver()->getObjValue();
                  fprintf( stderr, "STRONG BRANCHING RIGHT PROBE: j=%d(%s, %s),"
                        " DWM LB = %9.2f (%g)\n\n", i, (out == 1)? "out": "in",
                        (pref == 0)? "small": "large", rightGain, approx );
                  if (rightGain >= ub - 1.0 + IntegerTolerance)
                     fathomRight = true;
               }

               // stop if fathomed in any side
               if (fathomLeft || fathomRight)
                  return;

               // update the branching choice
               double minGain = leftGain;
               if (minGain > rightGain) minGain = rightGain;
               if (minGain > bestMinGain)
               {
                  bestMinGain = minGain;
                  bestJob = i;
                  bestIsOut = (out == 1);
                  bestPref = pref;
               }
            }
         }
      }
   }

   // initialize the boolean matrix of arcs
   arcs.clear();
   arcs.resize(inst->jobs()+1);
   for (int i = 0; i <= inst->jobs(); i++)
      arcs[i].resize(inst->jobs()+1, false);

   // choose the criterion to be considered
   if (bestPref == 0)
      costs = &costsSmall[0];
   else
      costs = &costsLarge[0];

   // select the arcs to be fixed on branching
   double bestLeft = 0.0;
   double bestRight = 0.0;
   int k = 0;
   for (int j = 0; j <= inst->jobs(); j++)
   {
      if (j == bestJob) continue;
      if (!bestIsOut)
      {
         adjArcs[k].i = j;
         adjArcs[k].j = bestJob;
      }
      else
      {
         adjArcs[k].i = bestJob;
         adjArcs[k].j = j;
      }
      k++;
   }
   sort(adjArcs.begin(), adjArcs.end(), ArcComparison(costs));
   for (k = 0; k < inst->jobs(); k++)
   {
      if (bestLeft <= bestRight)
      {
         bestLeft += costs[adjArcs[k].i][adjArcs[k].j];
         arcs[adjArcs[k].i][adjArcs[k].j] = true;
      }
      else
         bestRight += costs[adjArcs[k].i][adjArcs[k].j];
   }
   fprintf( stderr, "Branching values: left = %g, right = %g\n", bestLeft, bestRight );
}

#else // STRONG_BRANCHING

void BBNode::setBranchingArcs( std::vector< std::vector<bool> >& arcs,
         bool& fathomLeft, bool& fathomRight )
{
   fathomLeft = fathomRight = false;

   // calculate the cost of each arc
   std::vector< std::vector<double> > costs( inst->jobs() + 1 );
   for (int i = 0; i <= inst->jobs(); i++)
      costs[i].resize( inst->jobs() + 1, 1E-5 );
   double costJob = FracCostAddBranch * ub / double(inst->jobs());
   for (int a = 0; a < sol.numArcsCap; a++)
   {
      // add the cost of the current capacitated arc to the corresponding arc
      int i = sol.arcsCap[a].i;
      int j = sol.arcsCap[a].j;
      int t = sol.arcsCap[a].d;
      costs[i][j] += (costJob + inst->getCost(j, t)) * sol.arcsCap[a].value;
   }

   // for each vertex...
   double bestMinGain = 0.0;
   int bestJob = 0;
   bool bestIsOut = false;
   std::vector<CUTS_ArcVariable> adjArcs(inst->jobs());
   for (int i = 1; i <= inst->jobs(); i++)
   {
      for (int out = 0; out < 2; out++)
      {
         // build an array of adjacent arcs sorted non-increasing by costs
         int k = 0;
         for (int j = 0; j <= inst->jobs(); j++)
         {
            if (j == i) continue;
            if (out == 0)
            {
               adjArcs[k].i = j;
               adjArcs[k].j = i;
            }
            else
            {
               adjArcs[k].i = i;
               adjArcs[k].j = j;
            }
            k++;
         }
         sort(adjArcs.begin(), adjArcs.end(), ArcComparison(&costs[0]));

         // divide the costs between left and right branching
         // -> sorting helps the balance being not so bad although suboptimal
         double left = 0.0;
         double right = 0.0;
         for (k = 0; k < inst->jobs(); k++)
         {
            if (left <= right)
               left += costs[adjArcs[k].i][adjArcs[k].j];
            else
               right += costs[adjArcs[k].i][adjArcs[k].j];
         }
         double minGain = left;
         if (minGain > right) minGain = right;

         // update the branching choice
         if (minGain > bestMinGain)
         {
            bestMinGain = minGain;
            bestJob = i;
            bestIsOut = (out == 1);
         }
      }
   }

   // initialize the boolean matrix of arcs
   arcs.clear();
   arcs.resize(inst->jobs()+1);
   for (int i = 0; i <= inst->jobs(); i++)
      arcs[i].resize(inst->jobs()+1, false);

   // select the arcs to be fixed on branching
   double bestLeft = 0.0;
   double bestRight = 0.0;
   int k = 0;
   for (int j = 0; j <= inst->jobs(); j++)
   {
      if (j == bestJob) continue;
      if (!bestIsOut)
      {
         adjArcs[k].i = j;
         adjArcs[k].j = bestJob;
      }
      else
      {
         adjArcs[k].i = bestJob;
         adjArcs[k].j = j;
      }
      k++;
   }
   sort(adjArcs.begin(), adjArcs.end(), ArcComparison(&costs[0]));
   for (k = 0; k < inst->jobs(); k++)
   {
      if (bestLeft <= bestRight)
      {
         bestLeft += costs[adjArcs[k].i][adjArcs[k].j];
         arcs[adjArcs[k].i][adjArcs[k].j] = true;
      }
      else
         bestRight += costs[adjArcs[k].i][adjArcs[k].j];
   }
   fprintf( stderr, "Branching values: left = %g, right = %g\n", bestLeft, bestRight );
}

#endif // STRONG_BRANCHING else
#endif // DANIEL BRANCHING else
