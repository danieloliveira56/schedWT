////////////////////////////////////////////
// File: Constants.h                      //
////////////////////////////////////////////
// File with the definitions of the       //
// constants to the branch and cut and    //
// price (exceptions and other constants) //
////////////////////////////////////////////
// Authors:       Marcelo Ladeira Reis    //
//                Ricardo Fukasawa        //
////////////////////////////////////////////
// Created on: 02/10/2003                 //
// Last modified: 02/10/2003              //
////////////////////////////////////////////

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

//---------------------------------------------------------------------------//

// Exceptions numbers

// Memory cannot be allocated
#define LOWMEMORYEXCEPTION 1

// Some memory was not deallocated
#define MEMNOTDEALLOCEXCEPTION 2

// Cplex environment coud not be initialized
#define SOLVERNOTINITIALIZEDEXCEPTION 3

// Could not open a file
#define FILENOTOPENEXCEPTION 4

// LP variable of solver is not initialized
#define LPVARIABLENOTINITIALIZED 5

// The LP has no solution
#define LPHASNOSOLUTIONEXCEPTION 6

// Exception while trying to retrieve the objective function val
#define LPSOLVEROBJFUNCTIONEXCEPTION 7

// Exception while trying to retrieve the solution
#define LPSOLVERSOLUTIONEXCEPTION 8

// Exception while trying to retrieve the dual solution
#define LPSOLVERDUALSOLEXCEPTION 9

// Trying to access an invalid position of an array
#define ARRAYINVALIDPOSITION 10

// Master is not initialized
#define MASTERNOTINITEXCEPTION 11

// Node is not initialized
#define NODENOTINITEXCEPTION 12

// Exception while trying to retrieve the rhs of lp
#define LPSOLVERRHSEXCEPTION 13

// Exception while trying to change a coefficient of a variable
#define LPSOLVERCHGCOEFFEXCEPTION 14

// Exception while trying to get a coefficient of a variable
#define LPSOLVERGETCOEFFEXCEPTION 15

// Exception while trying to create new columns
#define LPSOLVERNEWCOLSEXCEPTION 16

// Exception while trying to create new rows
#define LPSOLVERNEWROWSEXCEPTION 17

// Number of columns of the basis is wrong
#define BASISCOLSIZEEXCEPTION 18

// Number of rows of the basis is wrong
#define BASISROWSIZEEXCEPTION 19

// Cannot get the basis of the lp
#define GETBASISEXCEPTION 20

// Cannot copy the basis of the lp
#define COPYBASISEXCEPTION 21

// Number of parameters passed to the program is invalid
#define PROMPTINVALIDPARAMNUMEXCEPTION 22

// Invalid problem type
#define PROMPTINVALIDPTYPEEXCEPTION 23

// Lower bound of the son lower than the father
#define LOWSONLOWERBOUND 24

// Solver error while trying to delete rows
#define SOLVERDELROWSEXCEPTION 25

// Solver error while trying to delete cols
#define SOLVERDELCOLSEXCEPTION 26

// Solver error while trying to set an integer parameter
#define SOLVERSETINTPARAMEXCEPTION 27

// Solver error while trying to solve MIP
#define SOLVERMIPEXCEPTION 28

// Solver error while trying to get the MIP solution
#define SOLVERMIPXEXCEPTION 29

// Solver error while trying to get the MIP obj value
#define SOLVERMIPOBJEXCEPTION 30

// Invalid command line argument
#define INVALIDINPUTOPTION 31

// Invalid value of k (k-cycles)
#define PROMPTINVALIDCAPEXCEPTION 32

// Invalid value of kcycle when setting/querying column
#define PROMPTINVALIDHQ 33



#define PROMPTINVALIDPARAMPMIN 34

#define PROMPTINVALIDPARAMPMAX 35

#define PROMPTINVALIDPARAMSBCGITER 36

// FUKA - DUALSTAB - BEGIN
#define LPSOLVERCHGUBEXCEPTION 37
// FUKA - DUALSTAB - END


// FUKA2 - BEGIN
#define SBNOBRANCHSETS 38
// FUKA2 - END

#define LPSOLVERGETCOEFEXCEPTION 39


//---------------------------------------------------------------------------//
// DEFAULT INPUT PARAMETERS FOR THE COMMAND LINE
//---------------------------------------------------------------------------//
#define DEFAULTKCYCLE 3

#define DEFAULTPMIN   5

#define DEFAULTPMAX   8

#define DEFAULTSBCGITERS 20000

#define DEFAULTMAXECCDEPTH 5

//---------------------------------------------------------------------------//

// Tolerance default values
#define OBJEPS      1e-3
#define INTEPS      1e-6
#define CUTEPS      1e-3
#define CGEPS       1e-6
#ifdef ACVRP_PROBLEM
#define CGEPS_MULT  1e-5    // value of CGEPS when costs as multiplied by 1000
#endif
#define ARTIFEPS    1e-8    // artificial variables (non-configurable)
#define DUALEPS     1e-11   // dual variables (non-configurable)
#define PRIMALEPS   1e-11   // primal variables (non-configurable)

//---------------------------------------------------------------------------//

// Default value of infinity
#define INFINITY_VALUE 10000000

//---------------------------------------------------------------------------//

// Pool size default values for clean operation
#define MAX_COL_POOL_SIZE 5000
#define MAX_CUT_POOL_SIZE 1000

//---------------------------------------------------------------------------//

/**
 * Node status for a node of the branch-and-cut-and-price.
 */
enum NodeStatus{ FEASIBLE, FEASIBLE_IMPROVED, BRANCH, FATHOMED, INFEASIBLE };

//---------------------------------------------------------------------------//

/** Defines if the code is to be compiled with extra checking code to
 * debug it or not */
#define CMSTDEBUG 1

#define MEMDEBUG 0

//---------------------------------------------------------------------------//

#define VAREPS      1e-8    // primal variables (non-configurable)
#define DUALEPS     1e-11   // dual variables (non-configurable)

#endif
