#ifndef _CMSTDEBUG_H
#define _CMSTDEBUG_H

#include <stdio.h>

/** Debug functions should be called with a parenthesis on the string to be
   printed. This is because MSVC 6 does not support macros with a variable
   number of arguments */

#ifndef __func__
#define __func__ ""
#endif

#define __CMST_PRINT_LOC__(F) fprintf(((F==0)?stdout:F), " in (%s:%d)\n", __FILE__, __LINE__)

#define CMST_EXIT(A,REST) {\
 printf REST ; \
 __CMST_PRINT_LOC__(stdout); \
 exit(A); }
 
#define CMST_GOTO(A,REST) {\
 printf REST ; \
 __CMST_PRINT_LOC__(stdout); \
 goto A;}

#define CMST_IF_EXIT(A,B,REST) {\
 if(A) {\
 printf REST ; \
 __CMST_PRINT_LOC__(stdout); \
 exit(B);}}

#define CMST_CHECKRVAL(A,B) {\
 if(A) {\
   __CMST_PRINT_LOC__(stdout); \
   return B; } }

#define CMST_IF_WARNING(A, REST) {\
  if(A) {\
     printf REST ; \
      __CMST_PRINT_LOC__(stdout); \
      }}

#define CMST_IF_GOTO(A,B,REST) {\
 if(A) CMST_GOTO(B, REST) }

#define CMST_THROW(A,REST) {\
   printf REST ; \
   __CMST_PRINT_LOC__(stdout); \
   throw(A); }

#define CMST_IF_THROW(A,B,REST) {\
if(A) CMST_THROW(B,REST) ; }
  

#endif
