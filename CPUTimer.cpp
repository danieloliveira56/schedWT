#include "CPUTimer.h"
#include <time.h>

#include <stdio.h>

 
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
CPUTimer::CPUTimer()
{
  CPUCurrSecs         = 0;
  CPUTotalSecs        = 0;
  CronoCurrSecs       = 0;
  CronoTotalSecs      = 0;

  started             = false;

   //---------------------------------------------------------------------
   hFile = GetCurrentProcess();

   GetProcessTimes( hFile,
                  &FT_CreationTime,
                  &FT_ExitTime,
                  &FT_KernelTime,
                  &FT_UserTime );

   FileTimeToLocalFileTime( &FT_CreationTime,
                          &FT_Aux1Time );
}



//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
double CPUTimer::getCPUCurrSecs()
{
   return CPUCurrSecs;
}


//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
double CPUTimer::getCPUTotalSecs()
{
   return CPUTotalSecs;
}


//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
double CPUTimer::getCronoCurrSecs()
{
   return CronoCurrSecs;
}


//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
double CPUTimer::getCronoTotalSecs()
{
   return CronoTotalSecs;
}


//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
bool CPUTimer::start()
{
   bool status = true;

   CPUCurrSecs   = 0;
   CronoCurrSecs = 0;

   //---------------------------------------------------------------------
   GetProcessTimes( hFile,
                  &FT_CreationTime,
                  &FT_ExitTime,
                  &FT_KernelTime,
                  &FT_UserTime );
   FileTimeToLocalFileTime( &FT_KernelTime, &FT_Aux1Time );
   cpu_ti  = *(ULONG64 *) &FT_Aux1Time;
   FileTimeToLocalFileTime( &FT_UserTime, &FT_Aux1Time );
   cpu_ti += *(ULONG64 *) &FT_Aux1Time;

   GetSystemTimeAsFileTime( &FT_Aux1Time );
   FileTimeToLocalFileTime( &FT_Aux1Time, &FT_CronoStartTime );
   FileTimeToSystemTime( &FT_CronoStartTime, &ST_CronoStartTime );
   crono_ti = *(ULONG64 *) &FT_CronoStartTime;
   //---------------------------------------------------------------------

   gottime = false;
   started = status;

   return( status );
}


//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
bool CPUTimer::stop()
{
   bool status = true;

   if (started)
   {
      // Tempo de CPU
      // ----------------------------------------------------------------------------
      GetProcessTimes( hFile,
                       &FT_CreationTime,
                       &FT_ExitTime,
                       &FT_KernelTime,
                       &FT_UserTime );
      FileTimeToLocalFileTime( &FT_KernelTime, &FT_Aux1Time );
      cpu_tf  = *(ULONG64 *) &FT_Aux1Time;
      FileTimeToLocalFileTime( &FT_UserTime, &FT_Aux1Time );
      cpu_tf += *(ULONG64 *) &FT_Aux1Time;

      aux  = (LONG64) (cpu_tf - cpu_ti);
      FT_Aux1Time = *(FILETIME *) &aux;
      FileTimeToSystemTime( &FT_Aux1Time, &ST_CPUCurrTime );

      CPUCurrSecs   = (double) (aux / 10000000U) + ((aux % 10000000U) * 0.0000001);
      CPUTotalSecs += CPUCurrSecs;

      secs  = (LONG64) (CPUTotalSecs * 10000000U);
      FT_Aux1Time = *(FILETIME *) &secs;
      FileTimeToSystemTime( &FT_Aux1Time, &ST_CPUTotalTime );
      // ----------------------------------------------------------------------------


      // Tempo cronológico
      // ----------------------------------------------------------------------------
      GetSystemTimeAsFileTime( &FT_Aux1Time );
      FileTimeToLocalFileTime( &FT_Aux1Time, &FT_CronoStopTime );
      FileTimeToSystemTime( &FT_CronoStopTime, &ST_CronoStopTime );
      crono_tf = *(ULONG64 *) &FT_CronoStopTime;

      aux = (LONG64) (crono_tf - crono_ti);
      CronoCurrSecs   = (double) (aux / 10000000U) + ((aux % 10000000U) * 0.0000001);
      CronoTotalSecs += CronoCurrSecs;

      secs  = (LONG64) (CronoTotalSecs * 10000000U);
      FT_Aux1Time = *(FILETIME *) &secs;
      FileTimeToSystemTime( &FT_Aux1Time, &ST_CronoTotalTime );

      gottime = true;
      // ----------------------------------------------------------------------------
   }
   else
   {
      fprintf(stderr,"CPUTimer::stop(): called without calling CPUTimer::start() first!\n");
      status = false;
   }
   started = false;

   return status;
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
void CPUTimer::reset()
{
   started = false;

   CPUCurrSecs    = 0;
   CPUTotalSecs   = 0;
   CronoCurrSecs  = 0;
   CronoTotalSecs = 0;

   memset( &ST_CronoStartTime, 0, sizeof(SYSTEMTIME) );
   memset( &ST_CronoStopTime,  0, sizeof(SYSTEMTIME) );
   memset( &ST_CronoTotalTime, 0, sizeof(SYSTEMTIME) );
}

void CPUTimer::operator += ( CPUTimer  t )
{
   CPUCurrSecs += t.getCPUCurrSecs();
   CPUTotalSecs += t.getCPUTotalSecs();
   CronoCurrSecs += t.getCronoCurrSecs();
   CronoTotalSecs += t.getCronoTotalSecs();
}

