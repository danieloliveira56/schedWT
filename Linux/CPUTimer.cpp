#include "CPUTimer.h"
#include <time.h>
#include <sys/resource.h>

#include <stdio.h>

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
CPUTimer::CPUTimer()
{
  CPUCurrSecs         = 0;
  CPUTotalSecs        = 0;
  CronoCurrSecs       = 0;
  CronoTotalSecs      = 0;

  started             = false;
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
  CPUTStart = zeit();
  CronoTStart = real_zeit();
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
    CPUTStop = zeit();
    CronoTStop = real_zeit();
    CPUTotalSecs += CPUTStop - CPUTStart;
    CronoTotalSecs += CronoTStop - CronoTStart;
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
}

void CPUTimer::operator += ( CPUTimer  t )
{
  CPUCurrSecs += t.getCPUCurrSecs();
  CPUTotalSecs += t.getCPUTotalSecs();
  CronoCurrSecs += t.getCronoCurrSecs();
  CronoTotalSecs += t.getCronoTotalSecs();

}

double CPUTimer::zeit (void)
{
   struct rusage ru;

   getrusage (RUSAGE_SELF, &ru);

   return ((double) ru.ru_utime.tv_sec) +
      ((double) ru.ru_utime.tv_usec) / 1000000.0;
}

double CPUTimer::real_zeit (void)
{
   return (double) time (0);
}
