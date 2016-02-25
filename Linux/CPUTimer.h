#ifndef CPUTIME_H
#define CPUTIME_H


#include <stdlib.h>

class CPUTimer
{
  public:
    CPUTimer();



    // Retorna o tempo (em segs e msegs) de CPU cronometrado para uma rotina.
    // Se apenas uma cronometragem foi realizada, entao os valores retornados
    // por getCPUCurrSecs() e getCPUTtotalSecs sao iguais.
    double     getCPUCurrSecs();

    // Retorna o tempo total (em segs e msegs) de CPU cronometrado para uma rotina
    double     getCPUTotalSecs();

    // Retorna o tempo (em segs e msegs) de execucao cronometrado para uma rotina.
    // Se apenas uma cronometragem foi realizada, entao os valores retornados
    // por getCPUCurrSecs() e getCPUTtotalSecs sao iguais.
    double     getCronoCurrSecs();

    // Retorna o tempo total (em segs e msegs) de execucao cronometrado para uma rotina.
    double     getCronoTotalSecs();

    // Inicia a cronometragem (tempo de execucao e de CPU) de uma rotina
    bool       start();

    // Encerra a cronometragem (tempo de execucao e de CPU) de uma rotina
    bool       stop();

    // Prepara o ambiente de cronometragem para ser utilizado em outra rotina
    void       reset();


    // Operator to add cputimers
    void operator +=( CPUTimer t );

    inline void increaseCPUTotalSecs( double s){CPUTotalSecs += s; };

    bool         started;             // is the timer started?
  private:

    double CPUTStart;            // the start time
    double CPUTStop;             // the stop time
    double CronoTStart;
    double CronoTStop;
    double zeit();
    double real_zeit();

    double       CPUCurrSecs;         // tempo de cpu cronometrado para uma rotina (segs e msegs)
    double       CPUTotalSecs;        // total do tempo de cpu cronometrado para uma rotina (segs e msegs)

    double       CronoCurrSecs;       // tempo de execucao cronometrado para uma rotina (segs e msegs)
    double       CronoTotalSecs;      // total do tempo de execucao cronometrado para uma rotina (segs e msegs)

    bool         gottime;             // do we have a measured time we can return?
};
#endif

