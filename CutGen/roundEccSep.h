#ifndef __ROUND_ECC_SEP_H__
#define __ROUND_ECC_SEP_H__

#include <map>
#include <vector>

#define ECC_VIOLATED_EPS      0.05

#define COL_EXPAND_EPS        1e-6

#define LONG_DOUBLE  double

struct Rational
{
  int num;
  int den;
};

struct RationalCmp
{
  bool operator()( const Rational& r1, const Rational& r2 ) const
  {
    return ( ((long long)(r1.num) * (long long)(r2.den)) <
             ((long long)(r2.num) * (long long)(r1.den)) );
  }
};

typedef std::map<Rational, LONG_DOUBLE, RationalCmp> RationalMap;

struct YZChange
{
   // changed values for y variables
   std::vector<int> demY;
   std::vector<double> valY;

   // changed values for z variables
   std::vector<int> demZ;
   std::vector<double> valZ;

   // previous set demand
   int prevSetDem;
};

void RECC_DoRoundings( YZChange& chg, RationalMap& IncBeforeMap,
      RationalMap& IncAfterMap, int Cap, int SetDemand, int& BestNum,
      int& BestDenom );

#endif // __ROUND_ECC_SEP_H__
