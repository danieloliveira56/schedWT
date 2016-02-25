#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "roundEccSep.h"

#define MAX_CAP_HECC 3000

void RECC_DoRoundings( YZChange& chg, RationalMap& IncBeforeMap,
      RationalMap& IncAfterMap, int Cap, int SetDemand, int& BestNum,
      int& BestDenom )
{
   int i,j,q;
   LONG_DOUBLE LHS;
   double Viol,BestViol,Eps;
   double BestLHS = 0.0;
   double BestRHS = 0.0;
   int Num = 0;
   int Denom = 1;
   int RHS;

   RationalMap::iterator it, it2;

   Eps = COL_EXPAND_EPS;

   //fprintf(stderr,"BEGIN DoRoundings\n");
   //fprintf(stderr,"  Cap=%d, SetDemand=%d\n",Cap,SetDemand);
   //fprintf(stderr,"Eps = %lf\n",Eps);

   // update the contribution from the RHS through all multipliers
   for (i = 0; i < chg.prevSetDem; i++)
   {
      if (i * Cap > chg.prevSetDem * MAX_CAP_HECC) break;
      Rational r;
      if (i == 0)
      {
         r.num = 0;
         r.den = chg.prevSetDem;
      }
      else
      {
         r.num = i;
         r.den = chg.prevSetDem;
      }

      /* Update the j'th LHS increase */
      it = IncAfterMap.find(r);
      if (it == IncAfterMap.end())
         IncAfterMap[r] = 1.0;
      else
      {
         it->second += (LONG_DOUBLE)1.0;
         if (fabs(it->second) < Eps)
            IncAfterMap.erase(it);
      }
   }
   for (i = 0; i < SetDemand; i++)
   {
      if (i * Cap > SetDemand * MAX_CAP_HECC) break;
      Rational r;
      if (i == 0)
      {
         r.num = 0;
         r.den = SetDemand;
      }
      else
      {
         r.num = i;
         r.den = SetDemand;
      }

      /* Update the j'th LHS increase */
      it = IncAfterMap.find(r);
      if (it == IncAfterMap.end())
         IncAfterMap[r] = -1.0;
      else
      {
         it->second -= (LONG_DOUBLE)1.0;
         if (fabs(it->second) < Eps)
            IncAfterMap.erase(it);
      }
   }

   // update the LHS changes due to the y variables
   for (j = 0; j < (int)chg.demY.size(); j++)
   {
      q = chg.demY[j];

      // update the contribution from Y[q] through all multipliers
      for (i = 0; i < q; i++)
      {
         // obtain the corresponding multiplier
         if (i * Cap > q * MAX_CAP_HECC) break;
         Rational r;
         if (i == 0)
         {
            r.num = 0;
            r.den = q;
         }
         else
         {
            r.num = i;
            r.den = q;
         }

         // update the map of LHS increases
         it = IncAfterMap.find(r);
         if (it == IncAfterMap.end())
         {
            IncAfterMap[r] = chg.valY[j];
            //fprintf(stderr,"LHSInc(%d) += Y(%d) = %lf => LHSInc(%d) = %lf\n",i,q,
            //    chg.valY[j],i,IncAfterMap[r]);
         }
         else
         {
            it->second += (LONG_DOUBLE)chg.valY[j];
            //fprintf(stderr,"LHSInc(%d) += Y(%d) = %lf => LHSInc(%d) = %lf\n",i,q,
            //    chg.valY[j],i,IncAfterMap[r]);
            if (fabs(it->second) < Eps)
               IncAfterMap.erase(it);
         }
      }
   }

   // update the LHS changes due to the z variables
   for (j = 0; j < (int)chg.demZ.size(); j++)
   {
      q = chg.demZ[j];

      // update the contribution from Z[q] through all multipliers
      for (i = 1; i < q; i++)
      {
         // obtain the corresponding multiplier
         if (i * Cap > q * MAX_CAP_HECC) break;
         Rational r;
         r.num = i;
         r.den = q;

         // update the map of LHS increases
         it = IncBeforeMap.find(r);
         if (it == IncBeforeMap.end())
         {
            IncBeforeMap[r] = -chg.valZ[j];
            //fprintf(stderr,"LHSInc(%d) -= Z(%d) = %lf => LHSInc(%d) = %lf\n",j,q,
            //    -chg.valZ[j],j,IncBeforeMap[r]);
         }
         else
         {
            it->second -= (LONG_DOUBLE)chg.valZ[j];
            //fprintf(stderr,"LHSInc(%d) -= Z(%d) = %lf => LHSInc(%d) = %lf\n",j,q,
            //    -chg.valZ[j],j,IncBeforeMap[r]);
            if (fabs(it->second) < Eps)
               IncBeforeMap.erase(it);
         }
      }
   }

   /* Find the best multiplier */
   BestViol = -2.0*Cap;
   BestNum = 0;
   BestDenom = 1;

   LHS = 0.0;
   it = IncBeforeMap.begin();
   it2 = IncAfterMap.begin();
   bool changed = false;
   int prevNum = -1;
   int prevDen = 1;
   while ((it != IncBeforeMap.end()) || (it2 != IncAfterMap.end()))
   {
      bool incAfter = false;
      if (!changed)
      {
         if (it2 == IncAfterMap.end())
         {
            LHS += it->second;
            changed = true;
            Num = it->first.num;
            Denom = it->first.den;
            it++;
         }
         else if (it != IncBeforeMap.end())
         {
            long long m1 = (long long)(it->first.num) * (long long)(it2->first.den);
            long long m2 = (long long)(it2->first.num) * (long long)(it->first.den);
            if (m1 <= m2)
            {
               changed = true;
               Num = it->first.num;
               Denom = it->first.den;
               LHS += it->second;
               it++;
            }
            if (m2 <= m1)
               incAfter = true;
         }
         else
            incAfter = true;
      }

      if (changed)
      {
         RHS = ((long long)(SetDemand) * (long long)(Num)) / (long long)(Denom);
         if ((((long long)(SetDemand) * (long long)(Num)) % (long long)(Denom)) > 0) RHS++;

         if ((long long)Num * (long long)prevDen <=
               (long long)Denom * (long long)prevNum)
            fprintf( stderr, "ERROR IN THE RECC MULTIPLIERS (%d/%d <= %d/%d).\n",
                  Num, Denom, prevNum, prevDen );
         prevNum = Num;
         prevDen = Denom;

         //fprintf(stderr,"Dem=%d, r=%d/%d => LHS-RHS=%lf, RHS=%d\n",
         //      SetDemand,Num,Denom,LHS,RHS);

         Viol = -LHS;
         //if ((Viol/(1.0 * RHS)) > BestViol)
         //if ((Viol/(0.1 * SetDemand + 1.0 * RHS)) > BestViol)
         if (Viol > BestViol)
         //if ((Viol * RHS) > BestViol)
         {
            //BestViol = Viol/(1.0 * RHS);
            //BestViol = Viol/(0.1 * SetDemand + 1.0 * RHS);
            BestViol = Viol;
            //BestViol = Viol * RHS;
            BestNum = Num;
            BestDenom = Denom;

            BestLHS = LHS + (1.0 * RHS);
            BestRHS = 1.0 * RHS;

            // check if the vilation is acceptable
            if (Viol > 1.01)
            {
               fprintf(stderr,"Invalid violation = %lf:\n", BestViol);
               fprintf(stderr,"  Dem=%d, r=%d/%d => LHS=%lf, RHS=%g\n",
                     SetDemand,BestNum,BestDenom,BestLHS,BestRHS);
               fprintf(stderr,"next inc before = %lf: r=%d/%d\n",
                     it->second, it->first.num, it->first.den);
               fprintf(stderr,"next inc after = %lf: r=%d/%d\n",
                     it2->second, it2->first.num, it2->first.den);
               throw(-93772);
            }
         }

         changed = false;
      }

      if (incAfter)
      {
         changed = true;
         Num = it2->first.num;
         Denom = it2->first.den;
         LHS += it2->second;
         it2++;

         // add a sufficiently small fraction to the multiplier
         int maxDen = Cap;
         int maxDen2 = maxDen;
         if (SetDemand > maxDen2) maxDen2 = SetDemand;
         i = ((long long)(maxDen)*(long long)(maxDen2))/(long long)(Denom) + 1;
         Num = Num*i + 1;
         Denom *= i;
      }
   }

   //if (BestViol > 0.0)
   //   fprintf(stderr,"ECC: r=%d/%d, Dem=%d, LHS=%lf, RHS=%g, Viol=%lf\n",
   //         BestNum,BestDenom,SetDemand,BestLHS,BestRHS,BestViol);
}
