#ifndef INSTANCE_HPP_INCLUDED
#define INSTANCE_HPP_INCLUDED

#include <algorithm>

class Instance
{
public:
   Instance( const char *fileName );

   const char* getName() const { return &(name[0]); }

   int jobs() const { return n; }

   int machines() const { return m; }

   const int *ptime() const { return p; }

   // maximum time with some activity in some machine
   int T() const { return _T; }

   // cost of job j finishing at time t
   int getCost ( const int j, const int t ) const { return std::max( 0, w[j]*(t-d[j])); }

   int Tmin() const { return _Tmin; }

   int psum() const { return _psum; }

   int pmax() const { return _pmax; }

   int *weight() { return w; }

   int *duedate() { return d; }
   
   const char* getFileName() { return pfname; }

   virtual ~Instance();
private:

   const char* pfname;
   
   char name[256];

   int _psum;

   int _T;

   int _Tmin;

   int _pmax;

   int n;
   int m;

   int *p;
   int *w;
   int *d;
};

#endif // INSTANCE_HPP_INCLUDED
