#include "Instance.hpp"

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#define STR_SIZE 256

#define ALLOW_LOW_MACHINE_LOAD

#ifdef _WIN32

char *strtok_r(char *s1, const char *s2, char **s3)
{
    char* s4;
    if (s1 != NULL) s4 = s1;
    else s4 = *s3;
    *s3 = strstr(s4, s2);
    if (*s3 != NULL)
    {
        **s3 = '\0';
        *s3 += strlen(s2);
    }
    return s4;
}

#endif


void readIntVector( FILE *file, int vector[], const int n );
char *getInstName(char *destination, const char *fileWithPath);
char *applyInversion(char* str);

Instance::Instance( const char *fileName ) :
   pfname(fileName),
   _psum(0),
   _pmax(0)
{
   char s1[STR_SIZE], s2[STR_SIZE], s3[STR_SIZE];
   char *ptrSave;

   getInstName( &(name[0]), fileName );
   //printf("instName: %s\n", &(instName[0]));

   // sample instance file name
   // wt40-2m-1.txt
   // wtJobs-machines-instanceNumber
   strcpy( &(s1[0]) , strtok_r( &(name[0]), "-", &ptrSave ) );
   strcpy( &(s2[0]) , strtok_r( NULL, "-", &ptrSave ) );
   strcpy( &(s3[0]) , strtok_r( NULL, "-", &ptrSave ) );

   // getting n. of jobs
   if (strlen(s1)<3)
   {
      fprintf( stderr, "Instance name should start with wt.\n" );
      exit(1);
   }

   this->n = atoi( (&(s1[0]) + 2) );

   if (strlen(s2)<2 )
   {
      fprintf( stderr, "Instance name should contain number of machines.\n" );
      exit(1);
   }

   s2[ strlen(s2)-1 ] = '\0';

   this->m = atoi( s2 );

   //printf(" %s %s %s \n", &(s1[0]), &(s2[0]), &(s3[0]) );

   this->p = new int[ n+1 ];
   this->w = new int[ n+1 ];
   this->d = new int[ n+1 ];

   FILE *f = fopen( fileName, "r" );
   if (!f)
   {
      fprintf( stderr, "Could not open instance file." );
      exit(1);
   }

   // first read must be discarded
   if (fscanf( f, "%d", &(p[0]) )!=1)
   {
      fprintf(stderr, "Wrong file format.\n");
      exit(1);
   }

   p[0] = 0;
   w[0] = 0;
   d[0] = 0;
   readIntVector( f, &p[1], n );
   readIntVector( f, &w[1], n );
   readIntVector( f, &d[1], n );

   for ( int i=1 ; ( i<=n ) ; i++ )
      _psum += p[i];

   for ( int i=1 ; ( i<=n ) ; i++ )
      if (p[i]>_pmax)
         _pmax = p[i];

   if ( m == 1 )
   {
      _T    = _psum;
      _Tmin = _psum;
   }
   else
   {
      _T = (int) (floor( ((double)(_psum-_pmax)) / ((double)m) )) + _pmax;
#ifdef ALLOW_LOW_MACHINE_LOAD
      if (n < 70)    // only for not too big instances
         _Tmin = 1;
      else
         _Tmin = (int) (ceil( ((double)(_psum-(m-1)*_pmax)) / ((double)m) ));
#else
      _Tmin = (int) (ceil( ((double)(_psum-(m-1)*_pmax)) / ((double)m) ));
#endif
   }

   fprintf( stderr, "T in [%d, %d]\n", _Tmin, _T );

   //for ( int i=1 ; i<=n ; i++ )
   //   printf( "d %d\n", (d[i]) );

   fclose( f );

   char nameDetails[128];
   sprintf( &(nameDetails[0]), "_%dm_%s", machines(), &(s3[0]) );
   strcat( &(name[0]), &(nameDetails[0]) );
   strcat( &(name[0]), ".txt" );

   fprintf( stderr,"Instance read. Dimensions: n: %d m: %d\n", this->n, this->m );
}

char *getInstName(char *destination, const char *fileWithPath) {
	int i,pos;
	// Returning till found a slash
	for ( i=(strlen(fileWithPath)-1),pos=0 ; (i>=0) ; i--,pos++ ) {
		if (fileWithPath[i]=='/') break;
		if (fileWithPath[i]=='\\') break;
		destination[pos] = fileWithPath[i];
	}
	destination[pos]='\0';
	applyInversion(destination);
	// Now, removing the .
	for ( i=0 ; i<(int)strlen(destination) ;i++) {
		if (destination[i]=='.') break;
	}
	destination[i]='\0';
	return destination;
}

char* applyInversion(char* str) {
	int head, tail;
	char t;
	for ( head=0,tail=(strlen(str)-1) ; head<tail ; head++,tail-- ) {
		t = str[head];
		str[head] = str[tail];
		str[tail] = t;
	}
	return str;
}

void readIntVector( FILE *file, int vector[], const int n )
{
   for ( int i=0 ; (i<n) ; i++ )
   {
      if (fscanf( file, "%d", &(vector[i]) ) != 1)
      {
         fprintf( stderr, "Could not read integer.\n");
         exit(1);
      }
   }
}

Instance::~Instance()
{
   delete[] p;
   delete[] w;
   delete[] d;
}
