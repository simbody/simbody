/*
   CDSArray.C
   my own Array template class - same as CDSVector minus
   the mathematical operations
*/

#include "cdsMath.h"
#include "cdsVector.h"
#include "sthead.h"

#include <cdsIostream.h>
#include <cdsAlloc.h>

#include <stdlib.h>
#include <string.h>
#include <subVector.h>
#include <standAloneOps.h>

int
CDSVector_test()
{
 int exit=0;
 cout << "testing CDSVector...";


 CDSVector<double> v0(0); //zero sized vector
 if ( v0.size() != 0 ) {
   cerr << "error in zero-sized vector\n";
   exit=1;
 }


 CDSVector<double> v1;                 //test all constructors
 CDSVector<double> v2(2);   v2(0) = 1; v2(1) = 2;
 CDSVector<double,1> v3(1); v3(1) = 3;
 CDSVector<double> v4(v2);

 v4.resize(3); v4(2) = 4;
 
 if ( v3.size() != 1  || v3.offset() != 1 ) {
   cerr << "error in size or offset of v3\n";
   exit=1;
 }

 double a[] = { 5., 6. };
 v4.copy(a,2);

 if ( v4(0)!=5. || v4(1)!=6. ) {
   cerr << "error in copy: " << v4 << " != " << "{ 5, 6 }\n";
   exit=1;
 }

 v4 = v2;

 if ( v4(0)!=1. || v4(1)!=2. ) {
   cerr << "error in operator=: " << v4 << " != " << "{ 1, 2 }\n";
   exit=1;
 }

 {
   CDSVector<double> v5; v5 = v2.vector();

   if ( v5(0)!=1. || v5(1)!=2. ) {
     cerr << "error in operator=: " << v5 << " != " << "{ 1, 2 }\n";
     exit=1;
   }
 }

 v4.set(0.);

 if ( v4(0)!=0. || v4(1)!=0. ) {
   cerr << "error in set: " << v4 << " != " << "{ 0, 0 }\n";
   exit=1;
 }

 //
 // math operations
 //
 v2 *= 2.0;
 if ( v2(0)!=2 || v2(1)!=4 ) {
   cerr << "error in scale: " << v2 << " != " << "{ 2, 4 }\n";
   exit=1;
 }

 v2+=v2;
 if ( v2(0)!=4 || v2(1)!=8 ) {
   cerr << "error in operator+=: " << v2 << " != " << "{ 3, 5 }\n";
   exit=1;
 }

 v2-=v2;
 if ( v2(0)!=0 || v2(1)!=0 ) {
   cerr << "error in operator-=: " << v2 << " != " << "{ 0, 0 }\n";
   exit=1;
 }

 CDSVector<double,1> v5(3,4);
 if ( v5(1)!=4 || v5(2)!=4 || v5(3)!=4 ) {
   cerr << "error in constructor with initializer: " << v5 
	<< " != " << "{ 4, 4, 4 }\n";
   exit=1;
 }

 CDSVector<double,1> v6;
 if ( v6.size() != 0 ) {
   cerr << "error in default constructor: v6.size() != 0\n";
   exit=1;
 }

 CDSVector<double> bV = blockVec(v2,v5);
 if ( bV(0)!=0 || bV(1)!=0 || bV(2)!=4 || bV(3)!=4 || bV(4)!=4 ) {
   cerr << "error in blockVec: " << bV
	<< " != " << "{ 0, 0, 4, 4, 4 }\n";
   exit=1;
 }

 CDSVector<double,1> sV( SubVector< CDSVector<double> >(bV,1,2) );
 if ( sV(1)!=0 || sV(2)!=4 ) {
   cerr << "error in subVec: " << sV
	<< " != " << "{ 0, 4 }\n";
   exit=1;
 }

 v2.set(3.0);
 v4.set(1.0);
 v2 = v2 + 0.5*sV + v4;
 if ( v2(0)!=4 || v2(1)!=6 ) {
   cerr << "error in operator+ or operator*: " << v2
	<< " != " << "{ 4, 6 }\n";
   exit=1;
 }
 
   

 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* test */


