#include <fixedVector.h>
#include <fixedMatrix.h>
#include <subVector.h>
#include <sthead.h>

//#include <stream.h>
#include <cdsIostream.h>
#include <cdsIomanip.h>
#include <math.h> /* for sqrt */
//#include <cdsMath.h>

#include <stdlib.h>
#include <string.h>
#include <cdsSStream.h>

int
FixedVector_test()
{
 int exit=0;
 cout << "testing FixedVectorBase...";

 FixedVector<int,2> l;              //test all constructors
 FixedVector<int,2> l1( 9 );
 const int ar[2] = { 4, 5 };
 FixedVector<int,2> l2( ar );
 FixedVector<int,2> l3(l2);

 l = l2;

 if (l1(0) != 9 || l1(1) != 9) {
   cerr << "index failed: l1 != { 9, 9}\n";
   exit=1;
 }
 
 if (l2(0) !=4 || l2(1) !=5) {
   cerr << "index failed: l2 != { 4, 5}\n";
   exit=1;
 }

 if (l3(0) !=4 || l3(1) !=5) {
   cerr << "index failed: l3 != { 4, 5}\n";
   exit=1;
 }

 if (l(0) !=4 || l(1) !=5) {
   cerr << "index failed: l != { 4, 5}\n";
   exit=1;
 }

 l+=l1;

 if (l(0) !=13 || l(1) !=14) {
   cerr << "operator+= failed: l != { 13, 14}\n";
   exit=1;
 }

 l-=l2;

 if (l(0) !=9 || l(1) !=9) {
   cerr << "operator-= failed: l != { 9, 9}\n";
   exit=1;
 }

 l/=2;

 if (l(0) !=4 || l(1) !=4) {
   cerr << "operator/= failed: l != { 4, 4}\n";
   exit=1;
 }

 l*=3;

 if (l(0) !=12 || l(1) !=12) {
   cerr << "operator*= failed: l != { 12, 12}\n";
   exit=1;
 }

 int a[] = {1, 2};
 l.copyArray(a);

 if (l(0) !=1 || l(1) !=2) {
   cerr << "copyArray failed: l != { 1, 2}\n";
   exit=1;
 }

 OStringStream os; os << l << ends;
 if ( strcmp(os.str(),"{ 1, 2 }") ) {
   cerr << "operator<< error: " << os.str() << " != { 1, 2 }\n";
   exit=1;
 }

 {
   OStringStream os2; os2 << "{ 3, 4}" << ends;
   bool fail=0;
   try {
     istream is(os2.rdbuf()); is >> l;
   }
   catch (...) {
     fail=1;
   }
   if (fail || l(0) !=3 || l(1) !=4) {
     cerr << "operator>> error: " << os2.str() << " != { 3, 4 }\n";
     exit=1;
   }
 }

 FixedVector<int,2,1>& nl = (FixedVector<int,2,1>&)l;
 if (nl(1) !=3 || nl(2) !=4) {
   cerr << "error in conversion to diff. offset: " << nl << " != { 3, 4 }\n";
   exit=1;
 }

 // ADD: tests for blockVec and subVec(vec)
 {
   int a[] = { 1, 2, 
	       3, 4 };
   FixedMatrix<int,2,2> m(a);
   FixedVector<int,2,1> v = FixedVector<int,2,1>::subCol(m,0,0,1);
   if (v(1) != 1 || v(2) != 3) {
     cerr << "error in subVec(matrix): " << v << " != { 1, 3 }\n";
     exit=1;
   }
 }
   
 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* test */
