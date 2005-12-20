

#include <fixedSymMatrix.h>
#include <fixedVector.h>
#include <cdsIostream.h>
#include <sthead.h>
#include <cdsMath.h>
#include <cdsString.h>

#include <stdlib.h>
#include <string.h>
#include <cdsSStream.h>
#include <matrixTools.h>

using MatrixTools::inverse;
//using CDS::orthoTransform;

int
FixedSymMatrix_test()
{
 int exit=0;
 cout << "testing FixedSymMatrix...";

 FixedSymMatrix<double,3> l;              //test all constructors
 FixedSymMatrix<double,2> l1(9);
 double ar[] = { 4, 
                 2, 5, 
		 6, 7, 8 }; // from row-major array (lower- tri)
 FixedSymMatrix<double,3> l2(ar);
 FixedSymMatrix<double,3> l3(l2);

 l = l2;

 if (l1(0,0) != 9 || l1(0,1) != 9 || 
     l1(1,0) != 9 || l1(1,1) != 9 ) {
   cerr << "index failed: l1 != { {9, 9},{9, 9} }\n";
   exit=1;
 }
 
 if (l2(0,0) !=4 || l2(0,1) !=2 || l2(0,2) !=6 ||
     l2(1,0) !=2 || l2(1,1) !=5 || l2(1,2) !=7 || 
     l2(2,0) !=6 || l2(2,1) !=7 || l2(2,2) != 8) {
   cerr << "index failed: l2 = " << l2 << '\n';
   exit=1;
 }

 if (l3(0,0) !=4 || l3(0,1) !=2 || l3(0,2) !=6 ||
     l3(1,0) !=2 || l3(1,1) !=5 || l3(1,2) !=7 ||
     l3(2,0) !=6 || l3(2,1) !=7 || l3(2,2) != 8) {
   cerr << "index failed: l3 != { {4, 2}, {5, 7}}\n";
   exit=1;
 }

 if (l(0,0) !=4 || l(0,1) !=2 || l(0,2) !=6 ||
     l(1,0) !=2 || l(1,1) !=5 || l(1,2) !=7 ||
     l(2,0) !=6 || l(2,1) !=7 || l(2,2) != 8) {
   cerr << "index failed: l != { {4, 2}, {5, 7}}\n";
   exit=1;
 }

 l3.set( 1 );
 l += l3;

 if (l(0,0) !=5 || l(0,1) !=3 || l(0,2) !=7 ||
     l(1,0) !=3 || l(1,1) !=6 || l(1,2) !=8 ||
     l(2,0) !=7 || l(2,1) !=8 || l(2,2) != 9) {
   cerr << "operator+= failed: l != { {5, 3}, {6, 8}}\n";
   exit=1;
 }

 l-=l2;

 if (l(0,0) !=1 || l(0,1) !=1 || l(0,2) !=1 ||
     l(1,0) !=1 || l(1,1) !=1 || l(1,2) !=1 ||
     l(2,0) !=1 || l(2,1) !=1 || l(2,2) != 1) {
   cerr << "operator-= failed: l != { {1, 1}, {1, 1}}\n";
   exit=1;
 }

 l(0,1) = 3;
 l.setDiag(2);

 if (l(0,0) !=2 || l(0,1) !=3 || l(0,2) !=1 ||
     l(1,0) !=3 || l(1,1) !=2 || l(1,2) !=1 ||
     l(2,0) !=1 || l(2,1) !=1 || l(2,2) != 2) {
   cerr << "setDiag(2) failed: " << l << " != { {2, 3}, {1, 2}}\n";
   exit=1;
 }

// l2 = inverse(l);
//
// if (l2(0,0) !=2 || l2(0,1) !=-3 ||
//     l2(1,0) !=-1 || l2(1,1) !=2) {
//   cerr << "inverse failed: " << l2 << "!= { {2, -3}, {-1, 2}}\n";
//   exit=1;
// }

 OStringStream os; os << l << ends;
 if ( strcmp(os.str(),"{ { 2, 3, 1 }, { 3, 2, 1 }, { 1, 1, 2 } }") ) {
   cerr << "operator<< error: " << os.str() 
	<< " != "
	<< "{ { 2, 3, 1 }, { 3, 2, 1 }, { 1, 1, 2 } }"
	<< '\n';
   exit=1;
 }

 // FIX:
// OStringStream os2; os2 << "{{1,2, 3},{ 2, 4, 5}, {3, 5, 7} }" << ends;
// istream is(os2.rdbuf()); is >> l;
// if (l(0,0) !=1 || l(0,1) !=2 || l(0,2) !=3 ||
//     l(1,0) !=2 || l(1,1) !=4 || l(1,2) !=5 ||
//     l(2,0) !=3 || l(2,1) !=5 || l(2,2) != 7) {
//   cerr << "operator>> error: " << os2.str() << " != {{1, 2},{ 3, 4 }}\n";
//   exit=1;
// }

// l3 = l * l2;
// if (l3(0,0) != 0 || l3(0,1) != 1 ||
//     l3(1,0) != 2 || l3(1,1) !=-1) {
//   cerr << "operator* failed: " << l3 << "!= { {0, 1}, {2, -1}}\n";
//   exit=1;
// }
//
// FixedVector<double,2> v1;
// FixedVector<double,2> v2(1.0); v2(1) = -1;
//
// v1 = l3 * v2;
// if (v1(0) != -1 || v1(1) != 3 ) {
//   cerr << "operator*(m,v) failed: " << v1 << "!= { -1,3 }\n";
//   exit=1;
// }

// FixedVector<double,2> v3 = v1 + v2;
// if (v3(0) != 0 || v3(1) != 2 ) {
//   cerr << "operator+(v,v) failed: " << v3 << "!= { 0,2 }\n";
//   exit=1;
// }

 

 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* test */

