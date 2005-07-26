

#include <fixedMatrix.h>
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

int
FixedMatrix_test()
{
 int exit=0;
 cout << "testing FixedMatrix...";

 FixedMatrix<double,2> l;              //test all constructors
 FixedMatrix<double,2,3> l1(9);
 double ar[] = { 4, 2, 5, 7 }; //constructor from row-major array
 FixedMatrix<double,2> l2(ar);
 FixedMatrix<double,2> l3(l2);

 l = l2;

 if (l1(0,0) != 9 || l1(0,1) != 9 || l1(0,2) != 9 || 
     l1(1,0) != 9 || l1(1,1) != 9 || l1(1,2) != 9 ) {
   cerr << "index failed: l1 != { {9, 9},{9, 9},{9, 9},}\n";
   exit=1;
 }
 
 if (l2(0,0) !=4 || l2(0,1) !=2 ||
     l2(1,0) !=5 || l2(1,1) !=7) {
   cerr << "index failed: l2 != { {4, 2}, {5, 7}}\n";
   exit=1;
 }

 if (l3(0,0) !=4 || l3(0,1) !=2 ||
     l3(1,0) !=5 || l3(1,1) !=7) {
   cerr << "index failed: l3 != { {4, 2}, {5, 7}}\n";
   exit=1;
 }

 if (l(0,0) !=4 || l(0,1) !=2 ||
     l(1,0) !=5 || l(1,1) !=7) {
   cerr << "index failed: l != { {4, 2}, {5, 7}}\n";
   exit=1;
 }

 l3.set( 1 );
 l += l3;

 if (l(0,0) !=5 || l(0,1) !=3 ||
     l(1,0) !=6 || l(1,1) !=8) {
   cerr << "operator+= failed: l != { {5, 3}, {6, 8}}\n";
   exit=1;
 }

 l-=l2;

 if (l(0,0) !=1 || l(0,1) !=1 ||
     l(1,0) !=1 || l(1,1) !=1) {
   cerr << "operator-= failed: l != { {1, 1}, {1, 1}}\n";
   exit=1;
 }

// // FixedMatrix<double,2,2,1,3>& nl = (FixedMatrix<double,2,2,1,3>&)l;
// FixedMatrix<double,2,2,1,3>& nl = l;
// if (nl(1,3) !=1 || nl(1,4) !=1 ||
//     nl(2,3) !=1 || nl(2,4) !=1) {
//   cerr << "error in conversion to diff. offset: " << nl 
//	  << " != {{1,1},{1,1}}\n";
//   exit=1;
// }

 l(0,1) = 3;
 l.setDiag(2);

 if (l(0,0) !=2 || l(0,1) !=3 ||
     l(1,0) !=1 || l(1,1) !=2) {
   cerr << "setDiag(2) failed: " << l << " != { {2, 3}, {1, 2}}\n";
   exit=1;
 }

 l2 = transpose(l);

 if (l2(0,0) !=2 || l2(0,1) !=1 ||
     l2(1,0) !=3 || l2(1,1) !=2) {
   cerr << "transpose failed: l != { {2, 1}, {3, 2}}\n";
   exit=1;
 }

 l2 = inverse(l);

 if (l2(0,0) !=2 || l2(0,1) !=-3 ||
     l2(1,0) !=-1 || l2(1,1) !=2) {
   cerr << "inverse failed: " << l2 << "!= { {2, -3}, {-1, 2}}\n";
   exit=1;
 }

 OStringStream os; os << l << ends;
 if ( strcmp(os.str(),"{ { 2, 3 }, { 1, 2 } }") ) {
   cerr << "operator<< error: " << os.str() << " != { { 2, 3 }, { 1, 2 } }\n";
   exit=1;
 }

 OStringStream os2; os2 << "{{1,2},{ 3, 4}}" << ends;
 istream is(os2.rdbuf()); is >> l;
 if (l(0,0) !=1 || l(0,1) !=2 ||
     l(1,0) !=3 || l(1,1) !=4) {
   cerr << "operator>> error: " << os2.str() << " != {{1, 2},{ 3, 4 }}\n";
   exit=1;
 }

 l3 = l * l2;
 if (l3(0,0) != 0 || l3(0,1) != 1 ||
     l3(1,0) != 2 || l3(1,1) !=-1) {
   cerr << "operator* failed: " << l3 << "!= { {0, 1}, {2, -1}}\n";
   exit=1;
 }

 FixedVector<double,2> v1;
 FixedVector<double,2> v2(1.0); v2(1) = -1;

 v1 = l3 * v2;
 if (v1(0) != -1 || v1(1) != 3 ) {
   cerr << "operator*(m,v) failed: " << v1 << "!= { -1,3 }\n";
   exit=1;
 }

 FixedVector<double,2> v3 = v1 + v2;
 if (v3(0) != 0 || v3(1) != 2 ) {
   cerr << "operator+(v,v) failed: " << v3 << "!= { 0,2 }\n";
   exit=1;
 }

 { // test rectangular matrices
   FixedMatrix<double,2,1> m1; m1(0,0) = 1; m1(1,0) = 2;
   
   FixedMatrix<double,1,2> m2 = transpose(m1);
   if (m2(0,0) != 1 || m2(0,1) != 2) {
     cerr << "transpose failed: " << m2 << "!= { 1, 2 }\n";
     exit=1;
   }
   
   FixedVector<double,2> v1; v1(0)=3 ; v1(1)=4;
   FixedVector<double,1> v2 = m2 * v1;   

   if ( v2(0) != 11 ) {
     cerr << "operator*(m,v) failed: " << v2 << "!= { 11 }\n";
     exit=1;
   }

   FixedMatrix<double,2,2> A(0.0);   // test blockMat
   FixedMatrix<double,2,1> B(1.0);
   FixedMatrix<double,1,2> C(2.0);
   FixedMatrix<double,1,1> D(3.0);
   FixedMatrix<double,3,3,1,1> M = blockMat22(A,B,
					      C,D);
   FixedMatrix<double,2,3,1,1> M2 = blockMat12(A,B);
   FixedMatrix<double,3,2,1,1> M3 = blockMat21(A,C);

   if ( M(1,1) != 0.0 || M(1,2) != 0.0 || M(1,3) != 1.0 ||
        M(2,1) != 0.0 || M(2,2) != 0.0 || M(2,3) != 1.0 ||
        M(3,1) != 2.0 || M(3,2) != 2.0 || M(3,3) != 3.0 ) {
     cerr << "blockMat22 failed: " << M << "!= "
	  << "{{ 0, 0, 1},{0,0,1},{2,2,3}}\n";
     exit=1;
   }

   if ( M2(1,1) != 0.0 || M2(1,2) != 0.0 || M2(1,3) != 1.0 ||
        M2(2,1) != 0.0 || M2(2,2) != 0.0 || M2(2,3) != 1.0  ) {
     cerr << "blockMat12 failed: " << M2 << "!= "
	  << "{{ 0, 0, 1},{0,0,1}}\n";
     exit=1;
   }
   
   if ( M3(1,1) != 0.0 || M3(1,2) != 0.0 || 
        M3(2,1) != 0.0 || M3(2,2) != 0.0 || 
        M3(3,1) != 2.0 || M3(3,2) != 2.0    ) {
     cerr << "blockMat21 failed: " << M3 << "!= "
	  << "{{ 0, 0},{0,0},{2,2}}\n";
     exit=1;
   }
 }
 { //test orthoTransform
   int a[] = {1, 2};       FixedMatrix<int,1,2> S(a);
   int b[] = {1, 2, 3, 4}; FixedMatrix<int,2,2> m(b);
   FixedMatrix<int,1,1> n = orthoTransform(m,S);
   if ( n(0,0) != 27 ) {
     cerr << "orthoTransform failed: " << n << "!= "
	  << "{{ 27 }}\n";
     exit=1;
   }
 }
 

 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* test */


