
#include <subMatrix.h>
#include <cdsExcept.h>

#include <stdlib.h>
#include <string.h>
#include <cdsSStream.h>
#include <fixedMatrix.h>

int
SubMatrix_test()
{
 int exit=0;
 cout << "testing SubMatrix...";

 int fMa[] = {1,2,3, 4,5,6, 7,8,9};
 FixedMatrix<int,3,3> fixedMat(fMa);

 SubMatrix< FixedMatrix<int,3,3> > sM(fixedMat,1,1,2,2);

 if ( sM(0,0)!=5 || sM(0,1)!=6 || sM(1,0)!=8 || sM(1,1)!=9 ) {
   cerr << "indexing failed.\n";
   exit=1;
 }

 SubMatrix< FixedMatrix<int,3,3> > sM2(fixedMat,0,1,2,2);
 sM2 = sM;

 if ( sM2(0,0)!=5 || sM2(0,1)!=6 || sM2(1,0)!=8 || sM2(1,1)!=9 ) {
   cerr << "operator=(submatrix) failed: " << sM2 << "\n";
   exit=1;
 }

 int fM2a[] = {4,3,2,1};
 FixedMatrix<int,2,2> fixedMat2(fM2a);
 sM2 = fixedMat2;
 
 if ( sM2(0,0)!=4 || sM2(0,1)!=3 || sM2(1,0)!=2 || sM2(1,1)!=1 ) {
   cerr << "operator=(matrix) failed.\n";
   exit=1;
 }

 SubMatrix< FixedMatrix<int,2,2> > sMp(fixedMat2,0,0,2,2);
 sMp += sMp;
 if ( sMp(0,0)!=8 || sMp(0,1)!=6 || sMp(1,0)!=4 || sMp(1,1)!=2 ) {
   cerr << "operator+= failed: " << sMp << "\n";
   exit=1;
 }

 sM(0,0) = 3; sM(0,1) = 4; sM(1,0) = 5; sM(1,1) = 6; 
 FixedMatrix<int,2,2> fM3 = sM * fixedMat2;
 
 if ( fM3(0,0)!=40 || fM3(0,1)!=26  || fM3(1,0)!=64 || fM3(1,1)!=42 ) {
   cerr << "operator*(submatrix,matrix) failed: " << fM3 << "\n";
   exit=1;
 }

 fM3 = fixedMat2 * sM;
 if ( fM3(0,0)!=54 || fM3(0,1)!=68 || fM3(1,0)!=22|| fM3(1,1)!=28) {
   cerr << "operator*(matrix,submatrix) failed.\n";
   exit=1;
 }

 fM3 = fixedMat2 + sM;
 if ( fM3(0,0)!=11 || fM3(0,1)!=10 || fM3(1,0)!=9|| fM3(1,1)!=8) {
   cerr << "operator+(matrix,submatrix) failed.\n";
   exit=1;
 }

 fM3 = sM + fixedMat2;
 if ( fM3(0,0)!=11 || fM3(0,1)!=10 || fM3(1,0)!=9|| fM3(1,1)!=8) {
   cerr << "operator+(submatrix,matrix) failed.\n";
   exit=1;
 }

 fM3 = sM - fixedMat2;
 if ( fM3(0,0)!=-5 || fM3(0,1)!=-2 || fM3(1,0)!=1 || fM3(1,1)!=4 ) {
   cerr << "operator-(submatrix,matrix) failed.\n";
   exit=1;
 }

 fM3 = fixedMat2 - sM;
 if ( fM3(0,0)!=5 || fM3(0,1)!=2 || fM3(1,0)!=-1 || fM3(1,1)!=-4 ) {
   cerr << "operator-(matrix,submatrix) failed.\n";
   exit=1;
 }

 {
   FixedVector<int,3> v; v(0)=5; v(1)=6; v(2)=7; 
   SubMatrix< FixedMatrix<int,3,3> > sM(fixedMat,0,1,3,1);
   sM.assignFromVector(v);
   
   for (int i=0 ; i<3 ; i++)
     if ( fixedMat(i,1) != v(i) ) {
       cerr << "assignFromVector (" << i << "): "
	    << fixedMat(i,1) << " != " << v(i) << '\n';
       exit=1;
     }
   //FixedVector<int,3> v2( CDS::Vector<CDSVector<int> >(sM.convertToVector()) );
   FixedVector<int,3> v2( sM.convertToVector().vector() );
   for (int i=0 ; i<3 ; i++)
     if ( fixedMat(i,1) != v2(i) ) {
       cerr << "convertToVector (" << i << "): "
	    << fixedMat(i,1) << " != " << v2(i) << '\n';
       exit=1;
     }

   SubMatrix< FixedMatrix<int,3,3> > sM2(fixedMat,2,0,1,3);
   sM2.assignFromVector(v);
   
   for (int i=0 ; i<3 ; i++)
     if ( fixedMat(2,i) != v(i) ) {
       cerr << "assignFromVector (" << i << "): "
	    << fixedMat(2,i) << " != " << v(i) << '\n';
       exit=1;
     }

   FixedVector<int,3> v3 = sM2.convertToVector().vector();
   for (int i=0 ; i<3 ; i++)
     if ( fixedMat(2,i) != v3(i) ) {
       cerr << "convertToVector (" << i << "): "
	    << fixedMat(2,i) << " != " << v3(i) << '\n';
       exit=1;
     }

 }

       

 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* test */

