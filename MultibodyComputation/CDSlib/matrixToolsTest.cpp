
#include "matrixTools.h"
#include <lapack.h>
#include <cdsAuto_arr.h>
#include <cdsExcept.h>
#include <cdsMath.h>
#include <cdsSStream.h>
#include <fixedSymMatrix.h>

#include "fixedVector.h"
#include "fixedMatrix.h"

using CDS::Complex;


static double
dot(const FixedVector<double,3>& v1,
    const CDSVector<double>& v2)
{  double ret=0; 
 for (int i=0;i<v1.size();i++) ret += v1(i)*v2(i);
 return ret;
} /* dot */

typedef Complex<double> DComplex;

static DComplex
dot(const FixedVector<DComplex,3>& v1,
    const CDSVector<DComplex>& v2)
{ 
 DComplex ret=0; 
 for (int i=0;i<v1.size();i++) ret += v1(i)*v2(i);
 return ret;
} /* dot */


int
MatrixTools::matrixToolsTest()
{
 int exit=0;
 double tol=1e-8;
 cout << "testing MatrixTools...";

 double ar[] = { 4, 2, 5, 7 }; 
 FixedMatrix<double,2,2> m(ar);

 FixedMatrix<double,2,2> t = transpose(m);
 if ( t(0,0) != 4 ||
      t(0,1) != 5 ||
      t(1,0) != 2 ||
      t(1,1) != 7  ) {
   cerr << "transpose failed: t != [ 4 5 ]\n"
	<< "                       [ 2 7 ]\n"
	<< " t = " << t << '\n';
   exit=1;
 }

 {
   FixedMatrix<double,2,2> mI = inverse(m);

   double ma[] = {  .3888888889,   -.1111111111,
		   -.2777777778,    .2222222222 };
   FixedMatrix<double,2,2> mIc(ma);
   
   for (int i=0 ; i<2 ; i++)
     for (int j=0 ; j<2 ; j++)
       if ( fabs(mI(i,j)-mIc(i,j))>tol ) {
	 cerr << "inverse error in mI(" << i << ',' << j << "): " << mI(i,j)
	      << " != " << mIc(i,j) << '\n';
	 exit=1;
       }
 }

 {
   double ar[] = { 4, 2, 5, 7, 3, 1 }; 
   FixedMatrix<double,2,3> m(ar);
   SVDResults< double > svdRet = svd(m);

   if ( svdRet.info != 0 ) {
     cerr << "svd info: " << svdRet.info << " != 0\n";
     exit=1;
   }
     

   if ( fabs(svdRet.sigma(0)-9.572002168)>tol ||
	fabs(svdRet.sigma(1)-3.518064026)>tol  ) {
     cerr << "svd sigma " << svdRet.sigma 
	  << "!= {9.572002168, 3.518064026}\n";
     exit=1;
   }
   
   double ua[] = {  -.6416135664,  -.7670280512,
		    -.7670280512,   .6416135664};
   FixedMatrix<double,2,2> u(ua);
   for (int i=0 ; i<2 ; i++)
     for (int j=0 ; j<2 ; j++)
       if ( fabs(u(i,j)-svdRet.u(i,j))>tol ) {
	 cerr << "svd error in u(" << i << ',' << j << "): " << u(i,j)
	      << " != " << svdRet.u(i,j) << '\n';
	 exit=1;
       }
   

   //note that the sign of the 3rd column is indeterminate.
   double va[] = {-.8290481432,    .4045357758,   -.3860440161 ,
		  -.3744578431,    .1110794442,    .9205664999 ,
		  -.4152836380,   -.9077511571,   -.05939138709};

   FixedMatrix<double,3,3> v(va);
   for (int i=0 ; i<3 ; i++)
     for (int j=0 ; j<3 ; j++)
       if ( fabs(v(i,j)-svdRet.vT(j,i))>tol ) {
	 cerr << "svd error in v(" << i << ',' << j << "): " << v(i,j)
	      << " != " << svdRet.vT(j,i) << '\n';
	 exit=1;
       }

   //CDSMatrix<double> sMat(2,3,0.);
   //sMat(0,0) = svdRet.sigma(0);
   //sMat(1,1) = svdRet.sigma(1);
   //CDSMatrix<double> result = svdRet.u * sMat * svdRet.vT;
   //cout << result << '\n';
       
 }

 {
   double ar[] = { 1.,
		   2, 4,
		   3, 5, 6 };
   FixedSymMatrix<double,3> m(ar);

   EigenResults<SymmetricMatrix<double> > ret = eigen(m);

   if ( ret.info != 0 ) {
     cerr << "sym eigen info: " << ret.info << " != 0\n";
     exit=1;
   }

   if ( (fabs(ret.eigenPairs[0].value+.51572947158925714020)>tol) ||
	(fabs(ret.eigenPairs[1].value-.17091518882717945220)>tol) ||
	(fabs(ret.eigenPairs[2].value-11.344814282762077688)>tol) ) {
     cerr << "eigen: eigenvalue error. eigenvalues: " 
	  << ret.eigenPairs[0].value << " "
	  << ret.eigenPairs[1].value << " "
	  << ret.eigenPairs[2].value << '\n';
     exit=1;
   }

	
   double ar0[] = {-.73697622909957824233,
		   -.32798527760568176797,  .59100904850610352543};
   FixedVector<double,3> v0(ar0);

   if ( fabs( (1-fabs(dot(v0,ret.eigenPairs[0].vector))))>tol ) {
     cerr << "eigen: eigenvector error for 0: " 
	  << ret.eigenPairs[0].vector << '\n';
     exit=1;
   }

   double ar1[] = {.59100904850610352539, -.73697622909957824241,
		   .32798527760568176779};
   FixedVector<double,3> v1(ar1);

   if ( fabs( (1-fabs(dot(v1,ret.eigenPairs[1].vector))))>tol ) {
     cerr << "eigen: eigenvector error for 1: " 
	  << ret.eigenPairs[1].vector << '\n';
     exit=1;
   }

   double ar2[] = {.32798527760568177039,   
		   .59100904850610352377,  .73697622909957824250};
   FixedVector<double,3> v2(ar2);

   if ( fabs( (1-fabs(dot(v2,ret.eigenPairs[2].vector))))>tol ) {
     cerr << "eigen: eigenvector error for 2: " 
	  << ret.eigenPairs[2].vector << '\n';
     exit=1;
   }



 }

 {
   double ar[] = { 4, 2, 5, 7 }; //constructor from row-major array
   FixedMatrix<double,2,2> fl2(ar);
   CDSMatrix<double,1,2> m(fl2);
   CDSVector<double> v = getColumn(m,2);
   if ( v(0) != 4 || v(1) != 5 ) {
     cerr << "getColumn failed: " << v << " != " << "{4,5}\n";
     exit=1;
   }
   v = getRow(m,2);
   if ( v(0) != 5 || v(1) != 7 ) {
     cerr << "getRow failed: " << v << " != " << "{5,7}\n";
     exit=1;
   }
   FixedVector<double,2,3> v2; v2(3) = 42; v2(4) = -3;
   setColumn(m,3,v2);
   if ( m(1,3) != 42 || m(2,3) != -3 ) {
     cerr << "setColumn failed: " << m << " != " << "{42,-3}\n";
     exit=1;
   }
   v2(3) = 51; v2(4) = 61;
   setRow(m,2,v2);
   if ( m(2,2) != 51 || m(2,3) != 61 ) {
     cerr << "setRow failed: " << m << " != " << "{51,61}\n";
     exit=1;
   }
 }


 using CDS::norm;
 using CDS::exp;
 using CDS::conj;
 {

   double ar[] = { 1, 2, 3,
		   4, 5, 6,
		   7, 8, 9 };
   FixedMatrix<double,3,3> m(ar);

   EigenResults<FullMatrix<double> > ret = eigen(m);

   if ( ret.info != 0 ) {
     cerr << "full eigen info: " << ret.info << " != 0\n";
     exit=1;
   }

   if ( norm(ret.eigenPairs[2].value)>tol ||
	norm(ret.eigenPairs[0].value-
	     Complex<double>(16.116843969807042990,0))>tol ||
	norm(ret.eigenPairs[1].value-
	     Complex<double>(-1.1168439698070429898,0))>tol ) {
     cerr << "full eigen: eigenvalue error. eigenvalues: " 
	  << ret.eigenPairs[0].value << " "
	  << ret.eigenPairs[1].value << " "
	  << ret.eigenPairs[2].value << '\n';
     exit=1;
   }

   FixedVector<DComplex,3> v0;
   v0[0] = DComplex(.23197068724628586166,0);
   v0[1] = DComplex(.52532209330123369344,0);
   v0[2] = DComplex(.81867349935618152522,0);

   if ( fabs( (1-norm(dot(v0,ret.eigenPairs[0].vector))))>tol ) {
     cerr << "eigen: eigenvector error for 0: " 
	  << ret.eigenPairs[0].vector << '\n';
     exit=1;
   }

   FixedVector<DComplex,3> v1;
   v1[0] = DComplex(.78583023874206710169,0);
   v1[1] = DComplex(.086751339256628453617,0);
   v1[2] = DComplex(-.61232756022881019446,0);

   if ( fabs( (1-norm(dot(v1,ret.eigenPairs[1].vector))))>tol ) {
     cerr << "eigen: eigenvector error for 1: " 
	  << ret.eigenPairs[1].vector << '\n';
     exit=1;
   }

   FixedVector<DComplex,3> v2;
   v2[0] = DComplex( .40824829046386301637);
   v2[1] = DComplex(-.81649658092772603274);
   v2[2] = DComplex( .40824829046386301637);

   if ( fabs( (1-norm(dot(v2,ret.eigenPairs[2].vector))))>tol ) {
     cerr << "eigen: eigenvector error for 2: " 
	  << ret.eigenPairs[2].vector << '\n';
     exit=1;
   }




 }
 {
   // test for complex eigenvalues
   double ar[] = 
     {.40450849718747371205 , .86602540378443864675 , .29389262614623656458,
      -.70062926922203677235 , .50000000000000000000 , -.50903696045512718345,
      -.58778525229247312917 , 0. , .80901699437494742410};
   FixedMatrix<double,3,3> m(ar);

   EigenResults<FullMatrix<double> > ret = eigen(m);

   if ( ret.info != 0 ) {
     cerr << "full eigen info: " << ret.info << " != 0\n";
     exit=1;
   }

   if ( fabs(ret.eigenPairs[2].value.re-1)>tol ||
	CDS::norm(ret.eigenPairs[0].value-
		  DComplex(.35676274578121056808,
			   .93419502419069398562))>tol ||
	norm(ret.eigenPairs[1].value-
	     DComplex(.35676274578121056808,
		      -.93419502419069398562))>tol ) {
     cerr << "full eigen: eigenvalue error. eigenvalues: " 
	  << ret.eigenPairs[0].value << " "
	  << ret.eigenPairs[1].value << " "
	  << ret.eigenPairs[2].value << '\n';
     exit=1;
   }

   FixedVector<DComplex,3> v0;
   v0[0] = DComplex(1.0590380375758614244, -.11709778340913245969);
   v0[1] = DComplex(-.041012847258213331484, .97547006019991125602);
   v0[2] = DComplex(.32102121300891695584, .51092514975162756950);

   
   FixedVector<DComplex,3> v0p = ret.eigenPairs[0].vector.vector();
   v0 *= (v0p[0]/v0[0]); //fix phase

   if ( norm(::norm(v0-v0p)) > tol ) {
     cerr << "eigen: eigenvector error for 0: " 
	  << ret.eigenPairs[0].vector << " != " 
	  << v0 << '\n';
     exit=1;
   }

   FixedVector<DComplex,3> v1;
   v1[0] = DComplex(1.0590380375758614244, .11709778340913245969);
   v1[1] = DComplex(-.041012847258213331484, -.97547006019991125602);
   v1[2] = DComplex(.32102121300891695584, -.51092514975162756950);
   //v1 *= exp( DComplex(0,atan2(-v1[0].im,v1[0].re)) ); //set phase
   //v1 /= sqrt(norm(v1*conj(v1.vector())));

   FixedVector<DComplex,3> v1p = ret.eigenPairs[1].vector.vector();
   v1 *= (v1p[0]/v1[0]); //fix phase

   if ( norm(::norm(v1-v1p)) > tol ) {
     cerr << "eigen: eigenvector error for 1: " 
	  << ret.eigenPairs[1].vector << " != " 
	  << v1 << '\n';
     exit=1;
   }

   FixedVector<DComplex,3> v2;
   v2[0] = DComplex(-.24835132946212412380);
   v2[1] = DComplex(-.43015712075567639685);
   v2[2] = DComplex( .76434679812116689667);
   v2 /= sqrt(norm(v2*conj(v2.vector())));

   if ( fabs( (1-norm(dot(v2,ret.eigenPairs[2].vector))))>tol ) {
     cerr << "eigen: eigenvector error for 2: " 
	  << ret.eigenPairs[2].vector << " != " 
	  << v2 << '\n';
     exit=1;
   }




 }


 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* matrixToolsTest */

