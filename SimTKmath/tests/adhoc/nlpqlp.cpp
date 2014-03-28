

#include "SimTKmath.h"
using SimTK::Real;
using SimTK::Vector;
using SimTK::Matrix;
using SimTK::isNaN;
using SimTK::isFinite;
using SimTK::Optimizer;
using SimTK::OptimizerSystem;

#include <cstdlib>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <limits>

#include "nlpqlp.h"

static int pow_ii(int *ap, int *bp) {
    int pow, x, n;

    pow = 1;
    x = *ap;
    n = *bp;

    if(n < 0)
	    { }
    else if(n > 0)
	    for( ; ; )
		    {
		    if(n & 01)
			    pow *= x;
		    if(n >>= 1)
			    x *= x;
		    else
			    break;
		    }
    return(pow);
}

static const Real log10e = 0.43429448190325182765;
Real  d_lg10(Real *x)
{

return( log10e * log(*x) );
}
static Real d_sign(Real *a, Real *b)
{
   Real x;
   x = (*a >= 0 ? *a : - *a);
   return( *b >= 0 ? x : -x);
}
static Real pow_di(Real *ap, int *bp)
{
Real pow, x;
int n;

pow = 1;
x = *ap;
n = *bp;

if(n != 0)
	{
	if(n < 0)
		{
		n = -n;
		x = 1/x;
		}
	for( ; ; )
		{
		if(n & 01)
			pow *= x;
		if(n >>= 1)
			x *= x;
		else
			break;
		}
	}
return(pow);
}

int (*current_test)(int*);

typedef struct l20_values {
             bool lex;
        int nex;
        Real fex, *xex;
}*L20_vals; 

typedef struct l1_values {
        int n,nili,ninl,neli,nenl;
}*L1_vals; 

int tp1_(int*),tp2_(int*),tp3_(int*),tp4_(int*),tp5_(int*),tp6_(int*),tp7_(int*),tp8_(int*),tp9_(int*),tp10_(int*);
int tp11_(int*),tp12_(int*),tp13_(int*),tp14_(int*),tp15_(int*),tp16_(int*),tp17_(int*),tp18_(int*),tp19_(int*),tp20_(int*);
int tp21_(int*),tp22_(int*),tp23_(int*),tp24_(int*),tp25_(int*),tp26_(int*),tp27_(int*),tp28_(int*),tp29_(int*),tp30_(int*);
int tp31_(int*),tp32_(int*),tp33_(int*),tp34_(int*),tp35_(int*),tp36_(int*),tp37_(int*),tp38_(int*),tp39_(int*),tp40_(int*);
int tp41_(int*),tp42_(int*),tp43_(int*),tp44_(int*),tp45_(int*),tp46_(int*),tp47_(int*),tp48_(int*),tp49_(int*),tp50_(int*);
int tp51_(int*),tp52_(int*),tp53_(int*),tp54_(int*),tp55_(int*),tp56_(int*),tp57_(int*),tp58_(int*),tp59_(int*),tp60_(int*);
int tp61_(int*),tp62_(int*),tp63_(int*),tp64_(int*),tp65_(int*),tp66_(int*),tp67_(int*),tp68_(int*),tp69_(int*),tp70_(int*);
int tp71_(int*),tp72_(int*),tp73_(int*),tp74_(int*),tp75_(int*),tp76_(int*),tp77_(int*),tp78_(int*),tp79_(int*),tp80_(int*);
int tp81_(int*),tp82_(int*),tp83_(int*),tp84_(int*),tp85_(int*),tp86_(int*),tp87_(int*),tp88_(int*),tp89_(int*),tp90_(int*);
int tp91_(int*),tp92_(int*),tp93_(int*),tp94_(int*),tp95_(int*),tp96_(int*),tp97_(int*),tp98_(int*),tp99_(int*),tp100_(int*);
int tp101_(int*),tp102_(int*),tp103_(int*),tp104_(int*),tp105_(int*),tp106_(int*),tp107_(int*),tp108_(int*),tp109_(int*),tp110_(int*);
int tp111_(int*),tp112_(int*),tp113_(int*),tp114_(int*),             tp116_(int*),tp117_(int*),tp118_(int*),tp119_(int*);
int tp201_(int*),tp202_(int*),tp203_(int*),tp204_(int*),tp205_(int*),tp206_(int*),tp207_(int*),tp208_(int*),tp209_(int*),tp210_(int*);
int tp211_(int*),tp212_(int*),tp213_(int*),tp214_(int*),tp215_(int*),tp216_(int*),tp217_(int*),tp218_(int*),tp219_(int*),tp220_(int*);
int tp221_(int*),tp222_(int*),tp223_(int*),tp224_(int*),tp225_(int*),tp226_(int*),tp227_(int*),tp228_(int*),tp229_(int*),tp230_(int*);
int tp231_(int*),tp232_(int*),tp233_(int*),tp234_(int*),tp235_(int*),tp236_(int*),tp237_(int*),tp238_(int*),tp239_(int*),tp240_(int*);
int tp241_(int*),tp242_(int*),tp243_(int*),tp244_(int*),tp245_(int*),tp246_(int*),tp247_(int*),tp248_(int*),tp249_(int*),tp250_(int*);
int tp251_(int*),tp252_(int*),tp253_(int*),tp254_(int*),tp255_(int*),tp256_(int*),tp257_(int*),tp258_(int*),tp259_(int*),tp260_(int*);
int tp261_(int*),tp262_(int*),tp263_(int*),tp264_(int*),tp265_(int*),tp266_(int*),tp267_(int*),tp268_(int*),tp269_(int*),tp270_(int*);
int tp271_(int*),tp272_(int*),tp273_(int*),tp274_(int*),tp275_(int*),tp276_(int*),tp277_(int*),tp278_(int*),tp279_(int*),tp280_(int*);
int tp281_(int*),tp282_(int*),tp283_(int*),tp284_(int*),tp285_(int*),tp286_(int*),tp287_(int*),tp288_(int*),tp289_(int*),tp290_(int*);
int tp291_(int*),tp292_(int*),tp293_(int*),tp294_(int*),tp295_(int*),tp296_(int*),tp297_(int*),tp298_(int*),tp299_(int*),tp300_(int*);
int tp301_(int*),tp302_(int*),tp303_(int*),tp304_(int*),tp305_(int*),tp306_(int*),tp307_(int*),tp308_(int*),tp309_(int*),tp310_(int*);
int tp311_(int*),tp312_(int*),tp313_(int*),tp314_(int*),tp315_(int*),tp316_(int*),tp317_(int*),tp318_(int*),tp319_(int*),tp320_(int*);
int tp321_(int*),tp322_(int*),tp323_(int*),tp324_(int*),tp325_(int*),tp326_(int*),tp327_(int*),tp328_(int*),tp329_(int*),tp330_(int*);
int tp331_(int*),tp332_(int*),tp333_(int*),tp334_(int*),tp335_(int*),tp336_(int*),tp337_(int*),tp338_(int*),tp339_(int*),tp340_(int*);
int tp341_(int*),tp342_(int*),tp343_(int*),tp344_(int*),tp345_(int*),tp346_(int*),tp347_(int*),tp348_(int*),tp349_(int*),tp350_(int*);
int tp351_(int*),tp352_(int*),tp353_(int*),tp354_(int*),tp355_(int*),tp356_(int*),tp357_(int*),tp358_(int*),tp359_(int*),tp360_(int*);
int tp361_(int*),tp362_(int*),tp363_(int*);

static int (*(tests[]))(int*) = {  tp1_,  tp2_,  tp3_,  tp4_,  tp5_,  tp6_,  tp7_,  tp8_,  tp9_, tp10_, 
                           tp11_, tp12_, tp13_, tp14_, tp15_, tp16_, tp17_, tp18_, tp19_, tp20_,
                           tp21_, tp22_, tp23_, tp24_, tp25_, tp26_, tp27_, tp28_, tp29_, tp30_,
                           tp31_, tp32_, tp33_, tp34_, tp35_, tp36_, tp37_, tp38_, tp39_, tp40_, 
                           tp41_, tp42_, tp43_, tp44_, tp45_, tp46_, tp47_, tp48_, tp49_, tp50_,
                           tp51_, tp52_, tp53_, tp54_, tp55_, tp56_, tp57_, tp58_, tp59_, tp60_,
                           tp61_, tp62_, tp63_, tp64_, tp65_, tp66_, tp67_, tp68_, tp69_, tp70_, 
                           tp71_, tp72_, tp73_, tp74_, tp75_, tp76_, tp77_, tp78_, tp79_, tp80_, 
                           tp81_,        tp83_, tp84_, tp85_, tp86_, tp87_, tp88_, tp89_, tp90_, 
                           tp91_, tp92_, tp93_,        tp95_, tp96_, tp97_, tp98_, tp99_, tp100_, 
                           tp101_, tp102_, tp103_, tp104_, tp105_, tp106_, tp107_, tp108_, tp109_, tp110_,
                           tp111_, tp112_, tp113_, tp114_,         tp116_, tp117_, tp118_, tp119_, 
                           tp201_, tp202_, tp203_, tp204_, tp205_, tp206_, tp207_, tp208_, tp209_, tp210_,
                           tp211_, tp212_, tp213_, tp214_, tp215_, tp216_, tp217_, tp218_, tp219_, tp220_,
                           tp221_, tp222_, tp223_, tp224_, tp225_, tp226_, tp227_, tp228_, tp229_, tp230_,
                           tp231_, tp232_, tp233_, tp234_, tp235_, tp236_, tp237_, tp238_, tp239_, tp240_, 
                           tp241_, tp242_, tp243_, tp244_, tp245_, tp246_, tp247_, tp248_, tp249_, tp250_,
                           tp251_, tp252_, tp253_, tp254_, tp255_, tp256_, tp257_, tp258_, tp259_, tp260_,
                           tp261_, tp262_, tp263_, tp264_, tp265_, tp266_, tp267_, tp268_, tp269_, tp270_, 
                           tp271_, tp272_, tp273_, tp274_, tp275_, tp276_, tp277_, tp278_, tp279_, tp280_, 
                           tp281_, tp282_, tp283_, tp284_, tp285_, tp286_, tp287_, tp288_, tp289_, tp290_, 
                           tp291_, tp292_, tp293_, tp294_, tp295_, tp296_, tp297_, tp298_, tp299_, tp300_, 
                           tp301_, tp302_, tp303_, tp304_, tp305_, tp306_, tp307_, tp308_, tp309_, tp310_,
                           tp311_, tp312_, tp313_, tp314_, tp315_, tp316_, tp317_, tp318_, tp319_, tp320_,
                           tp321_, tp322_, tp323_, tp324_, tp325_, tp326_, tp327_, tp328_, tp329_, tp330_,
                           tp331_, tp332_, tp333_, tp334_, tp335_, tp336_, tp337_, tp338_, tp339_, tp340_, 
                           tp341_, tp342_, tp343_, tp344_, tp345_, tp346_, tp347_, tp348_, tp349_, tp350_,
                           tp351_, tp352_, tp353_, tp354_, tp355_, tp356_, tp357_, tp358_, tp359_, tp360_,
                           tp361_, tp362_, tp363_
                         };

static const char *test_names[] = {     "tp1_",  "tp2_",  "tp3_",  "tp4_",  "tp5_",  "tp6_",  "tp7_",  "tp8_",  "tp9_", "tp10_", 
                           "tp11_", "tp12_", "tp13_", "tp14_", "tp15_", "tp16_", "tp17_", "tp18_", "tp19_", "tp20_",
                           "tp21_", "tp22_", "tp23_", "tp24_", "tp25_", "tp26_", "tp27_", "tp28_", "tp29_", "tp30_",
                           "tp31_", "tp32_", "tp33_", "tp34_", "tp35_", "tp36_", "tp37_", "tp38_", "tp39_", "tp40_", 
                           "tp41_", "tp42_", "tp43_", "tp44_", "tp45_", "tp46_", "tp47_", "tp48_", "tp49_", "tp50_",
                           "tp51_", "tp52_", "tp53_", "tp54_", "tp55_", "tp56_", "tp57_", "tp58_", "tp59_", "tp60_",
                           "tp61_", "tp62_", "tp63_", "tp64_", "tp65_", "tp66_", "tp67_", "tp68_", "tp69_", "tp70_", 
                           "tp71_", "tp72_", "tp73_", "tp74_", "tp75_", "tp76_", "tp77_", "tp78_", "tp79_", "tp80_", 
                           "tp81_",        "tp83_", "tp84_", "tp85_", "tp86_", "tp87_", "tp88_", "tp89_", "tp90_", 
                           "tp91_", "tp92_", "tp93_",        "tp95_", "tp96_", "tp97_", "tp98_", "tp99_", "tp100_", 
                           "tp101_", "tp102_", "tp103_", "tp104_", "tp105_", "tp106_", "tp107_", "tp108_", "tp109_", "tp110_",
                           "tp111_", "tp112_", "tp113_", "tp114_",         "tp116_", "tp117_", "tp118_", "tp119_", 
                           "tp201_", "tp202_", "tp203_", "tp204_", "tp205_", "tp206_", "tp207_", "tp208_", "tp209_", "tp210_",
                           "tp211_", "tp212_", "tp213_", "tp214_", "tp215_", "tp216_", "tp217_", "tp218_", "tp219_", "tp220_",
                           "tp221_", "tp222_", "tp223_", "tp224_", "tp225_", "tp226_", "tp227_", "tp228_", "tp229_", "tp230_",
                           "tp231_", "tp232_", "tp233_", "tp234_", "tp235_", "tp236_", "tp237_", "tp238_", "tp239_", "tp240_", 
                           "tp241_", "tp242_", "tp243_", "tp244_", "tp245_", "tp246_", "tp247_", "tp248_", "tp249_", "tp250_",
                           "tp251_", "tp252_", "tp253_", "tp254_", "tp255_", "tp256_", "tp257_", "tp258_", "tp259_", "tp260_",
                           "tp261_", "tp262_", "tp263_", "tp264_", "tp265_", "tp266_", "tp267_", "tp268_", "tp269_", "tp270_", 
                           "tp271_", "tp272_", "tp273_", "tp274_", "tp275_", "tp276_", "tp277_", "tp278_", "tp279_", "tp280_", 
                           "tp281_", "tp282_", "tp283_", "tp284_", "tp285_", "tp286_", "tp287_", "tp288_", "tp289_", "tp290_", 
                           "tp291_", "tp292_", "tp293_", "tp294_", "tp295_", "tp296_", "tp297_", "tp298_", "tp299_", "tp300_", 
                           "tp301_", "tp302_", "tp303_", "tp304_", "tp305_", "tp306_", "tp307_", "tp308_", "tp309_", "tp310_",
                           "tp311_", "tp312_", "tp313_", "tp314_", "tp315_", "tp316_", "tp317_", "tp318_", "tp319_", "tp320_",
                           "tp321_", "tp322_", "tp323_", "tp324_", "tp325_", "tp326_", "tp327_", "tp328_", "tp329_", "tp330_",
                           "tp331_", "tp332_", "tp333_", "tp334_", "tp335_", "tp336_", "tp337_", "tp338_", "tp339_", "tp340_", 
                           "tp341_", "tp342_", "tp343_", "tp344_", "tp345_", "tp346_", "tp347_", "tp348_", "tp349_", "tp350_",
                           "tp351_", "tp352_", "tp353_", "tp354_", "tp355_", "tp356_", "tp357_", "tp358_", "tp359_", "tp360_",
                           "tp361_", "tp362_", "tp363_"
                         };
extern L1 l1_;
extern L2 l2_;
extern L3 l3_;
extern L4 l4_;
extern L5 l5_;
extern L6 l6_;
extern L9 l9_;
extern L10 l10_;
extern L11 l11_;
extern L12 l12_;
extern L13 l13_;
extern L14 l14_;
extern L20 l20_;
void init_common_blocks() {
  int *iptr,i;
  Real *rptr;
  bool *bptr;

  iptr = (int*)&l1_;  for(i=0;i<4;i++) *iptr++= std::numeric_limits<int>::quiet_NaN(); 
  rptr = (SimTK::Real*)&l2_;  for(i=0;i<100;i++) *rptr++= std::numeric_limits<double>::quiet_NaN();
  rptr = (SimTK::Real*)&l3_;  for(i=0;i<45;i++) *rptr++= std::numeric_limits<double>::quiet_NaN();
  rptr = (SimTK::Real*)&l4_;  for(i=0;i<100;i++) *rptr++= std::numeric_limits<double>::quiet_NaN();
  rptr = (SimTK::Real*)&l5_;  for(i=0;i<435;i++) *rptr++= std::numeric_limits<double>::quiet_NaN();
  rptr = (SimTK::Real*)&l6_;  for(i=0;i<1;i++) *rptr++= std::numeric_limits<double>::quiet_NaN();
  bptr = (bool*)&l9_;  for(i=0;i<45;i++) *bptr++= false;
  bptr = (bool*)&l10_;  for(i=0;i<38;i++) *bptr++= false;
  bptr = (bool*)&l11_;  for(i=0;i<100;i++) *bptr++= false;
  bptr = (bool*)&l12_;  for(i=0;i<100;i++) *bptr++= false;
  rptr = (SimTK::Real*)&l13_;  for(i=0;i<100;i++) *rptr++= std::numeric_limits<double>::quiet_NaN();
  rptr = (SimTK::Real*)&l14_;  for(i=0;i<100;i++) *rptr++= std::numeric_limits<double>::quiet_NaN();


  return;
}

/* Function Implementations */
class ProblemSystem : public OptimizerSystem {
public:

    int m;

    int objectiveFunc(  const Vector &coefficients, bool new_coefficients, Real& f ) const {
        int mode = 2;
        int i;

        Real *nx = (Real *)&l2_;
        for(i=0;i<getNumParameters();i++) nx[i] = coefficients(i);
        current_test(&mode);

        f = l6_.fx;
        if( !isFinite(f) ) {
            printf(" objective is not finite \n");
            exit(0);
        }
        /*
        printf("f = %f   x=",f);
        for(i=0;i<numParameters;i++)printf(" %f",coefficients(i)); printf("\n");
        */

        return (0);
    }

    int gradientFunc( const Vector &coefficients, bool new_coefficients, Vector &gradient ) const {
        int i,mode = 3;
        int fmode = 2;

        Real *nx = (Real *)&l2_;
        for(i=0;i<getNumParameters();i++) nx[i] = coefficients(i);

        current_test(&fmode);   // bug in tp260_ assumes objective has been calculated prior to gradient calc.
        current_test(&mode);

        Real *g = (Real*)&l4_;
        for(i=0;i<getNumParameters();i++)  {
            if( !isFinite(g[i]) ) {
                printf(" gradient is not finite i=%d g[i]=%f\n",i,g[i]);
                exit(0);
            }
            gradient(i) = g[i];
        }

        return(0);
    }

    int constraintFunc( const Vector &coefficients, bool new_coefficients, Vector &constraints)  const {
        int i,mode = 4;

        Real *nx = (Real *)&l2_;
        for(i=0;i<getNumParameters();i++) nx[i] = coefficients[i]; 
        current_test(&mode);

        Real *cv = (Real*)&l3_;
        for(i=0;i<m;i++)  {
            if( !isFinite(cv[i]) ) {
                printf(" constraint is not finite \n");
                exit(0);
            }
            constraints(i) = cv[i];
        }

        return(0);
    }

    int constraintJacobian( const Vector& coefficients, bool new_coefficients, Matrix& jac)  const {
        int i,j;
        int mode = 5;


        Real *nx = (Real *)&l2_;
        for(i=0;i<getNumParameters();i++) nx[i] = coefficients(i); 

        current_test(&mode);
        Real *jv = (Real*)&l5_;
        for(j=0;j<getNumEqualityConstraints();j++) {
            for(i=0;i<getNumParameters();i++) {
                if( !isFinite(jv[i*m+j]) ) {
                    printf(" jacobian(%d,%d) is not finite index=%d\n",j,i,i*getNumEqualityConstraints()+j);
                    exit(0);
                }
                jac(j,i) = jv[i*m+j];   // original
                //          jac(j,i) = jv[i+getNumParameters()*j];
            }
        }
        return (0);
    }

    ProblemSystem( const int numParams, const int numEqualityConstraints, const int numInEqualityConstraints)
    :   OptimizerSystem( numParams ) { 
        setNumEqualityConstraints( numEqualityConstraints );
        setNumInequalityConstraints( numInEqualityConstraints );

        m = numEqualityConstraints+ numInEqualityConstraints;
    } 

    ProblemSystem( const int numParams) : OptimizerSystem( numParams ) {}

}; 

/* Main Program */
int main()
{
  int j, num_tests=0,num_passed=0;
  bool run_test( int );
  void init_common_blocks();
  int failures [] = { 27, 74, 75, 76, 97, 99, 103, 107, 112, 114, 133, 135, 136, 
                      162, 164, 176, 182, 184,  197, 222, 223, 
                      224, 241, 248, 265,  277 };
   int fail64bit[] = { 76, 99, 103, 136, 164, 184 };


   for( j=0;j<279;j++) {
//   for( j=0;j<20;j++) {
//   for( j=83;j<86;j++) {
//   for( l=0;l<26;l++) {
//   for( l=0;l<1;l++) {
//      j = failures[l]-1;
   
      num_tests++;
      init_common_blocks();
      if( run_test( j ) ) num_passed++;
   }

   printf(" %d PASSED 	out of  %d tests   pass rate = %f %%\n",num_passed, num_tests, 100.0*(float)num_passed/(float)num_tests);
   return 0;
}

bool run_test(int j ) {

  L20_vals expect;
  L1_vals ndim;

  Real *initial_values,f;
  int mode = 1;
  int n,m,i;
  bool *index1, *index2;
  
  current_test = tests[j];
  current_test(&mode);

  ndim = (L1_vals)&l1_;
  n = ndim->n;
  m = ndim->nili + ndim->ninl + ndim->neli + ndim->nenl;

  ProblemSystem sys( n,  ndim->neli + ndim->nenl, ndim->nili + ndim->ninl );

  index1 = (bool *)&l9_;
  for(i=0;i<m;i++)  index1[i] = true;

  /* initialize the L10 common block for which constaint gradients are computed for mode 5. */
  index2 = (bool *)&l10_;
  for(i=0;i<m;i++)  index2[i] = true;


  Vector results(n);
  Vector lower_bounds(n);
  Vector upper_bounds(n);
  initial_values = (Real *)&l2_;
  for(i=0;i<n;i++) {
     results(i) = initial_values[i];
  }

  /* set the values for the variable bounds */

  bool *lxl = (bool *)&l11_;
  bool *lxu = (bool *)&l12_;
  Real *xl = (Real *)&l13_;
  Real *xu = (Real *)&l14_;
  bool has_limits = false;
  for (i=0; i<n; i++) {
     if( lxl[i] ) { 
        has_limits = true;
        lower_bounds(i) = xl[i];
     } else {
        lower_bounds(i) = -2e19;
     }
     if( lxu[i] ) { 
        has_limits = true;
        upper_bounds(i) = xu[i];
     } else {
        upper_bounds(i) = 2e19;
     }
  }

  if( has_limits ) {
     sys.setParameterLimits( lower_bounds, upper_bounds );
  }
  /* set the number of constraints and allocate space for the bounds */

  
  bool converged = true;
  try {

    Optimizer opt( sys );

//    opt.setConvergenceTolerance( .0001 );
    opt.setConvergenceTolerance( 1e-3 );

//    opt.useNumericalGradient( true );
    opt.useNumericalJacobian( true );
//    opt.setDiagnosticsLevel( 7 );

    // some tests do nor supply gradients
    if( j == 248-1 ||
        j == 264-1 ||
        j == 265-1 ||
        j == 272-1 ||
        j == 273-1 ||
        j == 277-1 ||
        j == 278-1 ||
        j == 279-1 ) {
         opt.useNumericalGradient( true );
         opt.useNumericalJacobian( true );
    }

    opt.setDiagnosticsLevel( 0 );

    printf("%d %s ",j+1, test_names[j]);
    /* compute  optimization */
    f = opt.optimize( results );
  }
  catch (const std::exception&) {
    converged = false;
  }

  if (converged ) {
    expect = (L20_vals)&l20_;
    printf("f(x*) = %e  expected = %e  diff=%e", f,expect->fex,f - expect->fex);
    if( fabs(f - expect->fex) < 1.0e-4 ) {
         printf("\n");
    } else {    
       printf(" *** %d \n",j+1);
    }
  } else {
    printf("#### Failed to converge #### \n" );
  }

  return( converged );
 } /* end of test */


/* prob.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

L1 l1_;
L2 l2_;
L3 l3_;
L4 l4_;
L5 l5_;
L6 l6_;
L9 l9_;
L10 l10_;
L11 l11_;
L12 l12_;
L13 l13_;
L14 l14_;
L20 l20_;


/* Table of constant values */

static Real c_b74 = .33333333333333331;
static Real c_b305 = 2.;
static Real c_b306 = -.33333333333333331;
static Real c_b308 = -.5;
static Real c_b310 = -.91666666666666663;
static Real c_b312 = -.25;
static Real c_b590 = .5;
static int c_n1 = -1;
static Real c_b933 = .66666666666666663;
static Real c_b934 = .25;
static Real c_b940 = 1.5;
static Real c_b949 = .75;
static Real c_b993 = .67;
static Real c_b997 = -.33;
static Real c_b998 = -.67;
static Real c_b1002 = -1.67;
static Real c_b1008 = -.71;
static Real c_b1009 = -1.3;
static Real c_b1015 = -1.71;
static Real c_b1016 = -2.3;
static Real c_b1157 = .2;
static Real c_b1523 = 3.;
static int c__1 = 1;
static int c__2 = 2;
static Real c_b1920 = -.66666666666666663;
static Real c_b2003 = -1.;
static Real c_b2368 = .35;
static Real c_b2380 = .33333333;
static Real c_b2383 = .782;
static Real c_b2387 = .546;
static Real c_b2390 = .3;
static Real c_b2391 = .467;
static Real c_b2746 = .666;
static Real c_b2748 = 1.6;
static Real c_b2752 = 3.55;
static Real c_b2753 = 5.66;
static Real c_b2754 = 1.2;
static Real c_b2757 = 2.42;
static Real c_b2758 = 4.17;
static Real c_b2761 = .803;
static Real c_b2782 = .334;
static Real c_b2789 = .6;
static Real c_b2800 = 4.66;
static Real c_b2801 = 2.55;
static Real c_b2813 = 1.42;
static Real c_b2816 = 3.17;
static Real c_b2823 = .197;
static Real c_b3046 = .827;
static Real c_b3047 = .182;
static Real c_b3048 = 10.;
static Real c_b3049 = 3.378;
static Real c_b3050 = .126;
static Real c_b3053 = 20.;
static Real c_b3054 = .656;

/* Initialized data */

static struct {
    Real e_1[2];
    Real fill_2[201];
    int e_3;
    } b_ = { .3, 7.263e-4, {0}, 98 };


/* Subroutine */ int tp1_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -2.;
    l2_1.x[1] = 1.;
    l11_1.lxl[0] = false;
    l11_1.lxl[1] = true;
    l12_1.lxu[0] = false;
    l12_1.lxu[1] = false;
    l13_1.xl[1] = -1.5;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 1.;
    l20_1.fex = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = d__1 * d__1 * 100. + d__3 * d__3;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (l2_1.x[1] - d__1 * d__1) * 200.;
    l4_1.gf[0] = (l2_1.x[0] * (l4_1.gf[1] - 1.) + 1.) * -2.;
labelL4:
    return 0;
} /* tp1_ */


/* Subroutine */ int tp2_(int *mode)
{
    /* System generated locals */
    Real r__1;
    Real d__1, d__2, d__3;

    /* Local variables */
    static Real w1;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -2.;
    l2_1.x[1] = 1.;
    l11_1.lxl[0] = false;
    l11_1.lxl[1] = true;
    l12_1.lxu[0] = false;
    l12_1.lxu[1] = false;
    l13_1.xl[1] = 1.5;
    l20_1.lex = true;
    l20_1.nex = 1;
    w1 = std::sqrt(.49833333333333335);
/* Computing 3rd power */
    r__1 = w1;
    l20_1.xex[0] = w1 * 2. * std::cos(std::acos(.0025 / (r__1 * (r__1 * r__1))) / 3.);
    l20_1.xex[1] = 1.5;
/* Computing 2nd power */
    d__2 = l20_1.xex[0];
/* Computing 2nd power */
    d__1 = l20_1.xex[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l20_1.xex[0];
    l20_1.fex = d__1 * d__1 * 100. + d__3 * d__3;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = d__1 * d__1 * 100. + d__3 * d__3;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (l2_1.x[1] - d__1 * d__1) * 200.;
    l4_1.gf[0] = (l2_1.x[0] * (l4_1.gf[1] - 1.) + 1.) * -2.;
labelL4:
    return 0;
} /* tp2_ */


/* Subroutine */ int tp3_(int *mode)
{
    /* System generated locals */
    Real d__1;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 10.;
    l2_1.x[1] = 1.;
    l11_1.lxl[0] = false;
    l11_1.lxl[1] = true;
    l12_1.lxu[0] = false;
    l12_1.lxu[1] = false;
    l13_1.xl[1] = 0.;
    l20_1.lex = true;
    l20_1.xex[0] = 0.;
    l20_1.xex[1] = 0.;
    l20_1.fex = 0.;
    l20_1.nex = 1;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[1] - l2_1.x[0];
    l6_1.fx = l2_1.x[1] + d__1 * d__1 * 1e-5;
    return 0;
labelL3:
    l4_1.gf[0] = (l2_1.x[1] - l2_1.x[0]) * -2. * 1e-5;
    l4_1.gf[1] = 1. - l4_1.gf[0];
labelL4:
    return 0;
} /* tp3_ */


/* Subroutine */ int tp4_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 1.125;
    l2_1.x[1] = .125;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = true;
    }
    l13_1.xl[0] = 1.;
    l13_1.xl[1] = 0.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 0.;
    l20_1.fex = 2.6666666666666665;
    l4_1.gf[1] = 1.;
    return 0;
labelL2:
/* Computing 3rd power */
    d__1 = l2_1.x[0] + 1.;
    l6_1.fx = d__1 * (d__1 * d__1) / 3. + l2_1.x[1];
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0] + 1.;
    l4_1.gf[0] = d__1 * d__1;
labelL4:
    return 0;
} /* tp4_ */


/* Subroutine */ int tp5_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real a;
    static int i__;
    static Real v1, v2;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 0.;
    l2_1.x[1] = 0.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = true;
/* labelL6: */
	l12_1.lxu[i__ - 1] = true;
    }
    l13_1.xl[0] = -1.5;
    l13_1.xl[1] = -3.;
    l14_1.xu[0] = 4.;
    l14_1.xu[1] = 3.;
    a = std::atan(1.) * 4.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = .5 - a / 3.;
    l20_1.xex[1] = l20_1.xex[0] - 1.;
    l20_1.fex = -std::sqrt(3.) / 2. - a / 3.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - l2_1.x[1];
    l6_1.fx = std::sin(l2_1.x[0] + l2_1.x[1]) + d__1 * d__1 - l2_1.x[0] * 1.5 + 
	    l2_1.x[1] * 2.5 + 1.;
    return 0;
labelL3:
    v1 = std::cos(l2_1.x[0] + l2_1.x[1]);
    v2 = (l2_1.x[0] - l2_1.x[1]) * 2.;
    l4_1.gf[0] = v1 + v2 - 1.5;
    l4_1.gf[1] = v1 - v2 + 2.5;
labelL4:
    return 0;
} /* tp5_ */


/* Subroutine */ int tp6_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    l2_1.x[0] = -1.2;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 1.;
    l20_1.fex = 0.;
    l5_1.gg[1] = 10.;
    l4_1.gf[1] = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = 1. - l2_1.x[0];
    l6_1.fx = d__1 * d__1;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2. - 2.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_1.g[0] = (l2_1.x[1] - d__1 * d__1) * 10.;
    }
    return 0;
labelL5:
    if (l10_2.index2[0]) {
	l5_1.gg[0] = l2_1.x[0] * -20.;
    }
    return 0;
} /* tp6_ */


/* Subroutine */ int tp7_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    l2_1.x[0] = 2.;
    l2_1.x[1] = 2.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.fex = -std::sqrt(3.);
    l20_1.xex[0] = 0.;
    l20_1.xex[1] = -l20_1.fex;
    l20_1.nex = 1;
    l4_1.gf[1] = -1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l6_1.fx = std::log(d__1 * d__1 + 1.) - l2_1.x[1];
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[0] = l2_1.x[0] * 2. / (d__1 * d__1 + 1.);
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__2 = l2_1.x[0];
/* Computing 2nd power */
	d__1 = d__2 * d__2 + 1.;
/* Computing 2nd power */
	d__3 = l2_1.x[1];
	l3_1.g[0] = d__1 * d__1 + d__3 * d__3 - 4.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l5_1.gg[0] = l2_1.x[0] * 4. * (d__1 * d__1 + 1.);
    l5_1.gg[1] = l2_1.x[1] * 2.;
L7:
    return 0;
} /* tp7_ */


/* Subroutine */ int tp8_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real a, b;
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL3;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 2;
    l2_1.x[0] = 2.;
    l2_1.x[1] = 1.;
    a = std::sqrt((std::sqrt(301.) + 25.) / 2.);
    b = std::sqrt((25. - std::sqrt(301.)) / 2.);
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_2.lex = true;
    l20_2.nex = 4;
    l20_2.xex[0] = a;
    l20_2.xex[1] = 9. / a;
    l20_2.xex[4] = b;
    l20_2.xex[5] = 9. / b;
    for (i__ = 3; i__ <= 7; i__ += 4) {
	for (j = 1; j <= 2; ++j) {
/* L30: */
	    l20_2.xex[i__ + j - 2] = -l20_2.xex[i__ + j - 4];
	}
    }
    l20_2.fex = -1.;
    l4_1.gf[0] = 0.;
    l4_1.gf[1] = 0.;
    l6_1.fx = -1.;
labelL3:
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_2.g[0] = d__1 * d__1 + d__2 * d__2 - 25.;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = l2_1.x[0] * l2_1.x[1] - 9.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_2.gg[0] = l2_1.x[0] * 2.;
    l5_2.gg[2] = l2_1.x[1] * 2.;
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    l5_2.gg[1] = l2_1.x[1];
    l5_2.gg[3] = l2_1.x[0];
L8:
    return 0;
} /* tp8_ */


/* Subroutine */ int tp9_(int *mode)
{
    /* Local variables */
    static int i__;
    static Real v, v1, v2, v3, v4;

    v = std::atan(1.) * 4.;
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    l2_1.x[0] = 0.;
    l2_1.x[1] = 0.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = -1;
    l20_1.fex = -.5;
    l20_1.xex[0] = -3.;
    l20_1.xex[1] = -4.;
    l5_1.gg[0] = 4.;
    l5_1.gg[1] = -3.;
    return 0;
labelL2:
    l6_1.fx = std::sin(v * l2_1.x[0] / 12.) * std::cos(v * l2_1.x[1] / 16.);
    return 0;
labelL3:
    v3 = v / 12.;
    v4 = v / 16.;
    v1 = v3 * l2_1.x[0];
    v2 = v4 * l2_1.x[1];
    l4_1.gf[0] = v3 * std::cos(v1) * std::cos(v2);
    l4_1.gf[1] = -v4 * std::sin(v1) * std:: sin(v2);
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = l2_1.x[0] * 4. - l2_1.x[1] * 3.;
    }
labelL5:
    return 0;
} /* tp9_ */


/* Subroutine */ int tp10_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -10.;
    l2_1.x[1] = 10.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = 0.;
    l20_1.xex[1] = 1.;
    l20_1.fex = -1.;
    l4_1.gf[0] = 1.;
    l4_1.gf[1] = -1.;
    return 0;
labelL2:
    l6_1.fx = l2_1.x[0] - l2_1.x[1];
labelL3:
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_1.g[0] = d__1 * d__1 * -3. + l2_1.x[0] * 2. * l2_1.x[1] - d__2 * 
		d__2 + 1.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_1.gg[0] = l2_1.x[0] * -6. + l2_1.x[1] * 2.;
    l5_1.gg[1] = (l2_1.x[0] - l2_1.x[1]) * 2.;
L7:
    return 0;
} /* tp10_ */


/* Subroutine */ int tp11_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;
    static Real aw, aex, qaw;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l2_1.x[0] = 4.9;
    l2_1.x[1] = .1;
    l20_1.lex = true;
    l20_1.nex = 1;
    aex = std::sqrt(6.) * 7.5;
/* Computing 2nd power */
    d__2 = aex;
    d__1 = std::sqrt(d__2 * d__2 + 1.) + aex;
    aw = pow_dd(&d__1, &c_b74);
/* Computing 2nd power */
    d__1 = aw;
    qaw = d__1 * d__1;
    l20_1.xex[0] = (aw - 1. / aw) / std::sqrt(6.);
    l20_1.xex[1] = (qaw - 2. + 1. / qaw) / 6.;
/* Computing 2nd power */
    d__1 = l20_1.xex[0] - 5.;
/* Computing 2nd power */
    d__2 = l20_1.xex[1];
    l20_1.fex = d__1 * d__1 + d__2 * d__2 - 25.;
    l5_1.gg[1] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 5.;
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 - 25.;
    return 0;
labelL3:
    l4_1.gf[0] = (l2_1.x[0] - 5.) * 2.;
    l4_1.gf[1] = l2_1.x[1] * 2.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_1.g[0] = -(d__1 * d__1) + l2_1.x[1];
    }
    return 0;
labelL5:
    if (l10_2.index2[0]) {
	l5_1.gg[0] = l2_1.x[0] * -2.;
    }
    return 0;
} /* tp11_ */


/* Subroutine */ int tp12_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 0.;
    l2_1.x[1] = 0.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = 2.;
    l20_1.xex[1] = 3.;
    l20_1.fex = -30.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l6_1.fx = d__1 * d__1 * .5 + d__2 * d__2 - l2_1.x[0] * l2_1.x[1] - l2_1.x[
	    0] * 7. - l2_1.x[1] * 7.;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] - l2_1.x[1] - 7.;
    l4_1.gf[1] = l2_1.x[1] * 2. - l2_1.x[0] - 7.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_1.g[0] = 25. - d__1 * d__1 * 4. - d__2 * d__2;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_1.gg[0] = l2_1.x[0] * -8.;
    l5_1.gg[1] = l2_1.x[1] * -2.;
L7:
    return 0;
} /* tp12_ */


/* Subroutine */ int tp13_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -2.;
    l2_1.x[1] = -2.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
	l11_1.lxl[i__ - 1] = true;
/* labelL6: */
	l13_1.xl[i__ - 1] = 0.;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 0.;
    l20_1.fex = 1.;
    l5_1.gg[1] = -1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 2.;
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = (l2_1.x[0] - 2.) * 2.;
    l4_1.gf[1] = l2_1.x[1] * 2.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 3rd power */
	d__1 = 1. - l2_1.x[0];
	l3_1.g[0] = d__1 * (d__1 * d__1) - l2_1.x[1];
    }
    return 0;
labelL5:
    if (l10_2.index2[0]) {
/* Computing 2nd power */
	d__1 = 1. - l2_1.x[0];
	l5_1.gg[0] = d__1 * d__1 * -3.;
    }
    return 0;
} /* tp13_ */


/* Subroutine */ int tp14_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;
    static Real w7;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    l2_1.x[0] = 2.;
    l2_1.x[1] = 2.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    w7 = std::sqrt(7.);
    l20_1.xex[0] = (w7 - 1.) * .5;
    l20_1.xex[1] = (w7 + 1.) * .25;
    l20_1.fex = 9. - w7 * 23. / 8.;
    l5_2.gg[1] = 1.;
    l5_2.gg[3] = -2.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 2.;
/* Computing 2nd power */
    d__2 = l2_1.x[1] - 1.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = (l2_1.x[0] - 2.) * 2.;
    l4_1.gf[1] = (l2_1.x[1] - 1.) * 2.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_2.g[0] = 1. - d__1 * d__1 * .25 - d__2 * d__2;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = l2_1.x[0] - l2_1.x[1] * 2. + 1.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_2.gg[0] = -l2_1.x[0] * .5;
    l5_2.gg[2] = l2_1.x[1] * -2.;
L7:
    return 0;
} /* tp14_ */


/* Subroutine */ int tp15_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -2.;
    l2_1.x[1] = 1.;
    l11_1.lxl[0] = false;
    l11_1.lxl[1] = false;
    l12_1.lxu[0] = true;
    l12_1.lxu[1] = false;
    l14_1.xu[0] = .5;
    l20_1.lex = true;
    l20_1.xex[0] = .5;
    l20_1.xex[1] = (float)2.;
    l20_1.fex = (float)3.065;
    l20_1.nex = 1;
    l5_2.gg[1] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = d__1 * d__1 + d__3 * d__3 * (float).01;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (l2_1.x[1] - d__1 * d__1) * (float)2.;
    l4_1.gf[0] = (l2_1.x[0] * (l4_1.gf[1] - 1.) + 1.) * -.02;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_1.x[0] * l2_1.x[1] - 1.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[1];
	l3_2.g[1] = d__1 * d__1 + l2_1.x[0];
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_2.gg[0] = l2_1.x[1];
    l5_2.gg[2] = l2_1.x[0];
L7:
    if (l10_3.index2[1]) {
	l5_2.gg[3] = l2_1.x[1] * 2.;
    }
    return 0;
} /* tp15_ */


/* Subroutine */ int tp16_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -2.;
    l2_1.x[1] = 1.;
    l11_1.lxl[0] = true;
    l11_1.lxl[1] = false;
    l12_1.lxu[0] = true;
    l12_1.lxu[1] = true;
    l13_1.xl[0] = -2.;
    l14_1.xu[0] = .5;
    l14_1.xu[1] = 1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = .5;
    l20_1.xex[1] = .25;
    l20_1.fex = .25;
    l5_2.gg[0] = 1.;
    l5_2.gg[3] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = d__1 * d__1 * 100. + d__3 * d__3;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (l2_1.x[1] - d__1 * d__1) * 200.;
    l4_1.gf[0] = (l2_1.x[0] * (l4_1.gf[1] - 1.) + 1.) * -2.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[1];
	l3_2.g[0] = d__1 * d__1 + l2_1.x[0];
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_2.g[1] = d__1 * d__1 + l2_1.x[1];
    }
    return 0;
labelL5:
    if (l10_3.index2[0]) {
	l5_2.gg[2] = l2_1.x[1] * 2.;
    }
    if (l10_3.index2[1]) {
	l5_2.gg[1] = l2_1.x[0] * 2.;
    }
    return 0;
} /* tp16_ */


/* Subroutine */ int tp17_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = (float)-2.;
    l2_1.x[1] = 1.;
    l11_1.lxl[0] = true;
    l11_1.lxl[1] = false;
    l12_1.lxu[0] = true;
    l12_1.lxu[1] = true;
    l13_1.xl[0] = (float)-2.;
    l14_1.xu[0] = .5;
    l14_1.xu[1] = 1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = 0.;
    l20_1.xex[1] = 0.;
    l20_1.fex = .01;
    l5_2.gg[0] = -1.;
    l5_2.gg[3] = -1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = (d__1 * d__1 * 100. + d__3 * d__3) * (float).01;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (l2_1.x[1] - d__1 * d__1) * 200. * (float).01;
    l4_1.gf[0] = (l2_1.x[0] * (l4_1.gf[1] - 1.) + 1.) * -2. * (float).01;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[1];
	l3_2.g[0] = d__1 * d__1 - l2_1.x[0];
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_2.g[1] = d__1 * d__1 - l2_1.x[1];
    }
    return 0;
labelL5:
    if (l10_3.index2[0]) {
	l5_2.gg[2] = l2_1.x[1] * 2.;
    }
    if (l10_3.index2[1]) {
	l5_2.gg[1] = l2_1.x[0] * 2.;
    }
    return 0;
} /* tp17_ */


/* Subroutine */ int tp18_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 2.;
    l2_1.x[1] = 2.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = true;
	l12_1.lxu[i__ - 1] = true;
/* labelL6: */
	l14_1.xu[i__ - 1] = 50.;
    }
    l13_1.xl[0] = 2.;
    l13_1.xl[1] = 0.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = std::sqrt(250.);
    l20_1.xex[1] = l20_1.xex[0] * .1;
    l20_1.fex = 5.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l6_1.fx = d__1 * d__1 * .01 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * .02;
    l4_1.gf[1] = l2_1.x[1] * 2.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_1.x[0] * l2_1.x[1] - 25.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_2.g[1] = d__1 * d__1 + d__2 * d__2 - 25.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_2.gg[0] = l2_1.x[1];
    l5_2.gg[2] = l2_1.x[0];
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    l5_2.gg[1] = l2_1.x[0] * 2.;
    l5_2.gg[3] = l2_1.x[1] * 2.;
L8:
    return 0;
} /* tp18_ */


/* Subroutine */ int tp19_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real saex;
    static int i__;
    static Real aex;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 20.1;
    l2_1.x[1] = 5.84;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = true;
	l12_1.lxu[i__ - 1] = true;
/* labelL6: */
	l14_1.xu[i__ - 1] = 100.;
    }
    l13_1.xl[0] = 13.;
    l13_1.xl[1] = 0.;
    l20_1.lex = true;
    l20_1.nex = 1;
    saex = 17.280975;
    aex = std::sqrt(saex);
    l20_1.xex[0] = 14.095;
    l20_1.xex[1] = 5. - aex;
/* Computing 3rd power */
    d__1 = aex + 15.;
    l20_1.fex = (68.669157374999998 - d__1 * (d__1 * d__1)) * (float)1e-4;
    return 0;
labelL2:
/* Computing 3rd power */
    d__1 = l2_1.x[0] - 10.;
/* Computing 3rd power */
    d__2 = l2_1.x[1] - 20.;
    l6_1.fx = (d__1 * (d__1 * d__1) + d__2 * (d__2 * d__2)) * (float)1e-4;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 10.;
    l4_1.gf[0] = d__1 * d__1 * 3. * (float)1e-4;
/* Computing 2nd power */
    d__1 = l2_1.x[1] - 20.;
    l4_1.gf[1] = d__1 * d__1 * 3. * (float)1e-4;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0] - 5.;
/* Computing 2nd power */
	d__2 = l2_1.x[1] - 5.;
	l3_2.g[0] = d__1 * d__1 + d__2 * d__2 - 100.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0] - 6.;
/* Computing 2nd power */
	d__2 = l2_1.x[1] - 5.;
	l3_2.g[1] = 82.81 - d__1 * d__1 - d__2 * d__2;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_2.gg[0] = (l2_1.x[0] - 5.) * 2.;
    l5_2.gg[2] = (l2_1.x[1] - 5.) * 2.;
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    l5_2.gg[1] = (l2_1.x[0] - 6.) * -2.;
    l5_2.gg[3] = (l2_1.x[1] - 5.) * -2.;
L8:
    return 0;
} /* tp19_ */


/* Subroutine */ int tp20_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 3;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 1.;
    l2_1.x[1] = 1.;
    l11_1.lxl[0] = true;
    l11_1.lxl[1] = false;
    l12_1.lxu[0] = true;
    l12_1.lxu[1] = false;
    l13_1.xl[0] = -.5;
    l14_1.xu[0] = .5;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = .5;
    l20_1.xex[1] = std::sqrt(3.) * .5;
    l20_1.fex = 81.5 - std::sqrt(3.) * 25.;
    l5_3.gg[0] = 1.;
    l5_3.gg[4] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = d__1 * d__1 * 100. + d__3 * d__3;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (l2_1.x[1] - d__1 * d__1) * 200. * (float).01;
    l4_1.gf[0] = (l2_1.x[0] * (l4_1.gf[1] - 1.) + 1.) * -2. * (float).01;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[1];
	l3_3.g[0] = d__1 * d__1 + l2_1.x[0];
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_3.g[1] = d__1 * d__1 + l2_1.x[1];
    }
    if (l9_4.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_3.g[2] = d__1 * d__1 + d__2 * d__2 - 1.;
    }
    return 0;
labelL5:
    if (l10_4.index2[0]) {
	l5_3.gg[3] = l2_1.x[1] * 2.;
    }
    if (l10_4.index2[1]) {
	l5_3.gg[1] = l2_1.x[0] * 2.;
    }
    if (! l10_4.index2[2]) {
	goto labelL9;
    }
    l5_3.gg[2] = l2_1.x[0] * 2.;
    l5_3.gg[5] = l2_1.x[1] * 2.;
labelL9:
    return 0;
} /* tp20_ */


/* Subroutine */ int tp21_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 1;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -1.;
    l2_1.x[1] = -1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = true;
	l12_1.lxu[i__ - 1] = true;
/* labelL6: */
	l14_1.xu[i__ - 1] = 50.;
    }
    l13_1.xl[0] = 2.;
    l13_1.xl[1] = -50.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = 2.;
    l20_1.xex[1] = 0.;
    l20_1.fex = -.99959999999999993;
    l5_1.gg[0] = 10.;
    l5_1.gg[1] = -1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l6_1.fx = (d__1 * d__1 * (float).01 + d__2 * d__2 - 100.) * (float).01;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * (float).02 * (float).01;
    l4_1.gf[1] = l2_1.x[1] * (float)2. * (float).01;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = l2_1.x[0] * 10. - l2_1.x[1] - 10.;
    }
labelL5:
    return 0;
} /* tp21_ */


/* Subroutine */ int tp22_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 1;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 2.;
    l2_1.x[1] = 2.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 1.;
    l20_1.fex = 1.;
    l5_2.gg[0] = -1.;
    l5_2.gg[2] = -1.;
    l5_2.gg[3] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 2.;
/* Computing 2nd power */
    d__2 = l2_1.x[1] - 1.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = (l2_1.x[0] - 2.) * 2.;
    l4_1.gf[1] = (l2_1.x[1] - 1.) * 2.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = 2. - l2_1.x[0] - l2_1.x[1];
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_2.g[1] = l2_1.x[1] - d__1 * d__1;
    }
    return 0;
labelL5:
    if (l10_3.index2[1]) {
	l5_2.gg[1] = l2_1.x[0] * -2.;
    }
    return 0;
} /* tp22_ */


/* Subroutine */ int tp23_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 1;
    l1_1.ninl = 4;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 3.;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = true;
	l12_1.lxu[i__ - 1] = true;
	l13_1.xl[i__ - 1] = -50.;
/* labelL6: */
	l14_1.xu[i__ - 1] = 50.;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 1.;
    l20_1.fex = 2.;
    l5_4.gg[0] = 1.;
    l5_4.gg[5] = 1.;
    l5_4.gg[8] = -1.;
    l5_4.gg[4] = -1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2.;
    l4_1.gf[1] = l2_1.x[1] * 2.;
    return 0;
labelL4:
    if (l9_5.index1[0]) {
	l3_4.g[0] = l2_1.x[0] + l2_1.x[1] - 1.;
    }
    if (l9_5.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_4.g[1] = d__1 * d__1 + d__2 * d__2 - 1.;
    }
    if (l9_5.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_4.g[2] = d__1 * d__1 * 9. + d__2 * d__2 - 9.;
    }
    if (l9_5.index1[3]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_4.g[3] = d__1 * d__1 - l2_1.x[1];
    }
    if (l9_5.index1[4]) {
/* Computing 2nd power */
	d__1 = l2_1.x[1];
	l3_4.g[4] = d__1 * d__1 - l2_1.x[0];
    }
    return 0;
labelL5:
    if (! l10_5.index2[1]) {
	goto L8;
    }
    l5_4.gg[1] = l2_1.x[0] * 2.;
    l5_4.gg[6] = l2_1.x[1] * 2.;
L8:
    if (! l10_5.index2[2]) {
	goto labelL9;
    }
    l5_4.gg[2] = l2_1.x[0] * 18.;
    l5_4.gg[7] = l2_1.x[1] * 2.;
labelL9:
    if (l10_5.index2[3]) {
	l5_4.gg[3] = l2_1.x[0] * 2.;
    }
    if (l10_5.index2[4]) {
	l5_4.gg[9] = l2_1.x[1] * 2.;
    }
    return 0;
} /* tp23_ */


/* Subroutine */ int tp24_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real a;
    static int i__;

    a = std::sqrt(3.);
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 3;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 1.;
    l2_1.x[1] = .5;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = true;
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l13_1.xl[i__ - 1] = 0.;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.xex[0] = 3.;
    l20_1.xex[1] = a;
    l20_1.fex = -1.;
    l5_3.gg[0] = 1. / a;
    l5_3.gg[3] = -1.;
    l5_3.gg[1] = 1.;
    l5_3.gg[4] = a;
    l5_3.gg[2] = -1.;
    l5_3.gg[5] = -a;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 3.;
/* Computing 3rd power */
    d__2 = l2_1.x[1];
    l6_1.fx = (d__1 * d__1 - 9.) * (d__2 * (d__2 * d__2)) / (a * 27.);
    return 0;
labelL3:
/* Computing 3rd power */
    d__1 = l2_1.x[1];
    l4_1.gf[0] = (l2_1.x[0] - 3.) * 2. * (d__1 * (d__1 * d__1)) / (a * 27.);
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 3.;
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l4_1.gf[1] = (d__1 * d__1 - 9.) * (d__2 * d__2) / (a * 9.);
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_1.x[0] / a - l2_1.x[1];
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_1.x[0] + l2_1.x[1] * a;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = 6. - l2_1.x[1] * a - l2_1.x[0];
    }
labelL5:
    return 0;
} /* tp24_ */


/* Subroutine */ int tp25_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real a[99], b[99];
    static int i__, j;
    static Real s, t, u[99], v1, v2, da[297]	/* was [99][3] */, 
	    v11, v22;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_2.x[0] = 100.;
    l2_2.x[1] = 12.5;
    l2_2.x[2] = 3.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = true;
/* labelL6: */
	l12_2.lxu[i__ - 1] = true;
    }
    l13_2.xl[0] = .1;
    l13_2.xl[1] = 1e-5;
    l13_2.xl[2] = 1e-5;
    l14_2.xu[0] = 100.;
    l14_2.xu[1] = 25.6;
    l14_2.xu[2] = 5.;
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.xex[0] = 50.;
    l20_3.xex[1] = 25.;
    l20_3.xex[2] = 3.;
    l20_3.fex = 0.;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 99; ++i__) {
	v1 = .66666666666666663;
	d__1 = std::log((Real) i__ * .01) * -50.;
	u[i__ - 1] = pow_dd(&d__1, &v1) + 25.;
	v1 = u[i__ - 1] - l2_2.x[1];
	if (v1 < 0.) {
	    goto L7;
	}
	v11 = -pow_dd(&v1, &l2_2.x[2]) / l2_2.x[0];
	b[i__ - 1] = std::exp(v11);
/* L30: */
	a[i__ - 1] = b[i__ - 1] - (Real) i__ * .01;
    }
    t = 0.;
    for (i__ = 1; i__ <= 99; ++i__) {
/* L31: */
/* Computing 2nd power */
	d__1 = a[i__ - 1];
	t += d__1 * d__1;
    }
    l6_1.fx = t;
    return 0;
L7:
    s = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__1 = l2_2.x[i__ - 1] - 5.;
	s += d__1 * d__1;
    }
    l6_1.fx = s;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 99; ++i__) {
	v1 = .66666666666666663;
	d__1 = std::log((Real) i__ * .01) * -50.;
	u[i__ - 1] = pow_dd(&d__1, &v1) + 25.;
	v2 = u[i__ - 1] - l2_2.x[1];
	if (v2 <= 0.) {
	    goto labelL9;
	}
	v22 = -pow_dd(&v2, &l2_2.x[2]) / l2_2.x[0];
	if (v22 > -150.) {
	    goto L42;
	}
	b[i__ - 1] = 0.;
	goto L43;
L42:
	b[i__ - 1] = std::exp(v22);
L43:
	a[i__ - 1] = b[i__ - 1] - (Real) i__ * .01;
/* Computing 2nd power */
	d__1 = l2_2.x[0];
	da[i__ - 1] = pow_dd(&v2, &l2_2.x[2]) / (d__1 * d__1) * b[i__ - 1];
	d__1 = l2_2.x[2] - 1.;
	da[i__ + 98] = l2_2.x[2] * pow_dd(&v2, &d__1) / l2_2.x[0] * b[i__ - 1]
		;
/* L36: */
	da[i__ + 197] = -pow_dd(&v2, &l2_2.x[2]) / l2_2.x[0] * std::log(v2) * b[
		i__ - 1];
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	t = 0.;
	for (j = 1; j <= 99; ++j) {
/* L33: */
	    t += a[j - 1] * 2. * da[j + i__ * 99 - 100];
	}
/* L34: */
	l4_2.gf[i__ - 1] = t;
    }
    return 0;
labelL9:
    for (i__ = 1; i__ <= 3; ++i__) {
/* labelL10: */
	l4_2.gf[i__ - 1] = (l2_2.x[i__ - 1] - 5.) * 2.;
    }
labelL4:
    return 0;
} /* tp25_ */


/* Subroutine */ int tp26_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real a;
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    l2_2.x[0] = -2.6;
    l2_2.x[1] = 2.;
    l2_2.x[2] = 2.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l20_4.lex = true;
    l20_4.nex = 2;
    l20_4.xex[0] = 1.;
    l20_4.xex[1] = 1.;
    l20_4.xex[2] = 1.;
    a = std::sqrt(1.287037037037037);
    d__1 = a - 1.1296296296296295;
    d__2 = a + 1.1296296296296295;
    l20_4.xex[3] = pow_dd(&d__1, &c_b74) - pow_dd(&d__2, &c_b74) - 
	    .66666666666666663;
    l20_4.xex[4] = l20_4.xex[3];
    l20_4.xex[5] = l20_4.xex[3];
    l20_4.fex = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0] - l2_2.x[1];
/* Computing 4th power */
    d__2 = l2_2.x[1] - l2_2.x[2], d__2 *= d__2;
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_2.gf[0] = (l2_2.x[0] - l2_2.x[1]) * 2.;
/* Computing 3rd power */
    d__1 = l2_2.x[1] - l2_2.x[2];
    l4_2.gf[2] = d__1 * (d__1 * d__1) * -4.;
    l4_2.gf[1] = -l4_2.gf[0] - l4_2.gf[2];
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[1];
/* Computing 4th power */
	d__2 = l2_2.x[2], d__2 *= d__2;
	l3_1.g[0] = l2_2.x[0] * (d__1 * d__1 + 1.) + d__2 * d__2 - 3.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
/* Computing 2nd power */
    d__1 = l2_2.x[1];
    l5_5.gg[0] = d__1 * d__1 + 1.;
    l5_5.gg[1] = l2_2.x[0] * 2. * l2_2.x[1];
/* Computing 3rd power */
    d__1 = l2_2.x[2];
    l5_5.gg[2] = d__1 * (d__1 * d__1) * 4.;
L7:
    return 0;
} /* tp26_ */


/* Subroutine */ int tp27_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 2.;
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.xex[0] = -1.;
    l20_3.xex[1] = 1.;
    l20_3.xex[2] = 0.;
    l20_3.fex = (float)4.;
    l4_2.gf[2] = 0.;
    l5_5.gg[0] = 1.;
    l5_5.gg[1] = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0] - 1.;
/* Computing 2nd power */
    d__3 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1] - d__3 * d__3;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 * (float)100.;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
    l4_2.gf[1] = (l2_2.x[1] - d__1 * d__1) * 200.;
    l4_2.gf[0] = (l2_2.x[0] * (.01 - l4_2.gf[1]) - .01) * 200.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[2];
	l3_1.g[0] = l2_2.x[0] + d__1 * d__1 + 1.;
    }
    return 0;
labelL5:
    if (l10_2.index2[0]) {
	l5_5.gg[2] = l2_2.x[2] * 2.;
    }
    return 0;
} /* tp27_ */


/* Subroutine */ int tp28_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    l2_2.x[0] = -4.;
    l2_2.x[1] = 1.;
    l2_2.x[2] = 1.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.xex[0] = .5;
    l20_3.xex[1] = -.5;
    l20_3.xex[2] = .5;
    l20_3.fex = 0.;
    l5_5.gg[0] = 1.;
    l5_5.gg[1] = 2.;
    l5_5.gg[2] = 3.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0] + l2_2.x[1];
/* Computing 2nd power */
    d__2 = l2_2.x[1] + l2_2.x[2];
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_2.gf[0] = (l2_2.x[0] + l2_2.x[1]) * 2.;
    l4_2.gf[2] = (l2_2.x[1] + l2_2.x[2]) * 2.;
    l4_2.gf[1] = l4_2.gf[0] + l4_2.gf[2];
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = l2_2.x[0] + l2_2.x[1] * 2. + l2_2.x[2] * 3. - 1.;
    }
labelL5:
    return 0;
} /* tp28_ */


/* Subroutine */ int tp29_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 1.;
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l20_5.lex = true;
    l20_5.nex = 4;
    l20_5.xex[0] = 4.;
    l20_5.xex[1] = std::sqrt(2.) * 2.;
    l20_5.xex[2] = 2.;
    l20_5.xex[3] = l20_5.xex[0];
    l20_5.xex[4] = -l20_5.xex[1];
    l20_5.xex[5] = -l20_5.xex[2];
    l20_5.xex[6] = -l20_5.xex[0];
    l20_5.xex[7] = l20_5.xex[1];
    l20_5.xex[8] = -l20_5.xex[2];
    l20_5.xex[9] = -l20_5.xex[0];
    l20_5.xex[10] = -l20_5.xex[1];
    l20_5.xex[11] = l20_5.xex[2];
    l20_5.fex = std::sqrt(2.) * -16.;
    return 0;
labelL2:
    l6_1.fx = -l2_2.x[0] * l2_2.x[1] * l2_2.x[2];
    return 0;
labelL3:
    l4_2.gf[0] = -l2_2.x[1] * l2_2.x[2];
    l4_2.gf[1] = -l2_2.x[0] * l2_2.x[2];
    l4_2.gf[2] = -l2_2.x[0] * l2_2.x[1];
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[1];
/* Computing 2nd power */
	d__3 = l2_2.x[2];
	l3_1.g[0] = 48. - d__1 * d__1 - d__2 * d__2 * 2. - d__3 * d__3 * 4.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_5.gg[0] = l2_2.x[0] * -2.;
    l5_5.gg[1] = l2_2.x[1] * -4.;
    l5_5.gg[2] = l2_2.x[2] * -8.;
L7:
    return 0;
} /* tp29_ */


/* Subroutine */ int tp30_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 1.;
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = true;
/* labelL6: */
	l14_2.xu[i__ - 1] = 10.;
    }
    l13_2.xl[0] = 1.;
    l13_2.xl[1] = -10.;
    l13_2.xl[2] = -10.;
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.xex[0] = 1.;
    l20_3.xex[1] = 0.;
    l20_3.xex[2] = 0.;
    l20_3.fex = 1.;
    l5_5.gg[2] = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1];
/* Computing 2nd power */
    d__3 = l2_2.x[2];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return 0;
labelL3:
    l4_2.gf[0] = l2_2.x[0] * 2.;
    l4_2.gf[1] = l2_2.x[1] * 2.;
    l4_2.gf[2] = l2_2.x[2] * 2.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[1];
	l3_1.g[0] = d__1 * d__1 + d__2 * d__2 - 1.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_5.gg[0] = l2_2.x[0] * 2.;
    l5_5.gg[1] = l2_2.x[1] * 2.;
L7:
    return 0;
} /* tp30_ */


/* Subroutine */ int tp31_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 1.;
	l11_2.lxl[i__ - 1] = true;
/* labelL6: */
	l12_2.lxu[i__ - 1] = true;
    }
    l13_2.xl[0] = -10.;
    l13_2.xl[1] = 1.;
    l13_2.xl[2] = -10.;
    l14_2.xu[0] = 10.;
    l14_2.xu[1] = 10.;
    l14_2.xu[2] = 1.;
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.xex[0] = 1. / std::sqrt(3.);
    l20_3.xex[1] = std::sqrt(3.);
    l20_3.xex[2] = 0.;
    l20_3.fex = 6.;
    l5_5.gg[2] = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1];
/* Computing 2nd power */
    d__3 = l2_2.x[2];
    l6_1.fx = d__1 * d__1 * 9. + d__2 * d__2 + d__3 * d__3 * 9.;
    return 0;
labelL3:
    l4_2.gf[0] = l2_2.x[0] * 18.;
    l4_2.gf[1] = l2_2.x[1] * 2.;
    l4_2.gf[2] = l2_2.x[2] * 18.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = l2_2.x[0] * l2_2.x[1] - 1.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_5.gg[0] = l2_2.x[1];
    l5_5.gg[1] = l2_2.x[0];
L7:
    return 0;
} /* tp31_ */


/* Subroutine */ int tp32_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    l2_2.x[0] = .1;
    l2_2.x[1] = .7;
    l2_2.x[2] = .2;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = false;
/* labelL6: */
	l13_2.xl[i__ - 1] = 0.;
    }
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.xex[0] = 0.;
    l20_3.xex[1] = 0.;
    l20_3.xex[2] = 1.;
    l20_3.fex = 1.;
    l5_3.gg[2] = 6.;
    l5_3.gg[4] = 4.;
    l5_3.gg[1] = -1.;
    l5_3.gg[3] = -1.;
    l5_3.gg[5] = -1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0] + l2_2.x[1] * 3. + l2_2.x[2];
/* Computing 2nd power */
    d__2 = l2_2.x[0] - l2_2.x[1];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 * 4.;
    return 0;
labelL3:
    l4_2.gf[0] = l2_2.x[0] * 10. - l2_2.x[1] * 2. + l2_2.x[2] * 2.;
    l4_2.gf[1] = l2_2.x[0] * -2. + l2_2.x[1] * 26. + l2_2.x[2] * 6.;
    l4_2.gf[2] = (l2_2.x[0] + l2_2.x[1] * 3. + l2_2.x[2]) * 2.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 3rd power */
	d__1 = l2_2.x[0];
	l3_2.g[0] = -(d__1 * (d__1 * d__1)) + l2_2.x[1] * 6. + l2_2.x[2] * 4. 
		- 3.;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = 1. - l2_2.x[0] - l2_2.x[1] - l2_2.x[2];
    }
    return 0;
labelL5:
    if (l10_3.index2[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
	l5_3.gg[0] = d__1 * d__1 * -3.;
    }
    return 0;
} /* tp32_ */


/* Subroutine */ int tp33_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;


    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_2.x[0] = 0.;
    l2_2.x[1] = 0.;
    l2_2.x[2] = 3.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = true;
/* labelL6: */
	l13_2.xl[i__ - 1] = 0.;
    }
    l12_2.lxu[0] = false;
    l12_2.lxu[1] = false;
    l12_2.lxu[2] = true;
    l14_2.xu[2] = 5.;
    l20_3.lex = true;
    l20_3.nex = 2;
    l20_3.xex[0] = 0.;
    l20_3.xex[1] = std::sqrt(2.);
    l20_3.xex[2] = std::sqrt(2.);
    l20_3.fex = std::sqrt(2.) - (float)6.;
    l4_2.gf[1] = 0.;
    l4_2.gf[2] = 1.;
    return 0;
labelL2:
    l6_1.fx = (l2_2.x[0] - 1.) * (l2_2.x[0] - 2.) * (l2_2.x[0] - 3.) + l2_2.x[
	    2];
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
    l4_2.gf[0] = d__1 * d__1 * 3. - l2_2.x[0] * 12. + 11.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[2];
/* Computing 2nd power */
	d__2 = l2_2.x[0];
/* Computing 2nd power */
	d__3 = l2_2.x[1];
	l3_2.g[0] = d__1 * d__1 - d__2 * d__2 - d__3 * d__3;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[1];
/* Computing 2nd power */
	d__3 = l2_2.x[2];
	l3_2.g[1] = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 - 4.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_3.gg[0] = l2_2.x[0] * -2.;
    l5_3.gg[2] = l2_2.x[1] * -2.;
    l5_3.gg[4] = l2_2.x[2] * 2.;
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    l5_3.gg[1] = l2_2.x[0] * 2.;
    l5_3.gg[3] = l2_2.x[1] * 2.;
    l5_3.gg[5] = l2_2.x[2] * 2.;
L8:
    return 0;
} /* tp33_ */


/* Subroutine */ int tp34_(int *mode)
{
    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_2.x[0] = 0.;
    l2_2.x[1] = 1.05;
    l2_2.x[2] = 2.9;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = true;
/* labelL6: */
	l13_2.xl[i__ - 1] = 0.;
    }
    l14_2.xu[0] = 100.;
    l14_2.xu[1] = 100.;
    l14_2.xu[2] = 10.;
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.xex[0] = std::log(std::log(10.));
    l20_3.xex[1] = std::log(10.);
    l20_3.xex[2] = 10.;
    l20_3.fex = -l20_3.xex[0];
    l4_2.gf[0] = -1.;
    l4_2.gf[1] = 0.;
    l4_2.gf[2] = 0.;
    l5_3.gg[2] = 1.;
    l5_3.gg[4] = 0.;
    l5_3.gg[1] = 0.;
    l5_3.gg[5] = 1.;
    return 0;
labelL2:
    l6_1.fx = -l2_2.x[0];
labelL3:
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_2.x[1] - std::exp(l2_2.x[0]);
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = l2_2.x[2] - std::exp(l2_2.x[1]);
    }
    return 0;
labelL5:
    if (l10_3.index2[0]) {
	l5_3.gg[0] = -std::exp(l2_2.x[0]);
    }
    if (l10_3.index2[1]) {
	l5_3.gg[3] = -std::exp(l2_2.x[1]);
    }
    return 0;
} /* tp34_ */


/* Subroutine */ int tp35_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 1;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = .5;
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = false;
/* labelL6: */
	l13_2.xl[i__ - 1] = 0.;
    }
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.xex[0] = 1.3333333333333333;
    l20_3.xex[1] = .77777777777777779;
    l20_3.xex[2] = .44444444444444442;
    l20_3.fex = .1111111111111111;
    l5_5.gg[0] = -1.;
    l5_5.gg[1] = -1.;
    l5_5.gg[2] = -2.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1];
/* Computing 2nd power */
    d__3 = l2_2.x[2];
    l6_1.fx = 9. - l2_2.x[0] * 8. - l2_2.x[1] * 6. - l2_2.x[2] * 4. + d__1 * 
	    d__1 * 2. + d__2 * d__2 * 2. + d__3 * d__3 + l2_2.x[0] * 2. * 
	    l2_2.x[1] + l2_2.x[0] * 2. * l2_2.x[2];
    return 0;
labelL3:
    l4_2.gf[0] = l2_2.x[0] * 4. - 8. + l2_2.x[1] * 2. + l2_2.x[2] * 2.;
    l4_2.gf[1] = l2_2.x[1] * 4. - 6. + l2_2.x[0] * 2.;
    l4_2.gf[2] = l2_2.x[2] * 2. - 4. + l2_2.x[0] * 2.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = -l2_2.x[0] - l2_2.x[1] - l2_2.x[2] * 2. + 3.;
    }
labelL5:
    return 0;
} /* tp35_ */


/* Subroutine */ int tp36_(int *mode)
{
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 1;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 10.;
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = true;
/* labelL6: */
	l13_2.xl[i__ - 1] = 0.;
    }
    l14_2.xu[0] = 20.;
    l14_2.xu[1] = 11.;
    l14_2.xu[2] = 42.;
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.xex[0] = 20.;
    l20_3.xex[1] = 11.;
    l20_3.xex[2] = 15.;
    l20_3.fex = -3300.;
    l5_5.gg[0] = -1.;
    l5_5.gg[1] = -2.;
    l5_5.gg[2] = -2.;
    return 0;
labelL2:
    l6_1.fx = -l2_2.x[0] * l2_2.x[1] * l2_2.x[2];
    return 0;
labelL3:
    l4_2.gf[0] = -l2_2.x[1] * l2_2.x[2];
    l4_2.gf[1] = -l2_2.x[0] * l2_2.x[2];
    l4_2.gf[2] = -l2_2.x[0] * l2_2.x[1];
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = 72. - l2_2.x[0] - l2_2.x[1] * 2. - l2_2.x[2] * 2.;
    }
    return 0;
labelL5:
    return 0;
} /* tp36_ */


/* Subroutine */ int tp37_(int *mode)
{
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 2;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 10.;
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = true;
	l14_2.xu[i__ - 1] = 42.;
/* labelL6: */
	l13_2.xl[i__ - 1] = 0.;
    }
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.xex[0] = 24.;
    l20_3.xex[1] = 12.;
    l20_3.xex[2] = 12.;
    l20_3.fex = -3456.;
    l5_3.gg[0] = -1.;
    l5_3.gg[2] = -2.;
    l5_3.gg[4] = -2.;
    l5_3.gg[1] = 1.;
    l5_3.gg[3] = 2.;
    l5_3.gg[5] = 2.;
    return 0;
labelL2:
    l6_1.fx = -l2_2.x[0] * l2_2.x[1] * l2_2.x[2];
    return 0;
labelL3:
    l4_2.gf[0] = -l2_2.x[1] * l2_2.x[2];
    l4_2.gf[1] = -l2_2.x[0] * l2_2.x[2];
    l4_2.gf[2] = -l2_2.x[0] * l2_2.x[1];
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = 72. - l2_2.x[0] - l2_2.x[1] * 2. - l2_2.x[2] * 2.;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = l2_2.x[0] + l2_2.x[1] * 2. + l2_2.x[2] * 2.;
    }
labelL5:
    return 0;
} /* tp37_ */


/* Subroutine */ int tp38_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_3.x[0] = -3.;
    l2_3.x[1] = -1.;
    l2_3.x[2] = -3.;
    l2_3.x[3] = -1.;
    for (i__ = 1; i__ <= 4; ++i__) {
	l11_3.lxl[i__ - 1] = true;
	l12_3.lxu[i__ - 1] = true;
	l13_3.xl[i__ - 1] = -10.;
/* labelL6: */
	l14_3.xu[i__ - 1] = 10.;
    }
    l20_6.lex = true;
    l20_6.nex = 1;
    for (i__ = 1; i__ <= 4; ++i__) {
/* L30: */
	l20_6.xex[i__ - 1] = 1.;
    }
    l20_6.fex = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_3.x[0];
/* Computing 2nd power */
    d__1 = l2_3.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_3.x[0];
/* Computing 2nd power */
    d__5 = l2_3.x[2];
/* Computing 2nd power */
    d__4 = l2_3.x[3] - d__5 * d__5;
/* Computing 2nd power */
    d__6 = 1. - l2_3.x[2];
/* Computing 2nd power */
    d__7 = l2_3.x[1] - 1.;
/* Computing 2nd power */
    d__8 = l2_3.x[3] - 1.;
    l6_1.fx = (d__1 * d__1 * 100. + d__3 * d__3 + d__4 * d__4 * 90. + d__6 * 
	    d__6 + (d__7 * d__7 + d__8 * d__8) * 10.1 + (l2_3.x[1] - 1.) * 
	    19.8 * (l2_3.x[3] - 1.)) * 1e-5;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l4_3.gf[0] = (l2_3.x[0] * -400. * (l2_3.x[1] - d__1 * d__1) - (1. - 
	    l2_3.x[0]) * 2.) * 1e-5;
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l4_3.gf[1] = ((l2_3.x[1] - d__1 * d__1) * 200. + (l2_3.x[1] - 1.) * 20.2 
	    + (l2_3.x[3] - 1.) * 19.8) * 1e-5;
/* Computing 2nd power */
    d__1 = l2_3.x[2];
    l4_3.gf[2] = (l2_3.x[2] * -360. * (l2_3.x[3] - d__1 * d__1) - (1. - 
	    l2_3.x[2]) * 2.) * 1e-5;
/* Computing 2nd power */
    d__1 = l2_3.x[2];
    l4_3.gf[3] = ((l2_3.x[3] - d__1 * d__1) * 180. + (l2_3.x[3] - 1.) * 20.2 
	    + (l2_3.x[1] - 1.) * 19.8) * 1e-5;
labelL4:
    return 0;
} /* tp38_ */


/* Subroutine */ int tp39_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 2;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 2.;
	l11_3.lxl[i__ - 1] = false;
/* labelL6: */
	l12_3.lxu[i__ - 1] = false;
    }
    l20_6.nex = 1;
    l20_6.lex = true;
    l20_6.xex[0] = 1.;
    l20_6.xex[1] = 1.;
    l20_6.xex[2] = 0.;
    l20_6.xex[3] = 0.;
    l20_6.fex = -1.;
    l4_3.gf[0] = -1.;
    l4_3.gf[1] = 0.;
    l4_3.gf[2] = 0.;
    l4_3.gf[3] = 0.;
    l5_6.gg[2] = 1.;
    l5_6.gg[6] = 0.;
    l5_6.gg[3] = -1.;
    l5_6.gg[5] = 0.;
    return 0;
labelL2:
    l6_1.fx = -l2_3.x[0];
labelL3:
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 3rd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[2];
	l3_2.g[0] = l2_3.x[1] - d__1 * (d__1 * d__1) - d__2 * d__2;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[3];
	l3_2.g[1] = d__1 * d__1 - l2_3.x[1] - d__2 * d__2;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l5_6.gg[0] = d__1 * d__1 * -3.;
    l5_6.gg[4] = l2_3.x[2] * -2.;
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    l5_6.gg[1] = l2_3.x[0] * 2.;
    l5_6.gg[7] = l2_3.x[3] * -2.;
L8:
    return 0;
} /* tp39_ */


/* Subroutine */ int tp40_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 3;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = .8;
	l11_3.lxl[i__ - 1] = false;
/* labelL6: */
	l12_3.lxu[i__ - 1] = false;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 4; ++j) {
/* L15: */
	    l5_7.gg[i__ + j * 3 - 4] = 0.;
	}
    }
    l5_7.gg[7] = -1.;
    l5_7.gg[5] = -1.;
    l20_2.lex = true;
    l20_2.xex[0] = pow_dd(&c_b305, &c_b306);
    l20_2.xex[1] = pow_dd(&c_b305, &c_b308);
    l20_2.xex[2] = pow_dd(&c_b305, &c_b310);
    l20_2.xex[3] = pow_dd(&c_b305, &c_b312);
    l20_2.xex[4] = l20_2.xex[0];
    l20_2.xex[5] = l20_2.xex[1];
    l20_2.xex[6] = -l20_2.xex[2];
    l20_2.xex[7] = -l20_2.xex[3];
    l20_2.fex = -.25;
    l20_2.nex = 2;
    return 0;
labelL2:
    l6_1.fx = -l2_3.x[0] * l2_3.x[1] * l2_3.x[2] * l2_3.x[3];
    return 0;
labelL3:
    l4_3.gf[0] = -l2_3.x[1] * l2_3.x[2] * l2_3.x[3];
    l4_3.gf[1] = -l2_3.x[0] * l2_3.x[2] * l2_3.x[3];
    l4_3.gf[2] = -l2_3.x[0] * l2_3.x[1] * l2_3.x[3];
    l4_3.gf[3] = -l2_3.x[0] * l2_3.x[1] * l2_3.x[2];
    return 0;
labelL4:
    if (l9_4.index1[0]) {
/* Computing 3rd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[1];
	l3_3.g[0] = d__1 * (d__1 * d__1) + d__2 * d__2 - 1.;
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_3.x[0];
	l3_3.g[1] = d__1 * d__1 * l2_3.x[3] - l2_3.x[2];
    }
    if (l9_4.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_3.x[3];
	l3_3.g[2] = d__1 * d__1 - l2_3.x[1];
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto L7;
    }
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l5_7.gg[0] = d__1 * d__1 * 3.;
    l5_7.gg[3] = l2_3.x[1] * 2.;
L7:
    if (! l10_4.index2[1]) {
	goto L8;
    }
    l5_7.gg[1] = l2_3.x[0] * 2. * l2_3.x[3];
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l5_7.gg[10] = d__1 * d__1;
L8:
    if (! l10_4.index2[2]) {
	goto labelL9;
    }
    l5_7.gg[11] = l2_3.x[3] * 2.;
labelL9:
    return 0;
} /* tp40_ */


/* Subroutine */ int tp41_(int *mode)
{
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 2.;
	l11_3.lxl[i__ - 1] = true;
	l12_3.lxu[i__ - 1] = true;
/* labelL6: */
	l13_3.xl[i__ - 1] = 0.;
    }
    l14_3.xu[0] = 1.;
    l14_3.xu[1] = 1.;
    l14_3.xu[2] = 1.;
    l14_3.xu[3] = 2.;
    l4_3.gf[3] = 0.;
    l5_2.gg[0] = 1.;
    l5_2.gg[1] = 2.;
    l5_2.gg[2] = 2.;
    l5_2.gg[3] = -1.;
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.xex[0] = .66666666666666663;
    l20_6.xex[1] = .33333333333333331;
    l20_6.xex[2] = l20_6.xex[1];
    l20_6.xex[3] = 2.;
    l20_6.fex = 1.9259259259259258;
    return 0;
labelL2:
    l6_1.fx = 2. - l2_3.x[0] * l2_3.x[1] * l2_3.x[2];
    return 0;
labelL3:
    l4_3.gf[0] = -l2_3.x[1] * l2_3.x[2];
    l4_3.gf[1] = -l2_3.x[0] * l2_3.x[2];
    l4_3.gf[2] = -l2_3.x[0] * l2_3.x[1];
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = l2_3.x[0] + l2_3.x[1] * 2. + l2_3.x[2] * 2. - l2_3.x[3];
    }
labelL5:
    return 0;
} /* tp41_ */


/* Subroutine */ int tp42_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 1;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 1.;
	l12_3.lxu[i__ - 1] = false;
/* labelL6: */
	l11_3.lxl[i__ - 1] = false;
    }
    l5_6.gg[0] = 1.;
    l5_6.gg[2] = 0.;
    l5_6.gg[4] = 0.;
    l5_6.gg[6] = 0.;
    l5_6.gg[1] = 0.;
    l5_6.gg[3] = 0.;
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.xex[0] = 2.;
    l20_6.xex[1] = 2.;
    l20_6.xex[2] = std::sqrt(.72);
    l20_6.xex[3] = std::sqrt(1.28);
    l20_6.fex = 28. - std::sqrt(2.) * 10.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_3.x[0] - 1.;
/* Computing 2nd power */
    d__2 = l2_3.x[1] - 2.;
/* Computing 2nd power */
    d__3 = l2_3.x[2] - 3.;
/* Computing 2nd power */
    d__4 = l2_3.x[3] - 4.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 4; ++i__) {
/* L21: */
	l4_3.gf[i__ - 1] = (l2_3.x[i__ - 1] - (Real) i__) * 2.;
    }
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_3.x[0] - 2.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_3.x[2];
/* Computing 2nd power */
	d__2 = l2_3.x[3];
	l3_2.g[1] = d__1 * d__1 + d__2 * d__2 - 2.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    l5_6.gg[5] = l2_3.x[2] * 2.;
    l5_6.gg[7] = l2_3.x[3] * 2.;
L8:
    return 0;
} /* tp42_ */


/* Subroutine */ int tp43_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 3;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 0.;
	l11_3.lxl[i__ - 1] = false;
/* labelL6: */
	l12_3.lxu[i__ - 1] = false;
    }
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.xex[0] = 0.;
    l20_6.xex[1] = 1.;
    l20_6.xex[2] = 2.;
    l20_6.xex[3] = -1.;
    l20_6.fex = -44.;
    l5_7.gg[11] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_3.x[0];
/* Computing 2nd power */
    d__2 = l2_3.x[1];
/* Computing 2nd power */
    d__3 = l2_3.x[2];
/* Computing 2nd power */
    d__4 = l2_3.x[3];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 * 2. + d__4 * d__4 - 
	    l2_3.x[0] * 5. - l2_3.x[1] * 5. - l2_3.x[2] * 21. + l2_3.x[3] * 
	    7.;
    return 0;
labelL3:
    l4_3.gf[0] = l2_3.x[0] * 2. - 5.;
    l4_3.gf[1] = l2_3.x[1] * 2. - 5.;
    l4_3.gf[2] = l2_3.x[2] * 4. - 21.;
    l4_3.gf[3] = l2_3.x[3] * 2. + 7.;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[1];
/* Computing 2nd power */
	d__3 = l2_3.x[2];
/* Computing 2nd power */
	d__4 = l2_3.x[3];
	l3_3.g[0] = -(d__1 * d__1) - d__2 * d__2 - d__3 * d__3 - d__4 * d__4 
		- l2_3.x[0] + l2_3.x[1] - l2_3.x[2] + l2_3.x[3] + 8.;
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[1];
/* Computing 2nd power */
	d__3 = l2_3.x[2];
/* Computing 2nd power */
	d__4 = l2_3.x[3];
	l3_3.g[1] = -(d__1 * d__1) - d__2 * d__2 * 2. - d__3 * d__3 - d__4 * 
		d__4 * 2. + l2_3.x[0] + l2_3.x[3] + 10.;
    }
    if (l9_4.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[1];
/* Computing 2nd power */
	d__3 = l2_3.x[2];
	l3_3.g[2] = d__1 * d__1 * -2. - d__2 * d__2 - d__3 * d__3 - l2_3.x[0] 
		* 2. + l2_3.x[1] + l2_3.x[3] + 5.;
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto L7;
    }
    l5_7.gg[0] = l2_3.x[0] * -2. - 1.;
    l5_7.gg[3] = l2_3.x[1] * -2. + 1.;
    l5_7.gg[6] = l2_3.x[2] * -2. - 1.;
    l5_7.gg[9] = l2_3.x[3] * -2. + 1.;
L7:
    if (! l10_4.index2[1]) {
	goto L8;
    }
    l5_7.gg[1] = l2_3.x[0] * -2. + 1.;
    l5_7.gg[4] = l2_3.x[1] * -4.;
    l5_7.gg[7] = l2_3.x[2] * -2.;
    l5_7.gg[10] = l2_3.x[3] * -4. + 1.;
L8:
    if (! l10_4.index2[2]) {
	goto labelL9;
    }
    l5_7.gg[2] = l2_3.x[0] * -4. - 2.;
    l5_7.gg[5] = l2_3.x[1] * -2. + 1.;
    l5_7.gg[8] = l2_3.x[2] * -2.;
labelL9:
    return 0;
} /* tp43_ */


/* Subroutine */ int tp44_(int *mode)
{
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 6;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 0.;
	l13_3.xl[i__ - 1] = 0.;
	l11_3.lxl[i__ - 1] = true;
/* labelL6: */
	l12_3.lxu[i__ - 1] = false;
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	for (j = 1; j <= 4; ++j) {
/* L15: */
	    l5_8.gg[i__ + j * 6 - 7] = 0.;
	}
    }
    l5_8.gg[0] = -1.;
    l5_8.gg[6] = -2.;
    l5_8.gg[1] = -4.;
    l5_8.gg[7] = -1.;
    l5_8.gg[2] = -3.;
    l5_8.gg[8] = -4.;
    l5_8.gg[15] = -2.;
    l5_8.gg[21] = -1.;
    l5_8.gg[16] = -1.;
    l5_8.gg[22] = -2.;
    l5_8.gg[17] = -1.;
    l5_8.gg[23] = -1.;
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.xex[0] = 0.;
    l20_6.xex[1] = 3.;
    l20_6.xex[2] = 0.;
    l20_6.xex[3] = 4.;
    l20_6.fex = -15.;
    return 0;
labelL2:
    l6_1.fx = l2_3.x[0] - l2_3.x[1] - l2_3.x[2] - l2_3.x[0] * l2_3.x[2] + 
	    l2_3.x[0] * l2_3.x[3] + l2_3.x[1] * l2_3.x[2] - l2_3.x[1] * 
	    l2_3.x[3];
    return 0;
labelL3:
    l4_3.gf[0] = 1. - l2_3.x[2] + l2_3.x[3];
    l4_3.gf[1] = l2_3.x[2] - 1. - l2_3.x[3];
    l4_3.gf[2] = -1. - l2_3.x[0] + l2_3.x[1];
    l4_3.gf[3] = l2_3.x[0] - l2_3.x[1];
    return 0;
labelL4:
    if (l9_6.index1[0]) {
	l3_5.g[0] = 8. - l2_3.x[0] - l2_3.x[1] * 2.;
    }
    if (l9_6.index1[1]) {
	l3_5.g[1] = 12. - l2_3.x[0] * 4. - l2_3.x[1];
    }
    if (l9_6.index1[2]) {
	l3_5.g[2] = 12. - l2_3.x[0] * 3. - l2_3.x[1] * 4.;
    }
    if (l9_6.index1[3]) {
	l3_5.g[3] = 8. - l2_3.x[2] * 2. - l2_3.x[3];
    }
    if (l9_6.index1[4]) {
	l3_5.g[4] = 8. - l2_3.x[2] - l2_3.x[3] * 2.;
    }
    if (l9_6.index1[5]) {
	l3_5.g[5] = 5. - l2_3.x[2] - l2_3.x[3];
    }
labelL5:
    return 0;
} /* tp44_ */


/* Subroutine */ int tp45_(int *mode)
{
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 5; ++i__) {
	l2_4.x[i__ - 1] = 2.;
	l11_4.lxl[i__ - 1] = true;
	l12_4.lxu[i__ - 1] = true;
	l13_4.xl[i__ - 1] = 0.;
/* labelL6: */
	l14_4.xu[i__ - 1] = (Real) i__;
    }
    l20_7.lex = true;
    l20_7.nex = 1;
    for (i__ = 1; i__ <= 5; ++i__) {
/* L30: */
	l20_7.xex[i__ - 1] = (Real) i__;
    }
    l20_7.fex = 1.;
    return 0;
labelL2:
    l6_1.fx = 2. - l2_4.x[0] * l2_4.x[1] * l2_4.x[2] * l2_4.x[3] * l2_4.x[4] /
	     120.;
    return 0;
labelL3:
    l4_4.gf[0] = -l2_4.x[1] * l2_4.x[2] * l2_4.x[3] * l2_4.x[4] / 120.;
    l4_4.gf[1] = -l2_4.x[0] * l2_4.x[2] * l2_4.x[3] * l2_4.x[4] / 120.;
    l4_4.gf[2] = -l2_4.x[0] * l2_4.x[1] * l2_4.x[3] * l2_4.x[4] / 120.;
    l4_4.gf[3] = -l2_4.x[0] * l2_4.x[1] * l2_4.x[2] * l2_4.x[4] / 120.;
    l4_4.gf[4] = -l2_4.x[0] * l2_4.x[1] * l2_4.x[2] * l2_4.x[3] / 120.;
labelL4:
    return 0;
} /* tp45_ */


/* Subroutine */ int tp46_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 2;
    l2_4.x[0] = std::sqrt(2.) * .5;
    l2_4.x[1] = 1.75;
    l2_4.x[2] = .5;
    l2_4.x[3] = 2.;
    l2_4.x[4] = 2.;
    for (i__ = 1; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = false;
/* labelL6: */
	l12_4.lxu[i__ - 1] = false;
    }
    l5_4.gg[2] = 0.;
    l5_4.gg[4] = 0.;
    l5_4.gg[1] = 0.;
    l5_4.gg[3] = 1.;
    l5_4.gg[9] = 0.;
    l20_7.lex = true;
    l20_7.nex = 1;
    for (i__ = 1; i__ <= 5; ++i__) {
/* L30: */
	l20_7.xex[i__ - 1] = 1.;
    }
    l20_7.fex = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_4.x[0] - l2_4.x[1];
/* Computing 2nd power */
    d__2 = l2_4.x[2] - 1.;
/* Computing 4th power */
    d__3 = l2_4.x[3] - 1., d__3 *= d__3;
/* Computing 6th power */
    d__4 = l2_4.x[4] - 1., d__4 *= d__4;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * (d__4 * d__4);
    return 0;
labelL3:
    l4_4.gf[0] = (l2_4.x[0] - l2_4.x[1]) * 2.;
    l4_4.gf[1] = -l4_4.gf[0];
    l4_4.gf[2] = (l2_4.x[2] - 1.) * 2.;
/* Computing 3rd power */
    d__1 = l2_4.x[3] - 1.;
    l4_4.gf[3] = d__1 * (d__1 * d__1) * 4.;
/* Computing 5th power */
    d__1 = l2_4.x[4] - 1., d__2 = d__1, d__1 *= d__1;
    l4_4.gf[4] = d__2 * (d__1 * d__1) * 6.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_4.x[0];
	l3_2.g[0] = d__1 * d__1 * l2_4.x[3] + std::sin(l2_4.x[3] - l2_4.x[4]) - 1.;
    }
    if (l9_3.index1[1]) {
/* Computing 4th power */
	d__1 = l2_4.x[2], d__1 *= d__1;
/* Computing 2nd power */
	d__2 = l2_4.x[3];
	l3_2.g[1] = l2_4.x[1] + d__1 * d__1 * (d__2 * d__2) - 2.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_4.gg[0] = l2_4.x[0] * 2. * l2_4.x[3];
    l5_4.gg[8] = -std::cos(l2_4.x[3] - l2_4.x[4]);
/* Computing 2nd power */
    d__1 = l2_4.x[0];
    l5_4.gg[6] = d__1 * d__1 - l5_4.gg[8];
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
/* Computing 3rd power */
    d__1 = l2_4.x[2];
/* Computing 2nd power */
    d__2 = l2_4.x[3];
    l5_4.gg[5] = d__1 * (d__1 * d__1) * 4. * (d__2 * d__2);
/* Computing 4th power */
    d__1 = l2_4.x[2], d__1 *= d__1;
    l5_4.gg[7] = d__1 * d__1 * 2. * l2_4.x[3];
L8:
    return 0;
} /* tp46_ */


/* Subroutine */ int tp47_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__, j;
    static Real v1, v2, v3, v4;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 3;
    l2_4.x[0] = 2.;
    l2_4.x[1] = std::sqrt(2.);
    l2_4.x[2] = -1.;
    l2_4.x[3] = 2. - l2_4.x[1];
    l2_4.x[4] = .5;
    for (i__ = 1; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = false;
/* labelL6: */
	l12_4.lxu[i__ - 1] = false;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 5; ++j) {
/* L15: */
	    l5_9.gg[i__ + j * 3 - 4] = 0.;
	}
    }
    l5_9.gg[0] = 1.;
    l5_9.gg[4] = 1.;
    l5_9.gg[10] = 1.;
    l20_7.lex = true;
    l20_7.nex = 1;
    for (i__ = 1; i__ <= 5; ++i__) {
/* L30: */
	l20_7.xex[i__ - 1] = 1.;
    }
    l20_7.fex = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_4.x[0] - l2_4.x[1];
/* Computing 2nd power */
    d__2 = l2_4.x[1] - l2_4.x[2];
/* Computing 4th power */
    d__3 = l2_4.x[2] - l2_4.x[3], d__3 *= d__3;
/* Computing 4th power */
    d__4 = l2_4.x[3] - l2_4.x[4], d__4 *= d__4;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4;
    return 0;
labelL3:
    v1 = (l2_4.x[0] - l2_4.x[1]) * 2.;
    v2 = (l2_4.x[1] - l2_4.x[2]) * 2.;
/* Computing 3rd power */
    d__1 = l2_4.x[2] - l2_4.x[3];
    v3 = d__1 * (d__1 * d__1) * 4.;
/* Computing 3rd power */
    d__1 = l2_4.x[3] - l2_4.x[4];
    v4 = d__1 * (d__1 * d__1) * 4.;
    l4_4.gf[0] = v1;
    l4_4.gf[1] = -v1 + v2;
    l4_4.gf[2] = -v2 + v3;
    l4_4.gf[3] = -v3 + v4;
    l4_4.gf[4] = -v4;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_4.x[1];
/* Computing 3rd power */
	d__2 = l2_4.x[2];
	l3_3.g[0] = l2_4.x[0] + d__1 * d__1 + d__2 * (d__2 * d__2) - 3.;
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_4.x[2];
	l3_3.g[1] = l2_4.x[1] - d__1 * d__1 + l2_4.x[3] - 1.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = l2_4.x[0] * l2_4.x[4] - 1.;
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto L7;
    }
    l5_9.gg[3] = l2_4.x[1] * 2.;
/* Computing 2nd power */
    d__1 = l2_4.x[2];
    l5_9.gg[6] = d__1 * d__1 * 3.;
L7:
    if (l10_4.index2[1]) {
	l5_9.gg[7] = l2_4.x[2] * -2.;
    }
/* L8: */
    if (! l10_4.index2[2]) {
	goto labelL9;
    }
    l5_9.gg[2] = l2_4.x[4];
    l5_9.gg[14] = l2_4.x[0];
labelL9:
    return 0;
} /* tp47_ */


/* Subroutine */ int tp48_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 2;
    l1_1.nenl = 0;
    l2_4.x[0] = 3.;
    l2_4.x[1] = 5.;
    l2_4.x[2] = -3.;
    l2_4.x[3] = 2.;
    l2_4.x[4] = -2.;
    for (i__ = 1; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = false;
/* labelL6: */
	l12_4.lxu[i__ - 1] = false;
    }
    l5_4.gg[1] = 0.;
    l5_4.gg[3] = 0.;
    for (i__ = 1; i__ <= 5; ++i__) {
/* labelL20: */
	l5_4.gg[(i__ << 1) - 2] = 1.;
    }
    l5_4.gg[5] = 1.;
    l5_4.gg[7] = -2.;
    l5_4.gg[9] = -2.;
    l20_7.lex = true;
    l20_7.nex = 1;
    for (i__ = 1; i__ <= 5; ++i__) {
/* L30: */
	l20_7.xex[i__ - 1] = 1.;
    }
    l20_7.fex = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_4.x[0] - 1.;
/* Computing 2nd power */
    d__2 = l2_4.x[1] - l2_4.x[2];
/* Computing 2nd power */
    d__3 = l2_4.x[3] - l2_4.x[4];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return 0;
labelL3:
    l4_4.gf[0] = (l2_4.x[0] - 1.) * 2.;
    l4_4.gf[1] = (l2_4.x[1] - l2_4.x[2]) * 2.;
    l4_4.gf[2] = -l4_4.gf[1];
    l4_4.gf[3] = (l2_4.x[3] - l2_4.x[4]) * 2.;
    l4_4.gf[4] = -l4_4.gf[3];
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_4.x[0] + l2_4.x[1] + l2_4.x[2] + l2_4.x[3] + l2_4.x[4] 
		- 5.;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = l2_4.x[2] - (l2_4.x[3] + l2_4.x[4]) * 2. + 3.;
    }
labelL5:
    return 0;
} /* tp48_ */


/* Subroutine */ int tp49_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 2;
    l1_1.nenl = 0;
    l2_4.x[0] = 10.;
    l2_4.x[1] = 7.;
    l2_4.x[2] = 2.;
    l2_4.x[3] = -3.;
    l2_4.x[4] = .8;
    for (i__ = 1; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = false;
/* labelL6: */
	l12_4.lxu[i__ - 1] = false;
    }
    l5_4.gg[8] = 0.;
    l5_4.gg[1] = 0.;
    l5_4.gg[3] = 0.;
    l5_4.gg[7] = 0.;
    l5_4.gg[0] = 1.;
    l5_4.gg[2] = 1.;
    l5_4.gg[4] = 1.;
    l5_4.gg[6] = 4.;
    l5_4.gg[5] = 1.;
    l5_4.gg[9] = 5.;
    l20_7.lex = true;
    l20_7.nex = 1;
    for (i__ = 1; i__ <= 5; ++i__) {
/* L30: */
	l20_7.xex[i__ - 1] = 1.;
    }
    l20_7.fex = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_4.x[0] - l2_4.x[1];
/* Computing 2nd power */
    d__2 = l2_4.x[2] - 1.;
/* Computing 4th power */
    d__3 = l2_4.x[3] - 1., d__3 *= d__3;
/* Computing 6th power */
    d__4 = l2_4.x[4] - 1., d__4 *= d__4;
    l6_1.fx = (d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * (d__4 * d__4))
	     * .001;
    return 0;
labelL3:
    l4_4.gf[0] = (l2_4.x[0] - l2_4.x[1]) * 2. * .001;
    l4_4.gf[1] = -l4_4.gf[0];
    l4_4.gf[2] = (l2_4.x[2] - 1.) * 2. * .001;
/* Computing 3rd power */
    d__1 = l2_4.x[3] - 1.;
    l4_4.gf[3] = d__1 * (d__1 * d__1) * 4. * .001;
/* Computing 5th power */
    d__1 = l2_4.x[4] - 1., d__2 = d__1, d__1 *= d__1;
    l4_4.gf[4] = d__2 * (d__1 * d__1) * 6. * .001;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_4.x[0] + l2_4.x[1] + l2_4.x[2] + l2_4.x[3] * 4. - 7.;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = l2_4.x[2] + l2_4.x[4] * 5. - 6.;
    }
labelL5:
    return 0;
} /* tp49_ */


/* Subroutine */ int tp50_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__, j;
    static Real v1, v2, v3, v4;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 3;
    l1_1.nenl = 0;
    l2_4.x[0] = 35.;
    l2_4.x[1] = -31.;
    l2_4.x[2] = 11.;
    l2_4.x[3] = 5.;
    l2_4.x[4] = -5.;
    for (i__ = 1; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = false;
/* labelL6: */
	l12_4.lxu[i__ - 1] = false;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 5; ++j) {
/* L15: */
	    l5_9.gg[i__ + j * 3 - 4] = 0.;
	}
    }
    l5_9.gg[0] = 1.;
    l5_9.gg[3] = 2.;
    l5_9.gg[6] = 3.;
    l5_9.gg[4] = 1.;
    l5_9.gg[7] = 2.;
    l5_9.gg[10] = 3.;
    l5_9.gg[8] = 1.;
    l5_9.gg[11] = 2.;
    l5_9.gg[14] = 3.;
    l20_7.lex = true;
    l20_7.nex = 1;
    for (i__ = 1; i__ <= 5; ++i__) {
/* L30: */
	l20_7.xex[i__ - 1] = 1.;
    }
    l20_7.fex = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_4.x[0] - l2_4.x[1];
/* Computing 2nd power */
    d__2 = l2_4.x[1] - l2_4.x[2];
/* Computing 4th power */
    d__3 = l2_4.x[2] - l2_4.x[3], d__3 *= d__3;
/* Computing 4th power */
    d__4 = l2_4.x[3] - l2_4.x[4], d__4 *= d__4;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4;
    return 0;
labelL3:
    v1 = (l2_4.x[0] - l2_4.x[1]) * 2.;
    v2 = (l2_4.x[1] - l2_4.x[2]) * 2.;
/* Computing 3rd power */
    d__1 = l2_4.x[2] - l2_4.x[3];
    v3 = d__1 * (d__1 * d__1) * 4.;
/* Computing 3rd power */
    d__1 = l2_4.x[3] - l2_4.x[4];
    v4 = d__1 * (d__1 * d__1) * 4.;
    l4_4.gf[0] = v1;
    l4_4.gf[1] = -v1 + v2;
    l4_4.gf[2] = -v2 + v3;
    l4_4.gf[3] = -v3 + v4;
    l4_4.gf[4] = -v4;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_4.x[0] + l2_4.x[1] * 2. + l2_4.x[2] * 3. - 6.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_4.x[1] + l2_4.x[2] * 2. + l2_4.x[3] * 3. - 6.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = l2_4.x[2] + l2_4.x[3] * 2. + l2_4.x[4] * 3. - 6.;
    }
labelL5:
    return 0;
} /* tp50_ */


/* Subroutine */ int tp51_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 3;
    l1_1.nenl = 0;
    l2_4.x[0] = 2.5;
    l2_4.x[1] = .5;
    l2_4.x[2] = 2.;
    l2_4.x[3] = -1.;
    l2_4.x[4] = .5;
    for (i__ = 1; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = false;
/* labelL6: */
	l12_4.lxu[i__ - 1] = false;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 5; ++j) {
/* L15: */
	    l5_9.gg[i__ + j * 3 - 4] = 0.;
	}
    }
    l5_9.gg[0] = 1.;
    l5_9.gg[3] = 3.;
    l5_9.gg[7] = 1.;
    l5_9.gg[10] = 1.;
    l5_9.gg[13] = -2.;
    l5_9.gg[5] = 1.;
    l5_9.gg[14] = -1.;
    l20_7.lex = true;
    l20_7.nex = 1;
    for (i__ = 1; i__ <= 5; ++i__) {
/* L30: */
	l20_7.xex[i__ - 1] = 1.;
    }
    l20_7.fex = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_4.x[0] - l2_4.x[1];
/* Computing 2nd power */
    d__2 = l2_4.x[1] + l2_4.x[2] - 2.;
/* Computing 2nd power */
    d__3 = l2_4.x[3] - 1.;
/* Computing 2nd power */
    d__4 = l2_4.x[4] - 1.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4;
    return 0;
labelL3:
    l4_4.gf[0] = (l2_4.x[0] - l2_4.x[1]) * 2.;
    l4_4.gf[2] = (l2_4.x[1] + l2_4.x[2] - 2.) * 2.;
    l4_4.gf[1] = l4_4.gf[2] - l4_4.gf[0];
    l4_4.gf[3] = (l2_4.x[3] - 1.) * 2.;
    l4_4.gf[4] = (l2_4.x[4] - 1.) * 2.;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_4.x[0] + l2_4.x[1] * 3. - 4.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_4.x[2] + l2_4.x[3] - l2_4.x[4] * 2.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = l2_4.x[1] - l2_4.x[4];
    }
labelL5:
    return 0;
} /* tp51_ */


/* Subroutine */ int tp52_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 3;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 5; ++i__) {
	l2_4.x[i__ - 1] = 2.;
	l11_4.lxl[i__ - 1] = false;
/* labelL6: */
	l12_4.lxu[i__ - 1] = false;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 5; ++j) {
/* L15: */
	    l5_9.gg[i__ + j * 3 - 4] = 0.;
	}
    }
    l5_9.gg[0] = 1.;
    l5_9.gg[3] = 3.;
    l5_9.gg[7] = 1.;
    l5_9.gg[10] = 1.;
    l5_9.gg[13] = -2.;
    l5_9.gg[5] = 1.;
    l5_9.gg[14] = -1.;
    l20_7.lex = true;
    l20_7.nex = 1;
    l20_7.xex[0] = -.094555873925501438;
    l20_7.xex[1] = .03151862464183381;
    l20_7.xex[2] = .51575931232091687;
    l20_7.xex[3] = -.45272206303724927;
    l20_7.xex[4] = l20_7.xex[1];
    l20_7.fex = 5.3266475644699138;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_4.x[0] * 4. - l2_4.x[1];
/* Computing 2nd power */
    d__2 = l2_4.x[1] + l2_4.x[2] - 2.;
/* Computing 2nd power */
    d__3 = l2_4.x[3] - 1.;
/* Computing 2nd power */
    d__4 = l2_4.x[4] - 1.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4;
    return 0;
labelL3:
    l4_4.gf[0] = (l2_4.x[0] * 4. - l2_4.x[1]) * 8.;
    l4_4.gf[2] = (l2_4.x[1] + l2_4.x[2] - 2.) * 2.;
    l4_4.gf[1] = l4_4.gf[0] * -.25 + l4_4.gf[2];
    l4_4.gf[3] = (l2_4.x[3] - 1.) * 2.;
    l4_4.gf[4] = (l2_4.x[4] - 1.) * 2.;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_4.x[0] + l2_4.x[1] * 3.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_4.x[2] + l2_4.x[3] - l2_4.x[4] * 2.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = l2_4.x[1] - l2_4.x[4];
    }
labelL5:
    return 0;
} /* tp52_ */


/* Subroutine */ int tp53_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 3;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 5; ++i__) {
	l2_4.x[i__ - 1] = 2.;
	l11_4.lxl[i__ - 1] = true;
	l12_4.lxu[i__ - 1] = true;
	l13_4.xl[i__ - 1] = -10.;
/* labelL6: */
	l14_4.xu[i__ - 1] = 10.;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 5; ++j) {
/* L15: */
	    l5_9.gg[i__ + j * 3 - 4] = 0.;
	}
    }
    l5_9.gg[0] = 1.;
    l5_9.gg[3] = 3.;
    l5_9.gg[7] = 1.;
    l5_9.gg[10] = 1.;
    l5_9.gg[13] = -2.;
    l5_9.gg[5] = 1.;
    l5_9.gg[14] = -1.;
    l20_7.lex = true;
    l20_7.nex = 1;
    l20_7.xex[0] = -.76744186046511631;
    l20_7.xex[1] = .2558139534883721;
    l20_7.xex[2] = .62790697674418605;
    l20_7.xex[3] = -.11627906976744186;
    l20_7.xex[4] = .2558139534883721;
    l20_7.fex = 4.0930232558139537;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_4.x[0] - l2_4.x[1];
/* Computing 2nd power */
    d__2 = l2_4.x[1] + l2_4.x[2] - 2.;
/* Computing 2nd power */
    d__3 = l2_4.x[3] - 1.;
/* Computing 2nd power */
    d__4 = l2_4.x[4] - 1.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4;
    return 0;
labelL3:
    l4_4.gf[0] = (l2_4.x[0] - l2_4.x[1]) * 2.;
    l4_4.gf[2] = (l2_4.x[1] + l2_4.x[2] - 2.) * 2.;
    l4_4.gf[1] = l4_4.gf[2] - l4_4.gf[0];
    l4_4.gf[3] = (l2_4.x[3] - 1.) * 2.;
    l4_4.gf[4] = (l2_4.x[4] - 1.) * 2.;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_4.x[0] + l2_4.x[1] * 3.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_4.x[2] + l2_4.x[3] - l2_4.x[4] * 2.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = l2_4.x[1] - l2_4.x[4];
    }
labelL5:
    return 0;
} /* tp53_ */


/* Subroutine */ int tp54_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static int i__;
    static Real q, v1, v2, v3, v4, v5, v6, v7, v8, v9, dq[6];

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 6;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    l2_5.x[0] = 6e3;
    l2_5.x[1] = 1.5;
    l2_5.x[2] = 4e6;
    l2_5.x[3] = 2.;
    l2_5.x[4] = .003;
    l2_5.x[5] = 5e7;
    for (i__ = 1; i__ <= 6; ++i__) {
	l11_5.lxl[i__ - 1] = true;
/* labelL6: */
	l12_5.lxu[i__ - 1] = true;
    }
    l13_5.xl[0] = 0.;
    l13_5.xl[1] = -10.;
    l13_5.xl[2] = 0.;
    l13_5.xl[3] = 0.;
    l13_5.xl[4] = -1.;
    l13_5.xl[5] = 0.;
    l14_5.xu[0] = 2e4;
    l14_5.xu[1] = 10.;
    l14_5.xu[2] = 1e7;
    l14_5.xu[3] = 20.;
    l14_5.xu[4] = 1.;
    l14_5.xu[5] = 2e8;
    l5_3.gg[0] = 1.;
    l5_3.gg[1] = 4e3;
    for (i__ = 3; i__ <= 6; ++i__) {
/* L30: */
	l5_3.gg[i__ - 1] = (float)0.;
    }
    l20_4.lex = true;
    l20_4.nex = 1;
    l20_4.xex[0] = 13085.714285714286;
    l20_4.xex[1] = 1.1285714285714286;
    l20_4.xex[2] = 2e6;
    l20_4.xex[3] = 10.;
    l20_4.xex[4] = .001;
    l20_4.xex[5] = 1e8;
    l20_4.fex = -std::exp(-.096428571428571433);
    return 0;
labelL2:
    v1 = l2_5.x[0] - 1e4;
    v2 = l2_5.x[1] - 1.;
    v3 = l2_5.x[2] - 2e6;
    v4 = l2_5.x[3] - 10.;
    v5 = l2_5.x[4] - .001;
    v6 = l2_5.x[5] - 1e8;
    v7 = 1.0416666666666667;
    v8 = 2.0408163265306122e-14;
    v9 = 4.0816326530612245e-14;
/* Computing 2nd power */
    d__1 = v1;
/* Computing 2nd power */
    d__2 = v2;
/* Computing 2nd power */
    d__3 = l2_5.x[2] - 2e6;
/* Computing 2nd power */
    d__4 = l2_5.x[3] - 10.;
/* Computing 2nd power */
    d__5 = l2_5.x[4] - .001;
/* Computing 2nd power */
    d__6 = l2_5.x[5] - 1e8;
    q = (d__1 * d__1 * 1.5625e-8 + v1 * 5e-5 * v2 + d__2 * d__2) * v7 + d__3 *
	     d__3 * v8 + d__4 * d__4 * 4e-4 + d__5 * d__5 * 400. + d__6 * 
	    d__6 * 4e-18;
/*      Q = ((X(1)-1.0D6)**2/6.4D+7 + (X(1)-1.0D+4)*(X(2)-1.0D0)/2.0D4 */
/*     /       + (X(2)-1.0D0)**2)*(X(3)-2.0D6)**2/(0.96*4.9D13) */
/*     /       + (X(4)-1.0D1)**2/2.5D3 + (X(5)-1.0D-3)**2/2.5D-3 */
/*     /       + (X(6)-1.0D8)**2/2.5D17 */
    l6_1.fx = -std::exp(q * -.5);
    return 0;
labelL3:
    v1 = l2_5.x[0] - 1e4;
    v2 = l2_5.x[1] - 1.;
    v3 = l2_5.x[2] - 2e6;
    v4 = l2_5.x[3] - 10.;
    v5 = l2_5.x[4] - .001;
    v6 = l2_5.x[5] - 1e8;
    v7 = 1.0416666666666667;
    v8 = 2.0408163265306122e-14;
    v9 = 4.0816326530612245e-14;
/* Computing 2nd power */
    d__1 = v1;
/* Computing 2nd power */
    d__2 = v2;
/* Computing 2nd power */
    d__3 = v3;
/* Computing 2nd power */
    d__4 = v4;
/* Computing 2nd power */
    d__5 = v5;
/* Computing 2nd power */
    d__6 = v6;
    q = (d__1 * d__1 * 1.5625e-8 + v1 * 5e-5 * v2 + d__2 * d__2) * v7 + d__3 *
	     d__3 * v8 + d__4 * d__4 * 4e-4 + d__5 * d__5 * 400. + d__6 * 
	    d__6 * 4e-18;
    dq[0] = (v1 * 3.125e-8 + v2 * 5e-5) * v7;
    dq[1] = (v1 * 5e-5 + v2 * 2.) * v7;
    dq[2] = v3 * v9;
    dq[3] = v4 * 8e-4;
    dq[4] = v5 * 800.;
    dq[5] = v6 * 8e-18;
    for (i__ = 1; i__ <= 6; ++i__) {
/* L31: */
	l4_5.gf[i__ - 1] = std::exp(q * -.5) * .5 * dq[i__ - 1];
    }
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = l2_5.x[0] + l2_5.x[1] * 4e3 - 17600.;
    }
labelL5:
    return 0;
} /* tp54_ */


/* Subroutine */ int tp55_(int *mode)
{
    /* Local variables */
    static int i__, j;
    static Real v1;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 6;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 6;
    l1_1.nenl = 0;
    l2_5.x[0] = 1.;
    l2_5.x[1] = 2.;
    l2_5.x[2] = 0.;
    l2_5.x[3] = 0.;
    l2_5.x[4] = 0.;
    l2_5.x[5] = 2.;
    for (i__ = 1; i__ <= 6; ++i__) {
	l11_5.lxl[i__ - 1] = true;
	l12_5.lxu[i__ - 1] = false;
/* labelL6: */
	l13_5.xl[i__ - 1] = 0.;
    }
    l12_5.lxu[0] = true;
    l12_5.lxu[3] = true;
    l14_5.xu[0] = 1.;
    l14_5.xu[3] = 1.;
    for (i__ = 1; i__ <= 6; ++i__) {
	for (j = 1; j <= 6; ++j) {
/* L15: */
	    l5_10.gg[i__ + j * 6 - 7] = 0.;
	}
    }
    l4_5.gf[1] = 2.;
    l4_5.gf[2] = 0.;
    l4_5.gf[4] = 4.;
    l4_5.gf[5] = 0.;
    l5_10.gg[0] = 1.;
    l5_10.gg[6] = 2.;
    l5_10.gg[24] = 5.;
    l5_10.gg[1] = 1.;
    l5_10.gg[7] = 1.;
    l5_10.gg[13] = 1.;
    l5_10.gg[20] = 1.;
    l5_10.gg[26] = 1.;
    l5_10.gg[32] = 1.;
    l5_10.gg[3] = 1.;
    l5_10.gg[21] = 1.;
    l5_10.gg[10] = 1.;
    l5_10.gg[28] = 1.;
    l5_10.gg[17] = 1.;
    l5_10.gg[35] = 1.;
    l20_4.lex = true;
    l20_4.nex = 1;
    l20_4.xex[0] = 0.;
    l20_4.xex[1] = 1.3333333333333333;
    l20_4.xex[2] = 1.6666666666666667;
    l20_4.xex[3] = 1.;
    l20_4.xex[4] = .66666666666666663;
    l20_4.xex[5] = .33333333333333331;
    l20_4.fex = 6.333333333333333;
    return 0;
labelL2:
    l6_1.fx = l2_5.x[0] + l2_5.x[1] * 2. + l2_5.x[4] * 4. + std::exp(l2_5.x[0] * 
	    l2_5.x[3]);
    return 0;
labelL3:
    v1 = std::exp(l2_5.x[0] * l2_5.x[3]);
    l4_5.gf[0] = l2_5.x[3] * v1 + 1.;
    l4_5.gf[3] = l2_5.x[0] * v1;
    return 0;
labelL4:
    if (l9_6.index1[0]) {
	l3_5.g[0] = l2_5.x[0] + l2_5.x[1] * 2. + l2_5.x[4] * 5. - 6.;
    }
    if (l9_6.index1[1]) {
	l3_5.g[1] = l2_5.x[0] + l2_5.x[1] + l2_5.x[2] - 3.;
    }
    if (l9_6.index1[2]) {
	l3_5.g[2] = l2_5.x[3] + l2_5.x[4] + l2_5.x[5] - 2.;
    }
    if (l9_6.index1[3]) {
	l3_5.g[3] = l2_5.x[0] + l2_5.x[3] - 1.;
    }
    if (l9_6.index1[4]) {
	l3_5.g[4] = l2_5.x[1] + l2_5.x[4] - 2.;
    }
    if (l9_6.index1[5]) {
	l3_5.g[5] = l2_5.x[2] + l2_5.x[5] - 2.;
    }
labelL5:
    return 0;
} /* tp55_ */


/* Subroutine */ int tp56_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 7;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 4;
    l2_6.x[0] = 1.;
    l2_6.x[1] = 1.;
    l2_6.x[2] = 1.;
    for (i__ = 1; i__ <= 7; ++i__) {
	l11_6.lxl[i__ - 1] = false;
/* labelL6: */
	l12_6.lxu[i__ - 1] = false;
    }
    l2_6.x[3] = std::asin(std::sqrt(.23809523809523808));
    l2_6.x[4] = l2_6.x[3];
    l2_6.x[5] = l2_6.x[3];
    l2_6.x[6] = std::asin(std::sqrt(.69444444444444442));
    for (i__ = 4; i__ <= 7; ++i__) {
/* L30: */
	l4_6.gf[i__ - 1] = 0.;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	for (j = 1; j <= 7; ++j) {
/* L15: */
	    l5_11.gg[i__ + (j << 2) - 5] = 0.;
	}
    }
    l5_11.gg[0] = 1.;
    l5_11.gg[5] = 1.;
    l5_11.gg[10] = 1.;
    l5_11.gg[3] = 1.;
    l5_11.gg[7] = 2.;
    l5_11.gg[11] = 2.;
    l20_8.lex = true;
    l20_8.nex = -1;
    l20_8.xex[0] = 2.4;
    l20_8.xex[1] = 1.2;
    l20_8.xex[2] = 1.2;
    l20_8.xex[3] = std::asin(std::sqrt(.5714285714285714));
    l20_8.xex[4] = std::asin(std::sqrt(.2857142857142857));
    l20_8.xex[5] = l20_8.xex[4];
    l20_8.xex[6] = std::atan(1.) * 2.;
    l20_8.fex = -3.456;
    return 0;
labelL2:
    l6_1.fx = -l2_6.x[0] * l2_6.x[1] * l2_6.x[2];
    return 0;
labelL3:
    l4_6.gf[0] = -l2_6.x[1] * l2_6.x[2];
    l4_6.gf[1] = -l2_6.x[0] * l2_6.x[2];
    l4_6.gf[2] = -l2_6.x[0] * l2_6.x[1];
    return 0;
labelL4:
    if (l9_7.index1[0]) {
/* Computing 2nd power */
	d__1 = std::sin(l2_6.x[3]);
	l3_6.g[0] = l2_6.x[0] - d__1 * d__1 * 4.2;
    }
    if (l9_7.index1[1]) {
/* Computing 2nd power */
	d__1 = std::sin(l2_6.x[4]);
	l3_6.g[1] = l2_6.x[1] - d__1 * d__1 * 4.2;
    }
    if (l9_7.index1[2]) {
/* Computing 2nd power */
	d__1 = std::sin(l2_6.x[5]);
	l3_6.g[2] = l2_6.x[2] - d__1 * d__1 * 4.2;
    }
    if (l9_7.index1[3]) {
/* Computing 2nd power */
	d__1 = std::sin(l2_6.x[6]);
	l3_6.g[3] = l2_6.x[0] + l2_6.x[1] * 2. + l2_6.x[2] * 2. - d__1 * d__1 
		* 7.2;
    }
    return 0;
labelL5:
    if (l10_7.index2[0]) {
	l5_11.gg[12] = std::sin(l2_6.x[3]) * -8.4 * std::cos(l2_6.x[3]);
    }
    if (l10_7.index2[1]) {
	l5_11.gg[17] = std::sin(l2_6.x[4]) * -8.4 * std::cos(l2_6.x[4]);
    }
    if (l10_7.index2[2]) {
	l5_11.gg[22] = std::sin(l2_6.x[5]) * -8.4 * std::cos(l2_6.x[5]);
    }
    if (l10_7.index2[3]) {
	l5_11.gg[27] = std::sin(l2_6.x[6]) * -14.4 * std::cos(l2_6.x[6]);
    }
    return 0;
} /* tp56_ */


/* Subroutine */ int tp57_(int *mode)
{
    /* System generated locals */
    int i__1;
    Real d__1;

    /* Local variables */
    static Real a[44], b[44], f[44];
    static int i__, j;
    static Real s[2], t, v1, df[88]	/* was [44][2] */;

    if (*mode - 2 >= 0) {
	goto L18;
    } else {
	goto labelL1;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = .42;
    l2_1.x[1] = 5.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = true;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l13_1.xl[0] = .4;
    l13_1.xl[1] = -4.;
    l20_1.lex = false;
    l20_1.xex[0] = .419952674511;
    l20_1.xex[1] = 1.2848456293;
    l20_1.fex = 2.8459669721299998;
    return 0;
L18:
    for (i__ = 1; i__ <= 2; ++i__) {
	a[i__ - 1] = 8.;
	a[i__ + 15] = 18.;
	a[i__ + 29] = 28.;
	a[i__ + 34] = 32.;
	a[i__ + 37] = 36.;
	a[i__ + 39] = 38.;
	b[i__ - 1] = .49;
	b[i__ + 5] = .46;
	b[i__ + 10] = .43;
	b[i__ + 13] = .43;
	b[i__ + 17] = .42;
	b[i__ + 20] = .41;
	b[i__ + 24] = .4;
	b[i__ + 28] = .41;
	b[i__ + 35] = .4;
	b[i__ + 39] = .4;
/* labelL20: */
	b[i__ + 41] = .39;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	a[i__ + 9] = 14.;
	a[i__ + 12] = 16.;
	a[i__ + 17] = 20.;
	a[i__ + 20] = 22.;
	a[i__ + 23] = 24.;
	a[i__ + 26] = 26.;
	a[i__ + 31] = 30.;
/* L21: */
	b[i__ + 30] = .4;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	a[i__ + 1] = 10.;
/* L22: */
	a[i__ + 5] = 12.;
    }
    a[37] = 34.;
    a[42] = 40.;
    a[43] = 42.;
    b[2] = .48;
    b[3] = .47;
    b[4] = .48;
    b[5] = .47;
    b[8] = .45;
    b[9] = .43;
    b[10] = .45;
    b[13] = .44;
    b[16] = .46;
    b[17] = .45;
    b[20] = .43;
    b[23] = .4;
    b[24] = .42;
    b[27] = .41;
    b[28] = .4;
    b[34] = .38;
    b[35] = .41;
    b[38] = .41;
    b[39] = .38;
    if ((i__1 = *mode - 4) < 0) {
	goto L17;
    } else if (i__1 == 0) {
	goto labelL4;
    } else {
	goto labelL5;
    }
L17:
    for (i__ = 1; i__ <= 44; ++i__) {
/* L30: */
	f[i__ - 1] = b[i__ - 1] - l2_1.x[0] - (.49 - l2_1.x[0]) * std::exp(-l2_1.x[
		1] * (a[i__ - 1] - 8.));
    }
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
    }
labelL2:
    t = 0.;
    for (i__ = 1; i__ <= 44; ++i__) {
/* L19: */
/* Computing 2nd power */
	d__1 = f[i__ - 1];
	t += d__1 * d__1;
    }
    l6_1.fx = t * (float)100.;
    return 0;
labelL3:
    s[0] = 0.;
    s[1] = 0.;
    for (i__ = 1; i__ <= 44; ++i__) {
	v1 = std::exp(-l2_1.x[1] * (a[i__ - 1] - 8.));
	df[i__ - 1] = v1 - 1.;
	df[i__ + 43] = (a[i__ - 1] - 8.) * (.49 - l2_1.x[0]) * v1;
	for (j = 1; j <= 2; ++j) {
/* L32: */
	    s[j - 1] += f[i__ - 1] * 2. * df[i__ + j * 44 - 45];
	}
/* L31: */
    }
    l4_1.gf[0] = s[0] * (float)100.;
    l4_1.gf[1] = s[1] * (float)100.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = -l2_1.x[0] * l2_1.x[1] + l2_1.x[1] * .49 - .09;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_1.gg[0] = -l2_1.x[1];
    l5_1.gg[1] = -l2_1.x[0] + .49;
L7:
    return 0;
} /* tp57_ */


/* Subroutine */ int tp58_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 3;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -2.;
    l2_1.x[1] = 1.;
    l11_1.lxl[0] = true;
    l12_1.lxu[0] = true;
    l11_1.lxl[1] = false;
    l12_1.lxu[1] = false;
    l13_1.xl[0] = -.5;
    l14_1.xu[0] = .5;
    l5_3.gg[0] = -1.;
    l5_3.gg[4] = -1.;
    l20_1.lex = false;
    l20_1.xex[0] = -.786150483331;
    l20_1.xex[1] = .618034533851;
    l20_1.fex = 3.19033354957;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = d__1 * d__1 * 100. + d__3 * d__3;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (l2_1.x[1] - d__1 * d__1) * 200.;
    l4_1.gf[0] = (l2_1.x[0] * (l4_1.gf[1] - 1.) + 1.) * -2.;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[1];
	l3_3.g[0] = d__1 * d__1 - l2_1.x[0];
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_3.g[1] = d__1 * d__1 - l2_1.x[1];
    }
    if (l9_4.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_3.g[2] = d__1 * d__1 + d__2 * d__2 - 1.;
    }
    return 0;
labelL5:
    if (l10_4.index2[0]) {
	l5_3.gg[3] = l2_1.x[1] * 2.;
    }
    if (l10_4.index2[1]) {
	l5_3.gg[1] = l2_1.x[0] * 2.;
    }
    if (! l10_4.index2[2]) {
	goto labelL9;
    }
    l5_3.gg[2] = l2_1.x[0] * 2.;
    l5_3.gg[5] = l2_1.x[1] * 2.;
labelL9:
    return 0;
} /* tp58_ */


/* Subroutine */ int tp59_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;
    static Real x11, x12, x13, x14, x21, x22, x23, x24, xx12, xx21, 
	    xx31;
    static Real xx11;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 3;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 90.;
    l2_1.x[1] = 10.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = true;
	l12_1.lxu[i__ - 1] = true;
/* labelL6: */
	l13_1.xl[i__ - 1] = 0.;
    }
    l14_1.xu[0] = 75.;
    l14_1.xu[1] = 65.;
    l5_3.gg[4] = 1.;
    l5_3.gg[2] = -5.;
    l20_1.lex = false;
/*     XEX(1) = 0.463995762710D+2 */
/*     XEX(2) = 0.522196899513D+2 */
/*     FEX = -0.675456604292D+1 */
    l20_1.xex[0] = 13.5501042366;
    l20_1.xex[1] = 51.6601812877;
    l20_1.fex = -7.80422632408;
    return 0;
labelL2:
    x11 = l2_1.x[0];
    x12 = x11 * x11;
    x13 = x12 * x11;
    x14 = x13 * x11;
    x21 = l2_1.x[1];
    x22 = x21 * x21;
    x23 = x22 * x21;
    x24 = x23 * x21;
    l6_1.fx = x11 * 3.8112 - 75.196 - x12 * .12694 + x13 * .0020567 - x14 * 
	    1.0345e-5 + x21 * 6.8306 - x11 * .030234 * x21 + x12 * .00128134 *
	     x21 - x13 * 3.5256e-5 * x21 + x14 * 2.266e-7 * x21 - x22 * 
	    .25645 + x23 * .0034604 - x24 * 1.3514e-5 + 28.106 / (x21 + 1.) + 
	    x12 * 5.2375e-6 * x22 + x13 * 6.3e-8 * x22 - x13 * 7e-10 * x23 - 
	    x11 * 3.4054e-4 * x22 + x11 * 1.6638e-6 * x23 + std::exp(x11 * 5e-4 * 
	    x21) * 2.8673;
    return 0;
labelL3:
    x11 = l2_1.x[0];
    x12 = x11 * x11;
    x13 = x12 * x11;
    x14 = x13 * x11;
    x21 = l2_1.x[1];
    x22 = x21 * x21;
    x23 = x22 * x21;
    xx11 = x11 * x21;
    xx12 = x11 * x22;
    xx21 = x12 * x21;
    xx31 = x13 * x21;
    l4_1.gf[0] = 3.8112 - x11 * .25388 + x12 * .0061701 - x13 * 4.138e-5 - 
	    x21 * .030234 + xx11 * .00256268 - xx21 * 1.05768e-4 + xx31 * 
	    9.064e-7 + xx12 * 1.0475e-5 + x12 * 1.89e-7 * x22 - x12 * 2.1e-9 *
	     x23 - x22 * 3.4054e-4 + x23 * 1.6638e-6 + x21 * .00143365 * std::exp(
	    xx11 * 5e-4);
/* Computing 2nd power */
    d__1 = x21 + 1.;
    l4_1.gf[1] = 6.8306 - x11 * .030234 + x12 * .00128134 - x13 * 3.5256e-5 + 
	    x14 * 2.266e-7 - x21 * .5129 + x22 * .0103812 - x23 * 5.4056e-5 - 
	    28.106 / (d__1 * d__1) + xx21 * 1.0475e-5 + xx31 * 1.26e-7 - x13 *
	     2.1e-9 * x22 - xx11 * 6.8108e-4 + xx12 * 4.9914e-6 + x11 * 
	    .00143365 * std::exp(xx11 * 5e-4);
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_1.x[0] * l2_1.x[1] - 700.;
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_3.g[1] = l2_1.x[1] - d__1 * d__1 * .008;
    }
    if (l9_4.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_1.x[1] - 50.;
	l3_3.g[2] = d__1 * d__1 - (l2_1.x[0] - 55.) * 5.;
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto L7;
    }
    l5_3.gg[0] = l2_1.x[1];
    l5_3.gg[3] = l2_1.x[0];
L7:
    if (l10_4.index2[1]) {
	l5_3.gg[1] = l2_1.x[0] * -.016;
    }
    if (l10_4.index2[2]) {
	l5_3.gg[5] = (l2_1.x[1] - 50.) * 2.;
    }
    return 0;
} /* tp59_ */


/* Subroutine */ int tp60_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;
    static Real v1;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 2.;
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = true;
	l13_2.xl[i__ - 1] = -10.;
/* labelL6: */
	l14_2.xu[i__ - 1] = 10.;
    }
    l20_3.lex = false;
    l20_3.xex[0] = 1.10485902423;
    l20_3.xex[1] = 1.19667419413;
    l20_3.xex[2] = 1.53526225739;
    l20_3.fex = .0325682002513;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0] - 1.;
/* Computing 2nd power */
    d__2 = l2_2.x[0] - l2_2.x[1];
/* Computing 4th power */
    d__3 = l2_2.x[1] - l2_2.x[2], d__3 *= d__3;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return 0;
labelL3:
    v1 = (l2_2.x[0] - l2_2.x[1]) * 2.;
    l4_2.gf[0] = (l2_2.x[0] - 1.) * 2. + v1;
/* Computing 3rd power */
    d__1 = l2_2.x[1] - l2_2.x[2];
    l4_2.gf[2] = d__1 * (d__1 * d__1) * -4.;
    l4_2.gf[1] = -l4_2.gf[2] - v1;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[1];
/* Computing 4th power */
	d__2 = l2_2.x[2], d__2 *= d__2;
	l3_1.g[0] = l2_2.x[0] * (d__1 * d__1 + 1.) + d__2 * d__2 - 4. - std::sqrt( 2.) * 3.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
/* Computing 2nd power */
    d__1 = l2_2.x[1];
    l5_5.gg[0] = d__1 * d__1 + 1.;
    l5_5.gg[1] = l2_2.x[0] * 2. * l2_2.x[1];
/* Computing 3rd power */
    d__1 = l2_2.x[2];
    l5_5.gg[2] = d__1 * (d__1 * d__1) * 4.;
L7:
    return 0;
} /* tp60_ */


/* Subroutine */ int tp61_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 2;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 0.;
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l20_3.lex = false;
    l20_3.xex[0] = 5.32677015744;
    l20_3.xex[1] = -2.11899863998;
    l20_3.xex[2] = 3.21046423906;
    l20_3.fex = -143.646142201;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1];
/* Computing 2nd power */
    d__3 = l2_2.x[2];
    l6_1.fx = d__1 * d__1 * 4. + d__2 * d__2 * 2. + d__3 * d__3 * 2. - l2_2.x[
	    0] * 33. + l2_2.x[1] * 16. - l2_2.x[2] * 24.;
    return 0;
labelL3:
    l4_2.gf[0] = l2_2.x[0] * 8. - 33.;
    l4_2.gf[1] = l2_2.x[1] * 4. + 16.;
    l4_2.gf[2] = l2_2.x[2] * 4. - 24.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[1];
	l3_2.g[0] = l2_2.x[0] * 3. - d__1 * d__1 * 2. - 7.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_2.x[2];
	l3_2.g[1] = l2_2.x[0] * 4. - d__1 * d__1 - 11.;
    }
    return 0;
labelL5:
    if (l10_3.index2[0]) {
	l5_3.gg[2] = l2_2.x[1] * -4.;
    }
    if (l10_3.index2[1]) {
	l5_3.gg[5] = l2_2.x[2] * -2.;
    }
    l5_3.gg[0] = 3.;
    l5_3.gg[4] = 0.;
    l5_3.gg[1] = 4.;
    l5_3.gg[3] = 0.;
    return 0;
} /* tp61_ */


/* Subroutine */ int tp62_(int *mode)
{
    /* System generated locals */
    int i__1;
    Real d__1;

    /* Local variables */
    static int i__;
    static Real s, b1, b2, b3, c1, c2, c3, v1, v2, v3, v4, v5, v6, v7, 
	    rb1, rb2, rb3, rc1, rc2, rc3;

    if (*mode - 2 >= 0) {
	goto L17;
    } else {
	goto labelL1;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    l2_2.x[0] = .7;
    l2_2.x[1] = .2;
    l2_2.x[2] = .1;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = true;
	l13_2.xl[i__ - 1] = 0.;
	l14_2.xu[i__ - 1] = 1.;
/* labelL6: */
	l5_5.gg[i__ - 1] = 1.;
    }
    l20_3.lex = false;
    l20_3.xex[0] = .61781329821;
    l20_3.xex[1] = .328202155786;
    l20_3.xex[2] = .0539845460119;
    l20_3.fex = -26272.5144873;
    return 0;
L17:
    if ((i__1 = *mode - 4) < 0) {
	goto L18;
    } else if (i__1 == 0) {
	goto labelL4;
    } else {
	goto labelL5;
    }
L18:
    b3 = l2_2.x[2] + .03;
    c3 = l2_2.x[2] * .13 + .03;
    b2 = b3 + l2_2.x[1];
    c2 = b3 + l2_2.x[1] * .07;
    b1 = b2 + l2_2.x[0];
    c1 = b2 + l2_2.x[0] * .09;
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
    }
labelL2:
    v5 = b1 / c1;
    v6 = b2 / c2;
    v7 = b3 / c3;
    if (v5 <= 0. || v6 <= 0. || v7 <= 0.) {
	goto L7;
    }
    l6_1.fx = (std::log(v5) * 255. + std::log(v6) * 280. + std::log(v7) * 290.) * -32.174;
    return 0;
L7:
    s = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__1 = l2_2.x[i__ - 1] - 5.;
	s += d__1 * d__1;
    }
    l6_1.fx = s + 1e3 - 26700.;
    return 0;
labelL3:
    rb1 = 1. / b1;
    rb2 = 1. / b2;
    rb3 = 1. / b3;
    rc1 = 1. / c1;
    rc2 = 1. / c2;
    rc3 = 1. / c3;
    v1 = -8204.369999999999;
    v2 = -9008.7199999999993;
    v3 = -9330.4599999999991;
    v4 = v1 * (rb1 - rc1);
    l4_2.gf[0] = v1 * (rb1 - rc1 * .09);
    l4_2.gf[1] = v4 + v2 * (rb2 - rc2 * .07);
    l4_2.gf[2] = v4 + v2 * (rb2 - rc2) + v3 * (rb3 - rc3 * .13);
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = l2_2.x[0] + l2_2.x[1] + l2_2.x[2] - 1.;
    }
labelL5:
    return 0;
} /* tp62_ */


/* Subroutine */ int tp63_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 1;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 2.;
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = false;
/* labelL6: */
	l13_2.xl[i__ - 1] = 0.;
    }
    l5_3.gg[0] = 8.;
    l5_3.gg[2] = 14.;
    l5_3.gg[4] = 7.;
    l20_3.lex = false;
    l20_3.xex[0] = 3.51211841492;
    l20_3.xex[1] = .216988174172;
    l20_3.xex[2] = 3.55217403459;
    l20_3.fex = 961.715172127;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1];
/* Computing 2nd power */
    d__3 = l2_2.x[2];
    l6_1.fx = 1e3 - d__1 * d__1 - d__2 * d__2 * 2. - d__3 * d__3 - l2_2.x[0] *
	     l2_2.x[1] - l2_2.x[0] * l2_2.x[2];
    return 0;
labelL3:
    l4_2.gf[0] = l2_2.x[0] * -2. - l2_2.x[1] - l2_2.x[2];
    l4_2.gf[1] = l2_2.x[1] * -4. - l2_2.x[0];
    l4_2.gf[2] = l2_2.x[2] * -2. - l2_2.x[0];
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_2.x[0] * 8. + l2_2.x[1] * 14. + l2_2.x[2] * 7. - 56.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[1];
/* Computing 2nd power */
	d__3 = l2_2.x[2];
	l3_2.g[1] = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 - 25.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    l5_3.gg[1] = l2_2.x[0] * 2.;
    l5_3.gg[3] = l2_2.x[1] * 2.;
    l5_3.gg[5] = l2_2.x[2] * 2.;
L8:
    return 0;
} /* tp63_ */


/* Subroutine */ int tp64_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 1.;
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = false;
/* labelL6: */
	l13_2.xl[i__ - 1] = 1e-5;
    }
    l20_3.lex = false;
    l20_3.xex[0] = 108.734717597;
    l20_3.xex[1] = 85.1261394257;
    l20_3.xex[2] = 204.324707858;
    l20_3.fex = 6299.84242821;
    return 0;
labelL2:
    l6_1.fx = l2_2.x[0] * 5. + 5e4 / l2_2.x[0] + l2_2.x[1] * 20. + 7.2e4 / 
	    l2_2.x[1] + l2_2.x[2] * 10. + 1.44e5 / l2_2.x[2];
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
    l4_2.gf[0] = 5. - 5e4 / (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_2.x[1];
    l4_2.gf[1] = 20. - 7.2e4 / (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_2.x[2];
    l4_2.gf[2] = 10. - 1.44e5 / (d__1 * d__1);
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = 1. - 4. / l2_2.x[0] - 32. / l2_2.x[1] - 120. / l2_2.x[2];
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
/* Computing 2nd power */
    d__1 = l2_2.x[0];
    l5_5.gg[0] = 4. / (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_2.x[1];
    l5_5.gg[1] = 32. / (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_2.x[2];
    l5_5.gg[2] = 120. / (d__1 * d__1);
L7:
    return 0;
} /* tp64_ */


/* Subroutine */ int tp65_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;
    static Real v1, v2;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_2.x[0] = -5.;
    l2_2.x[1] = 5.;
    l2_2.x[2] = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = true;
/* labelL6: */
	l12_2.lxu[i__ - 1] = true;
    }
    l13_2.xl[0] = -4.5;
    l13_2.xl[1] = -4.5;
    l13_2.xl[2] = -5.;
    l14_2.xu[0] = 4.5;
    l14_2.xu[1] = 4.5;
    l14_2.xu[2] = 5.;
    l20_3.lex = false;
    l20_3.xex[0] = 3.65046182158;
    l20_3.xex[1] = 3.6504616894;
    l20_3.xex[2] = 4.62041750754;
    l20_3.fex = .953528856757;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0] - l2_2.x[1];
/* Computing 2nd power */
    d__2 = (l2_2.x[0] + l2_2.x[1] - 10.) / 3.;
/* Computing 2nd power */
    d__3 = l2_2.x[2] - 5.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return 0;
labelL3:
    v1 = (l2_2.x[0] - l2_2.x[1]) * 2.;
    v2 = (l2_2.x[0] + l2_2.x[1] - 10.) * 2. / 9.;
    l4_2.gf[0] = v1 + v2;
    l4_2.gf[1] = -v1 + v2;
    l4_2.gf[2] = (l2_2.x[2] - 5.) * 2.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[1];
/* Computing 2nd power */
	d__3 = l2_2.x[2];
	l3_1.g[0] = 48. - d__1 * d__1 - d__2 * d__2 - d__3 * d__3;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_5.gg[0] = l2_2.x[0] * -2.;
    l5_5.gg[1] = l2_2.x[1] * -2.;
    l5_5.gg[2] = l2_2.x[2] * -2.;
L7:
    return 0;
} /* tp65_ */


/* Subroutine */ int tp66_(int *mode)
{
    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_2.x[0] = 0.;
    l2_2.x[1] = 1.05;
    l2_2.x[2] = 2.9;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = true;
	l13_2.xl[i__ - 1] = 0.;
/* labelL6: */
	l14_2.xu[i__ - 1] = 100.;
    }
    l14_2.xu[2] = 10.;
    l4_2.gf[0] = -.8;
    l4_2.gf[1] = 0.;
    l4_2.gf[2] = .2;
    l5_3.gg[2] = 1.;
    l5_3.gg[4] = 0.;
    l5_3.gg[1] = 0.;
    l5_3.gg[5] = 1.;
    l20_3.lex = false;
    l20_3.xex[0] = .184126487951;
    l20_3.xex[1] = 1.20216787321;
    l20_3.xex[2] = 3.32732232258;
    l20_3.fex = .518163274159;
    return 0;
labelL2:
    l6_1.fx = l2_2.x[2] * .2 - l2_2.x[0] * .8;
labelL3:
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_2.x[1] - std::exp(l2_2.x[0]);
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = l2_2.x[2] - std::exp(l2_2.x[1]);
    }
    return 0;
labelL5:
    if (l10_3.index2[0]) {
	l5_3.gg[0] = -std::exp(l2_2.x[0]);
    }
    if (l10_3.index2[1]) {
	l5_3.gg[3] = -std::exp(l2_2.x[1]);
    }
    return 0;
} /* tp66_ */


/* Subroutine */ int tp67_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real a[3];
    static int i__, j;
    static Real y[8], v1, v2, v3, dy[24]	/* was [8][3] */, rx, y2c, 
	    y4c, dy2c[3], dy4c[3];

    if (*mode - 2 >= 0) {
	goto L17;
    } else {
	goto labelL1;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 14;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_2.x[0] = 1745.;
    l2_2.x[1] = 1.2e4;
    l2_2.x[2] = 110.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = true;
/* labelL6: */
	l13_2.xl[i__ - 1] = 1e-5;
    }
    l14_2.xu[0] = 2e3;
    l14_2.xu[1] = 1.6e4;
    l14_2.xu[2] = 120.;
    l20_3.lex = false;
    l20_3.xex[0] = 1728.37128614;
    l20_3.xex[1] = 1.6e4;
    l20_3.xex[2] = 98.1415140238;
    l20_3.fex = -1162.03650728;
    return 0;
L17:
    rx = 1. / l2_2.x[0];
    y[1] = l2_2.x[0] * 1.6;
    dy[1] = 1.6;
    dy[9] = 0.;
    dy[17] = 0.;
L100:
    y[2] = y[1] * 1.22 - l2_2.x[0];
    dy[2] = dy[1] * 1.22 - 1.;
    dy[10] = dy[9] * 1.22;
    dy[18] = dy[17] * 1.22;
    y[5] = (l2_2.x[1] + y[2]) * rx;
/* Computing 2nd power */
    d__1 = rx;
    dy[5] = (l2_2.x[0] * dy[2] - l2_2.x[1] - y[2]) * (d__1 * d__1);
    dy[13] = (dy[10] + 1.) * rx;
    dy[21] = dy[2] * rx;
    v1 = l2_2.x[0] * .01 * (13.167 - y[5] * 1.3333999999999999);
    v2 = ((13.167 - y[5] * .6667) * y[5] + 112.) * .01;
    y2c = l2_2.x[0] * v2;
    dy2c[0] = v2 + v1 * dy[5];
    dy2c[1] = v1 * dy[13];
    dy2c[2] = v1 * dy[21];
    if ((d__1 = y2c - y[1], std::abs(d__1)) - .001 <= 0.) {
	goto L102;
    } else {
	goto L101;
    }
L101:
    for (i__ = 1; i__ <= 3; ++i__) {
/* L103: */
	dy[(i__ << 3) - 7] = dy2c[i__ - 1];
    }
    y[1] = y2c;
    goto L100;
L102:
    y[3] = 93.;
    for (i__ = 1; i__ <= 3; ++i__) {
/* L104: */
	dy[(i__ << 3) - 5] = 0.;
    }
L105:
/* Computing 2nd power */
    d__1 = y[5];
    y[4] = y[5] * 1.098 + 86.35 - d__1 * d__1 * .038 + (y[3] - 89.) * .325;
    y[7] = y[4] * 3. - 133.;
    y[6] = 35.82 - y[7] * .222;
    for (i__ = 1; i__ <= 3; ++i__) {
	dy[(i__ << 3) - 4] = dy[(i__ << 3) - 3] * 1.098 - y[5] * .076 * dy[(
		i__ << 3) - 3] + dy[(i__ << 3) - 5] * .325;
	dy[(i__ << 3) - 1] = dy[(i__ << 3) - 4] * 3.;
/* L106: */
	dy[(i__ << 3) - 2] = dy[(i__ << 3) - 1] * -.222;
    }
    v3 = 1. / (y[1] * y[6] + l2_2.x[2] * 1e3);
    y4c = l2_2.x[2] * 9.8e4 * v3;
    for (i__ = 1; i__ <= 2; ++i__) {
/* L107: */
/* Computing 2nd power */
	d__1 = y[1] * y[6] + l2_2.x[2] * 1e3;
	dy4c[i__ - 1] = l2_2.x[2] * -9.8e4 * (y[1] * dy[(i__ << 3) - 2] + y[6]
		 * dy[(i__ << 3) - 7]) / (d__1 * d__1);
    }
/* Computing 2nd power */
    d__1 = v3;
    dy4c[2] = (y[1] * y[6] - l2_2.x[2] * (y[1] * dy[22] + y[6] * dy[17])) * 
	    9.8e4 * (d__1 * d__1);
    if ((d__1 = y4c - y[3], std::abs(d__1)) - 1e-4 <= 0.) {
	goto L109;
    } else {
	goto L108;
    }
L108:
    y[3] = y4c;
    for (i__ = 1; i__ <= 3; ++i__) {
/* L110: */
	dy[(i__ << 3) - 5] = dy4c[i__ - 1];
    }
    goto L105;
L109:
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL2:
    l6_1.fx = -(y[1] * .063 * y[4] - l2_2.x[0] * 5.04 - y[2] * 3.36 - l2_2.x[
	    1] * .035 - l2_2.x[2] * 10.);
    return 0;
labelL3:
    for (i__ = 1; i__ <= 3; ++i__) {
/* L120: */
	a[i__ - 1] = (dy[(i__ << 3) - 7] * y[4] + dy[(i__ << 3) - 4] * y[1]) *
		 -.063 + dy[(i__ << 3) - 6] * 3.36;
    }
    l4_2.gf[0] = a[0] + 5.04;
    l4_2.gf[1] = a[1] + .035;
    l4_2.gf[2] = a[2] + 10.;
    return 0;
labelL4:
    if (l9_8.index1[0]) {
	l3_7.g[0] = y[1];
    }
    if (l9_8.index1[1]) {
	l3_7.g[1] = y[2];
    }
    if (l9_8.index1[2]) {
	l3_7.g[2] = y[3] - 85.;
    }
    if (l9_8.index1[3]) {
	l3_7.g[3] = y[4] - 90.;
    }
    if (l9_8.index1[4]) {
	l3_7.g[4] = y[5] - 3.;
    }
    if (l9_8.index1[5]) {
	l3_7.g[5] = y[6] - .01;
    }
    if (l9_8.index1[6]) {
	l3_7.g[6] = y[7] - 145.;
    }
    if (l9_8.index1[7]) {
	l3_7.g[7] = 5e3 - y[1];
    }
    if (l9_8.index1[8]) {
	l3_7.g[8] = 2e3 - y[2];
    }
    if (l9_8.index1[9]) {
	l3_7.g[9] = 93. - y[3];
    }
    if (l9_8.index1[10]) {
	l3_7.g[10] = 95. - y[4];
    }
    if (l9_8.index1[11]) {
	l3_7.g[11] = 12. - y[5];
    }
    if (l9_8.index1[12]) {
	l3_7.g[12] = 4. - y[6];
    }
    if (l9_8.index1[13]) {
	l3_7.g[13] = 162. - y[7];
    }
    return 0;
labelL5:
    for (j = 1; j <= 7; ++j) {
	if (! l10_8.index2[j - 1]) {
	    goto L131;
	}
	for (i__ = 1; i__ <= 3; ++i__) {
/* L130: */
	    l5_12.gg[j + i__ * 14 - 15] = dy[j + 1 + (i__ << 3) - 9];
	}
L131:
	if (! l10_8.index2[j + 6]) {
	    goto L133;
	}
	for (i__ = 1; i__ <= 3; ++i__) {
/* L132: */
	    l5_12.gg[j + 7 + i__ * 14 - 15] = -dy[j + 1 + (i__ << 3) - 9];
	}
L133:
	;
    }
    return 0;
} /* tp67_ */


/* Subroutine */ int tp68_0_(int n__, int  *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real a[2], b[2], d__[2];
    static int i__;
    static Real z__[2], h1, h2, h3, v1, v2, v3, v4, v5;
    int mdnord_(Real*, Real*);
    static int kn1;

    switch(n__) {
	case 1: goto L_tp69;
	}

    kn1 = 1;
    l20_6.xex[0] = .0678587452312;
    l20_6.xex[1] = 3.64617174165;
    l20_6.xex[2] = 2.66175189694e-4;
    l20_6.xex[3] = .894862212037;
    l20_6.fex = -.920425020704;
    goto labelL9;

L_tp69:
    kn1 = 2;
    l20_6.xex[0] = .029371418083;
    l20_6.xex[1] = 1.19025343488;
    l20_6.xex[2] = .233946796758;
    l20_6.xex[3] = .791667815694;
    l20_6.fex = -956.712887064;
labelL9:
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 2;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_3.x[i__ - 1] = 1.;
	l2_3.x[i__ + 1] = 1.;
	l11_3.lxl[i__ - 1] = true;
	l11_3.lxl[i__ + 1] = true;
	l12_3.lxu[i__ - 1] = true;
	l12_3.lxu[i__ + 1] = true;
	l13_3.xl[i__ + 1] = 0.;
	l14_3.xu[i__ - 1] = 100.;
/* labelL6: */
	l14_3.xu[i__ + 1] = 2.;
    }
    l13_3.xl[0] = 1e-4;
    l13_3.xl[1] = 0.;
    a[0] = 1e-4;
    a[1] = .1;
    b[0] = 1.;
    b[1] = 1e3;
    d__[0] = 1.;
    d__[1] = 1.;
    z__[0] = 24.;
    z__[1] = 4.;
    l5_6.gg[0] = 0.;
    l5_6.gg[6] = 0.;
    l5_6.gg[4] = 1.;
    l5_6.gg[1] = 0.;
    l5_6.gg[5] = 0.;
    l5_6.gg[7] = 1.;
    l4_3.gf[1] = 0.;
    l20_6.lex = false;
    return 0;
labelL2:
    v1 = std::exp(l2_3.x[0]) - 1.;
    l6_1.fx = (a[kn1 - 1] * z__[kn1 - 1] - l2_3.x[3] * (b[kn1 - 1] * v1 - 
	    l2_3.x[2]) / (v1 + l2_3.x[3])) / l2_3.x[0];
    return 0;
labelL3:
    v1 = std::exp(l2_3.x[0]);
    v2 = v1 - 1.;
    v3 = 1. / (v2 + l2_3.x[3]);
    v4 = 1. / l2_3.x[0];
    v5 = (b[kn1 - 1] * v2 - l2_3.x[2]) * v4;
    l4_3.gf[0] = -((v1 * (l2_3.x[3] * b[kn1 - 1] + l2_3.x[2]) * v3 - v5) * 
	    l2_3.x[3] * v3 + a[kn1 - 1] * z__[kn1 - 1] * v4) * v4;
    l4_3.gf[2] = l2_3.x[3] * v4 * v3;
/* Computing 2nd power */
    d__1 = v3;
    l4_3.gf[3] = -v5 * v2 * (d__1 * d__1);
    return 0;
labelL4:
    if (! l9_3.index1[0]) {
	goto L30;
    }
    d__1 = -l2_3.x[1];
    mdnord_(&d__1, &h1);
    l3_2.g[0] = l2_3.x[2] - h1 * 2.;
L30:
    if (! l9_3.index1[1]) {
	return 0;
    }
    d__1 = -l2_3.x[1] + std::sqrt(z__[kn1 - 1]);
    mdnord_(&d__1, &h1);
    d__1 = -l2_3.x[1] - std::sqrt(z__[kn1 - 1]);
    mdnord_(&d__1, &h2);
    l3_2.g[1] = l2_3.x[3] - h1 - h2;
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
/* Computing 2nd power */
    d__1 = l2_3.x[1];
    l5_6.gg[2] = std::exp(d__1 * d__1 * -.5) * 2. / std::sqrt(std::atan(1.) * 8.);
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    h1 = -l2_3.x[1] - d__[kn1 - 1] * std::sqrt(z__[kn1 - 1]);
    h2 = -l2_3.x[1] + d__[kn1 - 1] * std::sqrt(z__[kn1 - 1]);
    h3 = 1. / std::sqrt(std::atan(1.) * 8.);
/* Computing 2nd power */
    d__1 = h1;
/* Computing 2nd power */
    d__2 = h2;
    l5_6.gg[3] = (std::exp(d__1 * d__1 * -.5) + std::exp(d__2 * d__2 * -.5)) * h3;
L8:
    return 0;
} /* tp68_ */

/* Subroutine */ int tp68_(int *mode)
{
    return tp68_0_(0, mode);
    }

/* Subroutine */ int tp69_(int *mode)
{
    return tp68_0_(1, mode);
    }


/* Subroutine */ int tp70_(int *mode)
{
    /* System generated locals */
    int i__1;
    Real d__1;

    /* Local variables */
    static Real b, c__[19], f[19];
    static int i__, j;
    static Real s, t, h1, h2, h3, h4, h5, h6, h7, h8, h9, u1[19], u2[19]
	    , v3[19], v4[19], z1, z2, v8[19], v9[19], z5, z6, z7, h10, df[76]	
	    /* was [19][4] */, h30, h40, h11, h12, h13, h14, h15, h16, h17, 
	    h18, v10, v11, v12, h19, h20, yc[19], h21, h22, h23, yo[19];
    static bool log__;
    static Real sum;

    log__ = false;
    if (*mode - 2 >= 0) {
	goto L17;
    } else {
	goto labelL1;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_3.x[0] = 2.;
    l2_3.x[1] = 4.;
    l2_3.x[2] = .04;
    l2_3.x[3] = 2.;
    for (i__ = 1; i__ <= 4; ++i__) {
	l12_3.lxu[i__ - 1] = true;
	l11_3.lxl[i__ - 1] = true;
	l13_3.xl[i__ - 1] = 1e-5;
/* labelL6: */
	l14_3.xu[i__ - 1] = 100.;
    }
    l14_3.xu[2] = 1.;
    l5_2.gg[0] = 0.;
    l5_2.gg[1] = 0.;
    l20_6.xex[0] = 12.2769537557;
    l20_6.xex[1] = 4.63178796886;
    l20_6.xex[2] = .312862470961;
    l20_6.xex[3] = 2.0292903157;
    l20_6.fex = .00749846356143;
    l20_6.lex = false;
    return 0;
L17:
    c__[0] = .1;
    for (i__ = 2; i__ <= 19; ++i__) {
/* labelL20: */
	c__[i__ - 1] = (Real) (i__ - 1);
    }
    yo[0] = .00189;
    yo[1] = .1038;
    yo[2] = .268;
    yo[3] = .506;
    yo[4] = .577;
    yo[5] = .604;
    yo[6] = .725;
    yo[7] = .898;
    yo[8] = .947;
    yo[9] = .845;
    yo[10] = .702;
    yo[11] = .528;
    yo[12] = .385;
    yo[13] = .257;
    yo[14] = .159;
    yo[15] = .0869;
    yo[16] = .0453;
    yo[17] = .01509;
    yo[18] = .00189;
    if ((i__1 = *mode - 4) < 0) {
	goto L18;
    } else if (i__1 == 0) {
	goto labelL4;
    } else {
	goto labelL5;
    }
L18:
    b = l2_3.x[2] + (1. - l2_3.x[2]) * l2_3.x[3];
    h1 = l2_3.x[0] - 1.;
    h2 = l2_3.x[1] - 1.;
    h3 = .13058239749281797;
    h5 = b * h3;
    h4 = h5 / l2_3.x[3];
    h6 = l2_3.x[0] * 12. / (l2_3.x[0] * 12. + 1.);
    h7 = l2_3.x[1] * 12. / (l2_3.x[1] * 12. + 1.);
    h30 = 0.;
    h40 = 0.;
    v10 = l2_3.x[1] / 6.2832;
    v11 = b / l2_3.x[3];
    v12 = l2_3.x[0] / 6.2832;
    if (b < 0. || v10 < 0. || v11 < 0. || v12 < 0.) {
	log__ = true;
    }
    if (log__ && *mode == 2) {
	goto L8;
    }
    if (log__ && *mode == 3) {
	goto labelL9;
    }
    z1 = l2_3.x[2] * pow_dd(&b, &l2_3.x[1]);
    z2 = pow_dd(&v10, &c_b590);
    z5 = 1. - l2_3.x[2];
    z6 = pow_dd(&v11, l2_3.x);
    z7 = pow_dd(&v12, &c_b590);
    for (i__ = 1; i__ <= 19; ++i__) {
	d__1 = c__[i__ - 1] * h3;
	v3[i__ - 1] = pow_dd(&d__1, &h2);
	v4[i__ - 1] = std::exp(l2_3.x[1] * (1. - c__[i__ - 1] * h5));
	d__1 = c__[i__ - 1] * h3;
	v8[i__ - 1] = pow_dd(&d__1, &h1);
	v9[i__ - 1] = std::exp(l2_3.x[0] * (1. - c__[i__ - 1] * h4));
	u1[i__ - 1] = z1 * z2 * v3[i__ - 1] * v4[i__ - 1] * h7;
	u2[i__ - 1] = z5 * z6 * z7 * v8[i__ - 1] * v9[i__ - 1] * h6;
	yc[i__ - 1] = u1[i__ - 1] + u2[i__ - 1];
/* L30: */
	f[i__ - 1] = yc[i__ - 1] - yo[i__ - 1];
    }
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
    }
labelL2:
    t = 0.;
    for (i__ = 1; i__ <= 19; ++i__) {
/* L31: */
/* Computing 2nd power */
	d__1 = f[i__ - 1];
	t += d__1 * d__1;
    }
    l6_1.fx = t;
    return 0;
labelL3:
    h8 = l2_3.x[3] - 1.;
    h9 = l2_3.x[2] - 1.;
    h10 = 1. / b;
    h11 = 1. / l2_3.x[3];
    h12 = (1. / (l2_3.x[0] * 12. + 1.) + .5) / l2_3.x[0] + 1.;
    h13 = (1. / (l2_3.x[1] * 12. + 1.) + .5) / l2_3.x[1] + 1.;
    h16 = l2_3.x[1] * h8;
    h17 = l2_3.x[0] * h8;
    h18 = 1. / l2_3.x[2] - h16 * h10;
    h19 = 1. / h9 - h17 * h10;
    h16 *= h3;
    h17 = h17 * h11 * h3;
    h20 = l2_3.x[1] * h9;
/* Computing 2nd power */
    d__1 = h11;
    h21 = l2_3.x[0] * l2_3.x[2] * (d__1 * d__1);
    h22 = h20 * h3;
    h23 = h21 * h3;
    h20 *= h10;
    h21 = h21 * h10 * l2_3.x[3];
    for (j = 1; j <= 19; ++j) {
	h14 = h4 * c__[j - 1];
	h15 = h5 * c__[j - 1];
	if (h14 > 0.) {
	    goto labelL10;
	}
	l4_3.gf[0] = (l2_3.x[0] - 5.) * 2.;
	h30 = 1.;
	goto labelL11;
labelL10:
	df[j - 1] = u2[j - 1] * (h12 - h14 + std::log(h14));
labelL11:
	if (h15 > 0.) {
	    goto labelL12;
	}
	l4_3.gf[1] = (l2_3.x[1] - 5.) * 2.;
	h40 = 1.;
	goto labelL13;
labelL12:
	df[j + 18] = u1[j - 1] * (h13 - h15 + std::log(h15));
labelL13:
	df[j + 37] = u1[j - 1] * (h18 + c__[j - 1] * h16) + u2[j - 1] * (h19 
		+ c__[j - 1] * h17);
/* L33: */
	df[j + 56] = u1[j - 1] * (h22 * c__[j - 1] - h20) + u2[j - 1] * (h23 *
		 c__[j - 1] - h21);
    }
    if (h30 == 1. || h40 == 1.) {
	goto labelL14;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	s = 0.;
	for (j = 1; j <= 19; ++j) {
	    s += f[j - 1] * 2. * df[j + i__ * 19 - 20];
/* L36: */
	}
	l4_3.gf[i__ - 1] = s;
/* L37: */
    }
    return 0;
labelL14:
    for (i__ = 1; i__ <= 4; ++i__) {
	if (i__ == 1 && h30 == 1.) {
	    goto L38;
	}
	if (i__ == 2 && h40 == (float).1) {
	    goto L38;
	}
	s = 0.;
	for (j = 1; j <= 19; ++j) {
/* L39: */
	    s += f[j - 1] * 2. * df[j + i__ * 19 - 20];
	}
	l4_3.gf[i__ - 1] = s;
L38:
	;
    }
    h30 = 0.;
    h40 = 0.;
    return 0;
L8:
    log__ = false;
    sum = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
/* L40: */
/* Computing 2nd power */
	d__1 = l2_3.x[i__ - 1] - 5.;
	sum += d__1 * d__1;
    }
    l6_1.fx = sum + 1e3 + .0075;
    return 0;
labelL9:
    log__ = false;
    for (i__ = 1; i__ <= 4; ++i__) {
/* L41: */
	l4_3.gf[i__ - 1] = (l2_3.x[i__ - 1] - 5.) * 2.;
    }
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = l2_3.x[2] + (1. - l2_3.x[2]) * l2_3.x[3];
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_2.gg[2] = 1. - l2_3.x[3];
    l5_2.gg[3] = -l2_3.x[2] + 1.;
L7:
    return 0;
} /* tp70_ */


/* Subroutine */ int tp71_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 4; ++i__) {
	l11_3.lxl[i__ - 1] = true;
	l12_3.lxu[i__ - 1] = true;
	l13_3.xl[i__ - 1] = 1.;
/* labelL6: */
	l14_3.xu[i__ - 1] = 5.;
    }
    l2_3.x[0] = 1.;
    l2_3.x[1] = 5.;
    l2_3.x[2] = 5.;
    l2_3.x[3] = 1.;
    l20_6.xex[0] = 1.;
    l20_6.xex[1] = 4.74299937545;
    l20_6.xex[2] = 3.82115032617;
    l20_6.xex[3] = 1.37940824585;
    l20_6.fex = 17.0140172895;
    l20_6.lex = false;
    return 0;
labelL2:
    l6_1.fx = l2_3.x[0] * l2_3.x[3] * (l2_3.x[0] + l2_3.x[1] + l2_3.x[2]) + 
	    l2_3.x[2];
    return 0;
labelL3:
    l4_3.gf[0] = l2_3.x[3] * (l2_3.x[0] * 2. + l2_3.x[1] + l2_3.x[2]);
    l4_3.gf[1] = l2_3.x[0] * l2_3.x[3];
    l4_3.gf[2] = l4_3.gf[1] + 1.;
    l4_3.gf[3] = l2_3.x[0] * (l2_3.x[0] + l2_3.x[1] + l2_3.x[2]);
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = (l2_3.x[0] * l2_3.x[1] * l2_3.x[2] * l2_3.x[3] - 25.) / 
		25.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[1];
/* Computing 2nd power */
	d__3 = l2_3.x[2];
/* Computing 2nd power */
	d__4 = l2_3.x[3];
	l3_2.g[1] = (d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4 - 
		40.) / 40.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_6.gg[0] = l2_3.x[1] * l2_3.x[2] * l2_3.x[3] / 25.;
    l5_6.gg[2] = l2_3.x[0] * l2_3.x[2] * l2_3.x[3] / 25.;
    l5_6.gg[4] = l2_3.x[0] * l2_3.x[1] * l2_3.x[3] / 25.;
    l5_6.gg[6] = l2_3.x[0] * l2_3.x[1] * l2_3.x[2] / 25.;
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
/* labelL20: */
	l5_6.gg[(i__ << 1) - 1] = l2_3.x[i__ - 1] * 2. / 40.;
    }
L8:
    return 0;
} /* tp71_ */


/* Subroutine */ int tp72_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 1.;
	l13_3.xl[i__ - 1] = .001;
	l14_3.xu[i__ - 1] = (5. - (Real) i__) * 1e5;
	l11_3.lxl[i__ - 1] = true;
/* labelL6: */
	l12_3.lxu[i__ - 1] = true;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
/* L30: */
	l4_3.gf[i__ - 1] = 1.;
    }
    l20_6.lex = false;
    l20_6.xex[0] = 193.407050141;
    l20_6.xex[1] = 179.547504555;
    l20_6.xex[2] = 185.018587841;
    l20_6.xex[3] = 168.706233485;
    l20_6.fex = 727.679376021;
    return 0;
labelL2:
    l6_1.fx = l2_3.x[0] + 1. + l2_3.x[1] + l2_3.x[2] + l2_3.x[3];
labelL3:
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = -4. / l2_3.x[0] - 2.25 / l2_3.x[1] - 1. / l2_3.x[2] - .25 
		/ l2_3.x[3] + .0401;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = -.16 / l2_3.x[0] - .36 / l2_3.x[1] - (1. / l2_3.x[2] + 1. 
		/ l2_3.x[3]) * .64 + .010085;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l5_6.gg[0] = 4. / (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_3.x[1];
    l5_6.gg[2] = 2.25 / (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_3.x[2];
    l5_6.gg[4] = 1. / (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_3.x[3];
    l5_6.gg[6] = .25 / (d__1 * d__1);
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l5_6.gg[1] = .16 / (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_3.x[1];
    l5_6.gg[3] = .36 / (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_3.x[2];
    l5_6.gg[5] = .64 / (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_3.x[3];
    l5_6.gg[7] = .64 / (d__1 * d__1);
L8:
    return 0;
} /* tp72_ */


/* Subroutine */ int tp73_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5;


    /* Local variables */
    static Real a[4];
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 1;
    l1_1.ninl = 1;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 1.;
	l13_3.xl[i__ - 1] = 0.;
	l11_3.lxl[i__ - 1] = true;
/* labelL6: */
	l12_3.lxu[i__ - 1] = false;
    }
    l4_3.gf[0] = 24.55;
    l4_3.gf[1] = 26.75;
    l4_3.gf[2] = 39.;
    l4_3.gf[3] = 40.5;
    l5_7.gg[0] = 2.3;
    l5_7.gg[3] = 5.6;
    l5_7.gg[6] = 11.1;
    l5_7.gg[9] = 1.3;
    for (i__ = 1; i__ <= 4; ++i__) {
/* L31: */
	l5_7.gg[i__ * 3 - 1] = 1.;
    }
    l20_6.lex = false;
    l20_6.xex[0] = .635521568605;
    l20_6.xex[1] = -1.1786227376e-12;
    l20_6.xex[2] = .312701880754;
    l20_6.xex[3] = .0517765506011;
    l20_6.fex = 29.8943781573;
    return 0;
labelL2:
    l6_1.fx = l2_3.x[0] * 24.55 + l2_3.x[1] * 26.75 + l2_3.x[2] * 39. + 
	    l2_3.x[3] * 40.5;
labelL3:
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_3.x[0] * 2.3 + l2_3.x[1] * 5.6 + l2_3.x[2] * 11.1 + 
		l2_3.x[3] * 1.3 - 5.;
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__2 = l2_3.x[0];
/* Computing 2nd power */
	d__3 = l2_3.x[1];
/* Computing 2nd power */
	d__4 = l2_3.x[2];
/* Computing 2nd power */
	d__5 = l2_3.x[3];
	d__1 = d__2 * d__2 * .28 + d__3 * d__3 * .19 + d__4 * d__4 * 20.5 + 
		d__5 * d__5 * .62;
	l3_3.g[1] = l2_3.x[0] * 12. + l2_3.x[1] * 11.9 + l2_3.x[2] * 41.8 + 
		l2_3.x[3] * 52.1 - pow_dd(&d__1, &c_b590) * 1.645 - 21.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = l2_3.x[0] + l2_3.x[1] + l2_3.x[2] + l2_3.x[3] - 1.;
    }
    return 0;
labelL5:
    if (! l10_4.index2[1]) {
	goto L8;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
/* L30: */
/* Computing 2nd power */
	d__2 = l2_3.x[0];
/* Computing 2nd power */
	d__3 = l2_3.x[1];
/* Computing 2nd power */
	d__4 = l2_3.x[2];
/* Computing 2nd power */
	d__5 = l2_3.x[3];
	d__1 = d__2 * d__2 * .28 + d__3 * d__3 * .19 + d__4 * d__4 * 20.5 + 
		d__5 * d__5 * .62;
	a[i__ - 1] = l2_3.x[i__ - 1] * 1.645 * pow_dd(&d__1, &c_b308);
    }
    l5_7.gg[1] = 12. - a[0] * .28;
    l5_7.gg[4] = 11.9 - a[1] * .19;
    l5_7.gg[7] = 41.8 - a[2] * 20.5;
    l5_7.gg[10] = 52.1 - a[3] * .62;
L8:
    return 0;
} /* tp73_ */


/* Subroutine */ int tp74_0_(int n__, int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real a[2];
    static int i__;
    static Real v1, v2;
    static int kn1;

    switch(n__) {
	case 1: goto L_tp75;
	}

    l20_6.xex[0] = 679.945319802;
    l20_6.xex[1] = 1026.06713256;
    l20_6.xex[2] = .11887636449;
    l20_6.xex[3] = -.39623355318;
    l20_6.fex = 5126.49810934;
    kn1 = 1;
    goto L7;

L_tp75:
    kn1 = 2;
    l20_6.xex[0] = 776.159220293;
    l20_6.xex[1] = 925.194939196;
    l20_6.xex[2] = .0511087936804;
    l20_6.xex[3] = -.428891137432;
    l20_6.fex = 5174.41288686;
L7:
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 2;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 3;
    a[0] = .55;
    a[1] = .48;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 0.;
	l11_3.lxl[i__ - 1] = true;
/* labelL6: */
	l12_3.lxu[i__ - 1] = true;
    }
    l13_3.xl[0] = 0.;
    l13_3.xl[1] = 0.;
    l13_3.xl[2] = -a[kn1 - 1];
    l13_3.xl[3] = -a[kn1 - 1];
    l14_3.xu[0] = 1200.;
    l14_3.xu[1] = 1200.;
    l14_3.xu[2] = a[kn1 - 1];
    l14_3.xu[3] = a[kn1 - 1];
    l4_3.gf[2] = 0.;
    l4_3.gf[3] = 0.;
    l5_13.gg[0] = 0.;
    l5_13.gg[5] = 0.;
    l5_13.gg[10] = -1.;
    l5_13.gg[15] = 1.;
    l5_13.gg[1] = 0.;
    l5_13.gg[6] = 0.;
    l5_13.gg[11] = 1.;
    l5_13.gg[16] = -1.;
    l5_13.gg[2] = -1.;
    l5_13.gg[7] = 0.;
    l5_13.gg[3] = 0.;
    l5_13.gg[8] = -1.;
    l5_13.gg[4] = 0.;
    l5_13.gg[9] = 0.;
    l20_6.lex = false;
    return 0;
labelL2:
/* Computing 3rd power */
    d__1 = l2_3.x[0];
/* Computing 3rd power */
    d__2 = l2_3.x[1];
    l6_1.fx = l2_3.x[0] * 3. + d__1 * (d__1 * d__1) * 1e-6 + l2_3.x[1] * 2. + 
	    d__2 * (d__2 * d__2) * 6.666666666666666e-7;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l4_3.gf[0] = d__1 * d__1 * 3e-6 + 3.;
/* Computing 2nd power */
    d__1 = l2_3.x[1];
    l4_3.gf[1] = d__1 * d__1 * 2e-6 + 2.;
    return 0;
labelL4:
    if (l9_5.index1[0]) {
	l3_4.g[0] = l2_3.x[3] - l2_3.x[2] + a[kn1 - 1];
    }
    if (l9_5.index1[1]) {
	l3_4.g[1] = l2_3.x[2] - l2_3.x[3] + a[kn1 - 1];
    }
    if (l9_5.index1[2]) {
	l3_4.g[2] = (sin(-l2_3.x[2] - .25) + std::sin(-l2_3.x[3] - .25)) * 1e3 + 
		894.8 - l2_3.x[0];
    }
    if (l9_5.index1[3]) {
	l3_4.g[3] = (sin(l2_3.x[2] - .25) + std::sin(l2_3.x[2] - l2_3.x[3] - .25)) 
		* 1e3 + 894.8 - l2_3.x[1];
    }
    if (l9_5.index1[4]) {
	l3_4.g[4] = (sin(l2_3.x[3] - .25) + std::sin(l2_3.x[3] - l2_3.x[2] - .25)) 
		* 1e3 + 1294.8;
    }
    return 0;
labelL5:
    if (! l10_5.index2[2]) {
	goto labelL9;
    }
    l5_13.gg[12] = std::cos(-l2_3.x[2] - .25) * -1e3;
    l5_13.gg[17] = std::cos(-l2_3.x[3] - .25) * -1e3;
labelL9:
    if (! l10_5.index2[3]) {
	goto labelL10;
    }
    v1 = std::cos(l2_3.x[2] - l2_3.x[3] - .25);
    l5_13.gg[13] = (std::cos(l2_3.x[2] - .25) + v1) * 1e3;
    l5_13.gg[18] = v1 * -1e3;
labelL10:
    if (! l10_5.index2[4]) {
	goto labelL11;
    }
    v2 = std::cos(l2_3.x[3] - l2_3.x[2] - .25);
    l5_13.gg[14] = v2 * -1e3;
    l5_13.gg[19] = (std::cos(l2_3.x[3] - .25) + v2) * 1e3;
labelL11:
    return 0;
} /* tp74_ */

/* Subroutine */ int tp74_(int *mode)
{
    return tp74_0_(0, mode);
    }

/* Subroutine */ int tp75_(int *mode)
{
    return tp74_0_(1, mode);
    }


/* Subroutine */ int tp76_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 3;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	l11_3.lxl[i__ - 1] = true;
	l12_3.lxu[i__ - 1] = false;
	l2_3.x[i__ - 1] = .5;
/* labelL6: */
	l13_3.xl[i__ - 1] = 0.;
    }
    l5_7.gg[0] = -1.;
    l5_7.gg[3] = -2.;
    l5_7.gg[6] = -1.;
    l5_7.gg[9] = -1.;
    l5_7.gg[1] = -3.;
    l5_7.gg[4] = -1.;
    l5_7.gg[7] = -2.;
    l5_7.gg[10] = 1.;
    l5_7.gg[2] = 0.;
    l5_7.gg[5] = 1.;
    l5_7.gg[8] = 4.;
    l5_7.gg[11] = 0.;
    l20_6.lex = false;
    l20_6.xex[0] = .272727272717;
    l20_6.xex[1] = 2.09090909094;
    l20_6.xex[2] = -2.63371889808e-11;
    l20_6.xex[3] = .545454545496;
    l20_6.fex = -4.68181818182;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_3.x[0];
/* Computing 2nd power */
    d__2 = l2_3.x[2];
/* Computing 2nd power */
    d__3 = l2_3.x[1];
/* Computing 2nd power */
    d__4 = l2_3.x[3];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + (d__3 * d__3 + d__4 * d__4) * .5 - 
	    l2_3.x[0] * l2_3.x[2] + l2_3.x[2] * l2_3.x[3] - l2_3.x[0] - 
	    l2_3.x[1] * 3. + l2_3.x[2] - l2_3.x[3];
    return 0;
labelL3:
    l4_3.gf[0] = l2_3.x[0] * 2. - l2_3.x[2] - 1.;
    l4_3.gf[1] = l2_3.x[1] - 3.;
    l4_3.gf[2] = l2_3.x[2] * 2. - l2_3.x[0] + l2_3.x[3] + 1.;
    l4_3.gf[3] = l2_3.x[3] + l2_3.x[2] - 1.;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = -l2_3.x[0] - l2_3.x[1] * 2. - l2_3.x[2] - l2_3.x[3] + 5.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_3.x[0] * -3. - l2_3.x[1] - l2_3.x[2] * 2. + l2_3.x[3] 
		+ 4.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = l2_3.x[1] + l2_3.x[2] * 4. - 1.5;
    }
labelL5:
    return 0;
} /* tp76_ */


/* Subroutine */ int tp77_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static int i__;
    static Real v1;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 2;
    for (i__ = 1; i__ <= 5; ++i__) {
	l2_4.x[i__ - 1] = 2.;
	l11_4.lxl[i__ - 1] = false;
/* labelL6: */
	l12_4.lxu[i__ - 1] = false;
    }
    l5_4.gg[2] = 0.;
    l5_4.gg[4] = 0.;
    l5_4.gg[1] = 0.;
    l5_4.gg[3] = 1.;
    l5_4.gg[9] = 0.;
    l20_7.lex = false;
    l20_7.xex[0] = 1.16617219726;
    l20_7.xex[1] = 1.18211138813;
    l20_7.xex[2] = 1.38025704044;
    l20_7.xex[3] = 1.50603627961;
    l20_7.xex[4] = .610920257517;
    l20_7.fex = .241505128786;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_4.x[0] - 1.;
/* Computing 2nd power */
    d__2 = l2_4.x[0] - l2_4.x[1];
/* Computing 2nd power */
    d__3 = l2_4.x[2] - 1.;
/* Computing 4th power */
    d__4 = l2_4.x[3] - 1., d__4 *= d__4;
/* Computing 6th power */
    d__5 = l2_4.x[4] - 1., d__5 *= d__5;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4 + d__5 * (
	    d__5 * d__5);
    return 0;
labelL3:
    l4_4.gf[0] = (l2_4.x[0] * 2. - l2_4.x[1] - 1.) * 2.;
    l4_4.gf[1] = (l2_4.x[0] - l2_4.x[1]) * -2.;
    l4_4.gf[2] = (l2_4.x[2] - 1.) * 2.;
/* Computing 3rd power */
    d__1 = l2_4.x[3] - 1.;
    l4_4.gf[3] = d__1 * (d__1 * d__1) * 4.;
/* Computing 5th power */
    d__1 = l2_4.x[4] - 1., d__2 = d__1, d__1 *= d__1;
    l4_4.gf[4] = d__2 * (d__1 * d__1) * 6.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_4.x[0];
	l3_2.g[0] = d__1 * d__1 * l2_4.x[3] + std::sin(l2_4.x[3] - l2_4.x[4]) - 
		std::sqrt(2.) * 2.;
    }
    if (l9_3.index1[1]) {
/* Computing 4th power */
	d__1 = l2_4.x[2], d__1 *= d__1;
/* Computing 2nd power */
	d__2 = l2_4.x[3];
	l3_2.g[1] = l2_4.x[1] + d__1 * d__1 * (d__2 * d__2) - 8. - std::sqrt(2.);
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    v1 = std::cos(l2_4.x[3] - l2_4.x[4]);
    l5_4.gg[0] = l2_4.x[0] * 2. * l2_4.x[3];
/* Computing 2nd power */
    d__1 = l2_4.x[0];
    l5_4.gg[6] = d__1 * d__1 + v1;
    l5_4.gg[8] = -v1;
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
/* Computing 3rd power */
    d__1 = l2_4.x[2];
/* Computing 2nd power */
    d__2 = l2_4.x[3];
    l5_4.gg[5] = d__1 * (d__1 * d__1) * 4. * (d__2 * d__2);
/* Computing 4th power */
    d__1 = l2_4.x[2], d__1 *= d__1;
    l5_4.gg[7] = d__1 * d__1 * 2. * l2_4.x[3];
L8:
    return 0;
} /* tp77_ */


/* Subroutine */ int tp78_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 3;
    for (i__ = 1; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = false;
/* labelL6: */
	l12_4.lxu[i__ - 1] = false;
    }
    l2_4.x[0] = -2.;
    l2_4.x[1] = 1.5;
    l2_4.x[2] = 2.;
    l2_4.x[3] = -1.;
    l2_4.x[4] = -1.;
    l5_9.gg[1] = 0.;
    l5_9.gg[8] = 0.;
    l5_9.gg[11] = 0.;
    l5_9.gg[14] = 0.;
    l20_7.lex = false;
    l20_7.xex[0] = -1.7171423423;
    l20_7.xex[1] = 1.59570826805;
    l20_7.xex[2] = 1.82724803488;
    l20_7.xex[3] = -.7636429466;
    l20_7.xex[4] = -.763643482853;
    l20_7.fex = -2.91970040911;
    return 0;
labelL2:
    l6_1.fx = l2_4.x[0] * l2_4.x[1] * l2_4.x[2] * l2_4.x[3] * l2_4.x[4];
    return 0;
labelL3:
    l4_4.gf[0] = l2_4.x[1] * l2_4.x[2] * l2_4.x[3] * l2_4.x[4];
    l4_4.gf[1] = l2_4.x[0] * l2_4.x[2] * l2_4.x[3] * l2_4.x[4];
    l4_4.gf[2] = l2_4.x[0] * l2_4.x[1] * l2_4.x[3] * l2_4.x[4];
    l4_4.gf[3] = l2_4.x[0] * l2_4.x[1] * l2_4.x[2] * l2_4.x[4];
    l4_4.gf[4] = l2_4.x[0] * l2_4.x[1] * l2_4.x[2] * l2_4.x[3];
    return 0;
labelL4:
    if (l9_4.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_4.x[0];
/* Computing 2nd power */
	d__2 = l2_4.x[1];
/* Computing 2nd power */
	d__3 = l2_4.x[2];
/* Computing 2nd power */
	d__4 = l2_4.x[3];
/* Computing 2nd power */
	d__5 = l2_4.x[4];
	l3_3.g[0] = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4 + 
		d__5 * d__5 - 10.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_4.x[1] * l2_4.x[2] - l2_4.x[3] * 5. * l2_4.x[4];
    }
    if (l9_4.index1[2]) {
/* Computing 3rd power */
	d__1 = l2_4.x[0];
/* Computing 3rd power */
	d__2 = l2_4.x[1];
	l3_3.g[2] = d__1 * (d__1 * d__1) + d__2 * (d__2 * d__2) + 1.;
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto L7;
    }
    for (i__ = 1; i__ <= 5; ++i__) {
/* L30: */
	l5_9.gg[i__ * 3 - 3] = l2_4.x[i__ - 1] * 2.;
    }
L7:
    if (! l10_4.index2[1]) {
	goto L8;
    }
    l5_9.gg[4] = l2_4.x[2];
    l5_9.gg[7] = l2_4.x[1];
    l5_9.gg[10] = l2_4.x[4] * -5.;
    l5_9.gg[13] = l2_4.x[3] * -5.;
L8:
    if (! l10_4.index2[2]) {
	goto labelL9;
    }
/* Computing 2nd power */
    d__1 = l2_4.x[0];
    l5_9.gg[2] = d__1 * d__1 * 3.;
/* Computing 2nd power */
    d__1 = l2_4.x[1];
    l5_9.gg[5] = d__1 * d__1 * 3.;
labelL9:
    return 0;
} /* tp78_ */


/* Subroutine */ int tp79_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static int i__;
    static Real v1, v2, v3, v4;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 3;
    for (i__ = 1; i__ <= 5; ++i__) {
	l2_4.x[i__ - 1] = 2.;
	l11_4.lxl[i__ - 1] = false;
/* labelL6: */
	l12_4.lxu[i__ - 1] = false;
    }
    l5_9.gg[0] = 1.;
    l5_9.gg[9] = 0.;
    l5_9.gg[12] = 0.;
    l5_9.gg[1] = 0.;
    l5_9.gg[4] = 1.;
    l5_9.gg[10] = 1.;
    l5_9.gg[13] = 0.;
    l5_9.gg[5] = 0.;
    l5_9.gg[8] = 0.;
    l5_9.gg[11] = 0.;
    l20_7.lex = false;
    l20_7.xex[0] = 1.19112745626;
    l20_7.xex[1] = 1.36260316492;
    l20_7.xex[2] = 1.4728179315;
    l20_7.xex[3] = 1.63501661894;
    l20_7.xex[4] = 1.67908143619;
    l20_7.fex = .0787768208538;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_4.x[0] - 1.;
/* Computing 2nd power */
    d__2 = l2_4.x[0] - l2_4.x[1];
/* Computing 2nd power */
    d__3 = l2_4.x[1] - l2_4.x[2];
/* Computing 4th power */
    d__4 = l2_4.x[2] - l2_4.x[3], d__4 *= d__4;
/* Computing 4th power */
    d__5 = l2_4.x[3] - l2_4.x[4], d__5 *= d__5;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4 + d__5 * 
	    d__5;
    return 0;
labelL3:
    v1 = l2_4.x[0] - l2_4.x[1];
    v2 = l2_4.x[1] - l2_4.x[2];
    v3 = l2_4.x[2] - l2_4.x[3];
    v4 = l2_4.x[3] - l2_4.x[4];
    l4_4.gf[0] = (l2_4.x[0] - 1. + v1) * 2.;
    l4_4.gf[1] = (v2 - v1) * 2.;
/* Computing 3rd power */
    d__1 = v3;
    l4_4.gf[2] = v2 * -2. + d__1 * (d__1 * d__1) * 4.;
/* Computing 3rd power */
    d__1 = v4;
/* Computing 3rd power */
    d__2 = v3;
    l4_4.gf[3] = (d__1 * (d__1 * d__1) - d__2 * (d__2 * d__2)) * 4.;
/* Computing 3rd power */
    d__1 = v4;
    l4_4.gf[4] = d__1 * (d__1 * d__1) * -4.;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_4.x[1];
/* Computing 3rd power */
	d__2 = l2_4.x[2];
	l3_3.g[0] = l2_4.x[0] + d__1 * d__1 + d__2 * (d__2 * d__2) - 2. - 
		std::sqrt(2.) * 3.;
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_4.x[2];
	l3_3.g[1] = l2_4.x[1] - d__1 * d__1 + l2_4.x[3] + 2. - std::sqrt(2.) * 2.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = l2_4.x[0] * l2_4.x[4] - 2.;
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto L7;
    }
    l5_9.gg[3] = l2_4.x[1] * 2.;
/* Computing 2nd power */
    d__1 = l2_4.x[2];
    l5_9.gg[6] = d__1 * d__1 * 3.;
L7:
    if (l10_4.index2[1]) {
	l5_9.gg[7] = l2_4.x[2] * -2.;
    }
    if (! l10_4.index2[2]) {
	goto labelL9;
    }
    l5_9.gg[2] = l2_4.x[4];
    l5_9.gg[14] = l2_4.x[0];
labelL9:
    return 0;
} /* tp79_ */


/* Subroutine */ int tp80_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static int i__;
    static Real t, v1, v2;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 3;
    for (i__ = 1; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = true;
/* labelL6: */
	l12_4.lxu[i__ - 1] = true;
    }
    l2_4.x[0] = -2.;
    l2_4.x[1] = 2.;
    l2_4.x[2] = 2.;
    l2_4.x[3] = -1.;
    l2_4.x[4] = -1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l13_4.xl[i__ - 1] = -2.3;
/* labelL20: */
	l14_4.xu[i__ - 1] = 2.3;
    }
    for (i__ = 3; i__ <= 5; ++i__) {
	l13_4.xl[i__ - 1] = -3.2;
/* L21: */
	l14_4.xu[i__ - 1] = 3.2;
    }
    l5_9.gg[1] = 0.;
    l5_9.gg[8] = 0.;
    l5_9.gg[11] = 0.;
    l5_9.gg[14] = 0.;
    l20_7.lex = false;
    l20_7.xex[0] = -1.71714294417;
    l20_7.xex[1] = 1.59570896503;
    l20_7.xex[2] = 1.82724691654;
    l20_7.xex[3] = -.763641279311;
    l20_7.xex[4] = -.763645016315;
    l20_7.fex = .0539498477624;
    return 0;
labelL2:
    l6_1.fx = std::exp(l2_4.x[0] * l2_4.x[1] * l2_4.x[2] * l2_4.x[3] * l2_4.x[4]);
    return 0;
labelL3:
    v1 = l2_4.x[3] * l2_4.x[4];
    v2 = l2_4.x[0] * l2_4.x[1];
    t = std::exp(v2 * l2_4.x[2] * v1);
    l4_4.gf[0] = l2_4.x[1] * l2_4.x[2] * v1 * t;
    l4_4.gf[1] = l2_4.x[0] * l2_4.x[2] * v1 * t;
    l4_4.gf[2] = v1 * v2 * t;
    l4_4.gf[3] = v2 * l2_4.x[2] * l2_4.x[4] * t;
    l4_4.gf[4] = v2 * l2_4.x[2] * l2_4.x[3] * t;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_4.x[0];
/* Computing 2nd power */
	d__2 = l2_4.x[1];
/* Computing 2nd power */
	d__3 = l2_4.x[2];
/* Computing 2nd power */
	d__4 = l2_4.x[3];
/* Computing 2nd power */
	d__5 = l2_4.x[4];
	l3_3.g[0] = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4 + 
		d__5 * d__5 - 10.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_4.x[1] * l2_4.x[2] - l2_4.x[3] * 5. * l2_4.x[4];
    }
    if (l9_4.index1[2]) {
/* Computing 3rd power */
	d__1 = l2_4.x[0];
/* Computing 3rd power */
	d__2 = l2_4.x[1];
	l3_3.g[2] = d__1 * (d__1 * d__1) + d__2 * (d__2 * d__2) + 1.;
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto L7;
    }
    for (i__ = 1; i__ <= 5; ++i__) {
/* L30: */
	l5_9.gg[i__ * 3 - 3] = l2_4.x[i__ - 1] * 2.;
    }
L7:
    if (! l10_4.index2[1]) {
	goto L8;
    }
    l5_9.gg[4] = l2_4.x[2];
    l5_9.gg[7] = l2_4.x[1];
    l5_9.gg[10] = l2_4.x[4] * -5.;
    l5_9.gg[13] = l2_4.x[3] * -5.;
L8:
    if (! l10_4.index2[2]) {
	goto labelL9;
    }
/* Computing 2nd power */
    d__1 = l2_4.x[0];
    l5_9.gg[2] = d__1 * d__1 * 3.;
/* Computing 2nd power */
    d__1 = l2_4.x[1];
    l5_9.gg[5] = d__1 * d__1 * 3.;
labelL9:
    return 0;
} /* tp80_ */


/* Subroutine */ int tp81_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static int i__;
    static Real t, v1, v2, v3;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 3;
    for (i__ = 1; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = true;
/* labelL6: */
	l12_4.lxu[i__ - 1] = true;
    }
    for (i__ = 1; i__ <= 2; ++i__) {
	l13_4.xl[i__ - 1] = -2.3;
/* labelL10: */
	l14_4.xu[i__ - 1] = 2.3;
    }
    for (i__ = 3; i__ <= 5; ++i__) {
	l13_4.xl[i__ - 1] = -3.2;
/* labelL20: */
	l14_4.xu[i__ - 1] = 3.2;
    }
    l2_4.x[0] = -2.;
    l2_4.x[1] = 2.;
    l2_4.x[2] = 2.;
    l2_4.x[3] = -1.;
    l2_4.x[4] = -1.;
    l5_9.gg[1] = 0.;
    l5_9.gg[8] = 0.;
    l5_9.gg[11] = 0.;
    l5_9.gg[14] = 0.;
    l20_7.lex = false;
    l20_7.xex[0] = -1.71714240091;
    l20_7.xex[1] = 1.59570833592;
    l20_7.xex[2] = 1.82724792592;
    l20_7.xex[3] = -.763647440817;
    l20_7.xex[4] = -.763638975604;
    l20_7.fex = .0539498477749;
    return 0;
labelL2:
/* Computing 3rd power */
    d__2 = l2_4.x[0];
/* Computing 3rd power */
    d__3 = l2_4.x[1];
/* Computing 2nd power */
    d__1 = d__2 * (d__2 * d__2) + d__3 * (d__3 * d__3) + 1.;
    l6_1.fx = std::exp(l2_4.x[0] * l2_4.x[1] * l2_4.x[2] * l2_4.x[3] * l2_4.x[4]) 
	    - d__1 * d__1 * .5;
    return 0;
labelL3:
/* Computing 3rd power */
    d__1 = l2_4.x[0];
/* Computing 3rd power */
    d__2 = l2_4.x[1];
    v1 = d__1 * (d__1 * d__1) + d__2 * (d__2 * d__2) + 1.;
    v2 = l2_4.x[0] * l2_4.x[1];
    v3 = l2_4.x[3] * l2_4.x[4];
    t = std::exp(v2 * v3 * l2_4.x[2]);
/* Computing 2nd power */
    d__1 = l2_4.x[0];
    l4_4.gf[0] = l2_4.x[1] * l2_4.x[2] * v3 * t - d__1 * d__1 * 3. * v1;
/* Computing 2nd power */
    d__1 = l2_4.x[1];
    l4_4.gf[1] = l2_4.x[0] * l2_4.x[2] * v3 * t - d__1 * d__1 * 3. * v1;
    l4_4.gf[2] = v2 * v3 * t;
    l4_4.gf[3] = v2 * l2_4.x[2] * l2_4.x[4] * t;
    l4_4.gf[4] = v2 * l2_4.x[2] * l2_4.x[3] * t;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_4.x[0];
/* Computing 2nd power */
	d__2 = l2_4.x[1];
/* Computing 2nd power */
	d__3 = l2_4.x[2];
/* Computing 2nd power */
	d__4 = l2_4.x[3];
/* Computing 2nd power */
	d__5 = l2_4.x[4];
	l3_3.g[0] = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4 + 
		d__5 * d__5 - 10.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_4.x[1] * l2_4.x[2] - l2_4.x[3] * 5. * l2_4.x[4];
    }
    if (l9_4.index1[2]) {
/* Computing 3rd power */
	d__1 = l2_4.x[0];
/* Computing 3rd power */
	d__2 = l2_4.x[1];
	l3_3.g[2] = d__1 * (d__1 * d__1) + d__2 * (d__2 * d__2) + 1.;
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto L7;
    }
    for (i__ = 1; i__ <= 5; ++i__) {
/* L30: */
	l5_9.gg[i__ * 3 - 3] = l2_4.x[i__ - 1] * 2.;
    }
L7:
    if (! l10_4.index2[1]) {
	goto L8;
    }
    l5_9.gg[4] = l2_4.x[2];
    l5_9.gg[7] = l2_4.x[1];
    l5_9.gg[10] = l2_4.x[4] * -5.;
    l5_9.gg[13] = l2_4.x[3] * -5.;
L8:
    if (! l10_4.index2[2]) {
	goto labelL9;
    }
/* Computing 2nd power */
    d__1 = l2_4.x[0];
    l5_9.gg[2] = d__1 * d__1 * 3.;
/* Computing 2nd power */
    d__1 = l2_4.x[1];
    l5_9.gg[5] = d__1 * d__1 * 3.;
labelL9:
    return 0;
} /* tp81_ */


/* Subroutine */ int tp83_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real a, b, c__, d__;
    static int i__;
    static Real a1, a2, a3, a4, a5, a6, a7, a8, a9, v1, v2, v3, a10, 
	    a11, a12;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 6;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_4.x[0] = 78.;
    l2_4.x[1] = 33.;
    l2_4.x[2] = 27.;
    l2_4.x[3] = 27.;
    l2_4.x[4] = 27.;
    for (i__ = 1; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = true;
/* labelL6: */
	l12_4.lxu[i__ - 1] = true;
    }
    l13_4.xl[0] = 78.;
    l13_4.xl[1] = 33.;
    l14_4.xu[0] = 102.;
    l14_4.xu[1] = 45.;
    for (i__ = 3; i__ <= 5; ++i__) {
	l13_4.xl[i__ - 1] = 27.;
/* L31: */
	l14_4.xu[i__ - 1] = 45.;
    }
    a = 5.3578547;
    b = .8356891;
    c__ = 37.293239;
    d__ = 40792.141;
    a1 = 85.334407;
    a2 = .0056858;
    a3 = 6.262e-4;
    a4 = .0022053;
    a5 = 80.51249;
    a6 = .0071317;
    a7 = .0029955;
    a8 = .0021813;
    a9 = 9.300961;
    a10 = .0047026;
    a11 = .0012547;
    a12 = .0019085;
    l4_4.gf[1] = 0.;
    l4_4.gf[3] = 0.;
    l5_14.gg[19] = 0.;
    l5_14.gg[8] = 0.;
    l5_14.gg[22] = 0.;
    l5_14.gg[11] = 0.;
    l20_7.lex = false;
    l20_7.xex[0] = 78.;
    l20_7.xex[1] = 33.;
    l20_7.xex[2] = 29.9952560253;
    l20_7.xex[3] = 45.;
    l20_7.xex[4] = 36.7758129081;
    l20_7.fex = -30665.5386717;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_4.x[2];
    l6_1.fx = a * (d__1 * d__1) + b * l2_4.x[0] * l2_4.x[4] + c__ * l2_4.x[0] 
	    - d__;
    return 0;
labelL3:
    l4_4.gf[0] = b * l2_4.x[4] + c__;
    l4_4.gf[2] = a * 2. * l2_4.x[2];
    l4_4.gf[4] = b * l2_4.x[0];
    return 0;
labelL4:
    if (! (l9_6.index1[0] || l9_6.index1[3])) {
	goto L41;
    }
    v1 = a1 + a2 * l2_4.x[1] * l2_4.x[4] + a3 * l2_4.x[0] * l2_4.x[3] - a4 * 
	    l2_4.x[2] * l2_4.x[4];
    if (l9_6.index1[0]) {
	l3_5.g[0] = v1;
    }
    if (l9_6.index1[3]) {
	l3_5.g[3] = 92. - v1;
    }
L41:
    if (! (l9_6.index1[1] || l9_6.index1[4])) {
	goto L42;
    }
/* Computing 2nd power */
    d__1 = l2_4.x[2];
    v2 = a5 + a6 * l2_4.x[1] * l2_4.x[4] + a7 * l2_4.x[0] * l2_4.x[1] + a8 * (
	    d__1 * d__1) - 90.;
    if (l9_6.index1[1]) {
	l3_5.g[1] = v2;
    }
    if (l9_6.index1[4]) {
	l3_5.g[4] = 20. - v2;
    }
L42:
    if (! (l9_6.index1[2] || l9_6.index1[5])) {
	goto L43;
    }
    v3 = a9 + a10 * l2_4.x[2] * l2_4.x[4] + a11 * l2_4.x[0] * l2_4.x[2] + a12 
	    * l2_4.x[2] * l2_4.x[3] - 20.;
    if (l9_6.index1[2]) {
	l3_5.g[2] = v3;
    }
    if (l9_6.index1[5]) {
	l3_5.g[5] = 5. - v3;
    }
L43:
    return 0;
labelL5:
    if (! l10_6.index2[0]) {
	goto L7;
    }
    l5_14.gg[0] = a3 * l2_4.x[3];
    l5_14.gg[6] = a2 * l2_4.x[4];
    l5_14.gg[12] = -a4 * l2_4.x[4];
    l5_14.gg[18] = a3 * l2_4.x[0];
    l5_14.gg[24] = a2 * l2_4.x[1] - a4 * l2_4.x[2];
L7:
    if (! l10_6.index2[1]) {
	goto L8;
    }
    l5_14.gg[1] = a7 * l2_4.x[1];
    l5_14.gg[7] = a6 * l2_4.x[4] + a7 * l2_4.x[0];
    l5_14.gg[13] = a8 * 2. * l2_4.x[2];
    l5_14.gg[25] = a6 * l2_4.x[1];
L8:
    if (! l10_6.index2[2]) {
	goto labelL9;
    }
    l5_14.gg[2] = a11 * l2_4.x[2];
    l5_14.gg[14] = a10 * l2_4.x[4] + a11 * l2_4.x[0] + a12 * l2_4.x[3];
    l5_14.gg[20] = a12 * l2_4.x[2];
    l5_14.gg[26] = a10 * l2_4.x[2];
labelL9:
    if (! l10_6.index2[3]) {
	goto labelL10;
    }
    l5_14.gg[3] = -a3 * l2_4.x[3];
    l5_14.gg[9] = -a2 * l2_4.x[4];
    l5_14.gg[15] = a4 * l2_4.x[4];
    l5_14.gg[21] = -a3 * l2_4.x[0];
    l5_14.gg[27] = -a2 * l2_4.x[1] + a4 * l2_4.x[2];
labelL10:
    if (! l10_6.index2[4]) {
	goto labelL11;
    }
    l5_14.gg[4] = -a7 * l2_4.x[1];
    l5_14.gg[10] = -a6 * l2_4.x[4] - a7 * l2_4.x[0];
    l5_14.gg[16] = a8 * -2. * l2_4.x[2];
    l5_14.gg[28] = -a6 * l2_4.x[1];
labelL11:
    if (! l10_6.index2[5]) {
	goto labelL12;
    }
    l5_14.gg[5] = -a11 * l2_4.x[2];
    l5_14.gg[17] = -a10 * l2_4.x[4] - a11 * l2_4.x[0] - a12 * l2_4.x[3];
    l5_14.gg[23] = -a12 * l2_4.x[2];
    l5_14.gg[29] = -a10 * l2_4.x[2];
labelL12:
    return 0;
} /* tp83_ */


/* Subroutine */ int tp84_(int *mode)
{
    static Real a[21], b[3];
    static int i__, j, i1;
    static Real v1;

    a[0] = -24345.;
    a[1] = -8720288.849;
    a[2] = 150512.5253;
    a[3] = -156.6950325;
    a[4] = 476470.3222;
    a[5] = 729482.8271;
    a[6] = -145421.402;
    a[7] = 2931.1506;
    a[8] = -40.427932;
    a[9] = 5106.192;
    a[10] = 15711.36;
    a[11] = -155011.1084;
    a[12] = 4360.53352;
    a[13] = 12.9492344;
    a[14] = 10236.884;
    a[15] = 13176.786;
    a[16] = -326669.5104;
    a[17] = 7390.68412;
    a[18] = -27.8986976;
    a[19] = 16643.076;
    a[20] = 30988.146;
    b[0] = 2.94e5;
    b[1] = 2.94e5;
    b[2] = 277200.;
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 6;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_4.x[0] = 2.52;
    l2_4.x[1] = 2.;
    l2_4.x[2] = 37.5;
    l2_4.x[3] = 9.25;
    l2_4.x[4] = 6.8;
    for (i__ = 1; i__ <= 5; ++i__) {
	l12_4.lxu[i__ - 1] = true;
/* labelL6: */
	l11_4.lxl[i__ - 1] = true;
    }
    l13_4.xl[0] = 0.;
    l13_4.xl[1] = 1.2;
    l13_4.xl[2] = 20.;
    l13_4.xl[3] = 9.;
    l13_4.xl[4] = 6.5;
    l14_4.xu[0] = 1e3;
    l14_4.xu[1] = 2.4;
    l14_4.xu[2] = 60.;
    l14_4.xu[3] = 9.3;
    l14_4.xu[4] = 7.;
    l20_7.lex = false;
    l20_7.xex[0] = 4.53743097466;
    l20_7.xex[1] = 2.40000000002;
    l20_7.xex[2] = 60.;
    l20_7.xex[3] = 9.29999999999;
    l20_7.xex[4] = 7.;
    l20_7.fex = -52.803351330600002;
    return 0;
labelL2:
    l6_1.fx = -(a[0] + l2_4.x[0] * (a[1] + a[2] * l2_4.x[1] + a[3] * l2_4.x[2]
	     + a[4] * l2_4.x[3] + a[5] * l2_4.x[4]));
    l6_1.fx *= 1e-5;
    return 0;
labelL3:
    l4_4.gf[0] = -(a[1] + a[2] * l2_4.x[1] + a[3] * l2_4.x[2] + a[4] * l2_4.x[
	    3] + a[5] * l2_4.x[4]) * 1e-5;
    for (i__ = 2; i__ <= 5; ++i__) {
/* L30: */
	l4_4.gf[i__ - 1] = -a[i__] * l2_4.x[0] * 1e-5;
    }
    return 0;
labelL4:
    for (i__ = 1; i__ <= 3; ++i__) {
	if (! (l9_6.index1[i__ - 1] || l9_6.index1[i__ + 2])) {
	    goto L80;
	}
	i1 = i__ * 5;
	v1 = l2_4.x[0] * (a[i1 + 1] + a[i1 + 2] * l2_4.x[1] + a[i1 + 3] * 
		l2_4.x[2] + a[i1 + 4] * l2_4.x[3] + a[i1 + 5] * l2_4.x[4]);
	if (l9_6.index1[i__ - 1]) {
	    l3_5.g[i__ - 1] = v1;
	}
	if (l9_6.index1[i__ + 2]) {
	    l3_5.g[i__ + 2] = b[i__ - 1] - v1;
	}
L80:
	;
    }
    return 0;
labelL5:
    for (i__ = 1; i__ <= 3; ++i__) {
	if (! (l10_6.index2[i__ - 1] || l10_6.index2[i__ + 2])) {
	    goto L90;
	}
	i1 = i__ * 5 + 1;
	if (! l10_6.index2[i__ - 1]) {
	    goto L95;
	}
	l5_14.gg[i__ - 1] = a[i1] + a[i1 + 1] * l2_4.x[1] + a[i1 + 2] * 
		l2_4.x[2] + a[i1 + 3] * l2_4.x[3] + a[i1 + 4] * l2_4.x[4];
	for (j = 2; j <= 5; ++j) {
/* L91: */
	    l5_14.gg[i__ + j * 6 - 7] = a[i1 + j - 1] * l2_4.x[0];
	}
	if (! l10_6.index2[i__ + 2]) {
	    goto L90;
	}
	for (j = 1; j <= 5; ++j) {
/* L92: */
	    l5_14.gg[i__ + 3 + j * 6 - 7] = -l5_14.gg[i__ + j * 6 - 7];
	}
	goto L90;
L95:
	l5_14.gg[i__ + 2] = -(a[i1] + a[i1 + 1] * l2_4.x[1] + a[i1 + 2] * 
		l2_4.x[2] + a[i1 + 3] * l2_4.x[3] + a[i1 + 4] * l2_4.x[4]);
	for (j = 2; j <= 5; ++j) {
/* L96: */
	    l5_14.gg[i__ + 3 + j * 6 - 7] = -a[i1 + j - 1] * l2_4.x[0];
	}
L90:
	;
    }
    return 0;
} /* tp84_ */


/* Subroutine */ int tp85_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real a[17], b[17], c__[17];
    static int i__, j;
    static Real y[17], v1, v2, v3, v4, v5, v6, v7, dc[85]	/* was [17][5]
	     */, dy[85]	/* was [17][5] */;

    if (*mode - 2 >= 0) {
	goto L17;
    } else {
	goto labelL1;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 3;
    l1_1.ninl = 35;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_4.x[0] = 900.;
    l2_4.x[1] = 80.;
    l2_4.x[2] = 115.;
    l2_4.x[3] = 267.;
    l2_4.x[4] = 27.;
    for (i__ = 1; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = true;
/* labelL6: */
	l12_4.lxu[i__ - 1] = true;
    }
    l13_4.xl[0] = 704.4148;
    l13_4.xl[1] = 68.6;
    l13_4.xl[2] = 0.;
    l13_4.xl[3] = 193.;
    l13_4.xl[4] = 25.;
    l14_4.xu[0] = 906.3855;
    l14_4.xu[1] = 288.88;
    l14_4.xu[2] = 134.75;
    l14_4.xu[3] = 287.0966;
    l14_4.xu[4] = 84.1988;
    a[1] = 17.505;
    a[2] = 11.275;
    a[3] = 214.228;
    a[4] = 7.458;
    a[5] = .961;
    a[6] = 1.612;
    a[7] = .146;
    a[8] = 107.99;
    a[9] = 922.693;
    a[10] = 926.832;
    a[11] = 18.766;
    a[12] = 1072.163;
    a[13] = 8961.448;
    a[14] = .063;
    a[15] = 71084.33;
    a[16] = 2802713.;
    b[1] = 1053.6667;
    b[2] = 35.03;
    b[3] = 665.585;
    b[4] = 584.463;
    b[5] = 265.916;
    b[6] = 7.046;
    b[7] = .222;
    b[8] = 273.366;
    b[9] = 1286.105;
    b[10] = 1444.046;
    b[11] = 537.141;
    b[12] = 3247.039;
    b[13] = 26844.086;
    b[14] = .386;
    b[15] = 1.4e5;
    b[16] = 12146108.;
    for (i__ = 1; i__ <= 5; ++i__) {
	l5_15.gg[i__ * 38 - 38] = 0.;
	l5_15.gg[i__ * 38 - 37] = 0.;
	l5_15.gg[i__ * 38 - 36] = 0.;
	dc[i__ * 17 - 17] = 0.;
	dc[i__ * 17 - 13] = 0.;
	dc[i__ * 17 - 8] = 0.;
/* L61: */
    }
    l5_15.gg[38] = 1.5;
    l5_15.gg[76] = -1.;
    l5_15.gg[39] = 1.;
    l5_15.gg[77] = 1.;
    l5_15.gg[40] = -1.;
    l5_15.gg[78] = -1.;
    dy[0] = 0.;
    dy[17] = 1.;
    dy[34] = 1.;
    dy[51] = 0.;
    dy[68] = 0.;
    dc[51] = .024;
    dc[21] = 100.;
    l20_7.lex = false;
    l20_7.xex[0] = 705.180328772;
    l20_7.xex[1] = 68.6000529425;
    l20_7.xex[2] = 102.900013236;
    l20_7.xex[3] = 282.324998587;
    l20_7.xex[4] = 37.5850413432;
    l20_7.fex = -1.90513375046;
    return 0;
L17:
    y[0] = l2_4.x[1] + l2_4.x[2] + 41.6;
    c__[0] = l2_4.x[3] * .024 - 4.62;
    y[1] = 12.5 / c__[0] + 12.;
    v3 = y[1] * l2_4.x[0];
    c__[1] = (l2_4.x[0] * 3.535e-4 + .5311) * l2_4.x[0] + v3 * .08705;
    c__[2] = l2_4.x[0] * .052 + 78. + v3 * .002377;
    y[2] = c__[1] / c__[2];
    y[3] = y[2] * 19.;
    v1 = l2_4.x[0] - y[2];
    c__[3] = (v1 * .1956 / l2_4.x[1] + .04782) * v1 + y[3] * .6376 + y[2] * 
	    1.594;
    c__[4] = l2_4.x[1] * 100.;
    c__[5] = v1 - y[3];
    c__[6] = .95 - c__[3] / c__[4];
    y[4] = c__[5] * c__[6];
    v2 = y[4] + y[3];
    y[5] = v1 - v2;
    c__[7] = v2 * .995;
    y[6] = c__[7] / y[0];
    y[7] = c__[7] / 3798.;
    c__[8] = y[6] - y[6] * .0663 / y[7] - .3153;
    y[8] = 96.82 / c__[8] + y[0] * .321;
    y[9] = y[4] * 1.29 + y[3] * 1.258 + y[2] * 2.29 + y[5] * 1.71;
    y[10] = l2_4.x[0] * 1.71 - y[3] * .452 + y[2] * .58;
    c__[9] = .016349860428020738;
    c__[10] = v3 * 1.74125;
    c__[11] = y[9] * .995 + 1998.;
    y[11] = c__[9] * l2_4.x[0] + c__[10] / c__[11];
    y[12] = c__[11] - y[1] * 1.75;
    v4 = y[8] + l2_4.x[4];
    y[13] = l2_4.x[1] * 64.4 + 3623. + l2_4.x[2] * 58.4 + 146312. / v4;
    c__[12] = y[9] * .995 + l2_4.x[1] * 60.8 + l2_4.x[3] * 48. - y[13] * 
	    .1121 - 5095.;
    y[14] = y[12] / c__[12];
    y[15] = 1.48e5 - y[14] * 3.31e5 + y[12] * 40. - y[14] * 61. * y[12];
    c__[13] = y[9] * 2324. - y[1] * 2.874e7;
    y[16] = 1.413e7 - y[9] * 1328. - y[10] * 531. + c__[13] / c__[11];
    c__[14] = y[12] / y[14] - y[12] / .52;
    c__[15] = 1.104 - y[14] * .72;
    c__[16] = v4;
    if (*mode == 3 || *mode == 5) {
	goto L71;
    }
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
L71:
    for (i__ = 1; i__ <= 5; ++i__) {
/* L30: */
/* Computing 2nd power */
	d__1 = c__[0];
	dy[i__ * 17 - 16] = dc[i__ * 17 - 17] * -12.5 / (d__1 * d__1);
    }
    v5 = y[1] + l2_4.x[0] * dy[1];
    dc[1] = l2_4.x[0] * 7.07e-4 + .5311 + v5 * .08705;
    dc[2] = v5 * .002377 + .052;
    for (i__ = 2; i__ <= 5; ++i__) {
	v6 = l2_4.x[0] * dy[i__ * 17 - 16];
	dc[i__ * 17 - 16] = v6 * .08705;
/* L32: */
	dc[i__ * 17 - 15] = v6 * .002377;
    }
    for (i__ = 1; i__ <= 5; ++i__) {
/* Computing 2nd power */
	d__1 = c__[2];
	dy[i__ * 17 - 15] = (c__[2] * dc[i__ * 17 - 16] - c__[1] * dc[i__ * 
		17 - 15]) / (d__1 * d__1);
/* L33: */
	dy[i__ * 17 - 14] = dy[i__ * 17 - 15] * 19.;
    }
    dc[3] = (v1 * .3912 / l2_4.x[1] + .04782) * (1. - dy[2]) + dy[3] * .6376 
	    + dy[2] * 1.594;
/* Computing 2nd power */
    d__1 = l2_4.x[1];
    dc[20] = v1 * -.1956 * (v1 + l2_4.x[1] * 2. * dy[19]) / (d__1 * d__1) + 
	    dy[19] * 1.54618 + dy[20] * .6376;
    for (i__ = 3; i__ <= 5; ++i__) {
/* L34: */
	dc[i__ * 17 - 14] = (1.54618 - v1 * .3912 / l2_4.x[1]) * dy[i__ * 17 
		- 15] + dy[i__ * 17 - 14] * .6376;
    }
    dc[5] = 1. - dy[2] - dy[3];
    for (i__ = 2; i__ <= 5; ++i__) {
/* L35: */
	dc[i__ * 17 - 12] = -dy[i__ * 17 - 15] - dy[i__ * 17 - 14];
    }
    for (i__ = 1; i__ <= 5; ++i__) {
/* Computing 2nd power */
	d__1 = c__[4];
	dc[i__ * 17 - 11] = -(c__[4] * dc[i__ * 17 - 14] - c__[3] * dc[i__ * 
		17 - 13]) / (d__1 * d__1);
/* L36: */
	dy[i__ * 17 - 13] = c__[5] * dc[i__ * 17 - 11] + c__[6] * dc[i__ * 17 
		- 12];
    }
    for (i__ = 1; i__ <= 5; ++i__) {
/* L37: */
	dy[i__ * 17 - 12] = -dy[i__ * 17 - 13] - dy[i__ * 17 - 14] - dy[i__ * 
		17 - 15];
    }
    dy[5] += 1.;
    for (i__ = 1; i__ <= 5; ++i__) {
	dc[i__ * 17 - 10] = (dy[i__ * 17 - 13] + dy[i__ * 17 - 14]) * .995;
/* Computing 2nd power */
	d__1 = y[0];
	dy[i__ * 17 - 11] = (y[0] * dc[i__ * 17 - 10] - c__[7] * dy[i__ * 17 
		- 17]) / (d__1 * d__1);
	dy[i__ * 17 - 10] = dc[i__ * 17 - 10] / 3798.;
/* Computing 2nd power */
	d__1 = y[7];
	dc[i__ * 17 - 9] = dy[i__ * 17 - 11] - (y[7] * dy[i__ * 17 - 11] - y[
		6] * dy[i__ * 17 - 10]) * .0663 / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = c__[8];
	dy[i__ * 17 - 9] = dc[i__ * 17 - 9] * -96.82 / (d__1 * d__1) + dy[i__ 
		* 17 - 17] * .321;
/* L38: */
	dy[i__ * 17 - 8] = dy[i__ * 17 - 13] * 1.29 + dy[i__ * 17 - 14] * 
		1.258 + dy[i__ * 17 - 15] * 2.29 + dy[i__ * 17 - 12] * 1.71;
    }
    dy[10] = 1.71 - dy[3] * .452 + dy[2] * .58;
    for (i__ = 2; i__ <= 5; ++i__) {
/* L39: */
	dy[i__ * 17 - 7] = dy[i__ * 17 - 14] * -.452 + dy[i__ * 17 - 15] * 
		.58;
    }
    dc[10] = (y[1] + l2_4.x[0] * dy[1]) * 1.74125;
    for (i__ = 2; i__ <= 5; ++i__) {
/* L40: */
	dc[i__ * 17 - 7] = l2_4.x[0] * 1.74125 * dy[i__ * 17 - 16];
    }
    for (i__ = 1; i__ <= 5; ++i__) {
/* L41: */
	dc[i__ * 17 - 6] = dy[i__ * 17 - 8] * .995;
    }
/* Computing 2nd power */
    d__1 = c__[11];
    dy[11] = c__[9] + l2_4.x[0] * dc[9] + (c__[11] * dc[10] - c__[10] * dc[11]
	    ) / (d__1 * d__1);
    for (i__ = 2; i__ <= 5; ++i__) {
/* L42: */
/* Computing 2nd power */
	d__1 = c__[11];
	dy[i__ * 17 - 6] = (c__[11] * dc[i__ * 17 - 7] - c__[10] * dc[i__ * 
		17 - 6]) / (d__1 * d__1);
    }
    for (i__ = 1; i__ <= 5; ++i__) {
/* L43: */
	dy[i__ * 17 - 5] = dc[i__ * 17 - 6] - dy[i__ * 17 - 16] * 1.75;
    }
/* Computing 2nd power */
    d__1 = v4;
    v7 = -146312. / (d__1 * d__1);
    dy[13] = v7 * dy[8];
    dy[30] = v7 * dy[25] + 64.4;
    dy[47] = v7 * dy[42] + 58.4;
    dy[64] = v7 * dy[59];
    dy[81] = v7 * (dy[76] + 1.);
    for (i__ = 1; i__ <= 5; ++i__) {
/* L44: */
	dc[i__ * 17 - 5] = dy[i__ * 17 - 8] * .995 - dy[i__ * 17 - 4] * .1121;
    }
    dc[29] += 60.8;
    dc[63] += 48.;
    for (i__ = 1; i__ <= 5; ++i__) {
/* Computing 2nd power */
	d__1 = c__[12];
	dy[i__ * 17 - 3] = (c__[12] * dy[i__ * 17 - 5] - y[12] * dc[i__ * 17 
		- 5]) / (d__1 * d__1);
	dy[i__ * 17 - 2] = dy[i__ * 17 - 3] * -3.31e5 + dy[i__ * 17 - 5] * 
		40. - (y[14] * dy[i__ * 17 - 5] + y[12] * dy[i__ * 17 - 3]) * 
		61.;
	dc[i__ * 17 - 4] = dy[i__ * 17 - 8] * 2324. - dy[i__ * 17 - 16] * 
		2.874e7;
/* Computing 2nd power */
	d__1 = c__[11];
	dy[i__ * 17 - 1] = dy[i__ * 17 - 8] * -1328. - dy[i__ * 17 - 7] * 
		531. + (c__[11] * dc[i__ * 17 - 4] - c__[13] * dc[i__ * 17 - 
		6]) / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = y[14];
	dc[i__ * 17 - 3] = (y[14] * dy[i__ * 17 - 5] - y[12] * dy[i__ * 17 - 
		3]) / (d__1 * d__1) - dy[i__ * 17 - 5] / .52;
/* L45: */
	dc[i__ * 17 - 2] = dy[i__ * 17 - 3] * -.72;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
/* L46: */
	dc[i__ * 17 - 1] = dy[i__ * 17 - 9];
    }
    dc[84] = dy[76] + 1.;
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL2:
    l6_1.fx = -(y[16] * 5.843e-7 - y[13] * 1.17e-4 - .1365 - y[12] * 2.358e-5 
	    - y[15] * 1.502e-6 - y[11] * .0321 - y[4] * .004324 - c__[14] * 
	    1e-4 / c__[15] - y[1] * 37.48 / c__[11]);
    return 0;
labelL3:
    for (i__ = 1; i__ <= 5; ++i__) {
/* Computing 2nd power */
	d__1 = c__[15];
/* Computing 2nd power */
	d__2 = c__[11];
	l4_4.gf[i__ - 1] = dy[i__ * 17 - 1] * -5.843e-7 + dy[i__ * 17 - 4] * 
		1.17e-4 + dy[i__ * 17 - 5] * 2.358e-5 + dy[i__ * 17 - 2] * 
		1.502e-6 + dy[i__ * 17 - 6] * .0321 + dy[i__ * 17 - 13] * 
		.004324 + (c__[15] * dc[i__ * 17 - 3] - c__[14] * dc[i__ * 17 
		- 2]) * 1e-4 / (d__1 * d__1) + (c__[11] * dy[i__ * 17 - 16] - 
		y[1] * dc[i__ * 17 - 6]) * 37.48 / (d__2 * d__2);
/* L47: */
    }
    return 0;
labelL4:
    if (l9_9.index1[0]) {
	l3_8.g[0] = l2_4.x[1] * 1.5 - l2_4.x[2];
    }
    if (l9_9.index1[1]) {
	l3_8.g[1] = y[0] - 213.1;
    }
    if (l9_9.index1[2]) {
	l3_8.g[2] = 405.23 - y[0];
    }
    for (i__ = 1; i__ <= 16; ++i__) {
	if (l9_9.index1[i__ + 2]) {
	    l3_8.g[i__ + 2] = y[i__] - a[i__];
	}
	if (l9_9.index1[i__ + 18]) {
	    l3_8.g[i__ + 18] = b[i__] - y[i__];
	}
/* L50: */
    }
    if (l9_9.index1[35]) {
	l3_8.g[35] = y[3] - y[4] * .28 / .72;
    }
    if (l9_9.index1[36]) {
	l3_8.g[36] = 21. - y[1] * 3496. / c__[11];
    }
    if (l9_9.index1[37]) {
	l3_8.g[37] = 62212. / c__[16] - 110.6 - y[0];
    }
    return 0;
labelL5:
    for (i__ = 1; i__ <= 16; ++i__) {
	if (! l10_9.index2[i__ + 2]) {
	    goto L52;
	}
	for (j = 1; j <= 5; ++j) {
/* L51: */
	    l5_15.gg[i__ + 3 + j * 38 - 39] = dy[i__ + 1 + j * 17 - 18];
	}
L52:
	if (! l10_9.index2[i__ + 18]) {
	    goto L54;
	}
	for (j = 1; j <= 5; ++j) {
/* L53: */
	    l5_15.gg[i__ + 19 + j * 38 - 39] = -dy[i__ + 1 + j * 17 - 18];
	}
L54:
	;
    }
    if (! l10_9.index2[35]) {
	goto L56;
    }
    for (j = 1; j <= 5; ++j) {
/* L55: */
	l5_15.gg[j * 38 - 3] = dy[j * 17 - 14] - dy[j * 17 - 13] * .28 / .72;
    }
L56:
    if (! l10_9.index2[36]) {
	goto L58;
    }
    for (j = 1; j <= 5; ++j) {
/* L57: */
/* Computing 2nd power */
	d__1 = c__[11];
	l5_15.gg[j * 38 - 2] = (c__[11] * dy[j * 17 - 16] - y[1] * dc[j * 17 
		- 6]) * -3496. / (d__1 * d__1);
    }
L58:
    if (! l10_9.index2[37]) {
	goto L60;
    }
    for (j = 1; j <= 5; ++j) {
/* L59: */
/* Computing 2nd power */
	d__1 = c__[16];
	l5_15.gg[j * 38 - 1] = dc[j * 17 - 1] * -62212. / (d__1 * d__1) - dy[
		j * 17 - 17];
    }
L60:
    return 0;
} /* tp85_ */


/* Subroutine */ int tp86_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real a[50]	/* was [10][5] */, b[10], c__[25]	/* 
	    was [5][5] */, d__[5], e[5];
    static int i__, j;
    static Real t;
    static int i1;
    static Real t1;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 10;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
/* L26: */
	l2_4.x[i__ - 1] = 0.;
    }
    l2_4.x[4] = 1.;
    for (i__ = 1; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = true;
	l12_4.lxu[i__ - 1] = false;
/* labelL6: */
	l13_4.xl[i__ - 1] = 0.;
    }
    e[0] = -15.;
    e[1] = -27.;
    e[2] = -36.;
    e[3] = -18.;
    e[4] = -12.;
    c__[0] = 30.;
    c__[5] = -20.;
    c__[10] = -10.;
    c__[15] = 32.;
    c__[20] = -10.;
    c__[6] = 39.;
    c__[11] = -6.;
    c__[16] = -31.;
    c__[21] = 32.;
    c__[12] = 10.;
    c__[17] = -6.;
    c__[22] = -10.;
    c__[18] = 39.;
    c__[23] = -20.;
    c__[24] = 30.;
    for (i__ = 1; i__ <= 4; ++i__) {
	i1 = i__ + 1;
	for (j = i1; j <= 5; ++j) {
/* L27: */
	    c__[j + i__ * 5 - 6] = c__[i__ + j * 5 - 6];
	}
    }
    d__[0] = 4.;
    d__[1] = 8.;
    d__[2] = 10.;
    d__[3] = 6.;
    d__[4] = 2.;
    a[0] = -16.;
    a[10] = 2.;
    a[20] = 0.;
    a[30] = 1.;
    a[40] = 0.;
    a[1] = 0.;
    a[11] = -2.;
    a[21] = 0.;
    a[31] = .4;
    a[41] = 2.;
    a[2] = -3.5;
    a[12] = 0.;
    a[22] = 2.;
    a[32] = 0.;
    a[42] = 0.;
    a[3] = 0.;
    a[13] = -2.;
    a[23] = 0.;
    a[33] = -4.;
    a[43] = -1.;
    a[4] = 0.;
    a[14] = -9.;
    a[24] = -2.;
    a[34] = 1.;
    a[44] = -2.8;
    a[5] = 2.;
    a[15] = 0.;
    a[25] = -4.;
    a[35] = 0.;
    a[45] = 0.;
    a[7] = -1.;
    a[17] = -2.;
    a[27] = -3.;
    a[37] = -2.;
    a[47] = -1.;
    for (i__ = 1; i__ <= 5; ++i__) {
	a[i__ * 10 - 4] = -1.;
	a[i__ * 10 - 2] = (Real) i__;
/* L29: */
	a[i__ * 10 - 1] = 1.;
    }
    b[0] = -40.;
    b[1] = -2.;
    b[2] = -.25;
    b[3] = -4.;
    b[4] = -4.;
    b[5] = -1.;
    b[6] = -40.;
    b[7] = -60.;
    b[8] = 5.;
    b[9] = 1.;
    for (i__ = 1; i__ <= 10; ++i__) {
	for (j = 1; j <= 5; ++j) {
/* L25: */
	    l5_16.gg[i__ + j * 10 - 11] = a[i__ + j * 10 - 11];
	}
    }
    l20_7.lex = false;
    l20_7.xex[0] = .299999999948;
    l20_7.xex[1] = .333467606492;
    l20_7.xex[2] = .400000000107;
    l20_7.xex[3] = .42831010474;
    l20_7.xex[4] = .223964873676;
    l20_7.fex = -32.3486789716;
    return 0;
labelL2:
    t = 0.;
    for (i__ = 1; i__ <= 5; ++i__) {
	t1 = 0.;
	for (j = 1; j <= 5; ++j) {
/* L21: */
	    t1 += c__[j + i__ * 5 - 6] * l2_4.x[i__ - 1] * l2_4.x[j - 1];
	}
/* labelL20: */
/* Computing 3rd power */
	d__1 = l2_4.x[i__ - 1];
	t = t + e[i__ - 1] * l2_4.x[i__ - 1] + d__[i__ - 1] * (d__1 * (d__1 * 
		d__1)) + t1;
    }
    l6_1.fx = t;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 5; ++i__) {
	t = 0.;
	for (j = 1; j <= 5; ++j) {
/* L22: */
	    t += (c__[i__ + j * 5 - 6] + c__[j + i__ * 5 - 6]) * l2_4.x[j - 1]
		    ;
	}
/* L23: */
/* Computing 2nd power */
	d__1 = l2_4.x[i__ - 1];
	l4_4.gf[i__ - 1] = e[i__ - 1] + t + d__[i__ - 1] * 3. * (d__1 * d__1);
    }
    return 0;
labelL4:
    for (i__ = 1; i__ <= 10; ++i__) {
/* L24: */
	if (l9_10.index1[i__ - 1]) {
	    l3_9.g[i__ - 1] = a[i__ - 1] * l2_4.x[0] + a[i__ + 9] * l2_4.x[1] 
		    + a[i__ + 19] * l2_4.x[2] + a[i__ + 29] * l2_4.x[3] + a[
		    i__ + 39] * l2_4.x[4] - b[i__ - 1];
	}
    }
labelL5:
    return 0;
} /* tp86_ */


/* Subroutine */ int tp87_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real a, b, c__, d__, e;
    static int i__;
    static Real f1, f2, v1, v2, v3, v4, v5;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 6;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 4;
    l2_5.x[0] = 390.;
    l2_5.x[1] = 1e3;
    l2_5.x[2] = 419.5;
    l2_5.x[3] = 340.5;
    l2_5.x[4] = 198.175;
    l2_5.x[5] = .5;
    for (i__ = 1; i__ <= 6; ++i__) {
	l11_5.lxl[i__ - 1] = true;
/* labelL6: */
	l12_5.lxu[i__ - 1] = true;
    }
    l13_5.xl[0] = 0.;
    l13_5.xl[1] = 0.;
    l13_5.xl[2] = 340.;
    l13_5.xl[3] = 340.;
    l13_5.xl[4] = -1e3;
    l13_5.xl[5] = 0.;
    l14_5.xu[0] = 400.;
    l14_5.xu[1] = 1e3;
    l14_5.xu[2] = 420.;
    l14_5.xu[3] = 420.;
    l14_5.xu[4] = 1e3;
    l14_5.xu[5] = .5236;
    a = 131.078;
    b = 1.48477;
    c__ = .90798;
    d__ =std::cos(1.47588);
    e = std::sin(1.47588);
    for (i__ = 3; i__ <= 6; ++i__) {
/* L70: */
	l4_5.gf[i__ - 1] = 0.;
    }
    l5_8.gg[0] = -1.;
    l5_8.gg[4] = 0.;
    l5_8.gg[16] = 0.;
    l5_8.gg[1] = 0.;
    l5_8.gg[5] = -1.;
    l5_8.gg[17] = 0.;
    l5_8.gg[2] = 0.;
    l5_8.gg[6] = 0.;
    l5_8.gg[18] = -1.;
    l5_8.gg[3] = 0.;
    l5_8.gg[7] = 0.;
    l5_8.gg[19] = 0.;
    l20_4.lex = false;
    l20_4.xex[0] = 107.811937779;
    l20_4.xex[1] = 196.318606955;
    l20_4.xex[2] = 373.830728516;
    l20_4.xex[3] = 420.;
    l20_4.xex[4] = 21.3071293896;
    l20_4.xex[5] = .153291953422;
    l20_4.fex = 8927.59773493;
    return 0;
labelL2:
    if (l2_5.x[0] - 300. >= 0.) {
	goto L32;
    } else {
	goto L31;
    }
L31:
    f1 = l2_5.x[0] * 30.;
    goto L33;
L32:
    f1 = l2_5.x[0] * 31.;
L33:
    if (l2_5.x[1] - 100. >= 0.) {
	goto L35;
    } else {
	goto L34;
    }
L34:
    f2 = l2_5.x[1] * 28.;
    goto L46;
L35:
    if (l2_5.x[1] - 200. >= 0.) {
	goto L37;
    } else {
	goto L36;
    }
L36:
    f2 = l2_5.x[1] * 29.;
    goto L46;
L37:
    f2 = l2_5.x[1] * 30.;
L46:
    l6_1.fx = f1 + f2;
    return 0;
labelL3:
    if (l2_5.x[0] - 300. >= 0.) {
	goto L39;
    } else {
	goto L38;
    }
L38:
    l4_5.gf[0] = 30.;
    goto L40;
L39:
    l4_5.gf[0] = 31.;
L40:
    if (l2_5.x[1] - 100. >= 0.) {
	goto L42;
    } else {
	goto L41;
    }
L41:
    l4_5.gf[1] = 28.;
    goto L45;
L42:
    if (l2_5.x[1] - 200. >= 0.) {
	goto L44;
    } else {
	goto L43;
    }
L43:
    l4_5.gf[1] = 29.;
    goto L45;
L44:
    l4_5.gf[1] = 30.;
L45:
    return 0;
labelL4:
    if (l9_7.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_5.x[2];
	l3_6.g[0] = -l2_5.x[0] + 300. - l2_5.x[2] * l2_5.x[3] / a *std::cos(b - 
		l2_5.x[5]) + c__ * (d__1 * d__1) / a * d__;
    }
    if (l9_7.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_5.x[3];
	l3_6.g[1] = -l2_5.x[1] - l2_5.x[2] * l2_5.x[3] / a *std::cos(b + l2_5.x[5]
		) + c__ * (d__1 * d__1) / a * d__;
    }
    if (l9_7.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_5.x[3];
	l3_6.g[2] = -l2_5.x[4] - l2_5.x[2] * l2_5.x[3] / a * std::sin(b + l2_5.x[5]
		) + c__ * (d__1 * d__1) / a * e;
    }
    if (l9_7.index1[3]) {
/* Computing 2nd power */
	d__1 = l2_5.x[2];
	l3_6.g[3] = 200. - l2_5.x[2] * l2_5.x[3] / a * std::sin(b - l2_5.x[5]) + 
		c__ * (d__1 * d__1) / a * e;
    }
    return 0;
labelL5:
    v1 = 1. / a;
    if (! (l10_7.index2[0] || l10_7.index2[3])) {
	goto L8;
    }
    v2 = b - l2_5.x[5];
    v3 =std::cos(v2) * v1;
    v4 = std::sin(v2) * v1;
    v5 = c__ * 2. * l2_5.x[2] * v1;
    if (! l10_7.index2[0]) {
	goto L7;
    }
    l5_8.gg[8] = -l2_5.x[3] * v3 + v5 * d__;
    l5_8.gg[12] = -l2_5.x[2] * v3;
    l5_8.gg[20] = -l2_5.x[2] * l2_5.x[3] * v4;
L7:
    if (! l10_7.index2[3]) {
	goto L8;
    }
    l5_8.gg[11] = -l2_5.x[3] * v4 + v5 * e;
    l5_8.gg[15] = -l2_5.x[2] * v4;
    l5_8.gg[23] = l2_5.x[2] * l2_5.x[3] * v3;
L8:
    if (! (l10_7.index2[1] || l10_7.index2[2])) {
	goto labelL10;
    }
    v2 = b + l2_5.x[5];
    v3 =std::cos(v2) * v1;
    v4 = std::sin(v2) * v1;
    v5 = c__ * 2. * l2_5.x[3] * v1;
    if (! l10_7.index2[1]) {
	goto labelL9;
    }
    l5_8.gg[9] = -l2_5.x[3] * v3;
    l5_8.gg[13] = -l2_5.x[2] * v3 + v5 * d__;
    l5_8.gg[21] = l2_5.x[2] * l2_5.x[3] * v4;
labelL9:
    if (! l10_7.index2[2]) {
	goto labelL10;
    }
    l5_8.gg[10] = -l2_5.x[3] * v4;
    l5_8.gg[14] = -l2_5.x[2] * v4 + v5 * e;
    l5_8.gg[22] = -l2_5.x[2] * l2_5.x[3] * v3;
labelL10:
    return 0;
} /* tp87_ */


/* Subroutine */ int tp88_0_(int n__, int *mode)
{
    /* System generated locals */
    int i__1, i__2;
    Real d__1;

    int pow_ii(int*, int*);

    /* Local variables */
    static Real a[30];
    static int i__, j;
    static Real t[6], u, w, z__, a1, intko;
    static int n1;
    static Real v1, v2, v3, pi, dv[6];
    Real gleich_(Real *);
    static Real dz[6];
    static Real dcosko[30];
    static Real ep1;
    static int kn1;
    static Real mue[30], rho[30];

    switch(n__) {
	case 1: goto L_tp89;
	case 2: goto L_tp90;
	case 3: goto L_tp91;
	case 4: goto L_tp92;
	}

    kn1 = 1;
    l20_4.xex[0] = 1.0743187294;
    l20_4.xex[1] = -.456613707247;
    l20_4.fex = 1.36265680997;
    goto L7;

L_tp89:
    kn1 = 2;
    l20_4.xex[0] = 1.07431872754;
    l20_4.xex[1] = -.456613706239;
    l20_4.xex[2] = 3.00836097604e-11;
    l20_4.fex = 1.36265680508;
    goto L7;

L_tp90:
    kn1 = 3;
    l20_4.xex[0] = .708479399007;
    l20_4.xex[1] = 2.37919269592e-5;
    l20_4.xex[2] = .807599939006;
    l20_4.xex[3] = -.456613723294;
    l20_4.fex = 1.36265681317;
    goto L7;

L_tp91:
    kn1 = 4;
    l20_4.xex[0] = .701892928031;
    l20_4.xex[1] = 2.21084326516e-12;
    l20_4.xex[2] = .813330836201;
    l20_4.xex[3] = .456613707134;
    l20_4.xex[4] = 8.99937588382e-12;
    l20_4.fex = 1.3626568091;
    goto L7;

L_tp92:
    kn1 = 5;
    l20_4.xex[0] = .494144465323;
    l20_4.xex[1] = -1.03530473697e-5;
    l20_4.xex[2] = .61495083955;
    l20_4.xex[3] = -2.42186612731e-6;
    l20_4.xex[4] = .729258528936;
    l20_4.xex[5] = -.456613099133;
    l20_4.fex = 1.36265681213;
L7:
    if (*mode - 2 >= 0) {
	goto L17;
    } else {
	goto labelL1;
    }
labelL1:
    l1_1.n = kn1 + 1;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_5.x[(i__ << 1) - 2] = .5;
/* labelL11: */
	l2_5.x[(i__ << 1) - 1] = -.5;
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	l11_5.lxl[i__ - 1] = true;
	l12_5.lxu[i__ - 1] = true;
	l13_5.xl[i__ - 1] = -10.;
/* labelL6: */
	l14_5.xu[i__ - 1] = 10.;
    }
    pi = std::atan(1.) * 4.;
    for (i__ = 1; i__ <= 30; ++i__) {
	z__ = pi * (Real) (i__ - 1);
	mue[i__ - 1] = gleich_(&z__);
	v1 = std::sin(mue[i__ - 1]);
	v2 =std::cos(mue[i__ - 1]);
/* Computing 2nd power */
	d__1 = mue[i__ - 1];
	dcosko[i__ - 1] = (v1 / mue[i__ - 1] - v2) / (d__1 * d__1);
/* labelL10: */
	a[i__ - 1] = v1 * 2. / (mue[i__ - 1] + v1 * v2);
    }
    intko = .13333333333333333;
    l20_4.lex = false;
    return 0;
L17:
    if (*mode - 4 >= 0) {
	goto L18;
    } else {
	goto L19;
    }
L18:
    n1 = l1_1.n - 1;
/* Computing 2nd power */
    d__1 = l2_5.x[l1_1.n - 1];
    t[l1_1.n - 1] = d__1 * d__1;
    i__1 = n1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__1 = l2_5.x[l1_1.n - i__ - 1];
	t[l1_1.n - i__ - 1] = t[l1_1.n - i__] + d__1 * d__1;
    }
    v1 = 0.;
    for (j = 1; j <= 30; ++j) {
	v2 = mue[j - 1];
/* Computing 2nd power */
	d__1 = v2;
	v3 = -(d__1 * d__1);
	rho[j - 1] = (Real) pow_ii(&c_n1, &l1_1.n);
	i__1 = n1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ep1 = 0.;
	    a1 = v3 * t[l1_1.n + 1 - i__ - 1];
	    if (a1 > -100.) {
		ep1 = std::exp(a1);
	    }
/* labelL14: */
	    i__2 = l1_1.n - i__;
	    rho[j - 1] += (Real) pow_ii(&c_n1, &i__2) * 2. * ep1;
	}
	ep1 = 0.;
	a1 = v3 * t[0];
	if (a1 > -100.) {
	    ep1 = std::exp(a1);
	}
	rho[j - 1] = (rho[j - 1] + ep1) / v3;
/* labelL13: */
	v1 -= v3 * a[j - 1] * rho[j - 1] * (v2 * std::sin(v2) * rho[j - 1] - 
		dcosko[j - 1] * 2.);
    }
L19:
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL2:
    u = 0.;
    i__2 = l1_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* labelL20: */
/* Computing 2nd power */
	d__1 = l2_5.x[i__ - 1];
	u += d__1 * d__1;
    }
    l6_1.fx = u;
    return 0;
labelL3:
    i__2 = l1_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L21: */
	l4_5.gf[i__ - 1] = l2_5.x[i__ - 1] * 2.;
    }
    return 0;
labelL4:
    l3_1.g[0] = 1e-4 - v1 - intko;
    return 0;
labelL5:
    i__2 = l1_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L22: */
	dv[i__ - 1] = 0.;
    }
    for (j = 1; j <= 30; ++j) {
	w = mue[j - 1];
/* Computing 2nd power */
	d__1 = w;
	v1 = d__1 * d__1 * a[j - 1] * (w * std::sin(w) * rho[j - 1] - dcosko[j - 1]
		);
	ep1 = 0.;
/* Computing 2nd power */
	d__1 = mue[j - 1];
	a1 = -(d__1 * d__1) * t[0];
	if (a1 > -100.) {
	    ep1 = std::exp(a1);
	}
	dz[0] = ep1;
	dv[0] += dz[0] * v1;
	i__2 = l1_1.n;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    ep1 = 0.;
/* Computing 2nd power */
	    d__1 = mue[j - 1];
	    a1 = -(d__1 * d__1) * t[i__ - 1];
	    if (a1 > -100.) {
		ep1 = std::exp(a1);
	    }
	    i__1 = i__ + 1;
	    dz[i__ - 1] = dz[i__ - 2] + (Real) pow_ii(&c_n1, &i__1) * 
		    2. * ep1;
/* L23: */
	    dv[i__ - 1] += dz[i__ - 1] * v1;
	}
/* L25: */
    }
    i__2 = l1_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L24: */
	l5_3.gg[i__ - 1] = dv[i__ - 1] * -4. * l2_5.x[i__ - 1];
    }
    return 0;
} /* tp88_ */

/* Subroutine */ int tp88_(int *mode)
{
    return tp88_0_(0, mode);
    }

/* Subroutine */ int tp89_(int *mode)
{
    return tp88_0_(1, mode);
    }

/* Subroutine */ int tp90_(int *mode)
{
    return tp88_0_(2, mode);
    }

/* Subroutine */ int tp91_(int *mode)
{
    return tp88_0_(3, mode);
    }

/* Subroutine */ int tp92_(int *mode)
{
    return tp88_0_(4, mode);
    }


/* Subroutine */ int tp93_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;
    static Real v1, v2, v3, v4, v5, v6, v7, v8, v9;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 6;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_5.x[0] = 5.54;
    l2_5.x[1] = 4.4;
    l2_5.x[2] = 12.02;
    l2_5.x[3] = 11.82;
    l2_5.x[4] = .702;
    l2_5.x[5] = .852;
    for (i__ = 1; i__ <= 6; ++i__) {
	l11_5.lxl[i__ - 1] = true;
	l12_5.lxu[i__ - 1] = false;
/* labelL6: */
	l13_5.xl[i__ - 1] = 0.;
    }
    l20_4.lex = false;
    l20_4.xex[0] = 5.33266639884;
    l20_4.xex[1] = 4.65674439073;
    l20_4.xex[2] = 10.4329901123;
    l20_4.xex[3] = 12.0823085893;
    l20_4.xex[4] = .752607369745;
    l20_4.xex[5] = .87865083685;
    l20_4.fex = 135.075961229;
    return 0;
labelL2:
    v1 = l2_5.x[0] + l2_5.x[1] + l2_5.x[2];
    v2 = l2_5.x[0] + l2_5.x[1] * 1.57 + l2_5.x[3];
    v3 = l2_5.x[0] * l2_5.x[3];
    v4 = l2_5.x[2] * l2_5.x[1];
/* Computing 2nd power */
    d__1 = l2_5.x[4];
/* Computing 2nd power */
    d__2 = l2_5.x[5];
    l6_1.fx = v3 * .0204 * v1 + v4 * .0187 * v2 + v3 * .0607 * v1 * (d__1 * 
	    d__1) + v4 * .0437 * v2 * (d__2 * d__2);
    return 0;
labelL3:
    v1 = l2_5.x[0] * l2_5.x[3];
    v2 = l2_5.x[1] * l2_5.x[2];
    v3 = l2_5.x[1] + l2_5.x[2];
    v4 = l2_5.x[0] + l2_5.x[3];
    v5 = v3 + l2_5.x[0];
    v6 = l2_5.x[1] * 1.57 + v4;
/* Computing 2nd power */
    d__1 = l2_5.x[4];
    v7 = l2_5.x[3] * (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_5.x[5];
    v8 = l2_5.x[2] * (d__1 * d__1);
    v9 = l2_5.x[0] * .0607 * v7;
    l4_5.gf[0] = v1 * .0408 + l2_5.x[3] * .0204 * v3 + v2 * .0187 + v9 * 2. + 
	    v7 * .0607 * v3 + l2_5.x[1] * .0437 * v8;
    l4_5.gf[1] = v1 * .0204 + v2 * .058718 + l2_5.x[2] * .0187 * v4 + v9 + 
	    l2_5.x[1] * .137218 * v8 + v8 * .0437 * v4;
/* Computing 2nd power */
    d__1 = l2_5.x[5];
    l4_5.gf[2] = v1 * .0204 + l2_5.x[1] * .0187 * v6 + v9 + l2_5.x[1] * .0437 
	    * (d__1 * d__1) * v6;
/* Computing 2nd power */
    d__1 = l2_5.x[4];
    l4_5.gf[3] = l2_5.x[0] * .0204 * v5 + v2 * .0187 + l2_5.x[1] * .0437 * v8 
	    + l2_5.x[0] * .0607 * (d__1 * d__1) * v5;
    l4_5.gf[4] = v1 * .1214 * l2_5.x[4] * v5;
    l4_5.gf[5] = l2_5.x[5] * .0874 * v2 * v6;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_5.x[0] * .001 * l2_5.x[1] * l2_5.x[2] * l2_5.x[3] * 
		l2_5.x[4] * l2_5.x[5] - 2.07;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_5.x[4];
/* Computing 2nd power */
	d__2 = l2_5.x[5];
	l3_2.g[1] = 1. - l2_5.x[0] * 6.2e-4 * l2_5.x[3] * (d__1 * d__1) * (
		l2_5.x[0] + l2_5.x[1] + l2_5.x[2]) - l2_5.x[1] * 5.8e-4 * 
		l2_5.x[2] * (d__2 * d__2) * (l2_5.x[0] + l2_5.x[1] * 1.57 + 
		l2_5.x[3]);
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    v1 = l2_5.x[0] * l2_5.x[1] * l2_5.x[2] * .001;
    v2 = l2_5.x[3] * l2_5.x[4] * l2_5.x[5] * .001;
    l5_7.gg[0] = l2_5.x[1] * l2_5.x[2] * v2;
    l5_7.gg[2] = l2_5.x[0] * l2_5.x[2] * v2;
    l5_7.gg[4] = l2_5.x[0] * l2_5.x[1] * v2;
    l5_7.gg[6] = l2_5.x[4] * l2_5.x[5] * v1;
    l5_7.gg[8] = l2_5.x[3] * l2_5.x[5] * v1;
    l5_7.gg[10] = l2_5.x[3] * l2_5.x[4] * v1;
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
/* Computing 2nd power */
    d__1 = l2_5.x[4];
    v1 = -(d__1 * d__1) * 6.2e-4;
/* Computing 2nd power */
    d__1 = l2_5.x[5];
    v2 = -(d__1 * d__1) * 5.8e-4;
    v3 = l2_5.x[0] + l2_5.x[1] + l2_5.x[2];
    v4 = l2_5.x[0] + l2_5.x[1] * 1.57 + l2_5.x[3];
    v5 = v1 * v3;
    v6 = v2 * v4;
    v7 = v1 * l2_5.x[0] * l2_5.x[3];
    v8 = v2 * l2_5.x[1] * l2_5.x[2];
    l5_7.gg[1] = v7 + v5 * l2_5.x[3] + v8;
    l5_7.gg[3] = v7 + v6 * l2_5.x[2] + v8 * 1.57;
    l5_7.gg[5] = v7 + v6 * l2_5.x[1];
    l5_7.gg[7] = v5 * l2_5.x[0] + v8;
    l5_7.gg[9] = l2_5.x[0] * -.00124 * l2_5.x[3] * l2_5.x[4] * v3;
    l5_7.gg[11] = l2_5.x[1] * -.00116 * l2_5.x[2] * l2_5.x[5] * v4;
L8:
    return 0;
} /* tp93_ */


/* Subroutine */ int tp95_0_(int n__, int *mode)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static Real b[4];
    static int i__, kn1;

    switch(n__) {
	case 1: goto L_tp96;
	case 2: goto L_tp97;
	case 3: goto L_tp98;
	}

    kn1 = 1;
    l20_4.xex[0] = -4.76149332788e-12;
    l20_4.xex[1] = -3.55239427962e-11;
    l20_4.xex[2] = -7.02611041315e-11;
    l20_4.xex[3] = -1.71856469485e-11;
    l20_4.xex[4] = -7.48993551642e-11;
    l20_4.xex[5] = .00332330328254;
    l20_4.fex = .0156195144282;
    goto labelL11;

L_tp96:
    kn1 = 2;
    l20_4.xex[0] = -5.19722825686e-12;
    l20_4.xex[1] = -3.87748184662e-11;
    l20_4.xex[2] = -7.66908552858e-11;
    l20_4.xex[3] = -1.87583442974e-11;
    l20_4.xex[4] = -8.17535626869e-11;
    l20_4.xex[5] = .00332330328612;
    l20_4.fex = .0156195134384;
    goto labelL11;

L_tp97:
    kn1 = 3;
    l20_4.xex[0] = .268564912352;
    l20_4.xex[1] = 0.;
    l20_4.xex[2] = 0.;
    l20_4.xex[3] = 0.;
    l20_4.xex[4] = .028;
    l20_4.xex[5] = .0134000000001;
    l20_4.fex = 3.13580912311;
    goto labelL11;

L_tp98:
    kn1 = 4;
    l20_4.xex[0] = .268564912323;
    l20_4.xex[1] = 0.;
    l20_4.xex[2] = 0.;
    l20_4.xex[3] = 0.;
    l20_4.xex[4] = .028;
    l20_4.xex[5] = .0134000000001;
    l20_4.fex = 3.13580912299;
labelL11:
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 6;
    l1_1.nili = 0;
    l1_1.ninl = 4;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 6; ++i__) {
	l2_5.x[i__ - 1] = 0.;
	l11_5.lxl[i__ - 1] = true;
	l12_5.lxu[i__ - 1] = true;
/* labelL6: */
	l13_5.xl[i__ - 1] = 0.;
    }
    l14_5.xu[0] = .31;
    l14_5.xu[1] = .046;
    l14_5.xu[2] = .068;
    l14_5.xu[3] = .042;
    l14_5.xu[4] = .028;
    l14_5.xu[5] = .0134;
    l4_5.gf[0] = 4.3;
    l4_5.gf[1] = 31.8;
    l4_5.gf[2] = 63.3;
    l4_5.gf[3] = 15.8;
    l4_5.gf[4] = 68.5;
    l4_5.gf[5] = 4.7;
    l5_8.gg[4] = 38.2;
    l5_8.gg[5] = 36.8;
    l5_8.gg[2] = 0.;
    l5_8.gg[6] = -273.;
    l5_8.gg[10] = 0.;
    l5_8.gg[22] = 0.;
    l5_8.gg[7] = -311.;
    l5_8.gg[11] = 0.;
    l5_8.gg[15] = 587.;
    l5_8.gg[19] = 391.;
    if (kn1 - 2 <= 0) {
	goto labelL20;
    } else {
	goto L21;
    }
labelL20:
    b[0] = 4.97;
    b[1] = -1.88;
    goto L22;
L21:
    b[0] = 32.97;
    b[1] = 25.12;
L22:
    if ((i__1 = kn1 - 3) < 0) {
	goto L25;
    } else if (i__1 == 0) {
	goto L24;
    } else {
	goto L23;
    }
L23:
    b[2] = -124.08;
    b[3] = -173.02;
    goto L27;
L24:
    b[2] = -29.08;
    b[3] = -78.02;
    goto L27;
L25:
    if (kn1 - 2 >= 0) {
	goto L26;
    } else {
	goto L24;
    }
L26:
    b[2] = -69.08;
    b[3] = -118.02;
L27:
    l20_4.lex = false;
    return 0;
labelL2:
    l6_1.fx = l2_5.x[0] * 4.3 + l2_5.x[1] * 31.8 + l2_5.x[2] * 63.3 + l2_5.x[
	    3] * 15.8 + l2_5.x[4] * 68.5 + l2_5.x[5] * 4.7;
labelL3:
    return 0;
labelL4:
    if (l9_7.index1[0]) {
	l3_6.g[0] = l2_5.x[0] * (float)17.1 + l2_5.x[1] * (float)38.2 + 
		l2_5.x[2] * (float)204.2 + l2_5.x[3] * (float)212.3 + l2_5.x[
		4] * (float)623.4 + l2_5.x[5] * 1495.5 - l2_5.x[0] * (float)
		169. * l2_5.x[2] - l2_5.x[2] * 3580. * l2_5.x[4] - l2_5.x[3] *
		 3810. * l2_5.x[4] - l2_5.x[3] * 18500. * l2_5.x[5] - l2_5.x[
		4] * 24300. * l2_5.x[5] - b[0];
    }
    if (l9_7.index1[1]) {
	l3_6.g[1] = l2_5.x[0] * (float)17.9 + l2_5.x[1] * (float)36.8 + 
		l2_5.x[2] * (float)113.9 + l2_5.x[3] * (float)169.7 + l2_5.x[
		4] * (float)337.8 + l2_5.x[5] * 1385.2 - l2_5.x[0] * (float)
		139. * l2_5.x[2] - l2_5.x[3] * 2450. * l2_5.x[4] - l2_5.x[3] *
		 16600. * l2_5.x[5] - l2_5.x[4] * 17200. * l2_5.x[5] - b[1];
    }
    if (l9_7.index1[2]) {
	l3_6.g[2] = l2_5.x[1] * (float)-273. - l2_5.x[3] * (float)70. - 
		l2_5.x[4] * (float)819. + l2_5.x[3] * 2.6e4 * l2_5.x[4] - b[2]
		;
    }
    if (l9_7.index1[3]) {
	l3_6.g[3] = l2_5.x[0] * 159.9 - l2_5.x[1] * 311. + l2_5.x[3] * 587. + 
		l2_5.x[4] * (float)391. + l2_5.x[5] * 2198. - l2_5.x[0] * 
		1.4e4 * l2_5.x[5] - b[3];
    }
    return 0;
labelL5:
    if (! l10_7.index2[0]) {
	goto L7;
    }
    l5_8.gg[0] = (float)17.1 - l2_5.x[2] * (float)169.;
    l5_8.gg[8] = (float)204.2 - l2_5.x[0] * (float)169. - l2_5.x[4] * 3580.;
    l5_8.gg[12] = (float)212.3 - l2_5.x[4] * 3810. - l2_5.x[5] * 18500.;
    l5_8.gg[16] = (float)623.4 - l2_5.x[2] * 3580. - l2_5.x[3] * 3810. - 
	    l2_5.x[5] * 24300.;
    l5_8.gg[20] = 1495.5 - l2_5.x[3] * 18500. - l2_5.x[4] * 24300.;
L7:
    if (! l10_7.index2[1]) {
	goto L8;
    }
    l5_8.gg[1] = (float)17.9 - l2_5.x[2] * (float)139.;
    l5_8.gg[9] = (float)113.9 - l2_5.x[0] * (float)139.;
    l5_8.gg[13] = (float)169.7 - l2_5.x[4] * 2450. - l2_5.x[5] * 16600.;
    l5_8.gg[17] = (float)337.8 - l2_5.x[3] * 2450. - l2_5.x[5] * 17200.;
    l5_8.gg[21] = 1385.2 - l2_5.x[3] * 16600. - l2_5.x[4] * 17200.;
L8:
    if (! l10_7.index2[2]) {
	goto labelL9;
    }
    l5_8.gg[14] = l2_5.x[4] * 2.6e4 - 70.;
    l5_8.gg[18] = l2_5.x[3] * 2.6e4 - 819.;
labelL9:
    if (! l10_7.index2[3]) {
	goto labelL10;
    }
    l5_8.gg[3] = (float)159.9 - l2_5.x[5] * 1.4e4;
    l5_8.gg[23] = 2198. - l2_5.x[0] * 1.4e4;
labelL10:
    return 0;
} /* tp95_ */

/* Subroutine */ int tp95_(int *mode)
{
    return tp95_0_(0, mode);
    }

/* Subroutine */ int tp96_(int *mode)
{
    return tp95_0_(1, mode);
    }

/* Subroutine */ int tp97_(int *mode)
{
    return tp95_0_(2, mode);
    }

/* Subroutine */ int tp98_(int *mode)
{
    return tp95_0_(3, mode);
    }


/* Subroutine */ int tp99_(int *mode)
{
    /* System generated locals */
    int i__1;
    Real d__1;

    /* Local variables */
    static Real a[8];
    static int i__, j;
    static Real p[8], q[8], r__[8], s[8], t[8];
    static int i1;
    static Real v1, v2, v3, v4, dp[56]	/* was [8][7] */, dq[56]	
	    /* was [8][7] */, dr[56]	/* was [8][7] */, ds[56]	/* 
	    was [8][7] */;

    if ((i__1 = *mode - 2) < 0) {
	goto labelL1;
    } else if (i__1 == 0) {
	goto L18;
    } else {
	goto L17;
    }
labelL1:
    l1_1.n = 7;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 2;
    for (i__ = 1; i__ <= 7; ++i__) {
	l2_6.x[i__ - 1] = .5;
	l11_6.lxl[i__ - 1] = true;
	l12_6.lxu[i__ - 1] = true;
	l13_6.xl[i__ - 1] = 0.;
/* labelL6: */
	l14_6.xu[i__ - 1] = 1.58;
    }
    a[0] = 0.;
    a[1] = 50.;
    a[2] = 50.;
    a[3] = 75.;
    a[4] = 75.;
    a[5] = 75.;
    a[6] = 100.;
    a[7] = 100.;
    t[0] = 0.;
    t[1] = 25.;
    t[2] = 50.;
    t[3] = 100.;
    t[4] = 150.;
    t[5] = 200.;
    t[6] = 290.;
    t[7] = 380.;
    p[0] = 0.;
    q[0] = 0.;
    r__[0] = 0.;
    s[0] = 0.;
    for (j = 1; j <= 7; ++j) {
	i__1 = j;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dp[i__ + (j << 3) - 9] = 0.;
	    dq[i__ + (j << 3) - 9] = 0.;
	    dr[i__ + (j << 3) - 9] = 0.;
/* L31: */
	    ds[i__ + (j << 3) - 9] = 0.;
	}
    }
    l20_8.lex = false;
    l20_8.xex[0] = .542460319142;
    l20_8.xex[1] = .529015870015;
    l20_8.xex[2] = .508450583169;
    l20_8.xex[3] = .480269265187;
    l20_8.xex[4] = .451235157238;
    l20_8.xex[5] = .409187805755;
    l20_8.xex[6] = .352784693565;
    l20_8.fex = -831079891.516;
    return 0;
L17:
    if (*mode - 4 <= 0) {
	goto L18;
    } else {
	goto L19;
    }
L18:
    for (i__ = 2; i__ <= 8; ++i__) {
	i1 = i__ - 1;
	v1 = a[i__ - 1] * std::sin(l2_6.x[i1 - 1]) - 32.;
	v2 = a[i__ - 1] *std::cos(l2_6.x[i1 - 1]);
	v3 = t[i__ - 1] - t[i1 - 1];
/* Computing 2nd power */
	d__1 = v3;
	v4 = d__1 * d__1 * .5;
	p[i__ - 1] = v2 * v4 + v3 * r__[i1 - 1] + p[i1 - 1];
	q[i__ - 1] = v1 * v4 + v3 * s[i1 - 1] + q[i1 - 1];
	r__[i__ - 1] = v2 * v3 + r__[i1 - 1];
/* L30: */
	s[i__ - 1] = v1 * v3 + s[i1 - 1];
    }
    if (*mode - 3 != 0) {
	goto L40;
    } else {
	goto L19;
    }
L19:
    for (i__ = 2; i__ <= 8; ++i__) {
	for (j = 1; j <= 7; ++j) {
	    if ((i__1 = j - i__ + 1) < 0) {
		goto L33;
	    } else if (i__1 == 0) {
		goto L32;
	    } else {
		goto L34;
	    }
L32:
	    i1 = i__ - 1;
	    v1 = a[i__ - 1] * std::sin(l2_6.x[i1 - 1]);
	    v2 = a[i__ - 1] *std::cos(l2_6.x[i1 - 1]);
	    v3 = t[i__ - 1] - t[i1 - 1];
/* Computing 2nd power */
	    d__1 = v3;
	    v4 = d__1 * d__1 * .5;
	    dp[i__ + (i1 << 3) - 9] = -v1 * v4 + v3 * dr[i1 + (i1 << 3) - 9] 
		    + dp[i1 + (i1 << 3) - 9];
	    dq[i__ + (i1 << 3) - 9] = v2 * v4 + v3 * ds[i1 + (i1 << 3) - 9] + 
		    dq[i1 + (i1 << 3) - 9];
	    dr[i__ + (i1 << 3) - 9] = -v1 * v3 + dr[i1 + (i1 << 3) - 9];
	    ds[i__ + (i1 << 3) - 9] = v2 * v3 + ds[i1 + (i1 << 3) - 9];
	    goto L34;
L33:
	    i1 = i__ - 1;
	    v1 = t[i__ - 1] - t[i1 - 1];
	    dp[i__ + (j << 3) - 9] = v1 * dr[i1 + (j << 3) - 9] + dp[i1 + (j 
		    << 3) - 9];
	    dq[i__ + (j << 3) - 9] = v1 * ds[i1 + (j << 3) - 9] + dq[i1 + (j 
		    << 3) - 9];
	    dr[i__ + (j << 3) - 9] = dr[i1 + (j << 3) - 9];
	    ds[i__ + (j << 3) - 9] = ds[i1 + (j << 3) - 9];
L34:
	    ;
	}
    }
L40:
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL2:
/* Computing 2nd power */
    d__1 = r__[7];
    l6_1.fx = -(d__1 * d__1);
    return 0;
labelL3:
    for (i__ = 1; i__ <= 7; ++i__) {
/* L35: */
	l4_6.gf[i__ - 1] = r__[7] * -2. * dr[(i__ << 3) - 1];
    }
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = q[7] - 1e5;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = s[7] - 1e3;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    for (i__ = 1; i__ <= 7; ++i__) {
/* L36: */
	l5_17.gg[(i__ << 1) - 2] = dq[(i__ << 3) - 1];
    }
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    for (i__ = 1; i__ <= 7; ++i__) {
/* L37: */
	l5_17.gg[(i__ << 1) - 1] = ds[(i__ << 3) - 1];
    }
L8:
    return 0;
} /* tp99_ */


/* Subroutine */ int tp100_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5, d__6, d__7;

    /* Local variables */
    static int i__;
    static Real v1, v2;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 7;
    l1_1.nili = 0;
    l1_1.ninl = 4;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 7; ++i__) {
	l11_6.lxl[i__ - 1] = false;
/* labelL6: */
	l12_6.lxu[i__ - 1] = false;
    }
    l2_6.x[0] = 1.;
    l2_6.x[1] = 2.;
    l2_6.x[2] = 0.;
    l2_6.x[3] = 4.;
    l2_6.x[4] = 0.;
    l2_6.x[5] = 1.;
    l2_6.x[6] = 1.;
    l5_11.gg[8] = -1.;
    l5_11.gg[16] = -5.;
    l5_11.gg[20] = 0.;
    l5_11.gg[24] = 0.;
    l5_11.gg[1] = -7.;
    l5_11.gg[5] = -3.;
    l5_11.gg[13] = -1.;
    l5_11.gg[17] = 1.;
    l5_11.gg[21] = 0.;
    l5_11.gg[25] = 0.;
    l5_11.gg[2] = -23.;
    l5_11.gg[10] = 0.;
    l5_11.gg[14] = 0.;
    l5_11.gg[18] = 0.;
    l5_11.gg[26] = 8.;
    l5_11.gg[15] = 0.;
    l5_11.gg[19] = 0.;
    l5_11.gg[23] = -5.;
    l5_11.gg[27] = 11.;
    l20_8.lex = false;
    l20_8.xex[0] = 2.33049937431;
    l20_8.xex[1] = 1.95137237315;
    l20_8.xex[2] = -.477541392625;
    l20_8.xex[3] = 4.36572623462;
    l20_8.xex[4] = -.624486970475;
    l20_8.xex[5] = 1.03813101881;
    l20_8.xex[6] = 1.59422671137;
    l20_8.fex = 680.630057275;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_6.x[0] - 10.;
/* Computing 2nd power */
    d__2 = l2_6.x[1] - 12.;
/* Computing 4th power */
    d__3 = l2_6.x[2], d__3 *= d__3;
/* Computing 2nd power */
    d__4 = l2_6.x[3] - 11.;
/* Computing 6th power */
    d__5 = l2_6.x[4], d__5 *= d__5;
/* Computing 2nd power */
    d__6 = l2_6.x[5];
/* Computing 4th power */
    d__7 = l2_6.x[6], d__7 *= d__7;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 * 5. + d__3 * d__3 + d__4 * d__4 * 3. 
	    + d__5 * (d__5 * d__5) * 10. + d__6 * d__6 * 7. + d__7 * d__7 - 
	    l2_6.x[5] * 4. * l2_6.x[6] - l2_6.x[5] * 10. - l2_6.x[6] * 8.;
    return 0;
labelL3:
    l4_6.gf[0] = (l2_6.x[0] - 10.) * 2.;
    l4_6.gf[1] = (l2_6.x[1] - 12.) * 10.;
/* Computing 3rd power */
    d__1 = l2_6.x[2];
    l4_6.gf[2] = d__1 * (d__1 * d__1) * 4.;
    l4_6.gf[3] = (l2_6.x[3] - 11.) * 6.;
/* Computing 5th power */
    d__1 = l2_6.x[4], d__2 = d__1, d__1 *= d__1;
    l4_6.gf[4] = d__2 * (d__1 * d__1) * 60.;
    l4_6.gf[5] = l2_6.x[5] * 14. - l2_6.x[6] * 4. - 10.;
/* Computing 3rd power */
    d__1 = l2_6.x[6];
    l4_6.gf[6] = d__1 * (d__1 * d__1) * 4. - l2_6.x[5] * 4. - 8.;
    return 0;
labelL4:
/* Computing 2nd power */
    d__1 = l2_6.x[0];
    v1 = d__1 * d__1 * 2.;
/* Computing 2nd power */
    d__1 = l2_6.x[1];
    v2 = d__1 * d__1;
    if (l9_7.index1[0]) {
/* Computing 2nd power */
	d__1 = v2;
/* Computing 2nd power */
	d__2 = l2_6.x[3];
	l3_6.g[0] = -v1 - d__1 * d__1 * 3. - l2_6.x[2] - d__2 * d__2 * 4. - 
		l2_6.x[4] * 5. + 127.;
    }
    if (l9_7.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_6.x[2];
	l3_6.g[1] = l2_6.x[0] * -7. - l2_6.x[1] * 3. - d__1 * d__1 * 10. - 
		l2_6.x[3] + l2_6.x[4] + 282.;
    }
    if (l9_7.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_6.x[5];
	l3_6.g[2] = l2_6.x[0] * -23. - v2 - d__1 * d__1 * 6. + l2_6.x[6] * 8. 
		+ 196.;
    }
    if (l9_7.index1[3]) {
/* Computing 2nd power */
	d__1 = l2_6.x[2];
	l3_6.g[3] = v1 * -2. - v2 + l2_6.x[0] * 3. * l2_6.x[1] - d__1 * d__1 *
		 2. - l2_6.x[5] * 5. + l2_6.x[6] * 11.;
    }
    return 0;
labelL5:
    if (! l10_7.index2[0]) {
	goto L7;
    }
    l5_11.gg[0] = l2_6.x[0] * -4.;
/* Computing 3rd power */
    d__1 = l2_6.x[1];
    l5_11.gg[4] = d__1 * (d__1 * d__1) * -12.;
    l5_11.gg[12] = l2_6.x[3] * -8.;
L7:
    if (l10_7.index2[1]) {
	l5_11.gg[9] = l2_6.x[2] * -20.;
    }
    if (! l10_7.index2[2]) {
	goto labelL9;
    }
    l5_11.gg[6] = l2_6.x[1] * -2.;
    l5_11.gg[22] = l2_6.x[5] * -12.;
labelL9:
    if (! l10_7.index2[3]) {
	goto labelL10;
    }
    l5_11.gg[3] = l2_6.x[0] * -8. + l2_6.x[1] * 3.;
    l5_11.gg[7] = l2_6.x[1] * -2. + l2_6.x[0] * 3.;
    l5_11.gg[11] = l2_6.x[2] * -4.;
labelL10:
    return 0;
} /* tp100_ */


/* Subroutine */ int tp101_0_(int n__, int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12, d__13;

    /* Local variables */
    static Real fmin[3], a[3];
    static int i__, k, m;
    static Real v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, 
	    v14, v15, gv[7], v16, v17, v18, v19, v20, v21, v22, v23, v24, v25,
	     v26, v27, v28, v29, v30, v31, v32, v33, v34, v35, v36, v37, v38, 
	    v39, v40, v41, v42, v43, v44, v45;
    static int kn1;
    static Real sum;

    switch(n__) {
	case 1: goto L_tp102;
	case 2: goto L_tp103;
	}

    kn1 = 1;
    l20_8.xex[0] = 2.85615855584;
    l20_8.xex[1] = .610823030755;
    l20_8.xex[2] = 2.15081256203;
    l20_8.xex[3] = 4.71287370945;
    l20_8.xex[4] = .999487540961;
    l20_8.xex[5] = 1.34750750498;
    l20_8.xex[6] = .0316527664991;
    l20_8.fex = 1809.76476556;
    goto labelL13;

L_tp102:
    kn1 = 2;
    l20_8.xex[0] = 3.89625319099;
    l20_8.xex[1] = .809358760118;
    l20_8.xex[2] = 2.66438599373;
    l20_8.xex[3] = 4.30091287458;
    l20_8.xex[4] = .853554935267;
    l20_8.xex[5] = 1.09528744459;
    l20_8.xex[6] = .0273104596581;
    l20_8.fex = 911.880571336;
    goto labelL13;

L_tp103:
    kn1 = 3;
    l20_8.xex[0] = 4.39410451026;
    l20_8.xex[1] = .854468738817;
    l20_8.xex[2] = 2.8432303138;
    l20_8.xex[3] = 3.39997866779;
    l20_8.xex[4] = .722926133025;
    l20_8.xex[5] = .87040638184;
    l20_8.xex[6] = .0246388263302;
    l20_8.fex = 543.667958424;
labelL13:
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 7;
    m = kn1;
    a[0] = -.25;
    a[1] = .125;
    a[2] = .5;
    l1_1.nili = 0;
    l1_1.ninl = 6;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 7; ++i__) {
	l2_6.x[i__ - 1] = 6.;
	l13_6.xl[i__ - 1] = .1;
	l14_6.xu[i__ - 1] = 10.;
	l11_6.lxl[i__ - 1] = true;
/* labelL6: */
	l12_6.lxu[i__ - 1] = true;
    }
    l13_6.xl[6] = .01;
    l5_12.gg[24] = 0.;
    l5_12.gg[37] = 0.;
    l5_12.gg[20] = 0.;
    l5_12.gg[33] = 0.;
    l20_8.lex = false;
    return 0;
labelL2:
    for (k = 1; k <= 7; ++k) {
/* L200: */
	if (l2_6.x[k - 1] < 0.) {
	    goto labelL14;
	}
    }
/* Computing 2nd power */
    d__1 = l2_6.x[3];
/* Computing 3rd power */
    d__2 = l2_6.x[5];
/* Computing 2nd power */
    d__3 = l2_6.x[1];
/* Computing 2nd power */
    d__4 = l2_6.x[0];
/* Computing 2nd power */
    d__5 = l2_6.x[4];
/* Computing 2nd power */
    d__6 = l2_6.x[0];
/* Computing 2nd power */
    d__7 = l2_6.x[1];
/* Computing 2nd power */
    d__8 = l2_6.x[5];
    l6_1.fx = l2_6.x[0] * 10. * (d__1 * d__1) * pow_dd(&l2_6.x[6], &a[m - 1]) 
	    / (l2_6.x[1] * (d__2 * (d__2 * d__2))) + l2_6.x[2] * 15. * l2_6.x[
	    3] / (l2_6.x[0] * (d__3 * d__3) * l2_6.x[4] * pow_dd(&l2_6.x[6], &
	    c_b590)) + l2_6.x[1] * 20. * l2_6.x[5] / (d__4 * d__4 * l2_6.x[3] 
	    * (d__5 * d__5)) + d__6 * d__6 * 25. * (d__7 * d__7) * pow_dd(&
	    l2_6.x[4], &c_b590) * l2_6.x[6] / (l2_6.x[2] * (d__8 * d__8));
    return 0;
labelL14:
    sum = 0.;
    for (i__ = 1; i__ <= 7; ++i__) {
/* L40: */
/* Computing 2nd power */
	d__1 = l2_6.x[i__ - 1] - 5.;
	sum += d__1 * d__1;
    }
    fmin[0] = 1800.;
    fmin[1] = 910.;
    fmin[2] = 540.;
    l6_1.fx = sum + 1e3 + fmin[kn1 - 1];
    return 0;
labelL3:
    for (k = 1; k <= 7; ++k) {
/* L201: */
	if (l2_6.x[k - 1] < 0.) {
	    goto L15;
	}
    }
/* Computing 2nd power */
    d__1 = l2_6.x[3];
    v1 = d__1 * d__1 * 10.;
    v2 = pow_dd(&l2_6.x[6], &a[m - 1]);
/* Computing 3rd power */
    d__1 = l2_6.x[5];
    v3 = l2_6.x[1] * (d__1 * (d__1 * d__1));
    v4 = l2_6.x[2] * 15. * l2_6.x[3];
/* Computing 2nd power */
    d__1 = l2_6.x[1];
    v5 = l2_6.x[0] * (d__1 * d__1) * l2_6.x[4] * pow_dd(&l2_6.x[6], &c_b590);
    v6 = l2_6.x[1] * 20. * l2_6.x[5];
/* Computing 2nd power */
    d__1 = l2_6.x[0];
/* Computing 2nd power */
    d__2 = l2_6.x[4];
    v7 = d__1 * d__1 * l2_6.x[3] * (d__2 * d__2);
    v8 = l2_6.x[0] * 25. * l2_6.x[1] * pow_dd(&l2_6.x[4], &c_b590) * l2_6.x[6]
	    ;
/* Computing 2nd power */
    d__1 = l2_6.x[5];
    v9 = l2_6.x[2] * (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_6.x[0];
/* Computing 2nd power */
    d__2 = l2_6.x[1];
    v10 = d__1 * d__1 * 12.5 * (d__2 * d__2) * l2_6.x[6];
    v11 = pow_dd(&l2_6.x[4], &c_b590);
    l4_6.gf[0] = v1 * v2 / v3 - v4 / (l2_6.x[0] * v5) - v6 * 2. / (l2_6.x[0] *
	     v7) + l2_6.x[1] * 2. * v8 / v9;
    l4_6.gf[1] = -v1 * l2_6.x[0] * v2 / (l2_6.x[1] * v3) - v4 * 2. / (l2_6.x[
	    1] * v5) + l2_6.x[5] * 20. / v7 + l2_6.x[0] * 2. * v8 / v9;
    l4_6.gf[2] = l2_6.x[3] * 15. / v5 - l2_6.x[0] * l2_6.x[1] * v8 / (l2_6.x[
	    2] * v9);
    l4_6.gf[3] = l2_6.x[0] * 20. * l2_6.x[3] * v2 / v3 + l2_6.x[2] * 15. / v5 
	    - v6 / (l2_6.x[3] * v7);
    l4_6.gf[4] = -v4 / (l2_6.x[4] * v5) - v6 * 2. / (l2_6.x[4] * v7) + v10 / (
	    v9 * v11);
    l4_6.gf[5] = v1 * -3. * l2_6.x[0] * v2 / (l2_6.x[5] * v3) + l2_6.x[1] * 
	    20. / v7 - v10 * 4. * v11 / (l2_6.x[5] * v9);
    d__1 = a[m - 1] - 1.;
    l4_6.gf[6] = a[m - 1] * v1 * l2_6.x[0] * pow_dd(&l2_6.x[6], &d__1) / v3 - 
	    v4 * .5 / (v5 * l2_6.x[6]) + v8 * l2_6.x[0] * l2_6.x[1] / (l2_6.x[
	    6] * v9);
    return 0;
L15:
    for (i__ = 1; i__ <= 7; ++i__) {
/* L50: */
	l4_6.gf[i__ - 1] = (l2_6.x[i__ - 1] - 5.) * 2.;
    }
    return 0;
labelL4:
    if (l9_6.index1[0]) {
	d__1 = std::abs(l2_6.x[0]);
/* Computing 2nd power */
	d__2 = l2_6.x[5];
/* Computing 3rd power */
	d__3 = l2_6.x[0];
	d__4 = std::abs(l2_6.x[6]);
/* Computing 2nd power */
	d__5 = l2_6.x[2];
	d__6 = std::abs(l2_6.x[5]);
	d__7 = std::abs(l2_6.x[6]);
	d__8 = std::abs(l2_6.x[3]);
	l3_5.g[0] = 1. - pow_dd(&d__1, &c_b590) * .5 * l2_6.x[6] / (l2_6.x[2] 
		* (d__2 * d__2)) - d__3 * (d__3 * d__3) * .7 * l2_6.x[1] * 
		l2_6.x[5] * pow_dd(&d__4, &c_b590) / (d__5 * d__5) - l2_6.x[2]
		 * .2 * pow_dd(&d__6, &c_b933) * pow_dd(&d__7, &c_b934) / (
		l2_6.x[1] * pow_dd(&d__8, &c_b590));
    }
    if (l9_6.index1[1]) {
	d__1 = std::abs(l2_6.x[0]);
/* Computing 2nd power */
	d__2 = l2_6.x[5];
	d__3 = std::abs(l2_6.x[1]);
	d__4 = std::abs(l2_6.x[5]);
/* Computing 2nd power */
	d__5 = l2_6.x[3];
	l3_5.g[1] = 1. - l2_6.x[1] * 1.3 * l2_6.x[5] / (pow_dd(&d__1, &c_b590)
		 * l2_6.x[2] * l2_6.x[4]) - l2_6.x[2] * .8 * (d__2 * d__2) / (
		l2_6.x[3] * l2_6.x[4]) - pow_dd(&d__3, &c_b590) * 3.1 * 
		pow_dd(&d__4, &c_b74) / (l2_6.x[0] * (d__5 * d__5) * l2_6.x[4]
		);
    }
    if (l9_6.index1[2]) {
	d__2 = std::abs(l2_6.x[6]);
	d__3 = std::abs(l2_6.x[2]);
	d__4 = (d__1 = l2_6.x[2] * l2_6.x[6], std::abs(d__1));
	d__5 = std::abs(l2_6.x[2]);
/* Computing 2nd power */
	d__6 = l2_6.x[1];
	l3_5.g[2] = 1. - l2_6.x[0] * 2. * l2_6.x[4] * pow_dd(&d__2, &c_b74) / 
		(pow_dd(&d__3, &c_b940) * l2_6.x[5]) - l2_6.x[1] * .1 * 
		l2_6.x[4] / (pow_dd(&d__4, &c_b590) * l2_6.x[5]) - l2_6.x[1] *
		 pow_dd(&d__5, &c_b590) * l2_6.x[4] / l2_6.x[0] - l2_6.x[2] * 
		.65 * l2_6.x[4] * l2_6.x[6] / (d__6 * d__6 * l2_6.x[5]);
    }
    if (l9_6.index1[3]) {
	d__1 = std::abs(l2_6.x[4]);
	d__2 = std::abs(l2_6.x[6]);
/* Computing 2nd power */
	d__3 = l2_6.x[0];
	d__4 = std::abs(l2_6.x[0]);
/* Computing 2nd power */
	d__5 = l2_6.x[1];
	d__6 = std::abs(l2_6.x[3]);
	d__7 = std::abs(l2_6.x[6]);
	d__8 = std::abs(l2_6.x[4]);
	d__9 = std::abs(l2_6.x[6]);
/* Computing 3rd power */
	d__10 = l2_6.x[0];
/* Computing 2nd power */
	d__11 = l2_6.x[1];
	d__12 = std::abs(l2_6.x[6]);
/* Computing 2nd power */
	d__13 = l2_6.x[2];
	l3_5.g[3] = 1. - l2_6.x[1] * .2 * pow_dd(&d__1, &c_b590) * pow_dd(&
		d__2, &c_b74) / (d__3 * d__3 * l2_6.x[3]) - pow_dd(&d__4, &
		c_b590) * .3 * (d__5 * d__5) * l2_6.x[2] * pow_dd(&d__6, &
		c_b74) * pow_dd(&d__7, &c_b934) / pow_dd(&d__8, &c_b933) - 
		l2_6.x[2] * .4 * l2_6.x[4] * pow_dd(&d__9, &c_b949) / (d__10 *
		 (d__10 * d__10) * (d__11 * d__11)) - l2_6.x[3] * .5 * pow_dd(
		&d__12, &c_b590) / (d__13 * d__13);
    }
    if (l9_6.index1[4]) {
/* Computing 2nd power */
	d__1 = l2_6.x[3];
	d__2 = std::abs(l2_6.x[6]);
/* Computing 3rd power */
	d__3 = l2_6.x[5];
/* Computing 2nd power */
	d__4 = l2_6.x[1];
	d__5 = std::abs(l2_6.x[6]);
/* Computing 2nd power */
	d__6 = l2_6.x[0];
/* Computing 2nd power */
	d__7 = l2_6.x[4];
/* Computing 2nd power */
	d__8 = l2_6.x[0];
/* Computing 2nd power */
	d__9 = l2_6.x[1];
	d__10 = std::abs(l2_6.x[4]);
/* Computing 2nd power */
	d__11 = l2_6.x[5];
	l3_5.g[4] = l2_6.x[0] * 10. * (d__1 * d__1) * pow_dd(&d__2, &a[m - 1])
		 / (l2_6.x[1] * (d__3 * (d__3 * d__3))) + l2_6.x[2] * 15. * 
		l2_6.x[3] / (l2_6.x[0] * (d__4 * d__4) * l2_6.x[4] * pow_dd(&
		d__5, &c_b590)) + l2_6.x[1] * 20. * l2_6.x[5] / (d__6 * d__6 *
		 l2_6.x[3] * (d__7 * d__7)) + d__8 * d__8 * 25. * (d__9 * 
		d__9) * pow_dd(&d__10, &c_b590) * l2_6.x[6] / (l2_6.x[2] * (
		d__11 * d__11)) - 100.;
    }
    if (l9_6.index1[5]) {
/* Computing 2nd power */
	d__1 = l2_6.x[3];
	d__2 = std::abs(l2_6.x[6]);
/* Computing 3rd power */
	d__3 = l2_6.x[5];
/* Computing 2nd power */
	d__4 = l2_6.x[1];
	d__5 = std::abs(l2_6.x[6]);
/* Computing 2nd power */
	d__6 = l2_6.x[0];
/* Computing 2nd power */
	d__7 = l2_6.x[4];
/* Computing 2nd power */
	d__8 = l2_6.x[0];
/* Computing 2nd power */
	d__9 = l2_6.x[1];
	d__10 = std::abs(l2_6.x[4]);
/* Computing 2nd power */
	d__11 = l2_6.x[5];
	l3_5.g[5] = -(l2_6.x[0] * 10. * (d__1 * d__1) * pow_dd(&d__2, &a[m - 
		1]) / (l2_6.x[1] * (d__3 * (d__3 * d__3))) + l2_6.x[2] * 15. *
		 l2_6.x[3] / (l2_6.x[0] * (d__4 * d__4) * l2_6.x[4] * pow_dd(&
		d__5, &c_b590)) + l2_6.x[1] * 20. * l2_6.x[5] / (d__6 * d__6 *
		 l2_6.x[3] * (d__7 * d__7)) + d__8 * d__8 * 25. * (d__9 * 
		d__9) * pow_dd(&d__10, &c_b590) * l2_6.x[6] / (l2_6.x[2] * (
		d__11 * d__11))) + 3e3;
    }
    return 0;
labelL5:
    if (! l10_6.index2[0]) {
	goto L7;
    }
    d__1 = std::abs(l2_6.x[0]);
    v1 = pow_dd(&d__1, &c_b590);
/* Computing 3rd power */
    d__1 = l2_6.x[0];
    v2 = d__1 * (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_6.x[2];
    v4 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = v4;
    v5 = d__1 * d__1;
    d__1 = std::abs(l2_6.x[3]);
    v6 = pow_dd(&d__1, &c_b590);
/* Computing 2nd power */
    d__1 = l2_6.x[5];
    v7 = d__1 * d__1;
    d__1 = std::abs(l2_6.x[5]);
    v8 = pow_dd(&d__1, &c_b933);
    d__1 = std::abs(l2_6.x[6]);
    v9 = pow_dd(&d__1, &c_b590);
    d__1 = std::abs(l2_6.x[6]);
    v10 = pow_dd(&d__1, &c_b934);
/* Computing 2nd power */
    d__1 = l2_6.x[0];
    l5_12.gg[0] = l2_6.x[6] * -.25 / (v1 * l2_6.x[2] * v7) - d__1 * d__1 * 
	    2.1 * l2_6.x[1] * l2_6.x[5] * v9 / v4;
/* Computing 2nd power */
    d__1 = l2_6.x[1];
    l5_12.gg[6] = v2 * -.7 * l2_6.x[5] * v9 / v4 + l2_6.x[2] * .2 * v8 * v10 /
	     (d__1 * d__1 * v6);
    l5_12.gg[12] = v1 * .5 * l2_6.x[6] / (v4 * v7) + v2 * 1.4 * l2_6.x[1] * 
	    l2_6.x[5] * v9 / (l2_6.x[2] * v4) - v8 * .2 * v10 / (l2_6.x[1] * 
	    v6);
    l5_12.gg[18] = l2_6.x[2] * .1 * v8 * v10 / (l2_6.x[1] * l2_6.x[3] * v6);
    d__1 = std::abs(l2_6.x[5]);
    l5_12.gg[30] = v1 * l2_6.x[6] / (l2_6.x[2] * v7 * l2_6.x[5]) - v2 * .7 * 
	    l2_6.x[1] * v9 / v4 - l2_6.x[2] * .13333333333333333 * v10 / (
	    l2_6.x[1] * v6 * pow_dd(&d__1, &c_b74));
    l5_12.gg[36] = v1 * -.5 / (l2_6.x[2] * v7) - v2 * .35 * l2_6.x[1] * 
	    l2_6.x[5] / (v4 * v9) - l2_6.x[2] * .05 * v8 / (l2_6.x[1] * v6 * 
	    v9 * v10);
L7:
    if (! l10_6.index2[1]) {
	goto L8;
    }
    d__1 = std::abs(l2_6.x[0]);
    v11 = pow_dd(&d__1, &c_b590);
    d__1 = std::abs(l2_6.x[1]);
    v12 = pow_dd(&d__1, &c_b590);
/* Computing 2nd power */
    d__1 = l2_6.x[3];
    v13 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = l2_6.x[4];
    v14 = d__1 * d__1;
    d__1 = std::abs(l2_6.x[5]);
    v15 = pow_dd(&d__1, &c_b74);
/* Computing 2nd power */
    d__1 = l2_6.x[5];
    v16 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = l2_6.x[0];
    l5_12.gg[1] = l2_6.x[1] * .65 * l2_6.x[5] / (l2_6.x[0] * v11 * l2_6.x[2] *
	     l2_6.x[4]) + v12 * 3.1 * v15 / (d__1 * d__1 * v13 * l2_6.x[4]);
    l5_12.gg[7] = l2_6.x[5] * -1.3 / (v11 * l2_6.x[2] * l2_6.x[4]) - v15 * 
	    1.55 / (l2_6.x[0] * v12 * v13 * l2_6.x[4]);
/* Computing 2nd power */
    d__1 = l2_6.x[2];
    l5_12.gg[13] = l2_6.x[1] * 1.3 * l2_6.x[5] / (v11 * (d__1 * d__1) * 
	    l2_6.x[4]) - v16 * .8 / (l2_6.x[3] * l2_6.x[4]);
    l5_12.gg[19] = l2_6.x[2] * .8 * v16 / (v13 * l2_6.x[4]) + v12 * 6.2 * v15 
	    / (l2_6.x[0] * v13 * l2_6.x[3] * l2_6.x[4]);
    l5_12.gg[25] = l2_6.x[1] * 1.3 * l2_6.x[5] / (v11 * l2_6.x[2] * v14) + 
	    l2_6.x[2] * .8 * v16 / (l2_6.x[3] * v14) + v12 * 3.1 * v15 / (
	    l2_6.x[0] * v13 * v14);
/* Computing 2nd power */
    d__1 = v15;
    l5_12.gg[31] = l2_6.x[1] * -1.3 / (v11 * l2_6.x[2] * l2_6.x[4]) - l2_6.x[
	    2] * 1.6 * l2_6.x[5] / (l2_6.x[3] * l2_6.x[4]) - v12 * 
	    1.0333333333333334 / (l2_6.x[0] * v13 * l2_6.x[4] * (d__1 * d__1))
	    ;
L8:
    if (! l10_6.index2[2]) {
	goto labelL9;
    }
/* Computing 2nd power */
    d__1 = l2_6.x[1];
    v17 = d__1 * d__1;
    d__1 = std::abs(l2_6.x[2]);
    v18 = pow_dd(&d__1, &c_b590);
    v19 = v18 * l2_6.x[2];
/* Computing 2nd power */
    d__1 = l2_6.x[5];
    v20 = d__1 * d__1;
    d__1 = std::abs(l2_6.x[6]);
    v21 = pow_dd(&d__1, &c_b74);
    d__1 = std::abs(l2_6.x[6]);
    v22 = pow_dd(&d__1, &c_b590);
/* Computing 2nd power */
    d__1 = l2_6.x[0];
    l5_12.gg[2] = l2_6.x[4] * -2. * v21 / (v19 * l2_6.x[5]) + l2_6.x[1] * v18 
	    * l2_6.x[4] / (d__1 * d__1);
    l5_12.gg[8] = -v18 * l2_6.x[4] / l2_6.x[0] + l2_6.x[2] * 1.3 * l2_6.x[4] *
	     l2_6.x[6] / (v17 * l2_6.x[1] * l2_6.x[5]) - l2_6.x[4] * .1 / (
	    v18 * v22 * l2_6.x[5]);
    l5_12.gg[14] = l2_6.x[0] * 3. * l2_6.x[4] * v21 / (l2_6.x[2] * v19 * 
	    l2_6.x[5]) + l2_6.x[1] * .05 * l2_6.x[4] / (v19 * l2_6.x[5] * v22)
	     - l2_6.x[1] * .5 * l2_6.x[4] / (l2_6.x[0] * v18) - l2_6.x[4] * 
	    .65 * l2_6.x[6] / (v17 * l2_6.x[5]);
    l5_12.gg[26] = l2_6.x[0] * -2. * v21 / (v19 * l2_6.x[5]) - l2_6.x[1] * .1 
	    / (v18 * l2_6.x[5] * v22) - l2_6.x[1] * v18 / l2_6.x[0] - l2_6.x[
	    2] * .65 * l2_6.x[6] / (v17 * l2_6.x[5]);
    l5_12.gg[32] = l2_6.x[0] * 2. * l2_6.x[4] * v21 / (v19 * v20) + l2_6.x[1] 
	    * .1 * l2_6.x[4] / (v18 * v20 * v22) + l2_6.x[2] * .65 * l2_6.x[4]
	     * l2_6.x[6] / (v17 * v20);
/* Computing 2nd power */
    d__1 = v21;
    l5_12.gg[38] = l2_6.x[0] * -.66666666666666663 * l2_6.x[4] / (v19 * 
	    l2_6.x[5] * (d__1 * d__1)) + l2_6.x[1] * .05 * l2_6.x[4] / (v18 * 
	    l2_6.x[5] * v22 * l2_6.x[6]) - l2_6.x[2] * .65 * l2_6.x[4] / (v17 
	    * l2_6.x[5]);
labelL9:
    if (! l10_6.index2[3]) {
	goto labelL10;
    }
    d__1 = std::abs(l2_6.x[0]);
    v23 = pow_dd(&d__1, &c_b590);
/* Computing 2nd power */
    d__1 = l2_6.x[0];
    v24 = d__1 * d__1;
    v25 = v24 * l2_6.x[0];
/* Computing 2nd power */
    d__1 = l2_6.x[1];
    v26 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = l2_6.x[2];
    v27 = d__1 * d__1;
    d__1 = std::abs(l2_6.x[3]);
    v28 = pow_dd(&d__1, &c_b74);
    d__1 = std::abs(l2_6.x[4]);
    v29 = pow_dd(&d__1, &c_b933);
    d__1 = std::abs(l2_6.x[4]);
    v30 = pow_dd(&d__1, &c_b590);
    d__1 = std::abs(l2_6.x[6]);
    v31 = pow_dd(&d__1, &c_b934);
/* Computing 2nd power */
    d__1 = v31;
    v32 = d__1 * d__1;
    v33 = v31 * v32;
    d__1 = std::abs(l2_6.x[6]);
    v34 = pow_dd(&d__1, &c_b74);
/* Computing 2nd power */
    d__1 = v24;
    l5_12.gg[3] = l2_6.x[1] * .4 * v30 * v34 / (v25 * l2_6.x[3]) - v26 * .15 *
	     l2_6.x[2] * v28 * v31 / (v23 * v29) + l2_6.x[2] * 1.2 * l2_6.x[4]
	     * v33 / (d__1 * d__1 * v26);
    l5_12.gg[9] = v30 * -.2 * v34 / (v24 * l2_6.x[3]) - v23 * .6 * l2_6.x[1] *
	     l2_6.x[2] * v28 * v31 / v29 + l2_6.x[2] * .8 * l2_6.x[4] * v33 / 
	    (v25 * v26 * l2_6.x[1]);
    l5_12.gg[15] = v23 * -.3 * v26 * v28 * v31 / v29 - l2_6.x[4] * .4 * v33 / 
	    (v25 * v26) + l2_6.x[3] * v32 / (v27 * l2_6.x[2]);
/* Computing 2nd power */
    d__1 = l2_6.x[3];
/* Computing 2nd power */
    d__2 = v28;
    l5_12.gg[21] = l2_6.x[1] * .2 * v30 * v34 / (v24 * (d__1 * d__1)) - v23 * 
	    .1 * v26 * l2_6.x[2] * v31 / (d__2 * d__2 * v29) - v32 * .5 / v27;
    l5_12.gg[27] = l2_6.x[1] * -.1 * v34 / (v24 * l2_6.x[3] * v30) + v23 * .2 
	    * v26 * l2_6.x[2] * v28 * v31 / (l2_6.x[4] * v29) - l2_6.x[2] * 
	    .4 * v33 / (v25 * v26);
/* Computing 2nd power */
    d__1 = v34;
    l5_12.gg[39] = l2_6.x[1] * -.066666666666666666 * v30 / (v24 * l2_6.x[3] *
	     (d__1 * d__1)) - v23 * .075 * v26 * l2_6.x[2] * v28 / (v29 * v33)
	     - l2_6.x[2] * .3 * l2_6.x[4] / (v25 * v26 * v31) - l2_6.x[3] * 
	    .25 / (v27 * v32);
labelL10:
    if (! l10_6.index2[4] && ! l10_6.index2[5]) {
	goto labelL12;
    }
/* Computing 2nd power */
    d__1 = l2_6.x[3];
    v35 = d__1 * d__1 * 10.;
    d__1 = std::abs(l2_6.x[6]);
    v36 = pow_dd(&d__1, &a[m - 1]);
/* Computing 3rd power */
    d__1 = l2_6.x[5];
    v37 = l2_6.x[1] * (d__1 * (d__1 * d__1));
    v38 = l2_6.x[2] * 15. * l2_6.x[3];
/* Computing 2nd power */
    d__1 = l2_6.x[1];
    d__2 = std::abs(l2_6.x[6]);
    v39 = l2_6.x[0] * (d__1 * d__1) * l2_6.x[4] * pow_dd(&d__2, &c_b590);
    v40 = l2_6.x[1] * 20. * l2_6.x[5];
/* Computing 2nd power */
    d__1 = l2_6.x[0];
/* Computing 2nd power */
    d__2 = l2_6.x[4];
    v41 = d__1 * d__1 * l2_6.x[3] * (d__2 * d__2);
    d__1 = std::abs(l2_6.x[4]);
    v42 = l2_6.x[0] * 25. * l2_6.x[1] * pow_dd(&d__1, &c_b590) * l2_6.x[6];
/* Computing 2nd power */
    d__1 = l2_6.x[5];
    v43 = l2_6.x[2] * (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_6.x[0];
/* Computing 2nd power */
    d__2 = l2_6.x[1];
    v44 = d__1 * d__1 * 12.5 * (d__2 * d__2) * l2_6.x[6];
    d__1 = std::abs(l2_6.x[4]);
    v45 = pow_dd(&d__1, &c_b590);
    gv[0] = v35 * v36 / v37 - v38 / (l2_6.x[0] * v39) - v40 * 2. / (l2_6.x[0] 
	    * v41) + l2_6.x[1] * 2. * v42 / v43;
    gv[1] = -v35 * l2_6.x[0] * v36 / (l2_6.x[1] * v37) - v38 * 2. / (l2_6.x[1]
	     * v39) + l2_6.x[5] * 20. / v41 + l2_6.x[0] * 2. * v42 / v43;
    gv[2] = l2_6.x[3] * 15. / v39 - l2_6.x[0] * l2_6.x[1] * v42 / (l2_6.x[2] *
	     v43);
    gv[3] = l2_6.x[0] * 20. * l2_6.x[3] * v36 / v37 + l2_6.x[2] * 15. / v39 - 
	    v40 / (l2_6.x[3] * v41);
    gv[4] = -v38 / (l2_6.x[4] * v39) - v40 * 2. / (l2_6.x[4] * v41) + v44 / (
	    v43 * v45);
    gv[5] = v35 * -3. * l2_6.x[0] * v36 / (l2_6.x[5] * v37) + l2_6.x[1] * 20. 
	    / v41 - v44 * 4. * v45 / (l2_6.x[5] * v43);
    d__1 = std::abs(l2_6.x[6]);
    d__2 = a[m - 1] - 1.;
    gv[6] = a[m - 1] * v35 * l2_6.x[0] * pow_dd(&d__1, &d__2) / v37 - v38 * 
	    .5 / (v39 * l2_6.x[6]) + v42 * l2_6.x[0] * l2_6.x[1] / (l2_6.x[6] 
	    * v43);
    if (! l10_6.index2[4]) {
	goto labelL11;
    }
    for (i__ = 1; i__ <= 7; ++i__) {
/* labelL20: */
	l5_12.gg[i__ * 6 - 2] = gv[i__ - 1];
    }
labelL11:
    if (! l10_6.index2[5]) {
	goto labelL12;
    }
    for (i__ = 1; i__ <= 7; ++i__) {
/* L30: */
	l5_12.gg[i__ * 6 - 1] = -gv[i__ - 1];
    }
labelL12:
    return 0;
} /* tp101_ */

/* Subroutine */ int tp101_(int *mode)
{
    return tp101_0_(0, mode);
    }

/* Subroutine */ int tp102_(int *mode)
{
    return tp101_0_(1, mode);
    }

/* Subroutine */ int tp103_(int *mode)
{
    return tp101_0_(2, mode);
    }


/* Subroutine */ int tp104_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static Real a;
    static int i__, j;
    static Real s, v1, v2, bx;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 8;
    l1_1.nili = 0;
    l1_1.ninl = 6;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 8; ++i__) {
	l11_7.lxl[i__ - 1] = true;
	l12_7.lxu[i__ - 1] = true;
	l13_7.xl[i__ - 1] = .1;
	l14_7.xu[i__ - 1] = 10.;
/* L33: */
    }
    a = .0588;
    l2_7.x[0] = 6.;
    l2_7.x[1] = 3.;
    l2_7.x[2] = .4;
    l2_7.x[3] = .2;
    l2_7.x[4] = 6.;
    l2_7.x[5] = 6.;
    l2_7.x[6] = 1.;
    l2_7.x[7] = .5;
    l4_7.gf[2] = 0.;
    l4_7.gf[3] = 0.;
    l4_7.gf[4] = 0.;
    l4_7.gf[5] = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
	for (j = 1; j <= 8; ++j) {
/* L34: */
	    l5_18.gg[i__ + j * 6 - 7] = 0.;
	}
    }
    l5_18.gg[0] = -.1;
    l5_18.gg[1] = -.1;
    l5_18.gg[7] = -.1;
    l20_2.lex = false;
    l20_2.xex[0] = 6.46511402796;
    l20_2.xex[1] = 2.23270864907;
    l20_2.xex[2] = .667397491303;
    l20_2.xex[3] = .595756422907;
    l20_2.xex[4] = 5.93267567811;
    l20_2.xex[5] = 5.52723456506;
    l20_2.xex[6] = 1.01332200907;
    l20_2.xex[7] = .400668229166;
    l20_2.fex = 3.95116343955;
    return 0;
labelL2:
    if (l2_7.x[0] / l2_7.x[6] < 0. || l2_7.x[1] / l2_7.x[7] < 0.) {
	goto labelL13;
    }
    d__1 = l2_7.x[0] / l2_7.x[6];
    d__2 = l2_7.x[1] / l2_7.x[7];
    l6_1.fx = (pow_dd(&d__1, &c_b993) + pow_dd(&d__2, &c_b993)) * .4 + 10. - 
	    l2_7.x[0] - l2_7.x[1];
    return 0;
labelL13:
    s = 0.;
    for (i__ = 1; i__ <= 8; ++i__) {
/* labelL14: */
/* Computing 2nd power */
	d__1 = l2_7.x[i__ - 1] - 5.;
	s += d__1 * d__1;
    }
    l6_1.fx = s + 1e3 + 3.9;
    return 0;
labelL3:
    if (l2_7.x[0] < 0. || l2_7.x[1] < 0. || l2_7.x[6] < 0. || l2_7.x[7] < 0.) 
	    {
	goto L15;
    }
    l4_7.gf[0] = pow_dd(l2_7.x, &c_b997) * .268 * pow_dd(&l2_7.x[6], &c_b998) 
	    - 1.;
    l4_7.gf[1] = pow_dd(&l2_7.x[1], &c_b997) * .268 * pow_dd(&l2_7.x[7], &
	    c_b998) - 1.;
    l4_7.gf[6] = pow_dd(l2_7.x, &c_b993) * -.268 * pow_dd(&l2_7.x[6], &
	    c_b1002);
    l4_7.gf[7] = pow_dd(&l2_7.x[1], &c_b993) * -.268 * pow_dd(&l2_7.x[7], &
	    c_b1002);
    return 0;
L15:
    for (i__ = 1; i__ <= 2; ++i__) {
	l4_7.gf[i__ - 1] = (l2_7.x[i__ - 1] - 5.) * 2.;
/* L16: */
	l4_7.gf[i__ + 5] = (l2_7.x[i__ + 5] - 5.) * 2.;
    }
    return 0;
labelL4:
    d__3 = (d__1 = l2_7.x[0] / l2_7.x[6], std::abs(d__1));
    d__4 = (d__2 = l2_7.x[1] / l2_7.x[7], std::abs(d__2));
    bx = (pow_dd(&d__3, &c_b993) + pow_dd(&d__4, &c_b993)) * .4 + 10. - 
	    l2_7.x[0] - l2_7.x[1];
    if (l9_6.index1[0]) {
	l3_5.g[0] = -a * l2_7.x[4] * l2_7.x[6] - l2_7.x[0] * .1 + 1.;
    }
    if (l9_6.index1[1]) {
	l3_5.g[1] = -a * l2_7.x[5] * l2_7.x[7] - l2_7.x[0] * .1 - l2_7.x[1] * 
		.1 + 1.;
    }
    if (l9_6.index1[2]) {
	d__1 = std::abs(l2_7.x[2]);
	d__2 = std::abs(l2_7.x[2]);
	l3_5.g[2] = (l2_7.x[2] * -4. - pow_dd(&d__1, &c_b1008) * 2.) / l2_7.x[
		4] - a * pow_dd(&d__2, &c_b1009) * l2_7.x[6] + 1.;
    }
    if (l9_6.index1[3]) {
	d__1 = std::abs(l2_7.x[3]);
	d__2 = std::abs(l2_7.x[3]);
	l3_5.g[3] = (l2_7.x[3] * -4. - pow_dd(&d__1, &c_b1008) * 2.) / l2_7.x[
		5] - a * pow_dd(&d__2, &c_b1009) * l2_7.x[7] + 1.;
    }
    if (l9_6.index1[4]) {
	l3_5.g[4] = bx - 1.;
    }
    if (l9_6.index1[5]) {
	l3_5.g[5] = 4.2 - bx;
    }
    return 0;
labelL5:
    if (! l10_6.index2[0]) {
	goto L7;
    }
    l5_18.gg[24] = -a * l2_7.x[6];
    l5_18.gg[36] = -a * l2_7.x[4];
L7:
    if (! l10_6.index2[1]) {
	goto L8;
    }
    l5_18.gg[31] = -a * l2_7.x[7];
    l5_18.gg[43] = -a * l2_7.x[5];
L8:
    if (! l10_6.index2[2]) {
	goto labelL9;
    }
/* Computing 2nd power */
    d__1 = l2_7.x[4];
    v1 = d__1 * d__1;
    d__1 = std::abs(l2_7.x[2]);
    d__2 = std::abs(l2_7.x[2]);
    l5_18.gg[14] = (pow_dd(&d__1, &c_b1015) * 1.42 - 4.) / l2_7.x[4] + a * 
	    1.3 * pow_dd(&d__2, &c_b1016) * l2_7.x[6];
    d__1 = std::abs(l2_7.x[2]);
    l5_18.gg[26] = (l2_7.x[2] * 4. + pow_dd(&d__1, &c_b1008) * 2.) / v1;
    d__1 = std::abs(l2_7.x[2]);
    l5_18.gg[38] = -a * pow_dd(&d__1, &c_b1009);
labelL9:
    if (! l10_6.index2[3]) {
	goto labelL10;
    }
/* Computing 2nd power */
    d__1 = l2_7.x[5];
    v2 = d__1 * d__1;
    d__1 = std::abs(l2_7.x[3]);
    d__2 = std::abs(l2_7.x[3]);
    l5_18.gg[21] = (pow_dd(&d__1, &c_b1015) * 1.42 - 4.) / l2_7.x[5] + a * 
	    1.3 * pow_dd(&d__2, &c_b1016) * l2_7.x[7];
    d__1 = std::abs(l2_7.x[3]);
    l5_18.gg[33] = (l2_7.x[3] * 4. + pow_dd(&d__1, &c_b1008) * 2.) / v2;
    d__1 = std::abs(l2_7.x[3]);
    l5_18.gg[45] = -a * pow_dd(&d__1, &c_b1009);
labelL10:
    if (! l10_6.index2[4] && ! l10_6.index2[5]) {
	goto labelL12;
    }
    d__1 = std::abs(l2_7.x[0]);
    d__2 = std::abs(l2_7.x[6]);
    l4_7.gf[0] = pow_dd(&d__1, &c_b997) * .268 * pow_dd(&d__2, &c_b998) - 1.;
    d__1 = std::abs(l2_7.x[1]);
    d__2 = std::abs(l2_7.x[7]);
    l4_7.gf[1] = pow_dd(&d__1, &c_b997) * .268 * pow_dd(&d__2, &c_b998) - 1.;
    d__1 = std::abs(l2_7.x[0]);
    d__2 = std::abs(l2_7.x[6]);
    l4_7.gf[6] = pow_dd(&d__1, &c_b993) * -.268 * pow_dd(&d__2, &c_b1002);
    d__1 = std::abs(l2_7.x[1]);
    d__2 = std::abs(l2_7.x[7]);
    l4_7.gf[7] = pow_dd(&d__1, &c_b993) * -.268 * pow_dd(&d__2, &c_b1002);
    if (! l10_6.index2[4]) {
	goto labelL11;
    }
    for (i__ = 1; i__ <= 8; ++i__) {
/* L31: */
	l5_18.gg[i__ * 6 - 2] = l4_7.gf[i__ - 1];
    }
labelL11:
    if (! l10_6.index2[5]) {
	goto labelL12;
    }
    for (i__ = 1; i__ <= 8; ++i__) {
/* L32: */
	l5_18.gg[i__ * 6 - 1] = -l4_7.gf[i__ - 1];
    }
labelL12:
    return 0;
} /* tp104_ */


/* Subroutine */ int tp105_(int *mode)
{
    /* System generated locals */
    int i__1;
    Real d__1;

    /* Local variables */
    static Real a[235], b[235], c__[235];
    static int i__, j;
    static Real s, v, y[235], t1, v0, v1, v2, v3, v4, v5, v6, v7, v8, 
	    v9, da[1880]	/* was [235][8] */, db[1880]	/* was [235][
	    8] */, dc[1880]	/* was [235][8] */, v10, v11;
    static Real sum;

    if (*mode - 1 <= 0) {
	goto labelL1;
    } else {
	goto labelL20;
    }
labelL1:
    l1_1.n = 8;
    l1_1.nili = 1;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 8; ++i__) {
	l11_7.lxl[i__ - 1] = true;
/* L59: */
	l12_7.lxu[i__ - 1] = true;
    }
    l13_7.xl[0] = .001;
    l13_7.xl[1] = .001;
    l13_7.xl[2] = 100.;
    l13_7.xl[3] = 130.;
    l13_7.xl[4] = 170.;
    l14_7.xu[0] = .499;
    l14_7.xu[1] = .499;
    l14_7.xu[2] = 180.;
    l14_7.xu[3] = 210.;
    l14_7.xu[4] = 240.;
    for (i__ = 6; i__ <= 8; ++i__) {
	l13_7.xl[i__ - 1] = 5.;
/* L62: */
	l14_7.xu[i__ - 1] = 25.;
    }
    l2_7.x[0] = .1;
    l2_7.x[1] = .2;
    l2_7.x[2] = 100.;
    l2_7.x[3] = 125.;
    l2_7.x[4] = 175.;
    l2_7.x[5] = 11.2;
    l2_7.x[6] = 13.2;
    l2_7.x[7] = 15.8;
    y[0] = 95.;
    y[1] = 105.;
    for (i__ = 3; i__ <= 6; ++i__) {
/* L30: */
	y[i__ - 1] = 110.;
    }
    for (i__ = 7; i__ <= 10; ++i__) {
/* L31: */
	y[i__ - 1] = 115.;
    }
    for (i__ = 11; i__ <= 25; ++i__) {
/* L32: */
	y[i__ - 1] = 120.;
    }
    for (i__ = 26; i__ <= 40; ++i__) {
/* L33: */
	y[i__ - 1] = 125.;
    }
    for (i__ = 41; i__ <= 55; ++i__) {
/* L34: */
	y[i__ - 1] = 130.;
    }
    for (i__ = 56; i__ <= 68; ++i__) {
/* L35: */
	y[i__ - 1] = 135.;
    }
    for (i__ = 69; i__ <= 89; ++i__) {
/* L36: */
	y[i__ - 1] = 140.;
    }
    for (i__ = 90; i__ <= 101; ++i__) {
/* L37: */
	y[i__ - 1] = 145.;
    }
    for (i__ = 102; i__ <= 118; ++i__) {
/* L38: */
	y[i__ - 1] = 150.;
    }
    for (i__ = 119; i__ <= 122; ++i__) {
/* L39: */
	y[i__ - 1] = 155.;
    }
    for (i__ = 123; i__ <= 142; ++i__) {
/* L40: */
	y[i__ - 1] = 160.;
    }
    for (i__ = 143; i__ <= 150; ++i__) {
/* L41: */
	y[i__ - 1] = 165.;
    }
    for (i__ = 151; i__ <= 167; ++i__) {
/* L42: */
	y[i__ - 1] = 170.;
    }
    for (i__ = 168; i__ <= 175; ++i__) {
/* L43: */
	y[i__ - 1] = 175.;
    }
    for (i__ = 176; i__ <= 181; ++i__) {
/* L44: */
	y[i__ - 1] = 180.;
    }
    for (i__ = 182; i__ <= 187; ++i__) {
/* L45: */
	y[i__ - 1] = 185.;
    }
    for (i__ = 188; i__ <= 194; ++i__) {
/* L46: */
	y[i__ - 1] = 190.;
    }
    for (i__ = 195; i__ <= 198; ++i__) {
/* L47: */
	y[i__ - 1] = 195.;
    }
    for (i__ = 199; i__ <= 201; ++i__) {
/* L48: */
	y[i__ - 1] = 200.;
    }
    for (i__ = 202; i__ <= 204; ++i__) {
/* L49: */
	y[i__ - 1] = 205.;
    }
    for (i__ = 205; i__ <= 212; ++i__) {
/* L50: */
	y[i__ - 1] = 210.;
    }
    y[212] = 215.;
    for (i__ = 214; i__ <= 219; ++i__) {
/* L51: */
	y[i__ - 1] = 220.;
    }
    for (i__ = 220; i__ <= 224; ++i__) {
/* L52: */
	y[i__ - 1] = 230.;
    }
    y[224] = 235.;
    for (i__ = 226; i__ <= 232; ++i__) {
/* L53: */
	y[i__ - 1] = 240.;
    }
    y[232] = 245.;
    y[233] = 260.;
    y[234] = 260.;
    l5_6.gg[0] = -1.;
    l5_6.gg[1] = -1.;
    for (i__ = 3; i__ <= 8; ++i__) {
/* L58: */
	l5_6.gg[i__ - 1] = 0.;
    }
    l20_2.lex = false;
    l20_2.xex[0] = .412892753597;
    l20_2.xex[1] = .403352658261;
    l20_2.xex[2] = 131.261311486;
    l20_2.xex[3] = 164.313514476;
    l20_2.xex[4] = 217.422221771;
    l20_2.xex[5] = 12.2801780396;
    l20_2.xex[6] = 15.7716989473;
    l20_2.xex[7] = 20.7468249193;
    l20_2.fex = 1138.4162396;
    return 0;
labelL20:
    if ((i__1 = *mode - 4) < 0) {
	goto L21;
    } else if (i__1 == 0) {
	goto labelL4;
    } else {
	goto labelL5;
    }
L21:
    s = 0.;
    v = 1. / std::sqrt(std::atan(1.) * 8.);
    v1 = l2_7.x[0] / l2_7.x[5];
    v2 = l2_7.x[1] / l2_7.x[6];
    v3 = (1. - l2_7.x[0] - l2_7.x[1]) / l2_7.x[7];
/* Computing 2nd power */
    d__1 = l2_7.x[5];
    v4 = 1. / (d__1 * d__1 * 2.);
/* Computing 2nd power */
    d__1 = l2_7.x[6];
    v5 = 1. / (d__1 * d__1 * 2.);
/* Computing 2nd power */
    d__1 = l2_7.x[7];
    v6 = 1. / (d__1 * d__1 * 2.);
    for (i__ = 1; i__ <= 235; ++i__) {
/* Computing 2nd power */
	d__1 = y[i__ - 1] - l2_7.x[2];
	a[i__ - 1] = v1 * std::exp(-(d__1 * d__1) * v4);
/* Computing 2nd power */
	d__1 = y[i__ - 1] - l2_7.x[3];
	b[i__ - 1] = v2 * std::exp(-(d__1 * d__1) * v5);
/* Computing 2nd power */
	d__1 = y[i__ - 1] - l2_7.x[4];
	c__[i__ - 1] = v3 * std::exp(-(d__1 * d__1) * v6);
	v11 = (a[i__ - 1] + b[i__ - 1] + c__[i__ - 1]) * v;
	if (v11 <= 0.) {
	    goto L70;
	}
/* L54: */
	s += std::log(v11);
    }
    if (*mode == 3) {
	goto labelL3;
    }
/* labelL2: */
    l6_1.fx = -s;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 235; ++i__) {
	for (j = 1; j <= 8; ++j) {
	    da[i__ + j * 235 - 236] = 0.;
	    db[i__ + j * 235 - 236] = 0.;
/* L60: */
	    dc[i__ + j * 235 - 236] = 0.;
	}
    }
    for (i__ = 1; i__ <= 235; ++i__) {
/* Computing 2nd power */
	d__1 = l2_7.x[5];
	v0 = d__1 * d__1;
	v2 = y[i__ - 1] - l2_7.x[2];
/* Computing 2nd power */
	d__1 = v2;
	v1 = std::exp(-(d__1 * d__1) / (v0 * 2.));
	da[i__ - 1] = v1 / l2_7.x[5];
/* Computing 3rd power */
	d__1 = l2_7.x[5];
	da[i__ + 469] = l2_7.x[0] * v2 / (d__1 * (d__1 * d__1)) * v1;
/* Computing 2nd power */
	d__1 = v2;
	da[i__ + 1174] = l2_7.x[0] / v0 * (d__1 * d__1 / v0 - 1.) * v1;
/* Computing 2nd power */
	d__1 = l2_7.x[6];
	v3 = d__1 * d__1;
	v4 = y[i__ - 1] - l2_7.x[3];
/* Computing 2nd power */
	d__1 = v4;
	v5 = std::exp(-(d__1 * d__1) / (v3 * 2.));
	db[i__ + 234] = v5 / l2_7.x[6];
/* Computing 3rd power */
	d__1 = l2_7.x[6];
	db[i__ + 704] = l2_7.x[1] * v4 / (d__1 * (d__1 * d__1)) * v5;
/* Computing 2nd power */
	d__1 = v4;
	db[i__ + 1409] = l2_7.x[1] / v3 * (d__1 * d__1 / v3 - 1.) * v5;
/* Computing 2nd power */
	d__1 = l2_7.x[7];
	v7 = d__1 * d__1;
	v9 = y[i__ - 1] - l2_7.x[4];
/* Computing 2nd power */
	d__1 = v9;
	v8 = std::exp(-(d__1 * d__1) / (v7 * 2.));
	v10 = 1. - l2_7.x[0] - l2_7.x[1];
	dc[i__ - 1] = -v8 / l2_7.x[7];
	dc[i__ + 234] = dc[i__ - 1];
/* Computing 3rd power */
	d__1 = l2_7.x[7];
	dc[i__ + 939] = v10 * v9 / (d__1 * (d__1 * d__1)) * v8;
/* Computing 2nd power */
	d__1 = v9;
	dc[i__ + 1644] = v10 / v7 * (d__1 * d__1 / v7 - 1.) * v8;
/* L55: */
    }
    for (j = 1; j <= 8; ++j) {
	t1 = 0.;
	for (i__ = 1; i__ <= 235; ++i__) {
/* L56: */
	    t1 += (da[i__ + j * 235 - 236] + db[i__ + j * 235 - 236] + dc[i__ 
		    + j * 235 - 236]) / (a[i__ - 1] + b[i__ - 1] + c__[i__ - 
		    1]);
	}
	l4_7.gf[j - 1] = -t1;
/* L57: */
    }
    return 0;
L70:
    for (i__ = 1; i__ <= 8; ++i__) {
	sum = (float)0.;
/* L71: */
/* Computing 2nd power */
	d__1 = l2_7.x[i__ - 1] - 5.;
	sum += d__1 * d__1;
    }
    l6_1.fx = sum + 2090.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = 1. - l2_7.x[0] - l2_7.x[1];
    }
labelL5:
    return 0;
} /* tp105_ */


/* Subroutine */ int tp106_(int *mode)
{
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 8;
    l1_1.nili = 3;
    l1_1.ninl = 3;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_7.x[i__ - 1] = 5e3;
	l13_7.xl[i__ - 1] = 1e3;
/* L23: */
	l14_7.xu[i__ - 1] = 1e4;
    }
    l13_7.xl[0] = 100.;
    for (i__ = 4; i__ <= 8; ++i__) {
	l13_7.xl[i__ - 1] = 10.;
/* L24: */
	l14_7.xu[i__ - 1] = 1e3;
    }
    l2_7.x[3] = 200.;
    l2_7.x[4] = 350.;
    l2_7.x[5] = 150.;
    l2_7.x[6] = 225.;
    l2_7.x[7] = 425.;
    for (i__ = 1; i__ <= 8; ++i__) {
	l11_7.lxl[i__ - 1] = true;
/* labelL6: */
	l12_7.lxu[i__ - 1] = true;
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	for (j = 1; j <= 8; ++j) {
/* L22: */
	    l5_18.gg[i__ + j * 6 - 7] = 0.;
	}
    }
    l4_7.gf[0] = 1.;
    l4_7.gf[1] = 1.;
    l4_7.gf[2] = 1.;
    for (i__ = 4; i__ <= 8; ++i__) {
/* labelL20: */
	l4_7.gf[i__ - 1] = 0.;
    }
    l5_18.gg[18] = -.0025;
    l5_18.gg[30] = -.0025;
    l5_18.gg[25] = -.0025;
    l5_18.gg[37] = -.0025;
    l5_18.gg[19] = .0025;
    l5_18.gg[26] = .01;
    l5_18.gg[44] = -.01;
    l5_18.gg[21] = -.0083333252000000017;
    l5_18.gg[28] = -.012500000000000001;
    l20_2.lex = false;
    l20_2.xex[0] = 579.316670388;
    l20_2.xex[1] = 1359.94292292;
    l20_2.xex[2] = 5110.07132983;
    l20_2.xex[3] = 182.017426603;
    l20_2.xex[4] = 295.598526195;
    l20_2.xex[5] = 217.979874037;
    l20_2.xex[6] = 286.41620107;
    l20_2.xex[7] = 395.597851351;
    l20_2.fex = 7049.33092308;
    return 0;
labelL2:
    l6_1.fx = l2_7.x[0] + l2_7.x[1] + l2_7.x[2];
labelL3:
    return 0;
labelL4:
    if (l9_6.index1[0]) {
	l3_5.g[0] = (l2_7.x[3] + l2_7.x[5]) * -.0025 + 1.;
    }
    if (l9_6.index1[1]) {
	l3_5.g[1] = (l2_7.x[4] + l2_7.x[6] - l2_7.x[3]) * -.0025 + 1.;
    }
    if (l9_6.index1[2]) {
	l3_5.g[2] = (l2_7.x[7] - l2_7.x[4]) * -.01 + 1.;
    }
    if (l9_6.index1[3]) {
	l3_5.g[3] = (l2_7.x[3] * -833.33252 - l2_7.x[0] * 100. + 83333.333 + 
		l2_7.x[0] * l2_7.x[5]) * 1e-5;
    }
    if (l9_6.index1[4]) {
	l3_5.g[4] = (l2_7.x[4] * -1250. - l2_7.x[1] * l2_7.x[3] + l2_7.x[3] * 
		1250. + l2_7.x[1] * l2_7.x[6]) * 1e-5;
    }
    if (l9_6.index1[5]) {
	l3_5.g[5] = (-1.25e6 - l2_7.x[2] * l2_7.x[4] + l2_7.x[4] * 2500. + 
		l2_7.x[2] * l2_7.x[7]) * 1e-5;
    }
    return 0;
labelL5:
    if (! l10_6.index2[3]) {
	goto labelL10;
    }
    l5_18.gg[3] = (l2_7.x[5] - 100.) * 1e-5;
    l5_18.gg[33] = l2_7.x[0] * 1e-5;
labelL10:
    if (! l10_6.index2[4]) {
	goto labelL11;
    }
    l5_18.gg[10] = (-l2_7.x[3] + l2_7.x[6]) * 1e-5;
    l5_18.gg[22] = (-l2_7.x[1] + 1250.) * 1e-5;
    l5_18.gg[40] = l2_7.x[1] * 1e-5;
labelL11:
    if (! l10_6.index2[5]) {
	goto labelL12;
    }
    l5_18.gg[17] = (-l2_7.x[4] + l2_7.x[7]) * 1e-5;
    l5_18.gg[29] = (-l2_7.x[2] + 2500.) * 1e-5;
    l5_18.gg[47] = l2_7.x[2] * 1e-5;
labelL12:
    return 0;
} /* tp106_ */


/* Subroutine */ int tp107_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real c__, d__;
    static int i__, j;
    static Real v1, v2, v3, v4, v5, v6, v7, v8, v9, y1, y2, y3, y4, y5, 
	    y6, v10, v11, v12, v13, v14, v15;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 9;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 6;
    l2_8.x[0] = .8;
    l2_8.x[1] = .8;
    l2_8.x[2] = .2;
    l2_8.x[3] = .2;
    l2_8.x[4] = 1.0454;
    l2_8.x[5] = 1.0454;
    l2_8.x[6] = 1.0454;
    l2_8.x[7] = 0.;
    l2_8.x[8] = 0.;
    v1 = .96460459183673464;
    c__ = v1 * std::sin(.25);
    d__ = v1 *std::cos(.25);
    for (i__ = 1; i__ <= 6; ++i__) {
	for (j = 1; j <= 9; ++j) {
/* labelL20: */
	    l5_19.gg[i__ + j * 6 - 7] = 0.;
	}
    }
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_8.lxl[i__ - 1] = true;
	l12_8.lxu[i__ - 1] = false;
/* L21: */
	l13_8.xl[i__ - 1] = 0.;
    }
    for (i__ = 3; i__ <= 4; ++i__) {
	l11_8.lxl[i__ - 1] = false;
/* L22: */
	l12_8.lxu[i__ - 1] = false;
    }
    for (i__ = 5; i__ <= 7; ++i__) {
	l11_8.lxl[i__ - 1] = true;
	l12_8.lxu[i__ - 1] = true;
	l13_8.xl[i__ - 1] = .90909;
/* L23: */
	l14_8.xu[i__ - 1] = 1.0909;
    }
    for (i__ = 8; i__ <= 9; ++i__) {
	l11_8.lxl[i__ - 1] = false;
/* L24: */
	l12_8.lxu[i__ - 1] = false;
    }
    for (i__ = 3; i__ <= 9; ++i__) {
/* L25: */
	l4_8.gf[i__ - 1] = 0.;
    }
    l5_19.gg[0] = -1.;
    l5_19.gg[7] = -1.;
    l5_19.gg[15] = -1.;
    l5_19.gg[22] = -1.;
    l20_9.lex = false;
    l20_9.xex[0] = .667009506909;
    l20_9.xex[1] = 1.02238816675;
    l20_9.xex[2] = .228287932605;
    l20_9.xex[3] = .184821729352;
    l20_9.xex[4] = 1.09090000001;
    l20_9.xex[5] = 1.09090000001;
    l20_9.xex[6] = 1.06903593236;
    l20_9.xex[7] = .106612642267;
    l20_9.xex[8] = -.338786658776;
    l20_9.fex = 5055.01180339;
    return 0;
labelL2:
/* Computing 3rd power */
    d__1 = l2_8.x[0];
/* Computing 3rd power */
    d__2 = l2_8.x[1];
    l6_1.fx = l2_8.x[0] * 3e3 + d__1 * (d__1 * d__1) * 1e3 + l2_8.x[1] * 2e3 
	    + d__2 * (d__2 * d__2) * 666.667;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_8.x[0];
    l4_8.gf[0] = d__1 * d__1 * 3e3 + 3e3;
/* Computing 2nd power */
    d__1 = l2_8.x[1];
    l4_8.gf[1] = d__1 * d__1 * 2000.001 + 2e3;
    return 0;
labelL4:
    y1 = std::sin(l2_8.x[7]);
    y2 =std::cos(l2_8.x[7]);
    y3 = std::sin(l2_8.x[8]);
    y4 =std::cos(l2_8.x[8]);
    y5 = std::sin(l2_8.x[7] - l2_8.x[8]);
    y6 =std::cos(l2_8.x[7] - l2_8.x[8]);
    if (*mode == 5) {
	goto labelL5;
    }
    if (l9_6.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_8.x[4];
	l3_5.g[0] = .4 - l2_8.x[0] + c__ * 2. * (d__1 * d__1) + l2_8.x[4] * 
		l2_8.x[5] * (-d__ * y1 - c__ * y2) + l2_8.x[4] * l2_8.x[6] * (
		-d__ * y3 - c__ * y4);
    }
    if (l9_6.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_8.x[5];
	l3_5.g[1] = .4 - l2_8.x[1] + c__ * 2. * (d__1 * d__1) + l2_8.x[4] * 
		l2_8.x[5] * (d__ * y1 - c__ * y2) + l2_8.x[5] * l2_8.x[6] * (
		d__ * y5 - c__ * y6);
    }
    if (l9_6.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_8.x[6];
	l3_5.g[2] = c__ * 2. * (d__1 * d__1) + .8 + l2_8.x[4] * l2_8.x[6] * (
		d__ * y3 - c__ * y4) + l2_8.x[5] * l2_8.x[6] * (-d__ * y5 - 
		c__ * y6);
    }
    if (l9_6.index1[3]) {
/* Computing 2nd power */
	d__1 = l2_8.x[4];
	l3_5.g[3] = .2 - l2_8.x[2] + d__ * 2. * (d__1 * d__1) - l2_8.x[4] * 
		l2_8.x[5] * (-c__ * y1 + d__ * y2) - l2_8.x[4] * l2_8.x[6] * (
		-c__ * y3 + d__ * y4);
    }
    if (l9_6.index1[4]) {
/* Computing 2nd power */
	d__1 = l2_8.x[5];
	l3_5.g[4] = .2 - l2_8.x[3] + d__ * 2. * (d__1 * d__1) - l2_8.x[4] * 
		l2_8.x[5] * (c__ * y1 + d__ * y2) - l2_8.x[5] * l2_8.x[6] * (
		c__ * y5 + d__ * y6);
    }
    if (l9_6.index1[5]) {
/* Computing 2nd power */
	d__1 = l2_8.x[6];
	l3_5.g[5] = d__ * 2. * (d__1 * d__1) - .337 - l2_8.x[4] * l2_8.x[6] * 
		(c__ * y3 + d__ * y4) - l2_8.x[5] * l2_8.x[6] * (-c__ * y5 + 
		d__ * y6);
    }
    return 0;
labelL5:
    if (! l10_6.index2[0]) {
	goto L7;
    }
    v1 = -d__ * y1 - c__ * y2;
    v2 = -d__ * y3 - c__ * y4;
    l5_19.gg[24] = c__ * 4. * l2_8.x[4] + l2_8.x[5] * v1 + l2_8.x[6] * v2;
    l5_19.gg[30] = l2_8.x[4] * v1;
    l5_19.gg[36] = l2_8.x[4] * v2;
    l5_19.gg[42] = l2_8.x[4] * l2_8.x[5] * (-d__ * y2 + c__ * y1);
    l5_19.gg[48] = l2_8.x[4] * l2_8.x[6] * (-d__ * y4 + c__ * y3);
L7:
    if (! l10_6.index2[1]) {
	goto L8;
    }
    v2 = d__ * y1 - c__ * y2;
    v3 = d__ * y6 + c__ * y5;
    v4 = d__ * y5 - c__ * y6;
    l5_19.gg[25] = l2_8.x[5] * v2;
    l5_19.gg[31] = c__ * 4. * l2_8.x[5] + l2_8.x[4] * v2 + l2_8.x[6] * v4;
    l5_19.gg[37] = l2_8.x[5] * v4;
    l5_19.gg[43] = l2_8.x[4] * l2_8.x[5] * (d__ * y2 + c__ * y1) + l2_8.x[5] *
	     l2_8.x[6] * v3;
    l5_19.gg[49] = -l2_8.x[5] * l2_8.x[6] * v3;
L8:
    if (! l10_6.index2[2]) {
	goto labelL9;
    }
    v5 = d__ * y3 - c__ * y4;
    v6 = -d__ * y5 - c__ * y6;
    v7 = -d__ * y6 + c__ * y5;
    l5_19.gg[26] = l2_8.x[6] * v5;
    l5_19.gg[32] = l2_8.x[6] * v6;
    l5_19.gg[38] = c__ * 4. * l2_8.x[6] + l2_8.x[4] * v5 + l2_8.x[5] * v6;
    l5_19.gg[44] = l2_8.x[5] * l2_8.x[6] * v7;
    l5_19.gg[50] = l2_8.x[4] * l2_8.x[6] * (d__ * y4 + c__ * y3) - l2_8.x[5] *
	     l2_8.x[6] * v7;
labelL9:
    if (! l10_6.index2[3]) {
	goto labelL10;
    }
    v8 = -c__ * y1 + d__ * y2;
    v9 = -c__ * y3 + d__ * y4;
    l5_19.gg[27] = d__ * 4. * l2_8.x[4] - l2_8.x[5] * v8 - l2_8.x[6] * v9;
    l5_19.gg[33] = -l2_8.x[4] * v8;
    l5_19.gg[39] = -l2_8.x[4] * v9;
    l5_19.gg[45] = l2_8.x[4] * l2_8.x[5] * (c__ * y2 + d__ * y1);
    l5_19.gg[51] = l2_8.x[4] * l2_8.x[6] * (c__ * y4 + d__ * y3);
labelL10:
    if (! l10_6.index2[4]) {
	goto labelL11;
    }
    v10 = c__ * y1 + d__ * y2;
    v11 = c__ * y5 + d__ * y6;
    v12 = (c__ * y6 - d__ * y5) * l2_8.x[5];
    l5_19.gg[28] = -l2_8.x[5] * v10;
    l5_19.gg[34] = d__ * 4. * l2_8.x[5] - l2_8.x[4] * v10 - l2_8.x[6] * v11;
    l5_19.gg[40] = -l2_8.x[5] * v11;
    l5_19.gg[46] = -l2_8.x[4] * l2_8.x[5] * (c__ * y2 - d__ * y1) - l2_8.x[6] 
	    * v12;
    l5_19.gg[52] = l2_8.x[6] * v12;
labelL11:
    if (! l10_6.index2[5]) {
	goto labelL12;
    }
    v13 = c__ * y3 + d__ * y4;
    v14 = -c__ * y5 + d__ * y6;
    v15 = (c__ * y6 + d__ * y5) * l2_8.x[5] * l2_8.x[6];
    l5_19.gg[29] = -l2_8.x[6] * v13;
    l5_19.gg[35] = -l2_8.x[6] * v14;
    l5_19.gg[41] = d__ * 4. * l2_8.x[6] - l2_8.x[4] * v13 - l2_8.x[5] * v14;
    l5_19.gg[47] = v15;
    l5_19.gg[53] = -l2_8.x[4] * l2_8.x[6] * (c__ * y4 - d__ * y3) - v15;
labelL12:
    return 0;
} /* tp107_ */


/* Subroutine */ int tp108_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 9;
    l1_1.nili = 0;
    l1_1.ninl = 13;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (j = 1; j <= 9; ++j) {
/* L21: */
	l2_8.x[j - 1] = 1.;
    }
    for (i__ = 1; i__ <= 9; ++i__) {
/*      XL(I)=0.0D0 */
	l11_8.lxl[i__ - 1] = false;
/*      XU(I)=1.0D0 */
/* L22: */
	l12_8.lxu[i__ - 1] = false;
    }
    l11_8.lxl[8] = true;
    l13_8.xl[8] = 0.;
    l20_9.lex = false;
    l20_9.xex[0] = .884129216724;
    l20_9.xex[1] = .467242472598;
    l20_9.xex[2] = .0374207573677;
    l20_9.xex[3] = .99929959821;
    l20_9.xex[4] = .884129216724;
    l20_9.xex[5] = .467242472594;
    l20_9.xex[6] = .0374207573618;
    l20_9.xex[7] = .99929959821;
    l20_9.xex[8] = 2.61984643608e-20;
    l20_9.fex = -.866025403841;
    return 0;
labelL2:
    l6_1.fx = (l2_8.x[0] * l2_8.x[3] - l2_8.x[1] * l2_8.x[2] + l2_8.x[2] * 
	    l2_8.x[8] - l2_8.x[4] * l2_8.x[8] + l2_8.x[4] * l2_8.x[7] - 
	    l2_8.x[5] * l2_8.x[6]) * -.5;
    return 0;
labelL3:
    l4_8.gf[0] = l2_8.x[3] * -.5;
    l4_8.gf[1] = l2_8.x[2] * .5;
    l4_8.gf[2] = (l2_8.x[1] - l2_8.x[8]) * .5;
    l4_8.gf[3] = l2_8.x[0] * -.5;
    l4_8.gf[4] = (l2_8.x[8] - l2_8.x[7]) * .5;
    l4_8.gf[5] = l2_8.x[6] * .5;
    l4_8.gf[6] = l2_8.x[5] * .5;
    l4_8.gf[7] = l2_8.x[4] * -.5;
    l4_8.gf[8] = -l4_8.gf[7] - l4_8.gf[1];
    return 0;
labelL4:
    if (l9_11.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_8.x[2];
/* Computing 2nd power */
	d__2 = l2_8.x[3];
	l3_10.g[0] = 1. - d__1 * d__1 - d__2 * d__2;
    }
    if (l9_11.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_8.x[8];
	l3_10.g[1] = 1. - d__1 * d__1;
    }
    if (l9_11.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_8.x[4];
/* Computing 2nd power */
	d__2 = l2_8.x[5];
	l3_10.g[2] = 1. - d__1 * d__1 - d__2 * d__2;
    }
    if (l9_11.index1[3]) {
/* Computing 2nd power */
	d__1 = l2_8.x[0];
/* Computing 2nd power */
	d__2 = l2_8.x[1] - l2_8.x[8];
	l3_10.g[3] = 1. - d__1 * d__1 - d__2 * d__2;
    }
    if (l9_11.index1[4]) {
/* Computing 2nd power */
	d__1 = l2_8.x[0] - l2_8.x[4];
/* Computing 2nd power */
	d__2 = l2_8.x[1] - l2_8.x[5];
	l3_10.g[4] = 1. - d__1 * d__1 - d__2 * d__2;
    }
    if (l9_11.index1[5]) {
/* Computing 2nd power */
	d__1 = l2_8.x[0] - l2_8.x[6];
/* Computing 2nd power */
	d__2 = l2_8.x[1] - l2_8.x[7];
	l3_10.g[5] = 1. - d__1 * d__1 - d__2 * d__2;
    }
    if (l9_11.index1[6]) {
/* Computing 2nd power */
	d__1 = l2_8.x[2] - l2_8.x[4];
/* Computing 2nd power */
	d__2 = l2_8.x[3] - l2_8.x[5];
	l3_10.g[6] = 1. - d__1 * d__1 - d__2 * d__2;
    }
    if (l9_11.index1[7]) {
/* Computing 2nd power */
	d__1 = l2_8.x[2] - l2_8.x[6];
/* Computing 2nd power */
	d__2 = l2_8.x[3] - l2_8.x[7];
	l3_10.g[7] = 1. - d__1 * d__1 - d__2 * d__2;
    }
    if (l9_11.index1[8]) {
/* Computing 2nd power */
	d__1 = l2_8.x[6];
/* Computing 2nd power */
	d__2 = l2_8.x[7] - l2_8.x[8];
	l3_10.g[8] = 1. - d__1 * d__1 - d__2 * d__2;
    }
    if (l9_11.index1[9]) {
	l3_10.g[9] = l2_8.x[0] * l2_8.x[3] - l2_8.x[1] * l2_8.x[2];
    }
    if (l9_11.index1[10]) {
	l3_10.g[10] = l2_8.x[2] * l2_8.x[8];
    }
    if (l9_11.index1[11]) {
	l3_10.g[11] = -l2_8.x[4] * l2_8.x[8];
    }
    if (l9_11.index1[12]) {
	l3_10.g[12] = l2_8.x[4] * l2_8.x[7] - l2_8.x[5] * l2_8.x[6];
    }
    return 0;
labelL5:
    for (i__ = 1; i__ <= 13; ++i__) {
	for (j = 1; j <= 9; ++j) {
/* L23: */
	    l5_20.gg[i__ + j * 13 - 14] = 0.;
	}
    }
    if (! l10_11.index2[0]) {
	goto L7;
    }
    l5_20.gg[26] = l2_8.x[2] * -2.;
    l5_20.gg[39] = l2_8.x[3] * -2.;
L7:
    if (! l10_11.index2[1]) {
	goto L8;
    }
    l5_20.gg[105] = l2_8.x[8] * -2.;
L8:
    if (! l10_11.index2[2]) {
	goto labelL9;
    }
    l5_20.gg[54] = l2_8.x[4] * -2.;
    l5_20.gg[67] = l2_8.x[5] * -2.;
labelL9:
    if (! l10_11.index2[3]) {
	goto labelL10;
    }
    l5_20.gg[3] = l2_8.x[0] * -2.;
    l5_20.gg[16] = (l2_8.x[1] - l2_8.x[8]) * -2.;
    l5_20.gg[107] = -l5_20.gg[16];
labelL10:
    if (! l10_11.index2[4]) {
	goto labelL11;
    }
    l5_20.gg[4] = (l2_8.x[0] - l2_8.x[4]) * -2.;
    l5_20.gg[17] = (l2_8.x[1] - l2_8.x[5]) * -2.;
    l5_20.gg[56] = -l5_20.gg[4];
    l5_20.gg[69] = -l5_20.gg[17];
labelL11:
    if (! l10_11.index2[5]) {
	goto labelL12;
    }
    l5_20.gg[5] = (l2_8.x[0] - l2_8.x[6]) * -2.;
    l5_20.gg[18] = (l2_8.x[1] - l2_8.x[7]) * -2.;
    l5_20.gg[83] = -l5_20.gg[5];
    l5_20.gg[96] = -l5_20.gg[18];
labelL12:
    if (! l10_11.index2[6]) {
	goto labelL13;
    }
    l5_20.gg[32] = (l2_8.x[2] - l2_8.x[4]) * -2.;
    l5_20.gg[45] = (l2_8.x[3] - l2_8.x[5]) * -2.;
    l5_20.gg[58] = -l5_20.gg[32];
    l5_20.gg[71] = -l5_20.gg[45];
labelL13:
    if (! l10_11.index2[7]) {
	goto labelL14;
    }
    l5_20.gg[33] = (l2_8.x[2] - l2_8.x[6]) * -2.;
    l5_20.gg[46] = (l2_8.x[3] - l2_8.x[7]) * -2.;
    l5_20.gg[85] = -l5_20.gg[33];
    l5_20.gg[98] = -l5_20.gg[46];
labelL14:
    if (! l10_11.index2[8]) {
	goto L15;
    }
    l5_20.gg[86] = l2_8.x[6] * -2.;
    l5_20.gg[99] = (l2_8.x[7] - l2_8.x[8]) * -2.;
    l5_20.gg[112] = -l5_20.gg[99];
L15:
    if (! l10_11.index2[9]) {
	goto L16;
    }
    l5_20.gg[9] = l2_8.x[3];
    l5_20.gg[22] = -l2_8.x[2];
    l5_20.gg[35] = -l2_8.x[1];
    l5_20.gg[48] = l2_8.x[0];
L16:
    if (! l10_11.index2[10]) {
	goto L17;
    }
    l5_20.gg[36] = l2_8.x[8];
    l5_20.gg[114] = l2_8.x[2];
L17:
    if (! l10_11.index2[11]) {
	goto L18;
    }
    l5_20.gg[63] = -l2_8.x[8];
    l5_20.gg[115] = -l2_8.x[4];
L18:
    if (! l10_11.index2[12]) {
	goto L19;
    }
    l5_20.gg[64] = l2_8.x[7];
    l5_20.gg[77] = -l2_8.x[6];
    l5_20.gg[90] = -l2_8.x[5];
    l5_20.gg[103] = l2_8.x[4];
L19:
    return 0;
} /* tp108_ */


/* Subroutine */ int tp109_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real a, b, c__;
    static int i__, j;
    static Real v1, v2, v3, v4, v5, v6, v7, v8, v9, ra, v10, v11, v12, 
	    v13, v14, v15, v16, v17, v18, v19, v20, hv1;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 9;
    l1_1.nili = 2;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 6;
    a = 50.176;
    ra = 1. / a;
    b = std::sin(.25);
    c__ =std::cos(.25);
    for (i__ = 1; i__ <= 2; ++i__) {
	l13_8.xl[i__ - 1] = 0.;
	l11_8.lxl[i__ - 1] = true;
	l12_8.lxu[i__ - 1] = false;
	l13_8.xl[i__ + 1] = -.55;
	l14_8.xu[i__ + 1] = .55;
	l13_8.xl[i__ + 6] = -400.;
	l14_8.xu[i__ + 6] = 800.;
	l13_8.xl[i__ + 3] = 196.;
	l14_8.xu[i__ + 3] = 252.;
/* labelL20: */
    }
    l13_8.xl[6] = 196.;
    l14_8.xu[6] = 252.;
    for (i__ = 3; i__ <= 9; ++i__) {
	l11_8.lxl[i__ - 1] = true;
	l12_8.lxu[i__ - 1] = true;
/* L21: */
    }
    for (j = 1; j <= 9; ++j) {
	l2_8.x[j - 1] = 0.;
	for (i__ = 1; i__ <= 10; ++i__) {
	    l5_21.gg[i__ + j * 10 - 11] = 0.;
/* L22: */
	}
    }
    for (i__ = 3; i__ <= 9; ++i__) {
/* L30: */
	l4_8.gf[i__ - 1] = 0.;
    }
    l5_21.gg[20] = -1.;
    l5_21.gg[30] = 1.;
    l5_21.gg[21] = 1.;
    l5_21.gg[31] = -1.;
    l5_21.gg[4] = -1.;
    l5_21.gg[15] = -1.;
    l5_21.gg[77] = 1.;
    l5_21.gg[88] = 1.;
    l20_9.lex = false;
    l20_9.xex[0] = 674.888100445;
    l20_9.xex[1] = 1134.1703947;
    l20_9.xex[2] = .133569060261;
    l20_9.xex[3] = -.371152592466;
    l20_9.xex[4] = 252.;
    l20_9.xex[5] = 252.;
    l20_9.xex[6] = 201.464535316;
    l20_9.xex[7] = 426.660777226;
    l20_9.xex[8] = 368.494083867;
    l20_9.fex = 5362.06927538;
    return 0;
labelL2:
/* Computing 3rd power */
    d__1 = l2_8.x[0];
/* Computing 3rd power */
    d__2 = l2_8.x[1];
    l6_1.fx = l2_8.x[0] * 3. + d__1 * (d__1 * d__1) * 1e-6 + d__2 * (d__2 * 
	    d__2) * 5.22074e-7 + l2_8.x[1] * 2.;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_8.x[0];
    l4_8.gf[0] = d__1 * d__1 * 3e-6 + 3.;
/* Computing 2nd power */
    d__1 = l2_8.x[1];
    l4_8.gf[1] = d__1 * d__1 * 1.566222e-6 + 2.;
    return 0;
labelL4:
    if (l9_10.index1[0]) {
	l3_9.g[0] = l2_8.x[3] - l2_8.x[2] + .55;
    }
    if (l9_10.index1[1]) {
	l3_9.g[1] = l2_8.x[2] - l2_8.x[3] + .55;
    }
    if (l9_10.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_8.x[0];
/* Computing 2nd power */
	d__2 = l2_8.x[7];
	l3_9.g[2] = 2.25e6 - d__1 * d__1 - d__2 * d__2;
    }
    if (l9_10.index1[3]) {
/* Computing 2nd power */
	d__1 = l2_8.x[1];
/* Computing 2nd power */
	d__2 = l2_8.x[8];
	l3_9.g[3] = 2.25e6 - d__1 * d__1 - d__2 * d__2;
    }
    if (l9_10.index1[4]) {
/* Computing 2nd power */
	d__1 = l2_8.x[4];
	l3_9.g[4] = (l2_8.x[4] * l2_8.x[5] * std::sin(-l2_8.x[2] - .25) + l2_8.x[4]
		 * l2_8.x[6] * std::sin(-l2_8.x[3] - .25) + d__1 * d__1 * 2. * b) *
		 ra + 400. - l2_8.x[0];
    }
    if (l9_10.index1[5]) {
/* Computing 2nd power */
	d__1 = l2_8.x[5];
	l3_9.g[5] = (l2_8.x[4] * l2_8.x[5] * std::sin(l2_8.x[2] - .25) + l2_8.x[5] 
		* l2_8.x[6] * std::sin(l2_8.x[2] - l2_8.x[3] - .25) + d__1 * d__1 *
		 2. * b) * ra + 400. - l2_8.x[1];
    }
    if (l9_10.index1[6]) {
/* Computing 2nd power */
	d__1 = l2_8.x[6];
	l3_9.g[6] = (l2_8.x[4] * l2_8.x[6] * std::sin(l2_8.x[3] - .25) + l2_8.x[5] 
		* l2_8.x[6] * std::sin(l2_8.x[3] - l2_8.x[2] - .25) + d__1 * d__1 *
		 2. * b) * ra + 881.779;
    }
    if (l9_10.index1[7]) {
/* Computing 2nd power */
	d__1 = l2_8.x[4];
/* Computing 2nd power */
	d__2 = l2_8.x[4];
	l3_9.g[7] = l2_8.x[7] + (l2_8.x[4] * l2_8.x[5] *std::cos(-l2_8.x[2] - .25)
		 + l2_8.x[4] * l2_8.x[6] *std::cos(-l2_8.x[3] - .25) - d__1 * 
		d__1 * 2. * c__) * ra + d__2 * d__2 * 7.533e-4 - 200.;
    }
    if (l9_10.index1[8]) {
/* Computing 2nd power */
	d__1 = l2_8.x[5];
/* Computing 2nd power */
	d__2 = l2_8.x[5];
	l3_9.g[8] = l2_8.x[8] + (l2_8.x[4] * l2_8.x[5] *std::cos(l2_8.x[2] - .25) 
		+ l2_8.x[6] * l2_8.x[5] *std::cos(l2_8.x[2] - l2_8.x[3] - .25) - 
		d__1 * d__1 * 2. * c__) * ra + d__2 * d__2 * 7.533e-4 - 200.;
    }
    if (l9_10.index1[9]) {
/* Computing 2nd power */
	d__1 = l2_8.x[6];
/* Computing 2nd power */
	d__2 = l2_8.x[6];
	l3_9.g[9] = (l2_8.x[4] * l2_8.x[6] *std::cos(l2_8.x[3] - .25) + l2_8.x[5] 
		* l2_8.x[6] *std::cos(l2_8.x[3] - l2_8.x[2] - .25) - d__1 * d__1 *
		 2. * c__) * ra + d__2 * d__2 * 7.533e-4 - 22.938;
    }
    return 0;
labelL5:
/* L8: */
    if (! l10_10.index2[2]) {
	goto labelL9;
    }
    l5_21.gg[2] = l2_8.x[0] * -2.;
    l5_21.gg[72] = l2_8.x[7] * -2.;
labelL9:
    if (! l10_10.index2[3]) {
	goto labelL10;
    }
    l5_21.gg[13] = l2_8.x[1] * -2.;
    l5_21.gg[83] = l2_8.x[8] * -2.;
labelL10:
    if (! l10_10.index2[4]) {
	goto labelL11;
    }
    v1 = std::sin(-l2_8.x[2] - .25);
    v2 = std::sin(-l2_8.x[3] - .25);
    v3 = l2_8.x[4] * ra;
    l5_21.gg[24] = -l2_8.x[5] * v3 *std::cos(-l2_8.x[2] - .25);
    l5_21.gg[34] = -l2_8.x[6] * v3 *std::cos(-l2_8.x[3] - .25);
    l5_21.gg[44] = (l2_8.x[5] * v1 + l2_8.x[6] * v2 + l2_8.x[4] * 4. * b) * 
	    ra;
    l5_21.gg[54] = v3 * v1;
    l5_21.gg[64] = v3 * v2;
labelL11:
    if (! l10_10.index2[5]) {
	goto labelL12;
    }
    hv1 = l2_8.x[2] - l2_8.x[3] - .25;
    v3 =std::cos(hv1);
    v4 = std::sin(l2_8.x[2] - .25);
    v5 = l2_8.x[5] * ra;
    v6 = std::sin(hv1);
    l5_21.gg[25] = l2_8.x[4] * v5 *std::cos(l2_8.x[2] - .25) + l2_8.x[6] * v5 * 
	    v3;
    l5_21.gg[35] = -l2_8.x[6] * v5 * v3;
    l5_21.gg[45] = v5 * v4;
    l5_21.gg[55] = (l2_8.x[4] * v4 + l2_8.x[6] * v6) * ra + v5 * 4. * b;
    l5_21.gg[65] = v5 * v6;
labelL12:
    if (! l10_10.index2[6]) {
	goto labelL13;
    }
    hv1 = l2_8.x[3] - l2_8.x[2] - .25;
    v7 = l2_8.x[6] * ra;
    v8 =std::cos(hv1);
    v9 = std::sin(l2_8.x[3] - .25);
    v10 = std::sin(hv1);
    l5_21.gg[26] = -l2_8.x[5] * v7 * v8;
    l5_21.gg[36] = l2_8.x[4] * v7 *std::cos(l2_8.x[3] - .25) + l2_8.x[5] * v7 * 
	    v8;
    l5_21.gg[46] = v7 * v9;
    l5_21.gg[56] = v7 * v10;
    l5_21.gg[66] = (l2_8.x[4] * v9 + l2_8.x[5] * v10) * ra + v7 * 4. * b;
labelL13:
    if (! l10_10.index2[7]) {
	goto labelL14;
    }
    v11 = l2_8.x[4] * ra;
    v12 =std::cos(-l2_8.x[2] - .25) * ra;
    v13 =std::cos(-l2_8.x[3] - .25) * ra;
    l5_21.gg[27] = l2_8.x[5] * v11 * std::sin(-l2_8.x[2] - .25);
    l5_21.gg[37] = l2_8.x[6] * v11 * std::sin(-l2_8.x[3] - .25);
    l5_21.gg[47] = l2_8.x[5] * v12 + l2_8.x[6] * v13 - v11 * 4. * c__ + 
	    l2_8.x[4] * .0015066;
    l5_21.gg[57] = l2_8.x[4] * v12;
    l5_21.gg[67] = l2_8.x[4] * v13;
labelL14:
    if (! l10_10.index2[8]) {
	goto L15;
    }
    hv1 = l2_8.x[2] - l2_8.x[3] - .25;
    v14 = std::sin(hv1) * l2_8.x[5] * ra;
    v15 =std::cos(l2_8.x[2] - .25) * ra;
    v16 =std::cos(hv1) * ra;
    l5_21.gg[28] = -l2_8.x[4] * l2_8.x[5] * std::sin(l2_8.x[2] - .25) * ra - 
	    l2_8.x[6] * v14;
    l5_21.gg[38] = l2_8.x[6] * v14;
    l5_21.gg[48] = l2_8.x[5] * v15;
    l5_21.gg[58] = l2_8.x[4] * v15 + l2_8.x[6] * v16 - l2_8.x[5] * 4. * c__ * 
	    ra + l2_8.x[5] * .0015066;
    l5_21.gg[68] = l2_8.x[5] * v16;
L15:
    if (! l10_10.index2[9]) {
	goto L16;
    }
    hv1 = l2_8.x[3] - l2_8.x[2] - .25;
    v17 = std::sin(hv1) * l2_8.x[5] * ra;
    v18 =std::cos(l2_8.x[3] - .25) * ra;
    v19 =std::cos(hv1) * ra;
    v20 = l2_8.x[6] * ra;
    l5_21.gg[29] = l2_8.x[6] * v17;
    l5_21.gg[39] = -l2_8.x[4] * v20 * std::sin(l2_8.x[3] - .25) - l2_8.x[6] * v17;
    l5_21.gg[49] = l2_8.x[6] * v18;
    l5_21.gg[59] = l2_8.x[6] * v19;
    l5_21.gg[69] = l2_8.x[4] * v18 + l2_8.x[5] * v19 - v20 * 4. * c__ + 
	    l2_8.x[6] * .0015066;
L16:
    return 0;
} /* tp109_ */


/* Subroutine */ int tp110_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    Real d_sign(Real*, Real*);
    static int i__;
    static Real s, t, u, sum;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 10;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 10; ++i__) {
	l11_9.lxl[i__ - 1] = true;
	l12_9.lxu[i__ - 1] = true;
	l13_9.xl[i__ - 1] = 2.001;
	l14_9.xu[i__ - 1] = 9.999;
	l2_9.x[i__ - 1] = 9.;
/* labelL20: */
    }
    l20_10.lex = false;
    l20_10.xex[0] = 9.35025654733;
    l20_10.xex[1] = 9.35025654733;
    l20_10.xex[2] = 9.35025654733;
    l20_10.xex[3] = 9.35025654733;
    l20_10.xex[4] = 9.35025654733;
    l20_10.xex[5] = 9.35025654733;
    l20_10.xex[6] = 9.35025654733;
    l20_10.xex[7] = 9.35025654733;
    l20_10.xex[8] = 9.35025654733;
    l20_10.xex[9] = 9.35025654733;
    l20_10.fex = -45.7784697153;
    return 0;
labelL2:
    t = 1.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L30: */
	t *= l2_9.x[i__ - 1];
    }
    d__1 = std::abs(t);
    s = pow_dd(&d__1, &c_b1157);
    s = d_sign(&s, &t);
    if (*mode == 3) {
	goto labelL3;
    }
    u = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
	if (l2_9.x[i__ - 1] - 2. <= 0. || 10. - l2_9.x[i__ - 1] <= 0.) {
	    goto L33;
	}
/* L31: */
/* Computing 2nd power */
	d__1 = std::log(l2_9.x[i__ - 1] - 2.);
/* Computing 2nd power */
	d__2 = std::log(10. - l2_9.x[i__ - 1]);
	u = u + d__1 * d__1 + d__2 * d__2;
    }
    l6_1.fx = u - s;
    return 0;
L33:
    sum = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L34: */
/* Computing 2nd power */
	d__1 = l2_9.x[i__ - 1] - 5.;
	sum += d__1 * d__1;
    }
    l6_1.fx = sum + 1e3 - 45.8;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 10; ++i__) {
	if (l2_9.x[i__ - 1] - 2. <= 0. || 10. - l2_9.x[i__ - 1] <= 0.) {
	    goto L35;
	}
	l4_9.gf[i__ - 1] = (std::log(l2_9.x[i__ - 1] - 2.) / (l2_9.x[i__ - 1] - 2.)
		 - log(10. - l2_9.x[i__ - 1]) / (10. - l2_9.x[i__ - 1])) * 2. 
		- s / l2_9.x[i__ - 1] * .2;
	goto L32;
L35:
	l4_9.gf[i__ - 1] = (l2_9.x[i__ - 1] - 5.) * 2.;
L32:
	;
    }
labelL4:
    return 0;
} /* tp110_ */


/* Subroutine */ int tp111_(int *mode)
{
    /* Local variables */
    static Real c__[10];
    static int i__, j;
    static Real s, t;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 10;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 3;
    for (j = 1; j <= 10; ++j) {
	l2_9.x[j - 1] = -2.3;
	for (i__ = 1; i__ <= 3; ++i__) {
/* labelL20: */
	    l5_14.gg[i__ + j * 3 - 4] = 0.;
	}
    }
    c__[0] = -6.089;
    c__[1] = -17.164;
    c__[2] = -34.054;
    c__[3] = -5.914;
    c__[4] = -24.721;
    c__[5] = -14.986;
    c__[6] = -24.1;
    c__[7] = -10.708;
    c__[8] = -26.662;
    c__[9] = -22.179;
    for (i__ = 1; i__ <= 10; ++i__) {
	l13_9.xl[i__ - 1] = -100.;
	l14_9.xu[i__ - 1] = 100.;
	l11_9.lxl[i__ - 1] = true;
/* labelL6: */
	l12_9.lxu[i__ - 1] = true;
    }
    l20_10.lex = false;
    l20_10.xex[0] = -3.20121253241;
    l20_10.xex[1] = -1.91205959435;
    l20_10.xex[2] = -.244441308369;
    l20_10.xex[3] = -6.53748856532;
    l20_10.xex[4] = -.723152425984;
    l20_10.xex[5] = -7.26773826993;
    l20_10.xex[6] = -3.59671064233;
    l20_10.xex[7] = -4.01776873216;
    l20_10.xex[8] = -3.28746169619;
    l20_10.xex[9] = -2.33558183059;
    l20_10.fex = -47.7610902637;
    return 0;
labelL2:
    t = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L30: */
	t += std::exp(l2_9.x[i__ - 1]);
    }
    if (*mode == 3) {
	goto labelL3;
    }
    s = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L31: */
	s += std::exp(l2_9.x[i__ - 1]) * (c__[i__ - 1] + l2_9.x[i__ - 1] - std::log(t));
    }
    l6_1.fx = s;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 10; ++i__) {
/* L33: */
	l4_9.gf[i__ - 1] = std::exp(l2_9.x[i__ - 1]) * (c__[i__ - 1] + l2_9.x[i__ 
		- 1] - std::log(t));
    }
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = std::exp(l2_9.x[0]) + std::exp(l2_9.x[1]) * 2. + std::exp(l2_9.x[2]) * 
		2. + std::exp(l2_9.x[5]) + std::exp(l2_9.x[9]) - 2.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = std::exp(l2_9.x[3]) + std::exp(l2_9.x[4]) * 2. + std::exp(l2_9.x[5]) + 
		std::exp(l2_9.x[6]) - 1.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = std::exp(l2_9.x[2]) + std::exp(l2_9.x[6]) + std::exp(l2_9.x[7]) + std::exp(
		l2_9.x[8]) * 2. + std::exp(l2_9.x[9]) - 1.;
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto L7;
    }
    l5_14.gg[0] = std::exp(l2_9.x[0]);
    l5_14.gg[3] = std::exp(l2_9.x[1]) * 2.;
    l5_14.gg[6] = std::exp(l2_9.x[2]) * 2.;
    l5_14.gg[15] = std::exp(l2_9.x[5]);
    l5_14.gg[27] = std::exp(l2_9.x[9]);
L7:
    if (! l10_4.index2[1]) {
	goto L8;
    }
    l5_14.gg[10] = std::exp(l2_9.x[3]);
    l5_14.gg[13] = std::exp(l2_9.x[4]) * 2.;
    l5_14.gg[16] = std::exp(l2_9.x[5]);
    l5_14.gg[19] = std::exp(l2_9.x[6]);
L8:
    if (! l10_4.index2[2]) {
	goto labelL9;
    }
    l5_14.gg[8] = std::exp(l2_9.x[2]);
    l5_14.gg[20] = std::exp(l2_9.x[6]);
    l5_14.gg[23] = std::exp(l2_9.x[7]);
    l5_14.gg[26] = std::exp(l2_9.x[8]) * 2.;
    l5_14.gg[29] = std::exp(l2_9.x[9]);
labelL9:
    return 0;
} /* tp111_ */


/* Subroutine */ int tp112_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real c__[10];
    static int i__, j;
    static Real s, t, dlogt;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 10;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 3;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 10; ++i__) {
	l13_9.xl[i__ - 1] = 1e-4;
	l11_9.lxl[i__ - 1] = true;
/* labelL6: */
	l12_9.lxu[i__ - 1] = false;
    }
    c__[0] = -6.089;
    c__[1] = -17.164;
    c__[2] = -34.054;
    c__[3] = -5.914;
    c__[4] = -24.721;
    c__[5] = -14.986;
    c__[6] = -24.1;
    c__[7] = -10.708;
    c__[8] = -26.662;
    c__[9] = -22.179;
    for (j = 1; j <= 10; ++j) {
	l2_9.x[j - 1] = .1;
	for (i__ = 1; i__ <= 3; ++i__) {
/* labelL20: */
	    l5_14.gg[i__ + j * 3 - 4] = 0.;
	}
    }
    l5_14.gg[0] = 1.;
    l5_14.gg[3] = 2.;
    l5_14.gg[6] = 2.;
    l5_14.gg[15] = 1.;
    l5_14.gg[27] = 1.;
    l5_14.gg[10] = 1.;
    l5_14.gg[13] = 2.;
    l5_14.gg[16] = 1.;
    l5_14.gg[19] = 1.;
    l5_14.gg[8] = 1.;
    l5_14.gg[20] = 1.;
    l5_14.gg[23] = 1.;
    l5_14.gg[26] = 2.;
    l5_14.gg[29] = 1.;
    l20_10.lex = false;
    l20_10.xex[0] = .0177354776881;
    l20_10.xex[1] = .0820018011109;
    l20_10.xex[2] = .88256455892;
    l20_10.xex[3] = 7.23325625629e-4;
    l20_10.xex[4] = .490785079062;
    l20_10.xex[5] = 4.33546900325e-4;
    l20_10.xex[6] = .0172729773078;
    l20_10.xex[7] = .00776563912291;
    l20_10.xex[8] = .0198492864597;
    l20_10.xex[9] = .0526982611793;
    l20_10.fex = -.47761086;
    return 0;
labelL2:
    t = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L30: */
	t += l2_9.x[i__ - 1];
    }
    if (*mode == 3) {
	goto labelL3;
    }
    if (t < 1e-5) {
	goto L34;
    }
    dlogt = std::log(t);
    s = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
	if (l2_9.x[i__ - 1] < 0.) {
	    goto L34;
	}
/* L31: */
	s += l2_9.x[i__ - 1] * (c__[i__ - 1] + std::log(l2_9.x[i__ - 1]) - dlogt);
    }
    l6_1.fx = s * (float).01;
    return 0;
L34:
    s = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L35: */
	if (l2_9.x[i__ - 1] < 0.) {
/* Computing 2nd power */
	    d__1 = l2_9.x[i__ - 1] - 5.;
	    s += d__1 * d__1;
	}
    }
    l6_1.fx = (s + 1e3 - 47.8) * (float).01;
    return 0;
labelL3:
    if (t < 1e-5) {
	goto L36;
    }
    dlogt = std::log(t);
    for (i__ = 1; i__ <= 10; ++i__) {
	if (l2_9.x[i__ - 1] < 0.) {
	    goto L36;
	}
/* L33: */
	l4_9.gf[i__ - 1] = (c__[i__ - 1] + std::log(l2_9.x[i__ - 1]) - dlogt) * (
		float).01;
    }
    return 0;
L36:
    for (i__ = 1; i__ <= 10; ++i__) {
	l4_9.gf[i__ - 1] = 0.;
/* L37: */
	if (l2_9.x[i__ - 1] < 0.) {
	    l4_9.gf[i__ - 1] = (l2_9.x[i__ - 1] - 5.) * 2. * (float).01;
	}
    }
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_9.x[0] + l2_9.x[1] * 2. + l2_9.x[2] * 2. + l2_9.x[5] + 
		l2_9.x[9] - (float)2.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_9.x[3] + l2_9.x[4] * 2. + l2_9.x[5] + l2_9.x[6] - (
		float)1.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = l2_9.x[2] + l2_9.x[6] + l2_9.x[7] + l2_9.x[8] * 2. + 
		l2_9.x[9] - (float)1.;
    }
labelL5:
    return 0;
} /* tp112_ */


/* Subroutine */ int tp113_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 10;
    l1_1.nili = 3;
    l1_1.ninl = 5;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 8; ++i__) {
	for (j = 1; j <= 10; ++j) {
/* labelL20: */
	    l5_22.gg[i__ + (j << 3) - 9] = 0.;
	}
    }
    l2_9.x[0] = 2.;
    l2_9.x[1] = 3.;
    l2_9.x[2] = 5.;
    l2_9.x[3] = 5.;
    l2_9.x[4] = 1.;
    l2_9.x[5] = 2.;
    l2_9.x[6] = 7.;
    l2_9.x[7] = 3.;
    l2_9.x[8] = 6.;
    l2_9.x[9] = 10.;
    for (i__ = 1; i__ <= 10; ++i__) {
	l11_9.lxl[i__ - 1] = false;
/* labelL6: */
	l12_9.lxu[i__ - 1] = false;
    }
    l5_22.gg[0] = -4.;
    l5_22.gg[8] = -5.;
    l5_22.gg[48] = 3.;
    l5_22.gg[56] = -9.;
    l5_22.gg[1] = -10.;
    l5_22.gg[9] = 8.;
    l5_22.gg[49] = 17.;
    l5_22.gg[57] = -2.;
    l5_22.gg[2] = 8.;
    l5_22.gg[10] = -2.;
    l5_22.gg[66] = -5.;
    l5_22.gg[74] = 2.;
    l5_22.gg[27] = 7.;
    l5_22.gg[12] = -8.;
    l5_22.gg[28] = 2.;
    l5_22.gg[45] = 1.;
    l5_22.gg[38] = -14.;
    l5_22.gg[46] = 6.;
    l5_22.gg[7] = 3.;
    l5_22.gg[15] = -6.;
    l5_22.gg[79] = 7.;
    l20_10.lex = false;
    l20_10.xex[0] = 2.17199637118;
    l20_10.xex[1] = 2.36368297378;
    l20_10.xex[2] = 8.77392573847;
    l20_10.xex[3] = 5.09598448797;
    l20_10.xex[4] = .990654764992;
    l20_10.xex[5] = 1.43057397893;
    l20_10.xex[6] = 1.32164420805;
    l20_10.xex[7] = 9.82872580801;
    l20_10.xex[8] = 8.28009167017;
    l20_10.xex[9] = 8.37592666387;
    l20_10.fex = 24.3062090641;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_9.x[0];
/* Computing 2nd power */
    d__2 = l2_9.x[1];
/* Computing 2nd power */
    d__3 = l2_9.x[2] - 10.;
/* Computing 2nd power */
    d__4 = l2_9.x[3] - 5.;
/* Computing 2nd power */
    d__5 = l2_9.x[4] - 3.;
/* Computing 2nd power */
    d__6 = l2_9.x[5] - 1.;
/* Computing 2nd power */
    d__7 = l2_9.x[6];
/* Computing 2nd power */
    d__8 = l2_9.x[7] - 11.;
/* Computing 2nd power */
    d__9 = l2_9.x[8] - 10.;
/* Computing 2nd power */
    d__10 = l2_9.x[9] - 7.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + l2_9.x[0] * l2_9.x[1] - l2_9.x[0] * 
	    14. - l2_9.x[1] * 16. + d__3 * d__3 + d__4 * d__4 * 4. + d__5 * 
	    d__5 + d__6 * d__6 * 2. + d__7 * d__7 * 5. + d__8 * d__8 * 7. + 
	    d__9 * d__9 * 2. + d__10 * d__10 + 45.;
    return 0;
labelL3:
    l4_9.gf[0] = l2_9.x[0] * 2. + l2_9.x[1] - 14.;
    l4_9.gf[1] = l2_9.x[1] * 2. + l2_9.x[0] - 16.;
    l4_9.gf[2] = (l2_9.x[2] - 10.) * 2.;
    l4_9.gf[3] = (l2_9.x[3] - 5.) * 8.;
    l4_9.gf[4] = (l2_9.x[4] - 3.) * 2.;
    l4_9.gf[5] = (l2_9.x[5] - 1.) * 4.;
    l4_9.gf[6] = l2_9.x[6] * 10.;
    l4_9.gf[7] = (l2_9.x[7] - 11.) * 14.;
    l4_9.gf[8] = (l2_9.x[8] - 10.) * 4.;
    l4_9.gf[9] = (l2_9.x[9] - 7.) * 2.;
    return 0;
labelL4:
    if (l9_12.index1[0]) {
	l3_11.g[0] = l2_9.x[0] * -4. - l2_9.x[1] * 5. + l2_9.x[6] * 3. - 
		l2_9.x[7] * 9. + 105.;
    }
    if (l9_12.index1[1]) {
	l3_11.g[1] = l2_9.x[0] * -10. + l2_9.x[1] * 8. + l2_9.x[6] * 17. - 
		l2_9.x[7] * 2.;
    }
    if (l9_12.index1[2]) {
	l3_11.g[2] = l2_9.x[0] * 8. - l2_9.x[1] * 2. - l2_9.x[8] * 5. + 
		l2_9.x[9] * 2. + 12.;
    }
    if (l9_12.index1[3]) {
/* Computing 2nd power */
	d__1 = l2_9.x[0] - 2.;
/* Computing 2nd power */
	d__2 = l2_9.x[1] - 3.;
/* Computing 2nd power */
	d__3 = l2_9.x[2];
	l3_11.g[3] = d__1 * d__1 * -3. - d__2 * d__2 * 4. - d__3 * d__3 * 2. 
		+ l2_9.x[3] * 7. + 120.;
    }
    if (l9_12.index1[4]) {
/* Computing 2nd power */
	d__1 = l2_9.x[0];
/* Computing 2nd power */
	d__2 = l2_9.x[2] - 6.;
	l3_11.g[4] = d__1 * d__1 * -5. - l2_9.x[1] * 8. - d__2 * d__2 + 
		l2_9.x[3] * 2. + 40.;
    }
    if (l9_12.index1[5]) {
/* Computing 2nd power */
	d__1 = l2_9.x[0] - 8.;
/* Computing 2nd power */
	d__2 = l2_9.x[1] - 4.;
/* Computing 2nd power */
	d__3 = l2_9.x[4];
	l3_11.g[5] = d__1 * d__1 * -.5 - d__2 * d__2 * 2. - d__3 * d__3 * 3. 
		+ l2_9.x[5] + 30.;
    }
    if (l9_12.index1[6]) {
/* Computing 2nd power */
	d__1 = l2_9.x[0];
/* Computing 2nd power */
	d__2 = l2_9.x[1] - 2.;
	l3_11.g[6] = -(d__1 * d__1) - d__2 * d__2 * 2. + l2_9.x[0] * 2. * 
		l2_9.x[1] - l2_9.x[4] * 14. + l2_9.x[5] * 6.;
    }
    if (l9_12.index1[7]) {
/* Computing 2nd power */
	d__1 = l2_9.x[8] - 8.;
	l3_11.g[7] = l2_9.x[0] * 3. - l2_9.x[1] * 6. - d__1 * d__1 * 12. + 
		l2_9.x[9] * 7.;
    }
    return 0;
labelL5:
    if (! l10_12.index2[3]) {
	goto labelL10;
    }
    l5_22.gg[3] = (l2_9.x[0] - 2.) * -6.;
    l5_22.gg[11] = (l2_9.x[1] - 3.) * -8.;
    l5_22.gg[19] = l2_9.x[2] * -4.;
labelL10:
    if (! l10_12.index2[4]) {
	goto labelL11;
    }
    l5_22.gg[4] = l2_9.x[0] * -10.;
    l5_22.gg[20] = (l2_9.x[2] - 6.) * -2.;
labelL11:
    if (! l10_12.index2[5]) {
	goto labelL12;
    }
    l5_22.gg[5] = 8. - l2_9.x[0];
    l5_22.gg[13] = (l2_9.x[1] - 4.) * -4.;
    l5_22.gg[37] = l2_9.x[4] * -6.;
labelL12:
    if (! l10_12.index2[6]) {
	goto labelL13;
    }
    l5_22.gg[6] = l2_9.x[0] * -2. + l2_9.x[1] * 2.;
    l5_22.gg[14] = (l2_9.x[1] - 2.) * -4. + l2_9.x[0] * 2.;
labelL13:
    if (l10_12.index2[7]) {
	l5_22.gg[71] = (l2_9.x[8] - 8.) * -24.;
    }
    return 0;
} /* tp113_ */


/* Subroutine */ int tp114_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__, j;
    static Real v1, v2, v3;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 10;
    l1_1.nili = 4;
    l1_1.ninl = 4;
    l1_1.neli = 1;
    l1_1.nenl = 2;
    l2_9.x[0] = 1745.;
    l2_9.x[1] = 1.2e4;
    l2_9.x[2] = 110.;
    l2_9.x[3] = 3048.;
    l2_9.x[4] = 1974.;
    l2_9.x[5] = 89.2;
    l2_9.x[6] = 92.8;
    l2_9.x[7] = 8.;
    l2_9.x[8] = 3.6;
    l2_9.x[9] = 145.;
    for (i__ = 1; i__ <= 5; ++i__) {
/* L25: */
	l13_9.xl[i__ - 1] = 1e-5;
    }
    l13_9.xl[5] = 85.;
    l13_9.xl[6] = 90.;
    l13_9.xl[7] = 3.;
    l13_9.xl[8] = 1.2;
    l13_9.xl[9] = 145.;
    l14_9.xu[0] = 2e3;
    l14_9.xu[1] = 1.6e4;
    l14_9.xu[2] = 120.;
    l14_9.xu[3] = 5e3;
    l14_9.xu[4] = 2e3;
    l14_9.xu[5] = 93.;
    l14_9.xu[6] = 95.;
    l14_9.xu[7] = 12.;
    l14_9.xu[8] = 4.;
    l14_9.xu[9] = 162.;
    for (i__ = 1; i__ <= 10; ++i__) {
	l11_9.lxl[i__ - 1] = true;
	l12_9.lxu[i__ - 1] = true;
	for (j = 1; j <= 11; ++j) {
/* labelL6: */
	    l5_23.gg[j + i__ * 11 - 12] = 0.;
	}
    }
    l4_9.gf[0] = 5.04e-4;
    l4_9.gf[1] = 3.5000000000000004e-6;
    l4_9.gf[2] = .001;
    l4_9.gf[4] = 3.3599999999999998e-4;
    l4_9.gf[5] = 0.;
    l4_9.gf[7] = 0.;
    l4_9.gf[8] = 0.;
    l4_9.gf[9] = 0.;
    l5_23.gg[88] = -.9;
    l5_23.gg[99] = -.222;
    l5_23.gg[67] = 3.;
    l5_23.gg[100] = -.99;
    l5_23.gg[90] = 1.1111111111111112;
    l5_23.gg[101] = .222;
    l5_23.gg[69] = -3.;
    l5_23.gg[102] = 1.0101010101010102;
    l5_23.gg[37] = -.99;
    l5_23.gg[60] = .325;
    l5_23.gg[71] = -.99;
    l5_23.gg[39] = 1.0101010101010102;
    l5_23.gg[62] = -.325;
    l5_23.gg[73] = 1.0101010101010102;
    l5_23.gg[8] = -1.;
    l5_23.gg[41] = 1.22;
    l5_23.gg[52] = -1.;
    l5_23.gg[64] = -1.;
    l5_23.gg[87] = -1.;
    l20_10.lex = false;
    l20_10.xex[0] = 1698.09564792;
    l20_10.xex[1] = 15818.7256985;
    l20_10.xex[2] = 54.1022827849;
    l20_10.xex[3] = 3031.22594099;
    l20_10.xex[4] = 2e3;
    l20_10.xex[5] = 90.1153669236;
    l20_10.xex[6] = 95.;
    l20_10.xex[7] = 10.4933580864;
    l20_10.xex[8] = 1.5616363638;
    l20_10.xex[9] = 153.535353535;
    l20_10.fex = -17.688069634399998;
    return 0;
labelL2:
    l6_1.fx = l2_9.x[0] * 5.04 + l2_9.x[1] * .035 + l2_9.x[2] * 10. + l2_9.x[
	    4] * 3.36 - l2_9.x[3] * .063 * l2_9.x[6];
    l6_1.fx *= .01;
labelL3:
    l4_9.gf[3] = l2_9.x[6] * -.063 * .01;
    l4_9.gf[6] = l2_9.x[3] * -.063 * .01;
    return 0;
labelL4:
    if (l9_13.index1[0]) {
	l3_12.g[0] = 35.82 - l2_9.x[9] * .222 - l2_9.x[8] * .9;
    }
    if (l9_13.index1[1]) {
	l3_12.g[1] = l2_9.x[6] * 3. - 133. - l2_9.x[9] * .99;
    }
    if (l9_13.index1[2]) {
	l3_12.g[2] = l2_9.x[9] * .222 - 35.82 + l2_9.x[8] * 
		1.1111111111111112;
    }
    if (l9_13.index1[3]) {
	l3_12.g[3] = 133. - l2_9.x[6] * 3. + l2_9.x[9] / .99;
    }
    if (l9_13.index1[4]) {
/* Computing 2nd power */
	d__1 = l2_9.x[7];
	l3_12.g[4] = l2_9.x[0] * 1.12 + l2_9.x[0] * .13167 * l2_9.x[7] - 
		l2_9.x[0] * .00667 * (d__1 * d__1) - l2_9.x[3] * .99;
    }
    if (l9_13.index1[5]) {
/* Computing 2nd power */
	d__1 = l2_9.x[7];
	l3_12.g[5] = l2_9.x[7] * 1.098 + 57.425 - d__1 * d__1 * .038 + l2_9.x[
		5] * .325 - l2_9.x[6] * .99;
    }
    if (l9_13.index1[6]) {
/* Computing 2nd power */
	d__1 = l2_9.x[7];
	l3_12.g[6] = l2_9.x[0] * -1.12 - l2_9.x[0] * .13167 * l2_9.x[7] + 
		l2_9.x[0] * .00667 * (d__1 * d__1) + l2_9.x[3] / .99;
    }
    if (l9_13.index1[7]) {
/* Computing 2nd power */
	d__1 = l2_9.x[7];
	l3_12.g[7] = -57.425 - l2_9.x[7] * 1.098 + d__1 * d__1 * .038 - 
		l2_9.x[5] * .325 + l2_9.x[6] / .99;
    }
    if (l9_13.index1[8]) {
	l3_12.g[8] = l2_9.x[3] * 1.22 - l2_9.x[0] - l2_9.x[4];
    }
    if (l9_13.index1[9]) {
	l3_12.g[9] = l2_9.x[2] * 9.8e4 / (l2_9.x[3] * l2_9.x[8] + l2_9.x[2] * 
		1e3) - l2_9.x[5];
    }
    if (l9_13.index1[10]) {
	l3_12.g[10] = (l2_9.x[1] + l2_9.x[4]) / l2_9.x[0] - l2_9.x[7];
    }
    return 0;
labelL5:
    if (! l10_13.index2[4]) {
	goto labelL11;
    }
/* Computing 2nd power */
    d__1 = l2_9.x[7];
    l5_23.gg[4] = l2_9.x[7] * .13167 + 1.12 - d__1 * d__1 * .00667;
    l5_23.gg[81] = l2_9.x[0] * .13167 - l2_9.x[0] * .013339999999999999 * 
	    l2_9.x[7];
labelL11:
    if (l10_13.index2[5]) {
	l5_23.gg[82] = 1.098 - l2_9.x[7] * .076;
    }
    if (! l10_13.index2[6]) {
	goto labelL13;
    }
    l5_23.gg[83] = l2_9.x[0] * -.13167 + l2_9.x[0] * .013339999999999999 * 
	    l2_9.x[7];
/* Computing 2nd power */
    d__1 = l2_9.x[7];
    l5_23.gg[6] = l2_9.x[7] * -.13167 + d__1 * d__1 * .00667 - 1.12;
labelL13:
    if (l10_13.index2[7]) {
	l5_23.gg[84] = l2_9.x[7] * .076 - 1.098;
    }
    if (! l10_13.index2[9]) {
	goto L16;
    }
/* Computing 2nd power */
    d__1 = l2_9.x[3] * l2_9.x[8] + l2_9.x[2] * 1e3;
    v1 = d__1 * d__1;
    v2 = l2_9.x[8] * 9.8e4;
    v3 = v2 / v1;
    l5_23.gg[31] = l2_9.x[3] * v3;
    l5_23.gg[42] = -l2_9.x[2] * v3;
    l5_23.gg[97] = l2_9.x[2] * -9.8e4 * l2_9.x[3] / v1;
L16:
    if (! l10_13.index2[10]) {
	goto L17;
    }
/* Computing 2nd power */
    d__1 = l2_9.x[0];
    l5_23.gg[10] = -(l2_9.x[1] + l2_9.x[4]) / (d__1 * d__1);
    l5_23.gg[21] = 1. / l2_9.x[0];
    l5_23.gg[54] = l5_23.gg[21];
L17:
    return 0;
} /* tp114_ */


/* Subroutine */ int tp116_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 13;
    l1_1.nili = 5;
    l1_1.ninl = 10;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_10.x[0] = .5;
    l2_10.x[1] = .8;
    l2_10.x[2] = .9;
    l2_10.x[3] = .1;
    l2_10.x[4] = .14;
    l2_10.x[5] = .5;
    l2_10.x[6] = 489.;
    l2_10.x[7] = 80.;
    l2_10.x[8] = 650.;
    l2_10.x[9] = 450.;
    l2_10.x[10] = 150.;
    l2_10.x[11] = 150.;
    l2_10.x[12] = 150.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* labelL6: */
	l13_10.xl[i__ - 1] = .1;
    }
    l13_10.xl[3] = 1e-4;
    l13_10.xl[8] = 500.;
    l13_10.xl[10] = 1.;
    l13_10.xl[11] = 1e-4;
    l13_10.xl[12] = 1e-4;
    for (i__ = 1; i__ <= 3; ++i__) {
	l14_10.xu[i__ - 1] = 1.;
	l14_10.xu[i__ + 5] = 1e3;
/* L7: */
	l14_10.xu[i__ + 9] = 150.;
    }
    l14_10.xu[3] = .1;
    l14_10.xu[4] = .9;
    l14_10.xu[5] = .9;
    l14_10.xu[9] = 500.;
    for (i__ = 1; i__ <= 13; ++i__) {
	l11_10.lxl[i__ - 1] = true;
	l12_10.lxu[i__ - 1] = true;
	for (j = 1; j <= 15; ++j) {
/* L32: */
	    l5_24.gg[j + i__ * 15 - 16] = 0.;
	}
    }
    for (i__ = 1; i__ <= 10; ++i__) {
/* L30: */
	l4_10.gf[i__ - 1] = 0.;
    }
    for (i__ = 11; i__ <= 13; ++i__) {
/* L31: */
	l4_10.gf[i__ - 1] = 1.;
    }
    l5_24.gg[15] = -1.;
    l5_24.gg[30] = 1.;
    l5_24.gg[1] = -1.;
    l5_24.gg[16] = 1.;
    l5_24.gg[92] = -.002;
    l5_24.gg[107] = .002;
    l5_24.gg[153] = 1.;
    l5_24.gg[168] = 1.;
    l5_24.gg[183] = 1.;
    l5_24.gg[154] = -1.;
    l5_24.gg[169] = -1.;
    l5_24.gg[184] = -1.;
    l5_24.gg[185] = 1.;
    l5_24.gg[163] = 1.;
    l5_24.gg[179] = 1.;
    l20_11.lex = false;
    l20_11.xex[0] = .803770278595;
    l20_11.xex[1] = .89998603367;
    l20_11.xex[2] = .970972419495;
    l20_11.xex[3] = .0999995162129;
    l20_11.xex[4] = .190815447786;
    l20_11.xex[5] = .460571745738;
    l20_11.xex[6] = 574.080310673;
    l20_11.xex[7] = 74.0804261398;
    l20_11.xex[8] = 500.016155317;
    l20_11.xex[9] = .0999999999985;
    l20_11.xex[10] = 20.2341325935;
    l20_11.xex[11] = 77.3475459898;
    l20_11.xex[12] = .00673039736648;
    l20_11.fex = 97.5884089805;
    return 0;
labelL2:
    l6_1.fx = l2_10.x[10] + l2_10.x[11] + l2_10.x[12];
labelL3:
    return 0;
labelL4:
    if (l9_14.index1[0]) {
	l3_13.g[0] = l2_10.x[2] - l2_10.x[1];
    }
    if (l9_14.index1[1]) {
	l3_13.g[1] = l2_10.x[1] - l2_10.x[0];
    }
    if (l9_14.index1[2]) {
	l3_13.g[2] = 1. - (l2_10.x[6] - l2_10.x[7]) * .002;
    }
    if (l9_14.index1[3]) {
	l3_13.g[3] = l2_10.x[10] + l2_10.x[11] + l2_10.x[12] - 50.;
    }
    if (l9_14.index1[4]) {
	l3_13.g[4] = 250. - l2_10.x[10] - l2_10.x[11] - l2_10.x[12];
    }
    if (l9_14.index1[5]) {
	l3_13.g[5] = l2_10.x[12] - l2_10.x[9] * 1.262626 + l2_10.x[2] * 
		1.231059 * l2_10.x[9];
    }
    if (l9_14.index1[6]) {
/* Computing 2nd power */
	d__1 = l2_10.x[1];
	l3_13.g[6] = l2_10.x[4] - l2_10.x[1] * .03475 - l2_10.x[1] * .975 * 
		l2_10.x[4] + d__1 * d__1 * .00975;
    }
    if (l9_14.index1[7]) {
/* Computing 2nd power */
	d__1 = l2_10.x[2];
	l3_13.g[7] = l2_10.x[5] - l2_10.x[2] * .03475 - l2_10.x[2] * .975 * 
		l2_10.x[5] + d__1 * d__1 * .00975;
    }
    if (l9_14.index1[8]) {
	l3_13.g[8] = l2_10.x[4] * l2_10.x[6] - l2_10.x[0] * l2_10.x[7] - 
		l2_10.x[3] * l2_10.x[6] + l2_10.x[3] * l2_10.x[7];
    }
    if (l9_14.index1[9]) {
	l3_13.g[9] = (l2_10.x[1] * l2_10.x[8] + l2_10.x[4] * l2_10.x[7] - 
		l2_10.x[0] * l2_10.x[7] - l2_10.x[5] * l2_10.x[8]) * -.002 - 
		l2_10.x[5] - l2_10.x[4] + 1.;
    }
    if (l9_14.index1[10]) {
	l3_13.g[10] = l2_10.x[1] * l2_10.x[8] - l2_10.x[2] * l2_10.x[9] - 
		l2_10.x[5] * l2_10.x[8] - (l2_10.x[1] - l2_10.x[5]) * 500. + 
		l2_10.x[1] * l2_10.x[9];
    }
    if (l9_14.index1[11]) {
	l3_13.g[11] = l2_10.x[1] - .9 - (l2_10.x[1] * l2_10.x[9] - l2_10.x[2] 
		* l2_10.x[9]) * .002;
    }
    if (l9_14.index1[12]) {
/* Computing 2nd power */
	d__1 = l2_10.x[0];
	l3_13.g[12] = l2_10.x[3] - l2_10.x[0] * .03475 - l2_10.x[0] * .975 * 
		l2_10.x[3] + d__1 * d__1 * .00975;
    }
    if (l9_14.index1[13]) {
	l3_13.g[13] = l2_10.x[10] - l2_10.x[7] * 1.262626 + l2_10.x[0] * 
		1.231059 * l2_10.x[7];
    }
    if (l9_14.index1[14]) {
	l3_13.g[14] = l2_10.x[11] - l2_10.x[8] * 1.262626 + l2_10.x[1] * 
		1.231059 * l2_10.x[8];
    }
    return 0;
labelL5:
    if (! l10_14.index2[5]) {
	goto labelL12;
    }
    l5_24.gg[35] = l2_10.x[9] * 1.231059;
    l5_24.gg[140] = l2_10.x[2] * 1.231059 - 1.262626;
labelL12:
    if (! l10_14.index2[6]) {
	goto labelL13;
    }
    l5_24.gg[21] = -.03475 - l2_10.x[4] * .975 + l2_10.x[1] * .0195;
    l5_24.gg[66] = 1. - l2_10.x[1] * .975;
labelL13:
    if (! l10_14.index2[7]) {
	goto labelL14;
    }
    l5_24.gg[37] = -.03475 - l2_10.x[5] * .975 + l2_10.x[2] * .0195;
    l5_24.gg[82] = 1. - l2_10.x[2] * .975;
labelL14:
    if (! l10_14.index2[8]) {
	goto L15;
    }
    l5_24.gg[8] = -l2_10.x[7];
    l5_24.gg[53] = -l2_10.x[6] + l2_10.x[7];
    l5_24.gg[68] = l2_10.x[6];
    l5_24.gg[98] = l2_10.x[4] - l2_10.x[3];
    l5_24.gg[113] = -l2_10.x[0] + l2_10.x[3];
L15:
    if (! l10_14.index2[9]) {
	goto L16;
    }
    l5_24.gg[9] = l2_10.x[7] * .002;
    l5_24.gg[24] = l2_10.x[8] * -.002;
    l5_24.gg[69] = l2_10.x[7] * -.002 - 1.;
    l5_24.gg[84] = l2_10.x[8] * .002 - 1.;
    l5_24.gg[114] = (l2_10.x[4] - l2_10.x[0]) * -.002;
    l5_24.gg[129] = (l2_10.x[1] - l2_10.x[5]) * -.002;
L16:
    if (! l10_14.index2[10]) {
	goto L17;
    }
    l5_24.gg[25] = l2_10.x[8] - 500. + l2_10.x[9];
    l5_24.gg[40] = -l2_10.x[9];
    l5_24.gg[85] = -l2_10.x[8] + 500.;
    l5_24.gg[130] = l2_10.x[1] - l2_10.x[5];
    l5_24.gg[145] = -l2_10.x[2] + l2_10.x[1];
L17:
    if (! l10_14.index2[11]) {
	goto L18;
    }
    l5_24.gg[26] = 1. - l2_10.x[9] * .002;
    l5_24.gg[41] = l2_10.x[9] * .002;
    l5_24.gg[146] = (l2_10.x[2] - l2_10.x[1]) * .002;
L18:
    if (! l10_14.index2[12]) {
	goto L19;
    }
    l5_24.gg[12] = -.03475 - l2_10.x[3] * .975 + l2_10.x[0] * .0195;
    l5_24.gg[57] = 1. - l2_10.x[0] * .975;
L19:
    if (! l10_14.index2[13]) {
	goto labelL20;
    }
    l5_24.gg[13] = l2_10.x[7] * 1.231059;
    l5_24.gg[118] = l2_10.x[0] * 1.231059 - 1.262626;
labelL20:
    if (! l10_14.index2[14]) {
	goto L21;
    }
    l5_24.gg[29] = l2_10.x[8] * 1.231059;
    l5_24.gg[134] = l2_10.x[1] * 1.231059 - 1.262626;
L21:
    return 0;
} /* tp116_ */


/* Subroutine */ int tp117_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real a[50]	/* was [10][5] */, b[10], c__[25]	/* 
	    was [5][5] */, d__[5], e[5];
    static int i__, j;
    static Real t1, t2, t3, t4[5], t5[5], t6[5];

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL2;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 15;
    l1_1.nili = 0;
    l1_1.ninl = 5;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 15; ++i__) {
/* labelL20: */
	l2_11.x[i__ - 1] = .001;
    }
    l2_11.x[6] = 60.;
    for (i__ = 1; i__ <= 15; ++i__) {
	l13_11.xl[i__ - 1] = 0.;
	l11_11.lxl[i__ - 1] = true;
	l12_11.lxu[i__ - 1] = false;
	for (j = 1; j <= 5; ++j) {
/* L22: */
	    l5_25.gg[j + i__ * 5 - 6] = 0.;
	}
    }
    e[0] = -15.;
    e[1] = -27.;
    e[2] = -36.;
    e[3] = -18.;
    e[4] = -12.;
    c__[0] = 30.;
    c__[5] = -20.;
    c__[10] = -10.;
    c__[15] = 32.;
    c__[20] = -10.;
    c__[6] = 39.;
    c__[11] = -6.;
    c__[16] = -31.;
    c__[21] = 32.;
    c__[12] = 10.;
    c__[17] = -6.;
    c__[22] = -10.;
    c__[18] = 39.;
    c__[23] = -20.;
    c__[24] = 30.;
    for (i__ = 1; i__ <= 5; ++i__) {
	for (j = 1; j <= 5; ++j) {
/* L70: */
	    c__[j + i__ * 5 - 6] = c__[i__ + j * 5 - 6];
	}
    }
    d__[0] = 4.;
    d__[1] = 8.;
    d__[2] = 10.;
    d__[3] = 6.;
    d__[4] = 2.;
    for (i__ = 1; i__ <= 6; ++i__) {
	for (j = 1; j <= 5; ++j) {
/* L72: */
	    a[i__ + j * 10 - 11] = 0.;
	}
    }
    a[0] = -16.;
    a[10] = 2.;
    a[30] = 1.;
    a[11] = -2.;
    a[31] = .4;
    a[41] = 2.;
    a[2] = -3.5;
    a[22] = 2.;
    a[13] = -2.;
    a[33] = -4.;
    a[43] = -1.;
    a[14] = -9.;
    a[24] = -2.;
    a[34] = 1.;
    a[44] = -2.8;
    a[5] = 2.;
    a[25] = -4.;
    a[7] = -1.;
    a[17] = -2.;
    a[27] = -3.;
    a[37] = -2.;
    a[47] = -1.;
    for (i__ = 1; i__ <= 5; ++i__) {
	a[i__ * 10 - 4] = -1.;
	a[i__ * 10 - 2] = (Real) i__;
/* L73: */
	a[i__ * 10 - 1] = 1.;
    }
    b[0] = -40.;
    b[1] = -2.;
    b[2] = -.25;
    b[3] = -4.;
    b[4] = -4.;
    b[5] = -1.;
    b[6] = -40.;
    b[7] = -60.;
    b[8] = 5.;
    b[9] = 1.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L35: */
	l4_11.gf[i__ - 1] = -b[i__ - 1];
    }
    for (i__ = 1; i__ <= 10; ++i__) {
	for (j = 1; j <= 5; ++j) {
/* L40: */
	    l5_25.gg[j + i__ * 5 - 6] = -a[i__ + j * 10 - 11];
	}
    }
    for (i__ = 1; i__ <= 5; ++i__) {
	for (j = 1; j <= 5; ++j) {
	    if (i__ - j != 0) {
		goto L45;
	    } else {
		goto L42;
	    }
L45:
	    l5_25.gg[j + (i__ + 10) * 5 - 6] = c__[i__ + j * 5 - 6] * 2.;
L42:
	    ;
	}
/* L41: */
    }
    l20_12.lex = false;
    l20_12.xex[0] = 0.;
    l20_12.xex[1] = 0.;
    l20_12.xex[2] = 5.1741363068;
    l20_12.xex[3] = 0.;
    l20_12.xex[4] = 3.06109271525;
    l20_12.xex[5] = 11.839676029;
    l20_12.xex[6] = 0.;
    l20_12.xex[7] = 0.;
    l20_12.xex[8] = .103907059194;
    l20_12.xex[9] = 0.;
    l20_12.xex[10] = .299992902601;
    l20_12.xex[11] = .333470928832;
    l20_12.xex[12] = .399990975915;
    l20_12.xex[13] = .428314541579;
    l20_12.xex[14] = .223960749729;
    l20_12.fex = 32.3486789791;
    return 0;
labelL2:
    t1 = 0.;
    t2 = 0.;
    for (j = 1; j <= 5; ++j) {
/* Computing 3rd power */
	d__1 = l2_11.x[j + 9];
	t2 += d__[j - 1] * (d__1 * (d__1 * d__1));
	for (i__ = 1; i__ <= 5; ++i__) {
	    t1 += c__[i__ + j * 5 - 6] * l2_11.x[i__ + 9] * l2_11.x[j + 9];
/* L30: */
	}
    }
    t3 = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L31: */
	t3 += b[i__ - 1] * l2_11.x[i__ - 1];
    }
    for (j = 1; j <= 5; ++j) {
	t4[j - 1] = 0.;
	t5[j - 1] = 0.;
	for (i__ = 1; i__ <= 5; ++i__) {
/* L32: */
	    t4[j - 1] += c__[i__ + j * 5 - 6] * l2_11.x[i__ + 9];
	}
	for (i__ = 1; i__ <= 10; ++i__) {
/* L33: */
	    t5[j - 1] += a[i__ + j * 10 - 11] * l2_11.x[i__ - 1];
	}
/* L34: */
    }
    if (*mode == 4) {
	goto labelL4;
    }
    l6_1.fx = -(t3 - t1 - t2 * 2.);
    return 0;
labelL3:
    for (i__ = 1; i__ <= 5; ++i__) {
	t6[i__ - 1] = 0.;
	for (j = 1; j <= 5; ++j) {
/* L37: */
	    t6[i__ - 1] += (c__[i__ + j * 5 - 6] + c__[j + i__ * 5 - 6]) * 
		    l2_11.x[j + 9];
	}
    }
    for (i__ = 1; i__ <= 5; ++i__) {
/* L36: */
/* Computing 2nd power */
	d__1 = l2_11.x[i__ + 9];
	l4_11.gf[i__ + 9] = t6[i__ - 1] + d__[i__ - 1] * 6. * (d__1 * d__1);
    }
    return 0;
labelL4:
    for (j = 1; j <= 5; ++j) {
	if (l9_5.index1[j - 1]) {
/* Computing 2nd power */
	    d__1 = l2_11.x[j + 9];
	    l3_4.g[j - 1] = t4[j - 1] * 2. + d__[j - 1] * 3. * (d__1 * d__1) 
		    + e[j - 1] - t5[j - 1];
	}
/* L38: */
    }
    return 0;
labelL5:
    for (j = 1; j <= 5; ++j) {
	if (! l10_5.index2[j - 1]) {
	    goto L39;
	}
	l5_25.gg[j + (j + 10) * 5 - 6] = c__[j + j * 5 - 6] * 2. + d__[j - 1] 
		* 6. * l2_11.x[j + 9];
L39:
	;
    }
    return 0;
} /* tp117_ */

/* Subroutine */ int tp118_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__, j, k, m;
    static Real t;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 15;
    l1_1.nili = 29;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 15; ++i__) {
/* labelL6: */
	l2_11.x[i__ - 1] = 20.;
    }
    l2_11.x[1] = 55.;
    l2_11.x[2] = 15.;
    l2_11.x[4] = 60.;
    l2_11.x[7] = 60.;
    l2_11.x[10] = 60.;
    l2_11.x[13] = 60.;
    l13_11.xl[0] = 8.;
    l13_11.xl[1] = 43.;
    l13_11.xl[2] = 3.;
    l14_11.xu[0] = 21.;
    l14_11.xu[1] = 57.;
    l14_11.xu[2] = 16.;
    for (i__ = 1; i__ <= 4; ++i__) {
	l13_11.xl[i__ * 3] = 0.;
	l13_11.xl[i__ * 3 + 1] = 0.;
	l13_11.xl[i__ * 3 + 2] = 0.;
	l14_11.xu[i__ * 3] = 90.;
	l14_11.xu[i__ * 3 + 1] = 120.;
	l14_11.xu[i__ * 3 + 2] = 60.;
/* L22: */
    }
    for (i__ = 1; i__ <= 15; ++i__) {
	l11_11.lxl[i__ - 1] = true;
	l12_11.lxu[i__ - 1] = true;
	for (j = 1; j <= 29; ++j) {
/* L25: */
	    l5_26.gg[j + i__ * 29 - 30] = 0.;
	}
    }
    for (k = 1; k <= 4; ++k) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    l5_26.gg[k + (i__ << 2) - 4 + (k * 3 + i__) * 29 - 30] = 1.;
	    l5_26.gg[k + (i__ << 2) - 4 + (k * 3 + i__ - 3) * 29 - 30] = -1.;
	    l5_26.gg[k + 8 + (i__ << 2) + (k * 3 + i__) * 29 - 30] = -1.;
	    l5_26.gg[k + 8 + (i__ << 2) + (k * 3 + i__ - 3) * 29 - 30] = 1.;
/* labelL20: */
	}
    }
    for (k = 1; k <= 5; ++k) {
	for (i__ = 1; i__ <= 3; ++i__) {
/* L21: */
	    l5_26.gg[k + 24 + (k * 3 - 3 + i__) * 29 - 30] = 1.;
	}
    }
    l20_12.lex = false;
    l20_12.xex[0] = 8.;
    l20_12.xex[1] = 49.;
    l20_12.xex[2] = 3.;
    l20_12.xex[3] = 1.;
    l20_12.xex[4] = 56.;
    l20_12.xex[5] = 0.;
    l20_12.xex[6] = .999999999545;
    l20_12.xex[7] = 63.;
    l20_12.xex[8] = 6.;
    l20_12.xex[9] = 2.99999999965;
    l20_12.xex[10] = 70.;
    l20_12.xex[11] = 12.;
    l20_12.xex[12] = 4.99999999971;
    l20_12.xex[13] = 77.;
    l20_12.xex[14] = 18.;
    l20_12.fex = 664.820449993;
    return 0;
labelL2:
    t = 0.;
    for (m = 1; m <= 5; ++m) {
	i__ = m - 1;
/* L30: */
/* Computing 2nd power */
	d__1 = l2_11.x[i__ * 3];
/* Computing 2nd power */
	d__2 = l2_11.x[i__ * 3 + 1];
/* Computing 2nd power */
	d__3 = l2_11.x[i__ * 3 + 2];
	t = t + l2_11.x[i__ * 3] * 2.3 + d__1 * d__1 * 1e-4 + l2_11.x[i__ * 3 
		+ 1] * 1.7 + d__2 * d__2 * 1e-4 + l2_11.x[i__ * 3 + 2] * 2.2 
		+ d__3 * d__3 * 1.5e-4;
    }
    l6_1.fx = t;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 5; ++i__) {
	l4_11.gf[i__ * 3 - 3] = l2_11.x[i__ * 3 - 3] * 2e-4 + 2.3;
	l4_11.gf[i__ * 3 - 2] = l2_11.x[i__ * 3 - 2] * 2e-4 + 1.7;
/* L31: */
	l4_11.gf[i__ * 3 - 1] = l2_11.x[i__ * 3 - 1] * 3e-4 + 2.2;
    }
    return 0;
labelL4:
    for (i__ = 1; i__ <= 4; ++i__) {
	if (l9_15.index1[i__ - 1]) {
	    l3_14.g[i__ - 1] = l2_11.x[i__ * 3] - l2_11.x[i__ * 3 - 3] + 7.;
	}
	if (l9_15.index1[i__ + 3]) {
	    l3_14.g[i__ + 3] = l2_11.x[i__ * 3 + 1] - l2_11.x[i__ * 3 - 2] + 
		    7.;
	}
	if (l9_15.index1[i__ + 7]) {
	    l3_14.g[i__ + 7] = l2_11.x[i__ * 3 + 2] - l2_11.x[i__ * 3 - 1] + 
		    7.;
	}
	if (l9_15.index1[i__ + 11]) {
	    l3_14.g[i__ + 11] = l2_11.x[i__ * 3 - 3] - l2_11.x[i__ * 3] + 6.;
	}
	if (l9_15.index1[i__ + 15]) {
	    l3_14.g[i__ + 15] = l2_11.x[i__ * 3 - 2] - l2_11.x[i__ * 3 + 1] + 
		    7.;
	}
/* L32: */
	if (l9_15.index1[i__ + 19]) {
	    l3_14.g[i__ + 19] = l2_11.x[i__ * 3 - 1] - l2_11.x[i__ * 3 + 2] + 
		    6.;
	}
    }
    if (l9_15.index1[24]) {
	l3_14.g[24] = l2_11.x[0] + l2_11.x[1] + l2_11.x[2] - 60.;
    }
    if (l9_15.index1[25]) {
	l3_14.g[25] = l2_11.x[3] + l2_11.x[4] + l2_11.x[5] - 50.;
    }
    if (l9_15.index1[26]) {
	l3_14.g[26] = l2_11.x[6] + l2_11.x[7] + l2_11.x[8] - 70.;
    }
    if (l9_15.index1[27]) {
	l3_14.g[27] = l2_11.x[9] + l2_11.x[10] + l2_11.x[11] - 85.;
    }
    if (l9_15.index1[28]) {
	l3_14.g[28] = l2_11.x[12] + l2_11.x[13] + l2_11.x[14] - 100.;
    }
labelL5:
    return 0;
} /* tp118_ */

/* Subroutine */ int tp119_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real a[256]	/* was [16][16] */, b[128]	/* was [8][16]
	     */, c__[8];
    static int i__, j;
    static Real s[16], t;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 16;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 8;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 16; ++i__) {
	for (j = 1; j <= 16; ++j) {
/* L21: */
	    a[i__ + (j << 4) - 17] = 0.;
	}
	a[i__ + (i__ << 4) - 17] = 1.;
	for (j = 1; j <= 8; ++j) {
/* L81: */
	    b[j + (i__ << 3) - 9] = 0.;
	}
/* L80: */
    }
    for (i__ = 1; i__ <= 8; ++i__) {
/* L22: */
	b[i__ + ((i__ + 8) << 3) - 9] = 1.;
    }
    a[48] = 1.;
    a[96] = 1.;
    a[112] = 1.;
    a[240] = 1.;
    a[33] = 1.;
    a[97] = 1.;
    a[145] = 1.;
    a[98] = 1.;
    a[130] = 1.;
    a[146] = 1.;
    a[210] = 1.;
    a[99] = 1.;
    a[163] = 1.;
    a[227] = 1.;
    a[84] = 1.;
    a[148] = 1.;
    a[180] = 1.;
    a[244] = 1.;
    a[117] = 1.;
    a[229] = 1.;
    a[166] = 1.;
    a[198] = 1.;
    a[151] = 1.;
    a[231] = 1.;
    a[184] = 1.;
    a[248] = 1.;
    a[217] = 1.;
    a[202] = 1.;
    a[219] = 1.;
    a[220] = 1.;
    b[0] = .22;
    b[8] = .2;
    b[16] = .19;
    b[24] = .25;
    b[32] = .15;
    b[40] = .11;
    b[48] = .12;
    b[56] = .13;
    b[1] = -1.46;
    b[17] = -1.3;
    b[25] = 1.82;
    b[33] = -1.15;
    b[49] = .8;
    b[2] = 1.29;
    b[10] = -.89;
    b[34] = -1.16;
    b[42] = -.96;
    b[58] = -.49;
    b[3] = -1.1;
    b[11] = -1.06;
    b[19] = .95;
    b[27] = -.54;
    b[43] = -1.78;
    b[51] = -.41;
    b[28] = -1.43;
    b[36] = 1.51;
    b[44] = .59;
    b[52] = -.33;
    b[60] = -.43;
    b[13] = -1.72;
    b[21] = -.33;
    b[37] = 1.62;
    b[45] = 1.24;
    b[53] = .21;
    b[61] = -.26;
    b[6] = 1.12;
    b[30] = .31;
    b[54] = 1.12;
    b[70] = -.36;
    b[15] = .45;
    b[23] = .26;
    b[31] = -1.1;
    b[39] = .58;
    b[55] = -1.03;
    b[63] = .1;
    c__[0] = 2.5;
    c__[1] = 1.1;
    c__[2] = -3.1;
    c__[3] = -3.5;
    c__[4] = 1.3;
    c__[5] = 2.1;
    c__[6] = 2.3;
    c__[7] = -1.5;
    for (i__ = 1; i__ <= 16; ++i__) {
	l11_12.lxl[i__ - 1] = true;
	l12_12.lxu[i__ - 1] = true;
	l2_12.x[i__ - 1] = 10.;
	l13_12.xl[i__ - 1] = 0.;
	l14_12.xu[i__ - 1] = 5.;
	for (j = 1; j <= 8; ++j) {
/* labelL20: */
	    l5_27.gg[j + (i__ << 3) - 9] = b[j + (i__ << 3) - 9];
	}
    }
    l20_13.lex = false;
    l20_13.xex[0] = .0398473514099;
    l20_13.xex[1] = .791983155694;
    l20_13.xex[2] = .202870330224;
    l20_13.xex[3] = .844357916347;
    l20_13.xex[4] = 1.26990645286;
    l20_13.xex[5] = .934738707827;
    l20_13.xex[6] = 1.68196196924;
    l20_13.xex[7] = .15530087749;
    l20_13.xex[8] = 1.56787033356;
    l20_13.xex[9] = -3.59021173251e-12;
    l20_13.xex[10] = -6.12900888082e-12;
    l20_13.xex[11] = -8.86794857449e-13;
    l20_13.xex[12] = .660204066;
    l20_13.xex[13] = -2.54340725727e-12;
    l20_13.xex[14] = .674255926901;
    l20_13.xex[15] = -1.10433723798e-11;
    l20_13.fex = 244.899697515;
    return 0;
labelL2:
    t = 0.;
    for (i__ = 1; i__ <= 16; ++i__) {
	for (j = 1; j <= 16; ++j) {
/* L30: */
/* Computing 2nd power */
	    d__1 = l2_12.x[i__ - 1];
/* Computing 2nd power */
	    d__2 = l2_12.x[j - 1];
	    t += a[i__ + (j << 4) - 17] * (d__1 * d__1 + l2_12.x[i__ - 1] + 
		    1.) * (d__2 * d__2 + l2_12.x[j - 1] + 1.);
	}
    }
    l6_1.fx = t;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 16; ++i__) {
	s[i__ - 1] = 0.;
	for (j = 1; j <= 16; ++j) {
/* L32: */
/* Computing 2nd power */
	    d__1 = l2_12.x[j - 1];
	    s[i__ - 1] += (a[i__ + (j << 4) - 17] + a[j + (i__ << 4) - 17]) * 
		    (d__1 * d__1 + l2_12.x[j - 1] + 1.) * (l2_12.x[i__ - 1] * 
		    2. + 1.);
	}
/* L31: */
	l4_12.gf[i__ - 1] = s[i__ - 1];
    }
    return 0;
labelL4:
    for (i__ = 1; i__ <= 8; ++i__) {
	if (! l9_12.index1[i__ - 1]) {
	    goto L33;
	}
	s[i__ - 1] = 0.;
	for (j = 1; j <= 16; ++j) {
/* L34: */
	    s[i__ - 1] += b[i__ + (j << 3) - 9] * l2_12.x[j - 1];
	}
	l3_11.g[i__ - 1] = s[i__ - 1] - c__[i__ - 1];
L33:
	;
    }
labelL5:
    return 0;
} /* tp119_ */


/* Subroutine */ int tp201_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 8.;
    l2_1.x[1] = 9.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 5.;
    l20_1.xex[1] = 6.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 5.;
/* Computing 2nd power */
    d__2 = l2_1.x[1] - 6.;
    l6_1.fx = d__1 * d__1 * 4. + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = (l2_1.x[0] - 5.) * 8.;
    l4_1.gf[1] = (l2_1.x[1] - 6.) * 2.;
labelL4:
    return 0;
} /* tp201_ */


/* Subroutine */ int tp202_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 15.;
    l2_1.x[1] = -2.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = true;
/* labelL6: */
	l11_1.lxl[i__ - 1] = true;
    }
    l14_1.xu[0] = (float)20.;
    l14_1.xu[1] = (float)5.;
    l13_1.xl[0] = (float)1.;
    l13_1.xl[1] = (float)-5.;
    l20_6.lex = true;
    l20_6.nex = 2;
    l20_6.fex = 0.;
    l20_6.xex[0] = 5.;
    l20_6.xex[1] = 4.;
    l20_6.xex[2] = 48.98425;
    l20_6.xex[3] = -.89681;
    l15_1.lsum = 2;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[1];
/* Computing 3rd power */
    d__2 = l2_1.x[1];
    l16_1.f[0] = l2_1.x[0] - 13. - l2_1.x[1] * 2. + d__1 * d__1 * 5. - d__2 * 
	    (d__2 * d__2);
/* Computing 2nd power */
    d__1 = l2_1.x[1];
/* Computing 3rd power */
    d__2 = l2_1.x[1];
    l16_1.f[1] = l2_1.x[0] - 29. - l2_1.x[1] * 14. + d__1 * d__1 + d__2 * (
	    d__2 * d__2);
/* Computing 2nd power */
    d__1 = l16_1.f[0];
/* Computing 2nd power */
    d__2 = l16_1.f[1];
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[1];
/* Computing 3rd power */
    d__2 = l2_1.x[1];
    l16_1.f[0] = l2_1.x[0] - 13. - l2_1.x[1] * 2. + d__1 * d__1 * 5. - d__2 * 
	    (d__2 * d__2);
/* Computing 2nd power */
    d__1 = l2_1.x[1];
/* Computing 3rd power */
    d__2 = l2_1.x[1];
    l16_1.f[1] = l2_1.x[0] - 29. - l2_1.x[1] * 14. + d__1 * d__1 + d__2 * (
	    d__2 * d__2);
    l17_1.df[0] = 1.;
/* Computing 2nd power */
    d__1 = l2_1.x[1];
    l17_1.df[2] = l2_1.x[1] * 10. - 2. - d__1 * d__1 * 3.;
    l17_1.df[1] = 1.;
/* Computing 2nd power */
    d__1 = l2_1.x[1];
    l17_1.df[3] = l2_1.x[1] * 2. - 14. + d__1 * d__1 * 3.;
    for (i__ = 1; i__ <= 2; ++i__) {
/* L7: */
	l4_1.gf[i__ - 1] = l16_1.f[0] * 2. * l17_1.df[(i__ << 1) - 2] + 
		l16_1.f[1] * 2. * l17_1.df[(i__ << 1) - 1];
    }
labelL4:
    return 0;
} /* tp202_ */


/* Subroutine */ int tp203_(int *mode)
{
    /* Initialized data */

    static Real c__[3] = { 1.5,2.25,2.625 };

    /* System generated locals */
    Real d__1;

    /* Builtin functions */
    Real  pow_di(Real* , int*);

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 2.;
    l2_1.x[1] = .2;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 3.;
    l20_1.xex[1] = .5;
    l15_1.lsum = 3;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l16_2.f[i__ - 1] = c__[i__ - 1] - l2_1.x[0] * (1. - pow_di(&l2_1.x[1],
		 &i__));
/* L7: */
/* Computing 2nd power */
	d__1 = l16_2.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 3; ++i__) {
/* L8: */
	l16_2.f[i__ - 1] = c__[i__ - 1] - l2_1.x[0] * (1. - pow_di(&l2_1.x[1],
		 &i__));
    }
    l17_2.df[0] = l2_1.x[1] - 1.;
    l17_2.df[3] = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1];
    l17_2.df[1] = d__1 * d__1 - 1.;
    l17_2.df[4] = l2_1.x[0] * 2. * l2_1.x[1];
/* Computing 3rd power */
    d__1 = l2_1.x[1];
    l17_2.df[2] = d__1 * (d__1 * d__1) - 1.;
/* Computing 2nd power */
    d__1 = l2_1.x[1];
    l17_2.df[5] = l2_1.x[0] * 3. * (d__1 * d__1);
    for (i__ = 1; i__ <= 2; ++i__) {
/* labelL9: */
	l4_1.gf[i__ - 1] = l16_2.f[0] * 2. * l17_2.df[i__ * 3 - 3] + l16_2.f[
		1] * 2. * l17_2.df[i__ * 3 - 2] + l16_2.f[2] * 2. * l17_2.df[
		i__ * 3 - 1];
    }
labelL4:
    return 0;
} /* tp203_ */


/* Subroutine */ int tp204_(int *mode)
{
    /* Initialized data */

    static Real a[3] = { .13294,-.244378,.325895 };
    static Real d__[3] = { 2.5074,-1.36401,1.02282 };
    static Real h__[6]	/* was [3][2] */ = { -.564255,-.404979,
	    -.0735084,.392417,.927589,.535493 };

    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static Real prod;
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = .1;
    l2_1.x[1] = .1;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = .183601;
    l20_1.xex[0] = 0.;
    l20_1.xex[1] = 0.;
    l15_1.lsum = 3;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 3; ++i__) {
	prod = h__[i__ - 1] * l2_1.x[0] + h__[i__ + 2] * l2_1.x[1];
/* labelL10: */
/* Computing 2nd power */
	d__1 = prod;
	l16_2.f[i__ - 1] = a[i__ - 1] + prod + d__1 * d__1 * .5 * d__[i__ - 1]
		;
    }
/* Computing 2nd power */
    d__1 = l16_2.f[0];
/* Computing 2nd power */
    d__2 = l16_2.f[1];
/* Computing 2nd power */
    d__3 = l16_2.f[2];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 3; ++i__) {
	prod = h__[i__ - 1] * l2_1.x[0] + h__[i__ + 2] * l2_1.x[1];
/* labelL11: */
/* Computing 2nd power */
	d__1 = prod;
	l16_2.f[i__ - 1] = a[i__ - 1] + prod + d__1 * d__1 * .5 * d__[i__ - 1]
		;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
/* L7: */
	l17_2.df[i__ - 1] = h__[i__ - 1] + (h__[i__ - 1] * l2_1.x[0] + h__[
		i__ + 2] * l2_1.x[1]) * h__[i__ - 1] * d__[i__ - 1];
    }
    for (i__ = 1; i__ <= 3; ++i__) {
/* L8: */
	l17_2.df[i__ + 2] = h__[i__ + 2] + (h__[i__ - 1] * l2_1.x[0] + h__[
		i__ + 2] * l2_1.x[1]) * h__[i__ + 2] * d__[i__ - 1];
    }
    for (i__ = 1; i__ <= 2; ++i__) {
/* labelL9: */
	l4_1.gf[i__ - 1] = (l16_2.f[0] * l17_2.df[i__ * 3 - 3] + l16_2.f[1] * 
		l17_2.df[i__ * 3 - 2] + l16_2.f[2] * l17_2.df[i__ * 3 - 1]) * 
		2.;
    }
labelL4:
    return 0;
} /* tp204_ */


/* Subroutine */ int tp205_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 0.;
    l2_1.x[1] = 0.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 3.;
    l20_1.xex[1] = .5;
    l15_1.lsum = 3;
    return 0;
labelL2:
    l16_2.f[0] = 1.5 - l2_1.x[0] * (1. - l2_1.x[1]);
/* Computing 2nd power */
    d__1 = l2_1.x[1];
    l16_2.f[1] = 2.25 - l2_1.x[0] * (1. - d__1 * d__1);
/* Computing 3rd power */
    d__1 = l2_1.x[1];
    l16_2.f[2] = 2.625 - l2_1.x[0] * (1. - d__1 * (d__1 * d__1));
/* Computing 2nd power */
    d__1 = l16_2.f[0];
/* Computing 2nd power */
    d__2 = l16_2.f[1];
/* Computing 2nd power */
    d__3 = l16_2.f[2];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return 0;
labelL3:
    l16_2.f[0] = 1.5 - l2_1.x[0] * (1. - l2_1.x[1]);
/* Computing 2nd power */
    d__1 = l2_1.x[1];
    l16_2.f[1] = 2.25 - l2_1.x[0] * (1. - d__1 * d__1);
/* Computing 3rd power */
    d__1 = l2_1.x[1];
    l16_2.f[2] = 2.625 - l2_1.x[0] * (1. - d__1 * (d__1 * d__1));
    l17_2.df[0] = l2_1.x[1] - 1.;
    l17_2.df[3] = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1];
    l17_2.df[1] = d__1 * d__1 - 1.;
    l17_2.df[4] = l2_1.x[0] * 2. * l2_1.x[1];
/* Computing 3rd power */
    d__1 = l2_1.x[1];
    l17_2.df[2] = d__1 * (d__1 * d__1) - 1.;
/* Computing 2nd power */
    d__1 = l2_1.x[1];
    l17_2.df[5] = l2_1.x[0] * 3. * (d__1 * d__1);
    for (i__ = 1; i__ <= 2; ++i__) {
/* L7: */
	l4_1.gf[i__ - 1] = l16_2.f[0] * 2. * l17_2.df[i__ * 3 - 3] + l16_2.f[
		1] * 2. * l17_2.df[i__ * 3 - 2] + l16_2.f[2] * 2. * l17_2.df[
		i__ * 3 - 1];
    }
labelL4:
    return 0;
} /* tp205_ */


/* Subroutine */ int tp206_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -1.2;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = d__1 * d__1 + d__3 * d__3 * 100.;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[0] = (l2_1.x[1] - d__1 * d__1) * -2. * 2. * l2_1.x[0] - (1. - 
	    l2_1.x[0]) * 200.;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (l2_1.x[1] - d__1 * d__1) * 2.;
labelL4:
    return 0;
} /* tp206_ */


/* Subroutine */ int tp207_(int *mode)
{
    /* Initialized data */

    static Real c__ = 1.;

    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int lsum, i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -1.2;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 1.;
    lsum = 2;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = c__ * (d__1 * d__1) + d__3 * d__3;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[0] = c__ * -4. * (l2_1.x[1] - d__1 * d__1) * l2_1.x[0] - (1. - 
	    l2_1.x[0]) * 2.;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = c__ * 2. * (l2_1.x[1] - d__1 * d__1);
labelL4:
    return 0;
} /* tp207_ */


/* Subroutine */ int tp208_(int *mode)
{
    /* Initialized data */

    static Real c__ = 100.;

    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -1.2;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = c__ * (d__1 * d__1) + d__3 * d__3;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[0] = c__ * -4. * (l2_1.x[1] - d__1 * d__1) * l2_1.x[0] - (1. - 
	    l2_1.x[0]) * 2.;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = c__ * 2. * (l2_1.x[1] - d__1 * d__1);
labelL4:
    return 0;
} /* tp208_ */

/* Subroutine */ int tp209_(int *mode)
{
    /* Initialized data */

    static Real c__ = 1e4;

    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -1.2;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = c__ * (d__1 * d__1) + d__3 * d__3;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[0] = c__ * -4. * (l2_1.x[1] - d__1 * d__1) * l2_1.x[0] - (1. - 
	    l2_1.x[0]) * 2.;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = c__ * 2. * (l2_1.x[1] - d__1 * d__1);
labelL4:
    return 0;
} /* tp209_ */

/* Subroutine */ int tp210_(int *mode)
{
    /* Initialized data */

    static Real c__ = 1e6;

    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -1.2;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = (c__ * (d__1 * d__1) + d__3 * d__3) / c__;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[0] = (c__ * -4. * (l2_1.x[1] - d__1 * d__1) * l2_1.x[0] - (1. - 
	    l2_1.x[0]) * 2.) / c__;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = c__ * 2. * (l2_1.x[1] - d__1 * d__1) / c__;
labelL4:
    return 0;
} /* tp210_ */

/* Subroutine */ int tp211_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -1.2;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 1.;
    return 0;
labelL2:
/* Computing 3rd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * (d__2 * d__2);
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = d__1 * d__1 * 100. + d__3 * d__3;
    return 0;
labelL3:
/* Computing 3rd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[0];
    l4_1.gf[0] = (l2_1.x[1] - d__1 * (d__1 * d__1)) * -200. * 3. * (d__2 * 
	    d__2) - (1. - l2_1.x[0]) * 2.;
/* Computing 3rd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (l2_1.x[1] - d__1 * (d__1 * d__1)) * 200.;
labelL4:
    return 0;
} /* tp211_ */

/* Subroutine */ int tp212_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 2.;
    l2_1.x[1] = 0.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 0.;
    l20_1.xex[1] = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = (l2_1.x[0] + l2_1.x[1]) * 4.;
/* Computing 2nd power */
    d__3 = l2_1.x[0] - 2.;
/* Computing 2nd power */
    d__4 = l2_1.x[1];
/* Computing 2nd power */
    d__2 = (l2_1.x[0] + l2_1.x[1]) * 4. + (l2_1.x[0] - l2_1.x[1]) * (d__3 * 
	    d__3 + d__4 * d__4 - 1.);
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 2.;
/* Computing 2nd power */
    d__2 = l2_1.x[1];
/* Computing 2nd power */
    d__3 = l2_1.x[0] - 2.;
/* Computing 2nd power */
    d__4 = l2_1.x[1];
    l4_1.gf[0] = (l2_1.x[0] + l2_1.x[1]) * 32. + ((l2_1.x[0] + l2_1.x[1]) * 
	    4. + (l2_1.x[0] - l2_1.x[1]) * (d__1 * d__1 + d__2 * d__2 - 1.)) *
	     2. * (d__3 * d__3 + d__4 * d__4 - 1. + 4. + (l2_1.x[0] - l2_1.x[
	    1]) * 2. * (l2_1.x[0] - 2.));
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 2.;
/* Computing 2nd power */
    d__2 = l2_1.x[1];
/* Computing 2nd power */
    d__3 = l2_1.x[0] - 2.;
/* Computing 2nd power */
    d__4 = l2_1.x[1];
    l4_1.gf[1] = (l2_1.x[0] + l2_1.x[1]) * 32. + ((l2_1.x[0] + l2_1.x[1]) * 
	    4. + (l2_1.x[0] - l2_1.x[1]) * (d__1 * d__1 + d__2 * d__2 - 1.)) *
	     2. * (4. - d__3 * d__3 + d__4 * d__4 - 1. + (l2_1.x[0] - l2_1.x[
	    1]) * 2. * l2_1.x[1]);
labelL4:
    return 0;
} /* tp212_ */


/* Subroutine */ int tp213_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 3.;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0] - l2_1.x[1];
/* Computing 2nd power */
    d__3 = l2_1.x[0] - 1.;
/* Computing 4th power */
    d__1 = d__2 * d__2 * 10. + d__3 * d__3, d__1 *= d__1;
    l6_1.fx = d__1 * d__1 * 1e-6;
    return 0;
labelL3:
/* Computing 2nd power */
    d__2 = l2_1.x[0] - l2_1.x[1];
/* Computing 2nd power */
    d__3 = l2_1.x[0] - 1.;
/* Computing 3rd power */
    d__1 = d__2 * d__2 * 10. + d__3 * d__3;
    l4_1.gf[0] = d__1 * (d__1 * d__1) * 4. * ((l2_1.x[0] - l2_1.x[1]) * 20. + 
	    (l2_1.x[0] - 1.) * 2.) * 1e-6;
/* Computing 2nd power */
    d__2 = l2_1.x[0] - l2_1.x[1];
/* Computing 2nd power */
    d__3 = l2_1.x[0] - 1.;
/* Computing 3rd power */
    d__1 = d__2 * d__2 * 10. + d__3 * d__3;
    l4_1.gf[1] = d__1 * (d__1 * d__1) * 4. * 20. * (l2_1.x[1] - l2_1.x[0]) * 
	    1e-6;
labelL4:
    return 0;
} /* tp213_ */


/* Subroutine */ int tp214_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -1.2;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0] - l2_1.x[1];
/* Computing 2nd power */
    d__3 = l2_1.x[0] - 1.;
    d__1 = d__2 * d__2 * 10. + d__3 * d__3;
    l6_1.fx = pow_dd(&d__1, &c_b934) * (float)100.;
    return 0;
labelL3:
/* Computing 2nd power */
    d__2 = l2_1.x[0] - l2_1.x[1];
/* Computing 2nd power */
    d__3 = l2_1.x[0] - 1.;
    d__1 = d__2 * d__2 * 10. + d__3 * d__3;
    l4_1.gf[0] = .25 / pow_dd(&d__1, &c_b949) * (l2_1.x[0] * 22. - l2_1.x[1] *
	     20. - 2.) * (float)100.;
/* Computing 2nd power */
    d__2 = l2_1.x[0] - l2_1.x[1];
/* Computing 2nd power */
    d__3 = l2_1.x[0] - 1.;
    d__1 = d__2 * d__2 * 10. + d__3 * d__3;
    l4_1.gf[1] = .25 / pow_dd(&d__1, &c_b949) * 20. * (l2_1.x[1] - l2_1.x[0]) 
	    * (float)100.;
labelL4:
    return 0;
} /* tp214_ */


/* Subroutine */ int tp215_(int *mode)
{
    /* System generated locals */
    Real d__1;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 1.;
    l2_1.x[1] = 1.;
    l12_1.lxu[0] = false;
    l12_1.lxu[1] = false;
    l11_1.lxl[0] = true;
    l11_1.lxl[1] = false;
    l13_1.xl[0] = 0.;
    l5_1.gg[1] = 1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 0.;
    l20_1.xex[1] = 0.;
    return 0;
labelL2:
    l6_1.fx = l2_1.x[1];
    return 0;
labelL3:
    l4_1.gf[0] = 0.;
    l4_1.gf[1] = 1.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_1.g[0] = l2_1.x[1] - d__1 * d__1;
    }
    return 0;
labelL5:
    if (l10_2.index2[0]) {
	l5_1.gg[0] = l2_1.x[0] * -2.;
    }
    return 0;
} /* tp215_ */


/* Subroutine */ int tp216_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    l2_1.x[0] = -1.2;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l14_1.xu[i__ - 1] = (float)10.;
	l12_1.lxu[i__ - 1] = true;
	l13_1.xl[i__ - 1] = (float)-3.;
/* labelL6: */
	l11_1.lxl[i__ - 1] = true;
    }
    l5_1.gg[1] = -2.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 1.;
    l20_1.xex[0] = 2.;
    l20_1.xex[1] = 4.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = d__2 * d__2 - l2_1.x[1];
/* Computing 2nd power */
    d__3 = l2_1.x[0] - 1.;
    l6_1.fx = d__1 * d__1 * 100. + d__3 * d__3;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[0] = (d__1 * d__1 - l2_1.x[1]) * 400. * l2_1.x[0] + (l2_1.x[0] - 
	    1.) * 2.;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (d__1 * d__1 - l2_1.x[1]) * -200.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = l2_1.x[0] * (l2_1.x[0] - 4.) - l2_1.x[1] * 2. + 12.;
    }
    return 0;
labelL5:
    if (l10_2.index2[0]) {
	l5_1.gg[0] = l2_1.x[0] * 2. - 4.;
    }
    return 0;
} /* tp216_ */


/* Subroutine */ int tp217_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 1;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    l2_1.x[0] = 10.;
    l2_1.x[1] = 10.;
    l12_1.lxu[0] = false;
    l12_1.lxu[1] = false;
    l11_1.lxl[0] = true;
    l11_1.lxl[1] = false;
    l13_1.xl[0] = 0.;
    l5_2.gg[0] = 1.;
    l5_2.gg[2] = -2.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = -.8;
    l20_1.xex[0] = .6;
    l20_1.xex[1] = .8;
    return 0;
labelL2:
    l6_1.fx = -l2_1.x[1];
    return 0;
labelL3:
    l4_1.gf[0] = 0.;
    l4_1.gf[1] = -1.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_1.x[0] + 1. - l2_1.x[1] * 2.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_2.g[1] = d__1 * d__1 + d__2 * d__2 - 1.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    l5_2.gg[1] = l2_1.x[0] * 2.;
    l5_2.gg[3] = l2_1.x[1] * 2.;
L8:
    return 0;
} /* tp217_ */


/* Subroutine */ int tp218_(int *mode)
{
    /* System generated locals */
    Real d__1;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 9.;
    l2_1.x[1] = 100.;
    l12_1.lxu[0] = false;
    l12_1.lxu[1] = false;
    l11_1.lxl[0] = false;
    l11_1.lxl[1] = true;
    l13_1.xl[1] = 0.;
    l5_1.gg[1] = 1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = 0.;
    l20_1.xex[1] = 0.;
    return 0;
labelL2:
    l6_1.fx = l2_1.x[1];
    return 0;
labelL3:
    l4_1.gf[0] = 0.;
    l4_1.gf[1] = 1.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_1.g[0] = l2_1.x[1] - d__1 * d__1;
    }
    return 0;
labelL5:
    if (l10_2.index2[0]) {
	l5_1.gg[0] = l2_1.x[0] * -2.;
    }
    return 0;
} /* tp218_ */


/* Subroutine */ int tp219_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 2;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 10.;
	l12_3.lxu[i__ - 1] = false;
/* L8: */
	l11_3.lxl[i__ - 1] = false;
    }
    l5_6.gg[2] = 1.;
    l5_6.gg[6] = 0.;
    l5_6.gg[3] = -1.;
    l5_6.gg[5] = 0.;
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.fex = -1.;
    l20_6.xex[0] = 1.;
    l20_6.xex[1] = 1.;
    l20_6.xex[2] = 0.;
    l20_6.xex[3] = 0.;
    return 0;
labelL2:
    l6_1.fx = -l2_3.x[0];
    return 0;
labelL3:
    l4_3.gf[0] = -1.0;
    l4_3.gf[1] = 0.0;
    l4_3.gf[2] = 0.0;
    l4_3.gf[3] = 0.0;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 3rd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[2];
	l3_2.g[0] = l2_3.x[1] - d__1 * (d__1 * d__1) - d__2 * d__2;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[3];
	l3_2.g[1] = d__1 * d__1 - l2_3.x[1] - d__2 * d__2;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto labelL6;
    }
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l5_6.gg[0] = d__1 * d__1 * -3.;
    l5_6.gg[4] = l2_3.x[2] * -2.;
labelL6:
    if (! l10_3.index2[1]) {
	goto L7;
    }
    l5_6.gg[1] = l2_3.x[0] * 2.;
    l5_6.gg[7] = l2_3.x[3] * -2.;
L7:
    return 0;
} /* tp219_ */


/* Subroutine */ int tp220_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = (float)2.5e4;
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = true;
    }
    l13_1.xl[0] = (float)1.;
    l13_1.xl[1] = (float)0.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 1.;
    l20_1.xex[0] = (float)1.;
    l20_1.xex[1] = (float)0.;
    return 0;
labelL2:
    l6_1.fx = l2_1.x[0];
    return 0;
labelL3:
    l4_1.gf[0] = (float)1.;
    l4_1.gf[1] = (float)0.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 3rd power */
	d__1 = l2_1.x[0] - (float)1.;
	l3_1.g[0] = d__1 * (d__1 * d__1) - l2_1.x[1];
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
/* Computing 2nd power */
    d__1 = l2_1.x[0] - (float)1.;
    l5_1.gg[0] = d__1 * d__1 * (float)3.;
    l5_1.gg[1] = (float)-1.;
L7:
    return 0;
} /* tp220_ */


/* Subroutine */ int tp221_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = .25;
	l12_1.lxu[i__ - 1] = true;
	l14_1.xu[i__ - 1] = (float)1.;
	l11_1.lxl[i__ - 1] = true;
/* labelL6: */
	l13_1.xl[i__ - 1] = (float)0.;
    }
    l5_1.gg[1] = (float)-1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = -1.;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 0.;
    return 0;
labelL2:
    l6_1.fx = -l2_1.x[0];
    return 0;
labelL3:
    l4_1.gf[0] = -1.;
    l4_1.gf[1] = 0.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 3rd power */
	d__1 = l2_1.x[0] - 1.;
	l3_1.g[0] = -(d__1 * (d__1 * d__1)) - l2_1.x[1];
    }
    return 0;
labelL5:
    if (l10_2.index2[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0] - 1.;
	l5_1.gg[0] = d__1 * d__1 * -3.;
    }
    return 0;
} /* tp221_ */


/* Subroutine */ int tp222_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 1.3;
    l2_1.x[1] = .2;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
	l11_1.lxl[i__ - 1] = true;
/* labelL6: */
	l13_1.xl[i__ - 1] = 0.;
    }
    l5_1.gg[1] = -1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = -1.5;
    l20_1.xex[0] = 1.5;
    l20_1.xex[1] = 0.;
    return 0;
labelL2:
    l6_1.fx = -l2_1.x[0];
    return 0;
labelL3:
    l4_1.gf[0] = -1.;
    l4_1.gf[1] = 0.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 3rd power */
	d__1 = l2_1.x[0] - 1.;
	l3_1.g[0] = .125 - d__1 * (d__1 * d__1) - l2_1.x[1];
    }
    return 0;
labelL5:
    if (l10_2.index2[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0] - 1.;
	l5_1.gg[0] = d__1 * d__1 * -3.;
    }
    return 0;
} /* tp222_ */

/* Subroutine */ int tp223_(int *mode)
{
    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = .1;
    l2_1.x[1] = 3.3;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = true;
	l11_1.lxl[i__ - 1] = true;
	l14_1.xu[i__ - 1] = 10.;
/* labelL6: */
	l13_1.xl[i__ - 1] = 0.;
    }
    l14_1.xu[0] = (float)1.;
    l5_2.gg[2] = 0.;
    l5_2.gg[3] = 1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = -std::log(std::log(10.));
    l20_1.xex[0] = std::log(std::log(10.));
    l20_1.xex[1] = 10.;
    return 0;
labelL2:
    l6_1.fx = -l2_1.x[0];
    return 0;
labelL3:
    l4_1.gf[0] = -1.;
    l4_1.gf[1] = 0.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = std::exp(std::exp(l2_1.x[0]));
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = l2_1.x[1] - std::exp(std::exp(l2_1.x[0]));
    }
    return 0;
labelL5:
    if (l10_3.index2[0]) {
	l5_2.gg[0] = std::exp(l2_1.x[0]) * std::exp(std::exp(l2_1.x[0]));
    }
    if (l10_3.index2[1]) {
	l5_2.gg[1] = -std::exp(l2_1.x[0]) * std::exp(std::exp(l2_1.x[0]));
    }
    return 0;
} /* tp223_ */

/* Subroutine */ int tp224_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 4;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = .1;
	l12_1.lxu[i__ - 1] = true;
	l11_1.lxl[i__ - 1] = true;
	l14_1.xu[i__ - 1] = 6.;
/* labelL6: */
	l13_1.xl[i__ - 1] = 0.;
    }
    l5_6.gg[0] = 1.;
    l5_6.gg[4] = 3.;
    l5_6.gg[1] = -1.;
    l5_6.gg[5] = -3.;
    l5_6.gg[2] = 1.;
    l5_6.gg[6] = 1.;
    l5_6.gg[3] = -1.;
    l5_6.gg[7] = -1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = -304.;
    l20_1.xex[0] = 4.;
    l20_1.xex[1] = 4.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l6_1.fx = d__1 * d__1 * 2. + d__2 * d__2 - l2_1.x[0] * 48. - l2_1.x[1] * 
	    40.;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 4. - 48.;
    l4_1.gf[1] = l2_1.x[1] * 2. - 40.;
    return 0;
labelL4:
    if (l9_7.index1[0]) {
	l3_6.g[0] = l2_1.x[0] + l2_1.x[1] * 3.;
    }
    if (l9_7.index1[1]) {
	l3_6.g[1] = 18. - l2_1.x[0] - l2_1.x[1] * 3.;
    }
    if (l9_7.index1[2]) {
	l3_6.g[2] = l2_1.x[0] + l2_1.x[1];
    }
    if (l9_7.index1[3]) {
	l3_6.g[3] = 8. - l2_1.x[0] - l2_1.x[1];
    }
labelL5:
    return 0;
} /* tp224_ */

/* Subroutine */ int tp225_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 5;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 3.;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l5_4.gg[0] = 1.;
    l5_4.gg[5] = 1.;
    l5_4.gg[8] = -1.;
    l5_4.gg[4] = -1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 2.;
    l20_1.xex[0] = 1.;
    l20_1.xex[1] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2.;
    l4_1.gf[1] = l2_1.x[1] * 2.;
    return 0;
labelL4:
    if (l9_5.index1[0]) {
	l3_4.g[0] = l2_1.x[0] + l2_1.x[1] - 1.;
    }
    if (l9_5.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_4.g[1] = d__1 * d__1 + d__2 * d__2 - 1.;
    }
    if (l9_5.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_4.g[2] = d__1 * d__1 * 9. + d__2 * d__2 - 9.;
    }
    if (l9_5.index1[3]) {
	l3_4.g[3] = pow_dd(l2_1.x, &c_b305) - l2_1.x[1];
    }
    if (l9_5.index1[4]) {
	l3_4.g[4] = pow_dd(&l2_1.x[1], &c_b305) - l2_1.x[0];
    }
    return 0;
labelL5:
    if (! l10_5.index2[1]) {
	goto L7;
    }
    l5_4.gg[1] = l2_1.x[0] * 2.;
    l5_4.gg[6] = l2_1.x[1] * 2.;
L7:
    if (! l10_5.index2[2]) {
	goto L8;
    }
    l5_4.gg[2] = l2_1.x[0] * 18.;
    l5_4.gg[7] = l2_1.x[1] * 2.;
L8:
    if (l10_5.index2[3]) {
	l5_4.gg[3] = l2_1.x[0] * 2.;
    }
    if (l10_5.index2[4]) {
	l5_4.gg[9] = l2_1.x[1] * 2.;
    }
    return 0;
} /* tp225_ */

/* Subroutine */ int tp226_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = .8;
    l2_1.x[1] = .05;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
	l11_1.lxl[i__ - 1] = true;
	l13_1.xl[i__ - 1] = 0.;
/* labelL6: */
	l20_1.xex[i__ - 1] = 1. / std::sqrt(2.);
    }
    l20_1.lex = true;
    l20_1.fex = -.5;
    return 0;
labelL2:
    l6_1.fx = -l2_1.x[0] * l2_1.x[1];
    return 0;
labelL3:
    l4_1.gf[0] = -l2_1.x[1];
    l4_1.gf[1] = -l2_1.x[0];
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_2.g[0] = d__1 * d__1 + d__2 * d__2;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_2.g[1] = 1. - d__1 * d__1 - d__2 * d__2;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_2.gg[0] = l2_1.x[0] * 2.;
    l5_2.gg[2] = l2_1.x[1] * 2.;
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    l5_2.gg[1] = l2_1.x[0] * -2.;
    l5_2.gg[3] = l2_1.x[1] * -2.;
L8:
    return 0;
} /* tp226_ */

/* Subroutine */ int tp227_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = .5;
	l12_1.lxu[i__ - 1] = false;
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l20_1.xex[i__ - 1] = 1.;
    }
    l5_2.gg[2] = 1.;
    l5_2.gg[1] = 1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 2.;
/* Computing 2nd power */
    d__2 = l2_1.x[1] - 1.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = (l2_1.x[0] - 2.) * 2.;
    l4_1.gf[1] = (l2_1.x[1] - 1.) * 2.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_2.g[0] = -(d__1 * d__1) + l2_1.x[1];
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[1];
	l3_2.g[1] = l2_1.x[0] - d__1 * d__1;
    }
    return 0;
labelL5:
    if (l10_3.index2[0]) {
	l5_2.gg[0] = l2_1.x[0] * -2.;
    }
    if (l10_3.index2[1]) {
	l5_2.gg[3] = l2_1.x[1] * -2.;
    }
    return 0;
} /* tp227_ */

/* Subroutine */ int tp228_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 1;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 0.;
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.xex[0] = 0.;
    l20_1.xex[1] = -3.;
    l5_2.gg[0] = -1.;
    l5_2.gg[2] = -1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = -3.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l6_1.fx = d__1 * d__1 + l2_1.x[1];
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2.;
    l4_1.gf[1] = 1.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = -l2_1.x[0] - l2_1.x[1] + 1.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_2.g[1] = -(d__1 * d__1 + d__2 * d__2) + 9.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[1]) {
	goto L7;
    }
    l5_2.gg[1] = l2_1.x[0] * -2.;
    l5_2.gg[3] = l2_1.x[1] * -2.;
L7:
    return 0;
} /* tp228_ */

/* Subroutine */ int tp229_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -1.2;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = true;
	l11_1.lxl[i__ - 1] = true;
	l14_1.xu[i__ - 1] = 2.;
	l13_1.xl[i__ - 1] = -2.;
/* labelL6: */
	l20_1.xex[i__ - 1] = 1.;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = d__1 * d__1 * 100. + d__3 * d__3;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[0] = l2_1.x[0] * -400. * (l2_1.x[1] - d__1 * d__1) - (1. - l2_1.x[
	    0]) * 2.;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (l2_1.x[1] - d__1 * d__1) * 200.;
labelL4:
    return 0;
} /* tp229_ */

/* Subroutine */ int tp230_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 0.;
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.xex[0] = .5;
    l20_1.xex[1] = .375;
    l5_2.gg[2] = 1.;
    l5_2.gg[3] = 1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = .375;
    return 0;
labelL2:
    l6_1.fx = l2_1.x[1];
    return 0;
labelL3:
    l4_1.gf[0] = 0.;
    l4_1.gf[1] = 1.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 3rd power */
	d__2 = l2_1.x[0];
	l3_2.g[0] = d__1 * d__1 * -2. + d__2 * (d__2 * d__2) + l2_1.x[1];
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = 1. - l2_1.x[0];
/* Computing 3rd power */
	d__2 = 1. - l2_1.x[0];
	l3_2.g[1] = d__1 * d__1 * -2. + d__2 * (d__2 * d__2) + l2_1.x[1];
    }
    return 0;
labelL5:
    if (l10_3.index2[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l5_2.gg[0] = l2_1.x[0] * -4. + d__1 * d__1 * 3.;
    }
    if (l10_3.index2[1]) {
/* Computing 2nd power */
	d__1 = 1. - l2_1.x[0];
	l5_2.gg[1] = (1. - l2_1.x[0]) * 4. - d__1 * d__1 * 3.;
    }
    return 0;
} /* tp230_ */

/* Subroutine */ int tp231_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 2;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -1.2;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l20_1.xex[i__ - 1] = 1.;
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l5_2.gg[0] = .33333333333333331;
    l5_2.gg[2] = 1.;
    l5_2.gg[1] = -.33333333333333331;
    l5_2.gg[3] = 1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = d__1 * d__1 * 100. + d__3 * d__3;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[0] = l2_1.x[0] * -400. * (l2_1.x[1] - d__1 * d__1) - (1. - l2_1.x[
	    0]) * 2.;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (l2_1.x[1] - d__1 * d__1) * 200.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_1.x[0] / 3. + l2_1.x[1] + .1;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = -l2_1.x[0] / 3. + l2_1.x[1] + .1;
    }
labelL5:
    return 0;
} /* tp231_ */

/* Subroutine */ int tp232_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;
    static Real hv;

    hv = std::sqrt(3.);
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 3;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 2.;
    l2_1.x[1] = .5;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
	l11_1.lxl[i__ - 1] = true;
/* labelL6: */
	l13_1.xl[i__ - 1] = 0.;
    }
    l5_3.gg[0] = 1 / hv;
    l5_3.gg[3] = -1.;
    l5_3.gg[1] = 1.;
    l5_3.gg[4] = hv;
    l5_3.gg[2] = -1.;
    l5_3.gg[5] = -hv;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = (float)-1.;
    l20_1.xex[0] = 3.;
    l20_1.xex[1] = hv;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 3.;
/* Computing 3rd power */
    d__2 = l2_1.x[1];
    l6_1.fx = -(9. - d__1 * d__1) * (d__2 * (d__2 * d__2)) / (hv * 27.);
    return 0;
labelL3:
    l4_1.gf[0] = (l2_1.x[0] - 3.) * 2. * pow_dd(&l2_1.x[1], &c_b1523) / (hv * 
	    27.);
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 3.;
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l4_1.gf[1] = -(9. - d__1 * d__1) * (d__2 * d__2) / (hv * 9.);
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_1.x[0] / hv - l2_1.x[1];
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_1.x[0] + hv * l2_1.x[1];
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = 6. - l2_1.x[0] - hv * l2_1.x[1];
    }
labelL5:
    return 0;
} /* tp232_ */

/* Subroutine */ int tp233_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 1.2;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l20_1.xex[i__ - 1] = 1.;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = l2_1.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0];
    l6_1.fx = d__1 * d__1 * 100. + d__3 * d__3;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[0] = l2_1.x[0] * -400. * (l2_1.x[1] - d__1 * d__1) - (1. - l2_1.x[
	    0]) * 2.;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (l2_1.x[1] - d__1 * d__1) * 200.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_1.g[0] = d__1 * d__1 + d__2 * d__2 - .25;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_1.gg[0] = l2_1.x[0] * 2.;
    l5_1.gg[1] = l2_1.x[1] * 2.;
L7:
    return 0;
} /* tp233_ */


/* Subroutine */ int tp234_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 1.;
	l12_1.lxu[i__ - 1] = true;
	l11_1.lxl[i__ - 1] = true;
	l13_1.xl[i__ - 1] = .2;
	l14_1.xu[i__ - 1] = 2.;
/* labelL6: */
	l20_1.xex[i__ - 1] = .2;
    }
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = -.8;
    return 0;
labelL2:
/* Computing 4th power */
    d__1 = l2_1.x[1] - l2_1.x[0], d__1 *= d__1;
    l6_1.fx = d__1 * d__1 - (1. - l2_1.x[0]);
    return 0;
labelL3:
/* Computing 3rd power */
    d__1 = l2_1.x[1] - l2_1.x[0];
    l4_1.gf[0] = d__1 * (d__1 * d__1) * -4. + 1.;
/* Computing 3rd power */
    d__1 = l2_1.x[1] - l2_1.x[0];
    l4_1.gf[1] = d__1 * (d__1 * d__1) * 4.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_1.g[0] = -(d__1 * d__1) - d__2 * d__2 + 1.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_1.gg[0] = l2_1.x[0] * -2.;
    l5_1.gg[1] = l2_1.x[1] * -2.;
L7:
    return 0;
} /* tp234_ */


/* Subroutine */ int tp235_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    l2_2.x[0] = -2.;
    l2_2.x[1] = 3.;
    l2_2.x[2] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_2.lxu[i__ - 1] = false;
/* labelL6: */
	l11_2.lxl[i__ - 1] = false;
    }
    l5_5.gg[0] = 1.;
    l5_5.gg[1] = 0.;
    l20_3.xex[0] = -1.;
    l20_3.xex[1] = 1.;
    l20_3.xex[2] = 0.;
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.fex = .04;
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_2.x[0];
/* Computing 2nd power */
    d__1 = l2_2.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = l2_2.x[0] - 1.;
    l6_1.fx = d__1 * d__1 + d__3 * d__3 * .01;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
    l4_2.gf[0] = l2_2.x[0] * -4. * (l2_2.x[1] - d__1 * d__1) + (l2_2.x[0] - 
	    1.) * .02;
/* Computing 2nd power */
    d__1 = l2_2.x[0];
    l4_2.gf[1] = (l2_2.x[1] - d__1 * d__1) * 2.;
    l4_2.gf[2] = 0.0;

    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[2];
	l3_1.g[0] = l2_2.x[0] + d__1 * d__1 + 1.;
    }
    return 0;
labelL5:
    if (l10_2.index2[0]) {
	l5_5.gg[2] = l2_2.x[2] * 2.;
    }
    return 0;
} /* tp235_ */


/* Subroutine */ int tp236239_(int *imode)
{
    /* Initialized data */

    static Real b[20] = { 75.1963666677,-3.8112755343,.1269366345,
	    -.0020567665,1.0345e-5,-6.8306567613,.0302344793,-.0012813448,
	    3.52559e-5,-2.266e-7,.2564581253,-.003460403,1.35139e-5,
	    -28.1064434908,-5.2375e-6,-6.3e-9,7e-10,3.405462e-4,-1.6638e-6,
	    -2.8673112392 };

    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12, d__13, d__14, d__15, d__16, d__17;


    switch ((int)*imode) {
	case 1:  goto labelL2;
	case 2:  goto labelL3;
    }
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 3rd power */
    d__2 = l2_1.x[0];
/* Computing 4th power */
    d__3 = l2_1.x[0], d__3 *= d__3;
/* Computing 2nd power */
    d__4 = l2_1.x[0];
/* Computing 3rd power */
    d__5 = l2_1.x[0];
/* Computing 4th power */
    d__6 = l2_1.x[0], d__6 *= d__6;
/* Computing 2nd power */
    d__7 = l2_1.x[1];
/* Computing 3rd power */
    d__8 = l2_1.x[1];
/* Computing 4th power */
    d__9 = l2_1.x[1], d__9 *= d__9;
/* Computing 2nd power */
    d__10 = l2_1.x[0];
/* Computing 2nd power */
    d__11 = l2_1.x[1];
/* Computing 3rd power */
    d__12 = l2_1.x[0];
/* Computing 2nd power */
    d__13 = l2_1.x[1];
/* Computing 3rd power */
    d__14 = l2_1.x[0];
/* Computing 3rd power */
    d__15 = l2_1.x[1];
/* Computing 2nd power */
    d__16 = l2_1.x[1];
/* Computing 3rd power */
    d__17 = l2_1.x[1];
    l6_1.fx = b[0] + b[1] * l2_1.x[0] + b[2] * (d__1 * d__1) + b[3] * (d__2 * 
	    (d__2 * d__2)) + b[4] * (d__3 * d__3) + b[5] * l2_1.x[1] + b[6] * 
	    l2_1.x[0] * l2_1.x[1] + b[7] * (d__4 * d__4) * l2_1.x[1] + b[8] * 
	    (d__5 * (d__5 * d__5)) * l2_1.x[1] + b[9] * (d__6 * d__6) * 
	    l2_1.x[1] + b[10] * (d__7 * d__7) + b[11] * (d__8 * (d__8 * d__8))
	     + b[12] * (d__9 * d__9) + b[13] * (1. / (l2_1.x[1] + 1.)) + b[14]
	     * (d__10 * d__10) * (d__11 * d__11) + b[15] * (d__12 * (d__12 * 
	    d__12)) * (d__13 * d__13) + b[16] * (d__14 * (d__14 * d__14)) * (
	    d__15 * (d__15 * d__15)) + b[17] * l2_1.x[0] * (d__16 * d__16) + 
	    b[18] * l2_1.x[0] * (d__17 * (d__17 * d__17)) + b[19] * std::exp(
	    l2_1.x[0] * 5e-4 * l2_1.x[1]);
    l6_1.fx = -l6_1.fx;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 3rd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__3 = l2_1.x[0];
/* Computing 3rd power */
    d__4 = l2_1.x[0];
/* Computing 2nd power */
    d__5 = l2_1.x[1];
/* Computing 2nd power */
    d__6 = l2_1.x[0];
/* Computing 2nd power */
    d__7 = l2_1.x[1];
/* Computing 2nd power */
    d__8 = l2_1.x[0];
/* Computing 3rd power */
    d__9 = l2_1.x[1];
/* Computing 2nd power */
    d__10 = l2_1.x[1];
/* Computing 3rd power */
    d__11 = l2_1.x[1];
    l4_1.gf[0] = b[1] + b[2] * 2. * l2_1.x[0] + b[3] * 3. * (d__1 * d__1) + b[
	    4] * 4. * (d__2 * (d__2 * d__2)) + b[6] * l2_1.x[1] + b[7] * 2. * 
	    l2_1.x[0] * l2_1.x[1] + b[8] * 3. * (d__3 * d__3) * l2_1.x[1] + b[
	    9] * 4. * (d__4 * (d__4 * d__4)) * l2_1.x[1] + b[14] * 2. * 
	    l2_1.x[0] * (d__5 * d__5) + b[15] * 3. * (d__6 * d__6) * (d__7 * 
	    d__7) + b[16] * 3. * (d__8 * d__8) * (d__9 * (d__9 * d__9)) + b[
	    17] * (d__10 * d__10) + b[18] * (d__11 * (d__11 * d__11)) + b[19] 
	    * std::exp(l2_1.x[0] * 5e-4 * l2_1.x[1]) * (l2_1.x[1] * 5e-4);
    l4_1.gf[0] = -l4_1.gf[0];
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 3rd power */
    d__2 = l2_1.x[0];
/* Computing 4th power */
    d__3 = l2_1.x[0], d__3 *= d__3;
/* Computing 2nd power */
    d__4 = l2_1.x[1];
/* Computing 3rd power */
    d__5 = l2_1.x[1];
/* Computing 2nd power */
    d__6 = l2_1.x[1] + 1.;
/* Computing 2nd power */
    d__7 = l2_1.x[0];
/* Computing 3rd power */
    d__8 = l2_1.x[0];
/* Computing 3rd power */
    d__9 = l2_1.x[0];
/* Computing 2nd power */
    d__10 = l2_1.x[1];
/* Computing 2nd power */
    d__11 = l2_1.x[1];
    l4_1.gf[1] = b[5] + b[6] * l2_1.x[0] + b[7] * (d__1 * d__1) + b[8] * (
	    d__2 * (d__2 * d__2)) + b[9] * (d__3 * d__3) + b[10] * 2. * 
	    l2_1.x[1] + b[11] * 3. * (d__4 * d__4) + b[12] * 4. * (d__5 * (
	    d__5 * d__5)) + b[13] * (-1. / (d__6 * d__6)) + b[14] * (d__7 * 
	    d__7) * 2 * l2_1.x[1] + b[15] * (d__8 * (d__8 * d__8)) * 2 * 
	    l2_1.x[1] + b[16] * (d__9 * (d__9 * d__9)) * 3 * (d__10 * d__10) 
	    + b[17] * l2_1.x[0] * 2. * l2_1.x[1] + b[18] * l2_1.x[0] * 3. * (
	    d__11 * d__11) + b[19] * std::exp(l2_1.x[0] * 5e-4 * l2_1.x[1]) * (
	    l2_1.x[0] * 5e-4);
    l4_1.gf[1] = -l4_1.gf[1];
    return 0;
} /* tp236239_ */


/* Subroutine */ int tp236_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;
    int tp236239_(int *);

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 10.;
    l2_1.x[1] = 10.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = true;
	l11_1.lxl[i__ - 1] = true;
/* labelL6: */
	l13_1.xl[i__ - 1] = 0.;
    }
    l14_1.xu[0] = 75.;
    l14_1.xu[1] = 65.;
    l5_2.gg[3] = 1.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = -58.903436;
    l20_1.xex[0] = 75.;
    l20_1.xex[1] = 65.;
    return 0;
labelL2:
    tp236239_(&c__1);
    return 0;
labelL3:
    tp236239_(&c__2);
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_1.x[0] * l2_1.x[1] - 700.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0] / 25.;
	l3_2.g[1] = l2_1.x[1] - d__1 * d__1 * 5.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_2.gg[0] = l2_1.x[1];
    l5_2.gg[2] = l2_1.x[0];
L7:
    if (l10_3.index2[1]) {
	l5_2.gg[1] = l2_1.x[0] * -.016;
    }
    return 0;
} /* tp236_ */


/* Subroutine */ int tp237_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    int tp236239_(int *);

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 3;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 10.;
    l2_1.x[1] = 10.;
    l12_1.lxu[0] = true;
    l12_1.lxu[1] = true;
    l11_1.lxl[0] = true;
    l11_1.lxl[1] = false;
    l14_1.xu[0] = 75.;
    l14_1.xu[1] = 65.;
    l13_1.xl[0] = 54.;
    l5_3.gg[4] = 1.;
    l5_3.gg[2] = -5.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = -58.903436;
    l20_1.xex[0] = 75.;
    l20_1.xex[1] = 65.;
    return 0;
labelL2:
    tp236239_(&c__1);
    return 0;
labelL3:
    tp236239_(&c__2);
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_1.x[0] * l2_1.x[1] - 700.;
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0] / 25.;
	l3_3.g[1] = l2_1.x[1] - d__1 * d__1 * 5.;
    }
    if (l9_4.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_1.x[1] - 50.;
	l3_3.g[2] = d__1 * d__1 - (l2_1.x[0] - 55.) * 5.;
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto labelL6;
    }
    l5_3.gg[0] = l2_1.x[1];
    l5_3.gg[3] = l2_1.x[0];
labelL6:
    if (l10_4.index2[1]) {
	l5_3.gg[1] = l2_1.x[0] * -.016;
    }
    if (l10_4.index2[2]) {
	l5_3.gg[5] = l2_1.x[1] * 2. - 100.;
    }
    return 0;
} /* tp237_ */


/* Subroutine */ int tp238_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;
    int tp236239_(int *);

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 3;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 10.;
    l2_1.x[1] = 10.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = true;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l14_1.xu[0] = 75.;
    l14_1.xu[1] = 65.;
    l5_3.gg[4] = 1.;
    l5_3.gg[2] = -5.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = -58.903436;
    l20_1.xex[0] = 75.;
    l20_1.xex[1] = 65.;
    return 0;
labelL2:
    tp236239_(&c__1);
    return 0;
labelL3:
    tp236239_(&c__2);
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_1.x[0] * l2_1.x[1] - 700.;
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0] / 25.;
	l3_3.g[1] = l2_1.x[1] - d__1 * d__1 * 5.;
    }
    if (l9_4.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_1.x[1] - 50.;
	l3_3.g[2] = d__1 * d__1 - (l2_1.x[0] - 55.) * 5.;
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto L7;
    }
    l5_3.gg[0] = l2_1.x[1];
    l5_3.gg[3] = l2_1.x[0];
L7:
    if (l10_4.index2[1]) {
	l5_3.gg[1] = l2_1.x[0] * -.016;
    }
    if (l10_4.index2[2]) {
	l5_3.gg[5] = l2_1.x[1] * 2. - 100.;
    }
    return 0;
} /* tp238_ */


/* Subroutine */ int tp239_(int *mode)
{
    static int i__;
    int tp236239_(int *);

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 10.;
    l2_1.x[1] = 10.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = true;
	l11_1.lxl[i__ - 1] = true;
/* labelL6: */
	l13_1.xl[i__ - 1] = 0.;
    }
    l14_1.xu[0] = 75.;
    l14_1.xu[1] = 65.;
    l20_1.lex = true;
    l20_1.nex = 1;
    l20_1.fex = (float)-58.903436;
    l20_1.xex[0] = (float)75.;
    l20_1.xex[1] = (float)65.;
    return 0;
labelL2:
    tp236239_(&c__1);
    return 0;
labelL3:
    tp236239_(&c__2);
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = l2_1.x[0] * l2_1.x[1] - 700.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_1.gg[0] = l2_1.x[1];
    l5_1.gg[1] = l2_1.x[0];
L7:
    return 0;
} /* tp239_ */


/* Subroutine */ int tp240_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_2.x[0] = 100.;
    l2_2.x[1] = -1.;
    l2_2.x[2] = 2.5;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = false;
	l12_2.lxu[i__ - 1] = false;
/* labelL6: */
	l20_3.xex[i__ - 1] = 0.;
    }
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.fex = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0] - l2_2.x[1] + l2_2.x[2];
/* Computing 2nd power */
    d__2 = -l2_2.x[0] + l2_2.x[1] + l2_2.x[2];
/* Computing 2nd power */
    d__3 = l2_2.x[0] + l2_2.x[1] - l2_2.x[2];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
labelL3:
    l4_2.gf[0] = (l2_2.x[0] - l2_2.x[1] + l2_2.x[2]) * 2. - (l2_2.x[1] - 
	    l2_2.x[0] + l2_2.x[2]) * 2. + (l2_2.x[0] + l2_2.x[1] - l2_2.x[2]) 
	    * 2.;
    l4_2.gf[1] = (l2_2.x[1] - l2_2.x[0] + l2_2.x[2]) * 2. - (l2_2.x[0] - 
	    l2_2.x[1] + l2_2.x[2]) * 2. + (l2_2.x[0] + l2_2.x[1] - l2_2.x[2]) 
	    * 2.;
    l4_2.gf[2] = (l2_2.x[0] - l2_2.x[1] + l2_2.x[2]) * 2. + (l2_2.x[1] - 
	    l2_2.x[0] + l2_2.x[2]) * 2. - (l2_2.x[0] + l2_2.x[1] - l2_2.x[2]) 
	    * 2.;
labelL4:
    return 0;
} /* tp240_ */

/* Subroutine */ int tp241_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l15_1.lsum = 5;
    l2_2.x[0] = 1.;
    l2_2.x[1] = 2.;
    l2_2.x[2] = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.fex = 0.;
    l20_3.xex[0] = 0.;
    l20_3.xex[1] = 0.;
    l20_3.xex[2] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1];
/* Computing 2nd power */
    d__3 = l2_2.x[2];
    l16_3.f[0] = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 - 1.;
/* Computing 2nd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1];
/* Computing 2nd power */
    d__3 = l2_2.x[2] - 2.;
    l16_3.f[1] = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 - 1.;
    l16_3.f[2] = l2_2.x[0] + l2_2.x[1] + l2_2.x[2] - 1.;
    l16_3.f[3] = l2_2.x[0] + l2_2.x[1] - l2_2.x[2] + 1.;
/* Computing 3rd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1];
/* Computing 2nd power */
    d__3 = l2_2.x[2] * 5. - l2_2.x[0] + 1.;
    l16_3.f[4] = d__1 * (d__1 * d__1) + d__2 * d__2 * 3. + d__3 * d__3 - 36.;
    if (*mode == 3) {
	goto labelL3;
    }
/* Computing 2nd power */
    d__1 = l16_3.f[0];
/* Computing 2nd power */
    d__2 = l16_3.f[1];
/* Computing 2nd power */
    d__3 = l16_3.f[2];
/* Computing 2nd power */
    d__4 = l16_3.f[3];
/* Computing 2nd power */
    d__5 = l16_3.f[4];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4 + d__5 * 
	    d__5;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 3; ++i__) {
	l17_3.df[i__ * 5 - 5] = l2_2.x[i__ - 1] * 2.;
/* L7: */
	l17_3.df[i__ * 5 - 3] = 1.;
    }
    l17_3.df[1] = l17_3.df[0];
    l17_3.df[6] = l17_3.df[5];
    l17_3.df[11] = (l2_2.x[2] - 2.) * 2.;
    l17_3.df[3] = 1.;
    l17_3.df[8] = 1.;
    l17_3.df[13] = -1.;
/* Computing 2nd power */
    d__1 = l2_2.x[0];
    l17_3.df[4] = d__1 * d__1 * 3. - (l2_2.x[2] * 5. - l2_2.x[0] + 1.) * 2.;
    l17_3.df[9] = l2_2.x[1] * 6.;
    l17_3.df[14] = (l2_2.x[2] * 5. - l2_2.x[0] + 1.) * 10.;
    l4_2.gf[0] = 0.;
    l4_2.gf[1] = 0.;
    l4_2.gf[2] = 0.;
    for (i__ = 1; i__ <= 5; ++i__) {
	l4_2.gf[0] += l16_3.f[i__ - 1] * l17_3.df[i__ - 1] * 2.;
	l4_2.gf[1] += l16_3.f[i__ - 1] * l17_3.df[i__ + 4] * 2.;
/* L8: */
	l4_2.gf[2] += l16_3.f[i__ - 1] * l17_3.df[i__ + 9] * 2.;
    }
labelL4:
    return 0;
} /* tp241_ */

/* Subroutine */ int tp242_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;
    static Real ti;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_2.x[0] = 2.5;
    l2_2.x[1] = 10.;
    l2_2.x[2] = l2_2.x[1];
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = true;
	l13_2.xl[i__ - 1] = 0.;
/* labelL6: */
	l14_2.xu[i__ - 1] = 10.;
    }
    l15_1.lsum = 10;
    l20_3.lex = true;
    l20_3.nex = -1;
    l20_3.fex = 0.;
    l20_3.xex[0] = 1.;
    l20_3.xex[1] = 10.;
    l20_3.xex[2] = 1.;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
	ti = (Real) i__ * .1;
/* L7: */
	l16_4.f[i__ - 1] = std::exp(-l2_2.x[0] * ti) - std::exp(-l2_2.x[1] * ti) - 
		l2_2.x[2] * (std::exp(-ti) - std::exp(ti * -10.));
    }
    if (*mode == 3) {
	goto labelL3;
    }
    for (i__ = 1; i__ <= 10; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__1 = l16_4.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    l4_2.gf[0] = 0.;
    l4_2.gf[1] = 0.;
    l4_2.gf[2] = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
	ti = (Real) i__ * .1;
	l16_4.f[i__ - 1] = std::exp(-l2_2.x[0] * ti) - std::exp(-l2_2.x[1] * ti) - 
		l2_2.x[2] * (std::exp(-ti) - std::exp(ti * -10.));
	l17_4.df[i__ - 1] = -ti * std::exp(-l2_2.x[0] * ti);
	l17_4.df[i__ + 9] = ti * std::exp(-l2_2.x[1] * ti);
	l17_4.df[i__ + 19] = std::exp(ti * -10.) - std::exp(-ti);
	l4_2.gf[0] += l16_4.f[i__ - 1] * l17_4.df[i__ - 1] * 2.;
	l4_2.gf[1] += l16_4.f[i__ - 1] * l17_4.df[i__ + 9] * 2.;
/* labelL9: */
	l4_2.gf[2] += l16_4.f[i__ - 1] * l17_4.df[i__ + 19] * 2.;
    }
labelL4:
    return 0;
} /* tp242_ */

/* Subroutine */ int tp243_(int *mode)
{
    /* Initialized data */

    static Real a[4] = { .14272,-.184981,-.521869,-.685306 };
    static Real d__[4] = { 1.75168,-1.35195,-.479048,-.3648 };
    static Real b[9]	/* was [3][3] */ = { 2.95137,4.87407,-2.0506,
	    4.87407,9.39321,-3.93185,-2.0506,-3.93189,2.64745 };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real dxbx[3], g[12]	/* was [4][3] */;
    static int i__;
    static Real xbx;

    g[0] = -.564255;
    g[4] = .392417;
    g[8] = -.404979;
    g[1] = .927589;
    g[5] = -.0735083;
    g[9] = .535493;
    g[2] = .658799;
    g[6] = -.636666;
    g[10] = -.681091;
    g[3] = -.869487;
    g[7] = .586387;
    g[11] = .289826;
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l15_1.lsum = 4;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = false;
	l12_2.lxu[i__ - 1] = false;
	l2_2.x[i__ - 1] = .1;
/* labelL6: */
	l20_3.xex[i__ - 1] = 0.;
    }
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.fex = .7966;
    return 0;
labelL2:
    l6_1.fx = 0.;
    xbx = (l2_2.x[0] * b[0] + l2_2.x[1] * b[1] + l2_2.x[2] * b[2]) * l2_2.x[0]
	     + (l2_2.x[0] * b[3] + l2_2.x[1] * b[4] + l2_2.x[2] * b[5]) * 
	    l2_2.x[1] + (l2_2.x[0] * b[6] + l2_2.x[1] * b[7] + l2_2.x[2] * b[
	    8]) * l2_2.x[2];
    for (i__ = 1; i__ <= 4; ++i__) {
/* L7: */
	l16_5.f[i__ - 1] = a[i__ - 1] + g[i__ - 1] * l2_2.x[0] + g[i__ + 3] * 
		l2_2.x[1] + g[i__ + 7] * l2_2.x[2] + xbx * .5 * d__[i__ - 1];
    }
    if (*mode == 3) {
	goto labelL3;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
/* labelL10: */
/* Computing 2nd power */
	d__1 = l16_5.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    dxbx[0] = (l2_2.x[0] * b[0] + l2_2.x[1] * b[1] + l2_2.x[2] * b[2]) * 2.;
    dxbx[1] = (l2_2.x[0] * b[3] + l2_2.x[1] * b[4] + l2_2.x[2] * b[5]) * 2.;
    dxbx[2] = (l2_2.x[0] * b[6] + l2_2.x[1] * b[7] + l2_2.x[2] * b[8]) * 2.;
    for (i__ = 1; i__ <= 3; ++i__) {
/* labelL9: */
	l4_2.gf[i__ - 1] = 0.;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	l17_5.df[i__ - 1] = g[i__ - 1] + dxbx[0] * d__[i__ - 1] * .5;
	l17_5.df[i__ + 3] = g[i__ + 3] + dxbx[1] * d__[i__ - 1] * .5;
	l17_5.df[i__ + 7] = g[i__ + 7] + dxbx[2] * d__[i__ - 1] * .5;
	l4_2.gf[0] += l16_5.f[i__ - 1] * l17_5.df[i__ - 1] * 2.;
	l4_2.gf[1] += l16_5.f[i__ - 1] * l17_5.df[i__ + 3] * 2.;
/* L8: */
	l4_2.gf[2] += l16_5.f[i__ - 1] * l17_5.df[i__ + 7] * 2.;
    }
labelL4:
    return 0;
} /* tp243_ */


/* Subroutine */ int tp244_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;
    static Real yi, zi;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l15_1.lsum = 10;
    l20_3.nex = 1;
    l2_2.x[0] = 1.;
    l2_2.x[1] = 2.;
    l2_2.x[2] = 1.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l13_2.xl[i__ - 1] = (float)0.;
	l14_2.xu[i__ - 1] = 1e10;
	l11_2.lxl[i__ - 1] = true;
/* labelL6: */
	l12_2.lxu[i__ - 1] = true;
    }
    l20_3.lex = true;
    l20_3.fex = 0.;
    l20_3.xex[0] = 1.;
    l20_3.xex[1] = 10.;
    l20_3.xex[2] = 5.;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
	zi = (Real) i__ * .1;
	yi = std::exp(-zi) - std::exp(zi * -10.) * 5.;
/* L7: */
	l16_4.f[i__ - 1] = std::exp(-l2_2.x[0] * zi) - l2_2.x[2] * std::exp(-l2_2.x[1] *
		 zi) - yi;
    }
    if (*mode == 3) {
	goto labelL3;
    }
    for (i__ = 1; i__ <= 8; ++i__) {
/* labelL10: */
/* Computing 2nd power */
	d__1 = l16_4.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 3; ++i__) {
/* L8: */
	l4_2.gf[i__ - 1] = 0.;
    }
    for (i__ = 1; i__ <= 10; ++i__) {
	zi = (Real) i__ * .1;
	yi = std::exp(-zi) - std::exp(zi * -10.) * 5.;
	l17_4.df[i__ - 1] = -zi * std::exp(-l2_2.x[0] * zi);
	l17_4.df[i__ + 9] = zi * l2_2.x[2] * std::exp(-l2_2.x[1] * zi);
	l17_4.df[i__ + 19] = -exp(-l2_2.x[1] * zi);
	l4_2.gf[0] += l16_4.f[i__ - 1] * l17_4.df[i__ - 1] * 2.;
	l4_2.gf[1] += l16_4.f[i__ - 1] * l17_4.df[i__ + 9] * 2.;
/* labelL9: */
	l4_2.gf[2] += l16_4.f[i__ - 1] * l17_4.df[i__ + 19] * 2.;
    }
labelL4:
    return 0;
} /* tp244_ */


/* Subroutine */ int tp245_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;
    static Real di;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l20_3.nex = -1;
    l2_2.x[0] = 0.;
    l2_2.x[1] = 10.;
    l2_2.x[2] = 20.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l13_13.xl[i__ - 1] = (float)0.;
	l11_2.lxl[i__ - 1] = true;
	l14_13.xu[i__ - 1] = (float)20.;
/* labelL6: */
	l12_2.lxu[i__ - 1] = true;
    }
    l14_13.xu[0] = (float)12.;
    l14_13.xu[1] = (float)12.;
    l20_3.lex = true;
    l20_3.fex = 0.;
    l20_3.xex[0] = 1.;
    l20_3.xex[1] = 10.;
    l20_3.xex[2] = 1.;
    l15_1.lsum = 10;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
	di = (Real) i__;
/* L7: */
	l16_4.f[i__ - 1] = std::exp(-di * l2_2.x[0] / 10.) - std::exp(-di * l2_2.x[1] / 
		10.) - l2_2.x[2] * (exp(-di / 10.) - std::exp(-di));
    }
    if (*mode == 3) {
	goto labelL3;
    }
    for (i__ = 1; i__ <= 10; ++i__) {
/* labelL10: */
/* Computing 2nd power */
	d__1 = l16_4.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 3; ++i__) {
/* L8: */
	l4_2.gf[i__ - 1] = 0.;
    }
    for (i__ = 1; i__ <= 10; ++i__) {
	di = (Real) i__;
	l17_4.df[i__ - 1] = -di / 10. * std::exp(-di * l2_2.x[0] / 10.);
	l17_4.df[i__ + 9] = di / 10. * std::exp(-di * l2_2.x[1] / 10.);
	l17_4.df[i__ + 19] = std::exp(-di) - std::exp(-di / 10.);
	l4_2.gf[0] += l16_4.f[i__ - 1] * l17_4.df[i__ - 1] * 2.;
	l4_2.gf[1] += l16_4.f[i__ - 1] * l17_4.df[i__ + 9] * 2.;
/* labelL9: */
	l4_2.gf[2] += l16_4.f[i__ - 1] * l17_4.df[i__ + 19] * 2.;
    }
labelL4:
    return 0;
} /* tp245_ */


/* Subroutine */ int tp246_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_2.x[0] = -1.2;
    l2_2.x[1] = 2.;
    l2_2.x[2] = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* L7: */
	    l17_6.df[i__ + j * 3 - 4] = 0.;
	}
	l11_2.lxl[i__ - 1] = false;
	l12_2.lxu[i__ - 1] = false;
/* labelL6: */
	l20_3.xex[i__ - 1] = 1.;
    }
    l17_6.df[6] = 10.;
    l17_6.df[1] = -1.;
    l17_6.df[5] = -1.;
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.fex = 0.;
    l15_1.lsum = 3;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = (l2_2.x[0] + l2_2.x[1]) / 2.;
    l16_2.f[0] = (l2_2.x[2] - d__1 * d__1) * 10.;
    l16_2.f[1] = 1. - l2_2.x[0];
    l16_2.f[2] = 1. - l2_2.x[1];
    if (*mode == 3) {
	goto labelL3;
    }
/* Computing 2nd power */
    d__1 = l16_2.f[0];
/* Computing 2nd power */
    d__2 = l16_2.f[1];
/* Computing 2nd power */
    d__3 = l16_2.f[2];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return 0;
labelL3:
    l17_6.df[0] = (l2_2.x[0] + l2_2.x[1]) * -10.;
    l17_6.df[3] = (l2_2.x[0] + l2_2.x[1]) * -10.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l4_2.gf[i__ - 1] = 0.;
	for (j = 1; j <= 3; ++j) {
/* L8: */
	    l4_2.gf[i__ - 1] += l16_2.f[j - 1] * 2. * l17_6.df[j + i__ * 3 - 
		    4];
	}
    }
labelL4:
    return 0;
} /* tp246_ */


/* Subroutine */ int tp247_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static int i__;
    static Real theta, dtheta[3], xpi;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l20_3.nex = 1;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l11_2.lxl[2] = true;
    l12_2.lxu[2] = true;
    l11_2.lxl[0] = true;
    l13_2.xl[0] = .1;
    l13_2.xl[2] = -2.5;
    l14_2.xu[2] = 7.5;
    l20_3.lex = true;
    l2_2.x[0] = -1.;
    l2_2.x[1] = 0.;
    l2_2.x[2] = 0.;
    l20_3.fex = 0.;
    l20_3.xex[0] = 1.;
    l20_3.xex[1] = 0.;
    l20_3.xex[2] = 0.;
    return 0;
labelL2:
    xpi = std::asin(1.) * 2.;
    theta = 1. / (xpi * 2.) * std::atan(l2_2.x[1] / l2_2.x[0]);
    if (l2_2.x[0] < 0.) {
	theta += .5;
    }
/* Computing 2nd power */
    d__1 = l2_2.x[2] - theta * 10.;
/* Computing 2nd power */
    d__3 = l2_2.x[0];
/* Computing 2nd power */
    d__4 = l2_2.x[1];
/* Computing 2nd power */
    d__2 = std::sqrt(d__3 * d__3 + d__4 * d__4) - 1.;
/* Computing 2nd power */
    d__5 = l2_2.x[2];
    l6_1.fx = (d__1 * d__1 + d__2 * d__2) * 100. + d__5 * d__5;
    return 0;
labelL3:
    xpi = std::asin(1.) * 2.;
    theta = 1. / (xpi * 2.) * std::atan(l2_2.x[1] / l2_2.x[0]);
/* Computing 2nd power */
    d__1 = l2_2.x[1] / l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[0];
    dtheta[0] = -l2_2.x[1] / ((d__1 * d__1 + 1.) * (d__2 * d__2));
/* Computing 2nd power */
    d__1 = l2_2.x[1] / l2_2.x[0];
    dtheta[1] = 1. / ((d__1 * d__1 + 1.) * l2_2.x[0]);
    dtheta[2] = 0.;
    if (l2_2.x[0] < 0.) {
	theta += .5;
    }
/* Computing 2nd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1];
/* Computing 2nd power */
    d__3 = l2_2.x[0];
/* Computing 2nd power */
    d__4 = l2_2.x[1];
    l4_2.gf[0] = ((l2_2.x[2] - theta * 10.) * 20. * dtheta[0] + (sqrt(d__1 * 
	    d__1 + d__2 * d__2) - 1.) * 2. / std::sqrt(d__3 * d__3 + d__4 * d__4) *
	     l2_2.x[0]) * 100.;
/* Computing 2nd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1];
/* Computing 2nd power */
    d__3 = l2_2.x[0];
/* Computing 2nd power */
    d__4 = l2_2.x[1];
    l4_2.gf[1] = ((l2_2.x[2] - theta * 10.) * 20. * dtheta[1] + (std::sqrt(d__1 * 
	    d__1 + d__2 * d__2) - 1.) * 2. / std::sqrt(d__3 * d__3 + d__4 * d__4) *
	     l2_2.x[1]) * 100.;
    l4_2.gf[2] = (l2_2.x[2] - theta * 10.) * 2. * 100. + l2_2.x[2] * 2.;
labelL4:
    return 0;
} /* tp247_ */

/* Subroutine */ int tp248_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 1;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    l2_2.x[0] = -.1;
    l2_2.x[1] = -1.;
    l2_2.x[2] = .1;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l5_3.gg[0] = 1.;
    l5_3.gg[2] = -2.;
    l5_3.gg[4] = 0.;
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.fex = -.8;
    l20_3.xex[0] = .6;
    l20_3.xex[1] = .8;
    l20_3.xex[2] = 0.;
    l4_2.gf[0] = 0.;
    l4_2.gf[1] = -1.;
    l4_2.gf[2] = 0.;
    return 0;
labelL2:
    l6_1.fx = -l2_2.x[1];
labelL3:
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = 1. - l2_2.x[1] * 2. + l2_2.x[0];
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[1];
/* Computing 2nd power */
	d__3 = l2_2.x[2];
	l3_2.g[1] = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 - 1.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
/* labelL9: */
	l5_3.gg[(i__ << 1) - 1] = l2_2.x[i__ - 1] * 2.;
    }
L8:
    return 0;
} /* tp248_ */

/* Subroutine */ int tp249_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 1.;
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l11_2.lxl[0] = true;
    l5_5.gg[2] = 0.;
    l13_2.xl[0] = 1.;
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.fex = 1.;
    l20_3.xex[0] = 1.;
    l20_3.xex[1] = 0.;
    l20_3.xex[2] = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1];
/* Computing 2nd power */
    d__3 = l2_2.x[2];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 3; ++i__) {
/* L7: */
	l4_2.gf[i__ - 1] = l2_2.x[i__ - 1] * 2.;
    }
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[1];
	l3_1.g[0] = d__1 * d__1 + d__2 * d__2 - 1.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L8;
    }
    l5_5.gg[0] = l2_2.x[0] * 2.;
    l5_5.gg[1] = l2_2.x[1] * 2.;
L8:
    return 0;
} /* tp249_ */

/* Subroutine */ int tp250_(int *mode)
{
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 2;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l20_3.nex = 1;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = true;
	l13_2.xl[i__ - 1] = 0.;
/* labelL6: */
	l2_2.x[i__ - 1] = 10.;
    }
    l14_2.xu[0] = 20.;
    l14_2.xu[1] = 11.;
    l14_2.xu[2] = 42.;
    l20_3.lex = true;
    l20_3.fex = -3300.;
    l20_3.xex[0] = 20.;
    l20_3.xex[1] = 11.;
    l20_3.xex[2] = 15.;
    l5_3.gg[0] = 1.;
    l5_3.gg[2] = 2.;
    l5_3.gg[4] = 2.;
    l5_3.gg[1] = -1.;
    l5_3.gg[3] = -2.;
    l5_3.gg[5] = -2.;
    return 0;
labelL2:
    l6_1.fx = -l2_2.x[0] * l2_2.x[1] * l2_2.x[2];
    return 0;
labelL3:
    l4_2.gf[0] = -l2_2.x[1] * l2_2.x[2];
    l4_2.gf[1] = -l2_2.x[0] * l2_2.x[2];
    l4_2.gf[2] = -l2_2.x[0] * l2_2.x[1];
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_2.x[0] + l2_2.x[1] * 2. + l2_2.x[2] * 2.;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = 72. - l2_2.x[0] - l2_2.x[1] * 2. - l2_2.x[2] * 2.;
    }
labelL5:
    return 0;
} /* tp250_ */

/* Subroutine */ int tp251_(int *mode)
{
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 1;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 10.;
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = true;
	l13_2.xl[i__ - 1] = 0.;
/* labelL6: */
	l14_2.xu[i__ - 1] = 42.;
    }
    l5_5.gg[0] = -1.;
    l5_5.gg[1] = -2.;
    l5_5.gg[2] = -2.;
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.fex = -3456.;
    l20_3.xex[0] = 24.;
    l20_3.xex[1] = 12.;
    l20_3.xex[2] = 12.;
    return 0;
labelL2:
    l6_1.fx = -l2_2.x[0] * l2_2.x[1] * l2_2.x[2];
    return 0;
labelL3:
    l4_2.gf[0] = -l2_2.x[1] * l2_2.x[2];
    l4_2.gf[1] = -l2_2.x[0] * l2_2.x[2];
    l4_2.gf[2] = -l2_2.x[0] * l2_2.x[1];
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = 72. - l2_2.x[0] - l2_2.x[1] * 2. - l2_2.x[2] * 2.;
    }
labelL5:
    return 0;
} /* tp251_ */

/* Subroutine */ int tp252_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    l2_2.x[0] = -1.;
    l2_2.x[1] = 2.;
    l2_2.x[2] = 2.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l12_2.lxu[0] = true;
    l14_2.xu[0] = -1.;
    l5_5.gg[0] = 1.;
    l5_5.gg[1] = 0.;
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.fex = .04;
    l20_3.xex[0] = -1.;
    l20_3.xex[1] = 1.;
    l20_3.xex[2] = 0.;
    l4_2.gf[2] = 0.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0] - 1.;
/* Computing 2nd power */
    d__3 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1] - d__3 * d__3;
    l6_1.fx = d__1 * d__1 * .01 + d__2 * d__2;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
    l4_2.gf[0] = (l2_2.x[0] - 1.) * .02 - (l2_2.x[1] - d__1 * d__1) * 4. * 
	    l2_2.x[0];
/* Computing 2nd power */
    d__1 = l2_2.x[0];
    l4_2.gf[1] = (l2_2.x[1] - d__1 * d__1) * 2.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[2];
	l3_1.g[0] = l2_2.x[0] + d__1 * d__1 + 1.;
    }
    return 0;
labelL5:
    if (l10_2.index2[0]) {
	l5_5.gg[2] = l2_2.x[2] * 2.;
    }
    return 0;
} /* tp252_ */


/* Subroutine */ int tp253_(int *mode)
{
    /* Initialized data */

    static Real a[24]	/* was [3][8] */ = { 0.,0.,0.,10.,0.,0.,10.,
	    10.,0.,0.,10.,0.,0.,0.,10.,10.,0.,10.,10.,10.,10.,0.,10.,10. };

    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__, j;

/*      DATA ((A(I,J),J=1,8),I=1,3) */
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 1;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_2.x[0] = 0.;
    l2_2.x[1] = 2.;
    l2_2.x[2] = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l12_2.lxu[i__ - 1] = false;
	l11_2.lxl[i__ - 1] = true;
/* labelL6: */
	l13_2.xl[i__ - 1] = 0.;
    }
    l5_5.gg[0] = -3.;
    l5_5.gg[1] = 0.;
    l5_5.gg[2] = -3.;
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.fex = 69.282032;
    for (i__ = 1; i__ <= 3; ++i__) {
/* L7: */
	l20_3.xex[i__ - 1] = (float).5;
    }
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (j = 1; j <= 8; ++j) {
/* L8: */
/* Computing 2nd power */
	d__1 = a[j * 3 - 3] - l2_2.x[0];
/* Computing 2nd power */
	d__2 = a[j * 3 - 2] - l2_2.x[1];
/* Computing 2nd power */
	d__3 = a[j * 3 - 1] - l2_2.x[2];
	l6_1.fx += std::sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 3; ++i__) {
	l4_2.gf[i__ - 1] = 0.;
	for (j = 1; j <= 8; ++j) {
/* labelL9: */
/* Computing 2nd power */
	    d__1 = a[j * 3 - 3] - l2_2.x[0];
/* Computing 2nd power */
	    d__2 = a[j * 3 - 2] - l2_2.x[1];
/* Computing 2nd power */
	    d__3 = a[j * 3 - 1] - l2_2.x[2];
	    l4_2.gf[i__ - 1] += (l2_2.x[i__ - 1] - a[i__ + j * 3 - 4]) / std::sqrt(
		    d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	}
    }
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = 30. - l2_2.x[0] * 3. - l2_2.x[2] * 3.;
    }
labelL5:
    return 0;
} /* tp253_ */


/* Subroutine */ int tp254_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;
    Real d_lg10(Real *);

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 2;
    l2_2.x[0] = 1.;
    l2_2.x[1] = 1.;
    l2_2.x[2] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_2.lxu[i__ - 1] = false;
/* labelL6: */
	l11_2.lxl[i__ - 1] = false;
    }
    l12_2.lxu[2] = false;
    l11_2.lxl[2] = true;
    l13_2.xl[2] = 1.;
    l5_3.gg[0] = 0.;
    l5_3.gg[3] = 0.;
    l5_3.gg[5] = 1.;
    l20_3.lex = true;
    l20_3.nex = 1;
    l20_3.fex = -std::sqrt(3.);
    l20_3.xex[0] = 0.;
    l20_3.xex[1] = std::sqrt(3.);
    l20_3.xex[2] = 1.;
    l4_2.gf[0] = 0.;
    l4_2.gf[1] = -1.;
    return 0;
labelL2:
    l6_1.fx = d_lg10(&l2_2.x[2]) - l2_2.x[1];
    return 0;
labelL3:
    l4_2.gf[2] = 1. / (l2_2.x[2] * std::log(10.));
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[1];
/* Computing 2nd power */
	d__2 = l2_2.x[2];
	l3_2.g[0] = d__1 * d__1 + d__2 * d__2 - 4.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
	l3_2.g[1] = l2_2.x[2] - 1. - d__1 * d__1;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_3.gg[2] = l2_2.x[1] * 2.;
    l5_3.gg[4] = l2_2.x[2] * 2.;
L7:
    if (l10_3.index2[1]) {
	l5_3.gg[1] = l2_2.x[0] * -2.;
    }
    return 0;
} /* tp254_ */


/* Subroutine */ int tp255_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_3.x[0] = -3.;
    l2_3.x[1] = 1.;
    l2_3.x[2] = -3.;
    l2_3.x[3] = 1.;
    for (i__ = 1; i__ <= 4; ++i__) {
	l14_3.xu[i__ - 1] = (float)10.;
	l12_3.lxu[i__ - 1] = true;
	l13_3.xl[i__ - 1] = (float)-10.;
/* labelL6: */
	l11_3.lxl[i__ - 1] = true;
    }
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.fex = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
/* L7: */
	l20_6.xex[i__ - 1] = 1.;
    }
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_3.x[0];
/* Computing 2nd power */
    d__2 = (float)1. - l2_3.x[0];
/* Computing 2nd power */
    d__3 = l2_3.x[2];
/* Computing 2nd power */
    d__4 = (float)1. - l2_3.x[2];
/* Computing 2nd power */
    d__5 = l2_3.x[1] - (float)1.;
/* Computing 2nd power */
    d__6 = l2_3.x[3] - (float)1.;
    l6_1.fx = (l2_3.x[1] - d__1 * d__1) * (float)100. + d__2 * d__2 + (l2_3.x[
	    3] - d__3 * d__3) * (float)90. + d__4 * d__4 + (d__5 * d__5 + 
	    d__6 * d__6) * (float)10.1 + (l2_3.x[1] - (float)1.) * (float)
	    19.8 * (l2_3.x[3] - (float)1.);
/* Computing 2nd power */
    d__1 = l6_1.fx;
    l6_1.fx = d__1 * d__1 * (float).5;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_3.x[0];
/* Computing 2nd power */
    d__2 = (float)1. - l2_3.x[0];
/* Computing 2nd power */
    d__3 = l2_3.x[2];
/* Computing 2nd power */
    d__4 = (float)1. - l2_3.x[2];
/* Computing 2nd power */
    d__5 = l2_3.x[1] - (float)1.;
/* Computing 2nd power */
    d__6 = l2_3.x[3] - (float)1.;
    l6_1.fx = (l2_3.x[1] - d__1 * d__1) * (float)100. + d__2 * d__2 + (l2_3.x[
	    3] - d__3 * d__3) * (float)90. + d__4 * d__4 + (d__5 * d__5 + 
	    d__6 * d__6) * (float)10.1 + (l2_3.x[1] - (float)1.) * (float)
	    19.8 * (l2_3.x[3] - (float)1.);
    l4_3.gf[0] = l6_1.fx * (l2_3.x[0] * -198. - 2.);
    l4_3.gf[1] = l6_1.fx * (l2_3.x[1] * 20.2 + l2_3.x[3] * 19.8 + 60.);
    l4_3.gf[2] = l6_1.fx * (l2_3.x[2] * -178. - 2.);
    l4_3.gf[3] = l6_1.fx * (l2_3.x[1] * 19.8 + l2_3.x[3] * 20.2 + 50.);
labelL4:
    return 0;
} /* tp255_ */


/* Subroutine */ int tp256_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_3.x[0] = 3.;
    l2_3.x[1] = -1.;
    l2_3.x[2] = 0.;
    l2_3.x[3] = 1.;
    for (i__ = 1; i__ <= 4; ++i__) {
	l12_3.lxu[i__ - 1] = false;
/* labelL6: */
	l11_3.lxl[i__ - 1] = false;
    }
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.fex = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
/* L7: */
	l20_6.xex[i__ - 1] = 0.;
    }
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_3.x[0] + l2_3.x[1] * 10.;
/* Computing 2nd power */
    d__2 = l2_3.x[2] - l2_3.x[3];
/* Computing 4th power */
    d__3 = l2_3.x[1] - l2_3.x[2] * 2., d__3 *= d__3;
/* Computing 4th power */
    d__4 = l2_3.x[0] - l2_3.x[3], d__4 *= d__4;
    l6_1.fx = (d__1 * d__1 + d__2 * d__2 * 5. + d__3 * d__3 + d__4 * d__4 * 
	    10.) * 1e-4;
    return 0;
labelL3:
/* Computing 3rd power */
    d__1 = l2_3.x[0] - l2_3.x[3];
    l4_3.gf[0] = ((l2_3.x[0] + l2_3.x[1] * 10.) * 2. + d__1 * (d__1 * d__1) * 
	    40.) * 1e-4;
/* Computing 3rd power */
    d__1 = l2_3.x[1] - l2_3.x[2] * 2.;
    l4_3.gf[1] = ((l2_3.x[0] + l2_3.x[1] * 10.) * 20. + d__1 * (d__1 * d__1) *
	     4.) * 1e-4;
/* Computing 3rd power */
    d__1 = l2_3.x[1] - l2_3.x[2] * 2.;
    l4_3.gf[2] = ((l2_3.x[2] - l2_3.x[3]) * 10. - d__1 * (d__1 * d__1) * 8.) *
	     1e-4;
/* Computing 3rd power */
    d__1 = l2_3.x[0] - l2_3.x[3];
    l4_3.gf[3] = ((l2_3.x[2] - l2_3.x[3]) * -10. - d__1 * (d__1 * d__1) * 40.)
	     * 1e-4;
labelL4:
    return 0;
} /* tp256_ */


/* Subroutine */ int tp257_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_3.x[0] = -3.;
    l2_3.x[1] = -1.;
    l2_3.x[2] = -3.;
    l2_3.x[3] = -1.;
    for (i__ = 1; i__ <= 4; ++i__) {
/* labelL6: */
	l12_3.lxu[i__ - 1] = false;
    }
    l11_3.lxl[0] = true;
    l11_3.lxl[1] = false;
    l11_3.lxl[2] = true;
    l11_3.lxl[3] = false;
    l13_3.xl[0] = 0.;
    l13_3.xl[2] = 0.;
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.fex = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
/* L7: */
	l20_6.xex[i__ - 1] = 1.;
    }
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_3.x[0];
/* Computing 2nd power */
    d__1 = d__2 * d__2 - l2_3.x[1];
/* Computing 2nd power */
    d__3 = l2_3.x[0] - 1.;
/* Computing 2nd power */
    d__5 = l2_3.x[2];
/* Computing 2nd power */
    d__4 = d__5 * d__5 - l2_3.x[3];
/* Computing 2nd power */
    d__6 = l2_3.x[2] - 1.;
/* Computing 2nd power */
    d__7 = l2_3.x[1] - 1.;
/* Computing 2nd power */
    d__8 = l2_3.x[3] - 1.;
    l6_1.fx = d__1 * d__1 * 100. + d__3 * d__3 + d__4 * d__4 * 90. + d__6 * 
	    d__6 + (d__7 * d__7 + d__8 * d__8) * 10.1 + (l2_3.x[0] - 1.) * 
	    19.8 * (l2_3.x[3] - 1.);
    return 0;
labelL3:
/* Computing 3rd power */
    d__1 = l2_3.x[0];
    l4_3.gf[0] = (d__1 * (d__1 * d__1) - l2_3.x[0] * l2_3.x[1]) * 400. + 
	    l2_3.x[0] * 2. + l2_3.x[3] * 19.8 - 21.8;
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l4_3.gf[1] = d__1 * d__1 * -200. + l2_3.x[1] * 220.2 - 20.2;
/* Computing 3rd power */
    d__1 = l2_3.x[2];
    l4_3.gf[2] = (d__1 * (d__1 * d__1) - l2_3.x[2] * l2_3.x[3]) * 360. + 
	    l2_3.x[2] * 2. - 2.;
/* Computing 2nd power */
    d__1 = l2_3.x[2];
    l4_3.gf[3] = d__1 * d__1 * -180. + l2_3.x[3] * 200.2 + l2_3.x[0] * 19.8 - 
	    40.;
labelL4:
    return 0;
} /* tp257_ */


/* Subroutine */ int tp258_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_3.x[0] = -3.;
    l2_3.x[1] = -1.;
    l2_3.x[2] = -3.;
    l2_3.x[3] = -1.;
    for (i__ = 1; i__ <= 4; ++i__) {
	l12_3.lxu[i__ - 1] = false;
/* labelL6: */
	l11_3.lxl[i__ - 1] = false;
    }
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.fex = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
/* L7: */
	l20_6.xex[i__ - 1] = 1.;
    }
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_3.x[0];
/* Computing 2nd power */
    d__1 = l2_3.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_3.x[0];
/* Computing 2nd power */
    d__5 = l2_3.x[2];
/* Computing 2nd power */
    d__4 = l2_3.x[3] - d__5 * d__5;
/* Computing 2nd power */
    d__6 = 1. - l2_3.x[2];
/* Computing 2nd power */
    d__7 = l2_3.x[1] - 1.;
/* Computing 2nd power */
    d__8 = l2_3.x[3] - 1.;
    l6_1.fx = (d__1 * d__1 * 100. + d__3 * d__3 + d__4 * d__4 * 90. + d__6 * 
	    d__6 + (d__7 * d__7 + d__8 * d__8) * 10.1 + (l2_3.x[1] - 1.) * 
	    19.8 * (l2_3.x[3] - 1.)) * 1e-4;
    return 0;
labelL3:
/* Computing 3rd power */
    d__1 = l2_3.x[0];
    l4_3.gf[0] = ((d__1 * (d__1 * d__1) - l2_3.x[0] * l2_3.x[1]) * 400. + 
	    l2_3.x[0] * 2. - 2.) * 1e-4;
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l4_3.gf[1] = (d__1 * d__1 * -200. + l2_3.x[1] * 220.2 + l2_3.x[3] * 19.8 
	    - 40.) * 1e-4;
/* Computing 3rd power */
    d__1 = l2_3.x[2];
    l4_3.gf[2] = ((d__1 * (d__1 * d__1) - l2_3.x[2] * l2_3.x[3]) * 360. + 
	    l2_3.x[2] * 2. - 2.) * 1e-4;
/* Computing 2nd power */
    d__1 = l2_3.x[2];
    l4_3.gf[3] = (d__1 * d__1 * -180. + l2_3.x[3] * 200.2 + l2_3.x[1] * 19.8 
	    - 40.) * 1e-4;
labelL4:
    return 0;
} /* tp258_ */


/* Subroutine */ int tp259_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 0.;
	l12_3.lxu[i__ - 1] = false;
/* labelL6: */
	l11_3.lxl[i__ - 1] = false;
    }
    l12_3.lxu[3] = true;
    l14_3.xu[3] = 1.;
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.fex = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
/* L7: */
	l20_6.xex[i__ - 1] = 1.;
    }
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_3.x[0];
/* Computing 2nd power */
    d__1 = l2_3.x[1] - d__2 * d__2;
/* Computing 2nd power */
    d__3 = 1. - l2_3.x[0];
/* Computing 2nd power */
    d__5 = l2_3.x[2];
/* Computing 2nd power */
    d__4 = l2_3.x[3] - d__5 * d__5;
/* Computing 3rd power */
    d__6 = 1. - l2_3.x[2];
/* Computing 2nd power */
    d__7 = l2_3.x[1] - 1.;
/* Computing 2nd power */
    d__8 = l2_3.x[3] - 1.;
    l6_1.fx = d__1 * d__1 * 100. + d__3 * d__3 + d__4 * d__4 * 90. + d__6 * (
	    d__6 * d__6) + d__7 * d__7 * 10.1 + d__8 * d__8 + (l2_3.x[1] - 1.)
	     * 19.8 * (l2_3.x[3] - 1.);
    return 0;
labelL3:
/* Computing 3rd power */
    d__1 = l2_3.x[0];
    l4_3.gf[0] = (d__1 * (d__1 * d__1) - l2_3.x[0] * l2_3.x[1]) * 400. + 
	    l2_3.x[0] * 2. - 2.;
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l4_3.gf[1] = d__1 * d__1 * -200. + l2_3.x[1] * 220.2 + l2_3.x[3] * 19.8 - 
	    40.;
/* Computing 3rd power */
    d__1 = l2_3.x[2];
/* Computing 2nd power */
    d__2 = 1. - l2_3.x[2];
    l4_3.gf[2] = (d__1 * (d__1 * d__1) - l2_3.x[2] * l2_3.x[3]) * 360. - d__2 
	    * d__2 * 3.;
/* Computing 2nd power */
    d__1 = l2_3.x[2];
    l4_3.gf[3] = d__1 * d__1 * -180. + l2_3.x[3] * 182. + l2_3.x[1] * 19.8 - 
	    21.8;
labelL4:
    return 0;
} /* tp259_ */


/* Subroutine */ int tp260_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_3.x[0] = -3.;
    l2_3.x[1] = -1.;
    l2_3.x[2] = -3.;
    l2_3.x[3] = -1.;
    for (i__ = 1; i__ <= 4; ++i__) {
	l12_3.lxu[i__ - 1] = false;
/* labelL6: */
	l11_3.lxl[i__ - 1] = false;
    }
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.fex = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
/* L7: */
	l20_6.xex[i__ - 1] = 1.;
    }
    l15_1.lsum = 7;
    for (i__ = 1; i__ <= 7; ++i__) {
	for (j = 1; j <= 4; ++j) {
/* L8: */
	    l17_7.df[i__ + j * 7 - 8] = 0.;
	}
    }
    l17_7.df[7] = 10.;
    l17_7.df[1] = -1.;
    l17_7.df[23] = std::sqrt(90.);
    l17_7.df[17] = -1.;
    l17_7.df[11] = std::sqrt(9.9);
    l17_7.df[25] = std::sqrt(9.9);
    l17_7.df[12] = std::sqrt(.2);
    l17_7.df[27] = std::sqrt(.2);
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l16_6.f[0] = (l2_3.x[1] - d__1 * d__1) * 10.;
    l16_6.f[1] = 1. - l2_3.x[0];
/* Computing 2nd power */
    d__1 = l2_3.x[2];
    l16_6.f[2] = std::sqrt(90.) * (l2_3.x[3] - d__1 * d__1);
    l16_6.f[3] = 1. - l2_3.x[2];
    l16_6.f[4] = std::sqrt(9.9) * (l2_3.x[1] - 1. + (l2_3.x[3] - 1.));
    l16_6.f[5] = std::sqrt(.2) * (l2_3.x[1] - 1.);
    l16_6.f[6] = std::sqrt(.2) * (l2_3.x[3] - 1.);
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 7; ++i__) {
/* labelL9: */
/* Computing 2nd power */
	d__1 = l16_6.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    l6_1.fx *= 1e-4;
    return 0;
labelL3:
    l17_7.df[0] = l2_3.x[0] * -20. * 1e-4;
    l17_7.df[16] = -std::sqrt(90.) * (float)2. * l2_3.x[2] * 1e-4;
    for (i__ = 1; i__ <= 4; ++i__) {
	l4_3.gf[i__ - 1] = 0.;
	for (j = 1; j <= 7; ++j) {
/* labelL10: */
	    l4_3.gf[i__ - 1] += l16_6.f[j - 1] * 2. * l17_7.df[j + i__ * 7 - 
		    8];
	}
    }
labelL4:
    return 0;
} /* tp260_ */


/* Subroutine */ int tp261_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static Real a, b, c__;
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 4;
    l15_1.lsum = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	l13_14.xl[i__ - 1] = (float)0.;
	l11_3.lxl[i__ - 1] = true;
	l14_14.xu[i__ - 1] = (float)1e5;
	l12_3.lxu[i__ - 1] = true;
	l2_3.x[i__ - 1] = 0.;
	for (j = 1; j <= 5; ++j) {
/* labelL6: */
	    l17_8.df[j + i__ * 5 - 6] = 0.;
	}
    }
    l17_8.df[19] = 1.;
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.fex = 0.;
    l20_6.xex[0] = 0.;
    for (i__ = 2; i__ <= 4; ++i__) {
/* L7: */
	l20_6.xex[i__ - 1] = 1.;
    }
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = std::exp(l2_3.x[0]) - l2_3.x[1];
    l16_3.f[0] = d__1 * d__1;
/* Computing 3rd power */
    d__1 = l2_3.x[1] - l2_3.x[2];
    l16_3.f[1] = d__1 * (d__1 * d__1) * 10.;
/* Computing 2nd power */
    d__1 = std::tan(l2_3.x[2] - l2_3.x[3]);
    l16_3.f[2] = d__1 * d__1;
/* Computing 4th power */
    d__1 = l2_3.x[0], d__1 *= d__1;
    l16_3.f[3] = d__1 * d__1;
    l16_3.f[4] = l2_3.x[3] - 1.;
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 5; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__1 = l16_3.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    a = std::exp(l2_3.x[0]) - l2_3.x[1];
    b = std::tan(l2_3.x[2] - l2_3.x[3]);
/* Computing 2nd power */
    d__1 =std::cos(l2_3.x[2] - l2_3.x[3]);
    c__ = b / (d__1 * d__1);
    l17_8.df[0] = std::exp(l2_3.x[0]) * 2. * a;
    l17_8.df[5] = a * -2.;
/* Computing 2nd power */
    d__1 = l2_3.x[1] - l2_3.x[2];
    l17_8.df[6] = d__1 * d__1 * 30.;
    l17_8.df[11] = -l17_8.df[6];
    l17_8.df[12] = c__ * 2.;
    l17_8.df[17] = -l17_8.df[12];
/* Computing 3rd power */
    d__1 = l2_3.x[0];
    l17_8.df[3] = d__1 * (d__1 * d__1) * 4.;
/* Computing 3rd power */
    d__1 = a;
/* Computing 7th power */
    d__2 = l2_3.x[0], d__3 = d__2, d__2 *= d__2, d__3 *= d__2;
    l4_3.gf[0] = std::exp(l2_3.x[0]) * 4. * (d__1 * (d__1 * d__1)) + d__3 * (d__2 *
	     d__2) * 8.;
/* Computing 3rd power */
    d__1 = a;
/* Computing 5th power */
    d__2 = l2_3.x[1] - l2_3.x[2], d__3 = d__2, d__2 *= d__2;
    l4_3.gf[1] = d__1 * (d__1 * d__1) * -4. + d__3 * (d__2 * d__2) * 600.;
/* Computing 5th power */
    d__1 = l2_3.x[1] - l2_3.x[2], d__2 = d__1, d__1 *= d__1;
    l4_3.gf[2] = b * 4. * b * c__ - d__2 * (d__1 * d__1) * 600.;
    l4_3.gf[3] = b * -4. * b * c__ + (l2_3.x[3] - 1.) * 2.;
labelL4:
    return 0;
} /* tp261_ */


/* Subroutine */ int tp262_(int *mode)
{
    /* Initialized data */

    static Real hgg[16]	/* was [4][4] */ = { -1.,-.2,-2.,1.,-1.,-.5,
	    -1.,1.,-1.,-1.,-.5,1.,-1.,-2.,-.2,-2. };

    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL3;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 3;
    l1_1.ninl = 0;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 1.;
	l11_3.lxl[i__ - 1] = true;
	l12_3.lxu[i__ - 1] = false;
	l13_3.xl[i__ - 1] = 0.;
	for (j = 1; j <= 4; ++j) {
	    l5_28.gg[i__ + (j << 2) - 5] = hgg[i__ + (j << 2) - 5];
/* labelL6: */
	}
    }
    for (i__ = 1; i__ <= 3; i__ += 2) {
	l4_3.gf[i__ - 1] = -.5;
/* L7: */
	l4_3.gf[i__] = -1.;
    }
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.fex = -10.;
    l20_6.xex[0] = 0.;
    l20_6.xex[1] = 8.6666666666666661;
    l20_6.xex[2] = 0.;
    l20_6.xex[3] = 1.3333333333333333;
    return 0;
labelL2:
    l6_1.fx = l2_3.x[0] * -.5 - l2_3.x[1] - l2_3.x[2] * .5 - l2_3.x[3];
labelL3:
    return 0;
labelL4:
    if (l9_7.index1[0]) {
	l3_6.g[0] = 10. - l2_3.x[0] - l2_3.x[1] - l2_3.x[2] - l2_3.x[3];
    }
    if (l9_7.index1[1]) {
	l3_6.g[1] = 10. - l2_3.x[0] * .2 - l2_3.x[1] * .5 - l2_3.x[2] - 
		l2_3.x[3] * 2.;
    }
    if (l9_7.index1[2]) {
	l3_6.g[2] = 10. - l2_3.x[0] * 2. - l2_3.x[1] - l2_3.x[2] * .5 - 
		l2_3.x[3] * .2;
    }
    if (l9_7.index1[3]) {
	l3_6.g[3] = l2_3.x[0] + l2_3.x[1] + l2_3.x[2] - l2_3.x[3] * 2. - 6.;
    }
    return 0;
} /* tp262_ */


/* Subroutine */ int tp263_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 2;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 10.;
	l11_3.lxl[i__ - 1] = false;
/* labelL6: */
	l12_3.lxu[i__ - 1] = false;
    }
    l5_28.gg[4] = 1.;
    l5_28.gg[8] = 0.;
    l5_28.gg[5] = -1.;
    l5_28.gg[9] = 0.;
    l5_28.gg[6] = 1.;
    l5_28.gg[7] = -1.;
    l5_28.gg[11] = 0.;
    l4_3.gf[0] = -1.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l5_28.gg[i__ + 11] = 0.;
/* L7: */
	l4_3.gf[i__] = 0.;
    }
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.fex = -1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l20_6.xex[i__ - 1] = 1.;
/* L8: */
	l20_6.xex[i__ + 1] = 0.;
    }
    return 0;
labelL2:
    l6_1.fx = -l2_3.x[0];
labelL3:
    return 0;
labelL4:
    if (l9_7.index1[0]) {
/* Computing 3rd power */
	d__1 = l2_3.x[0];
	l3_6.g[0] = l2_3.x[1] - d__1 * (d__1 * d__1);
    }
    if (l9_7.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_3.x[0];
	l3_6.g[1] = d__1 * d__1 - l2_3.x[1];
    }
    if (l9_7.index1[2]) {
/* Computing 3rd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[2];
	l3_6.g[2] = l2_3.x[1] - d__1 * (d__1 * d__1) - d__2 * d__2;
    }
    if (l9_7.index1[3]) {
/* Computing 2nd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[3];
	l3_6.g[3] = d__1 * d__1 - l2_3.x[1] - d__2 * d__2;
    }
    return 0;
labelL5:
    if (l10_7.index2[0]) {
/* Computing 2nd power */
	d__1 = l2_3.x[0];
	l5_28.gg[0] = d__1 * d__1 * -3.;
    }
    if (l10_7.index2[1]) {
	l5_28.gg[1] = l2_3.x[0] * 2.;
    }
    if (! l10_7.index2[2]) {
	goto labelL9;
    }
/* Computing 2nd power */
    d__1 = l2_3.x[0];
    l5_28.gg[2] = d__1 * d__1 * -3.;
    l5_28.gg[10] = l2_3.x[2] * -2.;
labelL9:
    if (! l10_7.index2[3]) {
	goto labelL10;
    }
    l5_28.gg[3] = l2_3.x[0] * 2.;
    l5_28.gg[15] = l2_3.x[3] * -2.;
labelL10:
    return 0;
} /* tp263_ */


/* Subroutine */ int tp264_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 3;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 0.;
	l11_3.lxl[i__ - 1] = false;
/* labelL6: */
	l12_3.lxu[i__ - 1] = false;
    }
    l5_7.gg[11] = 1.;
    l20_6.lex = true;
    l20_6.nex = 1;
    l20_6.fex = -.44;
    l20_6.xex[0] = 0.;
    l20_6.xex[1] = 1.;
    l20_6.xex[2] = 2.;
    l20_6.xex[3] = -1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_3.x[0];
/* Computing 2nd power */
    d__2 = l2_3.x[1];
/* Computing 2nd power */
    d__3 = l2_3.x[2];
/* Computing 2nd power */
    d__4 = l2_3.x[3];
    l6_1.fx = (d__1 * d__1 + d__2 * d__2 + d__3 * d__3 * 2. + d__4 * d__4 - 
	    l2_3.x[0] * 5. - l2_3.x[1] * 5. - l2_3.x[2] * 21. + l2_3.x[3] * 
	    7.) * (float).01;
    return 0;
labelL3:
    l4_3.gf[0] = (l2_3.x[0] * 2. - 5.) * (float).01;
    l4_3.gf[1] = (l2_3.x[1] * 2. - 5.) * (float).01;
    l4_3.gf[2] = (l2_3.x[2] * 4. - 21.) * (float).01;
    l4_3.gf[3] = (l2_3.x[3] * 2. + 7.) * (float).01;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[1];
/* Computing 2nd power */
	d__3 = l2_3.x[2];
/* Computing 2nd power */
	d__4 = l2_3.x[3];
	l3_3.g[0] = 8. - d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * 
		d__4 - l2_3.x[0] + l2_3.x[1] - l2_3.x[2] + l2_3.x[3];
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[1];
/* Computing 2nd power */
	d__3 = l2_3.x[2];
/* Computing 2nd power */
	d__4 = l2_3.x[3];
	l3_3.g[1] = 9. - d__1 * d__1 - d__2 * d__2 * 2. - d__3 * d__3 - d__4 *
		 d__4 * 2. + l2_3.x[0] + l2_3.x[3];
    }
    if (l9_4.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_3.x[0];
/* Computing 2nd power */
	d__2 = l2_3.x[1];
/* Computing 2nd power */
	d__3 = l2_3.x[2];
	l3_3.g[2] = 5. - d__1 * d__1 * 2. - d__2 * d__2 - d__3 * d__3 - 
		l2_3.x[0] * 2. + l2_3.x[1] + l2_3.x[3];
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto labelL9;
    }
    l5_7.gg[0] = -1. - l2_3.x[0] * 2.;
    l5_7.gg[3] = 1. - l2_3.x[1] * 2.;
    l5_7.gg[6] = -1. - l2_3.x[2] * 2.;
    l5_7.gg[9] = 1. - l2_3.x[3] * 2.;
labelL9:
    if (! l10_4.index2[1]) {
	goto labelL10;
    }
    l5_7.gg[1] = 1. - l2_3.x[0] * 2.;
    l5_7.gg[4] = l2_3.x[1] * -4.;
    l5_7.gg[7] = l2_3.x[2] * -2.;
    l5_7.gg[10] = -1. - l2_3.x[3] * 4.;
labelL10:
    if (! l10_4.index2[2]) {
	goto labelL11;
    }
    l5_7.gg[2] = -2. - l2_3.x[0] * 4.;
    l5_7.gg[5] = 1. - l2_3.x[1] * 2.;
    l5_7.gg[8] = l2_3.x[2] * -2.;
labelL11:
    return 0;
} /* tp264_ */


/* Subroutine */ int tp265_(int *mode)
{
    /* Initialized data */

    static Real hgg[8]	/* was [2][4] */ = { 1.,0.,1.,0.,0.,1.,0.,1. }
	    ;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 2;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_3.x[i__ - 1] = 0.;
	l11_3.lxl[i__ - 1] = true;
	l13_3.xl[i__ - 1] = 0.;
	l12_3.lxu[i__ - 1] = false;
	l14_3.xu[i__ - 1] = 1.;
	for (j = 1; j <= 2; ++j) {
	    l5_6.gg[j + (i__ << 1) - 3] = hgg[j + (i__ << 1) - 3];
/* labelL6: */
	}
    }
    l20_2.lex = true;
    l20_2.nex = 2;
    l20_2.fex = .97474658;
    for (i__ = 1; i__ <= 3; i__ += 2) {
	l20_2.xex[i__ - 1] = 1.;
/* L7: */
	l20_2.xex[i__] = 0.;
    }
    for (i__ = 5; i__ <= 8; ++i__) {
/* L8: */
	l20_2.xex[i__ - 1] = l20_2.xex[9 - i__ - 1];
    }
    return 0;
labelL2:
    l6_1.fx = 2.;
    for (i__ = 1; i__ <= 2; ++i__) {
/* labelL9: */
	l6_1.fx -= std::exp(l2_3.x[i__ - 1] * -10. * std::exp(-l2_3.x[i__ + 1]));
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 2; ++i__) {
	l4_3.gf[i__ - 1] = std::exp(l2_3.x[i__ - 1] * -10. * std::exp(-l2_3.x[i__ + 1]) 
		- l2_3.x[i__ + 1]) * 10.;
/* labelL10: */
	l4_3.gf[i__ + 1] = -l2_3.x[i__ - 1] * l4_3.gf[i__ - 1];
    }
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_3.x[0] + l2_3.x[1] - 1.;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = l2_3.x[2] + l2_3.x[3] - 1.;
    }
labelL5:
    return 0;
} /* tp265_ */


/* Subroutine */ int tp266_(int *mode)
{
    /* Initialized data */

    static Real a[10] = { .0426149,.0352053,.0878058,.0330812,.0580924,
	    .649704,.344144,-.627443,.001828,-.224783 };
    static Real d__[10] = { 2.34659,2.84048,1.13888,3.02286,1.72139,
	    .153917,.290577,-.159378,54.691,-.444873 };
    static Real c__[50]	/* was [10][5] */ = { -.564255,.535493,
	    .586387,.608734,.774227,-.435033,.759468,-.152448,-.821772,
	    .819831,.392417,.658799,.289826,.984915,.325421,-.688583,-.627795,
	    -.546437,-.53412,-.910632,-.404979,-.0636666,.854402,.375699,
	    -.151719,.222278,.0403142,.484134,-.798498,-.480344,.927589,
	    -.681091,.789312,.239547,.448051,-.524653,.724666,.353951,
	    -.658572,-.871758,-.0735083,-.869487,.949721,.463136,.149926,
	    .413248,-.0182537,.887866,.662362,-.978666 };
    static Real b[25]	/* was [5][5] */ = { .354033,-.0230349,
	    -.211938,-.0554288,.220429,-.0230349,.29135,-.00180333,-.111141,
	    .0485461,-.211938,-.00180333,-.815808,-.133538,-.38067,-.0554288,
	    -.111141,-.133538,.389198,-.131586,.220429,.0485461,-.38067,
	    -.131586,.534706 };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int ninl;
    Real tp266a_(Real*, Real*, Real*, Real*, int*, Real*);
    static int i__, k, l;
    static Real hf;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_2.n = 5;
    l15_1.lsum = 10;
    l1_2.nili = 0;
    ninl = 0;
    l1_2.neli = 0;
    l1_2.nenl = 0;
    for (i__ = 1; i__ <= 5; ++i__) {
	l2_4.x[i__ - 1] = .1;
	l11_4.lxl[i__ - 1] = false;
	l12_4.lxu[i__ - 1] = false;
/* labelL6: */
	l20_7.xex[i__ - 1] = 0.;
    }
    l20_7.nex = 1;
    l20_7.lex = true;
/*      FEX=0.0 */
/*      DO I=1,10 */
/*         FEX=FEX+A(I)**2 */
/*      ENDDO */
    l20_7.fex = 1e4;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 10; ++i__) {
/* labelL20: */
	l16_4.f[i__ - 1] = tp266a_(&a[i__ - 1], b, c__, &d__[i__ - 1], &i__, 
		l2_4.x);
    }
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L7: */
/* Computing 2nd power */
	d__1 = l16_4.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    l6_1.fx *= 1e4;
    return 0;
labelL3:
    for (k = 1; k <= 5; ++k) {
	l4_4.gf[k - 1] = 0.;
	for (i__ = 1; i__ <= 10; ++i__) {
	    hf = 0.;
	    for (l = 1; l <= 5; ++l) {
/* labelL10: */
		hf += (b[k + l * 5 - 6] + b[l + k * 5 - 6]) * l2_4.x[l - 1];
	    }
	    l17_9.df[i__ + k * 10 - 11] = c__[i__ + k * 10 - 11] + d__[i__ - 
		    1] * .5 * hf;
	    l4_4.gf[k - 1] += l16_4.f[i__ - 1] * 2. * l17_9.df[i__ + k * 10 - 
		    11];
/* L8: */
	    l4_4.gf[k - 1] *= 1e4;
	}
    }
labelL4:
    return 0;
} /* tp266_ */

Real tp266a_(Real *ai, Real *b, Real *c__, Real *di, int *i__, Real *x)
{
    /* System generated locals */
    Real ret_val;

    /* Local variables */
    static int k, l;
    static Real hf;


/*  FUNKTION ZUR BERECHNUNG VON */
/*      F(I,X)=(A+C*X+0.5*X'*B*X*D)(I) */
/*  D.H. DAS I-TE ELEMENT DES REAL*8 VEKTORS F(10) FUER I=1,..,10 */

    /* Parameter adjustments */
    --x;
    c__ -= 11;
    b -= 6;

    /* Function Body */
    ret_val = *ai;
    for (k = 1; k <= 5; ++k) {
	hf = 0.;
	for (l = 1; l <= 5; ++l) {
/* labelL20: */
	    hf += b[k + l * 5] * x[l];
	}
	ret_val += x[k] * (c__[*i__ + k * 10] + *di * .5 * hf);
/* labelL10: */
    }
    return ret_val;
} /* tp266a_ */


/* Subroutine */ int tp267_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real h__;
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 5; ++i__) {
	l2_4.x[i__ - 1] = 2.;
	l11_4.lxl[i__ - 1] = false;
	l14_4.xu[i__ - 1] = (float)15.;
/* labelL6: */
	l12_4.lxu[i__ - 1] = true;
    }
    l11_4.lxl[0] = true;
    l11_4.lxl[1] = true;
    l11_4.lxl[4] = true;
    l13_4.xl[0] = (float)0.;
    l13_4.xl[1] = (float)0.;
    l13_4.xl[4] = (float)0.;
    l15_1.lsum = 11;
    l20_10.lex = true;
    l20_10.nex = 2;
    l20_10.fex = 0.;
    l20_10.xex[0] = 1.;
    l20_10.xex[1] = 10.;
    l20_10.xex[2] = 1.;
    l20_10.xex[3] = 5.;
    l20_10.xex[4] = 4.;
    l20_10.xex[5] = 10.;
    l20_10.xex[6] = 1.;
    l20_10.xex[7] = -5.;
    l20_10.xex[8] = -1.;
    l20_10.xex[9] = 4.;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 11; ++i__) {
	h__ = (Real) i__ * .1;
/*      HF(1)=X(3)*DEXP(-X(1)*H) */
/*      HF(2)=X(4)*DEXP(-X(2)*H) */
/*      HF(3)=3.D+0*DEXP(-X(5)*H) */
/*      HF(4)=-DEXP(-H)+5.D+0*DEXP(-1.D+1*H)-3.D+0*DEXP(-4.D+0*H) */
/*      F(I)=HF(1)-HF(2)+HF(3)+HF(4) */
/* labelL20: */
	l16_7.f[i__ - 1] = l2_4.x[2] * std::exp(-l2_4.x[0] * h__) - l2_4.x[3] * 
		exp(-l2_4.x[1] * h__) + std::exp(-l2_4.x[4] * h__) * 3. - (exp(
		-h__) - std::exp(h__ * -10.) * 5. + std::exp(h__ * -4.) * 3.);
    }
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 11; ++i__) {
/* L7: */
/* Computing 2nd power */
	d__1 = l16_7.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (j = 1; j <= 11; ++j) {
	h__ = (Real) j * .1;
	l17_10.df[j + 21] = std::exp(-l2_4.x[0] * h__);
	l17_10.df[j + 32] = -exp(-l2_4.x[1] * h__);
	l17_10.df[j - 1] = -h__ * l2_4.x[2] * l17_10.df[j + 21];
	l17_10.df[j + 10] = -h__ * l2_4.x[3] * l17_10.df[j + 32];
/* L8: */
	l17_10.df[j + 43] = -h__ * 3. * std::exp(-l2_4.x[4] * h__);
    }
    for (i__ = 1; i__ <= 5; ++i__) {
	l4_4.gf[i__ - 1] = 0.;
	for (j = 1; j <= 11; ++j) {
/* labelL13: */
	    l4_4.gf[i__ - 1] += l16_7.f[j - 1] * 2. * l17_10.df[j + i__ * 11 
		    - 12];
	}
    }
labelL4:
    return 0;
} /* tp267_ */


/* Subroutine */ int tp268_(int *mode)
{
    /* Initialized data */

    static int hgg[25]	/* was [5][5] */ = { -1,10,-8,8,-4,-1,10,1,-1,
	    -2,-1,-3,-2,2,3,-1,5,-5,5,-5,-1,4,3,-3,1 };
    static int dd[25]	/* was [5][5] */ = { 10197,-12454,-1013,1948,
	    329,-12454,20909,-1733,-4914,-186,-1013,-1733,1755,1089,-174,1948,
	    -4914,1089,1515,-22,329,-186,-174,-22,27 };
    static int ddvekt[5] = { -9170,17099,-2271,-4336,-43 };
    static int dvdv = 14463;

    static int ninl, i__, j;
    static Real hf;

/* KONSTANTE DATEN DER ZIELFUNKTION: */
/*         DD=D'*D */
/*         DDVEKT=DVEKT'*D */
/*         DVDV=DVEKT'*DVEKT=14463 */
/* MIT */
/*          -                      - */
/*          |  -74   80   18  -11  -4  | */
/*          |   14  -69   21   28   0  | */
/*     D= |   66  -72   -5    7   1  | */
/*          |  -12   66  -30  -23   3  | */
/*          |    3    8   -7   -4   1  | */
/*          |    4  -12    4    4   0  | */
/*          -                       - */

/*     DVEKT=( 51, -61, -56, 69, 10, -12 )' */
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_2.n = 5;
    l1_2.nili = 5;
    ninl = 0;
    l1_2.neli = 0;
    l1_2.nenl = 0;
    for (i__ = 1; i__ <= 5; ++i__) {
	l2_4.x[i__ - 1] = 1.;
	l11_4.lxl[i__ - 1] = false;
	l12_4.lxu[i__ - 1] = false;
	for (j = 1; j <= 5; ++j) {
/* labelL6: */
	    l5_29.gg[i__ + j * 5 - 6] = (Real) hgg[i__ + j * 5 - 6];
	}
    }
    l20_7.lex = true;
    l20_7.nex = 1;
    l20_7.fex = 0.;
    l20_7.xex[0] = 1.;
    l20_7.xex[1] = 2.;
    l20_7.xex[2] = -1.;
    l20_7.xex[3] = 3.;
    l20_7.xex[4] = -4.;
    return 0;
labelL2:
    l6_1.fx = (Real) dvdv;
    for (i__ = 1; i__ <= 5; ++i__) {
	hf = 0.;
	for (j = 1; j <= 5; ++j) {
/* L8: */
	    hf += (Real) dd[i__ + j * 5 - 6] * l2_4.x[j - 1];
	}
/* L7: */
	l6_1.fx += l2_4.x[i__ - 1] * (hf - (Real) ddvekt[i__ - 1] * 2.);
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 5; ++i__) {
	l4_4.gf[i__ - 1] = (Real) ddvekt[i__ - 1] * -2.;
	for (j = 1; j <= 5; ++j) {
/* labelL9: */
	    l4_4.gf[i__ - 1] += (Real) (dd[i__ + j * 5 - 6] + dd[j + 
		    i__ * 5 - 6]) * l2_4.x[j - 1];
	}
    }
    return 0;
labelL4:
    if (l9_5.index1[0]) {
	l3_4.g[0] = -l2_4.x[0] - l2_4.x[1] - l2_4.x[2] - l2_4.x[3] - l2_4.x[4]
		 + 5;
    }
    if (l9_5.index1[1]) {
	l3_4.g[1] = l2_4.x[0] * 10. + l2_4.x[1] * 10. - l2_4.x[2] * 3. + 
		l2_4.x[3] * 5. + l2_4.x[4] * 4. - 20.;
    }
    if (l9_5.index1[2]) {
	l3_4.g[2] = l2_4.x[0] * -8. + l2_4.x[1] - l2_4.x[2] * 2. - l2_4.x[3] *
		 5. + l2_4.x[4] * 3. + 40.;
    }
    if (l9_5.index1[3]) {
	l3_4.g[3] = l2_4.x[0] * 8. - l2_4.x[1] + l2_4.x[2] * 2. + l2_4.x[3] * 
		5. - l2_4.x[4] * 3. - 11.;
    }
    if (l9_5.index1[4]) {
	l3_4.g[4] = l2_4.x[0] * -4. - l2_4.x[1] * 2. + l2_4.x[2] * 3. - 
		l2_4.x[3] * 5. + l2_4.x[4] + 30.;
    }
labelL5:
    return 0;
} /* tp268_ */


/* Subroutine */ int tp269_(int *mode)
{
    /* Initialized data */

    static Real hgg[15]	/* was [3][5] */ = { 1.,0.,0.,3.,0.,1.,0.,1.,
	    0.,0.,1.,0.,0.,-2.,-1. };
    static Real hdf[20]	/* was [4][5] */ = { 1.,0.,0.,0.,-1.,1.,0.,0.,
	    0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1. };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l15_1.lsum = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 3;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 5; ++i__) {
	l2_4.x[i__ - 1] = 2.;
	l11_4.lxl[i__ - 1] = false;
	l12_4.lxu[i__ - 1] = false;
	for (j = 1; j <= 3; ++j) {
	    l5_9.gg[j + i__ * 3 - 4] = hgg[j + i__ * 3 - 4];
/* L7: */
	    l17_8.df[j + (i__ << 2) - 5] = hdf[j + (i__ << 2) - 5];
	}
/* labelL6: */
	l17_8.df[(i__ << 2) - 1] = hdf[(i__ << 2) - 1];
    }
    l20_7.lex = true;
    l20_7.nex = 1;
    l20_7.fex = 4.0930232558139537;
    l20_7.xex[0] = -.76744186046511631;
    l20_7.xex[1] = .2558139534883721;
    l20_7.xex[2] = .62790697674418605;
    l20_7.xex[3] = -.11627906976744186;
    l20_7.xex[4] = .2558139534883721;
    return 0;
labelL2:
    l16_3.f[0] = l2_4.x[0] - l2_4.x[1];
    l16_3.f[1] = l2_4.x[1] + l2_4.x[2] - 2.;
    l16_3.f[2] = l2_4.x[3] - 1.;
    l16_3.f[3] = l2_4.x[4] - 1.;
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 4; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__1 = l16_3.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    l4_4.gf[0] = (l2_4.x[0] - l2_4.x[1]) * 2.;
    l4_4.gf[1] = (l2_4.x[1] * 2. + l2_4.x[2] - l2_4.x[0] - 2.) * 2.;
    l4_4.gf[2] = (l2_4.x[1] + l2_4.x[2] - 2.) * 2.;
    l4_4.gf[3] = (l2_4.x[3] - 1.) * 2.;
    l4_4.gf[4] = (l2_4.x[4] - 1.) * 2.;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_4.x[0] + l2_4.x[1] * 3.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_4.x[2] + l2_4.x[3] - l2_4.x[4] * 2.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = l2_4.x[1] - l2_4.x[4];
    }
labelL5:
    return 0;
} /* tp269_ */

/* Subroutine */ int tp270_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 4; ++i__) {
	l2_4.x[i__ - 1] = (Real) i__ + .1;
	l20_7.xex[i__ - 1] = (Real) i__;
	l11_4.lxl[i__ - 1] = true;
	l12_4.lxu[i__ - 1] = false;
/* labelL6: */
	l13_4.xl[i__ - 1] = (Real) i__;
    }
    l11_4.lxl[4] = false;
    l12_4.lxu[4] = false;
    l2_4.x[4] = -1.;
    l20_7.xex[4] = 2.;
    l20_7.nex = 1;
    l20_7.lex = true;
    l20_7.fex = -1.;
    return 0;
labelL2:
/* Computing 4th power */
    d__1 = l2_4.x[4], d__1 *= d__1;
/* Computing 3rd power */
    d__2 = l2_4.x[4];
/* Computing 2nd power */
    d__3 = l2_4.x[4];
    l6_1.fx = l2_4.x[0] * l2_4.x[1] * l2_4.x[2] * l2_4.x[3] - l2_4.x[0] * 3. *
	     l2_4.x[1] * l2_4.x[3] - l2_4.x[0] * 4. * l2_4.x[1] * l2_4.x[2] + 
	    l2_4.x[0] * 12. * l2_4.x[1] - l2_4.x[1] * l2_4.x[2] * l2_4.x[3] + 
	    l2_4.x[1] * 3. * l2_4.x[3] + l2_4.x[1] * 4. * l2_4.x[2] - l2_4.x[
	    1] * 12. - l2_4.x[0] * 2. * l2_4.x[2] * l2_4.x[3] + l2_4.x[0] * 
	    6. * l2_4.x[3] + l2_4.x[0] * 8. * l2_4.x[2] - l2_4.x[0] * 24. + 
	    l2_4.x[2] * 2. * l2_4.x[3] - l2_4.x[3] * 6. - l2_4.x[2] * 8. + 
	    24. + d__1 * d__1 * 1.5 - d__2 * (d__2 * d__2) * 5.75 + d__3 * 
	    d__3 * 5.25;
    return 0;
labelL3:
    l4_4.gf[0] = l2_4.x[1] * l2_4.x[2] * l2_4.x[3] - l2_4.x[1] * 3. * l2_4.x[
	    3] - l2_4.x[1] * 4. * l2_4.x[2] + l2_4.x[1] * 12. - l2_4.x[2] * 
	    2. * l2_4.x[3] + l2_4.x[3] * 6. + l2_4.x[2] * 8. - 24.;
    l4_4.gf[1] = l2_4.x[0] * l2_4.x[2] * l2_4.x[3] - l2_4.x[0] * 3. * l2_4.x[
	    3] - l2_4.x[0] * 4. * l2_4.x[2] + l2_4.x[0] * 12. - l2_4.x[2] * 
	    l2_4.x[3] + l2_4.x[3] * 3. + l2_4.x[2] * 4. - 12.;
    l4_4.gf[2] = l2_4.x[0] * l2_4.x[1] * l2_4.x[3] - l2_4.x[0] * 4. * l2_4.x[
	    1] - l2_4.x[1] * l2_4.x[3] + l2_4.x[1] * 4. - l2_4.x[0] * 2. * 
	    l2_4.x[3] + l2_4.x[0] * 8. + l2_4.x[3] * 2. - 8.;
    l4_4.gf[3] = l2_4.x[0] * l2_4.x[1] * l2_4.x[2] - l2_4.x[0] * 3. * l2_4.x[
	    1] - l2_4.x[1] * l2_4.x[2] + l2_4.x[1] * 3. - l2_4.x[0] * 2. * 
	    l2_4.x[2] + l2_4.x[0] * 6. + l2_4.x[2] * 2. - 6.;
/* Computing 2nd power */
    d__1 = l2_4.x[4];
/* Computing 3rd power */
    d__2 = l2_4.x[4];
    l4_4.gf[4] = l2_4.x[4] * 10.5 - d__1 * d__1 * 17.25 + d__2 * (d__2 * d__2)
	     * 6.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_4.x[0];
/* Computing 2nd power */
	d__2 = l2_4.x[1];
/* Computing 2nd power */
	d__3 = l2_4.x[2];
/* Computing 2nd power */
	d__4 = l2_4.x[3];
/* Computing 2nd power */
	d__5 = l2_4.x[4];
	l3_1.g[0] = 34. - d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * 
		d__4 - d__5 * d__5;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    for (i__ = 1; i__ <= 5; ++i__) {
/* L8: */
	l5_30.gg[i__ - 1] = l2_4.x[i__ - 1] * -2.;
    }
L7:
    return 0;
} /* tp270_ */

/* Subroutine */ int tp271_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 6;
    l15_1.lsum = 6;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 6; ++i__) {
	l11_5.lxl[i__ - 1] = false;
	l12_5.lxu[i__ - 1] = false;
	l2_5.x[i__ - 1] = 0.;
/* labelL6: */
	l20_4.xex[i__ - 1] = 1.;
    }
    l20_4.lex = true;
    l20_4.fex = 0.;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 6; ++i__) {
/* labelL10: */
	l16_8.f[i__ - 1] = std::sqrt((Real) ((16 - i__) * 10)) * (l2_5.x[i__ 
		- 1] - 1.);
    }
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 6; ++i__) {
/* L7: */
/* Computing 2nd power */
	d__1 = l16_8.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 6; ++i__) {
	l4_5.gf[i__ - 1] = (Real) ((16 - i__) * 20) * (l2_5.x[i__ - 1] 
		- 1.);
	for (j = 1; j <= 6; ++j) {
/* labelL9: */
	    l17_11.df[i__ + j * 6 - 7] = 0.;
	}
/* L8: */
	l17_11.df[i__ + i__ * 6 - 7] = std::sqrt((Real) ((16 - i__) * 10));
    }
labelL4:
    return 0;
} /* tp271_ */

/* Subroutine */ int tp272_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real h__;
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 6;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 6; ++i__) {
	l2_5.x[i__ - 1] = 1.;
	l13_5.xl[i__ - 1] = (float)0.;
	l11_5.lxl[i__ - 1] = true;
/* labelL6: */
	l12_5.lxu[i__ - 1] = false;
    }
    l2_5.x[1] = 2.;
    l15_1.lsum = 13;
    l20_4.lex = true;
    l20_4.nex = 1;
    l20_4.fex = 0.;
    l20_4.xex[0] = 1.;
    l20_4.xex[1] = 10.;
    l20_4.xex[2] = 4.;
    l20_4.xex[3] = 1.;
    l20_4.xex[4] = 5.;
    l20_4.xex[5] = 3.;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 13; ++i__) {
	h__ = (Real) i__ * .1;
/* labelL20: */
	l16_9.f[i__ - 1] = l2_5.x[3] * std::exp(-l2_5.x[0] * h__) - l2_5.x[4] * 
		exp(-l2_5.x[1] * h__) + l2_5.x[5] * std::exp(-l2_5.x[2] * h__) - 
		exp(-h__) + std::exp(h__ * -10.) * 5. - std::exp(h__ * -4.) * 3.;
    }
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* L7: */
/* Computing 2nd power */
	d__1 = l16_9.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (j = 1; j <= 13; ++j) {
	h__ = (Real) j * .1;
	l17_12.df[j + 38] = std::exp(-l2_5.x[0] * h__);
	l17_12.df[j + 51] = -exp(-l2_5.x[1] * h__);
	l17_12.df[j + 64] = std::exp(-l2_5.x[2] * h__);
	for (i__ = 1; i__ <= 3; ++i__) {
/* L8: */
	    l17_12.df[j + i__ * 13 - 14] = -h__ * l2_5.x[i__ + 2] * l17_12.df[
		    j + (i__ + 3) * 13 - 14];
	}
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	l4_5.gf[i__ - 1] = 0.;
	for (j = 1; j <= 13; ++j) {
/* labelL13: */
	    l4_5.gf[i__ - 1] += l17_12.df[j + i__ * 13 - 14] * 2. * l16_9.f[j 
		    - 1];
	}
    }
labelL4:
    return 0;
} /* tp272_ */

/* Subroutine */ int tp273_(int *mode)
{
    static int ninl;
    Real tp273a_(Real *);
    static int i__;
    static Real hx;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_2.n = 6;
    l1_2.nili = 0;
    ninl = 0;
    l1_2.neli = 0;
    l1_2.nenl = 0;
    for (i__ = 1; i__ <= 6; ++i__) {
	l2_5.x[i__ - 1] = 0.;
	l20_4.xex[i__ - 1] = 1.;
	l11_5.lxl[i__ - 1] = false;
/* labelL6: */
	l12_5.lxu[i__ - 1] = false;
    }
    l20_4.lex = true;
    l20_4.nex = 1;
    l20_4.fex = 0.;
    return 0;
labelL2:
    hx = tp273a_(l2_5.x);
    l6_1.fx = hx * 10. * (hx + 1.);
    return 0;
labelL3:
    hx = tp273a_(l2_5.x);
    for (i__ = 1; i__ <= 6; ++i__) {
/* L7: */
	l4_5.gf[i__ - 1] = (16. - (Real) i__) * 20. * (l2_5.x[i__ - 1] 
		- 1.) * (hx * 2. + 1.);
    }
labelL4:
    return 0;
} /* tp273_ */

Real tp273a_(Real *x)
{
    /* System generated locals */
    Real ret_val, d__1;

    /* Local variables */
    static int i__;


/*  BERECHNUNG VON H(X) IN TP273 */
/*  HX=SUM( (16-I)*(X(I)-1)**2,I=1,..,6) */

    /* Parameter adjustments */
    --x;

    /* Function Body */
    ret_val = 0.;
    for (i__ = 1; i__ <= 6; ++i__) {
/* labelL10: */
/* Computing 2nd power */
	d__1 = x[i__] - 1.;
	ret_val += (16. - (Real) i__) * (d__1 * d__1);
    }
    return ret_val;
} /* tp273a_ */

/* Subroutine */ int tp274_0_(int n__, int *mode)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static Real a[36]	/* was [6][6] */;
    static int i__, j;

    switch(n__) {
	case 1: goto L_tp275;
	case 2: goto L_tp276;
	}

    l1_1.n = 2;
    goto labelL10;

L_tp275:
    l1_1.n = 4;
    goto labelL10;

L_tp276:
    l1_1.n = 6;
labelL10:
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l2_5.x[i__ - 1] = -4. / (Real) i__;
	i__2 = l1_1.n;
	for (j = 1; j <= i__2; ++j) {
/* labelL11: */
	    a[i__ + j * 6 - 7] = 1. / (Real) (i__ + j - 1);
	}
	l20_4.xex[i__ - 1] = 0.;
	l11_5.lxl[i__ - 1] = false;
/* labelL6: */
	l12_5.lxu[i__ - 1] = false;
    }
    l20_4.lex = true;
    l20_4.nex = 1;
    l20_4.fex = 0.;
    return 0;
labelL2:
    l6_1.fx = 0.;
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = l1_1.n;
	for (j = 1; j <= i__2; ++j) {
/* L7: */
	    l6_1.fx += a[i__ + j * 6 - 7] * l2_5.x[i__ - 1] * l2_5.x[j - 1];
	}
    }
    return 0;
labelL3:
    i__2 = l1_1.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	l4_5.gf[i__ - 1] = 0.;
	i__1 = l1_1.n;
	for (j = 1; j <= i__1; ++j) {
/* L8: */
	    l4_5.gf[i__ - 1] += l2_5.x[j - 1] * (a[i__ + j * 6 - 7] + a[j + 
		    i__ * 6 - 7]);
	}
    }
labelL4:
    return 0;
} /* tp274_ */

/* Subroutine */ int tp274_(int *mode)
{
    return tp274_0_(0, mode);
    }

/* Subroutine */ int tp275_(int *mode)
{
    return tp274_0_(1, mode);
    }

/* Subroutine */ int tp276_(int *mode)
{
    return tp274_0_(2, mode);
    }

/* Subroutine */ int tp277_0_(int n__, int *mode)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static Real h__;
    static int i__, j;

    switch(n__) {
	case 1: goto L_tp278;
	case 2: goto L_tp279;
	case 3: goto L_tp280;
	}

    l1_1.n = 4;
    goto labelL10;

L_tp278:
    l1_1.n = 6;
    goto labelL10;

L_tp279:
    l1_1.n = 8;
    goto labelL10;

L_tp280:
    l1_1.n = 10;
labelL10:
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.nili = l1_1.n;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l20_10.fex = 0.;
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l2_9.x[i__ - 1] = 0.;
	i__2 = l1_1.n;
	for (j = 1; j <= i__2; ++j) {
	    l20_10.fex += 1. / (Real) (i__ + j - 1);
/* L16: */
	    l5_31.gg[j + (i__ * l1_1.n) - (l1_1.n+1)] = 1. / (Real) (i__ + j - 1);
//printf("jac[%d] = %f \n",j + (i__ * l1_1.n) - (l1_1.n+1),l5_31.gg[j + (i__ * l1_1.n) - (l1_1.n+1)] );
	}
	l11_9.lxl[i__ - 1] = true;
	l12_9.lxu[i__ - 1] = false;
	l13_9.xl[i__ - 1] = 0.;
/* labelL6: */
	l20_10.xex[i__ - 1] = 1.;
    }
    l20_10.lex = true;
    l20_10.nex = 1;
    return 0;
labelL2:
    l6_1.fx = 0.;
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__ = 0.;
	i__2 = l1_1.n;
	for (j = 1; j <= i__2; ++j) {
/* L17: */
	    h__ += 1. / (Real) (i__ + j - 1);
	}
/* L7: */
	l6_1.fx += h__ * l2_9.x[i__ - 1];
    }
    return 0;
labelL3:
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__ = 0.;
	i__2 = l1_1.n;
	for (j = 1; j <= i__2; ++j) {
/* L18: */
	    h__ += 1. / (Real) (i__ + j - 1);
	}
/* L8: */
	l4_9.gf[i__ - 1] = h__;
    }
    return 0;
labelL4:
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__ = 0.;
	i__2 = l1_1.n;
	for (j = 1; j <= i__2; ++j) {
/* L19: */
	    h__ += (l2_9.x[j - 1] - 1.) / (Real) (i__ + j - 1);
	}
/* labelL9: */
	if (l9_10.index1[i__ - 1]) {
	    l3_9.g[i__ - 1] = h__;
	}
    }
labelL5:
    return 0;
} /* tp277_ */

/* Subroutine */ int tp277_(int *mode)
{
    int tp277_0_(int, int *);
    return tp277_0_(0, mode);
    }

/* Subroutine */ int tp278_(int *mode)
{
    int tp277_0_(int, int *);
    return tp277_0_(1, mode);
    }

/* Subroutine */ int tp279_(int *mode)
{
    int tp277_0_(int, int *);
    return tp277_0_(2, mode);
    }

/* Subroutine */ int tp280_(int *mode)
{
    int tp277_0_(int, int *);
    return tp277_0_(3, mode);
    }


/* Subroutine */ int tp281_(int *mode)
{
    /* System generated locals */
    int i__1;
    Real d__1;

    /* Local variables */
    static Real h__;
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 10;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 10; ++i__) {
	l2_9.x[i__ - 1] = 0.;
	l11_9.lxl[i__ - 1] = false;
	l12_9.lxu[i__ - 1] = false;
/* labelL6: */
	l20_10.xex[i__ - 1] = 1.;
    }
    l20_10.lex = true;
    l20_10.nex = 1;
    l20_10.fex = 0.;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L7: */
/* Computing 3rd power */
	i__1 = i__;
/* Computing 2nd power */
	d__1 = l2_9.x[i__ - 1] - 1.;
	l6_1.fx += (Real) (i__1 * (i__1 * i__1)) * (d__1 * d__1);
    }
    l6_1.fx = pow_dd(&l6_1.fx, &c_b74);
    return 0;
labelL3:
    h__ = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L8: */
/* Computing 3rd power */
	i__1 = i__;
/* Computing 2nd power */
	d__1 = l2_9.x[i__ - 1] - 1.;
	h__ += (Real) (i__1 * (i__1 * i__1)) * (d__1 * d__1);
    }
    for (i__ = 1; i__ <= 10; ++i__) {
/* labelL13: */
/* Computing 3rd power */
	i__1 = i__;
	l4_9.gf[i__ - 1] = (Real) (i__1 * (i__1 * i__1)) * 
		.66666666666666663 * (l2_9.x[i__ - 1] - 1.) * pow_dd(&h__, &
		c_b1920);
   if( !isFinite(l4_9.gf[i__ - 1])) l4_9.gf[i__ - 1] = 0.0; 
    }
labelL4:
    return 0;
} /* tp281_ */


/* Subroutine */ int tp282_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real h__;
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 10;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 10; ++i__) {
	l2_9.x[i__ - 1] = 0.;
	l11_9.lxl[i__ - 1] = false;
	l12_9.lxu[i__ - 1] = false;
/* labelL6: */
	l20_10.xex[i__ - 1] = 1.;
    }
    l2_9.x[0] = -1.2;
    l20_10.lex = true;
    l20_10.nex = 1;
    l20_10.fex = 0.;
    l15_1.lsum = 11;
    return 0;
labelL2:
    l16_7.f[9] = l2_9.x[0] - 1.;
    l16_7.f[10] = l2_9.x[9] - 1.;
    for (i__ = 1; i__ <= 9; ++i__) {
/* labelL20: */
/* Computing 2nd power */
	d__1 = l2_9.x[i__ - 1];
	l16_7.f[i__ - 1] = std::sqrt((10. - (Real) i__) * 10.) * (d__1 * 
		d__1 - l2_9.x[i__]);
    }
    if (*mode == 3) {
	goto labelL3;
    }
/* Computing 2nd power */
    d__1 = l16_7.f[9];
/* Computing 2nd power */
    d__2 = l16_7.f[10];
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    for (i__ = 1; i__ <= 9; ++i__) {
/* L7: */
/* Computing 2nd power */
	d__1 = l16_7.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 11; ++i__) {
	for (j = 1; j <= 10; ++j) {
/* L8: */
	    l17_13.df[i__ + j * 11 - 12] = 0.;
	}
    }
    for (i__ = 1; i__ <= 9; ++i__) {
	h__ = std::sqrt((10. - (Real) i__) * 10.);
	l17_13.df[i__ + i__ * 11 - 12] = h__ * 2. * l2_9.x[i__ - 1];
/* labelL13: */
	l17_13.df[i__ + (i__ + 1) * 11 - 12] = -h__;
    }
    l17_13.df[9] = 1.;
    l17_13.df[109] = 1.;
    for (i__ = 1; i__ <= 10; ++i__) {
	l4_9.gf[i__ - 1] = 0.;
	for (j = 1; j <= 11; ++j) {
/* L18: */
	    l4_9.gf[i__ - 1] += l17_13.df[j + i__ * 11 - 12] * 2. * l16_7.f[j 
		    - 1];
	}
    }
labelL4:
    return 0;
} /* tp282_ */


/* Subroutine */ int tp283_(int *mode)
{
    /* System generated locals */
    int i__1;
    Real d__1;

    /* Local variables */
    static Real h__;
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 10;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 10; ++i__) {
	l2_9.x[i__ - 1] = 0.;
	l11_9.lxl[i__ - 1] = false;
	l12_9.lxu[i__ - 1] = false;
/* labelL6: */
	l20_10.xex[i__ - 1] = 1.;
    }
    l20_10.lex = true;
    l20_10.nex = 1;
    l20_10.fex = 0.;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L7: */
/* Computing 3rd power */
	i__1 = i__;
/* Computing 2nd power */
	d__1 = l2_9.x[i__ - 1] - 1.;
	l6_1.fx += (Real) (i__1 * (i__1 * i__1)) * (d__1 * d__1);
    }
/* Computing 3rd power */
    d__1 = l6_1.fx;
    l6_1.fx = d__1 * (d__1 * d__1) * 1e-5;
    return 0;
labelL3:
    h__ = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L8: */
/* Computing 3rd power */
	i__1 = i__;
/* Computing 2nd power */
	d__1 = l2_9.x[i__ - 1] - 1.;
	h__ += (Real) (i__1 * (i__1 * i__1)) * (d__1 * d__1);
    }
    for (i__ = 1; i__ <= 10; ++i__) {
/* labelL13: */
/* Computing 3rd power */
	i__1 = i__;
/* Computing 2nd power */
	d__1 = h__;
	l4_9.gf[i__ - 1] = (Real) (i__1 * (i__1 * i__1)) * (float)6. * (
		l2_9.x[i__ - 1] - 1.) * (d__1 * d__1) * 1e-5;
    }
labelL4:
    return 0;
} /* tp283_ */


/* Subroutine */ int tp284_(int *mode)
{
    /* Initialized data */

    static int c__[15] = { 20,40,400,20,80,20,40,140,380,280,80,40,140,40,
	    120 };
    static int b[10] = { 385,470,560,565,645,430,485,455,390,460 };
    static int a[150]	/* was [10][15] */ = { 100,90,70,50,50,40,30,
	    20,10,5,100,100,50,0,10,0,60,30,70,10,10,10,0,0,70,50,30,40,10,
	    100,5,35,55,65,60,95,90,25,35,5,10,20,25,35,45,50,0,40,25,20,0,5,
	    100,100,45,35,30,25,65,5,0,0,40,35,0,10,5,15,0,10,25,35,50,60,35,
	    60,25,10,30,35,0,55,0,0,65,0,0,80,0,95,10,25,30,15,5,45,70,20,0,
	    70,55,20,60,0,75,15,20,30,25,20,5,0,10,75,100,20,25,30,0,10,45,40,
	    30,35,75,0,70,5,15,35,20,25,0,30,10,5,15,65,50,10,0,10,40,65,0,5,
	    15,20,55,30 };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__, j;
    static Real sum;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 15;
    l1_1.nili = 0;
    l1_1.ninl = 10;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l20_12.fex = 0.;
    for (i__ = 1; i__ <= 15; ++i__) {
	l2_11.x[i__ - 1] = 0.;
	l20_12.xex[i__ - 1] = 1.;
	l20_12.fex -= (Real) c__[i__ - 1];
	l12_11.lxu[i__ - 1] = false;
/* labelL6: */
	l11_11.lxl[i__ - 1] = false;
    }
    l20_12.lex = true;
    l20_12.nex = 1;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 15; ++i__) {
/* L7: */
	l6_1.fx -= (Real) c__[i__ - 1] * l2_11.x[i__ - 1];
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 15; ++i__) {
/* L8: */
	l4_11.gf[i__ - 1] = -((Real) c__[i__ - 1]);
    }
    return 0;
labelL4:
    for (i__ = 1; i__ <= 10; ++i__) {
	sum = 0.;
	for (j = 1; j <= 15; ++j) {
/* labelL9: */
/* Computing 2nd power */
	    d__1 = l2_11.x[j - 1];
	    sum += (Real) a[i__ + j * 10 - 11] * (d__1 * d__1);
	}
/* labelL10: */
	if (l9_10.index1[i__ - 1]) {
	    l3_9.g[i__ - 1] = (Real) b[i__ - 1] - sum;
	}
    }
    return 0;
labelL5:
    for (j = 1; j <= 10; ++j) {
	for (i__ = 1; i__ <= 15; ++i__) {
	    if (! l10_10.index2[j - 1]) {
		goto labelL12;
	    }
/* labelL11: */
	    l5_32.gg[(i__ - 1) * 10 + j - 1] = (Real) a[j + i__ * 10 - 
		    11] * -2. * l2_11.x[i__ - 1];
	}
labelL12:
	;
    }
    return 0;
} /* tp284_ */

/* Subroutine */ int tp285_(int *mode)
{
    /* Initialized data */

    static int c__[15] = { 486,640,758,776,477,707,175,619,627,614,475,
	    377,524,468,529 };
    static int b[10] = { 385,470,560,565,645,430,485,455,390,460 };
    static int a[150]	/* was [10][15] */ = { 100,90,70,50,50,40,30,
	    20,10,5,100,100,50,0,10,0,60,30,70,10,10,10,0,0,70,50,30,40,10,
	    100,5,35,55,65,60,95,90,25,35,5,10,20,25,35,45,50,0,40,25,20,0,5,
	    100,100,45,35,30,25,65,5,0,0,40,35,0,10,5,15,0,10,25,35,50,60,35,
	    60,25,10,30,35,0,55,0,0,65,0,0,80,0,95,10,25,30,15,5,45,70,20,0,
	    70,55,20,60,0,75,15,20,30,25,20,5,0,10,75,100,20,25,30,0,10,45,40,
	    30,35,75,0,70,5,15,35,20,25,0,30,10,5,15,65,50,10,0,10,40,65,0,5,
	    15,20,55,30 };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__, j;
    static Real sum;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 15;
    l1_1.nili = 0;
    l1_1.ninl = 10;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l20_12.fex = 0.;
    for (i__ = 1; i__ <= 15; ++i__) {
	l2_11.x[i__ - 1] = 0.;
	l20_12.xex[i__ - 1] = 1.;
	l20_12.fex -= (Real) c__[i__ - 1];
	l12_11.lxu[i__ - 1] = false;
/* labelL6: */
	l11_11.lxl[i__ - 1] = false;
    }
    l20_12.lex = true;
    l20_12.nex = 1;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 15; ++i__) {
/* L7: */
	l6_1.fx -= (Real) c__[i__ - 1] * l2_11.x[i__ - 1];
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 15; ++i__) {
/* L8: */
	l4_11.gf[i__ - 1] = -((Real) c__[i__ - 1]);
    }
    return 0;
labelL4:
    for (i__ = 1; i__ <= 10; ++i__) {
	sum = 0.;
	for (j = 1; j <= 15; ++j) {
/* labelL9: */
/* Computing 2nd power */
	    d__1 = l2_11.x[j - 1];
	    sum += (Real) a[i__ + j * 10 - 11] * (d__1 * d__1);
	}
/* labelL10: */
	if (l9_10.index1[i__ - 1]) {
	    l3_9.g[i__ - 1] = (Real) b[i__ - 1] - sum;
	}
    }
    return 0;
labelL5:
    for (j = 1; j <= 10; ++j) {
	for (i__ = 1; i__ <= 15; ++i__) {
	    if (! l10_10.index2[j - 1]) {
		goto labelL12;
	    }
/* labelL11: */
	    l5_32.gg[(i__ - 1) * 10 + j - 1] = (Real) a[j + i__ * 10 - 
		    11] * -2. * l2_11.x[i__ - 1];
	}
labelL12:
	;
    }
    return 0;
} /* tp285_ */

/* Subroutine */ int tp286_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL3;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 20;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 20; ++i__) {
	l2_13.x[i__ - 1] = -1.2;
	l20_14.xex[i__ - 1] = 1.;
	l12_13.lxu[i__ - 1] = false;
/* labelL6: */
	l11_13.lxl[i__ - 1] = false;
    }
    for (i__ = 11; i__ <= 20; ++i__) {
/* L7: */
	l2_13.x[i__ - 1] = 1.;
    }
    l20_14.lex = true;
    l20_14.nex = 1;
    l20_14.fex = 0.;
    l15_1.lsum = 20;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 20; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__1 = l16_10.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 10; ++i__) {
	l16_10.f[i__ - 1] = l2_13.x[i__ - 1] - 1.;
/* labelL12: */
/* Computing 2nd power */
	d__1 = l2_13.x[i__ - 1];
	l16_10.f[i__ + 9] = (d__1 * d__1 - l2_13.x[i__ + 9]) * 10.;
    }
    if (*mode == 2) {
	goto labelL2;
    }
    for (i__ = 1; i__ <= 10; ++i__) {
	for (j = 1; j <= 20; ++j) {
	    l17_14.df[i__ + j * 198 - 199] = 0.;
	    if (i__ == j) {
		l17_14.df[i__ + j * 198 - 199] = 1.;
	    }
	    l17_14.df[i__ + 10 + j * 198 - 199] = 0.;
	    if (i__ == j) {
		l17_14.df[i__ + 10 + j * 198 - 199] = l2_13.x[i__ - 1] * 20.;
	    }
/* labelL9: */
	    if (j == i__ + 10) {
		l17_14.df[i__ + 10 + j * 198 - 199] = -10.;
	    }
	}
    }
    for (j = 1; j <= 20; ++j) {
	l4_13.gf[j - 1] = 0.;
	for (i__ = 1; i__ <= 20; ++i__) {
/* labelL11: */
	    l4_13.gf[j - 1] += l16_10.f[i__ - 1] * 2. * l17_14.df[i__ + j * 
		    198 - 199];
	}
    }
labelL4:
    return 0;
} /* tp286_ */

/* Subroutine */ int tp287_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 20;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 5; ++i__) {
	l2_13.x[i__ - 1] = -3.;
	l2_13.x[i__ + 4] = -1.;
	l2_13.x[i__ + 9] = -3.;
/* labelL6: */
	l2_13.x[i__ + 14] = -1.;
    }
    for (i__ = 1; i__ <= 20; ++i__) {
	l12_13.lxu[i__ - 1] = false;
	l11_13.lxl[i__ - 1] = false;
/* L7: */
	l20_14.xex[i__ - 1] = 1.;
    }
    l20_14.lex = true;
    l20_14.nex = 1;
    l20_14.fex = 0.;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 5; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__2 = l2_13.x[i__ - 1];
/* Computing 2nd power */
	d__1 = d__2 * d__2 - l2_13.x[i__ + 4];
/* Computing 2nd power */
	d__3 = l2_13.x[i__ - 1] - 1.;
/* Computing 2nd power */
	d__5 = l2_13.x[i__ + 9];
/* Computing 2nd power */
	d__4 = d__5 * d__5 - l2_13.x[i__ + 14];
/* Computing 2nd power */
	d__6 = l2_13.x[i__ + 9] - 1.;
/* Computing 2nd power */
	d__7 = l2_13.x[i__ + 4] - 1.;
/* Computing 2nd power */
	d__8 = l2_13.x[i__ + 14] - 1.;
	l6_1.fx = l6_1.fx + d__1 * d__1 * 100. + d__3 * d__3 + d__4 * d__4 * 
		90. + d__6 * d__6 + (d__7 * d__7 + d__8 * d__8) * 10.1 + (
		l2_13.x[i__ + 4] - 1.) * 19.8 * (l2_13.x[i__ + 14] - 1.);
    }
    l6_1.fx *= 1e-5;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 5; ++i__) {
/* Computing 2nd power */
	d__1 = l2_13.x[i__ - 1];
	l4_13.gf[i__ - 1] = (d__1 * d__1 - l2_13.x[i__ + 4]) * 400. * l2_13.x[
		i__ - 1] + (l2_13.x[i__ - 1] - 1.) * 2.;
	l4_13.gf[i__ - 1] *= 1e-5;
/* Computing 2nd power */
	d__1 = l2_13.x[i__ - 1];
	l4_13.gf[i__ + 4] = (d__1 * d__1 - l2_13.x[i__ + 4]) * -200. + (
		l2_13.x[i__ + 4] - 1.) * 20.2 + (l2_13.x[i__ + 14] - 1.) * 
		19.8;
	l4_13.gf[i__ + 4] *= 1e-5;
/* Computing 2nd power */
	d__1 = l2_13.x[i__ + 9];
	l4_13.gf[i__ + 9] = l2_13.x[i__ + 9] * 360. * (d__1 * d__1 - l2_13.x[
		i__ + 14]) + (l2_13.x[i__ + 9] - 1.) * 2.;
	l4_13.gf[i__ + 9] *= 1e-5;
/* Computing 2nd power */
	d__1 = l2_13.x[i__ + 9];
	l4_13.gf[i__ + 14] = (d__1 * d__1 - l2_13.x[i__ + 14]) * -180. + (
		l2_13.x[i__ + 14] - 1.) * 20.2 + (l2_13.x[i__ + 4] - 1.) * 
		19.8;
/* labelL9: */
	l4_13.gf[i__ + 14] *= 1e-5;
    }
labelL4:
    return 0;
} /* tp287_ */

/* Subroutine */ int tp288_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL3;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 20;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 5; ++i__) {
	l2_13.x[i__ - 1] = 3.;
	l2_13.x[i__ + 4] = -1.;
	l2_13.x[i__ + 9] = 0.;
/* labelL6: */
	l2_13.x[i__ + 14] = 1.;
    }
    for (i__ = 1; i__ <= 20; ++i__) {
	l20_14.xex[i__ - 1] = 0.;
	l11_13.lxl[i__ - 1] = false;
/* L7: */
	l12_13.lxu[i__ - 1] = false;
    }
    l20_14.fex = 0.;
    l20_14.lex = true;
    l20_14.nex = 1;
    l15_1.lsum = 20;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 20; ++i__) {
/* labelL9: */
/* Computing 2nd power */
	d__1 = l16_10.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 5; ++i__) {
	l16_10.f[i__ - 1] = l2_13.x[i__ - 1] + l2_13.x[i__ + 4] * 10.;
	l16_10.f[i__ + 4] = std::sqrt(5.) * (l2_13.x[i__ + 9] - l2_13.x[i__ + 14]);
/* Computing 2nd power */
	d__1 = l2_13.x[i__ + 4] - l2_13.x[i__ + 9] * 2.;
	l16_10.f[i__ + 9] = d__1 * d__1;
/* L8: */
/* Computing 2nd power */
	d__1 = l2_13.x[i__ - 1] - l2_13.x[i__ + 14];
	l16_10.f[i__ + 14] = std::sqrt(10.) * (d__1 * d__1);
    }
    if (*mode == 2) {
	goto labelL2;
    }
    for (i__ = 1; i__ <= 5; ++i__) {
	for (j = 1; j <= 20; ++j) {
	    l17_14.df[i__ + j * 198 - 199] = 0.;
	    if (j == i__) {
		l17_14.df[i__ + j * 198 - 199] = 1.;
	    }
	    if (j == i__ + 5) {
		l17_14.df[i__ + j * 198 - 199] = 10.;
	    }
	    l17_14.df[i__ + 5 + j * 198 - 199] = 0.;
	    if (j == i__ + 10) {
		l17_14.df[i__ + 5 + j * 198 - 199] = std::sqrt(5.);
	    }
	    if (j == i__ + 15) {
		l17_14.df[i__ + 5 + j * 198 - 199] = -std::sqrt(5.);
	    }
	    l17_14.df[i__ + 10 + j * 198 - 199] = 0.;
	    if (j == i__ + 5) {
		l17_14.df[i__ + 10 + j * 198 - 199] = (l2_13.x[i__ + 4] - 
			l2_13.x[i__ + 9] * 2.) * 2.;
	    }
	    if (j == i__ + 10) {
		l17_14.df[i__ + 10 + j * 198 - 199] = (l2_13.x[i__ + 4] - 
			l2_13.x[i__ + 9] * 2.) * -4.;
	    }
	    l17_14.df[i__ + 15 + j * 198 - 199] = 0.;
	    if (j == i__) {
		l17_14.df[i__ + 15 + j * 198 - 199] = std::sqrt(10.) * 2. * (
			l2_13.x[i__ - 1] - l2_13.x[i__ + 14]);
	    }
/* labelL11: */
	    if (j == i__ + 15) {
		l17_14.df[i__ + 15 + j * 198 - 199] = -std::sqrt(10.) * 2. * (
			l2_13.x[i__ - 1] - l2_13.x[i__ + 14]);
	    }
	}
    }
    for (j = 1; j <= 20; ++j) {
	l4_13.gf[j - 1] = 0.;
	for (i__ = 1; i__ <= 20; ++i__) {
/* labelL10: */
	    l4_13.gf[j - 1] += l16_10.f[i__ - 1] * 2. * l17_14.df[i__ + j * 
		    198 - 199];
	}
    }
labelL4:
    return 0;
} /* tp288_ */

/* Subroutine */ int tp289_(int *mode)
{
    /* System generated locals */
    Real d__1;
    Real pow_di( Real*, int*);

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 30;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 30; ++i__) {
	l2_14.x[i__ - 1] = pow_di(&c_b2003, &i__) * ((Real) i__ / 30. + 
		1.);
	l20_15.xex[i__ - 1] = 0.;
	l12_14.lxu[i__ - 1] = false;
/* labelL6: */
	l11_14.lxl[i__ - 1] = false;
    }
    l20_15.lex = true;
    l20_15.nex = 1;
    l20_15.fex = 0.;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 30; ++i__) {
/* L7: */
/* Computing 2nd power */
	d__1 = l2_14.x[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    l6_1.fx = 1. - std::exp(-l6_1.fx / 60.);
    if (*mode == 2) {
	return 0;
    }
    for (i__ = 1; i__ <= 30; ++i__) {
/* L8: */
	l4_14.gf[i__ - 1] = (l6_1.fx - 1.) * (l2_14.x[i__ - 1] * -2. / 60.);
    }
labelL4:
    return 0;
} /* tp289_ */

/* Subroutine */ int tp290_0_(int n__, int *mode)
{
    /* System generated locals */
    int i__1;
    Real d__1;

    /* Local variables */
    static int i__, j;

    switch(n__) {
	case 1: goto L_tp291;
	case 2: goto L_tp292;
	case 3: goto L_tp293;
	}

    l1_1.n = 2;
    goto labelL10;

L_tp291:
    l1_1.n = 10;
    goto labelL10;

L_tp292:
    l1_1.n = 30;
    goto labelL10;

L_tp293:
    l1_1.n = 50;
labelL10:
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL3;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l2_15.x[i__ - 1] = 1.;
	l20_16.xex[i__ - 1] = 0.;
	l12_15.lxu[i__ - 1] = false;
/* labelL6: */
	l11_15.lxl[i__ - 1] = false;
    }
    l20_16.fex = 0.;
    l20_16.lex = true;
    l20_16.nex = 1;
    l15_1.lsum = 1;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l16_11.f[0];
    l6_1.fx = d__1 * d__1;
    return 0;
labelL3:
    l16_11.f[0] = 0.;
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L7: */
/* Computing 2nd power */
	d__1 = l2_15.x[i__ - 1];
	l16_11.f[0] += (Real) i__ * (d__1 * d__1);
    }
    if (*mode == 2) {
	goto labelL2;
    }
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L8: */
	l17_9.df[i__ - 1] = (Real) i__ * 2. * l2_15.x[i__ - 1];
    }
    i__1 = l1_1.n;
    for (j = 1; j <= i__1; ++j) {
/* labelL9: */
	l4_15.gf[j - 1] = l16_11.f[0] * 2. * l17_9.df[j - 1];
    }
labelL4:
    return 0;
} /* tp290_ */

/* Subroutine */ int tp290_(int *mode)
{
    return tp290_0_(0, mode);
    }

/* Subroutine */ int tp291_(int *mode)
{
    return tp290_0_(1, mode);
    }

/* Subroutine */ int tp292_(int *mode)
{
    return tp290_0_(2, mode);
    }

/* Subroutine */ int tp293_(int *mode)
{
    return tp290_0_(3, mode);
    }


/* Subroutine */ int tp294_0_(int n__, int *mode)
{
    /* System generated locals */
    int i__1, i__2;
    Real d__1, d__2;

    /* Local variables */
    static int i__, j, k;

    switch(n__) {
	case 1: goto L_tp295;
	case 2: goto L_tp296;
	case 3: goto L_tp297;
	case 4: goto L_tp298;
	case 5: goto L_tp299;
	}

    l1_1.n = 6;
    goto labelL20;

L_tp295:
    l1_1.n = 10;
    goto labelL20;

L_tp296:
    l1_1.n = 16;
    goto labelL20;

L_tp297:
    l1_1.n = 30;
    goto labelL20;

L_tp298:
    l1_1.n = 50;
    goto labelL20;

L_tp299:
    l1_1.n = 100;
labelL20:
    k = l1_1.n - 1;
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL3;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l2_16.x[i__ - 1] = -1.2;
	l20_17.xex[i__ - 1] = 1.;
	l12_16.lxu[i__ - 1] = false;
/* labelL6: */
	l11_16.lxl[i__ - 1] = false;
    }
    i__1 = l1_1.n / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L7: */
	l2_16.x[(i__ << 1) - 1] = 1.;
    }
    l20_17.fex = 0.;
    l20_17.lex = true;
    l20_17.nex = 1;
    l15_1.lsum = k << 1;
    return 0;
labelL2:
    l6_1.fx = 0.;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__1 = l16_12.f[i__ - 1];
/* Computing 2nd power */
	d__2 = l16_12.f[i__ + k - 1];
	l6_1.fx = l6_1.fx + d__1 * d__1 + d__2 * d__2;
    }
    l6_1.fx *= 1e-4;
    return 0;
labelL3:
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = l2_16.x[i__ - 1];
	l16_12.f[i__ - 1] = (l2_16.x[i__] - d__1 * d__1) * 10.;
/* labelL12: */
	l16_12.f[i__ + k - 1] = 1. - l2_16.x[i__ - 1];
    }
    if (*mode == 2) {
	goto labelL2;
    }
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = l1_1.n;
	for (j = 1; j <= i__2; ++j) {
	    l17_14.df[i__ + j * 198 - 199] = 0.;
	    if (j == i__) {
		l17_14.df[i__ + j * 198 - 199] = l2_16.x[i__ - 1] * -20.;
	    }
	    if (j == i__ + 1) {
		l17_14.df[i__ + j * 198 - 199] = 10.;
	    }
	    l17_14.df[i__ + k + j * 198 - 199] = 0.;
/* labelL9: */
	    if (j == i__) {
		l17_14.df[i__ + k + j * 198 - 199] = -1.;
	    }
	}
    }
    i__2 = l1_1.n;
    for (j = 1; j <= i__2; ++j) {
	l4_16.gf[j - 1] = 0.;
	i__1 = l15_1.lsum;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* labelL13: */
	    l4_16.gf[j - 1] += l16_12.f[i__ - 1] * 2. * l17_14.df[i__ + j * 
		    198 - 199];
	}
	l4_16.gf[j - 1] *= 1e-4;
/* labelL11: */
    }
labelL4:
    return 0;
} /* tp294_ */

/* Subroutine */ int tp294_(int *mode)
{
    return tp294_0_(0, mode);
    }

/* Subroutine */ int tp295_(int *mode)
{
    return tp294_0_(1, mode);
    }

/* Subroutine */ int tp296_(int *mode)
{
    return tp294_0_(2, mode);
    }

/* Subroutine */ int tp297_(int *mode)
{
    return tp294_0_(3, mode);
    }

/* Subroutine */ int tp298_(int *mode)
{
    return tp294_0_(4, mode);
    }

/* Subroutine */ int tp299_(int *mode)
{
    return tp294_0_(5, mode);
    }


/* Subroutine */ int tp300_0_(int n__, int *mode)
{
    /* System generated locals */
    int i__1;
    Real d__1;

    /* Local variables */
    static int i__;

    switch(n__) {
	case 1: goto L_tp301;
	case 2: goto L_tp302;
	}

    l1_1.n = 20;
    l20_17.fex = -20.;
    goto labelL10;

L_tp301:
    l1_1.n = 50;
    l20_17.fex = -50.;
    goto labelL10;

L_tp302:
    l1_1.n = 100;
    l20_17.fex = -100.;
labelL10:
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l2_16.x[i__ - 1] = 0.;
	l20_17.xex[i__ - 1] = (Real) (l1_1.n - i__) + 1.;
	l12_16.lxu[i__ - 1] = false;
/* labelL6: */
	l11_16.lxl[i__ - 1] = false;
    }
    l20_17.lex = true;
    l20_17.nex = 1;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_16.x[0];
    l6_1.fx = d__1 * d__1 - l2_16.x[0] * 2.;
    i__1 = l1_1.n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L7: */
/* Computing 2nd power */
	d__1 = l2_16.x[i__ - 1];
	l6_1.fx = l6_1.fx + d__1 * d__1 * 2. - l2_16.x[i__ - 2] * 2. * 
		l2_16.x[i__ - 1];
    }
    return 0;
labelL3:
    l4_16.gf[0] = l2_16.x[0] * 2. - l2_16.x[1] * 2. - 2.;
    i__1 = l1_1.n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L8: */
	l4_16.gf[i__ - 1] = l2_16.x[i__ - 1] * 4. - (l2_16.x[i__ - 2] + 
		l2_16.x[i__]) * 2.;
    }
    l4_16.gf[l1_1.n - 1] = l2_16.x[l1_1.n - 1] * 4. - l2_16.x[l1_1.n - 2] * 
	    2.;
labelL4:
    return 0;
} /* tp300_ */

/* Subroutine */ int tp300_(int *mode)
{
    return tp300_0_(0, mode);
    }

/* Subroutine */ int tp301_(int *mode)
{
    return tp300_0_(1, mode);
    }

/* Subroutine */ int tp302_(int *mode)
{
    return tp300_0_(2, mode);
    }


/* Subroutine */ int tp303_0_(int n__, int *mode)
{
    /* System generated locals */
    int i__1;
    Real d__1, d__2;

    /* Local variables */
    static int i__;
    static Real pom;

    switch(n__) {
	case 1: goto L_tp304;
	case 2: goto L_tp305;
	}

    l1_1.n = 20;
    goto labelL10;

L_tp304:
    l1_1.n = 50;
    goto labelL10;

L_tp305:
    l1_1.n = 100;
labelL10:
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l2_16.x[i__ - 1] = .1;
	l20_17.xex[i__ - 1] = 0.;
	l12_16.lxu[i__ - 1] = false;
/* labelL6: */
	l11_16.lxl[i__ - 1] = false;
    }
    l20_17.lex = true;
    l20_17.nex = 1;
    l20_17.fex = 0.;
    return 0;
labelL2:
    pom = 0.;
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L7: */
	pom += (Real) i__ * .5 * l2_16.x[i__ - 1];
    }
/* Computing 2nd power */
    d__1 = pom;
/* Computing 4th power */
    d__2 = pom, d__2 *= d__2;
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__1 = l2_16.x[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* labelL9: */
/* Computing 3rd power */
	d__1 = pom;
	l4_16.gf[i__ - 1] = l2_16.x[i__ - 1] * 2. + pom * (Real) i__ + (
		Real) i__ * 2. * (d__1 * (d__1 * d__1));
    }
labelL4:
    return 0;
} /* tp303_ */

/* Subroutine */ int tp303_(int *mode)
{
    return tp303_0_(0, mode);
    }

/* Subroutine */ int tp304_(int *mode)
{
    return tp303_0_(1, mode);
    }

/* Subroutine */ int tp305_(int *mode)
{
    return tp303_0_(2, mode);
    }


/* Subroutine */ int tp306_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real a, b;
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 1.;
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_6.lex = false;
    l20_6.nex = 2;
    l20_6.fex = -1.1036;
    l20_6.xex[0] = 0.;
    l20_6.xex[1] = 1.;
/*      XEX(3)=.0D+0 */
/*      XEX(4)=-.1D+1 */
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l6_1.fx = -exp(-l2_1.x[0] - l2_1.x[1]) * (d__1 * d__1 * 2. + d__2 * d__2 *
	     3.);
    return 0;
labelL3:
    a = std::exp(-l2_1.x[0] - l2_1.x[1]);
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    b = d__1 * d__1 * 2. + d__2 * d__2 * 3.;
    l4_1.gf[0] = a * (b - l2_1.x[0] * 4.);
    l4_1.gf[1] = a * (b - l2_1.x[1] * 6.);
labelL4:
    return 0;
} /* tp306_ */


/* Subroutine */ int tp307_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;
    static Real wi, xi1, xi2;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = (float).3;
    l2_1.x[1] = (float).4;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = true;
	l13_1.xl[i__ - 1] = (float)0.;
	l12_1.lxu[i__ - 1] = true;
	l14_20.xu[i__ - 1] = (float).26;
/* labelL6: */
	l20_1.xex[i__ - 1] = .25783;
    }
    l20_1.nex = 1;
    l20_1.lex = false;
/*      FEX=0.12436*0.001 */
    l20_1.fex = (float)0.;
    l15_1.lsum = 10;
    return 0;
labelL2:
    l6_1.fx = (float)0.;
    goto labelL9;
labelL3:
    l4_1.gf[0] = (float)0.;
    l4_1.gf[1] = (float)0.;
labelL9:
    for (i__ = 1; i__ <= 10; ++i__) {
	wi = (Real) i__;
	xi1 = wi * l2_1.x[0];
	if (xi1 > 20.) {
	    xi1 = (float)0.;
	}
	xi2 = wi * l2_1.x[1];
	if (xi2 > 20.) {
	    xi2 = (float)0.;
	}
	l16_4.f[i__ - 1] = wi * (float)2. + (float)2. - std::exp(xi1) - std::exp(xi2);
	if (*mode == 2) {
	    goto L8;
	}
	l17_8.df[i__ - 1] = -wi * std::exp(xi1);
	l17_8.df[i__ + 9] = -wi * std::exp(xi2);
	l4_1.gf[0] += l16_4.f[i__ - 1] * (float)2. * l17_8.df[i__ - 1];
	l4_1.gf[1] += l16_4.f[i__ - 1] * (float)2. * l17_8.df[i__ + 9];
	goto L7;
L8:
/* Computing 2nd power */
	d__1 = l16_4.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
	l6_1.fx *= (float).001;
L7:
	;
    }
    l4_1.gf[0] *= (float).001;
    l4_1.gf[1] *= (float).001;
labelL4:
    return 0;
} /* tp307_ */


/* Subroutine */ int tp308_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 3.;
    l2_1.x[1] = .1;
    l12_1.lxu[0] = false;
    l12_1.lxu[1] = false;
    l11_1.lxl[0] = false;
    l11_1.lxl[1] = false;
    l20_1.lex = false;
    l15_1.lsum = 3;
    l20_1.nex = 1;
    l20_1.fex = .77319906;
    l20_1.xex[0] = -.15543724;
    l20_1.xex[1] = .69456378;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l16_2.f[0] = d__1 * d__1 + d__2 * d__2 + l2_1.x[0] * l2_1.x[1];
    l16_2.f[1] = std::sin(l2_1.x[0]);
    l16_2.f[2] =std::cos(l2_1.x[1]);
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
/* labelL5: */
/* Computing 2nd power */
	d__1 = l16_2.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    l17_2.df[0] = l2_1.x[0] * 2. + l2_1.x[1];
    l17_2.df[3] = l2_1.x[1] * 2. + l2_1.x[0];
    l17_2.df[1] = std::sin(l2_1.x[0]) * 2. *std::cos(l2_1.x[0]);
    l17_2.df[4] = 0.;
    l17_2.df[2] = 0.;
    l17_2.df[5] =std::cos(l2_1.x[1]) * -2. * std::sin(l2_1.x[1]);
    l4_1.gf[0] = 0.;
    l4_1.gf[1] = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l4_1.gf[0] += l16_2.f[i__ - 1] * 2. * l17_2.df[i__ - 1];
/* labelL6: */
	l4_1.gf[1] += l16_2.f[i__ - 1] * 2. * l17_2.df[i__ + 2];
    }
labelL4:
    return 0;
} /* tp308_ */

/* Subroutine */ int tp309_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 0.;
	l11_1.lxl[i__ - 1] = false;
/* labelL5: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.xex[0] = 3.4826826;
    l20_1.xex[1] = 3.9;
    l20_1.fex = -3.9871708;
labelL4:
    return 0;
labelL2:
/* Computing 4th power */
    d__1 = l2_1.x[0], d__1 *= d__1;
/* Computing 3rd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__3 = l2_1.x[0];
/* Computing 2nd power */
    d__4 = l2_1.x[1] - 3.9;
    l6_1.fx = d__1 * d__1 * 1.41 - d__2 * (d__2 * d__2) * 12.76 + d__3 * d__3 
	    * 39.91 - l2_1.x[0] * 51.93 + 24.37 + d__4 * d__4;
    return 0;
labelL3:
/* Computing 3rd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[0];
    l4_1.gf[0] = d__1 * (d__1 * d__1) * 5.64 - d__2 * d__2 * 38.28 + l2_1.x[0]
	     * 79.82 - 51.93;
    l4_1.gf[1] = l2_1.x[1] * 2. - 7.8;
    return 0;
} /* tp309_ */

/* Subroutine */ int tp310_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static Real a, b, c__;
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -1.2;
    l2_1.x[1] = 1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
/* L7: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = false;
    l20_1.nex = -1;
    l20_1.fex = 0.;
    l20_1.xex[0] = -1.2;
    l20_1.xex[1] = 1.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] * l2_1.x[1];
/* Computing 2nd power */
    d__2 = 1. - l2_1.x[0];
/* Computing 5th power */
    d__4 = 1. - l2_1.x[0], d__5 = d__4, d__4 *= d__4;
/* Computing 2nd power */
    d__3 = 1. - l2_1.x[0] - l2_1.x[1] * (d__5 * (d__4 * d__4));
    l6_1.fx = d__1 * d__1 * (d__2 * d__2) * (d__3 * d__3);
    return 0;
labelL3:
    a = l2_1.x[0] * l2_1.x[1];
    b = 1. - l2_1.x[0];
/* Computing 5th power */
    d__1 = b, d__2 = d__1, d__1 *= d__1;
    c__ = b - l2_1.x[1] * (d__2 * (d__1 * d__1));
/* Computing 4th power */
    d__1 = b, d__1 *= d__1;
    l4_1.gf[0] = a * 2. * b * c__ * (l2_1.x[1] - 1. - l2_1.x[1] * 5. * (d__1 *
	     d__1));
/* Computing 5th power */
    d__1 = b, d__2 = d__1, d__1 *= d__1;
    l4_1.gf[1] = a * 2. * b * c__ * (l2_1.x[0] - d__2 * (d__1 * d__1));
labelL4:
    return 0;
} /* tp310_ */

/* Subroutine */ int tp311_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 1.;
	l11_1.lxl[i__ - 1] = false;
/* labelL5: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_6.lex = false;
    l20_6.nex = 2;
    l20_6.xex[0] = 3.;
    l20_6.xex[1] = 2.;
    l20_6.xex[2] = 3.58443;
    l20_6.xex[3] = -1.84813;
    l20_6.fex = 0.;
labelL4:
    return 0;
labelL2:
/* Computing 2nd power */
    d__2 = l2_1.x[0];
/* Computing 2nd power */
    d__1 = d__2 * d__2 + l2_1.x[1] - 11.;
/* Computing 2nd power */
    d__4 = l2_1.x[1];
/* Computing 2nd power */
    d__3 = l2_1.x[0] + d__4 * d__4 - 7.;
    l6_1.fx = d__1 * d__1 + d__3 * d__3;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l4_1.gf[0] = l2_1.x[0] * 4. * (d__1 * d__1 + l2_1.x[1] - 11.) + (l2_1.x[0]
	     + d__2 * d__2 - 7.) * 2.;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l4_1.gf[1] = (d__1 * d__1 + l2_1.x[1] - 11.) * 2. + l2_1.x[1] * 4. * (
	    l2_1.x[0] + d__2 * d__2 - 7.);
    return 0;
} /* tp311_ */

/* Subroutine */ int tp312_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real a, b;
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 1.;
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 0.;
    l20_1.xex[0] = -21.026652;
    l20_1.xex[1] = -36.760009;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    a = d__1 * d__1 + l2_1.x[1] * 12. - 1.;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    b = (d__1 * d__1 + d__2 * d__2) * 49. + l2_1.x[0] * 84. + l2_1.x[1] * 
	    2324. - 681.;
/* Computing 2nd power */
    d__1 = a;
/* Computing 2nd power */
    d__2 = b;
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    l6_1.fx *= 1e-4;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    a = d__1 * d__1 + l2_1.x[1] * 12. - 1.;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    b = (d__1 * d__1 + d__2 * d__2) * 49. + l2_1.x[0] * 84. + l2_1.x[1] * 
	    2324. - 681.;
    l4_1.gf[0] = (l2_1.x[0] * 2. * a + b * (l2_1.x[0] * 98. + 84.)) * 2. * 
	    1e-4;
    l4_1.gf[1] = (a * 12. + b * (l2_1.x[1] * 98. + 2324.)) * 2. * 1e-4;
labelL4:
    return 0;
} /* tp312_ */

/* Subroutine */ int tp313_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;
    static Real xh;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 0.;
    l2_1.x[1] = -1.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = false;
/* labelL5: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = .199786;
    l20_1.xex[0] = 3.;
    l20_1.xex[1] = 2.850214;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 3.;
    l6_1.fx = d__1 * d__1 * 1e-4 - (l2_1.x[1] - l2_1.x[0]) + std::exp((l2_1.x[1] - 
	    l2_1.x[0]) * 20.);
    return 0;
labelL3:
    xh = std::exp((l2_1.x[1] - l2_1.x[0]) * 20.) * 20.;
    l4_1.gf[0] = (l2_1.x[0] - 3.) * 2e-4 + 1. - xh;
    l4_1.gf[1] = xh - 1.;
labelL4:
    return 0;
} /* tp313_ */

/* Subroutine */ int tp314_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static Real a, b;
    static int i__;
    static Real g1, h1;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 2.;
	l12_1.lxu[i__ - 1] = false;
/* labelL6: */
	l11_1.lxl[i__ - 1] = false;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = .16904;
    l20_1.xex[0] = 1.789039;
    l20_1.xex[1] = 1.3740024;
    return 0;
labelL2:
    a = l2_1.x[0] - 2.;
    b = l2_1.x[1] - 1.;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    g1 = d__1 * d__1 / -4. - d__2 * d__2 + 1.;
    h1 = l2_1.x[0] - l2_1.x[1] * 2. + 1.;
/* Computing 2nd power */
    d__1 = a;
/* Computing 2nd power */
    d__2 = b;
/* Computing 2nd power */
    d__3 = h1;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + .04 / g1 + d__3 * d__3 / .2;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    g1 = d__1 * d__1 / -4. - d__2 * d__2 + 1.;
    h1 = l2_1.x[0] - l2_1.x[1] * 2. + 1.;
/* Computing 2nd power */
    d__1 = g1;
    l4_1.gf[0] = (l2_1.x[0] - 2. + l2_1.x[0] * .01 / (d__1 * d__1) + h1 * 5.) 
	    * 2.;
/* Computing 2nd power */
    d__1 = g1;
    l4_1.gf[1] = (l2_1.x[1] - 1. + l2_1.x[1] * .04 / (d__1 * d__1) - h1 * 10.)
	     * 2.;
labelL4:
    return 0;
} /* tp314_ */

/* Subroutine */ int tp315_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 1;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = -.1;
    l2_1.x[1] = -.9;
    l12_1.lxu[0] = false;
    l12_1.lxu[1] = false;
    l11_1.lxl[0] = false;
    l11_1.lxl[1] = false;
    l5_3.gg[0] = 1.;
    l5_3.gg[3] = -2.;
    l4_1.gf[0] = 0.;
    l4_1.gf[1] = -1.;
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = -.8;
    l20_1.xex[0] = .6;
    l20_1.xex[1] = .8;
    return 0;
labelL2:
    l6_1.fx = -l2_1.x[1];
labelL3:
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = 1. - l2_1.x[1] * 2. + l2_1.x[0];
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_3.g[1] = d__1 * d__1 + d__2 * d__2;
    }
    if (l9_4.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_3.g[2] = 1. - d__1 * d__1 - d__2 * d__2;
    }
    return 0;
labelL5:
    if (! l10_4.index2[1]) {
	goto labelL6;
    }
    l5_3.gg[1] = l2_1.x[0] * 2.;
    l5_3.gg[4] = l2_1.x[1] * 2.;
labelL6:
    if (! l10_4.index2[2]) {
	goto L7;
    }
    l5_3.gg[2] = l2_1.x[0] * -2.;
    l5_3.gg[5] = l2_1.x[1] * -2.;
L7:
    return 0;
} /* tp315_ */

/* Subroutine */ int tp316_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 0.;
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 334.31458;
    l20_1.xex[0] = 7.0710678;
    l20_1.xex[1] = -7.0710678;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 20.;
/* Computing 2nd power */
    d__2 = l2_1.x[1] + 20.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2. - 40.;
    l4_1.gf[1] = l2_1.x[1] * 2. + 40.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_1.g[0] = d__1 * d__1 * .01 + d__2 * d__2 * .01 - 1.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	return 0;
    }
    l5_1.gg[0] = l2_1.x[0] * .02;
    l5_1.gg[1] = l2_1.x[1] * .02;
    return 0;
} /* tp316_ */

/* Subroutine */ int tp317_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 0.;
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 372.46661;
    l20_1.xex[0] = 7.3519262;
    l20_1.xex[1] = -5.422866;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 20.;
/* Computing 2nd power */
    d__2 = l2_1.x[1] + 20.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2. - 40.;
    l4_1.gf[1] = l2_1.x[1] * 2. + 40.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_1.g[0] = d__1 * d__1 * .01 + d__2 * d__2 / 64. - 1.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	return 0;
    }
    l5_1.gg[0] = l2_1.x[0] * .02;
    l5_1.gg[1] = l2_1.x[1] * 2. / 64.;
    return 0;
} /* tp317_ */

/* Subroutine */ int tp318_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 0.;
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 412.75005;
    l20_1.xex[0] = 7.8091266;
    l20_1.xex[1] = -3.7478414;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 20.;
/* Computing 2nd power */
    d__2 = l2_1.x[1] + 20.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2. - 40.;
    l4_1.gf[1] = l2_1.x[1] * 2. + 40.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_1.g[0] = d__1 * d__1 * .01 + d__2 * d__2 / 36. - 1.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	return 0;
    }
    l5_1.gg[0] = l2_1.x[0] * .02;
    l5_1.gg[1] = l2_1.x[1] * 2. / 36.;
    return 0;
} /* tp318_ */

/* Subroutine */ int tp319_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 0.;
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 452.4044;
    l20_1.xex[0] = 8.4922857;
    l20_1.xex[1] = -2.1121017;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 20.;
/* Computing 2nd power */
    d__2 = l2_1.x[1] + 20.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2. - 40.;
    l4_1.gf[1] = l2_1.x[1] * 2. + 40.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_1.g[0] = d__1 * d__1 * .01 + d__2 * d__2 / 16. - 1.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	return 0;
    }
    l5_1.gg[0] = l2_1.x[0] * .02;
    l5_1.gg[1] = l2_1.x[1] * 2. / 16.;
    return 0;
} /* tp319_ */

/* Subroutine */ int tp320_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 0.;
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 485.53146;
    l20_1.xex[0] = 9.39525;
    l20_1.xex[1] = -.68459019;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 20.;
/* Computing 2nd power */
    d__2 = l2_1.x[1] + 20.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2. - 40.;
    l4_1.gf[1] = l2_1.x[1] * 2. + 40.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_1.g[0] = d__1 * d__1 * .01 + d__2 * d__2 / 4. - 1.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	return 0;
    }
    l5_1.gg[0] = l2_1.x[0] * .02;
    l5_1.gg[1] = l2_1.x[1] * .5;
    return 0;
} /* tp320_ */

/* Subroutine */ int tp321_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 0.;
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 496.11237;
    l20_1.xex[0] = 9.8160292;
    l20_1.xex[1] = -.19093377;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 20.;
/* Computing 2nd power */
    d__2 = l2_1.x[1] + 20.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2. - 40.;
    l4_1.gf[1] = l2_1.x[1] * 2. + 40.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_1.g[0] = d__1 * d__1 * .01 + d__2 * d__2 - 1.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	return 0;
    }
    l5_1.gg[0] = l2_1.x[0] * .02;
    l5_1.gg[1] = l2_1.x[1] * 2.;
    return 0;
} /* tp321_ */


/* Subroutine */ int tp322_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = (float)1e-4;
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 499.96001;
    l20_1.xex[0] = 9.9980018;
    l20_1.xex[1] = -.0019990011;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - (float)20.;
/* Computing 2nd power */
    d__2 = l2_1.x[1] + (float)20.;
    l6_1.fx = d__1 * d__1 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2. - 40.;
    l4_1.gf[1] = l2_1.x[1] * 2. + 40.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_1.g[0] = d__1 * d__1 * .01 + d__2 * d__2 * 100. - 1.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	return 0;
    }
    l5_1.gg[0] = l2_1.x[0] * .02;
    l5_1.gg[1] = l2_1.x[1] * 200.;
    return 0;
} /* tp322_ */

/* Subroutine */ int tp323_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 1;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 0.;
    l2_1.x[1] = 1.;
    l11_1.lxl[0] = true;
    l13_1.xl[0] = 0.;
    l11_1.lxl[1] = true;
    l13_1.xl[1] = 0.;
    l12_1.lxu[0] = false;
    l12_1.lxu[1] = false;
    l5_2.gg[0] = 1.;
    l5_2.gg[2] = -1.;
    l5_2.gg[3] = 1.;
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 3.7989446;
    l20_1.xex[0] = .55357378;
    l20_1.xex[1] = 1.3064439;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 - l2_1.x[0] * 4. + 4.;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2. - 4.;
    l4_1.gf[1] = l2_1.x[1] * 2.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_1.x[0] - l2_1.x[1] + 2.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_2.g[1] = d__1 * d__1 * -1. + l2_1.x[1] - 1.;
    }
    return 0;
labelL5:
    if (l10_3.index2[1]) {
	l5_2.gg[1] = l2_1.x[0] * -2.;
    }
    return 0;
} /* tp323_ */

/* Subroutine */ int tp324_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 2.;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l11_1.lxl[0] = true;
    l11_1.lxl[1] = false;
    l13_1.xl[0] = 2.;
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 5.;
    l20_1.xex[0] = 15.811389;
    l20_1.xex[1] = 1.5811387;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l6_1.fx = d__1 * d__1 * .01 + d__2 * d__2;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * .02;
    l4_1.gf[1] = l2_1.x[1] * 2.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_1.x[0] * l2_1.x[1] - 25.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_2.g[1] = d__1 * d__1 + d__2 * d__2 - 25.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_2.gg[0] = l2_1.x[1];
    l5_2.gg[2] = l2_1.x[0];
L7:
    if (! l10_3.index2[1]) {
	return 0;
    }
    l5_2.gg[1] = l2_1.x[0] * 2.;
    l5_2.gg[3] = l2_1.x[1] * 2.;
    return 0;
} /* tp324_ */

/* Subroutine */ int tp325_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 1;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    l2_1.x[0] = -3.;
    l2_1.x[1] = 0.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l5_3.gg[i__ * 3 - 3] = -1.;
	l11_1.lxl[i__ - 1] = false;
/* labelL6: */
	l12_1.lxu[i__ - 1] = false;
    }
    l5_3.gg[1] = -1.;
    l4_1.gf[1] = 1.;
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 3.7913415;
    l20_1.xex[0] = -2.3722813;
    l20_1.xex[1] = -1.8363772;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l6_1.fx = d__1 * d__1 + l2_1.x[1];
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2.;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = -(l2_1.x[0] + l2_1.x[1]) + 1.;
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[1];
	l3_3.g[1] = -(l2_1.x[0] + d__1 * d__1) + 1.;
    }
    if (l9_4.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
/* Computing 2nd power */
	d__2 = l2_1.x[1];
	l3_3.g[2] = d__1 * d__1 + d__2 * d__2 - 9.;
    }
    return 0;
labelL5:
    if (l10_4.index2[1]) {
	l5_3.gg[4] = l2_1.x[1] * -2.;
    }
    if (! l10_4.index2[2]) {
	return 0;
    }
    l5_3.gg[2] = l2_1.x[0] * 2.;
    l5_3.gg[5] = l2_1.x[1] * 2.;
    return 0;
} /* tp325_ */

/* Subroutine */ int tp326_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = 4.;
    l2_1.x[1] = 3.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_1.lxl[i__ - 1] = true;
	l12_1.lxu[i__ - 1] = true;
	l14_1.xu[i__ - 1] = (float)10.;
/* labelL6: */
	l13_1.xl[i__ - 1] = 0.;
    }
    l5_2.gg[2] = -4.;
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = -79.807821;
    l20_1.xex[0] = 5.2396091;
    l20_1.xex[1] = 3.7460378;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
    l6_1.fx = d__1 * d__1 + d__2 * d__2 - l2_1.x[0] * 16. - l2_1.x[1] * 10.;
    return 0;
labelL3:
    l4_1.gf[0] = l2_1.x[0] * 2. - 16.;
    l4_1.gf[1] = l2_1.x[1] * 2. - 10.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0];
	l3_2.g[0] = 11. - d__1 * d__1 + l2_1.x[0] * 6. - l2_1.x[1] * 4.;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = l2_1.x[0] * l2_1.x[1] - l2_1.x[1] * 3. - std::exp(l2_1.x[0] - 
		3.) + 1.;
    }
    return 0;
labelL5:
    if (l10_3.index2[0]) {
	l5_2.gg[0] = l2_1.x[0] * -2. + 6.;
    }
    if (! l10_3.index2[1]) {
	return 0;
    }
    l5_2.gg[1] = l2_1.x[1] - std::exp(l2_1.x[0] - 3.);
    l5_2.gg[3] = l2_1.x[0] - 3.;
    return 0;
} /* tp326_ */


/* Subroutine */ int tp327_(int *mode)
{
    /* Initialized data */

    static Real y[44] = { .49,.49,.48,.47,.48,.47,.46,.46,.45,.43,.45,
	    .43,.43,.44,.43,.43,.46,.45,.42,.42,.43,.41,.41,.4,.42,.4,.4,.41,
	    .4,.41,.41,.4,.4,.4,.38,.41,.4,.4,.41,.38,.4,.4,.39,.39 };
    static Real z__[44] = { 8.,8.,10.,10.,10.,10.,12.,12.,12.,12.,14.,
	    14.,14.,16.,16.,16.,18.,18.,20.,20.,20.,22.,22.,22.,24.,24.,24.,
	    26.,26.,26.,28.,28.,30.,30.,30.,32.,32.,34.,36.,36.,38.,38.,40.,
	    42. };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = .42;
    l2_1.x[1] = 5.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l12_1.lxu[i__ - 1] = false;
	l11_1.lxl[i__ - 1] = true;
/* labelL12: */
	l13_1.xl[i__ - 1] = .4;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = .030646306;
    l20_1.xex[0] = .42190424;
    l20_1.xex[1] = 5.0000526;
    l15_1.lsum = 44;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 44; ++i__) {
/* labelL6: */
	l16_13.f[i__ - 1] = y[i__ - 1] - l2_1.x[0] - (.49 - l2_1.x[0]) * std::exp(
		-l2_1.x[1] * (z__[i__ - 1] - 8.));
    }
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 44; ++i__) {
/* L7: */
/* Computing 2nd power */
	d__1 = l16_13.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    l4_1.gf[0] = 0.;
    l4_1.gf[1] = 0.;
    for (i__ = 1; i__ <= 44; ++i__) {
	l17_15.df[i__ - 1] = std::exp(-l2_1.x[1] * (z__[i__ - 1] - 8.)) - 1.;
	l17_15.df[i__ + 43] = (.49 - l2_1.x[0]) * std::exp(-l2_1.x[1] * (z__[i__ - 
		1] - 8.)) * (z__[i__ - 1] - 8.);
	l4_1.gf[0] += l17_15.df[i__ - 1] * l16_13.f[i__ - 1] * 2.;
/* L8: */
	l4_1.gf[1] += l17_15.df[i__ + 43] * l16_13.f[i__ - 1] * 2.;
    }
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = -.09 - l2_1.x[0] * l2_1.x[1] + l2_1.x[1] * .49;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	return 0;
    }
    l5_1.gg[0] = -l2_1.x[1];
    l5_1.gg[1] = .49 - l2_1.x[0];
    return 0;
} /* tp327_ */

/* Subroutine */ int tp328_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static Real a, b;
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = .5;
	l12_1.lxu[i__ - 1] = true;
	l14_1.xu[i__ - 1] = 3.;
	l11_1.lxl[i__ - 1] = true;
/* labelL6: */
	l13_1.xl[i__ - 1] = 1.;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 1.744152;
    l20_1.xex[0] = 1.743439;
    l20_1.xex[1] = 2.0297056;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_1.x[1];
/* Computing 2nd power */
    d__2 = l2_1.x[0];
    a = (d__1 * d__1 + 1.) / (d__2 * d__2);
/* Computing 2nd power */
    d__1 = l2_1.x[0] * l2_1.x[1];
/* Computing 4th power */
    d__2 = l2_1.x[0] * l2_1.x[1], d__2 *= d__2;
    b = (d__1 * d__1 + 100.) / (d__2 * d__2);
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l6_1.fx = (d__1 * d__1 + 12. + a + b) / 10.;
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[1];
/* Computing 3rd power */
    d__2 = l2_1.x[0];
    a = (d__1 * d__1 + 1.) / (d__2 * (d__2 * d__2));
/* Computing 3rd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
/* Computing 5th power */
    d__3 = l2_1.x[0], d__4 = d__3, d__3 *= d__3;
/* Computing 4th power */
    d__5 = l2_1.x[1], d__5 *= d__5;
    b = 1. / (d__1 * (d__1 * d__1) * (d__2 * d__2)) + 200. / (d__4 * (d__3 * 
	    d__3) * (d__5 * d__5));
    l4_1.gf[0] = (l2_1.x[0] - a - b) / 5.;
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 3rd power */
    d__2 = l2_1.x[1];
/* Computing 4th power */
    d__3 = l2_1.x[0], d__3 *= d__3;
/* Computing 5th power */
    d__4 = l2_1.x[1], d__5 = d__4, d__4 *= d__4;
    a = 1. / (d__1 * d__1 * (d__2 * (d__2 * d__2))) + 200. / (d__3 * d__3 * (
	    d__5 * (d__4 * d__4)));
/* Computing 2nd power */
    d__1 = l2_1.x[0];
    l4_1.gf[1] = (l2_1.x[1] / (d__1 * d__1) - a) / 5.;
labelL4:
    return 0;
} /* tp328_ */

/* Subroutine */ int tp329_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 3;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l11_1.lxl[0] = true;
    l11_1.lxl[1] = true;
    l12_1.lxu[0] = true;
    l12_1.lxu[1] = true;
    l2_1.x[0] = 14.35;
    l2_1.x[1] = 8.6;
    l13_1.xl[0] = 13.;
    l13_1.xl[1] = 0.;
    l14_1.xu[0] = 16.;
    l14_1.xu[1] = 15.;
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = -6961.8139;
    l20_1.xex[0] = 14.095;
    l20_1.xex[1] = .84296079;
    return 0;
labelL2:
/* Computing 3rd power */
    d__1 = l2_1.x[0] - 10.;
/* Computing 3rd power */
    d__2 = l2_1.x[1] - 20.;
    l6_1.fx = d__1 * (d__1 * d__1) + d__2 * (d__2 * d__2);
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0] - 10.;
    l4_1.gf[0] = d__1 * d__1 * 3.;
/* Computing 2nd power */
    d__1 = l2_1.x[1] - 20.;
    l4_1.gf[1] = d__1 * d__1 * 3.;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0] - 5.;
/* Computing 2nd power */
	d__2 = l2_1.x[1] - 5.;
	l3_3.g[0] = d__1 * d__1 + d__2 * d__2 - 100.;
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0] - 6.;
/* Computing 2nd power */
	d__2 = l2_1.x[1] - 5.;
	l3_3.g[1] = d__1 * d__1 + d__2 * d__2;
    }
    if (l9_4.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_1.x[0] - 6.;
/* Computing 2nd power */
	d__2 = l2_1.x[1] - 5.;
	l3_3.g[2] = 82.81 - d__1 * d__1 - d__2 * d__2;
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto labelL6;
    }
    l5_3.gg[0] = l2_1.x[0] * 2. - 10.;
    l5_3.gg[3] = l2_1.x[1] * 2. - 10.;
labelL6:
    if (! l10_4.index2[1]) {
	goto L7;
    }
    l5_3.gg[1] = l2_1.x[0] * 2. - 12.;
    l5_3.gg[4] = l2_1.x[1] * 2. - 10.;
L7:
    if (! l10_4.index2[2]) {
	return 0;
    }
    l5_3.gg[2] = l2_1.x[0] * -2. + 12.;
    l5_3.gg[5] = l2_1.x[1] * -2. + 10.;
    return 0;
} /* tp329_ */

/* Subroutine */ int tp330_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = 2.5;
	l11_1.lxl[i__ - 1] = true;
	l12_1.lxu[i__ - 1] = true;
	l13_1.xl[i__ - 1] = 0.;
/* labelL6: */
	l14_1.xu[i__ - 1] = 5.;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 1.6205833;
    l20_1.xex[0] = 1.2866773;
    l20_1.xex[1] = .53046181;
    return 0;
labelL2:
/* Computing 3rd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
/* Computing 3rd power */
    d__3 = l2_1.x[1];
    l6_1.fx = d__1 * (d__1 * d__1) * .044 / (d__2 * d__2) + 1. / l2_1.x[0] + 
	    l2_1.x[0] * .0592 / (d__3 * (d__3 * d__3));
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_1.x[0];
/* Computing 2nd power */
    d__2 = l2_1.x[1];
/* Computing 2nd power */
    d__3 = 1 / l2_1.x[0];
/* Computing 3rd power */
    d__4 = l2_1.x[1];
    l4_1.gf[0] = d__1 * d__1 * .132 / (d__2 * d__2) - d__3 * d__3 + .0592 / (
	    d__4 * (d__4 * d__4));
/* Computing 3rd power */
    d__1 = l2_1.x[0];
/* Computing 3rd power */
    d__2 = l2_1.x[1];
/* Computing 4th power */
    d__3 = l2_1.x[1], d__3 *= d__3;
    l4_1.gf[1] = d__1 * (d__1 * d__1) * -.088 / (d__2 * (d__2 * d__2)) - 
	    l2_1.x[0] * .1776 / (d__3 * d__3);
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 3rd power */
	d__1 = l2_1.x[1];
	l3_1.g[0] = 1. - d__1 * (d__1 * d__1) * 8.62 / l2_1.x[0];
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	return 0;
    }
/* Computing 3rd power */
    d__1 = l2_1.x[1];
/* Computing 2nd power */
    d__2 = l2_1.x[0];
    l5_1.gg[0] = d__1 * (d__1 * d__1) * 8.62 / (d__2 * d__2);
/* Computing 2nd power */
    d__1 = l2_1.x[1];
    l5_1.gg[1] = d__1 * d__1 * -25.86 / l2_1.x[0];
    return 0;
} /* tp330_ */


/* Subroutine */ int tp331_(int *mode)
{

    /* Local variables */
    static Real a, b, c__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 1;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_1.x[0] = (float).5;
    l2_1.x[1] = (float).1;
    l12_1.lxu[0] = true;
    l12_1.lxu[1] = true;
    l14_1.xu[0] = (float).7;
    l14_1.xu[1] = (float).2;
    l11_1.lxl[0] = true;
    l13_1.xl[0] = (float).3;
    l11_1.lxl[1] = true;
    l13_1.xl[1] = (float).1;
    l5_1.gg[0] = -1.;
    l5_1.gg[1] = -1.;
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 4.258;
    l20_1.xex[0] = .6175;
    l20_1.xex[1] = .1039;
    return 0;
labelL2:
    l6_1.fx = std::log(std::log(l2_1.x[1]) * (float)2. / std::log(l2_1.x[0] + l2_1.x[1])) / 
	    l2_1.x[0];
    return 0;
labelL3:
    a = l2_1.x[0] + l2_1.x[1];
    b = std::log(a);
    c__ = std::log(l2_1.x[1]) * (float)2.;
    l4_1.gf[0] = -1. / l2_1.x[0] * (std::log(c__ / b) / l2_1.x[0] + 1. / (b * a));
    l4_1.gf[1] = (b * 2. / l2_1.x[1] - c__ / a) / (c__ * b * l2_1.x[0]);
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_1.g[0] = (float)1. - l2_1.x[0] - l2_1.x[1];
    }
labelL5:
    return 0;
} /* tp331_ */


/* Subroutine */ int tp332_(int *mode)
{
    /* Initialized data */

    static Real pi = 3.141592653;

    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real pbig, a, b, c__;
    static int i__;
    static Real tr, pangle, pim, xxx, yyy;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL3;
    }
labelL1:
    l1_1.n = 2;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 2; ++i__) {
	l2_1.x[i__ - 1] = .75;
	l12_1.lxu[i__ - 1] = true;
	l14_1.xu[i__ - 1] = 1.5;
	l11_1.lxl[i__ - 1] = true;
/* labelL6: */
	l13_1.xl[i__ - 1] = 0.;
    }
    l20_1.lex = false;
    l20_1.nex = 1;
    l20_1.fex = 114.95015;
    l20_1.xex[0] = .91139872;
    l20_1.xex[1] = .029280207;
labelL3:
    return 0;
labelL2:
    pim = pi / 3.6;
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 100; ++i__) {
	tr = pi * (((Real) i__ - 1.) / 180. + .33333333333333331);
	a = std::log(tr);
	b = std::sin(tr);
	c__ =std::cos(tr);
	xxx = (a + l2_1.x[1]) * b + l2_1.x[0] * c__;
	yyy = (a + l2_1.x[1]) * c__ - l2_1.x[0] * b;
/* L7: */
/* Computing 2nd power */
	d__1 = xxx;
/* Computing 2nd power */
	d__2 = yyy;
	l6_1.fx += pim * (d__1 * d__1 + d__2 * d__2);
    }
    return 0;
labelL4:
    pbig = -360.;
    pim = 180. / pi;
    for (i__ = 1; i__ <= 100; ++i__) {
	tr = pi * (((Real) i__ - 1.) / 180. + .33333333333333331);
	a = 1. / tr - l2_1.x[0];
	b = std::log(tr) + l2_1.x[1];
	pangle = pim * std::atan((d__1 = a / b, std::abs(d__1)));
/* L8: */
	if (pangle > pbig) {
	    pbig = pangle;
	}
    }
    if (l9_3.index1[0]) {
	l3_2.g[0] = 30. - pbig;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = pbig + 30.;
    }
    return 0;
} /* tp332_ */


/* Subroutine */ int tp333_(int *mode)
{
    /* Initialized data */

    static Real a[8] = { 4.,5.75,7.5,24.,32.,48.,72.,96. };
    static Real y[8] = { 72.1,65.6,55.9,17.1,9.8,4.5,1.3,.6 };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_2.x[0] = 30.;
    l2_2.x[1] = .04;
    l2_2.x[2] = 3.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l12_2.lxu[i__ - 1] = false;
/* labelL6: */
	l11_2.lxl[i__ - 1] = false;
    }
    l11_2.lxl[1] = true;
    l13_2.xl[1] = (float)0.;
    l12_2.lxu[1] = true;
    l14_2.xu[1] = (float).07;
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.fex = (float)43.200000000000003;
    l20_3.xex[0] = 89.902;
    l20_3.xex[1] = .06699;
    l20_3.xex[2] = .47809;
    l15_1.lsum = 8;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 8; ++i__) {
/* L7: */
	l16_14.f[i__ - 1] = (y[i__ - 1] - l2_2.x[0] * std::exp(-l2_2.x[1] * a[i__ 
		- 1]) - l2_2.x[2]) / y[i__ - 1];
    }
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 8; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__1 = l16_14.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    l6_1.fx *= (float)1e3;
    return 0;
labelL3:
    l4_2.gf[0] = 0.;
    l4_2.gf[1] = 0.;
    l4_2.gf[2] = 0.;
    for (i__ = 1; i__ <= 8; ++i__) {
	l17_16.df[i__ - 1] = -exp(-l2_2.x[1] * a[i__ - 1]) / y[i__ - 1];
	l17_16.df[i__ + 7] = l2_2.x[0] * a[i__ - 1] * std::exp(-l2_2.x[1] * a[i__ 
		- 1]) / y[i__ - 1];
	l17_16.df[i__ + 15] = -1. / y[i__ - 1];
	l4_2.gf[0] += l17_16.df[i__ - 1] * l16_14.f[i__ - 1] * 2.;
	l4_2.gf[1] += l17_16.df[i__ + 7] * l16_14.f[i__ - 1] * 2.;
/* labelL9: */
	l4_2.gf[2] += l17_16.df[i__ + 15] * l16_14.f[i__ - 1] * 2.;
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	l4_2.gf[i__ - 1] *= (float)1e3;
    }
labelL4:
    return 0;
} /* tp333_ */


/* Subroutine */ int tp334_(int *mode)
{
    /* Initialized data */

    static Real y[15] = { .14,.18,.22,.25,.29,.32,.35,.39,.37,.58,.73,
	    .96,1.34,2.1,4.39 };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;
    static Real ui, vi, wi;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l15_1.lsum = 15;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 1.;
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = .082481481;
    l20_3.xex[1] = 1.1354102;
    l20_3.xex[2] = 2.3413942;
    l20_3.fex = .0082149184;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 15; ++i__) {
	ui = (Real) i__;
	vi = (Real) (16 - i__);
	wi = std::min(ui,vi);
/* L7: */
	l16_15.f[i__ - 1] = y[i__ - 1] - (l2_2.x[0] + i__ / (l2_2.x[1] * vi + 
		l2_2.x[2] * wi));
    }
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 15; ++i__) {
/* labelL10: */
/* Computing 2nd power */
	d__1 = l16_15.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 3; ++i__) {
/* L8: */
	l4_2.gf[i__ - 1] = 0.;
    }
    for (i__ = 1; i__ <= 15; ++i__) {
	ui = (Real) i__;
	vi = (Real) (16 - i__);
	wi = std::min(ui,vi);
	l17_17.df[i__ - 1] = -1.;
/* Computing 2nd power */
	d__1 = l2_2.x[1] * vi + l2_2.x[2] * wi;
	l17_17.df[i__ + 14] = ui * vi / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = l2_2.x[1] * vi + l2_2.x[2] * wi;
	l17_17.df[i__ + 29] = ui * wi / (d__1 * d__1);
	l4_2.gf[0] += l17_17.df[i__ - 1] * l16_15.f[i__ - 1] * 2.;
	l4_2.gf[1] += l17_17.df[i__ + 14] * l16_15.f[i__ - 1] * 2.;
/* labelL9: */
	l4_2.gf[2] += l17_17.df[i__ + 29] * l16_15.f[i__ - 1] * 2.;
    }
labelL4:
    return 0;
} /* tp334_ */

/* Subroutine */ int tp335_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 2;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 1.;
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l5_3.gg[4] = -1.;
    l5_3.gg[5] = 1.;
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = 2.0309475e-6;
    l20_3.xex[1] = .0044721349;
    l20_3.xex[2] = .0020000032;
    l20_3.fex = -.004472137;
    return 0;
labelL2:
    l6_1.fx = -(l2_2.x[0] * .001 + l2_2.x[1]);
    return 0;
labelL3:
    l4_2.gf[0] = -.001;
    l4_2.gf[1] = -1.;
    l4_2.gf[2] = 0.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[1];
	l3_2.g[0] = d__1 * d__1 * 1e3 + d__2 * d__2 * 100. - l2_2.x[2];
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[1];
	l3_2.g[1] = d__1 * d__1 * 100. + d__2 * d__2 * 400. + l2_2.x[2] - .01;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_3.gg[0] = l2_2.x[0] * 2e3;
    l5_3.gg[2] = l2_2.x[1] * 200.;
    if (! l10_3.index2[1]) {
	goto L7;
    }
    l5_3.gg[1] = l2_2.x[0] * 200.;
    l5_3.gg[3] = l2_2.x[1] * 800.;
L7:
    return 0;
} /* tp335_ */

/* Subroutine */ int tp336_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 1;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 0.;
	l12_2.lxu[i__ - 1] = false;
/* labelL6: */
	l11_2.lxl[i__ - 1] = false;
    }
    l5_3.gg[0] = 5.;
    l5_3.gg[2] = 5.;
    l5_3.gg[4] = -3.;
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = .53459441;
    l20_3.xex[1] = .53397092;
    l20_3.xex[2] = -.21905778;
    l20_3.fex = -.33789573;
    return 0;
labelL2:
    l6_1.fx = l2_2.x[0] * 7. - l2_2.x[1] * 6. + l2_2.x[2] * 4.;
    return 0;
labelL3:
    l4_2.gf[0] = 7.;
    l4_2.gf[1] = -6.;
    l4_2.gf[2] = 4.;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_2.x[0] * 5. + l2_2.x[1] * 5. - l2_2.x[2] * 3. - 6.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[1];
/* Computing 2nd power */
	d__3 = l2_2.x[2];
	l3_2.g[1] = d__1 * d__1 + d__2 * d__2 * 2. + d__3 * d__3 * 3. - 1.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    l5_3.gg[1] = l2_2.x[0] * 2.;
    l5_3.gg[3] = l2_2.x[1] * 4.;
    l5_3.gg[5] = l2_2.x[2] * 6.;
L8:
    return 0;
} /* tp336_ */

/* Subroutine */ int tp337_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
/* labelL6: */
	l2_2.x[i__ - 1] = 1.;
    }
    l12_2.lxu[0] = false;
    l11_2.lxl[0] = false;
    l12_2.lxu[1] = false;
    l11_2.lxl[1] = true;
    l12_2.lxu[2] = true;
    l11_2.lxl[2] = false;
    l13_2.xl[1] = 1.;
    l14_2.xu[2] = 1.;
    l5_5.gg[2] = 0.;
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = .57735194;
    l20_3.xex[1] = 1.7320458;
    l20_3.xex[2] = -2.0256839e-6;
    l20_3.fex = 6.;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1];
/* Computing 2nd power */
    d__3 = l2_2.x[2];
    l6_1.fx = d__1 * d__1 * 9. + d__2 * d__2 + d__3 * d__3 * 9.;
    return 0;
labelL3:
    l4_2.gf[0] = l2_2.x[0] * 18.;
    l4_2.gf[1] = l2_2.x[1] * 2.;
    l4_2.gf[2] = l2_2.x[2] * 18.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = l2_2.x[0] * l2_2.x[1] - 1.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_5.gg[0] = l2_2.x[1];
    l5_5.gg[1] = l2_2.x[0];
L7:
    return 0;
} /* tp337_ */

/* Subroutine */ int tp338_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 1;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 0.;
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l5_3.gg[0] = .5;
    l5_3.gg[2] = 1.;
    l5_3.gg[4] = 1.;
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = .36689438;
    l20_3.xex[1] = 2.2437202;
    l20_3.xex[2] = -1.4271674;
    l20_3.fex = -7.2056984;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[1];
/* Computing 2nd power */
    d__3 = l2_2.x[2];
    l6_1.fx = -(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    return 0;
labelL3:
    for (i__ = 1; i__ <= 3; ++i__) {
/* L7: */
	l4_2.gf[i__ - 1] = l2_2.x[i__ - 1] * -2.;
    }
    return 0;
labelL4:
    if (l9_3.index1[0]) {
	l3_2.g[0] = l2_2.x[0] * .5 + l2_2.x[1] + l2_2.x[2] - 1.;
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[1];
/* Computing 2nd power */
	d__3 = l2_2.x[2];
	l3_2.g[1] = d__1 * d__1 + d__2 * d__2 * .66666666666666663 + d__3 * 
		d__3 * .25 - 4.;
    }
    return 0;
labelL5:
    if (! l10_3.index2[1]) {
	goto L8;
    }
    l5_3.gg[1] = l2_2.x[0] * 2.;
    l5_3.gg[3] = l2_2.x[1] * 1.3333333333333333;
    l5_3.gg[5] = l2_2.x[2] * .5;
L8:
    return 0;
} /* tp338_ */

/* Subroutine */ int tp339_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 1.;
	l12_2.lxu[i__ - 1] = false;
	l11_2.lxl[i__ - 1] = true;
/* labelL6: */
	l13_2.xl[i__ - 1] = 0.;
    }
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = 2.3797626;
    l20_3.xex[1] = .31622787;
    l20_3.xex[2] = 1.9429359;
    l20_3.fex = 3.3616797;
    return 0;
labelL2:
    l6_1.fx = .2 / (l2_2.x[0] * l2_2.x[1] * l2_2.x[2]) + 4. / l2_2.x[0] + 3. /
	     l2_2.x[2];
    return 0;
labelL3:
/* Computing 2nd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[0];
    l4_2.gf[0] = -.2 / (l2_2.x[1] * l2_2.x[2] * (d__1 * d__1)) - 4. / (d__2 * 
	    d__2);
/* Computing 2nd power */
    d__1 = l2_2.x[1];
    l4_2.gf[1] = -.2 / (l2_2.x[0] * l2_2.x[2] * (d__1 * d__1));
/* Computing 2nd power */
    d__1 = l2_2.x[2];
/* Computing 2nd power */
    d__2 = l2_2.x[2];
    l4_2.gf[2] = -.2 / (l2_2.x[0] * l2_2.x[1] * (d__1 * d__1)) - 3. / (d__2 * 
	    d__2);
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = 10. - l2_2.x[0] * 2. * l2_2.x[2] - l2_2.x[0] * l2_2.x[1];
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L8;
    }
    l5_5.gg[0] = l2_2.x[2] * -2. - l2_2.x[1];
    l5_5.gg[1] = -l2_2.x[0];
    l5_5.gg[2] = l2_2.x[0] * -2.;
L8:
    return 0;
} /* tp339_ */

/* Subroutine */ int tp340_(int *mode)
{
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 1;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 1.;
	l12_2.lxu[i__ - 1] = false;
/* labelL6: */
	l11_2.lxl[i__ - 1] = false;
    }
    l12_2.lxu[0] = true;
    l14_2.xu[0] = 1.;
    l5_5.gg[0] = -1.;
    l5_5.gg[1] = -2.;
    l5_5.gg[2] = -2.;
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = .60000009;
    l20_3.xex[1] = .29999998;
    l20_3.xex[2] = .29999998;
    l20_3.fex = -.054;
    return 0;
labelL2:
    l6_1.fx = -l2_2.x[0] * l2_2.x[1] * l2_2.x[2];
    return 0;
labelL3:
    l4_2.gf[0] = -l2_2.x[1] * l2_2.x[2];
    l4_2.gf[1] = -l2_2.x[0] * l2_2.x[2];
    l4_2.gf[2] = -l2_2.x[0] * l2_2.x[1];
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = 1.8 - l2_2.x[0] - l2_2.x[1] * 2. - l2_2.x[2] * 2.;
    }
labelL5:
    return 0;
} /* tp340_ */

/* Subroutine */ int tp341_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 1.;
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = false;
/* labelL6: */
	l13_2.xl[i__ - 1] = 0.;
    }
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = 4.;
    l20_3.xex[1] = 2.8284271;
    l20_3.xex[2] = 2.;
    l20_3.fex = -22.627417;
    return 0;
labelL2:
    l6_1.fx = -l2_2.x[0] * l2_2.x[1] * l2_2.x[2];
    return 0;
labelL3:
    l4_2.gf[0] = -l2_2.x[1] * l2_2.x[2];
    l4_2.gf[1] = -l2_2.x[0] * l2_2.x[2];
    l4_2.gf[2] = -l2_2.x[0] * l2_2.x[1];
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[1];
/* Computing 2nd power */
	d__3 = l2_2.x[2];
	l3_1.g[0] = d__1 * d__1 * -1. - d__2 * d__2 * 2. - d__3 * d__3 * 4. + 
		48.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_5.gg[0] = l2_2.x[0] * -2.;
    l5_5.gg[1] = l2_2.x[1] * -4.;
    l5_5.gg[2] = l2_2.x[2] * -8.;
L7:
    return 0;
} /* tp341_ */

/* Subroutine */ int tp342_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 1.;
	l12_2.lxu[i__ - 1] = false;
	l11_2.lxl[i__ - 1] = true;
/* labelL6: */
	l13_2.xl[i__ - 1] = 0.;
    }
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = 4.;
    l20_3.xex[1] = 2.8284271;
    l20_3.xex[2] = 2.;
    l20_3.fex = -22.627417;
    return 0;
labelL2:
    l6_1.fx = -l2_2.x[0] * l2_2.x[1] * l2_2.x[2];
    return 0;
labelL3:
    l4_2.gf[0] = -l2_2.x[1] * l2_2.x[2];
    l4_2.gf[1] = -l2_2.x[0] * l2_2.x[2];
    l4_2.gf[2] = -l2_2.x[0] * l2_2.x[1];
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[1];
/* Computing 2nd power */
	d__3 = l2_2.x[2];
	l3_1.g[0] = 48. - d__1 * d__1 - d__2 * d__2 * 2. - d__3 * d__3 * 4.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L8;
    }
    l5_5.gg[0] = l2_2.x[0] * -2.;
    l5_5.gg[1] = l2_2.x[1] * -4.;
    l5_5.gg[2] = l2_2.x[2] * -8.;
L8:
    return 0;
} /* tp342_ */

/* Subroutine */ int tp343_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l12_2.lxu[i__ - 1] = true;
	l11_2.lxl[i__ - 1] = true;
/* labelL6: */
	l13_2.xl[i__ - 1] = 0.;
    }
    l2_2.x[0] = 22.3;
    l2_2.x[1] = .5;
    l2_2.x[2] = 125.;
    l14_2.xu[0] = 36.;
    l14_2.xu[1] = 5.;
    l14_2.xu[2] = 125.;
    l5_3.gg[4] = 0.;
    l5_3.gg[3] = 0.;
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = 16.508383;
    l20_3.xex[1] = 2.4768216;
    l20_3.xex[2] = 123.99452;
    l20_3.fex = -5.6847825;
    return 0;
labelL2:
/* Computing 4th power */
    d__1 = l2_2.x[0], d__1 *= d__1;
/* Computing 2nd power */
    d__2 = l2_2.x[2];
    l6_1.fx = d__1 * d__1 * -.0201 * l2_2.x[1] * (d__2 * d__2) * 1e-7;
    return 0;
labelL3:
/* Computing 3rd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[2];
    l4_2.gf[0] = d__1 * (d__1 * d__1) * -.0804 * l2_2.x[1] * (d__2 * d__2) * 
	    1e-7;
/* Computing 4th power */
    d__1 = l2_2.x[0], d__1 *= d__1;
/* Computing 2nd power */
    d__2 = l2_2.x[2];
    l4_2.gf[1] = d__1 * d__1 * -.0201 * (d__2 * d__2) * 1e-7;
/* Computing 4th power */
    d__1 = l2_2.x[0], d__1 *= d__1;
    l4_2.gf[2] = d__1 * d__1 * -.0402 * l2_2.x[1] * l2_2.x[2] * 1e-7;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
	l3_2.g[0] = 675. - d__1 * d__1 * l2_2.x[1];
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[2];
	l3_2.g[1] = .419 - d__1 * d__1 * (d__2 * d__2) * 1e-7;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_3.gg[0] = l2_2.x[0] * -2. * l2_2.x[1];
/* Computing 2nd power */
    d__1 = l2_2.x[0];
    l5_3.gg[2] = -(d__1 * d__1);
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
/* Computing 2nd power */
    d__1 = l2_2.x[2];
    l5_3.gg[1] = l2_2.x[0] * -2e-7 * (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_2.x[0];
    l5_3.gg[5] = d__1 * d__1 * -2e-7 * l2_2.x[2];
L8:
    return 0;
} /* tp343_ */

/* Subroutine */ int tp344_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 2.;
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = false;
    }
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = 1.104859;
    l20_3.xex[1] = 1.1966742;
    l20_3.xex[2] = 1.5352623;
    l20_3.fex = .0325682;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0] - 1.;
/* Computing 2nd power */
    d__2 = l2_2.x[0] - l2_2.x[1];
/* Computing 4th power */
    d__3 = l2_2.x[1] - l2_2.x[2], d__3 *= d__3;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return 0;
labelL3:
    l4_2.gf[0] = (l2_2.x[0] - 1.) * 2. + (l2_2.x[0] - l2_2.x[1]) * 2.;
/* Computing 3rd power */
    d__1 = l2_2.x[1] - l2_2.x[2];
    l4_2.gf[1] = (l2_2.x[0] - l2_2.x[1]) * -2. + d__1 * (d__1 * d__1) * 4.;
/* Computing 3rd power */
    d__1 = l2_2.x[1] - l2_2.x[2];
    l4_2.gf[2] = d__1 * (d__1 * d__1) * -4.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[1];
/* Computing 4th power */
	d__2 = l2_2.x[2], d__2 *= d__2;
	l3_1.g[0] = l2_2.x[0] * (d__1 * d__1 + 1.) + d__2 * d__2 - 4. - std::sqrt(
		2.) * 3.;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
/* Computing 2nd power */
    d__1 = l2_2.x[1];
    l5_5.gg[0] = d__1 * d__1 + 1.;
    l5_5.gg[1] = l2_2.x[0] * 2. * l2_2.x[1];
/* Computing 3rd power */
    d__1 = l2_2.x[2];
    l5_5.gg[2] = d__1 * (d__1 * d__1) * 4.;
L7:
    return 0;
} /* tp344_ */

/* Subroutine */ int tp345_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 3; ++i__) {
	l2_2.x[i__ - 1] = 0.;
	l12_2.lxu[i__ - 1] = false;
/* labelL6: */
	l11_2.lxl[i__ - 1] = false;
    }
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = 1.1048584;
    l20_3.xex[1] = 1.1966752;
    l20_3.xex[2] = 1.5352622;
    l20_3.fex = .0325682;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_2.x[0] - 1.;
/* Computing 2nd power */
    d__2 = l2_2.x[0] - l2_2.x[1];
/* Computing 4th power */
    d__3 = l2_2.x[1] - l2_2.x[2], d__3 *= d__3;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    return 0;
labelL3:
    l4_2.gf[0] = (l2_2.x[0] - 1) * 2. + (l2_2.x[0] - l2_2.x[1]) * 2.;
/* Computing 3rd power */
    d__1 = l2_2.x[1] - l2_2.x[2];
    l4_2.gf[1] = (l2_2.x[0] - l2_2.x[1]) * -2. + d__1 * (d__1 * d__1) * 4.;
/* Computing 3rd power */
    d__1 = l2_2.x[1] - l2_2.x[2];
    l4_2.gf[2] = d__1 * (d__1 * d__1) * -4.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[1];
/* Computing 4th power */
	d__2 = l2_2.x[2], d__2 *= d__2;
	l3_1.g[0] = l2_2.x[0] * (d__1 * d__1 + 1.) + d__2 * d__2 - 4. - std::sqrt(
		18.);
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L8;
    }
/* Computing 2nd power */
    d__1 = l2_2.x[1];
    l5_5.gg[0] = d__1 * d__1 + 1.;
    l5_5.gg[1] = l2_2.x[0] * 2. * l2_2.x[1];
/* Computing 3rd power */
    d__1 = l2_2.x[2];
    l5_5.gg[2] = d__1 * (d__1 * d__1) * 4.;
L8:
    return 0;
} /* tp345_ */

/* Subroutine */ int tp346_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 3; ++i__) {
	l12_2.lxu[i__ - 1] = true;
	l11_2.lxl[i__ - 1] = true;
/* labelL6: */
	l13_2.xl[i__ - 1] = 0.;
    }
    l2_2.x[0] = 22.3;
    l2_2.x[1] = .5;
    l2_2.x[2] = 125.;
    l14_2.xu[0] = 36.;
    l14_2.xu[1] = 5.;
    l14_2.xu[2] = 125.;
    l5_3.gg[4] = 0.;
    l5_3.gg[3] = 0.;
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = 16.508383;
    l20_3.xex[1] = 2.4768216;
    l20_3.xex[2] = 123.99452;
    l20_3.fex = -5.6847825;
    return 0;
labelL2:
/* Computing 4th power */
    d__1 = l2_2.x[0], d__1 *= d__1;
/* Computing 2nd power */
    d__2 = l2_2.x[2];
    l6_1.fx = d__1 * d__1 * -.0201 * l2_2.x[1] * (d__2 * d__2) * 1e-7;
    return 0;
labelL3:
/* Computing 3rd power */
    d__1 = l2_2.x[0];
/* Computing 2nd power */
    d__2 = l2_2.x[2];
    l4_2.gf[0] = d__1 * (d__1 * d__1) * -.0804 * l2_2.x[1] * (d__2 * d__2) * 
	    1e-7;
/* Computing 4th power */
    d__1 = l2_2.x[0], d__1 *= d__1;
/* Computing 2nd power */
    d__2 = l2_2.x[2];
    l4_2.gf[1] = d__1 * d__1 * -.0201 * (d__2 * d__2) * 1e-7;
/* Computing 4th power */
    d__1 = l2_2.x[0], d__1 *= d__1;
    l4_2.gf[2] = d__1 * d__1 * -.0402 * l2_2.x[1] * l2_2.x[2] * 1e-7;
    return 0;
labelL4:
    if (l9_3.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
	l3_2.g[0] = 675. - d__1 * d__1 * l2_2.x[1];
    }
    if (l9_3.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_2.x[0];
/* Computing 2nd power */
	d__2 = l2_2.x[2];
	l3_2.g[1] = .419 - d__1 * d__1 * (d__2 * d__2) * 1e-7;
    }
    return 0;
labelL5:
    if (! l10_3.index2[0]) {
	goto L7;
    }
    l5_3.gg[0] = l2_2.x[0] * -2. * l2_2.x[1];
/* Computing 2nd power */
    d__1 = l2_2.x[0];
    l5_3.gg[2] = -(d__1 * d__1);
L7:
    if (! l10_3.index2[1]) {
	goto L8;
    }
/* Computing 2nd power */
    d__1 = l2_2.x[2];
    l5_3.gg[1] = l2_2.x[0] * -2e-7 * (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_2.x[0];
    l5_3.gg[5] = d__1 * d__1 * -2e-7 * l2_2.x[2];
L8:
    return 0;
} /* tp346_ */

/* Subroutine */ int tp347_(int *mode)
{
    /* Initialized data */

    static Real a[3] = { 8204.37,9008.72,9330.46 };

    /* Local variables */
    static int ninl;
    static Real h__[8];
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_3.n = 3;
    l1_3.nili = 0;
    ninl = 0;
    l1_3.neli = 1;
    l1_3.nenl = 0;
    l2_2.x[0] = .7;
    l2_2.x[1] = .2;
    l2_2.x[2] = .1;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = true;
	l12_2.lxu[i__ - 1] = true;
	l13_2.xl[i__ - 1] = 0.;
	l14_2.xu[i__ - 1] = 1.;
/* L7: */
	l5_5.gg[i__ - 1] = 1.;
    }
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = 0.;
    l20_3.xex[1] = 0.;
    l20_3.xex[2] = 1.;
    l20_3.fex = 17374.625;
    return 0;
labelL2:
    h__[0] = l2_2.x[0] + l2_2.x[1] + l2_2.x[2] + .03;
    h__[1] = l2_2.x[0] * .09 + l2_2.x[1] + l2_2.x[2] + .03;
    h__[2] = h__[0] * h__[1];
    h__[3] = l2_2.x[1] + l2_2.x[2] + .03;
    h__[4] = l2_2.x[1] * .07 + l2_2.x[2] + .03;
    h__[5] = h__[3] * h__[4];
    h__[6] = l2_2.x[2] + .03;
    h__[7] = l2_2.x[2] * .13 + .03;
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = a[0] * std::log(h__[0] / h__[1]) + a[1] * std::log(h__[3] / h__[4]) + a[2]
	     * std::log(h__[6] / h__[7]);
    return 0;
labelL3:
    l4_2.gf[0] = a[0] * (h__[1] - h__[0] * .09) / h__[2];
    l4_2.gf[1] = a[0] * (h__[1] - h__[0]) / h__[2] + a[1] * (h__[4] - h__[3] *
	     .07) / h__[5];
    l4_2.gf[2] = a[0] * (h__[1] - h__[0]) / h__[2] + a[1] * (h__[4] - h__[3]) 
	    / h__[5] + a[2] * (h__[7] - h__[6] * .13) / (h__[6] * h__[7]);
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = l2_2.x[0] + l2_2.x[1] + l2_2.x[2] - 1.;
    }
labelL5:
    return 0;
} /* tp347_ */


/* Subroutine */ int tp348_(int *mode)
{
    /* Initialized data */

    static Real rho = .0747;
    static Real xmu = .0443;
    static Real cp = .24;
    static Real pr = .709;
    static Real pi = 3.14159;
    static Real d__ = .525;
    static Real tin = 75.;
    static Real tsurf = 45.;
    static Real h__ = 13.13;
    static Real w = 3.166;
    static Real rhoc = 559.;
    static Real rhoa = 169.;

    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real etaf, delp, etas, xval;
    static int i__;
    static Real q, costf, costm, h1, xmdot, costt, ac, af, gi, at, re, 
	    ho, hef;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_2.lxl[i__ - 1] = false;
/* labelL6: */
	l12_2.lxu[i__ - 1] = true;
    }
    l11_2.lxl[1] = true;
    l13_2.xl[1] = 13.13;
    l14_2.xu[0] = .044;
    l14_2.xu[1] = 24.;
    l14_2.xu[2] = 600.;
    l2_2.x[0] = .1;
    l2_2.x[1] = 18.;
    l2_2.x[2] = 144.;
    l20_3.lex = false;
    l20_3.nex = 1;
    l20_3.xex[0] = .044;
    l20_3.xex[1] = 24.;
    l20_3.xex[2] = 85.607576;
    l20_3.fex = 36.97084;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = d__;
    af = l2_2.x[1] / l2_2.x[0] * 2. * (w * h__ - pi * 30. * (d__1 * d__1) / 
	    4.) / 144.;
    at = pi * 30. * d__ * l2_2.x[1] / 144.;
    ac = (h__ * l2_2.x[1] - d__ * 10. * l2_2.x[1] - l2_2.x[1] / l2_2.x[0] * 
	    .006 * h__) / 144.;
    gi = rho * l2_2.x[2] * (h__ * l2_2.x[1]) / (ac * 144.) * 60.;
    re = gi * 1.083 / (xmu * 12.);
    if (re < 1e-10) {
	re = 1e-10;
    }
    ho = gi * .195 * cp / (pow_dd(&pr, &c_b993) * pow_dd(&re, &c_b2368));
    xmdot = rho * l2_2.x[2] * h__ * l2_2.x[1] / 144. * 60.;
/* Computing 2nd power */
    d__1 = gi;
    delp = 1.833e-6 / rho * (d__1 * d__1) * 3. * (af / ac * pow_dd(&re, &
	    c_b308) + at * .1 / ac);
    if (ho < 1e-10) {
	ho = 1e-10;
    }
    xval = std::sqrt(ho) * .0732;
    etaf = std::tanh(xval) / xval;
    etas = 1. - af / (af + at) * (1. - etaf);
    hef = 1. - std::exp(-etas * ho * (af + at) / (xmdot * cp));
    q = hef * (tin - tsurf) * xmdot * cp;
    if (*mode == 4) {
	goto L7;
    }
    h1 = delp / rho * xmdot / 1.98e6;
    if (h1 < 1e-10) {
	h1 = 1e-10;
    }
    costm = std::sqrt(h1) / .0718 + 4.;
/* Computing 2nd power */
    d__1 = d__;
/* Computing 2nd power */
    d__2 = d__ - .036;
    costt = l2_2.x[1] * 30.300000000000001 * pi / 4. * (d__1 * d__1 - d__2 * 
	    d__2);
    costf = h__ * .47 * w * .006 * rhoa / 1728. * l2_2.x[1] / l2_2.x[0];
    costt = costt * rhoc / 1728.;
    l6_1.fx = costm + costt + costf;
labelL3:
    return 0;
labelL4:
    if (! l9_2.index1[0]) {
	goto labelL5;
    }
    goto labelL2;
L7:
    l3_1.g[0] = 6e3 - q;
labelL5:
    return 0;
} /* tp348_ */


/* Subroutine */ int tp349_(int *mode)
{
    /* Initialized data */

    static Real p2 = 10.;
    static Real c1f = .075;
    static Real c2f = .025;
    static Real h1 = 8e3;
    static Real h2 = 8e3;
    static Real e1 = 1e3;
    static Real e2 = 1e3;
    static Real cpp = 3.6938503;
    static Real p = 20.;

    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static Real area, heat, argu, wate, cost, vest, temp1, temp2, temp3;
    static int i__;
    static Real v, c0, c1, c2, xlmtd, c3, press, c4, c5, c6, c7, p1, 
	    a11, a12, a13, a22, a23, a31, a32, a33, ut, xk1, xk2, hea, dia, 
	    are, phi[9];

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL3;
    }
labelL1:
    l1_1.n = 3;
    l1_1.nili = 0;
    l1_1.ninl = 9;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_2.x[0] = 5e3;
    l2_2.x[1] = 200.;
    l2_2.x[2] = 100.;
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_2.lxl[i__ - 1] = true;
/* labelL6: */
	l12_2.lxu[i__ - 1] = true;
    }
    l11_2.lxl[2] = false;
    l12_2.lxu[2] = false;
    l13_2.xl[0] = 1e3;
    l13_2.xl[1] = 100.;
    l14_2.xu[0] = 8e3;
    l14_2.xu[1] = 500.;
    l20_3.lex = false;
    l20_3.fex = -4.1489499;
    l20_3.xex[0] = 7828.7954;
    l20_3.xex[1] = 188.81406;
    l20_3.xex[2] = 113.81406;
    return 0;
labelL2:
    p1 = 100.;
    for (i__ = 1; i__ <= 3; ++i__) {
/* L35: */
	if (l2_2.x[i__ - 1] < 1e-6) {
	    l2_2.x[i__ - 1] = 1e-6;
	}
    }
    xk1 = p1 * std::exp(-e1 / (l2_2.x[1] + 460.));
    xk2 = p2 * std::exp(-e2 / (l2_2.x[1] + 460.));
    v = p * l2_2.x[0] / (xk2 * (l2_2.x[0] * c2f - p));
    c1 = (l2_2.x[0] * c1f - p) / (l2_2.x[0] + v * xk1);
    ut = l2_2.x[1] * .0452 + 43.;
L39:
    argu = (l2_2.x[1] - l2_2.x[2] - 75.) / (l2_2.x[1] - 100.);
    if (argu == 0.) {
	goto L48;
    }
    xlmtd = (25. - l2_2.x[2]) / std::log((std::abs(argu)));
    heat = l2_2.x[0] * cpp * (100. - l2_2.x[1]) + xk1 * (l2_2.x[0] * c1f - p) 
	    * v * h1 / (l2_2.x[0] + v * xk1) + p * h2;
    area = heat / (ut * xlmtd);
    are = std::abs(area);
    hea = std::abs(heat);
    d__1 = v / 12.72;
    dia = pow_dd(&d__1, &c_b2380);
    if (l2_2.x[1] < 200.) {
	goto L40;
    }
/* Computing 3rd power */
    d__1 = l2_2.x[1];
    press = d__1 * (d__1 * d__1) * 3.3e-6 + 23.6;
    goto L41;
L48:
    l2_2.x[1] *= 1.0001;
    goto L39;
L40:
    press = 50.;
L41:
/* Computing 3rd power */
    d__1 = dia;
/* Computing 2nd power */
    d__2 = dia;
/* Computing 2nd power */
    d__3 = dia;
    wate = (d__1 * (d__1 * d__1) * .0909 + d__2 * d__2 * .482) * press + d__3 
	    * d__3 * 36.6 + dia * 160.5;
    c1 = pow_dd(&wate, &c_b2383) * 4.8;
    if (l2_2.x[1] < 200.) {
	goto L42;
    }
/* Computing 2nd power */
    d__1 = dia;
    c2 = (l2_2.x[1] * .0133 + 17.2) * (d__1 * d__1);
    goto L43;
L42:
    c2 = 0.;
L43:
    if (press < 150.) {
	goto L44;
    }
/* Computing 3rd power */
    d__1 = l2_2.x[1];
    c3 = pow_dd(&are, &c_b2387) * 270. * (d__1 * (d__1 * d__1) * 1.68e-7 + 
	    .962);
    goto L45;
L44:
    c3 = pow_dd(&are, &c_b2387) * 270.;
L45:
    c4 = dia * 140. + 1400.;
    d__1 = v * .05;
    c5 = pow_dd(&d__1, &c_b2390) * 875.;
/* Computing 3rd power */
    d__2 = l2_2.x[1];
    d__1 = d__2 * (d__2 * d__2) * 4.59e-11 + 6.95e-4 + l2_2.x[0];
    c6 = pow_dd(&d__1, &c_b2391) * 812.;
    if (l2_2.x[1] < 250.) {
	goto L46;
    }
    d__1 = hea * 298. / l2_2.x[2];
    c7 = pow_dd(&d__1, &c_b2391) * 1291.;
    goto L47;
L46:
    d__1 = hea * 298. / l2_2.x[2];
    c7 = pow_dd(&d__1, &c_b2391) * 812.;
L47:
    cost = c1 + c2 + c3 + c4 + c5 + c6 + c7;
    vest = cost * 5.;
/* Computing 3rd power */
    d__1 = l2_2.x[1];
    c0 = vest * .18 + 2.2e4 + v * 3.1 + (d__1 * (d__1 * d__1) * 4.59e-11 + 
	    6.95e-4) * l2_2.x[0] * 61.1 + heat * .00115 + heat * 6.92 + 
	    l2_2.x[0] * 574. * (c1f - c1) + 114800.;
    l6_1.fx = (6.88e5 - c0) / (vest * 2.) * -.001;
labelL3:
    return 0;
labelL4:
    p1 = 100.;
    for (i__ = 1; i__ <= 3; ++i__) {
/* L37: */
	if (l2_2.x[i__ - 1] < 1e-6) {
	    l2_2.x[i__ - 1] = 1e-6;
	}
    }
    xk1 = p1 * std::exp(-e1 / (l2_2.x[1] + 460.));
    xk2 = p2 * std::exp(-e2 / (l2_2.x[1] + 460.));
    v = p * l2_2.x[0] / (xk2 * (l2_2.x[0] * c2f - p));
    c1 = (l2_2.x[0] * c1f - p) / (l2_2.x[0] + v * xk1);
    ut = l2_2.x[1] * .0452 + 43.;
L36:
    argu = (l2_2.x[1] - l2_2.x[2] - 75.) / (l2_2.x[1] - 100.);
    if (argu == 0.) {
	goto L49;
    }
    xlmtd = (25. - l2_2.x[2]) / std::log((std::abs(argu)));
    heat = l2_2.x[0] * cpp * (100. - l2_2.x[1]) + xk1 * (l2_2.x[0] * c1f - p) 
	    * v * h1 / (l2_2.x[0] + v * xk1) + p * h2;
    area = heat / (ut * xlmtd);
    d__1 = v / 12.72;
    dia = pow_dd(&d__1, &c_b2380);
    if (l2_2.x[1] < 200.) {
	goto L50;
    }
/* Computing 3rd power */
    d__1 = l2_2.x[1];
    press = d__1 * (d__1 * d__1) * 3.3e-6 + 23.6;
    goto L51;
L49:
    l2_2.x[1] *= 1.0001;
    goto L36;
L50:
    press = 50.;
L51:
    phi[0] = dia - 1.25;
    phi[1] = 9.67 - dia;
    phi[2] = area - 50.;
    phi[3] = 4e3 - area;
    a11 = xk1 + l2_2.x[0] / v;
    a12 = xk2;
/* Computing 2nd power */
    d__1 = l2_2.x[1] + 460.;
/* Computing 2nd power */
    d__2 = l2_2.x[1] + 460.;
    a13 = (l2_2.x[0] * c1f - press) * xk1 * e1 / ((l2_2.x[0] + v * xk1) * (
	    d__1 * d__1)) + press * e2 / (v * (d__2 * d__2));
    a22 = xk2 + l2_2.x[0] / v;
/* Computing 2nd power */
    d__1 = l2_2.x[1] + 460.;
    a23 = press * e2 / (v * (d__1 * d__1));
    a31 = -h1 * xk1 / cpp;
    a32 = -h2 * xk2 / cpp;
/* Computing 2nd power */
    d__1 = l2_2.x[1] + 460.;
/* Computing 2nd power */
    d__2 = l2_2.x[1] + 460.;
    a33 = l2_2.x[0] / v + ut * area / (v * cpp) - (l2_2.x[0] * c1f - press) * 
	    xk1 * e1 * h1 / ((l2_2.x[0] + v * xk1) * cpp * (d__1 * d__1)) - 
	    press * e2 * h2 / (v * cpp * (d__2 * d__2));
    temp1 = a11 + a22 + a33;
    phi[4] = temp1;
    temp2 = a11 * a22 + a22 * a33 + a33 * a11 - a13 * a31 - a23 * a32;
    phi[5] = temp2;
    temp3 = a11 * a22 * a33 + a12 * a23 * a31 - a13 * a31 * a22 - a23 * a32 * 
	    a11;
    phi[6] = temp3;
    phi[7] = temp1 * temp2 - temp3;
    phi[8] = heat;
    for (i__ = 1; i__ <= 9; ++i__) {
/* L7: */
	if (l9_16.index1[i__ - 1]) {
	    l3_15.g[i__ - 1] = phi[i__ - 1];
	}
    }
    return 0;
} /* tp349_ */


/* Subroutine */ int tp350_(int *mode)
{
    /* Initialized data */

    static Real y[11] = { .1957,.1947,.1735,.16,.0844,.0627,.0456,.0342,
	    .0323,.0235,.0246 };
    static Real u[11] = { 4.,2.,1.,.5,.25,.167,.125,.1,.0833,.0714,
	    .0625 };

    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static Real h__[11];
    static int i__, j;

    for (i__ = 1; i__ <= 11; ++i__) {
/* labelL12: */
/* Computing 2nd power */
	d__1 = u[i__ - 1];
	h__[i__ - 1] = d__1 * d__1 + l2_3.x[2] * u[i__ - 1] + l2_3.x[3];
    }
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_3.x[0] = .25;
    l2_3.x[1] = .39;
    l2_3.x[2] = .415;
    l2_3.x[3] = .39;
    for (i__ = 1; i__ <= 4; ++i__) {
	l12_3.lxu[i__ - 1] = false;
/* labelL6: */
	l11_3.lxl[i__ - 1] = false;
    }
    l15_1.lsum = 11;
    l20_6.lex = false;
    l20_6.nex = 1;
    l20_6.fex = 3.0750561000000003;
    l20_6.xex[0] = .19280644;
    l20_6.xex[1] = .19126279;
    l20_6.xex[2] = .12305098;
    l20_6.xex[3] = .13605235;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 11; ++i__) {
/* labelL20: */
/* Computing 2nd power */
	d__1 = u[i__ - 1];
	l16_7.f[i__ - 1] = y[i__ - 1] - l2_3.x[0] / h__[i__ - 1] * (d__1 * 
		d__1 + l2_3.x[1] * u[i__ - 1]);
    }
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 11; ++i__) {
/* L7: */
/* Computing 2nd power */
	d__1 = l16_7.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    l6_1.fx *= 1e4;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 11; ++i__) {
/* Computing 2nd power */
	d__1 = u[i__ - 1];
	l17_18.df[i__ - 1] = (-(d__1 * d__1) - l2_3.x[1] * u[i__ - 1]) / h__[
		i__ - 1];
	l17_18.df[i__ + 10] = -l2_3.x[0] * u[i__ - 1] / h__[i__ - 1];
/* Computing 2nd power */
	d__1 = u[i__ - 1];
/* Computing 2nd power */
	d__2 = h__[i__ - 1];
	l17_18.df[i__ + 21] = l2_3.x[0] * u[i__ - 1] * (d__1 * d__1 + l2_3.x[
		1] * u[i__ - 1]) / (d__2 * d__2);
/* L8: */
/* Computing 2nd power */
	d__1 = u[i__ - 1];
/* Computing 2nd power */
	d__2 = h__[i__ - 1];
	l17_18.df[i__ + 32] = l2_3.x[0] * (d__1 * d__1 + l2_3.x[1] * u[i__ - 
		1]) / (d__2 * d__2);
    }
    for (j = 1; j <= 4; ++j) {
	l4_3.gf[j - 1] = 0.;
	for (i__ = 1; i__ <= 11; ++i__) {
/* labelL10: */
	    l4_3.gf[j - 1] += l16_7.f[i__ - 1] * 2. * l17_18.df[i__ + j * 11 
		    - 12];
	}
	l4_3.gf[j - 1] *= 1e4;
/* labelL9: */
    }
labelL4:
    return 0;
} /* tp350_ */


/* Subroutine */ int tp351_(int *mode)
{
    /* Initialized data */

    static Real a[7] = { 0.,4.28e-4,.001,.00161,.00209,.00348,.00525 };
    static Real b[7] = { 7.391,11.18,16.44,16.2,22.2,24.02,31.32 };

    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__, j;
    static Real xh1, xh2, xh3, xh4, xh5;

/* Computing 2nd power */
    d__1 = l2_3.x[0];
    xh1 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = l2_3.x[1];
    xh2 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = l2_3.x[2];
    xh3 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = l2_3.x[3];
    xh4 = d__1 * d__1;
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_3.x[0] = 2.7;
    l2_3.x[1] = 90.;
    l2_3.x[2] = 1500.;
    l2_3.x[3] = 10.;
    for (i__ = 1; i__ <= 4; ++i__) {
	l12_3.lxu[i__ - 1] = false;
/* labelL6: */
	l11_3.lxl[i__ - 1] = false;
    }
    l20_6.lex = false;
    l20_6.nex = 1;
    l20_6.fex = 318.58082;
    l20_6.xex[0] = 2.714;
    l20_6.xex[1] = 140.4;
    l20_6.xex[2] = 1707.;
    l20_6.xex[3] = 31.51;
    l15_1.lsum = 7;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 7; ++i__) {
/* labelL20: */
	l16_6.f[i__ - 1] = ((xh1 + a[i__ - 1] * xh2 + a[i__ - 1] * a[i__ - 1] 
		* xh3) / (a[i__ - 1] * xh4 + 1.) - b[i__ - 1]) / b[i__ - 1] * 
		100.;
    }
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 7; ++i__) {
/* L7: */
/* Computing 2nd power */
	d__1 = l16_6.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (j = 1; j <= 4; ++j) {
/* L8: */
	l4_3.gf[j - 1] = 0.;
    }
    for (i__ = 1; i__ <= 7; ++i__) {
	xh5 = (xh4 * a[i__ - 1] + 1.) * b[i__ - 1];
	l17_7.df[i__ - 1] = l2_3.x[0] * 200. / xh5;
	l17_7.df[i__ + 6] = l2_3.x[1] * 200. * a[i__ - 1] / xh5;
/* Computing 2nd power */
	d__1 = a[i__ - 1];
	l17_7.df[i__ + 13] = l2_3.x[2] * 200. * (d__1 * d__1) / xh5;
/* Computing 2nd power */
	d__1 = a[i__ - 1];
/* Computing 2nd power */
	d__2 = xh5;
	l17_7.df[i__ + 20] = l2_3.x[3] * -200. * a[i__ - 1] * b[i__ - 1] * (
		xh1 + xh2 * a[i__ - 1] + xh3 * (d__1 * d__1)) / (d__2 * d__2);
	for (j = 1; j <= 4; ++j) {
/* labelL10: */
	    l4_3.gf[j - 1] += l16_6.f[i__ - 1] * 2. * l17_7.df[i__ + j * 7 - 
		    8];
	}
    }
labelL4:
    return 0;
} /* tp351_ */

/* Subroutine */ int tp352_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__, j;
    static Real ti;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_3.x[0] = 25.;
    l2_3.x[1] = 5.;
    l2_3.x[2] = -5.;
    l2_3.x[3] = -1.;
    for (i__ = 1; i__ <= 4; ++i__) {
	l12_3.lxu[i__ - 1] = false;
/* labelL6: */
	l11_3.lxl[i__ - 1] = false;
    }
    l20_6.lex = false;
    l20_6.nex = 1;
    l20_6.fex = 903.23433;
    l20_6.xex[0] = -10.223574;
    l20_6.xex[1] = 11.908429;
    l20_6.xex[2] = -.45804134;
    l20_6.xex[3] = .58031996;
    l15_1.lsum = 40;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 20; ++i__) {
	ti = i__ * .2;
	l16_16.f[i__ - 1] = l2_3.x[0] + l2_3.x[1] * ti - std::exp(ti);
/* labelL20: */
	l16_16.f[i__ + 19] = l2_3.x[2] + l2_3.x[3] * std::sin(ti) -std::cos(ti);
    }
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 20; ++i__) {
/* L7: */
/* Computing 2nd power */
	d__1 = l16_16.f[i__ - 1];
/* Computing 2nd power */
	d__2 = l16_16.f[i__ + 19];
	l6_1.fx = l6_1.fx + d__1 * d__1 + d__2 * d__2;
    }
    return 0;
labelL3:
    for (j = 1; j <= 4; ++j) {
/* L8: */
	l4_3.gf[j - 1] = 0.;
    }
    for (i__ = 1; i__ <= 20; ++i__) {
	ti = i__ * .2;
	l17_19.df[i__ - 1] = 1.;
	l17_19.df[i__ + 19] = 0.;
	l17_19.df[i__ + 39] = ti;
	l17_19.df[i__ + 59] = 0.;
	l17_19.df[i__ + 79] = 0.;
	l17_19.df[i__ + 99] = 1.;
	l17_19.df[i__ + 119] = 0.;
	l17_19.df[i__ + 139] = std::sin(ti);
	for (j = 1; j <= 4; ++j) {
/* labelL10: */
	    l4_3.gf[j - 1] = l4_3.gf[j - 1] + l16_16.f[i__ - 1] * 2. * 
		    l17_19.df[i__ + j * 40 - 41] + l16_16.f[i__ + 19] * 2. * 
		    l17_19.df[i__ + 20 + j * 40 - 41];
	}
    }
labelL4:
    return 0;
} /* tp352_ */

/* Subroutine */ int tp353_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;
    static Real q;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 1;
    l1_1.ninl = 1;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    l2_3.x[0] = 0.;
    l2_3.x[1] = 0.;
    l2_3.x[2] = .4;
    l2_3.x[3] = .6;
    l5_7.gg[0] = 2.3;
    l5_7.gg[3] = 5.6;
    l5_7.gg[6] = 11.1;
    l5_7.gg[9] = 1.3;
    for (i__ = 1; i__ <= 4; ++i__) {
	l5_7.gg[i__ * 3 - 1] = 1.;
	l11_3.lxl[i__ - 1] = true;
	l12_3.lxu[i__ - 1] = false;
/* labelL6: */
	l13_3.xl[i__ - 1] = 0.;
    }
    l20_6.lex = false;
    l20_6.nex = 1;
    l20_6.fex = -39.933673;
    l20_6.xex[0] = 0.;
    l20_6.xex[1] = 0.;
    l20_6.xex[2] = .37755102;
    l20_6.xex[3] = .62244898;
    return 0;
labelL2:
    l6_1.fx = l2_3.x[0] * 24.55 + l2_3.x[1] * 26.75 + l2_3.x[2] * 39. + 
	    l2_3.x[3] * 40.5;
    l6_1.fx = -l6_1.fx;
    return 0;
labelL3:
    l4_3.gf[0] = -24.55;
    l4_3.gf[1] = -26.75;
    l4_3.gf[2] = -39.;
    l4_3.gf[3] = -40.5;
    return 0;
labelL4:
/* Computing 2nd power */
    d__1 = l2_3.x[0] * .53;
/* Computing 2nd power */
    d__2 = l2_3.x[1] * .44;
/* Computing 2nd power */
    d__3 = l2_3.x[2] * 4.5;
/* Computing 2nd power */
    d__4 = l2_3.x[3] * .79;
    q = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4;
    if (*mode == 5) {
	goto labelL5;
    }
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_3.x[0] * 2.3 + l2_3.x[1] * 5.6 + l2_3.x[2] * 11.1 + 
		l2_3.x[3] * 1.3 - 5.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_3.x[0] * 12. + l2_3.x[1] * 11.9 + l2_3.x[2] * 41.8 + 
		l2_3.x[3] * 52.1 - std::sqrt(q) * 1.645 - 12.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = l2_3.x[0] + l2_3.x[1] + l2_3.x[2] + l2_3.x[3] - 1.;
    }
    return 0;
labelL5:
    if (! l10_4.index2[1]) {
	goto L7;
    }
    l5_7.gg[1] = 12. - l2_3.x[0] * .46208050000000006 / std::sqrt(q);
    l5_7.gg[4] = 11.9 - l2_3.x[1] * .31847199999999998 / std::sqrt(q);
    l5_7.gg[7] = 41.8 - l2_3.x[2] * 33.311250000000001 / std::sqrt(q);
    l5_7.gg[10] = 52.1 - l2_3.x[3] * 1.0266445000000002 / std::sqrt(q);
L7:
    return 0;
} /* tp353_ */

/* Subroutine */ int tp354_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 1;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_3.x[0] = 3.;
    l2_3.x[1] = -1.;
    l2_3.x[2] = 0.;
    l2_3.x[3] = 1.;
    for (i__ = 1; i__ <= 4; ++i__) {
	l5_2.gg[i__ - 1] = 1.;
	l12_3.lxu[i__ - 1] = true;
	l11_3.lxl[i__ - 1] = false;
/* labelL6: */
	l14_3.xu[i__ - 1] = 20.;
    }
    l20_6.lex = false;
    l20_6.nex = 1;
    l20_6.fex = .11378386;
    l20_6.xex[0] = .50336521;
    l20_6.xex[1] = -.04556907;
    l20_6.xex[2] = .23580504;
    l20_6.xex[3] = .30639882;
    return 0;
labelL2:
/* Computing 2nd power */
    d__1 = l2_3.x[0] + l2_3.x[1] * 10.;
/* Computing 2nd power */
    d__2 = l2_3.x[2] - l2_3.x[3];
/* Computing 4th power */
    d__3 = l2_3.x[1] - l2_3.x[2] * 2., d__3 *= d__3;
/* Computing 4th power */
    d__4 = l2_3.x[0] - l2_3.x[3], d__4 *= d__4;
    l6_1.fx = d__1 * d__1 + d__2 * d__2 * 5. + d__3 * d__3 + d__4 * d__4 * 
	    10.;
    return 0;
labelL3:
/* Computing 3rd power */
    d__1 = l2_3.x[0] - l2_3.x[3];
    l4_3.gf[0] = l2_3.x[0] * 2. + l2_3.x[1] * 20. + d__1 * (d__1 * d__1) * 
	    40.;
/* Computing 3rd power */
    d__1 = l2_3.x[1] - l2_3.x[2] * 2.;
    l4_3.gf[1] = l2_3.x[0] * 20. + l2_3.x[1] * 200. + d__1 * (d__1 * d__1) * 
	    4.;
/* Computing 3rd power */
    d__1 = l2_3.x[1] - l2_3.x[2] * 2.;
    l4_3.gf[2] = l2_3.x[2] * 10. - l2_3.x[3] * 10. - d__1 * (d__1 * d__1) * 
	    8.;
/* Computing 3rd power */
    d__1 = l2_3.x[0] - l2_3.x[3];
    l4_3.gf[3] = l2_3.x[2] * -10. + l2_3.x[3] * 10. - d__1 * (d__1 * d__1) * 
	    40.;
    return 0;
labelL4:
    if (l9_2.index1[0]) {
	l3_1.g[0] = l2_3.x[0] + l2_3.x[1] + l2_3.x[2] + l2_3.x[3] - 1.;
    }
labelL5:
    return 0;
} /* tp354_ */


/* Subroutine */ int tp355_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;
    static Real r__[4], h0, h1, h2, h3, h4, h5, h6;

    h0 = l2_3.x[1] * l2_3.x[3];
    r__[0] = 11. - l2_3.x[0] * l2_3.x[3] - h0 + l2_3.x[2] * l2_3.x[3];
    r__[1] = l2_3.x[0] + l2_3.x[1] * 10. - l2_3.x[2] + l2_3.x[3] + h0 * (
	    l2_3.x[2] - l2_3.x[0]);
    r__[2] = 11. - l2_3.x[0] * 4. * l2_3.x[3] - h0 * 4. + l2_3.x[2] * l2_3.x[
	    3];
    r__[3] = l2_3.x[0] * 2. + l2_3.x[1] * 20. - l2_3.x[2] * .5 + l2_3.x[3] * 
	    2. + h0 * 2. * (l2_3.x[2] - l2_3.x[0] * 4.);
    h1 = r__[0] * 2.;
    h2 = r__[1] * 2.;
    h3 = r__[2] * 2.;
    h4 = r__[3] * 2.;
    h5 = h1 * l2_3.x[3];
    h6 = h2 * (1. - h0);
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 4; ++i__) {
	l12_3.lxu[i__ - 1] = false;
	l11_3.lxl[i__ - 1] = true;
/* labelL6: */
	l2_3.x[i__ - 1] = 0.;
    }
    l13_3.xl[0] = .1;
    l13_3.xl[1] = .1;
    l13_3.xl[2] = 0.;
    l13_3.xl[3] = 0.;
    l20_6.lex = false;
    l20_6.nex = 1;
    l20_6.fex = 69.675463;
    l20_6.xex[0] = 1.916633;
    l20_6.xex[1] = .1;
    l20_6.xex[2] = 0.;
    l20_6.xex[3] = 1.9718118;
    return 0;
labelL2:
    l6_1.fx = r__[0] * r__[0] + r__[1] * r__[1];
    return 0;
labelL3:
    l4_3.gf[0] = -h5 + h6;
    l4_3.gf[1] = -h5 + h2 * (l2_3.x[3] * (l2_3.x[2] - l2_3.x[0]) + 10.);
    l4_3.gf[2] = h5 - h6;
    l4_3.gf[3] = h1 * (-l2_3.x[0] - l2_3.x[1] + l2_3.x[2]) + h2 * (l2_3.x[1] *
	     (l2_3.x[2] - l2_3.x[0]) + 1.);
    return 0;
labelL4:
    if (l9_2.index1[0]) {
/* Computing 2nd power */
	d__1 = r__[0];
/* Computing 2nd power */
	d__2 = r__[1];
/* Computing 2nd power */
	d__3 = r__[2];
/* Computing 2nd power */
	d__4 = r__[3];
	l3_1.g[0] = d__1 * d__1 + d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
    }
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto L7;
    }
    l5_2.gg[0] = -h5 + h6 - h3 * -4. * l2_3.x[3] - h4 * (2. - h0 * 8.);
    l5_2.gg[1] = -h5 + h2 * (l2_3.x[3] * (l2_3.x[2] - l2_3.x[0]) + 10.) - h3 *
	     -4. * l2_3.x[3] - h4 * (l2_3.x[3] * 2. * (l2_3.x[2] - l2_3.x[0] *
	     4.) + 20.);
    l5_2.gg[2] = h5 - h6 - h3 * l2_3.x[3] - h4 * (h0 * 2. - .5);
    l5_2.gg[3] = h1 * (-l2_3.x[0] - l2_3.x[1] + l2_3.x[2]) + h2 * (l2_3.x[1] *
	     (l2_3.x[2] - l2_3.x[0]) + 1.) - h3 * (l2_3.x[0] * -4. - l2_3.x[1]
	     * 4. + l2_3.x[2]) - h4 * (l2_3.x[1] * 2. * (l2_3.x[2] - l2_3.x[0]
	     * 4.) + 2.);
L7:
    return 0;
} /* tp355_ */


/* Subroutine */ int tp356_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real eidc, load, cosa, sigd, eitc, e;
    static int i__;
    static Real j, l, m, r__, t, reidc, reitc, t1, t2, fh, ei, gh, gj, 
	    pc, td, wp, del, phi[5], sig;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL3;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 1;
    l1_1.ninl = 4;
    l1_1.nenl = 0;
    l1_1.nenl = 0;
    l2_3.x[0] = 1.;
    l2_3.x[1] = 7.;
    l2_3.x[2] = 8.;
    l2_3.x[3] = 1.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l12_3.lxu[i__ - 1] = false;
/* labelL6: */
	l11_3.lxl[i__ - 1] = true;
    }
    l11_3.lxl[3] = false;
    l12_3.lxu[3] = false;
    l13_3.xl[0] = .125;
    l13_3.xl[1] = 0.;
    l13_3.xl[2] = 0.;
    l20_6.lex = false;
    l20_6.nex = 1;
    l20_6.fex = 2.3811648;
    l20_6.xex[0] = .24436898;
    l20_6.xex[1] = 6.2187934;
    l20_6.xex[2] = 8.2914714;
    l20_6.xex[3] = .24436898;
    return 0;
labelL2:
    l6_1.fx = l2_3.x[0] * 1.10471 * l2_3.x[0] * l2_3.x[1] + l2_3.x[2] * 
	    .04811 * l2_3.x[3] * (l2_3.x[1] + 14.);
labelL3:
    return 0;
labelL4:
    l = 14.;
    load = 6e3;
    td = 13600.;
    sigd = 3e4;
    fh = load;
    t1 = fh / (l2_3.x[0] * 1.414 * l2_3.x[1]);
    m = fh * (l + l2_3.x[1] / 2.);
/* Computing 2nd power */
    d__1 = (l2_3.x[2] + l2_3.x[0]) / 2.;
    r__ = std::sqrt(l2_3.x[1] * l2_3.x[1] / 4. + d__1 * d__1);
/* Computing 2nd power */
    d__1 = (l2_3.x[2] + l2_3.x[0]) / 2.;
    j = l2_3.x[0] * .707 * l2_3.x[1] * (l2_3.x[1] * l2_3.x[1] / 12. + d__1 * 
	    d__1) * 2.;
    t2 = m * r__ / j;
    cosa = l2_3.x[1] / (r__ * 2.);
    wp = (d__1 = t1 * t1 + t1 * 2. * t2 * cosa + t2 * t2, std::abs(d__1));
    if (wp > (float)0.) {
	t = std::sqrt(wp);
    } else {
	t = (float)0.;
    }
    sig = fh * 6. * l / (l2_3.x[3] * l2_3.x[2] * l2_3.x[2]);
    phi[0] = (td - t) / 1e4;
    phi[1] = (sigd - sig) / 1e4;
    phi[2] = l2_3.x[3] - l2_3.x[0];
    e = 3e7;
    ei = e * l2_3.x[2] * l2_3.x[3] * l2_3.x[3] * l2_3.x[3] / 12.;
    gh = 1.2e7;
    gj = gh * l2_3.x[2] * l2_3.x[3] * l2_3.x[3] * l2_3.x[3] / 3.;
    eitc = ei * gj;
    eidc = ei / gj;
    if (eitc > (float)0.) {
	reitc = std::sqrt(eitc);
    } else {
	reitc = (float)0.;
    }
    if (eidc > (float)0.) {
	reidc = std::sqrt(eidc);
    } else {
	reidc = (float)0.;
    }
    pc = reitc * 4.013 * (1. - l2_3.x[2] / (l * 2.) * reidc) / (l * l);
    phi[3] = (pc - fh) / 1e4;
    del = fh * 4. * l * l * l / (e * l2_3.x[3] * l2_3.x[2] * l2_3.x[2] * 
	    l2_3.x[2]);
    phi[4] = .25 - del;
    if (l9_5.index1[0]) {
	l3_4.g[0] = phi[2];
    }
    if (l9_5.index1[1]) {
	l3_4.g[1] = phi[0];
    }
    if (l9_5.index1[2]) {
	l3_4.g[2] = phi[1];
    }
    if (l9_5.index1[3]) {
	l3_4.g[3] = phi[3];
    }
    if (l9_5.index1[4]) {
	l3_4.g[4] = phi[4];
    }
    return 0;
} /* tp356_ */


/* Subroutine */ int tp357_(int *mode)
{
    /* Initialized data */

    static Real xpt[36] = { 113.,110.1,106.2,101.3,95.4,88.8,81.6,74.,
	    66.1,58.4,51.,44.3,38.7,34.5,32.4,32.9,36.4,42.8,50.9,59.,65.8,
	    71.5,76.5,81.1,85.6,90.2,94.6,98.9,103.,106.7,109.9,112.5,114.4,
	    115.5,115.7,114.9 };
    static Real ypt[36] = { 40.2,46.8,53.3,59.4,65.,69.9,73.9,76.9,78.9,
	    79.8,79.7,78.5,76.5,73.6,70.2,66.,60.9,54.3,45.8,36.1,26.5,18.1,
	    11.4,6.2,2.6,.3,-.7,-.6,.7,3.1,6.4,10.5,15.5,21.,27.1,33.6 };
    static Real p0 = 90.;
    static Real q0 = 0.;
    static Real r0 = 0.;
    static Real s0 = 0.;

    /* System generated locals */
    Real d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static Real test, test1, a, b, c__;
    static int i__;
    static Real j;
    static int k;
    static Real alpha, calcx, calcy, p1, q1, r1, s1, ca, cp, sa, ph, 
	    dalpha, pi, qi, ri, si, sp, sql, sum, aabb;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL3;
	case 5:  goto labelL3;
    }
labelL1:
    l1_1.n = 4;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.nenl = 0;
    l1_1.nenl = 0;
    l2_3.x[0] = 136.;
    l2_3.x[1] = 0.;
    l2_3.x[2] = 74.8;
    l2_3.x[3] = 75.7;
    for (i__ = 1; i__ <= 4; ++i__) {
	l11_3.lxl[i__ - 1] = true;
/* labelL6: */
	l12_3.lxu[i__ - 1] = true;
    }
    l20_6.lex = false;
    l20_6.nex = 1;
    l20_6.fex = .3584566;
    l20_6.xex[0] = 136.00762;
    l20_6.xex[1] = .031371415;
    l20_6.xex[2] = 73.59439;
    l20_6.xex[3] = 72.187426;
    for (k = 1; k <= 4; ++k) {
/* L7: */
	l13_3.xl[k - 1] = 0.;
    }
    l14_3.xu[0] = 150.;
    l14_3.xu[1] = 50.;
    l14_3.xu[2] = 100.;
    l14_3.xu[3] = 100.;
    return 0;
labelL2:
    dalpha = .17452927777777777;
    sum = 0.;
    p1 = l2_3.x[0];
    q1 = l2_3.x[1];
    r1 = l2_3.x[2];
    s1 = l2_3.x[3];
    for (i__ = 2; i__ <= 36; ++i__) {
	alpha = dalpha * (i__ - 1);
	ca =std::cos(alpha);
	sa = std::sin(alpha);
	pi = p1 * ca - q1 * sa + p0 * (1. - ca) + q0 * sa;
	qi = p1 * sa + q1 * ca + q0 * (1. - ca) - p0 * sa;
	a = r0 * s1 - s0 * r1 - q1 * r0 + p1 * s0 + pi * q1 - p1 * qi + qi * 
		r1 - pi * s1;
	b = -r0 * r1 - s0 * s1 + p1 * r0 + q1 * s0 - p1 * pi - q1 * qi + pi * 
		r1 + qi * s1;
	c__ = -r1 * r0 - s1 * s0 + pi * r0 + qi * s0 + p1 * r1 + q1 * s1 - (
		p1 * p1 + q1 * q1 + pi * pi + qi * qi) / 2.;
	aabb = a * a + b * b;
	if (aabb < 1e-30) {
	    goto L50;
	}
	test = c__ / std::sqrt(aabb);
	if (std::abs(test) > 1.) {
	    goto L51;
	}
	j = 1.;
L52:
	ph = std::asin(test) - std::atan(b / a);
L55:
	sp = std::sin(ph);
	cp =std::cos(ph);
	ri = r1 * cp - s1 * sp + pi - p1 * cp + q1 * sp;
	si = r1 * sp + s1 * cp + qi - p1 * sp - q1 * cp;
/* Computing 2nd power */
	d__1 = r1 - r0;
/* Computing 2nd power */
	d__2 = s1 - s0;
	test1 = d__1 * d__1 + d__2 * d__2;
	if (test1 < 1e-10) {
	    test1 = 1e-10;
	}
/* Computing 2nd power */
	d__2 = ri - r0;
/* Computing 2nd power */
	d__3 = si - s0;
	if ((d__1 = (test1 - d__2 * d__2 - d__3 * d__3) / test1, std::abs(d__1)) < 
		.001) {
	    goto L53;
	}
	if (j == 2.) {
	    goto L51;
	}
	test = -test;
	j = 2.;
	goto L52;
L50:
	ph = -std::atan(b / a);
	goto L55;
L51:
	l6_1.fx = 1e20;
	return 0;
L53:
	calcx = xpt[0] * cp - ypt[0] * sp + pi - p1 * cp + q1 * sp;
	calcy = xpt[0] * sp + ypt[0] * cp + qi - p1 * sp - q1 * cp;
/* L54: */
/* Computing 2nd power */
	d__1 = calcx - xpt[i__ - 1];
/* Computing 2nd power */
	d__2 = calcy - ypt[i__ - 1];
	sum = sum + d__1 * d__1 + d__2 * d__2;
    }
/* Computing 2nd power */
    d__1 = r1 - r0;
/* Computing 2nd power */
    d__2 = s1 - s0;
/* Computing 2nd power */
    d__3 = r1 - p1;
/* Computing 2nd power */
    d__4 = s1 - q1;
/* Computing 2nd power */
    d__5 = p1 - p0;
/* Computing 2nd power */
    d__6 = q1 - q0;
    sql = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4 + d__5 * d__5 
	    + d__6 * d__6;
    l6_1.fx = sum / 100. + sql / 62500.;
labelL3:
    return 0;
} /* tp357_ */


/* Subroutine */ int tp358_(int *mode)
{
    /* Initialized data */

    static Real y[33] = { .844,.908,.932,.936,.925,.908,.881,.85,.818,
	    .784,.751,.718,.685,.658,.628,.603,.58,.558,.538,.522,.506,.49,
	    .478,.467,.457,.448,.438,.431,.424,.42,.414,.411,.406 };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__, j;
    static Real ti;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.nenl = 0;
    l1_1.nenl = 0;
    l2_4.x[0] = .5;
    l2_4.x[1] = 1.5;
    l2_4.x[2] = -1.;
    l2_4.x[3] = .01;
    l2_4.x[4] = .02;
    for (i__ = 1; i__ <= 5; ++i__) {
	l12_4.lxu[i__ - 1] = true;
/* labelL6: */
	l11_4.lxl[i__ - 1] = true;
    }
    l13_4.xl[0] = -.5;
    l14_4.xu[0] = .5;
    l13_4.xl[1] = 1.5;
    l14_4.xu[1] = 2.5;
    l13_4.xl[2] = -2.;
    l14_4.xu[2] = -1.;
    l13_4.xl[3] = .001;
    l14_4.xu[3] = .1;
    l13_4.xl[4] = .001;
    l14_4.xu[4] = .1;
    l20_7.lex = false;
    l20_7.nex = 1;
    l20_7.fex = 5.46e-5;
    l20_7.xex[0] = .3754;
    l20_7.xex[1] = 1.9358;
    l20_7.xex[2] = -1.4647;
    l20_7.xex[3] = .01287;
    l20_7.xex[4] = .02212;
    l15_1.lsum = 33;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 33; ++i__) {
	ti = (Real) (i__ - 1) * 10.;
/* labelL20: */
	l16_17.f[i__ - 1] = y[i__ - 1] - (l2_4.x[0] + l2_4.x[1] * std::exp(-l2_4.x[
		3] * ti) + l2_4.x[2] * std::exp(-l2_4.x[4] * ti));
    }
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 33; ++i__) {
/* L7: */
/* Computing 2nd power */
	d__1 = l16_17.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (j = 1; j <= 5; ++j) {
/* L8: */
	l4_4.gf[j - 1] = 0.;
    }
    for (i__ = 1; i__ <= 33; ++i__) {
	ti = (Real) (i__ - 1) * 10.;
	l17_20.df[i__ - 1] = -1.;
	l17_20.df[i__ + 32] = -exp(-l2_4.x[3] * ti);
	l17_20.df[i__ + 65] = -exp(-l2_4.x[4] * ti);
	l17_20.df[i__ + 98] = l2_4.x[1] * std::exp(-l2_4.x[3] * ti) * ti;
	l17_20.df[i__ + 131] = l2_4.x[2] * std::exp(-l2_4.x[4] * ti) * ti;
	for (j = 1; j <= 5; ++j) {
/* labelL10: */
	    l4_4.gf[j - 1] += l16_17.f[i__ - 1] * 2. * l17_20.df[i__ + j * 33 
		    - 34];
	}
    }
labelL4:
    return 0;
} /* tp358_ */


/* Subroutine */ int tp359_(int *mode)
{
    /* Initialized data */

    static Real b[5] = { -145421.402,2931.1506,-40.427932,5106.192,
	    15711.36 };
    static Real c__[5] = { -155011.1084,4360.53352,12.9492344,10236.884,
	    13176.786 };
    static Real d__[5] = { -326669.5104,7390.68412,-27.8986976,
	    16643.076,30988.146 };
    static Real a[5] = { -8720288.849,150512.5253,-156.6950325,
	    476470.3222,729482.8271 };

    static Real h__[6];
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 14;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_4.x[0] = 2.52;
    l2_4.x[1] = 5.04;
    l2_4.x[2] = 94.5;
    l2_4.x[3] = 23.31;
    l2_4.x[4] = 17.136;
    for (i__ = 1; i__ <= 5; ++i__) {
/* labelL10: */
	l4_4.gf[i__ - 1] = -a[i__ - 1];
    }
    for (j = 1; j <= 8; ++j) {
	for (i__ = 1; i__ <= 5; ++i__) {
/* labelL6: */
	    l5_33.gg[j + i__ * 14 - 15] = 0.;
	}
    }
    for (i__ = 1; i__ <= 4; ++i__) {
	l5_33.gg[(i__ << 1) - 1 + (i__ + 1) * 14 - 15] = -1.;
/* L7: */
	l5_33.gg[(i__ << 1) + (i__ + 1) * 14 - 15] = 1.;
    }
    l5_33.gg[0] = 2.4;
    l5_33.gg[1] = -1.2;
    l5_33.gg[2] = 60.;
    l5_33.gg[3] = -20.;
    l5_33.gg[4] = 9.3;
    l5_33.gg[5] = -9.;
    l5_33.gg[6] = 7.;
    l5_33.gg[7] = -6.5;
    for (i__ = 1; i__ <= 5; ++i__) {
	l5_33.gg[i__ * 14 - 6] = b[i__ - 1];
	l5_33.gg[i__ * 14 - 5] = c__[i__ - 1];
	l5_33.gg[i__ * 14 - 4] = d__[i__ - 1];
	l5_33.gg[i__ * 14 - 3] = -b[i__ - 1];
	l5_33.gg[i__ * 14 - 2] = -c__[i__ - 1];
	l5_33.gg[i__ * 14 - 1] = -d__[i__ - 1];
	l12_4.lxu[i__ - 1] = false;
/* L8: */
	l11_4.lxl[i__ - 1] = false;
    }
    l20_7.lex = false;
    l20_7.nex = 1;
    l20_7.fex = -5280416.8;
    l20_7.xex[0] = 4.5374;
    l20_7.xex[1] = 10.889;
    l20_7.xex[2] = 272.24;
    l20_7.xex[3] = 42.198;
    l20_7.xex[4] = 31.762;
    return 0;
labelL2:
    l6_1.fx = -24345.;
    for (i__ = 1; i__ <= 5; ++i__) {
/* labelL9: */
	l6_1.fx += a[i__ - 1] * l2_4.x[i__ - 1];
    }
    l6_1.fx = -l6_1.fx;
labelL3:
    return 0;
labelL4:
    if (l9_8.index1[0]) {
	l3_7.g[0] = l2_4.x[0] * 2.4 - l2_4.x[1];
    }
    if (l9_8.index1[1]) {
	l3_7.g[1] = l2_4.x[0] * -1.2 + l2_4.x[1];
    }
    if (l9_8.index1[2]) {
	l3_7.g[2] = l2_4.x[0] * 60. - l2_4.x[2];
    }
    if (l9_8.index1[3]) {
	l3_7.g[3] = l2_4.x[0] * -20. + l2_4.x[2];
    }
    if (l9_8.index1[4]) {
	l3_7.g[4] = l2_4.x[0] * 9.3 - l2_4.x[3];
    }
    if (l9_8.index1[5]) {
	l3_7.g[5] = l2_4.x[0] * -9. + l2_4.x[3];
    }
    if (l9_8.index1[6]) {
	l3_7.g[6] = l2_4.x[0] * 7. - l2_4.x[4];
    }
    if (l9_8.index1[7]) {
	l3_7.g[7] = l2_4.x[0] * -6.5 + l2_4.x[4];
    }
    for (i__ = 1; i__ <= 3; ++i__) {
/* labelL11: */
	h__[i__ - 1] = 0.;
    }
    h__[3] = 2.94e5;
    h__[4] = 2.94e5;
    h__[5] = 277200.;
    for (i__ = 1; i__ <= 5; ++i__) {
	h__[0] += b[i__ - 1] * l2_4.x[i__ - 1];
	h__[1] += c__[i__ - 1] * l2_4.x[i__ - 1];
	h__[2] += d__[i__ - 1] * l2_4.x[i__ - 1];
	h__[3] -= b[i__ - 1] * l2_4.x[i__ - 1];
	h__[4] -= c__[i__ - 1] * l2_4.x[i__ - 1];
/* labelL12: */
	h__[5] -= d__[i__ - 1] * l2_4.x[i__ - 1];
    }
    for (i__ = 9; i__ <= 14; ++i__) {
/* labelL13: */
	if (l9_8.index1[i__ - 1]) {
	    l3_7.g[i__ - 1] = h__[i__ - 9];
	}
    }
labelL5:
    return 0;
} /* tp359_ */


/* Subroutine */ int tp360_(int *mode)
{
    /* Initialized data */

    static Real c__[10] = { -8720288.849,150512.5233,-156.6950325,
	    476470.3222,729482.8271,-326669.5104,7390.68412,-27.8986976,
	    16643.076,30988.146 };

    static Real h__;
    static int i__;
    static Real hh[5];

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_4.x[0] = 2.52;
    l2_4.x[1] = 2.;
    l2_4.x[2] = 37.5;
    l2_4.x[3] = 9.25;
    l2_4.x[4] = 6.8;
    l11_4.lxl[0] = true;
    l12_4.lxu[0] = false;
    for (i__ = 2; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = true;
/* labelL6: */
	l12_4.lxu[i__ - 1] = true;
    }
    l13_4.xl[0] = 0.;
    l13_4.xl[1] = 1.2;
    l13_4.xl[2] = 20.;
    l13_4.xl[3] = 9.;
    l13_4.xl[4] = 6.5;
    l14_4.xu[1] = 2.4;
    l14_4.xu[2] = 60.;
    l14_4.xu[3] = 9.3;
    l14_4.xu[4] = 7.;
    l20_7.lex = false;
    l20_7.nex = 1;
    l20_7.fex = -5280335.1;
    l20_7.xex[0] = 4.537431;
    l20_7.xex[1] = 2.4;
    l20_7.xex[2] = 60.;
    l20_7.xex[3] = 9.3;
    l20_7.xex[4] = 7.;
    return 0;
labelL2:
    l6_1.fx = (-c__[0] - c__[1] * l2_4.x[1] - c__[2] * l2_4.x[2] - c__[3] * 
	    l2_4.x[3] - c__[4] * l2_4.x[4]) * l2_4.x[0] + 24345.;
    return 0;
labelL3:
    l4_4.gf[0] = -c__[0] - c__[1] * l2_4.x[1] - c__[2] * l2_4.x[2] - c__[3] * 
	    l2_4.x[3] - c__[4] * l2_4.x[4];
    for (i__ = 2; i__ <= 5; ++i__) {
/* L7: */
	l4_4.gf[i__ - 1] = -c__[i__ - 1] * l2_4.x[0];
    }
    return 0;
labelL4:
    h__ = (c__[5] + c__[6] * l2_4.x[1] + c__[7] * l2_4.x[2] + c__[8] * l2_4.x[
	    3] + c__[9] * l2_4.x[4]) * l2_4.x[0];
    if (l9_3.index1[0]) {
	l3_2.g[0] = h__;
    }
    if (l9_3.index1[1]) {
	l3_2.g[1] = 277200. - h__;
    }
    return 0;
labelL5:
    hh[0] = c__[5] + c__[6] * l2_4.x[1] + c__[7] * l2_4.x[2] + c__[8] * 
	    l2_4.x[3] + c__[9] * l2_4.x[4];
    for (i__ = 2; i__ <= 5; ++i__) {
/* L8: */
	hh[i__ - 1] = c__[i__ + 4] * l2_4.x[0];
    }
    if (! l10_3.index2[0]) {
	goto labelL11;
    }
    for (i__ = 1; i__ <= 5; ++i__) {
/* labelL9: */
	l5_4.gg[(i__ << 1) - 2] = hh[i__ - 1];
    }
labelL11:
    if (! l10_3.index2[1]) {
	goto labelL12;
    }
    for (i__ = 1; i__ <= 5; ++i__) {
/* labelL10: */
	l5_4.gg[(i__ << 1) - 1] = -hh[i__ - 1];
    }
labelL12:
    return 0;
} /* tp360_ */


/* Subroutine */ int tp361_(int *mode)
{
    /* Initialized data */

    static Real a[5] = { -8720288.849,150512.5253,-156.6950325,
	    476470.3222,729482.8271 };
    static Real b[5] = { -145421.402,2931.1506,-40.427932,5106.192,
	    15711.36 };
    static Real c__[5] = { -155011.1084,4360.53352,12.9492344,10236.884,
	    13176.786 };
    static Real d__[5] = { -326669.5104,7390.68412,-27.8986976,
	    16643.076,30988.146 };

    static Real h__[3];
    static int i__, j;
    static Real hh[15]	/* was [3][5] */;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 6;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_4.x[0] = 2.52;
    l2_4.x[1] = 2.;
    l2_4.x[2] = 37.5;
    l2_4.x[3] = 9.25;
    l2_4.x[4] = 6.8;
    l11_4.lxl[4] = false;
    l12_4.lxu[0] = false;
    for (i__ = 1; i__ <= 4; ++i__) {
	l11_4.lxl[i__ - 1] = true;
/* labelL6: */
	l12_4.lxu[i__] = true;
    }
    l13_4.xl[0] = 0.;
    l13_4.xl[1] = 1.2;
    l13_4.xl[2] = 20.;
    l13_4.xl[3] = 9.;
    l14_4.xu[1] = 2.4;
    l14_4.xu[2] = 60.;
    l14_4.xu[3] = 9.3;
    l14_4.xu[4] = 7.;
    l20_7.lex = false;
    l20_7.nex = 1;
    l20_7.fex = -776412.12;
    l20_7.xex[0] = .68128605;
    l20_7.xex[1] = .24;
    l20_7.xex[2] = 20.;
    l20_7.xex[3] = 9.3;
    l20_7.xex[4] = 7.;
    return 0;
labelL2:
    l6_1.fx = a[0];
    for (i__ = 2; i__ <= 5; ++i__) {
/* L7: */
	l6_1.fx += a[i__ - 1] * l2_4.x[i__ - 1];
    }
    l6_1.fx = l2_4.x[0] * l6_1.fx - 24345.;
    l6_1.fx = -l6_1.fx;
    return 0;
labelL3:
    l4_4.gf[0] = a[0];
    for (i__ = 2; i__ <= 5; ++i__) {
	l4_4.gf[0] += a[i__ - 1] * l2_4.x[i__ - 1];
/* L8: */
	l4_4.gf[i__ - 1] = a[i__ - 1] * l2_4.x[0];
    }
    for (i__ = 1; i__ <= 5; ++i__) {
/* labelL20: */
	l4_4.gf[i__ - 1] = -l4_4.gf[i__ - 1];
    }
    return 0;
labelL4:
    h__[0] = b[0];
    h__[1] = c__[0];
    h__[2] = d__[0];
    for (i__ = 2; i__ <= 5; ++i__) {
	h__[0] += b[i__ - 1] * l2_4.x[i__ - 1];
	h__[1] += c__[i__ - 1] * l2_4.x[i__ - 1];
/* labelL9: */
	h__[2] += d__[i__ - 1] * l2_4.x[i__ - 1];
    }
    for (i__ = 1; i__ <= 3; ++i__) {
/* labelL10: */
	h__[i__ - 1] = l2_4.x[0] * h__[i__ - 1];
    }
    for (i__ = 1; i__ <= 3; ++i__) {
/* labelL11: */
	if (l9_6.index1[i__ - 1]) {
	    l3_5.g[i__ - 1] = h__[i__ - 1];
	}
    }
    if (l9_6.index1[3]) {
	l3_5.g[3] = 29400. - h__[0];
    }
    if (l9_6.index1[4]) {
	l3_5.g[4] = 29400. - h__[1];
    }
    if (l9_6.index1[5]) {
	l3_5.g[5] = 277200. - h__[2];
    }
    return 0;
labelL5:
    hh[0] = b[0];
    hh[1] = c__[0];
    hh[2] = d__[0];
    for (i__ = 2; i__ <= 5; ++i__) {
	hh[0] += b[i__ - 1] * l2_4.x[i__ - 1];
	hh[1] += c__[i__ - 1] * l2_4.x[i__ - 1];
	hh[2] += d__[i__ - 1] * l2_4.x[i__ - 1];
	hh[i__ * 3 - 3] = b[i__ - 1] * l2_4.x[0];
	hh[i__ * 3 - 2] = c__[i__ - 1] * l2_4.x[0];
/* labelL12: */
	hh[i__ * 3 - 1] = d__[i__ - 1] * l2_4.x[0];
    }
    for (j = 1; j <= 3; ++j) {
	if (! l10_6.index2[j - 1]) {
	    goto labelL13;
	}
	for (i__ = 1; i__ <= 5; ++i__) {
/* L130: */
	    l5_14.gg[j + i__ * 6 - 7] = hh[j + i__ * 3 - 4];
	}
labelL13:
	;
    }
    for (j = 4; j <= 6; ++j) {
	if (! l10_6.index2[j - 1]) {
	    goto labelL14;
	}
	for (i__ = 1; i__ <= 5; ++i__) {
/* L140: */
	    l5_14.gg[j + i__ * 6 - 7] = -hh[j - 3 + i__ * 3 - 4];
	}
labelL14:
	;
    }
    return 0;
} /* tp361_ */


/* Subroutine */ int tp362_(int *mode)
{
    Real tp362a_(Real *);
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 4;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_4.x[0] = 15.;
    l2_4.x[1] = 9.05;
    l2_4.x[2] = 6.14;
    l2_4.x[3] = 4.55;
    l2_4.x[4] = 3.61;
    for (i__ = 2; i__ <= 4; ++i__) {
	l11_4.lxl[i__ - 1] = false;
/* labelL14: */
	l12_4.lxu[i__ - 1] = false;
    }
    l11_4.lxl[0] = true;
    l12_4.lxu[0] = true;
    l11_4.lxl[4] = true;
    l12_4.lxu[4] = false;
    l13_4.xl[0] = 15.;
    l14_4.xu[0] = 20.;
    l13_4.xl[4] = 2.;
    l20_7.lex = false;
    l20_7.nex = 1;
    l20_7.fex = .26229998;
    l20_7.xex[0] = 15.050962;
    l20_7.xex[1] = 8.8751199;
    l20_7.xex[2] = 5.908823;
    l20_7.xex[3] = 4.860481;
    l20_7.xex[4] = 4.399269;
    return 0;
labelL2:
    l6_1.fx = tp362a_(l2_4.x);
    l6_1.fx = l6_1.fx;
labelL3:
    return 0;
labelL4:
    for (i__ = 1; i__ <= 4; ++i__) {
/* L41: */
	if (l9_7.index1[i__ - 1]) {
	    l3_6.g[i__ - 1] = l2_4.x[i__ - 1] - l2_4.x[i__];
	}
    }
labelL5:
    return 0;
} /* tp362_ */

Real tp362a_(Real *x)
{
    /* Initialized data */

    static Real rad = 1.085;
    static Real con1 = 1.466667;
    static Real con2 = 12.90842;
    static Real rpmin = 600.;
    static Real rpmax = 5700.;
    static Real ei = .6;
    static Real vi = 98.;
    static Real dt = .01;
    static Real vmax = 100.;
    static Real v0 = 5.;
    static Real tshift = .25;
    static Real tmax = 100.;

    /* System generated locals */
    Real ret_val, d__1, d__2;

    /* Local variables */
    static int i__;
    static Real t, v, force;
    static int it;
    static Real tt, torque, acc, rpm, acc0;

    /* Parameter adjustments */
    --x;

    /* Function Body */
/* labelL13: */
    it = 0;
    acc = 0.;
    v = v0;
    i__ = 1;
L302:
/* Computing 2nd power */
    d__1 = v;
    force = d__1 * d__1 * .0239 + 31.2;
L301:
    rpm = v * con2 * x[i__];
    if (rpm < rpmin) {
	goto L300;
    }
    if (rpm >= rpmax) {
	goto L305;
    }
    if (rpm >= 600. && rpm < 1900.) {
/* Computing 3rd power */
	d__1 = rpm;
/* Computing 2nd power */
	d__2 = rpm;
	torque = d__1 * (d__1 * d__1) * 3.846154e-8 - d__2 * d__2 * 
		2.108974359e-4 + rpm * .42455128205133 - 187.11538461540295;
    }
    if (rpm >= 1900. && rpm < 3e3) {
/* Computing 3rd power */
	d__1 = rpm;
/* Computing 2nd power */
	d__2 = rpm;
	torque = d__1 * (d__1 * d__1) * -4.92424e-9 + d__2 * d__2 * 
		1.867424242e-5 + rpm * .01229545454547 + 64.999999999986;
    }
    if (rpm >= 3e3 && rpm < 4500.) {
/* Computing 3rd power */
	d__1 = rpm;
/* Computing 2nd power */
	d__2 = rpm;
	torque = d__1 * (d__1 * d__1) * -2.6667e-10 + d__2 * d__2 * 3e-6 - 
		rpm * .01263333333336 + 155.10000000002947;
    }
    if (rpm >= 4500. && rpm < 5600.) {
/* Computing 3rd power */
	d__1 = rpm;
/* Computing 2nd power */
	d__2 = rpm;
	torque = d__1 * (d__1 * d__1) * -6.64141e-9 + d__2 * d__2 * 
		8.337626263e-5 - rpm * .34351868688129 + 597.3636363847145;
    }
    if (rpm >= 5600. && rpm < 6e3) {
/* Computing 3rd power */
	d__1 = rpm;
/* Computing 2nd power */
	d__2 = rpm;
	torque = d__1 * (d__1 * d__1) * -2.539683e-8 + d__2 * d__2 * 
		3.8158730157e-4 - rpm * 1.9223492062348 + 3380.66666645715304;
    }
    acc0 = acc;
/* Computing 2nd power */
    d__1 = x[i__];
    acc = rad * (x[i__] * torque - force * rad) / (ei * (d__1 * d__1) + vi);
    ++it;
    t = dt * it;
    v += (acc0 + acc) / 2. * dt / con1;
    if (t > tmax) {
	goto L311;
    }
    if (v >= vmax) {
	goto L311;
    }
    goto L302;
L300:
    ret_val = tmax;
    return ret_val;
L305:
    ++i__;
    if (i__ > 5) {
	goto L311;
    }
    if (t == 0.) {
	goto L301;
    }
    tt = t + tshift;
L306:
/* Computing 2nd power */
    d__1 = rad;
    acc = -force * (d__1 * d__1) / vi;
    ++it;
    t = dt * it;
    v += acc * dt / con1;
    if (t < tt) {
	goto L307;
    }
    goto L302;
L307:
/* Computing 2nd power */
    d__1 = v;
    force = d__1 * d__1 * .0293 + 31.2;
    goto L306;
L311:
    ret_val = t / 100.;
    return ret_val;
} /* tp362a_ */


/* Subroutine */ int tp363_(int *mode)
{
    Real tp363a_(Real *);
    int tp363b_(bool *, Real *, Real *);
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 5;
    l1_1.nili = 0;
    l1_1.ninl = 3;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_4.x[0] = -.3359769;
    l2_4.x[1] = -1.432398;
    l2_4.x[2] = 0.;
    l2_4.x[3] = 4.;
    l2_4.x[4] = 9.;
    for (i__ = 1; i__ <= 5; ++i__) {
	l11_4.lxl[i__ - 1] = true;
	l12_4.lxu[i__ - 1] = true;
/* labelL11: */
	l14_4.xu[i__ - 1] = 10.;
    }
    l13_4.xl[0] = -10.;
    l13_4.xl[1] = -10.;
    l13_4.xl[2] = -10.;
    l13_4.xl[3] = .1;
    l13_4.xl[4] = 1.;
    l20_7.lex = false;
    l20_7.nex = 1;
    l20_7.xex[0] = .19852438;
    l20_7.xex[1] = -3.01059794;
    l20_7.xex[2] = -.0530266138;
    l20_7.xex[3] = 2.83165094;
    l20_7.xex[4] = 10.;
    l20_7.fex = (float)-5.55840576;
    return 0;
labelL2:
    l6_1.fx = tp363a_(l2_4.x);
labelL3:
    return 0;
labelL4:
    tp363b_(l9_4.index1, l2_4.x, l3_3.g);
labelL5:
    return 0;
} /* tp363_ */

Real tp363a_(Real *x)
{
    /* Initialized data */

    static Real xi = 1.;
    static Real x6 = 1.;
    static Real x11 = 1.;

    /* System generated locals */
    int i__1;
    Real ret_val;

    /* Local variables */
    int tp363c_(Real *x, Real *thick, Real *dthick, Real *c__, int *kkk, Real *c0, Real *xi); 
    int tp363d_(Real *, Real *);
    static int i__;
    static Real dr, xc[100];

    /* Parameter adjustments */
    --x;

    /* Function Body */
    b_1.w2 = x6 * 100.;
    dr = x[5] - xi;
    xc[0] = xi;
    i__1 = b_1.kkk;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L60: */
	xc[i__] = xc[i__ - 1] + dr / (Real) b_1.kkk;
    }
    tp363c_(xc, b_1.thick, b_1.dthick, &x[1], &b_1.kkk, &x11, &xi);
    tp363d_(xc, &ret_val);
    ret_val = -ret_val / 1e6;
    return ret_val;
} /* tp363a_ */

/* Subroutine */ int tp363b_(bool *index1, Real *c__, Real *g)
{
    /* Initialized data */

    static Real yi = 0.;
    static Real xi = 1.;
    static Real eps = .1;
    static Real c6 = 1.;

    /* System generated locals */
    int i__1;
    Real d__1, d__2, d__3;

    /* Local variables */
    static Real smax, root, stot[100];
    int tp363c_(Real *, Real *, Real *, Real *, int *, Real *, Real *); 
    int tp363e_(Real *, Real *, Real *, Real *, int *, Real *, Real *, Real*); 
    int tp363f_(Real *, Real *, Real *, Real *, Real *, Real *, Real *, Real *, Real *, Real *, Real *, Real *, int *); 
    int tp363g_(Real *, Real *, Real *, Real *, int *, Real *);
    static int i__;
    static Real v;
    static Real x[100], y[100], y1[100];
    static int ii;
    static Real dr;
    static int nn;
    static Real xl, xr, y1i, fxl, fxr, rst[100], tst[100];

    /* Parameter adjustments */
    --g;
    --c__;
    --index1;

    /* Function Body */
    b_1.xmu = (float).3;
    b_1.rho = 7.263e-4;
    ii = 99;
    b_1.kkk = 98;
    b_1.w2 = c6 * 100.;
    dr = c__[5] - xi;
    x[0] = xi;
    i__1 = b_1.kkk;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L60: */
	x[i__] = x[i__ - 1] + dr / (Real) b_1.kkk;
    }
    tp363c_(x, b_1.thick, b_1.dthick, &c__[1], &b_1.kkk, &c__[4], &xi);
/* labelL5: */
    y1i = 1e5;
/* labelL2: */
    tp363e_(&y1i, &yi, &xi, &c__[5], &ii, y, y1, x);
    fxl = -y[ii - 1];
    xl = y1i;
    yi = 0.;
    y1i = 2.503e6;
    tp363e_(&y1i, &yi, &xi, &c__[5], &ii, y, y1, x);
    fxr = -y[ii - 1];
    xr = y1i;
    tp363f_(&xl, &xr, &fxl, &fxr, &eps, y, y1, x, &xi, &c__[5], &yi, &root, &ii);
    smax = 0.;
    st_1.value2 = 0.;
    i__1 = ii;
    for (nn = 1; nn <= i__1; ++nn) {
	rst[nn - 1] = y[nn - 1] / (b_1.thick[nn - 1] * x[nn - 1]);
/* Computing 2nd power */
	d__1 = b_1.w2;
/* Computing 2nd power */
	d__2 = x[nn - 1];
	tst[nn - 1] = (y1[nn - 1] + b_1.thick[nn - 1] * b_1.rho * (d__1 * 
		d__1) * (d__2 * d__2)) / b_1.thick[nn - 1];
/* Computing 2nd power */
	d__1 = rst[nn - 1] - tst[nn - 1];
/* Computing 2nd power */
	d__2 = rst[nn - 1];
/* Computing 2nd power */
	d__3 = tst[nn - 1];
	stot[nn - 1] = std::sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	if (stot[nn - 1] > smax) {
	    smax = stot[nn - 1];
	}
/* Computing 2nd power */
	d__1 = 3e4 - stot[nn - 1];
	st_1.value2 += d__1 * d__1;
/* L300: */
    }
    st_1.value2 = std::sqrt(st_1.value2);
    if (index1[1]) {
	g[1] = (3e4 - smax) / 1e3;
    }
    if (index1[2]) {
	g[2] = 5. - tfn1_1.tmax;
    }
    tp363g_(&c__[4], &c__[1], &v, b_1.thick, &b_1.kkk, x);
    if (index1[3]) {
	g[3] = (625. - v) / 10.;
    }
/* L900: */
    return 0;
} /* tp363b_ */

/* Subroutine */ int tp363c_(Real *x, Real *thick, Real *dthick, Real *c__, int *kkk, Real *c0, Real *xi)
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int nfst, i__, lm;
    static Real xl;
    static int jkl;

    /* Parameter adjustments */
    --c__;
    --dthick;
    --thick;
    --x;

    /* Function Body */
    thick[1] = *c0;
    xl = x[*kkk + 1] - x[1];
    nfst = 2;
    tfn1_1.tmax = thick[1];
    i__1 = *kkk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	thick[i__ + 1] = *c0 + c__[1] * (x[i__ + 1] - *xi);
	i__2 = nfst;
	for (lm = 1; lm <= i__2; ++lm) {
	    jkl = lm + 1;
/* labelL9: */
	    thick[i__ + 1] += c__[jkl] * std::sin(((Real) jkl * (float)2. - (
		    float)3.) * 3.1415926535897932 * (x[i__] - x[1]) / xl);
	}
	if (thick[i__ + 1] > tfn1_1.tmax) {
	    tfn1_1.tmax = thick[i__ + 1];
	}
/* labelL10: */
	dthick[i__] = (thick[i__ + 1] - thick[i__]) / (x[i__ + 1] - x[i__]);
    }
    return 0;
} /* tp363c_ */

/* Subroutine */ int tp363d_(Real *x, Real *xke)
{
    /* System generated locals */
    int i__1;
    Real d__1;

    /* Local variables */
    static int i__;
    static Real const__;

    /* Parameter adjustments */
    --x;

    /* Function Body */
/* Computing 2nd power */
    d__1 = b_2.w;
    const__ = b_2.rho * 3.1415926535897932 * (d__1 * d__1);
    *xke = 0.;
    i__1 = b_2.kkk;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* labelL10: */
/* Computing 3rd power */
	d__1 = x[i__];
	*xke += d__1 * (d__1 * d__1) * b_2.thick[i__ - 1] * (x[i__ + 1] - x[
		i__]);
    }
    *xke *= const__;
    return 0;
} /* tp363d_ */

/* Subroutine */ int tp363e_(Real *y1i, Real *yi, Real *xi, Real *xf, int *ii, Real *y, Real *y1, Real *x)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    Real tp363h_(Real *, Real *, Real *, int *);
    static Real h__;
    static int j;
    static Real m0, m1, m2, m3;
    static int kk, ll;
    static Real xr, yr, y1r;

    /* Parameter adjustments */
    --x;
    --y1;
    --y;

    /* Function Body */
    x[1] = *xi;
    y[1] = *yi;
    y1[1] = *y1i;
    h__ = (*xf - *xi) / (Real) (*ii - 1);
    kk = *ii - 1;
    i__1 = kk;
    for (j = 1; j <= i__1; ++j) {
	ll = j;
	xr = x[j];
	yr = y[j];
	y1r = y1[j];
	m0 = h__ * tp363h_(&xr, &yr, &y1r, &ll);
	xr = x[j] + h__ / 2.;
	yr = y[j] + h__ * y1[j] / 2.;
	y1r = y1[j] + m0 / 2.;
	m1 = h__ * tp363h_(&xr, &yr, &y1r, &ll);
	yr += h__ * m0 / 4.;
	y1r = y1[j] + m1 / 2.;
	m2 = h__ * tp363h_(&xr, &yr, &y1r, &ll);
	xr = x[j] + h__;
	yr = y[j] + h__ * y1[j] + h__ * m1 / 2.;
	y1r = y1[j] + m2;
	m3 = h__ * tp363h_(&xr, &yr, &y1r, &ll);
	y[j + 1] = y[j] + h__ * y1[j] + h__ / 6. * (m0 + m1 + m2);
	y1[j + 1] = y1[j] + (m0 + m1 * 2. + m2 * 2. + m3) / 6.;
/* labelL10: */
	x[j + 1] = x[j] + h__;
    }
    return 0;
} /* tp363e_ */

/* Subroutine */ int tp363f_(Real *xl, Real *xr, Real *fxl, Real *fxr, Real *eps, Real *y, Real *y1, Real *x, Real *xi, Real *xf, Real *yi, 
	Real *root, int *ii)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real xapp;
    int tp363e_(Real *, Real *, Real *, Real *, int *, Real *, Real *, Real *);
    static Real value, fxapp, xsave;

    /* Parameter adjustments */
    --x;
    --y1;
    --y;

    /* Function Body */
    xsave = *xl;
L105:
    xapp = *xl + *fxl * (*xr - *xl) / (*fxl - *fxr);
    tp363e_(&xapp, yi, xi, xf, ii, &y[1], &y1[1], &x[1]);
    fxapp = -y[*ii];
    if ((d__1 = xapp - xsave, std::abs(d__1)) / xapp <= *eps) {
	goto L250;
    }
    value = fxapp * *fxl;
    if (value < 0.) {
	goto L110;
    }
    *xl = xapp;
    xsave = *xl;
    *fxl = fxapp;
    goto L105;
L110:
    *xr = xapp;
    xsave = *xr;
    *fxr = fxapp;
    goto L105;
L250:
    *root = xapp;
    return 0;
} /* tp363f_ */

/* Subroutine */ int tp363g_(Real *c0, Real *c__, Real *v, Real *thick, int *kkk, Real *x)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;
    static Real deltx, r1, r2, r3, pi;
    static int lmn;

    /* Parameter adjustments */
    --x;
    --thick;
    --c__;

    /* Function Body */
    *v = 0.;
    pi = 3.141592654;
    deltx = (x[*kkk + 1] - x[1]) / *kkk;
    lmn = *kkk - 1;
    i__1 = lmn;
    for (i__ = 1; i__ <= i__1; i__ += 2) {
	r1 = (x[i__ + 1] + x[i__]) / 2.;
	r2 = (x[i__ + 1] + x[i__ + 2]) / 2.;
	r3 = (r1 + r2) / 2.;
/* labelL10: */
	*v += pi * 2. * deltx / 3. * (thick[i__] * r1 + thick[i__ + 1] * 4. * 
		r3 + thick[i__ + 2] * r2);
    }
    return 0;
} /* tp363g_ */

Real tp363h_(Real *xr, Real *yr, Real *y1r, int *i__)
{
    /* System generated locals */
    Real ret_val, d__1, d__2;

/* Computing 2nd power */
    d__1 = *xr;
/* Computing 2nd power */
    d__2 = b_1.w2;
    ret_val = (1. / b_1.thick[*i__ - 1] * b_1.dthick[*i__ - 1] - 1. / *xr) * *
	    y1r + (1. / (d__1 * d__1) - b_1.xmu / (*xr * b_1.thick[*i__ - 1]) 
	    * b_1.dthick[*i__ - 1]) * *yr - (b_1.xmu + 3.) * b_1.rho * (d__2 *
	     d__2) * b_1.thick[*i__ - 1] * *xr;
    return ret_val;
} /* tp363h_ */


/* Subroutine */ int tp364_(int *mode)
{
    /* Local variables */
    Real tp364a_(Real *);
    static int i__;
    static Real xmu1, xmu2;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 6;
    l1_1.nili = 2;
    l1_1.ninl = 2;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_5.x[0] = 1.;
    l2_5.x[1] = 4.5;
    l2_5.x[2] = 4.;
    l2_5.x[3] = 5.;
    l2_5.x[4] = 3.;
    l2_5.x[5] = 3.;
    for (i__ = 1; i__ <= 4; ++i__) {
/* labelL11: */
	l11_5.lxl[i__ - 1] = true;
    }
    l11_5.lxl[4] = false;
    l11_5.lxl[5] = false;
    l12_5.lxu[4] = false;
    l12_5.lxu[5] = false;
    l12_5.lxu[0] = true;
    l12_5.lxu[3] = true;
    l12_5.lxu[1] = false;
    l12_5.lxu[2] = false;
    l13_5.xl[0] = .5;
    l13_5.xl[1] = 0.;
    l13_5.xl[2] = 0.;
    l13_5.xl[3] = 2.;
    l14_5.xu[0] = 3.;
    l14_5.xu[3] = 10.;
    l20_4.lex = false;
    l20_4.nex = 1;
    l20_4.fex = (float).0606002;
    l20_4.xex[0] = (float).99616882;
    l20_4.xex[1] = 4.1960616;
    l20_4.xex[2] = 2.9771652;
    l20_4.xex[3] = 3.9631949;
    l20_4.xex[4] = 1.6536702;
    l20_4.xex[5] = 1.2543998;
    return 0;
labelL2:
    l6_1.fx = tp364a_(l2_5.x);
labelL3:
    return 0;
labelL4:
    xmu1 = .7853981633;
    xmu2 = 2.356194491;
    if (l9_7.index1[0]) {
	l3_6.g[0] = -l2_5.x[0] + l2_5.x[1] + l2_5.x[2] - l2_5.x[3];
    }
    if (l9_7.index1[1]) {
	l3_6.g[1] = -l2_5.x[0] - l2_5.x[1] + l2_5.x[2] + l2_5.x[3];
    }
    if (l9_7.index1[2]) {
	l3_6.g[2] = -l2_5.x[1] * l2_5.x[1] - l2_5.x[2] * l2_5.x[2] + (l2_5.x[
		3] - l2_5.x[0]) * (l2_5.x[3] - l2_5.x[0]) + l2_5.x[1] * 2. * 
		l2_5.x[2] *std::cos(xmu1);
    }
    if (l9_7.index1[3]) {
	l3_6.g[3] = l2_5.x[1] * l2_5.x[1] + l2_5.x[2] * l2_5.x[2] - (l2_5.x[3]
		 + l2_5.x[0]) * (l2_5.x[3] + l2_5.x[0]) - l2_5.x[1] * 2. * 
		l2_5.x[2] *std::cos(xmu2);
    }
labelL5:
    return 0;
} /* tp364_ */

Real tp364a_(Real *x)
{
    /* System generated locals */
    Real ret_val, d__1, d__2;

    /* Local variables */
    static Real xinc, coss, sins, cosy, siny;
    int tp364b_(Real*, Real*, Real*), tp364c_(Real*, Real*, Real*);
    static int i__;
    static Real x1[31], y1[31], pi, wp, x1a[31], y1a[31], phi[31];

    /* Parameter adjustments */
    --x;

    /* Function Body */
    pi = 3.141592654;
    xinc = pi * 2. / 30.;
    for (i__ = 1; i__ <= 31; ++i__) {
/* labelL1: */
	phi[i__ - 1] = xinc * (Real) (i__ - 1);
    }
    tp364b_(phi, x1, y1);
    ret_val = 0.;
    for (i__ = 1; i__ <= 31; ++i__) {
	tp364c_(&x[1], &phi[i__ - 1], &coss);
	wp = (d__1 = 1. - coss * coss, std::abs(d__1));
	if (wp > (float)0.) {
	    sins = std::sqrt(wp);
	} else {
	    sins = (float)0.;
	}
	cosy = (x[4] + x[3] * coss - x[1] *std::cos(phi[i__ - 1])) / x[2];
	siny = (x[3] * sins - x[1] * std::sin(phi[i__ - 1])) / x[2];
	x1a[i__ - 1] = x[1] *std::cos(phi[i__ - 1]) + x[5] * cosy - x[6] * siny;
	y1a[i__ - 1] = x[1] * std::sin(phi[i__ - 1]) + x[5] * siny + x[6] * cosy;
/* labelL2: */
/* Computing 2nd power */
	d__1 = x1a[i__ - 1] - x1[i__ - 1];
/* Computing 2nd power */
	d__2 = y1a[i__ - 1] - y1[i__ - 1];
	ret_val = ret_val + d__1 * d__1 + d__2 * d__2;
    }
    wp = ret_val / 31.;
    if (wp > (float)0.) {
	ret_val = std::sqrt(wp);
    } else {
	ret_val = (float)0.;
    }
    return ret_val;
} /* tp364a_ */

/* Subroutine */ int tp364b_(Real *phi, Real *x1, Real *y1)
{
    /* Local variables */
    static int i__;
    static Real pi;

    /* Parameter adjustments */
    --y1;
    --x1;
    --phi;

    /* Function Body */
    pi = 3.141592654;
    for (i__ = 1; i__ <= 31; ++i__) {
	x1[i__] = std::sin(pi * 2. * ((pi - phi[i__]) / (pi * 2.) - .16)) + .4;
/* labelL1: */
	y1[i__] = std::sin(pi - phi[i__]) * .9 + 2.;
    }
    return 0;
} /* tp364b_ */

/* Subroutine */ int tp364c_(Real *x, Real *phi, Real *w)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real term, a, b, c__, k, l, m, pi;

    /* Parameter adjustments */
    --x;

    /* Function Body */
    pi = 3.141592654;
    m = x[1] * 2. * x[3] * std::sin(*phi);
    l = x[3] * 2. * x[4] - x[1] * 2. * x[3] *std::cos(*phi);
    k = x[1] * x[1] - x[2] * x[2] + x[3] * x[3] + x[4] * x[4] - x[4] * 2. * x[
	    1] *std::cos(*phi);
    a = l * l + m * m;
    b = k * 2. * l;
    c__ = k * k - m * m;
    term = std::sqrt((d__1 = b * b - a * 4. * c__, std::abs(d__1)));
    if (pi - *phi < 0.) {
	term = -term;
    }
    *w = (-b + term) / (a * 2.);
    return 0;
} /* tp364c_ */


/* Subroutine */ int tp365_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;
    static Real p, q;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 7;
    l1_1.nili = 0;
    l1_1.ninl = 5;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_6.x[0] = 3.;
    l2_6.x[1] = 0.;
    l2_6.x[2] = 2.;
    l2_6.x[3] = -1.5;
    l2_6.x[4] = 1.5;
    l2_6.x[5] = 5.;
    l2_6.x[6] = 0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	l11_6.lxl[(i__ << 1) - 1] = false;
/* labelL11: */
	l11_6.lxl[(i__ << 1) - 2] = true;
    }
    l11_6.lxl[6] = true;
    for (i__ = 1; i__ <= 7; ++i__) {
/* labelL12: */
	l12_6.lxu[i__ - 1] = false;
    }
    l13_6.xl[0] = 0.;
    l13_6.xl[2] = 0.;
    l13_6.xl[4] = 1.;
    l13_6.xl[6] = 1.;
    l20_8.lex = false;
    l20_8.nex = 1;
    l20_8.fex = 23.313708;
    l20_8.xex[0] = 4.8284266;
    l20_8.xex[1] = 4.7529555e-6;
    l20_8.xex[2] = 4.8284276;
    l20_8.xex[3] = 1.0000024;
    l20_8.xex[4] = 2.4142144;
    l20_8.xex[5] = 2.4142151;
    l20_8.xex[6] = 1.;
    return 0;
labelL2:
    l6_1.fx = l2_6.x[0] * l2_6.x[2];
labelL3:
    return 0;
labelL4:
/* Computing 2nd power */
    d__1 = l2_6.x[1];
/* Computing 2nd power */
    d__2 = l2_6.x[2];
    p = std::sqrt(d__1 * d__1 + d__2 * d__2);
/* Computing 2nd power */
    d__1 = l2_6.x[2];
/* Computing 2nd power */
    d__2 = l2_6.x[1] - l2_6.x[0];
    q = std::sqrt(d__1 * d__1 + d__2 * d__2);
    if (l9_5.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_6.x[3] - l2_6.x[5];
/* Computing 2nd power */
	d__2 = l2_6.x[4] - l2_6.x[6];
	l3_4.g[0] = d__1 * d__1 + d__2 * d__2 - 4.;
    }
    if (l9_5.index1[1]) {
	l3_4.g[1] = (l2_6.x[2] * l2_6.x[3] - l2_6.x[1] * l2_6.x[4]) / p - 1.;
    }
    if (l9_5.index1[2]) {
	l3_4.g[2] = (l2_6.x[2] * l2_6.x[5] - l2_6.x[1] * l2_6.x[6]) / p - 1.;
    }
    if (l9_5.index1[3]) {
	l3_4.g[3] = (l2_6.x[0] * l2_6.x[2] + (l2_6.x[1] - l2_6.x[0]) * l2_6.x[
		4] - l2_6.x[2] * l2_6.x[3]) / q - 1.;
    }
    if (l9_5.index1[4]) {
	l3_4.g[4] = (l2_6.x[0] * l2_6.x[2] + (l2_6.x[1] - l2_6.x[0]) * l2_6.x[
		6] - l2_6.x[2] * l2_6.x[5]) / q - 1.;
    }
labelL5:
    return 0;
} /* tp365_ */


/* Subroutine */ int tp366_(int *mode)
{
    /* Initialized data */

    static Real c__[38] = { 5.9553571e-4,.88392857,-.1175625,1.1088,
	    .1303533,-.0066033,6.6173269e-4,.017239878,-.0056595559,
	    -.019120592,56.85075,1.08702,.32175,-.03762,.006198,2462.3121,
	    -25.125634,161.18996,5e3,-489510.,44.333333,.33,.022556,-.007595,
	    6.1e-4,-5e-4,.819672,.819672,24500.,-250.,.010204082,1.2244898e-5,
	    6.25e-5,6.25e-5,-7.625e-5,1.22,1.,-1. };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 7;
    l1_1.nili = 0;
    l1_1.ninl = 14;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_6.x[0] = 1745.;
    l2_6.x[1] = 110.;
    l2_6.x[2] = 3048.;
    l2_6.x[3] = 89.;
    l2_6.x[4] = 92.8;
    l2_6.x[5] = 8.;
    l2_6.x[6] = 145.;
    for (i__ = 1; i__ <= 7; ++i__) {
	l11_6.lxl[i__ - 1] = true;
/* labelL11: */
	l12_6.lxu[i__ - 1] = true;
    }
    l13_6.xl[0] = 1.;
    l13_6.xl[1] = 1.;
    l13_6.xl[2] = 1.;
    l13_6.xl[3] = 85.;
    l13_6.xl[4] = 90.;
    l13_6.xl[5] = 3.;
    l13_6.xl[6] = 145.;
    l14_6.xu[0] = 2e3;
    l14_6.xu[1] = 120.;
    l14_6.xu[2] = 5e3;
    l14_6.xu[3] = 93.;
    l14_6.xu[4] = 95.;
    l14_6.xu[5] = 12.;
    l14_6.xu[6] = 162.;
    l20_8.lex = false;
    l20_8.nex = 1;
    l20_8.fex = 704.3056;
    l20_8.xex[0] = 905.40351;
    l20_8.xex[1] = 36.394998;
    l20_8.xex[2] = 2381.4783;
    l20_8.xex[3] = 88.987691;
    l20_8.xex[4] = 95.;
    l20_8.xex[5] = 12.;
    l20_8.xex[6] = 153.53535;
    return 0;
labelL2:
    l6_1.fx = l2_6.x[0] * 1.715 + l2_6.x[0] * .035 * l2_6.x[5] + l2_6.x[2] * 
	    4.0565 + l2_6.x[1] * 10. + 3e3 - l2_6.x[2] * .063 * l2_6.x[4];
labelL3:
    return 0;
labelL4:
    if (l9_8.index1[0]) {
/* Computing 2nd power */
	d__1 = l2_6.x[5];
	l3_7.g[0] = 1. - c__[0] * (d__1 * d__1) - c__[1] * l2_6.x[2] / l2_6.x[
		0] - c__[2] * l2_6.x[5];
    }
    if (l9_8.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_6.x[5];
	l3_7.g[1] = 1. - c__[3] * l2_6.x[0] / l2_6.x[2] - c__[4] * l2_6.x[0] /
		 l2_6.x[2] * l2_6.x[5] - c__[5] * l2_6.x[0] / l2_6.x[2] * (
		d__1 * d__1);
    }
    if (l9_8.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_6.x[5];
	l3_7.g[2] = 1. - c__[6] * (d__1 * d__1) - c__[7] * l2_6.x[4] - c__[8] 
		* l2_6.x[3] - c__[9] * l2_6.x[5];
    }
    if (l9_8.index1[3]) {
/* Computing 2nd power */
	d__1 = l2_6.x[5];
	l3_7.g[3] = 1. - c__[10] / l2_6.x[4] - c__[11] / l2_6.x[4] * l2_6.x[5]
		 - c__[12] * l2_6.x[3] / l2_6.x[4] - c__[13] / l2_6.x[4] * (
		d__1 * d__1);
    }
    if (l9_8.index1[4]) {
	l3_7.g[4] = 1. - c__[14] * l2_6.x[6] - c__[15] * l2_6.x[1] / l2_6.x[2]
		 / l2_6.x[3] - c__[16] * l2_6.x[1] / l2_6.x[2];
    }
    if (l9_8.index1[5]) {
	l3_7.g[5] = 1. - c__[17] / l2_6.x[6] - c__[18] * l2_6.x[1] / l2_6.x[2]
		 / l2_6.x[6] - c__[19] * l2_6.x[1] / l2_6.x[2] / l2_6.x[3] / 
		l2_6.x[6];
    }
    if (l9_8.index1[6]) {
	l3_7.g[6] = 1. - c__[20] / l2_6.x[4] - c__[21] * l2_6.x[6] / l2_6.x[4]
		;
    }
    if (l9_8.index1[7]) {
	l3_7.g[7] = 1. - c__[22] * l2_6.x[4] - c__[23] * l2_6.x[6];
    }
    if (l9_8.index1[8]) {
	l3_7.g[8] = 1. - c__[24] * l2_6.x[2] - c__[25] * l2_6.x[0];
    }
    if (l9_8.index1[9]) {
	l3_7.g[9] = 1. - c__[26] * l2_6.x[0] / l2_6.x[2] - c__[27] / l2_6.x[2]
		;
    }
    if (l9_8.index1[10]) {
	l3_7.g[10] = 1. - c__[28] * l2_6.x[1] / l2_6.x[2] / l2_6.x[3] - c__[
		29] * l2_6.x[1] / l2_6.x[2];
    }
    if (l9_8.index1[11]) {
	l3_7.g[11] = 1. - c__[30] * l2_6.x[3] - c__[31] / l2_6.x[1] * l2_6.x[
		2] * l2_6.x[3];
    }
    if (l9_8.index1[12]) {
	l3_7.g[12] = 1. - c__[32] * l2_6.x[0] * l2_6.x[5] - c__[33] * l2_6.x[
		0] - c__[34] * l2_6.x[2];
    }
    if (l9_8.index1[13]) {
	l3_7.g[13] = 1. - c__[35] / l2_6.x[0] * l2_6.x[2] - c__[36] / l2_6.x[
		0] - c__[37] * l2_6.x[5];
    }
labelL5:
    return 0;
} /* tp366_ */


/* Subroutine */ int tp367_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 7;
    l1_1.nili = 2;
    l1_1.ninl = 1;
    l1_1.neli = 1;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 7; ++i__) {
	l2_6.x[i__ - 1] = .1;
	l11_6.lxl[i__ - 1] = true;
	l12_6.lxu[i__ - 1] = false;
/* labelL11: */
	l13_6.xl[i__ - 1] = 0.;
    }
    l20_8.lex = false;
    l20_8.nex = 1;
    l20_8.fex = -37.41296;
    l20_8.xex[0] = 1.4688103;
    l20_8.xex[1] = 1.9839711;
    l20_8.xex[2] = .35187754;
    l20_8.xex[3] = 1.1953411;
    l20_8.xex[4] = .56940029;
    l20_8.xex[5] = .78474478;
    l20_8.xex[6] = 1.4121216;
    for (i__ = 1; i__ <= 7; ++i__) {
/* labelL12: */
	l5_34.gg[i__ * 5 - 5] = -1.;
    }
    for (i__ = 1; i__ <= 4; ++i__) {
/* labelL13: */
	l5_34.gg[i__ * 5 - 4] = -1.;
    }
    for (i__ = 5; i__ <= 7; ++i__) {
/* labelL14: */
	l5_34.gg[i__ * 5 - 4] = 0.;
    }
    l5_34.gg[2] = -1.;
    l5_34.gg[7] = 0.;
    l5_34.gg[12] = -1.;
    l5_34.gg[17] = 0.;
    l5_34.gg[22] = -1.;
    l5_34.gg[3] = 0.;
    l5_34.gg[8] = 0.;
    l5_34.gg[13] = 0.;
    l5_34.gg[18] = 2.;
    l5_34.gg[23] = 1.;
    l5_34.gg[28] = .8;
    l5_34.gg[33] = 1.;
    l5_34.gg[4] = 0.;
    l5_34.gg[19] = 0.;
    l5_34.gg[34] = 0.;
    l4_6.gf[1] = -5.;
    l4_6.gf[3] = -6.;
    return 0;
labelL2:
    l6_1.fx = l2_6.x[0] * -5. - l2_6.x[1] * 5. - l2_6.x[2] * 4. - l2_6.x[0] * 
	    l2_6.x[2] - l2_6.x[3] * 6. - l2_6.x[4] * 5. / (l2_6.x[4] + 1.) - 
	    l2_6.x[5] * 8. / (l2_6.x[5] + 1.) - (1. - std::exp(-l2_6.x[6]) * 2. + 
	    std::exp(l2_6.x[6] * -2.)) * 10.;
    return 0;
labelL3:
    l4_6.gf[0] = -5. - l2_6.x[2];
    l4_6.gf[2] = -4. - l2_6.x[0];
/* Computing 2nd power */
    d__1 = l2_6.x[4] + 1.;
    l4_6.gf[4] = -5. / (d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_6.x[5] + 1.;
    l4_6.gf[5] = -8. / (d__1 * d__1);
    l4_6.gf[6] = (exp(-l2_6.x[6]) - std::exp(l2_6.x[6] * -2.)) * -20.;
    return 0;
labelL4:
    if (l9_5.index1[0]) {
	l3_4.g[0] = 10. - l2_6.x[0] - l2_6.x[1] - l2_6.x[2] - l2_6.x[3] - 
		l2_6.x[4] - l2_6.x[5] - l2_6.x[6];
    }
    if (l9_5.index1[1]) {
	l3_4.g[1] = 5. - l2_6.x[0] - l2_6.x[1] - l2_6.x[2] - l2_6.x[3];
    }
    if (l9_5.index1[2]) {
/* Computing 2nd power */
	d__1 = l2_6.x[5];
/* Computing 2nd power */
	d__2 = l2_6.x[6];
	l3_4.g[2] = 5. - l2_6.x[0] - l2_6.x[2] - l2_6.x[4] - d__1 * d__1 - 
		d__2 * d__2;
    }
    if (l9_5.index1[3]) {
	l3_4.g[3] = l2_6.x[3] * 2. + l2_6.x[4] + l2_6.x[5] * .8 + l2_6.x[6] - 
		5.;
    }
    if (l9_5.index1[4]) {
/* Computing 2nd power */
	d__1 = l2_6.x[1];
/* Computing 2nd power */
	d__2 = l2_6.x[2];
/* Computing 2nd power */
	d__3 = l2_6.x[4];
/* Computing 2nd power */
	d__4 = l2_6.x[5];
	l3_4.g[4] = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4 - 
		5.;
    }
    return 0;
labelL5:
    if (! l10_5.index2[2]) {
	goto L51;
    }
    l5_34.gg[27] = l2_6.x[5] * -2.;
    l5_34.gg[32] = l2_6.x[6] * -2.;
L51:
    if (! l10_5.index2[4]) {
	goto L52;
    }
    l5_34.gg[9] = l2_6.x[1] * 2.;
    l5_34.gg[14] = l2_6.x[2] * 2.;
    l5_34.gg[24] = l2_6.x[4] * 2.;
    l5_34.gg[29] = l2_6.x[5] * 2.;
L52:
    return 0;
} /* tp367_ */


/* Subroutine */ int tp368_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;
    static Real s2, s3, s4;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 8;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 8; ++i__) {
	l2_7.x[i__ - 1] = (float)1. - (float)1. / (Real) i__;
	l11_7.lxl[i__ - 1] = true;
	l12_7.lxu[i__ - 1] = true;
	l13_7.xl[i__ - 1] = 0.;
/* labelL11: */
	l14_7.xu[i__ - 1] = 1.;
    }
    l2_7.x[6] = (float).7;
    l2_7.x[7] = (float).7;
    l20_2.lex = false;
    l20_2.nex = 1;
    l20_2.fex = -.74997564;
    l20_2.xex[0] = .49834105;
    l20_2.xex[1] = .4997795;
    l20_2.xex[2] = .50201378;
    l20_2.xex[3] = .50378302;
    l20_2.xex[4] = .50263008;
    l20_2.xex[5] = .50232579;
    l20_2.xex[6] = 1.;
    l20_2.xex[7] = 1.;
    l20_2.lex = true;
    l20_2.nex = 1;
    l20_2.fex = (float)-.75;
    l20_2.xex[0] = (float).5;
    l20_2.xex[1] = (float).5;
    l20_2.xex[2] = (float).5;
    l20_2.xex[3] = (float).5;
    l20_2.xex[4] = (float).5;
    l20_2.xex[5] = (float).5;
    l20_2.xex[6] = 1.;
    l20_2.xex[7] = 1.;
    return 0;
labelL2:
    s2 = 0.;
    s3 = 0.;
    s4 = 0.;
    for (i__ = 1; i__ <= 8; ++i__) {
/* Computing 2nd power */
	d__1 = l2_7.x[i__ - 1];
	s2 += d__1 * d__1;
/* Computing 3rd power */
	d__1 = l2_7.x[i__ - 1];
	s3 += d__1 * (d__1 * d__1);
/* labelL10: */
/* Computing 4th power */
	d__1 = l2_7.x[i__ - 1], d__1 *= d__1;
	s4 += d__1 * d__1;
    }
    if (*mode == 3) {
	goto labelL3;
    }
/* Computing 2nd power */
    d__1 = s3;
    l6_1.fx = -s2 * s4 + d__1 * d__1;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 8; ++i__) {
/* L31: */
/* Computing 3rd power */
	d__1 = l2_7.x[i__ - 1];
/* Computing 2nd power */
	d__2 = l2_7.x[i__ - 1];
	l4_7.gf[i__ - 1] = l2_7.x[i__ - 1] * -2. * s4 - d__1 * (d__1 * d__1) *
		 4. * s2 + d__2 * d__2 * 6. * s3;
    }
labelL4:
    return 0;
} /* tp368_ */


/* Subroutine */ int tp369_(int *mode)
{
    /* Initialized data */

    static Real c__[16] = { 833.33252,100.,-83333.333,1250.,1.,-1250.,
	    1.25e6,1.,-2500.,.0025,.0025,.0025,.0025,-.0025,.01,-.01 };

    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 8;
    l1_1.nili = 3;
    l1_1.ninl = 3;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_7.x[0] = 5e3;
    l2_7.x[1] = 5e3;
    l2_7.x[2] = 5e3;
    l2_7.x[3] = 200.;
    l2_7.x[4] = 350.;
    l2_7.x[5] = 150.;
    l2_7.x[6] = 225.;
    l2_7.x[7] = 425.;
    for (i__ = 1; i__ <= 8; ++i__) {
	l12_7.lxu[i__ - 1] = true;
/* labelL11: */
	l11_7.lxl[i__ - 1] = true;
    }
    l13_7.xl[0] = 100.;
    l13_7.xl[1] = 1e3;
    l13_7.xl[2] = 1e3;
    for (i__ = 1; i__ <= 3; ++i__) {
/* labelL12: */
	l14_7.xu[i__ - 1] = 1e4;
    }
    for (i__ = 4; i__ <= 8; ++i__) {
	l13_7.xl[i__ - 1] = 10.;
/* labelL13: */
	l14_7.xu[i__ - 1] = 1e3;
    }
    l20_2.lex = false;
    l20_2.nex = 1;
    l20_2.fex = 7049.248;
    l20_2.xex[0] = 579.30657;
    l20_2.xex[1] = 1359.9705;
    l20_2.xex[2] = 5109.9709;
    l20_2.xex[3] = 182.01769;
    l20_2.xex[4] = 295.60116;
    l20_2.xex[5] = 217.98231;
    l20_2.xex[6] = 286.41653;
    l20_2.xex[7] = 395.60116;
    return 0;
labelL2:
    l6_1.fx = l2_7.x[0] + l2_7.x[1] + l2_7.x[2];
labelL3:
    return 0;
labelL4:
    if (l9_6.index1[0]) {
	l3_5.g[0] = (float)1. - c__[9] * l2_7.x[3] - c__[10] * l2_7.x[5];
    }
    if (l9_6.index1[1]) {
	l3_5.g[1] = (float)1. - c__[11] * l2_7.x[4] - c__[12] * l2_7.x[6] - 
		c__[13] * l2_7.x[3];
    }
    if (l9_6.index1[2]) {
	l3_5.g[2] = (float)1. - c__[14] * l2_7.x[7] - c__[15] * l2_7.x[4];
    }
    if (l9_6.index1[3]) {
	l3_5.g[3] = (float)1. - c__[0] / l2_7.x[0] * l2_7.x[3] / l2_7.x[5] - 
		c__[1] / l2_7.x[5] - c__[2] / l2_7.x[0] / l2_7.x[5];
    }
    if (l9_6.index1[4]) {
	l3_5.g[4] = (float)1. - c__[3] / l2_7.x[1] * l2_7.x[4] / l2_7.x[6] - 
		c__[4] * l2_7.x[3] / l2_7.x[6] - c__[5] / l2_7.x[1] * l2_7.x[
		3] / l2_7.x[6];
    }
    if (l9_6.index1[5]) {
	l3_5.g[5] = (float)1. - c__[6] / l2_7.x[2] / l2_7.x[7] - c__[7] * 
		l2_7.x[4] / l2_7.x[7] - c__[8] / l2_7.x[2] * l2_7.x[4] / 
		l2_7.x[7];
    }
labelL5:
    return 0;
} /* tp369_ */


/* Subroutine */ int tp370_0_(int n__, int *mode)
{
    /* System generated locals */
    int i__1, i__2;
    Real d__1;

    /* Builtin functions */
    Real pow_di(Real*, int *);

    /* Local variables */
    static int i__, j;
    static Real scale;
    static Real basis, sum, sum1, sum2[31];

    switch(n__) {
	case 1: goto L_tp371;
	}

    l1_1.n = 6;
    scale = (float)1.;
    l20_9.fex = .00228767005355;
    l20_9.xex[0] = -.015725;
    l20_9.xex[1] = 1.012435;
    l20_9.xex[2] = -.232992;
    l20_9.xex[3] = -1.26043;
    l20_9.xex[4] = -1.513729;
    l20_9.xex[5] = .992996;
    goto labelL10;

L_tp371:
    l1_1.n = 9;
    scale = (float)1e5;
    l20_9.fex = scale * 1.3997660138e-6;
    l20_9.xex[0] = 1.5e-5;
    l20_9.xex[1] = .99979;
    l20_9.xex[2] = -.014764;
    l20_9.xex[3] = -.146342;
    l20_9.xex[4] = 1.000821;
    l20_9.xex[5] = -2.617731;
    l20_9.xex[6] = 4.104403;
    l20_9.xex[7] = -2.143612;
    l20_9.xex[8] = 1.052627;
labelL10:
    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL4;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    i__1 = l1_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l2_8.x[i__ - 1] = 0.;
	l11_8.lxl[i__ - 1] = false;
/* labelL11: */
	l12_8.lxu[i__ - 1] = false;
    }
    l15_1.lsum = 31;
    l20_9.lex = false;
    l20_9.nex = 1;
    return 0;
labelL2:
    l16_18.f[0] = l2_8.x[0];
/* Computing 2nd power */
    d__1 = l2_8.x[0];
    l16_18.f[1] = l2_8.x[1] - d__1 * d__1 - 1.;
    for (i__ = 2; i__ <= 30; ++i__) {
	basis = (Real) (i__ - 1) / 29.;
	sum1 = 0.;
	i__1 = l1_1.n;
	for (j = 2; j <= i__1; ++j) {
/* L21: */
	    i__2 = j - 2;
	    sum1 += l2_8.x[j - 1] * (Real) (j - 1) * pow_di(&basis, &
		    i__2);
	}
	sum = 0.;
	i__2 = l1_1.n;
	for (j = 1; j <= i__2; ++j) {
/* L23: */
	    i__1 = j - 1;
	    sum += l2_8.x[j - 1] * pow_di(&basis, &i__1);
	}
	sum2[i__] = sum;
/* labelL20: */
/* Computing 2nd power */
	d__1 = sum2[i__];
	l16_18.f[i__] = sum1 - d__1 * d__1 - 1.;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 31; ++i__) {
/* L22: */
/* Computing 2nd power */
	d__1 = l16_18.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    l6_1.fx *= scale;
labelL4:
    return 0;
} /* tp370_ */

/* Subroutine */ int tp370_(int *mode)
{
    return tp370_0_(0, mode);
    }

/* Subroutine */ int tp371_(int *mode)
{
    return tp370_0_(1, mode);
    }


/* Subroutine */ int tp372_(int *mode)
{
    /* System generated locals */
    Real d__1;


    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 9;
    l1_1.nili = 0;
    l1_1.ninl = 12;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_8.x[0] = 300.;
    l2_8.x[1] = -100.;
    l2_8.x[2] = -.1997;
    l2_8.x[3] = -127.;
    l2_8.x[4] = -151.;
    l2_8.x[5] = 379.;
    l2_8.x[6] = 421.;
    l2_8.x[7] = 460.;
    l2_8.x[8] = 426.;
    for (i__ = 1; i__ <= 12; ++i__) {
	for (j = 4; j <= 9; ++j) {
/* labelL12: */
	    l5_36.gg[i__ + j * 12 - 13] = 0.;
	}
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	l5_36.gg[i__ - 1] = 1.;
/* L18: */
	l5_36.gg[i__ + (i__ + 3) * 12 - 13] = 1.;
    }
    for (i__ = 7; i__ <= 12; ++i__) {
	l5_36.gg[i__ - 1] = -1.;
/* labelL13: */
	l5_36.gg[i__ + (i__ - 3) * 12 - 13] = 1.;
    }
    l4_8.gf[0] = 0.;
    l4_8.gf[1] = 0.;
    l4_8.gf[2] = 0.;
    for (i__ = 4; i__ <= 9; ++i__) {
	l11_8.lxl[i__ - 1] = true;
	l12_8.lxu[i__ - 1] = false;
/* labelL14: */
	l13_8.xl[i__ - 1] = 0.;
    }
    for (i__ = 1; i__ <= 2; ++i__) {
	l11_8.lxl[i__ - 1] = false;
/* L15: */
	l12_8.lxu[i__ - 1] = false;
    }
    l13_8.xl[2] = (float)-1.;
    l14_8.xu[2] = (float)0.;
    l11_8.lxl[2] = true;
    l12_8.lxu[2] = true;
    l20_9.xex[0] = 523.30555;
    l20_9.xex[1] = -156.94787;
    l20_9.xex[2] = -.19966457;
    l20_9.xex[3] = 29.608067;
    l20_9.xex[4] = 86.615521;
    l20_9.xex[5] = 47.326718;
    l20_9.xex[6] = 26.235604;
    l20_9.xex[7] = 22.915985;
    l20_9.xex[8] = 39.470742;
    l20_9.lex = false;
    l20_9.nex = 1;
    l20_9.fex = 13390.093;
    l15_1.lsum = 6;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 4; i__ <= 9; ++i__) {
	l16_8.f[i__ - 4] = l2_8.x[i__ - 1];
/* L21: */
/* Computing 2nd power */
	d__1 = l2_8.x[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (i__ = 4; i__ <= 9; ++i__) {
/* L35: */
	l16_8.f[i__ - 4] = l2_8.x[i__ - 1];
    }
    for (i__ = 4; i__ <= 9; ++i__) {
/* L31: */
	l4_8.gf[i__ - 1] = l2_8.x[i__ - 1] * 2.;
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	for (j = 1; j <= 9; ++j) {
/* L33: */
	    l17_22.df[i__ + j * 6 - 7] = 0.;
	}
    }
    for (i__ = 1; i__ <= 6; ++i__) {
/* L34: */
	l17_22.df[i__ + (i__ + 3) * 6 - 7] = 1.;
    }
    return 0;
labelL4:
    if (l9_17.index1[0]) {
	l3_16.g[0] = l2_8.x[0] + l2_8.x[1] * std::exp(l2_8.x[2] * -5.) + l2_8.x[3] 
		- 127.;
    }
    if (l9_17.index1[1]) {
	l3_16.g[1] = l2_8.x[0] + l2_8.x[1] * std::exp(l2_8.x[2] * -3.) + l2_8.x[4] 
		- 151.;
    }
    if (l9_17.index1[2]) {
	l3_16.g[2] = l2_8.x[0] + l2_8.x[1] * std::exp(-l2_8.x[2]) + l2_8.x[5] - 
		379.;
    }
    if (l9_17.index1[3]) {
	l3_16.g[3] = l2_8.x[0] + l2_8.x[1] * std::exp(l2_8.x[2]) + l2_8.x[6] - 
		421.;
    }
    if (l9_17.index1[4]) {
	l3_16.g[4] = l2_8.x[0] + l2_8.x[1] * std::exp(l2_8.x[2] * 3.) + l2_8.x[7] 
		- 460.;
    }
    if (l9_17.index1[5]) {
	l3_16.g[5] = l2_8.x[0] + l2_8.x[1] * std::exp(l2_8.x[2] * 5.) + l2_8.x[8] 
		- 426.;
    }
    if (l9_17.index1[6]) {
	l3_16.g[6] = -l2_8.x[0] - l2_8.x[1] * std::exp(l2_8.x[2] * -5.) + l2_8.x[3]
		 + 127.;
    }
    if (l9_17.index1[7]) {
	l3_16.g[7] = -l2_8.x[0] - l2_8.x[1] * std::exp(l2_8.x[2] * -3.) + l2_8.x[4]
		 + 151.;
    }
    if (l9_17.index1[8]) {
	l3_16.g[8] = -l2_8.x[0] - l2_8.x[1] * std::exp(-l2_8.x[2]) + l2_8.x[5] + 
		379.;
    }
    if (l9_17.index1[9]) {
	l3_16.g[9] = -l2_8.x[0] - l2_8.x[1] * std::exp(l2_8.x[2]) + l2_8.x[6] + 
		421.;
    }
    if (l9_17.index1[10]) {
	l3_16.g[10] = -l2_8.x[0] - l2_8.x[1] * std::exp(l2_8.x[2] * 3.) + l2_8.x[7]
		 + 460.;
    }
    if (l9_17.index1[11]) {
	l3_16.g[11] = -l2_8.x[0] - l2_8.x[1] * std::exp(l2_8.x[2] * 5.) + l2_8.x[8]
		 + 426.;
    }
    return 0;
labelL5:
    for (i__ = 1; i__ <= 6; ++i__) {
	if (! l10_16.index2[i__ + 5]) {
	    goto L51;
	}
	l5_36.gg[i__ + 17] = -exp((Real) ((i__ << 1) - 7) * l2_8.x[2]);
	l5_36.gg[i__ + 29] = -l2_8.x[1] * (Real) ((i__ << 1) - 7) * 
		l5_36.gg[i__ + 17];
L51:
	;
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	if (! l10_16.index2[i__ - 1]) {
	    goto L52;
	}
	l5_36.gg[i__ + 11] = std::exp((Real) ((i__ << 1) - 7) * l2_8.x[2]);
	l5_36.gg[i__ + 23] = l2_8.x[1] * (Real) ((i__ << 1) - 7) * 
		l5_36.gg[i__ + 11];
L52:
	;
    }
    return 0;
} /* tp372_ */


/* Subroutine */ int tp373_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 6;
    l1_1.n = 9;
    for (i__ = 1; i__ <= 9; ++i__) {
	l11_8.lxl[i__ - 1] = false;
/* labelL11: */
	l12_8.lxu[i__ - 1] = false;
    }
    l11_8.lxl[2] = true;
    l12_8.lxu[2] = true;
    l13_8.xl[2] = (float)-1.;
    l14_8.xu[2] = (float)0.;
    l2_8.x[0] = 300.;
    l2_8.x[1] = -100.;
    l2_8.x[2] = -.1997;
    l2_8.x[3] = -127.;
    l2_8.x[4] = -151.;
    l2_8.x[5] = 379.;
    l2_8.x[6] = 421.;
    l2_8.x[7] = 460.;
    l2_8.x[8] = 426.;
    l4_8.gf[0] = 0.;
    l4_8.gf[1] = 0.;
    l4_8.gf[2] = 0.;
    for (i__ = 1; i__ <= 6; ++i__) {
	for (j = 4; j <= 9; ++j) {
/* L17: */
	    l5_19.gg[i__ + j * 6 - 7] = 0.;
	}
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	l5_19.gg[i__ - 1] = 1.;
/* L18: */
	l5_19.gg[i__ + (i__ + 3) * 6 - 7] = 1.;
    }
    l20_9.lex = false;
    l20_9.nex = 1;
    l20_9.fex = 13390.093;
    l20_9.xex[0] = 523.30542;
    l20_9.xex[1] = -156.9477;
    l20_9.xex[2] = -.19966472;
    l20_9.xex[3] = 29.608061;
    l20_9.xex[4] = -86.615571;
    l20_9.xex[5] = 47.326669;
    l20_9.xex[6] = 26.235575;
    l20_9.xex[7] = 22.915982;
    l20_9.xex[8] = -39.470718;
    l15_1.lsum = 6;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 4; i__ <= 9; ++i__) {
	l16_8.f[i__ - 4] = l2_8.x[i__ - 1];
/* L21: */
/* Computing 2nd power */
	d__1 = l2_8.x[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (i__ = 4; i__ <= 9; ++i__) {
/* L35: */
	l16_8.f[i__ - 4] = l2_8.x[i__ - 1];
    }
    for (i__ = 4; i__ <= 9; ++i__) {
/* L31: */
	l4_8.gf[i__ - 1] = l2_8.x[i__ - 1] * 2.;
    }
    for (i__ = 1; i__ <= 6; ++i__) {
	for (j = 1; j <= 9; ++j) {
/* L33: */
	    l17_22.df[i__ + j * 6 - 7] = 0.;
	}
    }
    for (i__ = 1; i__ <= 6; ++i__) {
/* L34: */
	l17_22.df[i__ + (i__ + 3) * 6 - 7] = 1.;
    }
    return 0;
labelL4:
    if (l9_6.index1[0]) {
	l3_5.g[0] = l2_8.x[0] + l2_8.x[1] * std::exp(l2_8.x[2] * -5.) + l2_8.x[3] 
		- 127.;
    }
    if (l9_6.index1[1]) {
	l3_5.g[1] = l2_8.x[0] + l2_8.x[1] * std::exp(l2_8.x[2] * -3.) + l2_8.x[4] 
		- 151.;
    }
    if (l9_6.index1[2]) {
	l3_5.g[2] = l2_8.x[0] + l2_8.x[1] * std::exp(-l2_8.x[2]) + l2_8.x[5] - 
		379.;
    }
    if (l9_6.index1[3]) {
	l3_5.g[3] = l2_8.x[0] + l2_8.x[1] * std::exp(l2_8.x[2]) + l2_8.x[6] - 421.;
    }
    if (l9_6.index1[4]) {
	l3_5.g[4] = l2_8.x[0] + l2_8.x[1] * std::exp(l2_8.x[2] * 3.) + l2_8.x[7] - 
		460.;
    }
    if (l9_6.index1[5]) {
	l3_5.g[5] = l2_8.x[0] + l2_8.x[1] * std::exp(l2_8.x[2] * 5.) + l2_8.x[8] - 
		426.;
    }
    return 0;
labelL5:
    for (i__ = 1; i__ <= 6; ++i__) {
	if (! l10_6.index2[i__ - 1]) {
	    goto L52;
	}
	l5_19.gg[i__ + 5] = std::exp((Real) ((i__ << 1) - 7) * l2_8.x[2]);
	l5_19.gg[i__ + 11] = l2_8.x[1] * (Real) ((i__ << 1) - 7) * 
		l5_19.gg[i__ + 5];
L52:
	;
    }
    return 0;
} /* tp373_ */


/* Subroutine */ int tp374_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    Real tp374a_(Real *, Real*), tp374b_(Real *, Real *), tp374g_(Real *, Real *);
    static int i__, k;
    static Real z__, pi;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 10;
    l1_1.nili = 0;
    l1_1.ninl = 35;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 10; ++i__) {
	l2_9.x[i__ - 1] = .1;
	l12_9.lxu[i__ - 1] = false;
/* labelL6: */
	l11_9.lxl[i__ - 1] = false;
    }
    l20_14.lex = false;
    l20_14.nex = 2;
    l20_14.fex = .233264;
    l20_14.xex[0] = .218212;
    l20_14.xex[1] = .23264;
    l20_14.xex[2] = .278457;
    l20_14.xex[3] = .268125;
    l20_14.xex[4] = .21201;
    l20_14.xex[5] = .125918;
    l20_14.xex[6] = .034102;
    l20_14.xex[7] = -.026136;
    l20_14.xex[8] = -.142233;
    l20_14.xex[9] = .233264;
    l20_14.xex[10] = -.142233;
    l20_14.xex[11] = -.026136;
    l20_14.xex[12] = .034102;
    l20_14.xex[13] = .125918;
    l20_14.xex[14] = .21201;
    l20_14.xex[15] = .268125;
    l20_14.xex[16] = .278457;
    l20_14.xex[17] = .23264;
    l20_14.xex[18] = .218212;
    l20_14.xex[19] = .233264;
    for (i__ = 1; i__ <= 9; ++i__) {
/* L46: */
	l4_9.gf[i__ - 1] = 0.;
    }
    l4_9.gf[9] = 1.;
    return 0;
labelL2:
    l6_1.fx = l2_9.x[9];
labelL3:
    return 0;
labelL4:
    pi = std::atan(1.) * 4.;
    for (i__ = 1; i__ <= 10; ++i__) {
	z__ = pi / 4. * ((Real) (i__ - 1) * .1);
/* L8: */
	if (l9_18.index1[i__ - 1]) {
/* Computing 2nd power */
	    d__1 = 1. - l2_9.x[9];
	    l3_17.g[i__ - 1] = tp374g_(&z__, l2_9.x) - d__1 * d__1;
	}
    }
    for (i__ = 11; i__ <= 20; ++i__) {
	z__ = pi / 4. * ((Real) (i__ - 11) * .1);
/* labelL9: */
	if (l9_18.index1[i__ - 1]) {
/* Computing 2nd power */
	    d__1 = l2_9.x[9] + 1.;
	    l3_17.g[i__ - 1] = d__1 * d__1 - tp374g_(&z__, l2_9.x);
	}
    }
    for (i__ = 21; i__ <= 35; ++i__) {
	z__ = pi / 4. * ((Real) (i__ - 21) * .2 + 1.2);
/* labelL10: */
	if (l9_18.index1[i__ - 1]) {
/* Computing 2nd power */
	    d__1 = l2_9.x[9];
	    l3_17.g[i__ - 1] = d__1 * d__1 - tp374g_(&z__, l2_9.x);
	}
    }
    return 0;
labelL5:
    pi = std::atan(1.) * 4.;
    for (i__ = 1; i__ <= 10; ++i__) {
	if (! l10_17.index2[i__ - 1]) {
	    goto L50;
	}
	z__ = pi / 4. * ((Real) (i__ - 1) * .1);
	for (k = 1; k <= 9; ++k) {
	    l5_37.gg[i__ + k * 35 - 36] = (tp374a_(&z__, l2_9.x) *std::cos(k * 
		    z__) + tp374b_(&z__, l2_9.x) * std::sin(k * z__)) * 2.;
/* L51: */
	    l5_37.gg[i__ + 314] = (1. - l2_9.x[9]) * 2.;
	}
L50:
	;
    }
    for (i__ = 11; i__ <= 20; ++i__) {
	if (! l10_17.index2[i__ - 1]) {
	    goto L52;
	}
	z__ = pi / 4. * ((Real) (i__ - 11) * .1);
	for (k = 1; k <= 9; ++k) {
	    l5_37.gg[i__ + k * 35 - 36] = (tp374a_(&z__, l2_9.x) *std::cos(k * 
		    z__) + tp374b_(&z__, l2_9.x) * std::sin(k * z__)) * -2.;
/* L53: */
	    l5_37.gg[i__ + 314] = (l2_9.x[9] + 1.) * 2.;
	}
L52:
	;
    }
    for (i__ = 21; i__ <= 35; ++i__) {
	if (! l10_17.index2[i__ - 1]) {
	    goto L54;
	}
	z__ = pi / 4. * ((Real) (i__ - 21) * .2 + 1.2);
	for (k = 1; k <= 9; ++k) {
	    l5_37.gg[i__ + k * 35 - 36] = (tp374a_(&z__, l2_9.x) *std::cos(k * 
		    z__) + tp374b_(&z__, l2_9.x) * std::sin(k * z__)) * -2.;
/* L55: */
	    l5_37.gg[i__ + 314] = l2_9.x[9] * 2.;
	}
L54:
	;
    }
    return 0;
} /* tp374_ */

Real tp374g_(Real *a, Real *x)
{
    /* System generated locals */
    Real ret_val, d__1, d__2;

    /* Local variables */
    Real tp374a_(Real *, Real*), tp374b_(Real *, Real *);

    /* Parameter adjustments */
    --x;

    /* Function Body */
/* Computing 2nd power */
    d__1 = tp374a_(a, &x[1]);
/* Computing 2nd power */
    d__2 = tp374b_(a, &x[1]);
    ret_val = d__1 * d__1 + d__2 * d__2;
    return ret_val;
} /* tp374g_ */

Real tp374a_(Real *a, Real *x)
{
    /* System generated locals */
    Real ret_val;


    /* Local variables */
    static int k;

    /* Parameter adjustments */
    --x;

    /* Function Body */
    ret_val = 0.;
    for (k = 1; k <= 9; ++k) {
/* labelL10: */
	ret_val += x[k] *std::cos(k * *a);
    }
    return ret_val;
} /* tp374a_ */

Real tp374b_(Real *a, Real *x)
{
    /* System generated locals */
    Real ret_val;

    /* Local variables */
    static int k;

    /* Parameter adjustments */
    --x;

    /* Function Body */
    ret_val = 0.;
    for (k = 1; k <= 9; ++k) {
/* labelL10: */
	ret_val += x[k] * std::sin(k * *a);
    }
    return ret_val;
} /* tp374b_ */


/* Subroutine */ int tp375_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    Real tp375a_(int *, int *);
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 10;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 8;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 10; ++i__) {
	l2_9.x[i__ - 1] = 1.;
	l12_9.lxu[i__ - 1] = false;
/* labelL6: */
	l11_9.lxl[i__ - 1] = false;
    }
    l20_10.lex = false;
    l20_10.nex = 1;
    l20_10.fex = -15.16;
    for (i__ = 1; i__ <= 8; ++i__) {
/* L60: */
	l20_10.xex[i__ - 1] = .1064;
    }
    l20_10.xex[8] = 2.843;
    l20_10.xex[9] = -2.642;
    for (j = 1; j <= 8; ++j) {
	for (i__ = 1; i__ <= 10; ++i__) {
/* L7: */
	    l5_21.gg[j + i__ * 9 - 10] = 1.;
	}
    }
    for (i__ = 1; i__ <= 8; ++i__) {
/* L8: */
	l5_21.gg[i__ + i__ * 9 - 10] = .5;
    }
    l15_1.lsum = 10;
    for (i__ = 1; i__ <= 10; ++i__) {
	for (j = 1; j <= 10; ++j) {
	    l17_23.df[i__ + j * 10 - 11] = 0.;
/* L17: */
	    l17_23.df[i__ + i__ * 10 - 11] = -1.;
	}
    }
    return 0;
labelL2:
    for (i__ = 1; i__ <= 10; ++i__) {
/* L16: */
	l16_4.f[i__ - 1] = -l2_9.x[i__ - 1];
    }
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* labelL9: */
/* Computing 2nd power */
	d__1 = l2_9.x[i__ - 1];
	l6_1.fx -= d__1 * d__1;
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 10; ++i__) {
/* labelL10: */
	l4_9.gf[i__ - 1] = l2_9.x[i__ - 1] * -2.;
    }
labelL4:
    for (j = 1; j <= 8; ++j) {
	if (! l9_16.index1[j - 1]) {
	    goto labelL11;
	}
	l3_15.g[j - 1] = 0.;
	for (i__ = 1; i__ <= 10; ++i__) {
/* labelL12: */
	    l3_15.g[j - 1] += l2_9.x[i__ - 1] / tp375a_(&i__, &j);
	}
	l3_15.g[j - 1] += -1.;
labelL11:
	;
    }
    if (! l9_16.index1[8]) {
	goto L18;
    }
    l3_15.g[8] = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* labelL13: */
	l3_15.g[8] += pow_dd(&l2_9.x[i__ - 1], &c_b305) / ((Real) (i__ 
		- 1) / 3. + 1.);
    }
    l3_15.g[8] += -4.;
L18:
    return 0;
labelL5:
    if (! l10_18.index2[8]) {
	goto L15;
    }
    for (i__ = 1; i__ <= 10; ++i__) {
/* labelL14: */
	l5_21.gg[i__ * 9 - 1] = l2_9.x[i__ - 1] * 2. / ((Real) (i__ - 1)
		 / 3. + 1.);
    }
L15:
    return 0;
} /* tp375_ */

Real tp375a_(int *i__, int *j)
{
    /* System generated locals */
    Real ret_val;

    ret_val = 1.;
    if (*i__ == *j) {
	ret_val = 2.;
    }
    return ret_val;
} /* tp375a_ */


/* Subroutine */ int tp376_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 10;
    l1_1.nili = 0;
    l1_1.ninl = 14;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    l2_9.x[0] = 10.;
    l2_9.x[1] = .005;
    l2_9.x[2] = .0081;
    l2_9.x[3] = 100.;
    l2_9.x[4] = .0017;
    l2_9.x[5] = .0013;
    l2_9.x[6] = .0027;
    l2_9.x[7] = .002;
    l2_9.x[8] = .15;
    l2_9.x[9] = .105;
    for (i__ = 1; i__ <= 10; ++i__) {
	l12_9.lxu[i__ - 1] = true;
/* labelL6: */
	l11_9.lxl[i__ - 1] = true;
    }
    l13_9.xl[0] = 0.;
    l13_9.xl[1] = 0.;
    l13_9.xl[2] = 5e-5;
    l13_9.xl[3] = 10.;
    for (i__ = 5; i__ <= 10; ++i__) {
/* L7: */
	l13_9.xl[i__ - 1] = 5e-5;
    }
    l14_9.xu[0] = 10.;
    l14_9.xu[1] = .1;
    l14_9.xu[2] = .0081;
    l14_9.xu[3] = 1e3;
    l14_9.xu[4] = .0017;
    l14_9.xu[5] = .0013;
    l14_9.xu[6] = .0027;
    l14_9.xu[7] = .002;
    l14_9.xu[8] = 1.;
    l14_9.xu[9] = 1.;
    l5_32.gg[0] = 1.;
    l5_32.gg[15] = 0.;
    for (i__ = 5; i__ <= 10; ++i__) {
/* labelL9: */
	l5_32.gg[i__ * 15 - 15] = 0.;
    }
    l5_32.gg[1] = 1.;
    l5_32.gg[16] = 0.;
    l5_32.gg[31] = 0.;
    for (i__ = 6; i__ <= 8; ++i__) {
/* labelL10: */
	l5_32.gg[i__ * 15 - 14] = 0.;
    }
    l5_32.gg[136] = 0.;
    l5_32.gg[2] = 1.;
    l5_32.gg[17] = 0.;
    l5_32.gg[32] = 0.;
    l5_32.gg[62] = 0.;
    for (i__ = 7; i__ <= 9; ++i__) {
/* labelL11: */
	l5_32.gg[i__ * 15 - 13] = 0.;
    }
    l5_32.gg[3] = 1.;
    l5_32.gg[18] = 0.;
    l5_32.gg[33] = 0.;
    l5_32.gg[63] = 0.;
    l5_32.gg[78] = 0.;
    for (i__ = 8; i__ <= 10; ++i__) {
/* labelL12: */
	l5_32.gg[i__ * 15 - 12] = 0.;
    }
    l5_32.gg[4] = 1.;
    l5_32.gg[19] = 0.;
    l5_32.gg[34] = 0.;
    for (i__ = 5; i__ <= 7; ++i__) {
/* labelL13: */
	l5_32.gg[i__ * 15 - 11] = 0.;
    }
    l5_32.gg[124] = 0.;
    l5_32.gg[139] = 0.;
    l5_32.gg[5] = 0.;
    l5_32.gg[20] = 1e4;
    l5_32.gg[35] = 0.;
    for (i__ = 6; i__ <= 8; ++i__) {
/* labelL14: */
	l5_32.gg[i__ * 15 - 10] = 0.;
    }
    l5_32.gg[140] = 0.;
    l5_32.gg[6] = 0.;
    l5_32.gg[21] = 1e4;
    l5_32.gg[36] = 0.;
    l5_32.gg[66] = 0.;
    for (i__ = 7; i__ <= 9; ++i__) {
/* L15: */
	l5_32.gg[i__ * 15 - 9] = 0.;
    }
    l5_32.gg[7] = 0.;
    l5_32.gg[22] = 1e4;
    l5_32.gg[37] = 0.;
    l5_32.gg[67] = 0.;
    l5_32.gg[82] = 0.;
    for (i__ = 8; i__ <= 10; ++i__) {
/* L16: */
	l5_32.gg[i__ * 15 - 8] = 0.;
    }
    l5_32.gg[8] = 0.;
    l5_32.gg[23] = 1e4;
    l5_32.gg[38] = 0.;
    for (i__ = 5; i__ <= 7; ++i__) {
/* L17: */
	l5_32.gg[i__ * 15 - 7] = 0.;
    }
    l5_32.gg[128] = 0.;
    l5_32.gg[143] = 0.;
    l5_32.gg[9] = 0.;
    l5_32.gg[24] = 1e4;
    for (i__ = 5; i__ <= 10; ++i__) {
/* L18: */
	l5_32.gg[i__ * 15 - 6] = 0.;
    }
    l5_32.gg[10] = 0.;
    l5_32.gg[25] = 1e4;
    for (i__ = 5; i__ <= 10; ++i__) {
/* L19: */
	l5_32.gg[i__ * 15 - 5] = 0.;
    }
    l5_32.gg[11] = 0.;
    l5_32.gg[26] = 1e4;
    for (i__ = 5; i__ <= 10; ++i__) {
/* labelL20: */
	l5_32.gg[i__ * 15 - 4] = 0.;
    }
    l5_32.gg[12] = 0.;
    l5_32.gg[27] = 1e4;
    for (i__ = 5; i__ <= 10; ++i__) {
/* L21: */
	l5_32.gg[i__ * 15 - 3] = 0.;
    }
    l5_32.gg[13] = 0.;
    l5_32.gg[28] = 0.;
    l5_32.gg[103] = 0.;
    l5_32.gg[133] = 0.;
    l5_32.gg[148] = 0.;
    for (i__ = 1; i__ <= 8; ++i__) {
/* L8: */
	l5_32.gg[i__ * 15 - 1] = 0.;
    }
    l5_32.gg[134] = 1.;
    l5_32.gg[149] = 1.;
    l20_10.lex = false;
    l20_10.nex = 1;
    l20_10.fex = -4430.0879;
    l20_10.xex[0] = (float).14727222;
    l20_10.xex[1] = (float).1;
    l20_10.xex[2] = .0081;
    l20_10.xex[3] = 628.71731;
    l20_10.xex[4] = .0017;
    l20_10.xex[5] = .0011816143;
    l20_10.xex[6] = .0027;
    l20_10.xex[7] = .00135;
    l20_10.xex[8] = (float).15740741;
    l20_10.xex[9] = .097592593;
    for (i__ = 3; i__ <= 10; ++i__) {
/* L45: */
	l4_9.gf[i__ - 1] = 0.;
    }
    return 0;
labelL2:
    l6_1.fx = (l2_9.x[0] * .15 + l2_9.x[1] * 14. - .06) * 2e4 / (l2_9.x[0] + 
	    .002 + l2_9.x[1] * 60.);
    l6_1.fx = -l6_1.fx;
    return 0;
labelL3:
    d__1 = l2_9.x[0] + .002 + l2_9.x[1] * 60.;
    l4_9.gf[0] = ((l2_9.x[0] + .002 + l2_9.x[1] * 60.) * .15 - (l2_9.x[0] * 
	    .15 + l2_9.x[1] * 14. - .06)) * 2e4 / pow_dd(&d__1, &c_b305);
    d__1 = l2_9.x[0] + .002 + l2_9.x[1] * 60.;
    l4_9.gf[1] = ((l2_9.x[0] + .002 + l2_9.x[1] * 60.) * 14. - (l2_9.x[0] * 
	    .15 + l2_9.x[1] * 14. - .06) * 60.) * 2e4 / pow_dd(&d__1, &c_b305)
	    ;
    l4_9.gf[0] = -l4_9.gf[0];
    l4_9.gf[1] = -l4_9.gf[1];
    return 0;
labelL4:
    if (l9_14.index1[0]) {
	l3_13.g[0] = l2_9.x[0] - .75 / l2_9.x[2] / l2_9.x[3];
    }
    if (l9_14.index1[1]) {
	l3_13.g[1] = l2_9.x[0] - l2_9.x[8] / l2_9.x[4] / l2_9.x[3];
    }
    if (l9_14.index1[2]) {
	l3_13.g[2] = l2_9.x[0] - l2_9.x[9] / l2_9.x[5] / l2_9.x[3] - 10. / 
		l2_9.x[3];
    }
    if (l9_14.index1[3]) {
	l3_13.g[3] = l2_9.x[0] - .19 / l2_9.x[6] / l2_9.x[3] - 10. / l2_9.x[3]
		;
    }
    if (l9_14.index1[4]) {
	l3_13.g[4] = l2_9.x[0] - .125 / l2_9.x[7] / l2_9.x[3];
    }
    if (l9_14.index1[5]) {
	l3_13.g[5] = l2_9.x[1] * 1e4 - l2_9.x[8] * .00131 * pow_dd(&l2_9.x[4],
		 &c_b2746) * pow_dd(&l2_9.x[3], &c_b940);
    }
    if (l9_14.index1[6]) {
	l3_13.g[6] = l2_9.x[1] * 1e4 - l2_9.x[9] * .001038 * pow_dd(&l2_9.x[5]
		, &c_b2748) * pow_dd(&l2_9.x[3], &c_b1523);
    }
    if (l9_14.index1[7]) {
	l3_13.g[7] = l2_9.x[1] * 1e4 - pow_dd(&l2_9.x[6], &c_b2746) * 2.23e-4 
		* pow_dd(&l2_9.x[3], &c_b940);
    }
    if (l9_14.index1[8]) {
	l3_13.g[8] = l2_9.x[1] * 1e4 - pow_dd(&l2_9.x[7], &c_b2752) * 7.6e-5 *
		 pow_dd(&l2_9.x[3], &c_b2753);
    }
    if (l9_14.index1[9]) {
/* Computing 2nd power */
	d__1 = l2_9.x[3];
	l3_13.g[9] = l2_9.x[1] * 1e4 - pow_dd(&l2_9.x[2], &c_b2754) * 6.98e-4 
		* (d__1 * d__1);
    }
    if (l9_14.index1[10]) {
	l3_13.g[10] = l2_9.x[1] * 1e4 - pow_dd(&l2_9.x[2], &c_b2748) * 5e-5 * 
		pow_dd(&l2_9.x[3], &c_b1523);
    }
    if (l9_14.index1[11]) {
	l3_13.g[11] = l2_9.x[1] * 1e4 - pow_dd(&l2_9.x[2], &c_b2757) * 
		6.54e-6 * pow_dd(&l2_9.x[3], &c_b2758);
    }
    if (l9_14.index1[12]) {
	l3_13.g[12] = l2_9.x[1] * 1e4 - pow_dd(&l2_9.x[2], &c_b2746) * 
		2.57e-4 * pow_dd(&l2_9.x[3], &c_b940);
    }
    if (l9_14.index1[13]) {
	l3_13.g[13] = 30. - l2_9.x[4] * 2.003 * l2_9.x[3] - l2_9.x[5] * 1.885 
		* l2_9.x[3] - l2_9.x[7] * .184 * l2_9.x[3] - pow_dd(&l2_9.x[2]
		, &c_b2761) * 2. * l2_9.x[3];
    }
    if (l9_14.index1[14]) {
	l3_13.g[14] = l2_9.x[8] + l2_9.x[9] - .255;
    }
    return 0;
labelL5:
    if (! l10_14.index2[0]) {
	goto L50;
    }
    l5_32.gg[30] = .75 / pow_dd(&l2_9.x[2], &c_b305) / l2_9.x[3];
    l5_32.gg[45] = .75 / l2_9.x[2] / pow_dd(&l2_9.x[3], &c_b305);
L50:
    if (! l10_14.index2[1]) {
	goto L51;
    }
    l5_32.gg[46] = l2_9.x[8] / l2_9.x[4] / pow_dd(&l2_9.x[3], &c_b305);
    l5_32.gg[61] = l2_9.x[8] / pow_dd(&l2_9.x[4], &c_b305) / l2_9.x[3];
    l5_32.gg[121] = -1. / l2_9.x[4] / l2_9.x[3];
L51:
    if (! l10_14.index2[2]) {
	goto L52;
    }
    l5_32.gg[47] = l2_9.x[9] / l2_9.x[5] / pow_dd(&l2_9.x[3], &c_b305) + 10. /
	     pow_dd(&l2_9.x[3], &c_b305);
    l5_32.gg[77] = l2_9.x[9] / pow_dd(&l2_9.x[5], &c_b305) / l2_9.x[3];
    l5_32.gg[137] = -1. / l2_9.x[5] / l2_9.x[3];
L52:
    if (! l10_14.index2[3]) {
	goto L53;
    }
    l5_32.gg[48] = .19 / l2_9.x[6] / pow_dd(&l2_9.x[3], &c_b305) + 10. / 
	    pow_dd(&l2_9.x[3], &c_b305);
    l5_32.gg[93] = .19 / pow_dd(&l2_9.x[6], &c_b305) / l2_9.x[3];
L53:
    if (! l10_14.index2[4]) {
	goto L54;
    }
    l5_32.gg[49] = .125 / l2_9.x[7] / pow_dd(&l2_9.x[3], &c_b305);
    l5_32.gg[109] = .125 / pow_dd(&l2_9.x[7], &c_b305) / l2_9.x[3];
L54:
    if (! l10_14.index2[5]) {
	goto L55;
    }
    l5_32.gg[50] = l2_9.x[8] * -.0019649999999999997 * pow_dd(&l2_9.x[4], &
	    c_b2746) * pow_dd(&l2_9.x[3], &c_b590);
    l5_32.gg[65] = l2_9.x[8] * -8.7246000000000003e-4 / pow_dd(&l2_9.x[4], &
	    c_b2782) * pow_dd(&l2_9.x[3], &c_b940);
    l5_32.gg[125] = pow_dd(&l2_9.x[4], &c_b2746) * -.00131 * pow_dd(&l2_9.x[3]
	    , &c_b940);
L55:
    if (! l10_14.index2[6]) {
	goto L56;
    }
    l5_32.gg[51] = l2_9.x[9] * -.0031140000000000004 * pow_dd(&l2_9.x[5], &
	    c_b2748) * pow_dd(&l2_9.x[3], &c_b305);
    l5_32.gg[81] = l2_9.x[9] * -.0016608000000000003 * pow_dd(&l2_9.x[5], &
	    c_b2789) * pow_dd(&l2_9.x[3], &c_b1523);
    l5_32.gg[141] = pow_dd(&l2_9.x[5], &c_b2748) * -.001038 * pow_dd(&l2_9.x[
	    3], &c_b1523);
L56:
    if (! l10_14.index2[7]) {
	goto L57;
    }
    l5_32.gg[52] = pow_dd(&l2_9.x[6], &c_b2746) * -3.345e-4 * pow_dd(&l2_9.x[
	    3], &c_b590);
    l5_32.gg[97] = -1.4851800000000002e-4 / pow_dd(&l2_9.x[6], &c_b2782) * 
	    pow_dd(&l2_9.x[3], &c_b940);
L57:
    if (! l10_14.index2[8]) {
	goto L58;
    }
    l5_32.gg[53] = pow_dd(&l2_9.x[7], &c_b2752) * -4.3016000000000002e-4 * 
	    pow_dd(&l2_9.x[3], &c_b2800);
    l5_32.gg[113] = pow_dd(&l2_9.x[7], &c_b2801) * -2.698e-4 * pow_dd(&l2_9.x[
	    3], &c_b2753);
L58:
    if (! l10_14.index2[9]) {
	goto L59;
    }
    l5_32.gg[39] = pow_dd(&l2_9.x[2], &c_b1157) * -8.3760000000000008e-4 * 
	    pow_dd(&l2_9.x[3], &c_b305);
    l5_32.gg[54] = pow_dd(&l2_9.x[2], &c_b2754) * -.0013960000000000001 * 
	    l2_9.x[3];
L59:
    if (! l10_14.index2[10]) {
	goto L60;
    }
    l5_32.gg[40] = pow_dd(&l2_9.x[2], &c_b2789) * -8.0000000000000007e-5 * 
	    pow_dd(&l2_9.x[3], &c_b1523);
    l5_32.gg[55] = pow_dd(&l2_9.x[2], &c_b2748) * -1.5000000000000001e-4 * 
	    pow_dd(&l2_9.x[3], &c_b305);
L60:
    if (! l10_14.index2[11]) {
	goto L61;
    }
    l5_32.gg[41] = pow_dd(&l2_9.x[2], &c_b2813) * -1.5826799999999999e-5 * 
	    pow_dd(&l2_9.x[3], &c_b2758);
    l5_32.gg[56] = pow_dd(&l2_9.x[2], &c_b2757) * -2.7271799999999999e-5 * 
	    pow_dd(&l2_9.x[3], &c_b2816);
L61:
    if (! l10_14.index2[12]) {
	goto L62;
    }
    l5_32.gg[42] = -1.7116200000000001e-4 / pow_dd(&l2_9.x[2], &c_b2782) * 
	    pow_dd(&l2_9.x[3], &c_b940);
    l5_32.gg[57] = pow_dd(&l2_9.x[2], &c_b2746) * -3.8549999999999999e-4 * 
	    pow_dd(&l2_9.x[3], &c_b590);
L62:
    if (! l10_14.index2[13]) {
	goto L63;
    }
    l5_32.gg[43] = -1.6060000000000001 / pow_dd(&l2_9.x[2], &c_b2823) * 
	    l2_9.x[3];
    l5_32.gg[58] = l2_9.x[4] * -2.003 - l2_9.x[5] * 1.885 - l2_9.x[7] * .184 
	    - pow_dd(&l2_9.x[2], &c_b2761) * 2.;
    l5_32.gg[73] = l2_9.x[3] * -2.003;
    l5_32.gg[88] = l2_9.x[3] * -1.885;
    l5_32.gg[118] = l2_9.x[3] * -.184;
L63:
    return 0;
} /* tp376_ */


/* Subroutine */ int tp377_(int *mode)
{
    /* Initialized data */

    static Real a[10] = { -6.089,-17.164,-34.054,-5.914,-24.721,-14.986,
	    -24.1,-10.708,-26.662,-22.179 };

    /* Local variables */
    static int i__;
    static Real sum;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 10;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 3;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 10; ++i__) {
	l2_9.x[i__ - 1] = .1;
	l12_9.lxu[i__ - 1] = true;
	l11_9.lxl[i__ - 1] = true;
	l13_9.xl[i__ - 1] = 1e-5;
/* labelL6: */
	l14_9.xu[i__ - 1] = 10.;
    }
    l5_14.gg[0] = 1.;
    l5_14.gg[3] = -2.;
    l5_14.gg[6] = 2.;
    l5_14.gg[9] = 0.;
    l5_14.gg[12] = 0.;
    l5_14.gg[15] = 1.;
    for (i__ = 7; i__ <= 9; ++i__) {
/* L7: */
	l5_14.gg[i__ * 3 - 3] = 0.;
    }
    l5_14.gg[27] = 1.;
    for (i__ = 1; i__ <= 3; ++i__) {
/* L8: */
	l5_14.gg[i__ * 3 - 2] = 0.;
    }
    l5_14.gg[10] = 1.;
    l5_14.gg[13] = -2.;
    l5_14.gg[16] = 1.;
    l5_14.gg[19] = 1.;
    for (i__ = 8; i__ <= 10; ++i__) {
/* labelL9: */
	l5_14.gg[i__ * 3 - 2] = 0.;
    }
    l5_14.gg[2] = 0.;
    l5_14.gg[5] = 0.;
    l5_14.gg[8] = 1.;
    for (i__ = 4; i__ <= 6; ++i__) {
/* labelL10: */
	l5_14.gg[i__ * 3 - 1] = 0.;
    }
    l5_14.gg[20] = 1.;
    l5_14.gg[23] = 1.;
    l5_14.gg[26] = 2.;
    l5_14.gg[29] = 1.;
    l20_10.lex = false;
    l20_10.nex = 1;
    l20_10.fex = (float)-795.001;
    l20_10.xex[0] = (float)10.;
    l20_10.xex[1] = (float)10.;
    l20_10.xex[2] = (float)1.;
    l20_10.xex[3] = (float)10.;
    l20_10.xex[4] = (float)9.5;
    l20_10.xex[5] = (float)10.;
    l20_10.xex[6] = 1e-4;
    l20_10.xex[7] = 1e-4;
    l20_10.xex[8] = 1e-4;
    l20_10.xex[9] = 1e-4;
    return 0;
labelL2:
    l6_1.fx = 0.;
    sum = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* labelL11: */
	sum += l2_9.x[i__ - 1];
    }
    for (i__ = 1; i__ <= 10; ++i__) {
/* labelL12: */
	l6_1.fx += l2_9.x[i__ - 1] * (a[i__ - 1] + std::log(l2_9.x[i__ - 1] / sum))
		;
    }
labelL3:
    return 0;
    sum = 0.;
    for (i__ = 1; i__ <= 10; ++i__) {
/* L45: */
	sum += l2_9.x[i__ - 1];
    }
    for (i__ = 1; i__ <= 10; ++i__) {
/* L46: */
	l4_9.gf[i__ - 1] = a[i__ - 1] + std::log(l2_9.x[i__ - 1] / sum);
    }
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = l2_9.x[0] - l2_9.x[1] * 2. + l2_9.x[2] * 2. + l2_9.x[5] + 
		l2_9.x[9] - 2.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = l2_9.x[3] - l2_9.x[4] * 2. + l2_9.x[5] + l2_9.x[6] - 1.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = l2_9.x[2] + l2_9.x[6] + l2_9.x[7] + l2_9.x[8] * 2. + 
		l2_9.x[9] - 1.;
    }
labelL5:
    return 0;
} /* tp377_ */


/* Subroutine */ int tp378_(int *mode)
{
    /* Initialized data */

    static Real a[10] = { -6.089,-17.164,-34.054,-5.914,-24.721,-14.986,
	    -24.1,-10.708,-26.662,-22.179 };

    /* Local variables */
    static int i__, j;
    static Real con;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 10;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 3;
    for (i__ = 1; i__ <= 10; ++i__) {
	l2_9.x[i__ - 1] = -2.3;
	l14_9.xu[i__ - 1] = (float)-.1;
	l12_9.lxu[i__ - 1] = true;
	l13_9.xl[i__ - 1] = (float)-16.;
/* labelL6: */
	l11_9.lxl[i__ - 1] = true;
    }
    l20_10.lex = false;
    l20_10.nex = 1;
    l20_10.fex = -47.76;
    l20_10.xex[0] = -3.2024;
    l20_10.xex[1] = -1.9123;
    l20_10.xex[2] = -.2444;
    l20_10.xex[3] = -15.67;
    l20_10.xex[4] = -.7217;
    l20_10.xex[5] = -7.2736;
    l20_10.xex[6] = -3.5965;
    l20_10.xex[7] = -4.0206;
    l20_10.xex[8] = -3.2885;
    l20_10.xex[9] = -2.3344;
    l5_14.gg[9] = 0.;
    l5_14.gg[12] = 0.;
    l5_14.gg[18] = 0.;
    l5_14.gg[21] = 0.;
    l5_14.gg[24] = 0.;
    l5_14.gg[1] = 0.;
    l5_14.gg[4] = 0.;
    l5_14.gg[7] = 0.;
    l5_14.gg[22] = 0.;
    l5_14.gg[25] = 0.;
    l5_14.gg[28] = 0.;
    l5_14.gg[2] = 0.;
    l5_14.gg[5] = 0.;
    l5_14.gg[11] = 0.;
    l5_14.gg[14] = 0.;
    l5_14.gg[17] = 0.;
    return 0;
labelL2:
    l6_1.fx = 0.;
    con = 0.;
    for (j = 1; j <= 10; ++j) {
/* L7: */
	con += std::exp(l2_9.x[j - 1]);
    }
    con = std::log(con);
    for (i__ = 1; i__ <= 10; ++i__) {
/* L8: */
	l6_1.fx += std::exp(l2_9.x[i__ - 1]) * (a[i__ - 1] + l2_9.x[i__ - 1] - con)
		;
    }
labelL3:
    return 0;
    con = 0.;
    for (j = 1; j <= 10; ++j) {
/* L45: */
	con += std::exp(l2_9.x[j - 1]);
    }
    con = std::log(con);
    for (i__ = 1; i__ <= 10; ++i__) {
	l4_9.gf[i__ - 1] = std::exp(l2_9.x[i__ - 1]) * (a[i__ - 1] + l2_9.x[i__ - 
		1] - con);
/* L46: */
    }
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = std::exp(l2_9.x[0]) + std::exp(l2_9.x[1]) * 2. + std::exp(l2_9.x[2]) * 
		2. + std::exp(l2_9.x[5]) + std::exp(l2_9.x[9]) - 2.;
    }
    if (l9_4.index1[1]) {
	l3_3.g[1] = std::exp(l2_9.x[3]) + std::exp(l2_9.x[4]) * 2. + std::exp(l2_9.x[5]) + 
		exp(l2_9.x[6]) - 1.;
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = std::exp(l2_9.x[2]) + std::exp(l2_9.x[6]) + std::exp(l2_9.x[7]) + std::exp(
		l2_9.x[8]) * 2. + std::exp(l2_9.x[9]) - 1.;
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto L47;
    }
    l5_14.gg[0] = std::exp(l2_9.x[0]);
    l5_14.gg[3] = std::exp(l2_9.x[1]) * 2.;
    l5_14.gg[6] = std::exp(l2_9.x[2]) * 2.;
    l5_14.gg[15] = std::exp(l2_9.x[5]);
    l5_14.gg[27] = std::exp(l2_9.x[9]);
L47:
    if (! l10_4.index2[1]) {
	goto L48;
    }
    l5_14.gg[10] = std::exp(l2_9.x[3]);
    l5_14.gg[13] = std::exp(l2_9.x[4]) * 2.;
    l5_14.gg[16] = std::exp(l2_9.x[5]);
    l5_14.gg[19] = std::exp(l2_9.x[6]);
L48:
    if (! l10_4.index2[2]) {
	goto L49;
    }
    l5_14.gg[8] = std::exp(l2_9.x[2]);
    l5_14.gg[20] = std::exp(l2_9.x[6]);
    l5_14.gg[23] = std::exp(l2_9.x[7]);
    l5_14.gg[26] = std::exp(l2_9.x[8]) * 2.;
    l5_14.gg[29] = std::exp(l2_9.x[9]);
L49:
    return 0;
} /* tp378_ */

/* Subroutine */ int tp379_(int *mode)
{
    /* Initialized data */

    static Real y[65] = { 1.366,1.191,1.112,1.013,.991,.885,.831,.847,
	    .786,.725,.746,.679,.608,.655,.616,.606,.602,.626,.651,.724,.649,
	    .649,.694,.644,.624,.661,.612,.558,.533,.495,.5,.423,.395,.375,
	    .372,.391,.396,.405,.428,.429,.523,.562,.607,.653,.672,.708,.633,
	    .668,.645,.632,.591,.559,.597,.625,.739,.71,.729,.72,.636,.581,
	    .428,.292,.162,.098,.054 };

    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static int i__, j;
    static Real t;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL2;
	case 4:  goto labelL4;
	case 5:  goto labelL4;
    }
labelL1:
    l1_1.n = 11;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l15_1.lsum = 65;
    l2_17.x[0] = 1.3;
    l2_17.x[1] = .65;
    l2_17.x[2] = .65;
    l2_17.x[3] = .7;
    l2_17.x[4] = .6;
    l2_17.x[5] = 3.;
    l2_17.x[6] = 5.;
    l2_17.x[7] = 7.;
    l2_17.x[8] = 2.;
    l2_17.x[9] = 4.5;
    l2_17.x[10] = 5.5;
    for (i__ = 1; i__ <= 11; ++i__) {
	l12_17.lxu[i__ - 1] = false;
	l13_21.xl[i__ - 1] = (float)0.;
/* labelL6: */
	l11_17.lxl[i__ - 1] = true;
    }
    l20_18.lex = false;
    l20_18.nex = 1;
    l20_18.fex = .0401377;
    l20_18.xex[0] = 1.30997;
    l20_18.xex[1] = .431554;
    l20_18.xex[2] = .63366;
    l20_18.xex[3] = .59943;
    l20_18.xex[4] = .754183;
    l20_18.xex[5] = .904286;
    l20_18.xex[6] = 1.36581;
    l20_18.xex[7] = 4.82369;
    l20_18.xex[8] = 2.39868;
    l20_18.xex[9] = 4.56887;
    l20_18.xex[10] = 5.67534;
    return 0;
labelL2:
    for (i__ = 1; i__ <= 65; ++i__) {
	t = (Real) (i__ - 1) * .1;
/* L7: */
/* Computing 2nd power */
	d__1 = t - l2_17.x[8];
/* Computing 2nd power */
	d__2 = t - l2_17.x[9];
/* Computing 2nd power */
	d__3 = t - l2_17.x[10];
	l16_19.f[i__ - 1] = y[i__ - 1] - (l2_17.x[0] * std::exp(-l2_17.x[4] * t) + 
		l2_17.x[1] * std::exp(-l2_17.x[5] * (d__1 * d__1)) + l2_17.x[2] * 
		exp(-l2_17.x[6] * (d__2 * d__2)) + l2_17.x[3] * std::exp(-l2_17.x[
		7] * (d__3 * d__3)));
    }
    if (*mode == 3) {
	goto labelL3;
    }
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 65; ++i__) {
/* L70: */
/* Computing 2nd power */
	d__1 = l16_19.f[i__ - 1];
	l6_1.fx += d__1 * d__1;
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 65; ++i__) {
	t = (Real) (i__ - 1) * .1;
	l17_24.df[i__ - 1] = -exp(-l2_17.x[4] * t);
/* Computing 2nd power */
	d__1 = t - l2_17.x[8];
	l17_24.df[i__ + 64] = -exp(-l2_17.x[5] * (d__1 * d__1));
/* Computing 2nd power */
	d__1 = t - l2_17.x[9];
	l17_24.df[i__ + 129] = -exp(-l2_17.x[6] * (d__1 * d__1));
/* Computing 2nd power */
	d__1 = t - l2_17.x[10];
	l17_24.df[i__ + 194] = -exp(-l2_17.x[7] * (d__1 * d__1));
	l17_24.df[i__ + 259] = l2_17.x[0] * t * std::exp(-l2_17.x[4] * t);
/* Computing 2nd power */
	d__1 = t - l2_17.x[8];
/* Computing 2nd power */
	d__2 = t - l2_17.x[8];
	l17_24.df[i__ + 324] = l2_17.x[1] * (d__1 * d__1) * std::exp(-l2_17.x[5] * 
		(d__2 * d__2));
/* Computing 2nd power */
	d__1 = t - l2_17.x[9];
/* Computing 2nd power */
	d__2 = t - l2_17.x[9];
	l17_24.df[i__ + 389] = l2_17.x[2] * (d__1 * d__1) * std::exp(-l2_17.x[6] * 
		(d__2 * d__2));
/* Computing 2nd power */
	d__1 = t - l2_17.x[10];
/* Computing 2nd power */
	d__2 = t - l2_17.x[10];
	l17_24.df[i__ + 454] = l2_17.x[3] * (d__1 * d__1) * std::exp(-l2_17.x[7] * 
		(d__2 * d__2));
/* Computing 2nd power */
	d__1 = t - l2_17.x[8];
	l17_24.df[i__ + 519] = -l2_17.x[1] * l2_17.x[5] * 2. * (t - l2_17.x[8]
		) * std::exp(-l2_17.x[5] * (d__1 * d__1));
/* Computing 2nd power */
	d__1 = t - l2_17.x[9];
	l17_24.df[i__ + 584] = -l2_17.x[2] * l2_17.x[6] * 2. * (t - l2_17.x[9]
		) * std::exp(-l2_17.x[6] * (d__1 * d__1));
/* L8: */
/* Computing 2nd power */
	d__1 = t - l2_17.x[10];
	l17_24.df[i__ + 649] = -l2_17.x[3] * l2_17.x[7] * 2. * (t - l2_17.x[
		10]) * std::exp(-l2_17.x[7] * (d__1 * d__1));
    }
    for (i__ = 1; i__ <= 11; ++i__) {
/* L19: */
	l4_17.gf[i__ - 1] = 0.;
    }
    for (j = 1; j <= 11; ++j) {
	for (i__ = 1; i__ <= 65; ++i__) {
/* labelL20: */
	    l4_17.gf[j - 1] += l16_19.f[i__ - 1] * 2. * l17_24.df[i__ + j * 
		    65 - 66];
	}
    }
labelL4:
    return 0;
} /* tp379_ */


/* Subroutine */ int tp380_(int *mode)
{
    /* Initialized data */

    static Real a[11] = { -.00133172,-.002270927,-.00248546,-4.67,
	    -4.671973,-.00814,-.008092,-.005,-9.09e-4,-8.8e-4,-.00119 };
    static Real c__[30] = { .05367373,.021863746,.097733533,.0066940803,
	    1e-6,1e-5,1e-6,1e-10,1e-8,.01,1e-4,.10898645,1.6108052e-4,1e-23,
	    1.9304541e-6,.001,1e-6,1e-5,1e-6,1e-9,1e-9,.001,.001,.10898645,
	    1.6108052e-5,1e-23,1.9304541e-8,1e-5,1.1184059e-4,1e-4 };

    /* System generated locals */
    Real d__1, d__2, d__3;

    /* Local variables */
    static Real temp;
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 12;
    l1_1.nili = 0;
    l1_1.ninl = 3;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 12; ++i__) {
	l2_18.x[i__ - 1] = 4.;
	l12_18.lxu[i__ - 1] = true;
	l11_18.lxl[i__ - 1] = true;
	l13_22.xl[i__ - 1] = .1;
/* labelL6: */
	l14_22.xu[i__ - 1] = 100.;
    }
    l5_10.gg[0] = -c__[0];
    l5_10.gg[3] = -c__[1];
    l5_10.gg[6] = -c__[2];
    for (i__ = 6; i__ <= 12; ++i__) {
/* L7: */
	l5_10.gg[i__ * 3 - 3] = 0.;
    }
    l5_10.gg[1] = -c__[4];
    l5_10.gg[7] = -c__[6];
    l5_10.gg[22] = 0.;
    l5_10.gg[25] = 0.;
    l5_10.gg[31] = 0.;
    l5_10.gg[8] = -c__[18];
    l5_10.gg[17] = -c__[21];
    l5_10.gg[20] = 0.;
    l5_10.gg[23] = -c__[22];
    l5_10.gg[29] = 0.;
    l5_10.gg[32] = -c__[29];
    l5_10.gg[35] = 0.;
    l20_5.lex = false;
    l20_5.nex = 1;
/*      FEX=0.316859D+1 */
    l20_5.fex = 316822.15000000002;
    l20_5.xex[0] = 2.6631947068;
    l20_5.xex[1] = 4.517277762;
    l20_5.xex[2] = 7.133802907;
    l20_5.xex[3] = 2.237268448;
    l20_5.xex[4] = 4.07840382657;
    l20_5.xex[5] = 1.31827569;
    l20_5.xex[6] = 4.125187034;
    l20_5.xex[7] = 2.856195978;
    l20_5.xex[8] = 1.6765929748;
    l20_5.xex[9] = 2.1789111052;
    l20_5.xex[10] = 5.12343515;
    l20_5.xex[11] = 6.659338016;
    return 0;
labelL2:
    l6_1.fx = 1e5;
    for (i__ = 1; i__ <= 11; ++i__) {
	temp = l2_18.x[i__ - 1];
	if (l2_18.x[i__ - 1] < 1e-15) {
	    temp = 1e-15;
	}
/* labelL20: */
	l6_1.fx *= pow_dd(&temp, &a[i__ - 1]);
    }
    l6_1.fx *= 1e5;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 11; ++i__) {
	temp = l2_18.x[i__ - 1];
	if (l2_18.x[i__ - 1] < 1e-15) {
	    temp = 1e-15;
	}
	l4_18.gf[i__ - 1] = l6_1.fx * (a[i__ - 1] / temp);
/* labelL10: */
	l4_18.gf[i__ - 1] *= 1e5;
    }
    l4_18.gf[11] = 0.;
    return 0;
labelL4:
    if (l9_4.index1[0]) {
	l3_3.g[0] = 1. - c__[0] * l2_18.x[0] - c__[1] * l2_18.x[1] - c__[2] * 
		l2_18.x[2] - c__[3] * l2_18.x[3] * l2_18.x[4];
    }
    if (l9_4.index1[1]) {
/* Computing 2nd power */
	d__1 = l2_18.x[11];
	l3_3.g[1] = 1. - c__[4] * l2_18.x[0] - c__[5] * l2_18.x[1] - c__[6] * 
		l2_18.x[2] - c__[7] * l2_18.x[3] * l2_18.x[11] - c__[8] * 
		l2_18.x[4] / l2_18.x[11] - c__[9] * l2_18.x[5] / l2_18.x[11] 
		- c__[10] * l2_18.x[6] * l2_18.x[11] - c__[11] * l2_18.x[3] * 
		l2_18.x[4] - c__[12] * l2_18.x[1] * l2_18.x[4] / l2_18.x[11] 
		- c__[13] * l2_18.x[1] * l2_18.x[3] * l2_18.x[4] - c__[14] * 
		l2_18.x[1] / l2_18.x[3] * l2_18.x[4] / (d__1 * d__1) - c__[15]
		 * l2_18.x[9] / l2_18.x[11];
    }
    if (l9_4.index1[2]) {
	l3_3.g[2] = 1. - c__[16] * l2_18.x[0] - c__[17] * l2_18.x[1] - c__[18]
		 * l2_18.x[2] - c__[19] * l2_18.x[3] - c__[20] * l2_18.x[4] - 
		c__[21] * l2_18.x[5] - c__[22] * l2_18.x[7] - c__[23] * 
		l2_18.x[3] * l2_18.x[4] - c__[24] * l2_18.x[1] * l2_18.x[4] - 
		c__[25] * l2_18.x[1] * l2_18.x[3] * l2_18.x[4] - c__[26] * 
		l2_18.x[1] * l2_18.x[4] / l2_18.x[3] - c__[27] * l2_18.x[8] - 
		c__[28] * l2_18.x[0] * l2_18.x[8] - c__[29] * l2_18.x[10];
    }
    return 0;
labelL5:
    if (! l10_4.index2[0]) {
	goto L50;
    }
    l5_10.gg[9] = -c__[3] * l2_18.x[4];
    l5_10.gg[12] = -c__[3] * l2_18.x[3];
L50:
    if (! l10_4.index2[1]) {
	goto L51;
    }
/* Computing 2nd power */
    d__1 = l2_18.x[11];
    l5_10.gg[4] = -c__[5] - c__[12] * l2_18.x[4] / l2_18.x[11] - c__[13] * 
	    l2_18.x[3] * l2_18.x[4] - c__[14] / l2_18.x[3] * l2_18.x[4] / (
	    d__1 * d__1);
/* Computing 2nd power */
    d__1 = l2_18.x[3];
/* Computing 2nd power */
    d__2 = l2_18.x[11];
    l5_10.gg[10] = -c__[7] * l2_18.x[11] - c__[11] * l2_18.x[4] - c__[13] * 
	    l2_18.x[1] * l2_18.x[4] + c__[14] * l2_18.x[1] / (d__1 * d__1) * 
	    l2_18.x[4] / (d__2 * d__2);
/* Computing 2nd power */
    d__1 = l2_18.x[11];
    l5_10.gg[13] = -c__[8] / l2_18.x[11] - c__[11] * l2_18.x[3] - c__[12] * 
	    l2_18.x[1] / l2_18.x[11] - c__[13] * l2_18.x[1] * l2_18.x[3] - 
	    c__[14] * l2_18.x[1] / l2_18.x[3] / (d__1 * d__1);
    l5_10.gg[16] = -c__[9] / l2_18.x[11];
    l5_10.gg[19] = -c__[10] * l2_18.x[11];
    l5_10.gg[28] = -c__[15] / l2_18.x[11];
/* Computing 2nd power */
    d__1 = l2_18.x[11];
/* Computing 2nd power */
    d__2 = l2_18.x[11];
/* Computing 2nd power */
    d__3 = l2_18.x[11];
    l5_10.gg[34] = -c__[7] * l2_18.x[3] + c__[8] * l2_18.x[4] / (d__1 * d__1) 
	    + c__[9] * l2_18.x[5] / pow_dd(&l2_18.x[11], &c_b305) - c__[10] * 
	    l2_18.x[6] + c__[12] * l2_18.x[1] * l2_18.x[4] / (d__2 * d__2) + 
	    c__[14] * 2 * l2_18.x[1] / l2_18.x[3] * l2_18.x[4] / pow_dd(&
	    l2_18.x[11], &c_b1523) + c__[15] * l2_18.x[9] / (d__3 * d__3);
L51:
    if (! l10_4.index2[2]) {
	goto L52;
    }
    l5_10.gg[2] = -c__[16] - c__[28] * l2_18.x[8];
    l5_10.gg[5] = -c__[17] - c__[24] * l2_18.x[4] - c__[25] * l2_18.x[3] * 
	    l2_18.x[4] - c__[26] * l2_18.x[4] / l2_18.x[3];
/* Computing 2nd power */
    d__1 = l2_18.x[3];
    l5_10.gg[11] = -c__[19] - c__[23] * l2_18.x[4] - c__[25] * l2_18.x[1] * 
	    l2_18.x[4] + c__[26] * l2_18.x[1] * l2_18.x[4] / (d__1 * d__1);
    l5_10.gg[14] = -c__[20] - c__[23] * l2_18.x[3] - c__[24] * l2_18.x[1] - 
	    c__[25] * l2_18.x[1] * l2_18.x[3] - c__[26] * l2_18.x[1] / 
	    l2_18.x[3];
    l5_10.gg[26] = -c__[27] - c__[28] * l2_18.x[0];
L52:
    return 0;
} /* tp380_ */


/* Subroutine */ int tp381_(int *mode)
{
    /* Initialized data */

    static Real r__[13] = { .8,1.1,.85,3.45,2.,2.1,3.,.8,.45,.72,1.8,3.,
	    .6 };
    static Real s[13] = { 11.6,13.7,9.5,48.5,31.9,51.1,65.5,0.,0.,0.,
	    21.8,46.9,0. };
    static Real u[13] = { .05,.07,0.,.33,0.,1.27,1.27,23.35,35.84,.81,
	    1.79,7.34,0. };
    static Real v[13] = { .35,.37,.1,.62,0.,1.03,1.69,18.21,.01,.08,.31,
	    1.59,22.45 };

    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 13;
    l1_1.nili = 3;
    l1_1.ninl = 0;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 13; ++i__) {
	l2_10.x[i__ - 1] = .1;
	l12_10.lxu[i__ - 1] = false;
	l11_10.lxl[i__ - 1] = true;
/* labelL6: */
	l13_10.xl[i__ - 1] = 0.;
    }
    for (i__ = 1; i__ <= 13; ++i__) {
	l5_38.gg[(i__ << 2) - 4] = s[i__ - 1];
	l5_38.gg[(i__ << 2) - 3] = u[i__ - 1];
/* L7: */
	l5_38.gg[(i__ << 2) - 1] = 1.;
    }
    for (i__ = 1; i__ <= 13; ++i__) {
/* L8: */
	l5_38.gg[(i__ << 2) - 2] = v[i__ - 1];
    }
    l20_11.lex = false;
    l20_11.nex = 1;
    l20_11.fex = 1.0149;
    l20_11.xex[0] = .785586;
    l20_11.xex[1] = 0.;
    l20_11.xex[2] = 0.;
    l20_11.xex[3] = 0.;
    l20_11.xex[4] = 0.;
    l20_11.xex[5] = .173918;
    l20_11.xex[6] = 0.;
    l20_11.xex[7] = 0.;
    l20_11.xex[8] = .020643;
    l20_11.xex[9] = 0.;
    l20_11.xex[10] = 0.;
    l20_11.xex[11] = 0.;
    l20_11.xex[12] = .019853;
    for (i__ = 1; i__ <= 13; ++i__) {
/* labelL11: */
	l4_10.gf[i__ - 1] = r__[i__ - 1];
    }
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* labelL10: */
	l6_1.fx += r__[i__ - 1] * l2_10.x[i__ - 1];
    }
labelL3:
    return 0;
labelL4:
    if (! l9_7.index1[0]) {
	goto labelL12;
    }
    l3_6.g[0] = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* labelL13: */
	l3_6.g[0] += s[i__ - 1] * l2_10.x[i__ - 1];
    }
    l3_6.g[0] += -18.;
labelL12:
    if (! l9_7.index1[1]) {
	goto labelL14;
    }
    l3_6.g[1] = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* L15: */
	l3_6.g[1] += u[i__ - 1] * l2_10.x[i__ - 1];
    }
    l3_6.g[1] += -1.;
labelL14:
    if (! l9_7.index1[2]) {
	goto L16;
    }
    l3_6.g[2] = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* L17: */
	l3_6.g[2] += v[i__ - 1] * l2_10.x[i__ - 1];
    }
    l3_6.g[2] += -.9;
L16:
    if (! l9_7.index1[3]) {
	goto labelL5;
    }
    l3_6.g[3] = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* L19: */
	l3_6.g[3] += l2_10.x[i__ - 1];
    }
    l3_6.g[3] += -1.;
labelL5:
    return 0;
} /* tp381_ */

/* Subroutine */ int tp382_(int *mode)
{
    /* Initialized data */

    static Real r__[13] = { .8,1.1,.85,3.45,2.,2.1,3.,.8,.45,.72,1.8,3.,
	    .6 };
    static Real s[13] = { 11.6,13.7,9.5,48.5,31.9,51.1,65.5,0.,0.,0.,
	    21.8,46.9,0. };
    static Real z1[13] = { .4844,.3003,.1444,.0588,4.9863,.0653,21.0222,
	    0.,0.,0.,.297,9.2933,0. };
    static Real u[13] = { .05,.07,0.,.33,0.,1.27,1.27,23.35,35.84,.81,
	    1.79,7.34,0. };
    static Real z2[13] = { 1e-4,0.,0.,0.,0.,.004,.1404,1.3631,.5138,
	    .0289,.0097,.3893,0. };
    static Real v[13] = { .35,.37,.1,.62,0.,1.03,1.69,18.21,.01,.08,.31,
	    1.59,22.45 };
    static Real z3[13] = { .001,9e-4,1e-4,5e-4,0.,.0021,.0825,.2073,0.,
	    4e-4,5e-4,.0107,1.0206 };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real help;
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }

labelL1:
    l1_1.n = 13;
    l1_1.nili = 0;
    l1_1.ninl = 3;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 13; ++i__) {
	l2_10.x[i__ - 1] = .1;
	l12_10.lxu[i__ - 1] = false;
	l11_10.lxl[i__ - 1] = true;
	l13_10.xl[i__ - 1] = 0.;
/* labelL6: */
	l5_38.gg[(i__ << 2) - 1] = 1.;
    }
    l20_11.lex = false;
    l20_11.nex = 1;
    l20_11.fex = 1.03831;
    l20_11.xex[0] = .13205;
    l20_11.xex[1] = 0.;
    l20_11.xex[2] = 0.;
    l20_11.xex[3] = 0.;
    l20_11.xex[4] = 0.;
    l20_11.xex[5] = .32627;
    l20_11.xex[6] = 0.;
    l20_11.xex[7] = 0.;
    l20_11.xex[8] = .51668;
    l20_11.xex[9] = 0.;
    l20_11.xex[10] = 0.;
    l20_11.xex[11] = 0.;
    l20_11.xex[12] = .025004;
    for (i__ = 1; i__ <= 13; ++i__) {
/* L8: */
	l4_10.gf[i__ - 1] = r__[i__ - 1];
    }
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* L7: */
	l6_1.fx += r__[i__ - 1] * l2_10.x[i__ - 1];
    }
labelL3:
    return 0;
labelL4:
    if (! l9_7.index1[0]) {
	goto labelL9;
    }
    l3_6.g[0] = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* labelL10: */
/* Computing 2nd power */
	d__1 = l2_10.x[i__ - 1];
	l3_6.g[0] += z1[i__ - 1] * (d__1 * d__1);
    }
    l3_6.g[0] = -std::sqrt(l3_6.g[0]) * 1.645 - 18.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* labelL11: */
	l3_6.g[0] += s[i__ - 1] * l2_10.x[i__ - 1];
    }
labelL9:
    if (! l9_7.index1[1]) {
	goto labelL12;
    }
    l3_6.g[1] = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* labelL13: */
/* Computing 2nd power */
	d__1 = l2_10.x[i__ - 1];
	l3_6.g[1] += z2[i__ - 1] * (d__1 * d__1);
    }
    l3_6.g[1] = -std::sqrt(l3_6.g[1]) * 1.645 - 1.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* labelL14: */
	l3_6.g[1] += u[i__ - 1] * l2_10.x[i__ - 1];
    }
labelL12:
    if (! l9_7.index1[2]) {
	goto L15;
    }
    l3_6.g[2] = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* L16: */
/* Computing 2nd power */
	d__1 = l2_10.x[i__ - 1];
	l3_6.g[2] += z3[i__ - 1] * (d__1 * d__1);
    }
    l3_6.g[2] = -std::sqrt(l3_6.g[2]) * 1.645 - .9;
    for (i__ = 1; i__ <= 13; ++i__) {
/* L17: */
	l3_6.g[2] += v[i__ - 1] * l2_10.x[i__ - 1];
    }
L15:
    if (! l9_7.index1[3]) {
	goto L18;
    }
    l3_6.g[3] = -1.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* L19: */
	l3_6.g[3] += l2_10.x[i__ - 1];
    }
L18:
    return 0;
labelL5:
    if (! l10_7.index2[0]) {
	goto L27;
    }
    help = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* L21: */
/* Computing 2nd power */
	d__1 = l2_10.x[i__ - 1];
	help += z1[i__ - 1] * (d__1 * d__1);
    }
    help = -.82250000000000001 / std::sqrt(help);
    for (i__ = 1; i__ <= 13; ++i__) {
/* L22: */
	l5_38.gg[(i__ << 2) - 4] = s[i__ - 1] + help * 2. * z1[i__ - 1] * 
		l2_10.x[i__ - 1];
    }
L27:
    if (! l10_7.index2[1]) {
	goto L28;
    }
    help = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* L23: */
/* Computing 2nd power */
	d__1 = l2_10.x[i__ - 1];
	help += z2[i__ - 1] * (d__1 * d__1);
    }
    help = -.82250000000000001 / std::sqrt(help);
    for (i__ = 1; i__ <= 13; ++i__) {
/* L24: */
	l5_38.gg[(i__ << 2) - 3] = u[i__ - 1] + help * 2. * z2[i__ - 1] * 
		l2_10.x[i__ - 1];
    }
L28:
    if (! l10_7.index2[2]) {
	goto labelL20;
    }
    help = 0.;
    for (i__ = 1; i__ <= 13; ++i__) {
/* L25: */
/* Computing 2nd power */
	d__1 = l2_10.x[i__ - 1];
	help += z3[i__ - 1] * (d__1 * d__1);
    }
    help = -.82250000000000001 / std::sqrt(help);
    for (i__ = 1; i__ <= 13; ++i__) {
/* L26: */
	l5_38.gg[(i__ << 2) - 2] = v[i__ - 1] + help * 2. * z3[i__ - 1] * 
		l2_10.x[i__ - 1];
    }
labelL20:
    return 0;
} /* tp382_ */


/* Subroutine */ int tp383_(int *mode)
{
    /* Initialized data */

    static Real a[14] = { 12842.275,634.25,634.25,634.125,1268.,633.875,
	    633.75,1267.,760.05,633.25,1266.25,632.875,394.46,940.838 };
    static Real b[14] = { 25.,26.,26.,27.,28.,29.,30.,32.,33.,34.,35.,
	    37.,38.,36. };
    static Real c__[14] = { 5.47934,.83234,.94749,1.11082,2.64824,
	    1.55868,1.73215,3.90896,2.74284,2.60541,5.96184,3.29522,1.83517,
	    2.81372 };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 14;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 1;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 14; ++i__) {
	l2_19.x[i__ - 1] = .01;
	l12_19.lxu[i__ - 1] = true;
	l11_19.lxl[i__ - 1] = true;
	l13_23.xl[i__ - 1] = 0.;
	l14_23.xu[i__ - 1] = 1. / b[i__ - 1];
/* labelL6: */
	l5_17.gg[i__ - 1] = c__[i__ - 1];
    }
    l20_19.lex = false;
    l20_19.nex = 1;
    l20_19.fex = 7.28566;
    l20_19.xex[0] = .04;
    l20_19.xex[1] = .0382;
    l20_19.xex[2] = .0358;
    l20_19.xex[3] = .033;
    l20_19.xex[4] = .0303;
    l20_19.xex[5] = .0279;
    l20_19.xex[6] = .0265;
    l20_19.xex[7] = .0249;
    l20_19.xex[8] = .023;
    l20_19.xex[9] = .0216;
    l20_19.xex[10] = .0202;
    l20_19.xex[11] = .0192;
    l20_19.xex[12] = .0203;
    l20_19.xex[13] = .0253;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 14; ++i__) {
/* L7: */
	l6_1.fx += a[i__ - 1] / (l2_19.x[i__ - 1] + 1e-16);
    }
    l6_1.fx *= 1e-5;
    return 0;
labelL3:
    for (i__ = 1; i__ <= 14; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__1 = l2_19.x[i__ - 1];
	l4_19.gf[i__ - 1] = -a[i__ - 1] / (d__1 * d__1) * 1e-5;
    }
    return 0;
labelL4:
    if (! l9_2.index1[0]) {
	goto labelL5;
    }
    l3_1.g[0] = 0.;
    for (i__ = 1; i__ <= 14; ++i__) {
/* labelL10: */
	l3_1.g[0] += c__[i__ - 1] * l2_19.x[i__ - 1];
    }
    l3_1.g[0] += -1.;
labelL5:
    return 0;
} /* tp383_ */


/* Subroutine */ int tp384_(int *mode)
{
    /* Initialized data */

    static Real b[10] = { 385.,470.,560.,565.,645.,430.,485.,455.,390.,
	    860. };
    static Real a[150]	/* was [10][15] */ = { 100.,90.,70.,50.,50.,
	    40.,30.,20.,10.,5.,100.,100.,50.,0.,10.,0.,60.,30.,70.,10.,10.,
	    10.,0.,0.,70.,50.,30.,40.,10.,500.,5.,35.,55.,65.,60.,95.,90.,25.,
	    35.,5.,10.,20.,25.,35.,45.,50.,0.,40.,25.,20.,0.,5.,100.,100.,45.,
	    35.,30.,25.,65.,5.,0.,0.,40.,35.,0.,10.,5.,15.,0.,10.,25.,35.,50.,
	    60.,35.,60.,25.,10.,30.,35.,0.,55.,0.,0.,65.,0.,0.,80.,0.,95.,10.,
	    25.,30.,15.,5.,45.,70.,20.,0.,70.,55.,20.,60.,0.,75.,15.,20.,30.,
	    25.,20.,5.,0.,10.,75.,100.,20.,25.,30.,0.,10.,45.,40.,30.,35.,75.,
	    0.,70.,5.,15.,35.,20.,25.,0.,30.,10.,5.,15.,65.,50.,10.,0.,10.,
	    40.,65.,0.,5.,15.,20.,55.,30. };
    static Real d__[15] = { 486.,640.,758.,776.,477.,707.,175.,619.,
	    627.,614.,475.,377.,524.,468.,529. };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real c__;
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 15;
    l1_1.nili = 0;
    l1_1.ninl = 10;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 15; ++i__) {
	l2_11.x[i__ - 1] = 0.;
	l12_11.lxu[i__ - 1] = false;
/* labelL6: */
	l11_11.lxl[i__ - 1] = false;
    }
    l20_12.lex = false;
    l20_12.nex = 1;
    l20_12.fex = -8310.259;
    l20_12.xex[0] = .86095379;
    l20_12.xex[1] = .91736139;
    l20_12.xex[2] = .91973646;
    l20_12.xex[3] = .89600562;
    l20_12.xex[4] = 1.0372946;
    l20_12.xex[5] = .97308908;
    l20_12.xex[6] = .82243629;
    l20_12.xex[7] = 1.1987219;
    l20_12.xex[8] = 1.156335;
    l20_12.xex[9] = 1.1443868;
    l20_12.xex[10] = 1.0305681;
    l20_12.xex[11] = .90949479;
    l20_12.xex[12] = 1.082045;
    l20_12.xex[13] = .84682383;
    l20_12.xex[14] = 1.172372;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 15; ++i__) {
/* L15: */
	l6_1.fx -= d__[i__ - 1] * l2_11.x[i__ - 1];
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 15; ++i__) {
/* labelL11: */
	l4_11.gf[i__ - 1] = -d__[i__ - 1];
    }
    return 0;
labelL4:
    for (i__ = 1; i__ <= 10; ++i__) {
	if (! l9_10.index1[i__ - 1]) {
	    goto L7;
	}
	c__ = 0.;
	for (j = 1; j <= 15; ++j) {
/* labelL9: */
/* Computing 2nd power */
	    d__1 = l2_11.x[j - 1];
	    c__ += a[i__ + j * 10 - 11] * (d__1 * d__1);
	}
	l3_9.g[i__ - 1] = b[i__ - 1] - c__;
L7:
	;
    }
    return 0;
labelL5:
    for (i__ = 1; i__ <= 10; ++i__) {
	if (! l10_10.index2[i__ - 1]) {
	    goto labelL10;
	}
	for (j = 1; j <= 15; ++j) {
/* labelL13: */
	    l5_32.gg[i__ + j * 10 - 11] = a[i__ + j * 10 - 11] * -2. * 
		    l2_11.x[j - 1];
	}
labelL10:
	;
    }
    return 0;
} /* tp384_ */

/* Subroutine */ int tp385_(int *mode)
{
    /* Initialized data */

    static Real b[10] = { 385.,470.,560.,565.,645.,430.,485.,455.,890.,
	    460. };
    static Real a[150]	/* was [10][15] */ = { 100.,90.,70.,50.,50.,
	    40.,30.,20.,10.,5.,100.,100.,50.,0.,10.,0.,60.,30.,70.,10.,10.,
	    10.,0.,0.,70.,50.,30.,40.,10.,100.,5.,35.,55.,65.,60.,95.,90.,25.,
	    35.,5.,10.,20.,25.,35.,45.,50.,0.,40.,25.,20.,0.,5.,100.,100.,45.,
	    35.,30.,25.,65.,5.,0.,0.,40.,35.,0.,10.,5.,15.,0.,10.,25.,35.,50.,
	    60.,35.,60.,25.,10.,30.,35.,0.,55.,0.,0.,65.,0.,0.,80.,500.,95.,
	    10.,25.,30.,15.,5.,45.,70.,20.,0.,70.,55.,20.,60.,0.,75.,15.,20.,
	    30.,25.,20.,5.,0.,10.,75.,100.,20.,25.,30.,0.,10.,45.,40.,30.,35.,
	    75.,0.,70.,5.,15.,35.,20.,25.,0.,30.,10.,5.,15.,65.,50.,10.,0.,
	    10.,40.,65.,0.,5.,15.,20.,55.,30. };
    static Real d__[15] = { 486.,640.,758.,776.,477.,707.,175.,619.,
	    627.,614.,475.,377.,524.,468.,529. };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real c__;
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 15;
    l1_1.nili = 0;
    l1_1.ninl = 10;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 15; ++i__) {
	l2_11.x[i__ - 1] = 0.;
	l12_11.lxu[i__ - 1] = false;
/* labelL6: */
	l11_11.lxl[i__ - 1] = false;
    }
    l20_12.lex = false;
    l20_12.nex = 1;
    l20_12.fex = -8315.2859;
    l20_12.xex[0] = .81347016;
    l20_12.xex[1] = .11327964;
    l20_12.xex[2] = 1.0861184;
    l20_12.xex[3] = .99832977;
    l20_12.xex[4] = 1.0754862;
    l20_12.xex[5] = 1.0688758;
    l20_12.xex[6] = .6278155;
    l20_12.xex[7] = 1.092998;
    l20_12.xex[8] = .91363225;
    l20_12.xex[9] = .86191244;
    l20_12.xex[10] = 1.0047312;
    l20_12.xex[11] = .87742949;
    l20_12.xex[12] = .986715;
    l20_12.xex[13] = 1.0411266;
    l20_12.xex[14] = 1.1860995;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 15; ++i__) {
/* L15: */
	l6_1.fx -= d__[i__ - 1] * l2_11.x[i__ - 1];
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 15; ++i__) {
/* labelL11: */
	l4_11.gf[i__ - 1] = -d__[i__ - 1];
    }
    return 0;
labelL4:
    for (i__ = 1; i__ <= 10; ++i__) {
	if (! l9_10.index1[i__ - 1]) {
	    goto L7;
	}
	c__ = 0.;
	for (j = 1; j <= 15; ++j) {
/* labelL9: */
/* Computing 2nd power */
	    d__1 = l2_11.x[j - 1];
	    c__ += a[i__ + j * 10 - 11] * (d__1 * d__1);
	}
	l3_9.g[i__ - 1] = b[i__ - 1] - c__;
L7:
	;
    }
    return 0;
labelL5:
    for (i__ = 1; i__ <= 10; ++i__) {
	if (! l10_10.index2[i__ - 1]) {
	    goto labelL10;
	}
	for (j = 1; j <= 15; ++j) {
/* labelL13: */
	    l5_32.gg[i__ + j * 10 - 11] = a[i__ + j * 10 - 11] * -2. * 
		    l2_11.x[j - 1];
	}
labelL10:
	;
    }
    return 0;
} /* tp385_ */

/* Subroutine */ int tp386_(int *mode)
{
    /* Initialized data */

    static Real b[10] = { 385.,470.,560.,565.,645.,430.,485.,455.,390.,
	    460. };
    static Real a[150]	/* was [10][15] */ = { 100.,90.,70.,50.,50.,
	    40.,30.,20.,10.,5.,100.,100.,50.,0.,10.,0.,60.,30.,70.,10.,10.,
	    10.,0.,0.,70.,50.,30.,40.,10.,100.,5.,35.,55.,65.,60.,95.,90.,25.,
	    35.,5.,10.,20.,25.,35.,45.,50.,0.,40.,25.,20.,0.,5.,100.,100.,45.,
	    35.,30.,25.,65.,5.,0.,0.,40.,35.,0.,10.,5.,15.,0.,10.,25.,35.,50.,
	    60.,35.,60.,25.,10.,30.,35.,0.,55.,0.,0.,65.,0.,0.,80.,0.,95.,10.,
	    25.,30.,15.,5.,45.,70.,20.,0.,70.,55.,20.,60.,0.,75.,15.,20.,30.,
	    25.,20.,5.,0.,10.,75.,100.,20.,25.,30.,0.,10.,45.,40.,30.,35.,75.,
	    0.,70.,5.,15.,35.,20.,25.,0.,30.,10.,5.,15.,65.,50.,10.,0.,10.,
	    40.,65.,0.,5.,15.,20.,55.,30. };
    static Real d__[15] = { 486.,640.,758.,776.,477.,707.,175.,619.,
	    627.,614.,475.,377.,524.,468.,529. };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real c__;
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 15;
    l1_1.nili = 0;
    l1_1.ninl = 11;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 15; ++i__) {
	l2_11.x[i__ - 1] = 0.;
	l12_11.lxu[i__ - 1] = false;
/* labelL6: */
	l11_11.lxl[i__ - 1] = false;
    }
    l20_12.lex = false;
    l20_12.nex = 1;
    l20_12.fex = -8164.3688;
    l20_12.xex[0] = 1.0042725;
    l20_12.xex[1] = 1.0871174;
    l20_12.xex[2] = 1.10338;
    l20_12.xex[3] = 1.0307192;
    l20_12.xex[4] = .92857958;
    l20_12.xex[5] = 1.2568055;
    l20_12.xex[6] = .76058681;
    l20_12.xex[7] = .85688931;
    l20_12.xex[8] = 1.089778;
    l20_12.xex[9] = .98119425;
    l20_12.xex[10] = .85106387;
    l20_12.xex[11] = .96555941;
    l20_12.xex[12] = .9064419;
    l20_12.xex[13] = .83804049;
    l20_12.xex[14] = .80932365;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 15; ++i__) {
/* labelL20: */
	l6_1.fx -= d__[i__ - 1] * l2_11.x[i__ - 1];
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 15; ++i__) {
/* labelL11: */
	l4_11.gf[i__ - 1] = -d__[i__ - 1];
    }
    return 0;
labelL4:
    for (i__ = 1; i__ <= 10; ++i__) {
	if (! l9_13.index1[i__ - 1]) {
	    goto L7;
	}
	c__ = 0.;
	for (j = 1; j <= 15; ++j) {
/* labelL9: */
/* Computing 2nd power */
	    d__1 = l2_11.x[j - 1];
	    c__ += a[i__ + j * 10 - 11] * (d__1 * d__1);
	}
	l3_12.g[i__ - 1] = b[i__ - 1] - c__;
L7:
	;
    }
    if (! l9_13.index1[10]) {
	goto labelL12;
    }
    c__ = 0.;
    for (j = 1; j <= 15; ++j) {
/* labelL14: */
/* Computing 2nd power */
	d__1 = l2_11.x[j - 1] - 2.;
	c__ += (Real) j * (d__1 * d__1);
    }
    l3_12.g[10] = c__ / 2. - 70.;
labelL12:
    return 0;
labelL5:
    for (i__ = 1; i__ <= 10; ++i__) {
	if (! l10_13.index2[i__ - 1]) {
	    goto labelL10;
	}
	for (j = 1; j <= 15; ++j) {
/* labelL13: */
	    l5_39.gg[i__ + j * 11 - 12] = a[i__ + j * 10 - 11] * -2. * 
		    l2_11.x[j - 1];
	}
labelL10:
	;
    }
    if (! l10_13.index2[10]) {
	goto L15;
    }
    for (j = 1; j <= 15; ++j) {
/* L16: */
	l5_39.gg[j * 11 - 1] = (Real) j * (l2_11.x[j - 1] - 2.);
    }
L15:
    return 0;
} /* tp386_ */

/* Subroutine */ int tp387_(int *mode)
{
    /* Initialized data */

    static Real b[10] = { 385.,470.,560.,565.,645.,430.,485.,455.,390.,
	    460. };
    static Real a[150]	/* was [10][15] */ = { 100.,90.,70.,50.,50.,
	    40.,30.,20.,10.,5.,100.,100.,50.,0.,10.,0.,60.,30.,70.,10.,10.,
	    10.,0.,0.,70.,50.,30.,40.,10.,100.,5.,35.,55.,65.,60.,95.,90.,25.,
	    35.,5.,10.,20.,25.,35.,45.,50.,0.,40.,25.,20.,0.,5.,100.,100.,45.,
	    35.,30.,25.,65.,5.,0.,0.,40.,35.,0.,10.,5.,15.,0.,10.,25.,35.,50.,
	    60.,35.,60.,25.,10.,30.,35.,0.,55.,0.,0.,65.,0.,0.,80.,0.,95.,10.,
	    25.,30.,15.,5.,45.,70.,20.,0.,70.,55.,20.,60.,0.,75.,15.,20.,30.,
	    25.,20.,5.,0.,10.,75.,100.,20.,25.,30.,0.,10.,45.,40.,30.,35.,75.,
	    0.,70.,5.,15.,35.,20.,25.,0.,30.,10.,5.,15.,65.,50.,10.,0.,10.,
	    40.,65.,0.,5.,15.,20.,55.,30. };
    static Real d__[15] = { 486.,640.,758.,776.,477.,707.,175.,619.,
	    627.,614.,475.,377.,524.,468.,529. };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real c__;
    static int i__, j;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 15;
    l1_1.nili = 0;
    l1_1.ninl = 11;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 15; ++i__) {
	l2_11.x[i__ - 1] = 0.;
	l12_11.lxu[i__ - 1] = false;
/* labelL6: */
	l11_11.lxl[i__ - 1] = false;
    }
    l20_12.lex = false;
    l20_12.nex = 1;
    l20_12.fex = -8250.1417;
    l20_12.xex[0] = 1.0125415;
    l20_12.xex[1] = 1.0158505;
    l20_12.xex[2] = 1.0309039;
    l20_12.xex[3] = .99697018;
    l20_12.xex[4] = .98528372;
    l20_12.xex[5] = 1.0368532;
    l20_12.xex[6] = .99349349;
    l20_12.xex[7] = .9720116;
    l20_12.xex[8] = .99994095;
    l20_12.xex[9] = .99547294;
    l20_12.xex[10] = .9695385;
    l20_12.xex[11] = 1.0080569;
    l20_12.xex[12] = .98236999;
    l20_12.xex[13] = .99057993;
    l20_12.xex[14] = .97760168;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 15; ++i__) {
/* labelL20: */
	l6_1.fx -= d__[i__ - 1] * l2_11.x[i__ - 1];
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 15; ++i__) {
/* labelL11: */
	l4_11.gf[i__ - 1] = -d__[i__ - 1];
    }
    return 0;
labelL4:
    for (i__ = 1; i__ <= 10; ++i__) {
	if (! l9_13.index1[i__ - 1]) {
	    goto L7;
	}
	c__ = 0.;
	for (j = 1; j <= 15; ++j) {
/* labelL9: */
/* Computing 2nd power */
	    d__1 = l2_11.x[j - 1];
	    c__ += a[i__ + j * 10 - 11] * (d__1 * d__1);
	}
	l3_12.g[i__ - 1] = b[i__ - 1] - c__;
L7:
	;
    }
    if (! l9_13.index1[10]) {
	goto labelL12;
    }
    c__ = 0.;
    for (j = 1; j <= 15; ++j) {
/* labelL14: */
/* Computing 2nd power */
	d__1 = l2_11.x[j - 1] - 2.;
	c__ += (Real) j * (d__1 * d__1);
    }
    l3_12.g[10] = c__ / 2. - 61.;
labelL12:
    return 0;
labelL5:
    for (i__ = 1; i__ <= 10; ++i__) {
	if (! l10_13.index2[i__ - 1]) {
	    goto labelL10;
	}
	for (j = 1; j <= 15; ++j) {
/* labelL13: */
	    l5_39.gg[i__ + j * 11 - 12] = a[i__ + j * 10 - 11] * -2. * 
		    l2_11.x[j - 1];
	}
labelL10:
	;
    }
    if (! l10_13.index2[10]) {
	goto L15;
    }
    for (j = 1; j <= 15; ++j) {
/* L16: */
	l5_39.gg[j * 11 - 1] = (Real) j * (l2_11.x[j - 1] - 2.);
    }
L15:
    return 0;
} /* tp387_ */

/* Subroutine */ int tp388_(int *mode)
{
    /* Initialized data */

    static Real b[15] = { 385.,470.,560.,565.,645.,430.,485.,455.,390.,
	    460.,0.,70.,361.,265.,395. };
    static Real a[150]	/* was [10][15] */ = { 100.,90.,70.,50.,50.,
	    40.,30.,20.,10.,5.,100.,100.,50.,0.,10.,0.,60.,30.,70.,10.,10.,
	    10.,0.,0.,70.,50.,30.,40.,10.,100.,5.,35.,55.,65.,60.,95.,90.,25.,
	    35.,5.,10.,20.,25.,35.,45.,50.,0.,40.,25.,20.,0.,5.,100.,100.,45.,
	    35.,30.,25.,65.,5.,0.,0.,40.,35.,0.,10.,5.,15.,0.,10.,25.,35.,50.,
	    60.,35.,60.,25.,10.,30.,35.,0.,55.,0.,0.,65.,0.,0.,80.,0.,95.,10.,
	    25.,30.,15.,5.,45.,70.,20.,0.,70.,55.,20.,60.,0.,75.,15.,20.,30.,
	    25.,20.,5.,0.,10.,75.,100.,20.,25.,30.,0.,10.,45.,40.,30.,35.,75.,
	    0.,70.,5.,15.,35.,20.,25.,0.,30.,10.,5.,15.,65.,50.,10.,0.,10.,
	    40.,65.,0.,5.,15.,20.,55.,30. };
    static Real a1[60]	/* was [4][15] */ = { 1.,45.,53.,12.,2.,25.,
	    74.,43.,3.,35.,26.,51.,4.,85.,17.,39.,5.,40.,25.,58.,6.,73.,25.,
	    42.,7.,17.,26.,60.,8.,52.,24.,20.,9.,86.,85.,40.,10.,14.,35.,80.,
	    15.,30.,14.,75.,16.,50.,23.,85.,17.,40.,37.,95.,18.,70.,56.,23.,
	    19.,60.,10.,67. };
    static Real d__[15] = { 486.,640.,758.,776.,477.,707.,175.,619.,
	    627.,614.,475.,377.,524.,468.,529. };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real c__;
    static int i__, j, l;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 15;
    l1_1.nili = 4;
    l1_1.ninl = 11;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 15; ++i__) {
	l2_11.x[i__ - 1] = 0.;
	l12_11.lxu[i__ - 1] = false;
/* labelL6: */
	l11_11.lxl[i__ - 1] = false;
    }
    for (i__ = 12; i__ <= 15; ++i__) {
	for (j = 1; j <= 15; ++j) {
/* L19: */
	    l5_40.gg[i__ - 11 + j * 15 - 16] = -a1[i__ - 11 + (j << 2) - 5];
	}
    }
    l20_12.lex = false;
    l20_12.nex = 1;
    l20_12.fex = -5821.0842;
    l20_12.xex[0] = .62683876;
    l20_12.xex[1] = 1.4330999;
    l20_12.xex[2] = 1.4625963;
    l20_12.xex[3] = .73133338;
    l20_12.xex[4] = .7861424;
    l20_12.xex[5] = 1.2048598;
    l20_12.xex[6] = -1.1433978;
    l20_12.xex[7] = 1.0611103;
    l20_12.xex[8] = -.13389293;
    l20_12.xex[9] = 1.1820107;
    l20_12.xex[10] = .96917757;
    l20_12.xex[11] = -.84501289;
    l20_12.xex[12] = .48122454;
    l20_12.xex[13] = -.33986164;
    l20_12.xex[14] = .68589012;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 15; ++i__) {
/* labelL20: */
	l6_1.fx -= d__[i__ - 1] * l2_11.x[i__ - 1];
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 15; ++i__) {
/* labelL11: */
	l4_11.gf[i__ - 1] = -d__[i__ - 1];
    }
    return 0;
labelL4:
    for (i__ = 12; i__ <= 15; ++i__) {
	l = i__ - 11;
	if (! l9_14.index1[l - 1]) {
	    goto L7;
	}
	c__ = 0.;
	for (j = 1; j <= 15; ++j) {
/* labelL9: */
	    c__ += a1[l + (j << 2) - 5] * l2_11.x[j - 1];
	}
	l3_13.g[l - 1] = b[i__ - 1] - c__;
L7:
	;
    }
    for (i__ = 1; i__ <= 10; ++i__) {
	if (! l9_14.index1[i__ + 3]) {
	    goto labelL14;
	}
	c__ = 0.;
	for (j = 1; j <= 15; ++j) {
/* L16: */
/* Computing 2nd power */
	    d__1 = l2_11.x[j - 1];
	    c__ += a[i__ + j * 10 - 11] * (d__1 * d__1);
	}
	l3_13.g[i__ + 3] = b[i__ - 1] - c__;
labelL14:
	;
    }
    if (! l9_14.index1[14]) {
	goto L17;
    }
    c__ = 0.;
    for (j = 1; j <= 15; ++j) {
/* L18: */
/* Computing 2nd power */
	d__1 = l2_11.x[j - 1] - 2.;
	c__ += (Real) j * (d__1 * d__1);
    }
    l3_13.g[14] = c__ / 2. - 193.121;
L17:
    return 0;
labelL5:
    for (i__ = 1; i__ <= 10; ++i__) {
	if (! l10_14.index2[i__ + 3]) {
	    goto L22;
	}
	for (j = 1; j <= 15; ++j) {
/* L24: */
	    l5_40.gg[i__ + 4 + j * 15 - 16] = a[i__ + j * 10 - 11] * -2. * 
		    l2_11.x[j - 1];
	}
L22:
	;
    }
    if (! l10_14.index2[14]) {
	goto L25;
    }
    for (j = 1; j <= 15; ++j) {
/* L26: */
	l5_40.gg[j * 15 - 1] = (Real) j * (l2_11.x[j - 1] - 2.);
    }
L25:
    return 0;
} /* tp388_ */

/* Subroutine */ int tp389_(int *mode)
{
    /* Initialized data */

    static Real b[15] = { 385.,470.,560.,565.,645.,430.,485.,455.,390.,
	    460.,0.,70.,361.,265.,395. };
    static Real a[150]	/* was [10][15] */ = { 100.,90.,70.,50.,50.,
	    40.,30.,20.,10.,5.,100.,100.,50.,0.,10.,0.,60.,30.,70.,10.,10.,
	    10.,0.,0.,70.,50.,30.,40.,10.,100.,5.,35.,55.,65.,60.,95.,90.,25.,
	    35.,5.,10.,20.,25.,35.,45.,50.,0.,40.,25.,20.,0.,5.,100.,100.,45.,
	    35.,30.,25.,65.,5.,0.,0.,40.,35.,0.,10.,5.,15.,0.,10.,25.,35.,50.,
	    60.,35.,60.,25.,10.,30.,35.,0.,55.,0.,0.,65.,0.,0.,80.,0.,95.,10.,
	    25.,30.,15.,5.,45.,70.,20.,0.,70.,55.,20.,60.,0.,75.,15.,20.,30.,
	    25.,20.,5.,0.,10.,75.,100.,20.,25.,30.,0.,10.,45.,40.,30.,35.,75.,
	    0.,70.,5.,15.,35.,20.,25.,0.,30.,10.,5.,15.,65.,50.,10.,0.,10.,
	    40.,65.,0.,5.,15.,20.,55.,30. };
    static Real a1[60]	/* was [4][15] */ = { 1.,45.,53.,12.,2.,25.,
	    74.,43.,3.,35.,26.,51.,4.,85.,17.,39.,5.,40.,25.,58.,6.,73.,25.,
	    42.,7.,17.,26.,60.,8.,52.,24.,20.,9.,86.,85.,40.,10.,14.,35.,80.,
	    15.,30.,14.,75.,16.,50.,23.,85.,17.,40.,37.,95.,18.,70.,56.,23.,
	    19.,60.,10.,67. };
    static Real d__[15] = { 486.,640.,758.,776.,477.,707.,175.,619.,
	    627.,614.,475.,377.,524.,468.,529. };

    /* System generated locals */
    Real d__1;

    /* Local variables */
    static Real c__;
    static int i__, j, l;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 15;
    l1_1.nili = 4;
    l1_1.ninl = 11;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 15; ++i__) {
	l2_11.x[i__ - 1] = 0.;
	l12_11.lxu[i__ - 1] = false;
/* labelL6: */
	l11_11.lxl[i__ - 1] = false;
    }
    for (i__ = 12; i__ <= 15; ++i__) {
	for (j = 1; j <= 15; ++j) {
/* L19: */
	    l5_40.gg[i__ - 11 + j * 15 - 16] = -a1[i__ - 11 + (j << 2) - 5];
	}
    }
    l20_12.lex = false;
    l20_12.nex = 1;
    l20_12.fex = -5809.7197;
    l20_12.xex[0] = .67105172;
    l20_12.xex[1] = 1.38854;
    l20_12.xex[2] = 1.4676761;
    l20_12.xex[3] = .76023633;
    l20_12.xex[4] = .82935674;
    l20_12.xex[5] = 1.1638523;
    l20_12.xex[6] = -1.257829;
    l20_12.xex[7] = .98193399;
    l20_12.xex[8] = .068416463;
    l20_12.xex[9] = 1.1472773;
    l20_12.xex[10] = .98662969;
    l20_12.xex[11] = -.88834924;
    l20_12.xex[12] = .56465631;
    l20_12.xex[13] = -.58120082;
    l20_12.xex[14] = .72096897;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 15; ++i__) {
/* labelL20: */
	l6_1.fx -= d__[i__ - 1] * l2_11.x[i__ - 1];
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 15; ++i__) {
/* labelL11: */
	l4_11.gf[i__ - 1] = -d__[i__ - 1];
    }
    return 0;
labelL4:
    for (i__ = 12; i__ <= 15; ++i__) {
	l = i__ - 11;
	if (! l9_14.index1[l - 1]) {
	    goto L7;
	}
	c__ = 0.;
	for (j = 1; j <= 15; ++j) {
/* labelL9: */
	    c__ += a1[l + (j << 2) - 5] * l2_11.x[j - 1];
	}
	l3_13.g[l - 1] = b[i__ - 1] - c__;
L7:
	;
    }
    for (i__ = 1; i__ <= 10; ++i__) {
	if (! l9_14.index1[i__ + 3]) {
	    goto labelL14;
	}
	c__ = 0.;
	for (j = 1; j <= 15; ++j) {
/* L16: */
/* Computing 2nd power */
	    d__1 = l2_11.x[j - 1];
	    c__ += a[i__ + j * 10 - 11] * (d__1 * d__1);
	}
	l3_13.g[i__ + 3] = b[i__ - 1] - c__;
labelL14:
	;
    }
    if (! l9_14.index1[14]) {
	goto L17;
    }
    c__ = 0.;
    for (j = 1; j <= 15; ++j) {
/* L18: */
/* Computing 2nd power */
	d__1 = l2_11.x[j - 1] - 2.;
	c__ += (Real) j * (d__1 * d__1);
    }
    l3_13.g[14] = c__ / 2. - 200.;
L17:
    return 0;
labelL5:
    for (i__ = 1; i__ <= 10; ++i__) {
	if (! l10_14.index2[i__ + 3]) {
	    goto L22;
	}
	for (j = 1; j <= 15; ++j) {
/* L24: */
	    l5_40.gg[i__ + 4 + j * 15 - 16] = a[i__ + j * 10 - 11] * -2. * 
		    l2_11.x[j - 1];
	}
L22:
	;
    }
    if (! l10_14.index2[14]) {
	goto L25;
    }
    for (j = 1; j <= 15; ++j) {
/* L26: */
	l5_40.gg[j * 15 - 1] = (Real) j * (l2_11.x[j - 1] - 2.);
    }
L25:
    return 0;
} /* tp389_ */

/* Subroutine */ int tp390_(int *mode)
{
    /* System generated locals */
    Real d__1;

    /* Local variables */
    int tp390a_(Real *, Real *);
    static int i__;
    static Real zi1, zi2, zi3, psi[11];

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL3;
    }
labelL1:
    l1_1.n = 19;
    l1_1.nili = 1;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 11;
    l2_20.x[0] = .02;
    l2_20.x[1] = 4.;
    l2_20.x[2] = 100.;
    l2_20.x[3] = 100.;
    l2_20.x[4] = 15.;
    l2_20.x[5] = 15.;
    l2_20.x[6] = 100.;
    l2_20.x[7] = 1e3;
    l2_20.x[8] = 1e3;
    l2_20.x[9] = 1e3;
    l2_20.x[10] = 9e3;
    l2_20.x[11] = .001;
    l2_20.x[12] = .001;
    l2_20.x[13] = 1.;
    l2_20.x[14] = .001;
    l2_20.x[15] = .001;
    l2_20.x[16] = .1;
    l2_20.x[17] = 8e3;
    l2_20.x[18] = .001;
    for (i__ = 1; i__ <= 19; ++i__) {
	l12_20.lxu[i__ - 1] = true;
	l11_20.lxl[i__ - 1] = true;
	l14_24.xu[i__ - 1] = 1e5;
/* labelL6: */
	l13_24.xl[i__ - 1] = 1e-5;
    }
    for (i__ = 1; i__ <= 2; ++i__) {
	l14_24.xu[i__ - 1] = 50.;
/* L7: */
	l14_24.xu[i__ + 14] = 50.;
    }
    for (i__ = 3; i__ <= 6; ++i__) {
/* labelL9: */
	l14_24.xu[i__ - 1] = 100.;
    }
    for (i__ = 12; i__ <= 15; ++i__) {
/* labelL10: */
	l14_24.xu[i__ - 1] = 1.;
    }
    l20_20.lex = false;
    l20_20.nex = 1;
    l20_20.fex = 24.4724654;
    l20_20.xex[0] = .004473667;
    l20_20.xex[1] = 3.441565;
    l20_20.xex[2] = 99.34824;
    l20_20.xex[3] = 89.130035;
    l20_20.xex[4] = 15.279316;
    l20_20.xex[5] = 15.279316;
    l20_20.xex[6] = 94.726127;
    l20_20.xex[7] = 12304.197;
    l20_20.xex[8] = 12313.263;
    l20_20.xex[9] = 12313.263;
    l20_20.xex[10] = 95905.631;
    l20_20.xex[11] = 1e-5;
    l20_20.xex[12] = 1e-5;
    l20_20.xex[13] = .999989;
    l20_20.xex[14] = 1e-5;
    l20_20.xex[15] = 1e-5;
    l20_20.xex[16] = .1622235;
    l20_20.xex[17] = 8305.1515;
    l20_20.xex[18] = .0014797356;
    return 0;
labelL2:
    d__1 = l2_20.x[15] * 2268. * l2_20.x[0];
    zi1 = pow_dd(&d__1, &c_b3046) * 25.;
    zi2 = l2_20.x[16] * 1.75e5 + pow_dd(&l2_20.x[16], &c_b3047) * 36500.;
    zi3 = l2_20.x[17] * 12.6 + pow_dd(&c_b3048, &c_b3049) * 5.35 / pow_dd(&
	    l2_20.x[17], &c_b3050);
    l6_1.fx = (zi1 + zi2 + zi3 + 10950. + (l2_20.x[0] * (l2_20.x[12] - 
	    l2_20.x[13]) + l2_20.x[1] * (l2_20.x[11] + 1.) - (1. - l2_20.x[18]
	    ) * 3.) * 1150.) * 1.4;
    l6_1.fx /= 1e4;
labelL3:
    return 0;
labelL4:
    if (l9_17.index1[0]) {
	l3_16.g[0] = 1. - l2_20.x[12] - l2_20.x[13];
    }
    tp390a_(l2_20.x, psi);
    for (i__ = 2; i__ <= 12; ++i__) {
	if (! l9_17.index1[i__ - 1]) {
	    goto L8;
	}
	l3_16.g[i__ - 1] = psi[i__ - 2];
L8:
	;
    }
    return 0;
} /* tp390_ */

int tp390a_(Real *x, Real *psi)
{
    /* Local variables */
    static Real test, ak, ck, zj1, zj2, zj3, zj4, zj5, zk7, zj8, xz4, 
	    yz4, zj10, qz12;

    /* Parameter adjustments */
    --psi;
    --x;

    /* Function Body */
    ak = .64749999999999996 / pow_dd(&c_b3053, &c_b3054);
    xz4 = x[3] * std::exp(-ak * x[16]);
    zj1 = -(x[1] * x[13] * xz4 + x[19] * 300.);
    psi[1] = zj1 + x[1] * x[3] - x[2] * x[5] * x[12];
    yz4 = x[7] + (x[3] - xz4) * .5;
    zj2 = -x[13] * x[1] * yz4;
    psi[2] = zj2 + x[1] * x[7] - x[2] * x[9] * x[12];
    zj3 = (1. - x[19]) * -300. + x[6] * 3. * (1. - x[19]) - x[1] * x[14] * 
	    xz4;
    psi[3] = zj3 + x[2] * (x[4] - x[6]) + x[1] * x[6] * x[14];
    zj4 = x[11] * 3. * (1. - x[19]) + x[1] * x[14] * (x[11] - yz4);
    psi[4] = zj4 + x[2] * (x[8] - x[11]);
    zj5 = x[17] * (x[5] * .48 * x[9] / (x[5] + 100.));
    psi[5] = zj5 * -2. + x[2] * (x[4] - x[5]);
    psi[6] = zj5 + x[2] * (x[8] - x[9]) - x[9] * .048 * x[17];
    zk7 = x[1] * (1. - x[13] - x[14]);
    qz12 = x[1] * (1. - x[13] - x[14]) + x[2] * (1. - x[12]);
    psi[7] = -zk7 * xz4 + x[6] * qz12 - x[2] * x[5] * (1. - x[12]);
    zj8 = x[10] * qz12 - zk7 * yz4;
    psi[8] = zj8 - x[2] * x[9] * (1. - x[12]);
    psi[9] = (1. - x[15]) * 6. * (20. - x[6]) + x[11] * (x[2] - (1. - x[15]) *
	     3. - x[1] * x[14]) + x[19] * 3. * x[11] - x[10] * qz12;
    ck = .0013285402597402597;
    test = -ck * x[18] / qz12;
    if (test > 99.) {
	zj10 = std::sqrt((std::abs(x[10]))) * (float)-2.1 * std::exp(99.);
    }
    if (test < 99.) {
	zj10 = std::sqrt((std::abs(x[10]))) * (float)-2.1 * std::exp(-ck * x[18] / qz12);
    }
    psi[10] = zj10 + (20. - x[6]) * 2.;
    psi[11] = (1. - x[13]) * x[1] - x[12] * x[2] - x[19] * 3.;
    return 0;
} /* tp390a_ */

/* Subroutine */ int tp391_(int *mode)
{
    /* System generated locals */
    int i__1;
    Real d__1, d__2, d__3, d__4;

    /* Local variables */
    static Real wurz;
    static int i__, j;
    static Real sum;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL3;
	case 5:  goto labelL3;
    }
labelL1:
    l1_1.n = 30;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 30; ++i__) {
	sum = 0.;
	for (j = 1; j <= 30; ++j) {
	    if (j == i__) {
		goto L60;
	    }
	    wurz = std::sqrt((Real) i__ / (Real) j);
/* Computing 5th power */
	    d__1 = std::sin(std::log(wurz)), d__2 = d__1, d__1 *= d__1;
/* Computing 5th power */
	    d__3 =std::cos(std::log(wurz)), d__4 = d__3, d__3 *= d__3;
	    sum += wurz * (d__2 * (d__1 * d__1) + d__4 * (d__3 * d__3));
L60:
	    ;
	}
/* Computing 3rd power */
	i__1 = i__ - 15;
	l2_14.x[i__ - 1] = ((Real) (i__1 * (i__1 * i__1)) + sum) * 
		-2.8742711;
	l12_14.lxu[i__ - 1] = false;
/* labelL6: */
	l11_14.lxl[i__ - 1] = false;
    }
    l20_15.lex = false;
    l20_15.nex = 1;
    l20_15.fex = 0.;
    l20_15.xex[0] = 7898.8423;
    l20_15.xex[1] = 6316.093;
    l20_15.xex[2] = 4957.3076;
    l20_15.xex[3] = 3806.6336;
    l20_15.xex[4] = 2846.7099;
    l20_15.xex[5] = 2060.1112;
    l20_15.xex[6] = 1429.4605;
    l20_15.xex[7] = 937.41818;
    l20_15.xex[8] = 566.66665;
    l20_15.xex[9] = 299.90186;
    l20_15.xex[10] = 119.82826;
    l20_15.xex[11] = 9.1562939;
    l20_15.xex[12] = -49.399019;
    l20_15.xex[13] = -73.1188559;
    l20_15.xex[14] = -79.2811119;
    l20_15.xex[15] = -85.1607569;
    l20_15.xex[16] = -108.03012;
    l20_15.xex[17] = -165.15913;
    l20_15.xex[18] = -273.81552;
    l20_15.xex[19] = -451.26501;
    l20_15.xex[20] = -714.7715;
    l20_15.xex[21] = -1081.5972;
    l20_15.xex[22] = -1569.0029;
    l20_15.xex[23] = -2194.2479;
    l20_15.xex[24] = -2974.5902;
    l20_15.xex[25] = -3927.2868;
    l20_15.xex[26] = -5069.5936;
    l20_15.xex[27] = -6418.7656;
    l20_15.xex[28] = -7992.0568;
    l20_15.xex[29] = -9806.7206;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 30; ++i__) {
	sum = 0.;
	for (j = 1; j <= 30; ++j) {
	    if (j == i__) {
		goto L70;
	    }
/* Computing 2nd power */
	    d__1 = l2_14.x[j - 1];
	    wurz = std::sqrt(d__1 * d__1 + (Real) i__ / (Real) j);
/* Computing 5th power */
	    d__1 = std::sin(std::log(wurz)), d__2 = d__1, d__1 *= d__1;
/* Computing 5th power */
	    d__3 =std::cos(std::log(wurz)), d__4 = d__3, d__3 *= d__3;
	    sum += wurz * (d__2 * (d__1 * d__1) + d__4 * (d__3 * d__3));
L70:
	    ;
	}
/* L7: */
/* Computing 3rd power */
	i__1 = i__ - 15;
/* Computing 2nd power */
	d__1 = l2_14.x[i__ - 1] * 420. + (Real) (i__1 * (i__1 * i__1)) 
		+ sum;
	l6_1.fx += d__1 * d__1;
    }
labelL3:
    return 0;
} /* tp391_ */


/* Subroutine */ int tp392_(int *mode)
{
    /* Initialized data */

    static Real r1[15]	/* was [3][5] */ = { 1e3,520.,910.,1e3,520.,
	    910.,1e3,520.,1e3,1100.,600.,1e3,1100.,600.,1e3 };
    static Real r2[15]	/* was [3][5] */ = { .3,.1,.2,.3,.1,.2,.3,.1,
	    .2,.3,.1,.2,.3,.1,.2 };
    static Real ka[15]	/* was [3][5] */ = { 120.,65.,105.,150.,65.,
	    105.,150.,80.,120.,170.,80.,120.,170.,80.,120. };
    static Real k1[15]	/* was [3][5] */ = { 150.,75.,140.,150.,75.,
	    140.,150.,75.,140.,170.,90.,150.,170.,90.,150. };
    static Real kp[15]	/* was [3][5] */ = { 160.,75.,140.,160.,75.,
	    140.,160.,75.,140.,180.,90.,150.,180.,90.,150. };
    static Real k3[15]	/* was [3][5] */ = { .02,.01,.015,.2,.1,.15,
	    .25,.1,.15,.25,.15,.15,.25,.15,.15 };
    static Real kl1[15]	/* was [3][5] */ = { .005,.005,.005,.05,.05,
	    .05,.06,.06,.06,.06,.06,.06,.06,.06,.06 };
    static Real kl2[15]	/* was [3][5] */ = { 80.,45.,75.,80.,45.,75.,
	    100.,45.,90.,100.,50.,90.,100.,50.,90. };
    static Real h__[15]	/* was [3][5] */ = { 100.,280.,520.,180.,400.,
	    400.,220.,450.,500.,150.,450.,630.,100.,400.,600. };
    static Real t[9]	/* was [3][3] */ = { .6,.3,.36,.4,.1,.08,.1,
	    .12,.06 };
    static Real b[15]	/* was [3][5] */ = { 170.,170.,180.,170.,170.,
	    180.,170.,170.,180.,170.,170.,180.,170.,170.,180. };

    /* System generated locals */
    int i__1;
    Real d__1, d__2;

    /* Local variables */
    static int i__, j, k, l;
    static Real sum, sum1;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 30;
    l1_1.nili = 45;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 0;
    l2_14.x[0] = 80.;
    l2_14.x[1] = 100.;
    l2_14.x[2] = 400.;
    l2_14.x[3] = 100.;
    l2_14.x[4] = 200.;
    l2_14.x[5] = 200.;
    l2_14.x[6] = 100.;
    l2_14.x[7] = 250.;
    l2_14.x[8] = 400.;
    l2_14.x[9] = 50.;
    l2_14.x[10] = 200.;
    l2_14.x[11] = 500.;
    l2_14.x[12] = 50.;
    l2_14.x[13] = 200.;
    l2_14.x[14] = 500.;
    l2_14.x[15] = 100.;
    l2_14.x[16] = 120.;
    l2_14.x[17] = 410.;
    l2_14.x[18] = 120.;
    l2_14.x[19] = 250.;
    l2_14.x[20] = 250.;
    l2_14.x[21] = 150.;
    l2_14.x[22] = 300.;
    l2_14.x[23] = 410.;
    l2_14.x[24] = 600.;
    l2_14.x[25] = 250.;
    l2_14.x[26] = 510.;
    l2_14.x[27] = 100.;
    l2_14.x[28] = 250.;
    l2_14.x[29] = 510.;
    for (i__ = 1; i__ <= 30; ++i__) {
	l12_14.lxu[i__ - 1] = false;
	l11_14.lxl[i__ - 1] = true;
/* labelL6: */
	l13_17.xl[i__ - 1] = 0.;
    }
    l20_15.lex = false;
    l20_15.nex = 1;
    l20_15.fex = -1698.878;
    l20_15.xex[0] = 99.99;
    l20_15.xex[1] = 142.22;
    l20_15.xex[2] = 519.88;
    l20_15.xex[3] = 136.74;
    l20_15.xex[4] = 103.47;
    l20_15.xex[5] = 399.99;
    l20_15.xex[6] = 191.7;
    l20_15.xex[7] = 1.56;
    l20_15.xex[8] = 500.;
    l20_15.xex[9] = 143.43;
    l20_15.xex[10] = 82.39;
    l20_15.xex[11] = 629.82;
    l20_15.xex[12] = 99.92;
    l20_15.xex[13] = 125.22;
    l20_15.xex[14] = 600.;
    l20_15.xex[15] = 101.85;
    l20_15.xex[16] = 142.25;
    l20_15.xex[17] = 519.88;
    l20_15.xex[18] = 144.58;
    l20_15.xex[19] = 105.73;
    l20_15.xex[20] = 409.59;
    l20_15.xex[21] = 182.01;
    l20_15.xex[22] = 29.34;
    l20_15.xex[23] = 490.52;
    l20_15.xex[24] = 143.43;
    l20_15.xex[25] = 52.43;
    l20_15.xex[26] = 629.7;
    l20_15.xex[27] = 99.92;
    l20_15.xex[28] = 125.12;
    l20_15.xex[29] = 600.;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 5; ++i__) {
	sum = 0.;
	for (j = 1; j <= 3; ++j) {
	    sum1 = 0.;
	    i__1 = i__;
	    for (k = 1; k <= i__1; ++k) {
/* L72: */
		sum1 = sum1 + l2_14.x[j + 12 + k * 3 - 1] - l2_14.x[j - 3 + k 
			* 3 - 1];
	    }
/* L71: */
/* Computing 2nd power */
	    d__1 = l2_14.x[(i__ - 1) * 3 + j - 1];
/* Computing 2nd power */
	    d__2 = l2_14.x[i__ * 3 + 12 + j - 1] - l2_14.x[j + i__ * 3 - 4];
	    sum = sum + l2_14.x[(i__ - 1) * 3 + j - 1] * (r1[j + i__ * 3 - 4] 
		    - ka[j + i__ * 3 - 4]) - d__1 * d__1 * r2[j + i__ * 3 - 4]
		     - l2_14.x[i__ * 3 + 12 + j - 1] * (k1[j + i__ * 3 - 4] + 
		    kp[j + i__ * 3 - 4]) - d__2 * d__2 * (k3[j + i__ * 3 - 4] 
		    + kl1[j + i__ * 3 - 4]) - kl2[j + i__ * 3 - 4] * sum1;
	}
/* L70: */
	l6_1.fx -= sum;
    }
    l6_1.fx *= .001;
labelL3:
    return 0;
labelL4:
    for (i__ = 1; i__ <= 5; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    l = (i__ - 1) * 3 + j;
/* L8: */
	    if (l9_19.index1[l - 1]) {
		l3_18.g[l - 1] = h__[j + i__ * 3 - 4] - l2_14.x[l - 1];
	    }
	}
    }
    for (i__ = 1; i__ <= 5; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    l = (i__ - 1) * 3 + j + 15;
	    if (! l9_19.index1[l - 1]) {
		goto labelL9;
	    }
	    l3_18.g[l - 1] = b[j + i__ * 3 - 4];
	    for (k = 1; k <= 3; ++k) {
/* labelL10: */
		l3_18.g[l - 1] -= t[j + k * 3 - 4] * l2_14.x[i__ * 3 + 12 + k 
			- 1];
	    }
labelL9:
	    ;
	}
    }
    for (i__ = 1; i__ <= 5; ++i__) {
	for (j = 1; j <= 3; ++j) {
	    l = (i__ - 1) * 3 + j + 30;
	    if (! l9_19.index1[l - 1]) {
		goto labelL11;
	    }
	    l3_18.g[l - 1] = 0.;
	    i__1 = i__;
	    for (k = 1; k <= i__1; ++k) {
/* labelL12: */
		l3_18.g[l - 1] = l3_18.g[l - 1] + l2_14.x[k * 3 + 12 + j - 1] 
			- l2_14.x[j - 3 + k * 3 - 1];
	    }
labelL11:
	    ;
	}
    }
labelL5:
    return 0;
} /* tp392_ */


/* Subroutine */ int tp393_(int *mode)
{
    /* Local variables */
    int tp393b_(Real *, Real *);
    static Real c__, e;
    static int i__;
    static Real phi[1];

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 48;
    l1_1.nili = 0;
    l1_1.ninl = 1;
    l1_1.neli = 2;
    l1_1.nenl = 0;
    for (i__ = 1; i__ <= 24; ++i__) {
/* labelL6: */
	l2_21.x[i__ - 1] = 1.;
    }
    for (i__ = 25; i__ <= 30; ++i__) {
/* L7: */
	l2_21.x[i__ - 1] = 1.3;
    }
    for (i__ = 31; i__ <= 48; ++i__) {
/* L8: */
	l2_21.x[i__ - 1] = 1.;
    }
    for (i__ = 1; i__ <= 48; ++i__) {
	l11_21.lxl[i__ - 1] = true;
/* labelL9: */
	l13_25.xl[i__ - 1] = .002;
    }
    for (i__ = 1; i__ <= 24; ++i__) {
	l12_21.lxu[i__ - 1] = true;
/* labelL10: */
	l14_25.xu[i__ - 1] = 2.;
    }
    for (i__ = 25; i__ <= 48; ++i__) {
/* labelL11: */
	l12_21.lxu[i__ - 1] = false;
    }
    l20_21.lex = false;
    l20_21.nex = 1;
    l20_21.fex = .86337998;
    l20_21.xex[0] = 2.;
    l20_21.xex[1] = .002;
    l20_21.xex[2] = 2.;
    l20_21.xex[3] = .0339797;
    l20_21.xex[4] = .01657455;
    l20_21.xex[5] = 2.;
    l20_21.xex[6] = 1.8945347;
    l20_21.xex[7] = .002;
    l20_21.xex[8] = 2.;
    l20_21.xex[9] = .03424074;
    l20_21.xex[10] = .016670308;
    l20_21.xex[11] = 2.;
    l20_21.xex[12] = 2.;
    l20_21.xex[13] = .002;
    l20_21.xex[14] = 2.;
    l20_21.xex[15] = .002;
    l20_21.xex[16] = .002;
    l20_21.xex[17] = 1.988;
    l20_21.xex[18] = 2.;
    l20_21.xex[19] = .002;
    l20_21.xex[20] = 2.;
    l20_21.xex[21] = .002;
    l20_21.xex[22] = .002;
    l20_21.xex[23] = 2.;
    l20_21.xex[24] = 1.0159886;
    l20_21.xex[25] = .002;
    l20_21.xex[26] = 1.003163;
    l20_21.xex[27] = .002;
    l20_21.xex[28] = .002;
    l20_21.xex[29] = .999691944;
    l20_21.xex[30] = 1.11272844;
    l20_21.xex[31] = .002;
    l20_21.xex[32] = 1.1024463;
    l20_21.xex[33] = .002;
    l20_21.xex[34] = .002;
    l20_21.xex[35] = 1.1030764;
    l20_21.xex[36] = .92326572;
    l20_21.xex[37] = .9343325;
    l20_21.xex[38] = .92947437;
    l20_21.xex[39] = .91383802;
    l20_21.xex[40] = .90517162;
    l20_21.xex[41] = .89452569;
    l20_21.xex[42] = 1.174573;
    l20_21.xex[43] = .002;
    l20_21.xex[44] = 1.12080408;
    l20_21.xex[45] = .002;
    l20_21.xex[46] = .002;
    l20_21.xex[47] = 1.1163321536;
    return 0;
labelL2:
    e = 0.;
    for (i__ = 1; i__ <= 12; ++i__) {
	c__ = 1. - l2_21.x[i__ - 1];
/* L100: */
	e += c__ * 10. * c__;
    }
    for (i__ = 25; i__ <= 36; ++i__) {
	c__ = l2_21.x[i__ - 1] - 1.;
/* L120: */
	e += (c__ * 2. * (c__ + std::sqrt(c__ * c__ + .1)) + .1) * 1e3 / 4.;
    }
    for (i__ = 37; i__ <= 42; ++i__) {
	c__ = l2_21.x[i__ - 1] - 1.;
/* L140: */
	e += (c__ * 2. * (c__ + std::sqrt(c__ * c__ + .1)) + .1) * 2e3 / 4.;
    }
    for (i__ = 43; i__ <= 48; ++i__) {
/* L160: */
	e += l2_21.x[i__ - 1] * 100.;
    }
    l6_1.fx = e / 1e3;
labelL3:
    return 0;
labelL4:
    if (! l9_4.index1[0]) {
	goto labelL12;
    }
    tp393b_(l2_21.x, phi);
    l3_3.g[0] = phi[0];
labelL12:
    if (! l9_4.index1[1]) {
	goto labelL14;
    }
    l3_3.g[1] = 12.;
    for (i__ = 1; i__ <= 12; ++i__) {
/* labelL13: */
	l3_3.g[1] -= l2_21.x[i__ - 1];
    }
labelL14:
    if (! l9_4.index1[2]) {
	goto labelL5;
    }
    l3_3.g[2] = 12.;
    for (i__ = 1; i__ <= 12; ++i__) {
/* L15: */
	l3_3.g[2] -= l2_21.x[i__ + 11];
    }
labelL5:
    return 0;
} /* tp393_ */


/* Subroutine */ int tp393b_(Real *x, Real *phi)
{
    /* Initialized data */

    static Real a[18] = { .9,.8,1.1,1.,.7,1.1,1.,1.,1.1,.9,.8,1.2,.9,
	    1.2,1.2,1.,1.,.9 };

    static int i__;
    static Real r__, u[18];
    static int k1, k2, k3;
    static Real alp, sum;

    /* Parameter adjustments */
    --phi;
    --x;

    /* Function Body */
/*     1ST TIER OF GASFIERS */
    for (i__ = 1; i__ <= 6; ++i__) {
	k1 = i__ + 24;
	k2 = i__ + 42;
	k3 = i__ + 12;
	alp = x[k1] * x[k1] * a[i__ - 1] * 2. * x[k2] / (x[k2] + 1.) * x[k3];
/* labelL20: */
	u[i__ - 1] = x[i__] * x[i__] / (x[i__] + alp);
    }
/*     2ND TIER OF GASFIERS */
    for (i__ = 7; i__ <= 12; ++i__) {
	k1 = i__ + 24;
	k2 = i__ + 36;
	k3 = i__ + 12;
	alp = x[k1] * x[k1] * a[i__ - 1] * 2. * x[k2] / (x[k2] + 1.) * x[k3];
	sum = x[i__] + u[i__ - 7];
/* L40: */
	u[i__ - 1] = sum * sum / (sum + alp);
    }
/*     1ST TIER OF METHANATORS */
    for (i__ = 13; i__ <= 15; ++i__) {
	k1 = ((i__ - 10) << 1) + 1;
	k2 = i__ + 24;
	alp = x[k2] * x[k2] * a[i__ - 1];
	sum = u[k1 - 1] + u[k1];
/* L60: */
	u[i__ - 1] = sum * sum / (sum + alp);
    }
/*     2ND TIER OF METHANATORS */
    for (i__ = 16; i__ <= 18; ++i__) {
	k1 = i__ + 24;
	alp = x[k1] * x[k1] * a[i__ - 1];
	sum = u[i__ - 4];
/* L80: */
	u[i__ - 1] = sum * sum / (sum + alp);
    }
    r__ = u[15] + u[16] + u[17];
    phi[1] = 1.5 - r__;
    return 0;
} /* tp393b_ */


/* Subroutine */ int tp394_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 20;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 20; ++i__) {
	l2_13.x[i__ - 1] = 2.;
	l12_13.lxu[i__ - 1] = false;
/* labelL6: */
	l11_13.lxl[i__ - 1] = false;
    }
    l20_14.lex = false;
    l20_14.nex = 1;
    l20_14.fex = 1.9166667;
    l20_14.xex[0] = .9128716;
    l20_14.xex[1] = .4082468;
    l20_14.xex[2] = -1.6746493e-5;
    l20_14.xex[3] = -5.4074613e-6;
    l20_14.xex[4] = 1.9606096e-6;
    l20_14.xex[5] = -8.8626385e-6;
    l20_14.xex[6] = 8.1697576e-6;
    l20_14.xex[7] = -1.4386551e-5;
    l20_14.xex[8] = 2.18312e-5;
    l20_14.xex[9] = -1.3873341e-5;
    l20_14.xex[10] = 1.3498048e-5;
    l20_14.xex[11] = -3.9814429e-6;
    l20_14.xex[12] = -1.1023953e-5;
    l20_14.xex[13] = -1.280983e-5;
    l20_14.xex[14] = 7.9408513e-6;
    l20_14.xex[15] = 2.04589e-5;
    l20_14.xex[16] = 4.5644559e-6;
    l20_14.xex[17] = -9.4429887e-6;
    l20_14.xex[18] = -1.0142804e-5;
    l20_14.xex[19] = -1.3788343e-6;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 20; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__1 = l2_13.x[i__ - 1];
/* Computing 4th power */
	d__2 = l2_13.x[i__ - 1], d__2 *= d__2;
	l6_1.fx += (Real) i__ * (d__1 * d__1 + d__2 * d__2);
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 20; ++i__) {
/* labelL9: */
/* Computing 3rd power */
	d__1 = l2_13.x[i__ - 1];
	l4_13.gf[i__ - 1] = (Real) i__ * (l2_13.x[i__ - 1] * 2. + d__1 *
		 (d__1 * d__1) * 4.);
    }
    return 0;
labelL4:
    if (! l9_2.index1[0]) {
	goto labelL10;
    }
    l3_1.g[0] = 0.;
    for (i__ = 1; i__ <= 20; ++i__) {
/* labelL11: */
/* Computing 2nd power */
	d__1 = l2_13.x[i__ - 1];
	l3_1.g[0] += d__1 * d__1;
    }
    l3_1.g[0] += -1.;
labelL10:
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto labelL12;
    }
    for (i__ = 1; i__ <= 20; ++i__) {
/* labelL13: */
	l5_13.gg[i__ - 1] = l2_13.x[i__ - 1] * 2.;
    }
labelL12:
    return 0;
} /* tp394_ */

/* Subroutine */ int tp395_(int *mode)
{
    /* System generated locals */
    Real d__1, d__2;

    /* Local variables */
    static int i__;

    switch ((int)*mode) {
	case 1:  goto labelL1;
	case 2:  goto labelL2;
	case 3:  goto labelL3;
	case 4:  goto labelL4;
	case 5:  goto labelL5;
    }
labelL1:
    l1_1.n = 50;
    l1_1.nili = 0;
    l1_1.ninl = 0;
    l1_1.neli = 0;
    l1_1.nenl = 1;
    for (i__ = 1; i__ <= 50; ++i__) {
	l2_15.x[i__ - 1] = 2.;
	l12_15.lxu[i__ - 1] = false;
/* labelL6: */
	l11_15.lxl[i__ - 1] = false;
    }
    l20_16.lex = false;
    l20_16.nex = 1;
    l20_16.fex = 1.9166668;
    l20_16.xex[0] = .91285206;
    l20_16.xex[1] = .40829045;
    l20_16.xex[2] = -6.4969989e-6;
    l20_16.xex[3] = -9.9096716e-5;
    l20_16.xex[4] = 1.189129e-4;
    l20_16.xex[5] = -4.6486687e-5;
    l20_16.xex[6] = 5.7605078e-5;
    l20_16.xex[7] = -4.8016383e-5;
    l20_16.xex[8] = 2.5691371e-5;
    l20_16.xex[9] = 1.1670144e-5;
    l20_16.xex[10] = -3.0881321e-5;
    l20_16.xex[11] = 8.7202482e-6;
    l20_16.xex[12] = 1.998037e-5;
    l20_16.xex[13] = -1.2338706e-5;
    l20_16.xex[14] = -1.6390153e-5;
    l20_16.xex[15] = 7.3383634e-6;
    l20_16.xex[16] = 1.686298e-5;
    l20_16.xex[17] = 4.3922807e-6;
    l20_16.xex[18] = -5.8623189e-6;
    l20_16.xex[19] = -2.5188987e-6;
    l20_16.xex[20] = 4.5980202e-6;
    l20_16.xex[21] = 3.2507205e-6;
    l20_16.xex[22] = -6.6596023e-6;
    l20_16.xex[23] = -1.4419491e-5;
    l20_16.xex[24] = -1.2164937e-5;
    l20_16.xex[25] = -3.9129061e-6;
    l20_16.xex[26] = 9.8985037e-7;
    l20_16.xex[27] = 1.4776535e-7;
    l20_16.xex[28] = -6.8312704e-7;
    l20_16.xex[29] = 2.4242977e-6;
    l20_16.xex[30] = 5.3892372e-6;
    l20_16.xex[31] = 2.6662956e-6;
    l20_16.xex[32] = -2.928209e-6;
    l20_16.xex[33] = -3.8338271e-6;
    l20_16.xex[34] = 6.1198364e-7;
    l20_16.xex[35] = 4.367186e-6;
    l20_16.xex[36] = 4.1104627e-6;
    l20_16.xex[37] = 1.4549012e-6;
    l20_16.xex[38] = -1.2562117e-6;
    l20_16.xex[39] = -3.0092086e-6;
    l20_16.xex[40] = -3.8620459e-6;
    l20_16.xex[41] = -4.2627256e-6;
    l20_16.xex[42] = -4.5080325e-6;
    l20_16.xex[43] = -4.4852099e-6;
    l20_16.xex[44] = -3.7953194e-6;
    l20_16.xex[45] = -2.3440318e-6;
    l20_16.xex[46] = -7.4816106e-7;
    l20_16.xex[47] = -5.4626804e-8;
    l20_16.xex[48] = -1.0972677e-6;
    l20_16.xex[49] = -2.131277e-6;
    return 0;
labelL2:
    l6_1.fx = 0.;
    for (i__ = 1; i__ <= 50; ++i__) {
/* L8: */
/* Computing 2nd power */
	d__1 = l2_15.x[i__ - 1];
/* Computing 4th power */
	d__2 = l2_15.x[i__ - 1], d__2 *= d__2;
	l6_1.fx += (Real) i__ * (d__1 * d__1 + d__2 * d__2);
    }
    return 0;
labelL3:
    for (i__ = 1; i__ <= 50; ++i__) {
/* Computing 3rd power */
	d__1 = l2_15.x[i__ - 1];
	l4_15.gf[i__ - 1] = (Real) i__ * (l2_15.x[i__ - 1] * 2. + d__1 *
		 (d__1 * d__1) * 4.);
/* labelL9: */
	l4_15.gf[i__ - 1] = l4_15.gf[i__ - 1];
    }
    return 0;
labelL4:
    if (! l9_2.index1[0]) {
	goto labelL10;
    }
    l3_1.g[0] = 0.;
    for (i__ = 1; i__ <= 50; ++i__) {
/* labelL11: */
/* Computing 2nd power */
	d__1 = l2_15.x[i__ - 1];
	l3_1.g[0] += d__1 * d__1;
    }
    l3_1.g[0] += -1.;
labelL10:
    return 0;
labelL5:
    if (! l10_2.index2[0]) {
	goto labelL12;
    }
    for (i__ = 1; i__ <= 50; ++i__) {
/* labelL13: */
	l5_16.gg[i__ - 1] = l2_15.x[i__ - 1] * 2.;
    }
labelL12:
    return 0;
} /* tp395_ */


Real gleich_(Real *p)
{
    /* System generated locals */
    Real ret_val, d__1;

    /* Local variables */
    static Real a, f, y, eps;
    static Real eps2;

    eps = 1e-5;
    y = *p + 1.;
labelL2:
    f = y - *p - std::atan(1. / y);
    if (std::abs(f) <= eps) {
	goto labelL1;
    }
    a = y * y + 1.;
    a = (a + 1.) / a;
    y -= f / a;
    goto labelL2;
labelL1:
/* Computing 2nd power */
    d__1 = eps;
    eps2 = d__1 * d__1;
    if (y > eps2) {
	ret_val = y;
    } else {
	ret_val = eps2;
    }
    return ret_val;
} /* gleich_ */


/* Subroutine */ int mdnord_(Real *a, Real *b)
{
    Real norint_(Real *);

    *b = norint_(a);
    return 0;
} /* mdnord_ */


Real norint_(Real *x)
{
    /* Initialized data */

    static Real p1[9] = { 3723.5079815548067,7113.6632469540499,
	    6758.2169641104859,4032.2670108300497,1631.7602687537147,
	    456.26145870609263,86.082762211948595,10.064858974909542,
	    .56418958676181361 };
    static Real q1[10] = { 3723.5079815548065,11315.192081854405,
	    15802.535999402043,13349.346561284457,7542.4795101934758,
	    2968.0049014823087,817.62238630454408,153.07771075036222,
	    17.839498439139557,1. };
    static Real p2[6] = { 2.9788656263939929,7.4097406059647418,
	    6.1602098531096305,5.0190497267842675,1.275366644729966,
	    .56418958354775507 };
    static Real q2[7] = { 3.3690752069827528,9.6089653271927879,
	    17.081440747466004,12.048951927855129,9.3960340162350542,
	    2.260528520767327,1. };
    static Real sqrt2 = 1.41421356237390505;
    static Real rsqrtpi = .56418958354775629;

    /* System generated locals */
    Real ret_val, d__1;

    /* Local variables */
    static Real erfc, xabs, arg, erf, arg2;


/*   COMPUTES THE GAUSSIAN NORMAL DISTRIBUTION INTEGRAL */
/*   PRECISION ABOUT 16 DIGITS */

    xabs = std::abs(*x);
    if (xabs > .5) {
	if (xabs > 8.) {
	    if (xabs > 100.) {
		erfc = 0.;
	    } else {
		arg = xabs / sqrt2;
/* Computing 2nd power */
		d__1 = arg;
		erfc = (((((p2[5] * arg + p2[4]) * arg + p2[3]) * arg + p2[2])
			 * arg + p2[1]) * arg + p2[0]) / ((((((arg + q2[5]) * 
			arg + q2[4]) * arg + q2[3]) * arg + q2[2]) * arg + q2[
			1]) * arg + q2[0]) * std::exp(-(d__1 * d__1));
	    }
	} else {
	    arg = xabs / sqrt2;
/* Computing 2nd power */
	    d__1 = arg;
	    erfc = ((((((((p1[8] * arg + p1[7]) * arg + p1[6]) * arg + p1[5]) 
		    * arg + p1[4]) * arg + p1[3]) * arg + p1[2]) * arg + p1[1]
		    ) * arg + p1[0]) / (((((((((arg + q1[8]) * arg + q1[7]) * 
		    arg + q1[6]) * arg + q1[5]) * arg + q1[4]) * arg + q1[3]) 
		    * arg + q1[2]) * arg + q1[1]) * arg + q1[0]) * std::exp(-(d__1 
		    * d__1));
	}
	if (*x < 0.) {
	    ret_val = erfc * .5;
	    return ret_val;
	} else {
	    ret_val = (2. - erfc) * .5;
	    return ret_val;
	}
    } else {
	arg = xabs / sqrt2;
/* Computing 2nd power */
	d__1 = arg;
	arg2 = d__1 * d__1;
	erf = arg * 2. * rsqrtpi * ((((((((((arg2 / 210. - 
		.052631578947368418) * arg2 / 9. + .058823529411764705) * 
		arg2 / 8. - .066666666666666666) * arg2 / 7. + 
		.076923076923076927) * arg2 / 6. - .090909090909090912) * 
		arg2 / 5. + .1111111111111111) * arg2 / 4. - 
		.14285714285714285) * arg2 / 3. + .20000000000000001) * arg2 /
		 2. - .33333333333333331) * arg2 + 1.);
	if (*x >= 0.) {
	    ret_val = (erf + 1.) * .5;
	    return ret_val;
	} else {
	    ret_val = (1. - erf) * .5;
	    return ret_val;
	}
    }
    return ret_val;
} /* norint_ */
