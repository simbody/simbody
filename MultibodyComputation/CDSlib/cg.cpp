#include "linemin.h"
#include <cdsMath.h>
#include "cg.h"

#ifdef USE_CDS_NAMESPACE
using CDSMath::sq;
#endif

#include <fixedVector.h>
#include <cdsIostream.h>


typedef FixedVector<double,2> CDSVec2;
double costf(const CDSVec2& x) { return sq(x(0)-1)+(x(0)-1)*x(1)+sq(x(1)); }
void   gradf(const CDSVec2& x,
		   CDSVec2& grad)
{ grad(0) = 2*(x(0)-1) + x(1); grad(1) = (x(0)-1) + 2*x(1); }

class StopCondition {
  double  sqGradTol;
public:
  StopCondition(const double& gradTol) : sqGradTol(sq(gradTol)) {}
	
  bool operator()(      int     iter,
		  const double& cost,
		  const CDSVec2&   x,
		  const CDSVec2&   grad) {
//   cout << "iter: " << iter 
//	  << "  x: " << x
//	  << "  cost: " << cost 
//	  << "  grad: " << sqrt(abs2(grad)) << '\n';
   if (abs2(grad) <= sqGradTol )      return 1;
   return 0;
  }
}; /* StopCondition */
  
namespace CG {
 int test() {
   int exit=0;
   cout << "testing CG...";

   CDSVec2 x; x(0) = 3; x(1) = 3;
   double cost;
   cg(costf,gradf,x,cost,StopCondition(1e-16));
   if ( fabs(x(0)-1.0)>1e-12 || fabs(x(1)-0.0)>1e-12 ) {
     cerr<< "error: x (" << x << ") !=  { 1,0 } \n";
     exit=1;
   }
   cout << (exit?"failed":"ok") << '\n';
   return exit;
 } /* CG::test */
};

