
#include "dihedral.h"

#include "cdsVec3.h"
#include <atomSel.h>

Dihedral::Dihedral(const AtomSel& a0,
		   const AtomSel& a1,
		   const AtomSel& a2,
		   const AtomSel& a3)
  : a0(a0[0]), a1(a1[0]), a2(a2[0]), a3(a3[0])
{
 CDSString msg = "Dihedral: ";
 bool ok=1;

 if ( a0.size()!=1 ) {
   msg += "first atom selection must contain exactly one atom.\n";
   ok=0;
 }
 if ( a1.size()!=1 ) {
   msg += "second atom selection must contain exactly one atom.\n";
   ok=0;
 }
 if ( a2.size()!=1 ) {
   msg += "third atom selection must contain exactly one atom.\n";
   ok=0;
 }
 if ( a3.size()!=1 ) {
   msg += "forth atom selection must contain exactly one atom.\n";
   ok=0;
 }

 if (!ok) {
   //   cout << msg << endl;
   throw CDS::out_of_range(msg);
 } 

} /* constructor */

float_type
Dihedral::value()
{
 CDSVec3 v1 = a0.pos()-a1.pos();
 CDSVec3 v2 = a2.pos()-a1.pos();
 CDSVec3 v3 = a1.pos()-a2.pos();
 CDSVec3 v4 = a3.pos()-a2.pos();

 float_type dp = dot(unitVec(cross(v1,v2)), unitVec(cross(v3,v4)));
 CDSVec3 cp = cross(unitVec(cross(v1,v2)), unitVec(cross(v3,v4)));

 float_type val = acos(dp);

 if (dot(cp,v2)<0.)
   val *= -1;

 return val*180/PI;
} /* value */
