
#include "randomNum.h"

#include <math.h>

RandomNum::RandomNum() :
  seed_(314159)         // default used by xplor (in xplor.f)
{
}

static double D2P31M = 2147483647.0; // D2P31M=(2**31) - 1
static double D2P31  = 2147483648.0; // D2P31=(2**31)


void 
RandomNum::setSeed(int newVal)
{
 seed_ = newVal;
} /* setSeed */

int
RandomNum::seed() const
{
  return (int)seed_;
} /* seed */


double 
RandomNum::uniform()
{
 //
 // from IMSL's GGUBFS
 //
 //reference:
 //  LEWIS, P., GOODMAN, A. S., AND MILLER, J. M. 
 //  Pseudo-random number generator for the system/360. 
 //  IBM Syst. J. 8,2(1969), 300-312.
 //

 if (seed_ < 1)  
   seed_=314564;
  
 double newSeed = fmod(16807.0 * seed_, D2P31M);

 seed_ = newSeed;

 return (newSeed / D2P31);
} /* uniform */

