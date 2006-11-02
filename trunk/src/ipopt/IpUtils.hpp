// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpUtils.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPUTILS_HPP__
#define __IPUTILS_HPP__

// Standard Ip Include Files
#include "IpTypes.hpp"
#include "IpDebug.hpp"

namespace Ipopt
{

  inline ipfint Max(ipfint a, ipfint b)
  {
    return ((a) > (b) ? (a) : (b));
  }

  inline ipfint Min(ipfint a, ipfint b)
  {
    return ((a) < (b) ? (a) : (b));
  }

  inline Number Max(Number a, Number b)
  {
    return ((a) > (b) ? (a) : (b));
  }

  inline Number Max(Number a, Number b, Number c)
  {
    Number max = Max(a,b);
    max = Max(max, c);
    return max;
  }

  inline Number Max(Number a, Number b, Number c, Number d)
  {
    Number max = Max(a, b, c);
    max = Max(max, d);
    return max;
  }

  inline Number Min(Number a, Number b)
  {
    return ((a) < (b) ? (a) : (b));
  }

  inline Number Min(Number a, Number b, Number c)
  {
    Number min = Min(a,b);
    min = Min(min, c);
    return min;
  }

  inline Number Min(Number a, Number b, Number c, Number d)
  {
    Number min = Min(a, b, c);
    min = Min(min, d);
    return min;
  }

  /** Function returning true iff the argument is a valid double number
   *  (not NaN or Inf). */
  bool IsFiniteNumber(Number val);

} //namespace Ipopt

#endif
