// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpObserver.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpObserver.hpp"

namespace SimTKIpopt
{
#ifdef IP_DEBUG_OBSERVER
  const Index Observer::dbg_verbosity = 0;
  const Index Subject::dbg_verbosity = 0;
#endif

void preventNoSymbolsWarningInIpObserver() {}
} // namespace Ipopt
