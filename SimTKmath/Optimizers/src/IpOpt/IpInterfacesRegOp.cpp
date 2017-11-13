// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpInterfacesRegOp.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "IpInterfacesRegOp.hpp"
#include "IpRegOptions.hpp"
#include "IpIpoptApplication.hpp"
#include "IpTNLPAdapter.hpp"

namespace SimTKIpopt
{

  void RegisterOptions_Interfaces(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->SetRegisteringCategory("Uncategorized");
    IpoptApplication::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Uncategorized");
    TNLPAdapter::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Uncategorized");
  }

} // namespace Ipopt
