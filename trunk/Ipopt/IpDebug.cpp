// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpDebug.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpoptConfig.h"

#include "IpDebug.hpp"
#include "IpJournalist.hpp"

#ifdef IP_DEBUG

namespace Ipopt
{

  Index DebugJournalistWrapper::indentation_level_ = 0;
  Journalist* DebugJournalistWrapper::jrnl_ = NULL;

  DebugJournalistWrapper::DebugJournalistWrapper(
    std::string func_name,
    Index verbose_level)
      :
      func_name_(func_name),
      verbose_level_(verbose_level),
      method_owner_(NULL)
  {
    if (jrnl_==NULL) {
      verbose_level_ = 0;
      return;
    }
    DebugPrintf(1, "-> Calling to: %s\n", func_name_.c_str());
    if (verbose_level_>0) {
      indentation_level_++;
    }
  }

  DebugJournalistWrapper::DebugJournalistWrapper(
    std::string func_name, Index verbose_level,
    const void* const method_owner)
      :
      func_name_(func_name),
      verbose_level_(verbose_level),
      method_owner_(method_owner)
  {
    if (jrnl_==NULL) {
      verbose_level_ = 0;
      return;
    }
    DebugPrintf(1, "-> Calling to: %s in obj: 0x%x\n", func_name_.c_str(),
                method_owner_);
    if (verbose_level_>0) {
      indentation_level_++;
    }
  }

  DebugJournalistWrapper::~DebugJournalistWrapper()
  {
    if (verbose_level_>0) {
      indentation_level_--;
    }
    if (jrnl_) {
      if (method_owner_ == NULL) {
        DebugPrintf(1, "<- Returning from : %s\n", func_name_.c_str());
      }
      else {
        DebugPrintf(1, "<- Returning from : %s in obj: 0x%x\n",
                    func_name_.c_str(), method_owner_);
      }
    }
  }


  void DebugJournalistWrapper::SetJournalist(Journalist* jrnl)
  {

    if (jrnl == NULL) {
      jrnl_->Printf(J_ERROR, J_DBG,
                    "# Setting Journalist to NULL in DebugJournalistWrapper.\n");
      jrnl_ = NULL;
    }
    else if (!jrnl_) {
      jrnl_ = jrnl;
    }
  }


  void DebugJournalistWrapper::DebugPrintf(Index verbosity, const char* pformat, ...)
  {

    if (Verbosity() >= verbosity) {
      va_list(ap);

      va_start(ap, pformat);

      DBG_ASSERT(jrnl_);
      jrnl_->PrintfIndented(J_ERROR, J_DBG, indentation_level_, "# ");
      jrnl_->VPrintf(J_ERROR, J_DBG, pformat, ap);

      va_end(ap);
    }
  }

} // namespace Ipopt

#endif // #ifdef IP_DEBUG
