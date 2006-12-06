// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpOptionsList.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpoptConfig.h"
#include "IpOptionsList.hpp"

#ifdef HAVE_CCTYPE
# include <cctype>
#else
# ifdef HAVE_CTYPE_H
#  include <ctype.h>
# else
#  error "don't have header file for ctype"
# endif
#endif

#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

namespace Ipopt
{

  bool OptionsList::SetStringValue(const std::string& tag,
                                   const std::string& value,
                                   bool allow_clobber, /* = true */
                                   bool dont_print /* = false */)
  {
    if (IsValid(reg_options_)) {
      SmartPtr<const RegisteredOption> option = reg_options_->GetOption(tag);

      if (IsNull(option)) {
        if (IsValid(jnlst_)) {
          std::string msg = "Tried to set Option: " + tag;
          msg += ". It is not a valid option. Please check the list of available options.\n";
          jnlst_->Printf(J_ERROR, J_MAIN, msg.c_str());
        }
        //THROW_EXCEPTION(OPTION_INVALID, msg);
        return false;
      }

      if (option->Type() != OT_String) {
        if (IsValid(jnlst_)) {
          std::string msg = "Tried to set Option: " + tag;
          msg += ". It is a valid option, but it is of type ";
          if (option->Type() == OT_Number) {
            msg += " Number";
          }
          else if (option->Type() == OT_Integer) {
            msg += " Integer";
          }
          else {
            msg += " Unknown";
          }
          msg += ", not of type String. Please check the documentation for options.\n";
          jnlst_->Printf(J_ERROR, J_MAIN, msg.c_str());
          option->OutputDescription(*jnlst_);
        }
        //THROW_EXCEPTION(OPTION_INVALID, msg);
        return false;
      }

      if (!option->IsValidStringSetting(value)) {
        if (IsValid(jnlst_)) {
          std::string msg = "Setting: " + value;
          msg += " is not a valid setting for Option: ";
          msg += tag;
          msg += ". Check the option documentation.\n";
          jnlst_->Printf(J_ERROR, J_MAIN, msg.c_str());
          option->OutputDescription(*jnlst_);
        }
        //THROW_EXCEPTION(OPTION_INVALID, msg);
        return false;
      }
    }

    if (!will_allow_clobber(tag)) {
      if (IsValid(jnlst_)) {
        std::string msg = "WARNING: Tried to set option \"" + tag;
        msg += "\" to a value of \"" + value;
        msg += "\",\n         but the previous value is set to disallow clobbering.\n";
        msg += "         The setting will remain as: \"" + tag;
        msg += " " + options_[lowercase(tag)].GetValue();
        msg += "\"\n";
        jnlst_->Printf(J_WARNING, J_MAIN, msg.c_str());
      }
    }
    else {
      //    if (will_allow_clobber(tag)) {
      OptionsList::OptionValue optval(lowercase(value), allow_clobber,
                                      dont_print);
      options_[lowercase(tag)] = optval;
    }
    return true;

    //     std::string msg = "Option: \"" + tag;
    //     msg += " ";
    //     msg += value;
    //     msg += "\" not taken because a value of \n\"" ;
    //     msg += options_[lowercase(tag)].GetValue();
    //     msg += "\" already exists and is set to disallow clobbering.\n\n";
    //     jnlst_->Printf(J_ERROR, J_MAIN, msg.c_str());
    //     return false;
  }

  bool OptionsList::SetNumericValue(const std::string& tag, Number value,
                                    bool allow_clobber, /* = true */
                                    bool dont_print /* = false */)
  {
    char buffer[256];
    sprintf(buffer, "%g", value);

    if (IsValid(reg_options_)) {
      SmartPtr<const RegisteredOption> option = reg_options_->GetOption(tag);

      if (IsNull(option)) {
        if (IsValid(jnlst_)) {
          std::string msg = "Tried to set Option: " + tag;
          msg += ". It is not a valid option. Please check the list of available options.\n";
          jnlst_->Printf(J_ERROR, J_MAIN, msg.c_str());
        }
        //THROW_EXCEPTION(OPTION_INVALID, msg);
        return false;
      }

      if (option->Type() != OT_Number) {
        if (IsValid(jnlst_)) {
          std::string msg = "Tried to set Option: " + tag;
          msg += ". It is a valid option, but it is of type ";
          if (option->Type() == OT_String) {
            msg += " String";
          }
          else if (option->Type() == OT_Integer) {
            msg += " Integer";
          }
          else {
            msg += " Unknown";
          }
          msg += ", not of type Number. Please check the documentation for options.\n";
          jnlst_->Printf(J_ERROR, J_MAIN, msg.c_str());
          option->OutputDescription(*jnlst_);
        }
        //THROW_EXCEPTION(OPTION_INVALID, msg);
        return false;
      }

      if (!option->IsValidNumberSetting(value)) {
        if (IsValid(jnlst_)) {
          std::string msg = "Setting: ";
          msg += buffer;
          msg += " is not a valid setting for Option: ";
          msg += tag;
          msg += ". Check the option documentation.\n";
          jnlst_->Printf(J_ERROR, J_MAIN, msg.c_str());
          option->OutputDescription(*jnlst_);
        }
        //THROW_EXCEPTION(OPTION_INVALID, msg);
        return false;
      }
    }

    if (!will_allow_clobber(tag)) {
      if (IsValid(jnlst_)) {
        std::string msg = "WARNING: Tried to set option \"" + tag;
        msg += "\" to a value of \"";
        msg += buffer;
        msg += "\",\n         but the previous value is set to disallow clobbering.\n";
        msg += "         The setting will remain as: \"" + tag;
        msg += " " + options_[lowercase(tag)].GetValue();
        msg += "\"\n";
        jnlst_->Printf(J_WARNING, J_MAIN, msg.c_str());
      }
    }
    else {
      OptionsList::OptionValue optval(buffer, allow_clobber, dont_print);
      options_[lowercase(tag)] = optval;
    }
    return true;
  }

  bool OptionsList::SetIntegerValue(const std::string& tag, Index value,
                                    bool allow_clobber, /* = true */
                                    bool dont_print /* = false */)
  {
    char buffer[256];
    sprintf(buffer, "%d", value);

    if (IsValid(reg_options_)) {
      SmartPtr<const RegisteredOption> option = reg_options_->GetOption(tag);

      if (IsNull(option)) {
        std::string msg = "Tried to set Option: " + tag;
        msg += ". It is not a valid option. Please check the list of available options.\n";
        if (IsValid(jnlst_)) {
          jnlst_->Printf(J_ERROR, J_MAIN, msg.c_str());
        }
        //THROW_EXCEPTION(OPTION_INVALID, msg);
        return false;
      }

      if (option->Type() != OT_Integer) {
        if (IsValid(jnlst_)) {
          std::string msg = "Tried to set Option: " + tag;
          msg += ". It is a valid option, but it is of type ";
          if (option->Type() == OT_String) {
            msg += " String";
          }
          else if (option->Type() == OT_Number) {
            msg += " Number";
          }
          else {
            msg += " Unknown";
          }
          msg += ", not of type Integer. Please check the documentation for options.\n";
          jnlst_->Printf(J_ERROR, J_MAIN, msg.c_str());
          option->OutputDescription(*jnlst_);
        }
        //THROW_EXCEPTION(OPTION_INVALID, msg);
        return false;
      }

      if (!option->IsValidIntegerSetting(value)) {
        if (IsValid(jnlst_)) {
          std::string msg = "Setting: ";
          msg += buffer;
          msg += " is not a valid setting for Option: ";
          msg += tag;
          msg += ". Check the option documentation.\n";
          jnlst_->Printf(J_ERROR, J_MAIN, msg.c_str());
          option->OutputDescription(*jnlst_);
        }
        //THROW_EXCEPTION(OPTION_INVALID, msg);
        return false;
      }
    }

    if (!will_allow_clobber(tag)) {
      if (IsValid(jnlst_)) {
        std::string msg = "WARNING: Tried to set option \"" + tag;
        msg += "\" to a value of \"";
        msg += buffer;
        msg += "\",\n         but the previous value is set to disallow clobbering.\n";
        msg += "         The setting will remain as: \"" + tag;
        msg += " " + options_[lowercase(tag)].GetValue();
        msg += "\"\n";
        jnlst_->Printf(J_WARNING, J_MAIN, msg.c_str());
      }
    }
    else {
      //    if (will_allow_clobber(tag)) {
      OptionsList::OptionValue optval(buffer, allow_clobber, dont_print);
      options_[lowercase(tag)] = optval;
    }
    return true;
  }

  bool OptionsList::GetStringValue(const std::string& tag, std::string& value,
                                   const std::string& prefix) const
  {
    SmartPtr<const RegisteredOption> option = NULL;

    bool found = find_tag(tag, prefix, value);

    if (IsValid(reg_options_)) {
      option = reg_options_->GetOption(tag);
      if (IsNull(option)) {
        std::string msg = "IPOPT tried to get the value of Option: " + tag;
        msg += ". It is not a valid registered option.";
        THROW_EXCEPTION(OPTION_INVALID, msg);
      }

      if (option->Type() != OT_String) {
        std::string msg = "IPOPT tried to get the value of Option: " + tag;
        msg += ". It is a valid option, but it is of type ";
        if (option->Type() == OT_Integer) {
          msg += " Integer";
        }
        else if (option->Type() == OT_Number) {
          msg += " Number";
        }
        else {
          msg += " Unknown";
        }
        msg += ", not of type String. Please check the documentation for options.";
        if (IsValid(jnlst_)) {
          option->OutputDescription(*jnlst_);
        }
        THROW_EXCEPTION(OPTION_INVALID, msg);
      }

      if (found) {
        value = option->MapStringSetting(value);
      }
      else {
        value = option->DefaultString();
      }
    }

    return found;
  }

  bool OptionsList::GetEnumValue(const std::string& tag, Index& value,
                                 const std::string& prefix) const
  {
    std::string str;
    SmartPtr<const RegisteredOption> option = NULL;

    bool found = find_tag(tag, prefix, str);

    if (IsValid(reg_options_)) {
      option = reg_options_->GetOption(tag);
      if (IsNull(option)) {
        std::string msg = "IPOPT tried to get the value of Option: " + tag;
        msg += ". It is not a valid registered option.";
        THROW_EXCEPTION(OPTION_INVALID, msg);
      }

      if (option->Type() != OT_String) {
        std::string msg = "IPOPT tried to get the value of Option: " + tag;
        msg += ". It is a valid option, but it is of type ";
        if (option->Type() == OT_Integer) {
          msg += " Integer";
        }
        else if (option->Type() == OT_Number) {
          msg += " Number";
        }
        else {
          msg += " Unknown";
        }
        msg += ", not of type String. Please check the documentation for options.";
        if (IsValid(jnlst_)) {
          option->OutputDescription(*jnlst_);
        }
        THROW_EXCEPTION(OPTION_INVALID, msg);
      }

      if (found) {
        value = option->MapStringSettingToEnum(str);
      }
      else {
        value = option->DefaultStringAsEnum();
      }
    }

    return found;
  }

  bool OptionsList::GetBoolValue(const std::string& tag, bool& value,
                                 const std::string& prefix) const
  {
    std::string str;
    bool ret = GetStringValue(tag, str, prefix);
    if (str == "no" || str == "false" || str == "off") {
      value = false;
    }
    else if (str == "yes" || str == "true" || str == "on") {
      value = true;
    }
    else {
      THROW_EXCEPTION(OPTION_INVALID, "Tried to get a boolean from an option and failed.");
    }

    return ret;
  }

  bool OptionsList::GetNumericValue(const std::string& tag, Number& value,
                                    const std::string& prefix) const
  {
    SmartPtr<const RegisteredOption> option = NULL;

    if (IsValid(reg_options_)) {
      option = reg_options_->GetOption(tag);
      if (IsNull(option)) {
        std::string msg = "IPOPT tried to get the value of Option: " + tag;
        msg += ". It is not a valid registered option.";
        THROW_EXCEPTION(OPTION_INVALID, msg);
      }

      if (option->Type() != OT_Number) {
        std::string msg = "IPOPT tried to get the value of Option: " + tag;
        msg += ". It is a valid option, but it is of type ";
        if (option->Type() == OT_Integer) {
          msg += " Integer";
        }
        else if (option->Type() == OT_String) {
          msg += " String";
        }
        else {
          msg += " Unknown";
        }
        msg += ", not of type Number. Please check the documentation for options.";
        if (IsValid(jnlst_)) {
          option->OutputDescription(*jnlst_);
        }
        THROW_EXCEPTION(OPTION_INVALID, msg);
      }
    }

    std::string strvalue;
    if (find_tag(tag, prefix, strvalue)) {
      char* p_end;
      Number retval = strtod(strvalue.c_str(), &p_end);
      if (*p_end!='\0' && !isspace(*p_end)) {
        std::string msg = "Option \"" + tag +
                          "\": Double value expected, but non-numeric value \"" +
                          strvalue+"\" found.\n";
        THROW_EXCEPTION(OPTION_INVALID, msg);
      }
      value = retval;
      return true;
    }
    else if (IsValid(option)) {
      value = option->DefaultNumber();
      return false;
    }
    return false;
  }

  bool OptionsList::GetIntegerValue(const std::string& tag, Index& value,
                                    const std::string& prefix) const
  {
    SmartPtr<const RegisteredOption> option = NULL;

    if (IsValid(reg_options_)) {
      option = reg_options_->GetOption(tag);
      if (IsNull(option)) {
        std::string msg = "IPOPT tried to get the value of Option: " + tag;
        msg += ". It is not a valid registered option.";
        THROW_EXCEPTION(OPTION_INVALID, msg);
      }

      if (option->Type() != OT_Integer) {
        std::string msg = "IPOPT tried to get the value of Option: " + tag;
        msg += ". It is a valid option, but it is of type ";
        if (option->Type() == OT_Number) {
          msg += " Number";
        }
        else if (option->Type() == OT_String) {
          msg += " String";
        }
        else {
          msg += " Unknown";
        }
        msg += ", not of type Integer. Please check the documentation for options.";
        if (IsValid(jnlst_)) {
          option->OutputDescription(*jnlst_);
        }
        THROW_EXCEPTION(OPTION_INVALID, msg);
      }
    }

    std::string strvalue;
    if (find_tag(tag, prefix, strvalue)) {
      char* p_end;
      Index retval = strtol(strvalue.c_str(), &p_end, 10);
      if (*p_end!='\0' && !isspace(*p_end)) {
        std::string msg = "Option \"" + tag +
                          "\": Integer value expected, but non-integer value \"" +
                          strvalue+"\" found.\n";
        THROW_EXCEPTION(OPTION_INVALID, msg);
      }
      value = retval;
      return true;
    }
    else if (IsValid(option)) {
      value = option->DefaultInteger();
      return false;
    }

    return false;
  }

  const std::string& OptionsList::lowercase(const std::string tag) const
  {
    lowercase_buffer_ = tag;
    for(Index i=0; i<(Index)tag.length(); i++) {
      lowercase_buffer_[i] = tolower(tag[i]);
    }
    return lowercase_buffer_;
  }

  void OptionsList::PrintList(std::string& list) const
  {
    list.clear();
    char buffer[256];
    sprintf(buffer, "%40s   %-20s %s\n", "Name", "Value", "# times used");
    list += buffer;
    for(std::map< std::string, OptionValue >::const_iterator p = options_.begin();
        p != options_.end();
        p++ ) {
      sprintf(buffer, "%40s = %-20s %6d\n", p->first.c_str(),
              p->second.Value().c_str(), p->second.Counter());
      list += buffer;
    }
  }

  void OptionsList::PrintUserOptions(std::string& list) const
  {
    list.clear();
    char buffer[256];
    sprintf(buffer, "%40s   %-20s %s\n", "Name", "Value", "used");
    list += buffer;
    for(std::map< std::string, OptionValue >::const_iterator p = options_.begin();
        p != options_.end();
        p++ ) {
      if (!p->second.DontPrint()) {
        const char yes[] = "yes";
        const char no[] = "no";
        const char* used;
        if (p->second.Counter()>0) {
          used = yes;
        }
        else {
          used = no;
        }
        sprintf(buffer, "%40s = %-20s %4s\n", p->first.c_str(),
                p->second.Value().c_str(), used);
        list += buffer;
      }
    }
  }

  bool OptionsList::ReadFromStream(const Journalist& jnlst,
                                   std::istream& is)
  {
    jnlst.Printf(J_DETAILED, J_MAIN, "Start reading options from stream.\n");

    while (true) {
      std::string tag;
      std::string value;

      if (!readnexttoken(is, tag)) {
        // That's it - end of file reached.
        jnlst.Printf(J_DETAILED, J_MAIN,
                     "Finished reading options from file.\n");
        return true;
      }

      if (!readnexttoken(is, value)) {
        // Can't read value for a given tag
        jnlst.Printf(J_ERROR, J_MAIN,
                     "Error reading value for tag %s from file.\n",
                     tag.c_str());
        return false;
      }

      // Now add the value for the options list
      jnlst.Printf(J_DETAILED, J_MAIN,
                   "Adding option \"%s\" with value \"%s\" to OptionsList.\n",
                   tag.c_str(), value.c_str());

      if (IsValid(reg_options_)) {
        SmartPtr<const RegisteredOption> option = reg_options_->GetOption(tag);
        if (IsNull(option)) {
          std::string msg = "Read Option: ";
          msg += tag;
          msg += ". It is not a valid option. Check the list of available options.";
          THROW_EXCEPTION(OPTION_INVALID, msg);
        }

        if (option->Type() == OT_String) {
          bool result = SetStringValue(tag, value, false);
          ASSERT_EXCEPTION(result, OPTION_INVALID,
                           "Error setting string value read from option file.");
        }
        else if (option->Type() == OT_Number) {
          char* p_end;
          Number retval = strtod(value.c_str(), &p_end);
          if (*p_end!='\0' && !isspace(*p_end)) {
            std::string msg = "Option \"" + tag +
                              "\": Double value expected, but non-numeric option value \"" +
                              value + "\" found.\n";
            THROW_EXCEPTION(OPTION_INVALID, msg);
          }
          bool result = SetNumericValue(tag, retval, false);
          ASSERT_EXCEPTION(result, OPTION_INVALID,
                           "Error setting numeric value read from file.");
        }
        else if (option->Type() == OT_Integer) {
          char* p_end;
          Index retval = strtol(value.c_str(), &p_end, 10);
          if (*p_end!='\0' && !isspace(*p_end)) {
            std::string msg = "Option \"" + tag +
                              "\": Integer value expected, but non-integer option value \"" +
                              value + "\" found.\n";
            if (IsValid(jnlst_)) {
              option->OutputDescription(*jnlst_);
            }
            THROW_EXCEPTION(OPTION_INVALID, msg);
          }
          bool result = SetIntegerValue(tag, retval, false);
          ASSERT_EXCEPTION(result, OPTION_INVALID,
                           "Error setting integer value read from option file.");
        }
        else {
          DBG_ASSERT(false && "Option Type: Unknown");
        }
      }
      else {
        bool result = SetStringValue(tag, value, false);
        ASSERT_EXCEPTION(result, OPTION_INVALID,
                         "Error setting value read from option file.");
      }
    }
  }

  bool OptionsList::find_tag(const std::string& tag,
                             const std::string& prefix,
                             std::string& value) const
  {
    bool found=false;
    std::map< std::string, OptionValue >::const_iterator p;

    if (prefix != "") {
      p = options_.find(lowercase(prefix+tag));
      if (p != options_.end()) {
        found = true;
      }
    }

    if (!found) {
      p = options_.find(lowercase(tag));
      if (p != options_.end()) {
        found = true;
      }
    }

    if (found) {
      value = p->second.GetValue();
    }

    return found;
  }

  bool OptionsList::will_allow_clobber(const std::string& tag) const
  {
    bool allow_clobber=true;
    std::map< std::string, OptionValue >::const_iterator p;

    p = options_.find(lowercase(tag));
    if (p != options_.end()) {
      allow_clobber = p->second.AllowClobber();
    }

    return allow_clobber;
  }

  bool OptionsList::readnexttoken(std::istream& is, std::string& token)
  {
    token.clear();
    int c = is.get();

    // First get rid of all comments and white spaces
    while (!is.eof() && (isspace(c) || c=='#') ) {
      if (c=='#') {
        is.ignore(10000000, '\n');
      }
      c=is.get();
    }

    // Now read the token
    while (!is.eof() && !isspace(c)) {
      token += c;
      c = is.get();
    }

    return (!is.eof());
  }

} // namespace Ipopt

