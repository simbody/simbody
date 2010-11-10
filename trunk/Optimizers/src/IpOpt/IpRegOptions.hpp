// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpRegOptions.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-06-18

#ifndef __IPREGOPTIONS_HPP__
#define __IPREGOPTIONS_HPP__

#include "IpUtils.hpp"
#include "IpReferenced.hpp"
#include "IpException.hpp"
#include "IpSmartPtr.hpp"

#include <map>

namespace Ipopt
{

  enum RegisteredOptionType
  {
    OT_Number,
    OT_Integer,
    OT_String,
    OT_Unknown
  };

  /** Base class for registered options. The derived types are more
   *  specific to a string option or a Number (real) option, etc.
   */
  class RegisteredOption : public ReferencedObject
  {
  public:
    /** Constructors / Destructors */
    //@{
    RegisteredOption()
        :
        type_(OT_Unknown),
        has_lower_(false),
        has_upper_(false),
        counter_(0)
    {}

    RegisteredOption(const std::string& name,
                     const std::string& short_description,
                     const std::string& long_description,
                     const std::string& registering_category)
        :
        name_(name),
        short_description_(short_description),
        long_description_(long_description),
        registering_category_(registering_category),
        type_(OT_Unknown),
        has_lower_(false),
        has_upper_(false),
        counter_(next_counter_++)
    {}

    RegisteredOption(const RegisteredOption& copy)
        :
        name_(copy.name_),
        short_description_(copy.short_description_),
        long_description_(copy.long_description_),
        registering_category_(copy.registering_category_),
        type_(copy.type_),
        has_lower_(copy.has_lower_),
        lower_(copy.lower_),
        has_upper_(copy.has_upper_),
        upper_(copy.upper_),
        valid_strings_(copy.valid_strings_),
        counter_(copy.counter_)
    {}

    virtual ~RegisteredOption()
    {}
    //@}

    DECLARE_STD_EXCEPTION(ERROR_CONVERTING_STRING_TO_ENUM);

    /** Standard Get / Set Methods */
    //@{
    /** Get the option's name (tag in the input file) */
    const std::string& Name() const
    {
      return name_;
    }
    /** Set the option's name (tag in the input file) */
    void SetName(const std::string& name)
    {
      name_ = name;
    }
    /** Get the short description */
    const std::string& ShortDescription() const
    {
      return short_description_;
    }
    /** Get the long description */
    const std::string& LongDescription() const
    {
      return long_description_;
    }
    /** Set the short description */
    void SetShortDescription(const std::string& short_description)
    {
      short_description_ = short_description;
    }
    /** Set the long description */
    void SetLongDescription(const std::string& long_description)
    {
      long_description_ = long_description;
    }
    /** Get the registering class */
    const std::string& RegisteringCategory() const
    {
      return registering_category_;
    }
    /** Set the registering class */
    void SetRegisteringCategory(const std::string& registering_category)
    {
      registering_category_ = registering_category;
    }
    /** Get the Option's type */
    const RegisteredOptionType& Type() const
    {
      return type_;
    }
    /** Get the Option's type */
    void SetType(const RegisteredOptionType& type)
    {
      type_ = type;
    }
    /** Counter */
    Index Counter() const
    {
      return counter_;
    }
    //@}

    /** @name Get / Set methods valid for specific types - NOTE: the Type
     *  must be set before calling these methods.
     */
    //@{
    /** check if the option has a lower bound - can be called for
     *  OT_Number & OT_Integer*/
    const bool& HasLower() const
    {
      DBG_ASSERT(type_ == OT_Number || type_ == OT_Integer);
      return has_lower_;
    }
    /** check if the lower bound is strict - can be called for
    OT_Number */
    const bool& LowerStrict() const
    {
      DBG_ASSERT(type_ == OT_Number && has_lower_ == true);
      return lower_strict_;
    }
    /** get the Number version of the lower bound - can be called for
     *  OT_Number */
    Number LowerNumber() const
    {
      DBG_ASSERT(has_lower_ == true && type_ == OT_Number);
      return lower_;
    }
    /** set the Number version of the lower bound - can be called for
     *  OT_Number */
    void SetLowerNumber(const Number& lower, const bool& strict)
    {
      DBG_ASSERT(type_ == OT_Number);
      lower_ = lower;
      lower_strict_ = strict, has_lower_ = true;
    }
    /** get the Integer version of the lower bound can be called for
     *  OT_Integer*/
    Index LowerInteger() const
    {
      DBG_ASSERT(has_lower_ == true && type_ == OT_Integer);
      return (Index)lower_;
    }
    /** set the Integer version of the lower bound - can be called for
     *  OT_Integer */
    void SetLowerInteger(const Index& lower)
    {
      DBG_ASSERT(type_ == OT_Integer);
      lower_ = (Number)lower;
      has_lower_ = true;
    }
    /** check if the option has an upper bound - can be called for
     *  OT_Number & OT_Integer*/
    const bool& HasUpper() const
    {
      DBG_ASSERT(type_ == OT_Number || type_ == OT_Integer);
      return has_upper_;
    }
    /** check if the upper bound is strict - can be called for
     *  OT_Number */
    const bool& UpperStrict() const
    {
      DBG_ASSERT(type_ == OT_Number && has_upper_ == true);
      return upper_strict_;
    }
    /** get the Number version of the upper bound - can be called for
     *  OT_Number */
    Number UpperNumber()
    {
      DBG_ASSERT(has_upper_ == true && type_ == OT_Number);
      return upper_;
    }
    /** set the Number version of the upper bound - can be called for
     *  OT_Number */
    void SetUpperNumber(const Number& upper, const bool& strict)
    {
      DBG_ASSERT(type_ == OT_Number);
      upper_ = upper;
      upper_strict_ = strict;
      has_upper_ = true;
    }
    /** get the Integer version of the upper bound - can be called for
     *  OT_Integer*/
    Index UpperInteger() const
    {
      DBG_ASSERT(has_upper_ == true && type_ == OT_Integer);
      return (Index)upper_;
    }
    /** set the Integer version of the upper bound - can be called for
     *  OT_Integer */
    void SetUpperInteger(const Index& upper)
    {
      DBG_ASSERT(type_ == OT_Integer);
      upper_ = (Number)upper;
      has_upper_ = true;
    }
    /** method to add valid string entries - can be called for
     *  OT_String */
    void AddValidStringSetting(const std::string value,
                               const std::string description)
    {
      DBG_ASSERT(type_ == OT_String);
      valid_strings_.push_back(string_entry(value, description));
    }
    /** get the default as a Number - can be called for OT_Number */
    Number DefaultNumber() const
    {
      DBG_ASSERT(type_ == OT_Number);
      return default_number_;
    }
    /** Set the default as a Number - can be called for OT_Number */
    void SetDefaultNumber(const Number& default_value)
    {
      DBG_ASSERT(type_ == OT_Number);
      default_number_ = default_value;
    }
    /** get the default as an Integer - can be called for OT_Integer*/
    Index DefaultInteger() const
    {
      DBG_ASSERT(type_ == OT_Integer);
      return (Index)default_number_;
    }
    /** Set the default as an Integer - can be called for
    OT_Integer */
    void SetDefaultInteger(const Index& default_value)
    {
      DBG_ASSERT(type_ == OT_Integer);
      default_number_ = (Number)default_value;
    }
    /** get the default as a string - can be called for OT_String */
    std::string DefaultString() const
    {
      DBG_ASSERT(type_ == OT_String);
      return default_string_;
    }
    /** get the default as a string, but as the index of the string in
     *  the list - helps map from a string to an enum- can be called
     *  for OT_String */
    Index DefaultStringAsEnum() const
    {
      DBG_ASSERT(type_ == OT_String);
      return MapStringSettingToEnum(default_string_);
    }
    /** Set the default as a string - can be called for OT_String */
    void SetDefaultString(const std::string& default_value)
    {
      DBG_ASSERT(type_ == OT_String);
      default_string_ = default_value;
    }
    /** Check if the Number value is a valid setting - can be called
     *  for OT_Number */
    bool IsValidNumberSetting(const Number& value) const
    {
      DBG_ASSERT(type_ == OT_Number);
      if (has_lower_ && ((lower_strict_ == true && value <= lower_) ||
                         (lower_strict_ == false && value < lower_))) {
        return false;
      }
      if (has_upper_ && ((upper_strict_ == true && value >= upper_) ||
                         (upper_strict_ == false && value > upper_))) {
        return false;
      }
      return true;
    }
    /** Check if the Integer value is a valid setting - can be called
     *  for OT_Integer */
    bool IsValidIntegerSetting(const Index& value) const
    {
      DBG_ASSERT(type_ == OT_Integer);
      if (has_lower_ && value < lower_) {
        return false;
      }
      if (has_upper_ && value > upper_) {
        return false;
      }
      return true;
    }
    /** Check if the String value is a valid setting - can be called
     *  for OT_String */
    bool IsValidStringSetting(const std::string& value) const;

    /** Map a user setting (allowing any case) to the case used when
     *  the setting was registered.
     */
    std::string MapStringSetting(const std::string& value) const;

    /** Map a user setting (allowing any case) to the index of the
     *  matched setting in the list of string settings. Helps map a
     *  string setting to an enumeration.
     */
    Index MapStringSettingToEnum(const std::string& value) const;
    //@}

    /** output a description of the option */
    void OutputDescription(const Journalist& jnlst) const;
    /** output a more concise version */
    void OutputShortDescription(const Journalist& jnlst) const;
    /** output a latex version */
    void OutputLatexDescription(const Journalist& jnlst) const;

  private:
    std::string name_;
    std::string short_description_;
    std::string long_description_;
    std::string registering_category_;
    RegisteredOptionType type_;

    bool has_lower_;
    bool lower_strict_;
    Number lower_;
    bool has_upper_;
    bool upper_strict_;
    Number upper_;
    Number default_number_;

    void MakeValidLatexString(std::string source, std::string& dest) const;
    std::string MakeValidLatexNumber(Number value) const;

    /** Compare two strings and return true if they are equal (case
    insensitive comparison) */
    bool string_equal_insensitive(const std::string& s1,
                                  const std::string& s2) const;

    /** class to hold the valid string settings for a string option */
    class string_entry
    {
    public:
      string_entry(const std::string& value, const std::string& description)
          : value_(value), description_(description)
      {}
      std::string value_;
      std::string description_;
    };

    std::vector<string_entry> valid_strings_;
    std::string default_string_;

    /** Has the information as how many-th option this one was
     *  registered. */
    const Index counter_;

    static Index next_counter_;
  };

  /** Class for storing registered options. Used for validation and
   *  documentation.
   */
  class RegisteredOptions : public ReferencedObject
  {
  public:
    /** Constructors / Destructors */
    //@{
    /** Standard Constructor */
    RegisteredOptions()
        :
        current_registering_category_("Uncategorized")
    {}

    /** Standard Destructor */
    ~RegisteredOptions()
    {}
    //@}

    DECLARE_STD_EXCEPTION(OPTION_ALREADY_REGISTERED);

    /** Methods to interact with registered options */
    //@{
    /** set the registering class. All subsequent options will be
     *  added with the registered class */
    void SetRegisteringCategory(const std::string& registering_category)
    {
      current_registering_category_ = registering_category;
    }

    /** retrieve the value of the current registering category */
    std::string RegisteringCategory()
    {
      return current_registering_category_;
    }

    /** Add a Number option (with no restrictions) */
    void AddNumberOption(const std::string& name,
                         const std::string& short_description,
                         Number default_value,
                         const std::string& long_description="");
    /** Add a Number option (with a lower bound) */
    void AddLowerBoundedNumberOption(const std::string& name,
                                     const std::string& short_description,
                                     Number lower, bool strict,
                                     Number default_value,
                                     const std::string& long_description="");
    /** Add a Number option (with a upper bound) */
    void AddUpperBoundedNumberOption(const std::string& name,
                                     const std::string& short_description,
                                     Number upper, bool strict,
                                     Number default_value,
                                     const std::string& long_description="");
    /** Add a Number option (with a both bounds) */
    void AddBoundedNumberOption(const std::string& name,
                                const std::string& short_description,
                                Number lower, bool lower_strict,
                                Number upper, bool upper_strict,
                                Number default_value,
                                const std::string& long_description="");
    /** Add a Integer option (with no restrictions) */
    void AddIntegerOption(const std::string& name,
                          const std::string& short_description,
                          Index default_value,
                          const std::string& long_description="");
    /** Add a Integer option (with a lower bound) */
    void AddLowerBoundedIntegerOption(const std::string& name,
                                      const std::string& short_description,
                                      Index lower, Index default_value,
                                      const std::string& long_description="");
    /** Add a Integer option (with a upper bound) */
    void AddUpperBoundedIntegerOption(const std::string& name,
                                      const std::string& short_description,
                                      Index upper, Index default_value,
                                      const std::string& long_description="");
    /** Add a Integer option (with a both bounds) */
    void AddBoundedIntegerOption(const std::string& name,
                                 const std::string& short_description,
                                 Index lower, Index upper,
                                 Index default_value,
                                 const std::string& long_description="");

    /** Add a String option (with no restrictions) */
    void AddStringOption(const std::string& name,
                         const std::string& short_description,
                         const std::string& default_value,
                         const std::vector<std::string>& settings,
                         const std::vector<std::string>& descriptions,
                         const std::string& long_description="");
    /** Methods that make adding string options with only a few
     *  entries easier */
    void AddStringOption1(const std::string& name,
                          const std::string& short_description,
                          const std::string& default_value,
                          const std::string& setting1,
                          const std::string& description1,
                          const std::string& long_description="");
    void AddStringOption2(const std::string& name,
                          const std::string& short_description,
                          const std::string& default_value,
                          const std::string& setting1,
                          const std::string& description1,
                          const std::string& setting2,
                          const std::string& description2,
                          const std::string& long_description="");
    void AddStringOption3(const std::string& name,
                          const std::string& short_description,
                          const std::string& default_value,
                          const std::string& setting1,
                          const std::string& description1,
                          const std::string& setting2,
                          const std::string& description2,
                          const std::string& setting3,
                          const std::string& description3,
                          const std::string& long_description="");
    void AddStringOption4(const std::string& name,
                          const std::string& short_description,
                          const std::string& default_value,
                          const std::string& setting1,
                          const std::string& description1,
                          const std::string& setting2,
                          const std::string& description2,
                          const std::string& setting3,
                          const std::string& description3,
                          const std::string& setting4,
                          const std::string& description4,
                          const std::string& long_description="");
    void AddStringOption5(const std::string& name,
                          const std::string& short_description,
                          const std::string& default_value,
                          const std::string& setting1,
                          const std::string& description1,
                          const std::string& setting2,
                          const std::string& description2,
                          const std::string& setting3,
                          const std::string& description3,
                          const std::string& setting4,
                          const std::string& description4,
                          const std::string& setting5,
                          const std::string& description5,
                          const std::string& long_description="");
    void AddStringOption6(const std::string& name,
                          const std::string& short_description,
                          const std::string& default_value,
                          const std::string& setting1,
                          const std::string& description1,
                          const std::string& setting2,
                          const std::string& description2,
                          const std::string& setting3,
                          const std::string& description3,
                          const std::string& setting4,
                          const std::string& description4,
                          const std::string& setting5,
                          const std::string& description5,
                          const std::string& setting6,
                          const std::string& description6,
                          const std::string& long_description="");
    void AddStringOption7(const std::string& name,
                          const std::string& short_description,
                          const std::string& default_value,
                          const std::string& setting1,
                          const std::string& description1,
                          const std::string& setting2,
                          const std::string& description2,
                          const std::string& setting3,
                          const std::string& description3,
                          const std::string& setting4,
                          const std::string& description4,
                          const std::string& setting5,
                          const std::string& description5,
                          const std::string& setting6,
                          const std::string& description6,
                          const std::string& setting7,
                          const std::string& description7,
                          const std::string& long_description="");

    /** Get a registered option - this will return NULL if the option
     *  does not exist */
    SmartPtr<const RegisteredOption> GetOption(const std::string& name);

    /** Output documentation for the options - gives a description,
     *  etc. */
    void OutputOptionDocumentation(const Journalist& jnlst, std::list<std::string>& categories);

    /** Output documentation in Latex format to include in a latex file */
    void OutputLatexOptionDocumentation(const Journalist& jnlst, std::list<std::string>& categories);
    //@}

  private:
    std::string current_registering_category_;
    std::map<std::string, SmartPtr<RegisteredOption> > registered_options_;
  };
} // namespace Ipopt

#endif
