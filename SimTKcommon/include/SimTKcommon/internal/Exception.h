#ifndef SimTK_SimTKCOMMON_EXCEPTION_H_
#define SimTK_SimTKCOMMON_EXCEPTION_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

// Keeps MS VC++ 8 quiet about sprintf, strcpy, etc.
#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif

#include "SimTKcommon/internal/common.h"

#include <string>
#include <iostream>
#include <exception>
#include <cstdarg>
#include <cstdio>

namespace SimTK {

namespace Exception {
    
// SimTK::Exception::Base    
class Base : public std::exception {
public:
    explicit Base(const char* fn="<UNKNOWN>", int ln=0) 
      : fileName(fn), lineNo(ln) { } 
    virtual ~Base() throw() { }
    const std::string& getMessage()     const { return msg; }
    const std::string& getMessageText() const { return text; }

    // override virtual function from std::exception
    const char* what() const throw() {return getMessage().c_str();}
protected:
    void setMessage(const std::string& msgin) {
        text = msgin;
        msg = "SimTK Exception thrown at " + where() + ":\n  " + msgin;
    }
private:
    std::string    fileName;    // where the exception was thrown
    int            lineNo;    
    std::string    msg;        // a message formatted for display by catcher
    std::string text;      // the original passed-in text
    
    static std::string shortenFileName(const std::string& fn) 
    {   std::string::size_type pos = fn.find_last_of("/\\");
        if (pos+1>=fn.size()) pos=0;
        return std::string(fn,(int)(pos+1),(int)(fn.size()-(pos+1)));
    }
    
    std::string where() const {
        char buf[32];
        sprintf(buf,"%d",lineNo);
        return shortenFileName(fileName) + ":" + std::string(buf); 
    } 
};

/// This is for reporting internally-detected bugs only, not problems induced by 
/// confused users (that is, it is for confused developers instead). The exception
/// message accepts printf-style arguments and should contain lots of useful
/// information for developers. Don't throw 
/// this exception directly; use one of the family of SimTK_ASSERT and 
/// SimTK_ASSERT_ALWAYS macros.
class Assert : public Base {
public:
    Assert(const char* fn, int ln, const char* assertion, 
             const char* fmt ...) : Base(fn,ln)
    {
        char buf[1024];
        va_list args;
        va_start(args, fmt);
        vsprintf(buf, fmt, args);

        setMessage("Internal bug detected: " + std::string(buf)
                   + "\n  (Assertion '" + std::string(assertion) + "' failed).\n"
            "  Please file an Issue at https://github.com/simbody/simbody/issues.\n"
            "  Include the above information and anything else needed to reproduce the problem.");
        va_end(args);
    }
    virtual ~Assert() throw() { }
};

/// This is for reporting errors occurring during execution of SimTK core methods,
/// beyond those caused by mere improper API arguments, which should be reported with
/// APIArgcheck instead.  Nor is this intended for detection of internal
/// bugs; use Assert instead for that. It is expected that this error resulted from 
/// something the API user did, so the messages should be suitable for reporting to 
/// that programmer. The exception message accepts printf-style arguments and should 
/// contain lots of useful information for the API user. Don't throw this exception 
/// directly; use one of the family SimTK_ERRCHK and SimTK_ERRCHK_ALWAYS macros.
class ErrorCheck : public Base {
public:
    ErrorCheck(const char* fn, int ln, const char* assertion, 
           const char* whereChecked,    // e.g., ClassName::methodName()
           const char* fmt ...) : Base(fn,ln)
    {
        char buf[1024];
        va_list args;
        va_start(args, fmt);
        vsprintf(buf, fmt, args);

        setMessage("Error detected by Simbody method " 
            + std::string(whereChecked) + ": "
            + std::string(buf)
            + "\n  (Required condition '" + std::string(assertion) + "' was not met.)\n");
        va_end(args);
    }
    virtual ~ErrorCheck() throw() { }
};

/// This is for reporting problems detected by checking the caller's supplied arguments
/// to a SimTK API method. Messages should be suitable for SimTK API users. This is not
/// intended for detection of internal bugs where a SimTK developer passed bad arguments
/// to some internal routine -- use Assert instead for that. The exception message 
/// accepts printf-style arguments and should contain useful information for the API user. 
/// Don't throw this exception directly; use one of the family SimTK_APIARGCHECK and 
/// SimTK_APIARGCHECK_ALWAYS macros.
class APIArgcheckFailed : public Base {
public:
    APIArgcheckFailed(const char* fn, int ln, const char* assertion,
                      const char* className, const char* methodName,
                      const char* fmt ...) : Base(fn,ln)
    {
        char buf[1024];
        va_list args;
        va_start(args, fmt);
        vsprintf(buf, fmt, args);
        setMessage("Bad call to Simbody API method " 
                   + std::string(className) + "::" + std::string(methodName) + "(): "
                   + std::string(buf)
                   + "\n  (Required condition '" + std::string(assertion) + "' was not met.)");
        va_end(args);
    }
    virtual ~APIArgcheckFailed() throw() { }
};


class IndexOutOfRange : public Base {
public:
    IndexOutOfRange(const char* fn, int ln, const char* indexName,
                    long long lb, long long index, long long ub, const char* where)
      : Base(fn,ln)
    {
        char buf[1024];

        sprintf(buf, "Index out of range in %s: expected %lld <= %s < %lld but %s=%lld.",
            where,lb,indexName,ub,indexName,index);
        setMessage(std::string(buf));
    }
    virtual ~IndexOutOfRange() throw() { }
};

class SizeOutOfRange : public Base {
public:
    SizeOutOfRange(const char* fn, int ln, const char* szName,
                   unsigned long long sz, unsigned long long maxsz, const char* where)
      : Base(fn,ln)
    {
        char buf[1024];

        sprintf(buf, "Size out of range in %s: expected 0 <= %s <= %llu but %s=%llu.",
            where,szName,maxsz,szName,sz);
        setMessage(std::string(buf));
    }
    virtual ~SizeOutOfRange() throw() { }
};

class SizeWasNegative : public Base {
public:
    SizeWasNegative(const char* fn, int ln, const char* szName,
                   unsigned long long sz, const char* where)
      : Base(fn,ln)
    {
        char buf[1024];

        sprintf(buf, "Size argument was negative in %s: expected 0 <= %s but %s=%llu.",
            where,szName,szName,sz);
        setMessage(std::string(buf));
    }
    virtual ~SizeWasNegative() throw() { }
};

class ValueOutOfRange : public Base {
public:
    ValueOutOfRange(const char* fn, int ln, const char* valueName,
                    double lowerBound, double value, double upperBound, 
                    const char* where)
      : Base(fn,ln)
    {
        char buf[1024];

        sprintf(buf, "Value out of range in %s: expected %g <= %s <= %g but %s=%g.",
            where,lowerBound,valueName,upperBound,valueName,value);
        setMessage(std::string(buf));
    }
    virtual ~ValueOutOfRange() throw() { }
};

class ValueWasNegative : public Base {
public:
    ValueWasNegative(const char* fn, int ln, const char* valueName,
                     double value, const char* where)
      : Base(fn,ln)
    {
        char buf[1024];

        sprintf(buf, "Expected non-negative value for %s in %s but got %g.",
            valueName,where,value);
        setMessage(std::string(buf));
    }
    virtual ~ValueWasNegative() throw() { }
};

class UnimplementedMethod : public Base {
public:
    UnimplementedMethod(const char* fn, int ln, std::string methodName) 
    :   Base(fn,ln)
    { 
        setMessage("The method " + methodName
            + "is not yet implemented. Please post to the Simbody forum"
              " to find a workaround or request implementation.");
    }
    virtual ~UnimplementedMethod() throw() { }
};

class UnimplementedVirtualMethod : public Base {
public:
    UnimplementedVirtualMethod(const char* fn, int ln, 
        std::string baseClass, std::string methodName) 
        : Base(fn,ln)
    { 
        setMessage("The base class " + baseClass + 
            " dummy implementation of method " + methodName
            + "() was invoked because a derived class did not provide an implementation.");
    }
    virtual ~UnimplementedVirtualMethod() throw() { }
};

class IncompatibleValues : public Base {
public:
    IncompatibleValues(const char* fn, int ln, std::string src, std::string dest) : Base(fn,ln)
    {
        setMessage("Attempt to assign a Value<"+src+"> to a Value<"+dest+">");
    }
    virtual ~IncompatibleValues() throw() { }
};

class OperationNotAllowedOnView : public Base {
public:
    OperationNotAllowedOnView(const char* fn, int ln, const std::string& op) : Base(fn,ln)
    {
        setMessage("Operation '" + op + "' allowed only for owners, not views");
    }   
    virtual ~OperationNotAllowedOnView() throw() { }
};

class OperationNotAllowedOnOwner : public Base {
public:
    OperationNotAllowedOnOwner(const char* fn, int ln, const std::string& op) : Base(fn,ln)
    {
        setMessage("Operation '" + op + "' allowed only for views, not owners");
    }   
    virtual ~OperationNotAllowedOnOwner() throw() { }
};

class OperationNotAllowedOnNonconstReadOnlyView : public Base {
public:
    OperationNotAllowedOnNonconstReadOnlyView(const char* fn, int ln, const std::string& op) : Base(fn,ln)
    {
        setMessage("Operation '" + op + "' not allowed on non-const readonly view");
    }   
    virtual ~OperationNotAllowedOnNonconstReadOnlyView() throw() { }
};

// SimTK::Exception::Cant
class Cant : public Base {
public:
    Cant(const char* fn, int ln, const std::string& s) : Base(fn,ln)
    {
        setMessage("Can't perform operation: " + s);
    }    
    virtual ~Cant() throw() { }
};

} // namespace Exception
} // namespace SimTK

#define SimTK_THROW(exc) \
    throw exc(__FILE__, __LINE__)
#define SimTK_THROW1(exc,a1) \
    throw exc(__FILE__, __LINE__,a1)
#define SimTK_THROW2(exc,a1,a2) \
    throw exc(__FILE__, __LINE__,a1,a2)
#define SimTK_THROW3(exc,a1,a2,a3) \
    throw exc(__FILE__, __LINE__,a1,a2,a3)
#define SimTK_THROW4(exc,a1,a2,a3,a4) \
    throw exc(__FILE__, __LINE__,a1,a2,a3,a4)
#define SimTK_THROW5(exc,a1,a2,a3,a4,a5) \
    throw exc(__FILE__, __LINE__,a1,a2,a3,a4,a5)
#define SimTK_THROW6(exc,a1,a2,a3,a4,a5,a6) \
    throw exc(__FILE__, __LINE__,a1,a2,a3,a4,a5,a6)
#define SimTK_THROW7(exc,a1,a2,a3,a4,a5,a6,a7) \
    throw exc(__FILE__, __LINE__,a1,a2,a3,a4,a5,a6,a7)
#define SimTK_THROW8(exc,a1,a2,a3,a4,a5,a6,a7,a8) \
    throw exc(__FILE__, __LINE__,a1,a2,a3,a4,a5,a6,a7,a8)
#define SimTK_THROW9(exc,a1,a2,a3,a4,a5,a6,a7,a8,a9) \
    throw exc(__FILE__, __LINE__,a1,a2,a3,a4,a5,a6,a7,a8,a9)
#define SimTK_THROW10(exc,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10) \
    throw exc(__FILE__, __LINE__,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)

#endif // SimTK_SimTKCOMMON_EXCEPTION_H_

