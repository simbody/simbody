#ifndef SimTK_SimTKCOMMON_EXCEPTION_H_
#define SimTK_SimTKCOMMON_EXCEPTION_H_

/* Copyright (c) 2005 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/String.h"

#include <string>
#include <iostream>
#include <exception>
#include <cstdarg>

namespace SimTK {

namespace Exception {
	
// SimTK::Exception::Base	
class Base : public std::exception {
public:
	explicit Base(const char* fn="<UNKNOWN>", int ln=0) 
      : fileName(fn), lineNo(ln) { } 
	virtual ~Base() throw() { }
	const String& getMessage()     const { return msg; }
    const String& getMessageText() const { return text; }

    // override virtual function from std::exception
    const char* what() const throw() {return getMessage().c_str();}
protected:
	void setMessage(const String& msgin) {
        text = msgin;
        msg = "SimTK Exception thrown at " + where() + ":\n  " + msgin;
    }
private:
	String	fileName;	// where the exception was thrown
	int		lineNo;	
	String	msg;		// a message formatted for display by catcher
    String  text;      // the original passed-in text
    
    static String shortenFileName(const String& fn) 
    {   String::size_type pos = fn.find_last_of("/\\");
        if (pos+1>=fn.size()) pos=0;
        return String(fn,pos+1,fn.size()-(pos+1));
    }
	
	String where() const {
        char buf[32];
        sprintf(buf,"%d",lineNo);
        return shortenFileName(fileName) + ":" + String(buf); 
    } 
};

// This is for reporting internally-detected bugs only, not
// problems induced by confused users.
class Assert : public Base {
public:
    Assert(const char* fn, int ln, const char* assertion, 
           const char* fmt ...) : Base(fn,ln)
    {
        char buf[1024];
        va_list args;
        va_start(args, fmt);
        vsprintf(buf, fmt, args);

        setMessage("Internal SimTK bug detected: " + String(buf)
                   + " (Assertion '" + String(assertion) + "' failed). "
                     "Please report this to the appropriate authorities at SimTK.org.");
        va_end(args);
    }
private:
};


class APIArgcheckFailed : public Base {
public:
    APIArgcheckFailed(const char* fn, int ln, const char* className, const char* methodName,
                      const char* fmt ...) : Base(fn,ln)
    {
        char buf[1024];
        va_list args;
        va_start(args, fmt);
        vsprintf(buf, fmt, args);
        setMessage("Bad call to SimTK API method " 
                   + String(className) + "::" + String(methodName) + "(): "
                   + String(buf) + ".");
        va_end(args);
    }
private:
};


class IndexOutOfRange : public Base {
public:
    IndexOutOfRange(const char* fn, int ln, const char* indexName,
                    long lb, long index, long ub, const char* where)
      : Base(fn,ln)
    {
        char buf[1024];

        sprintf(buf, "Index out of range in %s: expected %ld <= %s < %ld but %s=%ld.",
            where,lb,indexName,ub,indexName,index);
        setMessage(String(buf));
    }
private:
};

class SizeOutOfRange : public Base {
public:
    SizeOutOfRange(const char* fn, int ln, const char* szName,
                   long sz, long maxsz, const char* where)
      : Base(fn,ln)
    {
        char buf[1024];

        sprintf(buf, "Size out of range in %s: expected 0 <= %s <= %ld but %s=%ld.",
            where,szName,maxsz,szName,sz);
        setMessage(String(buf));
    }
private:
};

class SizeWasNegative : public Base {
public:
    SizeWasNegative(const char* fn, int ln, const char* szName,
                   long sz, const char* where)
      : Base(fn,ln)
    {
        char buf[1024];

        sprintf(buf, "Size argument was negative in %s: expected 0 <= %s but %s=%ld.",
            where,szName,szName,sz);
        setMessage(String(buf));
    }
private:
};

class ValueOutOfRange : public Base {
public:
    ValueOutOfRange(const char* fn, int ln, const char* valueName,
                    Real lowerBound, Real value, Real upperBound, 
                    const char* where)
      : Base(fn,ln)
    {
        char buf[1024];

        sprintf(buf, "Value out of range in %s: expected %g <= %s <= %g but %s=%g.",
            where,lowerBound,valueName,upperBound,valueName,value);
        setMessage(String(buf));
    }
private:
};

class ValueWasNegative : public Base {
public:
    ValueWasNegative(const char* fn, int ln, const char* valueName,
                     Real value, const char* where)
      : Base(fn,ln)
    {
        char buf[1024];

        sprintf(buf, "Expected non-negative value for %s in %s but got %g.",
            valueName,where,value);
        setMessage(String(buf));
    }
private:
};

class UnimplementedVirtualMethod : public Base {
public:
    UnimplementedVirtualMethod(const char* fn, int ln, 
        String baseClass, String methodName) 
		: Base(fn,ln)
	{ 
		setMessage("The base class " + baseClass + 
            " dummy implementation of method " + methodName
            + "() was invoked because a derived class did not provide an implementation.");
	}
};

class IncompatibleValues : public Base {
public:
    IncompatibleValues(const char* fn, int ln, String src, String dest) : Base(fn,ln)
    {
        setMessage("Attempt to assign a Value<"+src+"> to a Value<"+dest+">");
    }
private:
};

class OperationNotAllowedOnView : public Base {
public:
    OperationNotAllowedOnView(const char* fn, int ln, const String& op) : Base(fn,ln)
    {
        setMessage("Operation '" + op + "' allowed only for owners, not views");
    }   
};

class OperationNotAllowedOnOwner : public Base {
public:
    OperationNotAllowedOnOwner(const char* fn, int ln, const String& op) : Base(fn,ln)
    {
        setMessage("Operation '" + op + "' allowed only for views, not owners");
    }   
};

class OperationNotAllowedOnNonconstReadOnlyView : public Base {
public:
    OperationNotAllowedOnNonconstReadOnlyView(const char* fn, int ln, const String& op) : Base(fn,ln)
    {
        setMessage("Operation '" + op + "' not allowed on non-const readonly view");
    }   
};

// SimTK::Exception::Cant
class Cant : public Base {
public:
	Cant(const char* fn, int ln, const String& s) : Base(fn,ln)
	{
		setMessage("Can't perform operation: " + s);
	}	
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

