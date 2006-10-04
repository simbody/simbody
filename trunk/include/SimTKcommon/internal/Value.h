#ifndef SimTK_SimTKCOMMON_VALUE_H_
#define SimTK_SimTKCOMMON_VALUE_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
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

#include "SimTKcommon/internal/List.h"
#include "SimTKcommon/internal/String.h"
#include "SimTKcommon/internal/Exception.h"

#include <limits>
#include <typeinfo>
#include <sstream>

namespace SimTK {

/**
 * Abstract base class representing an arbitrary value of self-describing type.
 */
class AbstractValue {
public:
	AbstractValue() { }
	virtual ~AbstractValue() { }

    virtual String      getTypeName() const = 0;    
	virtual String		getValueAsString() const = 0;
    virtual bool        isCompatible(const AbstractValue&) const = 0;
    virtual void compatibleAssign(const AbstractValue&) = 0;
        
    AbstractValue& operator=(const AbstractValue& v) { compatibleAssign(v); return *this; }
	
	virtual AbstractValue* clone() const = 0;
};

inline std::ostream& 
operator<<(std::ostream& o, const AbstractValue& v) { o << v.getValueAsString(); return o; }

/** 
 * Templatized version of the abstract class, providing generic type-specific
 * functionality that does not require specialization.
 */
template <class T> class ValueHelper : public AbstractValue {
public:
    ValueHelper() { } // contained value is default-constructed
    explicit ValueHelper(const T& t) : thing(t) { }
    // default copy, destructor

    // Define assignment explicitly here to avoid implicitly calling AbstractValue's
    // assignment operator.    
    ValueHelper& operator=(const ValueHelper& v) { thing = v.thing; return *this; }
    
    const T& get()      const { return thing; }
  
    // two ways to assign to a new value
    void     set(const T& t)  { thing = t; }    
    T&       upd()            { return thing; }

    bool isCompatible(const AbstractValue& v) const { return isA(v); }        
    void compatibleAssign(const AbstractValue& v) {
        if (!isA(v)) SimTK_THROW2(Exception::IncompatibleValues,v.getTypeName(),getTypeName());
        *this = downcast(v);
    }
    String getTypeName() const { return TypeInfo<T>::name(); }
    String getValueAsString() const 
    { std::ostringstream s; s << thing; return s.str(); }
    
    AbstractValue* clone() const { return new ValueHelper(*this); }
    SimTK_DOWNCAST(ValueHelper,AbstractValue);
protected:
    T thing;
};

/** 
 * A particular kind of AbstractValue, with automatic converstion
 * to the underlying type.
 */
template <class T> class Value : public ValueHelper<T> {
public:
    Value() { }
    explicit Value(const T& t) : ValueHelper<T>(t) { }
    // default copy, assignment, destructor
 
 	// "this->" kludge below appears to be required by gcc 3.4.4.   
    Value& operator=(const T& t) { this->thing=t; return *this; }
    operator const T&() const    { return this->thing; } // automatic conversion to T
    operator T&()                { return this->thing; } // automatic conversion to T
             
    SimTK_DOWNCAST2(Value,ValueHelper<T>,AbstractValue);
private:
    // NO DATA MEMBERS ALLOWED
};


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_VALUE_H_
