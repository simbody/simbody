#ifndef SimTK_SimTKCOMMON_VALUE_H_
#define SimTK_SimTKCOMMON_VALUE_H_

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

#include "SimTKcommon/internal/common.h"

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
 * functionality that does not require specialization, with automatic conversion
 * to the underlying type.
 *
 * Note that this class requires that type T supports an output operator "<<"
 * so that we can serialize abstract values.
 */
template <class T> class Value : public AbstractValue {
public:
    Value() { } // contained value is default-constructed
    explicit Value(const T& t) : thing(t) { }
    // default copy, destructor

    // Define assignment explicitly here to avoid implicitly calling AbstractValue's
    // assignment operator.    
    Value& operator=(const Value& v) { thing = v.thing; return *this; }
 
    Value& operator=(const T& t) { thing=t; return *this; }
    operator const T&() const    { return thing; } // automatic conversion to T
    operator T&()                { return thing; } // automatic conversion to T
    
    const T& get()      const { return thing; }
  
    // two ways to assign to a new value
    void     set(const T& t)  { thing = t; }    
    T&       upd()            { return thing; }

    bool isCompatible(const AbstractValue& v) const { return isA(v); }        
    void compatibleAssign(const AbstractValue& v) {
        if (!isA(v)) SimTK_THROW2(Exception::IncompatibleValues,v.getTypeName(),getTypeName());
        *this = downcast(v);
    }
    String getTypeName() const { return NiceTypeName<T>::name(); }
    // TODO: should have some general way to serialize these.
    String getValueAsString() const 
    { return "Value<" + getTypeName() + ">"; }
    
    AbstractValue* clone() const { return new Value(*this); }
    SimTK_DOWNCAST(Value,AbstractValue);
protected:
    T thing;
};



} // namespace SimTK

#endif // SimTK_SimTKCOMMON_VALUE_H_
