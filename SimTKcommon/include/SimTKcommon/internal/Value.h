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
 * Portions copyright (c) 2005-15 Stanford University and the Authors.        *
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
#include <utility>

namespace SimTK {

//==============================================================================
//                             ABSTRACT VALUE
//==============================================================================
/** Abstract base class representing an arbitrary value of unknown type.
This provides an ability to manipulate "values" abstractly without knowing the
specific type of that value. This is useful, for example, for discrete state
variables and Measures, where much of what we need to do with them is 
independent of their types. **/
class AbstractValue {
public:
    /** Create a deep copy of this object. **/
    virtual AbstractValue* clone() const = 0;

    /** Return a human-readable form of the object type that is stored in the
    concrete derived class underlying this %AbstractValue. **/
    virtual String getTypeName() const = 0;  

    /** Return a human-readable representation of the value stored in this
    %AbstractValue object. **/
    virtual String getValueAsString() const = 0;

    /** Check whether the `other` object contains a value that is assignment-
    compatible with this one. If so you can perform the assignment with the
    assignment operator or with compatibleAssign(). **/
    virtual bool isCompatible(const AbstractValue& other) const = 0;

    /** If the `source` contains a compatible value, assign a copy of
    that value into this object. Otherwise an exception is thrown. **/
    virtual void compatibleAssign(const AbstractValue& source) = 0;
    
    /** Invokes the compatibleAssign() method which will perform the assignment
    if the source object is compatible, or throw an exception otherwise. **/
    AbstractValue& operator=(const AbstractValue& v) 
    {   compatibleAssign(v); return *this; }

    virtual ~AbstractValue() {}   
};

/** Write a human-readable representation of an AbstractValue to an output
stream, using the getValueAsString() member function. 
@relates AbstractValue **/
inline std::ostream& operator<<(std::ostream& o, const AbstractValue& v) 
{   o << v.getValueAsString(); return o; }



//==============================================================================
//                               VALUE <T>
//==============================================================================
/** Concrete templatized class derived from AbstractValue, adding generic 
value type-specific functionality, with implicit conversion to the underlying 
type `T`. **/
template <class T> class Value : public AbstractValue {
public:
    /** Creates a `Value<T>` whose contained object of type `T` has been default 
    constructed. **/
    Value() {}

    /** Creates a `Value<T>` whose contained object is copy constructed from the
    given `value`. **/
    explicit Value(const T& value) : m_thing(value) {}

    /** Creates a `Value<T>` whose contained object is move constructed from the
    given `value`. If type `T` doesn't provide a move constructor then this
    is copy constructed. **/
    explicit Value(const T&& value) : m_thing(std::move(value)) {}

    ~Value() = default;

    /** Copy constructor invokes the type `T` copy constructor to copy the
    contained value. **/
    Value(const Value& value) : m_thing(value.m_thing) {}

    /** Move constructor invokes the type `T` move constructor, if available,
    to move the given `value` into this `Value<T>` object. **/
    Value(const Value&& value) : m_thing(std::move(value.m_thing)) {}

    /** The copy assignment here invokes type `T` copy assignment on the
    contained object, it does not invoke AbstractValue's copy assignment. **/
    Value& operator=(const Value& value) 
    {   m_thing = value.m_thing; return *this; }

    /** The move assignment here invokes type `T` move assignment on the
    contained object, it does not invoke AbstractValue's move assignment. If 
    type `T` doesn't provide move assignment then copy assignment is used
    instead. **/
    Value& operator=(const Value&& value) 
    {   m_thing = std::move(value.m_thing); return *this; }
 
    /** Assign a new value to the contained object, using the type `T`
    copy assignment operator. **/
    Value& operator=(const T& value) 
    {   m_thing = value; return *this; }

    /** Assign a new value to the contained object, using the type `T`
    move assignment operator, if available. **/
    Value& operator=(const T&& value) 
    {   m_thing = std::move(value); return *this; }
  
    /** Return a const reference to the object of type `T` that is contained
    in this `Value<T>` object. There is also an implicit conversion that
    performs the same function. **/
    const T& get() const {return m_thing;}

    /** Return a writable reference to the object of type `T` that is contained
    in this %Value object. There is also an implicit conversion that performs
    the same function. @see set() **/
    T& upd() {return m_thing;}
 
    /** Assign the contained object to the given `value` by invoking the
    type `T` copy assignment operator. @see upd() **/
    void set(const T& value)  {m_thing = value;}   

    /** Assign the contained object to the given `value` by invoking the
    type `T` move assignment operator, if available. **/
    void set(const T&& value)  {m_thing = std::move(value);}   

    /** Implicit conversion from a const reference to a `Value<T>` object to a
    const reference to the type `T` object it contains. This is identical to
    the get() method. **/
    operator const T&() const {return get();}

    /** Implicit conversion from a writable reference to a %Value object to a
    writable reference to the type `T` object it contains. This is identical
    to the upd() method. **/
    operator T&() {return upd();}

    /** Covariant implementation of the AbstractValue::clone() method. 
    (Covariant means that the return type has been changed to the derived 
    type.) **/
    Value* clone() const override {return new Value(*this);}

    /** Test whether a given AbstractValue is assignment-compatible with this
    %Value object. Currently this only returns true if the source is exactly
    the same type as this. **/
    bool isCompatible(const AbstractValue& value) const override 
    {   return isA(value); }

    /** If the given AbstractValue is assignment-compatible, perform the
    assignment. Otherwise, throw an exception (at least in a Debug build). **/
    void compatibleAssign(const AbstractValue& value) override {
        if (!isA(value)) 
            SimTK_THROW2(Exception::IncompatibleValues,value.getTypeName(),
                         getTypeName());
        *this = downcast(value);
    }

    /** Use NiceTypeName to produce a human-friendly representation of the 
    type `T`. **/
    String getTypeName() const override {return NiceTypeName<T>::namestr();}

    /** (Not implemented yet) Produce a human-friendly representation of the 
    contained value of type `T`. Currently just returns the type name. **/
    String getValueAsString() const 
    {   return "Value<" + getTypeName() + ">"; }
    
    /** Return true if the given AbstractValue is an object of this type
    `Value<T>`. **/ 
    static bool isA(const AbstractValue& value)
    {   return dynamic_cast<const Value*>(&value) != nullptr; }

    /** Downcast a const reference to an AbstractValue to a const reference
    to this type `Value<T>`. A std::bad_cast exception is thrown if the type
    is wrong, at least in Debug builds. **/
    static const Value& downcast(const AbstractValue& value)
    {   return SimTK_DYNAMIC_CAST_DEBUG<const Value&>(value); }

    /** Downcast a writable reference to an AbstractValue to a writable 
    reference to this type `Value<T>`. A std::bad_cast exception is thrown if 
    the type is wrong, at least in Debug builds. **/
    static Value& updDowncast(AbstractValue& value)
    {   return SimTK_DYNAMIC_CAST_DEBUG<Value&>(value); }

    /** Deprecated -- use updDowncast() instead. **/
    static Value& downcast(AbstractValue& value) {return updDowncast(value);}

private:
    T   m_thing;
};



} // namespace SimTK

#endif // SimTK_SimTKCOMMON_VALUE_H_
