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
 * Portions copyright (c) 2005-16 Stanford University and the Authors.        *
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
#include "SimTKcommon/internal/Xml.h"

#include <limits>
#include <typeinfo>
#include <type_traits>
#include <sstream>
#include <utility>
#include <unordered_map>
#include <functional>
#include <memory>

namespace SimTK {
// Forward declaration.
template<typename T> class Value;

//==============================================================================
//                             ABSTRACT VALUE
//==============================================================================
/** Abstract base class representing an arbitrary value of unknown type.
This provides an ability to manipulate "values" abstractly without knowing the
specific type of that value. This is useful, for example, for discrete state
variables and Measures, where much of what we need to do with them is 
independent of their types. Support is provided for serializing and 
deserializing %AbstractValue objects without every having to know the underlying
`Value<T>` object that is being dealt with. **/
class SimTK_SimTKCOMMON_EXPORT AbstractValue {
public:
    /** Create a deep copy of this object. **/
    virtual AbstractValue* clone() const = 0;

    /** Return a human-readable form of the object type that is stored in the
    concrete derived class underlying this %AbstractValue. This is the
    output from NiceTypeName's `namestr()` method instantiated for the concrete
    object type. **/
    virtual String getTypeName() const = 0;  


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

    /** Retrieve the original (type-erased) `thing` as read-only. The template
    argument must be exactly the non-reference type of the stored `thing`. **/
    template <typename T>
    const T& getValue() const {
      return Value<T>::downcast(*this).get();
    }

    /** Retrieve the original (type-erased) `thing` as read-write. The template
    argument must be exactly the non-reference type of the stored `thing`. **/
    template <typename T>
    T& updValue() {
      return Value<T>::updDowncast(*this).upd();
    }

    virtual ~AbstractValue() {}   

    /** Serialize the concrete `Value<T>` referenced by this %AbstractValue to
    an Xml::Element. The element is returned as a new orphan element with tag
    `<%Value>` containing a single subelement named `thing`. For simple types 
    the resulting XML looks like this: @verbatim
        <Value name="name" type="typename">
            <thing> <!-- a value of type typename --> </thing>
        </Value>
    @endverbatim while for more elaborate types @verbatim
        <Value name="name" type="typename">
            <TypeTag name="thing" version="1"> 
                <!-- a value of type typename --> 
            </TypeTag>
        </Value>
    @endverbatim where `TypeTag` is the arbitrary tag chosen by the contained
    type's serialization method and the value is arbitrarily complicated.
    The `name` attribute on the `<%Value>` element does not appear if the `name` 
    parameter is empty. A correct `Value<T>` object can then be
    deserialized from the element, provided that *typename* has been 
    registered with the %AbstractValue factory. **/
    virtual Xml::Element toXmlElement(const std::string& name) const = 0;

    /** This is the required signature for a method that can create a 
    concrete `Value<T>` from an XML element that specifies a particular type
    `T`, optionally checking for an expected name. **/
    using ValueFactory = 
        std::function<std::unique_ptr<AbstractValue>
                      (Xml::Element&,const std::string&)>;

    /** This is a factory for creating a type-erased %AbstractValue object from
    a serialized `Value<T>` of arbitrary type `T`. This assumes that the name
    used as the `type` attribute has been registered as the key for a 
    particular type-specific ValueFactory method. If you are expecting a
    particular name for this element you can provide it and an exception will
    be thrown if there is a mismatch. **/
    static std::unique_ptr<AbstractValue>
        createFromXmlElement(Xml::Element& elt,
                             const std::string& expectedName="");

    /** (Advanced) Register a (typeString,factoryMethod) pair to be used for 
    deserializing a `<%Value>` element that contains a value of the type 
    indicated by the key. This method is called automatically by `Value<T>`
    when first constructed with a type `T`. Users should not have to call 
    it. **/
    static const ValueFactory& registerValueFactory
       (const std::string& key, ValueFactory factory);

    /** <b>(Deprecated)</b> Use `toXmlElement().writeToString()` or the
    element's stream insertion `operator<<()` instead. **/
    DEPRECATED_14("use toXmlElement() instead")
    String getValueAsString() const
    {   String s; toXmlElement("").writeToString(s,true); return s; }
private:
    static std::unordered_map<std::string,ValueFactory> m_factory;
};

/** Write a human-readable representation of an AbstractValue to an output
stream, using the toXmlElement() member function to produce an unnamed XML
element with tag `<%Value>` and a `type` attribute providing a representation
of the underlying concrete type.
@relates AbstractValue **/
inline std::ostream& operator<<(std::ostream& o, const AbstractValue& v) 
{   o << v.toXmlElement(""); return o; }


//==============================================================================
//                               VALUE <T>
//==============================================================================
/** Concrete templatized class derived from AbstractValue, adding generic 
value type-specific functionality, with implicit conversion to the underlying 
type `T`. 

@tparam T   The underlying type. Must be an object that can reasonably be
            serialized -- not a function, pointer, reference, or void type.
            Must support default and copy construction and copy assignment.
**/
template <class T> 
class Value : public AbstractValue {

/** @cond **/ // Doxygen has trouble with this
// Don't allow types T that can't reasonably be copied or serialized.
static_assert(   std::is_object<T>::value 
              && !(   std::is_pointer<T>::value 
                   || std::is_member_pointer<T>::value
                   || std::is_null_pointer<T>::value),
    "Value<T>: T can't be a function, pointer, reference, or void type.");

static_assert(std::is_default_constructible<T>::value,
    "Value<T>: T must be default constructible.");
static_assert(std::is_copy_constructible<T>::value,
    "Value<T>: T must be copy constructible.");
static_assert(std::is_copy_assignable<T>::value,
    "Value<T>: T must be copy assignable.");
/** @endcond **/

public:
    /** Creates a `Value<T>` whose contained object of type `T` has been default 
    constructed. **/
    Value() 
    {   registerForDeserialization(); }

    /** Creates a `Value<T>` whose contained object is copy constructed from the
    given `value`. **/
    explicit Value(const T& value) : m_thing(value) 
    {   registerForDeserialization(); }

    /** Creates a `Value<T>` whose contained object is move constructed from the
    given `value`. If type `T` doesn't provide a move constructor then this
    is copy constructed. **/
    explicit Value(T&& value) : m_thing(moveIfMoveConstructible(value)) 
    {   registerForDeserialization(); }

    ~Value() = default;

    /** Copy constructor invokes the type `T` copy constructor to copy the
    contained value. **/
    Value(const Value& value) : m_thing(value.m_thing) 
    {   registerForDeserialization(); }

    /** Move constructor invokes the type `T` move constructor, if available,
    to move the given `value` into this `Value<T>` object, otherwise it uses
    the copy constructor. **/
    Value(Value&& value) : m_thing(moveIfMoveConstructible(value.m_thing)) 
    {   registerForDeserialization(); }

    /** The copy assignment here invokes type `T` copy assignment on the
    contained object, it does not invoke AbstractValue's copy assignment. **/
    Value& operator=(const Value& value) 
    {   m_thing = value.m_thing; return *this; }

    /** The move assignment here invokes type `T` move assignment on the
    contained object, it does not invoke AbstractValue's move assignment. If 
    type `T` doesn't provide move assignment then copy assignment is used
    instead. **/
    Value& operator=(Value&& value) 
    {   m_thing = moveIfMoveAssignable(value.m_thing); return *this; }
 
    /** Assign a new value to the contained object, using the type `T`
    copy assignment operator. **/
    Value& operator=(const T& value) 
    {   m_thing = value; return *this; }

    /** Assign a new value to the contained object, using the type `T`
    move assignment operator, if available, otherwise copy assignment. **/
    Value& operator=(T&& value) 
    {   m_thing = moveIfMoveAssignable(value); return *this; }
  
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
    type `T` move assignment operator, if available, otherwise copy 
    assignment. **/
    void set(T&& value)  {m_thing = moveIfMoveAssignable(value);}   

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
        *this = getDowncast(value);
    }

    /** Use NiceTypeName's `namestr()` method to produce a human-friendly 
    representation of the type `T`. **/
    String getTypeName() const override {return NiceTypeName<T>::namestr();}

    /** Serialize this `Value<T>` object to an XML element with tag `<Value>`, 
    with the given `name` as an attribute if it is not empty. A `type` 
    attribute contains a string uniquely identifying type `T`. Then there
    is a single subelement named "thing" containing a value of type `T`. **/
    Xml::Element toXmlElement(const std::string& name) const override {
        Xml::Element v("Value");
        v.setAttributeValue("type", NiceTypeName<T>::xmlstr());
        if (!name.empty()) v.setAttributeValue("name", name);
        v.appendNode(toXmlElementHelper(m_thing, "thing", true));
        return v;
    }

    /** Deserialize from a given XML element into the `Value<T>` object. The
    expected format is that produced by toXmlElement(). You can optionally
    specify a required name in which case the given element is checked to ensure
    that its name matches. **/
    void fromXmlElement(Xml::Element& e, const std::string& reqName="") {
        SimTK_ERRCHK1_ALWAYS(e.getElementTag()=="Value", 
            "Value<T>::fromXmlElement",
            "Expected tag 'Value' but got '%s'.", e.getElementTag().c_str());
        if (!reqName.empty()) {
            const String& name = e.getRequiredAttributeValue("name");
            SimTK_ERRCHK2_ALWAYS(name==reqName, "Value<T>::fromXmlElement",
                                 "Expected Value named '%s' but got '%s'.", 
                                 reqName.c_str(), name.c_str());
        }
        fromXmlElementHelper(m_thing, *e.element_begin(), "thing", true);
    }

    /** Given an XML element allegedly containing a serialized `Value<T>` where
    T is the identical type to the one instantiated here, create on the heap a 
    new `Value<T>` object whose underlying type T object has the same value
    as the serialized one had. The signature of this method is suitable for
    registration with AbstractValue as a ValueFactory to be used during
    deserialization. **/
    static std::unique_ptr<AbstractValue>
    createFromXmlElement(Xml::Element& e, const std::string& expectedName="") {
        std::unique_ptr<Value> vp(new Value());
        vp->fromXmlElement(e, expectedName);
        return std::unique_ptr<AbstractValue>(vp.release());
    }
    
    /** Return true if the given AbstractValue is an object of this type
    `Value<T>`. **/ 
    static bool isA(const AbstractValue& value)
    {   return dynamic_cast<const Value*>(&value) != nullptr; }

    /** Downcast a const reference to an AbstractValue to a const reference
    to this type `Value<T>`. A std::bad_cast exception is thrown if the type
    is wrong, at least in Debug builds. **/
    static const Value& getDowncast(const AbstractValue& value)
    {   return SimTK_DYNAMIC_CAST_DEBUG<const Value&>(value); }

    /** An acceptable abbreviation for getDowncast(), which returns a `const`
    reference, but not for updDowncast() which returns a writable reference. **/
    static const Value& downcast(const AbstractValue& value)
    {   return getDowncast(value); }

    /** Downcast a writable reference to an AbstractValue to a writable 
    reference to this type `Value<T>`. A std::bad_cast exception is thrown if 
    the type is wrong, at least in Debug builds. **/
    static Value& updDowncast(AbstractValue& value)
    {   return SimTK_DYNAMIC_CAST_DEBUG<Value&>(value); }

    /** <b>(Deprecated)</b> Use `updDowncast()` instead. **/
    DEPRECATED_14("use updDowncast() instead")
    static Value& downcast(AbstractValue& value) {return updDowncast(value);}

private:
    // Note that although this is invoked on every construction, its static
    // data is initialized only once so registration occurs only upon first
    // construction.
    void registerForDeserialization() const {
        static const ValueFactory& factory =
            AbstractValue::registerValueFactory(NiceTypeName<T>::xmlstr(),
                                                createFromXmlElement);
    }

    // Returns its argument as an rvalue reference if T has a move constructor,
    // otherwise as a const lvalue reference to force a copy.
    typename std::conditional<std::is_move_constructible<T>::value, 
                              T&&, const T&>::type
    moveIfMoveConstructible(T& val) {return std::move(val);}

    // Returns its argument as an rvalue reference if T has move assignment,
    // otherwise as a const lvalue reference to force a copy assignment.
    typename std::conditional<std::is_move_assignable<T>::value, 
                              T&&, const T&>::type
    moveIfMoveAssignable(T& val) {return std::move(val);}

private:
    T   m_thing;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_VALUE_H_
