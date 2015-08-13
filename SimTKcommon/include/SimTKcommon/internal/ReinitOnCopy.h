#ifndef SimTK_SimTKCOMMON_REINIT_ON_COPY_H_
#define SimTK_SimTKCOMMON_REINIT_ON_COPY_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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
#include <cassert>
#include <utility>
#include <memory>
#include <cstddef>
#include <type_traits>

namespace SimTK {

//==============================================================================
//                            REINIT ON COPY
//==============================================================================

// Only show the unspecialized helper class in Doxygen.

/** This helper class is used only by ReinitOnCopy and is specialized as 
necessary to support a variety of template types `T`. **/
template <class T, bool IsScalarType>
class ReinitOnCopyHelper {};

/** Ensures that a data member of type `T` is automatically reinitialized to a 
given initial value upon copy construction or copy assignment.\ This allows the
compiler-generated defaults for these to be used.

The template type `T` here is required to be `CopyAssignable` and 
`CopyConstructible`. There is space overhead here for one extra copy of `T`
used to hold the initial value. The default constructor is suppressed here; if 
you just need the data member reset to zero 
or to its default-constructed value, use class `ResetOnCopy<T>` instead. That 
class has zero overhead and accepts a wider range of template arguments.

Here are some usage examples:
@code
class Thing {
    // Able to use default constructors and assignments

    ReinitOnCopy<int>         m_which{-1};          // reinit to -1 on copy
    ReinitOnCopy<string>      m_name{"unknown"};    // back to "unknown" on copy
    ReinitOnCopy<const char*> m_desc{"none given"}; // similar

    // An example where ResetOnCopy is better.
    ResetOnCopy<unsigned>     m_count1;     // reset to 0 on copy; no overhead
    ReinitOnCopy<unsigned>    m_count2;     // error; no default construction
    ReinitOnCopy<unsigned>    m_count3{0};  // ok, but has space overhead

    // An example where ResetOnCopy is necessary (T doesn't allow copying).
    ResetOnCopy<std::unique_ptr<Something>> m_myThing;
};
@endcode

Other than copy behavior, an object of type `ReinitOnCopy<T>` behaves just like 
the underlying object of type `T`. It will implicitly convert to `T` when 
needed, and inherit constructors, assignment, and other methods from `T`. Move 
construction and move assignment behave as they would for `T`, and an assignment
from an object of type `T` to an object of type `ReinitOnCopy<T>` will invoke 
`T`'s ordinary copy assignment operator if there is one, and fail to compile if
an attempt is made to use a non-existent assignment operator.

@tparam  T  Template type that is a numeric, character, enum, or pointer 
            built-in type, or a class type that is CopyConstructible and 
            CopyAssignable, and has an accessible destructor. Array types are
            not allowed.

@see ResetOnCopy if you only need to reinitialize to the default value.
**/
template <class T> 
class ReinitOnCopy 
:   public ReinitOnCopyHelper<T, 
            std::is_scalar<typename std::remove_all_extents<T>::type>::value> {
    using Super = ReinitOnCopyHelper<T, 
            std::is_scalar<typename std::remove_all_extents<T>::type>::value>;
public:
    using Super::Super;
    using Super::operator=;
   
    /** Default constructor is deleted; use ResetOnCopy instead. **/
    ReinitOnCopy() = delete;

    /** Construct or implicitly convert from an object of type `T`.\ This sets
    both the current and remembered initial value to the given `value`. **/ 
    ReinitOnCopy(const T& value) : Super(value) {}

    /** Construct or implicitly convert from an rvalue object of type 
    `T`.\ This sets both the current and remembered initial value to the given
    `value`. **/ 
    ReinitOnCopy(T&& value) : Super(std::move(value)) {}

    /** Copy constructor sets the value and remembered initial value to the
    initial value in the source, using type `T`'s copy constructor. The 
    current value of the source is ignored. **/
    ReinitOnCopy(const ReinitOnCopy& source) : Super(source) {}

    /** Move constructor is simply a pass-through to the move constructor of
    the contained object for both the current and initial values. **/
    ReinitOnCopy(ReinitOnCopy&& source) : Super(std::move(source)) {} // default

    /** Copy assignment reinitializes this object to its original condition; the
    source argument is ignored. **/
    ReinitOnCopy& operator=(const ReinitOnCopy& ignored) 
    {   Super::operator=(ignored); return *this; }

    /** Move assignment uses type `T`'s move assignment for the current value
    but does not change the remembered initial value here. **/
    ReinitOnCopy& operator=(ReinitOnCopy&& source) 
    {   Super::operator=(std::move(source)); return *this; }

    /** Assignment from an object of type `T` uses `T`'s copy assignment
    operator; affects only the current value but does not change the remembered
    initial value. **/
    ReinitOnCopy& operator=(const T& value) 
    {   Super::operator=(value); return *this; }

    /** Assignment from an rvalue object of type `T` uses `T`'s move or copy
    assignment operator; affects only the current value but does not change the 
    remembered initial value. **/
    ReinitOnCopy& operator=(T&& value) 
    {   Super::operator=(std::move(value)); return *this; }

    /** Return a const reference to the contained object of type `T`. **/
    const T& getT() const {return Super::getT();}

    /** Return a writable reference to the contained object of type `T`. **/
    T& updT() {return Super::updT();}

    /** (Advanced) Return a const reference to the stored initial value. **/
    const T& getReinitT() const {return Super::getReinitT();}

    /** (Advanced) Return a writable reference to the stored initial 
    value.\ Use of this should be rare and restricted to constructors. **/
    T& updReinitT() {return Super::updReinitT();}
};

//==============================================================================
//                          REINIT ON COPY HELPERS
//==============================================================================


/** @cond **/ // hide helpers from doxygen
/** ReinitOnCopy helper class for built-in types (integral or floating point). 
These types are value initialized, so will be reset to zero. **/
template <class T>
class ReinitOnCopyHelper<T,true> {

    // TODO: should be possible to specialize this for arrays.
    static_assert(!std::is_array<T>::value, 
        "ReinitOnCopy<T> does not allow T to be an array. Use an array "
        "of ReinitOnCopy<E> instead, where E is the element type.");
public:
    // These three are just the defaults but for debugging it is helpful to
    // have them explicitly present to step into.

    // Default construction (note that members are "value initialized").
    ReinitOnCopyHelper() {}

    // Copy constructor sets the value and remembered initial value to the
    // initial value in the source, using type `T`'s copy constructor. The 
    // current value of the source is ignored.
    ReinitOnCopyHelper(const ReinitOnCopyHelper& source)
    :   ReinitOnCopyHelper(source.m_reinitValue) {}

    // Default move construction.
    ReinitOnCopyHelper(ReinitOnCopyHelper&& source) 
    :   m_value(std::move(source.m_value)), 
        m_reinitValue(std::move(source.m_reinitValue)) {}

    // Constructor from lvalue or rvalue `T` sets the value and remembered 
    // initial value to the given value, using type `T`'s copy constructor.
    explicit ReinitOnCopyHelper(const T& value) 
    :   m_value(value), m_reinitValue(value) {}

    // Move assignment moves the *value* from source to `this` but does not
    // move the recorded initial value which may differ.
    ReinitOnCopyHelper& operator=(ReinitOnCopyHelper&& source) 
    {   m_value = std::move(source.m_value); return *this; }

    // Copy assignment resets the current value to the remembered initial value
    // using type `T`'s copy assignment operator. The source is ignored.
    ReinitOnCopyHelper& operator=(const ReinitOnCopyHelper&) 
    {   m_value = m_reinitValue; return *this; }

    // Allow assignment from an lvalue or rvalue object of type T; affects only
    // the current value. Uses built-in type `T`'s copy assignment operator.
    ReinitOnCopyHelper& operator=(const T& value) 
    {   m_value = value; return *this; }

    // Implicit conversion to a const reference to the contained object of
    // type `T`.
    operator const T&() const {return getT();}

    // Implicit conversion to a writable reference to the contained object of
    // type `T`.
    operator T&() {return updT();}

    const T& getT() const {return m_value;}
    T&       updT()       {return m_value;}

    // Return a const reference to the stored initial value.
    const T& getReinitT() const {return m_reinitValue;}

    // Return a writable reference to the stored initial value. Use of this 
    // should be restricted to constructors.
    T& updReinitT() {return m_reinitValue;}

private:
    T           m_value{};          // note "value initialization"; i.e. zero
    const T     m_reinitValue{}; 
};


/** ReinitOnCopy helper class specialization for any type `T` that is not a
built-in ("scalar") type and that is `CopyConstructible` and `CopyAssignable`.
Those operators are used to reinitialize the object to a stored initial value
when copy constructor or copy assignment is performed. **/
template <class T>
class ReinitOnCopyHelper<T,false> : public T {

    // TODO: should be possible to specialize this for arrays.
    static_assert(!std::is_array<T>::value, 
        "ReinitOnCopy<T> does not allow T to be an array. Use an array "
        "of ReinitOnCopy<E> instead, where E is the element type.");

    static_assert(   std::is_copy_constructible<T>::value
                  && std::is_copy_assignable<T>::value,
        "ReinitOnCopy<T> requires type T to have a copy constructor and copy "
        "assignment operator. Use ResetOnCopy<T> instead to reinitialize to "
        "the default value.");

public:
    using T::T;
    using T::operator=;

    // Default construction (T and the stored initial value are default 
    // constructed).
    ReinitOnCopyHelper() : T(), m_reinitValue() {}

    // Move construction moves both the value and initial value from the 
    // source object. This is the same as default move construction.
    ReinitOnCopyHelper(ReinitOnCopyHelper&& source) 
    :   T(std::move(source)), m_reinitValue(std::move(source.m_reinitValue)) {}

    // Constructor from lvalue `T` sets the value and remembered initial value 
    // to the given value, using type `T`'s copy constructor.
    explicit ReinitOnCopyHelper(const T& value) 
    :   T(value), m_reinitValue(value) {}

    // Constructor from rvalue `T` sets the value and remembered initial value 
    // to the given value, using type `T`'s move and copy constructors.
    explicit ReinitOnCopyHelper(T&& value) 
    :   T(value), m_reinitValue(static_cast<const T&>(*this)) {} // careful!

    // Copy constructor sets the value and remembered initial value to the
    // initial value in the source, using type `T`'s copy constructor. The 
    // current value of the source is ignored.
    ReinitOnCopyHelper(const ReinitOnCopyHelper& source)
    :   ReinitOnCopyHelper(source.m_reinitValue) {}

    // Move assignment moves the *value* from source to `this` but does not
    // move the recorded initial value which may differ.
    ReinitOnCopyHelper& operator=(ReinitOnCopyHelper&& source) 
    {   T::operator=(std::move(source)); return *this; }

    // Copy assignment resets the current value to the remembered initial value
    // using type `T`'s copy assignment operator. The source is ignored.
    ReinitOnCopyHelper& operator=(const ReinitOnCopyHelper&) 
    {   T::operator=(m_reinitValue); return *this; }

    // Allow assignment from an lvalue object of type T; affects only the 
    // current value. Uses type `T`'s copy assignment operator.
    ReinitOnCopyHelper& operator=(const T& value) 
    {   T::operator=(value); return *this; }

    // Allow assignment from an rvalue object of type T; affects only the 
    // current value. Uses type `T`'s move assignment operator.
    ReinitOnCopyHelper& operator=(T&& value) 
    {   T::operator=(std::move(value)); return *this; }

    const T& getT() const {return static_cast<const T&>(*this);}
    T&       updT()       {return static_cast<T&>(*this);}

    // Return a const reference to the stored initial value.
    const T& getReinitT() const {return m_reinitValue;}

    // Return a writable reference to the stored initial value. Use of this 
    // should be restricted to constructors.
    T& updReinitT() {return m_reinitValue;}

private:
    // The remembered initial value is set to the already-constructed
    // base value in case we used one of the base constructors. This uses
    // the base class copy constructor.
    const T m_reinitValue {static_cast<const T&>(*this)}; 
};

/** @endcond **/

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_REINIT_ON_COPY_H_

