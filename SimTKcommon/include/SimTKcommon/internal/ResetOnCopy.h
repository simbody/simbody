#ifndef SimTK_SimTKCOMMON_RESET_ON_COPY_H_
#define SimTK_SimTKCOMMON_RESET_ON_COPY_H_

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
//                             RESET ON COPY
//==============================================================================

// Only show the unspecialized helper class in Doxygen.

/** This helper class is used only by ResetOnCopy and is specialized as
necessary to support a variety of template types `T`. **/
template <class T, bool IsScalarType>
class ResetOnCopyHelper {};

/** Ensures that a data member of type `T` is automatically reset to its default
value upon copy construction or copy assignment.\ This allows the
compiler-generated default copy methods to be used.

Here are some usage examples:
@code
class Thing {
    // Able to use default or initializing constructors, but note that any
    // initial values shown here disappear when this is copied; all members will
    // be default-initialized in the copy.

    ResetOnCopy<int>                     m_defint;
    ResetOnCopy<char>                    m_charZ = 'z';
    ResetOnCopy<string>                  m_defstr;
    ResetOnCopy<string>                  m_strHelloC = "hello";
    ResetOnCopy<string>                  m_strGoodbyeS = "goodbye"s;
    ResetOnCopy<short>                   m_shArr[3] = {9,8,7};
    ResetOnCopy<SubsystemIndex>          m_subIx{5};
    ResetOnCopy<std::vector<string>>     m_vstr {"one", "two", "three"};

    // Caution: not a smart pointer; will be cleared but not deleted.
    ResetOnCopy<const int*>              m_ptr;

    // Must use a smart pointer like this for automatic delete on copy.
    ResetOnCopy<unique_ptr<Array_<int>>> m_up{new Array_<int>({1,2,3})};
};
@endcode

Other than copy behavior, an object of type `ResetOnCopy<T>` behaves just like
the underlying object of type `T`. It will implicitly convert to `T` when
needed, and inherit constructors, assignment, and other methods from `T`. Move
construction and move assignment behave as they would for `T`, and an assignment
from an object of type `T` to an object of type `ResetOnCopy<T>` will invoke
`T`'s ordinary copy assignment operator if there is one, and fail to compile if
an attempt is made to use a non-existent assignment operator.

If `T` is a C++ built-in "scalar" type (arithmetic, character, or pointer type)
it will be reset to `T(0)` or `nullptr` when copy constructed or copy assigned.
We don't allow `T` to be an enumerated type here since resetting an enum
to zero isn't always reasonable; use `ReinitOnCopy<T>` instead and provide
an initializing enumerator. For class types `T`, copy construction and copy
assignment will use the default constructor. `ResetOnCopy<T>` adds
`CopyConstructible` and `CopyAssignable` concepts to move-only classes like
`std::unique_ptr`, but those operations just reset the object to its
default-constructed condition.

@tparam  T  Template type that is a numeric, character, or pointer built-in
            type, or a class type that is DefaultConstructible and
            `Destructible`. Enum and array types are not allowed.

@see ReinitOnCopy if you want to reset to a non-default value or use an enum.
**/
template <class T>
class ResetOnCopy : public ResetOnCopyHelper<T, std::is_scalar<T>::value> {
    using Super = ResetOnCopyHelper<T, std::is_scalar<T>::value>;

    /** @cond **/ // These confuse doxygen.
    // TODO: should be possible to specialize this for arrays.
    static_assert(!std::is_array<T>::value,
        "ResetOnCopy<T> does not allow T to be an array. Use an array "
        "of ResetOnCopy<E> instead, where E is the element type.");

    static_assert(!std::is_enum<T>::value,
        "ResetOnCopy<T> does not allow T to be an enum, because resetting to "
        "zero may be inappropriate. Use ReinitOnCopy<T> instead and "
        "provide the desired reinitialization enumerator value.");

    static_assert(std::is_default_constructible<T>::value,
        "ResetOnCopy<T> requires type T to have an accessible default "
        "constructor. Use ReinitOnCopy<T> instead to construct from an "
        "initial value.");

    static_assert(std::is_destructible<T>::value,
        "ResetOnCopy<T> requires type T to have an accessible destructor.");
    /** @endcond **/

public:
    using Super::Super;
    using Super::operator=;

    /** Default constructor performs zero-initialization for built-in types;
    default initialization for class types. **/
    ResetOnCopy() : Super() {} // same as default

    /** Construct or implicitly convert from an object of type `T` if there
    is a suitable copy constructor available. **/
    ResetOnCopy(const T& value) : Super(value) {}

    /** Construct or implicitly convert from an rvalue object of type `T` if
    thereis a suitable move constructor available, otherwise use the copy
    constructor if available, otherwise fails to compile. **/
    ResetOnCopy(T&& value) : Super(std::move(value)) {}

    /** Copy constructor behaves identically to the default constructor; the
    supplied source argument is ignored. **/
    ResetOnCopy(const ResetOnCopy&) : ResetOnCopy() {}

    /** Move constructor is simply a pass-through to the move constructor of
    the contained object so behaves normally. **/
    ResetOnCopy(ResetOnCopy&& source)
    :   Super(static_cast<Super&&>(source)) {} // default

    /** Copy assignment reinitializes this object to its default-constructed
    condition; the source argument is ignored. **/
    ResetOnCopy& operator=(const ResetOnCopy& ignored)
    {   Super::operator=(static_cast<const Super&>(ignored)); return *this; }

    /** Move assignment is simply a pass-through to the move assignment of the
    contained object so behaves normally. **/
    ResetOnCopy& operator=(ResetOnCopy&& source)
    {   Super::operator=(static_cast<Super&&>(source)); return *this; }

    /** Assignment from an object of type `T` uses `T`'s copy assignment
    operator if there is a suitable copy assignment operator available. **/
    ResetOnCopy& operator=(const T& value)
    {   Super::operator=(value); return *this; }

    /** Assignment from an rvalue object of type `T` uses `T`'s move assignment
    operator if available, otherwise uses copy assignment if available,
    otherwise fails to compile. **/
    ResetOnCopy& operator=(T&& value)
    {   Super::operator=(std::move(value)); return *this; }

    /** Return a const reference to the contained object of type `T`. **/
    const T& getT() const {return Super::getT();}

    /** Return a writable reference to the contained object of type `T`. **/
    T& updT() {return Super::updT();}
};

//==============================================================================
//                          RESET ON COPY HELPERS
//==============================================================================


/** @cond **/ // hide helpers from doxygen
/* ResetOnCopy helper class for built-in ("scalar") types (arithmetic,
character, and pointer types). We don't want to allow enums here because they
should always be initialized to one of their values, not necessarily zero.
These types are value initialized, so will be reset to zero. */
template <class T>
class ResetOnCopyHelper<T,true> {
public:
    // These three are just the defaults but for debugging it is helpful to
    // have them explicitly present to step into.

    // Default construction (note that m_value is "value initialized").
    ResetOnCopyHelper() {}

    // Default move construction.
    ResetOnCopyHelper(ResetOnCopyHelper&& source)
    :   m_value(std::move(source.m_value)) {}

    // Default move assignment.
    ResetOnCopyHelper& operator=(ResetOnCopyHelper&& source)
    {   m_value = std::move(source.m_value); return *this; }

    // Copy constructor isn't needed here since ResetOnCopy doesn't use it.
    ResetOnCopyHelper(const ResetOnCopyHelper&) = delete;

    // Constructor from lvalue or rvalue T sets an initial value but this
    // initial value will not be restored on copy.
    explicit ResetOnCopyHelper(const T& value) : m_value(value) {}

    // Make copy assignment produce the same result as if the target were
    // value initialized; the source is ignored.
    ResetOnCopyHelper& operator=(const ResetOnCopyHelper&)
    {   m_value = T{}; return *this; }

    // Allow assignment from an lvalue or rvalue object of type T.
    ResetOnCopyHelper& operator=(const T& value)
    {   m_value = value; return *this; }

    // Implicit conversion to a const reference to the contained object.
    operator const T&() const {return getT();}

    // Implicit conversion to a writable reference to the contained object.
    operator T&() {return updT();}

    const T& getT() const {return m_value;}
    T&       updT()       {return m_value;}

private:
    T   m_value{}; // value initialize; i.e., set to zero
};


/* ResetOnCopy helper class specialization for any type `T` that is not a
built-in ("scalar") type, as long as it is DefaultConstructible and
Destructible. This will add `CopyConstructible` and `CopyAssignable` concepts to
move-only classes like `std::unique`, but those operations just reset the object
to its default-constructed state. */
template <class T>
class ResetOnCopyHelper<T,false> : public T {
public:
    using T::T;
    using T::operator=;

    // These three are just the defaults but for debugging it is helpful to
    // have them explicitly present to step into.

    // Default construction (T is default constructed).
    ResetOnCopyHelper() : T() {}

    // Default move construction.
    ResetOnCopyHelper(ResetOnCopyHelper&& source) : T(std::move(source)) {}

    // Default move assignment.
    ResetOnCopyHelper& operator=(ResetOnCopyHelper&& source)
    {   T::operator=(std::move(source)); return *this; }

    // Copy constructor isn't needed here since ResetOnCopy doesn't use it.
    ResetOnCopyHelper(const ResetOnCopyHelper&) = delete;

    // Constructor from lvalue T sets an initial value but this initial value
    // will not be restored on copy. Won't compile if T doesn't have a copy
    // constructor.
    explicit ResetOnCopyHelper(const T& value) : T(value) {}

    // Constructor from rvalue T sets an initial value but this initial value
    // will not be restored on copy. Won't compile if T doesn't have a move
    // or copy constructor.
    explicit ResetOnCopyHelper(const T&& value) : T(std::move(value)) {}

    // Copy assignment just destroys and reconstructs this object using
    // the default constructor for type `T`.
    ResetOnCopyHelper& operator=(const ResetOnCopyHelper&) {
        T* thisT = static_cast<T*>(this); // upcast to base
        thisT->~T();        // destruct the current base object
        new (thisT) T();    // reconstruct it with default constructor
        return *this;
    }

    const T& getT() const {return static_cast<const T&>(*this);}
    T&       updT()       {return static_cast<T&>(*this);}
};
/** @endcond **/

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_RESET_ON_COPY_H_
