#ifndef SimTK_SimTKCOMMON_NULL_ON_COPY_UNIQUE_PTR_H_
#define SimTK_SimTKCOMMON_NULL_ON_COPY_UNIQUE_PTR_H_

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

namespace SimTK {

/** This is a smart pointer representing unique ownership like `std::unique_ptr`
but adds "reset" copy construction and copy assignment. Judicious use of this 
container will allow you to use compiler-generated copy constructors and copy 
assignment operators for classes which would otherwise have to implement their 
own in order to properly reset these pointer data members, which must not 
be copied.

The contained pointer is initialized to `nullptr` on construction, is properly
deleted on destruction, and the target is reset to to `nullptr` upon copy 
construction or copy assignment. Generally this is used when the owned object
has pointers into its owner, so should not be copied along with the owner 
since those pointers would point to the source object rather than the target.

The pointer *is* moved intact for move construction or move assignment. That
allows `std::vector<NullOnCopyUniquePtr<T>>` or 
`SimTK::Array_<NullOnCopyUniquePtr<T>>` to behave properly when their contents 
have to be moved for expansion (if copied they would end up null!).

This class derives publicly from `std::unique_ptr` and inherits its constructors
and other methods, so behaves identically except where noted.

@see ReferencePtr, ClonePtr, CloneOnWritePtr **/ 
template < class T, class D = std::default_delete<T> >
class NullOnCopyUniquePtr : public std::unique_ptr<T,D> {
    using Super = std::unique_ptr<T,D>;
public:
    //using Super::unique_ptr; // inherit most constructors (wait for VS2015)
                               // TODO: also fix has below for vs2015

    //--------------------------------------------------------------------------
    // Inheriting constructors doesn't work in VS2013. These can be removed
    // when we switch to VS2015.
    explicit NullOnCopyUniquePtr(typename Super::pointer p) NOEXCEPT_11 
    :   Super(p) {}
    explicit NullOnCopyUniquePtr(std::nullptr_t) NOEXCEPT_11 : Super() {}
    explicit NullOnCopyUniquePtr(Super&& up) NOEXCEPT_11 
    :   Super(std::move(up)) {}
    //--------------------------------------------------------------------------

    /** Default constructor creates a null pointer. **/
    NullOnCopyUniquePtr() NOEXCEPT_11 : Super() {}

    /** Copy constructor ignores source and creates a null pointer, 
    same as default construction. See class description for why. **/
    NullOnCopyUniquePtr(const NullOnCopyUniquePtr&) NOEXCEPT_11
    :   NullOnCopyUniquePtr() {}

    /** Move constructor transfers ownership (same as std::unique_ptr). **/
    NullOnCopyUniquePtr(NullOnCopyUniquePtr&& source) NOEXCEPT_11
    :   Super(std::move(source)) {}

    /** Copy assignment ignores source, leaves target null. See class
    description for why. **/
    NullOnCopyUniquePtr& operator=(const NullOnCopyUniquePtr&) NOEXCEPT_11 {
        Super::reset();
        return *this;
    }

    /** Move assignment transfers ownership (same as std::unique_ptr). **/
    NullOnCopyUniquePtr& operator=(NullOnCopyUniquePtr&& src) NOEXCEPT_11 {
        Super::operator=(std::move(src));
        return *this;
    }

    /** Allow move assignment from identical `std::unique_ptr`. **/
    NullOnCopyUniquePtr& operator=(Super&& up) NOEXCEPT_11 {
        Super::operator=(std::move(up));
        return *this;
    }

}; 



//==============================================================================
//                       SimTK namespace-scope functions
//==============================================================================
// These namespace-scope functions will be resolved by the compiler using
// "Koenig lookup" which examines the arguments' namespaces first.
// See Herb Sutter's discussion here: http://www.gotw.ca/publications/mill08.htm.

/** This is an overload of the STL std::swap() algorithm which uses the
cheap built-in swap() member of the NullOnCopyUniquePtr class. (This function
is defined in the `SimTK` namespace.)
@relates SimTK::NullOnCopyUniquePtr **/
template <class T, class D> inline void
swap(NullOnCopyUniquePtr<T,D>& p1, NullOnCopyUniquePtr<T,D>& p2) NOEXCEPT_11 {
    p1.swap(p2);
}

/** Output the system-dependent representation of the pointer contained
in a NullOnCopyUniquePtr object. This is equivalent to `os << p.get();`.
@relates NullOnCopyUniquePtr **/
template <class charT, class traits, class T, class D>
inline std::basic_ostream<charT,traits>& 
operator<<(std::basic_ostream<charT,traits>& os, 
           const NullOnCopyUniquePtr<T,D>&  p) 
{   os << p.get(); return os; }

} // namespace SimTK



//==============================================================================
//                             hash specialization
//==============================================================================
namespace std {
/** This is a specialization of std::hash<std::unique_ptr<T,D>> so that
SimTK::NullOnCopyUniquePtr hashes identically. This specialization is in the 
`std` namespace. 
@see SimTK::NullOnCopyUniquePtr **/
template <class T, class D>
struct hash<SimTK::NullOnCopyUniquePtr<T,D>> 
:   public hash<std::unique_ptr<T,D>> 
{
    //using hash<std::unique_ptr<T,D>>::hash; // inherit constructors (vs2015)
};
} // namespace std


#endif // SimTK_SimTKCOMMON_NULL_ON_COPY_UNIQUE_PTR_H_
