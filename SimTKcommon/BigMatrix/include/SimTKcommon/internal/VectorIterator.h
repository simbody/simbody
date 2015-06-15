#ifndef SimTK_SIMMATRIX_VECTORITERATOR_H_
#define SimTK_SIMMATRIX_VECTORITERATOR_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-15 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

/** @file
Defines an iterator that can be used to iterate over the elements of any
kind of Simbody vector. **/

#include <cstddef>

namespace SimTK {

//==============================================================================
//                            VECTOR ITERATOR
//==============================================================================
/** @brief This is an iterator for iterating over the elements of a Vector_ 
or Vec object.

@tparam ELT             The type of an element stored in the vector whose type
                        is given in \p VECTOR_CLASS.
@tparam VECTOR_CLASS    The type of vector container to iterate through. Its
                        element type must be \p ELT.

This random access iterator can be used with any container that supports
random-access indexing and a <code>size()</code> method. However, the intent is 
for use internally to allow writing a variety of vector math functions without
having to specialize them for the various flavors of vector we support.
**/
template <class ELT, class VECTOR_CLASS>
class VectorIterator {
public:
    typedef ELT         value_type;
    typedef ptrdiff_t   difference_type;
    typedef ELT&        reference;
    typedef ELT*        pointer;
    typedef std::random_access_iterator_tag iterator_category;

    /** Create an iterator for the supplied `vector` and set it to refer to the
    element at `index`. **/
    VectorIterator(VECTOR_CLASS& vector, ptrdiff_t index) 
    :   vectorp(&vector), index(index) {}

    // Default copy constructor, copy assignment, destructor.
    // No default constructor.

    ELT& operator*() {
        assert (index >= 0 && index < vectorp->size());
        return (*vectorp)[(int)index];
    }
    ELT& operator[](ptrdiff_t i) {
        assert (i >= 0 && i < vectorp->size());
        return (*vectorp)[(int)i];
    }
    VectorIterator& operator++() {
        assert (index < vectorp->size());
        ++index;
        return *this;
    }
    VectorIterator operator++(int) {
        assert (index < vectorp->size());
        VectorIterator current = *this;
        ++index;
        return current;
    }
    VectorIterator& operator--() {
        assert (index > 0);
        --index;
        return *this;
    }
    VectorIterator operator--(int) {
        assert (index > 0);
        VectorIterator current = *this;
        --index;
        return current;
    }
    VectorIterator& operator+=(ptrdiff_t n) {
        assert (0 <= index+n && index+n <= vectorp->size());
        index += n;
        return *this;
    }
    VectorIterator& operator-=(ptrdiff_t n) {
        assert (0 <= index-n && index-n <= vectorp->size());
        index -= n;
        return *this;
    }
    bool operator<(const VectorIterator& iter) const {
        return (index < iter.index);
    }
    bool operator>(const VectorIterator& iter) const {
        return (index > iter.index);
    }
    bool operator<=(const VectorIterator& iter) const {
        return (index <= iter.index);
    }
    bool operator>=(const VectorIterator& iter) const {
        return (index >= iter.index);
    }
    ptrdiff_t operator-(const VectorIterator& iter) const {
        return (index - iter.index);
    }
    VectorIterator operator-(ptrdiff_t n) const {
        return VectorIterator(*vectorp, index-n);
    }
    VectorIterator operator+(ptrdiff_t n) const {
        return VectorIterator(*vectorp, index+n);
    }
    bool operator==(const VectorIterator& iter) const {
        return (index == iter.index);
    }
    bool operator!=(const VectorIterator& iter) const {
        return (index != iter.index);
    }
private:
    VECTOR_CLASS*   vectorp;
    ptrdiff_t       index;
};

} //namespace SimTK

#endif // SimTK_SIMMATRIX_VECTORITERATOR_H_
