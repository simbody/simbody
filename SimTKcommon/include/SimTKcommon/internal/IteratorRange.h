#ifndef SimTK_SimTKCOMMON_ITERATOR_RANGE_H_
#define SimTK_SimTKCOMMON_ITERATOR_RANGE_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-16 Stanford University and the Authors.        *
 * Authors: Chris Dembia                                                      *
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

namespace SimTK {

/** Helper class to use range-based for loops with a pair of iterators. This
 * class should only be used when you're sure the iterators are valid. Don't
 * use this class directly; instead, use makeIteratorRange().
 *
 * @code
 * for (auto& x : makeIteratorRange(v.begin(), v.end())) {
 *     ...
 * }
 * @endcode
 * */
// http://stackoverflow.com/questions/6167598/why-was-pair-range-access-removed-from-c11
template <class Iterator>
class IteratorRange {
public:
    IteratorRange(Iterator first, Iterator last)
        : m_first(first), m_last(last) {}
    Iterator begin() const { return m_first; }
    Iterator end() const { return m_last; }
private:
    const Iterator m_first;
    const Iterator m_last;
};

/** Make an IteratorRange object to be used in a range-based for loop.
 * @relates IteratorRange
 */
template <class Iterator>
IteratorRange<Iterator> makeIteratorRange(Iterator first, Iterator last) {
    return IteratorRange<Iterator>(first, last);
}

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_ITERATOR_RANGE_H_
