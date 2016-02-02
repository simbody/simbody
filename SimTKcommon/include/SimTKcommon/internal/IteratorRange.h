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
 * Portions copyright (c) 2016 Stanford University and the Authors.           *
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
class should only be used when you're sure the iterators are valid. Don't
use this class directly; instead, use makeIteratorRange().

Here's an example of using iterators `first` and `last` to iterate over the
range `[first, last)` (that is, `last` won't be reached):
@code
std::vector<int> v {5, 10, 15, 20, 22};
auto first = std::lower_bound(v.begin(), v.end(), 10);
auto last = std::lower_bound(v.begin(), v.end(), 15); // actually points to 20.
for (auto& x : makeIteratorRange(first, last)) {
    ...
}
@endcode

You can also use this class with an std::pair of iterators, such as that
returned by std::multimap::equal_range(). We assume the first iterator in the
pair is the first iterator in the range, and the second iterator in the pair is
the last iterator in the range.
@code
std::multimap<std::string, int> map;
...
for (auto& x : makeIteratorRange(map.equal_range("some_key"))) {
    ...
}
@endcode
*/
// http://stackoverflow.com/questions/6167598/why-was-pair-range-access-removed-from-c11
template <class Iterator>
class IteratorRange {
public:
    /** This constructor allows you to iterate over the range `[first, last)`;
    this means `last` won't be reached. */
    IteratorRange(Iterator first, Iterator last)
        : m_first(first), m_last(last) {}
    explicit IteratorRange(const std::pair<Iterator, Iterator>& range)
        : m_first(range.first), m_last(range.second) {}
    Iterator begin() const { return m_first; }
    Iterator end() const { return m_last; }
private:
    const Iterator m_first;
    const Iterator m_last;
};

/** Make an IteratorRange object to be used in a range-based for loop, using
two iterators.
@relates IteratorRange
@see IteratorRange::IteratorRange()
*/
template <class Iterator>
IteratorRange<Iterator> makeIteratorRange(Iterator first, Iterator last) {
    return IteratorRange<Iterator>(first, last);
}
/** Make an IteratorRange object to be used in a range-based for loop, using
an std::pair of iterators.
@relates IteratorRange
*/
template <class Iterator>
IteratorRange<Iterator> makeIteratorRange(
        const std::pair<Iterator, Iterator>& range) {
    return IteratorRange<Iterator>(range);
}

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_ITERATOR_RANGE_H_
