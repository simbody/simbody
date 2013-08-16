#ifndef SimTK_SimTKCOMMON_STABLEARRAY_H_
#define SimTK_SimTKCOMMON_STABLEARRAY_H_

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
#include "SimTKcommon/internal/Array.h"

#include <cstddef>
#include <cassert>

namespace SimTK {

/**
 * StableArray<T> is like std::vector<T> (or SimTK::Array_<T>) but more stable 
 * in two ways:
 *  - the addresses of the inserted items never change, even if the array
 *     has to be resized, and
 *  - the index of an inserted item never changes either.
 *
 * The above means that once you insert an item (meaning that a copy of
 * it resides in the StableArray), you can save the address of that copy
 * and/or its index and be certain that adding or deleting other items
 * will leave those unaffected. Once an item has been erased, however, we
 * will feel free to reuse both the heap it once consumed and its index.
 *
 * As your punishment for the crime of enjoying this stability guarantee,
 * consecutive elements of a StableArray are not consecutive in memory.
 *
 * CAUTION: this is not suited for use across binary interfaces because
 * the implementation is fully exposed.
 */
template <class T> class StableArray {
public:
    StableArray() : nOccupiedSlots(0) { }

    // Allocate and fill every slot with the same value.
    explicit StableArray(size_t z, const T& ival=T()) : stuff(z), nOccupiedSlots(z) {
        for (size_t i=0; i<z; ++i) stuff[i] = new T(ival);
    }

    // Copy constructor must preserve slot numbers.
    StableArray(const StableArray& s) : stuff(s.size(),0), nOccupiedSlots(0) {
        for (size_t i=0; i<s.size(); ++i) 
            if (!s.empty(i)) initializeEmptyElement(i, s[i]);
        assert(nItems() == s.nItems());
    }

    // Assignment must preserve slot numbers.
    StableArray& operator=(const StableArray& s) {
        clear();
        stuff.resize(s.size(),0);
        for (size_t i=0; i<s.size(); ++i) 
            if (!s.empty(i)) initializeEmptyElement(i, s[i]);
        assert(nItems() == s.nItems());
        return *this;
    }

    ~StableArray() { clear(); }

    bool empty() const { return stuff.size()==0; }
    bool empty(size_t i) const {
        assert(i < stuff.size());
        return stuff[i] == 0;
    }
    size_t size()   const {return stuff.size();}
    size_t nItems() const {return nOccupiedSlots;}

    // This routine preserves as many items as possible and fills
    // in all empty slots with default values. The array will
    // thus have exactly newz slots with nItems==newz.
    // I'm not sure this is useful for anything.
    void resize(size_t newz, const T& ival=T()) {
        const size_t oldz = stuff.size();
        // If we're throwing anything away, destruct it.
        for (size_t i=newz; i < oldz; ++i)
            eraseElementIfNecessary(i);
        stuff.resize(newz,0);
        // Fill in all empty slots with default values.
        for (size_t i=0; i < newz; ++i)
            initializeElementIfNecessary(i,ival);
        assert(nItems() == newz);
    }

    void clear() {
        for (size_t i=0; i < stuff.size(); ++i)
            eraseElementIfNecessary(i);
        stuff.resize(0);
        assert(nItems()==0);
    }

    // You can push a new item onto the end, or put one in an
    // empty slot.
    void push_back(const T& t) {
        stuff.push_back(new T(t));
        ++nOccupiedSlots;
    }

    // Remove the last element and shrink the list by 1. If the second-to-the-last
    // was empty we'll reduce the length more, so that pop_back() guarantees either
    // that the the last element (back()) is not empty, or the entire list is empty.
    // Don't call this on an empty array.
    void pop_back() {
        assert(!empty() && stuff.back());
        eraseOccupiedElement(stuff.size()-1);   // last element *better* be occupied!

        // Skip over any exposed empties. No need to adjust count.
        do { stuff.pop_back(); } while (!stuff.empty() && !stuff.back());
    }

    void insertAt(size_t i, const T& t) {
        assert(i <= stuff.size());
        if (i == size()) push_back(t);
        else initializeEmptyElement(i,t);
    }

    size_t findFreeSlot() const {
        if (nItems() == size())
            return size();  // no room at the inn!
        for (size_t i=0; i<size(); ++i)
            if (empty(i)) return i;
        assert(false); // where's the empty slot???
        return size_t(-1);
    }

    // Look for the first occupied slot at or after i. Returns
    // size() (that is, one past the end) if none found. Use like this:
    //     for (size_t i=findNextItem(0); i < size(); i=findNextItem(i+1))
    //         do something with item[i] here
    size_t findNextItem(size_t i) {
        assert(i < stuff.size());
        for (; i < stuff.size() && !stuff[i]; ++i);
        return i;
    }

    // Insert the item in the first available slot, or grow the array
    // by one and stick it on the end if there are no free slots. The
    // slot in which it was placed is returned.
    size_t insert(const T& t) {
        const size_t free = findFreeSlot();
        insertAt(free, t);
        return free;
    }


    // Erasing an element slot if it isn't already empty. If we erase the
    // last element we don't have to leave a hole, and in fact we might
    // expose a hole in the second-to-the-last element too. In that 
    // case we erase by resizing away trailing detritus, otherwise we'll
    // make a hole.
    void erase(size_t i) {
        if (i == stuff.size()-1) pop_back();
        else eraseElementIfNecessary(i);
    }

    // This returns the first *occupied* element and blows up if there isn't one.
    const T& front() const {
        const size_t firstItem = findNextItem(0);
        assert(firstItem < stuff.size());
        return *stuff[firstItem];
    }
    T& front() {
        const size_t firstItem = findNextItem(0);
        assert(firstItem < stuff.size());
        return *stuff[firstItem];
    }

    // We don't need to search for back() because the rules here ensure that the
    // last element is not empty.
    const T& back()  const {assert(!empty() && stuff.back()); return *stuff.back();}
    T&       back()        {assert(!empty() && stuff.back()); return *stuff.back();}

    const T& operator[](size_t i) const {
        assert(i < stuff.size() && stuff[i]);
        return *stuff[i];
    }
    T& operator[](size_t i) {
        assert(i < stuff.size() && stuff[i]);
        return *stuff[i];
    }

private:
    size_t              nOccupiedSlots; // not counting empty slots
    Array_<T*,size_t>   stuff;

    // Note that this can leave empty slots at the end of the list which
    // is not a legitimate condition for the StableArray.

    void eraseOccupiedElement(size_t i) {
        assert(i < stuff.size() && stuff[i]);
        delete stuff[i]; stuff[i]=0; --nOccupiedSlots;
    }

    void initializeEmptyElement(size_t i, const T& t) {
        assert(i < stuff.size() && !stuff[i]);
        stuff[i] = new T(t); ++nOccupiedSlots;
    }

    void eraseElementIfNecessary(size_t i) {
        assert(i < stuff.size());
        if (stuff[i]) eraseOccupiedElement(i);
    }

    void initializeElementIfNecessary(size_t i, const T& t) {
        assert(i < stuff.size());
        if (!stuff[i]) initializeEmptyElement(i,t);
    }

};

} // namespace SimTK
  
#endif // SimTK_SimTKCOMMON_STABLEARRAY_H_
