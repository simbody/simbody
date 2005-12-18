#ifndef SIMTK_SIMBODY_COMMON_H_
#define SIMTK_SIMBODY_COMMON_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/** @file
 * Common include file for all Simbody modules.
 */

#include "simtk/SimTK.h"
#include "simmatrix/SmallMatrix.h"

#include <cassert>
#include <vector>

namespace simtk {

namespace Exception {

class APIMethodFailed : public Base {
public:
    APIMethodFailed(const char* fn, int ln, String method, String cause) : Base(fn,ln)
    {
        setMessage(method + " failed because:\n  " + cause);
    }
};


class FeatureAPIMethodFailed : public Base {
public:
    FeatureAPIMethodFailed(const char* fn, int ln, String fullFeatureName,
        String method, String argInfo, String cause) : Base(fn,ln)
    {
        setMessage("Feature('" + fullFeatureName + "')."
                   + method + "(" + argInfo + ") failed because:\n  " + cause);
    }
};

class PlacementAPIMethodFailed : public Base {
public:
    PlacementAPIMethodFailed(const char* fn, int ln,
        String method, String argInfo, String cause) : Base(fn,ln)
    {
        setMessage("Placement."
                   + method + "(" + argInfo + ") failed because:\n  " + cause);
    }
};

class AccessToInvalidPlacementValue : public Base {
public:
    AccessToInvalidPlacementValue(const char* fn, int ln, String featureName) : Base(fn,ln)
    {
        setMessage("An attempt was made to access an invalid value for Feature '" + featureName 
            + "'s Placement. Features must be realized before their values can be obtained.");
    }
private:
};

// This just reports rep-level bad things up to the API level with a helpful string.
class RepLevelException : public Base {
public:
    RepLevelException(const char* fn, int ln, String message) : Base(fn,ln)
    {
        setMessage(message);
    }
};

class EmptySubsystemPathname : public Base {
public:
    EmptySubsystemPathname(const char* fn, int ln) : Base(fn,ln)
    {
        setMessage("Subsystem pathname was empty.");
    }
private:
};

class IllegalSubsystemPathname : public Base {
public:
    IllegalSubsystemPathname(const char* fn, int ln, String pathname, String badseg) : Base(fn,ln)
    {
        setMessage("Subsystem pathname '" + pathname + "' is illegal at segment '" + badseg + "'.");
    }
private:
};

class SubsystemNameNotFound : public Base {
public:
    SubsystemNameNotFound(const char* fn, int ln, String subname, String parentname) : Base(fn,ln)
    {
        setMessage("Can't find any subsystem named '" + subname + "' in '" + parentname + "'.");
    }
private:
};

class FeatureHasAlreadyBeenPlaced : public Base {
public:
    FeatureHasAlreadyBeenPlaced(const char* fn, int ln, String featureName) : Base(fn,ln)
    {
        setMessage("Can't place Feature " + featureName + " because it already has a Placement.");
    }
private:
};

class FeatureCantBeUsedAsPlacement : public Base {
public:
    FeatureCantBeUsedAsPlacement(const char* fn, int ln, String featureName, String featureTypeName,
        String placementTypeNeeded) : Base(fn,ln)
    {
        setMessage("Can't use " + featureTypeName + " Feature '" + featureName 
                   + "' as a " + placementTypeNeeded + " Placement.");
    }
private:
};

class OnlyFeaturesHavePlacements : public Base {
public:
    OnlyFeaturesHavePlacements(const char* fn, int ln, String subsysName) : Base(fn,ln)
    {
        setMessage("An attempt was made to access a Placement for Subsystem " + subsysName 
                   + " which is not a Feature. Only Features have Placements.");
    }
private:
};

class ExpectedFeatureButGotSubsystem : public Base {
public:
    ExpectedFeatureButGotSubsystem(const char* fn, int ln, String subsysName, int index) : Base(fn,ln)
    {
        setMessage("Child Subsystem " + String(index) + " of Subsystem " + subsysName 
                   + " is not a Feature, but this operation expects a Feature.");
    }
private:
};

class ExpectedFeaturePrototypeButGotSubsystem : public Base {
public:
    ExpectedFeaturePrototypeButGotSubsystem(const char* fn, int ln, 
        String subsysName, String newName, String protoName) : Base(fn,ln)
    {
        setMessage("An attempt was made to add a Feature " + newName + " to Subsystem "
                   + subsysName + " using a prototype named " + protoName 
                   + " but the prototype was not a Feature.");
    }
private:
};

class NoFeatureLevelPlacementForThisKindOfFeature : public Base {
public:
    NoFeatureLevelPlacementForThisKindOfFeature(const char* fn, int ln, 
        String featureName, String featureTypeName) : Base(fn,ln)
    {
        setMessage("Can't use " + featureTypeName + " Feature '" + featureName 
                   + "' as a Placement because this kind of Feature has no"
                     " Feature-level Placement. Did you mean to use one of its Subfeatures?");
    }
private:
};

class IndexOutOfRangeForFeaturePlacementReference : public Base {
public:
    IndexOutOfRangeForFeaturePlacementReference(const char* fn, int ln, 
        String featureName, String featureTypeName, int index) : Base(fn,ln)
    {
        setMessage("Index " + String(index)
                   + " is out of range for a reference to the Placement of "
                   + featureTypeName + " Feature '" + featureName + "'.");
    }
private:
};

class PlacementCantBeUsedForThisFeature : public Base {
public:
    PlacementCantBeUsedForThisFeature(const char* fn, int ln,
                                      String placementTypeName,
                                      String featureName, String featureTypeName) : Base(fn,ln)
    {
        setMessage("Can't use " + placementTypeName + "Placement for "
                   + featureTypeName + "Feature '" + featureName 
                   + "'.");
    }
private:
};

class UnaryOperationNotAllowedForPlacementType : public Base {
public:
    UnaryOperationNotAllowedForPlacementType(const char* fn, int ln, 
                                             String opName, String placementTypeName) 
      : Base(fn,ln)
    {
        setMessage("Unary operator '" + opName + "' can't be used on a "
                   + placementTypeName + "Placement.");
    }
};

class InfixPlacementOperationNotAllowed : public Base {
public:
    InfixPlacementOperationNotAllowed(const char* fn, int ln, 
                                      String leftTypeName, String opName,
                                      String rightTypeName) 
      : Base(fn,ln)
    {
        setMessage("Operation (" + leftTypeName + "Placement " + opName + 
                   " " + rightTypeName + "Placement) is not supported.");
    }
};

class FeatureUsedAsFramePlacementMustBeOnFrame : public Base {
public:
    FeatureUsedAsFramePlacementMustBeOnFrame(const char* fn, int ln, 
        String featureName, String featureTypeName,
        String missingPlacementType) : Base(fn,ln)
    {
        setMessage("Can't use " + featureTypeName + " Feature '" + featureName 
         + "' as a Frame Placement because it doesn't have a parent Frame"
           " from which to inherit the " + missingPlacementType
           + " Placement.");
    }
private:
};

class PlacementMustBeLocal : public Base {
public:
    PlacementMustBeLocal(const char* fn, int ln, String method, 
                         String hostFeature, String offendingFeature) : Base(fn,ln)
    {
        setMessage(method + ": can't add placement expression to Feature '" 
            + hostFeature + "' because it references Feature '"
            + offendingFeature 
            + "' (and possibly others) which is not a descendent of '"
            + hostFeature + "'.");
    }
private:
};

class FeatureAndPlacementOnDifferentTrees : public Base {
public:
    FeatureAndPlacementOnDifferentTrees(const char* fn, int ln, 
                                        String hostFeature, String offendingFeature) : Base(fn,ln)
    {
        setMessage(
            "can't place Feature '" + hostFeature 
            + "' because the supplied placement references Feature '" + offendingFeature
            + "' (and possibly others) and there is no common ancestor.");
    }
private:
};

} // namespace simtk::Exception

/**
 * StableArray<T> is like std::vector<T> but more stable in two ways:
 * (1) the addresses of the inserted items never change, even if the array
 *     has to be resized, and
 * (2) the index of an inserted item never changes either.
 *
 * The above means that once you insert an item (meaning that a copy of
 * it resides in the StableArray), you can save the address of that copy
 * and/or its index and be certain that adding or deleting other items
 * will leave those unaffected. Once an item has been erased, however, we
 * will feel free to reuse both the heap it once consumed and its index.
 *
 * As your punishment for the crime of enjoying this stability guarantee,
 * consecutive elements of a StableArray are not consecutive in memory.
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
    size_t          nOccupiedSlots; // not counting empty slots
    std::vector<T*> stuff;

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


} // namespace simtk

#endif // SIMTK_SIMBODY_COMMON_H_
