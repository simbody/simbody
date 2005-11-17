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

class SubfeatureNameNotFound : public Base {
public:
    SubfeatureNameNotFound(const char* fn, int ln, String subname, String featurename) : Base(fn,ln)
    {
        setMessage("Can't find any subfeature named '" + subname + "' in '" + featurename + "'.");
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
 * StableArray<T> is like std::vector<T> but the addresses of the inserted items
 * never change, even if the array has to be resized. As a punishment
 * for this guarantee, consecutive elements of a StableArray are not consecutive
 * in memory.
 */
template <class T> class StableArray {
public:
    StableArray() { }
    explicit StableArray(size_t z, const T& ival=T()) {
        resize(z, ival);
    }
    StableArray(const StableArray& s) {
        resize(s.size());
        for (size_t i=0; i<s.size(); ++i)
            *stuff[i] = *s.stuff[i];
    }
    StableArray& operator=(const StableArray& s) {
        clear();
        resize(s.size());
        for (size_t i=0; i<s.size(); ++i)
            *stuff[i] = *s.stuff[i];
        return *this;
    }
    ~StableArray() { clear(); }

    bool   empty() const { return stuff.size()==0; }
    size_t size()  const { return stuff.size(); }
    void resize(size_t newz, const T& ival=T()) {
        const size_t oldz = stuff.size();
        // If we're throwing anything away, destruct it.
        for (size_t i=newz; i < oldz; ++i)
            {delete stuff[i]; stuff[i]=0;}
        stuff.resize(newz);
        // If we're adding anything new, initialize it.
        for (size_t i=oldz; i < newz; ++i)
            stuff[i] = new T(ival); 
    }
    void clear()               {resize(0);}
    void push_back(const T& t) {stuff.push_back(new T(t));}
    void pop_back()            {assert(!empty()); resize(size()-1);}

    const T& front() const {assert(!empty() && stuff[0]);        return *stuff[0];}
    T&       front()       {assert(!empty() && stuff[0]);        return *stuff[0];}
    const T& back()  const {assert(!empty() && stuff[size()-1]); return *stuff[size()-1];}
    T&       back()        {assert(!empty() && stuff[size()-1]); return *stuff[size()-1];}

    const T& operator[](size_t i) const {
        assert(i < stuff.size() && stuff[i]);
        return *stuff[i];
    }
    T& operator[](size_t i) {
        assert(i < stuff.size() && stuff[i]);
        return *stuff[i];
    }

private:
    std::vector<T*> stuff;
};


} // namespace simtk

#endif // SIMTK_SIMBODY_COMMON_H_
