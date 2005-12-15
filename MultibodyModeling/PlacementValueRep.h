#ifndef SIMTK_PLACEMENT_VALUE_REP_H_
#define SIMTK_PLACEMENT_VALUE_REP_H_

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

/**@file
 * The opaque, library side implementation of PlacementValues.
 */

#include "SimbodyCommon.h"

#include <string>
#include <cassert>
#include <sstream>

namespace simtk {

class Subsystem;
class Placement;

class PlacementValueRep {
public:
    PlacementValueRep() : valid(false), myHandle(0), owner(0), indexInOwner(-1), client(0) { }
    // warning: default copy & assignment are bitwise leaving bad pointers which must be corrected

    virtual ~PlacementValueRep() { }
    virtual PlacementValueRep* clone() const = 0;
    virtual std::string toString(const std::string& linePrefix) const = 0;

    // Create a copy of this PlacementValue using a new handle and
    // getting rid of the owner.
    void cloneUnownedWithNewHandle(PlacementValue& p) const {
        PlacementValueRep* pr = clone();
        pr->myHandle = &p;
        pr->owner = 0; pr->indexInOwner = -1;
        pr->client = 0;
        p.setRep(pr);
    }

    bool isValid() const { return valid; }
    void setValid(bool v) { valid=v; }

    void                  setMyHandle(PlacementValue& p) {myHandle = &p;}
    bool                  hasHandle()       const {return myHandle != 0;}
    const PlacementValue& getMyHandle()     const {assert(myHandle); return *myHandle;}
    PlacementValue&       updMyHandle()           {assert(myHandle); return *myHandle;} 

    void             setOwner(const Subsystem& s, int index) {owner = &s; indexInOwner=index;}
    bool             hasOwner()        const {return owner != 0;}
    const Subsystem& getOwner()        const {assert(owner);    return *owner;}
    int              getIndexInOwner() const {assert(owner);    return indexInOwner;}

    void             setClientPlacement(const Placement& p) {client = &p;}
    bool             hasClientPlacement()        const {return client != 0;}
    const Placement& getClientPlacement()        const {assert(client); return *client;}

    void checkPlacementValueConsistency(const Subsystem* expOwner, 
                                        int              expIndexInOwner,
                                        const Subsystem& expRoot) const;

private:
    bool valid;                    // Is the stored value (in the concrete Rep) meaningful?

    PlacementValue*  myHandle;     // the PlacementValue whose rep this is

    const Subsystem* owner;        // The Subsystem (if any) which owns this PlacementValue
    int              indexInOwner; // ... and the index in its placementValues list.


    const Placement* client;       // The Placement expression (if any) for which this is the Value.
                                   // That Placement's value pointer must point right
                                   // back here (through the PlacementValue handle).
};

template <class T> class PlacementValueRep_ : public PlacementValueRep {
public:
    PlacementValueRep_<T>() : PlacementValueRep() { }
    explicit PlacementValueRep_<T>(const T& v) : PlacementValueRep(), value(v) { }

    const PlacementValueRep_<T>& getMyHandle() const 
      { return PlacementValueRep_<T>::downcast(PlacementValueRep::getMyHandle()); }

    PlacementValueRep* clone() const { return new PlacementValueRep_<T>(*this); }
    std::string toString(const std::string&) const {
        std::stringstream s;
        s << TypeInfo<T>::name() << "(" << value << ")";   
        return s.str();
    }

    const T& getValue() const     { return value; }
    void     setValue(const T& v) {value=v; setValid(true);}

    SIMTK_DOWNCAST(PlacementValueRep_<T>, PlacementValueRep);
private:
    T value;
};

} // namespace simtk

#endif // SIMTK_PLACEMENT_VALUE_REP_H_
