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
#include "PlacementValue.h"

#include <string>
#include <cassert>
#include <sstream>

namespace simtk {

class Subsystem;
class Placement;

/**
 * This class holds a single PlacementValue within a Subsystem and handles the required
 * administrative overhead for it.
 */
class PlacementValueSlot {
public:
    PlacementValueSlot() 
        : value(), valid(false), owner(0), indexInOwner(-1), client(0) { }
    explicit PlacementValueSlot(const PlacementValue& v)
        : value(v), valid(false), owner(0), indexInOwner(-1), client(0) { }

    // Warning: we are intentionally using the default copy & assignment
    // which are *bitwise* leaving bad pointers which must be corrected afterwards.
    // Default destructor is OK.

    /// Obtain the PlacementValue stored in this slot. If it is not currently valid
    /// we'll throw an exception here.
    const PlacementValue& getValue() const {
        if (!isValid()) {
            throwInvalid(); 
            //NOTREACHED
        }
        return value;
    }

    /// Assign a new Value to this slot and mark it valid. The PlacementValue must
    /// be assignment compatible with the one that is already here.
    void setValue(const PlacementValue& v) {
        value = v;
        setValid(true);
    }

    /// If you want direct access to the PlacementValue to modify it, we'll mark
    /// it invalid here and you'll have to mark it valid yourself.
    PlacementValue& updValue() {
        setValid(false);
        return value;
    }

    bool isValid()        const {return valid;}
    void setValid(bool v)       {valid=v;}

    void                 setOwner(const Subsystem& s, int index) {owner = &s; indexInOwner=index;}
    bool                 hasOwner()        const {return owner != 0;}
    const Subsystem&     getOwner()        const {assert(owner); return *owner;}
    int                  getIndexInOwner() const {assert(owner); return indexInOwner;}

    void                 setClientPlacement(const Placement& p) {client = &p;}
    bool                 hasClientPlacement() const {return client != 0;}
    const Placement&     getClientPlacement() const {assert(client); return *client;}

    String toString(const String& linePrefix) const;
    void checkPlacementValueConsistency(const Subsystem* expOwner, 
                                        int              expIndexInOwner,
                                        const Subsystem& expRoot) const;
private:
    void throwInvalid() const;

private:
    PlacementValue       value;        // The stored value
    bool                 valid;        // ... and whether it is meaningful.

    const Subsystem*     owner;        // The Subsystem which owns this PlacementValueSlot
    int                  indexInOwner; // ... and the index in its placementValueSlots list.

    const Placement*     client;       // The PlacementSlot holding the expression for which
                                       // this slot holds the Value. That PlacementSlots's value
                                       // pointer must point right back here.
};

/**
 * This is the abstract implementation of the PlacementValue client-side class.
 * Derived concrete classes hold specific kinds of values (typically numeric but
 * that is not a requirement).
 */
class PlacementValueRep {
public:
    PlacementValueRep() : myHandle(0) { }
    virtual ~PlacementValueRep() { }
    virtual PlacementValueRep* clone() const = 0;
    virtual std::string toString(const std::string& linePrefix) const = 0;

    void                  setMyHandle(PlacementValue& p) {myHandle = &p;}
    bool                  hasHandle()       const {return myHandle != 0;}
    const PlacementValue& getMyHandle()     const {assert(myHandle); return *myHandle;}
    PlacementValue&       updMyHandle()           {assert(myHandle); return *myHandle;} 

private:
    PlacementValue*  myHandle;     // the PlacementValue whose rep this is
};

template <class T> class PlacementValueRep_ : public PlacementValueRep {
public:
    PlacementValueRep_<T>() : PlacementValueRep() { }
    // default destructor is fine

    /// Copy constructor copies the value but the resulting copy has no handle.
    PlacementValueRep_<T>(const PlacementValueRep_<T>& v)
        : PlacementValueRep(), value(v.getValue()) { }

    /// Copy assignment changes the value but not the handle.
    PlacementValueRep_<T>& operator=(const PlacementValueRep_<T>& v) {
        if (&v != this)
            value = v.getValue();
        return *this;
    }

    /// This constructor makes a handle-less PlacementValueRep holding the given value.
    explicit PlacementValueRep_<T>(const T& v) : PlacementValueRep(), value(v) { }

    /// Assign to a concrete value.
    PlacementValueRep_<T>& operator=(const T& v) {
        value = v;
        return *this;
    }

    const PlacementValueRep_<T>& getMyHandle() const 
      { return PlacementValueRep_<T>::downcast(PlacementValueRep::getMyHandle()); }

    /// Clone produces a handle-less PlacementValueRep containing a copy of
    /// the value stored here.
    PlacementValueRep* clone() const { return new PlacementValueRep_<T>(*this); }

    std::string toString(const std::string&) const {
        std::stringstream s;
        s << TypeInfo<T>::name() << "(" << value << ")";   
        return s.str();
    }

    const T& getValue() const     {return value;}
    T&       updValue()           {return value;}
    void     setValue(const T& v) {value=v;}

    SIMTK_DOWNCAST(PlacementValueRep_<T>, PlacementValueRep);
private:
    T value;
};

} // namespace simtk

#endif // SIMTK_PLACEMENT_VALUE_REP_H_
