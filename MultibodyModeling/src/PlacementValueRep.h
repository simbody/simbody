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

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/PlacementValue.h"

#include <string>
#include <cassert>
#include <sstream>

namespace simtk {

class Subsystem;
class Placement;


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
