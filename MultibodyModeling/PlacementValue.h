#ifndef SIMTK_PLACEMENT_VALUE_H_
#define SIMTK_PLACEMENT_VALUE_H_

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
 * Definitions for the API-level, client side PlacementValue classes.
 */

#include "SimbodyCommon.h"

#include <iostream>

namespace simtk {

class Subsystem;

/**
 * This class represents a PlacementValue of unknown type. This is ultimately
 * just an ordinary value object, such as a number, with no additional bookkeeping
 * done here. However, these are kept in PlacementValueSlots in Subsystems to serve
 * as controlled values for Placements, keeping track of current validity.
 */
class PlacementValue {
public:
    PlacementValue() : rep(0) { }
    PlacementValue(const PlacementValue&);
    PlacementValue& operator=(const PlacementValue&);
    ~PlacementValue();

    String toString(const String& linePrefix="") const;

protected:
    class PlacementValueRep* rep;
    friend class PlacementValueRep;

public:
    // internal use only
    explicit PlacementValue(class PlacementValueRep*);
    bool                     hasRep() const {return rep != 0;}
    const PlacementValueRep& getRep() const {assert(rep); return *rep;}
    PlacementValueRep&       updRep()       {assert(rep); return *rep;}
    void                     setRep(PlacementValueRep* p) {assert(!rep); rep=p;}
};
std::ostream& operator<<(std::ostream& o, const PlacementValue&);

/**
 * These are the concrete PlacementValue classes. They must be instantiated
 * in the library for every supported type; despite appearances this is not
 * extensible on the client side.
 */
template <class T> class PlacementValue_ : public PlacementValue {
public:
    PlacementValue_<T>(); // creates an invalid value of this type
    explicit PlacementValue_<T>(const T&);
    PlacementValue_<T>& operator=(const T&);
    const T& get() const;
    T& upd();
    void set(const T&);

    // implicit conversion to type T
    operator const T&() const;

    static bool                      isInstanceOf(const PlacementValue&);
    static const PlacementValue_<T>& downcast(const PlacementValue&);
    static PlacementValue_<T>&       downcast(PlacementValue&);
};

} // namespace simtk

#endif // SIMTK_PLACEMENT_VALUE_H_
