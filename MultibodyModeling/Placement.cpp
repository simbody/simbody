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
 * Implementation of Placement handles.
 */

#include "SimbodyCommon.h"
#include "Placement.h"
#include "PlacementRep.h"

#include <string>
#include <iostream> 
#include <sstream>
using std::endl;
using std::ostream;

namespace simtk {

    // PLACEMENT //
Placement::Placement(const Placement& src) { 
    rep = src.rep ? src.rep->clone(*this) : 0;
}
Placement& Placement::operator=(const Placement& src) {
    if (this != &src) {
        delete rep;
        rep = src.rep ? src.rep->clone(*this) : 0;
    }
    return *this;
}
Placement::~Placement() {
    if (rep==0) return;
    assert(&rep->getHandle() == this);
    delete rep; rep=0;
}
const Feature& Placement::getOwner() const {
    assert(rep && rep->hasOwner());
    assert(&rep->getHandle() == this);
    return rep->getOwner();
}
String Placement::toString(const String& linePrefix) const {
    std::stringstream s;
    s << "Placement ";
    if (!rep) {
        s << "at 0x" << this << " HAS NULL REP";
        return s.str();
    }
    if (&rep->getHandle() != this) {
        s << "at 0x" << this << " HAS MISMATCHED REP";
        return s.str();
    }
    s << rep->toString(linePrefix);
    return s.str();
}

std::ostream& operator<<(std::ostream& o, const Placement& p) {
    return o << p.toString() << endl;
}



    // STATION PLACEMENT //
StationPlacement::StationPlacement(const Vec3& v) {
    rep = new StationConstantPlacementRep(*this,v);
}
DirectionPlacement::DirectionPlacement(const Vec3& v) {
    rep = new DirectionConstantPlacementRep(*this,v);
}
} // namespace simtk
