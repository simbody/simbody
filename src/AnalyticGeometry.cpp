/* Copyright (c) 2006 Stanford University and Michael Sherman.
 * Contributors:
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
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/AnalyticGeometry.h"
#include "simbody/internal/DecorativeGeometry.h"

#include "AnalyticGeometryRep.h"

#include <cmath>

namespace SimTK {

    //////////////////////
    // AnalyticGeometry //
    //////////////////////


// This is an owner handle if there is no rep or if the rep points back
// to this handle.
bool AnalyticGeometry::isOwnerHandle() const {
    return rep==0 || rep->myHandle == this;
}

bool AnalyticGeometry::isEmptyHandle() const {return rep==0;}

AnalyticGeometry::~AnalyticGeometry() {
    if (isOwnerHandle())
        delete rep;
    rep = 0;
}

AnalyticGeometry::AnalyticGeometry(const AnalyticGeometry& src) : rep(0) {
    if (src.rep) {
        rep = src.rep->clone();
        rep->setMyHandle(*this);
    }
}

AnalyticGeometry& AnalyticGeometry::operator=(const AnalyticGeometry& src) {
    if (&src != this) {
        if (isOwnerHandle()) delete rep;
        rep = 0;
        if (src.rep) {
            rep = src.rep->clone();
            rep->setMyHandle(*this);
        }
    }
    return *this;
}

void AnalyticGeometry::setPlacement(const Transform& X_BG) {
    updRep().setPlacement(X_BG);
}
const Transform& AnalyticGeometry::getPlacement() const {
    return getRep().getPlacement();
}


DecorativeGeometry 
AnalyticGeometry::generateDecorativeGeometry() const {
    return DecorativeGeometry(*this);
}

    ///////////////////
    // AnalyticCurve //
    ///////////////////


    /////////////////////
    // AnalyticSurface //
    /////////////////////


    ////////////////////
    // AnalyticVolume //
    ////////////////////


    //////////////////
    // AnalyticLine //
    //////////////////

AnalyticLine::AnalyticLine(Real length) {
    rep = new AnalyticLineRep(length);
    rep->setMyHandle(*this);
}

} // namespace SimTK

