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

#include "simbody/SimbodyCommon.h"
#include "Placement.h"
#include "PlacementValue.h"
#include "PlacementValueRep.h"
#include "Feature.h"
#include "PlacementRep.h"
#include "FeatureRep.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace simtk {

    // PLACEMENT //

Placement::Placement(const Placement& src) : rep(0) { 
    if (src.rep) {
        rep = src.rep->clone();
        rep->setMyHandle(*this);
    }
}

Placement::Placement(PlacementRep* r) : rep(r) {
    assert(r && !r->hasHandle());
    r->setMyHandle(*this);
}

Placement& Placement::operator=(const Placement& src) {
    if (this != &src) {
        if (rep && (&rep->getMyHandle() == this)) delete rep; 
        rep=0;
        if (src.rep) {
            rep = src.rep->clone();
            rep->setMyHandle(*this);
        }
    }
    return *this;
}
Placement::~Placement() {
    // This will blow up if rep doesn't have a handle -- we shouldn't
    // be pointing to it in that case!
    if (rep && (&rep->getMyHandle() == this)) delete rep; 
    rep=0;
}

Placement::Placement(const Subsystem& s) : rep(0) {
    assert(Feature::isInstanceOf(s));
    const Feature& f = Feature::downcast(s);
    rep = f.getRep().createFeatureReference(*this);
}
Placement::Placement(const Subsystem& s, int i) : rep(0) {
    assert(Feature::isInstanceOf(s));
    const Feature& f = Feature::downcast(s);
    rep = f.getRep().createFeatureReference(*this, i);
}
Placement::Placement(const Real& r) : rep(0) {
    rep = new RealConstantPlacementRep(r);
    rep->setMyHandle(*this);
}
Placement::Placement(const Vec3& v) : rep(0) {
    rep = new Vec3ConstantPlacementRep(v);
    rep->setMyHandle(*this);
}
//Placement::Placement(const Mat33& m) : rep(0) {
//    rep = new Mat33ConstantPlacementRep(reinterpret_cast<RealPlacement&>(*this), r);
//}

// If this Placement expression depends on Features which need to be
// realized, we'll do that now. Then we will be able to do calcValue()
// when we want the value (that is, a Placement is an 'operator'
// rather than a 'response'.
void Placement::realize(Stage g) const {
    try {
        getRep().realize(g);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "realize", "", exc.getMessageText());
    }
}

PlacementValue Placement::calcValue() const {
    try {
        PlacementValue pv;
        getRep().evaluate(pv);
        return pv;
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "calcValue", "", exc.getMessageText());
    }
}

bool Placement::isConstant() const {
    return rep && rep->isConstant();
}

bool Placement::isFeatureReference() const {
    return rep && rep->isFeatureReference();
}

const Feature& 
Placement::getReferencedFeature() const {
    try {
        return getRep().getReferencedFeature();
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "getReferencedFeature", "", exc.getMessageText());
    }
}

String Placement::getPlacementTypeName() const {
    return rep ? rep->getPlacementTypeName() : "Void";
}

bool Placement::dependsOn(const Feature& f) const {
    return rep && rep->dependsOn(f);
}

String Placement::toString(const String& linePrefix) const {
    std::stringstream s;
    s << "Placement ";
    if (!rep) {
        s << "at 0x" << this << " HAS NULL REP";
        return s.str();
    }
    if (&rep->getMyHandle() != this) {
        s << "at 0x" << this << " HAS MISMATCHED REP";
        return s.str();
    }
    s << " " << rep->toString(linePrefix);
    return s.str();
}

std::ostream& operator<<(std::ostream& o, const Placement& p) {
    return o << p.toString() << std::endl;
}

// unary
Placement negate   (const Placement& p) {return p.getRep().genericNegate();}
Placement abs      (const Placement& p) {return p.getRep().genericAbs();}
Placement sqrt     (const Placement& p) {return p.getRep().genericSqrt();}
Placement exp      (const Placement& p) {return p.getRep().genericExp();}
Placement log      (const Placement& p) {return p.getRep().genericLog();}
Placement sin      (const Placement& p) {return p.getRep().genericSin();}
Placement cos      (const Placement& p) {return p.getRep().genericCos();}
Placement asin     (const Placement& p) {return p.getRep().genericAsin();}
Placement acos     (const Placement& p) {return p.getRep().genericAcos();}
Placement length   (const Placement& p) {return p.getRep().genericLength();}
Placement normalize(const Placement& p) {return p.getRep().genericNormalize();}

// binary
Placement add      (const Placement& l, const Placement& r) {return l.getRep().genericAdd(r);} 
Placement subtract (const Placement& l, const Placement& r) {return l.getRep().genericSub(r);}
Placement multiply (const Placement& l, const Placement& r) {return l.getRep().genericMul(r);}
Placement divide   (const Placement& l, const Placement& r) {return l.getRep().genericDvd(r);}
Placement distance (const Placement& l, const Placement& r) {return l.getRep().genericDistance(r);}
Placement angle    (const Placement& l, const Placement& r) {return l.getRep().genericAngle(r);}
Placement dot      (const Placement& l, const Placement& r) {return l.getRep().genericDotProduct(r);}
Placement cross    (const Placement& l, const Placement& r) {return l.getRep().genericCrossProduct(r);}

} // namespace simtk
