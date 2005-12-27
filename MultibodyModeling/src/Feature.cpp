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
 * Implementations of API-level Feature methods.
 */

#include "simbody/SimbodyCommon.h"
#include "SubsystemRep.h"
#include "Feature.h"
#include "FeatureRep.h"
#include "BasicFeatures.h"

#include <string>
#include <iostream> 
#include <sstream>


namespace simtk {

Feature::Feature(const Feature& src) : Subsystem(src) { }
Feature& Feature::operator=(const Feature& src)
  { Subsystem::operator=(src); return *this; }
Feature::~Feature() { }

/*static*/ bool             
Feature::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return FeatureRep::isA(f.getRep());
}
/*static*/ const Feature& 
Feature::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Feature&>(f);
}

/*static*/ Feature&       
Feature::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Feature&>(f);
}

static String featureHasNoRep(const Feature& f) {
    std::ostringstream s;
    s << "<FEATURE AT 0x" << &f << " WITH NULL REP>";
    return String(s.str());
}

String Feature::getFeatureTypeName() const {
    return hasRep() ? String(getRep().getFeatureTypeName()) 
                    : featureHasNoRep(*this);
}

const PlacementValue& Feature::getValue() const {
    try {
        return getRep().getPlacementSlot().getValueSlot().getValue();
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW4(Exception::FeatureAPIMethodFailed, getFullName(), 
            "getValue", "", exc.getMessageText());
    }
}

bool Feature::hasPlacement() const {
    return getRep().hasPlacement();
}
void Feature::place(const Placement& p) {
    try {
        updRep().place(p);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW2(Exception::APIMethodFailed, "Feature::place()", exc.getMessage());
    }
}
void Feature::replace(const Placement& p) {
    try {
        updRep().replace(p);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW2(Exception::APIMethodFailed, "Feature::replace()", exc.getMessage());
    }
}
const Placement& Feature::getPlacement() const {
    assert(hasPlacement());
    return getRep().getPlacementSlot().getPlacement();
}

} // namespace simtk
