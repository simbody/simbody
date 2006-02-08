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
 * Implementation of PlacementValue and PlacementValue_<T> API methods.
 */

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/Geometry.h"
#include "simbody/internal/Mechanics.h"
#include "simbody/internal/PlacementValue.h"

#include "PlacementValueRep.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace simtk {


/////////////////////////////////////////////////////////////
// PlacementValue definitions & instantiations.            //
// Note that every supported subtype must be instantiated. //
/////////////////////////////////////////////////////////////

    // PLACEMENT VALUE //

PlacementValue::PlacementValue(const PlacementValue& src) : rep(0) {
    if (src.rep) {
        rep = src.rep->clone();
        rep->setMyHandle(*this);
    }
}

PlacementValue::PlacementValue(PlacementValueRep* r) : rep(r) {
    assert(rep && !rep->hasHandle());
    rep->setMyHandle(*this);
}

PlacementValue&
PlacementValue::operator=(const PlacementValue& src) {
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

PlacementValue::~PlacementValue() {
    // This will blow up if rep doesn't have a handle -- we shouldn't
    // be pointing to it in that case!
    if (rep && (&rep->getMyHandle() == this)) delete rep; 
    rep=0;
}

String PlacementValue::toString(const String& linePrefix) const {
    std::stringstream s;
    s << "Placement Value ";
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

std::ostream& operator<<(std::ostream& o, const PlacementValue& v) {
    return o << v.toString();
}

    // PLACEMENT VALUE <T> //

template <class T>
PlacementValue_<T>::PlacementValue_<T>() : PlacementValue() { 
    rep = new PlacementValueRep_<T>(); 
    rep->setMyHandle(*this);
}

template <class T>
PlacementValue_<T>::PlacementValue_<T>(const T& v) : PlacementValue() {
    rep = new PlacementValueRep_<T>(v);
    rep->setMyHandle(*this);
}

template <class T>
PlacementValue_<T>& PlacementValue_<T>::operator=(const T& v) {
    if (rep) set(v);
    else {rep = new PlacementValueRep_<T>(v); rep->setMyHandle(*this);}
    return *this;
}

template <class T>
const T& PlacementValue_<T>::get() const {
    return PlacementValueRep_<T>::downcast(getRep()).getValue();
}

template <class T>
T& PlacementValue_<T>::upd() {
    return PlacementValueRep_<T>::downcast(updRep()).updValue();
}

template <class T>
PlacementValue_<T>::operator const T&() const { return get(); }

template <class T>
void PlacementValue_<T>::set(const T& v) { 
    PlacementValueRep_<T>::downcast(updRep()).setValue(v);
}

template <class T>
void PlacementValue_<T>::initializeToValueType(PlacementValue& pv) {
    if (!pv.hasRep()) {
        pv.setRep(new PlacementValueRep_<T>());
        pv.updRep().setMyHandle(pv);
    }
}

template <class T>
bool PlacementValue_<T>::isInstanceOf(const PlacementValue& pv) {
    return PlacementValueRep_<T>::isA(pv.getRep());
}

template <class T>
const PlacementValue_<T>& PlacementValue_<T>::downcast(const PlacementValue& pv) {
    assert(isInstanceOf(pv)); 
    return reinterpret_cast<const PlacementValue_<T>&>(pv);
}

template <class T>
PlacementValue_<T>& PlacementValue_<T>::downcast(PlacementValue& pv) { 
    assert(isInstanceOf(pv)); 
    return reinterpret_cast<PlacementValue_<T>&>(pv);
}

template class PlacementValue_<Real>;
template class PlacementValue_<Vec3>;
template class PlacementValue_<Mat33>;
template class PlacementValue_<TransformMat>;
template class PlacementValue_<InertiaMat>;
template class PlacementValue_<RotationMat>;
template class PlacementValue_<UnitVec3>;
} // namespace simtk
