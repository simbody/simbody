/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
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


/**@file
 *
 * Implementation of SimbodyForcesSubsystem and SimbodyForcesSubsystemRep.
 */

#include "Simbody.h"
#include "simbody/internal/ForceSubsystem.h"
#include "SimbodyForcesRep.h"


namespace SimTK {


    //////////////////////////
    // EmptyForcesSubsystem //
    //////////////////////////

EmptyForcesSubsystem::EmptyForcesSubsystem() : ForceSubsystem() {
    rep = new EmptyForcesSubsystemRep();
    rep->setMyHandle(*this);
}

    /////////////////////////////
    // TwoPointSpringSubsystem //
    /////////////////////////////


/*static*/ bool 
TwoPointSpringSubsystem::isInstanceOf(const ForceSubsystem& s) {
    return TwoPointSpringSubsystemRep::isA(s.getRep());
}
/*static*/ const TwoPointSpringSubsystem&
TwoPointSpringSubsystem::downcast(const ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const TwoPointSpringSubsystem&>(s);
}
/*static*/ TwoPointSpringSubsystem&
TwoPointSpringSubsystem::updDowncast(ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<TwoPointSpringSubsystem&>(s);
}

const TwoPointSpringSubsystemRep& 
TwoPointSpringSubsystem::getRep() const {
    return dynamic_cast<const TwoPointSpringSubsystemRep&>(*rep);
}
TwoPointSpringSubsystemRep&       
TwoPointSpringSubsystem::updRep() {
    return dynamic_cast<TwoPointSpringSubsystemRep&>(*rep);
}

TwoPointSpringSubsystem::TwoPointSpringSubsystem
    (int b1, const Vec3& s1, int b2, const Vec3& s2,
     const Real& k, const Real& x0) : ForceSubsystem() 
{
  rep = new TwoPointSpringSubsystemRep(b1,s1,b2,s2,k,x0);
  rep->setMyHandle(*this);
}


const Vec3& TwoPointSpringSubsystem::getGravity(const State& s) const {
    return getRep().getGravity(s);
}
Vec3& TwoPointSpringSubsystem::updGravity(State& s) const {
    return getRep().updGravity(s);
}
const Real& TwoPointSpringSubsystem::getDamping(const State& s) const {
    return getRep().getDamping(s);
}
Real& TwoPointSpringSubsystem::updDamping(State& s) const {
    return getRep().updDamping(s);
}
const Real& TwoPointSpringSubsystem::getStiffness(const State& s) const {
    return getRep().getStiffness(s);
}
Real& TwoPointSpringSubsystem::updStiffness(State& s) const {
    return getRep().updStiffness(s);
}

const Real& TwoPointSpringSubsystem::getNaturalLength(const State& s) const {
    return getRep().getNaturalLength(s);
}
Real& TwoPointSpringSubsystem::updNaturalLength(State& s) const {
    return getRep().updNaturalLength(s);
}

const Real& TwoPointSpringSubsystem::getPotentialEnergy(const State& s) const {
    return getRep().getPotentialEnergy(s);
}
const Vec3& TwoPointSpringSubsystem::getForceOnStation1(const State& s) const {
    return getRep().getForceOnStation1(s);
}

    ////////////////////////////
    // SimbodyForcesSubsystem //
    ////////////////////////////
/*
SimbodyForcesSubsystem::SimbodyForcesSubsystem() : MechanicalForcesSubsystem() {
    rep = new SimbodyForcesSubsystemRep();
    rep->setMyHandle(*this);
}

SimbodyForcesSubsystemRep& SimbodyForcesSubsystem::updRep() {
    return SimbodyForcesSubsystemRep::downcast(*rep);
}
const SimbodyForcesSubsystemRep& SimbodyForcesSubsystem::getRep() {
    return SimbodyForcesSubsystemRep::downcast(*rep);
}

const Real& SimbodyForcesSubsystem::getPotentialEnergy(const State& s) const {
    return getRep().getPotentialEnergy(s);
}
*/
} // namespace SimTK

