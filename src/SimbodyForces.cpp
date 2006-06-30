/* Portions copyright (c) 2006 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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

TwoPointSpringSubsystem::TwoPointSpringSubsystem() {
    rep = new TwoPointSpringSubsystemRep();
    rep->setMyHandle(*this);
}


int TwoPointSpringSubsystem::addLinearTwoPointSpring
   (int body1, const Vec3& s1,
    int body2, const Vec3& s2,
    const Real& stiffness,
    const Real& naturalLength) 
{
    return updRep().addLinearTwoPointSpring(body1,s1,body2,s2,stiffness,naturalLength);
}

int TwoPointSpringSubsystem::addGlobalMobilityDamping (const Real& dampingFactor) {
    return updRep().addGlobalMobilityDamping(dampingFactor);
}

} // namespace SimTK

