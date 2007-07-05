/* Portions copyright (c) 2007 Stanford University and Michael Sherman.
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
 * Private implementation of DecorationSubsystem.
 */

#include "SimTKsimbody.h"
#include "simbody/internal/Subsystem.h"
#include "simbody/internal/DecorationSubsystem.h"

#include "SubsystemRep.h"
#include "DecorationSubsystemRep.h"

namespace SimTK {


    //////////////////////////
    // DECORATION SUBSYSTEM //
    //////////////////////////

/*static*/ bool 
DecorationSubsystem::isInstanceOf(const Subsystem& s) {
    return DecorationSubsystemRep::isA(s.getRep());
}
/*static*/ const DecorationSubsystem&
DecorationSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const DecorationSubsystem&>(s);
}
/*static*/ DecorationSubsystem&
DecorationSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<DecorationSubsystem&>(s);
}

const DecorationSubsystemRep& 
DecorationSubsystem::getRep() const {
    assert(rep);
    return dynamic_cast<const DecorationSubsystemRep&>(*rep);
}
DecorationSubsystemRep&       
DecorationSubsystem::updRep() {
    assert(rep);
    return dynamic_cast<DecorationSubsystemRep&>(*rep);
}

// Create Subsystem but don't associate it with any System. This isn't much use except
// for making std::vector's, which require a default constructor to be available.
DecorationSubsystem::DecorationSubsystem()
  : Subsystem()
{
    rep = new DecorationSubsystemRep();
    rep->setMyHandle(*this);
}

DecorationSubsystem::DecorationSubsystem(MultibodySystem& mbs)
  : Subsystem() 
{
    rep = new DecorationSubsystemRep();
    rep->setMyHandle(*this);
    mbs.setDecorationSubsystem(*this);
}

void DecorationSubsystem::addBodyFixedDecoration
   (MobilizedBodyId body, const Transform& X_GD, const DecorativeGeometry& g) 
{
    updRep().addBodyFixedDecoration(body, X_GD, g);
}

void DecorationSubsystem::addRubberBandLine
   (MobilizedBodyId b1, const Vec3& station1,
    MobilizedBodyId b2, const Vec3& station2,
    const DecorativeLine& g)
{
    updRep().addRubberBandLine(b1,station1,b2,station2,g);
}


} // namespace SimTK

