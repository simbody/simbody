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
 * Implementation of ForceSubsystem, a still-abstract Subsystem.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"

#include "ForceSubsystemRep.h"

namespace SimTK {

    ////////////////////
    // ForceSubsystem //
    ////////////////////

// Default constructor is inline and creates an empty handle.
// Default copy & assignment just copy the parent class.
// Default destructor destructs the parent class.


// This is a PIMPL virtual method.
Real ForceSubsystem::calcPotentialEnergy(const State& s) const {
    return getRep().calcPotentialEnergy(s);
}

// This is a PIMPL virtual method.
void ForceSubsystem::addInForces(
    const State& s, const MatterSubsystem& matter,
    Vector_<SpatialVec>& rigidBodyForces,
    Vector_<Vec3>&       particleForces,
    Vector&              mobilityForces) const 
{
    getRep().addInForces(s, matter, 
                         rigidBodyForces, particleForces, mobilityForces);
}

void ForceSubsystem::setMatterSubsystemIndex(int subsys) {
    updRep().setMatterSubsystemIndex(subsys);
}
int ForceSubsystem::getMatterSubsystemIndex() const {
    return getRep().getMatterSubsystemIndex();
}

/*static*/ bool 
ForceSubsystem::isInstanceOf(const Subsystem& s) {
    return ForceSubsystemRep::isA(s.getRep());
}
/*static*/ const ForceSubsystem&
ForceSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const ForceSubsystem&>(s);
}
/*static*/ ForceSubsystem&
ForceSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<ForceSubsystem&>(s);
}


const ForceSubsystemRep& 
ForceSubsystem::getRep() const {
    return dynamic_cast<const ForceSubsystemRep&>(*rep);
}
ForceSubsystemRep&       
ForceSubsystem::updRep() {
    return dynamic_cast<ForceSubsystemRep&>(*rep);
}

} // namespace SimTK

