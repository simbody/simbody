/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */


/**@file
 *
 * Private implementation of DecorationSubsystem.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/DecorationSubsystem.h"

#include "DecorationSubsystemRep.h"

namespace SimTK {


    //////////////////////////
    // DECORATION SUBSYSTEM //
    //////////////////////////

/*static*/ bool 
DecorationSubsystem::isInstanceOf(const Subsystem& s) {
    return DecorationSubsystemGuts::isA(s.getSubsystemGuts());
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

const DecorationSubsystemGuts& 
DecorationSubsystem::getGuts() const {
    return dynamic_cast<const DecorationSubsystemGuts&>(getSubsystemGuts());
}
DecorationSubsystemGuts&       
DecorationSubsystem::updGuts() {
    return dynamic_cast<DecorationSubsystemGuts&>(updSubsystemGuts());
}

// Create Subsystem but don't associate it with any System. This isn't much use except
// for making std::vector's, which require a default constructor to be available.
DecorationSubsystem::DecorationSubsystem()
  : Subsystem()
{
     adoptSubsystemGuts(new DecorationSubsystemGuts());
}

DecorationSubsystem::DecorationSubsystem(MultibodySystem& mbs)
  : Subsystem() 
{
    adoptSubsystemGuts(new DecorationSubsystemGuts());
    mbs.setDecorationSubsystem(*this);
}

void DecorationSubsystem::addBodyFixedDecoration
   (MobilizedBodyId body, const Transform& X_GD, const DecorativeGeometry& g) 
{
    updGuts().addBodyFixedDecoration(body, X_GD, g);
}

void DecorationSubsystem::addRubberBandLine
   (MobilizedBodyId b1, const Vec3& station1,
    MobilizedBodyId b2, const Vec3& station2,
    const DecorativeLine& g)
{
    updGuts().addRubberBandLine(b1,station1,b2,station2,g);
}

    ///////////////////////////////
    // DECORATION SUBSYSTEM GUTS //
    ///////////////////////////////

// Return the MultibodySystem which owns this DecorationSubsystem.
const MultibodySystem& DecorationSubsystemGuts::getMultibodySystem() const {
    return MultibodySystem::downcast(getSystem());
}

int DecorationSubsystemGuts::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, Array<DecorativeGeometry>& geom) const
{
    switch(stage) {
    case Stage::Topology: {
        assert(subsystemTopologyHasBeenRealized());
        for (int i=0; i<(int)geometry.size(); ++i)
            geom.push_back(geometry[i]);
        break;
    }
    case Stage::Position: {
        assert(getStage(s) >= Stage::Position);
        const MultibodySystem&        mbs    = getMultibodySystem(); // my owner
        const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();
        for (int i=0; i<(int)rubberBandLines.size(); ++i) {
            const RubberBandLine& rb = rubberBandLines[i];
            geom.push_back(rb.line); // make a new copy
            DecorativeLine& line = DecorativeLine::updDowncast(geom.back()); // get access to copy
            line.setEndpoints(
                matter.getMobilizedBody(rb.body1).locateBodyPointOnGround(s,rb.station1),
                matter.getMobilizedBody(rb.body2).locateBodyPointOnGround(s,rb.station2));
        }
    }
    default: 
        assert(getStage(s) >= stage);
    }

    return 0;
}

} // namespace SimTK

