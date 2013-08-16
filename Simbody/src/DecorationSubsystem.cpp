/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
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
    return static_cast<const DecorationSubsystem&>(s);
}
/*static*/ DecorationSubsystem&
DecorationSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return static_cast<DecorationSubsystem&>(s);
}

const DecorationSubsystemGuts& 
DecorationSubsystem::getGuts() const {
    return SimTK_DYNAMIC_CAST_DEBUG<const DecorationSubsystemGuts&>(getSubsystemGuts());
}
DecorationSubsystemGuts&       
DecorationSubsystem::updGuts() {
    return SimTK_DYNAMIC_CAST_DEBUG<DecorationSubsystemGuts&>(updSubsystemGuts());
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
   (MobilizedBodyIndex body, const Transform& X_GD, const DecorativeGeometry& g) 
{
    updGuts().addBodyFixedDecoration(body, X_GD, g);
}

void DecorationSubsystem::addRubberBandLine
   (MobilizedBodyIndex b1, const Vec3& station1,
    MobilizedBodyIndex b2, const Vec3& station2,
    const DecorativeLine& g)
{
    updGuts().addRubberBandLine(b1,station1,b2,station2,g);
}

void DecorationSubsystem::addDecorationGenerator(Stage stage, DecorationGenerator* generator) {
    updGuts().addDecorationGenerator(stage, generator);
}

    ///////////////////////////////
    // DECORATION SUBSYSTEM GUTS //
    ///////////////////////////////

// Return the MultibodySystem which owns this DecorationSubsystem.
const MultibodySystem& DecorationSubsystemGuts::getMultibodySystem() const {
    return MultibodySystem::downcast(getSystem());
}

int DecorationSubsystemGuts::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
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
                matter.getMobilizedBody(rb.body1).findStationLocationInGround(s,rb.station1),
                matter.getMobilizedBody(rb.body2).findStationLocationInGround(s,rb.station2));
        }
    }
    default: 
        assert(getStage(s) >= stage);
    }
    for (int i = 0; i < (int) generators[stage].size(); i++)
        generators[stage][i]->generateDecorations(s, geom);

    return 0;
}

} // namespace SimTK

