#ifndef SimTK_SIMBODY_DECORATION_SUBSYSTEM_REP_H_
#define SimTK_SIMBODY_DECORATION_SUBSYSTEM_REP_H_

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
#include "simbody/internal/DecorationSubsystem.h"

namespace SimTK {

class DecorationSubsystemGuts : public Subsystem::Guts {

    struct RubberBandLine {
        RubberBandLine(MobilizedBodyIndex b1, const Vec3& s1,
                       MobilizedBodyIndex b2, const Vec3& s2,
                       const DecorativeLine& l)
          : body1(b1), body2(b2), station1(s1), station2(s2), line(l)
        {
        }
        MobilizedBodyIndex  body1, body2;
        Vec3 station1, station2;
        DecorativeLine line;
    };

public:
    DecorationSubsystemGuts()
      : Subsystem::Guts("DecorationSubsystem", "0.0.1"), generators(Stage::NValid)
    {
    }

    ~DecorationSubsystemGuts() {
        for (int i = 0; i < (int) generators.size(); i++)
            for (int j = 0; j < (int) generators[i].size(); j++)
                delete generators[i][j];
    }

    // Return the MultibodySystem which owns this DecorationSubsystem.
    const MultibodySystem& getMultibodySystem() const;

    // Add a permanent ("Topological") piece of geometry which is permanently fixed to a single body.
    // Thus the 3D polygonal representation can be precalculated once and for all at Stage::Topology,
    // then transformed at Stage::Position for display on a screen.
    //
    // This will make an internal copy of the supplied DecorativeGeometry. We'll save the
    // body Id in and apply the transform now to the saved copy, so that the geometry we return
    // later will be relative to the body frame only.
    void addBodyFixedDecoration(MobilizedBodyIndex body, const Transform& X_BD, const DecorativeGeometry& g)
    {
        invalidateSubsystemTopologyCache(); // this is a topological change
        geometry.push_back(g); // make a new copy
        DecorativeGeometry& myg = geometry.back();
        myg.setBodyId(body);
        myg.setTransform(X_BD*myg.getTransform());
    }

    // This will make an internal copy of the supplied DecorativeGeometry.
    void addRubberBandLine(MobilizedBodyIndex b1, const Vec3& station1, MobilizedBodyIndex b2, const Vec3& station2,
                           const DecorativeLine& g)
    {
        invalidateSubsystemTopologyCache(); // this is a topological change
        rubberBandLines.push_back(RubberBandLine(b1,station1,b2,station2,g));
        DecorativeLine& myg = rubberBandLines.back().line; // new copy
        myg.setBodyId(GroundIndex);    // make sure the generated geometry will display properly
        myg.setTransform(Transform());
    }

    void addDecorationGenerator(Stage stage, DecorationGenerator* generator) {
        generators[stage].push_back(generator);
    }

    DecorationSubsystemGuts* cloneImpl() const {
        return new DecorationSubsystemGuts(*this);
    }

    // Override virtual realize methods, although at the moment there is nothing here
    // so the default implementations would have been fine.

    int realizeSubsystemTopologyImpl(State& s) const {
        // No Topology-stage cache here
        return 0;
    }

    int realizeSubsystemModelImpl(State& s) const {
        // Sorry, no choices available at the moment.
        return 0;
    }

    int realizeSubsystemInstanceImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemTimeImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemPositionImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemVelocityImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemDynamicsImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemAccelerationImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemReportImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const;

    SimTK_DOWNCAST(DecorationSubsystemGuts, Subsystem::Guts);
private:
        // TOPOLOGY "STATE" VARIABLES
    Array_<DecorativeGeometry> geometry;
    Array_<RubberBandLine>     rubberBandLines;
    Array_<Array_<DecorationGenerator*> > generators;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_DECORATION_SUBSYSTEM_REP_H_







