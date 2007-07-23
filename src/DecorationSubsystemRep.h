#ifndef SimTK_DECORATION_SUBSYSTEM_REP_
#define SimTK_DECORATION_SUBSYSTEM_REP_

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

#include "SimTKcommon.h"
#include "SimTKcommon/internal/SubsystemGuts.h"

#include "simbody/internal/common.h"
#include "simbody/internal/DecorationSubsystem.h"

namespace SimTK {

class DecorationSubsystemGuts : public Subsystem::Guts {

    struct RubberBandLine {
        RubberBandLine(MobilizedBodyId b1, const Vec3& s1,
                       MobilizedBodyId b2, const Vec3& s2,
                       const DecorativeLine& l)
          : body1(b1), body2(b2), station1(s1), station2(s2), line(l)
        {
        }
        MobilizedBodyId  body1, body2;
        Vec3 station1, station2;
        DecorativeLine line;
    };

public:
    DecorationSubsystemGuts()
      : Subsystem::Guts("DecorationSubsystem", "0.0.1")
    {
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
    void addBodyFixedDecoration(MobilizedBodyId body, const Transform& X_BD, const DecorativeGeometry& g)
    {
        invalidateSubsystemTopologyCache(); // this is a topological change
        geometry.push_back(g); // make a new copy
        DecorativeGeometry& myg = geometry.back();
        myg.setBodyId(body);
        myg.setTransform(X_BD*myg.getTransform());
    }

    // This will make an internal copy of the supplied DecorativeGeometry.
    void addRubberBandLine(MobilizedBodyId b1, const Vec3& station1, MobilizedBodyId b2, const Vec3& station2,
                           const DecorativeLine& g)
    {
        invalidateSubsystemTopologyCache(); // this is a topological change
        rubberBandLines.push_back(RubberBandLine(b1,station1,b2,station2,g));
        DecorativeLine& myg = rubberBandLines.back().line; // new copy
        myg.setBodyId(GroundId);    // make sure the generated geometry will display properly
        myg.setTransform(Transform());
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
       (const State& s, Stage stage, Array<DecorativeGeometry>& geom) const;

    SimTK_DOWNCAST(DecorationSubsystemGuts, Subsystem::Guts);
private:
        // TOPOLOGY "STATE" VARIABLES
    std::vector<DecorativeGeometry> geometry;
    std::vector<RubberBandLine>     rubberBandLines;
};

} // namespace SimTK

#endif // SimTK_DECORATION_SUBSYSTEM_REP_







