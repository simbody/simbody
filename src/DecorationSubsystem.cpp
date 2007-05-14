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

namespace SimTK {

class DecorationSubsystemRep : public SubsystemRep {

    struct RubberBandLine {
        RubberBandLine(BodyId b1, const Vec3& s1,
                       BodyId b2, const Vec3& s2,
                       const DecorativeLine& l)
          : body1(b1), body2(b2), station1(s1), station2(s2), line(l)
        {
        }
        BodyId  body1, body2;
        Vec3 station1, station2;
        DecorativeLine line;
    };

public:
    DecorationSubsystemRep()
     : SubsystemRep("DecorationSubsystem", "0.0.1"), 
       built(false)
    {
    }

    // Return the MultibodySystem which owns this DecorationSubsystem.
    const MultibodySystem& getMultibodySystem() const {
        return MultibodySystem::downcast(getSystem());
    }

    // Add a permanent ("Topological") piece of geometry which is permanently fixed to a single body.
    // Thus the 3D polygonal representation can be precalculated once and for all at Stage::Topology,
    // then transformed at Stage::Position for display on a screen.
    //
    // This will make an internal copy of the supplied DecorativeGeometry. We'll save the
    // body Id in and apply the transform now to the saved copy, so that the geometry we return
    // later will be relative to the body frame only.
    void addBodyFixedDecoration(BodyId body, const Transform& X_BD, const DecorativeGeometry& g)
    {
        built = false;
        geometry.push_back(g); // make a new copy
        DecorativeGeometry& myg = geometry.back();
        myg.setBodyId(body);
        myg.setTransform(X_BD*myg.getTransform());
    }

    // This will make an internal copy of the supplied DecorativeGeometry.
    void addRubberBandLine(BodyId b1, const Vec3& station1, BodyId b2, const Vec3& station2,
                           const DecorativeLine& g)
    {
        built = false;
        rubberBandLines.push_back(RubberBandLine(b1,station1,b2,station2,g));
        DecorativeLine& myg = rubberBandLines.back().line; // new copy
        myg.setBodyId(GroundId);    // make sure the generated geometry will display properly
        myg.setTransform(Transform());
    }

    void calcDecorativeGeometryAndAppend(const State& s, Stage stage, Array<DecorativeGeometry>& geom) const {
        switch(stage) {
        case Stage::Topology: {
            assert(built);
            for (int i=0; i<(int)geometry.size(); ++i)
                geom.push_back(geometry[i]);
            break;
        }
        case Stage::Position: {
            assert(getStage(s) >= Stage::Position);
            const MultibodySystem& mbs    = getMultibodySystem(); // my owner
            const MatterSubsystem& matter = mbs.getMatterSubsystem();
            for (int i=0; i<(int)rubberBandLines.size(); ++i) {
                const RubberBandLine& rb = rubberBandLines[i];
                geom.push_back(rb.line); // make a new copy
                DecorativeLine& line = DecorativeLine::updDowncast(geom.back()); // get access to copy
                line.setEndpoints(
                    matter.locateBodyPointOnGround(s,rb.body1,rb.station1),
                    matter.locateBodyPointOnGround(s,rb.body2,rb.station2));
            }
        }
        default: 
            assert(getStage(s) >= stage);
        }
    }


    void realizeTopology(State& s) const {
        built = true;
    }

    void realizeModel(State& s) const {
        // Sorry, no choices available at the moment.
    }

    void realizeInstance(const State& s) const {
        // Nothing to compute here.
    }

    void realizeTime(const State& s) const {
        // Nothing to compute here.
    }

    void realizePosition(const State& s) const {
        // Nothing to compute here.
    }

    void realizeVelocity(const State& s) const {
        // Nothing to compute here.
    }

    void realizeDynamics(const State& s) const {
        // Nothing to compute here.
    }

    void realizeAcceleration(const State& s) const {
        // Nothing to compute here.
    }

    DecorationSubsystemRep* cloneSubsystemRep() const {
        return new DecorationSubsystemRep(*this);
    }

    SimTK_DOWNCAST(DecorationSubsystemRep, SubsystemRep);
private:
    // These must be filled in during realizeTopology and treated
    // as const thereafter. These are garbage unless built=true.
    mutable bool built;

    std::vector<DecorativeGeometry> geometry;
    std::vector<RubberBandLine>     rubberBandLines;
};


    /////////////////////////
    // DecorationSubsystem //
    /////////////////////////

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

DecorationSubsystem::DecorationSubsystem() {
    rep = new DecorationSubsystemRep();
    rep->setMyHandle(*this);
}

void DecorationSubsystem::addBodyFixedDecoration
   (BodyId body, const Transform& X_GD, const DecorativeGeometry& g) 
{
    updRep().addBodyFixedDecoration(body, X_GD, g);
}

void DecorationSubsystem::addRubberBandLine
   (BodyId b1, const Vec3& station1,
    BodyId b2, const Vec3& station2,
    const DecorativeLine& g)
{
    updRep().addRubberBandLine(b1,station1,b2,station2,g);
}


} // namespace SimTK

