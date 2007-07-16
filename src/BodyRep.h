#ifndef SimTK_SIMBODY_BODY_REP_H_
#define SimTK_SIMBODY_BODY_REP_H_

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
 * Private implementation of Body and its built-in subclasses.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Body.h"
#include "simbody/internal/DecorativeGeometry.h"

namespace SimTK {

    ///////////////
    // BODY REPS //
    ///////////////

class Body::BodyRep {
public:
    BodyRep() : myHandle(0) {
    }
    virtual ~BodyRep() { }
    virtual BodyRep* clone() const = 0;

    virtual const MassProperties& getDefaultRigidBodyMassProperties() const = 0;
    virtual void setDefaultRigidBodyMassProperties(const MassProperties&) = 0;

    // Add a permanent ("Topological") piece of geometry which is permanently fixed to the body.
    // Thus the 3D polygonal representation can be precalculated once and for all at Stage::Topology,
    // then transformed at Stage::Position for display on a screen.
    //
    // This will make an internal copy of the supplied DecorativeGeometry. We can't fill
    // in the MobilizedBodyId yet, but we apply the transform now to the saved copy so that
    // the geometry we return later will be relative to the body frame only.
    void addDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        geometry.push_back(g); // make a new copy
        DecorativeGeometry& myg = geometry.back();
        myg.setTransform(X_BD*myg.getTransform());
    }

    void appendDecorativeGeometry(MobilizedBodyId id, Array<DecorativeGeometry>& geom) const {
        for (int i=0; i<(int)geometry.size(); ++i) {
            geom.push_back(geometry[i]);
            geom.back().setBodyId(id);
        }
    }

    void setMyHandle(Body& h) {myHandle = &h;}
    const Body& getMyBodyHandle() const {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}
private:
    friend class Body;
    Body* myHandle;

        // TOPOLOGY "STATE" VARIABLES
    std::vector<DecorativeGeometry> geometry;
};

class Body::Ground::GroundRep : public Body::BodyRep {
public:
    GroundRep() : BodyRep(), infiniteMassProperties(NTraits<Real>::getInfinity(), 
                                                    Vec3(0),
                                                    NTraits<Real>::getInfinity()*Inertia(1))
    {
    }
    GroundRep* clone() const {
        return new GroundRep(*this);
    }
    const MassProperties& getDefaultRigidBodyMassProperties() const {
        static const MassProperties groundMass(NTraits<Real>::getInfinity(), 
                                               Vec3(0),
                                               NTraits<Real>::getInfinity()*Inertia(1));
        return infiniteMassProperties;
    }
    
    void setDefaultRigidBodyMassProperties(const MassProperties&) {
        SimTK_THROW1(Exception::Cant, "You can't change Ground's mass properties!");
    }

    const Body::Ground& getMyGroundBodyHandle() const {
        return reinterpret_cast<const Body::Ground&>(getMyBodyHandle());
    }

    SimTK_DOWNCAST(GroundRep, BodyRep);
private:
    const MassProperties infiniteMassProperties;
};

class Body::Massless::MasslessRep : public Body::BodyRep {
public:
    MasslessRep() : BodyRep(), zeroMassProperties(0, Vec3(0), Inertia(0)) {
    }
    MasslessRep* clone() const {
        return new MasslessRep(*this);
    }
    const MassProperties& getDefaultRigidBodyMassProperties() const {
        return zeroMassProperties;
    }
    
    void setDefaultRigidBodyMassProperties(const MassProperties&) {
        SimTK_THROW1(Exception::Cant, "You can't change a massless body's mass properties!");
    }

    const Body::Massless& getMyMasslessBodyHandle() const {
        return reinterpret_cast<const Body::Massless&>(getMyBodyHandle());
    }

    SimTK_DOWNCAST(MasslessRep, BodyRep);
private:
    const MassProperties zeroMassProperties;
};

class Body::Rigid::RigidRep : public Body::BodyRep {
public:
    RigidRep() : BodyRep(), defaultMassProperties(1,Vec3(0),Inertia(1)) {
    }
    explicit RigidRep(const MassProperties& m) : BodyRep(), defaultMassProperties(m) {
    }
    const MassProperties& getDefaultRigidBodyMassProperties() const {
        return defaultMassProperties;
    }
    void setDefaultRigidBodyMassProperties(const MassProperties& m) {
        defaultMassProperties = m;
    }
    RigidRep* clone() const {
        return new RigidRep(*this);
    }
    const Body::Rigid& getMyRigidBodyHandle() const {
        return reinterpret_cast<const Body::Rigid&>(getMyBodyHandle());
    }

    SimTK_DOWNCAST(RigidRep, BodyRep);
private:
    MassProperties defaultMassProperties;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_BODY_REP_H_