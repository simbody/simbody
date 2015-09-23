#ifndef SimTK_SIMBODY_BODY_REP_H_
#define SimTK_SIMBODY_BODY_REP_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-15 Stanford University and the Authors.        *
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
 * Private implementation of Body and its built-in subclasses.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Body.h"

namespace SimTK {

//==============================================================================
//                                 BODY REP
//==============================================================================
class Body::BodyRep {
public:
    BodyRep() : myHandle(0) {
    }
    virtual ~BodyRep() { }
    virtual BodyRep* clone() const = 0;

    virtual const MassProperties& getDefaultRigidBodyMassProperties() const = 0;
    virtual void setDefaultRigidBodyMassProperties(const MassProperties&) = 0;

    // Add a permanent ("Topological") piece of geometry which is permanently
    // fixed to the body. Thus the 3D polygonal representation can be
    // precalculated once and for all at Stage::Topology, then transformed at
    // Stage::Position for display on a screen.
    //
    // This will make an internal copy of the supplied DecorativeGeometry. We
    // can't fill in the MobilizedBodyIndex yet, but we apply the transform
    // now to the saved copy so that the geometry we return later will be
    // relative to the body frame only.
    //
    // Return an index that can be used to find this decoration later; the
    // index must be the same in all copies of this Body. The same index is
    // stored in the copied DecorativeGeometry object.
    int addDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        const int nxt = (int)decorations.size();
        decorations.push_back(g); // make a new copy
        DecorativeGeometry& myg = decorations.back();
        myg.setIndexOnBody(nxt);
        myg.setTransform(X_BD*myg.getTransform());
        return nxt;
    }

    void appendDecorativeGeometry(MobilizedBodyIndex id,
                                  Array_<DecorativeGeometry>& geom) const
    {
        for (int i=0; i<(int)decorations.size(); ++i) {
            geom.push_back(decorations[i]); // copies index and userRef
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
    Array_<std::pair<Transform,ContactSurface> >      surfaces;
    Array_<DecorativeGeometry>                        decorations;
};



//==============================================================================
//                                 RIGID REP
//==============================================================================
class Body::Rigid::RigidRep : public Body::BodyRep {
public:
    RigidRep()
    :   BodyRep(), defaultMassProperties(1,Vec3(0),Inertia(1))
    {
    }
    explicit RigidRep(const MassProperties& m) : BodyRep(), defaultMassProperties(m) {
    }
    const MassProperties& getDefaultRigidBodyMassProperties() const override {
        return defaultMassProperties;
    }
    void setDefaultRigidBodyMassProperties(const MassProperties& m) override {
        defaultMassProperties = m;
    }
    RigidRep* clone() const override {
        return new RigidRep(*this);
    }
    const Body::Rigid& getMyRigidBodyHandle() const {
        return static_cast<const Body::Rigid&>(getMyBodyHandle());
    }

    SimTK_DOWNCAST(RigidRep, BodyRep);
private:
    MassProperties defaultMassProperties;
};



//==============================================================================
//                                 GROUND REP
//==============================================================================
class Body::Ground::GroundRep : public Body::BodyRep {
public:
    GroundRep()
    :   BodyRep(), infiniteMassProperties(Infinity, Vec3(0), Inertia(Infinity))
    {
    }
    GroundRep* clone() const override {
        return new GroundRep(*this);
    }
    const MassProperties& getDefaultRigidBodyMassProperties() const override {
        return infiniteMassProperties;
    }

    void setDefaultRigidBodyMassProperties(const MassProperties&) override {
        SimTK_THROW1(Exception::Cant, "You can't change Ground's mass properties!");
    }

    const Body::Ground& getMyGroundBodyHandle() const {
        return static_cast<const Body::Ground&>(getMyBodyHandle());
    }

    SimTK_DOWNCAST(GroundRep, BodyRep);
private:
    const MassProperties infiniteMassProperties;
};



//==============================================================================
//                              MASSLESS REP
//==============================================================================
class Body::Massless::MasslessRep : public Body::BodyRep {
public:
    MasslessRep()
    :   BodyRep(), zeroMassProperties(0, Vec3(0), Inertia(0))
    {
    }
    MasslessRep* clone() const override {
        return new MasslessRep(*this);
    }
    const MassProperties& getDefaultRigidBodyMassProperties() const override {
        return zeroMassProperties;
    }

    void setDefaultRigidBodyMassProperties(const MassProperties&) override {
        SimTK_THROW1(Exception::Cant, "You can't change a massless body's mass properties!");
    }

    const Body::Massless& getMyMasslessBodyHandle() const {
        return static_cast<const Body::Massless&>(getMyBodyHandle());
    }

    SimTK_DOWNCAST(MasslessRep, BodyRep);
private:
    const MassProperties zeroMassProperties;
};


} // namespace SimTK

#endif // SimTK_SIMBODY_BODY_REP_H_



