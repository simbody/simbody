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

    void setMyHandle(Body& h) {myHandle = &h;}
    const Body& getMyBodyHandle() const {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}
private:
    friend class Body;
    Body* myHandle;
};

class Body::Ground::GroundRep : public Body::BodyRep {
public:
    GroundRep() : BodyRep() {
    }
    GroundRep* clone() const {
        return new GroundRep(*this);
    }
    const MassProperties& getDefaultRigidBodyMassProperties() const {
        static const MassProperties groundMass(NTraits<Real>::getInfinity(), 
                                               Vec3(0),
                                               NTraits<Real>::getInfinity()*Inertia(1));
        return groundMass;
    }

    const Body::Ground& getMyGroundBodyHandle() const {
        return reinterpret_cast<const Body::Ground&>(getMyBodyHandle());
    }

    SimTK_DOWNCAST(GroundRep, BodyRep);
private:
};

class Body::Rigid::RigidRep : public Body::BodyRep {
public:
    RigidRep() : BodyRep() {
    }
    RigidRep(const MassProperties& m) : BodyRep(), defaultMassProperties(m) {
    }
    const MassProperties& getDefaultRigidBodyMassProperties() const {
        return defaultMassProperties;
    }
    void setDefaultMassProperties(const MassProperties& m) {
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