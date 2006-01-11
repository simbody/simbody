#ifndef IVM_SIMBODY_INTERFACE_H_
#define IVM_SIMBODY_INTERFACE_H_

/** @file
 *
 * This is a Simbody-compatible interface to the IVM-derived rigid body
 * simulation toolset.
 */

#include "simmatrix/BigMatrix.h"
#include "simbody/Simbody.h"

using namespace simtk;

#include <iostream>

class IVMSimbodyInterface /* : public RigidBodyMechanicsResource(?) */ {
public:
    explicit IVMSimbodyInterface(const Multibody&, bool oldStyle=false);

    int getNBodies() const;
    int getNParameters() const;
    int getNQ() const;
    int getNU() const;

    void realizeParameters(const State&) const;
    void realizeConfiguration(const State&) const;
    void realizeMotion(const State&) const;
    void realizeReaction(const State&) const;

    void enforceConfigurationConstraints(State&) const;
    void enforceMotionConstraints(State&) const;

    const Vector& getQDot(const State&) const;
    Vector calcUDot(const State&, 
                    const Array<SpatialVec>& bodyForces,
                    const Vector&            hingeForces) const;

    // Initialize the two force arrays, resizing if necessary and setting all to zero.
    // WARNING: forces may be accumulated in whatever manner is best for the underlying
    // computations. Be sure to use the access routines below; do not attempt to
    // fill in these arrays yourself.
    void clearForces(Array<SpatialVec>& bodyForces, Vector& hingeForces) const {
        bodyForces.resize(getNBodies()); hingeForces.resize(getNU());
        bodyForces = SpatialVec(Vec3(0.)); hingeForces = 0.;
    }

    State getDefaultState() const;

    // Given a station on a body (measured from & expressed in the body's frame), and a
    // force vector (expressed in the body's frame), convert that to an appropriate
    // spatial force and update the appropriate entry(s) in the bodyForces array.
    //
    // This is a Configuration-stage operator; that is, the State must already have
    // been realized to at least the Configuration level. Of course it is possible that 
    // calculating the applied force must be done even later.
    void applyPointForce(const State&, const Body&, const Vec3& pt, const Vec3& frc, 
                         Array<SpatialVec>& bodyForces) const;

    // Similar operation for a torque.
    void applyBodyTorque(const State&, const Body&, const Vec3& trq, 
                         Array<SpatialVec>& bodyForces) const;

    // Given a uniform acceleration field g expressed in the ground frame, apply an
    // appropriate gravitational force to each body. This is a Configuration-stage
    // operator because we have to know how the bodies are oriented in order to
    // figure out in which direction gravity is tugging at them.
    void applyGravity(const State&, const Vec3& g, Array<SpatialVec>& bodyForces) const;

    // Apply a force along or around a particular generalized coordinate. This is
    // a Modeling-stage operator -- all we need to know is the allocation of
    // generalized coordinates to joints.
    void applyHingeForce(const State&, const Joint&, int axis, const Real& frc, 
                         Vector& hingeForces) const;

    // TODO: these should return cache references
    Frame         getBodyConfiguration(const State&, const Body&) const;
    SpatialVec getBodyVelocity     (const State&, const Body&) const;
private:
    class IVMSimbodyInterfaceRep* rep;
    friend class IVMSimbodyInterfaceRep;
};



#endif /* IVM_SIMBODY_INTERFACE_H_ */
