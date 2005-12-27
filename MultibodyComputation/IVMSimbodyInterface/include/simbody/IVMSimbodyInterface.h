#ifndef IVM_SIMBODY_INTERFACE_H_
#define IVM_SIMBODY_INTERFACE_H_

/** @file
 *
 * This is a Simbody-compatible interface to the IVM-derived rigid body
 * simulation toolset.
 */

#include "simmatrix/BigMatrix.h"
#include "simbody/Simbody.h"

using simtk::State;
using simtk::Vector;
using simtk::Multibody;
using simtk::SpatialVector;
using simtk::Frame;

#include <iostream>

class IVMSimbodyInterface /* : public RigidBodyMechanicsResource(?) */ {
public:
    IVMSimbodyInterface(const Multibody&);

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

    const Frame&         getBodyConfiguration(const State&, int body) const;
    const SpatialVector& getBodyVelocity     (const State&, int body) const;
    const SpatialVector& getBodyAcceleration (const State&, int body) const;

};



#endif /* IVM_SIMBODY_INTERFACE_H_ */
