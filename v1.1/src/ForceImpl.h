#ifndef SimTK_SIMBODY_FORCE_IMPL_H_
#define SimTK_SIMBODY_FORCE_IMPL_H_
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Force.h"

namespace SimTK {

class ForceImpl : public PIMPLImplementation<Force, ForceImpl> {
public:
    ForceImpl() {
    }
    ForceImpl(const ForceImpl& clone) {
        *this = clone;
    }
    virtual ~ForceImpl() {
    }
    virtual ForceImpl* clone() const = 0;
    virtual bool dependsOnlyOnPositions() const {
        return false;
    }
    ForceIndex getForceIndex() const {
        return index;
    }
    void setForceIndex(ForceIndex ind) {
        index = ind;
    }
    virtual void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const = 0;
private:
    ForceIndex index;
};

class Force::TwoPointLinearSpringImpl : public ForceImpl {
public:
    TwoPointLinearSpringImpl(const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real k, Real x0);
    TwoPointLinearSpringImpl* clone() const {
        return new TwoPointLinearSpringImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body1, body2;
    Vec3 station1, station2;
    Real k, x0;
};

class Force::TwoPointLinearDamperImpl : public ForceImpl {
public:
    TwoPointLinearDamperImpl(const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real damping);
    TwoPointLinearDamperImpl* clone() const {
        return new TwoPointLinearDamperImpl(*this);
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body1, body2;
    Vec3 station1, station2;
    Real damping;
};

class Force::TwoPointConstantForceImpl : public ForceImpl {
public:
    TwoPointConstantForceImpl(const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real force);
    TwoPointConstantForceImpl* clone() const {
        return new TwoPointConstantForceImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body1, body2;
    Vec3 station1, station2;
    Real force;
};

class Force::MobilityLinearSpringImpl : public ForceImpl {
public:
    MobilityLinearSpringImpl(const MobilizedBody& body, int coordinate, Real k, Real x0);
    MobilityLinearSpringImpl* clone() const {
        return new MobilityLinearSpringImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body;
    int coordinate;
    Real k, x0;
};

class Force::MobilityLinearDamperImpl : public ForceImpl {
public:
    MobilityLinearDamperImpl(const MobilizedBody& body, int coordinate, Real damping);
    MobilityLinearDamperImpl* clone() const {
        return new MobilityLinearDamperImpl(*this);
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body;
    int coordinate;
    Real damping;
};

class Force::MobilityConstantForceImpl : public ForceImpl {
public:
    MobilityConstantForceImpl(const MobilizedBody& body, int coordinate, Real force);
    MobilityConstantForceImpl* clone() const {
        return new MobilityConstantForceImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body;
    int coordinate;
    Real force;
};

class Force::ConstantForceImpl : public ForceImpl {
public:
    ConstantForceImpl(const MobilizedBody& body, const Vec3& station, const Vec3& force);
    ConstantForceImpl* clone() const {
        return new ConstantForceImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body;
    Vec3 station, force;
};

class Force::ConstantTorqueImpl : public ForceImpl {
public:
    ConstantTorqueImpl(const MobilizedBody& body, const Vec3& torque);
    ConstantTorqueImpl* clone() const {
        return new ConstantTorqueImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body;
    Vec3 torque;
};

class Force::GlobalDamperImpl : public ForceImpl {
public:
    GlobalDamperImpl(const SimbodyMatterSubsystem& matter, Real damping);
    GlobalDamperImpl* clone() const {
        return new GlobalDamperImpl(*this);
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const;
private:
    const SimbodyMatterSubsystem& matter;
    Real damping;
};

class Force::UniformGravityImpl : public ForceImpl {
public:
    UniformGravityImpl(const SimbodyMatterSubsystem& matter, const Vec3& g, Real zeroHeight);
    UniformGravityImpl* clone() const {
        return new UniformGravityImpl(*this);
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const;
private:
    const SimbodyMatterSubsystem& matter;
    Vec3 g;
    Real zeroHeight;
};

class Force::CustomImpl : public ForceImpl {
public:
    CustomImpl(Force::Custom::Implementation* implementation);
    CustomImpl* clone() const {
        return new CustomImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return implementation->dependsOnlyOnPositions();
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const;
    ~CustomImpl() {
        delete implementation;
    }
    const Force::Custom::Implementation& getImplementation() const {
        return *implementation;
    }
    Force::Custom::Implementation& updImplementation() {
        return *implementation;
    }
private:
    Force::Custom::Implementation* implementation;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_IMPL_H_
