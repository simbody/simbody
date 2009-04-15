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
 * Portions copyright (c) 2008-9 Stanford University and the Authors.         *
 * Authors: Peter Eastman, Michael Sherman                                    *
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
	ForceImpl() : forces(0) {
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
    ForceIndex getForceIndex() const {return index;}
	const GeneralForceSubsystem& getForceSubsystem() const {assert(forces); return *forces;}
	void setForceSubsystem(GeneralForceSubsystem& frcsub, ForceIndex ix) {
		forces = &frcsub;
		index  = ix;
	}
	void invalidateTopologyCache() const {
		if (forces) forces->invalidateSubsystemTopologyCache();
	}

    virtual void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const = 0;
    virtual Real calcPotentialEnergy(const State& state) const = 0;
    virtual void realizeTopology(State& state) const {
    }
    virtual void realizeModel(State& state) const {
    }
    virtual void realizeInstance(const State& state) const {
    }
    virtual void realizeTime(const State& state) const {
    }
    virtual void realizePosition(const State& state) const {
    }
    virtual void realizeVelocity(const State& state) const {
    }
    virtual void realizeDynamics(const State& state) const {
    }
    virtual void realizeAcceleration(const State& state) const {
    }
    virtual void realizeReport(const State& state) const {
    }
private:
    GeneralForceSubsystem* forces;	// just a reference; don't delete on destruction
    ForceIndex			   index;
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
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
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
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
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
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
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
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
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
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
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
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
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
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
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
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
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
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    Real damping;
};

class Force::ThermostatImpl : public ForceImpl {
public:
    ThermostatImpl(const SimbodyMatterSubsystem& matter, 
				   Real boltzmannsConstant, Real defBathTemp, Real defRelaxationTime)
    :   matter(matter), Kb(boltzmannsConstant), defaultNumChains(DefaultDefaultNumChains),
        defaultBathTemp(defBathTemp), defaultRelaxationTime(defRelaxationTime)
	{}

    ThermostatImpl* clone() const {return new ThermostatImpl(*this);}
    bool dependsOnlyOnPositions() const {return false;}

    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
				   Vector_<Vec3>& particleForces, Vector& mobilityForces) const;

	// Temperature does not contribute to potential energy.
	Real calcPotentialEnergy(const State& state) const {return 0;}

    void realizeTopology(State& state) const;
    void realizeModel(State& state) const;
    void realizeVelocity(const State& state) const;
    void realizeDynamics(const State& state) const;

	// Get the current number of chains.
	int getNumChains(const State& s) const {
		assert(dvNumChains.isValid());
		return Value<int>::downcast(getForceSubsystem().getDiscreteVariable(s, dvNumChains));
	}
	int& updNumChains(State& s) const {
		assert(dvNumChains.isValid());
		return Value<int>::updDowncast(getForceSubsystem().updDiscreteVariable(s, dvNumChains));
	}

	int getNumDOFs(const State& s) const;

	// Get the auxiliary continuous state index of the 0'th thermostat state variable z.
	ZIndex getZ0Index(const State& s) const {
		assert(cacheZ0Index.isValid());
		return Value<ZIndex>::downcast(getForceSubsystem().getCacheEntry(s, cacheZ0Index));
	}
	ZIndex& updZ0Index(State& s) const {
		assert(cacheZ0Index.isValid());
		return Value<ZIndex>::updDowncast(getForceSubsystem().updCacheEntry(s, cacheZ0Index));
	}

	// Get the current bath temperature (i.e., the desired temperature).
	Real getBathTemp(const State& s) const {
		assert(dvBathTemp.isValid());
		return Value<Real>::downcast(getForceSubsystem().getDiscreteVariable(s, dvBathTemp));
	}
	Real& updBathTemp(State& s) const {
		assert(dvBathTemp.isValid());
		return Value<Real>::updDowncast(getForceSubsystem().updDiscreteVariable(s, dvBathTemp));
	}

	// Get the current bath temperature (i.e., the desired temperature).
	Real getRelaxationTime(const State& s) const {
		assert(dvRelaxationTime.isValid());
		return Value<Real>::downcast(getForceSubsystem().getDiscreteVariable(s, dvRelaxationTime));
	}
	Real& updRelaxationTime(State& s) const {
		assert(dvRelaxationTime.isValid());
		return Value<Real>::updDowncast(getForceSubsystem().updDiscreteVariable(s, dvRelaxationTime));
	}

	// Get the calculated momentum M*u (after Stage::Velocity).
	const Vector& getMomentum(const State& s) const {
		assert(cacheMomentumIndex.isValid());
		return Value<Vector>::downcast(getForceSubsystem().getCacheEntry(s, cacheMomentumIndex));
	}
	Vector& updMomentum(const State& s) const {
		assert(cacheMomentumIndex.isValid());
		return Value<Vector>::updDowncast(getForceSubsystem().updCacheEntry(s, cacheMomentumIndex));
	}

	// Get the calculated system kinetic energy ~u*M*u/2 (after Stage::Velocity).
	const Real& getKE(const State& s) const {
		assert(cacheKEIndex.isValid());
		return Value<Real>::downcast(getForceSubsystem().getCacheEntry(s, cacheKEIndex));
	}
	Real& updKE(const State& s) const {
		assert(cacheKEIndex.isValid());
		return Value<Real>::updDowncast(getForceSubsystem().updCacheEntry(s, cacheKEIndex));
	}

	// Get the current value of one of the thermostat state variables.
	Real getZ(const State& s, int i) const {
		assert(0 <= i && i < 2*getNumChains(s));
		const ZIndex z0 = getZ0Index(s);
		return getForceSubsystem().getZ(s)[z0+i];
	}
	Real& updZ(State& s, int i) const {
		assert(0 <= i && i < 2*getNumChains(s));
		const ZIndex z0 = getZ0Index(s);
		return getForceSubsystem().updZ(s)[z0+i];
	}
	Real getZDot(const State& s, int i) const {
		assert(0 <= i && i < 2*getNumChains(s));
		const ZIndex z0 = getZ0Index(s);
		return getForceSubsystem().getZDot(s)[z0+i];
	}
	// State is const here because ZDot is a cache entry.
	Real& updZDot(const State& s, int i) const {
		assert(0 <= i && i < 2*getNumChains(s));
		const ZIndex z0 = getZ0Index(s);
		return getForceSubsystem().updZDot(s)[z0+i];
	}

	static const int DefaultDefaultNumChains = 3;
private:
    const SimbodyMatterSubsystem& matter;
    const Real  Kb;		// Boltzmann's constant in compatible units

	// Topology-stage "state" variables.
	int			defaultNumChains;		// # chains in a new State
	Real		defaultBathTemp;		// bath temperature
	Real		defaultRelaxationTime;	// relaxation time

	// These indices are Topology-stage "cache" variables.
	DiscreteVariableIndex dvNumChains;		// integer
	DiscreteVariableIndex dvBathTemp;		// Real
	DiscreteVariableIndex dvRelaxationTime;	// Real
	CacheEntryIndex		  cacheZ0Index;		// ZIndex
	CacheEntryIndex		  cacheMomentumIndex;	// M*u
	CacheEntryIndex		  cacheKEIndex;			// ~u*M*u/2

friend class Force::Thermostat;
};

class Force::UniformGravityImpl : public ForceImpl {
public:
    UniformGravityImpl(const SimbodyMatterSubsystem& matter, const Vec3& g, Real zeroHeight);
    UniformGravityImpl* clone() const {
        return new UniformGravityImpl(*this);
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
    Vec3 getGravity() const {
        return g;
    }
    void setGravity(const Vec3& gravity) {
        g = gravity;
        invalidateTopologyCache();
    }
    Real getZeroHeight() const {
        return zeroHeight;
    }
    void setZeroHeight(Real height) {
        zeroHeight = height;
        invalidateTopologyCache();
    }
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
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
    ~CustomImpl() {
        delete implementation;
    }
    const Force::Custom::Implementation& getImplementation() const {
        return *implementation;
    }
    Force::Custom::Implementation& updImplementation() {
        return *implementation;
    }
protected:
    void realizeTopology(State& state) const {
        implementation->realizeTopology(state);
    }
    void realizeModel(State& state) const {
        implementation->realizeModel(state);
    }
    void realizeInstance(const State& state) const {
        implementation->realizeInstance(state);
    }
    void realizeTime(const State& state) const {
        implementation->realizeTime(state);
    }
    void realizePosition(const State& state) const {
        implementation->realizePosition(state);
    }
    void realizeVelocity(const State& state) const {
        implementation->realizeVelocity(state);
    }
    void realizeDynamics(const State& state) const {
        implementation->realizeDynamics(state);
    }
    void realizeAcceleration(const State& state) const {
        implementation->realizeAcceleration(state);
    }
    void realizeReport(const State& state) const {
        implementation->realizeReport(state);
    }
private:
    Force::Custom::Implementation* implementation;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_IMPL_H_
