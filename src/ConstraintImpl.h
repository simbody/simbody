#ifndef SimTK_SIMBODY_CONSTRAINT_IMPL_H_
#define SimTK_SIMBODY_CONSTRAINT_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-9 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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

/**@file
 *
 * Private implementation of Constraint and its included subclasses which
 * represent the built-in constraint types.
 */

#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/SimbodyMatterSubtree.h"

#include <map>
#include <utility>  // std::pair
#include <iostream>
using std::cout; using std::endl;

class SimbodyMatterSubsystemRep;
class SBStateDigest;
class SBInstanceCache;
class SBModelCache;
class SBTreePositionCache;
class SBTreeVelocityCache;
class SBTreeAccelerationCache;

namespace SimTK {

class SimbodyMatterSubsystem;
class SimbodyMatterSubtree;
class MobilizedBody;

    /////////////////////
    // CONSTRAINT IMPL //
    /////////////////////

// This is what a Constraint handle points to.
class ConstraintImpl : public PIMPLImplementation<Constraint, ConstraintImpl> {
public:
    ConstraintImpl()
      : myMatterSubsystemRep(0), 
        defaultMp(0), defaultMv(0), defaultMa(0), defaultDisabled(false),
        myAncestorBodyIsNotGround(false)
    {
    }
    virtual ~ConstraintImpl() { }
    virtual ConstraintImpl* clone() const = 0;

    ConstraintImpl(int mp, int mv, int ma)
      : myMatterSubsystemRep(0), 
        defaultMp(mp), defaultMv(mv), defaultMa(ma), defaultDisabled(false),
        myAncestorBodyIsNotGround(false)
    {
    }

    void setDefaultNumConstraintEquations(int mp, int mv, int ma) {
        assert(mp >= 0 && mv >= 0 && ma >= 0);
        invalidateTopologyCache();
        defaultMp = mp;
        defaultMv = mv;
        defaultMa = ma;
    }

    void getDefaultNumConstraintEquations(int& mp, int& mv, int& ma) const {
        mp = defaultMp;
        mv = defaultMv;
        ma = defaultMa;
    }

    void setDisabledByDefault(bool shouldBeDisabled) {
        invalidateTopologyCache();
        defaultDisabled = shouldBeDisabled;
    }

    bool isDisabledByDefault() const {
        return defaultDisabled;
    }

    void setDisabled(State& s, bool shouldBeDisabled) const ;
    bool isDisabled(const State& s) const;

    typedef std::map<MobilizedBodyIndex,ConstrainedBodyIndex>       MobilizedBody2ConstrainedBodyMap;
    typedef std::map<MobilizedBodyIndex,ConstrainedMobilizerIndex>  MobilizedBody2ConstrainedMobilizerMap;

    // Call this during construction phase to add a body to the topological structure of
    // this Constraint. This body's mobilizer's mobilities are *not* part of the constraint; 
    // mobilizers must be added separately.
    ConstrainedBodyIndex addConstrainedBody(const MobilizedBody&);

    // Call this during construction phase to add a mobilizer to the topological structure of
    // this Constraint. All the coordinates q and mobilities u for this mobilizer are added also,
    // but we don't know how many of those there will be until Stage::Model.
    ConstrainedMobilizerIndex addConstrainedMobilizer(const MobilizedBody&);

    MobilizedBodyIndex getMobilizedBodyIndexOfConstrainedBody(ConstrainedBodyIndex c) const {
        assert(0 <= c && c < (int)myConstrainedBodies.size());
        return myConstrainedBodies[c];
    }
    MobilizedBodyIndex getMobilizedBodyIndexOfConstrainedMobilizer(ConstrainedMobilizerIndex c) const {
        assert(0 <= c && c < (int)myConstrainedMobilizers.size());
        return myConstrainedMobilizers[c];
    }

    //TODO: Constraint-local State allocation
    //int allocateDiscreteVariable(State& s, Stage g, AbstractValue* v) const;
    //int allocateCacheEntry(State& s, Stage g, AbstractValue* v) const;

    QIndex getQIndexOfConstrainedQ(const State& s, ConstrainedQIndex cqx) const;
    UIndex getUIndexOfConstrainedU(const State& s, ConstrainedUIndex cqx) const;

    void realizeTopology(State&)       const; // eventually calls realizeTopologyVirtual()
    void realizeModel   (State&)       const; // eventually calls realizeModelVirtual() 
    void realizeInstance(const State&) const; // eventually calls realizeInstanceVirtual() 

    // These are called in loops over all the Constraints from the SimbodyMatterSubsystem
    // realize() methods, which will have already deconstructed the State into
    // an SBStateDigest object.
    void realizeTime    (const SBStateDigest&) const; // eventually calls realizeTimeVirtual() 
    void realizePosition(const SBStateDigest&) const; // eventually calls realizePositionVirtual() 
    void realizeVelocity(const SBStateDigest&) const; // eventually calls realizeVelocityVirtual() 
    void realizeDynamics(const SBStateDigest&) const; // eventually calls realizeDynamicsVirtual() 
    void realizeAcceleration
                        (const SBStateDigest&) const; // eventually calls realizeAccelerationVirtual() 
    void realizeReport  (const State&) const; // eventually calls realizeReportVirtual() 

	// Given a state realized to Position stage, extract the position constraint errors
	// corresponding to this Constraint. The 'mp' argument is for sanity checking -- it
	// is an error if that isn't an exact match for the current number of holonomic
	// constraint equations generated by this Constraint. We expect that perr points
	// to an array of at least mp elements that we can write on.
	void getPositionErrors(const State& s, int mp, Real* perr) const;

	// Given a State realized to Velocity stage, extract the velocity constraint errors
	// corresponding to this Constraint. This includes velocity constraints which were
	// produced by differentiation of holonomic (position) constraints, and nonholonomic
	// constraints which are introduced at the velocity level. The 'mpv' argument is
	// for sanity checking -- it is an error if that isn't an exact match for the
	// current number of holonomic+nonholonomic (mp+mv) constraint equations generated
	// by this Constraint. We expect that pverr points to an array of at least mp+mv
	// elements that we can write on.
	void getVelocityErrors(const State& s, int mpv, Real* pverr) const;

	// Given a State realized to Acceleration stage, extract the accleration constraint errors
	// corresponding to this Constraint. This includes acceleration constraints which were
	// produced by twice differentiation of holonomic (position) constraints, and differentiation
	// of nonholonomic (velocity) constraints, and acceleration-only constraints which are
	// first introduced at the acceleration level. The 'mpva' argument is
	// for sanity checking -- it is an error if that isn't an exact match for the
	// current number of holonomic+nonholonomic+accelerationOnly (mp+mv+ma) constraint
	// equations generated by this Constraint. We expect that pvaerr points to an array
	// of at least mp+mv+ma elements that we can write on.
	void getAccelerationErrors(const State& s, int mpva, Real* pvaerr) const;

	// Given a State realized to Acceleration stage, extract the constraint multipliers lambda
	// corresponding to this constraint. This includes multipliers for all the holonomic,
	// nonholonomic, and acceleration-only constraints (but not quaternion constraints which
	// do not use multipliers). The 'mpva' argument is for sanity checking -- it is an error
	// if that isn't an exact match for the current number (mp+mv+ma) of constraint
	// equations generated by this Constraint. We expect that lambda points to an array
	// of at least mp+mv+ma elements that we can write on.
	void getMultipliers(const State& s, int mpva, Real* lambda) const;

	// Given a State realized to Position stage, and a set of m multipliers lambda, calculate in O(m) time
	// the constraint forces (body forces and torques and mobility forces) which would be generated
	// by those multipliers. You can restrict this to P,V,A subsets setting mp, mv, or ma to zero.
	void calcConstraintForcesFromMultipliers(const State& s, int mp, int mv, int ma, const Real* lambda,
							  Vector_<SpatialVec>& bodyForcesInA,		 // for each constrained body
							  Vector&              mobilityForces) const // for each constrained mobility
	{
		int actual_mp,actual_mv,actual_ma;
		getNumConstraintEquationsInUse(s, actual_mp, actual_mv, actual_ma);

		bodyForcesInA.resize(getNumConstrainedBodies());       bodyForcesInA  = SpatialVec(Vec3(0), Vec3(0));
		mobilityForces.resize(getNumConstrainedU(s));          mobilityForces = 0;

		if (mp) {
			assert(mp == actual_mp);
			applyPositionConstraintForces(s, mp, &lambda[0], bodyForcesInA, mobilityForces);
		}
		if (mv) {
			assert(mv == actual_mv);
			applyVelocityConstraintForces(s, mv, &lambda[mp], bodyForcesInA, mobilityForces);
		}
		if (ma) {
			assert(ma == actual_ma);
			applyAccelerationConstraintForces(s, ma, &lambda[mp+mv], bodyForcesInA, mobilityForces);
		}
	}

	// Given a State realized to Position stage, and a set of forces applied to the constrained
	// bodies and their mobilizers, convert these to an equivalent set of n generalized forces
	// applied at each of the participating mobilities, in O(n) time.
	void convertConstraintForcesToGeneralizedForces(const State& s,
	    const Vector_<SpatialVec>& bodyForcesInA, const Vector& mobilityForces,
		Vector& generalizedForces) const
	{
		// TODO
		assert(!"convertConstraintForcesToGeneralizedForces: not implemented yet");
	}

	// Calculate f = ~G*lambda in O(n+m) time. ~G=[~P ~V ~A] and you can work with any
	// subblock or combination by setting some of mp,mv,ma to zero. If nonzero they have
	// to match the actual number of holonomic, nonholonomic, and acceleration-only constraints.
	// Vector lambda (typically Lagrange multipliers but not necessarily) is segmented 
	// lambda=[mp|mv|ma] where some of the segments can be empty.
	void calcGTransposeLambda(const State& s, int mp, int mv, int ma, const Real* lambda,
							  Vector& f) const
	{
		Vector_<SpatialVec> bodyForcesInA;
		Vector              mobilityForces;
		calcConstraintForcesFromMultipliers(s, mp,mv,ma, lambda, bodyForcesInA, mobilityForces);
		convertConstraintForcesToGeneralizedForces(s, bodyForcesInA, mobilityForces, f);
	}

    // Find the indicated cache in the passed-in State. This requires that the 
    // cache entry has already been marked valid.

    const SBInstanceCache&		    getInstanceCache(const State&) const;
    const SBModelCache&             getModelCache(const State&) const;
    const SBTreePositionCache&		getTreePositionCache(const State&) const;
    const SBTreeVelocityCache&		getTreeVelocityCache(const State&) const;
    const SBTreeAccelerationCache&	getTreeAccelerationCache(const State&) const;


        // Methods for use with ConstrainedMobilizers.

    Real getOneQ(const State&, ConstrainedMobilizerIndex, MobilizerQIndex) const;
    Real getOneU(const State&, ConstrainedMobilizerIndex, MobilizerUIndex) const;

    //TODO: get rid of these annoying realization flags. State should have a "mark as valid"
    //feature for these built in cache entries rather than depending on stage.
    Real getOneQDot   (const State&, ConstrainedMobilizerIndex, MobilizerQIndex, bool realizingVelocity=false) const;
    Real getOneQDotDot(const State&, ConstrainedMobilizerIndex, MobilizerQIndex, bool realizingAcceleration=false) const;
    Real getOneUDot   (const State&, ConstrainedMobilizerIndex, MobilizerUIndex, bool realizingAcceleration=false) const;

    // Apply a generalized (mobility) force to a particular mobility of the given constrained body B,
    // adding it in to the appropriate slot of the mobilityForces vector.
    void addInOneMobilityForce(const State& s, ConstrainedMobilizerIndex M, MobilizerUIndex which,
                               Real f, Vector& mobilityForces) const 
    { 
        assert(mobilityForces.size() == getNumConstrainedU(s));
        assert(0 <= which && which < getNumConstrainedU(s, M));
        mobilityForces[getConstrainedUIndex(s,M,which)] += f;
    }

        // Methods for use with ConstrainedBodies.

    // These are used to retrieve the indicated values from the State cache, with all values
    // measured and expressed in the Ancestor (A) frame.
    const Transform&  getBodyTransform   (const State& s, ConstrainedBodyIndex B)     const; // X_AB
    const SpatialVec& getBodyVelocity    (const State& s, ConstrainedBodyIndex B)     const; // V_AB
    const SpatialVec& getBodyAcceleration(const State& s, ConstrainedBodyIndex B) const; // A_AB

    // Extract just the rotational quantities from the spatial quantities above.
    const Rotation& getBodyRotation       (const State& s, ConstrainedBodyIndex B)     const // R_AB
       {return getBodyTransform(s,B).R();}
    const Vec3& getBodyAngularVelocity    (const State& s, ConstrainedBodyIndex B)     const // w_AB
       {return getBodyVelocity(s,B)[0];}
    const Vec3& getBodyAngularAcceleration(const State& s, ConstrainedBodyIndex B) const // b_AB
       {return getBodyAcceleration(s,B)[0];}

    // Extract just the translational (linear) quantities from the spatial quantities above.
    const Vec3& getBodyOriginLocation    (const State& s, ConstrainedBodyIndex B)     const // p_AB
       {return getBodyTransform(s,B).p();}  
    const Vec3& getBodyOriginVelocity    (const State& s, ConstrainedBodyIndex B)     const // v_AB
       {return getBodyVelocity(s,B)[1];}     
    const Vec3& getBodyOriginAcceleration(const State& s, ConstrainedBodyIndex B) const // a_AB
       {return getBodyAcceleration(s,B)[1];} 

    Vec3 findStationLocation(const State& s, ConstrainedBodyIndex B, const Vec3& p_B) const {
        return getBodyTransform(s,B) * p_B; // re-measure and re-express
    }
    Vec3 findStationVelocity(const State& s, ConstrainedBodyIndex B, const Vec3& p_B) const {
        const Vec3        p_A  = getBodyRotation(s,B) * p_B; // rexpressed but not shifted
        const SpatialVec& V_AB = getBodyVelocity(s,B);
        return V_AB[1] + (V_AB[0] % p_A);
    }
    Vec3 findStationAcceleration(const State& s, ConstrainedBodyIndex B, const Vec3& p_B) const {
        const Vec3        p_A  = getBodyRotation(s,B) * p_B; // rexpressed but not shifted
        const Vec3&       w_AB = getBodyAngularVelocity(s,B);
        const SpatialVec& A_AB = getBodyAcceleration(s,B);
        const Vec3 a_A = A_AB[1] + (A_AB[0] % p_A) + w_AB % (w_AB % p_A); // careful: cross product is not associative
        return a_A;
    }

    // Apply an Ancestor-frame force to a B-frame station, updating the appropriate bodyForces entry.
    void addInStationForce(const State& s, ConstrainedBodyIndex B, const Vec3& p_B, 
                           const Vec3& forceInA, Vector_<SpatialVec>& bodyForcesInA) const 
    {
        assert(bodyForcesInA.size() == getNumConstrainedBodies());
        const Rotation& R_AB = getBodyRotation(s,B);
        bodyForcesInA[B] += SpatialVec((R_AB*p_B) % forceInA, forceInA); // rXf, f
    }

    // Apply an Ancestor-frame torque to body B, updating the appropriate bodyForces entry.
    void addInBodyTorque(const State& s, ConstrainedBodyIndex B,
                         const Vec3& torqueInA, Vector_<SpatialVec>& bodyForcesInA) const 
    {
        assert(bodyForcesInA.size() == getNumConstrainedBodies());
        bodyForcesInA[B][0] += torqueInA; // no force
    }


    // After realizeTopology() we can look at the values of modeling variables in the State.
    // A Constraint is free to use those in determining how many constraint equations of each
    // type to generate. The default implementation here doesn't look at the state but instead
    // returns the default numbers of equations supplied when the Constraint was constructed.
    void calcNumConstraintEquationsInUse(const State& s, int& mp, int& mv, int& ma) const {
        calcNumConstraintEquationsInUseVirtual(s,mp,mv,ma);
    }
    virtual void calcNumConstraintEquationsInUseVirtual(const State&, int& mp, int& mv, int& ma) const {
        mp = defaultMp; mv = defaultMv; ma = defaultMa;
    }

    void realizePositionErrors(const State& s, int mp,  Real* perr) const {
        realizePositionErrorsVirtual(s,mp,perr);
    }
    void realizePositionDotErrors(const State& s, int mp,  Real* pverr) const {
        realizePositionDotErrorsVirtual(s,mp,pverr);
    }
    void realizePositionDotDotErrors(const State& s, int mp,  Real* paerr) const {
        realizePositionDotDotErrorsVirtual(s,mp,paerr);
    }
    void applyPositionConstraintForces
       (const State& s, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const
    {
        applyPositionConstraintForcesVirtual(s,mp,multipliers,bodyForces,mobilityForces);
    }

    void realizeVelocityErrors(const State& s, int mv,  Real* verr) const {
        realizeVelocityErrorsVirtual(s,mv,verr);
    }
    void realizeVelocityDotErrors(const State& s, int mv,  Real* vaerr) const {
        realizeVelocityDotErrorsVirtual(s,mv,vaerr);
    }
    void applyVelocityConstraintForces
       (const State& s, int mv, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const
    {
        applyVelocityConstraintForcesVirtual(s,mv,multipliers,bodyForces,mobilityForces);
    }

    void realizeAccelerationErrors(const State& s, int ma,  Real* aerr) const {
        realizeAccelerationErrorsVirtual(s,ma,aerr);
    }
    void applyAccelerationConstraintForces
       (const State& s, int ma, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const
    {
        applyAccelerationConstraintForcesVirtual(s,ma,multipliers,bodyForces,mobilityForces);
    }

    void calcDecorativeGeometryAndAppend
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
    {
        // Let the individual constraint deal with any complicated stuff.
        calcDecorativeGeometryAndAppendVirtual(s,stage,geom);
    }

    virtual void realizeTopologyVirtual     (State&)        const {}
    virtual void realizeModelVirtual        (State&)        const {}
    virtual void realizeInstanceVirtual     (const State&)  const {}
    virtual void realizeTimeVirtual         (const State&)  const {}
    virtual void realizePositionVirtual     (const State&)  const {}
    virtual void realizeVelocityVirtual     (const State&)  const {}
    virtual void realizeDynamicsVirtual     (const State&)  const {}
    virtual void realizeAccelerationVirtual (const State&)  const {}
    virtual void realizeReportVirtual       (const State&)  const {}

    // These must be defined if there are any position (holonomic) constraints defined.
    virtual void realizePositionErrorsVirtual      (const State&, int mp,  Real* perr)  const;
    virtual void realizePositionDotErrorsVirtual   (const State&, int mp,  Real* pverr) const;
    virtual void realizePositionDotDotErrorsVirtual(const State&, int mp,  Real* paerr) const;
    virtual void applyPositionConstraintForcesVirtual
       (const State&, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const;

    // These must be defined if there are any velocity (nonholonomic) constraints defined.
    virtual void realizeVelocityErrorsVirtual   (const State&, int mv,  Real* verr)  const;
    virtual void realizeVelocityDotErrorsVirtual(const State&, int mv,  Real* vaerr) const;
    virtual void applyVelocityConstraintForcesVirtual
       (const State&, int mv, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const;

    // These must be defined if there are any acceleration-only constraints defined.
    virtual void realizeAccelerationErrorsVirtual(const State&, int ma,  Real* aerr) const;
    virtual void applyAccelerationConstraintForcesVirtual
       (const State&, int ma, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const;

    virtual void calcDecorativeGeometryAndAppendVirtual
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const {}


    void invalidateTopologyCache() const;
    bool subsystemTopologyHasBeenRealized() const;

    void setMyMatterSubsystem(SimbodyMatterSubsystem& matter,
                              ConstraintIndex id);

    const SimbodyMatterSubsystem& getMyMatterSubsystem() const;

    bool isInSubsystem() const {
        return myMatterSubsystemRep != 0;
    }

    // Is the supplied body in the same subsystem as this Constraint? (Returns false also if
    // either the Constraint or the MobilizedBody is not in a subsystem.)
    bool isInSameSubsystem(const MobilizedBody& body) const;

    int getNumConstrainedBodies() const {
        SimTK_ASSERT(subsystemTopologyHasBeenRealized(),
            "Number of constrained bodies is not available until Topology stage has been realized.");
        return (int)myConstrainedBodies.size();
    }
    int getNumConstrainedMobilizers() const {
        SimTK_ASSERT(subsystemTopologyHasBeenRealized(),
            "Number of constrained mobilizers is not available until Topology stage has been realized.");
        return (int)myConstrainedMobilizers.size();
    }

    const MobilizedBody& getMobilizedBodyFromConstrainedMobilizer(ConstrainedMobilizerIndex) const;
    const MobilizedBody& getMobilizedBodyFromConstrainedBody(ConstrainedBodyIndex) const;

    // Don't call this unless there is at least one Constrained Body.
    const MobilizedBody& getAncestorMobilizedBody() const;

	// Find out how many holonomic (position), nonholonomic (velocity),
	// and acceleration-only constraint equations are generated by this Constraint
    // as currently modeled. State must be realized to Stage::Model.
	void getNumConstraintEquationsInUse(const State&, int& mHolo, int& mNonholo, int& mAccOnly) const;

    // Find the first assigned slots for these constraint equations in the containing
    // SimbodyMatterSubsystem's QErr, UErr, and UDotErr/Multiplier arrays. There will be
    // (offset,length) slots at:
    //    (holo0,                                   mHolo)
    //    (totalNumHolo+nonholo0,                   mNonholo)
    //    (totalNumHolo+totalNumNonholo + accOnly0, mAccOnly)
    // Returns -1 if the Constraint has no constraint equations in the indicated category.
	// State must be realized to Stage::Model.
    void getConstraintEquationSlots(const State&, int& holo0, int& nonholo0, int& accOnly0) const;

    int getNumConstrainedQ(const State&) const;
    int getNumConstrainedU(const State&) const;
    int getNumConstrainedQ(const State&, ConstrainedMobilizerIndex) const;
    int getNumConstrainedU(const State&, ConstrainedMobilizerIndex) const;
    ConstrainedQIndex getConstrainedQIndex
       (const State&, ConstrainedMobilizerIndex, MobilizerQIndex which) const;
    ConstrainedUIndex getConstrainedUIndex
       (const State&, ConstrainedMobilizerIndex, MobilizerUIndex which) const;

    const SimbodyMatterSubsystemRep& getMyMatterSubsystemRep() const {
        SimTK_ASSERT(myMatterSubsystemRep,
            "Operation illegal on a Constraint that is not in a Subsystem.");
        return *myMatterSubsystemRep;
    }
    SimbodyMatterSubsystemRep& updMyMatterSubsystemRep() {
        SimTK_ASSERT(myMatterSubsystemRep,
            "Operation illegal on a Constraint that is not in a Subsystem.");
        return *myMatterSubsystemRep;
    }

    // Calculate the transform X_AB of each ConstrainedBody in its Ancestor frame, 
    // provided Ancestor!=Ground and A!=B. We expect a StateDigest of a State
    // that has been realized through Time stage, and a TreePositionCache in
    // which the mobilizer- and ground-frame position kinematics results have
    // already been calculated. We then fill in the missing ancestor-frame
    // results back into that same TreePositionCache. The TreePositionCache
    // in the StateDigest is ignored; it may or may not be the same as the
    // second argument depending on whether this is part of a realization or
    // part of an operator.
    void calcConstrainedBodyTransformInAncestor(const SBStateDigest&,   // in only
                                                SBTreePositionCache&    // in/out
                                                ) const;
    // Similarly we calculate V_AB and A_AB during realizeVelocity() and realizeAcceleration().
    // Similarly, here we expect a StateDigest realized through Position stage,
    // and a partly-filled-in VelocityCache where we'll put V_AB for the 
    // ConstrainedBodies.
    void calcConstrainedBodyVelocityInAncestor(const SBStateDigest&,    // in only
                                               SBTreeVelocityCache&     // in/out
                                               ) const;

    // Similarly, here we expect a StateDigest realized through Dynamics stage,
    // and a partly-filled-in AccelerationCache where we'll put A_AB for the 
    // ConstrainedBodies.
    void calcConstrainedBodyAccelerationInAncestor(const SBStateDigest&,    // in only
                                                   SBTreeAccelerationCache& // in/out
                                                   ) const;

private:
    friend class Constraint;

        // TOPOLOGY "STATE"

    // These data members are filled in once the Constraint is added to
    // a MatterSubsystem.
    SimbodyMatterSubsystemRep* myMatterSubsystemRep;
    ConstraintIndex            myConstraintIndex; // id within the matter subsystem

    // We'll keep the constrained bodies and constrained mobilizers each in two maps: 
    // one maps MobilizedBodyIndex->ConstrainedBody[Mobilizer]Index (O(log n) to look
    // up), and the other maps ConstrainedBody[Mobilizer]Index->MobilizedBodyIndex
    // (randomly addressable in constant time).
    MobilizedBody2ConstrainedBodyMap        myMobilizedBody2ConstrainedBodyMap;
    MobilizedBody2ConstrainedMobilizerMap   myMobilizedBody2ConstrainedMobilizerMap;

    Array_<MobilizedBodyIndex> myConstrainedBodies;     // index with ConstrainedBodyIndex
    Array_<MobilizedBodyIndex> myConstrainedMobilizers; // index with ConstrainedMobilizerIndex


    // These are the defaults for the number of position (holonomic) constraint equations,
    // the number of velocity (nonholonomic) constraint equations, and the number of
    // acceleration-only constraint equations.
    int defaultMp, defaultMv, defaultMa;

    // This says whether the Model-stage "disabled" flag for this Constraint should be initially
    // on or off. Most constraints are enabled by default.
    bool defaultDisabled;

        // TOPOLOGY "CACHE"

    // When topology is realized we study the constrained bodies to identify the
    // subtree of mobilized bodies which may be kinematically involved in satisfaction
    // of this Constraint. This requires finding the outmost common ancestor of 
    // the constrained bodies. All mobilized bodies on the paths inward from the
    // constrained bodies to the ancestor are included; nothing outboard of the
    // constrained bodies is included; and the ancestor is treated as ground so
    // that its mobilities are *not* included. The Ancestor may be one of the
    // Constrained Bodies or may be distinct.
    mutable SimbodyMatterSubtree mySubtree;

    // This is true only when (1) there are Constrained Bodies and (2) their Ancestor is
    // some MobilizedBody other than Ground.
    mutable bool myAncestorBodyIsNotGround;

    // When the Ancestor is not Ground, each of the Constrained Bodies (except
    // the Ancestor) needs some additional precalculated data in the State cache
    // so is assigned a MatterSubsystem-global AncestorConstrainedBodyPool slot.
    // If Ancestor isn't Ground there is an entry here for each ConstrainedBody
    // (including Ancestor if it is one), but the index is invalid for Ancestor's
    // entry. When Ancestor is Ground we don't allocate this array.
    mutable Array_<AncestorConstrainedBodyPoolIndex> myPoolIndex; // index with ConstrainedBodyIndex
};

    // ROD

class Constraint::RodImpl : public ConstraintImpl {
public:
    RodImpl() 
      : ConstraintImpl(1,0,0), defaultPoint1(0), defaultPoint2(0), defaultRodLength(1),
        pointRadius(-1) // this means "use default point radius"
    { 
        // Rod constructor sets all the data members here directly
    }
    RodImpl* clone() const { return new RodImpl(*this); }

    // Draw some end points and a rubber band line.
    void calcDecorativeGeometryAndAppendVirtual
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const;

    void setPointDisplayRadius(Real r) {
        // r == 0 means don't display point, r < 0 means use default which is some fraction of rod length
        invalidateTopologyCache();
        pointRadius= r > 0 ? r : 0;
    }
    Real getPointDisplayRadius() const {return pointRadius;}

    // Implementation of virtuals required for holonomic constraints.

	// TODO: this is badly scaled; consider using length instead of length^2
    // perr = (p^2 - d^2)/2
    void realizePositionErrorsVirtual(const State& s, int mp,  Real* perr) const {
        assert(mp==1 && perr);
        const Vec3 p1 = findStationLocation(s, B1, defaultPoint1); // meas from & expr in ancestor
        const Vec3 p2 = findStationLocation(s, B2, defaultPoint2);
        const Vec3 p = p2 - p1;
        //TODO: save p in state

        *perr = (dot(p, p) - square(defaultRodLength)) / 2;
    }

    // pverr = d/dt perr = pdot*p = v*p, where v=v2-v1 is relative velocity
    void realizePositionDotErrorsVirtual(const State& s, int mp,  Real* pverr) const {
        assert(mp==1 && pverr);
        //TODO: should be able to get p from State
        const Vec3 p1 = findStationLocation(s, B1, defaultPoint1); // meas from & expr in ancestor
        const Vec3 p2 = findStationLocation(s, B2, defaultPoint2);
        const Vec3 p = p2 - p1;

        const Vec3 v1 = findStationVelocity(s, B1, defaultPoint1); // meas & expr in ancestor
        const Vec3 v2 = findStationVelocity(s, B2, defaultPoint2);
        const Vec3 v = v2 - v1;
        *pverr = dot(v, p);
    }

    // paerr = d/dt verr = vdot*p + v*pdot =a*p+v*v, where a=a2-a1 is relative acceleration
    void realizePositionDotDotErrorsVirtual(const State& s, int mp,  Real* paerr) const {
        assert(mp==1 && paerr);
        //TODO: should be able to get p and v from State
        const Vec3 p1 = findStationLocation(s, B1, defaultPoint1); // meas from & expr in ancestor
        const Vec3 p2 = findStationLocation(s, B2, defaultPoint2);
        const Vec3 p = p2 - p1;
        const Vec3 v1 = findStationVelocity(s, B1, defaultPoint1); // meas & expr in ancestor
        const Vec3 v2 = findStationVelocity(s, B2, defaultPoint2);
        const Vec3 v = v2 - v1;

        const Vec3 a1 = findStationAcceleration(s, B1, defaultPoint1); // meas & expr in ancestor
        const Vec3 a2 = findStationAcceleration(s, B2, defaultPoint2);
        const Vec3 a = a2 - a1;

        *paerr = dot(a, p) + dot(v, v);
    }

    // Write this routine by inspection of the pverr routine, looking for terms involving
    // velocity. On point2 we see v2*p, on point1 we see -v1*p, so forces are m*p and -m*p,
    // respectively.
    void applyPositionConstraintForcesVirtual
       (const State& s, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForcesInA,
        Vector&              mobilityForces) const
    {
        assert(mp==1 && multipliers);
        const Real lambda = *multipliers;
        //TODO: should be able to get p from State
        const Vec3 p1 = findStationLocation(s, B1, defaultPoint1); // meas from & expr in ancestor
        const Vec3 p2 = findStationLocation(s, B2, defaultPoint2);
        const Vec3 p = p2 - p1;

        const Vec3 f2 = lambda * p;

        // The forces on either point have the same line of action because they are aligned
        // with the vector between the points. Applying the forces to any point along the line
        // would have the same effect (e.g., same point in space on both bodies) so this is
        // the same as an equal and opposite force applied to the same point and this constraint
        // will do no work even if the position or velocity constraints are not satisfied.
        addInStationForce(s, B2, defaultPoint2,  f2, bodyForcesInA);
        addInStationForce(s, B1, defaultPoint1, -f2, bodyForcesInA);
    }

    SimTK_DOWNCAST(RodImpl, ConstraintImpl);
private:
    friend class Constraint::Rod;

    ConstrainedBodyIndex B1, B2;

    Vec3            defaultPoint1; // on body 1, exp. in B1 frame
    Vec3            defaultPoint2; // on body 2, exp. in B2 frame
    Real            defaultRodLength;

    // This is just for visualization
    Real pointRadius;
};

    // POINT IN PLANE

class Constraint::PointInPlaneImpl : public ConstraintImpl {
public:
    PointInPlaneImpl()
      : ConstraintImpl(1,0,0), defaultPlaneNormal(), defaultPlaneHeight(0), defaultFollowerPoint(0),
        planeHalfWidth(1), pointRadius(0.05) 
    { }
    PointInPlaneImpl* clone() const { return new PointInPlaneImpl(*this); }

    void calcDecorativeGeometryAndAppendVirtual
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const;

    void setPlaneDisplayHalfWidth(Real h) {
        // h <= 0 means don't display plane
        invalidateTopologyCache();
        planeHalfWidth = h > 0 ? h : 0;
    }
    Real getPlaneDisplayHalfWidth() const {return planeHalfWidth;}

    void setPointDisplayRadius(Real r) {
        // r <= 0 means don't display point
        invalidateTopologyCache();
        pointRadius= r > 0 ? r : 0;
    }
    Real getPointDisplayRadius() const {return pointRadius;}

    // Implementation of virtuals required for holonomic constraints.

    // We have a point-in-plane connection between base body B, on which the plane is fixed,
	// and follower body F, on which the follower point S is fixed. All forces will be applied
	// at point S and the coincident material point C on B which is instantaneously at the same spatial
	// location as S. Then n is the plane normal (a constant unit vector in B), h is the
	// plane height measured from the B origin along n (a scalar constant).Point C's location in B
	// is given by the vector p_BC from B's origin to the current location of S, and expressed
	// in B. That vector expressed in A is p_BC_A (= p_AS-p_AB). We will express in the A frame but
    // differentiate in the B frame.
    //
    // Derivation:
    //   (1) Note that to take a derivative d/dt_B in a moving frame B, we can take the derivative d/dt_A
    //       and then add in the contribution d_A/dt_B from A's rotation in B (which is the angular
    //       velocity of A in B, w_BA=-w_AB).
    //   (2) p_CS = p_AS-p_AC = 0 by construction of C, but its derivative in A, 
    //       v_CS_A = d/dt_A p_CS != 0.
    //
    //    perr = p_CS * n + constant 
    //         = constant  (because P_CS==0 by construction)
    //
    //    verr = d/dt_B perr = d/dt_A perr + d_A/dt_B perr
    //         = [v_CS_A*n + p_CS * (w_AB X n)] + [(w_BA X p_CS) * n + p_CS * (w_BA X n)]
    //         = v_CS_A*n + p_CS * (w_AB X n) (because terms in 2nd [] cancel)
    //         = v_CS_A * n  (because p_CS==0 by construction)
    //
    //    aerr = d/dt_B verr = d/dt_A verr + d_A/dt_B verr
    //         = [a_CS_A*n + v_CS_A*(w_AB X n) + v_CS_A*(w_AB X n) + p_CS*(2 w_AB X(w_AB X n))]
    //           + [w_BAXv_CS_A*n + v_CS_A*w_BAXn]
    //         = (a_CS_A - 2 w_AB X v_CS_A) * n   (2nd bracket cancels, and p_CS==0)
    //
    // (The constant in perr is set so that S starts at the same height h as the plane.)
    //  
    // Then, from examination of verr noting that v_CS_A=v_AS-v_AC:
    //       ~v_AS*n                  (body F at point S) 
	//     - ~v_AC*n                  (body B at point C)
    // so we apply a forces lambda*n to F at S, -lambda*n to B at C.
    //
	//    --------------------------------
    //    perr = ~p_BS*n - h
	//    --------------------------------
    void realizePositionErrorsVirtual(const State& s, int mp,  Real* perr) const {
        assert(mp==1 && perr);

        const Transform& X_AB = getBodyTransform(s, planeBody);
        const Vec3       p_AS = findStationLocation(s, followerBody, defaultFollowerPoint);
        const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin and reexpress in B;
                                              // C is material point of B coincident with S

        // We'll calculate this scalar using B-frame vectors, but any frame would have done.
        *perr = dot(p_BC, defaultPlaneNormal) - defaultPlaneHeight;
    }

	//    --------------------------------
    //    verr = ~v_CS_A*n
	//    --------------------------------
    void realizePositionDotErrorsVirtual(const State& s, int mp,  Real* pverr) const {
        assert(mp==1 && pverr);
        //TODO: should be able to get p info from State
        const Transform& X_AB = getBodyTransform(s, planeBody);
        const Vec3       p_AS = findStationLocation(s, followerBody, defaultFollowerPoint);
        const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin and reexpress in B;
                                              // C is material point of B coincident with S
        const UnitVec3   n_A  = X_AB.R() * defaultPlaneNormal;

        const Vec3       v_AS = findStationVelocity(s, followerBody, defaultFollowerPoint);
        const Vec3       v_AC = findStationVelocity(s, planeBody, p_BC);

        // Calculate this scalar using A-frame vectors.
        *pverr = dot( v_AS-v_AC, n_A );
    }

	//    -------------------------------------
    //    aerr = ~(a_CS_A - 2 w_AB X v_CS_A) * n
	//    -------------------------------------
    void realizePositionDotDotErrorsVirtual(const State& s, int mp,  Real* paerr) const {
        assert(mp==1 && paerr);
        //TODO: should be able to get p and v info from State
        const Transform& X_AB = getBodyTransform(s, planeBody);
        const Vec3       p_AS = findStationLocation(s, followerBody, defaultFollowerPoint);
        const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin and reexpress in B;
                                              // C is material point of B coincident with S
        const UnitVec3   n_A  = X_AB.R() * defaultPlaneNormal;

        const Vec3&      w_AB = getBodyAngularVelocity(s, planeBody);
        const Vec3       v_AS = findStationVelocity(s, followerBody, defaultFollowerPoint);
        const Vec3       v_AC = findStationVelocity(s, planeBody, p_BC);

        const Vec3       a_AS = findStationAcceleration(s, followerBody, defaultFollowerPoint);;
        const Vec3       a_AC = findStationAcceleration(s, planeBody, p_BC);

        *paerr = dot( (a_AS-a_AC) - 2*w_AB % (v_AS-v_AC), n_A );
    }

	// apply f=lambda*n to the follower point S of body F,
	//       -f         to point C (coincident point) of body B
    void applyPositionConstraintForcesVirtual
       (const State& s, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForcesInA,
        Vector&              mobilityForces) const
    {
        assert(mp==1 && multipliers);
        const Real lambda = *multipliers;

        //TODO: should be able to get p info from State
        const Transform& X_AB    = getBodyTransform(s, planeBody);
        const Vec3&      p_FS    = defaultFollowerPoint; // measured & expressed in F
        const Vec3       p_AS    = findStationLocation(s, followerBody, defaultFollowerPoint);
        const Vec3       p_BC    = ~X_AB * p_AS;         // measured & expressed in B
		const Vec3       force_A = X_AB.R()*(lambda*defaultPlaneNormal);

        addInStationForce(s, followerBody, p_FS,  force_A, bodyForcesInA);
        addInStationForce(s, planeBody,    p_BC, -force_A, bodyForcesInA);
    }

    SimTK_DOWNCAST(PointInPlaneImpl, ConstraintImpl);
private:
    friend class Constraint::PointInPlane;

    ConstrainedBodyIndex planeBody;    // B1
    ConstrainedBodyIndex followerBody; // B2

    UnitVec3          defaultPlaneNormal;   // on body 1, exp. in B1 frame
    Real              defaultPlaneHeight;
    Vec3              defaultFollowerPoint; // on body 2, exp. in B2 frame

    // These are just for visualization
    Real planeHalfWidth;
    Real pointRadius;
};

    // POINT ON LINE

class Constraint::PointOnLineImpl : public ConstraintImpl {
public:
    PointOnLineImpl()
      : ConstraintImpl(2,0,0), defaultLineDirection(), defaultPointOnLine(), defaultFollowerPoint(0),
        lineHalfLength(1), pointRadius(0.05) 
    { }
    PointOnLineImpl* clone() const { return new PointOnLineImpl(*this); }

    void calcDecorativeGeometryAndAppendVirtual
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const;

    void setLineDisplayHalfLength(Real h) {
        // h <= 0 means don't display line
        invalidateTopologyCache();
        lineHalfLength = h > 0 ? h : 0;
    }
    Real getLineDisplayHalfLength() const {return lineHalfLength;}

    void setPointDisplayRadius(Real r) {
        // r <= 0 means don't display point
        invalidateTopologyCache();
        pointRadius= r > 0 ? r : 0;
    }
    Real getPointDisplayRadius() const {return pointRadius;}

    // Implementation of ContraintRep virtuals
    void realizeTopologyVirtual(State& s) const {
        x = defaultLineDirection.perp(); // x and y are mutable
        y = UnitVec3(defaultLineDirection % x);
    }

    // Implementation of virtuals required for holonomic constraints.

    // We have a point-on-line connection between base body B, on which the line is fixed,
	// and follower body F, on which the follower point S is fixed. All forces will be applied
	// at point S and the coincident material point C on B which is instantaneously at the same spatial
	// location as S. Then z is a unit vector in the direction of the line, and P is a point fixed
    // to B that the line passes through. We will enforce this using two point-on-plane constraints,
    // where the intersection of the two planes is the line. For that we need two plane normals
    // perpendicular to z. We'll use an arbitrary perpendicular x, then use y=z X x as the
    // other perpendicular. This establishes a right handed coordinate system where the line
    // is along the z axis, and we'll apply constraint forces in the x-y plane.
    //
    // See the point-in-plane constraint for details; here we're just picking x and y as
    // plane normals.

	//    --------------------------------
    //    perr = ~(p_BS-p_BP) * x
    //           ~(p_BS-p_BP) * y
	//    --------------------------------
    void realizePositionErrorsVirtual(const State& s, int mp,  Real* perr) const {
        assert(mp==2 && perr);

        const Transform& X_AB = getBodyTransform(s, lineBody);
        const Vec3       p_AS = findStationLocation(s, followerBody, defaultFollowerPoint);
        const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin and reexpress in B;
                                              // C is material point of B coincident with S
        const Vec3       p_PC = p_BC - defaultPointOnLine;

        // We'll calculate these two scalars using B-frame vectors, but any frame would have done.
        Vec2::updAs(perr) = Vec2(~p_PC*x, ~p_PC*y);
    }

	//    --------------------------------
    //    verr = ~v_CS_A*n
	//    --------------------------------
    void realizePositionDotErrorsVirtual(const State& s, int mp,  Real* pverr) const {
        assert(mp==2 && pverr);
        //TODO: should be able to get p info from State
        const Transform& X_AB = getBodyTransform(s, lineBody);
        const Vec3       p_AS = findStationLocation(s, followerBody, defaultFollowerPoint);
        const Vec3       p_BC = ~X_AB * p_AS;
        const Vec3       p_PC = p_BC - defaultPointOnLine;

        const Vec3       v_AS = findStationVelocity(s, followerBody, defaultFollowerPoint);
        const Vec3       v_AC = findStationVelocity(s, lineBody, p_BC);

        const Vec3       v_CS_B = ~X_AB.R()*(v_AS-v_AC); // reexpress in B

        // Calculate these scalar using B-frame vectors, but any frame would have done.
        Vec2::updAs(pverr) = Vec2(~v_CS_B*x, ~v_CS_B*y);
    }

	//    -------------------------------------
    //    aerr = ~(a_CS_A - 2 w_AB X v_CS_A) * n
	//    -------------------------------------
    void realizePositionDotDotErrorsVirtual(const State& s, int mp,  Real* paerr) const {
        assert(mp==2 && paerr);
        //TODO: should be able to get p and v info from State
        const Transform& X_AB = getBodyTransform(s, lineBody);
        const Vec3       p_AS = findStationLocation(s, followerBody, defaultFollowerPoint);
        const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin and reexpress in B;
                                              // C is material point of B coincident with S
        const Vec3       p_PC = p_BC - defaultPointOnLine;

        const Vec3&      w_AB = getBodyAngularVelocity(s, lineBody);
        const Vec3       v_AS = findStationVelocity(s, followerBody, defaultFollowerPoint);
        const Vec3       v_AC = findStationVelocity(s, lineBody, p_BC);

        const Vec3       a_AS = findStationAcceleration(s, followerBody, defaultFollowerPoint);
        const Vec3       a_AC = findStationAcceleration(s, lineBody, p_BC);
        const Vec3       a_CS_B = ~X_AB.R()*(a_AS-a_AC - 2 * w_AB % (v_AS-v_AC));

        // Calculate these scalar using B-frame vectors, but any frame would have done.
        Vec2::updAs(paerr) = Vec2( ~a_CS_B * x, ~a_CS_B * y );
    }

	// apply f=lambda0*x + lambda1*y to the follower point S of body F,
	//      -f                       to point C (coincident point) of body B
    void applyPositionConstraintForcesVirtual
       (const State& s, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForcesInA,
        Vector&              mobilityForces) const
    {
        assert(mp==2 && multipliers);
        const Vec2 lambda = Vec2::getAs(multipliers);

        //TODO: should be able to get p info from State
        const Transform& X_AB    = getBodyTransform(s, lineBody);
        const Vec3&      p_FS    = defaultFollowerPoint; // measured & expressed in F
        const Vec3       p_AS    = findStationLocation(s, followerBody, defaultFollowerPoint);
        const Vec3       p_BC    = ~X_AB * p_AS;         // measured & expressed in B

		const Vec3       force_B = lambda[0] * x + lambda[1] * y;
        const Vec3       force_A = X_AB.R() * force_B;

        addInStationForce(s, followerBody, p_FS,  force_A, bodyForcesInA);
        addInStationForce(s, lineBody,     p_BC, -force_A, bodyForcesInA);
    }

    SimTK_DOWNCAST(PointOnLineImpl, ConstraintImpl);
private:
    friend class Constraint::PointOnLine;

    ConstrainedBodyIndex lineBody;     // B
    ConstrainedBodyIndex followerBody; // F

    UnitVec3          defaultLineDirection;   // z on B, exp. in B frame
    Vec3              defaultPointOnLine;     // P on B, meas&exp in B frame
    Vec3              defaultFollowerPoint;   // S on F, meas&exp in F frame

    // These are just for visualization
    Real lineHalfLength;
    Real pointRadius;

    // TOPOLOGY CACHE (that is, calculated from construction data)
    mutable UnitVec3 x, y;
};

    // CONSTANT ANGLE

class Constraint::ConstantAngleImpl : public ConstraintImpl {
public:
    ConstantAngleImpl()
      : ConstraintImpl(1,0,0), defaultAxisB(), defaultAxisF(), defaultAngle(Pi/2),
        axisLength(1), axisThickness(1), cosineOfDefaultAngle(NaN)
    { }
    ConstantAngleImpl* clone() const { return new ConstantAngleImpl(*this); }

    void calcDecorativeGeometryAndAppendVirtual
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const;

    void setAxisLength(Real length) {
        // length <= 0 means don't display axis
        invalidateTopologyCache();
        axisLength = length > 0 ? length : 0;
    }
    Real getAxisLength() const {return axisLength;}

    void setAxisThickness(Real t) {
        // t <= 0 means don't display axis
        invalidateTopologyCache();
        axisThickness = t > 0 ? t : 0;
    }
    Real getAxisThickness() const {return axisThickness;}

    // Implementation of ContraintRep virtuals
    void realizeTopologyVirtual(State& s) const {
        cosineOfDefaultAngle = std::cos(defaultAngle);
    }


    // Implementation of virtuals required for holonomic constraints.

    // Let B=B1 be the "base" body onto which unit vector b is fixed, and F=B2 the "follower" 
    // body onto which unit vector f is fixed. The angle theta between these vectors is
    // given by cos(theta) = dot(b, f) with the axes expressed in a common basis.
    // This can range from 1 to -1, corresponding to angles 0 to 180 respectively.
    // We would like to enforce the constraint that cos(theta) is a constant. This can be done
    // with a single constraint equation as long as theta is sufficiently far away from 0 and
    // 180, with the numerically best performance at theta=90 degrees where cos(theta)==0.
    //
    // If you want to enforce that two axes are aligned with one another (that is, the angle
    // between them is 0 or 180), that takes *two* constraint equations since the only remaining
    // rotation is about the common axis.
    //
    // We will work in the A frame.
    //
    // ------------------------------
    // perr = ~b_A * f_A - cos(theta)
    // ------------------------------
    //
    // verr = d/dt perr (derivative taken in A)
    //      = ~b_A * (w_AF % f_A) + ~f_A * (w_AB % b_A)
    //      = ~w_AF * (f_A % b_A) - ~w_AB * (f_A % b_A)     (scalar triple product identity)
    // => ------------------------------
    // verr = ~(w_AF-w_AB) * (f_A % b_A)
    // ---------------------------------
    //
    // aerr = d/dt verr (derivative taken in A)
    //      = ~(b_AF-b_AB) * (f_A % b_A)
    //        + (w_AF-w_AB) * ((w_AF%f_A) % b_A)
    //        + (w_AF-w_AB) * (f_A % (w_AB%b_A))
    //      =   ~(b_AF-b_AB) * (f_A % b_A)
    //        + ~(w_AF-w_AB) * ((w_AF%f_A) % b_A) - (w_AB%b_A) % f_A)
    // => -----------------------------------------------------------
    // aerr =   ~(b_AF-b_AB) * (f_A % b_A)
    //        + ~(w_AF-w_AB) * ((w_AF%f_A) % b_A) - (w_AB%b_A) % f_A)
    // --------------------------------------------------------------
    //
    // Constraint torque can be determined by inspection of verr:
    //    lambda * (f_A % b_A) applied to body F
    //   -lambda * (f_A % b_A) applied to body B
    //

    // ------------------------------
    // perr = ~b_A * f_A - cos(theta)
    // ------------------------------
    void realizePositionErrorsVirtual(const State& s, int mp,  Real* perr) const {
        assert(mp==1 && perr);

        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const UnitVec3  b_A  = R_AB * defaultAxisB;
        const UnitVec3  f_A  = R_AF * defaultAxisF;

        *perr = dot(b_A, f_A) - cosineOfDefaultAngle;
    }

    // ----------------------------------
    // pverr = ~(w_AF-w_AB) * (f_A % b_A)
    // ----------------------------------
    void realizePositionDotErrorsVirtual(const State& s, int mp,  Real* pverr) const {
        assert(mp==1 && pverr);
        //TODO: should be able to get p info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const UnitVec3  b_A  = R_AB * defaultAxisB;
        const UnitVec3  f_A  = R_AF * defaultAxisF;
        const Vec3&     w_AB = getBodyAngularVelocity(s, B);
        const Vec3&     w_AF = getBodyAngularVelocity(s, F);

        *pverr = dot( w_AF-w_AB,  f_A % b_A );
    }

    // --------------------------------------------------------------
    // paerr =  ~(b_AF-b_AB) * (f_A % b_A)
    //        + ~(w_AF-w_AB) * ((w_AF%f_A) % b_A) - (w_AB%b_A) % f_A)
    // --------------------------------------------------------------
    void realizePositionDotDotErrorsVirtual(const State& s, int mp,  Real* paerr) const {
        assert(mp==1 && paerr);
        //TODO: should be able to get p and v info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const UnitVec3  b_A  = R_AB * defaultAxisB;
        const UnitVec3  f_A  = R_AF * defaultAxisF;
        const Vec3&     w_AB = getBodyAngularVelocity(s, B);
        const Vec3&     w_AF = getBodyAngularVelocity(s, F);
        const Vec3&     b_AB = getBodyAngularAcceleration(s, B);
        const Vec3&     b_AF = getBodyAngularAcceleration(s, F);

        *paerr =   dot( b_AF-b_AB, f_A % b_A )
                 + dot( w_AF-w_AB, (w_AF%f_A) % b_A - (w_AB%b_A) % f_A);
    }

    //    lambda * (f_A % b_A) applied to body F
    //   -lambda * (f_A % b_A) applied to body B
    void applyPositionConstraintForcesVirtual
       (const State& s, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForcesInA,
        Vector&              mobilityForces) const
    {
        assert(mp==1 && multipliers);
        const Real lambda = *multipliers;
        //TODO: should be able to get p info from State
        const Rotation&  R_AB = getBodyRotation(s, B);
        const Rotation&  R_AF = getBodyRotation(s, F);
        const UnitVec3   b_A = R_AB*defaultAxisB;
        const UnitVec3   f_A = R_AF*defaultAxisF;
        const Vec3       torque_F_A = lambda * (f_A % b_A); // on F, in A frame

        addInBodyTorque(s, F,  torque_F_A, bodyForcesInA);
        addInBodyTorque(s, B, -torque_F_A, bodyForcesInA);
    }

    SimTK_DOWNCAST(ConstantAngleImpl, ConstraintImpl);
private:
    friend class Constraint::ConstantAngle;

    ConstrainedBodyIndex B; // B1 is "base" body
    ConstrainedBodyIndex F; // B2 is "follower" body

    UnitVec3          defaultAxisB; // fixed to B, expressed in B frame
    UnitVec3          defaultAxisF; // fixed to F, expressed in F frame
    Real              defaultAngle; // required angle between axisB and axisF

    // These are just for visualization
    Real axisLength;
    Real axisThickness;

    // TOPOLOGY CACHE (that is, calculated from construction data)
    mutable Real cosineOfDefaultAngle;
};

    // BALL

class Constraint::BallImpl : public ConstraintImpl {
public:
    BallImpl() : ConstraintImpl(3,0,0), defaultPoint1(0), defaultPoint2(0), defaultRadius(0.1) { }
    BallImpl* clone() const { return new BallImpl(*this); }

    void calcDecorativeGeometryAndAppendVirtual
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const;

    void setDefaultRadius(Real r) {
        // r <= 0 means don't display
        invalidateTopologyCache();
        defaultRadius = r > 0 ? r : 0;
    }
    Real getDefaultRadius() const {return defaultRadius;}

    // Implementation of virtuals required for holonomic constraints.

    // We have a ball joint between base body B and follower body F, located at a point P fixed to B
    // and point S fixed on F. All forces will be applied at point S and the coincident material point 
    // C on B which is instantaneously at the same spatial location as S. We will work in the A frame.
    //
    //  First, find the material point C of B that is coincident
    //  in space with point S of F: p_BC = p_AS-p_AB. This vector
    //  is *constant* in the B frame because it is a material point,
    //  despite the fact that its definition involves a point which
    //  moves with respect to B. The velocity constraint is then
    //  very simple: the spatial velocity of point C of B should be
    //  the same as the spatial velocity of point S of F:
    //      verr = v_AS - v_AC = v_AS - (v_AB + w_AB X p_BC) = 0
    //  Integrating to get perr, we get
    //      perr = p_AS - p_AC + constant = 0
    //  But p_AC=p_AS by construction, so perr=constant=0.
    //  The constant is defined by the fact that we want material point
    //  C of B to be in the same spatial location as point P of B,
    //  so constant=p_BC-p_BP=0. Writing in the A frame we have:
    //      perr = p_AS-(p_AB+R_AB*p_BP) = 0 (a constant)
    //      verr = v_AS - (v_AB + w_AB X R_AB*p_BC)
    //      aerr = a_AS - (a_AB + b_AB X R_AB*p_BC + w_AB X (w_AB X R_AB*p_BC))
    //  apply +lambda to S of F, -lambda to C of B.
    //      
    //

    void realizePositionErrorsVirtual(const State& s, int mp,  Real* perr) const {
        assert(mp==3 && perr);

        const Vec3 p_AP = findStationLocation(s, B1, defaultPoint1);
        const Vec3 p_AS = findStationLocation(s, B2, defaultPoint2);

        // See above comments -- this is just the constant of integration; there is a missing (p_AS-p_AC)
        // term (always 0) here which is what we differentiate to get the verr equation.
        Vec3::updAs(perr) = p_AS - p_AP;
    }
 
    void realizePositionDotErrorsVirtual(const State& s, int mp,  Real* pverr) const {
        assert(mp==3 && pverr);
        //TODO: should be able to get p info from State
        const Transform&  X_AB   = getBodyTransform(s, B1);
        const Vec3        p_AS   = findStationLocation(s, B2, defaultPoint2);
        const Vec3        p_BC   = ~X_AB*p_AS; // C is a material point of body B

        const Vec3        v_AS    = findStationVelocity(s, B2, defaultPoint2);
        const Vec3        v_AC    = findStationVelocity(s, B1, p_BC);
        Vec3::updAs(pverr) = v_AS - v_AC;
    }

    void realizePositionDotDotErrorsVirtual(const State& s, int mp,  Real* paerr) const {
        assert(mp==3 && paerr);
        //TODO: should be able to get p and v info from State

        const Transform&  X_AB   = getBodyTransform(s, B1);
        const Vec3        p_AS   = findStationLocation(s, B2, defaultPoint2);
        const Vec3        p_BC   = ~X_AB*p_AS; // C is a material point of body B

        const Vec3        a_AS    = findStationAcceleration(s, B2, defaultPoint2);
        const Vec3        a_AC    = findStationAcceleration(s, B1, p_BC);
        Vec3::updAs(paerr) = a_AS - a_AC;
    }

    void applyPositionConstraintForcesVirtual
       (const State& s, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForcesInA,
        Vector&              mobilityForces) const
    {
        assert(mp==3 && multipliers);

        //TODO: should be able to get p info from State
        const Transform& X_AB  = getBodyTransform(s,B1);
        const Vec3&      p_FS  = defaultPoint2;
        const Vec3       p_AS  = findStationLocation(s, B2, p_FS);
        const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin and reexpress in B;
                                              // C is material point of B coincident with S

        const Vec3 force_A = Vec3::getAs(multipliers);

        // Multipliers are force to be applied to S on F, but
        // apply the -force not to point P of B, but to the point "C" of B
        // coincident with S, which won't be exactly the same place
        // as P if the position-level constraint isn't met exactly.

        addInStationForce(s, B2, p_FS,  force_A, bodyForcesInA);
        addInStationForce(s, B1, p_BC, -force_A, bodyForcesInA);
    }

    SimTK_DOWNCAST(BallImpl, ConstraintImpl);
private:
    friend class Constraint::Ball;

    ConstrainedBodyIndex B1;
    ConstrainedBodyIndex B2;

    Vec3            defaultPoint1; // on body 1, exp. in B1 frame
    Vec3            defaultPoint2; // on body 2, exp. in B2 frame
    Real            defaultRadius; // used for visualization only
};

    // CONSTANT ORIENTATION

class Constraint::ConstantOrientationImpl : public ConstraintImpl {
public:
    ConstantOrientationImpl()
      : ConstraintImpl(3,0,0), defaultRB(), defaultRF()
    { }
    ConstantOrientationImpl* clone() const { return new ConstantOrientationImpl(*this); }

    //TODO: visualization?


    // Implementation of virtuals required for holonomic constraints.

    // Let B=B1 be the "base" body onto which rotation matrix RB is fixed, and F=B2 the "follower" 
    // body onto which rotation matrix RF is fixed. We would like to enforce the constraint
    // that RB==RF when both are expressed in a common basis. (Remember that a rotation matrix
    // is just three axis vectors.)
    // 
    // Here the (redundant) assembly constraint is that all the axes are parallel, that is
    // RBx==RFx, RBy==RFy, and RBz==RFz. However, aligning two vectors takes *two* constraints
    // so that would be a total of 6 constraints, with only 3 independent.
    // The independent runtime constraints just enforce perpendicularity, but can be satisfied
    // in cases where some of the axes are antiparallel so are not suited for the initial assembly.
    // The runtime constraints are thus three "constant angle" constraints, where the angle
    // is always 90 degrees:
    //
    //    ~RFx * RBy = 0
    //    ~RFy * RBz = 0
    //    ~RFz * RBx = 0
    //
    // We'll work in A. See the "constant angle" constraint for details.
    //
    // -----------------
    // perr = ~RFx * RBy  (with all axes expressed in A)
    //        ~RFy * RBz
    //        ~RFz * RBx
    // -----------------
    //
    // ---------------------------------
    // verr = ~(w_AF-w_AB) * (RFx % RBy)
    //      = ~(w_AF-w_AB) * (RFy % RBz)
    //      = ~(w_AF-w_AB) * (RFz % RBx)
    // ---------------------------------
    //
    // -----------------------------------------------------------------------
    // aerr = ~(b_AF-b_AB) * (RFx % RBy)
    //                 + ~(w_AF-w_AB) * ((w_AF%RFx) % RBy) - (w_AB%RBy) % RFx)
    //        ~(b_AF-b_AB) * (RFy % RBz)
    //                 + ~(w_AF-w_AB) * ((w_AF%RFy) % RBz) - (w_AB%RBz) % RFy)
    //        ~(b_AF-b_AB) * (RFz % RBx)
    //                 + ~(w_AF-w_AB) * ((w_AF%RFz) % RBx) - (w_AB%RBx) % RFz)
    // -----------------------------------------------------------------------
    //
    // Constraint torque can be determined by inspection of verr:
    //    t_F =   lambda_x * (RFx % RBy)   (applied to body F)
    //          + lambda_y * (RFy % RBz)
    //          + lambda_z * (RFz % RBx)
    //    t_B = -t_F                       (applied to body B)
    //

    // -----------------
    // perr = ~RFx * RBy  (with all axes expressed in A)
    //        ~RFy * RBz
    //        ~RFz * RBx
    // -----------------
    void realizePositionErrorsVirtual(const State& s, int mp,  Real* perr) const {
        assert(mp==3 && perr);

        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const Rotation  RB = R_AB * defaultRB; // now expressed in A
        const Rotation  RF = R_AF * defaultRF;

        Vec3::updAs(perr) = Vec3(~RF.x()*RB.y(),
                                 ~RF.y()*RB.z(),
                                 ~RF.z()*RB.x());
    }

    // ----------------------------------
    // verr = ~(w_AF-w_AB) * (RFx % RBy)
    //      = ~(w_AF-w_AB) * (RFy % RBz)
    //      = ~(w_AF-w_AB) * (RFz % RBx)
    // ----------------------------------
    void realizePositionDotErrorsVirtual(const State& s, int mp,  Real* pverr) const {
        assert(mp==3 && pverr);
        //TODO: should be able to get p info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const Rotation  RB = R_AB * defaultRB; // now expressed in A
        const Rotation  RF = R_AF * defaultRF;

        const Vec3&     w_AB = getBodyAngularVelocity(s, B);
        const Vec3&     w_AF = getBodyAngularVelocity(s, F);
        const Vec3      w_BF = w_AF-w_AB; // in A

        Vec3::updAs(pverr) = Vec3( ~w_BF * (RF.x() % RB.y()),
                                   ~w_BF * (RF.y() % RB.z()),
                                   ~w_BF * (RF.z() % RB.x()) );
    }

    //------------------------------------------------------------------------
    // aerr = ~(b_AF-b_AB) * (RFx % RBy)
    //                 + ~(w_AF-w_AB) * ((w_AF%RFx) % RBy) - (w_AB%RBy) % RFx)
    //        ~(b_AF-b_AB) * (RFy % RBz)
    //                 + ~(w_AF-w_AB) * ((w_AF%RFy) % RBz) - (w_AB%RBz) % RFy)
    //        ~(b_AF-b_AB) * (RFz % RBx)
    //                 + ~(w_AF-w_AB) * ((w_AF%RFz) % RBx) - (w_AB%RBx) % RFz)
    //------------------------------------------------------------------------
    void realizePositionDotDotErrorsVirtual(const State& s, int mp,  Real* paerr) const {
        assert(mp==3 && paerr);
        //TODO: should be able to get p and v info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const Rotation  RB = R_AB * defaultRB; // now expressed in A
        const Rotation  RF = R_AF * defaultRF;

        const Vec3&     w_AB = getBodyAngularVelocity(s, B);
        const Vec3&     w_AF = getBodyAngularVelocity(s, F);
        const Vec3      w_BF = w_AF-w_AB; // in A

        const Vec3&     b_AB = getBodyAngularAcceleration(s, B);
        const Vec3&     b_AF = getBodyAngularAcceleration(s, F);
        const Vec3      b_BF = b_AF-b_AB; // in A

        Vec3::updAs(paerr) = 
             Vec3( dot( b_BF, RF.x() % RB.y() )
                        + dot( w_BF, (w_AF%RF.x()) % RB.y() - (w_AB%RB.y()) % RF.x()),
                   dot( b_BF, RF.y() % RB.z() )
                        + dot( w_BF, (w_AF%RF.y()) % RB.z() - (w_AB%RB.z()) % RF.y()),
                   dot( b_BF, RF.z() % RB.x() )
                        + dot( w_BF, (w_AF%RF.z()) % RB.x() - (w_AB%RB.x()) % RF.z()));
    }

    //    t_F =   lambda_x * (RFx % RBy)   (applied to body F)
    //          + lambda_y * (RFy % RBz)
    //          + lambda_z * (RFz % RBx)
    //    t_B = -t_F                       (applied to body B)
    void applyPositionConstraintForcesVirtual
       (const State& s, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForcesInA,
        Vector&              mobilityForces) const
    {
        assert(mp==3 && multipliers);
        const Vec3& lambda = Vec3::getAs(multipliers);

        //TODO: should be able to get p info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const Rotation  RB = R_AB * defaultRB; // now expressed in A
        const Rotation  RF = R_AF * defaultRF;

        const Vec3 torque_F_A =   lambda[0] * (RF.x() % RB.y())
                                + lambda[1] * (RF.y() % RB.z())
                                + lambda[2] * (RF.z() % RB.x());

        addInBodyTorque(s, F,  torque_F_A, bodyForcesInA);
        addInBodyTorque(s, B, -torque_F_A, bodyForcesInA);
    }

    SimTK_DOWNCAST(ConstantOrientationImpl, ConstraintImpl);
private:
    friend class Constraint::ConstantOrientation;

    ConstrainedBodyIndex B; // B1 is "base" body
    ConstrainedBodyIndex F; // B2 is "follower" body

    Rotation          defaultRB; // fixed to B, expressed in B frame; RB = R_B_RB
    Rotation          defaultRF; // fixed to F, expressed in F frame; RF = R_F_RF
};

    // WELD

class Constraint::WeldImpl : public ConstraintImpl {
    static Real getDefaultAxisDisplayLength() {return 1;}
    static Vec3 getDefaultFrameColor(int which) {
        return which==0 ? Blue : Purple;
    }
public:
    WeldImpl() 
      : ConstraintImpl(6,0,0), axisDisplayLength(-1), // means "use default axis length"
        frameBColor(-1), frameFColor(-1) // means "use default colors"
    {   // default Transforms are identity, i.e. body frames
    }
    WeldImpl* clone() const { return new WeldImpl(*this); }

    // Draw the two frames.
    void calcDecorativeGeometryAndAppendVirtual
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const;

    void setAxisDisplayLength(Real len) {
        // len == 0 means "don't display"
        // len < 0 means "use default"
        invalidateTopologyCache();
        axisDisplayLength = len >= 0 ? len : -1;
    }
    Real getAxisDisplayLength() const {
        return axisDisplayLength < 0 ? getDefaultAxisDisplayLength() : axisDisplayLength;
    }

    void setFrameColor(int which, const Vec3& color) {
        assert(which==0 || which==1);
        // color[0] < 0 means "use default color for this frame"
        invalidateTopologyCache();
        if (which==0) frameBColor = color[0] < 0 ? Vec3(-1) : color;
        else          frameFColor = color[0] < 0 ? Vec3(-1) : color;
    }
    Vec3 getFrameColor(int which) const {
        assert(which==0 || which==1);
        if (which==0) return frameBColor[0] < 0 ? getDefaultFrameColor(0) : frameBColor;
        else          return frameFColor[0] < 0 ? getDefaultFrameColor(1) : frameFColor;
    }

    // Implementation of virtuals required for holonomic constraints.

    // For theory, look at the ConstantOrientation (1st 3 equations) and 
    // Ball (last 3 equations) theory above. Otherwise just lay back and 
    // enjoy the ride.

    void realizePositionErrorsVirtual(const State& s, int mp,  Real* perr) const {
        assert(mp==6 && perr);

        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const Rotation  RB = R_AB * defaultFrameB.R(); // now expressed in A
        const Rotation  RF = R_AF * defaultFrameF.R();

        // Orientation error
        Vec3::updAs(perr) = Vec3(~RF.x()*RB.y(),
                                 ~RF.y()*RB.z(),
                                 ~RF.z()*RB.x());

        const Vec3 p_AF1 = findStationLocation(s, B, defaultFrameB.p());
        const Vec3 p_AF2 = findStationLocation(s, F, defaultFrameF.p());

        // position error
        Vec3::updAs(perr+3) = p_AF2 - p_AF1;
    }

    void realizePositionDotErrorsVirtual(const State& s, int mp,  Real* pverr) const {
        assert(mp==6 && pverr);
        //TODO: should be able to get p info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const Rotation  RB = R_AB * defaultFrameB.R(); // now expressed in A
        const Rotation  RF = R_AF * defaultFrameF.R();

        const Vec3&     w_AB = getBodyAngularVelocity(s, B);
        const Vec3&     w_AF = getBodyAngularVelocity(s, F);
        const Vec3      w_BF = w_AF-w_AB; // in A

        // orientation error
        Vec3::updAs(pverr) = Vec3( ~w_BF * (RF.x() % RB.y()),
                                   ~w_BF * (RF.y() % RB.z()),
                                   ~w_BF * (RF.z() % RB.x()) );

        //TODO: should be able to get p info from State
        const Transform&  X_AB   = getBodyTransform(s, B);
        const Vec3        p_AF2  = findStationLocation(s, F, defaultFrameF.p());
        const Vec3        p_BC   = ~X_AB*p_AF2; // C is a material point of body B

        const Vec3        v_AF2   = findStationVelocity(s, F, defaultFrameF.p());
        const Vec3        v_AC    = findStationVelocity(s, B, p_BC);
 
        // position error
        Vec3::updAs(pverr+3) = v_AF2 - v_AC;
    }

    void realizePositionDotDotErrorsVirtual(const State& s, int mp,  Real* paerr) const {
        assert(mp==6 && paerr);
        //TODO: should be able to get p and v info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const Rotation  RB = R_AB * defaultFrameB.R(); // now expressed in A
        const Rotation  RF = R_AF * defaultFrameF.R();

        const Vec3&     w_AB = getBodyAngularVelocity(s, B);
        const Vec3&     w_AF = getBodyAngularVelocity(s, F);
        const Vec3      w_BF = w_AF-w_AB; // in A

        const Vec3&     b_AB = getBodyAngularAcceleration(s, B);
        const Vec3&     b_AF = getBodyAngularAcceleration(s, F);
        const Vec3      b_BF = b_AF-b_AB; // in A

        // orientation error
        Vec3::updAs(paerr) = 
             Vec3( dot( b_BF, RF.x() % RB.y() )
                        + dot( w_BF, (w_AF%RF.x()) % RB.y() - (w_AB%RB.y()) % RF.x()),
                   dot( b_BF, RF.y() % RB.z() )
                        + dot( w_BF, (w_AF%RF.y()) % RB.z() - (w_AB%RB.z()) % RF.y()),
                   dot( b_BF, RF.z() % RB.x() )
                        + dot( w_BF, (w_AF%RF.z()) % RB.x() - (w_AB%RB.x()) % RF.z()));

        const Transform&  X_AB   = getBodyTransform(s, B);
        const Vec3        p_AF2  = findStationLocation(s, F, defaultFrameF.p());
        const Vec3        p_BC   = ~X_AB*p_AF2; // C is a material point of body B

        const Vec3        a_AF2  = findStationAcceleration(s, F, defaultFrameF.p());
        const Vec3        a_AC   = findStationAcceleration(s, B, p_BC);

        // position error
        Vec3::updAs(paerr+3) = a_AF2 - a_AC;
    }

    void applyPositionConstraintForcesVirtual
       (const State& s, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForcesInA,
        Vector&              mobilityForces) const
    {
        assert(mp==6 && multipliers);

        const Vec3& torques = Vec3::getAs(multipliers);
        const Vec3& force_A = Vec3::getAs(multipliers+3);

        //TODO: should be able to get p info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const Rotation  RB = R_AB * defaultFrameB.R(); // now expressed in A
        const Rotation  RF = R_AF * defaultFrameF.R();

        const Vec3 torque_F_A =   torques[0] * (RF.x() % RB.y())
                                + torques[1] * (RF.y() % RB.z())
                                + torques[2] * (RF.z() % RB.x());

        addInBodyTorque(s, F,  torque_F_A, bodyForcesInA);
        addInBodyTorque(s, B, -torque_F_A, bodyForcesInA);

        const Transform& X_AB  = getBodyTransform(s,B);
        const Vec3&      p_FF2 = defaultFrameF.p();
        const Vec3       p_AF2 = findStationLocation(s, F, p_FF2);
        const Vec3       p_BC = ~X_AB * p_AF2;

        addInStationForce(s, F, p_FF2, force_A, bodyForcesInA);
        addInStationForce(s, B, p_BC, -force_A, bodyForcesInA);
    }

    SimTK_DOWNCAST(WeldImpl, ConstraintImpl);
private:
    friend class Constraint::Weld;

    ConstrainedBodyIndex B; // aka "body 1"
    ConstrainedBodyIndex F; // aka "body 2"

    Transform       defaultFrameB; // on body 1, relative to B frame
    Transform       defaultFrameF; // on body 2, relative to F frame};

    // These are for visualization control only.
    Real axisDisplayLength; // for all 6 axes; <= 0 means "don't display"
    Vec3 frameBColor, frameFColor;
};

    // NO SLIP 1D

class Constraint::NoSlip1DImpl : public ConstraintImpl {
public:
    NoSlip1DImpl()
      : ConstraintImpl(0,1,0), defaultNoSlipDirection(), defaultContactPoint(0),
        directionLength(1), pointRadius(0.05) 
    { }
    NoSlip1DImpl* clone() const { return new NoSlip1DImpl(*this); }

    void calcDecorativeGeometryAndAppendVirtual
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const;

    void setDirectionDisplayLength(Real l) {
        // l <= 0 means don't display direction line
        invalidateTopologyCache();
        directionLength = l > 0 ? l : 0;
    }
    Real getDirectionDisplayLength() const {return directionLength;}

    void setPointDisplayRadius(Real r) {
        // r <= 0 means don't display point
        invalidateTopologyCache();
        pointRadius= r > 0 ? r : 0;
    }
    Real getPointDisplayRadius() const {return pointRadius;}

    // Implementation of virtuals required for nonholonomic constraints.

    // One non-holonomic constraint equation. There is a contact point P and a no-slip 
    // direction n fixed in a case body C. There are two moving bodies B0 and B1. The 
    // material point P0 of B0 and the material point P1 of B1 which are each coincident 
    // with the contact point P must have identical velocities in C, along the direction n.
    // This can be used to implement simple rolling contact between disks, such as occurs
    // in gear trains.
    //
    // There is no perr equation here since this is a non-holonomic (velocity) constraint.
    // In the C frame, the constraint we want is
    //    verr = ~(v_CP1 - v_CP0) * n_C
    // that is, the two contact points have no relative velocity in C along the normal.
    // We can calculate this in A instead since the velocities in C of each point will 
    // differ from their velocities in A by a constant (because they are both in the same
    // place in space). So:
    //    verr = ~(v_AP1 - v_AP0) * n_A
    // Differentiating material point velocities in A, we get the acceleration error
    //    aerr = ~(a_AP1 - a_AP0) * n_A + ~(v_AP1 - v_AP0) * (w_AC X n_A)
    //         = ~(a_AP1 - a_AP0 - w_AC X (v_AP1 - v_AP0)) * n_A
    // 
    void realizeVelocityErrorsVirtual(const State& s, int mv,  Real* verr) const {
        assert(mv==1 && verr);
        //TODO: should be able to get p info from State
        const Transform& X_AC  = getBodyTransform(s, caseBody);
        const Transform& X_AB0 = getBodyTransform(s, movingBody0);
        const Transform& X_AB1 = getBodyTransform(s, movingBody1);
        const Vec3       p_AP  = X_AC * defaultContactPoint; // P's location in A
        const Vec3       p_P0  = ~X_AB0 * p_AP;              // P0's station in B0
        const Vec3       p_P1  = ~X_AB1 * p_AP;              // P1's station in B1
        const UnitVec3   n_A   = X_AC.R() * defaultNoSlipDirection;

        const Vec3       v_AP0 = findStationVelocity(s, movingBody0, p_P0);
        const Vec3       v_AP1 = findStationVelocity(s, movingBody1, p_P1);

        // Calculate this scalar using A-frame vectors.
        *verr = ~(v_AP1-v_AP0) * n_A;
    }

    void realizeVelocityDotErrorsVirtual(const State& s, int mv,  Real* vaerr) const {
        assert(mv==1 && vaerr);
        //TODO: should be able to get p and v info from State
        const Transform& X_AC  = getBodyTransform(s, caseBody);
        const Transform& X_AB0 = getBodyTransform(s, movingBody0);
        const Transform& X_AB1 = getBodyTransform(s, movingBody1);
        const Vec3       p_AP  = X_AC * defaultContactPoint; // P's location in A
        const Vec3       p_P0  = ~X_AB0 * p_AP;              // P0's station in B0
        const Vec3       p_P1  = ~X_AB1 * p_AP;              // P1's station in B1
        const UnitVec3   n_A   = X_AC.R() * defaultNoSlipDirection;

        const Vec3       v_AP0 = findStationVelocity(s, movingBody0, p_P0);
        const Vec3       v_AP1 = findStationVelocity(s, movingBody1, p_P1);
        const Vec3&      w_AC  = getBodyAngularVelocity(s, caseBody);

        const Vec3       a_AP0 = findStationAcceleration(s, movingBody0, p_P0);
        const Vec3       a_AP1 = findStationAcceleration(s, movingBody1, p_P1);

        // Calculate this scalar using A-frame vectors.
        *vaerr = ~(a_AP1-a_AP0 - w_AC % (v_AP1-v_AP0)) * n_A;
    }

	// apply f=lambda*n to contact point P1 of B1,
	//      -f          to contact point P2 of B2
    void applyVelocityConstraintForcesVirtual
       (const State& s, int mv, const Real* multipliers,
        Vector_<SpatialVec>& bodyForcesInA,
        Vector&              mobilityForces) const
    {
        assert(mv==1 && multipliers);
        const Real lambda = *multipliers;

        //TODO: should be able to get p info from State
        const Transform& X_AC  = getBodyTransform(s, caseBody);
        const Transform& X_AB0 = getBodyTransform(s, movingBody0);
        const Transform& X_AB1 = getBodyTransform(s, movingBody1);
        const Vec3       p_AP  = X_AC * defaultContactPoint; // P's location in A
        const Vec3       p_P0  = ~X_AB0 * p_AP;              // P0's station in B0
        const Vec3       p_P1  = ~X_AB1 * p_AP;              // P1's station in B1

		const Vec3       force_A = X_AC.R()*(lambda*defaultNoSlipDirection);

        addInStationForce(s, movingBody1, p_P1,  force_A, bodyForcesInA);
        addInStationForce(s, movingBody0, p_P0, -force_A, bodyForcesInA);
    }

    SimTK_DOWNCAST(NoSlip1DImpl, ConstraintImpl);
private:
    friend class Constraint::NoSlip1D;

    ConstrainedBodyIndex caseBody;     // C
    ConstrainedBodyIndex movingBody0;  // B0
    ConstrainedBodyIndex movingBody1;  // B1

    UnitVec3          defaultNoSlipDirection;   // on body C, exp. in C frame
    Vec3              defaultContactPoint;      // on body C, exp. in C frame

    // These are just for visualization
    Real directionLength;
    Real pointRadius;
};

    // CONSTANT SPEED

class Constraint::ConstantSpeedImpl : public ConstraintImpl {
public:
    ConstantSpeedImpl()
      : ConstraintImpl(0,1,0), theMobilizer(), whichMobility(), prescribedSpeed(NaN)
    { }
    ConstantSpeedImpl* clone() const { return new ConstantSpeedImpl(*this); }

    // Implementation of virtuals required for nonholonomic constraints.

    // One non-holonomic (well, velocity-level) constraint equation.
    //    verr = u - s
    //    aerr = udot
    // 
    void realizeVelocityErrorsVirtual(const State& s, int mv,  Real* verr) const {
        assert(mv==1 && verr);
        *verr = getOneU(s, theMobilizer, whichMobility) - prescribedSpeed;
    }

    void realizeVelocityDotErrorsVirtual(const State& s, int mv,  Real* vaerr) const {
        assert(mv==1 && vaerr);
        *vaerr = getOneUDot(s, theMobilizer, whichMobility, true);
    }

	// apply generalized force lambda to the mobility
    void applyVelocityConstraintForcesVirtual
       (const State& s, int mv, const Real* multipliers,
        Vector_<SpatialVec>& bodyForcesInA,
        Vector&              mobilityForces) const
    {
        assert(mv==1 && multipliers);
        const Real lambda = *multipliers;
        addInOneMobilityForce(s, theMobilizer, whichMobility, lambda, mobilityForces);
    }

    SimTK_DOWNCAST(ConstantSpeedImpl, ConstraintImpl);
private:
    friend class Constraint::ConstantSpeed;

    ConstrainedMobilizerIndex   theMobilizer;
    MobilizerUIndex             whichMobility;
    Real                        prescribedSpeed;
};


    /////////////////////////////////////////////
    // CONSTRAINT::CUSTOM::IMPLEMENTATION IMPL //
    /////////////////////////////////////////////

// This class exists primarily to allow the Custom::Implementation class to keep
// a pointer to its handle class's CustomImpl class which is derived from ConstraintImpl
// which has all the goodies that are needed for defining a Constraint.
//
// At first this class is the owner of the CustomImpl. Then when this is put in a
// Custom handle, that handle takes over ownership of the CustomImpl and the 
// CustomImpl takes over ownership of this ImplementationImpl object.
class Constraint::Custom::ImplementationImpl 
  : public PIMPLImplementation<Implementation, ImplementationImpl>
{
public:
    // no default constructor
    explicit ImplementationImpl(CustomImpl* customImpl) : isOwner(true), builtInImpl(customImpl) { }
    inline ~ImplementationImpl(); // see below -- have to wait for CustomImpl's definition

    // Copying one of these just gives us a new one with a NULL CustomImpl pointer.
    ImplementationImpl(const ImplementationImpl& src) : isOwner(false), builtInImpl(0) { }

    ImplementationImpl* clone() const {return new ImplementationImpl(*this);}

    bool isOwnerOfCustomImpl() const {return builtInImpl && isOwner;}
    CustomImpl* removeOwnershipOfCustomImpl() {
        assert(isOwnerOfCustomImpl()); 
        isOwner=false; 
        return builtInImpl;
    }

    void setReferenceToCustomImpl(CustomImpl* cimpl) {
        assert(!builtInImpl); // you can only do this once
        isOwner=false;
        builtInImpl = cimpl;
    }

    bool hasCustomImpl() const {return builtInImpl != 0;}

    const CustomImpl& getCustomImpl() const {
        assert(builtInImpl);
        return *builtInImpl;
    }
    CustomImpl& updCustomImpl() {
        assert(builtInImpl);
        return *builtInImpl;
    }

private:
    bool isOwner;
    CustomImpl* builtInImpl; // just a reference; not owned

    // suppress assignment
    ImplementationImpl& operator=(const ImplementationImpl&);
};

    /////////////////////////////
    // CONSTRAINT::CUSTOM IMPL //
    /////////////////////////////

class Constraint::CustomImpl : public ConstraintImpl {
public:
    CustomImpl() : implementation(0) { }
    CustomImpl(int mp, int mv, int ma) : ConstraintImpl(mp,mv,ma), implementation(0) { }

    void takeOwnershipOfImplementation(Custom::Implementation* userImpl);

    explicit CustomImpl(Custom::Implementation* userImpl) : implementation(0) { 
        assert(userImpl);
        implementation = userImpl;
        implementation->updImpl().setReferenceToCustomImpl(this);
    }    

    // Copy constructor
    CustomImpl(const CustomImpl& src) : implementation(0) {
        if (src.implementation) {
            implementation = src.implementation->clone();
            implementation->updImpl().setReferenceToCustomImpl(this);
        }
    }
    
    ~CustomImpl() {
        delete implementation;
    }
    
    CustomImpl* clone() const { return new CustomImpl(*this); }

    const Custom::Implementation& getImplementation() const {
        assert(implementation);
        return *implementation;
    }

    Custom::Implementation& updImplementation() {
        assert(implementation);
        return *implementation;
    }

    // Forward all the virtuals to the Custom::Implementation virtuals.
    void realizeTopologyVirtual(State& s) const {getImplementation().realizeTopology(s);}
    void realizeModelVirtual   (State& s) const {getImplementation().realizeModel(s);}
    void realizeInstanceVirtual(const State& s) const {getImplementation().realizeInstance(s);}
    void realizeTimeVirtual    (const State& s) const {getImplementation().realizeTime(s);}
    void realizePositionVirtual(const State& s) const {getImplementation().realizePosition(s);}
    void realizeVelocityVirtual(const State& s) const {getImplementation().realizeVelocity(s);}
    void realizeDynamicsVirtual(const State& s) const {getImplementation().realizeDynamics(s);}
    void realizeAccelerationVirtual(const State& s) const {getImplementation().realizeAcceleration(s);}
    void realizeReportVirtual  (const State& s) const {getImplementation().realizeReport(s);}

    void realizePositionErrorsVirtual      (const State& s, int mp,  Real* perr) const
       {getImplementation().realizePositionErrors(s,mp,perr);}
    void realizePositionDotErrorsVirtual   (const State& s, int mp,  Real* pverr) const
       {getImplementation().realizePositionDotErrors(s,mp,pverr);}
    void realizePositionDotDotErrorsVirtual(const State& s, int mp,  Real* paerr) const
       {getImplementation().realizePositionDotDotErrors(s,mp,paerr);}
    void applyPositionConstraintForcesVirtual
           (const State& s, int mp, const Real* multipliers,
            Vector_<SpatialVec>& bodyForces,
            Vector&              mobilityForces) const
       {getImplementation().applyPositionConstraintForces(s,mp,multipliers,bodyForces,mobilityForces);}

    void realizeVelocityErrorsVirtual(const State& s, int mv,  Real* verr) const
       {getImplementation().realizeVelocityErrors(s,mv,verr);}
    void realizeVelocityDotErrorsVirtual(const State& s, int mv,  Real* vaerr) const
       {getImplementation().realizeVelocityDotErrors(s,mv,vaerr);}
    void applyVelocityConstraintForcesVirtual
           (const State& s, int mv, const Real* multipliers,
            Vector_<SpatialVec>& bodyForces,
            Vector&              mobilityForces) const
       {getImplementation().applyVelocityConstraintForces(s,mv,multipliers,bodyForces,mobilityForces);}

    void realizeAccelerationErrorsVirtual(const State& s, int ma,  Real* aerr) const
       {getImplementation().realizeAccelerationErrors(s,ma,aerr);}
    void applyAccelerationConstraintForcesVirtual
           (const State& s, int ma, const Real* multipliers,
            Vector_<SpatialVec>& bodyForces,
            Vector&              mobilityForces) const
       {getImplementation().applyAccelerationConstraintForces(s,ma,multipliers,bodyForces,mobilityForces);}

    void calcDecorativeGeometryAndAppendVirtual
          (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
       {getImplementation().calcDecorativeGeometryAndAppend(s,stage,geom);}

    SimTK_DOWNCAST(CustomImpl, ConstraintImpl);
private:
    friend class Constraint::Custom;

    Custom::Implementation* implementation;

    CustomImpl& operator=(const CustomImpl&); // suppress assignment
};

// Need definition for CustomImpl here in case we have to delete it.
inline Constraint::Custom::ImplementationImpl::~ImplementationImpl() {
    if (isOwner) 
        delete builtInImpl; 
    builtInImpl=0;
}

/////////////////////////////////////////
// CONSTRAINT::COORDINATE COUPLER IMPL //
/////////////////////////////////////////

class Constraint::CoordinateCouplerImpl : public Constraint::Custom::Implementation {
public:
    CoordinateCouplerImpl(SimbodyMatterSubsystem& matter, const Function* function, const Array_<MobilizedBodyIndex>& coordBody, const Array_<MobilizerQIndex>& coordIndex);
    
    ~CoordinateCouplerImpl() {
        if (--referenceCount[0] == 0) {
            delete function;
            delete[] referenceCount;
        }
    }
    
    Implementation* clone() const {
        referenceCount[0]++;
        return new CoordinateCouplerImpl(*this);
    }

    void realizePositionErrors(const State& s, int mp,  Real* perr) const;

    void realizePositionDotErrors(const State& s, int mp,  Real* pverr) const;

    void realizePositionDotDotErrors(const State& s, int mp,  Real* paerr) const;

    void applyPositionConstraintForces(const State& s, int mp, const Real* multipliers, Vector_<SpatialVec>& bodyForces, Vector& mobilityForces) const;

private:
    const Function* function;
    int* referenceCount;
    Array_<ConstrainedMobilizerIndex> coordBodies;
    Array_<MobilizerQIndex> coordIndices;
    mutable Vector temp;
};

////////////////////////////////////
// CONSTRAINT::SPEED COUPLER IMPL //
////////////////////////////////////

class Constraint::SpeedCouplerImpl : public Constraint::Custom::Implementation {
public:
    SpeedCouplerImpl(SimbodyMatterSubsystem& matter, const Function* function, const Array_<MobilizedBodyIndex>& speedBody, const Array_<MobilizerUIndex>& speedIndex,
            const Array_<MobilizedBodyIndex>& coordBody, const Array_<MobilizerQIndex>& coordIndex);
    
    ~SpeedCouplerImpl() {
        if (--referenceCount[0] == 0) {
            delete function;
            delete[] referenceCount;
        }
    }
    
    Implementation* clone() const {
        referenceCount[0]++;
        return new SpeedCouplerImpl(*this);
    }

    void realizeVelocityErrors(const State& s, int mv,  Real* verr) const;

    void realizeVelocityDotErrors(const State& s, int mv,  Real* vaerr) const;

    void applyVelocityConstraintForces(const State& s, int mv, const Real* multipliers, Vector_<SpatialVec>& bodyForces, Vector& mobilityForces) const;
    
private:
    void findArguments(const State& s) const {
        for (int i = 0; i < (int) speedBodies.size(); ++i)
            temp[i] = getOneU(s, speedBodies[i], speedIndices[i]);
        for (int i = 0; i < (int) coordBodies.size(); ++i)
            temp[i+speedBodies.size()] = getMatterSubsystem().getMobilizedBody(coordBodies[i]).getOneQ(s, coordIndices[i]);
    }

    const Function* function;
    int* referenceCount;
    Array_<ConstrainedMobilizerIndex> speedBodies;
    Array_<MobilizedBodyIndex> coordBodies;
    Array_<MobilizerUIndex> speedIndices;
    Array_<MobilizerQIndex> coordIndices;
    mutable Vector temp;
};

////////////////////////////////////////
// CONSTRAINT::PRESCRIBED MOTION IMPL //
////////////////////////////////////////

class Constraint::PrescribedMotionImpl : public Constraint::Custom::Implementation {
public:
    PrescribedMotionImpl(SimbodyMatterSubsystem& matter, const Function* function, MobilizedBodyIndex coordBody, MobilizerQIndex coordIndex);
    
    ~PrescribedMotionImpl() {
        if (--referenceCount[0] == 0) {
            delete function;
            delete[] referenceCount;
        }
    }
    
    Implementation* clone() const {
        referenceCount[0]++;
        return new PrescribedMotionImpl(*this);
    }

    void realizePositionErrors(const State& s, int mp,  Real* perr) const;

    void realizePositionDotErrors(const State& s, int mp,  Real* pverr) const;

    void realizePositionDotDotErrors(const State& s, int mp,  Real* paerr) const;

    void applyPositionConstraintForces(const State& s, int mp, const Real* multipliers, Vector_<SpatialVec>& bodyForces, Vector& mobilityForces) const;

private:
    const Function* function;
    int* referenceCount;
    ConstrainedMobilizerIndex coordBody;
    MobilizerQIndex coordIndex;
    mutable Vector temp;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_IMPL_H_



