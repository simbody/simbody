#ifndef SimTK_SIMBODY_CONSTRAINT_REP_H_
#define SimTK_SIMBODY_CONSTRAINT_REP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-8 Stanford University and the Authors.         *
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/SimbodyMatterSubsystem_Subtree.h"

#include <map>
#include <vector>
#include <utility>  // std::pair
#include <iostream>
using std::cout; using std::endl;

class ConstraintNode;
class SimbodyMatterSubsystemRep;

class SBModelCache;
class SBPositionCache;
class SBVelocityCache;
class SBAccelerationCache;

namespace SimTK {

    /////////////////////
    // CONSTRAINT REPS //
    /////////////////////

class Constraint::ConstraintRep {
public:
    ConstraintRep() : myHandle(0), myMatterSubsystemRep(0), 
        defaultMp(0), defaultMv(0), defaultMa(0),
        myConstraintNode(0) 
    {
    }

    ConstraintRep(int mp, int mv, int ma) : myHandle(0), myMatterSubsystemRep(0), 
        defaultMp(mp), defaultMv(mv), defaultMa(ma),
        myConstraintNode(0) 
    {
    }

    void setDefaultNumConstraints(int mp, int mv, int ma) {
        assert(mp >= 0 && mv >= 0 && ma >= 0);
        invalidateTopologyCache();
        defaultMp = mp;
        defaultMv = mv;
        defaultMa = ma;
    }

    typedef std::map<MobilizedBodyIndex,ConstrainedBodyIndex> Mobilized2ConstrainedMap;

    ConstrainedBodyIndex addConstrainedBody(const MobilizedBody& b) {
        assert(isInSameSubsystem(b));
        invalidateTopologyCache();

        const ConstrainedBodyIndex nextIx((int)myConstrainedBodies.size());

        // Add to the Mobilized->Constrained map and check for duplicates.
        std::pair<Mobilized2ConstrainedMap::iterator, bool> result;
        result = myMobilizedBodies.insert(
            Mobilized2ConstrainedMap::value_type(b.getMobilizedBodyIndex(), nextIx));
        assert(result.second); // can only add a body once

        // This is a new constrained body -- add it to the Constrained->Mobilized map too.
        myConstrainedBodies.push_back(b.getMobilizedBodyIndex());
        return nextIx;
    }

    MobilizedBodyIndex getMobilizedBodyIndexOfConstrainedBody(ConstrainedBodyIndex c) const {
        assert(0 <= c && c < (int)myConstrainedBodies.size());
        return myConstrainedBodies[c];
    }

    void realizeTopology(State&,int& nxtQErr, int& nxtUErr, int& nxtMult) const;



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
		getNumConstraintEquations(s, actual_mp, actual_mv, actual_ma);

		bodyForcesInA.resize(getNumConstrainedBodies());       bodyForcesInA  = SpatialVec(Vec3(0), Vec3(0));
        //TODO:
		//mobilityForces.resize(getNumConstrainedMobilities(s)); mobilityForces = 0;

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

    // Find the indicated cache in the passed-in State. This requires that realization has
    // been completed for the associated Stage. During realization, we will instead pass in
    // the appropriate cache entry rather than ask the State for it.
    const SBModelCache&			getModelCache(const State&) const;
    const SBPositionCache&		getPositionCache(const State&) const;
    const SBVelocityCache&		getVelocityCache(const State&) const;
    const SBAccelerationCache&	getAccelerationCache(const State&) const;

    // These are measured from and expressed in the ancestor (A) frame.

    //TODO: should precalculate in State, return reference
    // (Client "get" methods below should be changed to references also.) 

    // These are for use during realization of the associated stage.
    Transform  getBodyTransform   (const State& s, const SBPositionCache&, ConstrainedBodyIndex B) const; // X_AB
    SpatialVec getBodyVelocity    (const State& s, const SBVelocityCache&, ConstrainedBodyIndex B) const; // V_AB
    SpatialVec getBodyAcceleration(const State& s, const SBAccelerationCache&, ConstrainedBodyIndex B) const; // A_AB

    // These are for use when after realization of the associated stage has been completed.
    Transform  getBodyTransform(const State& s, ConstrainedBodyIndex B) const {
        return getBodyTransform(s, getPositionCache(s), B);
    }
    SpatialVec getBodyVelocity(const State& s, ConstrainedBodyIndex B) const {
        return getBodyVelocity(s, getVelocityCache(s), B);
    }
    SpatialVec getBodyAcceleration(const State& s, ConstrainedBodyIndex B) const {
        return getBodyAcceleration(s, getAccelerationCache(s), B);
    }

    // Extract just the rotational quantities from the spatial quantities above.
    //TODO: should be references (see above)
    const Rotation getBodyRotation           (const State& s, const SBPositionCache& pc, ConstrainedBodyIndex B)     const {return getBodyTransform(s,pc,B).R();}   // R_AB
    const Vec3     getBodyAngularVelocity    (const State& s, const SBVelocityCache& vc, ConstrainedBodyIndex B)     const {return getBodyVelocity(s,vc,B)[0];}     // w_AB
    const Vec3     getBodyAngularAcceleration(const State& s, const SBAccelerationCache& ac, ConstrainedBodyIndex B) const {return getBodyAcceleration(s,ac,B)[0];} // b_AB
    const Rotation getBodyRotation           (const State& s, ConstrainedBodyIndex B) const {return getBodyTransform(s,B).R();}   // R_AB
    const Vec3     getBodyAngularVelocity    (const State& s, ConstrainedBodyIndex B) const {return getBodyVelocity(s,B)[0];}     // w_AB
    const Vec3     getBodyAngularAcceleration(const State& s, ConstrainedBodyIndex B) const {return getBodyAcceleration(s,B)[0];} // b_AB

    // Extract just the translational (linear) quantities from the spatial quantities above.
    //TODO: should be references (see above)
    const Vec3 getBodyOriginLocation    (const State& s, const SBPositionCache& pc, ConstrainedBodyIndex B)     const {return getBodyTransform(s,pc,B).T();}   // p_AB
    const Vec3 getBodyOriginVelocity    (const State& s, const SBVelocityCache& vc, ConstrainedBodyIndex B)     const {return getBodyVelocity(s,vc,B)[1];}     // v_AB
    const Vec3 getBodyOriginAcceleration(const State& s, const SBAccelerationCache& ac, ConstrainedBodyIndex B) const {return getBodyAcceleration(s,ac,B)[1];} // a_AB
    const Vec3 getBodyOriginLocation    (const State& s, ConstrainedBodyIndex B) const {return getBodyTransform(s,B).T();}   // p_AB
    const Vec3 getBodyOriginVelocity    (const State& s, ConstrainedBodyIndex B) const {return getBodyVelocity(s,B)[1];}     // v_AB
    const Vec3 getBodyOriginAcceleration(const State& s, ConstrainedBodyIndex B) const {return getBodyAcceleration(s,B)[1];} // a_AB

    Vec3 calcStationLocation(const State& s, const SBPositionCache& pc, ConstrainedBodyIndex B, const Vec3& p_B) const {
        return getBodyTransform(s,pc,B) * p_B; // re-measure and re-express
    }
    Vec3 calcStationVelocity(const State& s, const SBVelocityCache& vc, ConstrainedBodyIndex B, const Vec3& p_B) const {
        const Vec3 p_A = getBodyRotation(s,B) * p_B; // rexpressed but not shifted
        const SpatialVec& V_AB = getBodyVelocity(s,vc,B);
        return V_AB[1] + (V_AB[0] % p_A);
    }
    Vec3 calcStationAcceleration(const State& s, const SBAccelerationCache& ac, ConstrainedBodyIndex B, const Vec3& p_B) const {
        const Vec3  p_A  = getBodyRotation(s,B) * p_B; // rexpressed but not shifted
        const Vec3& w_AB = getBodyVelocity(s,B)[0];
        const SpatialVec& A_AB = getBodyAcceleration(s,ac,B);
        const Vec3 a_A = A_AB[1] + (A_AB[0] % p_A) + w_AB % (w_AB % p_A); // careful: cross product is not associative
        return a_A;
    }

    // These are for use when after realization of the associated stage has been completed.
    Vec3 calcStationLocation(const State& s, ConstrainedBodyIndex B, const Vec3& p_B) const {
        return calcStationLocation(s, getPositionCache(s), B, p_B);
    }
    Vec3 calcStationVelocity(const State& s, ConstrainedBodyIndex B, const Vec3& p_B) const {
        return calcStationVelocity(s, getVelocityCache(s), B, p_B);
    }
    Vec3 calcStationAcceleration(const State& s, ConstrainedBodyIndex B, const Vec3& p_B) const {
        return calcStationAcceleration(s, getAccelerationCache(s), B, p_B);
    }

    // Apply an A-frame force to a B-frame station, updating the appropriate bodyForces entry.
    void addInStationForce(const State& s, ConstrainedBodyIndex B, const Vec3& p_B, 
                           const Vec3& forceInA, Vector_<SpatialVec>& bodyForcesInA) const 
    {
        assert(bodyForcesInA.size() == getNumConstrainedBodies());
        const Rotation& R_AB = getBodyRotation(s,B);
        bodyForcesInA[B] += SpatialVec((R_AB*p_B) % forceInA, forceInA); // rXf, f
    }

    // Apply an A-frame torque to body B, updating the appropriate bodyForces entry.
    void addInBodyTorque(const State& s, ConstrainedBodyIndex B, const Vec3& torqueInA,
                         Vector_<SpatialVec>& bodyForcesInA) const 
    {
        assert(bodyForcesInA.size() == getNumConstrainedBodies());
        bodyForcesInA[B][0] += torqueInA; // no force
    }

    // Apply a generalized (mobility) force to a particular mobility of the given constrained body B,
    // adding it in to the appropriate slot of the mobilityForces vector.
    void addInMobilityForce(const State& s, ConstrainedBodyIndex B, int which, Real f,
                            Vector& mobilityForces) const 
    { 
        assert(mobilityForces.size() == getNumConstrainedMobilities(s));
        assert(0 <= which && which < getNumConstrainedMobilities(s, B));
        mobilityForces[getConstrainedMobilityIndex(s,B,which)] += f;
    }

    virtual ~ConstraintRep();
    virtual ConstraintRep* clone() const = 0;

    // This creates a constraint node using the appropriate constraint type.
    virtual ConstraintNode* createConstraintNode() const = 0; 


    // After realizeTopology() we can look at the values of modeling variables in the State.
    // A Constraint is free to use those in determining how many constraint equations of each
    // type to generate. The default implementation here doesn't look at the state but instead
    // returns the default numbers of equations supplied when the Constraint was constructed.
    void calcNumConstraintEquations(const State& s, int& mp, int& mv, int& ma) const {
        calcNumConstraintEquationsVirtual(s,mp,mv,ma);
    }
    virtual void calcNumConstraintEquationsVirtual(const State&, int& mp, int& mv, int& ma) const {
        mp = defaultMp; mv = defaultMv; ma = defaultMa;
    }

    //NOTE: bodyForces and mobilityForces refer only to constrained bodies and their
    //associated mobilizers, not the system as a whole. They are initialized to zero
    //prior to the call so do not need to be set.

    //NOTE: each of these operators acts on the current state of this constraint's
    //Subtree, which may or may not be the same as that Subtree has in the global
    //State. This is controlled by the base class operator interface methods which
    //will call these only after setting the Subtree state properly.
    //TODO: Subtree

    void realizePositionErrors(const State& s, const SBPositionCache& pc, int mp,  Real* perr) const {
        realizePositionErrorsVirtual(s,pc,mp,perr);
    }
    void realizePositionDotErrors(const State& s, const SBVelocityCache& vc, int mp,  Real* pverr) const {
        realizePositionDotErrorsVirtual(s,vc,mp,pverr);
    }
    void realizePositionDotDotErrors(const State& s, const SBAccelerationCache& ac, int mp,  Real* paerr) const {
        realizePositionDotDotErrorsVirtual(s,ac,mp,paerr);
    }
    void applyPositionConstraintForces
       (const State& s, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const
    {
        applyPositionConstraintForcesVirtual(s,mp,multipliers,bodyForces,mobilityForces);
    }

    void realizeVelocityErrors(const State& s, const SBVelocityCache& vc, int mv,  Real* verr) const {
        realizeVelocityErrorsVirtual(s,vc,mv,verr);
    }
    void realizeVelocityDotErrors(const State& s, const SBAccelerationCache& ac, int mv,  Real* vaerr) const {
        realizeVelocityDotErrorsVirtual(s,ac,mv,vaerr);
    }
    void applyVelocityConstraintForces
       (const State& s, int mv, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const
    {
        applyVelocityConstraintForcesVirtual(s,mv,multipliers,bodyForces,mobilityForces);
    }

    void realizeAccelerationErrors(const State& s, const SBAccelerationCache& ac, int ma,  Real* aerr) const {
        realizeAccelerationErrorsVirtual(s,ac,ma,aerr);
    }
    void applyAccelerationConstraintForces
       (const State& s, int ma, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const
    {
        applyAccelerationConstraintForcesVirtual(s,ma,multipliers,bodyForces,mobilityForces);
    }

    virtual void realizeTopologyVirtual(State&) const { }
    virtual void realizeModelVirtual(State&) const { }
    virtual void realizeInstanceVirtual(const State&) const { }
    virtual void realizeTimeVirtual(const State&) const { }

    // These must be defined if there are any position (holonomic) constraints defined.
    virtual void realizePositionErrorsVirtual      (const State&, const SBPositionCache&, int mp,  Real* perr) const;
    virtual void realizePositionDotErrorsVirtual   (const State&, const SBVelocityCache&, int mp,  Real* pverr) const;
    virtual void realizePositionDotDotErrorsVirtual(const State&, const SBAccelerationCache&, int mp,  Real* paerr) const;
    virtual void applyPositionConstraintForcesVirtual
       (const State&, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const;

    // These must be defined if there are any velocity (nonholonomic) constraints defined.
    virtual void realizeVelocityErrorsVirtual   (const State&, const SBVelocityCache&, int mv,  Real* verr) const;
    virtual void realizeVelocityDotErrorsVirtual(const State&, const SBAccelerationCache&, int mv,  Real* vaerr) const;
    virtual void applyVelocityConstraintForcesVirtual
       (const State&, int mv, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const;

    // These must be defined if there are any acceleration-only constraints defined.
    virtual void realizeAccelerationErrorsVirtual(const State&, const SBAccelerationCache&, int ma,  Real* aerr) const;
    virtual void applyAccelerationConstraintForcesVirtual
       (const State&, int ma, const Real* multipliers,
        Vector_<SpatialVec>& bodyForces,
        Vector&              mobilityForces) const;


    virtual void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
    {
    }

    void calcDecorativeGeometryAndAppend
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
    {
        // Let the individual constraint deal with any complicated stuff.
        calcDecorativeGeometryAndAppendImpl(s,stage,geom);
    }

    void invalidateTopologyCache() const;
    bool subsystemTopologyHasBeenRealized() const;
    const ConstraintNode& getMyConstraintNode() const;

    void setMyMatterSubsystem(SimbodyMatterSubsystem& matter,
                              ConstraintIndex id);

    const SimbodyMatterSubsystem& getMyMatterSubsystem() const;

    bool isInSubsystem() const {
        return myMatterSubsystemRep != 0;
    }

    bool isInSameSubsystem(const MobilizedBody& body) const {
        return isInSubsystem() && body.isInSubsystem() 
            && getMyMatterSubsystem().isSameSubsystem(body.getMatterSubsystem());
    }

    int getNumConstrainedBodies() const {
        SimTK_ASSERT(subsystemTopologyHasBeenRealized(),
            "Number of constrained bodies is not available until Topology stage has been realized.");
        return (int)myConstrainedBodies.size();
    }

    const MobilizedBody& getConstrainedMobilizedBody(ConstrainedBodyIndex B) const;
    const MobilizedBody& getAncestorMobilizedBody() const;

	// Find out how many holonomic (position), nonholonomic (velocity),
	// and acceleration-only constraint equations are generated by this Constraint.
	// State must be realized to Stage::Model.
	void getNumConstraintEquations(const State& s, int& mp, int& mv, int& ma) const;

    int getNumConstrainedMobilities(const State& s) const {
        //TODO
        assert(!"Constraint::getNumConstrainedMobilities() not implemented yet.");
        return -1;
    }

    int getNumConstrainedMobilities(const State& s, ConstrainedBodyIndex B) const {
        //TODO
        assert(!"Constraint::getNumConstrainedMobilities(B) not implemented yet.");
        return -1;
    }

    int getConstrainedMobilityIndex(const State& s, ConstrainedBodyIndex B, int which) const {
        //TODO
        assert(!"Constraint::getConstrainedMobilityIndex(B) not implemented yet.");
        return -1;
    }

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

    void setMyHandle(Constraint& h) {myHandle = &h;}
    const Constraint& getMyHandle() const {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}

private:
    friend class Constraint;
    Constraint* myHandle;	// the owner handle of this rep

        // TOPOLOGY "STATE"

    // These data members are filled in once the Constraint is added to
    // a MatterSubsystem.
    SimbodyMatterSubsystemRep* myMatterSubsystemRep;
    ConstraintIndex            myConstraintIndex; // id within the matter subsystem

    // We'll keep the bodies in two maps: one maps MobilizedBodyIndex->ConstrainedBodyIndex
    // (O(log n) to look up), and the other maps ConstrainedBodyIndex->MobilizedBodyIndex
    // (randomly addressable in constant time).
    Mobilized2ConstrainedMap        myMobilizedBodies;
    std::vector<MobilizedBodyIndex> myConstrainedBodies; // index with ConstrainedBodyIndex

    // These are the defaults for the number of position (holonomic) constraint equations,
    // the number of velocity (nonholonomic) constraint equations, and the number of
    // acceleration-only constraint equations.
    int defaultMp, defaultMv, defaultMa;

        // TOPOLOGY "CACHE"

    // When topology is realized we study the constrained bodies to identify the
    // subtree of mobilized bodies which may be kinematically involved in satisfaction
    // of this Constraint. This requires finding the outmost common ancestor of 
    // the constrained bodies. All mobilized bodies on the paths inward from the
    // constrained bodies to the ancestor are included; nothing outboard of the
    // constrained bodies is included; and the ancestor is treated as ground so
    // that its mobilities are *not* included.
    mutable SimbodyMatterSubsystem::Subtree mySubtree;

    // A ConstraintNode is created for each Constraint during realizeTopology().
    // Think of it as the *computational* form of the Constraint; the Constraint
    // itself exists to make a nice API for the programmer. In fact is is possible
    // for several different Constraint classes to use the same underlying computational
    // form while providing a different API.
    //
    // This is a pointer to an abstract object whose heap space is *owned* by
    // this Constraint. Be sure to delete it upon destruction and whenever
    // topology is re-realized.
    mutable ConstraintNode* myConstraintNode;
};


class Constraint::Rod::RodRep : public Constraint::ConstraintRep {
public:
    RodRep() 
      : ConstraintRep(1,0,0), defaultPoint1(0), defaultPoint2(0), defaultRodLength(1)
    { 
        // Rod constructor sets all the data members here directly
    }
    RodRep* clone() const { return new RodRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    // Implementation of virtuals required for holonomic constraints.

    // perr = (p^2 - d^2)/2
    void realizePositionErrorsVirtual(const State& s, const SBPositionCache& pc, int mp,  Real* perr) const {
        assert(mp==1 && perr);
        const Vec3 p1 = calcStationLocation(s, pc, B1, defaultPoint1); // meas from & expr in ancestor
        const Vec3 p2 = calcStationLocation(s, pc, B2, defaultPoint2);
        const Vec3 p = p2 - p1;
        //TODO: save p in state

        *perr = (dot(p, p) - square(defaultRodLength)) / 2;
    }

    // pverr = d/dt perr = pdot*p = v*p, where v=v2-v1 is relative velocity
    void realizePositionDotErrorsVirtual(const State& s, const SBVelocityCache& vc, int mp,  Real* pverr) const {
        assert(mp==1 && pverr);
        //TODO: should be able to get p from State
        const Vec3 p1 = calcStationLocation(s, B1, defaultPoint1); // meas from & expr in ancestor
        const Vec3 p2 = calcStationLocation(s, B2, defaultPoint2);
        const Vec3 p = p2 - p1;

        const Vec3 v1 = calcStationVelocity(s, vc, B1, defaultPoint1); // meas & expr in ancestor
        const Vec3 v2 = calcStationVelocity(s, vc, B2, defaultPoint2);
        const Vec3 v = v2 - v1;
        *pverr = dot(v, p);
    }

    // paerr = d/dt verr = vdot*p + v*pdot =a*p+v*v, where a=a2-a1 is relative acceleration
    void realizePositionDotDotErrorsVirtual(const State& s, const SBAccelerationCache& ac, int mp,  Real* paerr) const {
        assert(mp==1 && paerr);
        //TODO: should be able to get p and v from State
        const Vec3 p1 = calcStationLocation(s, B1, defaultPoint1); // meas from & expr in ancestor
        const Vec3 p2 = calcStationLocation(s, B2, defaultPoint2);
        const Vec3 p = p2 - p1;
        const Vec3 v1 = calcStationVelocity(s, B1, defaultPoint1); // meas & expr in ancestor
        const Vec3 v2 = calcStationVelocity(s, B2, defaultPoint2);
        const Vec3 v = v2 - v1;

        const Vec3 a1 = calcStationAcceleration(s, ac, B1, defaultPoint1); // meas & expr in ancestor
        const Vec3 a2 = calcStationAcceleration(s, ac, B2, defaultPoint2);
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
        const Vec3 p1 = calcStationLocation(s, B1, defaultPoint1); // meas from & expr in ancestor
        const Vec3 p2 = calcStationLocation(s, B2, defaultPoint2);
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

    SimTK_DOWNCAST(RodRep, ConstraintRep);
private:
    friend class Constraint::Rod;

    ConstrainedBodyIndex B1, B2;

    Vec3            defaultPoint1; // on body 1, exp. in B1 frame
    Vec3            defaultPoint2; // on body 2, exp. in B2 frame
    Real            defaultRodLength;
};


class Constraint::PointInPlane::PointInPlaneRep : public Constraint::ConstraintRep {
public:
    PointInPlaneRep()
      : ConstraintRep(1,0,0), defaultPlaneNormal(), defaultPlaneHeight(0), defaultFollowerPoint(0),
        planeHalfWidth(1), pointRadius(0.05) 
    { }
    PointInPlaneRep* clone() const { return new PointInPlaneRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const;

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
	// at point S and the coincident point C on B which is instantaneously at the same spatial
	// location as S. Then n is the plane normal (a constant unit vector in B), h is the
	// plane height measured from the B origin along n (a scalar constant),
    // and h_C is the current height of material point C (on B) which is by construction
	// the same as the height of S (on F) over the plane (on B). Point C's location in B
	// is given by the vector p_BC from B's origin to the current location of S, and expressed
	// in B. That vector expressed in A is p_BC_A (= p_AS-p_AB).
	//
	// h_C is a calculated, state dependent quantity. The constraint we want to enforce is
	// that h_C(q)=h. We will work in the B frame where the plane normal n and its
	// height h are constant.
	//
	// Except for the projection along the normal, this constraint is identical to the ball
	// joint constraint (which has 3 constraints instead of the 1 here).
    //
    //    perr = h_C - h
    //         = ~p_BC*n - h
	//    --------------------------------
    //    perr = ~[R_BA*(p_AS-p_AB)]*n - h
	//    --------------------------------
    //
    //   (Below we're using the identity w_BA = -R_BA*w_AB, and the scalar triple
	//    product identity v*(w X z)=w*(z X v)=z*(v X w).)
    //
    //    verr = d/dt perr = ~v_BC * n (derivative taken in B)
    //         = ~[w_BA X R_BA*(p_AS-p_AB) 
    //                  + R_BA*(v_AS-v_AB)] * n
	//    -----------------------------------------
    //    verr = ~[R_BA*(         (v_AS-v_AB) 
    //                   - w_AB X (p_AS-p_AB))] * n
	//    -----------------------------------------
    //
    //    aerr = d/dt verr = ~a_BC * n (derivative taken in B)
    //         = ~[            b_BA X R_BA*(p_AS-p_AB) 
    //               + w_BA X (w_BA X R_BA*(p_AS-p_AB))
    //               +       2 w_BA X R_BA*(v_AS-v_AB)
    //               +                R_BA*(a_AS-a_AB) ] * n
    //
    //         = ~[              (-R_BA*b_AB) X R_BA*(p_AS-p_AB) 
    //               + R_BA*w_AB X (R_BA*w_AB X R_BA*(p_AS-p_AB))
    //               -            2 R_BA*w_AB X R_BA*(v_AS-v_AB)
    //               +                          R_BA*(a_AS-a_AB) ] * n
	//    ---------------------------------------------------------
    //    aerr = ~[R_BA*(         -b_AB X (p_AS-p_AB) 
    //                   + w_AB X (w_AB X (p_AS-p_AB))
    //                   -       2 w_AB X (v_AS-v_AB)
    //                   +                (a_AS-a_AB) ] * n
	//    ---------------------------------------------------------
    //  
    // Then, from examination of verr and some rearrangement,
	// we find velocities and angular velocities used like this:
    //       ~v_AS*n_A                  (body F at point S) 
	//     - ~v_AB*n_A                  (body B at B origin)
	//       ~w_AB*(n_A X (p_AS-p_AB))  (applied to body B)
    //
    // so we apply a forces lambda*n_A to S, -lambda*n_A to OB, and torque lambda*(n_A X (p_AS-p_AB)) to B.
    // More simply, just apply force lambda*n_A to S and -lambda*n_A to the point C of B which is
    // instantaneously coincident in space to S. That adds an rXF torque (p_AS-p_AB) X (-lambda*n_A).
    //

	//    --------------------------------
    //    perr = ~[R_BA*(p_AS-p_AB)]*n - h
	//    --------------------------------
    void realizePositionErrorsVirtual(const State& s, const SBPositionCache& pc, int mp,  Real* perr) const {
        assert(mp==1 && perr);

        const Rotation&  R_AB   = getBodyRotation(s, pc, planeBody);
        const Vec3&      p_AB   = getBodyOriginLocation(s, pc, planeBody);
        const Vec3       p_AS   = calcStationLocation(s, pc, followerBody, defaultFollowerPoint);
        const Vec3       p_BS_A = p_AS - p_AB;

        *perr = dot(~R_AB*p_BS_A, defaultPlaneNormal) - defaultPlaneHeight;
    }

	//    -----------------------------------------
    //    verr = ~[R_BA*(         (v_AS-v_AB) 
    //                   - w_AB X (p_AS-p_AB))] * n
	//    -----------------------------------------
    void realizePositionDotErrorsVirtual(const State& s, const SBVelocityCache& vc, int mp,  Real* pverr) const {
        assert(mp==1 && pverr);
        //TODO: should be able to get p info from State
        const Rotation&   R_AB   = getBodyRotation(s, planeBody);
        const Vec3&       p_AB   = getBodyOriginLocation(s, planeBody);
        const Vec3        p_AS   = calcStationLocation(s, followerBody, defaultFollowerPoint);
        const Vec3        p_BS_A = p_AS - p_AB;

        const Vec3&       w_AB    = getBodyAngularVelocity(s, vc, planeBody);
        const Vec3&       v_AB    = getBodyOriginVelocity(s, vc, planeBody);
        const Vec3        v_AS    = calcStationVelocity(s, vc, followerBody, defaultFollowerPoint);
        const Vec3        v_BS_A  = v_AS - v_AB;
        const Vec3        pverr_A = v_BS_A - w_AB % p_BS_A;

        *pverr = dot(~R_AB * pverr_A, defaultPlaneNormal);
    }

	//    ---------------------------------------------------------
    //    aerr = ~[R_BA*(         -b_AB X (p_AS-p_AB) 
    //                   + w_AB X (w_AB X (p_AS-p_AB))
    //                   -       2 w_AB X (v_AS-v_AB)
    //                   +                (a_AS-a_AB) ] * n
	//    ---------------------------------------------------------
    void realizePositionDotDotErrorsVirtual(const State& s, const SBAccelerationCache& ac, int mp,  Real* paerr) const {
        assert(mp==1 && paerr);
        //TODO: should be able to get p and v info from State

        const Rotation&   R_AB   = getBodyRotation(s, planeBody);
        const Vec3&       p_AB   = getBodyOriginLocation(s, planeBody);
        const Vec3        p_AS   = calcStationLocation(s, followerBody, defaultFollowerPoint);
        const Vec3        p_BS_A = p_AS - p_AB;

        const Vec3&       w_AB    = getBodyAngularVelocity(s, planeBody);
        const Vec3&       v_AB    = getBodyOriginVelocity(s, planeBody);
        const Vec3        v_AS    = calcStationVelocity(s, followerBody, defaultFollowerPoint);
        const Vec3        v_BS_A  = v_AS - v_AB;

        const Vec3&       b_AB    = getBodyAngularAcceleration(s, ac, planeBody);
        const Vec3&       a_AB    = getBodyOriginAcceleration(s, ac, planeBody);
        const Vec3        a_AS    = calcStationAcceleration(s, ac, followerBody, defaultFollowerPoint);
        const Vec3        a_BS_A  = a_AS - a_AB;

        const Vec3        paerr_A = a_BS_A - b_AB % p_BS_A
                                    - 2 * w_AB % v_BS_A
                                    + w_AB % (w_AB % p_BS_A);

        *paerr = dot(~R_AB * paerr_A, defaultPlaneNormal);
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
        const Rotation&  R_AB = getBodyRotation(s, planeBody);
        const Vec3&      p_AB = getBodyOriginLocation(s, planeBody);
        const Vec3&      p_FS = defaultFollowerPoint;
        const Vec3       p_AS = calcStationLocation(s, followerBody, defaultFollowerPoint);
        const Vec3       p_BC = ~R_AB * (p_AS-p_AB); // shift to B origin and reexpress in B;
                                                     // C is material point of B coincident with S

		const Vec3 force_B = lambda*defaultPlaneNormal;
		const Vec3 force_A = R_AB * force_B;

        addInStationForce(s, followerBody, p_FS,  force_A, bodyForcesInA);
        addInStationForce(s, planeBody,    p_BC, -force_A, bodyForcesInA);
    }

    SimTK_DOWNCAST(PointInPlaneRep, ConstraintRep);
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


class Constraint::ConstantAngle::ConstantAngleRep : public Constraint::ConstraintRep {
public:
    ConstantAngleRep()
      : ConstraintRep(1,0,0), defaultAxisB(), defaultAxisF(), defaultAngle(Pi/2),
        axisLength(1), axisThickness(1), cosineOfDefaultAngle(NaN)
    { }
    ConstantAngleRep* clone() const { return new ConstantAngleRep(*this); }

    ConstraintNode* createConstraintNode() const;

    void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const;

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
    //      = ~(b_AF-b_AB) * (f_A % b_A)
    //        + 2 (w_AF % w_AB) * (f_A % b_A)
    // => ------------------------------------------------
    // aerr = ~(b_AF - b_AB + 2 w_AF % w_AB) * (f_A % b_A)
    // ---------------------------------------------------
    //
    // Constraint torque can be determined by inspection of verr:
    //    lambda * (f_A % b_A) applied to body F
    //   -lambda * (f_A % b_A) applied to body B
    //

    // ------------------------------
    // perr = ~b_A * f_A - cos(theta)
    // ------------------------------
    void realizePositionErrorsVirtual(const State& s, const SBPositionCache& pc, int mp,  Real* perr) const {
        assert(mp==1 && perr);

        const Rotation& R_AB = getBodyRotation(s, pc, B);
        const Rotation& R_AF = getBodyRotation(s, pc, F);
        const UnitVec3  b_A  = R_AB * defaultAxisB;
        const UnitVec3  f_A  = R_AF * defaultAxisF;

        *perr = dot(b_A, f_A) - cosineOfDefaultAngle;
    }

    // ----------------------------------
    // pverr = ~(w_AF-w_AB) * (f_A % b_A)
    // ----------------------------------
    void realizePositionDotErrorsVirtual(const State& s, const SBVelocityCache& vc, int mp,  Real* pverr) const {
        assert(mp==1 && pverr);
        //TODO: should be able to get p info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const UnitVec3  b_A  = R_AB * defaultAxisB;
        const UnitVec3  f_A  = R_AF * defaultAxisF;
        const Vec3      w_AB = getBodyAngularVelocity(s, vc, B);
        const Vec3      w_AF = getBodyAngularVelocity(s, vc, F);

        *pverr = dot( w_AF-w_AB,  f_A % b_A );
    }

    // ----------------------------------------------------
    // paerr = ~(b_AF - b_AB + 2 w_AF % w_AB) * (f_A % b_A)
    // ----------------------------------------------------
    void realizePositionDotDotErrorsVirtual(const State& s, const SBAccelerationCache& ac, int mp,  Real* paerr) const {
        assert(mp==1 && paerr);
        //TODO: should be able to get p and v info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const UnitVec3  b_A  = R_AB * defaultAxisB;
        const UnitVec3  f_A  = R_AF * defaultAxisF;
        const Vec3      w_AB = getBodyAngularVelocity(s, B);
        const Vec3      w_AF = getBodyAngularVelocity(s, F);
        const Vec3      b_AB = getBodyAngularAcceleration(s, ac, B);
        const Vec3      b_AF = getBodyAngularAcceleration(s, ac, F);

        *paerr = dot( b_AF-b_AB + 2*(w_AF % w_AB), 
                      f_A % b_A );
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

    SimTK_DOWNCAST(ConstantAngleRep, ConstraintRep);
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


class Constraint::Ball::BallRep : public Constraint::ConstraintRep {
public:
    BallRep() : ConstraintRep(3,0,0), defaultPoint1(0), defaultPoint2(0), defaultRadius(0.1) { }
    BallRep* clone() const { return new BallRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const;

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

    void realizePositionErrorsVirtual(const State& s, const SBPositionCache& pc, int mp,  Real* perr) const {
        assert(mp==3 && perr);

        const Vec3 p_AP = calcStationLocation(s, pc, B1, defaultPoint1);
        const Vec3 p_AS = calcStationLocation(s, pc, B2, defaultPoint2);

        // See above comments -- this is just the constant of integration; there is a missing (p_AS-p_AC)
        // term (always 0) here which is what we differentiate to get the verr equation.
        Vec3::updAs(perr) = p_AS - p_AP;
    }
 
    void realizePositionDotErrorsVirtual(const State& s, const SBVelocityCache& vc, int mp,  Real* pverr) const {
        assert(mp==3 && pverr);
        //TODO: should be able to get p info from State
        const Transform&  X_AB   = getBodyTransform(s, B1);
        const Vec3        p_AS   = calcStationLocation(s, B2, defaultPoint2);
        const Vec3        p_BC   = ~X_AB*p_AS; // C is a material point of body B

        const Vec3        v_AS    = calcStationVelocity(s, vc, B2, defaultPoint2);
        const Vec3        v_AC    = calcStationVelocity(s, vc, B1, p_BC);
        Vec3::updAs(pverr) = v_AS - v_AC;
    }

    void realizePositionDotDotErrorsVirtual(const State& s, const SBAccelerationCache& ac, int mp,  Real* paerr) const {
        assert(mp==3 && paerr);
        //TODO: should be able to get p and v info from State

        const Transform&  X_AB   = getBodyTransform(s, B1);
        const Vec3        p_AS   = calcStationLocation(s, B2, defaultPoint2);
        const Vec3        p_BC   = ~X_AB*p_AS; // C is a material point of body B

        const Vec3        a_AS    = calcStationAcceleration(s, ac, B2, defaultPoint2);
        const Vec3        a_AC    = calcStationAcceleration(s, ac, B1, p_BC);
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
        const Vec3       p_AS  = calcStationLocation(s, B2, p_FS);
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

    SimTK_DOWNCAST(BallRep, ConstraintRep);
private:
    friend class Constraint::Ball;

    ConstrainedBodyIndex B1;
    ConstrainedBodyIndex B2;

    Vec3            defaultPoint1; // on body 1, exp. in B1 frame
    Vec3            defaultPoint2; // on body 2, exp. in B2 frame
    Real            defaultRadius; // used for visualization only
};

class Constraint::ConstantOrientation::ConstantOrientationRep : public Constraint::ConstraintRep {
public:
    ConstantOrientationRep()
      : ConstraintRep(1,0,0), defaultRB(), defaultRF(),
        axisLength(1), axisThickness(1)
    { }
    ConstantOrientationRep* clone() const { return new ConstantOrientationRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const;

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
    // ---------------------------------------------------
    // aerr = ~(b_AF - b_AB + 2 w_AF % w_AB) * (RFx % RBy)
    //      = ~(b_AF - b_AB + 2 w_AF % w_AB) * (RFy % RBz)
    //      = ~(b_AF - b_AB + 2 w_AF % w_AB) * (RFz % RBx)
    // ---------------------------------------------------
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
    void realizePositionErrorsVirtual(const State& s, const SBPositionCache& pc, int mp,  Real* perr) const {
        assert(mp==3 && perr);

        const Rotation& R_AB = getBodyRotation(s, pc, B);
        const Rotation& R_AF = getBodyRotation(s, pc, F);
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
    void realizePositionDotErrorsVirtual(const State& s, const SBVelocityCache& vc, int mp,  Real* pverr) const {
        assert(mp==3 && pverr);
        //TODO: should be able to get p info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const Rotation  RB = R_AB * defaultRB; // now expressed in A
        const Rotation  RF = R_AF * defaultRF;

        const Vec3&     w_AB = getBodyAngularVelocity(s, vc, B);
        const Vec3&     w_AF = getBodyAngularVelocity(s, vc, F);
        const Vec3      w_BF = w_AF-w_AB; // in A

        Vec3::updAs(pverr) = Vec3( ~w_BF * (RF.x() % RB.y()),
                                   ~w_BF * (RF.y() % RB.z()),
                                   ~w_BF * (RF.z() % RB.x()) );
    }

    // ----------------------------------------------------
    // aerr = ~(b_AF - b_AB + 2 w_AF % w_AB) * (RFx % RBy)
    //      = ~(b_AF - b_AB + 2 w_AF % w_AB) * (RFy % RBz)
    //      = ~(b_AF - b_AB + 2 w_AF % w_AB) * (RFz % RBx)
    // ----------------------------------------------------
    void realizePositionDotDotErrorsVirtual(const State& s, const SBAccelerationCache& ac, int mp,  Real* paerr) const {
        assert(mp==3 && paerr);
        //TODO: should be able to get p and v info from State
        const Rotation& R_AB = getBodyRotation(s, B);
        const Rotation& R_AF = getBodyRotation(s, F);
        const Rotation  RB = R_AB * defaultRB; // now expressed in A
        const Rotation  RF = R_AF * defaultRF;

        const Vec3&     w_AB = getBodyAngularVelocity(s, B);
        const Vec3&     w_AF = getBodyAngularVelocity(s, F);

        const Vec3&     b_AB = getBodyAngularAcceleration(s, ac, B);
        const Vec3&     b_AF = getBodyAngularAcceleration(s, ac, F);

        const Vec3      b_BF = (b_AF-b_AB + 2*(w_AF % w_AB)); // in A

        Vec3::updAs(paerr) = Vec3( ~b_BF * (RF.x() % RB.y()),
                                   ~b_BF * (RF.y() % RB.z()),
                                   ~b_BF * (RF.z() % RB.x()) );
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

    SimTK_DOWNCAST(ConstantOrientationRep, ConstraintRep);
private:
    friend class Constraint::ConstantOrientation;

    ConstrainedBodyIndex B; // B1 is "base" body
    ConstrainedBodyIndex F; // B2 is "follower" body

    Rotation          defaultRB; // fixed to B, expressed in B frame; RB = R_B_RB
    Rotation          defaultRF; // fixed to F, expressed in F frame; RF = R_F_RF

    // These are just for visualization
    Real axisLength;
    Real axisThickness;
};


class Constraint::Weld::WeldRep : public Constraint::ConstraintRep {
public:
    WeldRep() : ConstraintRep(6,0,0) { } // transforms are Identity
    WeldRep* clone() const { return new WeldRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    // Implementation of virtuals required for holonomic constraints.

    // TODO: THEORY GOES HERE

    void realizePositionErrorsVirtual(const State& s, const SBPositionCache& pc, int mp,  Real* perr) const {
        assert(mp==6 && perr);

        Vec6::updAs(perr) = 0; // orientation error, position error
    }

    // pverr = d/dt perr = 
    void realizePositionDotErrorsVirtual(const State& s, const SBVelocityCache& vc, int mp,  Real* pverr) const {
        assert(mp==6 && pverr);
        //TODO: should be able to get p info from State


        Vec6::updAs(pverr) = 0;
    }

    // paerr = d/dt verr = 
    void realizePositionDotDotErrorsVirtual(const State& s, const SBAccelerationCache& ac, int mp,  Real* paerr) const {
        assert(mp==6 && paerr);
        //TODO: should be able to get p and v info from State

        // XXX

        Vec6::updAs(paerr) = 0;
    }

    void applyPositionConstraintForcesVirtual
       (const State& s, int mp, const Real* multipliers,
        Vector_<SpatialVec>& bodyForcesInA,
        Vector&              mobilityForces) const
    {
        assert(mp==6 && multipliers);
        const Vec3 torque = Vec3::getAs(multipliers);
        const Vec3 force  = Vec3::getAs(multipliers + 3);
        //TODO: should be able to get p info from State

        const Transform& X_AB1  = getBodyTransform(s,B1);
        const Vec3&      p_B2F2 = defaultFrame2.T();
        const Vec3       p_AF2  = calcStationLocation(s, B2, p_B2F2);
        const Vec3       p_B1F2 = ~X_AB1 * p_AF2; // shift to B1 origin and reexpress in B1

        // Multipliers are torque and force to be applied to body2; apply
        // equal and opposite torque to body1, but apply the -force to the
        // point of body1 coincident with frame2's origin, which in general
        // won't be exactly the same as frame1's origin.

        addInBodyTorque(s, B2,  torque, bodyForcesInA);
        addInBodyTorque(s, B1, -torque, bodyForcesInA);

        addInStationForce(s, B2, p_B2F2,  force, bodyForcesInA);
        addInStationForce(s, B1, p_B1F2, -force, bodyForcesInA);
    }

    SimTK_DOWNCAST(WeldRep, ConstraintRep);
private:
    friend class Constraint::Weld;

    ConstrainedBodyIndex B1;
    ConstrainedBodyIndex B2;

    Transform       defaultFrame1; // on body 1, relative to B1 frame
    Transform       defaultFrame2; // on body 2, relative to B2 frame};
};


class Constraint::Custom::CustomRep : public Constraint::ConstraintRep {
public:
    CustomRep* clone() const { return new CustomRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    SimTK_DOWNCAST(CustomRep, ConstraintRep);
private:
    friend class Constraint::Custom;

    // TODO: God only knows what goes here!
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_REP_H_



