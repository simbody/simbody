#ifndef SimTK_SIMBODY_CONDITIONAL_CONSTRAINT_H_
#define SimTK_SIMBODY_CONDITIONAL_CONSTRAINT_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"

namespace SimTK {

class MobilizedBody;

/** TODO: Simbody model element representing a conditionally-enforced 
constraint.
**/
class SimTK_SIMBODY_EXPORT ConditionalConstraint {
public:

/** Given the specified minimum coefficient of restitution (COR), capture speed,
and the speed at which the minimum COR is attained, calculate an effective COR 
for a given impact speed. All speeds must be nonnegative. The COR will be zero
at or below capture speed, will be minCOR at or above minCORSpeed, and will 
rise linearly with decreasing impact speed between minCORSpeed and the capture
velocity. **/
static Real calcEffectiveCOR(Real minCOR, Real captureSpeed, Real minCORSpeed, 
                             Real impactSpeed)
{
    assert(0 <= minCOR && minCOR <= 1);
    assert(0 <= captureSpeed && captureSpeed <= minCORSpeed);
    assert(impactSpeed >= 0);

    if (impactSpeed <= captureSpeed) return 0;
    if (impactSpeed >= minCORSpeed)  return minCOR;
    // captureSpeed < impactSpeed < minCORSpeed
    const Real slope = (1-minCOR) / minCORSpeed;
    return 1 - slope*impactSpeed;
}

/** Given the coefficients of friction and slip-to-rolling transition 
speed, calculate the effective COF mu for a given slip velocity. Speeds
must be nonnegative. This utility calculates mu with an abrupt rise at
the transitionSpeed from the dynamic coefficient mu_d to the static 
coefficient mu_s. There is also a linear contribution mu_v*slipSpeed. **/
static Real calcEffectiveCOF(Real mu_s, Real mu_d, Real mu_v,
                             Real transitionSpeed, Real slipSpeed)
{
    assert(mu_s>=0 && mu_d>=0 && mu_v>=0);
    assert(mu_s>=mu_d);
    assert(transitionSpeed>=0 && slipSpeed>=0);
    const Real viscous = mu_v*slipSpeed; // typically zero
    return viscous + (slipSpeed <= transitionSpeed ? mu_s : mu_d);
}

};

class UnilateralContact; // normal + friction
class UnilateralSpeedConstraint;
class BoundedSpeedConstraint;

class ConstraintLimitedFriction;
class StateLimitedFriction;


/** A unilateral constraint uses a single constraint equation for which the
associated Lagrange multiplier lambda (or impulse pi) is restricted to be
nonpositive. Since in Simbody multipliers have the opposite sense from applied
forces, this is a requirement that the constraint force be nonnegative.

That is, we must satisfy the linear complementarity constraints 
<pre>
    lambda_i <= 0
    err_i    <= 0, err_i = A_i lambda - verr_i
    lambda_i*err_i = 0
</pre>
Note that we require constraints to be defined such that a constraint force 
-lambda generates a positive acceleration at this constraint. That is the same
as requiring that the diagonals of A=GM\~G are nonnegative.
**/
class SimTK_SIMBODY_EXPORT UnilateralConstraint {
public:
    UnilateralConstraint() {}
    virtual ~UnilateralConstraint() {}

    /** Get the multiplier associated with this constraint. **/
    MultiplierIndex getMultiplierIndex(const State& state) const;
private:
};

//==============================================================================
//                       UNILATERAL CONTACT CONSTRAINT
//==============================================================================
/** A unilateral contact constraint uses a single holonomic (position) 
constraint equation to prevent motion in one direction while leaving it 
unrestricted in the other direction. Examples are surface-surface contact, joint
stops, and inextensible ropes. These constraints are subject to violent impacts
that are treated with a coefficient of restitution that may be state dependent.

Some unilateral contacts may be associated with one or more friction elements
that are dependent on the normal force generated by the contact. Whenever the
unilateral contact is inactive (meaning its associated multiplier is zero), all
the associated friction elements are also inactive.
**/
class SimTK_SIMBODY_EXPORT UnilateralContact
:   public UnilateralConstraint {
public:
    UnilateralContact() {}

    /** Disable the normal and friction constraints if they were enabled. 
    Return true if we actually had to disable something. **/
    virtual bool disable(State& state) const = 0;

    /** Enable the normal and friction constraints if they were disabled. 
    Return true if we actually had to enable something. **/
    virtual bool enable(State& state) const = 0;

    /** Returns the effective coefficient of restitution (COR) for this contact,
    given an impact speed (a nonnegative scalar). For a given pair of contacting 
    materials this is typically a function of just the impact speed, but it
    could also depend on the time and configuration in \a state, which should 
    be realized through Stage::Position. The given default impact speed
    thresholds (also nonnegative) are used to calculate the COR unless this 
    Contact overrides those. **/
    virtual Real calcEffectiveCOR(const State& state,
                                  Real defaultCaptureSpeed,
                                  Real defaultMinCORSpeed,
                                  Real impactSpeed) const = 0;

    /** Return the position error (usually a signed distance function). **/
    virtual Real getPerr(const State& state) const = 0;
    /** Return the time derivative of the position error. **/
    virtual Real getVerr(const State& state) const = 0;
    /** Return the time derivative of the velocity error. **/
    virtual Real getAerr(const State& state) const = 0;
    /** Given the position constraint tolerance currently in use, is this 
    contact close enough to contacting that we should treat it as though
    it is in contact? Normally we just see if perr <= tol, but individual
    contacts can override this if they want to do some scaling. **/
    virtual bool isProximal(const State& state, Real ptol) const 
    {   return getPerr(state) <= ptol; }

    /** Return the multiplier index Simbody assigned for the unilateral 
    contact constraint (for contact, this is the normal constraint). If the
    constraint is not enabled, there is no multiplier and the returned index
    will be invalid. **/
    virtual MultiplierIndex 
    getContactMultiplierIndex(const State& state) const = 0;

    /** Returns \c true if there is a friction constraint associated with this
    contact constraint. If so, calcEffectiveCOF() must be overridden. **/
    virtual bool hasFriction(const State& state) const {return false;}

    /** Returns the effective coefficient of friction mu for this contact,
    given a relative slip speed (a nonnegative scalar). For a given pair of 
    contacting materials this is typically a function of just the slip speed, 
    but it could also depend on the time and configuration in \a state, which 
    should be realized through Stage::Position. The given default slip-to-roll
    transition speed threshold (also nonnegative) is used to calculate mu 
    unless this Contact overrides it with its own transition speed. **/
    virtual Real calcEffectiveCOF(const State& state,
                                  Real defaultTransitionSpeed,
                                  Real slipSpeed) const
    {   return NaN; }

    virtual Vec2 getSlipVelocity(const State& state) const 
    {   return Vec2(NaN); }

    /** If hasFriction(), this method returns the multipliers used for the
    x- and y-direction friction constraints. If no friction, or if this
    constraint is disabled, the returned values are invalid. **/
    virtual void
    getFrictionMultiplierIndices(const State&       state, 
                                 MultiplierIndex&   ix_x, 
                                 MultiplierIndex&   ix_y) const
    {   ix_x.invalidate(); ix_y.invalidate(); }

    /** TODO: kludge needed because we're misusing existing constraints. 
    This must be called while Stage::Position is valid. **/
    virtual Vec3 getPositionInfo(const State& state) const 
    {   return Vec3(NaN); }
    /** TODO: kludge to set instance parameters on internal constraints;
    this should be the same Vec3 you got from getPositionInfo(). **/
    virtual void setInstanceParameter(State& state, const Vec3& pos) const {}


    void setMyIndex(UnilateralContactIndex cx) {m_myIx = cx;}
    UnilateralContactIndex getMyIndex() const {return m_myIx;}
private:
    UnilateralContactIndex m_myIx;
};

//==============================================================================
//                       UNILATERAL SPEED CONSTRAINT
//==============================================================================
/** A unilateral speed constraint uses a single nonholonomic (velocity)
constraint equation to prevent relative slip in one direction but not in the 
other. Examples are ratchets and mechanical diodes.
**/
class SimTK_SIMBODY_EXPORT UnilateralSpeedConstraint
:   public UnilateralConstraint {
public:
    UnilateralSpeedConstraint() {}
    virtual ~UnilateralSpeedConstraint() {}

private:
};

//==============================================================================
//                         BOUNDED SPEED CONSTRAINT
//==============================================================================
/** A bounded speed constraint uses a single nonholonomic (velocity) constraint
equation to prevent relative slip provided it can do so while keeping the
multiplier value within a range given by a lower and upper bound. Outside that
range the connection will slip and the multiplier's value will be one of the
bounds, depending on the slip direction. An example is a torque-limited 
speed control motor. We enforce:
<pre>
    lower <= -lambda <= upper and verr=0
    or verr > 0 and -lambda=lower
    or verr < 0 and -lambda=upper
</pre>
The bounds (lower,upper) can be state dependent, for example, they may be
dependent on the current slip velocity. When lower=-upper, this is just 
a restriction on the magnitude |lambda|, like a friction constraint where the
normal force is known.

This constraint is workless when it is able to prevent slip with the multiplier
in range; it is maximally dissipative otherwise because the constraint force
opposes the slip velocity.
**/
class SimTK_SIMBODY_EXPORT BoundedSpeedConstraint {
public:
    BoundedSpeedConstraint() {}
    virtual ~BoundedSpeedConstraint() {}

    /** Return the currently effective lower and upper bounds on the
    associated multiplier as a Vec2(lower,upper). The bounds may depend on
    time, position, and velocity taken from the given \a state.
    **/
    virtual Vec2 calcEffectiveBounds(const State& state) const = 0;
private:
};

//==============================================================================
//                          STATE LIMITED FRICTION
//==============================================================================
class SimTK_SIMBODY_EXPORT StateLimitedFriction {
public:
    StateLimitedFriction() {}
    virtual ~StateLimitedFriction() {}

    /** Disable the friction constraints if they were enabled. Return true if 
    we actually had to disable something. **/
    virtual bool disable(State& state) const = 0;

    /** Enable the friction constraints if they were disabled. Return true if 
    we actually had to enable something. **/
    virtual bool enable(State& state) const = 0;

    /** Return the current value of the state-dependent normal force 
    magnitude that limits this friction element. **/
    virtual Real getNormalForceMagnitude(const State& state) const = 0;

    virtual Real calcEffectiveCOF(const State& state,
                                  Real defaultTransitionSpeed,
                                  Real slipSpeed) const = 0;

    virtual Real getSlipSpeed(const State& state) const = 0; 

    /** TODO: kludge needed because we're misusing existing constraints. 
    This must be called while Stage::Position is valid. **/
    virtual Vec3 getPositionInfo(const State& state) const 
    {   return Vec3(NaN); }
    /** TODO: kludge to set instance parameters on internal constraints;
    this should be the same Vec3 you got from getPositionInfo(). **/
    virtual void setInstanceParameter(State& state, const Vec3& pos) const {}

    void setMyIndex(StateLimitedFrictionIndex fx) {m_myIx = fx;}
    StateLimitedFrictionIndex getMyIndex() const {return m_myIx;}
private:
    StateLimitedFrictionIndex   m_myIx;
};


//==============================================================================
//                               HARD STOP
//==============================================================================
/** Limit on maximum and/or minimum value of a generalized coordinate. A
generalized force opposes further excursion of the coordinate, and a generalized
impulse is generated when the stop is hit with a non-zero velocity (an impact).
The specified material provides a coefficient of restitution (COR) that 
determines the rebound impulse that occurs as a result of impact. Typically the
COR is velocity-dependent. 

The same material is used for both the upper and
lower stops; if you want different materials use two HardStop objects each
restricting only one direction -- the other limit should be Infinity or 
-Infinity. **/
class SimTK_SIMBODY_EXPORT HardStop : public UnilateralContact {
public:
    HardStop(MobilizedBody& mobod, MobilizerQIndex which,
             Real defaultLowerLimit, Real defaultUpperLimit,
             Real minCOR);

private:
};

//==============================================================================
//                          POINT PLANE CONTACT
//==============================================================================
class SimTK_SIMBODY_EXPORT PointPlaneContact : public UnilateralContact {
public:
    PointPlaneContact(
        MobilizedBody& planeBodyB, const UnitVec3& normal_B, Real height,
        MobilizedBody& followerBodyF, const Vec3& point_F, 
        Real minCOR, Real mu_s, Real mu_d, Real mu_v);

    bool disable(State& state) const OVERRIDE_11 {
        if (m_ptInPlane.isDisabled(state)) return false;
        m_ptInPlane.disable(state);
        if (m_hasFriction) {m_noslipX.disable(state);m_noslipY.disable(state);}
        return true;
    }

    bool enable(State& state) const OVERRIDE_11 {
        if (!m_ptInPlane.isDisabled(state)) return false;
        m_ptInPlane.enable(state);
        if (m_hasFriction) {m_noslipX.enable(state);m_noslipY.enable(state);}
        return true;
    }

    // Currently have to fake the perr because the constraint might be
    // disabled in which case it won't calculate perr.
    Real getPerr(const State& state) const OVERRIDE_11;

    // We won't need to look at these except for proximal constraints which
    // will already have been enabled, so no need to fake.
    Real getVerr(const State& state) const OVERRIDE_11
    {   return m_ptInPlane.getVelocityError(state); }
    Real getAerr(const State& state) const OVERRIDE_11
    {   return m_ptInPlane.getAccelerationError(state); }


    Real calcEffectiveCOR(const State& state,
                          Real defaultCaptureSpeed,
                          Real defaultMinCORSpeed,
                          Real impactSpeed) const OVERRIDE_11 
    {
       return ConditionalConstraint::calcEffectiveCOR
               (m_minCOR, defaultCaptureSpeed, defaultMinCORSpeed,
                impactSpeed);
    }

    bool hasFriction(const State& state) const OVERRIDE_11
    {   return m_hasFriction; }

    Vec2 getSlipVelocity(const State& state) const  OVERRIDE_11 {
        return Vec2(m_noslipX.getVelocityError(state),
                    m_noslipY.getVelocityError(state));
    }

    Real calcEffectiveCOF(const State& state,
                          Real defaultTransitionSpeed,
                          Real slipSpeed) const OVERRIDE_11
    {
       return ConditionalConstraint::calcEffectiveCOF
               (m_mu_s, m_mu_d, m_mu_v, defaultTransitionSpeed, slipSpeed);
    }

    MultiplierIndex getContactMultiplierIndex(const State& s) const OVERRIDE_11;

    void getFrictionMultiplierIndices(const State&     s, 
                                      MultiplierIndex& ix_x, 
                                      MultiplierIndex& ix_y) const OVERRIDE_11;

    // Return the contact point in the plane body so we can make friction 
    // act there.
    Vec3 getPositionInfo(const State& state) const OVERRIDE_11;
    // Set the friction contact point to the contact point in the plane.
    void setInstanceParameter(State& state, const Vec3& pos) const OVERRIDE_11;

private:
    MobilizedBody&              m_planeBody;    // body P
    const Rotation              m_frame;        // z is normal; expressed in P
    const Real                  m_height;

    MobilizedBody&              m_follower;     // body F
    const Vec3                  m_point;        // measured & expressed in F

    Real                        m_minCOR;
    Real                        m_mu_s, m_mu_d, m_mu_v;

    bool                        m_hasFriction;

    Constraint::PointInPlane    m_ptInPlane;
    Constraint::NoSlip1D        m_noslipX, m_noslipY; // stiction
};


} // namespace SimTK

#endif // SimTK_SIMBODY_CONDITIONAL_CONSTRAINT_H_
