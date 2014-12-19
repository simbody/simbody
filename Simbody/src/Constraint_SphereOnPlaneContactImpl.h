#ifndef SimTK_SIMBODY_CONSTRAINT_SPHERE_ON_PLANE_CONTACT_IMPL_H_
#define SimTK_SIMBODY_CONSTRAINT_SPHERE_ON_PLANE_CONTACT_IMPL_H_

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

/**@file
Private implementation of Constraint::SphereOnPlaneContact. **/

#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"

#include "ConstraintImpl.h"

class SimbodyMatterSubsystemRep;

namespace SimTK {

class SimbodyMatterSubsystem;
class SimbodyMatterSubtree;
class MobilizedBody;

//==============================================================================
//                      SPHERE ON PLANE CONTACT IMPL
//==============================================================================
class Constraint::SphereOnPlaneContactImpl : public ConstraintImpl {
public:

// An object of this local type is used as the value of a Position-stage
// discrete state variable that allows plane and sphere parameters to be
// changed for a particular State.
struct Parameters {
    Parameters(const Transform& X_FP, const Vec3& p_BO, Real radius)
    :   m_X_FP(X_FP), m_p_BO(p_BO), m_radius(radius) {}
    Transform   m_X_FP;     // plane frame
    Vec3        m_p_BO;     // sphere center
    Real        m_radius;   // sphere radius
};

explicit SphereOnPlaneContactImpl(bool enforceRolling)
:   ConstraintImpl(1, enforceRolling?2:0, 0), 
    m_enforceRolling(enforceRolling), 
    m_def_X_FP(), m_def_p_BO(0), m_def_radius(NaN),
    m_planeHalfWidth(1)
{ }
SphereOnPlaneContactImpl* clone() const override
{   return new SphereOnPlaneContactImpl(*this); }

void calcDecorativeGeometryAndAppendVirtual
    (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
    override;

// Allocates the discrete state variable for the parameters.
void realizeTopologyVirtual(State& state) const override;

// Get the current value of the runtime-settable parameters from this state.
const Parameters& getParameters(const State& state) const;

// Get a writable reference to the value of the discrete variable; this 
// invalidated Position stage in the given state.
Parameters& updParameters(State& state) const;

void setPlaneDisplayHalfWidth(Real h) {
    // h <= 0 means don't display plane
    m_planeHalfWidth = h > 0 ? h : 0;
}
Real getPlaneDisplayHalfWidth() const {return m_planeHalfWidth;}

// Implementation of virtuals required for holonomic constraints.

/* Body B, the "ball" body has a sphere fixed to it with center O and radius r. 
Call the contact point at the bottom (w.r.t. plane normal) of the sphere C, 
with C=O-r*Pz. We'll define contact to occur at the contact point C.
The body to which the plane P is attached is called the "floor" body F.

The position constraint equation for contact enforces that point C must always
touch the plane P, that is, pz_PC=p_PC[2]=0, or equivalent

    (1) perr = pz_PC = p_PC.Pz = (p_PO-r*Pz).Pz = (p_FO-p_FP).Pz-r 
             = p_FO.Pz-(h+r)
             = (p_AO-p_AF).Pz-(h+r)
    (2) verr = d/dt_F perr = (d/dt_F (p_AO-p_AF)) . Pz 
                              + (p_AO-p_AF) . d/dt_F Pz  <-- this term is 0
             = [d/dt_A (p_AO-p_AF) - w_AF x (p_AO-p_AF)] . Pz
             = [(v_AO-v_AF) - w_AF x (p_AO-p_AF)] . Pz
             = v_FO . Pz
    (3) aerr = (d/dt)_F verr = a_FO . Pz
             = [a_AO-a_AF - b_AF x (p_AO-p_AF) - 2 w_AF x (v_AO-v_AF)
                          + w_AF x (w_AF x (p_AO-p_AF))] . Pz

To derive the final expression we have:
    a_FO = d/dt_F v_FO = d/dt_A v_FO - w_AF x v_FO
where
    d/dt_A v_FO = a_AO-a_AF - [b_AF x (p_AO-p_AF) + w_AF x (v_AO-v_AF)]
    w_AF x v_FO = w_AF x (v_AO-v_AF) - w_AF x (w_AF x (p_AO-p_AF))
Combining these gives the last expression.

Note that we took the time derivative in the F frame. But the scalars returned
by perr,verr, and aerr are measure numbers along Pz. The multiplier will 
consequently act along Pz.

Note that because the ball center and bottom point lie along the plane normal
and the force is applied along that line, it doesn't actually matter whether
we apply the force at the bottom point or at the center; the force and torque
applied at the rigid bodies' origins are the same either way so there is no
way to tell.

(See comments below for friction.)

    --------------------------------
    perr = p_FO_A.Pz_A - (h+r)
    --------------------------------
*/
void calcPositionErrorsVirtual      
   (const State&                                    s,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const override
{
    assert(allX_AB.size()==2 && constrainedQ.size()==0 && perr.size() == 1);

    const Parameters& params = getParameters(s);
    const Transform&  X_FP = params.m_X_FP;
    const Vec3&       p_BO = params.m_p_BO;
    const Real        r    = params.m_radius;

    const Vec3&     Po_F = X_FP.p();  // Plane origin in F
    const UnitVec3& Pz_F = X_FP.z();  // Plane normal direction in F
    const Real h = dot(Po_F, Pz_F);     // Height; could precalculate (5 flops)

    const Transform& X_AF = getBodyTransform(allX_AB, m_planeBody_F);
    const Vec3 p_AO =  findStationLocation(allX_AB, m_ballBody_B, p_BO); 
                                                                    // 18 flops
    const Vec3     p_FO_A = p_AO - X_AF.p();    //  3 flops
    const UnitVec3 Pz_A = X_AF.R() * Pz_F;      // 15 flops

    // Calculate this scalar using A-frame vectors; any frame would have done.
    perr[0] = ~p_FO_A * Pz_A - (h + r);  //  7 flops
}


/*  
    --------------------------------
    verr = v_FO_A.Pz_A
    -------------------------------- 
*/
void calcPositionDotErrorsVirtual      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const override
{
    assert(allV_AB.size()==2 && constrainedQDot.size()==0 && pverr.size() == 1);
    
    const Parameters& params = getParameters(s);
    const Transform&  X_FP = params.m_X_FP;
    const Vec3&       p_BO = params.m_p_BO;

    const Transform&  X_AF = getBodyTransformFromState(s, m_planeBody_F);
    const SpatialVec& V_AF = allV_AB[m_planeBody_F];
    const Vec3&       p_AF = X_AF.p();
    const Vec3&       w_AF = V_AF[0];
    const Vec3&       v_AF = V_AF[1];
    const UnitVec3    Pz_A = X_AF.R() * X_FP.z();                // 15 flops

    const Transform&  X_AB = getBodyTransformFromState(s, m_ballBody_B);
    const SpatialVec& V_AB = allV_AB[m_ballBody_B];
    const Vec3        p_BO_A = X_AB.R() * p_BO; // re-express in A (15 flops)

    Vec3 p_AO, v_AO;
    findStationInALocationVelocity(X_AB, V_AB, p_BO_A, p_AO, v_AO); // 15 flops

    const Vec3 p_FO_A     = p_AO - p_AF; // measure O from Fo, in A (3 flops)
    const Vec3 p_FO_A_dot = v_AO - v_AF; // derivative in A (3 flops)
    const Vec3 v_FO_A = p_FO_A_dot - w_AF % p_FO_A; // deriv in F (12 flops)

    // Calculate this scalar using A-frame vectors.
    pverr[0] = ~v_FO_A * Pz_A;                           // 5 flops
}

/*  
    --------------------------------
    aerr = a_FO_A.Pz_A
    -------------------------------- 
*/
void calcPositionDotDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const override
{
    assert(allA_AB.size()==2 && constrainedQDotDot.size()==0 && paerr.size()==1);
    
    const Parameters& params = getParameters(s);
    const Transform&  X_FP = params.m_X_FP;
    const Vec3&       p_BO = params.m_p_BO;

    const Transform&  X_AF = getBodyTransformFromState(s, m_planeBody_F);
    const SpatialVec& V_AF = getBodyVelocityFromState(s, m_planeBody_F);
    const SpatialVec& A_AF = allA_AB[m_planeBody_F];
    const Vec3&       p_AF = X_AF.p();
    const Vec3&       w_AF = V_AF[0];
    const Vec3&       v_AF = V_AF[1];
    const Vec3&       b_AF = A_AF[0];
    const Vec3&       a_AF = A_AF[1];
    const UnitVec3    Pz_A = X_AF.R() * X_FP.z();                // 15 flops

    const Transform&  X_AB = getBodyTransformFromState(s, m_ballBody_B);
    const SpatialVec& V_AB = getBodyVelocityFromState(s, m_ballBody_B);
    const SpatialVec& A_AB = allA_AB[m_ballBody_B];
    const Vec3        p_BO_A = X_AB.R() * p_BO; // re-express in A (15 flops)

    Vec3 p_AO, v_AO, a_AO;
    findStationInALocationVelocityAcceleration(X_AB, V_AB, A_AB, p_BO_A,
                                               p_AO, v_AO, a_AO); // 39 flops

    const Vec3 p_FO_A        = p_AO - p_AF; // measure O from Fo, in A (3 flops)
    const Vec3 p_FO_A_dot    = v_AO - v_AF; // deriv in A (3 flops)
    const Vec3 p_FO_A_dotdot = a_AO - a_AF; // 2nd deriv in A (3 flops)

    const Vec3 v_FO_A     = p_FO_A_dot    - w_AF % p_FO_A; // in F (12 flops)
    const Vec3 v_FO_A_dot = p_FO_A_dotdot -(b_AF % p_FO_A + w_AF % p_FO_A_dot);
                                                      // deriv in A (24 flops)
    const Vec3 a_FO_A = v_FO_A_dot - w_AF % v_FO_A; // deriv in F (12 flops)

    // Calculate this scalar using A-frame vectors.
    paerr[0] = ~a_FO_A * Pz_A;  // 5 flops
}

// apply f=lambda*Pz to the bottom point C of ball on body B,
//       -f         to point C (coincident point) of body F
void addInPositionConstraintForcesVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const override
{
    assert(multipliers.size()==1 && bodyForcesInA.size()==2 
           && qForces.size()==0);
    const Real lambda = multipliers[0];

    const Parameters& params = getParameters(s);
    const Transform&  X_FP = params.m_X_FP;
    const Vec3&       p_BO = params.m_p_BO;
    const Real        r    = params.m_radius;

    const Transform& X_AF = getBodyTransformFromState(s, m_planeBody_F);
    const Vec3&      p_AF = X_AF.p(); // p_AoFo_A
    const UnitVec3   Pz_A = X_AF.R() * X_FP.z();                   // 15 flops

    const Transform& X_AB = getBodyTransformFromState(s, m_ballBody_B);
    const Vec3&      p_AB = X_AB.p(); // p_AoBo_A

    const Vec3 p_BO_A = X_AB.R() * p_BO;   // re-express in A        (15 flops)
    const Vec3 p_BC_A = p_BO_A - r*Pz_A;   // bottom of ball         ( 6 flops)
    const Vec3 p_AC   = p_AB + p_BC_A;     // measure from Ao        ( 3 flops)
    const Vec3 p_FC_A = p_AC - p_AF;       // measure C from Fo, in A (3 flops)

    const Vec3 force_A = lambda * Pz_A;     // 3 flops
    // 30 flops for the next two calls
    addInStationInAForce(p_BC_A, force_A, bodyForcesInA[m_ballBody_B]);
    subInStationInAForce(p_FC_A, force_A, bodyForcesInA[m_planeBody_F]);
}


// Implementation of virtuals required for nonholonomic constraints.

/* The rolling friction constraint says that the velocity of the material 
point of B at the contact point C = O - r Pz 
(at the bottom of the ball with center O and radius r), measured 
with respect to the plane body F, must be zero in the plane's Px and Py 
directions. Note that the contact point is not a station of body B but moves 
in the B frame. We have
    p_BC = p_BO - r Pz
    p_FC = p_FO - r Pz
    verr = v_FC = d/dt_F p_FC 
         = v_FO - r w_FB x Pz
where   v_FO = (v_AO-v_AF) - w_AF x (p_AO-p_AF)
and     w_FB = (w_AB-w_AF)

You have to differentiate verr carefully. Because Pz is fixed in F, the 
result is not the acceleration of the material point at C, but rather of
the moving contact point C.

    aerr = d/dt_F verr = a_FO - r d/dt_F (w_FB x Pz)
    d/dt_F (w_FB x Pz) = (d/dt_F w_FB) x Pz + w_FB x (d/dt_F Pz)
                       = b_FB x Pz  [2nd term dropped out!]
so  aerr = a_FC = a_FO - r b_FB x Pz

where a_FO = d/dt_F v_FO
and   b_FB = d/dt_F w_FB = (b_AB-b_AF) - w_AF x (w_AB-w_AF)
           = (b_AB-b_AF) - w_AF x w_AB

d/dt_F v_FO = d/dt_A v_FO - w_AF x v_FO
d/dt_A v_FO = (a_AO-a_AF) - b_AF x (p_AO-p_AF) - w_AF x (v_AO-v_AF)
w_AF x v_FO = w_AF x (v_AO-v_AF) - w_AF x (w_AF x (p_AO-p_AF))

so a_FO = (a_AO-a_AF) - b_AF x (p_AO-p_AF) - 2 w_AF x (v_AO-v_AF)
          + w_AF x (w_AF x (p_AO-p_AF))
and
aerr = a_FC = a_FO - r [b_AB-b_AF - w_AF x w_AB] x Pz

    -------------------------------------------
    verr = [v_FC_A . Px_A]
           [v_FC_A . Py_A]
    -------------------------------------------
    ~ 110 flops TODO should precalculate where possible
*/
void calcVelocityErrorsVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr)   // mv of these
    const override
{
    assert(allV_AB.size()==2 && constrainedU.size()==0 && verr.size()==2);

    const Parameters& params = getParameters(s);
    const Transform&  X_FP = params.m_X_FP;
    const Vec3&       p_BO = params.m_p_BO;
    const Real        r    = params.m_radius;

    const Transform&  X_AF = getBodyTransformFromState(s, m_planeBody_F);
    const SpatialVec& V_AF = allV_AB[m_planeBody_F];
    const Vec3&       p_AF = X_AF.p();
    const Vec3&       w_AF = V_AF[0];
    const Vec3&       v_AF = V_AF[1];
    const Rotation    R_AP = X_AF.R() * X_FP.R(); //R_AP=[Px Py Pz]_A (45 flops)

    const Transform&  X_AB = getBodyTransformFromState(s, m_ballBody_B);
    const SpatialVec& V_AB = allV_AB[m_ballBody_B];
    const Vec3        p_BO_A = X_AB.R() * p_BO; // re-express in A (15 flops)
    const Vec3        p_BC_A = p_BO_A - r*R_AP.z(); // bottom of ball (6 flops)

    Vec3 p_AC, v_AC;
    findStationInALocationVelocity(X_AB, V_AB, p_BC_A, p_AC, v_AC); // 15 flops

    // Calculate relative velocity of C in F, expressed in A.
    const Vec3 p_FC_A     = p_AC - p_AF; // measure C from Fo, in A (3 flops)
    const Vec3 p_FC_A_dot = v_AC - v_AF; // derivative in A (3 flops)
    const Vec3 v_FC_A = p_FC_A_dot - w_AF % p_FC_A; // deriv in F (12 flops)

    // Calculate these scalars using A-frame vectors, but the results are
    // measure numbers in [Px Py].
    verr[0] = ~R_AP.x()*v_FC_A; // 5 flops
    verr[1] = ~R_AP.y()*v_FC_A; // 5 flops
}

/*
    -------------------------------------------
    a_FC_A = a_FO_A - r b_FB_A x Pz_A
    aerr = [a_FC_A . Px_A]
           [a_FC_A . Py_A]
    -------------------------------------------
    ~ 200 flops TODO should precalculate where possible
*/

void calcVelocityDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr)  // mv of these
    const override
{
    assert(allA_AB.size()==2 && constrainedUDot.size()==0 && vaerr.size()==2);

    const Parameters& params = getParameters(s);
    const Transform&  X_FP = params.m_X_FP;
    const Vec3&       p_BO = params.m_p_BO;
    const Real        r    = params.m_radius;

    const Transform&  X_AF = getBodyTransformFromState(s, m_planeBody_F);
    const SpatialVec& V_AF = getBodyVelocityFromState(s, m_planeBody_F);
    const SpatialVec& A_AF = allA_AB[m_planeBody_F];
    const Vec3&       p_AF = X_AF.p();
    const Vec3&       w_AF = V_AF[0];
    const Vec3&       v_AF = V_AF[1];
    const Vec3&       b_AF = A_AF[0];
    const Vec3&       a_AF = A_AF[1];
    const Rotation    R_AP = X_AF.R() * X_FP.R(); //R_AP=[Px Py Pz]_A (45 flops)

    const Transform&  X_AB = getBodyTransformFromState(s, m_ballBody_B);
    const SpatialVec& V_AB = getBodyVelocityFromState(s, m_ballBody_B);
    const SpatialVec& A_AB = allA_AB[m_ballBody_B];
    const Vec3&       w_AB = V_AB[0];
    const Vec3&       b_AB = A_AB[0];
    const Vec3        p_BO_A = X_AB.R() * p_BO; // re-express in A (15 flops)

    Vec3 p_AO, v_AO, a_AO;
    findStationInALocationVelocityAcceleration(X_AB, V_AB, A_AB, p_BO_A,
                                               p_AO, v_AO, a_AO); // 39 flops

    const Vec3 p_FO_A        = p_AO - p_AF; // measure O from Fo, in A (3 flops)
    const Vec3 p_FO_A_dot    = v_AO - v_AF; // deriv in A (3 flops)
    const Vec3 p_FO_A_dotdot = a_AO - a_AF; // 2nd deriv in A (3 flops)

    const Vec3 v_FO_A     = p_FO_A_dot    - w_AF % p_FO_A; // in F (12 flops)
    const Vec3 v_FO_A_dot = p_FO_A_dotdot -(b_AF % p_FO_A + w_AF % p_FO_A_dot);
                                                      // deriv in A (24 flops)
    const Vec3 a_FO_A = v_FO_A_dot - w_AF % v_FO_A; // now in F (12 flops)

    const Vec3 w_FB_A_dot = b_AB - b_AF;                        // 3 flops
    const Vec3 b_FB_A = w_FB_A_dot - w_AF % w_AB;               // 12 flops
    const Vec3 a_FC_A = a_FO_A - b_FB_A % (r*R_AP.z());         // 15 flops

    // Calculate these scalars using A-frame vectors, but the results are
    // measure numbers in [Px Py].
    vaerr[0] = ~R_AP.x()*a_FC_A; // 5 flops
    vaerr[1] = ~R_AP.y()*a_FC_A; // 5 flops
}

// apply f=lambda0*Px + lambda1*Py to the bottom point C of ball on B
//      -f           to point C (coincident point) of body F
void addInVelocityConstraintForcesVirtual
   (const State&                                    s,      // Stage::Velocity
    const Array_<Real>&                             multipliers, // mv of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedUIndex>&           mobilityForces) 
    const override
{
    assert(multipliers.size()==2 && mobilityForces.size()==0 
           && bodyForcesInA.size()==2);
    const Real lambda0 = multipliers[0], lambda1 = multipliers[1];

    const Parameters& params = getParameters(s);
    const Transform&  X_FP = params.m_X_FP;
    const Vec3&       p_BO = params.m_p_BO;
    const Real        r    = params.m_radius;

    const Transform& X_AF = getBodyTransformFromState(s, m_planeBody_F);
    const Vec3&      p_AF = X_AF.p(); // p_AoFo_A
    const Rotation   R_AP = X_AF.R() * X_FP.R(); //R_AP=[Px Py Pz]_A (45 flops)

    const Transform& X_AB = getBodyTransformFromState(s, m_ballBody_B);
    const Vec3&      p_AB = X_AB.p(); // p_AoBo_A

    const Vec3 p_BO_A = X_AB.R() * p_BO;   // re-express in A        (15 flops)
    const Vec3 p_BC_A = p_BO_A - r*R_AP.z(); // bottom of ball       ( 6 flops)
    const Vec3 p_AC   = p_AB + p_BC_A;     // measure from Ao        ( 3 flops)
    const Vec3 p_FC_A = p_AC - p_AF;       // measure C from Fo, in A (3 flops)

    const Vec3 force_A = lambda0*R_AP.x() + lambda1*R_AP.y(); // 9 flops

    // 30 flops for the next two calls
    addInStationInAForce(p_BC_A, force_A, bodyForcesInA[m_ballBody_B]);
    subInStationInAForce(p_FC_A, force_A, bodyForcesInA[m_planeBody_F]);
}

SimTK_DOWNCAST(SphereOnPlaneContactImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::SphereOnPlaneContact;

// These can't be changed after construction.
ConstrainedBodyIndex    m_planeBody_F;      // F (B1)
ConstrainedBodyIndex    m_ballBody_B;       // B (B2)
const bool              m_enforceRolling;   // if so, add 2 constraints

Transform               m_def_X_FP;     // default plane frame fixed in F
Vec3                    m_def_p_BO;     // default sphere ctr point fixed in B
Real                    m_def_radius;   // default sphere radius

// This is just for visualization
Real                    m_planeHalfWidth;

// This Position-stage variable holds the constraint parameters to be used
// for a particular State.
mutable DiscreteVariableIndex   parametersIx;
};



} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_SPHERE_ON_PLANE_CONTACT_IMPL_H_



