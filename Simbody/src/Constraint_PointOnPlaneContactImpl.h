#ifndef SimTK_SIMBODY_CONSTRAINT_POINT_ON_PLANE_CONTACT_IMPL_H_
#define SimTK_SIMBODY_CONSTRAINT_POINT_ON_PLANE_CONTACT_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.        *
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
Private implementation of Constraint::PointOnPlaneContact. **/

#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"

#include "ConstraintImpl.h"


namespace SimTK {

class SimbodyMatterSubsystem;
class SimbodyMatterSubsystemRep;
class SimbodyMatterSubtree;
class MobilizedBody;

//==============================================================================
//                      POINT ON PLANE CONTACT IMPL
//==============================================================================
class Constraint::PointOnPlaneContactImpl : public ConstraintImpl {
public:
PointOnPlaneContactImpl()
:   ConstraintImpl(1,2,0), m_X_SP(), m_p_BF(0), 
    m_planeHalfWidth(1), m_pointRadius(Real(0.05)) 
{ }
PointOnPlaneContactImpl* clone() const override 
{   return new PointOnPlaneContactImpl(*this); }

void calcDecorativeGeometryAndAppendVirtual
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
    override;

void setPlaneDisplayHalfWidth(Real h) {
    // h <= 0 means don't display plane
    invalidateTopologyCache();
    m_planeHalfWidth = h > 0 ? h : 0;
}
Real getPlaneDisplayHalfWidth() const {return m_planeHalfWidth;}

void setPointDisplayRadius(Real r) {
    // r <= 0 means don't display point
    invalidateTopologyCache();
    m_pointRadius= r > 0 ? r : 0;
}
Real getPointDisplayRadius() const {return m_pointRadius;}

// Implementation of virtuals required for the holonomic constraint.
// Here C is the station of S that is instantaneously coincident with F on B.

//    --------------------------------
//    perr = ~p_SF*Pz - h
//    --------------------------------
void calcPositionErrorsVirtual      
   (const State&                                    s,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const override
{
    assert(allX_AB.size()==2 && constrainedQ.size()==0 && perr.size() == 1);

    const Transform& X_AS = getBodyTransform(allX_AB, m_surfaceBody_S);
    const Vec3 p_AF = 
        findStationLocation(allX_AB, m_followerBody_B, m_p_BF); // 18 flops
    const Vec3 p_SC = ~X_AS * p_AF; // shift to So, reexpress in S (18 flops)

    // Calculate this scalar using S-frame vectors; any frame would have done.
    const Vec3& Po = m_X_SP.p(); // Plane origin in S
    perr[0] = dot(p_SC-Po, m_X_SP.z()); // 8 flops
}

//    -----------------------------------
//    verr = ~v_SF*Pz = ~(v_AF-v_AC)*Pz_A
//    -----------------------------------
void calcPositionDotErrorsVirtual      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const override 
{
    assert(allV_AB.size()==2 && constrainedQDot.size()==0 && pverr.size() == 1);

    const Transform&  X_AS = getBodyTransformFromState(s, m_surfaceBody_S);
    const SpatialVec& V_AS = allV_AB[m_surfaceBody_S];

    const Transform&  X_AB = getBodyTransformFromState(s, m_followerBody_B);
    const SpatialVec& V_AB = allV_AB[m_followerBody_B];

    const Vec3 p_BF_A = X_AB.R() * m_p_BF; // re-express in A (15 flops)
    Vec3 p_AF, v_AF;
    findStationInALocationVelocity(X_AB, V_AB, p_BF_A, p_AF, v_AF); // 15 flops

    const Vec3 p_SC_A = p_AF - X_AS.p(); // measure C from So, in A (3 flops)
    Vec3 p_AC, v_AC;
    findStationInALocationVelocity(X_AS, V_AS, p_SC_A, p_AC, v_AC); // 15 flops

    const UnitVec3 Pz_A  = X_AS.R() * m_X_SP.z();                // 15 flops

    // Calculate this scalar using A-frame vectors.
    pverr[0] = dot( v_AF-v_AC, Pz_A );                           // 8 flops
}

//    --------------------------------------------------------------
//    aerr = ~a_SF*Pz = ~((a_AF-a_AC) - 2 w_AS X (v_AF-v_AC)) * Pz_A
//    --------------------------------------------------------------
void calcPositionDotDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const override
{
    assert(allA_AB.size()==2 && constrainedQDotDot.size()==0 && paerr.size()==1);

    const Transform&  X_AS = getBodyTransformFromState(s, m_surfaceBody_S);
    const SpatialVec& V_AS = getBodyVelocityFromState(s, m_surfaceBody_S);
    const SpatialVec& A_AS = allA_AB[m_surfaceBody_S];

    const Transform&  X_AB = getBodyTransformFromState(s, m_followerBody_B);
    const SpatialVec& V_AB = getBodyVelocityFromState(s, m_followerBody_B);
    const SpatialVec& A_AB = allA_AB[m_followerBody_B];

    const Vec3 p_BF_A = X_AB.R() * m_p_BF; // re-express in A (15 flops)
    Vec3 p_AF, v_AF, a_AF;
    findStationInALocationVelocityAcceleration(X_AB, V_AB, A_AB, p_BF_A,
                                               p_AF, v_AF, a_AF); // 39 flops

    const Vec3 p_SC_A = p_AF - X_AS.p(); // measure C from So, in A  (3 flops)
    Vec3 p_AC, v_AC, a_AC;
    findStationInALocationVelocityAcceleration(X_AS, V_AS, A_AS, p_SC_A,
                                               p_AC, v_AC, a_AC); // 39 flops

    const UnitVec3 Pz_A  = X_AS.R() * m_X_SP.z();                 // 15 flops

    const Vec3& w_AS = V_AS[0];
    const Vec3 v_CF_A = v_AF-v_AC;          // 3 flops
    const Vec3 wXv2 = (2.*w_AS) % v_CF_A;   // 12 flops
    const Vec3 a_CF_A = a_AF-a_AC;          // 3 flops
    const Vec3 aRel_A = a_CF_A - wXv2; // relative accel of F and C  (3 flops)

    // Calculate this scalars using A-frame vectors, but the result is the
    // measure number along Pz.
    paerr[0] = ~Pz_A*aRel_A;                // 5 flops
}

// apply f=lambda*Pz to the follower point F of body B,
//      -f           to point C (coincident point) of body S
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

    const Transform& X_AS = getBodyTransformFromState(s, m_surfaceBody_S);
    const Transform& X_AB = getBodyTransformFromState(s, m_followerBody_B);

    const Vec3 p_BF_A = X_AB.R() * m_p_BF; // re-express in A (15 flops)
    const Vec3 p_AF   = X_AB.p() + p_BF_A; // shift to Ao     ( 3 flops)
    const Vec3 p_SC_A = p_AF - X_AS.p();   // measure C from So, in A (3 flops)

    const UnitVec3 Pz_A  = X_AS.R() * m_X_SP.z();                  // 15 flops

    const Vec3 force_A = lambda * Pz_A;

    addInStationInAForce(p_BF_A, force_A, bodyForcesInA[m_followerBody_B]);
    subInStationInAForce(p_SC_A, force_A, bodyForcesInA[m_surfaceBody_S]);
}

// Implementation of virtuals required for nonholonomic constraints.

//    -----------------------------------------
//    verr = [~Px] * v_SF = [~Px_A] * (v_AF-v_AC)
//           [~Py]          [~Py_A]
//    -----------------------------------------
// TODO: rework this using the more efficient API (see next method).
void calcVelocityErrorsVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr)   // mv of these
    const override
{
    assert(allV_AB.size()==2 && constrainedU.size()==0 && verr.size()==2);

    const Transform& X_AS = getBodyTransformFromState(s, m_surfaceBody_S);
    const Vec3 p_AF = 
        findStationLocationFromState(s, m_followerBody_B, m_p_BF); // 18 flops
    const Vec3 p_SC = ~X_AS * p_AF; // shift to So, reexpress in S (18 flops)
    const UnitVec3 Px_A  = X_AS.R() * m_X_SP.x();                  // 15 flops
    const UnitVec3 Py_A  = X_AS.R() * m_X_SP.y();                  // 15 flops
    
    // 54 flops for the two of these
    const Vec3 v_AF = findStationVelocity(s, allV_AB, m_followerBody_B, m_p_BF);
    const Vec3 v_AC = findStationVelocity(s, allV_AB, m_surfaceBody_S, p_SC);
    const Vec3 v_CF_A = v_AF-v_AC; // 3 flops

    // Calculate these scalars using A-frame vectors, but the results are
    // measure numbers in [Px Py].
    verr[0] = ~Px_A*v_CF_A; // 5 flops
    verr[1] = ~Py_A*v_CF_A; // 5 flops
}

//    ------------------------------------------------------------------
//    aerr = [~Px] * a_SF = [~Px_A] * ((a_AF-a_AC) - 2 w_AS X (v_AF-v_AC))
//           [~Py]          [~Py_A]
//    ------------------------------------------------------------------

void calcVelocityDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr)  // mv of these
    const override
{
    assert(allA_AB.size()==2 && constrainedUDot.size()==0 && vaerr.size()==2);

    const Transform&  X_AS = getBodyTransformFromState(s, m_surfaceBody_S);
    const SpatialVec& V_AS = getBodyVelocityFromState(s, m_surfaceBody_S);
    const SpatialVec& A_AS = allA_AB[m_surfaceBody_S];

    const Transform&  X_AB = getBodyTransformFromState(s, m_followerBody_B);
    const SpatialVec& V_AB = getBodyVelocityFromState(s, m_followerBody_B);
    const SpatialVec& A_AB = allA_AB[m_followerBody_B];

    const Vec3 p_BF_A = X_AB.R() * m_p_BF; // re-express in A (15 flops)
    Vec3 p_AF, v_AF, a_AF;
    findStationInALocationVelocityAcceleration(X_AB, V_AB, A_AB, p_BF_A,
                                               p_AF, v_AF, a_AF); // 39 flops

    const Vec3 p_SC_A = p_AF - X_AS.p(); // measure C from So, in A (3 flops)
    Vec3 p_AC, v_AC, a_AC;
    findStationInALocationVelocityAcceleration(X_AS, V_AS, A_AS, p_SC_A,
                                               p_AC, v_AC, a_AC); // 39 flops

    const UnitVec3 Px_A  = X_AS.R() * m_X_SP.x();                  // 15 flops
    const UnitVec3 Py_A  = X_AS.R() * m_X_SP.y();                  // 15 flops

    const Vec3& w_AS = V_AS[0];
    const Vec3 v_CF_A = v_AF-v_AC;          // 3 flops
    const Vec3 wXv2 = (2.*w_AS) % v_CF_A;   // 12 flops

    const Vec3 a_CF_A = a_AF-a_AC; // 3 flops
    const Vec3 aRel_A = a_CF_A - wXv2; // relative accel of F and C (3 flops)

    // Calculate these scalars using A-frame vectors, but the results are
    // measure numbers in [Px Py].
    vaerr[0] = ~Px_A*aRel_A; // 5 flops
    vaerr[1] = ~Py_A*aRel_A; // 5 flops
}

// apply f=lambda0*Px + lambda1*Py to the follower point F on B
//      -f           to point C (coincident point) of body S
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

    const Transform& X_AS = getBodyTransformFromState(s, m_surfaceBody_S);
    const Transform& X_AB = getBodyTransformFromState(s, m_followerBody_B);

    const Vec3 p_BF_A = X_AB.R() * m_p_BF; // re-express in A (15 flops)
    const Vec3 p_AF   = X_AB.p() + p_BF_A; // shift to Ao     ( 3 flops)
    const Vec3 p_SC_A = p_AF - X_AS.p();   // measure C from So, in A (3 flops)

    const UnitVec3 Px_A  = X_AS.R() * m_X_SP.x();                  // 15 flops
    const UnitVec3 Py_A  = X_AS.R() * m_X_SP.y();                  // 15 flops

    const Vec3 force_A = lambda0*Px_A + lambda1*Py_A; // 9 flops

    addInStationInAForce(p_BF_A, force_A, bodyForcesInA[m_followerBody_B]);
    subInStationInAForce(p_SC_A, force_A, bodyForcesInA[m_surfaceBody_S]);
}

SimTK_DOWNCAST(PointOnPlaneContactImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::PointOnPlaneContact;

ConstrainedBodyIndex    m_surfaceBody_S;    // S (B1)
ConstrainedBodyIndex    m_followerBody_B;   // B (B2)

Transform               m_X_SP; // default plane frame fixed in S
Vec3                    m_p_BF; // default follower point fixed in B

// These are just for visualization
Real                    m_planeHalfWidth;
Real                    m_pointRadius;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_POINT_ON_PLANE_CONTACT_IMPL_H_



