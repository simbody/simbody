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

#include "simbody/internal/common.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/Constraint.h"
#include "simbody/internal/ConditionalConstraint.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"


namespace SimTK {

//==============================================================================
//                           HARD STOP UPPER / LOWER
//==============================================================================
HardStopUpper::HardStopUpper(MobilizedBody& mobod, MobilizerQIndex whichQ,
                             Real defaultUpperLimit, Real minCOR)
:   UnilateralContact(-1), m_mobod(mobod), 
    m_defaultUpperLimit(defaultUpperLimit), m_minCOR(minCOR)
{
    SimTK_ERRCHK1_ALWAYS(0<=minCOR && minCOR<=1,
        "HardStopUpper()", "The coefficient of restitution must be between "
        "0 and 1 but was %g.", minCOR);

    // Set up the constraint.
    m_upper = Constraint::ConstantCoordinate(mobod, whichQ, defaultUpperLimit);
    m_upper.setIsConditional(true);
    m_upper.setDisabledByDefault(true);
}

HardStopLower::HardStopLower(MobilizedBody& mobod, MobilizerQIndex whichQ,
                             Real defaultLowerLimit, Real minCOR)
:   UnilateralContact(1), m_mobod(mobod), 
    m_defaultLowerLimit(defaultLowerLimit), m_minCOR(minCOR)
{
    SimTK_ERRCHK1_ALWAYS(0<=minCOR && minCOR<=1,
        "HardStopLower()", "The coefficient of restitution must be between "
        "0 and 1 but was %g.", minCOR);

    // Set up the constraint.
    m_lower = Constraint::ConstantCoordinate(mobod, whichQ, defaultLowerLimit);
    m_lower.setIsConditional(true);
    m_lower.setDisabledByDefault(true);
}

//------------------------------------------------------------------------------
//                            WHERE TO DISPLAY
//------------------------------------------------------------------------------
Vec3 HardStopUpper::whereToDisplay(const State& state) const {
    const Transform& X_BM = m_mobod.getOutboardFrame(state);
    return m_mobod.findStationLocationInGround(state,X_BM.p());
}
Vec3 HardStopLower::whereToDisplay(const State& state) const {
    const Transform& X_BM = m_mobod.getOutboardFrame(state);
    return m_mobod.findStationLocationInGround(state,X_BM.p());
}

//------------------------------------------------------------------------------
//                              GET PERR/VERR/AERR
//------------------------------------------------------------------------------
Real HardStopUpper::getPerr(const State& state) const
{   //return m_upper.getPositionError(state);
    const MobilizerQIndex whichQ = m_upper.getWhichQ();
    const Real p = m_upper.getPosition(state);
    const Real q = m_mobod.getOneQ(state, whichQ);
    return q - p; // positive when violated (q>p)
}
Real HardStopLower::getPerr(const State& state) const
{   //return m_lower.getPositionError(state);
    const MobilizerQIndex whichQ = m_lower.getWhichQ();
    const Real p = m_lower.getPosition(state);
    const Real q = m_mobod.getOneQ(state, whichQ);
    return q - p; // negative when violated (q<p)
}
Real HardStopUpper::getVerr(const State& state) const
{   return m_upper.getVelocityError(state); }

Real HardStopLower::getVerr(const State& state) const
{   return m_lower.getVelocityError(state); }

Real HardStopUpper::getAerr(const State& state) const
{   return m_upper.getAccelerationError(state); }

Real HardStopLower::getAerr(const State& state) const
{   return m_lower.getAccelerationError(state); }

//------------------------------------------------------------------------------
//                    GET CONTACT MULTIPLIER INDEX
//------------------------------------------------------------------------------
MultiplierIndex HardStopUpper::
getContactMultiplierIndex(const State& s) const {
    int mp, mv, ma;
    MultiplierIndex px0, vx0, ax0;
    m_upper.getNumConstraintEquationsInUse(s,mp,mv,ma);
    assert(mp==1 && mv==0 && ma==0); // don't call if not enabled
    m_upper.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    assert(px0.isValid() && !vx0.isValid() && !ax0.isValid());
    return px0;
}

MultiplierIndex HardStopLower::
getContactMultiplierIndex(const State& s) const {
    int mp, mv, ma;
    MultiplierIndex px0, vx0, ax0;
    m_lower.getNumConstraintEquationsInUse(s,mp,mv,ma);
    assert(mp==1 && mv==0 && ma==0); // don't call if not enabled
    m_lower.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    assert(px0.isValid() && !vx0.isValid() && !ax0.isValid());
    return px0;
}


//==============================================================================
//                                  ROPE
//==============================================================================
Rope::Rope(MobilizedBody& mobod1, const Vec3& point1, 
           MobilizedBody& mobod2, const Vec3& point2,
           Real defaultLengthLimit, Real minCOR)
:   UnilateralContact(-1), m_minCOR(minCOR)
{
    SimTK_ERRCHK1_ALWAYS(0<=minCOR && minCOR<=1,
        "Rope()", "The coefficient of restitution must be between "
        "0 and 1 but was %g.", minCOR);

    // Set up the constraint.
    m_rod = Constraint::Rod(mobod1, point1, mobod2, point2, defaultLengthLimit);
    m_rod.setIsConditional(true);
    m_rod.setDisabledByDefault(true);
}

//------------------------------------------------------------------------------
//                            WHERE TO DISPLAY
//------------------------------------------------------------------------------
Vec3 Rope::whereToDisplay(const State& state) const {
    const Vec3 p_GP = m_rod.getMobilizedBody1()
        .findStationLocationInGround(state, m_rod.getPointOnBody1(state));
    const Vec3 p_GQ = m_rod.getMobilizedBody2()
        .findStationLocationInGround(state, m_rod.getPointOnBody2(state));
    return (p_GP+p_GQ)/2;
}

//------------------------------------------------------------------------------
//                              GET PERR/VERR/AERR
//------------------------------------------------------------------------------

Real Rope::getPerr(const State& state) const
{   //return m_lower.getPositionError(state);
    return m_rod.findLengthViolation(state);
}

Real Rope::getVerr(const State& state) const
{   return m_rod.getVelocityError(state); }

Real Rope::getAerr(const State& state) const
{   return m_rod.getAccelerationError(state); }

//------------------------------------------------------------------------------
//                    GET CONTACT MULTIPLIER INDEX
//------------------------------------------------------------------------------
MultiplierIndex Rope::
getContactMultiplierIndex(const State& s) const {
    int mp, mv, ma;
    MultiplierIndex px0, vx0, ax0;
    m_rod.getNumConstraintEquationsInUse(s,mp,mv,ma);
    assert(mp==1 && mv==0 && ma==0); // don't call if not enabled
    m_rod.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    assert(px0.isValid() && !vx0.isValid() && !ax0.isValid());
    return px0;
}


//==============================================================================
//                    POINT PLANE FRICTIONLESS CONTACT
//==============================================================================
PointPlaneFrictionlessContact::PointPlaneFrictionlessContact
   (MobilizedBody& planeBodyB, const UnitVec3& normal_B, Real height,
    MobilizedBody& followerBodyF, const Vec3& point_F, Real minCOR)
:   m_planeBody(planeBodyB), m_frame(normal_B, ZAxis), m_height(height), 
    m_follower(followerBodyF), m_point(point_F), m_minCOR(minCOR)
{
    SimTK_ERRCHK1_ALWAYS(0<=minCOR && minCOR<=1,
        "PointPlaneFrictionlessContact()", 
        "The coefficient of restitution must be between 0 and 1 but was %g.", 
        minCOR);

    // Set up the contact constraint.
    m_ptInPlane = Constraint::PointInPlane
       (planeBodyB, normal_B, height, followerBodyF, point_F);
    m_ptInPlane.setIsConditional(true);
    m_ptInPlane.setDisabledByDefault(true);
}

//------------------------------------------------------------------------------
//                            WHERE TO DISPLAY
//------------------------------------------------------------------------------
Vec3 PointPlaneFrictionlessContact::whereToDisplay(const State& state) const {
    return m_follower.findStationLocationInGround(state,m_point);
}


//------------------------------------------------------------------------------
//                              GET PERR
//------------------------------------------------------------------------------
Real PointPlaneFrictionlessContact::getPerr(const State& state) const
{   //return m_ptInPlane.getPositionError(state); 
    const Vec3 p = m_follower.findStationLocationInAnotherBody
                                (state, m_point, m_planeBody);
    return ~p*m_frame.z() - m_height;
}

//------------------------------------------------------------------------------
//                    GET CONTACT MULTIPLIER INDEX
//------------------------------------------------------------------------------
MultiplierIndex PointPlaneFrictionlessContact::
getContactMultiplierIndex(const State& s) const {
    int mp, mv, ma;
    MultiplierIndex px0, vx0, ax0;
    m_ptInPlane.getNumConstraintEquationsInUse(s,mp,mv,ma);
    assert(mp==1 && mv==0 && ma==0); // don't call if not enabled
    m_ptInPlane.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    assert(px0.isValid() && !vx0.isValid() && !ax0.isValid());
    return px0;
}


//==============================================================================
//                           POINT PLANE CONTACT
//==============================================================================
PointPlaneContact::PointPlaneContact
   (MobilizedBody& planeBodyB, const UnitVec3& normal_B, Real height,
    MobilizedBody& followerBodyF, const Vec3& point_F, 
    Real minCOR, Real mu_s, Real mu_d, Real mu_v)
:   m_planeBody(planeBodyB), m_frame(normal_B, ZAxis), m_height(height), 
    m_follower(followerBodyF), m_point(point_F), m_minCOR(minCOR),
    m_mu_s(mu_s), m_mu_d(mu_d), m_mu_v(mu_v)
{
    SimTK_ERRCHK1_ALWAYS(0<=minCOR && minCOR<=1,
        "PointPlaneContact()", 
        "The coefficient of restitution must be between 0 and 1 but was %g.", 
        minCOR);
    SimTK_ERRCHK3_ALWAYS(mu_s >= 0 && mu_d >= 0 && mu_v >= 0,
        "PointPlaneContact()", 
        "All coefficients of friction must be nonnegative; got "
        "mu_s=%g, mu_d=%g, mu_v=%g.", mu_s, mu_d, mu_v);
    SimTK_ERRCHK2_ALWAYS(mu_d <= mu_s,
        "PointPlaneContact()", 
        "The dynamic coefficient of friction can't be larger than "
        "the static coefficient; got mu_s=%g, mu_d=%g.", mu_s, mu_d);

    // Set up the contact constraint.
    m_ptInPlane = Constraint::PointOnPlaneContact
       (planeBodyB, m_frame, followerBodyF, point_F);
    m_ptInPlane.setIsConditional(true);
    m_ptInPlane.setDisabledByDefault(true);
}

//------------------------------------------------------------------------------
//                            WHERE TO DISPLAY
//------------------------------------------------------------------------------
Vec3 PointPlaneContact::whereToDisplay(const State& state) const {
    return m_follower.findStationLocationInGround(state,m_point);
}


//------------------------------------------------------------------------------
//                              GET PERR
//------------------------------------------------------------------------------
Real PointPlaneContact::getPerr(const State& state) const
{   //return m_ptInPlane.getPositionError(state); 
    const Vec3 p = m_follower.findStationLocationInAnotherBody
                                (state, m_point, m_planeBody);
    return ~p*m_frame.z() - m_height;
}

//------------------------------------------------------------------------------
//                    GET CONTACT MULTIPLIER INDEX
//------------------------------------------------------------------------------
MultiplierIndex PointPlaneContact::
getContactMultiplierIndex(const State& s) const {
    int mp, mv, ma;
    MultiplierIndex px0, vx0, ax0;
    m_ptInPlane.getNumConstraintEquationsInUse(s,mp,mv,ma);
    assert(mp==1 && mv==2 && ma==0); // don't call if not enabled
    m_ptInPlane.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    assert(px0.isValid() && vx0.isValid() && !ax0.isValid());
    return px0;
}

//------------------------------------------------------------------------------
//                  GET FRICTION MULTIPLIER INDICES
//------------------------------------------------------------------------------
void PointPlaneContact::
getFrictionMultiplierIndices(const State&       s, 
                             MultiplierIndex&   ix_x, 
                             MultiplierIndex&   ix_y) const
{   ix_x.invalidate(); ix_y.invalidate(); 
    int mp, mv, ma;
    MultiplierIndex px0, vx0, ax0;
    m_ptInPlane.getNumConstraintEquationsInUse(s,mp,mv,ma);
    assert(mp==1 && mv==2 && ma==0); // don't call if not enabled
    m_ptInPlane.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    assert(px0.isValid() && vx0.isValid() && !ax0.isValid());
    ix_x = vx0;
    ix_y = MultiplierIndex(vx0+1);
}



//==============================================================================
//                           SPHERE PLANE CONTACT
//==============================================================================
SpherePlaneContact::SpherePlaneContact
   (MobilizedBody& planeBodyB, const UnitVec3& normal_B, Real height,
    MobilizedBody& followerBodyF, const Vec3& point_F, Real radius,
    Real minCOR, Real mu_s, Real mu_d, Real mu_v)
:   m_minCOR(minCOR), m_mu_s(mu_s), m_mu_d(mu_d), m_mu_v(mu_v)
{
    SimTK_ERRCHK1_ALWAYS(0<=minCOR && minCOR<=1,
        "SpherePlaneContact()", 
        "The coefficient of restitution must be between 0 and 1 but was %g.", 
        minCOR);
    SimTK_ERRCHK3_ALWAYS(mu_s >= 0 && mu_d >= 0 && mu_v >= 0,
        "SpherePlaneContact()", 
        "All coefficients of friction must be nonnegative; got "
        "mu_s=%g, mu_d=%g, mu_v=%g.", mu_s, mu_d, mu_v);
    SimTK_ERRCHK2_ALWAYS(mu_d <= mu_s,
        "SpherePlaneContact()", 
        "The dynamic coefficient of friction can't be larger than "
        "the static coefficient; got mu_s=%g, mu_d=%g.", mu_s, mu_d);

    // Set up the contact constraint.
    const Rotation frame(normal_B, ZAxis);
    m_sphereOnPlane = Constraint::SphereOnPlaneContact
       (planeBodyB, frame, followerBodyF, point_F, radius, true);
    m_sphereOnPlane.setIsConditional(true);
    m_sphereOnPlane.setDisabledByDefault(true);
}

//------------------------------------------------------------------------------
//                            WHERE TO DISPLAY
//------------------------------------------------------------------------------
Vec3 SpherePlaneContact::whereToDisplay(const State& state) const {
    return m_sphereOnPlane.findContactPointInG(state);
}


//------------------------------------------------------------------------------
//                              GET PERR
//------------------------------------------------------------------------------
Real SpherePlaneContact::getPerr(const State& state) const
{   //return m_ptInPlane.getPositionError(state); 
    return m_sphereOnPlane.findSeparation(state);
}

//------------------------------------------------------------------------------
//                    GET CONTACT MULTIPLIER INDEX
//------------------------------------------------------------------------------
MultiplierIndex SpherePlaneContact::
getContactMultiplierIndex(const State& s) const {
    int mp, mv, ma;
    MultiplierIndex px0, vx0, ax0;
    m_sphereOnPlane.getNumConstraintEquationsInUse(s,mp,mv,ma);
    assert(mp==1 && mv==2 && ma==0); // don't call if not enabled
    m_sphereOnPlane.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    assert(px0.isValid() && vx0.isValid() && !ax0.isValid());
    return px0;
}

//------------------------------------------------------------------------------
//                  GET FRICTION MULTIPLIER INDICES
//------------------------------------------------------------------------------
void SpherePlaneContact::
getFrictionMultiplierIndices(const State&       s, 
                             MultiplierIndex&   ix_x, 
                             MultiplierIndex&   ix_y) const
{   ix_x.invalidate(); ix_y.invalidate(); 
    int mp, mv, ma;
    MultiplierIndex px0, vx0, ax0;
    m_sphereOnPlane.getNumConstraintEquationsInUse(s,mp,mv,ma);
    assert(mp==1 && mv==2 && ma==0); // don't call if not enabled
    m_sphereOnPlane.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    assert(px0.isValid() && vx0.isValid() && !ax0.isValid());
    ix_x = vx0;
    ix_y = MultiplierIndex(vx0+1);
}



//==============================================================================
//                           SPHERE SPHERE CONTACT
//==============================================================================
SphereSphereContact::SphereSphereContact
   (MobilizedBody&      mobod_F, 
    const Vec3&         defaultCenterOnF, 
    Real                defaultRadiusOnF, 
    MobilizedBody&      mobod_B, 
    const Vec3&         defaultCenterOnB,
    Real                defaultRadiusOnB,
    Real minCOR, Real mu_s, Real mu_d, Real mu_v)
:   m_minCOR(minCOR), m_mu_s(mu_s), m_mu_d(mu_d), m_mu_v(mu_v)
{
    SimTK_ERRCHK1_ALWAYS(0<=minCOR && minCOR<=1,
        "SphereSphereContact()", 
        "The coefficient of restitution must be between 0 and 1 but was %g.", 
        minCOR);
    SimTK_ERRCHK3_ALWAYS(mu_s >= 0 && mu_d >= 0 && mu_v >= 0,
        "SphereSphereContact()", 
        "All coefficients of friction must be nonnegative; got "
        "mu_s=%g, mu_d=%g, mu_v=%g.", mu_s, mu_d, mu_v);
    SimTK_ERRCHK2_ALWAYS(mu_d <= mu_s,
        "SphereSphereContact()", 
        "The dynamic coefficient of friction can't be larger than "
        "the static coefficient; got mu_s=%g, mu_d=%g.", mu_s, mu_d);

    // Set up the contact constraint.
    m_sphereOnSphere = Constraint::SphereOnSphereContact
       (mobod_F, defaultCenterOnF, defaultRadiusOnF, 
       mobod_B, defaultCenterOnB, defaultRadiusOnB, true);
    m_sphereOnSphere.setIsConditional(true);
    m_sphereOnSphere.setDisabledByDefault(true);
}

//------------------------------------------------------------------------------
//                            WHERE TO DISPLAY
//------------------------------------------------------------------------------
Vec3 SphereSphereContact::whereToDisplay(const State& state) const {
    return m_sphereOnSphere.findContactFrameInG(state).p();
}


//------------------------------------------------------------------------------
//                              GET PERR
//------------------------------------------------------------------------------
Real SphereSphereContact::getPerr(const State& state) const
{   
    return m_sphereOnSphere.findSeparation(state);
}

//------------------------------------------------------------------------------
//                    GET CONTACT MULTIPLIER INDEX
//------------------------------------------------------------------------------
MultiplierIndex SphereSphereContact::
getContactMultiplierIndex(const State& s) const {
    int mp, mv, ma;
    MultiplierIndex px0, vx0, ax0;
    m_sphereOnSphere.getNumConstraintEquationsInUse(s,mp,mv,ma);
    assert(mp==1 && mv==2 && ma==0); // don't call if not enabled
    m_sphereOnSphere.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    assert(px0.isValid() && vx0.isValid() && !ax0.isValid());
    return px0;
}

//------------------------------------------------------------------------------
//                  GET FRICTION MULTIPLIER INDICES
//------------------------------------------------------------------------------
void SphereSphereContact::
getFrictionMultiplierIndices(const State&       s, 
                             MultiplierIndex&   ix_x, 
                             MultiplierIndex&   ix_y) const
{   ix_x.invalidate(); ix_y.invalidate(); 
    int mp, mv, ma;
    MultiplierIndex px0, vx0, ax0;
    m_sphereOnSphere.getNumConstraintEquationsInUse(s,mp,mv,ma);
    assert(mp==1 && mv==2 && ma==0); // don't call if not enabled
    m_sphereOnSphere.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    assert(px0.isValid() && vx0.isValid() && !ax0.isValid());
    ix_x = vx0;
    ix_y = MultiplierIndex(vx0+1);
}

} // namespace SimTK

