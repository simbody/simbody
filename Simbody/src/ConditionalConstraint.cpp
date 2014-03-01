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
//                           POINT PLANE CONTACT
//==============================================================================
PointPlaneContact::PointPlaneContact
   (MobilizedBody& planeBodyB, const UnitVec3& normal_B, Real height,
    MobilizedBody& followerBodyF, const Vec3& point_F, 
    Real minCOR, Real mu_s, Real mu_d, Real mu_v)
:   m_planeBody(planeBodyB), m_frame(normal_B, ZAxis), m_height(height), 
    m_follower(followerBodyF), m_point(point_F), m_minCOR(minCOR),
    m_mu_s(mu_s), m_mu_d(mu_d), m_mu_v(mu_v),
    m_hasFriction(m_mu_s>0 || m_mu_d>0 || m_mu_v>0)
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
    m_ptInPlane = Constraint::PointInPlane
       (planeBodyB, normal_B, height, followerBodyF, point_F);
    m_ptInPlane.setDisabledByDefault(true);

    // Create the friction constraints if needed.
    if (m_hasFriction) {
        m_noslipX = Constraint::NoSlip1D
           (m_planeBody, Vec3(0), // we'll change this point later 
            m_frame.x(), m_planeBody, m_follower);
        m_noslipY = Constraint::NoSlip1D
           (m_planeBody, Vec3(0), // we'll change this point later 
            m_frame.y(), m_planeBody, m_follower);
        m_noslipX.setDisabledByDefault(true);
        m_noslipY.setDisabledByDefault(true);
    }
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
    assert(mp==1 && mv==0 && ma==0); // don't call if not enabled
    m_ptInPlane.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    assert(px0.isValid() && !vx0.isValid() && !ax0.isValid());
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
    if (!m_hasFriction) return;
    int mp, mv, ma;
    MultiplierIndex px0, vx0, ax0;
    m_noslipX.getNumConstraintEquationsInUse(s,mp,mv,ma);
    assert(mp==0 && mv==1 && ma==0); // don't call if not enabled
    m_noslipX.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    assert(!px0.isValid() && vx0.isValid() && !ax0.isValid());
    ix_x = vx0;
    m_noslipY.getNumConstraintEquationsInUse(s,mp,mv,ma);
    assert(mp==0 && mv==1 && ma==0); // don't call if not enabled
    m_noslipY.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
    assert(!px0.isValid() && vx0.isValid() && !ax0.isValid());
    ix_y = vx0;  
}

//------------------------------------------------------------------------------
//                 GET POSITION INFO / SET INSTANCE PARAMETER
//------------------------------------------------------------------------------
// Return the contact point in the plane body so we can make friction 
// act there.
Vec3 PointPlaneContact::
getPositionInfo(const State& state) const {
    return m_follower.findStationLocationInAnotherBody
                                (state, m_point, m_planeBody); 
}

// Set the friction contact point to the contact point in the plane.
void PointPlaneContact::
setInstanceParameter(State& state, const Vec3& pos) const {
    if (m_hasFriction) {
        m_noslipX.setContactPoint(state, pos);
        m_noslipY.setContactPoint(state, pos);
    }
}


} // namespace SimTK

