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

/* Implementation of non-inline methods of the handle class
Constraint::SphereOnPlaneContact, and its implementation class
Constraint::SphereOnPlaneContactImpl. */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"
#include "simbody/internal/Constraint_SphereOnPlaneContact.h"

#include "Constraint_SphereOnPlaneContactImpl.h"
#include "SimbodyMatterSubsystemRep.h"

namespace SimTK {


//==============================================================================
//                         SPHERE ON PLANE CONTACT
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::SphereOnPlaneContact, 
                                        Constraint::SphereOnPlaneContactImpl, 
                                        Constraint);

Constraint::SphereOnPlaneContact::SphereOnPlaneContact
   (MobilizedBody&      planeBody, 
    const Transform&    defaultPlaneFrame, 
    MobilizedBody&      sphereBody, 
    const Vec3&         defaultSphereCenter,
    Real                defaultSphereRadius,
    bool                enforceRolling)
:   Constraint(new SphereOnPlaneContactImpl(enforceRolling))
{
    SimTK_ASSERT_ALWAYS(planeBody.isInSubsystem() && sphereBody.isInSubsystem(),
        "Constraint::SphereOnPlaneContact(): both bodies must already be in a "
        "SimbodyMatterSubsystem.");
    SimTK_ASSERT_ALWAYS(planeBody.isInSameSubsystem(sphereBody),
        "Constraint::SphereOnPlaneContact(): both bodies to be connected must be "
        "in the same SimbodyMatterSubsystem.");

    planeBody.updMatterSubsystem().adoptConstraint(*this);

    updImpl().m_planeBody_F     = updImpl().addConstrainedBody(planeBody);
    updImpl().m_ballBody_B      = updImpl().addConstrainedBody(sphereBody);
    updImpl().m_def_X_FP        = defaultPlaneFrame;
    updImpl().m_def_p_BO        = defaultSphereCenter;
    updImpl().m_def_radius      = defaultSphereRadius;
}

Constraint::SphereOnPlaneContact& Constraint::SphereOnPlaneContact::
setDefaultPlaneFrame(const Transform& defaultPlaneFrame) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_X_FP = defaultPlaneFrame;
    return *this;
}

Constraint::SphereOnPlaneContact& Constraint::SphereOnPlaneContact::
setDefaultSphereCenter(const Vec3& defaultSphereCenter) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_p_BO = defaultSphereCenter;
    return *this;
}

Constraint::SphereOnPlaneContact& Constraint::SphereOnPlaneContact::
setDefaultSphereRadius(Real defaultSphereRadius) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_radius = defaultSphereRadius;
    return *this;
}

const MobilizedBody& Constraint::SphereOnPlaneContact::
getPlaneMobilizedBody() const {
    const SphereOnPlaneContactImpl& impl = getImpl();
    return impl.getMobilizedBodyFromConstrainedBody(impl.m_planeBody_F);
}
const MobilizedBody& Constraint::SphereOnPlaneContact::
getSphereMobilizedBody() const {
    const SphereOnPlaneContactImpl& impl = getImpl();
    return impl.getMobilizedBodyFromConstrainedBody(impl.m_ballBody_B);
}

bool Constraint::SphereOnPlaneContact::isEnforcingRolling() const 
{   return getImpl().m_enforceRolling; }

const Transform& Constraint::SphereOnPlaneContact::getDefaultPlaneFrame() const {
    return getImpl().m_def_X_FP;
}

const Vec3& Constraint::SphereOnPlaneContact::getDefaultSphereCenter() const {
    return getImpl().m_def_p_BO;
}
Real Constraint::SphereOnPlaneContact::getDefaultSphereRadius() const {
    return getImpl().m_def_radius;
}

Constraint::SphereOnPlaneContact& Constraint::SphereOnPlaneContact::
setPlaneDisplayHalfWidth(Real h) {
    updImpl().setPlaneDisplayHalfWidth(h);
    return *this;
}

Real Constraint::SphereOnPlaneContact::getPlaneDisplayHalfWidth() const {
    return getImpl().getPlaneDisplayHalfWidth();
}

const Constraint::SphereOnPlaneContact& Constraint::SphereOnPlaneContact::
setPlaneFrame(State& state, const Transform& planeFrame) const {
    getImpl().updParameters(state).m_X_FP = planeFrame;
    return *this;
}

const Constraint::SphereOnPlaneContact& Constraint::SphereOnPlaneContact::
setSphereCenter(State& state, const Vec3& sphereCenter) const {
    getImpl().updParameters(state).m_p_BO = sphereCenter;
    return *this;
}

const Constraint::SphereOnPlaneContact& Constraint::SphereOnPlaneContact::
setSphereRadius(State& state, Real sphereRadius) const {
    getImpl().updParameters(state).m_radius = sphereRadius;
    return *this;
}

const Transform& Constraint::SphereOnPlaneContact::
getPlaneFrame(const State& state) const
{   return getImpl().getParameters(state).m_X_FP; }
const Vec3& Constraint::SphereOnPlaneContact::
getSphereCenter(const State& state) const
{   return getImpl().getParameters(state).m_p_BO; }
Real Constraint::SphereOnPlaneContact::
getSphereRadius(const State& state) const
{   return getImpl().getParameters(state).m_radius; }

Real Constraint::SphereOnPlaneContact::getPositionError(const State& s) const {
    Real perr;
    getImpl().getPositionErrors(s, 1, &perr);
    return perr;
}

Vec3 Constraint::SphereOnPlaneContact::getVelocityErrors(const State& s) const {
    const SphereOnPlaneContactImpl& impl = getImpl();
    Vec3 verr_PC; // result is velocity error in P frame 
    if (impl.m_enforceRolling) {
        Real verr[3];
        impl.getVelocityErrors(s, 3, verr);
        verr_PC = Vec3(verr[1],verr[2],verr[0]); // switch to x,y,z order
    } else {
        Real pverr;
        getImpl().getVelocityErrors(s, 1, &pverr);
        verr_PC = Vec3(0,0,pverr); // lone error is in z direction
    }
    return verr_PC;
}

Vec3 Constraint::SphereOnPlaneContact::getAccelerationErrors(const State& s) const {
    const SphereOnPlaneContactImpl& impl = getImpl();
    Vec3 aerr_PC; // result is acceleration error in P frame 
    if (impl.m_enforceRolling) {
        Real aerr[3];
        impl.getAccelerationErrors(s, 3, aerr);
        aerr_PC = Vec3(aerr[1],aerr[2],aerr[0]); // switch to x,y,z order
    } else {
        Real paerr;
        getImpl().getAccelerationErrors(s, 1, &paerr);
        aerr_PC = Vec3(0,0,paerr); // lone error is in z direction
    }
    return aerr_PC;
}

Vec3 Constraint::SphereOnPlaneContact::getMultipliers(const State& s) const {
    const SphereOnPlaneContactImpl& impl = getImpl();
    Vec3 lambda_PC; // result is -force on point F in P frame 
    if (impl.m_enforceRolling) {
        Real lambda[3];
        impl.getMultipliers(s, 3, lambda);
        lambda_PC = Vec3(lambda[1],lambda[2],lambda[0]); //switch to x,y,z order
    } else {
        Real lambda;
        getImpl().getMultipliers(s, 1, &lambda);
        lambda_PC = Vec3(0,0,lambda); // lone force is in z direction
    }
    return lambda_PC;
}

Vec3 Constraint::SphereOnPlaneContact::
findForceOnSphereInG(const State& state) const {
    const SphereOnPlaneContactImpl& impl = getImpl();
    if (impl.isDisabled(state)) 
        return Vec3(0);

    const Rotation& R_FP = impl.getParameters(state).m_X_FP.R();
    const MobilizedBody& bodyF =
        impl.getMobilizedBodyFromConstrainedBody(impl.m_planeBody_F);
    const Rotation& R_GF = bodyF.getBodyRotation(state);
    const Rotation R_GP = R_GF * R_FP; // orientation of P frame in G

    const Vec3 f_PC = -getMultipliers(state); // watch sign convention
    return R_GP * f_PC; // return f_GC
}


Vec3 Constraint::SphereOnPlaneContact::
findContactPointInG(const State& s) const {
    const SphereOnPlaneContactImpl& impl = getImpl();

    // Get the plane normal direction in G.
    const SphereOnPlaneContactImpl::Parameters& params = impl.getParameters(s);
    const UnitVec3&   Pz_F = params.m_X_FP.z();
    const Vec3&       p_BO = params.m_p_BO;
    const Real        r    = params.m_radius;

    const MobilizedBody& bodyF =
        impl.getMobilizedBodyFromConstrainedBody(impl.m_planeBody_F);
    const Rotation& R_GF = bodyF.getBodyRotation(s);
    const UnitVec3 Pz_G = R_GF * Pz_F;

    // Get the sphere center O in G.
    const MobilizedBody& bodyB =
        impl.getMobilizedBodyFromConstrainedBody(impl.m_ballBody_B);
    const Transform& X_GB = bodyB.getBodyTransform(s);
    const Vec3 p_GO = X_GB * p_BO;

    return p_GO - r * Pz_G;
}

Real Constraint::SphereOnPlaneContact::
findSeparation(const State& s) const {
    const SphereOnPlaneContactImpl& impl = getImpl();

    // The height of the sphere center should be its radius.
    const SphereOnPlaneContactImpl::Parameters& params = impl.getParameters(s);
    const Transform&  X_FP = params.m_X_FP;
    const Vec3&       p_BO = params.m_p_BO;
    const Real        r    = params.m_radius;

    const UnitVec3&   Pz_F = X_FP.z();
    const Vec3&       p_FP = X_FP.p();

    const MobilizedBody& bodyF =
        impl.getMobilizedBodyFromConstrainedBody(impl.m_planeBody_F);
    const MobilizedBody& bodyB =
        impl.getMobilizedBodyFromConstrainedBody(impl.m_ballBody_B);

    const Vec3 p_FO = bodyB.findStationLocationInAnotherBody(s,p_BO,bodyF);
    const Vec3 p_PO_F = p_FO - p_FP;
    return dot(p_PO_F, Pz_F) - r;
}


//==============================================================================
//                     SPHERE ON PLANE CONTACT IMPL
//==============================================================================

// The default plane and sphere parameters may be overridden by setting
// a discrete variable in the state. We allocate the state resources here.
void Constraint::SphereOnPlaneContactImpl::
realizeTopologyVirtual(State& state) const {
    parametersIx = getMyMatterSubsystemRep().
        allocateDiscreteVariable(state, Stage::Position, 
            new Value<Parameters>
               (Parameters(m_def_X_FP, m_def_p_BO, m_def_radius)));
}

// Return the pair of constrained station points, with the first expressed 
// in the body 1 frame and the second in the body 2 frame. Note that although
// these are used to define the position error, only the station on body 2
// is used to generate constraint forces; the point of body 1 that is 
// coincident with the body 2 point receives the equal and opposite force.
const Constraint::SphereOnPlaneContactImpl::Parameters& 
Constraint::SphereOnPlaneContactImpl::
getParameters(const State& state) const {
    return Value<Parameters>::downcast
       (getMyMatterSubsystemRep().getDiscreteVariable(state,parametersIx));
}

// Return a writable reference into the Instance-stage state variable 
// containing the pair of constrained station points, with the first expressed 
// in the body 1 frame and the second in the body 2 frame. Calling this
// method invalidates the Instance stage and above in the given state.
Constraint::SphereOnPlaneContactImpl::Parameters& 
Constraint::SphereOnPlaneContactImpl::
updParameters(State& state) const {
    return Value<Parameters>::updDowncast
       (getMyMatterSubsystemRep().updDiscreteVariable(state,parametersIx));
}

void Constraint::SphereOnPlaneContactImpl::
calcDecorativeGeometryAndAppendVirtual
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
    // We can't generate the artwork until we know the plane frame and the
    // sphere center and radius, which might not be until Position stage.
    if (   stage == Stage::Position 
        && getMyMatterSubsystemRep().getShowDefaultGeometry()) 
    {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        const Parameters& params = getParameters(s);
        const Transform& X_FP = params.m_X_FP;
        const Vec3&      p_BO = params.m_p_BO;
        const Real       r    = params.m_radius;

        // TODO: should be instance-stage data from State rather than 
        // topological data.
        // This makes z axis point along plane normal

        const MobilizedBodyIndex planeMBIx = 
            getMobilizedBodyIndexOfConstrainedBody(m_planeBody_F);
        const MobilizedBodyIndex ballMBIx = 
            getMobilizedBodyIndexOfConstrainedBody(m_ballBody_B);

        if (m_planeHalfWidth > 0) {
            // On the inboard body, draw a gray transparent rectangle, 
            // outlined in black lines.
            geom.push_back(DecorativeBrick
               (Vec3(m_planeHalfWidth,m_planeHalfWidth,r/10))
                .setColor(Gray)
                .setRepresentation(DecorativeGeometry::DrawSurface)
                .setOpacity(Real(0.3))
                .setBodyId(planeMBIx)
                .setTransform(X_FP));
            geom.push_back(DecorativeFrame(m_planeHalfWidth/5)
                           .setColor(Gray)
                           .setBodyId(planeMBIx)
                           .setTransform(X_FP));
        }

        // On the ball body draw an orange mesh sphere.
        geom.push_back(DecorativeSphere(r)
            .setColor(Orange)
            .setRepresentation(DecorativeGeometry::DrawWireframe)
            .setBodyId(ballMBIx)
            .setTransform(p_BO));
    }
}


} // namespace SimTK

