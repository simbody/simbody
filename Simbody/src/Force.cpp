/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
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

#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/Force.h"
#include "simbody/internal/Force_BuiltIns.h"

#include "ForceImpl.h"

namespace SimTK {

//------------------------------------------------------------------------------
//                                  FORCE
//------------------------------------------------------------------------------

void Force::setDisabledByDefault(bool shouldBeDisabled)
{   updImpl().setDisabledByDefault(shouldBeDisabled); }
bool Force::isDisabledByDefault() const
{   return getImpl().isDisabledByDefault(); }

void Force::disable(State& state) const
{   getForceSubsystem().setForceIsDisabled(state, getForceIndex(), true); }
void Force::enable(State& state) const
{   getForceSubsystem().setForceIsDisabled(state, getForceIndex(), false); }
bool Force::isDisabled(const State& state) const
{   return getForceSubsystem().isForceDisabled(state, getForceIndex()); }

void Force::calcForceContribution(const State&   state,
                           Vector_<SpatialVec>&  bodyForces,
                           Vector_<Vec3>&        particleForces,
                           Vector&               mobilityForces) const
{
    const MultibodySystem& mbs = getForceSubsystem().getMultibodySystem();
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    // Resize if necessary.
    bodyForces.resize(matter.getNumBodies());
    particleForces.resize(matter.getNumParticles());
    mobilityForces.resize(matter.getNumMobilities());

    // Set all forces to zero.
    bodyForces.setToZero();
    particleForces.setToZero();
    mobilityForces.setToZero();
    if (isDisabled(state)) return;

    // Add in force element contributions.
    getImpl().calcForce(state,bodyForces,particleForces,mobilityForces);
}

Real Force::calcPotentialEnergyContribution(const State& state) const {
    if (isDisabled(state)) return 0;
    return getImpl().calcPotentialEnergy(state);
}

const GeneralForceSubsystem& Force::getForceSubsystem() const
{   return getImpl().getForceSubsystem(); }
ForceIndex Force::getForceIndex() const
{   return getImpl().getForceIndex(); }

//-------------------------- TwoPointLinearSpring ------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::TwoPointLinearSpring, Force::TwoPointLinearSpringImpl, Force);

Force::TwoPointLinearSpring::TwoPointLinearSpring(GeneralForceSubsystem& forces, const MobilizedBody& body1, const Vec3& station1,
        const MobilizedBody& body2, const Vec3& station2, Real k, Real x0) : Force(new TwoPointLinearSpringImpl(
        body1, station1, body2, station2, k, x0)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::TwoPointLinearSpringImpl::TwoPointLinearSpringImpl(const MobilizedBody& body1, const Vec3& station1,
        const MobilizedBody& body2, const Vec3& station2, Real k, Real x0) : matter(body1.getMatterSubsystem()),
        body1(body1.getMobilizedBodyIndex()), station1(station1),
        body2(body2.getMobilizedBodyIndex()), station2(station2), k(k), x0(x0) {
}

void Force::TwoPointLinearSpringImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const Transform& X_GB1 = matter.getMobilizedBody(body1).getBodyTransform(state);
    const Transform& X_GB2 = matter.getMobilizedBody(body2).getBodyTransform(state);

    const Vec3 s1_G = X_GB1.R() * station1;
    const Vec3 s2_G = X_GB2.R() * station2;

    const Vec3 p1_G = X_GB1.p() + s1_G; // station measured from ground origin
    const Vec3 p2_G = X_GB2.p() + s2_G;

    const Vec3 r_G       = p2_G - p1_G; // vector from point1 to point2
    const Real d         = r_G.norm();  // distance between the points
    const Real stretch   = d - x0;      // + -> tension, - -> compression
    const Real frcScalar = k*stretch;   // k(x-x0)

    const Vec3 f1_G = (frcScalar/d) * r_G;
    bodyForces[body1] +=  SpatialVec(s1_G % f1_G, f1_G);
    bodyForces[body2] -=  SpatialVec(s2_G % f1_G, f1_G);
}

Real Force::TwoPointLinearSpringImpl::calcPotentialEnergy(const State& state) const {
    const Transform& X_GB1 = matter.getMobilizedBody(body1).getBodyTransform(state);
    const Transform& X_GB2 = matter.getMobilizedBody(body2).getBodyTransform(state);

    const Vec3 s1_G = X_GB1.R() * station1;
    const Vec3 s2_G = X_GB2.R() * station2;

    const Vec3 p1_G = X_GB1.p() + s1_G; // station measured from ground origin
    const Vec3 p2_G = X_GB2.p() + s2_G;

    const Vec3 r_G = p2_G - p1_G; // vector from point1 to point2
    const Real d   = r_G.norm();  // distance between the points
    const Real stretch   = d - x0; // + -> tension, - -> compression

    return k*stretch*stretch/2; // 1/2 k (x-x0)^2
}


void Force::TwoPointLinearSpringImpl::
calcDecorativeGeometryAndAppend(const State& state, Stage stage,
                                Array_<DecorativeGeometry>& geom) const {
    if (stage==Stage::Position && matter.getShowDefaultGeometry()) {
        const Transform& X_GB1 = matter.getMobilizedBody(body1)
                                       .getBodyTransform(state);
        const Transform& X_GB2 = matter.getMobilizedBody(body2)
                                       .getBodyTransform(state);

        const Vec3 s1_G = X_GB1.R() * station1;
        const Vec3 s2_G = X_GB2.R() * station2;

        const Vec3 p1_G = X_GB1.p() + s1_G; // measured from ground origin
        const Vec3 p2_G = X_GB2.p() + s2_G;

        geom.push_back(DecorativeLine(p1_G,p2_G).setColor(Orange));
    }
}


//-------------------------- TwoPointLinearDamper ------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::TwoPointLinearDamper, Force::TwoPointLinearDamperImpl, Force);

Force::TwoPointLinearDamper::TwoPointLinearDamper(GeneralForceSubsystem& forces, const MobilizedBody& body1, const Vec3& station1,
        const MobilizedBody& body2, const Vec3& station2, Real damping) : Force(new TwoPointLinearDamperImpl(
        body1, station1, body2, station2, damping)) {
    assert(damping >= 0);
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::TwoPointLinearDamperImpl::TwoPointLinearDamperImpl(const MobilizedBody& body1, const Vec3& station1,
        const MobilizedBody& body2, const Vec3& station2, Real damping) : matter(body1.getMatterSubsystem()),
        body1(body1.getMobilizedBodyIndex()), station1(station1),
        body2(body2.getMobilizedBodyIndex()), station2(station2), damping(damping) {
}

void Force::TwoPointLinearDamperImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const Transform& X_GB1 = matter.getMobilizedBody(body1).getBodyTransform(state);
    const Transform& X_GB2 = matter.getMobilizedBody(body2).getBodyTransform(state);

    const Vec3 s1_G = X_GB1.R() * station1;
    const Vec3 s2_G = X_GB2.R() * station2;

    const Vec3 p1_G = X_GB1.p() + s1_G; // station measured from ground origin
    const Vec3 p2_G = X_GB2.p() + s2_G;

    const Vec3 v1_G = matter.getMobilizedBody(body1).findStationVelocityInGround(state, station1);
    const Vec3 v2_G = matter.getMobilizedBody(body2).findStationVelocityInGround(state, station2);
    const Vec3 vRel = v2_G - v1_G; // relative velocity

    const UnitVec3 d(p2_G - p1_G); // direction from point1 to point2
    const Real frc = damping*dot(vRel, d); // c*v

    const Vec3 f1_G = frc*d;
    bodyForces[body1] +=  SpatialVec(s1_G % f1_G, f1_G);
    bodyForces[body2] -=  SpatialVec(s2_G % f1_G, f1_G);
}

Real Force::TwoPointLinearDamperImpl::calcPotentialEnergy(const State& state) const {
    return 0;
}


//-------------------------- TwoPointConstantForce -----------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::TwoPointConstantForce, Force::TwoPointConstantForceImpl, Force);

Force::TwoPointConstantForce::TwoPointConstantForce(GeneralForceSubsystem& forces, const MobilizedBody& body1, const Vec3& station1,
        const MobilizedBody& body2, const Vec3& station2, Real force) : Force(new TwoPointConstantForceImpl(
        body1, station1, body2, station2, force)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::TwoPointConstantForceImpl::TwoPointConstantForceImpl(const MobilizedBody& body1, const Vec3& station1,
        const MobilizedBody& body2, const Vec3& station2, Real force) : matter(body1.getMatterSubsystem()),
        body1(body1.getMobilizedBodyIndex()), station1(station1),
        body2(body2.getMobilizedBodyIndex()), station2(station2), force(force) {
}

void Force::TwoPointConstantForceImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const Transform& X_GB1 = matter.getMobilizedBody(body1).getBodyTransform(state);
    const Transform& X_GB2 = matter.getMobilizedBody(body2).getBodyTransform(state);

    const Vec3 s1_G = X_GB1.R() * station1;
    const Vec3 s2_G = X_GB2.R() * station2;

    const Vec3 p1_G = X_GB1.p() + s1_G; // station measured from ground origin
    const Vec3 p2_G = X_GB2.p() + s2_G;

    const Vec3 r_G = p2_G - p1_G; // vector from point1 to point2
    const Real x   = r_G.norm();  // distance between the points
    const UnitVec3 d(r_G/x, true);

    const Vec3 f2_G = force * d;
    bodyForces[body1] -=  SpatialVec(s1_G % f2_G, f2_G);
    bodyForces[body2] +=  SpatialVec(s2_G % f2_G, f2_G);
}

Real Force::TwoPointConstantForceImpl::calcPotentialEnergy(const State& state) const {
    return 0;
}


//--------------------------- MobilityLinearSpring -----------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::MobilityLinearSpring,
                                        Force::MobilityLinearSpringImpl,
                                        Force);

Force::MobilityLinearSpring::
MobilityLinearSpring(GeneralForceSubsystem&  forces,
                     const MobilizedBody&    mobod,
                     MobilizerQIndex         whichQ,
                     Real                    defaultStiffness,
                     Real                    defaultQZero)
:   Force(new MobilityLinearSpringImpl(mobod, whichQ,
                                       defaultStiffness, defaultQZero))
{
    SimTK_ERRCHK1_ALWAYS(defaultStiffness >= 0,
        "Force::MobilityLinearSpring::MobilityLinearSpring()",
        "Stiffness coefficient must be nonnegative "
        "(defaultStiffness=%g).", defaultStiffness);

    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::MobilityLinearSpring& Force::MobilityLinearSpring::
setDefaultStiffness(Real defaultStiffness) {
    SimTK_ERRCHK1_ALWAYS(defaultStiffness >= 0,
        "Force::MobilityLinearSpring::setDefaultStiffness()",
        "Stiffness coefficient must be nonnegative "
        "(defaultStiffness=%g).", defaultStiffness);

    MobilityLinearSpringImpl& impl = updImpl();
    if (impl.m_defaultStiffness != defaultStiffness) {
        impl.m_defaultStiffness = defaultStiffness;
        impl.invalidateTopologyCache();
    }
    return *this;
}


Force::MobilityLinearSpring& Force::MobilityLinearSpring::
setDefaultQZero(Real defaultQZero) {
    MobilityLinearSpringImpl& impl = updImpl();
    if (impl.m_defaultQZero != defaultQZero) {
        impl.m_defaultQZero = defaultQZero;
        impl.invalidateTopologyCache();
    }
    return *this;
}

Real Force::MobilityLinearSpring::
getDefaultStiffness() const
{   return getImpl().m_defaultStiffness; }

Real Force::MobilityLinearSpring::
getDefaultQZero() const
{   return getImpl().m_defaultQZero; }

const Force::MobilityLinearSpring& Force::MobilityLinearSpring::
setStiffness(State& state, Real stiffness) const {
    SimTK_ERRCHK1_ALWAYS(stiffness >= 0,
        "Force::MobilityLinearSpring::setStiffness()",
        "Stiffness coefficient must be nonnegative "
        "(stiffness=%g).", stiffness);

    getImpl().updParams(state).first = stiffness;
    return *this;
}

const Force::MobilityLinearSpring& Force::MobilityLinearSpring::
setQZero(State& state, Real qZero) const
{   getImpl().updParams(state).second = qZero; return *this; }

Real Force::MobilityLinearSpring::
getStiffness(const State& state) const
{   return getImpl().getParams(state).first; }

Real Force::MobilityLinearSpring::
getQZero(const State& state) const
{   return getImpl().getParams(state).second; }

Force::MobilityLinearSpringImpl::
MobilityLinearSpringImpl(const MobilizedBody&    mobod,
                         MobilizerQIndex         whichQ,
                         Real                    defaultStiffness,
                         Real                    defaultQZero)
:   m_matter(mobod.getMatterSubsystem()),
    m_mobodIx(mobod.getMobilizedBodyIndex()), m_whichQ(whichQ),
    m_defaultStiffness(defaultStiffness), m_defaultQZero(defaultQZero)
{
}

void Force::MobilityLinearSpringImpl::
calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
          Vector_<Vec3>& particleForces, Vector& mobilityForces) const
{
    const MobilizedBody& mb = m_matter.getMobilizedBody(m_mobodIx);
    const Real q = mb.getOneQ(state, m_whichQ);
    const std::pair<Real,Real>& params = getParams(state);
    const Real k = params.first, q0 = params.second;
    const Real frc = -k*(q-q0);
    // bug: this is depending on qdot=u
    mb.applyOneMobilityForce(state, MobilizerUIndex((int)m_whichQ),
                             frc, mobilityForces);
}

Real Force::MobilityLinearSpringImpl::
calcPotentialEnergy(const State& state) const {
    const MobilizedBody& mb = m_matter.getMobilizedBody(m_mobodIx);
    const Real q = mb.getOneQ(state, m_whichQ);
    const std::pair<Real,Real>& params = getParams(state);
    const Real k = params.first, q0 = params.second;
    const Real frc = -k*(q-q0);
    return k*square(q-q0)/2;
}



//--------------------------- MobilityLinearDamper -----------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::MobilityLinearDamper,
                                        Force::MobilityLinearDamperImpl,
                                        Force);

Force::MobilityLinearDamper::
MobilityLinearDamper(GeneralForceSubsystem& forces,
                     const MobilizedBody&   mobod,
                     MobilizerUIndex        whichU,
                     Real                   defaultDamping)
:   Force(new MobilityLinearDamperImpl(mobod, whichU, defaultDamping))
{
    SimTK_ERRCHK1_ALWAYS(defaultDamping >= 0,
        "Force::MobilityLinearDamper::MobilityLinearDamper()",
        "Damping coefficient must be nonnegative "
        "(defaultDamping=%g).", defaultDamping);

    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}


Force::MobilityLinearDamper& Force::MobilityLinearDamper::
setDefaultDamping(Real defaultDamping) {
    SimTK_ERRCHK1_ALWAYS(defaultDamping >= 0,
        "Force::MobilityLinearDamper::setDefaultDamping()",
        "Damping coefficient must be nonnegative "
        "(defaultDamping=%g).", defaultDamping);

    MobilityLinearDamperImpl& impl = updImpl();
    if (impl.m_defaultDamping != defaultDamping) {
        impl.m_defaultDamping = defaultDamping;
        impl.invalidateTopologyCache();
    }
    return *this;
}

Real Force::MobilityLinearDamper::
getDefaultDamping() const
{   return getImpl().m_defaultDamping; }

const Force::MobilityLinearDamper& Force::MobilityLinearDamper::
setDamping(State& state, Real damping) const
{
    SimTK_ERRCHK1_ALWAYS(damping >= 0,
        "Force::MobilityLinearDamper::setDamping()",
        "Damping coefficient must be nonnegative "
        "(damping=%g).", damping);

    getImpl().updDamping(state) = damping;
    return *this;
}

Real Force::MobilityLinearDamper::
getDamping(const State& state) const
{   return getImpl().getDamping(state); }


Force::MobilityLinearDamperImpl::
MobilityLinearDamperImpl(const MobilizedBody&   mobod,
                         MobilizerUIndex        whichU,
                         Real                   defaultDamping)
:   m_matter(mobod.getMatterSubsystem()),
    m_mobodIx(mobod.getMobilizedBodyIndex()), m_whichU(whichU),
    m_defaultDamping(defaultDamping)
{
}

void Force::MobilityLinearDamperImpl::
calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
          Vector_<Vec3>& particleForces, Vector& mobilityForces) const
{
    const MobilizedBody& mb = m_matter.getMobilizedBody(m_mobodIx);
    const Real u = mb.getOneU(state, m_whichU);
    const Real damping = getDamping(state);
    const Real frc = -damping*u;
    mb.applyOneMobilityForce(state, m_whichU, frc, mobilityForces);

    //Worthless loop to make the force calculations more expensive - make sure that
    //the parallelism is helping out
    volatile double junkInt = 0;
    for(double y = 0; y < 10000; y++)
        for(double x = 0; x < 99999; x++)
            junkInt = junkInt + .2 + 3.6 * 4 /2 + 3 % 2 *2321.234 -322488238382.124;
}

Real Force::MobilityLinearDamperImpl::
calcPotentialEnergy(const State& state) const {
    return 0;
}



//-------------------------- MobilityConstantForce -----------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::MobilityConstantForce,
                                        Force::MobilityConstantForceImpl,
                                        Force);

Force::MobilityConstantForce::MobilityConstantForce
   (GeneralForceSubsystem&  forces,
    const MobilizedBody&    mobod,
    MobilizerUIndex         whichU,
    Real                    defaultForce)
: Force(new MobilityConstantForceImpl(mobod, whichU, defaultForce)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::MobilityConstantForce& Force::MobilityConstantForce::
setDefaultForce(Real defaultForce) {
    MobilityConstantForceImpl& impl = updImpl();
    if (impl.m_defaultForce != defaultForce) {
        impl.m_defaultForce = defaultForce;
        impl.invalidateTopologyCache();
    }
    return *this;
}

Real Force::MobilityConstantForce::
getDefaultForce() const
{   return getImpl().m_defaultForce; }

void Force::MobilityConstantForce::
setForce(State& state, Real force) const {
    getImpl().updForce(state) = force;
}

Real Force::MobilityConstantForce::
getForce(const State& state) const
{   return getImpl().getForce(state); }

Force::MobilityConstantForceImpl::MobilityConstantForceImpl
   (const MobilizedBody&    mobod,
    MobilizerUIndex         whichU,
    Real                    defaultForce)
:   m_matter(mobod.getMatterSubsystem()),
    m_mobodIx(mobod.getMobilizedBodyIndex()), m_whichU(whichU),
    m_defaultForce(defaultForce)
{
}

void Force::MobilityConstantForceImpl::
calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
          Vector_<Vec3>& particleForces, Vector& mobilityForces) const
{
    const MobilizedBody& mb = m_matter.getMobilizedBody(m_mobodIx);
    mb.applyOneMobilityForce(state, m_whichU, getForce(state), mobilityForces);
}


//---------------------------- MobilityLinearStop ------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::MobilityLinearStop,
                                        Force::MobilityLinearStopImpl,
                                        Force);

Force::MobilityLinearStop::MobilityLinearStop
   (GeneralForceSubsystem&    forces,
    const MobilizedBody&      mobod,
    MobilizerQIndex           whichQ,
    Real                      defaultStiffness,
    Real                      defaultDissipation,
    Real                      defaultQLow,
    Real                      defaultQHigh)
:   Force(new MobilityLinearStopImpl(mobod,whichQ,
                                     defaultStiffness, defaultDissipation,
                                     defaultQLow, defaultQHigh))
{
    SimTK_ERRCHK2_ALWAYS(defaultStiffness >= 0 && defaultDissipation >= 0,
        "Force::MobilityLinearStop::MobilityLinearStop()",
        "Stiffness and dissipation coefficient must be nonnegative "
        "(defaultStiffness=%g, defaultDissipation=%g).",
        defaultStiffness, defaultDissipation);
    SimTK_ERRCHK2_ALWAYS(defaultQLow <= defaultQHigh,
        "Force::MobilityLinearStop::MobilityLinearStop()",
        "Lower bound can't be larger than upper bound "
        "(defaultQLow=%g, defaultQHigh=%g).",
        defaultQLow, defaultQHigh);

    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::MobilityLinearStop& Force::MobilityLinearStop::
setDefaultBounds(Real defaultQLow, Real defaultQHigh) {
    SimTK_ERRCHK2_ALWAYS(defaultQLow <= defaultQHigh,
        "Force::MobilityLinearStop::setDefaultBounds()",
        "Lower bound can't be larger than upper bound "
        "(defaultQLow=%g, defaultQHigh=%g).",
        defaultQLow, defaultQHigh);

    MobilityLinearStopImpl& impl = updImpl();
    if (impl.m_defQLow != defaultQLow || impl.m_defQHigh != defaultQHigh) {
        impl.m_defQLow = defaultQLow;
        impl.m_defQHigh = defaultQHigh;
        impl.invalidateTopologyCache();
    }
    return *this;
}

Force::MobilityLinearStop& Force::MobilityLinearStop::
setDefaultMaterialProperties(Real defaultStiffness, Real defaultDissipation) {
    SimTK_ERRCHK2_ALWAYS(defaultStiffness >= 0 && defaultDissipation >= 0,
        "Force::MobilityLinearStop::setDefaultMaterialProperties()",
        "Stiffness and dissipation coefficient must be nonnegative "
        "(defaultStiffness=%g, defaultDissipation=%g).",
        defaultStiffness, defaultDissipation);

    MobilityLinearStopImpl& impl = updImpl();
    if (   impl.m_defStiffness   != defaultStiffness
        || impl.m_defDissipation != defaultDissipation) {
        impl.m_defStiffness = defaultStiffness;
        impl.m_defDissipation = defaultDissipation;
        impl.invalidateTopologyCache();
    }
    return *this;
}

Real Force::MobilityLinearStop::getDefaultLowerBound() const
{   return getImpl().m_defQLow; }
Real Force::MobilityLinearStop::getDefaultUpperBound() const
{   return getImpl().m_defQHigh; }
Real Force::MobilityLinearStop::getDefaultStiffness() const
{   return getImpl().m_defStiffness; }
Real Force::MobilityLinearStop::getDefaultDissipation() const
{   return getImpl().m_defDissipation; }


void Force::MobilityLinearStop::
setBounds(State& state, Real qLow, Real qHigh) const {
    SimTK_ERRCHK2_ALWAYS(qLow <= qHigh,
        "Force::MobilityLinearStop::setBounds()",
        "Lower bound can't be larger than upper bound (qLow=%g, qHigh=%g).",
        qLow, qHigh);

    MobilityLinearStopImpl::Parameters& params =
        getImpl().updParameters(state); // invalidates Dynamics stage
    params.qLow = qLow; params.qHigh = qHigh;
}
void Force::MobilityLinearStop::
setMaterialProperties(State& state, Real stiffness, Real dissipation) const {
    SimTK_ERRCHK2_ALWAYS(stiffness >= 0 && dissipation >= 0,
        "Force::MobilityLinearStop::setMaterialProperties()",
        "Stiffness and dissipation coefficient must be nonnegative "
        "(stiffness=%g, dissipation=%g).",
        stiffness, dissipation);

    MobilityLinearStopImpl::Parameters& params =
        getImpl().updParameters(state); // invalidates Dynamics stage
    params.k = stiffness; params.d = dissipation;
}

Real Force::MobilityLinearStop::getLowerBound(const State& state) const
{   return getImpl().getParameters(state).qLow; }
Real Force::MobilityLinearStop::getUpperBound(const State& state) const
{   return getImpl().getParameters(state).qHigh; }
Real Force::MobilityLinearStop::getStiffness(const State& state) const
{   return getImpl().getParameters(state).k; }
Real Force::MobilityLinearStop::getDissipation(const State& state) const
{   return getImpl().getParameters(state).d; }


Force::MobilityLinearStopImpl::MobilityLinearStopImpl
   (const MobilizedBody&      mobod,
    MobilizerQIndex           whichQ,
    Real                      defaultStiffness,
    Real                      defaultDissipation,
    Real                      defaultQLow,
    Real                      defaultQHigh)
:   m_matter(mobod.getMatterSubsystem()),
    m_mobodIx(mobod.getMobilizedBodyIndex()), m_whichQ(whichQ),
    m_defStiffness(defaultStiffness), m_defDissipation(defaultDissipation),
    m_defQLow(defaultQLow), m_defQHigh(defaultQHigh)
{
}


void Force::MobilityLinearStopImpl::
calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
          Vector_<Vec3>& particleForces, Vector& mobilityForces) const
{
    const Parameters& param = getParameters(state);
    if (param.k == 0) return; // no stiffness, no force

    const MobilizedBody& mb = m_matter.getMobilizedBody(m_mobodIx);
    const Real q = mb.getOneQ(state, m_whichQ);

    // Don't ask for velocity-dependent qdot if there is no dissipation.
    const Real qdot = param.d != 0 ? mb.getOneQDot(state, m_whichQ)
                                   : Real(0);

    if (q > param.qHigh) {
        const Real x = q-param.qHigh;  // x > 0
        const Real fraw = param.k*x*(1+param.d*qdot); // should be >= 0
        mb.applyOneMobilityForce(state,
            MobilizerUIndex(m_whichQ), // TODO: only works qdot & u match
            std::min(Real(0), -fraw), mobilityForces);
    } else if (q < param.qLow) {
        const Real x = q-param.qLow;    // x < 0
        const Real fraw = param.k*x*(1-param.d*qdot); // should be <= 0
        mb.applyOneMobilityForce(state,
            MobilizerUIndex(m_whichQ), // TODO: only works qdot & u match
            std::max(Real(0), -fraw), mobilityForces);
    }
}

Real Force::MobilityLinearStopImpl::
calcPotentialEnergy(const State& state) const {
    const Parameters& param = getParameters(state);
    if (param.k == 0) return 0; // no stiffness, no energy stored

    const MobilizedBody& mb = m_matter.getMobilizedBody(m_mobodIx);
    const Real q = mb.getOneQ(state, m_whichQ);
    Real x = 0;
    if      (q > param.qHigh) x = q-param.qHigh;
    else if (q < param.qLow)  x = q-param.qLow;
    else return 0; // neither stop is engaged

    return param.k*x*x/2;  // 1/2 k x^2
}

//-------------------------- MobilityDiscreteForce -----------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::MobilityDiscreteForce,
                                        Force::MobilityDiscreteForceImpl,
                                        Force);

Force::MobilityDiscreteForce::MobilityDiscreteForce
   (GeneralForceSubsystem& forces, const MobilizedBody& mobod,
    MobilizerUIndex whichU, Real defaultForce)
:   Force(new MobilityDiscreteForceImpl(mobod, whichU, defaultForce)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::MobilityDiscreteForce& Force::MobilityDiscreteForce::
setDefaultMobilityForce(Real defaultForce) {
    updImpl().m_defaultVal = defaultForce;
    getImpl().invalidateTopologyCache();
    return *this;
}

Real Force::MobilityDiscreteForce::
getDefaultMobilityForce() const {
    return getImpl().m_defaultVal;
}

void Force::MobilityDiscreteForce::
setMobilityForce(State& state, Real f) const {
    getImpl().setMobilityForce(state, f);
}

Real Force::MobilityDiscreteForce::
getMobilityForce(const State& state) const {
    return getImpl().getMobilityForce(state);
}

Force::MobilityDiscreteForceImpl::MobilityDiscreteForceImpl
   (const MobilizedBody& mobod, MobilizerUIndex whichU, Real defaultForce)
:   m_matter(mobod.getMatterSubsystem()),
    m_mobodIx(mobod.getMobilizedBodyIndex()),
    m_whichU(whichU), m_defaultVal(defaultForce) {
}

void Force::MobilityDiscreteForceImpl::
setMobilityForce(State& state, Real f) const {
    const GeneralForceSubsystem& forces = getForceSubsystem();
    Real& fInState = Value<Real>::updDowncast
                            (forces.updDiscreteVariable(state, m_forceIx));
    fInState = f;
}

// Get the value of the generalized force to be applied.
Real Force::MobilityDiscreteForceImpl::
getMobilityForce(const State& state) const {
    const GeneralForceSubsystem& forces = getForceSubsystem();
    const Real& fInState = Value<Real>::downcast
                            (forces.getDiscreteVariable(state, m_forceIx));
    return fInState;
}

void Force::MobilityDiscreteForceImpl::
calcForce(  const State&         state,
            Vector_<SpatialVec>& /*bodyForces*/,
            Vector_<Vec3>&       /*particleForces*/,
            Vector&              mobilityForces) const
{
    const GeneralForceSubsystem& forces = getForceSubsystem();
    const Real f = Value<Real>::downcast
                            (forces.getDiscreteVariable(state, m_forceIx));
    const MobilizedBody& mobod = m_matter.getMobilizedBody(m_mobodIx);
    mobod.applyOneMobilityForce(state, m_whichU, f, mobilityForces);
}

void Force::MobilityDiscreteForceImpl::
realizeTopology(State& state) const {
    const GeneralForceSubsystem& forces = getForceSubsystem();
    m_forceIx = forces.allocateDiscreteVariable
        (state, Stage::Dynamics, new Value<Real>(m_defaultVal));
}




//------------------------------ DiscreteForces --------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::DiscreteForces,
                                        Force::DiscreteForcesImpl,
                                        Force);

Force::DiscreteForces::DiscreteForces
   (GeneralForceSubsystem& forces, const SimbodyMatterSubsystem& matter)
:   Force(new DiscreteForcesImpl(matter)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

void Force::DiscreteForces::
clearAllMobilityForces(State& state) const {
    getImpl().updAllMobilityForces(state).clear(); // set size to zero
}


void Force::DiscreteForces::
clearAllBodyForces(State& state) const {
    getImpl().updAllBodyForces(state).clear(); // set size to zero
}

void Force::DiscreteForces::
setOneMobilityForce(State& state, const MobilizedBody& mobod,
                    MobilizerUIndex whichU, Real f) const {
    Vector& mobForces = getImpl().updAllMobilityForces(state);
    if (mobForces.size() == 0) {
        mobForces.resize(state.getNU());
        mobForces.setToZero();
    }
    // Don't use "apply" here because that would add in the force.
    mobod.updOneFromUPartition(state, whichU, mobForces) = f;
}

Real Force::DiscreteForces::
getOneMobilityForce(const State& state, const MobilizedBody& mobod,
                    MobilizerUIndex whichU) const {
    const Vector& mobForces = getImpl().getAllMobilityForces(state);
    if (mobForces.size() == 0) {return 0;}

    return mobod.getOneFromUPartition(state, whichU, mobForces);
}


void Force::DiscreteForces::
setAllMobilityForces(State& state, const Vector& f) const {
    if (f.size()==0) {clearAllMobilityForces(state); return;}
    SimTK_ERRCHK2_ALWAYS(f.size() == state.getNU(),
        "Force::DiscreteForces::setAllMobilityForces()",
        "Mobility force vector f had wrong size %d; should have been %d.",
        f.size(), state.getNU());

    getImpl().updAllMobilityForces(state) = f;
}

const Vector& Force::DiscreteForces::
getAllMobilityForces(const State& state) const
{   return getImpl().getAllMobilityForces(state); }

void Force::DiscreteForces::
setOneBodyForce(State& state, const MobilizedBody& mobod,
                const SpatialVec& spatialForceInG) const {
    Vector_<SpatialVec>& bodyForces = getImpl().updAllBodyForces(state);
    if (bodyForces.size() == 0) {
        bodyForces.resize(getImpl().m_matter.getNumBodies());
        bodyForces.setToZero();
    }
    bodyForces[mobod.getMobilizedBodyIndex()] = spatialForceInG;
}


SpatialVec Force::DiscreteForces::
getOneBodyForce(const State& state, const MobilizedBody& mobod) const {
    const Vector_<SpatialVec>& bodyForces = getImpl().getAllBodyForces(state);
    if (bodyForces.size() == 0) {return SpatialVec(Vec3(0),Vec3(0));}

    return bodyForces[mobod.getMobilizedBodyIndex()];
}

const Vector_<SpatialVec>& Force::DiscreteForces::
getAllBodyForces(const State& state) const
{   return getImpl().getAllBodyForces(state); }

void Force::DiscreteForces::
setAllBodyForces(State& state, const Vector_<SpatialVec>& bodyForcesInG) const {
    if (bodyForcesInG.size()==0) {clearAllBodyForces(state); return;}
    const int numBodies = getImpl().m_matter.getNumBodies();
    SimTK_ERRCHK2_ALWAYS(bodyForcesInG.size() == numBodies,
      "Force::DiscreteForces::setAllBodyForces()",
      "Body force vector bodyForcesInG had wrong size %d; "
      "should have been %d (0th entry is for Ground).",
      bodyForcesInG.size(), numBodies);

    getImpl().updAllBodyForces(state) = bodyForcesInG;
}


void Force::DiscreteForces::
addForceToBodyPoint(State& state, const MobilizedBody& mobod,
                    const Vec3& pointInB, const Vec3& forceInG) const {
    Vector_<SpatialVec>& bodyForces = getImpl().updAllBodyForces(state);
    if (bodyForces.size() == 0) {
        bodyForces.resize(getImpl().m_matter.getNumBodies());
        bodyForces.setToZero();
    }
    mobod.applyForceToBodyPoint(state, pointInB, forceInG, bodyForces);
}

Force::DiscreteForcesImpl::DiscreteForcesImpl
   (const SimbodyMatterSubsystem& matter) : m_matter(matter) {}

const Vector& Force::DiscreteForcesImpl::
getAllMobilityForces(const State& state) const {
    const GeneralForceSubsystem& forces = getForceSubsystem();
    const Vector& fInState = Value<Vector>::downcast
                            (forces.getDiscreteVariable(state, m_mobForcesIx));
    return fInState;
}

Vector& Force::DiscreteForcesImpl::
updAllMobilityForces(State& state) const  {
    const GeneralForceSubsystem& forces = getForceSubsystem();
    Vector& fInState = Value<Vector>::updDowncast
                            (forces.updDiscreteVariable(state, m_mobForcesIx));
    return fInState;
}

const Vector_<SpatialVec>& Force::DiscreteForcesImpl::
getAllBodyForces(const State& state) const {
    const GeneralForceSubsystem& forces = getForceSubsystem();
    const Vector_<SpatialVec>& FInState = Value< Vector_<SpatialVec> >::downcast
                            (forces.getDiscreteVariable(state, m_bodyForcesIx));
    return FInState;
}
Vector_<SpatialVec>& Force::DiscreteForcesImpl::
updAllBodyForces(State& state) const {
    const GeneralForceSubsystem& forces = getForceSubsystem();
    Vector_<SpatialVec>& FInState = Value< Vector_<SpatialVec> >::updDowncast
                            (forces.updDiscreteVariable(state, m_bodyForcesIx));
    return FInState;
}

void Force::DiscreteForcesImpl::
calcForce(  const State&         state,
            Vector_<SpatialVec>& bodyForces,
            Vector_<Vec3>&       /*particleForces*/,
            Vector&              mobilityForces) const
{
    const GeneralForceSubsystem& forces = getForceSubsystem();
    const Vector& f = Value<Vector>::downcast
                            (forces.getDiscreteVariable(state, m_mobForcesIx));
    const Vector_<SpatialVec>& F = Value< Vector_<SpatialVec> >::downcast
                            (forces.getDiscreteVariable(state, m_bodyForcesIx));
    if (f.size()) mobilityForces += f;
    if (F.size()) bodyForces += F;
}

void Force::DiscreteForcesImpl::
realizeTopology(State& state) const {
    const GeneralForceSubsystem& forces = getForceSubsystem();
    m_mobForcesIx = forces.allocateDiscreteVariable
        (state, Stage::Dynamics, new Value<Vector>());
    m_bodyForcesIx = forces.allocateDiscreteVariable
        (state, Stage::Dynamics, new Value< Vector_<SpatialVec> >());
}



//------------------------------ ConstantForce ---------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::ConstantForce, Force::ConstantForceImpl, Force);

Force::ConstantForce::ConstantForce(GeneralForceSubsystem& forces, const MobilizedBody& body, const Vec3& station, const Vec3& force) :
        Force(new ConstantForceImpl(body, station, force)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::ConstantForceImpl::ConstantForceImpl(const MobilizedBody& body, const Vec3& station, const Vec3& force) :
        matter(body.getMatterSubsystem()), body(body.getMobilizedBodyIndex()), station(station), force(force) {
}

void Force::ConstantForceImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const Transform& X_GB = matter.getMobilizedBody(body).getBodyTransform(state);
    const Vec3 station_G = X_GB.R() * station;
    bodyForces[body] += SpatialVec(station_G % force, force);
}

Real Force::ConstantForceImpl::calcPotentialEnergy(const State& state) const {
    return 0;
}


//------------------------------ ConstantTorque --------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::ConstantTorque, Force::ConstantTorqueImpl, Force);

Force::ConstantTorque::ConstantTorque(GeneralForceSubsystem& forces, const MobilizedBody& body, const Vec3& torque) :
        Force(new ConstantTorqueImpl(body, torque)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::ConstantTorqueImpl::ConstantTorqueImpl(const MobilizedBody& body, const Vec3& torque) :
        matter(body.getMatterSubsystem()), body(body.getMobilizedBodyIndex()), torque(torque) {
}

void Force::ConstantTorqueImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    bodyForces[body][0] += torque;
}

Real Force::ConstantTorqueImpl::calcPotentialEnergy(const State& state) const {
    return 0;
}


//------------------------------- GlobalDamper ---------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::GlobalDamper, Force::GlobalDamperImpl, Force);

Force::GlobalDamper::GlobalDamper(GeneralForceSubsystem& forces, const SimbodyMatterSubsystem& matter,
        Real damping) : Force(new GlobalDamperImpl(matter, damping)) {
    assert(damping >= 0);
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Force::GlobalDamperImpl::GlobalDamperImpl(const SimbodyMatterSubsystem& matter, Real damping) : matter(matter), damping(damping) {
}

void Force::GlobalDamperImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    mobilityForces -= damping*matter.getU(state);
}

Real Force::GlobalDamperImpl::calcPotentialEnergy(const State& state) const {
    return 0;
}


//------------------------------ UniformGravity --------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::UniformGravity, Force::UniformGravityImpl, Force);

Force::UniformGravity::UniformGravity(GeneralForceSubsystem& forces, const SimbodyMatterSubsystem& matter,
        const Vec3& g, Real zeroHeight) : Force(new UniformGravityImpl(matter, g, zeroHeight)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

Vec3 Force::UniformGravity::getGravity() const {
    return getImpl().getGravity();
}

void Force::UniformGravity::setGravity(const Vec3& g) {
    updImpl().setGravity(g);
}

Real Force::UniformGravity::getZeroHeight() const {
    return getImpl().getZeroHeight();
}

void Force::UniformGravity::setZeroHeight(Real height) {
    updImpl().setZeroHeight(height);
}

Force::UniformGravityImpl::UniformGravityImpl(const SimbodyMatterSubsystem& matter, const Vec3& g, Real zeroHeight) : matter(matter), g(g), zeroHeight(zeroHeight) {
}

void Force::UniformGravityImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const int nBodies    = matter.getNumBodies();
    const int nParticles = matter.getNumParticles();

    if (nParticles) {
        const Vector& m = matter.getAllParticleMasses(state);
        const Vector_<Vec3>& loc_G = matter.getAllParticleLocations(state);
        for (int i=0; i < nParticles; ++i) {
            particleForces[i] += g * m[i];
        }
    }

    // no need to apply gravity to Ground!
    for (MobilizedBodyIndex i(1); i < nBodies; ++i) {
        const MassProperties& mprops = matter.getMobilizedBody(i).getBodyMassProperties(state);
        const Real&      m       = mprops.getMass();
        const Vec3&      com_B   = mprops.getMassCenter();
        const Transform& X_GB    = matter.getMobilizedBody(i).getBodyTransform(state);
        const Vec3       com_B_G = X_GB.R()*com_B;
        const Vec3       frc_G   = m*g;

        bodyForces[i] += SpatialVec(com_B_G % frc_G, frc_G);
    }
}

Real Force::UniformGravityImpl::calcPotentialEnergy(const State& state) const {
    const int nBodies    = matter.getNumBodies();
    const int nParticles = matter.getNumParticles();
    Real pe = 0.0;

    if (nParticles) {
        const Vector& m = matter.getAllParticleMasses(state);
        const Vector_<Vec3>& loc_G = matter.getAllParticleLocations(state);
        for (int i=0; i < nParticles; ++i) {
            pe -= m[i]*(~g*loc_G[i] + zeroHeight); // odd signs because height is in -g direction
        }
    }

    // no need to apply gravity to Ground!
    for (MobilizedBodyIndex i(1); i < nBodies; ++i) {
        const MassProperties& mprops = matter.getMobilizedBody(i).getBodyMassProperties(state);
        const Real&      m       = mprops.getMass();
        const Vec3&      com_B   = mprops.getMassCenter();
        const Transform& X_GB    = matter.getMobilizedBody(i).getBodyTransform(state);
        const Vec3       com_B_G = X_GB.R()*com_B;
        const Vec3       com_G   = X_GB.p() + com_B_G;

        pe -= m*(~g*com_G + zeroHeight); // odd signs because height is in -g direction
    }
    return pe;
}


//---------------------------------- Custom ------------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::Custom, Force::CustomImpl, Force);

Force::Custom::Custom(GeneralForceSubsystem& forces, Implementation* implementation) :
        Force(new CustomImpl(implementation)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

const Force::Custom::Implementation& Force::Custom::getImplementation() const {
    return getImpl().getImplementation();
}

Force::Custom::Implementation& Force::Custom::updImplementation() {
    return updImpl().updImplementation();
}


Force::CustomImpl::CustomImpl(Force::Custom::Implementation* implementation) : implementation(implementation) {
}

void Force::CustomImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    implementation->calcForce(state, bodyForces, particleForces, mobilityForces);
}

Real Force::CustomImpl::calcPotentialEnergy(const State& state) const {
    return implementation->calcPotentialEnergy(state);
}

} // namespace SimTK
