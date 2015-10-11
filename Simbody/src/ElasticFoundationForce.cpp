/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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
#include "simbody/internal/GeneralContactSubsystem.h"
#include "simbody/internal/MobilizedBody.h"
#include "ElasticFoundationForceImpl.h"
#include <map>
#include <set>

namespace SimTK {

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(ElasticFoundationForce, ElasticFoundationForceImpl, Force);

ElasticFoundationForce::ElasticFoundationForce(GeneralForceSubsystem& forces, GeneralContactSubsystem& contacts, ContactSetIndex set) :
        Force(new ElasticFoundationForceImpl(contacts, set)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

void ElasticFoundationForce::setBodyParameters
   (ContactSurfaceIndex bodyIndex, Real stiffness, Real dissipation,
    Real staticFriction, Real dynamicFriction, Real viscousFriction) {
    updImpl().setBodyParameters(bodyIndex, stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction);
}

Real ElasticFoundationForce::getTransitionVelocity() const {
    return getImpl().transitionVelocity;
}

void ElasticFoundationForce::setTransitionVelocity(Real v) {
    updImpl().transitionVelocity = v;
}

ElasticFoundationForceImpl::ElasticFoundationForceImpl
   (GeneralContactSubsystem& subsystem, ContactSetIndex set) :
        subsystem(subsystem), set(set), transitionVelocity(Real(0.01)) {
}

void ElasticFoundationForceImpl::setBodyParameters
   (ContactSurfaceIndex bodyIndex, Real stiffness, Real dissipation,
    Real staticFriction, Real dynamicFriction, Real viscousFriction) {
    SimTK_APIARGCHECK1(bodyIndex >= 0 && bodyIndex < subsystem.getNumBodies(set), "ElasticFoundationForceImpl", "setBodyParameters",
            "Illegal body index: %d", (int)bodyIndex);
    SimTK_APIARGCHECK1(subsystem.getBodyGeometry(set, bodyIndex).getTypeId()
                        == ContactGeometry::TriangleMesh::classTypeId(),
        "ElasticFoundationForceImpl", "setBodyParameters",
        "Body %d is not a triangle mesh", (int)bodyIndex);
    parameters[bodyIndex] =
        Parameters(stiffness, dissipation, staticFriction, dynamicFriction,
                   viscousFriction);
    const ContactGeometry::TriangleMesh& mesh =
        ContactGeometry::TriangleMesh::getAs
                (subsystem.getBodyGeometry(set, bodyIndex));
    Parameters& param = parameters[bodyIndex];
    param.springPosition.resize(mesh.getNumFaces());
    param.springNormal.resize(mesh.getNumFaces());
    param.springArea.resize(mesh.getNumFaces());
    Vec2 uv(Real(1./3.), Real(1./3.));
    for (int i = 0; i < (int) param.springPosition.size(); i++) {
        param.springPosition[i] =
           (mesh.getVertexPosition(mesh.getFaceVertex(i, 0))
            +mesh.getVertexPosition(mesh.getFaceVertex(i, 1))
            +mesh.getVertexPosition(mesh.getFaceVertex(i, 2)))/3;
        param.springNormal[i] = -mesh.findNormalAtPoint(i, uv);
        param.springArea[i] = mesh.getFaceArea(i);
    }
    subsystem.invalidateSubsystemTopologyCache();
}

void ElasticFoundationForceImpl::calcForce
   (const State& state, Vector_<SpatialVec>& bodyForces,
    Vector_<Vec3>& particleForces, Vector& mobilityForces) const
{
    const Array_<Contact>& contacts = subsystem.getContacts(state, set);
    Real& pe = Value<Real>::downcast
                (subsystem.updCacheEntry(state, energyCacheIndex));
    pe = 0.0;
    for (int i = 0; i < (int) contacts.size(); i++) {
        std::map<ContactSurfaceIndex, Parameters>::const_iterator iter1 =
            parameters.find(contacts[i].getSurface1());
        std::map<ContactSurfaceIndex, Parameters>::const_iterator iter2 =
            parameters.find(contacts[i].getSurface2());

        // If there are two meshes, scale each one's contributions by 50%.
        Real areaScale = (iter1==parameters.end() || iter2==parameters.end())
                         ? Real(1) : Real(0.5);

        if (iter1 != parameters.end()) {
            const TriangleMeshContact& contact =
                static_cast<const TriangleMeshContact&>(contacts[i]);
            processContact(state, contact.getSurface1(),
                contact.getSurface2(), iter1->second,
                contact.getSurface1Faces(), areaScale, bodyForces, pe);
        }

        if (iter2 != parameters.end()) {
            const TriangleMeshContact& contact =
                static_cast<const TriangleMeshContact&>(contacts[i]);
            processContact(state, contact.getSurface2(),
                contact.getSurface1(), iter2->second,
                contact.getSurface2Faces(), areaScale, bodyForces, pe);
        }
    }
}

void ElasticFoundationForceImpl::processContact
   (const State& state,
    ContactSurfaceIndex meshIndex, ContactSurfaceIndex otherBodyIndex,
    const Parameters& param, const std::set<int>& insideFaces,
    Real areaScale, Vector_<SpatialVec>& bodyForces, Real& pe) const
{
    const ContactGeometry& otherObject = subsystem.getBodyGeometry(set, otherBodyIndex);
    const MobilizedBody& body1 = subsystem.getBody(set, meshIndex);
    const MobilizedBody& body2 = subsystem.getBody(set, otherBodyIndex);
    const Transform t1g = body1.getBodyTransform(state)*subsystem.getBodyTransform(set, meshIndex); // mesh to ground
    const Transform t2g = body2.getBodyTransform(state)*subsystem.getBodyTransform(set, otherBodyIndex); // other object to ground
    const Transform t12 = ~t2g*t1g; // mesh to other object

    // Loop over all the springs, and evaluate the force from each one.

    for (std::set<int>::const_iterator iter = insideFaces.begin();
                                       iter != insideFaces.end(); ++iter) {
        int face = *iter;
        UnitVec3 normal;
        bool inside;
        Vec3 nearestPoint = otherObject.findNearestPoint(t12*param.springPosition[face], inside, normal);
        if (!inside)
            continue;

        // Find how much the spring is displaced.

        nearestPoint = t2g*nearestPoint;
        const Vec3 springPosInGround = t1g*param.springPosition[face];
        const Vec3 displacement = nearestPoint-springPosInGround;
        const Real distance = displacement.norm();
        if (distance == 0.0)
            continue;
        const Vec3 forceDir = displacement/distance;

        // Calculate the relative velocity of the two bodies at the contact point.

        const Vec3 station1 = body1.findStationAtGroundPoint(state, nearestPoint);
        const Vec3 station2 = body2.findStationAtGroundPoint(state, nearestPoint);
        const Vec3 v1 = body1.findStationVelocityInGround(state, station1);
        const Vec3 v2 = body2.findStationVelocityInGround(state, station2);
        const Vec3 v = v2-v1;
        const Real vnormal = dot(v, forceDir);
        const Vec3 vtangent = v-vnormal*forceDir;

        // Calculate the damping force.

        const Real area = areaScale * param.springArea[face];
        const Real f = param.stiffness*area*distance*(1+param.dissipation*vnormal);
        Vec3 force = (f > 0 ? f*forceDir : Vec3(0));

        // Calculate the friction force.

        const Real vslip = vtangent.norm();
        if (f > 0 && vslip != 0) {
            const Real vrel = vslip/transitionVelocity;
            const Real ffriction =
                f*(std::min(vrel, Real(1))
                 *(param.dynamicFriction+2*(param.staticFriction-param.dynamicFriction)
                 /(1+vrel*vrel))+param.viscousFriction*vslip);
            force += ffriction*vtangent/vslip;
        }

        body1.applyForceToBodyPoint(state, station1, force, bodyForces);
        body2.applyForceToBodyPoint(state, station2, -force, bodyForces);
        pe += param.stiffness*area*displacement.normSqr()/2;
    }
}

Real ElasticFoundationForceImpl::calcPotentialEnergy(const State& state) const {
    return Value<Real>::downcast
            (subsystem.getCacheEntry(state, energyCacheIndex));
}

void ElasticFoundationForceImpl::realizeTopology(State& state) const {
    energyCacheIndex = subsystem.allocateCacheEntry
                        (state, Stage::Dynamics, new Value<Real>());
}


} // namespace SimTK

