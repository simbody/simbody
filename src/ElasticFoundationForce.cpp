/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/Contact.h"
#include "simbody/internal/ContactGeometry.h"
#include "simbody/internal/GeneralContactSubsystem.h"
#include "simbody/internal/MobilizedBody.h"
#include "ElasticFoundationForceImpl.h"
#include <set>

using std::set;
using std::vector;

namespace SimTK {

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(ElasticFoundationForce, ElasticFoundationForceImpl, Force);

ElasticFoundationForce::ElasticFoundationForce(GeneralForceSubsystem& forces, GeneralContactSubsystem& contacts, ContactSetIndex set, int meshIndex) :
        Force(new ElasticFoundationForceImpl(contacts, set, meshIndex)) {
    updImpl().setForceIndex(forces.adoptForce(*this));
    const ContactGeometry::TriangleMesh& mesh = static_cast<const ContactGeometry::TriangleMesh&>(contacts.getBodyGeometry(set, meshIndex));
    updImpl().springPosition.resize(mesh.getNumFaces());
    updImpl().springNormal.resize(mesh.getNumFaces());
    updImpl().springArea.resize(mesh.getNumFaces());
    Vec2 uv(1.0/3.0, 1.0/3.0);
    for (int i = 0; i < getImpl().springPosition.size(); i++) {
        updImpl().springPosition[i] = (mesh.getVertexPosition(mesh.getFaceVertex(i, 0))+mesh.getVertexPosition(mesh.getFaceVertex(i, 1))+mesh.getVertexPosition(mesh.getFaceVertex(i, 2)))/3.0;
        updImpl().springNormal[i] = -mesh.findNormalAtPoint(i, uv);
        updImpl().springArea[i] = mesh.getFaceArea(i);
    }
}

Real ElasticFoundationForce::getYoungsModulus() const {
    return getImpl().young;
}

void ElasticFoundationForce::setYoungsModulus(Real v) {
    updImpl().young = v;
}

Real ElasticFoundationForce::getPoissonsRatio() const {
    return getImpl().poisson;
}

void ElasticFoundationForce::setPoissonsRatio(Real v) {
    SimTK_APIARGCHECK1(v > -1 && v < 0.5, "ElasticFoundationForce", "setPoissonsRatio", "Illegal value: %f", v);
    updImpl().poisson = v;
}

Real ElasticFoundationForce::getThickness() const {
    return getImpl().thickness;
}

void ElasticFoundationForce::setThickness(Real v) {
    updImpl().thickness = v;
}

Real ElasticFoundationForce::getTransitionVelocity() const {
    return getImpl().transitionVelocity;
}

void ElasticFoundationForce::setTransitionVelocity(Real v) {
    updImpl().transitionVelocity = v;
}

ElasticFoundationForceImpl::ElasticFoundationForceImpl(GeneralContactSubsystem& subsystem, ContactSetIndex set, int meshIndex) : 
        subsystem(subsystem), set(set), meshIndex(meshIndex), young(0), poisson(0), thickness(1), transitionVelocity(0.001), energyCacheIndex(-1) {
}

void ElasticFoundationForceImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const vector<Contact>& contacts = subsystem.getContacts(state, set);
    Real& pe = Value<Real>::downcast(state.updCacheEntry(subsystem.getMySubsystemIndex(), energyCacheIndex)).upd();
    pe = 0.0;
    for (int i = 0; i < contacts.size(); i++) {
        if (contacts[i].getFirstBody() != meshIndex && contacts[i].getSecondBody() != meshIndex)
            continue;
        const TriangleMeshContact& contact = static_cast<const TriangleMeshContact&>(contacts[i]);
        bool meshIsFirst = (contact.getFirstBody() == meshIndex);
        const std::set<int>& insideFaces = (meshIsFirst ? contact.getFirstBodyFaces() : contact.getSecondBodyFaces());
        int otherBodyIndex = meshIsFirst ? contact.getSecondBody() : contact.getFirstBody();
        const ContactGeometry& otherObject = subsystem.getBodyGeometry(set, otherBodyIndex);
        const MobilizedBody& body1 = subsystem.getBody(set, meshIndex);
        const MobilizedBody& body2 = subsystem.getBody(set, otherBodyIndex);
        const Transform t1g = body1.getBodyTransform(state)*subsystem.getBodyTransform(set, meshIndex); // mesh to ground
        const Transform t2g = body2.getBodyTransform(state)*subsystem.getBodyTransform(set, otherBodyIndex); // other object to ground
        const Transform t12 = ~t2g*t1g; // mesh to other object

        // Loop over all the springs, and evaluate the force from each one.
        
        Real k = young*(1-poisson)/((1+poisson)*(1-2*poisson)*thickness);
        for (std::set<int>::const_iterator iter = insideFaces.begin(); iter != insideFaces.end(); ++iter) {
            UnitVec3 normal;
            int face = *iter;
            
            bool inside;
            Vec3 nearestPoint = otherObject.findNearestPoint(t12*springPosition[face], inside, normal);
            if (!inside)
                continue;
            nearestPoint = t2g*nearestPoint;
            const Vec3 springPosInGround = t1g*springPosition[face];
            const Vec3 displacement = nearestPoint-springPosInGround;
            const Real distance = displacement.norm();
            Vec3 force = k*springArea[face]*displacement;
            const Vec3 station1 = body1.findStationAtGroundPoint(state, nearestPoint);
            const Vec3 station2 = body2.findStationAtGroundPoint(state, nearestPoint);
            body1.applyForceToBodyPoint(state, station1, force, bodyForces);
            body2.applyForceToBodyPoint(state, station2, -force, bodyForces);
            pe += 0.5*k*springArea[face]*distance*distance;
        }
        
//        const Parameters& param1 = getParameters(contacts[i].getFirstBody());
//        const Parameters& param2 = getParameters(contacts[i].getSecondBody());
//        
//        // Adjust the contact location based on the relative stiffness of the two materials.
//        
//        const Real s1 = param2.stiffness/(param1.stiffness+param2.stiffness);
//        const Real s2 = 1-s1;
//        const Real depth = contacts[i].getDepth();
//        const Vec3& normal = contacts[i].getNormal();
//        const Vec3 location = contacts[i].getLocation()+(depth*(0.5-s1))*normal;
//        
//        // Calculate the Hertz force.
//
//        const Real k = param1.stiffness*s1;
//        const Real c = param1.dissipation*s1 + param2.dissipation*s2;
//        const Real radius = contacts[i].getRadius();
//        const Real curvature = radius*radius/depth;
//        const Real fH = (4.0/3.0)*k*depth*std::sqrt(curvature*k*depth);
//        
//        // Calculate the relative velocity of the two bodies at the contact point.
//        
//        const MobilizedBody& body1 = subsystem.getBody(set, contacts[i].getFirstBody());
//        const MobilizedBody& body2 = subsystem.getBody(set, contacts[i].getSecondBody());
//        const Vec3 station1 = body1.findStationAtGroundPoint(state, location);
//        const Vec3 station2 = body2.findStationAtGroundPoint(state, location);
//        const Vec3 v1 = body1.findStationVelocityInGround(state, station1);
//        const Vec3 v2 = body2.findStationVelocityInGround(state, station2);
//        const Vec3 v = v1-v2;
//        const Real vnormal = dot(v, normal);
//        const Vec3 vtangent = v-vnormal*normal;
//        
//        // Calculate the Hunt-Crossley force.
//        
//        const Real f = fH*(1+1.5*c*vnormal);
//        Vec3 force = (f > 0 ? f*normal : Vec3(0));
//        
//        // Calculate the friction force.
//        
//        const Real vslip = vtangent.norm();
//        if (vslip != 0) {
//            const bool hasStatic = (param1.staticFriction != 0 || param2.staticFriction != 0);
//            const bool hasDynamic= (param1.dynamicFriction != 0 || param2.dynamicFriction != 0);
//            const bool hasViscous = (param1.viscousFriction != 0 || param2.viscousFriction != 0);
//            const Real us = hasStatic ? 2*param1.staticFriction*param2.staticFriction/(param1.staticFriction+param2.staticFriction) : 0;
//            const Real ud = hasDynamic ? 2*param1.dynamicFriction*param2.dynamicFriction/(param1.dynamicFriction+param2.dynamicFriction) : 0;
//            const Real uv = hasViscous ? 2*param1.viscousFriction*param2.viscousFriction/(param1.viscousFriction+param2.viscousFriction) : 0;
//            const Real vrel = vslip/getTransitionVelocity();
//            const Real ffriction = f*(std::min(vrel, 1.0)*(ud+2*(us-ud)/(1+vrel*vrel))+uv*vslip);
//            force += ffriction*vtangent/vslip;
//        }
//        
//        // Apply the force to the bodies.
//        
//        body1.applyForceToBodyPoint(state, station1, -force, bodyForces);
//        body2.applyForceToBodyPoint(state, station2, force, bodyForces);
    }
}

Real ElasticFoundationForceImpl::calcPotentialEnergy(const State& state) const {
    return Value<Real>::downcast(state.getCacheEntry(subsystem.getMySubsystemIndex(), energyCacheIndex)).get();
}

void ElasticFoundationForceImpl::realizeTopology(State& state) const {
        energyCacheIndex = state.allocateCacheEntry(subsystem.getMySubsystemIndex(), Stage::Dynamics, new Value<Real>());
}


} // namespace SimTK

