/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-19 Stanford University and the Authors.        *
 * Authors: Antoine Falisse, Gil Serrancoli                                   *
 * Contributors: Peter Eastman                                                *
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

#include "SimTKsimbody.h"

using namespace SimTK;
using namespace std;

const Real TOL = 1e-10;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

template <class T>
void assertEqual(T val1, T val2) {
    ASSERT(abs(val1-val2) < TOL);
}

template <int N>
void assertEqual(Vec<N> val1, Vec<N> val2) {
    for (int i = 0; i < N; ++i)
        ASSERT(abs(val1[i]-val2[i]) < TOL);
}

void testForces() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    const Vec3 gravity = Vec3(0, -9.8, 0);
    Force::UniformGravity(forces, matter, gravity, 0);
    const Real radius = 0.8;
    const Real k = 1500.0;
    const Real stiffness = 0.5*std::pow(k, 2.0/3.0);
    const Real dissipation = 0.5;
    const Real us = 1.0;
    const Real ud = 0.5;
    const Real uv = 0.1;
    const Real vt = 0.001;
    double eps = 1e-5;
    double bd = 300;
    double bv = 50;
    Random::Uniform random(0.0, 1.0);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Translation sphere(matter.updGround(),
        Transform(), body, Transform());
    SmoothSphereHalfplaneContact hc_smooth(forces);
    hc_smooth.setParameters(k,dissipation,us,ud,uv,vt);
    Vec3 normal(0,1,0);
    hc_smooth.setContactPlane(normal,.0);
    hc_smooth.setContactSphere(sphere);
    hc_smooth.setLocationContactSphere(Vec3(0));
    hc_smooth.setRadiusContactSphere(radius);
    State state = system.realizeTopology();
    // Position the sphere at a variety of positions and see if the normal
    // force and potential energy are correct (with horizontal ground plane)
    for (Real height = radius+0.2; height > 0; height -= 0.1) {
        sphere.setQToFitTranslation(state, Vec3(0, height, 0));
        system.realize(state, Stage::Dynamics);
        const Real depth = radius-height;
        Real f = (4./3.)*stiffness*std::pow(std::sqrt(depth*depth+eps),3./2.)
            *std::sqrt(radius*stiffness);
        Real f_smooth = f*(1./2.+(1./2.)*std::tanh(bd*depth));
        assertEqual(system.getRigidBodyForces(state, Stage::Dynamics)
            [sphere.getMobilizedBodyIndex()][1], gravity+Vec3(0, f_smooth, 0));
        assertEqual(hc_smooth.calcPotentialEnergyContribution(state),
            (2./ 5.)*f*depth);
    }

    // Now do it with a vertical velocity and see if the dissipation force is
    // correct (with horizontal ground plane)
    for (Real height = radius+0.2; height > 0; height -= 0.1) {
        sphere.setQToFitTranslation(state, Vec3(0, height, 0));
        const Real depth = radius-height;
        Real fh = (4./3.)*stiffness*std::pow(std::sqrt(depth*depth+eps),3./2.)
            *std::sqrt(radius*stiffness);
        Real fh_smooth = fh*(1./2.+(1./2.)*std::tanh(bd*depth));

        for (Real v = -1.0; v <= 1.0; v += 0.1) {
            sphere.setUToFitLinearVelocity(state, Vec3(0, -v, 0));
            system.realize(state, Stage::Dynamics);
            Real f = fh_smooth*(1.+(3./2.)*dissipation*v);
            Real f_smooth = f*(1./2.+(1./2.)
                *std::tanh(bv*(v+(2./(3.*dissipation)))));
            assertEqual(system.getRigidBodyForces(state, Stage::Dynamics)
                [sphere.getMobilizedBodyIndex()][1],
                gravity+Vec3(0, f_smooth, 0));
        }
    }

    // Now do it with a horizontal velocity and see if the friction force is
    // correct (with horizontal ground plane)
    Vector_<SpatialVec> expectedForce(matter.getNumBodies());
    for (Real height = radius+0.2; height > 0; height -= 0.1) {
        sphere.setQToFitTranslation(state, Vec3(0, height, 0));
        const Real depth = radius-height;
        Real fh = (4./3.)*stiffness*std::pow(std::sqrt(depth*depth+eps),3./2.)
            *std::sqrt(radius*stiffness);
        Real fh_smooth = fh*(1./2.+(1./2.)*std::tanh(bd*depth));

        for (Real v = -1.0; v <= 1.0; v += 0.1) {
            sphere.setUToFitLinearVelocity(state, Vec3(v, 0, 0));
            system.realize(state, Stage::Dynamics);
            Vec3 vec3v(v,0,0);
            Real vnormal = dot(vec3v, normal);
            Vec3 vtangent = vec3v - vnormal*normal;
            Real aux = pow(vtangent[0],2) + pow(vtangent[1],2) +
                pow(vtangent[2],2) + eps;
            Real vslip = pow(aux,1./2.);
            Real vrel = vslip / vt;
            Real ff_smooth_scalar = fh_smooth*(std::min(vrel,Real(1))*
                (ud+2*(us-ud)/(1+vrel*vrel))+uv*vslip);
            Vec3 ff_smooth = ff_smooth_scalar*(-vtangent) / vslip;
            const Vec3 totalForce = gravity + ff_smooth + fh_smooth*normal;
            expectedForce = SpatialVec(Vec3(0), Vec3(0));
            Vec3 contactPointInSphere = sphere.findStationAtGroundPoint(state,
                Vec3(0, -stiffness*depth/(stiffness+stiffness), 0));
            sphere.applyForceToBodyPoint(state, contactPointInSphere,
                totalForce, expectedForce);
            SpatialVec actualForce = system.getRigidBodyForces(state,
                Stage::Dynamics)[sphere.getMobilizedBodyIndex()];
            assertEqual(actualForce[0],
                expectedForce[sphere.getMobilizedBodyIndex()][0]);
            assertEqual(actualForce[1],
                expectedForce[sphere.getMobilizedBodyIndex()][1]);
        }
    }

    // Now do it with different ground angles and offsets
    Vector_<SpatialVec> expectedForce_contact(matter.getNumBodies());
    Vector_<SpatialVec> expectedForce_gravity(matter.getNumBodies());
    for (int i = 0; i < 20; i++) {
        for (Real offset = -0.1; offset <= 0.1; offset += 0.01) {
            Vec3 normal_rand((Real)(rand() % 100 - 50)  / RAND_MAX,
                (Real)(rand() % 100 - 50) / RAND_MAX,
                (Real)(rand() % 100 - 50) / RAND_MAX);
            normal_rand = normal_rand / normal_rand.norm();
            hc_smooth.setContactPlane(normal_rand, offset);
            Vec3 loc_inG(0.23, 0.10, 0.42); // random location of the sphere in
            // the space, in the ground frame
            sphere.setQToFitTranslation(state, loc_inG);
            Vec3 vel_inG(0.34, -0.65, -0.48); // random velocity of the sphere
            // in the space, in the ground frame
            sphere.setUToFitLinearVelocity(state, vel_inG);
            system.realize(state, Stage::Dynamics);

            Vec3 contactPoint = loc_inG - normal_rand*radius;
            const Real depth = -(dot(contactPoint,normal_rand)-offset);
            Real vnormal = dot(vel_inG, normal_rand);
            Vec3 vtangent = vel_inG - vnormal*normal_rand;
            Real IndentationVel = -vnormal;

            Real fh = (4./3.)*stiffness*
                std::pow(std::sqrt(depth*depth + eps), 3./2.)*
                std::sqrt(radius*stiffness);
            Real fhd = fh*(1.+(3./2.)*dissipation*IndentationVel);
            Real fh_smooth = fhd*(1./2.+(1./2.)*std::tanh(bd*depth))*
                (1./2.+(1./2.)*std::tanh(bv*(IndentationVel +
                (2./(3.*dissipation)))));;

            Real aux = pow(vtangent[0],2) + pow(vtangent[1],2) +
                pow(vtangent[2],2) + eps;
            Real vslip = pow(aux,1./2.);
            Real vrel = vslip/vt;
            Real ff_smooth_scalar = fh_smooth*(std::min(vrel, Real(1))*
                (ud+2*(us-ud)/(1+vrel*vrel))+uv*vslip);
            Vec3 ff_smooth = ff_smooth_scalar*(-vtangent) / vslip;
            const Vec3 totalForce = ff_smooth + fh_smooth*normal_rand;
            expectedForce_contact = SpatialVec(Vec3(0), Vec3(0));
            expectedForce_gravity = SpatialVec(Vec3(0), Vec3(0));
            expectedForce = SpatialVec(Vec3(0), Vec3(0));

            Vec3 contactPointInSphere = -normal_rand*radius +
                normal_rand*depth*(stiffness / (stiffness + stiffness));
            sphere.applyForceToBodyPoint(state, contactPointInSphere,
                totalForce, expectedForce_contact);
            sphere.applyForceToBodyPoint(state, Vec3(0,0,0),
                gravity, expectedForce_gravity);
            expectedForce = expectedForce_contact + expectedForce_gravity;

            SpatialVec actualForce = system.getRigidBodyForces(state,
                Stage::Dynamics)[sphere.getMobilizedBodyIndex()];
            assertEqual(actualForce[0],
                expectedForce[sphere.getMobilizedBodyIndex()][0]);
            assertEqual(actualForce[1],
                expectedForce[sphere.getMobilizedBodyIndex()][1]);
        }
    }
}

int main() {
    try {
        testForces();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
