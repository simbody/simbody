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
    const Real cf = 1e-5;
    const Real bd = 300;
    const Real bv = 50;

    Body::Rigid body1(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Translation sphere(matter.updGround(),
        Transform(), body1, Transform());

    Body::Rigid body2(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Free halfSpace(matter.updGround(),
        Transform(), body2, Transform());

    SmoothSphereHalfSpaceForce hc_smooth(forces);

    hc_smooth.setParameters(k,dissipation,us,ud,uv,vt,cf,bd,bv);
    hc_smooth.setContactSphereBody(sphere);
    hc_smooth.setContactSphereLocationInBody(Vec3(0));
    hc_smooth.setContactSphereRadius(radius);
    Transform testFrame(Rotation(-0.5*Pi, ZAxis), Vec3(0));
    hc_smooth.setContactHalfSpaceFrame(testFrame);
    hc_smooth.setContactHalfSpaceBody(halfSpace);
    State state = system.realizeTopology();
    // Position the sphere at a variety of positions and see if the normal
    // force and potential energy are correct.
    for (Real height = radius+0.2; height > 0; height -= 0.1) {
        sphere.setQToFitTranslation(state, Vec3(0, height, 0));
        system.realize(state, Stage::Dynamics);
        const Real depth = radius-height;
        Real f = (4./3.)*stiffness*std::pow(std::sqrt(depth*depth+cf),3./2.)
            *std::sqrt(radius*stiffness);
        Real f_smooth = f*(1./2.+(1./2.)*std::tanh(bd*depth));
        assertEqual(system.getRigidBodyForces(state, Stage::Dynamics)
            [sphere.getMobilizedBodyIndex()][1], gravity+Vec3(0, f_smooth, 0));
        assertEqual(hc_smooth.calcPotentialEnergyContribution(state),
            (2./5.)*f_smooth*depth);
    }

    // Now do it with a vertical velocity and see if the dissipation force is
    // correct.
    for (Real height = radius+0.2; height > 0; height -= 0.1) {
        sphere.setQToFitTranslation(state, Vec3(0, height, 0));
        const Real depth = radius-height;
        Real fh = (4./3.)*stiffness*std::pow(std::sqrt(depth*depth+cf),3./2.)
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
    // correct.
    Vector_<SpatialVec> expectedForce(matter.getNumBodies());
    for (Real height = radius+0.2; height > 0; height -= 0.1) {
        sphere.setQToFitTranslation(state, Vec3(0, height, 0));
        const Real depth = radius-height;
        Real fh = (4./3.)*stiffness*std::pow(std::sqrt(depth*depth+cf),3./2.)
            *std::sqrt(radius*stiffness);
        Real fh_smooth = fh*(1./2.+(1./2.)*std::tanh(bd*depth));

        for (Real v = -1.0; v <= 1.0; v += 0.1) {
            sphere.setUToFitLinearVelocity(state, Vec3(v, 0, 0));
            system.realize(state, Stage::Dynamics);
            Vec3 vec3v(v,0,0);
            UnitVec3 normal = (halfSpace.getBodyRotation(state)*testFrame.x());
            Real vnormal = dot(vec3v, normal);
            Vec3 vtangent = vec3v - vnormal*normal;
            Real aux = vtangent.normSqr() + cf;
            Real vslip = pow(aux,1./2.);
            Real vrel = vslip / vt;
            Real ff_smooth_scalar = fh_smooth*(std::min(vrel,Real(1))*
                (ud+2*(us-ud)/(1+vrel*vrel))+uv*vslip);
            Vec3 ff_smooth = ff_smooth_scalar*(vtangent) / vslip;
            const Vec3 totalForceOnSphere =
                gravity - ff_smooth - fh_smooth*normal;
            expectedForce = SpatialVec(Vec3(0), Vec3(0));
            Vec3 contactPointInSphere = sphere.findStationAtGroundPoint(state,
                Vec3(0, -stiffness*depth/(stiffness+stiffness), 0));
            sphere.applyForceToBodyPoint(state, contactPointInSphere,
                totalForceOnSphere, expectedForce);
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
