/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Ian Stavness, Michael Sherman                                     *
 * Contributors: Andreas Scholz                                               *
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
#include "simmath/internal/common.h"
#include "simmath/internal/Geodesic.h"
#include "simmath/internal/GeodesicIntegrator.h"
#include "simmath/internal/ContactGeometry.h"

#include "GeodesicEquations.h"

#include <iostream>
#include <cmath>

using namespace SimTK;

// Make sure these compile correctly.
template class GeodesicIntegrator<GeodesicOnImplicitSurface>;
template class GeodesicIntegrator<GeodesicOnParametricSurface>;

//==============================================================================
//                      GEODESIC ON IMPLICIT SURFACE
//==============================================================================

// See class header for documentation.
void GeodesicOnImplicitSurface::
calcDerivs(Real t, const Vec<N>& y, Vec<N>& ydot) const {
    const Vec3& p = getP(y);        // rename state variables
    const Vec3& v = getV(y);
    const Real& jr = getJRot(y);
    const Real& jt = getJTrans(y);

    // Evaluate the surface at p.
    const Vec3  g = geom.calcSurfaceGradient(p);
    const Mat33 H = geom.calcSurfaceHessian(p);
    Real Kg = geom.calcGaussianCurvature(g,H);

    const Real Gdotv = ~v*(H*v);
    const Real L = Gdotv/(~g*g);    // Lagrange multiplier

    // We have qdot = u; that part is easy.
    updQ(ydot) = getU(y);

    // These together are the udots.
    Vec3& a     = updV(ydot);          // d/dt v
    Real& jrdd  = updJRotDot(ydot);    // d/dt jdr
    Real& jtdd  = updJTransDot(ydot);  // d/dt jdt

    a    = -L*g;
    jrdd = -Kg*jr;
    jtdd = -Kg*jt;
}


// See class header for documentation.
void GeodesicOnImplicitSurface::
calcConstraintErrors(Real t, const Vec<N>& y, Vec<NC>& cerr) const {
    const Vec3& p = getP(y);
    const Vec3& v = getV(y);
    // This is the perr() equation that says the point must be on the surface.
    cerr[0] = geom.calcSurfaceValue(p);
    // These are the two verr() equations. The first is the derivative of
    // the above point-on-surface holonomic constraint above. The second is 
    // a nonholonomic velocity constraint restricting the velocity along 
    // the curve to be 1.
    cerr[1] = ~geom.calcSurfaceGradient(p)*v;
    cerr[2] = v.norm() - 1;
}

// Given a state y drive the infinity norm of the position and velocity 
// constraint errors to consTol or below by adjusting y.
bool GeodesicOnImplicitSurface::
projectIfNeeded(Real consTol, Real t, Vec<N>& y) const {
    const int MaxIter = 10;         // should take *far* fewer
    const Real OvershootFac = 0.1;  // try to do better than consTol
        
    const Real tryTol = consTol * OvershootFac;
    Vec3& p = updP(y); // aliases for the state variables
    Vec3& v = updV(y);

    // Fix the position constraint first. This is a Newton interation
    // that modifies only the point location to make sure it remains on
    // the surface. No position projection is done if we're already at
    // tryTol, which is a little tighter than the requested consTol.

    // NOTE: (sherm) I don't think this is exactly the right projection.
    // Here we project down the gradient, but the final result won't 
    // be exactly the nearest point on the surface if the gradient changes
    // direction on the way down. For correcting small errors this is
    // probably completely irrelevant since the starting and final gradient
    // directions will be the same.

    Real perr, ptolAchieved;
    int piters=0; 
    while (true) {
        perr = geom.calcSurfaceValue(p);
        ptolAchieved = std::abs(perr); 
        if (ptolAchieved <= tryTol || piters==MaxIter)
            break;

        ++piters;
        // We want a least squares solution dp to ~g*dp=perr which we
        // get using the pseudoinverse: dp=pinv(~g)*perr, where
        // pinv(~g) = g*inv(~g*g).
        const Vec3 g = geom.calcSurfaceGradient(p);
        const Vec3 pinvgt = g/(~g*g);
        const Vec3 dp = pinvgt*perr;

        p -= dp; // updates the state
    }


    // Now the velocities. There are two velocity constraints that have
    // to be satisfied simultaneously. They are (1) the time derivative of 
    // the perr equation which we just solved, and (2) the requirement that 
    // the velocity magnitude be 1. So verr=~[ ~g*v, |v|-1 ]. You might 
    // think these need to be solved simultaneously to find the least 
    // squares dv, but dv can be determined by two orthogonal projections.
    // The allowable velocity vectors form a unit circle whose normal is
    // in the gradient direction. The least squares dv is the shortest
    // vector from the end point of v to that cicle. To find the closest
    // point on the unit circle, first project the vector v onto the 
    // circle's plane by the shortest path (remove the normal component). 
    // Then stretch the result to unit length.
    // First we solve the linear least squares problem ~g*(v+dv0)=0 for
    // dv0, and set v0=v+dv0. Then set vfinal = v0/|v0|, giving dv=vfinal-v.

    // We're going to project velocities unconditionally because we
    // would have to evaluate the constraint anyway to see if it is
    // violated and that is most of the computation we need to fix it.

    const Vec3 g = geom.calcSurfaceGradient(p);
    const Vec3 pinvgt = g/(~g*g);
    const Real perrdot = ~g*v;

    const Vec3 dv0 = pinvgt*perrdot;
    const Vec3 v0 = v - dv0;    // fix direction
    v = v0/v0.norm();           // fix length; updates state

    const bool success = (ptolAchieved <= consTol);
    return success;
}