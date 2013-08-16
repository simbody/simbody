#ifndef SimTK_SIMMATH_GEODESIC_EQUATIONS_H_
#define SimTK_SIMMATH_GEODESIC_EQUATIONS_H_

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

// This is an internal header file; not part of the API. These are the 
// classes representing the differential equations that must be solved to
// calculate geodesics over general smooth surfaces.

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geodesic.h"
#include "simmath/internal/GeodesicIntegrator.h"
#include "simmath/internal/ContactGeometry.h"

#include "ContactGeometryImpl.h"

#include <iostream>
#include <cmath>

namespace SimTK {

//==============================================================================
//                      GEODESIC ON IMPLICIT SURFACE
//==============================================================================
/* This class satisfies the requirements for an Equations template argument
to the GeodesicIntegrator class. These are the differential-algebraic
equations that need to be solved to generate an arc length-parameterized
geodesic over an arbitrary implicit surface, starting at a given point and 
tangent direction. We also calculate the geodesic's Jacobi field[2,3] which 
provides scalars relating changes in direction or position at the start to 
their effects on the location of the end point.

Geodesic equations [1]
------------------
This system can be viewed as a 3d particle mass constrained to move along a 
surface with a unit velocity and no applied force except the constraint 
reaction force normal to the surface. With no applied force the particle traces
a geodesic along the surface.

The DAE for a generic multibody system is:
      (1) qdot               = Nu
      (2) M udot + ~G lambda = f
      (3) G udot             = b
      (4) perr(t,q)          = 0
      (5) verr(t,q,u)        = 0

Let   p be the 3d coordinates of the particle (=q)
      v = p'  = dp/ds the particle's velocity (=u)
      a = p'' = dv/ds the particle's acceleration (=udot)
      S(p) the implicit surface function
      g(p) = ~DS/Dp the gradient of the implicit surface function (3x1)
      H(p) = Dg/Dp the Hessian of the implicit surface function (3x3 symmetric)
      s is arc length along the geodesic curve

Notation
      D  = partial derivative
      ~A = transpose of matrix A
      A\ = inverse of matrix A

We express the implicit surface constraint and its derivatives with respect
to the independent variable s:
      (6)         S(p) = 0        (perr)
      (7)         ~g v = 0        (verr)
      (8) ~g a = -~g'v = - ~v H v (aerr)

There is also a non-holonomic (velocity) constraint required to maintain the
tangential velocity's magnitude to 1, ensuring that we get an arc length
parameterization of the geodesic:
      (9)          |v| = 1        (verr)

So in eqns. (2,3) the matrix M=I (unit mass particle), G = ~g and b = -vHv, and 
the equations of motion are:
              [  I  g ] [ a ]   [    0   ]
      (10)    [ ~g  0 ] [ L ] = [ - ~vHv ]
where L (the Lagrange multiplier) is proportional to the constraint force.
There is no multiplier associated with equation (9) since there are no
tangential forces.

Solving for L and then a gives:
       L = (~g g)\ ~vHv
       a = - g L


Jacobi field equations [2,3]
----------------------
There are two second order differential 
equations to solve for the rotational and translational Jacobi field 
amplitudes:
      (1) jr'' + Kg*jr = 0,   with jr(0)=0, jr'(0)=1
      (2) jt'' + Kg*jt = 0,   with jt(0)=1, jt'(0)=0
where jr=jr(s), jt=jt(s). Kg=Kg(s) is the Gaussian curvature of the surface 
evaluated at s along the curve.


Constraint projection
---------------------
Although the equations of motion satisfy the constraints at the acceleration
level, integration error will cause the solution to drift away from the 
position and velocity constraint manifolds. Unlike ref. [1] which uses
Baumgarte stabilization, we use the more tractable coordinate projection
method[4] to prevent this drift. This requires that at each step we perform
a least-squares projection of the state variables back onto the constraint
manifolds. Coordinate projection guarantees that the constraints are always
satisfied, and the solution to the differential equations is improved.

See the implementation of the projectIfNeeded() method for details.


References
----------
[1] De Sapio, V., Khatib, O., Delp, S. Least action principles and their 
application to constrained and task-level problems in robotics and 
biomechanics. Multibody System Dynamics 19(3):303 (2008), section 3.1.
[2] Do Carmo, M.P. Differential Geometry of Curves and Surfaces,
Chapter 5-5 Jacobi Fields and Conjugate Points. Prentice Hall (1976).
[3] For real understanding, see Andreas Scholz' master's thesis (2012).
[4] Eich, E. Convergence results for a coordinate projection method applied 
to mechanical systems with algebraic constraints. Siam Journal on Numerical 
Analysis 30(5):1467 (1993).
*/
class SimTK_SIMMATH_EXPORT GeodesicOnImplicitSurface {
public:
    // state y = q | u = px py pz jr jt | vx vy vz jrd jtd
    // NQ, NC are required by the GeodesicIntegrator
    enum { D = 3,       // 3 coordinates for a point on an implicit surface.
           NJ = 2,      // Number of Jacobi field equations.
           NQ = D+NJ,   // Number of 2nd order equations.
           N  = 2*NQ,   // Number of differential equations.
           NC = 3 };    // 3 constraints: point on surface, point velocity
                        // along surface, unit velocity
    
    GeodesicOnImplicitSurface(const ContactGeometryImpl& geom)
    :   geom(geom) {}

    // This method is required by the GeodesicIntegrator.
    // Calculate state derivatives ydot(t,y).
    // See above for documentation.
    void calcDerivs(Real t, const Vec<N>& y, Vec<N>& ydot) const;

    // This method is required by the GeodesicIntegrator.
    // Evaluate the position and velocity constraint errors cerr(t,y).
    // See above for documentation.
    void calcConstraintErrors(Real t, const Vec<N>& y, Vec<NC>& cerr) const;

    // This method is required by the GeodesicIntegrator.
    // Given a state y drive the infinity norm of the position and velocity 
    // constraint errors to consTol or below by adjusting y.
    bool projectIfNeeded(Real consTol, Real t, Vec<N>& y) const;

    // Utility routine for filling in the initial state given a starting
    // point and direction. Note that there is no guarantee that the resulting
    // state satisfies the constraints.
    static Vec<N> getInitialState(const Vec3& P, const UnitVec3& tP) {
        Vec<N> y;
        updP(y)      = P;   updV(y)         = tP.asVec3();
        updJRot(y)   = 0;   updJRotDot(y)   = 1;
        updJTrans(y) = 1;   updJTransDot(y) = 0;
        return y;
    }

    // These define the meanings of the state variables & derivatives.
    static const Vec<NQ>& getQ(const Vec<N>& y) {return Vec<NQ>::getAs(&y[0]);}
    static Vec<NQ>& updQ(Vec<N>& y) {return Vec<NQ>::updAs(&y[0]);}

    static const Vec<NQ>& getU(const Vec<N>& y) {return Vec<NQ>::getAs(&y[NQ]);}
    static Vec<NQ>& updU(Vec<N>& y) {return Vec<NQ>::updAs(&y[NQ]);}


    // Extract the point location from a full state y.
    static const Vec3& getP(const Vec<N>& y) {return Vec3::getAs(&y[0]);}
    static Vec3& updP(Vec<N>& y) {return Vec3::updAs(&y[0]);}
    // Extract the point velocity from a full state y.
    static const Vec3& getV(const Vec<N>& y) {return Vec3::getAs(&y[NQ]);}
    static Vec3& updV(Vec<N>& y) {return Vec3::updAs(&y[NQ]);}

    // Extract the value of the rotational Jacobi field from a state y.
    static const Real& getJRot(const Vec<N>& y) {return y[D];}
    static Real& updJRot(Vec<N>& y) {return y[D];}
    static const Real& getJRotDot(const Vec<N>& y) {return y[NQ+D];}
    static Real& updJRotDot(Vec<N>& y) {return y[NQ+D];}    
    // Extract the value of the translational Jacobi field from a state y.
    static const Real& getJTrans(const Vec<N>& y) {return y[D+1];}
    static Real& updJTrans(Vec<N>& y) {return y[D+1];}
    static const Real& getJTransDot(const Vec<N>& y) {return y[NQ+D+1];}
    static Real& updJTransDot(Vec<N>& y) {return y[NQ+D+1];}    

private:
    const ContactGeometryImpl& geom;       
};



//==============================================================================
//                      GEODESIC ON PARAMETRIC SURFACE
//==============================================================================

/* This class satisfies the requirements for an Equations template argument
to the GeodesicIntegrator class.
TODO!
*/
class GeodesicOnParametricSurface {
public:
    // state y = q | u = u v jr jt | ud vd jdr jdt
    enum { D = 2,       // 3 coordinates for a point on an implicit surface.
           NJ = 2,      // Number of Jacobi field equations.
           NQ = D+NJ,   // Number of 2nd order equations.
           N  = 2*NQ,   // Number of differential equations.
           NC = 1 };    // maintain unit velocity

    GeodesicOnParametricSurface(const ContactGeometryImpl& geom)
    :   geom(geom) {}

    void calcDerivs(Real t, const Vec<N>& y, Vec<N>& ydot) const {
        assert(!"not implemented");
    }

    void calcConstraintErrors(Real t, const Vec<N>& y, Vec<NC>& cerr) const {
        assert(!"not implemented");
    }     

    bool projectIfNeeded(Real consTol, Real t, Vec<N>& y) const {
        assert(!"not implemented");
        return false;
    }
private:
    const ContactGeometryImpl& geom;       
};


} // namespace SimTK

#endif // SimTK_SIMMATH_INLINE_INTEGRATOR_H_
