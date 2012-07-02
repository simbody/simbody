/* -------------------------------------------------------------------------- *
 *                         Simbody(tm)  ExampleWrapping                           *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Ian Stavness, Michael Sherman                                     *
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

/*
 * This example demonstrates cable path wrapping for one obstacle
 */

#include "Simbody.h"
#include "simmath/internal/Geodesic.h"
#include "simmath/internal/GeodesicGeometry.h"
#include "simmath/internal/ContactGeometry.h"

using namespace SimTK;
using std::cos;
using std::sin;
using std::cout;
using std::endl;

Real pauseDurationBetweenIterations = 1; // sec
Real estimatedPathErrorAccuracy = 1e-12;
Real ftol = 1e-12;
Real xtol = 1e-6;

/*
 * This class is used to calculate the path error for a single obstacle
 */
class PathError: public Differentiator::JacobianFunction {

public:
    PathError(int nf, int ny, GeodesicGeometry& geodgeom, Geodesic& geod,
            const Vec3& O, const Vec3& I) :
            Differentiator::JacobianFunction(nf, ny),
                    gg(geodgeom), geod(geod),
                    O(O), I(I) { }

    // x = ~[P, Q]
    int f(const Vector& x, Vector& fx) const  {
        Vec3 P(&x[0]);
        Vec3 Q(&x[3]);

        UnitVec3 r_OP(P-O);
        UnitVec3 r_QI(I-Q);

        geod.clear();
        gg.calcGeodesic(P, Q, r_OP, -r_QI, geod);

        Vec3 tP = geod.getTangents()[0];
        Vec3 tQ = geod.getTangents()[geod.getTangents().size()-1];
        Vec3 nP = gg.getGeom().calcSurfaceGradient((Vector)P);
        Vec3 nQ = gg.getGeom().calcSurfaceGradient((Vector)Q);
        Vec3 bP = cross(nP, tP);
        Vec3 bQ = cross(nQ, tQ);

        fx[0] = ~r_OP*nP;
        fx[1] = ~r_QI*nQ;
        fx[2] = ~r_OP*bP;
        fx[3] = ~r_QI*bQ;
        fx[4] = gg.getGeom().calcSurfaceValue((Vector)P);
        fx[5] = gg.getGeom().calcSurfaceValue((Vector)Q);

        return 0;
    }

    const Geodesic& getGeod() {
        return geod;
    }

private:
    GeodesicGeometry& gg;
    const Vec3& O;
    const Vec3& I;
    const Rotation R_SP;
    const Rotation R_SQ;

    // temporary variables
    Geodesic& geod;

}; // class PathError



int main() {
  try {

    // setup test problem
    double r = 1;
    double uP = -Pi/2;
    double vP = Pi/3;
    double uQ = 0;
    double vQ = 2;
    Vec3 O(-r, -r, 0.2);
    Vec3 I(r, r, -r);
    Vec3 P(r*cos(uP)*sin(vP), r*sin(uP)*sin(vP), r*cos(vP));
    Vec3 Q(r*cos(uQ)*sin(vQ), r*sin(uQ)*sin(vQ), r*cos(vQ));

    Vec3 r_OP = P-O;
    Vec3 r_IQ = Q-I;
    Vec3 tP = r_OP.normalize();
    Vec3 tQ = r_IQ.normalize();

    int n = 6; // problem size
    Vector x(n), dx(n), Fx(n), xold(n);
    Matrix J(n,n);

    ContactGeometry::Sphere sphere(r);
    GeodesicGeometry geodgeom(sphere);
    Geodesic geod;

    // Create a dummy MultibodySystem for visualization purposes
    MultibodySystem dummySystem;
    SimbodyMatterSubsystem matter(dummySystem);
    matter.updGround().addBodyDecoration(Transform(), DecorativeSphere(r)
            .setColor(Gray)
            .setOpacity(0.5));

    // Visualize with default options; ask for a report every 1/30 of a second
    // to match the Visualizer's default 30 frames per second rate.
    Visualizer viz(dummySystem);
    viz.setBackgroundType(Visualizer::SolidColor);
    dummySystem.addEventReporter(new Visualizer::Reporter(viz, 1./30));

    // add vizualization callbacks for geodesics, contact points, etc.
    viz.addDecorationGenerator(new GeodesicDecorator(geodgeom.getGeodP(), Red));
    viz.addDecorationGenerator(new GeodesicDecorator(geodgeom.getGeodQ(), Blue));
    //    viz.addDecorationGenerator(new GeodesicDecorator(geod, Orange));
    viz.addDecorationGenerator(new PlaneDecorator(geodgeom.getPlane(), Gray));
    viz.addDecorationGenerator(new PathDecorator(x, O, I, Green));
    dummySystem.realizeTopology();
    State dummyState = dummySystem.getDefaultState();

    // creat path error function
    PathError pathErrorFnc(n, n, geodgeom, geod, O, I);
    pathErrorFnc.setEstimatedAccuracy(estimatedPathErrorAccuracy);
    Differentiator diff(pathErrorFnc);

    // set initial conditions
    x[0]=P[0]; x[1]=P[1]; x[2]=P[2];
    x[3]=Q[0]; x[4]=Q[1]; x[5]=Q[2];

    pathErrorFnc.f(x, Fx);
    viz.report(dummyState);
    sleep(pauseDurationBetweenIterations);

    Real f = 0.5*~Fx*Fx;
    Real fold, lam = 1;


    int maxIter = 40;
    for (int i = 0; i < maxIter; ++i) {
        if (std::sqrt(f) < ftol) {
            std::cout << "path converged in " << i << " iterations" << std::endl;
            break;
        }
//        cout << "obstacle err = " << Fx << ", x = " << x << endl;

        diff.calcJacobian(x, Fx, J, Differentiator::ForwardDifference);
        fold = f;
        xold = x;
//        cout << "J = " << J << endl;
        dx = J.invert()*Fx;

        // backtracking
        lam = 1;
//        x = xold - lam*dx;
//        obstacleError.f(x, Fx);
        while (true) {
            x = xold - lam*dx;
            pathErrorFnc.f(x, Fx);
            f = 0.5*~Fx*Fx;
            if (f > fold) {
                lam = lam / 2;
            } else {
                break;
            }
        }

        viz.report(dummyState);
        sleep(pauseDurationBetweenIterations);

    }
    cout << "obstacle error = " << Fx << endl;

    cout << "num geod pts = " << geod.getPoints().size() << endl;


  } catch (const std::exception& e) {
    std::printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);

  } catch (...) {
    std::printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

    return 0;
}
