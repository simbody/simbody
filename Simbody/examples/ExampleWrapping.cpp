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

Real pauseBetweenPathIterations = 1; // sec
Real estimatedPathErrorAccuracy = 1e-12;

const Real vizInterval = Infinity; // set to 1/30. to vizualize shooting


class VizPeriodicReporter : public PeriodicEventReporter {
public:
    VizPeriodicReporter(const Visualizer& viz, const State& dummyState, Real interval) :
        PeriodicEventReporter(interval), viz(viz), dummyState(dummyState) {
    }

    void handleEvent(const State& state) const {
        viz.report(dummyState);
    }

private:
    const Visualizer& viz;
    const State& dummyState;
};

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
        Vec3 nP = gg.getGeom().calcSurfaceNormal((Vector)P);
        Vec3 nQ = gg.getGeom().calcSurfaceNormal((Vector)Q);
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

    // temporary variables
    Geodesic& geod;

}; // class PathError

/*
 * This class is used to calculate the path error for a single obstacle
 * using the split geodesic error
 */
class PathErrorSplit: public Differentiator::JacobianFunction {

public:
    PathErrorSplit(int nf, int ny, GeodesicGeometry& geodgeom, Geodesic& geod,
            const Vec3& O, const Vec3& I) :
            Differentiator::JacobianFunction(nf, ny),
                    gg(geodgeom), geod(geod),
                    O(O), I(I) { }

    // x = ~[P, Q]
    int f(const Vector& x, Vector& fx) const  {
        Vec3 P(&x[0]);
        Vec3 Q(&x[3]);

        // calculate plane bisecting P and Q, and use as termination condition for integrator
        UnitVec3 normal(Q-P);
        Real offset = (~(P+Q)*normal)/2 ;
        gg.setPlane(Plane(normal, offset));

        UnitVec3 r_OP(P-O);
        UnitVec3 r_QI(I-Q);

        geod.clear();
        Vec2 geodErr = gg.calcGeodError(P, Q, r_OP, -r_QI);

        Vec3 nP = gg.getGeom().calcSurfaceNormal((Vector)P);
        Vec3 nQ = gg.getGeom().calcSurfaceNormal((Vector)Q);

        fx[0] = ~r_OP*nP;
        fx[1] = ~r_QI*nQ;
        fx[2] = geodErr[0];
        fx[3] = geodErr[1];
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

    // temporary variables
    Geodesic& geod;

}; // class PathErrorSplit

static Real maxabs(Vector x) {
    Real maxVal = 0;
    for (int i = 0; i < x.size(); ++i) {
        if (std::abs(x[i]) > maxVal)
            maxVal = std::abs(x[i]);
    }
    return maxVal;
}

static Real maxabsdiff(Vector x, Vector xold) {
    ASSERT(x.size()==xold.size());
    Real maxVal = 0;
    for (int i = 0; i < x.size(); ++i) {
        if (std::abs(x[i]-xold[i])/std::max(x[i],1.0) > maxVal)
            maxVal = std::abs(x[i]-xold[i])/std::max(x[i],1.0);
    }
    return maxVal;
}


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

    ContactGeometry::Sphere geom(r);
//    r = 2;
//    Vec3 radii(1,2,3);
//    ContactGeometry::Ellipsoid geom(radii);
    GeodesicGeometry geodgeom(geom);
    Geodesic geod;

    // Create a dummy MultibodySystem for visualization purposes
    MultibodySystem dummySystem;
    SimbodyMatterSubsystem matter(dummySystem);
//    matter.updGround().addBodyDecoration(Transform(), DecorativeEllipsoid(radii)
    matter.updGround().addBodyDecoration(Transform(), DecorativeSphere(r)
            .setColor(Gray)
            .setOpacity(0.5)
            .setResolution(5));

    // Visualize with default options; ask for a report every 1/30 of a second
    // to match the Visualizer's default 30 frames per second rate.
    Visualizer viz(dummySystem);
    viz.setBackgroundType(Visualizer::SolidColor);
    dummySystem.addEventReporter(new Visualizer::Reporter(viz, 1./30));

    // add vizualization callbacks for geodesics, contact points, etc.
    viz.addDecorationGenerator(new GeodesicDecorator(geodgeom.getGeodP(), Red));
    viz.addDecorationGenerator(new GeodesicDecorator(geodgeom.getGeodQ(), Blue));
    viz.addDecorationGenerator(new GeodesicDecorator(geod, Orange));
    viz.addDecorationGenerator(new PlaneDecorator(geodgeom.getPlane(), Gray));
    viz.addDecorationGenerator(new PathDecorator(x, O, I, Green));
    dummySystem.realizeTopology();
    State dummyState = dummySystem.getDefaultState();


    // calculate the geodesic
    geodgeom.addVizReporter(new VizPeriodicReporter(viz, dummyState, vizInterval));
    viz.report(dummyState);

    // creat path error function
//    PathError pathErrorFnc(n, n, geodgeom, geod, O, I);
    PathErrorSplit pathErrorFnc(n, n, geodgeom, geod, O, I);
    pathErrorFnc.setEstimatedAccuracy(estimatedPathErrorAccuracy);
    Differentiator diff(pathErrorFnc);

    // set initial conditions
    x[0]=P[0]; x[1]=P[1]; x[2]=P[2];
    x[3]=Q[0]; x[4]=Q[1]; x[5]=Q[2];

    Real f, fold, lam;

    pathErrorFnc.f(x, Fx);
    viz.report(dummyState);
    usleep((useconds_t)(pauseBetweenPathIterations*1000000));

    f = std::sqrt(~Fx*Fx);
    for (int i = 0; i < maxNewtonIterations; ++i) {
        if (f < ftol) {
            std::cout << "path converged in " << i << " iterations" << std::endl;
//            cout << "obstacle err = " << Fx << ", x = " << x << endl;
            break;
        }

        diff.calcJacobian(x, Fx, J, Differentiator::ForwardDifference);
        dx = J.invert()*Fx;

        fold = f;
        xold = x;

        // backtracking
        lam = 1;
        while (true) {
            x = xold - lam*dx;
            pathErrorFnc.f(x, Fx);
            f = std::sqrt(~Fx*Fx);
            if (f > fold && lam > minlam) {
                lam = lam / 2;
            } else {
                break;
            }
        }
        if (maxabsdiff(x,xold) < xtol) {
            std::cout << "converged on step size after " << i << " iterations" << std::endl;
            std::cout << "error = " << Fx << std::endl;
            break;
        }
        viz.report(dummyState);
        usleep((useconds_t)(pauseBetweenPathIterations*1000000));

    }
    cout << "obstacle error = " << Fx << endl;

    cout << "num geodP pts = " << geodgeom.getGeodP().getPoints().size() << endl;


  } catch (const std::exception& e) {
    std::printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);

  } catch (...) {
    std::printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

  return 0;
}
