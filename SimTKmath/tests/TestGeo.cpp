/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-12 Stanford University and the Authors.        *
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

/* Tests for low-level geometric primitives and algorithms. */

#include "SimTKmath.h"
#include <vector>
#include <exception>

using namespace SimTK;
using namespace std;

static double Tol = NTraits<double>::getSignificant();
static float  fTol = NTraits<float>::getSignificant();


template <class P>
static void checkSphere(const Geo::Sphere_<P>& sph, const Array_<Vec<3,P> > pts) {
    const P radius = sph.getRadius();
    for (int i=0; i < (int)pts.size(); ++i) {
        const P dist = (pts[i]-sph.getCenter()).norm();
        SimTK_TEST(dist <= radius);
    }
}

static void addOctohedron(vector<Vec3>& vertices, vector<int>& faceIndices, 
                          Vec3 offset) {
    int start = (int)vertices.size();
    vertices.push_back(Vec3(0, 1, 0)+offset);
    vertices.push_back(Vec3(1, 0, 0)+offset);
    vertices.push_back(Vec3(0, 0, 1)+offset);
    vertices.push_back(Vec3(-1, 0, 0)+offset);
    vertices.push_back(Vec3(0, 0, -1)+offset);
    vertices.push_back(Vec3(0, -1, 0)+offset);
    int faces[8][3] = {{0, 2, 1}, {0, 3, 2}, {0, 4, 3}, {0, 1, 4}, 
                       {5, 1, 2}, {5, 2, 3}, {5, 3, 4}, {5, 4, 1}};
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]+start);
}

// From Peter Eastman's tri mesh tests.
void testTriMeshBoundingSphere() {
    Random::Uniform random(0, 10);
    Real ratio=0, worst=0, best=Infinity;
    const int NTrials = 100;

    for (int i = 0; i < NTrials; i++) {
        // Create a mesh consisting of a random number of octohedra at random 
        // places.
        
        vector<Vec3> vertices;
        vector<int> faceIndices;
        int numOctohedra = random.getIntValue()+1;
        for (int j = 0; j < numOctohedra; j++)
            addOctohedron(vertices, faceIndices, 
            Vec3(random.getValue(), random.getValue(), random.getValue()));
        ContactGeometry::TriangleMesh mesh(vertices, faceIndices);

        // Verify that all points are inside the bounding sphere.
        
        Vec3 center;
        Real radius;
        mesh.getBoundingSphere(center, radius);
        for (int j = 0; j < mesh.getNumVertices(); j++) {
            Real dist = (center-mesh.getVertexPosition(j)).norm();
            SimTK_TEST(dist <= radius);
        }

        
        // Make sure the bounding sphere is reasonably compact.      
        Vec3 boxRadius = 0.5*mesh.getOBBTreeNode().getBounds().getSize();
        SimTK_TEST(radius <= boxRadius.norm());

        // Compare with fast & crude Ritter sphere. On lucky occasions the
        // Ritter sphere can be just as good, and then roundoff might make
        // it trivially better but it shouldn't ever actually *be* better.
        Geo::Sphere ritter = Geo::Point::calcApproxBoundingSphere(vertices);
        SimTK_TEST(radius <= ritter.getRadius()*(1.01));

        const Real bsoas = cube(radius/ritter.getRadius());
        ratio += bsoas;
        if (bsoas > worst) worst=bsoas;
        if (bsoas < best)  best=bsoas;
    }
    ratio /= NTrials;
    printf("avg ratio=%g worst=%g best=%g\n", ratio, worst, best);
    SimTK_TEST(ratio <= Real(.85)); // volume ratio
    SimTK_TEST(worst <= Real(1.3));

}

// Generate many sets of random points, at difficult far-away places, then
// generate single- and double-precision bounding spheres and check them
// for admissibility (all points in) and optimality (no bigger than needed).
// For these the right answer is hard to come by so we can only check that
// the minimal spheres are never larger than Ritter spheres.
void testRandomPoints() {
    Random::Uniform random(0, 1000);
    Array_<Vec3> pts;
    Array_<fVec3> fpts;
    Real ratio=0, fratio=0, worst=0, fworst=0, best=Infinity, fbest=Infinity;
    const int NTrials = 10000;
    // TODO: At around 5,000,000 a case is generated where the "minimal" sphere
    // is more than 20% larger than the (bad) Ritter sphere.
    //const int NTrials = 10000000; 
    for (int trial=0; trial<NTrials; ++trial) {
        pts.clear(); fpts.clear();
        int numPoints = random.getIntValue()+1;
        Vec3 offs = Test::randDouble()*1000*Test::randVec3();
        fVec3 foffs((float)offs[0],(float)offs[1],(float)offs[2]);
        Real scale = offs.norm();

        for (int p=0; p<numPoints; ++p) {
            Vec3  pt(Test::randVec3());
            fVec3 fpt((float)pt[0],(float)pt[1],(float)pt[2]);
            pts.push_back(pt+offs); fpts.push_back(fpt+foffs);
        }

        Geo::Sphere bs = Geo::Point::calcBoundingSphere(pts);
        checkSphere(bs, pts);
        Geo::Sphere_<float> fbs = 
            Geo::Point_<float>::calcBoundingSphere(fpts);
        checkSphere(fbs, fpts);
        Geo::Sphere as = Geo::Point::calcApproxBoundingSphere(pts);
        checkSphere(as, pts);
        Geo::Sphere_<float> fas = 
            Geo::Point_<float>::calcApproxBoundingSphere(fpts);
        checkSphere(fas, fpts);

        const Real bsoas = cube(bs.getRadius()/as.getRadius());
        const Real fbsofas = cube((Real)(fbs.getRadius()/fas.getRadius()));
        ratio += bsoas; fratio += fbsofas;
        if (bsoas > worst) worst=bsoas;
        if (bsoas < best)  best=bsoas;
        if (fbsofas > fworst) fworst=fbsofas;
        if (fbsofas < fbest)  fbest=fbsofas;

        // The single and double precision spheres should be the same size
        // to within a small error. The Ritter sphere is more sensitive.
        const float frac = std::max((float)scale,1.f)*NTraits<float>::getSqrtEps();
        
        //if (!Test::numericallyEqual((float)bs.getRadius(), fbs.getRadius(), 1, frac))
        //    printf("bs=%g fbs=%g\n", bs.getRadius(), fbs.getRadius());
        SimTK_TEST_EQ_TOL((float)bs.getRadius(), fbs.getRadius(), 0.2f);
        SimTK_TEST_EQ_TOL((float)as.getRadius(), fas.getRadius(), 0.2f);

        // Compare Welzl spheres with fast & crude Ritter spheres. On lucky 
        // occasions the Ritter spheres can be just as good, and then roundoff
        // might make them trivially better but they shouldn't ever actually 
        // *be* better. TODO: check at the end. There are obscure cases 
        // where the minimal sphere is too big.
        //SimTK_TEST(bs.getRadius() <= as.getRadius()*(1.2));
        //SimTK_TEST(fbs.getRadius() <= fas.getRadius()*(1.2f));
    }
    ratio /= NTrials; fratio /= NTrials;
    printf("avg ratio=%g worst=%g best=%g\n", ratio, worst, best);
    printf("avg fratio=%g worst=%g best=%g\n", fratio, fworst, fbest);
    SimTK_TEST(ratio <= Real(.85)); // volume ratio
    SimTK_TEST(worst <= Real(1.3));
    SimTK_TEST(fratio <= Real(.85)); // volume ratio
    SimTK_TEST(fworst <= Real(1.3));}

void testCoplanarPoints() {
    // TODO
}

void testCosphericalPoints() {
    // TODO
}


void testCollinearPoints() {
    Random::Uniform random(0, 1000);
    Array_<Vec3> pts;
    Array_<fVec3> fpts;
    for (int trial=0; trial<1000; ++trial) {
        pts.clear(); fpts.clear();
        int numPoints = random.getIntValue()+1;
        Vec3 offs = Test::randDouble()*1000*Test::randVec3();
        UnitVec3 dir(Test::randVec3());
        fVec3 foffs((float)offs[0],(float)offs[1],(float)offs[2]);
        fUnitVec3 fdir((float)dir[0],(float)dir[1],(float)dir[2]);

        int minpos, maxpos;
        Real minval=Infinity, maxval=-Infinity;
        for (int p=0; p<numPoints; ++p) {
            Real pos = 100*Test::randDouble();
            if (pos > maxval) maxval=pos, maxpos=p;
            if (pos < minval) minval=pos, minpos=p;
            Vec3  pt(offs+pos*dir);
            fVec3 fpt((float)pt[0],(float)pt[1],(float)pt[2]);
            pts.push_back(pt); fpts.push_back(fpt);
        }

        Real radius = (pts[maxpos]-pts[minpos]).norm()/2;

        Geo::Sphere bs = Geo::Point::calcBoundingSphere(pts);
        checkSphere(bs, pts);
        Geo::Sphere_<float> fbs = 
            Geo::Point_<float>::calcBoundingSphere(fpts);
        checkSphere(fbs, fpts);
        Geo::Sphere as = Geo::Point::calcApproxBoundingSphere(pts);
        checkSphere(as, pts);
        Geo::Sphere_<float> fas = 
            Geo::Point_<float>::calcApproxBoundingSphere(fpts);
        checkSphere(fas, fpts);

        Real scale = std::max(std::max(max(offs.abs()), radius), One);
        float ftol = float(scale)*fTol;
        double tol = scale*Tol;


        SimTK_TEST_EQ_TOL(bs.getRadius(), radius, tol);
        SimTK_TEST_EQ_TOL(fbs.getRadius(), radius, ftol);
        SimTK_TEST_EQ_TOL(as.getRadius(), radius, tol);
        SimTK_TEST_EQ_TOL(fas.getRadius(), radius, ftol);

        // Repeat test with random noise added.
        for (int p=0; p<numPoints; ++p) {
            Vec3 noise(Test::randVec3());
            fVec3 fnoise((float)noise[0],(float)noise[1],(float)noise[2]);
            pts[p] += SignificantReal*noise;
            fpts[p] += NTraits<float>::getSignificant()*fnoise;
        }

        bs = Geo::Point::calcBoundingSphere(pts);
        checkSphere(bs, pts);
        fbs = Geo::Point_<float>::calcBoundingSphere(fpts);
        checkSphere(fbs, fpts);
        as = Geo::Point::calcApproxBoundingSphere(pts);
        checkSphere(as, pts);
        fas = Geo::Point_<float>::calcApproxBoundingSphere(fpts);
        checkSphere(fas, fpts);

        SimTK_TEST_EQ_TOL(bs.getRadius(), radius, tol);
        SimTK_TEST_EQ_TOL(fbs.getRadius(), radius, ftol);
        SimTK_TEST_EQ_TOL(as.getRadius(), radius, tol);
        SimTK_TEST_EQ_TOL(fas.getRadius(), radius, ftol);
    }
}

void testCollocatedPoints() {
    Random::Uniform random(0, 1000);
    Array_<Vec3> pts;
    Array_<fVec3> fpts;
    for (int trial=0; trial<1000; ++trial) {
        pts.clear(); fpts.clear();
        int numPoints = random.getIntValue()+1;
        Vec3 offs = Test::randDouble()*1000*Test::randVec3();
        fVec3 foffs((float)offs[0],(float)offs[1],(float)offs[2]);
        Real scale = offs.norm();

        for (int p=0; p<numPoints; ++p) {
            Vec3  pt(offs);
            fVec3 fpt((float)pt[0],(float)pt[1],(float)pt[2]);
            pts.push_back(pt); fpts.push_back(fpt);
        }

        Geo::Sphere bs = Geo::Point::calcBoundingSphere(pts);
        checkSphere(bs, pts);
        Geo::Sphere_<float> fbs = 
            Geo::Point_<float>::calcBoundingSphere(fpts);
        checkSphere(fbs, fpts);
        Geo::Sphere as = Geo::Point::calcApproxBoundingSphere(pts);
        checkSphere(as, pts);
        Geo::Sphere_<float> fas = 
            Geo::Point_<float>::calcApproxBoundingSphere(fpts);
        checkSphere(fas, fpts);

        SimTK_TEST(bs.getRadius() < SqrtEps)
        SimTK_TEST(fbs.getRadius() < NTraits<float>::getSqrtEps());
        SimTK_TEST(as.getRadius() < SqrtEps)
        SimTK_TEST(fas.getRadius() < NTraits<float>::getSqrtEps());

        // Repeat test with random noise added.
        for (int p=0; p<numPoints; ++p) {
            Vec3 noise(Test::randVec3());
            fVec3 fnoise((float)noise[0],(float)noise[1],(float)noise[2]);
            pts[p] += SignificantReal*noise;
            fpts[p] += NTraits<float>::getSignificant()*fnoise;
        }

        bs = Geo::Point::calcBoundingSphere(pts);
        checkSphere(bs, pts);
        fbs = Geo::Point_<float>::calcBoundingSphere(fpts);
        checkSphere(fbs, fpts);
        as = Geo::Point::calcApproxBoundingSphere(pts);
        checkSphere(as, pts);
        fas = Geo::Point_<float>::calcApproxBoundingSphere(fpts);
        checkSphere(fas, fpts);
    }
}

void testBox() {
    Geo::Box box(Vec3(3,4,2)); // half lengths
    SimTK_TEST(box.findVolume() == 8*24);
    SimTK_TEST(box.getOrderedAxis(0) == ZAxis); // smallest
    SimTK_TEST(box.getOrderedAxis(1) == XAxis); // medium
    SimTK_TEST(box.getOrderedAxis(2) == YAxis); // largest

    Geo::AlignedBox abox(Vec3(0), Vec3(1,2,3));
    abox.setCenter(Vec3(3,4,2)+Vec3(1,2,3)-Vec3(1e-6)); 
    SimTK_TEST(box.intersectsAlignedBox(abox));

    abox.setCenter(Vec3(3,4,2)+Vec3(1,2,3)+Vec3(1e-6)); 
    SimTK_TEST(!box.intersectsAlignedBox(abox));

    Geo::OrientedBox obox(Transform(), Vec3(1,2,3));
    SimTK_TEST(box.mayIntersectOrientedBox(obox)); // centers overlap
    SimTK_TEST(box.intersectsOrientedBox(obox));

    obox.setTransform(Vec3(10,0,0)); // x axis should separate
    SimTK_TEST(!box.mayIntersectOrientedBox(obox));
    SimTK_TEST(!box.intersectsOrientedBox(obox));

    obox.setTransform(Vec3(3.123,-1.3,.7)); // parallel boxes that intersect
    SimTK_TEST(box.mayIntersectOrientedBox(obox));
    SimTK_TEST(box.intersectsOrientedBox(obox));

    // Non-intersecting box for which no face will serve as separator.
    // In this case the fast method can't tell they are separated.
    obox.setTransform(Transform(
        Rotation(BodyRotationSequence, Pi/4, XAxis, Pi/8, YAxis, -Pi/4, ZAxis),
        Vec3(1.5, -5, 5.25)));
    SimTK_TEST(box.mayIntersectOrientedBox(obox));
    SimTK_TEST(!box.intersectsOrientedBox(obox));

    // This should make them intersect.
    obox.setTransform(Transform(
        Rotation(BodyRotationSequence, Pi/4, XAxis, Pi/8, YAxis, -Pi/4, ZAxis),
        Vec3(1.5, -5, 4.5)));
    SimTK_TEST(box.mayIntersectOrientedBox(obox));
    SimTK_TEST(box.intersectsOrientedBox(obox));

}

void testMiscGeo() {

    // TEST     GEO::POINT::POINTS ARE NUMERICALLY COINCIDENT
    Vec3 p1, p2, p3;
    p1 = Vec3(1,2,3);
    p2 = Vec3(2,3,4);
    p3 = p1;
    SimTK_TEST(!Geo::Point::pointsAreNumericallyCoincident(p1,p2));
    SimTK_TEST(Geo::Point::pointsAreNumericallyCoincident(p1,p3));
    p3 = p1 + Eps*p2;
    // Default tolerance is larger than Eps.
    SimTK_TEST(Geo::Point::pointsAreNumericallyCoincident(p1,p3));
    // This should be enough to separate them if they are near the origin.
    p3 = p1 + 10*Geo::getDefaultTol<Real>()*p2;
    SimTK_TEST(!Geo::Point::pointsAreNumericallyCoincident(p1,p3));
    // Shifting by 100 should mean that they are indistinguishable when
    // perturbed by 10*tol.
    p3 = 100.*p1 + 10*Geo::getDefaultTol<Real>()*p2;
    SimTK_TEST(Geo::Point::pointsAreNumericallyCoincident(100.*p1,p3));
    // But they should be distinguishable at 1000*tol.
    p3 = 100.*p1 + 1000*Geo::getDefaultTol<Real>()*p2;
    SimTK_TEST(!Geo::Point::pointsAreNumericallyCoincident(100.*p1,p3));
    // And again indistinguishable at a looser tolerance.
    SimTK_TEST(Geo::Point::pointsAreNumericallyCoincident(100.*p1,p3,
        10000.*Geo::getDefaultTol<Real>()));
}

int main() {
    SimTK_START_TEST("TestGeo");
        SimTK_SUBTEST(testMiscGeo);
        SimTK_SUBTEST(testBox);
        SimTK_SUBTEST(testTriMeshBoundingSphere);
        SimTK_SUBTEST(testRandomPoints);
        SimTK_SUBTEST(testCollinearPoints);
        SimTK_SUBTEST(testCollocatedPoints);
    SimTK_END_TEST();
}
