/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012-14 Stanford University and the Authors.        *
 * Authors: Ian Stavness                                                      *
 * Contributors: Michael Sherman, Andreas Scholz                              *
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
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Point.h"
#include "simmath/internal/Geo_Sphere.h"
#include "simmath/internal/ContactGeometry.h"

#include "ContactGeometryImpl.h"

#include <iostream>
#include <cmath>
#include <map>
#include <set>

using namespace SimTK;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::cout; using std::endl;


//==============================================================================
//                    CONTACT GEOMETRY :: TORUS & IMPL
//==============================================================================

ContactGeometry::Torus::Torus(Real torusRadius, Real tubeRadius)
:   ContactGeometry(new Torus::Impl(torusRadius, tubeRadius)) {}

/*static*/ ContactGeometryTypeId ContactGeometry::Torus::classTypeId()
{   return ContactGeometry::Torus::Impl::classTypeId(); }

Real ContactGeometry::Torus::getTorusRadius() const {
    return getImpl().getTorusRadius();
}

void ContactGeometry::Torus::setTorusRadius(Real radius) {
    updImpl().setTorusRadius(radius);
}

Real ContactGeometry::Torus::getTubeRadius() const {
    return getImpl().getTubeRadius();
}

void ContactGeometry::Torus::setTubeRadius(Real radius) {
    updImpl().setTubeRadius(radius);
}

const ContactGeometry::Torus::Impl& ContactGeometry::Torus::getImpl() const {
    assert(impl);
    return static_cast<const Torus::Impl&>(*impl);
}

ContactGeometry::Torus::Impl& ContactGeometry::Torus::updImpl() {
    assert(impl);
    return static_cast<Torus::Impl&>(*impl);
}

DecorativeGeometry ContactGeometry::Torus::Impl::createDecorativeGeometry() const {
    PolygonalMesh mesh;
    createPolygonalMesh(mesh);
    return DecorativeMesh(mesh);
}


Vec3 ContactGeometry::Torus::Impl::
findNearestPoint(const Vec3& Q, bool& inside, UnitVec3& normal) const {

    UnitVec3 Zdir(0,0,1);

    // find point P on circle in x-y plane that traces the centroid of the torus
    Vec3 P;
    Vec3 Qproj = Q - Zdir*(~Q*Zdir);
    Real normQproj = Qproj.norm();

    if (std::abs(normQproj) < SimTK::Eps) {
        // Q is along z-axis, therefore there is a locus of closest points,
        // we arbitrarily choose the one in the x-z plane
        P = Vec3(torusRadius,0,0);
    }
    else {
        P = Qproj/normQproj*torusRadius;
    }

    // find direction from centroid of tube to query point, and find near point Qhat
    UnitVec3 Qdir(Q-P);
    Vec3 Qhat = P + Qdir*tubeRadius;

    return Qhat;
}

//TODO
bool ContactGeometry::Torus::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction,
    Real& distance, UnitVec3& normal) const
{
    SimTK_ASSERT_ALWAYS(false, "ContactGeometry::Torus::Impl::intersectsRay unimplemented");
    return false;
}

void ContactGeometry::Torus::Impl::getBoundingSphere
    (Vec3& center, Real& radius) const {
    center = Vec3(0);
    radius = tubeRadius + torusRadius;
}

// Create a polygonal mesh for this torus using parameterization as follows:
// u = [0, 2*Pi] traces a circle in the x-y plane with radius torusRadius,
// which is the centroid of the torus. A point P on this circle is
// given by P = torusRadius*~[cos(u) sin(u) 0].
// v = [0, 2*Pi] traces a circle arond the cross-section (or tube) of the
// torus with radius tubeRadius, at a given u. A point Q on this circle
// is given by Q = (torusRadius + tubeRadius*cos(v))*e1 + tubeRadius*(~[0 0 1]*sin(v))
// where e1 = ~[sin(u) cos(u) 0]. The tube circle is in a plane spanned
// by e1 and the z-axis.
void ContactGeometry::Torus::Impl::createPolygonalMesh(PolygonalMesh& mesh) const {
    // TODO add resolution argument
    const int numSides = 12; //*resolution;
    const int numSlices = 36; //*resolution;

    // add vertices
    for (int i = 0; i < numSlices; ++i) {
      Real u = Real((i*2*SimTK_PI)/numSlices);
      UnitVec3 e1(std::sin(u), std::cos(u), 0); // torus circle aligned with z-axis (z-axis through hole)
      for (int j = 0; j < numSides; ++j) {
        Real v = Real((j*2*SimTK_PI)/numSides);
        Vec3 vtx = (torusRadius + tubeRadius*std::cos(v))*e1 + tubeRadius*std::sin(v)*Vec3(0,0,1); // use ZAXIS?
        mesh.addVertex(vtx);
      }
    }

    // add faces, be careful to wrap indices for the last slice
    int numVertices = mesh.getNumVertices();
//    cout << "num verts = " << numVertices << endl;
    for (int i = 0; i < numVertices; ++i) {
//      cout << "v" << i << ": " << mesh.getVertexPosition(i) << endl;
      // define counter-clockwise quad faces
      Array_<int> faceIndices;
      faceIndices.push_back(i); // u_i,v_i
      faceIndices.push_back((i+1)%numVertices); // u_i, v_i+1
      faceIndices.push_back((i+1+numSides)%numVertices); // u_i+1, v_i+1
      faceIndices.push_back((i+numSides)%numVertices); // u_i+1, v_i
      mesh.addFace(faceIndices);
    }

}

//TODO
void ContactGeometry::Torus::Impl::
calcCurvature(const Vec3& point, Vec2& curvature, Rotation& orientation) const {
    SimTK_ASSERT_ALWAYS(false, "ContactGeometry::Torus::Impl::calcCurvature unimplemented");
}

//TODO
Vec3  ContactGeometry::Torus::Impl::
calcSupportPoint(const UnitVec3& direction) const {
    SimTK_ASSERT_ALWAYS(false, "ContactGeometry::Torus::Impl::calcSupportPoint unimplemented");
    return Vec3(0);
}

Real TorusImplicitFunction::
calcValue(const Vector& x) const {
    return 1-(square(ownerp->getTorusRadius()-std::sqrt(x[0]*x[0]+x[1]*x[1]))+x[2]*x[2])/
            square(ownerp->getTubeRadius());
}

Real TorusImplicitFunction::
calcDerivative(const Array_<int>& derivComponents, const Vector& x) const {
    // first derivatives
    if (derivComponents.size() == 1) {
        if (derivComponents[0]<2) {
            Real sqrt_xy = std::sqrt(x[0]*x[0] + x[1]*x[1]);
            return 2*x[derivComponents[0]]*(ownerp->getTorusRadius() - sqrt_xy)/
                    (square(ownerp->getTubeRadius())*sqrt_xy);
        }
        else
            return -2*x[2]/square(ownerp->getTubeRadius());
    }

    // second derivatives
    if (derivComponents.size() == 2) {
        if (derivComponents[0] < 2) { // fx_ fy_
            if (derivComponents[1] < 2) {
                Real tubeRadiusSq = square(ownerp->getTubeRadius());
                Real xy = x[0]*x[0] + x[1]*x[1];
                Real sqrt_xy = std::sqrt(xy);
                Real den = tubeRadiusSq*xy*sqrt_xy;
                if (derivComponents[0]==derivComponents[1]) { // fxx or fyy
                    int idx = derivComponents[1]==0; // if 0 then idx=1, if 1 then idx=0
                    Real num = 2*ownerp->getTorusRadius()*x[idx]*x[idx];
                    return num/den - 2/tubeRadiusSq;
                }
                else { // fxy or fyx
                    return - 2*ownerp->getTorusRadius()*x[0]*x[1]/den;
                }
            }
            else // fxz = fyz = 0
                return 0;
        }
        else { // fz_
            if (derivComponents[1] == 2) // fzz
                return -2/square(ownerp->getTubeRadius());
            else // fzx = fzy = 0
                return 0;
        }
    }

    //TODO higher order derivatives
    SimTK_ASSERT1_ALWAYS(!"derivative not implemented",
        "Implicit Torus implements 1st&2nd derivs only but %d deriv requested.",
        derivComponents.size());
    return 0;
}

