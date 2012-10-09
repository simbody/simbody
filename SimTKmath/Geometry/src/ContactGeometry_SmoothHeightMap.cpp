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
 * Contributors: Matthew Millard                                              *
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
#include "simmath/internal/Geo_BicubicBezierPatch.h"
#include "simmath/internal/BicubicSurface.h"
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
//               CONTACT GEOMETRY :: SMOOTH HEIGHT MAP & IMPL
//==============================================================================

ContactGeometry::SmoothHeightMap::
SmoothHeightMap(const BicubicSurface& surface) 
:   ContactGeometry(new SmoothHeightMap::Impl(surface)) {}

/*static*/ ContactGeometryTypeId ContactGeometry::SmoothHeightMap::
classTypeId() 
{   return ContactGeometry::SmoothHeightMap::Impl::classTypeId(); }

const BicubicSurface& ContactGeometry::SmoothHeightMap::
getBicubicSurface() const {return getImpl().getBicubicSurface();}

const OBBTree& ContactGeometry::SmoothHeightMap::
getOBBTree() const {return getImpl().getOBBTree();}

const ContactGeometry::SmoothHeightMap::Impl& ContactGeometry::SmoothHeightMap::
getImpl() const {
    assert(impl);
    return static_cast<const SmoothHeightMap::Impl&>(*impl);
}

ContactGeometry::SmoothHeightMap::Impl& ContactGeometry::SmoothHeightMap::
updImpl() {
    assert(impl);
    return static_cast<SmoothHeightMap::Impl&>(*impl);
}

// This is the main constructor.
ContactGeometry::SmoothHeightMap::Impl::
Impl(const BicubicSurface& surface) 
:   surface(surface) { 
    implicitFunction.setOwner(*this); 

    createBoundingVolumes();

}

void ContactGeometry::SmoothHeightMap::Impl::
assignPatch(const Geo::BicubicBezierPatch& patch, 
            OBBNode& node, int depth,
            Array_<const Vec3*>* parentControlPoints) const 
{
    const Mat<4,4,Vec3>& nodeB = patch.getControlPoints();
    const Vec2& nodeB11 = nodeB(0,0).getSubVec<2>(0); // just x,y
    const Vec2& nodeB44 = nodeB(3,3).getSubVec<2>(0);
    node.centerUW = (nodeB11+nodeB44)/2;
    node.dims     = (nodeB11-nodeB44).abs()/2;

    // For now just split 4 ways; need to be done recursively based on
    // flatness of patch.
    node.children.resize(4);
    patch.split(0.5,0.5,node.children[0].patch, node.children[1].patch,
                        node.children[2].patch, node.children[3].patch);
    Array_<const Vec3*> myControlPoints;
    for (int c=0; c<4; ++c) {
        OBBNode& child = node.children[c];
        child.depth = depth+1;
        child.height = 0;
        child.box = child.patch.calcOrientedBoundingBox();
        const Mat<4,4,Vec3>& B = child.patch.getControlPoints();
        const Vec2& b11 = B(0,0).getSubVec<2>(0); // just x,y
        const Vec2& b44 = B(3,3).getSubVec<2>(0);
        child.centerUW = (b11+b44)/2;
        child.dims     = (b11-b44).abs()/2;
        for (int i=0; i<4; ++i) 
            for (int j=0; j<4; ++j) 
                myControlPoints.push_back(&B(i,j));
    }

    node.depth = depth;
    node.height = 1 + std::max(node.children[0].height,
                      std::max(node.children[1].height,
                      std::max(node.children[2].height,
                               node.children[3].height)));
    node.box = Geo::Point::calcOrientedBoundingBoxIndirect(myControlPoints);

    if (parentControlPoints)
        for (unsigned i=0; i<myControlPoints.size(); ++i)
            parentControlPoints->push_back(myControlPoints[i]);
}

void ContactGeometry::SmoothHeightMap::Impl::
splitPatches(int x0,int y0, int nx, int ny, 
             OBBNode& node, int depth,
             Array_<const Vec3*>* parentControlPoints) const {
    assert(nx>0 && ny>0 && depth>=0);


    node.x0=0; node.y0=0; node.nx=nx; node.ny=ny;
    if (nx==1 && ny==1) {
        assignPatch(surface.calcBezierPatch(x0,y0), node, depth,
                    parentControlPoints);
        return;
    } 

    // Add two children.
    node.children.resize(2);
    Array_<const Vec3*> myControlPoints;

    // Split on the long direction
    if (nx > ny) {
        splitPatches(x0,      y0, nx/2,    ny, node.children[0],
            depth+1, &myControlPoints);
        splitPatches(x0+nx/2, y0, nx-nx/2, ny, node.children[1],
            depth+1, &myControlPoints);
    } else {
        splitPatches(x0, y0,      nx, ny/2,    node.children[0],
            depth+1, &myControlPoints);
        splitPatches(x0, y0+ny/2, nx, ny-ny/2, node.children[1],
            depth+1, &myControlPoints);
    }
    node.depth = depth;
    node.height = 1 + std::max(node.children[0].height,
                               node.children[1].height);

    node.box = Geo::Point::calcOrientedBoundingBoxIndirect(myControlPoints);

    if (parentControlPoints) {
        for (unsigned i=0; i<myControlPoints.size(); ++i)
            parentControlPoints->push_back(myControlPoints[i]);
    }
}

void ContactGeometry::SmoothHeightMap::Impl::
createBoundingVolumes() {
    // Temporarily convert the surface into a set of Bezier patches (using
    // a lot more memory than the original).

    int nx,ny; surface.getNumPatches(nx,ny);
    OBBNode& root = obbTree.updRoot();
    splitPatches(0,0,nx,ny,root,0);


    // Create bounding sphere.
    // TODO: fake this using mesh; this needs to be done correctly instead
    // by the BicubicSurface itself. Using 5 subdivisions per patch.
    PolygonalMesh mesh = surface.createPolygonalMesh(5);

    // Collect all the vertices.
    const int n = mesh.getNumVertices();
    Array_<const Vec3*> points(n);
    for (int i=0; i<n; ++i)
        points[i] = &mesh.getVertexPosition(i);
    boundingSphere = Geo::Point::calcBoundingSphereIndirect(points);
    // Add 10% as a hack to make it less likely we'll miss part of the surface.
    boundingSphere.updRadius() *= 1.1;
}

DecorativeGeometry ContactGeometry::SmoothHeightMap::Impl::createDecorativeGeometry() const {
    return DecorativeMesh(surface.createPolygonalMesh());
}

Vec3 ContactGeometry::SmoothHeightMap::Impl::
findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    assert(false);
    return Vec3(NaN);
}

bool ContactGeometry::SmoothHeightMap::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction, 
    Real& distance, UnitVec3& normal) const 
{
    assert(false);
    return true;
}

Real SmoothHeightMapImplicitFunction::
calcValue(const Vector& p) const {
    const BicubicSurface&      surf = ownerp->getBicubicSurface();
    BicubicSurface::PatchHint& hint = ownerp->updHint();
    const Real z = surf.calcValue(Vec2(p[0],p[1]), hint);
    //TODO: this is negated from convention
    return z - p[2]; // negative outside, positive inside
}

// First deriv with respect to p[2] (z component) is -1 to match the above
// implicit function definition, all higher derivs are
// with respect to that component are 0.
Real SmoothHeightMapImplicitFunction::
calcDerivative(const Array_<int>& derivComponents, const Vector& p) const {
    if (derivComponents.empty()) return calcValue(p);
    if (derivComponents.size() == 1 && derivComponents[0]==2)
        return -1;
    for (unsigned i=0; i<derivComponents.size(); ++i)
        if (derivComponents[i]==2) return 0;

    // We're asking only for derivatives in x and y.
    const BicubicSurface&      surf = ownerp->getBicubicSurface();
    BicubicSurface::PatchHint& hint = ownerp->updHint();
    const Real d = surf.calcDerivative(derivComponents, Vec2(p[0],p[1]), hint);
    return d;
}


