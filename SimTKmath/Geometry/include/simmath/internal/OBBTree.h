#ifndef SimTK_SIMMATH_OBB_TREE_H_
#define SimTK_SIMMATH_OBB_TREE_H_

/* -------------------------------------------------------------------------- *
 *                        SimTK Simbody: SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
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

/** @file
Defines an oriented bounding box tree with generalized leaves. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Box.h"
#include "simmath/internal/Geo_BicubicBezierPatch.h"

#include <cassert>

namespace SimTK {

//==============================================================================
//                               OBB LEAF
//==============================================================================
/** TODO **/
class OBBLeaf {
public:
    virtual ~OBBLeaf() {}
};

//==============================================================================
//                               OBB NODE
//==============================================================================
/** TODO **/
class SimTK_SIMMATH_EXPORT OBBNode {
public:
    OBBNode() : contents(0) {clear();}
    ~OBBNode() {clear();}

    void clear() {
        delete contents; contents=0;
        x0=y0=nx=ny=-1;
        children.clear();
    }

    bool isLeaf() const {return children.empty();}
    int getNumChildren() const {return (int)children.size();}
    const OBBNode& getChild(int i) const {return children[i];}
    OBBNode& updChild(int i) {return children[i];}

    // A box enclosing the contents.
    Geo::OrientedBox    box;
    int                 depth;  // 0 is root
    int                 height; // a leaf is 0, node is max of children+1

    // A cone enclosing the entire range of normals.
    UnitVec3            normal; // central normal
    Real                coneHalfAngle;  // 0<=a<=pi, pi/2 makes a halfspace

    // An arbitrary point on the contained surface, used in 
    // distance queries where distance to box is min distance, distance
    // to point is max distance.
    Vec3                pointOnSurface;

    int                 x0,y0; // Range of patches in this node
    int                 nx,ny;

    Array_<OBBNode>     children;

    // If no children, leaf contents:
    OBBLeaf*            contents; // non-null only for leaf (NOT USED YET)
    Vec2                centerUW; // (u,w) parameters of patch center
    Vec2                dims;     // half-u, half-w sizes
    Geo::BicubicBezierPatch patch; // TODO: no need to keep this around
};


//==============================================================================
//                               OBB TREE
//==============================================================================
/** TODO **/
class SimTK_SIMMATH_EXPORT OBBTree {
public:
    const OBBNode& getRoot() const {return root;}
    OBBNode& updRoot() {return root;}
private:
    OBBNode root;
};



} // namespace SimTK

#endif // SimTK_SIMMATH_OBB_TREE_H_
