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

/** @file
Non-inline static methods from the Geo::Box_ class. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Box.h"
#include "simmath/internal/Geo_Sphere.h"

#include <cstdio>
#include <iostream>
using std::cout; using std::endl;

namespace SimTK {

//==============================================================================
//                               GEO :: BOX
//==============================================================================

// For testing whether this box intersects a given oriented box, we use the
// separating axis method first described in Gottschalk, S., Lin, MC,
// Manocha, D, "OBBTree: a hierarchical structure for rapid interference
// detection." Proceedings of the 23rd Annual Conference on Computer
// Graphics and Interactive Techniques, pp. 171-180, 1996. This is a series
// of 15 tests used to determine whether there exists a separating plane
// between the boxes. The possible plane normals are the 3 faces of each box,
// and the 9 possible cross products of the three edge directions from each
// box. In each case we project the boxes onto the plane normal lines
// and see whether the projections overlap. If they don't, we have found a
// separating plane and can return right away. This can be sped up by doing
// only the first six face tests and accepting the occasional false positive
// (allegedly less than 10%); we provide a method for doing that abbreviated
// version of the test as well as the full method.

// If the boxes have a parallel edge, they have to be handled specially because
// the cross product of those edges will be zero. Then the projections of both
// the box dimensions and the center-to-center distance will be zero, creating
// a false appearance of having found a separating plane. Christer Ericson
// advocates use of a "fudge factor" to ensure a conservative result, but David
// Eberly tackles it more directly by noting that if any edge is parallel then
// one of the six faces *must* be the separating plane if there is one. So if
// you've done six tests and not found a separating plane there is no need to
// check the edge cross products; the boxes do intersect. We follow Eberly's
// approach.



// First check only the B and O face normals as possible separating planes. If
// we find one of those works the boxes can't be intersecting, otherwise they
// might be. Ordering for speed:
//  - first check faces of this box B, because we're working in B frame
//    making these tests cheaper
//  - when working on a box, test smallest dimension first because that
//    covers the most space in one test
// For a few extra flops (1 per axis), we can also spot whether one box center
// is completely inside the other and if so report early that the boxes are
// definitely intersecting. Worst case here is about 74 flops. Minimum is 16
// flops if the first axis checked is a separating axis.

// Private helper method returns 0 if a separating plane is found (definitely
// no intersection), 2 if the centers overlap (definitely an intersection)
// and 1 otherwise (probably an intersection, but maybe not).

// It is convenient to calculate the absolute value of the Transform matrix
// because those are reused. We'll do that here and return the results in case
// the caller needs to use them. Note that the absolute value of a Rotation
// matrix is not necessarily a Rotation since the resulting vectors might not
// even be perpendicular.
template <class P> int Geo::Box_<P>::
intersectsOrientedBoxHelper(const OrientedBox_<P>&  O,
                            Mat33P&                 absR_BO,
                            Vec3P&                  absP_BO) const {
    const RotationP& R_BO = O.getTransform().R();
    const Vec3P&     p_BO = O.getTransform().p(); // center location

    absR_BO = R_BO.asMat33().abs(); // about 6 flops
    absP_BO = p_BO.abs();           // about 2 flops
    const Geo::Box_<P>& ob = O.getBox();
    const Vec3P&        oh = ob.getHalfLengths();

    bool centerOisInsideB = true; // set false when we see a counterexample

    // Try each face of B in increasing order of box dimension so that we
    // start with the skinniest slab of box B. Cost is 8 flops per axis.
    // "V" below is the 1st quadrant (+++) vertex of each box.
    for (int i=0; i<3; ++i) {
        const CoordinateAxis b = getOrderedAxis(i);
        const RealP rb = h[b];                  // proj. p_BV onto axis b
        const RealP ro = dot(oh, absR_BO[b]);   // proj. p_OV onto axis b (in O)
        const RealP d  = absP_BO[b];            // ctr-ctr distance along b
        if (d > rb+ro)
            return 0; // axis b is a separating plane normal
        if (d > rb)
            centerOisInsideB = false; // saw a counterexample
    }

    if (centerOisInsideB)
        return 2; // definitely does intersect!

    bool centerBisInsideO = true;

    // Try each face of O in increasing order of box dimension so that we
    // start with the skinniest slab of box O. Cost is 14 flops per axis.
    for (int i=0; i<3; ++i) {
        const CoordinateAxis o = ob.getOrderedAxis(i);
        // Note: round brackets indexing selects a column of R_BO1.
        const RealP rb = dot(h, absR_BO(o));    // proj. p_BV onto axis o (in B)
        const RealP ro = oh[o];                 // proj. of p_OV onto axis o
        const RealP d  = std::abs(dot(p_BO, R_BO(o))); //B-O dist along o (in B)
        if (d > rb+ro)
            return 0; // axis o is a separating plane normal
        if (d > ro)
            centerBisInsideO = false; // saw a counterexample
    }

    if (centerBisInsideO)
        return 2; // definitely does intersect!

    return 1; // probably intersect
}

// Worst case 74 flops.
template <class P> bool Geo::Box_<P>::
mayIntersectOrientedBox(const Geo::OrientedBox_<P>& ob) const {
    Mat33P absR_BO; Vec3P absP_BO;
    const int result = intersectsOrientedBoxHelper(ob, absR_BO, absP_BO);
    return result != 0;
}

// Worst case 74 + 9 + 9*12 = 191 flops.
template <class P> bool Geo::Box_<P>::
intersectsOrientedBox(const Geo::OrientedBox_<P>& ob) const {
    Mat33P absR_BO;
    Vec3P  absP_BO;
    const int result = intersectsOrientedBoxHelper(ob, absR_BO, absP_BO);
    if (result == 0) return false;
    else if (result == 2) return true;

    // Indeterminate -- we still don't know if there is an intersection.

    const RotationP& R_BO = ob.getTransform().R();
    const Vec3P&     p_BO = ob.getTransform().p(); // center location
    const Vec3P&     oh   = ob.getHalfLengths();

    // First we have to determine whether the boxes have a pair of edges that
    // are (almost) parallel. If so, one of the cross product tests below will
    // produce a near-zero separating vector candidate. Projections onto that
    // line will yield zero lengths and we'll end up comparing noise to noise
    // to see if we found a separating vector, producing unacceptable false
    // negatives. On the other hand, if there is a pair of parallel edges, then
    // one of the face normals we just tested would have to be a separating
    // plane if there is one (draw a picture). In that case we know that we
    // are intersecting and can return true now. This also serves as a fast
    // out for the common case of parallel boxes.

    // If there is a parallel edge, one of the entries for O's axes in |R_BO|
    // will be within noise of one of B's axes, that is, 100, 010, or 001.
    // David Eberly's algorithm looks for a "1" as cos(near zero) but that's a
    // bad idea because sin() and cos() are very flat near 1. It is better to
    // look for "0", as sin(near zero) because sin(theta)~=theta for small
    // angles. So we pick an angle tolerance that we declare is small enough to
    // be considered parallel, and then see if a column of |R_BO| has two
    // entries that size or smaller. It costs 9 flops to do this test.

    const RealP tol = Geo::getDefaultTol<P>(); // TODO: too tight?
    for (int i=0; i<3; ++i) {
        int numZeroes=0;
        const Vec3P& v = absR_BO(i); // get column; row would work too
        if (v[0] < tol) ++numZeroes; // 3 flops
        if (v[1] < tol) ++numZeroes;
        if (v[2] < tol) ++numZeroes;
        assert(numZeroes < 3); // can't happen in a Rotation!
        if (numZeroes==2)
            return true; // boxes intersect
    }

    // There are no parallel edges. Check the 9 remaining possibilites for
    // separating axes formed by the cross product of an edge from B and an
    // edge from O. B's edges are 100,010,001 so the cross products simplify
    // considerably: 100 X [o0 o1 o2] = [  0 -o2  o1]
    //               010 X [o0 o1 o2] = [ o2   0 -o0]
    //               001 X [o0 o1 o2] = [-o1  o0   0]
    // I don't know of any particularly good ordering for this step.

    RealP rb, ro; // Projections of box dimensions on separating vector.
    RealP d;      // Projection of center-center distance on separating vector.

    const Mat33P&   R    = R_BO.asMat33();       // abbreviations
    const Vec3P&    p    = p_BO;
    const Mat33P&   Rabs = absR_BO;

    // Each of these tests takes about 12 flops.
    rb = h[1] *Rabs(2, 0)+h[2] *Rabs(1, 0);         // b0 X oO
    ro = oh[1]*Rabs(0, 2)+oh[2]*Rabs(0, 1);
    d  = std::abs(p[2]*R(1, 0) - p[1]*R(2, 0));
    if (d > rb+ro) return false;

    rb = h[1] *Rabs(2, 1)+h[2] *Rabs(1, 1);         // b0 X o1
    ro = oh[0]*Rabs(0, 2)+oh[2]*Rabs(0, 0);
    d  = std::abs(p[2]*R(1, 1) - p[1]*R(2, 1));
    if (d > rb+ro) return false;

    rb = h[1] *Rabs(2, 2)+h[2] *Rabs(1, 2);         // b0 X o2
    ro = oh[0]*Rabs(0, 1)+oh[1]*Rabs(0, 0);
    d  = std::abs(p[2]*R(1, 2) - p[1]*R(2, 2));
    if (d > rb+ro) return false;

    rb = h[0] *Rabs(2, 0)+h[2] *Rabs(0, 0);         // b1 X oO
    ro = oh[1]*Rabs(1, 2)+oh[2]*Rabs(1, 1);
    d  = std::abs(p[0]*R(2, 0) - p[2]*R(0, 0));
    if (d > rb+ro) return false;

    rb = h[0] *Rabs(2, 1)+h[2] *Rabs(0, 1);         // b1 X o1
    ro = oh[0]*Rabs(1, 2)+oh[2]*Rabs(1, 0);
    d  = std::abs(p[0]*R(2, 1) - p[2]*R(0, 1));
    if (d > rb+ro) return false;

    rb = h[0] *Rabs(2, 2)+h[2] *Rabs(0, 2);         // b1 X o2
    ro = oh[0]*Rabs(1, 1)+oh[1]*Rabs(1, 0);
    d  = std::abs(p[0]*R(2, 2) - p[2]*R(0, 2));
    if (d > rb+ro) return false;

    rb = h[0] *Rabs(1, 0)+h[1] *Rabs(0, 0);         // b2 X oO
    ro = oh[1]*Rabs(2, 2)+oh[2]*Rabs(2, 1);
    d  = std::abs(p[1]*R(0, 0) - p[0]*R(1, 0));
    if (d > rb+ro) return false;

    rb = h[0] *Rabs(1, 1)+h[1] *Rabs(0, 1);         // b2 X o1
    ro = oh[0]*Rabs(2, 2)+oh[2]*Rabs(2, 0);
    d  = std::abs(p[1]*R(0, 1) - p[0]*R(1, 1));
    if (d > rb+ro) return false;

    rb = h[0] *Rabs(1, 2)+h[1] *Rabs(0, 2);         // b2 X o2
    ro = oh[0]*Rabs(2, 1)+oh[1]*Rabs(2, 0);
    d  = std::abs(p[1]*R(0, 2) - p[0]*R(1, 2));
    if (d > rb+ro) return false;

    return true; // The boxes intersect.
}



// Explicit instantiations for float and double.
template class Geo::Box_<float>;
template class Geo::Box_<double>;

//==============================================================================
//                           GEO :: ALIGNED BOX
//==============================================================================
// Explicit instantiations for float and double.
template class Geo::AlignedBox_<float>;
template class Geo::AlignedBox_<double>;

//==============================================================================
//                           GEO :: ORIENTED BOX
//==============================================================================



// Explicit instantiations for float and double.
template class Geo::OrientedBox_<float>;
template class Geo::OrientedBox_<double>;

}  // End of namespace SimTK
