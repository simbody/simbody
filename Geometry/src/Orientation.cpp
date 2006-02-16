/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * Implementations of non-inline methods of classes dealing
 * with orientation.
 */

#include "simtk/SimTK.h"
#include "simmatrix/SmallMatrix.h"
#include "simbody/internal/Orientation.h"

#include <iostream>

namespace simtk {

// Calculate a rotation matrix R_GF which defines the F
// coordinate frame by rotating the G frame z axis into alignment 
// with the passed-in zF vector (that is, zF is the z axis
// of R_GF expressed in G. The result is not unique.
// TODO: (sherm) I think this can be done with one sqrt and no trig functions.
// Just create a vector pointing in a substantially different
// direction than zF (e.g., move its largest coord), then cross it with zF to get the first 
// perpendicular. Normalize that, cross again and you have a frame.
// Must be careful about signs to get a right-handed set.
// TODO: (sherm) Uh, shouldn't the 3rd column of the result just be zF?
RotationMat::RotationMat(const UnitVec<1>& zF) {
    // Use the individual measure numbers (i.e., projections of zF
    // onto the G coordinate axes) to calculate spherical coordinates,
    // considering an equatorial x-y plane with z North.
    const Real ct=zF[2];    // cos(theta)
    const Real theta = std::acos( ct );             // zenith (90-elevation)
    const Real psi   = std::atan2( zF[0] , zF[1] ); // 90-azimuth
    const Real st=std::sin(theta), cp=std::cos(psi), sp=std::sin(psi);

    // This is a space fixed 1-2-3 sequence with angles
    // a1=-theta, a2=0, a3=-psi. That is, to get from G to F
    // first rotate by -theta around the G frame x axis, 
    // then rotate by -psi around the G frame z axis.

    BaseMat::operator=(BaseMat( cp, ct*sp, sp*st,
                               -sp, ct*cp, cp*st,
                                0.,   -st,    ct));
}

std::ostream& operator<<(std::ostream& o, const RotationMat& m) {
    return o << m.asMat33();
}
std::ostream& operator<<(std::ostream& o, const UnitVec3& v) {
    return o << v.asVec3();
}

std::ostream& operator<<(std::ostream& o, const TransformMat& x) {
    return o << "{" << x.R() << x.T() << "}";
}

} // namespace simtk

