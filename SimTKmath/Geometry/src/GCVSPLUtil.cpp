/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

#include "simmath/internal/GCVSPLUtil.h"

namespace SimTK {

void GCVSPLUtil::gcvspl(const Vector& x, const Vector& y, const Vector& wx, Real wy, int degree, int md, Real val, Vector& c, Vector& wk, int& ier) {
    gcvspl(x, reinterpret_cast<const Vector_<Vec1>&>(y), wx, Vec1(wy), degree, md, val, reinterpret_cast<Vector_<Vec1>&>(c), wk,  ier);
}

Real GCVSPLUtil::splder(int derivOrder, int degree, Real t, const Vector& x, const Vector& coeff) {
    return splder(derivOrder, degree, t, x, reinterpret_cast<const Vector_<Vec1>&>(coeff))[0];
}

} // namespace SimTK
