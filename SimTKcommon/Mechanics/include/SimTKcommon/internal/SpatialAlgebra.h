#ifndef SimTK_SIMMATRIX_SPATIAL_ALGEBRA_H_
#define SimTK_SIMMATRIX_SPATIAL_ALGEBRA_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
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
These are declarations for special matrices and vectors of use in implementing
Rodriguez and Jain's Spatial Operator Algebra. **/

#include "SimTKcommon/SmallMatrix.h"

#include <iostream>

namespace SimTK {


/**@defgroup SpatialAlgebraUtilities    Spatial Algebra Utilities
   @ingroup GlobalFunctions

These utility functions are used for manipulation of spatial quantities that 
are contained in SpatialVec or SpatialMat objects. These are intended for
expert use and are mostly used in the implemention of friendlier methods such
as those in MobilizedBody that are used to obtain various spatial quantities.

@note Although we use SpatialVec for both, there are two different spatial 
vector bases: one for motion quantities like velocities, accelerations, and 
momentum and another for forces and impulses; be sure to use the appropriate 
functions. Also, we use a pair of ordinary vectors, following Abhi Jain,
rather than the similar but subtly different Plucker basis vectors used by 
Roy Featherstone.

Spatial vectors are used for combined (rotational,translational) quantities. 
These include
@verbatim
     spatial velocity     = (angularVelocity,linearVelocity)
     spatial acceleration = (angularAcceleration,linearAcceleration)
     spatial force        = (moment,force)
@endverbatim

Spatial configuration (pose) has to be handled differently though since
orientation is not a vector quantity. We use the Transform class for this 
concept, which includes an orientation matrix and a translation vector.
@see Transform **/
/**@{**/

/** SpatialVec[0] is the rotational component; [1] is translational. **/
typedef Vec<2,   Vec3>  SpatialVec;
/** This is the type of a transposed SpatialVec. **/
typedef Row<2,   Row3>  SpatialRow;
/** This is used for primarily for spatial mass properties. **/
typedef Mat<2,2, Mat33> SpatialMat;


// Pre-declare methods here so that we can list them in whatever order we'd
// like them to appear in Doxygen.
inline SpatialVec findRelativeVelocity( const Transform&  X_FA,
                                        const SpatialVec& V_FA,
                                        const Transform&  X_FB,
                                        const SpatialVec& V_FB);
inline SpatialVec findRelativeVelocityInF(  const Vec3&       p_AB_F,
                                            const SpatialVec& V_FA,
                                            const SpatialVec& V_FB);

inline SpatialVec findRelativeAcceleration( const Transform&  X_FA,
                                            const SpatialVec& V_FA,
                                            const SpatialVec& A_FA,
                                            const Transform&  X_FB,
                                            const SpatialVec& V_FB,
                                            const SpatialVec& A_FB);
inline SpatialVec findRelativeAccelerationInF(  const Vec3&       p_AB_F,
                                                const SpatialVec& V_FA,
                                                const SpatialVec& A_FA,
                                                const SpatialVec& V_FB,
                                                const SpatialVec& A_FB);

inline SpatialVec reverseRelativeVelocity(const Transform&  X_AB,
                                          const SpatialVec& V_AB);
inline SpatialVec reverseRelativeVelocityInA(const Transform&  X_AB,
                                             const SpatialVec& V_AB);

inline SpatialVec shiftVelocityBy(const SpatialVec& V_AB, const Vec3& r_A);
inline SpatialVec shiftVelocityFromTo(const SpatialVec& V_A_BP, 
                                      const Vec3&       fromP_A,
                                      const Vec3&       toQ_A);

inline SpatialVec shiftForceBy(const SpatialVec& F_AP, const Vec3& r_A);
inline SpatialVec shiftForceFromTo(const SpatialVec& F_AP, 
                                   const Vec3&       fromP_A,
                                   const Vec3&       toQ_A);



//==============================================================================
//                           FIND RELATIVE VELOCITY
//==============================================================================
/** @brief Find the relative spatial velocity between two frames A and B whose 
individual spatial velocities are known with respect to a third frame F, with
the result returned in A.

@param[in]      X_FA
    The pose of frame A measured and expressed in frame F.
@param[in]      V_FA
    The spatial velocity of frame A measured and expressed in frame F.
@param[in]      X_FB
    The pose of frame B measured and expressed in frame F.
@param[in]      V_FB
    The spatial velocity of frame B measured and expressed in frame F.
@return V_AB, the relative spatial velocity of frame B in frame A, expressed 
    in A.

Given the spatial velocity V_FA of frame A in a reference frame F, and
the spatial velocity V_FB of frame B in F, and transforms giving the poses of
frames A and B in F, calculate the relative velocity V_AB of frame B in 
frame A, measured and expressed in A. Typical usage:
@code
    Transform  X_GA, X_GB;      // assume these are known from somewhere
    SpatialVec V_GA, V_GB;

    SpatialVec V_AB = findRelativeVelocity(X_GA, V_GA, 
                                           X_GB, V_GB);
@endcode
@note This returns the result expressed in A which is almost always what you
want; however, if you don't want it in that frame you can save 30 flops by
calling findRelativeVelocityInF() instead.

Cost is 51 flops. @see findRelativeVelocityInF() **/
inline SpatialVec findRelativeVelocity(const Transform&  X_FA,
                                       const SpatialVec& V_FA,
                                       const Transform&  X_FB,
                                       const SpatialVec& V_FB)
{
    const Vec3 p_AB_F = X_FB.p() - X_FA.p();                    //  3 flops
    return ~X_FA.R()*findRelativeVelocityInF(p_AB_F,V_FA,V_FB); // 48 flops
}



//==============================================================================
//                         FIND RELATIVE VELOCITY IN F
//==============================================================================
/** @brief Find the relative spatial velocity between two frames A and B whose 
individual spatial velocities are known in a third frame F, but leave the 
result in F.

@param[in]      p_AB_F
    The vector from the A frame origin OA to the B frame origin OB, but
    expressed in frame F.
@param[in]      V_FA
    The spatial velocity of frame A measured and expressed in frame F.
@param[in]      V_FB
    The spatial velocity of frame B measured and expressed in frame F.
@return V_AB_F, the relative spatial velocity of frame B in frame A, but still 
    expressed in F.

Typically the relative velocity of B in A would be returned in A; most
users will want to use findRelativeVelocity() instead which returns the result
in A. Use of this method saves the substantial cost of reexpressing the result,
so is useful in the rare case that you don't want the final result in A.
Example:
@code
    Transform  X_GA, X_GB;      // assume these are known from somewhere
    SpatialVec V_GA, V_GB;

    const Vec3 p_AB_G = X_GB.p() - X_GA.p();
    SpatialVec V_AB_G = findRelativeVelocityInF(p_AB_G, V_GA, V_GB);
@endcode
Cost is 18 flops. @see findRelativeVelocity() **/
inline SpatialVec findRelativeVelocityInF(const Vec3&       p_AB_F,
                                          const SpatialVec& V_FA,
                                          const SpatialVec& V_FB)
{
    // Relative angular velocity of B in A, expressed in F.
    const Vec3 w_AB_F     = V_FB[0] - V_FA[0];              // 3 flops
    // Relative linear velocity of B in A, taken and expressed in F.
    const Vec3 p_AB_F_dot = V_FB[1] - V_FA[1];              // 3 flops
    // Get linear velocity taken in A by removing the component due
    // to A's rotation in F (still expressed in F).
    const Vec3 v_AB_F = p_AB_F_dot - V_FA[0] % p_AB_F;      // 12 flops

    return SpatialVec(w_AB_F, v_AB_F);
}



//==============================================================================
//                        FIND RELATIVE ACCELERATION
//==============================================================================
/** @brief Find the relative spatial acceleration between two frames A and B 
whose individual spatial accelerations are known with respect to a third 
frame F, with the result returned in A.

@param[in]      X_FA
    The pose of frame A measured and expressed in frame F.
@param[in]      V_FA
    The spatial velocity of frame A measured and expressed in frame F.
@param[in]      A_FA
    The spatial acceleration of frame A measured and expressed in frame F.
@param[in]      X_FB
    The pose of frame B measured and expressed in frame F.
@param[in]      V_FB
    The spatial velocity of frame B measured and expressed in frame F.
@param[in]      A_FB
    The spatial acceleration of frame B measured and expressed in frame F.
@return A_AB, the relative spatial acceleration of frame B in frame A, 
    expressed in A.

Given the spatial acceleration A_FA of frame A in a reference frame F, and
the spatial acceleration A_FB of frame B in F, and corresonding pose and
velcoity information, calculate the relative acceleration A_AB of frame B in 
frame A, measured and expressed in A. Typical usage:
@code
    Transform  X_GA, X_GB;      // assume these are known from somewhere
    SpatialVec V_GA, V_GB;
    SpatialVec A_GA, A_GB;

    SpatialVec A_AB = findRelativeAcceleration(X_GA, V_GA, A_GA,
                                               X_GB, V_GB, A_GB);
@endcode
@note This returns the result expressed in A which is almost always what you
want; however, if you don't want it in that frame you can save 30 flops by
calling findRelativeAccelerationInF() instead.

Cost is 105 flops. @see findRelativeAccelerationInF() **/
inline SpatialVec findRelativeAcceleration( const Transform&  X_FA,
                                            const SpatialVec& V_FA,
                                            const SpatialVec& A_FA,
                                            const Transform&  X_FB,
                                            const SpatialVec& V_FB,
                                            const SpatialVec& A_FB)
{
    const Vec3 p_AB_F = X_FB.p() - X_FA.p();                        //  3 flops
    return ~X_FA.R() *                                              // 30 flops
           findRelativeAccelerationInF(p_AB_F,V_FA,A_FA,V_FB,A_FB); // 72 flops
}



//==============================================================================
//                       FIND RELATIVE ACCELERATION IN F
//==============================================================================
/** @brief Find the relative spatial acceleration between two frames A and B 
whose individual spatial acceleration are known in a third frame F, but leave 
the result in F.

@param[in]      p_AB_F
    The vector from the A frame origin OA to the B frame origin OB, but
    expressed in frame F.
@param[in]      V_FA
    The spatial velocity of frame A measured and expressed in frame F.
@param[in]      A_FA
    The spatial acceleration of frame A measured and expressed in frame F.
@param[in]      V_FB
    The spatial velocity of frame B measured and expressed in frame F.
@param[in]      A_FB
    The spatial acceleration of frame B measured and expressed in frame F.
@return A_AB_F, the relative spatial acceleration of frame B in frame A, but 
    still expressed in F.

Typically the relative acceleration of B in A would be returned in A; most
users will want to use findRelativeAcceleration() instead which returns the 
result in A. Use of this method saves the substantial cost of reexpressing the
result, so is useful in the rare case that you don't want the final result in A.
Example:
@code
    Transform  X_GA, X_GB;      // assume these are known from somewhere
    SpatialVec V_GA, V_GB;
    SpatialVec A_GA, A_GB;

    const Vec3 p_AB_G = X_GB.p() - X_GA.p();
    SpatialVec V_AB_G = findRelativeAccelerationInF(p_AB_G, V_GA, A_GA,
                                                            V_GB, A_GB);
@endcode
Cost is 72 flops. @see findRelativeAcceleration() **/
inline SpatialVec findRelativeAccelerationInF(  const Vec3&       p_AB_F,
                                                const SpatialVec& V_FA,
                                                const SpatialVec& A_FA,
                                                const SpatialVec& V_FB,
                                                const SpatialVec& A_FB)
{
    const Vec3& w_FA = V_FA[0];     // aliases for convenience
    const Vec3& w_FB = V_FB[0];
    const Vec3& b_FA = A_FA[0];
    const Vec3& b_FB = A_FB[0];

    const Vec3 p_AB_F_dot    = V_FB[1] - V_FA[1]; // d/dt p taken in F   (3 flops)
    const Vec3 p_AB_F_dotdot = A_FB[1] - A_FA[1]; // d^2/dt^2 taken in F (3 flops)

    const Vec3 w_AB_F =     // relative angvel of B in A, exp. in F
        w_FB - w_FA;        // (3 flops)
    const Vec3 v_AB_F =              // d/dt p taken in A, exp in F
        p_AB_F_dot - w_FA % p_AB_F;  // (12 flops)

    const Vec3 w_AB_F_dot = b_FB - b_FA; // d/dt of w_AB_F taken in F (3 flops)
    const Vec3 v_AB_F_dot =              // d/dt v_AB_F taken in F
        p_AB_F_dotdot - (b_FA % p_AB_F + w_FA % p_AB_F_dot); // (24 flops)
    
    // We have the derivative in F; change it to derivative in A by adding in 
    // contribution caused by motion of F in A, that is w_AF X w_AB_F. (Note 
    // that w_AF=-w_FA.)
    const Vec3 b_AB_F =             // ang. accel. of B in A, exp. in F
        w_AB_F_dot - w_FA % w_AB_F; // (12 flops)
    const Vec3 a_AB_F =             // taken in A, exp. in F
        v_AB_F_dot - w_FA % v_AB_F; // (12 flops)

    return SpatialVec(b_AB_F, a_AB_F); // taken in A, expressed in F
}



//==============================================================================
//                         REVERSE RELATIVE VELOCITY
//==============================================================================
/** @brief Given the relative velocity of frame B in frame A, reverse that to
give the relative velocity of frame A in B.

@param[in]      X_AB
    The pose of frame B in frame A, measured and expressed in A.
@param[in]      V_AB
    The relative spatial velocity of frame B in frame A, measured and 
    expressed in frame A.
@return V_BA, the relative spatial velocity of frame A in frame B, measured
    and expressed in B.

The input is expressed in the A frame; the result will be expressed in the B
frame instead. If you prefer that the result remain in the A frame you should 
call reverseRelativeVelocityInA() instead to avoid the extra cost of changing 
frames. Example:
@code
    Transform  X_AB;    // assume these are known from somewhere
    SpatialVec V_AB;

    SpatialVec V_BA = reverseRelativeVelocity(X_AB, V_AB);
@endcode
@note If the frame origins were in the same spatial location, then the result 
would just be the negative of the supplied velocity. However, since the linear
component of spatial velocity has to be measured at a point, and we're 
switching from measuring at a point coincident with B's origin OB to one
coincident with A's origin OA, there is going to be a change in the linear
part of the result. The angular velocity will just be negated, though, and
then reexpressed in B.

Cost is 51 flops. @see reverseRelativeVelocityInA() **/
inline SpatialVec reverseRelativeVelocity(const Transform&  X_AB,
                                          const SpatialVec& V_AB)
{
    // Reverse the velocity but with the result still expressed in A.
    const SpatialVec V_BA_A = reverseRelativeVelocityInA(X_AB,V_AB); 
                                                                // 21 flops
    // Then reexpress in B.
    return ~X_AB.R()*V_BA_A;                                    // 30 flops
}



//==============================================================================
//                       REVERSE RELATIVE VELOCITY IN A
//==============================================================================
/** @brief Given the relative velocity of frame B in frame A, reverse that to
give the relative velocity of frame A in B, but leave the result expressed
in frame A.

@param[in]      X_AB
    The pose of frame B in frame A, measured and expressed in A.
@param[in]      V_AB
    The relative spatial velocity of frame B in frame A, measured and 
    expressed in frame A.
@return V_BA_A, the relative velocity of frame A in frame B, but still 
    expressed in B.

The input V_AB is expressed in the A frame; you will almost always want the
output V_BA expressed in the B frame which is what the function 
reverseRelativeVelocity() does. However, if you're going to want it in
some other frame ultimately you may prefer to avoid the substantial cost of 
reexpressing it in B now, in which case this routine is useful.

See reverseRelativeVelocity() for more information about what this does. 
Example:
@code
    Transform  X_AB;    // assume these are known from somewhere
    SpatialVec V_AB;    // (expressed in A)

    // result is still expressed in A
    SpatialVec V_BA_A = reverseRelativeVelocityInA(X_AB, V_AB);
@endcode

Cost is 21 flops. @see reverseRelativeVelocity() **/
inline SpatialVec reverseRelativeVelocityInA(const Transform&  X_AB,
                                             const SpatialVec& V_AB)
{
    // Change the measurement point from a point coincident with OB
    // to a point coincident with OA, and negate since we want A's velocity
    // in B rather than the other way around.
    const SpatialVec V_BA_A = -shiftVelocityBy(V_AB, -X_AB.p()); // 21 flops
    return V_BA_A;
}



//==============================================================================
//                           SHIFT VELOCITY BY
//==============================================================================
/** @brief Shift a relative spatial velocity measured at some point to that
same relative spatial quantity but measured at a new point given by an offset
from the old one.

@param[in]      V_AB
    The relative spatial velocity of frame B in frame A, measured and 
    expressed in frame A.
@param[in]      r_A
    The vector offset, expressed in frame A, by which to change the point at 
    which the translational component of the relative spatial velocity is 
    measured.
@return V_A_BQ, the relative velocity of frame B in frame A, but measured at
    the point Q=Bo+r rather than at B's origin Bo.

Given the spatial velocity V_AB of frame B in A, measured at a point
coincident with B's origin Bo, change it to the spatial velocity V_A_BQ 
representing the same relationship but with the velocity measured at a new 
point Q=Bo+r for some position vector r. All vectors are measured and expressed
in frame A, including the vector r. Example:
@code
    SpatialVec V_AB;     // assume these are known from somewhere
    Vec3       offset_A; // Q = Bo + offset

    SpatialVec V_A_BQ = shiftVelocityBy(V_AB, offset_A);
@endcode

@note The shift in location leaves the relative angular velocity w the same but
results in the linear velocity changing by w X r.

Cost is 12 flops. @see shiftVelocityFromTo() **/
inline SpatialVec shiftVelocityBy(const SpatialVec& V_AB, const Vec3& r_A)
{   return SpatialVec( V_AB[0], V_AB[1] + V_AB[0] % r_A ); } // vp=v + wXr


//==============================================================================
//                          SHIFT VELOCITY FROM TO
//==============================================================================
/** @brief Shift a relative spatial velocity measured at some point P to that
same relative spatial quantity but measured at a new point Q given the points
P and Q.

@param[in]      V_A_BP
    The relative spatial velocity of frame B in frame A, measured and 
    expressed in frame A, with the linear component measured at a point P.
@param[in]      fromP_A
    The "from" point P at which the input linear velocity was measured, given
    as a vector from A's origin OA to the point P, expressed in A.
@param[in]      toQ_A
    The "to" point Q at which we want to re-measure the linear velocity, given
    as a vector from A's origin OA to the point Q, expressed in A.
@return V_A_BQ, the relative velocity of frame B in frame A, but measured at
    the point Q rather than at point P.

Given the spatial velocity V_A_BP of frame B in A, measured at a point P,
change it to the spatial velocity V_A_BQ representing the same relationship but
with the velocity measured at a new point Q. Example:
@code
    // assume these are known from somewhere
    Transform  X_AB;    // contains the vector from OA to OB  
    SpatialVec V_AB;    // linear velocity is measured at origin OB of B
    Vec3       p_AQ;    // vector from OA to some other point Q, in A

    SpatialVec V_A_BQ = shiftVelocityFromTo(V_AB, X_AB.p(), p_AQ);
@endcode

@note There is no way to know whether the supplied velocity was actually
measured at P; this method really just shifts the relative velocity by
the vector r=(to-from). Use it carefully.

Cost is 15 flops. @see shiftVelocityBy() **/
inline SpatialVec shiftVelocityFromTo(const SpatialVec& V_A_BP, 
                                      const Vec3&       fromP_A,
                                      const Vec3&       toQ_A)
{   return shiftVelocityBy(V_A_BP, toQ_A - fromP_A); }



//==============================================================================
//                         SHIFT ACCELERATION BY
//==============================================================================
/** @brief Shift a relative spatial acceleration measured at some point to that
same relative spatial quantity but measured at a new point given by an offset
from the old one.

@param[in]      A_AB
    The relative spatial acceleration of frame B in frame A, measured and 
    expressed in frame A.
@param[in]      w_AB
    The relative angular velocity of frame B in frame A, expressed in frame A.
@param[in]      r_A
    The vector offset, expressed in frame A, by which to change the point at 
    which the translational component of the relative spatial acceleration is 
    measured.
@return A_A_BQ, the relative acceleration of frame B in frame A, but measured at
    the point Q=Bo+r rather than at B's origin Bo.

Given the spatial acceleration A_AB and angular velocity w_AB of frame B in A, 
measured at a point coincident with B's origin Bo, change it to the spatial 
acceleration A_A_BQ representing the same relationship but with the acceleration
measured at a new point Q=Bo+r for some position vector r. All vectors are 
measured and expressed in frame A, including the vector r. Example:
@code
    SpatialVec A_AB;     // assume these are known from somewhere
    Vec3       w_AB;
    Vec3       offset_A; // Q = Bo + offset

    SpatialVec A_A_BQ = shiftAccelerationBy(A_AB, w_AB, offset_A);
@endcode

@note The shift in location leaves the relative angular acceleration b the same
but results in the linear acceleration changing by b X r + w X (w X r).

Cost is 33 flops. @see shiftAccelerationFromTo() **/
inline SpatialVec shiftAccelerationBy(const SpatialVec& A_AB, 
                                      const Vec3&       w_AB, 
                                      const Vec3&       r_A)
{   return SpatialVec( A_AB[0],   
                       A_AB[1] + A_AB[0] % r_A  + w_AB % (w_AB % r_A) ); } 



//==============================================================================
//                        SHIFT ACCELERATION FROM TO
//==============================================================================
/** @brief Shift a relative spatial acceleration measured at some point P to 
that same relative spatial quantity but measured at a new point Q given the 
points P and Q.

@param[in]      A_A_BP
    The relative spatial acceleration of frame B in frame A, measured and 
    expressed in frame A, with the linear component measured at a point P.
@param[in]      w_AB
    The relative angular velocity of frame B in frame A, expressed in frame A.
@param[in]      fromP_A
    The "from" point P at which the input linear acceleration was
    measured, given as a vector from A's origin Ao to the point P, 
    expressed in A.
@param[in]      toQ_A
    The "to" point Q at which we want to re-measure the linear acceleration, 
    given as a vector from A's origin Ao to the point Q, expressed 
    in A.
@return A_A_BQ, the relative acceleration of frame B in frame A, but measured at
    the point Q rather than at point P.

Given the spatial acceleration A_A_BP of frame B in A, measured at a point P,
change it to the spatial acceleration A_A_BQ representing the same relationship 
but with the acceleration measured at a new point Q. Example:
@code
    // assume these are known from somewhere
    Transform  X_AB;    // contains the vector from Ao to Bo  
    SpatialVec A_AB;    // linear acceleration is measured at origin Bo of B
    Vec3       w_AB;
    Vec3       p_AQ;    // vector from Ao to some other point Q, in A

    SpatialVec A_A_BQ = shiftAccelerationFromTo(A_AB, w_AB, X_AB.p(), p_AQ);
@endcode

@note There is no way to know whether the supplied acceleration was
actually measured at P; this method really just shifts the relative 
acceleration by the vector r=(to-from). Use it carefully.

Cost is 36 flops. @see shiftAccelerationBy() **/
inline SpatialVec shiftAccelerationFromTo(const SpatialVec& A_A_BP, 
                                          const Vec3&       w_AB,
                                          const Vec3&       fromP_A,
                                          const Vec3&       toQ_A)
{   return shiftAccelerationBy(A_A_BP, w_AB, toQ_A - fromP_A); }



//==============================================================================
//                              SHIFT FORCE BY
//==============================================================================
/** @brief Shift a spatial force applied at some point of a body to that
same spatial force applied at a new point given by an offset
from the old one.

@param[in]      F_AP
    A spatial force (moment and linear force), expressed in the A frame,
    whose translational component is applied at a point P.
@param[in]      r_A
    The vector offset, expressed in frame A, by which to change the point at 
    which the translational component of the input force is to be applied.
@return F_AQ, the same physical effect as the input but with the moment
    adjusted to reflect force application at point Q=P+r rather than at the
    original point P.

Given the spatial force F_AP including a pure moment m and a force vector f
applied at a point P, return the equivalent force F_AQ representing the same
physical quantity but as though the force were applied at a point Q=P+r
for some position vector r. All vectors are expressed in frame A. Example: 
@code
    SpatialVec F_AP;     // assume these are known from somewhere
    Vec3       offset_A; // Q = P + offset

    SpatialVec F_AQ = shiftForceBy(F_AP, offset_A);
@endcode

@note The shift in location leaves the force f the same but
results in an adjustment to the moment of -(r X f).

Cost is 12 flops. @see shiftForceFromTo() **/
inline SpatialVec shiftForceBy(const SpatialVec& F_AP, const Vec3& r_A)
{   return SpatialVec(F_AP[0] -  r_A % F_AP[1], F_AP[1]); } // mq = mp - r X f



//==============================================================================
//                           SHIFT FORCE FROM TO
//==============================================================================
/** @brief Shift a spatial force applied at some point P of a body to that
same spatial force applied at a new point Q, given P and Q.

@param[in]      F_AP
    A spatial force (moment and linear force), expressed in the A frame, whose 
    translational component is applied at a point P.
@param[in]      fromP_A
    The "from" point P at which the input force is applied, given as a 
    vector from A's origin OA to the point P, expressed in A.
@param[in]      toQ_A
    The "to" point Q to which we want to move the force application point, 
    given as a vector from A's origin OA to the point Q, expressed in A.
@return F_AQ, the same physical effect as the input but with the moment
    adjusted to reflect force application at point Q rather than at the 
    original point P.

Given the spatial force F_AP including a pure moment m and a force vector f
applied at a point P, return the equivalent force F_AQ representing the same
physical quantity but as though the force were applied at a new point Q. All 
vectors are expressed in frame A and points are measured from A's origin OA.
Example:
@code
    // assume these are known from somewhere
    SpatialVec F_AP;    // linear force is applied at point P
    Vec3       p_AP;    // vector from OA to P, in A
    Vec3       p_AQ;    // vector from OA to some other point Q, in A

    SpatialVec F_AQ = shiftForceFromTo(F_AP, p_AP, p_AQ);
@endcode

@note There is no way to know whether the supplied force was actually
applied at P; this method really just shifts the application point by
the vector r=(to-from). Use it carefully.

Cost is 15 flops. @see shiftForceBy() **/
inline SpatialVec shiftForceFromTo(const SpatialVec& F_AP, 
                                   const Vec3&       fromP_A,
                                   const Vec3&       toQ_A)
{   return shiftForceBy(F_AP, toQ_A - fromP_A); }

/**@}**/



//==============================================================================
//                                  PHI MATRIX
//==============================================================================
// support for efficient matrix multiplication involving the special phi
// matrix

class PhiMatrixTranspose;

class PhiMatrix {
public:
    typedef PhiMatrixTranspose TransposeType;

    PhiMatrix() { setToNaN(); }
    PhiMatrix(const Vec3& l) : l_(l) {}

    void setToZero() { l_ = 0; }
    void setToNaN()  { l_.setToNaN(); }

    SpatialMat toSpatialMat() const {
        return SpatialMat(Mat33(1), crossMat(l_),
                          Mat33(0),   Mat33(1));
    }

    const Vec3& l() const { return l_; }
private:
    Vec3 l_;
};

class PhiMatrixTranspose {
public:
    PhiMatrixTranspose(const PhiMatrix& phi) : phi(phi) {}

    SpatialMat toSpatialMat() const {
        return SpatialMat(   Mat33(1)    , Mat33(0),
                          crossMat(-l()) , Mat33(1));
    }

    const Vec3& l() const {return phi.l();}
private:
  const PhiMatrix& phi;
};

inline PhiMatrixTranspose
transpose(const PhiMatrix& phi)
{
    PhiMatrixTranspose ret(phi);
    return ret;
}

inline PhiMatrixTranspose
operator~(const PhiMatrix& phi) {return transpose(phi);}

inline SpatialVec
operator*(const PhiMatrix&  phi,
          const SpatialVec& v)
{
    return SpatialVec(v[0] + phi.l() % v[1], // 12 flops
                      v[1]);
}

inline SpatialMat
operator*(const PhiMatrix&  phi,
          const SpatialMat& m)
{
    const Mat33 x = crossMat(phi.l());  // 3 flops
    return SpatialMat( m(0,0) + x*m(1,0), m(0,1) + x*m(1,1), // 108 flops
                           m(1,0)       ,     m(1,1));
}

inline SpatialVec
operator*(const PhiMatrixTranspose& phiT,
          const SpatialVec&         v)
{
    return SpatialVec(v[0],
                      v[1] + v[0] % phiT.l());  // 12 flops
}


inline SpatialMat
operator*(const SpatialMat::THerm&  m,
          const PhiMatrixTranspose& phiT)
{
    const Mat33 x = crossMat(phiT.l()); // 3 flops
    return SpatialMat( m(0,0) - m(0,1) * x, m(0,1),     // 54 flops
                       m(1,0) - m(1,1) * x, m(1,1) );   // 54 flops
}

inline SpatialMat
operator*(const SpatialMat&         m,
          const PhiMatrixTranspose& phiT)
{
    const Mat33 x = crossMat(phiT.l()); // 3 flops
    return SpatialMat( m(0,0) - m(0,1) * x, m(0,1),     // 54 flops
                       m(1,0) - m(1,1) * x, m(1,1) );   // 54 flops
}

inline bool
operator==(const PhiMatrix& p1, const PhiMatrix& p2)
{
    return p1.l() == p2.l();
}

inline bool
operator==(const PhiMatrixTranspose& p1, const PhiMatrixTranspose& p2)
{
    return p1.l() == p2.l();
}
} // namespace SimTK

#endif // SimTK_SIMMATRIX_SPATIAL_ALGEBRA_H_
