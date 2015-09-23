/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-14 Stanford University and the Authors.        *
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

#include "SimTKmath.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/Assembler.h"
#include "simbody/internal/AssemblyCondition.h"
#include "simbody/internal/AssemblyCondition_Markers.h"
#include <map>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

//------------------------------------------------------------------------------
//                                  MARKERS
//------------------------------------------------------------------------------

Vec3 Markers::findCurrentMarkerLocation(MarkerIx mx) const {
    const SimbodyMatterSubsystem& matter = getMatterSubsystem();
    const Marker&                 marker = getMarker(mx);
    const MobilizedBody&          mobod  = matter.getMobilizedBody(marker.bodyB);
    const State&                  state  = getAssembler().getInternalState();
    const Transform&              X_GB   = mobod.getBodyTransform(state);
    return X_GB * marker.markerInB;
}

// goal = 1/2 sum( wi * ri^2 ) / sum(wi) for WRMS
int Markers::calcGoal(const State& state, Real& goal) const {
    const SimbodyMatterSubsystem& matter = getMatterSubsystem();
    goal = 0;
    // Loop over each body that has one or more active markers.
    Real wtot = 0;
    PerBodyMarkers::const_iterator bodyp = bodiesWithMarkers.begin();
    for (; bodyp != bodiesWithMarkers.end(); ++bodyp) {
        const MobilizedBodyIndex    mobodIx     = bodyp->first;
        const Array_<MarkerIx>&     bodyMarkers = bodyp->second;
        const MobilizedBody&        mobod = matter.getMobilizedBody(mobodIx);
        const Transform&            X_GB  = mobod.getBodyTransform(state);
        assert(bodyMarkers.size());
        // Loop over each marker on this body.
        for (unsigned m=0; m < bodyMarkers.size(); ++m) {
            const MarkerIx  mx = bodyMarkers[m];
            const Marker&   marker = markers[mx];
            assert(marker.bodyB == mobodIx); // better be on this body!
            const Vec3& location = observations[getObservationIxForMarker(mx)];
            if (location.isFinite()) { // skip NaNs
                goal += marker.weight
                        * (X_GB*marker.markerInB - location).normSqr();
                wtot += marker.weight;
            }
        }
    }

    goal /= (2*wtot);

    return 0;
}
// dgoal/dq = sum( wi * ri * dri/dq ) / sum(wi)
// This calculation is modeled after Peter Eastman's gradient implementation
// in ObservedPointFitter. It treats each marker position error as a potential
// energy function whose negative spatial gradient would be a spatial force F.
// We can then use Simbody's spatial force-to-generalized force method (using
// -F instead of F) to obtain the gradient in internal coordinates.
int Markers::calcGoalGradient(const State& state, Vector& gradient) const {
    const int np = getNumFreeQs();
    assert(gradient.size() == np);
    const SimbodyMatterSubsystem& matter = getMatterSubsystem();

    Vector_<SpatialVec> dEdR(matter.getNumBodies());
    dEdR = SpatialVec(Vec3(0), Vec3(0));
    // Loop over each body that has one or more active markers.
    Real wtot = 0;
    PerBodyMarkers::const_iterator bodyp = bodiesWithMarkers.begin();
    for (; bodyp != bodiesWithMarkers.end(); ++bodyp) {
        const MobilizedBodyIndex    mobodIx     = bodyp->first;
        const Array_<MarkerIx>&     bodyMarkers = bodyp->second;
        const MobilizedBody& mobod = matter.getMobilizedBody(mobodIx);
        const Transform& X_GB = mobod.getBodyTransform(state);
        assert(bodyMarkers.size());
        // Loop over each marker on this body.
        for (unsigned m=0; m < bodyMarkers.size(); ++m) {
            const MarkerIx  mx = bodyMarkers[m];
            const Marker&   marker = markers[mx];
            assert(marker.bodyB == mobodIx); // better be on this body!
            const Vec3& location = observations[getObservationIxForMarker(mx)];
            if (location.isFinite()) { // skip NaNs
                const Vec3 nf = marker.weight
                                * (X_GB*marker.markerInB - location);
                mobod.applyForceToBodyPoint(state, marker.markerInB, nf, dEdR);
                wtot += marker.weight;
            }
        }
    }
    // Convert spatial forces dEdR to generalized forces dEdU.
    Vector dEdU;
    matter.multiplyBySystemJacobianTranspose(state, dEdR, dEdU);

    dEdU /= wtot;

    const int nq = state.getNQ();
    if (np == nq) // gradient is full length
        matter.multiplyByNInv(state, true, dEdU, gradient);
    else { // calculate full gradient; extract the relevant parts
        Vector fullGradient(nq);
        matter.multiplyByNInv(state, true, dEdU, fullGradient);
        for (Assembler::FreeQIndex fx(0); fx < np; ++fx)
            gradient[fx] = fullGradient[getQIndexOfFreeQ(fx)];
    }


    return 0;
}

// TODO: We want the constraint version to minimize the same goal as above. But
// there can never be more than six independent constraints on the pose of
// a rigid body; this method should attempt to produce a minimal set so that
// the optimizer doesn't have to figure it out.
int Markers::calcErrors(const State& state, Vector& err) const
{   return AssemblyCondition::calcErrors(state,err); } //TODO


int Markers::calcErrorJacobian(const State& state, Matrix& jacobian) const
{   return AssemblyCondition::calcErrorJacobian(state,jacobian); } //TODO
int Markers::getNumErrors(const State& state) const
{   return AssemblyCondition::getNumErrors(state); } //TODO

// Run through all the Markers to find all the bodies that have at least one
// active marker. For each of those bodies, we collect all its markers so that
// we can process them all at once. Active markers are those whose weight is
// greater than zero. Also, if we haven't been given any observation<->marker
// correspondence, we're going to assume they map directly, with each
// ObservationIx the same as its MarkerIx.
int Markers::initializeCondition() const {
    // Fill in missing observation information if needed.
    if (observation2marker.empty()) {
        const Array_<MarkerIx> zeroLength; // gcc doesn't like this as a temp
        const_cast<Markers&>(*this).defineObservationOrder(zeroLength);
    }

    bodiesWithMarkers.clear();
    for (MarkerIx mx(0); mx < markers.size(); ++mx) {
        const Marker& marker = markers[mx];
        if (hasObservation(mx) && marker.weight > 0)
            bodiesWithMarkers[marker.bodyB].push_back(mx);
    }
    return 0;
}

// Throw away the bodiesWithMarkers map.
void Markers::uninitializeCondition() const {
    bodiesWithMarkers.clear();
}

