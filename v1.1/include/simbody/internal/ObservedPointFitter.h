#ifndef SimTK_SIMBODY_OBSERVED_POINT_FITTER_H_
#define SimTK_SIMBODY_OBSERVED_POINT_FITTER_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"

namespace SimTK {

/**
 * This class attempts to find the configuration of an internal coordinate model which best fits a set of
 * observed data.  The inputs to the algorithm are as follows:
 * 
 * - A MultibodySystem which describes the model to fit
 * - A set of points (called "stations") whose locations are defined relative to particular bodies
 * - The target location for each station, defined relative to ground
 * - (optional) A weight for each station, giving its relative importance for fitting
 * 
 * The output is a State giving the set of internal coordinates that best fit the stations to the target locations.
 */

class SimTK_SIMBODY_EXPORT ObservedPointFitter {
public:
    /**
     * Find the configuration of a MultibodySystem which best fits a set of target locations for stations.  This is identical to the other form of findBestFit(), but assumes every station
     * has a weight of 1.
     */

    static Real findBestFit(const MultibodySystem& system, State& state, const std::vector<MobilizedBodyIndex>& bodyIxs, const std::vector<std::vector<Vec3> >& stations, const std::vector<std::vector<Vec3> >& targetLocations, Real tolerance=0.001);

    /**
     * Find the configuration of a MultibodySystem which best fits a set of target locations for stations.
     * 
     * @param system      the MultibodySystem being analyzed
     * @param state       on exit, this State's Q vector contains the values which provide a best fit
     * @param bodyIxs     a list of MobilizedBodyIndexs corresponding to the bodies for which stations are defined
     * @param stations    the list of stations for each body.  stations[i][j] is the location of the j'th station for the body given by
     *                    bodyIxs[i], given in that body's reference frame.
     * @param targetLocations    the target locations for each body, given relative to ground.  targetLocations[i][j] is the target for stations[i][j].
     * @param weights     weights[i][j] is the weight to use for stations[i][j] when performing the fitting
     * @param tolerance   the distance tolerance within which the best fit should be found
     * @return the RMS distance of points in the best fit conformation from their target locations
     */

    static Real findBestFit(const MultibodySystem& system, State& state, const std::vector<MobilizedBodyIndex>& bodyIxs, const std::vector<std::vector<Vec3> >& stations, const std::vector<std::vector<Vec3> >& targetLocations, const std::vector<std::vector<Real> >& weights, Real tolerance=0.001);
private:
    static void createClonedSystem(const MultibodySystem& original, MultibodySystem& copy, const std::vector<MobilizedBodyIndex>& originalBodyIxs, std::vector<MobilizedBodyIndex>& copyBodyIxs);
    static void findUpstreamBodies(MobilizedBodyIndex currentBodyIx, const std::vector<int> numStations, const SimbodyMatterSubsystem& matter, std::vector<MobilizedBodyIndex>& bodyIxs, int requiredStations);
    static void findDownstreamBodies(MobilizedBodyIndex currentBodyIx, const std::vector<int> numStations, const std::vector<std::vector<MobilizedBodyIndex> > children, std::vector<MobilizedBodyIndex>& bodyIxs, int& requiredStations);
    static int findBodiesForClonedSystem(MobilizedBodyIndex primaryBodyIx, const std::vector<int> numStations, const SimbodyMatterSubsystem& matter, const std::vector<std::vector<MobilizedBodyIndex> > children, std::vector<MobilizedBodyIndex>& bodyIxs);
    class OptimizerFunction;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_OBSERVED_POINT_FITTER_H_
