#ifndef SimTK_SIMBODY_OBSERVED_POINT_FITTER_H_
#define SimTK_SIMBODY_OBSERVED_POINT_FITTER_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
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

    static Real findBestFit
       (const MultibodySystem&             system, 
        State&                             state, 
        const Array_<MobilizedBodyIndex>&  bodyIxs, 
        const Array_<Array_<Vec3> >&       stations, 
        const Array_<Array_<Vec3> >&       targetLocations, 
        Real                               tolerance=0.001);

    /** For compatibility with std::vector; requires extra copying. **/
    static Real findBestFit
       (const MultibodySystem&                  system, 
        State&                                  state, 
        const std::vector<MobilizedBodyIndex>&  bodyIxs, 
        const std::vector<std::vector<Vec3> >&  stations, 
        const std::vector<std::vector<Vec3> >&  targetLocations, 
        Real                                    tolerance=0.001) 
    {
        Array_<Array_<Vec3> > stationCopy(stations);
        Array_<Array_<Vec3> > targetCopy(targetLocations);
        return findBestFit(system,state,
                    ArrayViewConst_<MobilizedBodyIndex>(bodyIxs), // no copying here
                    stationCopy, targetCopy, tolerance);
    }

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

    static Real findBestFit
       (const MultibodySystem&             system, 
        State&                             state, 
        const Array_<MobilizedBodyIndex>&  bodyIxs, 
        const Array_<Array_<Vec3> >&       stations, 
        const Array_<Array_<Vec3> >&       targetLocations, 
        const Array_<Array_<Real> >&       weights, 
        Real                               tolerance=0.001);

    /** For compatibility with std::vector; requires extra copying. **/
    static Real findBestFit
       (const MultibodySystem&                  system, 
        State&                                  state, 
        const std::vector<MobilizedBodyIndex>&  bodyIxs, 
        const std::vector<std::vector<Vec3> >&  stations, 
        const std::vector<std::vector<Vec3> >&  targetLocations, 
        const std::vector<std::vector<Real> >&  weights, 
        Real                                    tolerance=0.001)
    {
        Array_<Array_<Vec3> > stationCopy(stations);
        Array_<Array_<Vec3> > targetCopy(targetLocations);
        Array_<Array_<Real> > weightCopy(weights);
        return findBestFit(system,state,
                    ArrayViewConst_<MobilizedBodyIndex>(bodyIxs), // no copying here
                    stationCopy, targetCopy, weightCopy,
                    tolerance);
    }


private:
    static void createClonedSystem(const MultibodySystem& original, MultibodySystem& copy, const Array_<MobilizedBodyIndex>& originalBodyIxs, Array_<MobilizedBodyIndex>& copyBodyIxs, bool& hasArtificialBaseBody);
    static void findUpstreamBodies(MobilizedBodyIndex currentBodyIx, const Array_<int> numStations, const SimbodyMatterSubsystem& matter, Array_<MobilizedBodyIndex>& bodyIxs, int requiredStations);
    static void findDownstreamBodies(MobilizedBodyIndex currentBodyIx, const Array_<int> numStations, const Array_<Array_<MobilizedBodyIndex> > children, Array_<MobilizedBodyIndex>& bodyIxs, int& requiredStations);
    static int findBodiesForClonedSystem(MobilizedBodyIndex primaryBodyIx, const Array_<int> numStations, const SimbodyMatterSubsystem& matter, const Array_<Array_<MobilizedBodyIndex> > children, Array_<MobilizedBodyIndex>& bodyIxs);
    class OptimizerFunction;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_OBSERVED_POINT_FITTER_H_
