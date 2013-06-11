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

#include "SimTKmath.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/LocalEnergyMinimizer.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include <vector>
#include <map>

using namespace SimTK;
using namespace std;

/**
 * This class defines the objective function which is passed to the Optimizer.
 */

class LocalEnergyMinimizer::OptimizerFunction : public OptimizerSystem {
public:
    // stateIn must be realized to Model stage already.
    OptimizerFunction(const MultibodySystem& system, const State& stateIn)
    :   OptimizerSystem(stateIn.getNQ()), system(system), state(stateIn) {
        state.updU() = 0;   // no velocities
        system.realize(state, Stage::Time); // we'll only change Position stage and above
        setNumEqualityConstraints(state.getNQErr());
    }
    int objectiveFunc(const Vector& parameters, bool new_parameters, Real& f) const {
        if (new_parameters)
            state.updQ() = parameters;
        system.realize(state, Stage::Dynamics);
        f = system.calcPotentialEnergy(state);
        return 0;
    }
    int gradientFunc(const Vector& parameters, bool new_parameters, Vector& gradient) const  {
        if (new_parameters)
            state.updQ() = parameters;
        system.realize(state, Stage::Dynamics);
        Vector_<SpatialVec> dEdR = system.getRigidBodyForces(state, Stage::Dynamics);
        const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
        Vector dEdU;
        // Convert spatial forces dEdR to generalized forces dEdU.
        matter.multiplyBySystemJacobianTranspose(state, dEdR, dEdU);
        dEdU -= system.getMobilityForces(state, Stage::Dynamics);
        matter.multiplyByNInv(state, true, -1*dEdU, gradient);
        return 0;
    }
    int constraintFunc(const Vector& parameters, bool new_parameters, Vector& constraints) const {
        state.updQ() = parameters;
        system.realize(state, Stage::Position);
        constraints = state.getQErr();
        return 0;
    }
    void optimize(Vector& q, Real tolerance) {
        Optimizer opt(*this);
        opt.useNumericalJacobian(true);
        opt.setConvergenceTolerance(tolerance);
        opt.optimize(q);
    }
private:
    const MultibodySystem& system;
    mutable State state;
};

void LocalEnergyMinimizer::minimizeEnergy(const MultibodySystem& system, State& state, Real tolerance) {
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    State tempState = state;
    matter.setUseEulerAngles(tempState, true);
    system.realizeModel(tempState);
    if (!matter.getUseEulerAngles(state))
        matter.convertToEulerAngles(state, tempState);
    OptimizerFunction optimizer(system, tempState);
    Vector q = tempState.getQ();
    optimizer.optimize(q, tolerance);
    if (matter.getUseEulerAngles(state))
        state.updQ() = q;
    else {
        tempState.updQ() = q;
        matter.convertToQuaternions(tempState, state);
    }
    system.realize(state, Stage::Dynamics);
}
