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
    OptimizerFunction(const MultibodySystem& system, const State& state) :
        OptimizerSystem(state.getNQ()), system(system), state(state) {
        setNumEqualityConstraints(state.getNQErr());
    }
    int objectiveFunc(const Vector& parameters, const bool new_parameters, Real& f) const {
        if (new_parameters)
            state.updQ() = parameters;
        system.realize(state, Stage::Dynamics);
        f = system.calcEnergy(state);
        return 0;
    }
    int gradientFunc(const Vector &parameters, const bool new_parameters, Vector &gradient) const  {
        if (new_parameters)
            state.updQ() = parameters;
        system.realize(state, Stage::Dynamics);
        Vector_<SpatialVec> dEdR = system.getRigidBodyForces(state, Stage::Dynamics);
        const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
        Vector dEdU;
        matter.calcInternalGradientFromSpatial(state, dEdR, dEdU);
        dEdU -= system.getMobilityForces(state, Stage::Dynamics);
        matter.multiplyByNInv(state, true, -1.0*dEdU, gradient);
        return 0;
    }
    int constraintFunc(const Vector& parameters, const bool new_parameters, Vector& constraints) const {
        state.updQ() = parameters;
        system.realize(state, Stage::Velocity);
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
}
