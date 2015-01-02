#ifndef SimTK_SIMBODY_ASSEMBLY_CONDITION_QVALUE_H_
#define SimTK_SIMBODY_ASSEMBLY_CONDITION_QVALUE_H_

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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Assembler.h"
#include "simbody/internal/AssemblyCondition.h"

namespace SimTK {


//------------------------------------------------------------------------------
//                                 Q VALUE
//------------------------------------------------------------------------------
/** This AssemblyCondition requests that a particular generalized coordinate
end up with a specified value. You can use this as a goal or a constraint,
depending on how serious you are about this requirement. **/
class QValue : public AssemblyCondition {
public:
    /** Construct an assembly condition that requests that the specified
    generalized coordinate be brought to the indicated value. The value 
    can be changed subsequently using setValue(). **/
    QValue(MobilizedBodyIndex mbx, MobilizerQIndex qx,
           Real value)
    :   AssemblyCondition("QValue"), 
        mobodIndex(mbx), qIndex(qx), value(value) {}

    /** Return the currently set value to be used for this generalized
    coordinate. **/
    Real getValue() const {return value;}
    /** Change the value to be used for this generalized coordinate; this
    can be done repeatedly during tracking to follow changing requirements. **/
    void setValue(Real newValue) {value=newValue;}

    // For constraint:
    int getNumEquations(const State&) const {return 1;}
    int calcErrors(const State& state, Vector& error) const override {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const MobilizedBody& mobod = matter.getMobilizedBody(mobodIndex);
        error.resize(1);
        error[0] = mobod.getOneQ(state, qIndex) - value;
        return 0;
    }
    // Error jacobian is a zero-row except for a 1 in this q's entry (if
    // this q is free).
    int calcErrorJacobian(const State& state, Matrix& J) const override {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const MobilizedBody& mobod = matter.getMobilizedBody(mobodIndex);
        J.resize(1, getNumFreeQs());
        J = 0; // will have at most one non-zero

        // Find the FreeQIndex corresponding to this q.
        const QIndex thisIx = QIndex(mobod.getFirstQIndex(state)+qIndex);
        const Assembler::FreeQIndex thisFreeIx = getFreeQIndexOfQ(thisIx);

        // If this q isn't free then there is no way to affect the error
        // so the Jacobian stays all-zero.
        if (thisFreeIx.isValid())
            J(0,thisFreeIx) = 1;

        return 0;
    }

    // For goal: goal = (q-value)^2 / 2 (the /2 is for gradient beauty)
    int calcGoal(const State& state, Real& goal) const override {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const MobilizedBody& mobod = matter.getMobilizedBody(mobodIndex);
        goal = square(mobod.getOneQ(state, qIndex) - value) / 2;
        return 0;
    }
    // Return a gradient with only this q's entry non-zero (if
    // this q is free).
    int calcGoalGradient(const State& state, Vector& grad) const override {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const MobilizedBody& mobod = matter.getMobilizedBody(mobodIndex);
        grad.resize(getNumFreeQs());
        grad = 0; // will have at most one non-zero

        // Find the FreeQIndex corresponding to this q.
        const QIndex thisIx = QIndex(mobod.getFirstQIndex(state)+qIndex);
        const Assembler::FreeQIndex thisFreeIx = getFreeQIndexOfQ(thisIx);

        // If this q isn't free then there is no way to affect the goal
        // so the gradient stays all-zero.
        if (thisFreeIx.isValid())
            grad[thisFreeIx] = mobod.getOneQ(state, qIndex) - value;

        return 0;
    }

private:
    MobilizedBodyIndex mobodIndex;
    MobilizerQIndex    qIndex;
    Real               value;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_ASSEMBLY_CONDITION_QVALUE_H_
