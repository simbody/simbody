/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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


/**@file
 *
 * Implementation of non-inline methods from the EventWitness classes.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Event.h"
#include "SimTKcommon/internal/EventTrigger.h"
#include "SimTKcommon/internal/EventWitness.h"

using namespace SimTK;


/*static*/ std::string EventWitness::toString(Range range) {
    switch(range) {
    case Bilateral:         return "Bilateral";
    case Unilateral:        return "Unilateral";
    default: return "BAD WITNESS Range " + std::to_string((int)range);
    }
}


/*static*/ std::string EventWitness::toString(Direction direction) {
    switch(direction) {
    case Rising:            return "Rising";
    case Falling:           return "Falling";
    case RisingAndFalling:  return "RisingAndFalling";
    default: return "BAD WITNESS Direction " + std::to_string((int)direction);
    }
}


/*static*/ std::string EventWitness::toString(Continuity continuity) {
    switch(continuity) {
    case Continuous:        return "Continuous";
    case Discontinuous:     return "Discontinuous";
    default: return "BAD WITNESS Continuity " + std::to_string((int)continuity);
    }
}

/*static*/ std::string EventWitness::toString(TransitionMask m) {
    if (m==NoTransition) return "NoEventTrigger";

    // Unmask one at a time.
    const TransitionMask allTransitions[] =
     { NegativeToZero,NegativeToPositive,ZeroToPositive,
       PositiveToNegative,PositiveToZero,ZeroToNegative, 
       NoTransition };
    const char* transitionNames[] =
     { "NegativeToZero","NegativeToPositive","ZeroToPositive",
       "PositiveToNegative","PositiveToZero","ZeroToNegative" };

    String s;
    for (int i=0; allTransitions[i] != NoTransition; ++i)
        if (m & allTransitions[i]) {
            if (s.size()) s += "|";
            s += transitionNames[i];
            m = TransitionMask((unsigned)m & ~((unsigned)allTransitions[i])); 
        }

    // should have accounted for everything by now
    if (m != NoTransition) {
        if (s.size()) s += " + ";
        s += "UNRECOGNIZED TRANSITION MASK GARBAGE ";
        s += String((unsigned)m, "0x%x");
    }
    return s;
}


