/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

#include "simbody/internal/common.h"
#include "simbody/internal/ImpulseSolver.h"

namespace SimTK {

// These static methods assume the "NA" value is -1 and the others count up
// contiguously from there.
const char* ImpulseSolver::getContactTypeName(ContactType ct) {
    static const char* nm[]={"TypeNA", "Observing", "Known", "Participating"};
    return TypeNA<=ct&&ct<=Participating ? nm[ct+1] : "UNKNOWNContactType";
}
const char* ImpulseSolver::getUniCondName(UniCond uc) {
    static const char* nm[]={"UniNA", "UniOff", "UniActive", "UniKnown"};
    return UniNA<=uc&&uc<=UniKnown ? nm[uc+1] : "UNKNOWNUniCond";
}
const char* ImpulseSolver::getFricCondName(FricCond fc) {
    static const char* nm[]={"FricNA", "FricOff", "Sliding", 
                            "Impending", "Rolling"};
    return FricNA<=fc&&fc<=Rolling ? nm[fc+1] : "UNKNOWNFricCond";
}
const char* ImpulseSolver::getBndCondName(BndCond bc) {
    static const char* nm[]={"BndNA", "SlipLow", "ImpendLow", "Engaged", 
                            "ImpendHigh", "SlipHigh"};
    return BndNA<=bc&&bc<=SlipHigh ? nm[bc+1] : "UNKNOWNBndCond";
}

void ImpulseSolver::
dumpUniContacts(const String& msg,
                const Array_<UniContactRT>& uniContacts) 
{
    printf("\n----------%s----------\n", msg.c_str());
    printf("Unilateral contact runtimes (n=%d):\n", (int)uniContacts.size());
    for (unsigned i=0; i < uniContacts.size(); ++i) {
        const UniContactRT& rt = uniContacts[i];
        printf("%d: UnilateralContactIndex %d, normal MultiplierIndex %d, "
               "sign=%g, hasFriction: %d\n", i, (int)rt.m_ucx, (int)rt.m_Nk,
               rt.m_sign, rt.hasFriction());
        if (rt.hasFriction()) {
            std::cout << "  friction MultIndices: " << rt.m_Fk   
                      << " effMu=" << rt.m_effMu << std::endl;
        }
        printf("  inputs: ContactType=%s, effCOR=%g\n",
               getContactTypeName(rt.m_type), rt.m_effCOR);
        printf("  outputs: normal cond=%s", 
               getUniCondName(rt.m_contactCond));
        if (rt.hasFriction()) {
            printf(" friction cond=%s slipV=%g %g, |slipV|=%g",
                   getFricCondName(rt.m_frictionCond),
                   rt.m_slipVel[0], rt.m_slipVel[1], rt.m_slipMag);
        }
        printf("\n");
    }
    printf("------------------------------\n\n");
}

} // namespace SimTK