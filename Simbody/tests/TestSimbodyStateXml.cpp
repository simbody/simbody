/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKsimbody                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-16 Stanford University and the Authors.        *
 * Authors: Chris Dembia                                                      *
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

#include "SimTKsimbody.h"

using namespace SimTK;

void testMotionLevelXml() {
    auto testMotionLevel = [](Motion::Level level) {
        Motion::Level ml0(level);
        Xml::Element mlXml = toXmlElementHelper(ml0, "", true);
        Motion::Level ml1;
        fromXmlElementHelper(ml1, mlXml, "", true);
        SimTK_TEST(ml1 == level);
    };
    
    testMotionLevel(Motion::Level::NoLevel);
    testMotionLevel(Motion::Level::Acceleration);
    testMotionLevel(Motion::Level::Velocity);
    testMotionLevel(Motion::Level::Position);
    

    // Test invalid level.
    {
        std::stringstream ss;
        ss.exceptions(std::stringstream::failbit);
        ss << "3";
        Motion::Level level;
        SimTK_TEST_MUST_THROW_EXC(ss >> level, std::stringstream::failure);
    }
    {
        Motion::Level ml0 = Motion::Level(3);
        Xml::Element mlXml = toXmlElementHelper(ml0, "", true);
        Motion::Level ml1;
        SimTK_TEST_MUST_THROW_EXC(fromXmlElementHelper(ml1, mlXml, "", true),
                Exception::ErrorCheck);
    }
}

int main() {
    SimTK_START_TEST("TestSimbodyStateXml");
        SimTK_SUBTEST(testMotionLevelXml);
    SimTK_END_TEST();
}
