/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
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

#include "SimTKcommon.h"

using namespace SimTK;

/* This function converts the provided instance of `T` to Xml, deserializes
it from Xml, ensures the deserialized instance is equal to the original
instance, reserializes the deserialized instance back to Xml and ensures
the original deserialization and the reserialization are equal.

One might set `checkInMemoryEquality` to false if `value` contains NaNs. */
template <typename T>
void testXmlRoundtrip(const T& value, bool checkInMemoryEquality = true,
                      std::string name = "") {
    Xml::Element xml = toXmlElementHelper(value, name, true);
    T deserialized;
    fromXmlElementHelper(deserialized, xml, name, true);
    if (checkInMemoryEquality) {
        try {
            SimTK_TEST(value == deserialized);
        }
        catch (Exception::Assert& e) {
            std::cout << "value: " << value << std::endl;
            std::cout << "serialized:\n" << xml << std::endl;
            std::cout << "deserialized: " << deserialized << std::endl;
            throw e;
        }
    }
    
    Xml::Element xmlReserialized = toXmlElementHelper(deserialized, name, true);
    String xmlStr, xmlReserializedStr;
    xml.writeToString(xmlStr);
    xmlReserialized.writeToString(xmlReserializedStr);
    try {
        SimTK_TEST(xmlStr == xmlReserializedStr);
    } catch (Exception::Assert& e) {
        std::cout << "value: " << value << std::endl;
        std::cout << "serialized:\n" << xml << std::endl;
        std::cout << "deserialized:" << deserialized << std::endl;
        std::cout << "reserialized:\n" << xmlReserialized << std::endl;
        throw e;
    }
        
}

void testStageXml() {
    auto testStageLevel = [](Stage::Level level) {
        Stage stage0(level);
        Xml::Element stageXml = toXmlElementHelper(stage0, "", true);
        Stage stage1;
        fromXmlElementHelper(stage1, stageXml, "", true);
        SimTK_TEST(stage1 == level);
    };
    
    testStageLevel(Stage::Empty);
    testStageLevel(Stage::Topology);
    testStageLevel(Stage::Model);
    testStageLevel(Stage::Instance);
    testStageLevel(Stage::Time);
    testStageLevel(Stage::Position);
    testStageLevel(Stage::Velocity);
    testStageLevel(Stage::Dynamics);
    testStageLevel(Stage::Acceleration);
    testStageLevel(Stage::Report);
    testStageLevel(Stage::Infinity);
    
    // Test invalid stage name.
    {
        std::stringstream ss;
        ss << "NotAStageLevelName";
        Stage stage;
        SimTK_TEST_MUST_THROW_EXC(ss >> stage, Exception::ErrorCheck);
    }
}

void testRotationXml() {
    testXmlRoundtrip(Rotation());
    testXmlRoundtrip(Rotation(BodyRotationSequence,
                              0.1, XAxis, 0.3, YAxis, 0.4, ZAxis));
    testXmlRoundtrip(Rotation(NaN, XAxis), false);
}

int main() {
    SimTK_START_TEST("TestStateXml");
        SimTK_SUBTEST(testStageXml);
        SimTK_SUBTEST(testRotationXml);
    SimTK_END_TEST();
}
