/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2016 Stanford University and the Authors.           *
 * Authors: Chris Dembia, Michael Sherman                                     *
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

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/Value.h"

#include <string>
#include <unordered_map>
#include <cstdio>

using namespace SimTK;

// This is filled in automatically upon first construction of an object of
// type Value<T> for a new T.
std::unordered_map<std::string,AbstractValue::ValueFactory> 
AbstractValue::m_factory;

/*static*/ std::unique_ptr<AbstractValue> 
AbstractValue::createFromXmlElement(Xml::Element&      elt,
                                    const std::string& expectedName)
{
    const String& key = elt.getRequiredAttributeValue("type");
    auto entry = m_factory.find(key);
    SimTK_ERRCHK1_ALWAYS(entry != m_factory.end(), 
        "AbstractValue::createFromXmlElement()",
        "No ValueFactory with key='%s' was registered.", key.c_str());

    return entry->second(elt, expectedName);
}


/*static*/ auto AbstractValue::
registerValueFactory(const std::string& key, ValueFactory factory) 
    -> const ValueFactory& {
    auto result = m_factory.emplace(key,factory);
    SimTK_ERRCHK1_ALWAYS(result.second, 
        "AbstractValue::registerValueFactory()",
        "A ValueFactory with key='%s' has already been registered.",
        key.c_str());
    return result.first->second;
}
