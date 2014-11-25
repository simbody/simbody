/* -------------------------------------------------------------------------- *
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

/* Implementation of URDFReader class. */

#include "Simbody.h"
#include "URDFReader.h"

#include <utility>
#include <string>
#include <vector>
#include <map>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;



//==============================================================================
//                    AUXILIARY FUNCTIONS FOR XML READING
//==============================================================================

// Convert the given pose in x,y,z,thetax,thetay,thetaz format to a Simbody
// Transform. The rotation angles are interpreted as a body-fixed sequence,
// meaning we rotation about x, then about the new y, then about the now twice-
// rotated z.
/*static*/ Transform URDF::origin2Transform(const Vec3& xyz, const Vec3& rpy) {
    Transform frame(Rotation(SimTK::BodyRotationSequence,
                             rpy[0], XAxis, rpy[1], YAxis, rpy[2], ZAxis),
                    xyz); 
    return frame;
}

// Convert a Simbody transform to a pose in x,y,z,thetax,thetay,thetaz format.
/*static*/ Vec6 URDF::transform2Pose(const Transform& X_AB) {
    Vec6 pose;
    pose.updSubVec<3>(0) = X_AB.p(); // position vector
    pose.updSubVec<3>(3) = X_AB.R().convertRotationToBodyFixedXYZ();
    return pose;
}

// If the given element contains an <origin> element, return it as a Transform.
// Otherwise return the identity Transform. If there is more than one <origin>
// element, only the first one is processed.
/*static*/ Transform URDF::getOrigin(Xml::Element element) {
    Xml::Element org = element.getOptionalElement("origin");
    if (!org.isValid()) return Transform();
    const Vec3 xyz = org.getOptionalAttributeValueAs<Vec3>("xyz",Vec3(0));
    const Vec3 rpy = org.getOptionalAttributeValueAs<Vec3>("rpy",Vec3(0));
    return origin2Transform(xyz,rpy);
}

// Look for an <inertial> element within the given element (which is probably
// a link but we don't care). If found, parse and transform to link (body)
// frame. Otherwise return unit mass and inertia. Result is returned as a
// Simbody MassProperties element containing mass, center of mass location,
// and unit inertia about body origin.
/*static*/ MassProperties URDF::getMassProperties(Xml::Element link) {
    Xml::Element inertial = link.getOptionalElement("inertial");
    if (!inertial.isValid())
        return MassProperties(1,Vec3(0),UnitInertia(1,1,1));
    Xml::Element massElt = inertial.getOptionalElement("mass");
    const Real mass = massElt.isValid() 
        ? massElt.getRequiredAttributeValueAs<Real>("value")
        : Real(1);
    Transform X_LI = getOrigin(inertial); // identity if not provided
    const Vec3 com_L = X_LI.p(); // vector from Lo to com, exp. in L
    Xml::Element inertia = inertial.getOptionalElement("inertia");
    if (mass==0 || !inertia.isValid())
        return MassProperties(mass,com_L,UnitInertia(1,1,1));
    // Get mass-weighted central inertia, expressed in I frame.
    Inertia Ic_I(inertia.getOptionalAttributeValueAs<Real>("ixx", 1.),
                 inertia.getOptionalAttributeValueAs<Real>("iyy", 1.),
                 inertia.getOptionalAttributeValueAs<Real>("izz", 1.),
                 inertia.getOptionalAttributeValueAs<Real>("ixy", 0.),
                 inertia.getOptionalAttributeValueAs<Real>("ixz", 0.),
                 inertia.getOptionalAttributeValueAs<Real>("iyz", 0.));
    // Re-express the central inertia from the I frame to the L frame.
    Inertia Ic_L = Ic_I.reexpress(~X_LI.R()); // Ic_L=R_LI*Ic_I*R_IL
    // Shift to L frame origin.
    Inertia Io_L = Ic_L.shiftFromMassCenter(-com_L, mass);
    return MassProperties(mass, com_L, Io_L); // converts to unit inertia
}

//==============================================================================
//                    IMPLEMENTATIONS OF URDF CLASSES
//==============================================================================

int URDFLinks::addLink(const URDFLinkInfo& info) {
    if (name2index.find(info.name) != name2index.end())
        throw std::runtime_error("URDFLinks::addLink(): Link name '" 
            + info.name + " was used more than once.");

    const int lx = (int)linkByIndex.size();
    linkByIndex.push_back(info);
    name2index[info.name] = lx;
    return lx;
}


int URDFJoints::addJoint(const URDFJointInfo& info) {
    if(name2index.find(info.name) != name2index.end())
        throw std::runtime_error(
            "URDFJoints::addJoint(): Joint name '" + info.name 
            + "' was used more than once.");

    const int jx = (int)jointByIndex.size();
    jointByIndex.push_back(info);
    name2index[info.name] = jx;
    return jx;
}



void URDFRobot::readRobot(Xml::Element robotElt) {
    name = robotElt.getRequiredAttributeValue("name");
    X_WM = URDF::getOrigin(robotElt);

    Array_<Xml::Element> linkElts = robotElt.getAllElements("link");
    Array_<Xml::Element> jointElts = robotElt.getAllElements("joint");
    printf("Reading robot '%s' with ground + %d links, %d joints.\n",
        name.c_str(), linkElts.size(), jointElts.size());

    // Create a World link, then read in the real links.
    URDFLinkInfo worldLink("world");
    worldLink.massProps = MassProperties(Infinity,Vec3(0),UnitInertia(1));
    links.addLink(worldLink);
    for (unsigned i=0; i < linkElts.size(); ++i) {
        Xml::Element elt = linkElts[i];
        URDFLinkInfo linfo(elt.getRequiredAttributeValue("name"));
        linfo.mustBeBaseLink = elt.getOptionalElementValueAs<bool>
                                                ("must_be_base_link", false);
        linfo.selfCollide = elt.getOptionalElementValueAs<bool>
                                                ("self_collide", false);
        linfo.element = elt;
        linfo.massProps = URDF::getMassProperties(elt);
        links.addLink(linfo);
    }

    // Read in the joints.
    for (unsigned i=0; i < jointElts.size(); ++i) {
        Xml::Element elt = jointElts[i];
        const std::string name = elt.getRequiredAttributeValue("name");
        const std::string type = elt.getRequiredAttributeValue("type");
        URDFJointInfo jinfo(name,type);
        jinfo.mustBreakLoopHere = elt.getOptionalElementValueAs<bool>
                                                ("must_be_loop_joint", false);
        jinfo.element = elt;

        std::string child, parent;
        Xml::Element cElt = elt.getRequiredElement("child");
        child = cElt.getRequiredAttributeValue("link");

        Xml::Element pElt = elt.getRequiredElement("parent");
        parent = pElt.getRequiredAttributeValue("link");

        assert(links.hasLink(parent) && links.hasLink(child));
        assert(child != parent);

        jinfo.parent = parent;
        jinfo.child  = child;

        jinfo.X_CB = Transform();          // joint frame on child body
        jinfo.X_PA = URDF::getOrigin(elt); // joint frame on parent body

        Xml::Element limit = elt.getOptionalElement("limit");
        if (limit.isValid()) {
            jinfo.effort = 
                limit.getOptionalAttributeValueAs<Real>("effort", Infinity);
            jinfo.lower = 
                limit.getOptionalAttributeValueAs<Real>("lower", -Infinity);
            jinfo.upper = 
                limit.getOptionalAttributeValueAs<Real>("upper", Infinity);
            jinfo.velocity = 
                limit.getOptionalAttributeValueAs<Real>("velocity", Infinity);

            if (jinfo.lower > jinfo.upper)
                throw std::runtime_error(
                    "URDFRobot::readRobot(): Joint name '" + name 
                    + "' had lower limit " + String(jinfo.lower)
                    + " > upper limit " + String(jinfo.upper) + ".");
        }

        joints.addJoint(jinfo);
    }
}

