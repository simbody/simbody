#ifndef SimTK_SIMBODY_EXAMPLE_ATLAS_URDF_READER_H_
#define SimTK_SIMBODY_EXAMPLE_ATLAS_URDF_READER_H_

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

/* This is a stripped-down .urdf robot file reader, with just enough 
functionality for this example. URDF is an XML format for ROS from the Open 
Source Robotics Foundation (http://osrfoundation.org). See
http://wiki.ros.org/urdf/XML for a description of the file format. Don't 
confuse this with Gazebo's .sdf format.

URDF file format notes
----------------------
A urdf file describes a robot in XML format using the following objects:
robot   Root element that contains links and joints.
link    Named mass- and geometry-carrying object (a body). There is a link 
        frame, and inertial, visual, and collision objects each with their own 
        frame given relative to the link frame. An initial pose for the link 
        frame is given, relative to the model frame.
joint   Connection between two links or between a link and the world.
        There is a parent and child link, given by name. There is a joint type 
        and then appropriate parameters depending on the type. There is a joint
        frame whose pose is given w.r.t. the parent frame. Unfortunately this
        must also be the child body frame.

The world frame is assumed to have x forward (roll axis), y left (pitch axis),
and z=x X y up (yaw axis).

An <origin> element (with an 'xyz' and 'rpy' attribute) is equivalent to a 
Simbody Transform, interpreted as x,y,z translation vector followed by a
roll-pitch-yaw (x-y-z) body-fixed Euler angle sequence given in radians.

An <inertial> element has <mass>, <origin>, <inertia> subelements. The origin
location is the link's COM; the inertia is given about the COM and oriented
as specified by the <origin>'s rpy orientation attribute.

A <visual> or <collision> element has <origin> and <geometry> subelements.
*/
#include "Simbody.h"

#include <utility>
#include <string>
#include <vector>
#include <map>
#include <iostream>

//==============================================================================
//                                URDF
//==============================================================================
// Some utility methods for URDF reading.
class URDF {
public:
    static SimTK::Transform 
        origin2Transform(const SimTK::Vec3& xyz, const SimTK::Vec3& rpy);
    static SimTK::Vec6 
        transform2Pose(const SimTK::Transform& X_AB);
    static SimTK::Transform 
        getOrigin(SimTK::Xml::Element element);
    static SimTK::MassProperties 
        getMassProperties(SimTK::Xml::Element link);
};

//==============================================================================
//                             URDF LINK INFO
//==============================================================================
// This is one link's information read from the .urdf input file and translated
// into Simbody's terminology and conventions. The mapping to Simbody 
// MobilizedBody is written here after we build the Simbody System.
// (Ignore master/slave stuff for URDF which can only do trees.)
class URDFLinkInfo {
public:
    explicit URDFLinkInfo(const std::string& name)
    :   name(name), mustBeBaseLink(false), selfCollide(false) {}

    // When a link is broken into several fragments (master and slaves), they
    // share the mass equally. Given the number of fragments, this returns the
    // appropriate mass properties to use for each fragment. Per Simbody's
    // convention, COM is measured from, and inertia taken about, the link 
    // origin and both are expressed in the link frame.
    SimTK::MassProperties getEffectiveMassProps(int numFragments) const {
        assert(numFragments > 0); // must be at least 1 for the master
        return SimTK::MassProperties(massProps.getMass()/numFragments,
                                     massProps.getMassCenter(),
                                     massProps.getUnitInertia());
    }

    SimTK::Xml::Element element;        // if this was in the Xml file
    std::string         name;
    bool                mustBeBaseLink;

    bool                selfCollide; // if true can collide with other links in
                                     // the same URDF model.

    // Mass properties converted to Simbody terms: com & inertia expressed in
    // body frame; inertia taken about body origin *not* body COM. This is the
    // full mass; if there are slaves the master and slaves will split the mass
    // properties evenly.
    SimTK::MassProperties                  massProps;

    // Which MobilizedBody corresponds to the master instance of this link.
    SimTK::MobilizedBody                   masterMobod;

    // If this link got split into a master and slaves, these are the 
    // MobilizedBodies used to mobilize the slaves.
    std::vector<SimTK::MobilizedBody>      slaveMobods;

    // And these are the Weld constraints used to attach slaves to master.
    std::vector<SimTK::Constraint::Weld>   slaveWelds;
};


//==============================================================================
//                            URDF JOINT INFO
//==============================================================================
// This is one joint's information read from the URDF input file and 
// translated into Simbody's terminology and conventions. The joint will
// typically have been modeled as a Mobilizer, but may have been modeled as
// a Constraint; in either case the corresponding Simbody element is written
// here after we build the Simbody System.
// (Ignore constraint for URDF.)
class URDFJointInfo {
public:
    URDFJointInfo(const std::string& name, const std::string& type) 
    :   name(name), type(type), mustBreakLoopHere(false),
        effort(SimTK::Infinity), lower(-SimTK::Infinity), 
        upper(SimTK::Infinity), velocity(SimTK::Infinity),
        isReversed(false)
    {}

    // These are set when we process the input.

    SimTK::Xml::Element     element; // if this was in the Xml file
    std::string             name, type, parent, child;
    bool                    mustBreakLoopHere;

    // Normally A=F, B=M. But if reversed, then B=F, A=M.
    SimTK::Transform        X_PA; // parent body frame to mobilizer frame
    SimTK::Transform        X_CB; // child body frame to mobilizer frame
    SimTK::Transform     defX_AB; // default mobilizer pose

    // Limits on force, position, velocity
    SimTK::Real            effort;         // |f| <= effort, effort >= 0
    SimTK::Real            lower, upper;   // lower <= upper
    SimTK::Real            velocity;       // |u| <= velocity, velocity >= 0

    // Members below here are set when we build the Simbody model.

    // How this joint was modeled in the Simbody System. We used either a
    // mobilizer or a constraint, but not both. The type of either one is the
    // same as the joint type above.
    SimTK::MobilizedBody    mobod;      // isValid() if we used a mobilizer
    bool                    isReversed; // if mobilizer, reverse parent&child?

    SimTK::Constraint       constraint; // isValid() if we used a constraint
};


//==============================================================================
//                               URDF LINKS
//==============================================================================
// This collects all the links for a particular model, and maintains a 
// name->URDFLinkInfo mapping. Be sure to add a world link first, with the
// appropriate pose measured from the model frame.
class URDFLinks {
public:
    // Add a new link to the collection and index it by name.
    int addLink(const URDFLinkInfo& info);
    // Return the number of links so far.
    int size() const {return (int)linkByIndex.size();}

    bool hasLink(const std::string& name) const 
    {   return name2index.find(name) != name2index.end(); }
    bool hasLink(int index) const
    {   assert(index>=0); return index < size(); }

    // Get link by name (const or writable).
    const URDFLinkInfo& getLink(const std::string& name) const 
    {   return getLink(getLinkIndex(name)); }
    URDFLinkInfo& updLink(const std::string& name) 
    {   return updLink(getLinkIndex(name)); }

    // Get link fast by index (const or writable).
    URDFLinkInfo& updLink(int index) 
    {   return linkByIndex[index]; }
    const URDFLinkInfo& getLink(int index) const
    {   return linkByIndex[index];}

private:
    int getLinkIndex(const std::string& name) const
    {   assert(hasLink(name)); return name2index.find(name)->second; }
    std::vector<URDFLinkInfo>   linkByIndex;
    std::map<std::string,int>   name2index;
};


//==============================================================================
//                             URDF JOINTS
//==============================================================================
// This collects all the joints for a particular model, and maintains a
// name->URDFJointInfo mapping.
class URDFJoints {
public:
    // Add a new joint to the collection, and return its assigned joint index,
    // starting from zero for the first joint and counting up from there.
    int addJoint(const URDFJointInfo& info);
    // Return the number of joints added so far.
    int size() const {return (int)jointByIndex.size();}

    bool hasJoint(const std::string& name) const 
    {   return name2index.find(name) != name2index.end(); }
    bool hasJoint(int index) const
    {   assert(index>=0); return index < size(); }

    // Get joint by name (const or writable).
    const URDFJointInfo& getJoint(const std::string& name) const 
    {   return getJoint(getJointIndex(name)); }
    URDFJointInfo& updJoint(const std::string& name) 
    {   return updJoint(getJointIndex(name)); }

    // Get joint fast by index (const or writable).
    const URDFJointInfo& getJoint(int index) const 
    {   return jointByIndex[index]; }
    URDFJointInfo& updJoint(int index) 
    {   return jointByIndex[index]; }

private:
    int getJointIndex(const std::string& name) const
    {   assert(hasJoint(name)); return name2index.find(name)->second; }

    std::vector<URDFJointInfo>  jointByIndex;
    std::map<std::string,int>   name2index;
};


//==============================================================================
//                               URDF ROBOT
//==============================================================================
// Contains robot-level information and the set of links and joints that
// were found in the input file for this robot.
class URDFRobot {
public:
    URDFRobot() : isStatic(false) {}
    void readRobot(SimTK::Xml::Element robotElt);

    std::string         name;
    SimTK::Transform    X_WM;     // model frame in the World frame
    bool                isStatic; // means all bodies are attached to Ground
    URDFLinks           links;
    URDFJoints          joints;

    // This is a grouping of contact surfaces that are invisible to one
    // another.
    SimTK::ContactCliqueId robotClique;
};

#endif // SimTK_SIMBODY_EXAMPLE_ATLAS_URDF_READER_H_