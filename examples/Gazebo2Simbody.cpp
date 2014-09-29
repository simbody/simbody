/* -------------------------------------------------------------------------- *
 *                   Simbody(tm) Example: Gazebo 2 Simbody                    *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

/* This example reads input in the form of a Gazebo ".sdf" file and generates
a roughly-equivalent Simbody model. Gazebo is a robot simulator from the Open 
Source Robotics Foundation (http://osrfoundation.org). This is not intended to
be comprehensive -- it is just a proof of concept.

We will construct the Simbody model here and maintain a mapping between the
Gazebo objects and the Simbody implementation of them to demonstrate how this
might be done in Gazebo.

Gazebo file format notes
------------------------
An sdf file describes a model in XML format using the following objects:
world   The Ground frame. May have associated geometry and lights. Contains one
        or more "model" objects. Normally the ground is the x-y plane, with
        +z the "up" direction and a ground plane at z=0. The world also acts
        as a link whose name is "world".
model   Named grouping of physical objects. Provides a model frame given
        relative to the world frame, however this frame is still fixed in the
        world. A model contains links and joints.
link    Named mass- and geometry-carrying object (a body). There is a link 
        frame, and inertial, visual, and collision objects each with their own 
        frame given relative to the link frame. An initial pose for the link 
        frame is given, relative to the model frame.
joint   Connection between two links or between a link and the world.
        There is a parent and child link, given by name, with "world"
        the name for the ground link. There is a joint type and then appropriate
        parameters depending on the type. There is a child frame and a
        parent frame. However, if only a single frame is given for the joint
        it is interpreted as the child frame and then the default link frame
        poses are used to find the coincident frame on the parent.

NOTE: for purposes of this example we have implemented some extensions:
- You can specify that a link must be a base link, i.e. must be connected
  directly to World. Use 
      <must_be_base_link> 1 </must_be_base_link> 
  in the <link> element.
- You can force a loop to be cut at a particular joint. Use
      <must_be_loop_joint> 1 </must_be_loop_joint> 
  in the <joint> element.
- We accept John Hsu's proposed modification allowing a joint to specify
  separate child and parent frames. Use 
      <child  link="childName" ><pose>...</pose></child> 
      <parent link="parentName"><pose>...</pose></parent> 
  to specify separate frame. (We also accept the old style.)

A <pose> element (containing six scalars) is equivalent to a Simbody Transform,
interpreted as x,y,z translation vector followed by an x-y-z body-fixed Euler 
angle sequence given in radians. Gazebo refers to this as pitch-roll-yaw.

A <visual> element has <geometry> and <material> subelements; here we'll use the
material only to pick a color.

A <collision> element has <geometry> and <surface> subelements, with the
latter giving the physical properties used in contact. We'll generate a gray
visual for these.

Generalized coordinate representation
-------------------------------------
We use Simbody's MultibodyGraphMaker class to take the Gazebo definition of
links and joints and construct a spanning tree whose root is
the World, with every link appearing exactly once. Any link that does not appear
in any joint will be given six degrees of freedom relative to the World via a
Simbody "Free" mobilizer. Note: If a disconnected link is a point mass (that is,
inertialess) it should be mobilized with a 3 dof "Translation" mobilizer but 
we're not doing that in this example. Also if there are disconnected branches we 
will pick one of the bodies as a "base" body and connect it to World by a Free 
(6dof) mobilizer.

MultibodyGraphMaker won't take the parent-child specification too seriously when 
building the tree since these can be arbitrary in the Gazebo description. To 
avoid confusion, we speak of the "inboard" and "outboard" bodies of a mobilizer;
usually inboard==parent and outboard==child but the relationship is opposite for 
a reverse mobilizer. For example, a three-body chain might be defined
     link1--j1--link2--j2--link3
     joint1 parent=link1 child=link2
     joint2 parent=link3 child=link2
There is no way to build a tree where inboard/outboard matches parent/child
for both joints in this example. MultibodyGraphMaker will tell us if we have
to reverse, which we'll do using Simbody's reverse mobilizer capability so
that the meaning of the joint is unchanged. 

After the tree is built, all the links will have been used but there may still
be some unused joints. Those joints will involve parent and child links both
of which are already in the tree. That means these joints form topological
loops in the multibody graph. MultibodyGraphMaker will normally handle that by 
splitting the child link to create a new "slave" body, mobilizing the slave
with the given joint, and then introducing a Weld constraint (-6 dofs) to hold 
the master link and its slaves together. Note that in this example the visual 
and contact geometry stays with the master; the mass and inertia are divided 
equally among the master and slaves. We'll also generate some "shadow" 
visualization for the slaves since they are interesting.

An alternative is used if you tell MultibodyGraphMaker that there is a good
constraint that can be used directly to replace a loop joint. Since Simbody's
Ball constraint (-3 dofs) is mostly indistinguishable from a Ball mobilizer
(+3 dofs), we'll substitute a Ball constraint when breaking a loop at a
ball joint, resulting in a smaller system overall. */
#include "Simbody.h"

#include <utility>
#include <string>
#include <vector>
#include <map>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

// Skip down to main() first -- these are just some helper classes to keep
// track of the Gazebo model and connect it to its Simbody implementation.

// Hacks for sherm's testing.
#define ANIMATE             // turn this off for timing tests
//#define ADD_JOINT_SPRINGS // makes all the revolute & prismatic joints stiff
//#define RAGDOLL2_ICS      // set initial conditions for ragdoll2 model
//#define USE_CONTACT_MESH  // use elastic foundation model for cylinder

//==============================================================================
//                            GAZEBO LINK INFO
//==============================================================================
// This is one link's information read from the Gazebo input file and translated
// into Simbody's terminology and conventions. The mapping to Simbody 
// MobilizedBody is written here after we build the Simbody System.
class GazeboLinkInfo {
public:
    explicit GazeboLinkInfo(const std::string& name)
    :   name(name), mustBeBaseLink(false), selfCollide(false) {}

    // When a link is broken into several fragments (master and slaves), they
    // share the mass equally. Given the number of fragments, this returns the
    // appropriate mass properties to use for each fragment. Per Simbody's
    // convention, COM is measured from, and inertia taken about, the link 
    // origin and both are expressed in the link frame.
    MassProperties getEffectiveMassProps(int numFragments) const {
        assert(numFragments > 0); // must be at least 1 for the master
        return MassProperties(massProps.getMass()/numFragments,
                              massProps.getMassCenter(),
                              massProps.getUnitInertia());
    }

    Xml::Element    element;        // if this was in the Xml file
    std::string     name;
    bool            mustBeBaseLink;

    bool            selfCollide; // if true can collide with other links in
                                 // the same Gazebo model.

    // Mass properties converted to Simbody terms: com & inertia expressed in
    // body frame; inertia taken about body origin *not* body COM. This is the
    // full mass; if there are slaves the master and slaves will split the mass
    // properties evenly.
    MassProperties                  massProps;
    Transform                       defX_ML;    // default pose in model frame

    // Which MobilizedBody corresponds to the master instance of this link.
    MobilizedBody                   masterMobod;

    // If this link got split into a master and slaves, these are the 
    // MobilizedBodies used to mobilize the slaves.
    std::vector<MobilizedBody>      slaveMobods;

    // And these are the Weld constraints used to attach slaves to master.
    std::vector<Constraint::Weld>   slaveWelds;
};


//==============================================================================
//                            GAZEBO JOINT INFO
//==============================================================================
// This is one joint's information read from the Gazebo input file and 
// translated into Simbody's terminology and conventions. The joint will
// typically have been modeled as a Mobilizer, but may have been modeled as
// a Constraint; in either case the corresponding Simbody element is written
// here after we build the Simbody System.
class GazeboJointInfo {
public:
    GazeboJointInfo(const std::string& name, const std::string& type) 
    :   name(name), type(type), mustBreakLoopHere(false), isReversed(false) {}

    // These are set when we process the input.

    Xml::Element    element; // if this was in the Xml file
    std::string     name, type, parent, child;
    bool            mustBreakLoopHere;

    // Normally A=F, B=M. But if reversed, then B=F, A=M.
    Transform       X_PA; // parent body frame to mobilizer frame
    Transform       X_CB; // child body frame to mobilizer frame
    Transform    defX_AB; // default mobilizer pose

    // Members below here are set when we build the Simbody model.

    // How this joint was modeled in the Simbody System. We used either a
    // mobilizer or a constraint, but not both. The type of either one is the
    // same as the joint type above.
    MobilizedBody   mobod;      // isValid() if we used a mobilizer
    bool            isReversed; // if mobilizer, did it reverse parent&child?

    Constraint      constraint; // isValid() if we used a constraint
};


//==============================================================================
//                               GAZEBO LINKS
//==============================================================================
// This collects all the links for a particular model, and maintains a 
// name->GazeboLinkInfo mapping. Be sure to add a world link first, with the
// appropriate pose measured from the model frame.
class GazeboLinks {
public:
    // Add a new link to the collection and index it by name.
    int addLink(const GazeboLinkInfo& info);
    // Return the number of links so far.
    int size() const {return (int)linkByIndex.size();}

    bool hasLink(const std::string& name) const 
    {   return name2index.find(name) != name2index.end(); }
    bool hasLink(int index) const
    {   assert(index>=0); return index < size(); }

    // Get link by name (const or writable).
    const GazeboLinkInfo& getLink(const std::string& name) const 
    {   return getLink(getLinkIndex(name)); }
    GazeboLinkInfo& updLink(const std::string& name) 
    {   return updLink(getLinkIndex(name)); }

    // Get link fast by index (const or writable).
    GazeboLinkInfo& updLink(int index) 
    {   return linkByIndex[index]; }
    const GazeboLinkInfo& getLink(int index) const
    {   return linkByIndex[index];}

private:
    int getLinkIndex(const std::string& name) const
    {   assert(hasLink(name)); return name2index.find(name)->second; }
    std::vector<GazeboLinkInfo>     linkByIndex;
    std::map<std::string,int>       name2index;
};


//==============================================================================
//                             GAZEBO JOINTS
//==============================================================================
// This collects all the joints for a particular model, and maintains a
// name->GazeboJointInfo mapping.
class GazeboJoints {
public:
    // Add a new joint to the collection, and return its assigned joint index,
    // starting from zero for the first joint and counting up from there.
    int addJoint(const GazeboJointInfo& info);
    // Return the number of joints added so far.
    int size() const {return (int)jointByIndex.size();}

    bool hasJoint(const std::string& name) const 
    {   return name2index.find(name) != name2index.end(); }
    bool hasJoint(int index) const
    {   assert(index>=0); return index < size(); }

    // Get joint by name (const or writable).
    const GazeboJointInfo& getJoint(const std::string& name) const 
    {   return getJoint(getJointIndex(name)); }
    GazeboJointInfo& updJoint(const std::string& name) 
    {   return updJoint(getJointIndex(name)); }

    // Get joint fast by index (const or writable).
    const GazeboJointInfo& getJoint(int index) const 
    {   return jointByIndex[index]; }
    GazeboJointInfo& updJoint(int index) 
    {   return jointByIndex[index]; }

private:
    int getJointIndex(const std::string& name) const
    {   assert(hasJoint(name)); return name2index.find(name)->second; }

    std::vector<GazeboJointInfo>    jointByIndex;
    std::map<std::string,int>       name2index;
};


//==============================================================================
//                              GAZEBO MODEL
//==============================================================================
// Contains model-level information and the set of links and joints that
// were found in the input file for this model.
class GazeboModel {
public:
    GazeboModel() : isStatic(false) {}
    void readModel(Xml::Element modelElt);

    std::string     name;
    Transform       X_WM;     // model frame in the World frame
    bool            isStatic; // means all bodies are attached to Ground
    GazeboLinks     links;
    GazeboJoints    joints;

    // This is a grouping of contact surfaces that are invisible to one
    // another. Gazebo's convention is that surfaces in the same model don't
    // collide unless they have been explicitly marked <self_collide>. The
    // Simbody equivalent is to put the unmarked ones into the same clique.
    ContactCliqueId modelClique;
};

//==============================================================================
//                                   MAIN
//==============================================================================



static void createMultibodyGraph(GazeboModel&           model,
                                 MultibodyGraphMaker&   mbgraph);

static void addModelToSimbodySystem(const MultibodyGraphMaker& mbgraph,
                                    GazeboModel&               model,
                                    MultibodySystem&           mbs,
                                    SimbodyMatterSubsystem&    matter,
                                    GeneralForceSubsystem&     forces,
                                    CompliantContactSubsystem& contact);

static void runSimulation(const MultibodySystem&          mbs,
                          const std::vector<GazeboModel>& gzModels);

int main(int argc, const char* argv[]) {
    std::string sdfFileName;
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " Gazebo_filename.sdf\n";
        std::cout << "Trying Gazebo_ragdoll.sdf by default.\n";
        sdfFileName = "Gazebo_ragdoll.sdf";
    } else
        sdfFileName = argv[1];

    try {
    //---------------------------- OPEN INPUT FILE -----------------------------
    // Attempt to identify this as a Gazebo input file, find the World
    // to get gravity and locate the Models.
    std::cout << "Working directory: " 
              << Pathname::getCurrentWorkingDirectory() << std::endl;
    std::cout << "Reading file: " << sdfFileName << std::endl;
    Xml::Document sdf(sdfFileName);

    if (sdf.getRootTag() != "sdf" && sdf.getRootTag() != "gazebo")
        throw std::runtime_error
           ("Expected to see document tag <sdf> or <gazebo> but saw <"
            + sdf.getRootTag() + "> instead.");

    // This is a Gazebo document.
    Xml::Element root = sdf.getRootElement();
    cout << "sdf version=" 
         << root.getOptionalAttributeValue("version", "unspecified") << endl;

    Xml::Element world = root.getRequiredElement("world");
    Xml::Element physics = world.getOptionalElement("physics");
    const Vec3 gravity = // default is std gravity at Earth's surface, m/s^2
        world.hasElement("gravity") 
        ? world.getRequiredElementValueAs<Vec3>("gravity")
        : (physics.hasElement("gravity") ? physics.getRequiredElementValueAs<Vec3>("gravity")
                                         : Vec3(0,0,-9.80665));

    Array_<Xml::Element> models = world.getAllElements("model");
    cout << "World '" << world.getOptionalAttributeValue("name","NONAME")
         << "' contains " << models.size() << " model(s).\n";
    cout << "Gravity=" << gravity << endl;

    if (models.empty()) {
        cout << "File contained no model -- nothing to do. Goodbye.\n";
        return 0;
    }
   
    //------------------------ CREATE SIMBODY SYSTEM ---------------------------
    // Create a Simbody System and populate it with Subsystems we'll need.
    MultibodySystem mbs; 
    SimbodyMatterSubsystem matter(mbs);
    GeneralForceSubsystem forces(mbs);
    ContactTrackerSubsystem tracker(mbs);
    CompliantContactSubsystem contact(mbs, tracker);
    // Tweak visualization to get Z up and to prevent Simbody from generating
    // its own sketchy visuals.
    mbs.setUpDirection(ZAxis);
    matter.setShowDefaultGeometry(false);
    // Set stiction max slip velocity to make it less stiff.
    contact.setTransitionVelocity(0.05);
    // Specify gravity (read in above from world).
    Force::UniformGravity(forces, matter, gravity);
    // Define a material to use for ground contact. This is not very stiff.
    ContactMaterial groundMaterial(1e6,   // stiffness
                             0.1,  // dissipation
                             0.7,   // mu_static
                             0.5,   // mu_dynamic
                             0.5);  // mu_viscous
    // Add a contact surface to represent the ground.
    // Half space normal is -x; must rotate about y to make it +z.
    matter.Ground().updBody().addContactSurface(Rotation(Pi/2,YAxis),
       ContactSurface(ContactGeometry::HalfSpace(), groundMaterial));
    // Draw world frame and model frame.
    matter.Ground().addBodyDecoration(Vec3(0), 
        DecorativeFrame(1).setColor(Green).setLineThickness(3)); // World

    std::vector<GazeboModel> gzModels(models.size());
    //------------------------ READ IN GAZEBO MODELS ---------------------------
    for (unsigned m=0; m < models.size(); ++m) {
        gzModels[m].readModel(models[m]); // TODO: ignore other models for now

        //----------------------- GENERATE MULTIBODY GRAPH -------------------------
        MultibodyGraphMaker mbgraph;
        createMultibodyGraph(gzModels[m], mbgraph);
        // Optional: dump the graph to stdout for debugging or curiosity.
        mbgraph.dumpGraph(std::cout);

        addModelToSimbodySystem(mbgraph, gzModels[m], 
                                mbs, matter, forces, contact);
    }

    //--------------------------- RUN A SIMULATION -----------------------------
    // The Simbody System has been built successfully. Now lets run it for
    // a while to see what it looks like.
    runSimulation(mbs, gzModels);

    } catch (const std::exception& e) {
        cout << "EXCEPTION: " << e.what() << "\n";
        return 1;
    }

    printf("DONE.\n");
    return 0;
}


//==============================================================================
//                           CREATE MULTIBODY GRAPH
//==============================================================================
// Define Gazebo joint types, then use links and joints in the given model
// to construct a reasonable spanning-tree-plus-constraints multibody graph
// to represent that model. An exception will be thrown if this fails.
// Note that this step is not Simbody dependent.
static void createMultibodyGraph(GazeboModel&           model,
                                 MultibodyGraphMaker&   mbgraph) {
    // Step 1: Tell MultibodyGraphMaker about joints it should know about.
    // Note: "weld" and "free" are always predefined at 0 and 6 dofs, resp.
    //                  Gazebo name  #dofs     Simbody equivalent
    mbgraph.addJointType("revolute",  1);   // Pin
    mbgraph.addJointType("revolute2", 2);   // ?
    mbgraph.addJointType("prismatic", 1);   // Slider
    mbgraph.addJointType("universal", 2);   // Universal
    mbgraph.addJointType("piston",    2);   // Cylinder
    // Simbody has a Ball constraint that is a good choice if you need to
    // break a loop at a ball joint.
    mbgraph.addJointType("ball", 3, true);  // Ball

    // Step 2: Tell it about all the links we read from the input file, 
    // starting with world, and provide a reference pointer.
    for (int lx=0; lx < model.links.size(); ++lx) {
        GazeboLinkInfo& link = model.links.updLink(lx);
        mbgraph.addBody(link.name, link.massProps.getMass(), 
                        link.mustBeBaseLink, &link);
    }

    // Step 3: Tell it about all the joints we read from the input file,
    // and provide a reference pointer.
    for (int jx=0; jx < model.joints.size(); ++jx) {
        GazeboJointInfo& joint = model.joints.updJoint(jx);
        mbgraph.addJoint(joint.name, joint.type, joint.parent, joint.child, 
                            joint.mustBreakLoopHere, &joint);
    }

    // Setp 4. Generate the multibody graph.
    mbgraph.generateGraph();
}



//==============================================================================
//                    AUXILIARY FUNCTIONS FOR XML READING
//==============================================================================

// Convert the given pose in x,y,z,thetax,thetay,thetaz format to a Simbody
// Transform. The rotation angles are interpreted as a body-fixed sequence,
// meaning we rotation about x, then about the new y, then about the now twice-
// rotated z.
static Transform pose2Transform(const Vec6& pose) {
    Transform frame(Rotation(SimTK::BodyRotationSequence,
                             pose[3], XAxis, pose[4], YAxis, pose[5], ZAxis),
                    pose.getSubVec<3>(0)); 
    return frame;
}

// Convert a Simbody transform to a pose in x,y,z,thetax,thetay,thetaz format.
static Vec6 transform2Pose(const Transform& X_AB) {
    Vec6 pose;
    pose.updSubVec<3>(0) = X_AB.p(); // position vector
    pose.updSubVec<3>(3) = X_AB.R().convertRotationToBodyFixedXYZ();
    return pose;
}

// If the given element contains a <pose> element, return it as a Transform.
// Otherwise return the identity Transform. If there is more than one <pose>
// element, only the first one is processed.
static Transform getPose(Xml::Element element) {
    const Vec6 pose = element.getOptionalElementValueAs<Vec6>("pose",Vec6(0));
    return pose2Transform(pose);
}

// Look for an <inertial> element within the given element (which is probably
// a link but we don't care). If found, parse and transform to link (body)
// frame. Otherwise return unit mass and inertia. Result is returned as a
// Simbody MassProperties element containing mass, center of mass location,
// and unit inertia about body origin.
static MassProperties getMassProperties(Xml::Element link) {
    Xml::Element inertial = link.getOptionalElement("inertial");
    if (!inertial.isValid())
        return MassProperties(1,Vec3(0),UnitInertia(1,1,1));
    const Real mass = inertial.getOptionalElementValueAs<Real>("mass", 1.);
    Transform X_LI = getPose(inertial); // identity if not provided
    const Vec3 com_L = X_LI.p(); // vector from Lo to com, exp. in L
    Xml::Element inertia = inertial.getOptionalElement("inertia");
    if (mass==0 || !inertia.isValid())
        return MassProperties(mass,com_L,UnitInertia(1,1,1));
    // Get mass-weighted central inertia, expressed in I frame.
    Inertia Ic_I(inertia.getOptionalElementValueAs<Real>("ixx", 1.),
                 inertia.getOptionalElementValueAs<Real>("iyy", 1.),
                 inertia.getOptionalElementValueAs<Real>("izz", 1.),
                 inertia.getOptionalElementValueAs<Real>("ixy", 0.),
                 inertia.getOptionalElementValueAs<Real>("ixz", 0.),
                 inertia.getOptionalElementValueAs<Real>("iyz", 0.));
    // Re-express the central inertia from the I frame to the L frame.
    Inertia Ic_L = Ic_I.reexpress(~X_LI.R()); // Ic_L=R_LI*Ic_I*R_IL
    // Shift to L frame origin.
    Inertia Io_L = Ic_L.shiftFromMassCenter(-com_L, mass);
    return MassProperties(mass, com_L, Io_L); // converts to unit inertia
}


//==============================================================================
//                       ADD MODEL TO SIMBODY SYSTEM
//==============================================================================
// Given a desired multibody graph, gravity, and the Gazebo model that was
// used to generate the graph, add elements to the Simbody System to represent
// it. There are many limitations here, especially in the handling of contact. 
// Any Gazebo features that we haven't modeled are just ignored.
// The GazeboModel is updated so that its links and joints have references to
// their corresponding Simbody elements.
// We set up some visualization here so we can see what's happening but this
// would not be needed in Gazebo since it does its own visualization.
static void addModelToSimbodySystem(const MultibodyGraphMaker& mbgraph,
                                    GazeboModel&               model,
                                    MultibodySystem&           mbs,
                                    SimbodyMatterSubsystem&    matter,
                                    GeneralForceSubsystem&     forces,
                                    CompliantContactSubsystem& contact) 
{
    // Define a material to use for contact. This is not very stiff.
    ContactMaterial material(1e6,   // stiffness
                             0.1,  // dissipation
                             0.7,   // mu_static
                             0.5,   // mu_dynamic
                             0.5);  // mu_viscous

    // Draw the model frame.
    matter.Ground().addBodyDecoration(model.X_WM, 
        DecorativeFrame(.75).setColor(Orange).setLineThickness(3));  // Model
    matter.Ground().addBodyDecoration(model.X_WM.p()+Vec3(0,0,.1),
        DecorativeText(model.name).setScale(.2)
        .setColor(model.isStatic?Green:Cyan));

    // Generate a contact clique we can put collision geometry in to prevent
    // self-collisions.
    model.modelClique = ContactSurface::createNewContactClique();

    // Record the MobilizedBody for the World link.
    model.links.updLink(0).masterMobod = matter.Ground();

    // Run through all the mobilizers in the multibody graph, adding a Simbody
    // MobilizedBody for each one. Also add visual and collision geometry to the
    // bodies when they are mobilized.
    for (int mobNum=0; mobNum < mbgraph.getNumMobilizers(); ++mobNum) {
        // Get a mobilizer from the graph, then extract its corresponding
        // joint and bodies. Note that these don't necessarily have equivalents
        // in the GazeboLink and GazeboJoint inputs.
        const MultibodyGraphMaker::Mobilizer& mob = mbgraph.getMobilizer(mobNum);
        const std::string& type = mob.getJointTypeName();

        // The inboard body always corresponds to one of the input links,
        // because a slave link is always the outboard body of a mobilizer.
        // The outboard body may be slave, but its master body is one of the
        // Gazebo input links.
        const bool isSlave = mob.isSlaveMobilizer();
        GazeboLinkInfo& gzInb  = *(GazeboLinkInfo*)mob.getInboardBodyRef();
        GazeboLinkInfo& gzOutb = *(GazeboLinkInfo*)mob.getOutboardMasterBodyRef();

        const MassProperties massProps = 
            gzOutb.getEffectiveMassProps(mob.getNumFragments());

        // This will reference the new mobilized body once we create it.
        MobilizedBody mobod; 

        if (model.isStatic) {
            mobod = matter.updGround();
        } else if (mob.isAddedBaseMobilizer()) {
            // There is no corresponding Gazebo joint for this mobilizer.
            // Create the joint and set its default position to be the default
            // pose of the base link relative to the Ground frame.
            assert(type=="free"); // May add more types later
            if (type == "free") {
                MobilizedBody::Free freeJoint(
                    gzInb.masterMobod,  Transform(),
                    massProps,    Transform());
                freeJoint.setDefaultTransform(~gzInb.defX_ML*gzOutb.defX_ML);
                mobod = freeJoint;
            }
        } else {
            // This mobilizer does correspond to one of the input joints.
            GazeboJointInfo& gzJoint = *(GazeboJointInfo*)mob.getJointRef();
            const bool isReversed = mob.isReversedFromJoint();

            // Find inboard and outboard frames for the mobilizer; these are
            // parent and child frames or the reverse.

            const Transform& X_IF0 = isReversed ? gzJoint.X_CB : gzJoint.X_PA;
            const Transform& X_OM0 = isReversed ? gzJoint.X_PA : gzJoint.X_CB;

            const MobilizedBody::Direction direction =
                isReversed ? MobilizedBody::Reverse : MobilizedBody::Forward;

            if (type == "free") {
                MobilizedBody::Free freeJoint(
                    gzInb.masterMobod,  X_IF0,
                    massProps,          X_OM0, 
                    direction);
                Transform defX_FM = isReversed ? Transform(~gzJoint.defX_AB)
                                               : gzJoint.defX_AB;
                freeJoint.setDefaultTransform(defX_FM);
                mobod = freeJoint;
            } else if (type == "revolute") {
                Xml::Element axisElt = gzJoint.element.getRequiredElement("axis");
                UnitVec3 axis = 
                    UnitVec3(axisElt.getRequiredElementValueAs<Vec3>("xyz")); 
                Rotation R_JZ(axis, ZAxis); // Simbody's pin is along Z
                Transform X_IF(X_IF0.R()*R_JZ, X_IF0.p());
                Transform X_OM(X_OM0.R()*R_JZ, X_OM0.p());
                MobilizedBody::Pin pinJoint(
                    gzInb.masterMobod,      X_IF,
                    massProps,              X_OM, 
                    direction);
                mobod = pinJoint;

                #ifdef ADD_JOINT_SPRINGS
                // KLUDGE add spring (stiffness proportional to mass)
                Force::MobilityLinearSpring(forces,mobod,0,
                                            30*massProps.getMass(),0);
                #endif
            } else if (type == "prismatic") {
                Xml::Element axisElt = gzJoint.element.getRequiredElement("axis");
                UnitVec3 axis = 
                    UnitVec3(axisElt.getRequiredElementValueAs<Vec3>("xyz")); 
                Rotation R_JX(axis, XAxis); // Simbody's slider is along X
                Transform X_IF(X_IF0.R()*R_JX, X_IF0.p());
                Transform X_OM(X_OM0.R()*R_JX, X_OM0.p());
                MobilizedBody::Slider sliderJoint(
                    gzInb.masterMobod,      X_IF,
                    massProps,              X_OM, 
                    direction);
                mobod = sliderJoint;

                #ifdef ADD_JOINT_SPRINGS
                // KLUDGE add spring (stiffness proportional to mass)
                Force::MobilityLinearSpring(forces,mobod,0,
                                            30*massProps.getMass(),0);
                #endif
            } else if (type == "ball") {
                MobilizedBody::Ball ballJoint(
                    gzInb.masterMobod,  X_IF0,
                    massProps,          X_OM0, 
                    direction);
                Rotation defR_FM = isReversed 
                    ? Rotation(~gzJoint.defX_AB.R())
                    : gzJoint.defX_AB.R();
                ballJoint.setDefaultRotation(defR_FM);
                mobod = ballJoint;
            } else if (type == "weld") {
                MobilizedBody::Weld weldJoint(
                    gzInb.masterMobod,  X_IF0,
                    massProps,          X_OM0);
                mobod = weldJoint;
            } 

            // Created a mobilizer that corresponds to gzJoint. Keep track.
            gzJoint.mobod = mobod;
            gzJoint.isReversed = isReversed;
        }

        // Link gzOutb has been mobilized; keep track for later.
        if (isSlave) gzOutb.slaveMobods.push_back(mobod);
        else gzOutb.masterMobod = mobod;

        // A mobilizer has been created; now add the visual and collision
        // geometry for the new mobilized body.

        Xml::Element master = gzOutb.element;
        Vec3 color = isSlave ? Red : Cyan;
        Real scale = isSlave ? 0.9 : 1.;
        if (master.isValid()) {
            // VISUAL
            Array_<Xml::Element> visuals = master.getAllElements("visual");
            for (unsigned i=0; i < visuals.size(); ++i)  {
                Transform X_LV = getPose(visuals[i]);
                Xml::Element geo = visuals[i].getRequiredElement("geometry");
                Xml::Element box = geo.getOptionalElement("box");
                if (box.isValid()) {
                    Vec3 sz = box.getRequiredElementValueAs<Vec3>("size");
                    mobod.addBodyDecoration(
                        X_LV, DecorativeBrick(sz/2).setOpacity(.5)
                                .setColor(color).setScale(scale));
                }
                Xml::Element sphere = geo.getOptionalElement("sphere");
                if (sphere.isValid()) {
                    Real r = sphere.getRequiredElementValueAs<Real>("radius");
                    mobod.addBodyDecoration(
                        X_LV, DecorativeSphere(r).setOpacity(.5)
                                .setColor(color).setScale(scale));
                }
                Xml::Element cyl = geo.getOptionalElement("cylinder");
                if (cyl.isValid()) {
                    Real r = cyl.getRequiredElementValueAs<Real>("radius");
                    Real l = cyl.getRequiredElementValueAs<Real>("length");
                    // Cylinder is along Z in Gazebo, Y in Simbody
                    Rotation YtoZ(Pi/2, XAxis);
                    mobod.addBodyDecoration(
                        Transform(X_LV.R()*YtoZ, X_LV.p()),
                        DecorativeCylinder(r, l/2).setOpacity(.5)
                                .setColor(color).setScale(scale));
                }

            } 

            // COLLISION
            Array_<Xml::Element> coll = master.getAllElements("collision");
            const Vec3 collColor = model.isStatic ? Green : Gray;
            for (unsigned i=0; i < coll.size(); ++i) {
                Transform X_LC = getPose(coll[i]);
                Xml::Element geo = coll[i].getRequiredElement("geometry");

                // Model sphere collision surface.
                Xml::Element sphere = geo.getOptionalElement("sphere");
                if (sphere.isValid()) {
                    Real r = sphere.getRequiredElementValueAs<Real>("radius");
                    mobod.addBodyDecoration(
                        X_LC, DecorativeSphere(r)
                                .setRepresentation(DecorativeGeometry::DrawWireframe)
                                .setColor(collColor));
                    ContactSurface surface(ContactGeometry::Sphere(r),
                                           material);
                    if (!gzOutb.selfCollide)
                        surface.joinClique(model.modelClique);
                    mobod.updBody().addContactSurface(X_LC, surface);
                }

                // Model cylinder collision surface (must fake with ellipsoid).
                Xml::Element cyl = geo.getOptionalElement("cylinder");
                if (cyl.isValid()) {
                    Real r   = cyl.getRequiredElementValueAs<Real>("radius");
                    Real len = cyl.getRequiredElementValueAs<Real>("length");
                    // Cylinder is along Z in Gazebo
#ifndef USE_CONTACT_MESH
                    Vec3 esz = Vec3(r,r,len/2); // Use ellipsoid instead
                    mobod.addBodyDecoration(X_LC, 
                        DecorativeEllipsoid(esz)
                            .setRepresentation(DecorativeGeometry::DrawWireframe)
                            .setColor(collColor));
                    ContactSurface surface(ContactGeometry::Ellipsoid(esz),
                                           material);
#else
                    const int resolution = 0; // chunky hexagonal shape
                    const PolygonalMesh mesh = PolygonalMesh::
                        createCylinderMesh(ZAxis,r,len/2,resolution);
                    const ContactGeometry::TriangleMesh triMesh(mesh);
    
                    mobod.addBodyDecoration(X_LC, 
                        DecorativeMesh(triMesh.createPolygonalMesh())
                        .setRepresentation(DecorativeGeometry::DrawWireframe)
                        .setColor(collColor));
                    ContactSurface surface(triMesh, material,1 /*Thickness*/);
#endif
                    if (!gzOutb.selfCollide)
                        surface.joinClique(model.modelClique);
                    mobod.updBody().addContactSurface(X_LC, surface);
                }


                // Model box collision surface (must fake with ellipsoid).
                Xml::Element box = geo.getOptionalElement("box");
                if (box.isValid()) {
                    Vec3 hsz = box.getRequiredElementValueAs<Vec3>("size")/2;
#ifndef USE_CONTACT_MESH
                    mobod.addBodyDecoration(X_LC, 
                        DecorativeEllipsoid(hsz) // use half dimensions
                            .setRepresentation(DecorativeGeometry::DrawWireframe)
                            .setColor(collColor));
                    ContactSurface surface(ContactGeometry::Ellipsoid(hsz),
                                           material);
#else
                    const int resolution = 20; // need dense to get near corners
                    const PolygonalMesh mesh = PolygonalMesh::
                        createBrickMesh(hsz,resolution);
                    const ContactGeometry::TriangleMesh triMesh(mesh);
    
                    mobod.addBodyDecoration(X_LC, 
                        DecorativeMesh(triMesh.createPolygonalMesh())
                        .setRepresentation(DecorativeGeometry::DrawWireframe)
                        .setColor(collColor));
                    ContactSurface surface(triMesh, material,1 /*Thickness*/);
#endif
                    if (!gzOutb.selfCollide)
                        surface.joinClique(model.modelClique);
                    mobod.updBody().addContactSurface(X_LC, surface);
                }
            }
        }
    }

    // Weld the slaves to their masters.
    for (int lx=0; lx < model.links.size(); ++lx) {
        GazeboLinkInfo& link = model.links.updLink(lx);
        if (link.slaveMobods.empty()) continue;
        for (unsigned i=0; i < link.slaveMobods.size(); ++i) {
            Constraint::Weld weld(link.masterMobod, link.slaveMobods[i]);
            link.slaveWelds.push_back(weld); // in case we want to know later
        }
    }

    // Add the loop joints if any.
    for (int lcx=0; lcx < mbgraph.getNumLoopConstraints(); ++lcx) {
        const MultibodyGraphMaker::LoopConstraint& loop =
            mbgraph.getLoopConstraint(lcx);

        GazeboJointInfo& joint  = *(GazeboJointInfo*)loop.getJointRef();
        GazeboLinkInfo&  parent = *(GazeboLinkInfo*) loop.getParentBodyRef();
        GazeboLinkInfo&  child  = *(GazeboLinkInfo*) loop.getChildBodyRef();

        if (joint.type == "weld") {
            Constraint::Weld weld(parent.masterMobod, joint.X_PA, 
                                  child.masterMobod,  joint.X_CB);
            joint.constraint = weld;
        } else if (joint.type == "ball") {
            Constraint::Ball ball(parent.masterMobod, joint.X_PA.p(), 
                                  child.masterMobod,  joint.X_CB.p());
            joint.constraint = ball;
        } else if (joint.type == "free") {
            // A "free" loop constraint is no constraint at all so we can
            // just ignore it. It might be more convenient if there were
            // a 0-constraint Constraint::Free, just as there is a 0-mobility
            // MobilizedBody::Weld.
        } else
            throw std::runtime_error(
                "Unrecognized loop constraint type '" + joint.type + "'.");
    }
}


//==============================================================================
//                            RUN SIMULATION
//==============================================================================
// Run a simulation and extract some data from it to show how that can be done.
static void runSimulation(const MultibodySystem&          mbs,
                          const std::vector<GazeboModel>& gzModels) {
    // Create a Visualizer so we can see what's what.
    Visualizer viz(mbs);
    viz.setShowSimTime(true);

    // Initialize the system and obtain the default state.
    State state = mbs.realizeTopology();

// The ragdoll2 model needs some reasonable initial conditions in order to
// assemble its loop joints properly.
#ifdef RAGDOLL2_ICS
    gzModels[0].joints.getJoint("l_upper_arm_joint").mobod.setOneQ(state,0,-.4);
    gzModels[0].joints.getJoint("r_upper_arm_joint").mobod.setOneQ(state,0,-.4);
    gzModels[0].joints.getJoint("torso_joint").mobod.setOneQ(state,0,.3);
#endif

    viz.report(state);  // Show the initial state.

    // If there are any constraints due to topological loops in the graph
    // then the system might not be assembled yet. Use Simbody's Assembler
    // to find a state that assembles the system.
    // Note: assembly mail fail if the system starts out too far from
    // assembled.

    printf("INITIAL STATE err=%g -- ENTER to assemble\n",
        state.getQErr().norm());
    getchar();
    Assembler assembler(mbs);
    Real assemblyTol = assembler.assemble(state);
    viz.report(state);
    printf("ASSEMBLED to err=%g -- ENTER to simulate\n", 
        state.getQErr().norm());
    getchar();

    // Assembled. Now simulate.

    const Real RunTime = 20;
    const Real ReportTime = 1./30.;
    //const Real Accuracy = 0.001; // 0.1% is Simbody default
    //const Real Accuracy = 0.05; // 5%
    const Real Accuracy = 0.1; // 1%
    //const Real Accuracy = 0.10; // 10%
    //const Real Accuracy = 0.20; // 20%

    // Use a low order integrator if you force the max step size to be small.
    //const Real MaxStepSize = .001;
    //RungeKutta2Integrator integ(mbs);
    const Real MaxStepSize = Infinity;
    RungeKutta3Integrator integ(mbs);
    //RungeKuttaMersonIntegrator integ(mbs);    // 4th order
    //CPodesIntegrator integ(mbs); // implicit integrator

    integ.setAccuracy(Accuracy);
    integ.setConstraintTolerance(std::min(1e-3, Accuracy/10)); 

    integ.initialize(state);
    viz.report(integ.getState());
    cout << "q=" << state.getQ() << endl;
    double startReal = realTime();
    double startCPU = cpuTime();
    while (integ.getTime() < RunTime) {
        Real nextReport = std::min(RunTime, integ.getTime() + ReportTime);
        do {
            integ.stepTo(nextReport, integ.getTime()+MaxStepSize);
        } while (integ.getTime() < nextReport);

        const State& state = integ.getState();
        const Real t = state.getTime();

        #ifdef ANIMATE // suppress for more accurate CPU time measurement
        viz.report(state);
        #endif

        // Show body origin locations and joint angles and rates at integer 
        // times to demo data extraction.
        if (std::abs(std::floor(t+ReportTime/2)-t) > ReportTime/2) 
            continue;

        for (unsigned m=0; m < gzModels.size(); ++m) {
            const GazeboModel& model = gzModels[m];
            printf("\nMODEL %s --------------------\n", model.name.c_str());
            printf("\n  LINKS t=%g\n", t);
            for (int i=0; i < model.links.size(); ++i) {
                const GazeboLinkInfo& link = model.links.getLink(i);
                const Vec3& loc = link.masterMobod.getBodyOriginLocation(state);
                printf("  %20s %10.3g %10.3g %10.3g\n", 
                    link.name.c_str(), loc[0], loc[1], loc[2]); 
            }
            printf("\n  JOINTS t=%g\n", t);
            for (int i=0; i < model.joints.size(); ++i) {
                const GazeboJointInfo& joint = model.joints.getJoint(i);
                if (joint.mobod.getNumU(state)) {
                    const Real q0 = joint.mobod.getOneQ(state, 0);
                    const Real u0 = joint.mobod.getOneU(state, 0);
                    printf("  %20s %10.3g %10.3g\n", 
                        joint.name.c_str(), q0, u0); 
                } else { // weld joint
                    printf("  %20s weld; no dofs\n", 
                        joint.name.c_str());
                }
            }
        }

    }
    viz.report(integ.getState()); // show final state

    // Finished simulating; dump out some stats. On Windows CPUtime is not
    // reliable if very little time is spent in this thread.
    const double timeInSec = realTime()-startReal;
    const double cpuInSec = cpuTime()-startCPU;
    const int evals = integ.getNumRealizations();
    cout << "Done -- took " << integ.getNumStepsTaken() << " steps in " <<
        timeInSec << "s for " << integ.getTime() << "s sim (avg step=" 
        << (1000*integ.getTime())/integ.getNumStepsTaken() << "ms) " 
        << (1000*integ.getTime())/evals << "sim ms/eval\n";
    cout << "CPUtime (not reliable when visualizing) " << cpuInSec << endl;

    printf("Used Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), 
        integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), 
        integ.getNumProjections());
}


//==============================================================================
//                    IMPLEMENTATIONS OF GAZEBO CLASSES
//==============================================================================

int GazeboLinks::addLink(const GazeboLinkInfo& info) {
    if (name2index.find(info.name) != name2index.end())
        throw std::runtime_error("GazeboLinks::addLink(): Link name '" 
            + info.name + " was used more than once.");

    const int lx = (int)linkByIndex.size();
    linkByIndex.push_back(info);
    name2index[info.name] = lx;
    return lx;
}


int GazeboJoints::addJoint(const GazeboJointInfo& info) {
    if(name2index.find(info.name) != name2index.end())
        throw std::runtime_error(
            "GazeboJoints::addJoint(): Joint name '" + info.name 
            + "' was used more than once.");

    const int jx = (int)jointByIndex.size();
    jointByIndex.push_back(info);
    name2index[info.name] = jx;
    return jx;
}



void GazeboModel::readModel(Xml::Element modelElt) {
    name = modelElt.getRequiredAttributeValue("name");
    X_WM = getPose(modelElt);
    isStatic = modelElt.getOptionalElementValueAs<bool>("static", false);

    Array_<Xml::Element> linkElts = modelElt.getAllElements("link");
    Array_<Xml::Element> jointElts = modelElt.getAllElements("joint");
    printf("Reading %s model '%s' with ground + %d links, %d joints.\n",
        isStatic ? "STATIC" : "DYNAMIC", name.c_str(), 
        linkElts.size(), jointElts.size());
    cout << "  Model frame X_WM as pose=" << transform2Pose(X_WM) << endl;

    // Create a World link, then read in the real links.
    GazeboLinkInfo worldLink("world");
    worldLink.defX_ML = ~X_WM; // different for each model in this world
    worldLink.massProps = MassProperties(Infinity,Vec3(0),UnitInertia(1));
    links.addLink(worldLink);
    for (unsigned i=0; i < linkElts.size(); ++i) {
        Xml::Element elt = linkElts[i];
        GazeboLinkInfo linfo(elt.getRequiredAttributeValue("name"));
        linfo.mustBeBaseLink = elt.getOptionalElementValueAs<bool>
                                                ("must_be_base_link", false);
        linfo.selfCollide = elt.getOptionalElementValueAs<bool>
                                                ("self_collide", false);
        linfo.element = elt;
        linfo.defX_ML = getPose(elt); // default link pose in Model frame
        linfo.massProps = getMassProperties(elt);
        links.addLink(linfo);
    }

    // Read in the joints.
    for (unsigned i=0; i < jointElts.size(); ++i) {
        Xml::Element elt = jointElts[i];
        const std::string name = elt.getRequiredAttributeValue("name");
        const std::string type = elt.getRequiredAttributeValue("type");
        GazeboJointInfo jinfo(name,type);
        jinfo.mustBreakLoopHere = elt.getOptionalElementValueAs<bool>
                                                ("must_be_loop_joint", false);
        jinfo.element = elt;

        // We'll accept either old style:
        //    <pose>pose on child</pose>
        //    <child>childName</child>
        //    <parent>parentName</parent>
        // or new style
        //    <child>
        //      <link_name>childName</linkName>
        //      <pose>childPose</pose>
        //    </child>
        //    <parent>
        //      <link_name>parentName</linkName>
        //      <pose>parentPose</pose>
        //    </parent>
        std::string child, parent;
        jinfo.X_CB = getPose(elt); // default joint frame on child body
        Xml::Element cElt = elt.getRequiredElement("child");
        if (cElt.hasElement("link_name")) {
            child = cElt.getRequiredElementValue("link_name");
            if (cElt.hasElement("pose")) 
                jinfo.X_CB = getPose(cElt);
            // else leave it as it was
        } else child = cElt.getValue(); // old style
        // At this point we have the child frame in X_CB

        Xml::Element pElt = elt.getRequiredElement("parent");
        bool gotParentPose = false;
        if (pElt.hasElement("link_name")) {
            parent = pElt.getRequiredElementValue("link_name");
            if (pElt.hasElement("pose")) {
                jinfo.X_PA = getPose(pElt);
                gotParentPose = true;
            }
        } else parent = pElt.getValue(); // old style

        assert(links.hasLink(parent) && links.hasLink(child));
        assert(child != parent);

        jinfo.parent = parent;
        jinfo.child  = child;
        if (!gotParentPose) {
            // Must derive joint frame on parent using the default body poses
            // TODO: this should be fixed in the sdf file
            const Transform X_MC = links.getLink(child).defX_ML;
            const Transform X_MP = links.getLink(parent).defX_ML;
            const Transform X_PC = ~X_MP*X_MC;
            jinfo.X_PA = X_PC*jinfo.X_CB; // i.e., A spatially coincident with B
        }

        joints.addJoint(jinfo);
    }
}

