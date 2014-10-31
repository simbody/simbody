/* -------------------------------------------------------------------------- *
 *             Simbody(tm) Example: Boston Dynamics Atlas robot               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors: Jack Wang, Chris Dembia, John Hsu                            *
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

#include "Simbody.h"
#include "URDFReader.h"
#include "Atlas.h"

#include <cstdio>
#include <iostream>

using namespace SimTK;

namespace {
// This local event handler is called periodically to sample the robot's 
// sensors.
class AtlasJointSampler : public PeriodicEventHandler {
public:
    AtlasJointSampler(const Atlas& realRobot, 
                      Real          interval) 
    :   PeriodicEventHandler(interval), m_realRobot(realRobot) {}

    // This method is called whenever this event occurs.
    void handleEvent(State& state, Real, bool&) const OVERRIDE_11;

private:
    const Atlas&                m_realRobot;
    Random::Gaussian            m_gaussian; // mean=0, stddev=1
};
}


static void createMultibodyGraph(URDFRobot&             robot,
                                 MultibodyGraphMaker&   mbgraph);

static void addRobotToSimbodySystem(const MultibodyGraphMaker& mbgraph,
                                    URDFRobot&                 robot,
                                    MultibodySystem&           mbs,
                                    SimbodyMatterSubsystem&    matter,
                                    GeneralForceSubsystem&     forces,
                                    CompliantContactSubsystem& contact);

//------------------------------------------------------------------------------
//                            ATLAS CONSTRUCTOR
//------------------------------------------------------------------------------
// Build a Simbody System of the Boston Dynamics Atlas robot.
Atlas::Atlas(const std::string& fileNameAndExt)
:   m_matter(*this), m_forces(*this), m_tracker(*this), 
    m_contact(*this, m_tracker),
    m_sampledAngles(*this, Stage::Dynamics, Vector()),
    m_sampledRates(*this, Stage::Dynamics, Vector()),
    m_sampledPelvisPose(*this, Stage::Dynamics, Transform()),
    m_sampledEndEffectorPos(*this, Stage::Dynamics, Vec3(0)),
    m_qNoise(*this,Stage::Dynamics,Zero), m_uNoise(*this,Stage::Dynamics,Zero),
    m_endEffectorLinkName("r_hand"), m_endEffectorStation(0,-.15,0)
{
    const std::string urdfPathName = "models/" + fileNameAndExt;

    setUpDirection(ZAxis);
    m_matter.setShowDefaultGeometry(false);

    // Set the sensor sampling rate. TODO: should be settable.
    addEventHandler(new AtlasJointSampler(*this, 0.002));  

    //--------------------------------------------------------------------------
    //                          Read the robot file
    //--------------------------------------------------------------------------
    std::cout << "Reading file: " << urdfPathName << std::endl;
    Xml::Document urdf(urdfPathName);

    if (urdf.getRootTag() != "robot")
        throw std::runtime_error
           ("Expected to see document tag <robot>but saw <"
            + urdf.getRootTag() + "> instead.");

    // This is a URDF robot document.
    Xml::Element root = urdf.getRootElement();
    std::cout << "Reading robot '" 
              << root.getOptionalAttributeValue("name","NONAME") << "'\n";

    m_urdfRobot.readRobot(root);
    createMultibodyGraph(m_urdfRobot, m_mbGraph);
    // Optional: dump the graph to stdout for debugging or curiosity.
    //m_mbGraph.dumpGraph(std::cout);


    //--------------------------------------------------------------------------
    //                          Gravity and Ground
    //--------------------------------------------------------------------------
    m_gravity = Force::Gravity(m_forces, m_matter, -SimTK::ZAxis, 9.80665);

    // Define a material to use for ground contact. This is not very stiff.
    ContactMaterial groundMaterial(1e6,   // stiffness
                             0.1,  // dissipation
                             0.7,   // mu_static
                             0.5,   // mu_dynamic
                             0.5);  // mu_viscous
    // Add a contact surface to represent the ground.
    // Half space normal is -x; must rotate about y to make it +z.
    m_matter.Ground().updBody().addContactSurface(Rotation(Pi/2,YAxis),
       ContactSurface(ContactGeometry::HalfSpace(), groundMaterial));
    // Draw world frame.
    m_matter.Ground().addBodyDecoration(Vec3(0), 
        DecorativeFrame(1).setColor(Green).setLineThickness(3)); // World

    //--------------------------------------------------------------------------
    //                          Build Simbody System
    //--------------------------------------------------------------------------
    addRobotToSimbodySystem(m_mbGraph, m_urdfRobot, 
                            *this, m_matter, m_forces, m_contact);
}


//------------------------------------------------------------------------------
//                 ATLAS JOINT SAMPLER :: HANDLE EVENT
//------------------------------------------------------------------------------
// This method is called whenever this event occurs.
void AtlasJointSampler::handleEvent
    (State& state, Real, bool&) const 
{
    const int nq = state.getNQ(), nu = state.getNU();
    const Vector& q = state.getQ();
    const Vector& u = state.getU();
    const Real qNoise = m_realRobot.getAngleNoise(state);
    const Real uNoise = m_realRobot.getRateNoise(state);
    Vector qSample(nq), uSample(nu);
    for (int i=0; i<nq; ++i)
        qSample[i] = q[i] + qNoise * m_gaussian.getValue();
    for (int i=0; i<nu; ++i)
        uSample[i] = u[i] + uNoise * m_gaussian.getValue();
    m_realRobot.setSampledAngles(state, qSample);
    m_realRobot.setSampledRates(state, uSample);
    m_realRobot.setSampledPelvisPose(state, 
        m_realRobot.getBody("pelvis").getBodyTransform(state));
    m_realRobot.setSampledEndEffectorPos(state, 
        m_realRobot.getActualEndEffectorPosition(state));
}


//==============================================================================
//                           CREATE MULTIBODY GRAPH
//==============================================================================
// Define URDF joint types, then use links and joints in the given model
// to construct a reasonable spanning-tree-plus-constraints multibody graph
// to represent that model. An exception will be thrown if this fails.
// Note that this step is not Simbody dependent.
static void createMultibodyGraph(URDFRobot&             robot,
                                 MultibodyGraphMaker&   mbgraph) {
    // Step 1: Tell MultibodyGraphMaker about joints it should know about.

    mbgraph.setWeldJointTypeName("fixed");      // 0 dofs, Weld
    mbgraph.setFreeJointTypeName("floating");   // 6 dofs, Free

    //                   URDF name  #dofs     Simbody equivalent
    mbgraph.addJointType("revolute",   1);   // Pin with limits
    mbgraph.addJointType("continuous", 1);   // Pin with no limits
    mbgraph.addJointType("planar",     3);   // Planar

    // Step 2: Tell it about all the links we read from the input file, 
    // starting with world, and provide a reference pointer.
    for (int lx=0; lx < robot.links.size(); ++lx) {
        URDFLinkInfo& link = robot.links.updLink(lx);
        mbgraph.addBody(link.name, link.massProps.getMass(), 
                        link.mustBeBaseLink, &link);
    }

    // Step 3: Tell it about all the joints we read from the input file,
    // and provide a reference pointer.
    for (int jx=0; jx < robot.joints.size(); ++jx) {
        URDFJointInfo& joint = robot.joints.updJoint(jx);
        mbgraph.addJoint(joint.name, joint.type, joint.parent, joint.child, 
                            joint.mustBreakLoopHere, &joint);
    }

    // Setp 4. Generate the multibody graph.
    mbgraph.generateGraph();
}

void Atlas::initialize(SimTK::State& state) {
    state = realizeTopology();
    realizeModel(state);
    const int nq=state.getNQ(), nu=state.getNU();
    m_sampledAngles.setValue(state, SimTK::Vector(nq, SimTK::Zero));
    m_sampledRates.setValue(state, SimTK::Vector(nu, SimTK::Zero));

    m_effortLimits.resize(nu); m_effortLimits =  Infinity;
    m_speedLimits.resize(nu);  m_speedLimits  =  Infinity;
    m_lowerLimits.resize(nq);  m_lowerLimits  = -Infinity;
    m_upperLimits.resize(nq);  m_upperLimits  =  Infinity;

    for (int j=0; j < (int)m_urdfRobot.joints.size(); ++j) {
        const URDFJointInfo& joint = m_urdfRobot.joints.getJoint(j);
        const MobilizedBody& mobod = joint.mobod;
        const int nu = mobod.getNumU(state), nq = mobod.getNumQ(state);
        const int u0 = mobod.getFirstUIndex(state);
        const int q0 = mobod.getFirstQIndex(state);
        for (int u=0; u < nu; ++u) {
            m_effortLimits[u0+u] = joint.effort;
            m_speedLimits[u0+u] = joint.velocity;
        }
        for (int q=0; q < nq; ++q) {
            m_lowerLimits[q0+q] = joint.lower;
            m_upperLimits[q0+q] = joint.upper;
        }
    }
}



//==============================================================================
//                       ADD ROBOT TO SIMBODY SYSTEM
//==============================================================================
// Given a desired multibody graph, gravity, and the URDF robot that was
// used to generate the graph, add elements to the Simbody System to represent
// it. There are many limitations here, especially in the handling of contact. 
// Any URDF features that we haven't roboted are just ignored.
// The URDFRobot is updated so that its links and joints have references to
// their corresponding Simbody elements.

static void addRobotToSimbodySystem(const MultibodyGraphMaker& mbgraph,
                                    URDFRobot&                 robot,
                                    MultibodySystem&           mbs,
                                    SimbodyMatterSubsystem&    matter,
                                    GeneralForceSubsystem&     forces,
                                    CompliantContactSubsystem& contact) 
{
    const std::string& freeJointName = mbgraph.getFreeJointTypeName();
    const std::string& weldJointName = mbgraph.getWeldJointTypeName();

    // Define a material to use for contact. This is not very stiff.
    ContactMaterial material(1e6,   // stiffness
                             0.1,  // dissipation
                             0.7,   // mu_static
                             0.5,   // mu_dynamic
                             0.5);  // mu_viscous

    // Draw the robot frame.
    matter.Ground().addBodyDecoration(robot.X_WM, 
        DecorativeFrame(.75).setColor(Orange).setLineThickness(3));  // Model
    matter.Ground().addBodyDecoration(robot.X_WM.p()+Vec3(0,0,.1),
        DecorativeText(robot.name).setScale(.2)
        .setColor(robot.isStatic?Green:Cyan));

    // Generate a contact clique we can put collision geometry in to prevent
    // self-collisions.
    robot.robotClique = ContactSurface::createNewContactClique();

    // Record the MobilizedBody for the World link.
    robot.links.updLink(0).masterMobod = matter.Ground();

    // Run through all the mobilizers in the multibody graph, adding a Simbody
    // MobilizedBody for each one. Also add visual and collision geometry to the
    // bodies when they are mobilized.
    for (int mobNum=0; mobNum < mbgraph.getNumMobilizers(); ++mobNum) {
        // Get a mobilizer from the graph, then extract its corresponding
        // joint and bodies. Note that these don't necessarily have equivalents
        // in the URDFLink and URDFJoint inputs.
        const MultibodyGraphMaker::Mobilizer& mob = mbgraph.getMobilizer(mobNum);
        const std::string& type = mob.getJointTypeName();

        // The inboard body always corresponds to one of the input links,
        // because a slave link is always the outboard body of a mobilizer.
        // The outboard body may be slave, but its master body is one of the
        // URDF input links.
        const bool isSlave = mob.isSlaveMobilizer();
        URDFLinkInfo& gzInb  = *(URDFLinkInfo*)mob.getInboardBodyRef();
        URDFLinkInfo& gzOutb = *(URDFLinkInfo*)mob.getOutboardMasterBodyRef();

        const MassProperties massProps = 
            gzOutb.getEffectiveMassProps(mob.getNumFragments());

        //std::cout << "link=" << gzOutb.name << std::endl;
        //std::cout << "     " << massProps << std::endl;

        // This will reference the new mobilized body once we create it.
        MobilizedBody mobod; 

        if (robot.isStatic) {
            mobod = matter.updGround();
        } else if (mob.isAddedBaseMobilizer()) {
            // There is no corresponding URDF joint for this mobilizer.
            // Create the joint and set its default position to be the default
            // pose of the base link relative to the Ground frame.
            assert(type==freeJointName); // May add more types later
            if (type == freeJointName) {
                MobilizedBody::Free freeJoint(
                    gzInb.masterMobod,  Transform(),
                    massProps,    Transform());
                freeJoint.setDefaultTransform(Transform());
                mobod = freeJoint;
            }
        } else {
            // This mobilizer does correspond to one of the input joints.
            URDFJointInfo& gzJoint = *(URDFJointInfo*)mob.getJointRef();
            const bool isReversed = mob.isReversedFromJoint();

            // Find inboard and outboard frames for the mobilizer; these are
            // parent and child frames or the reverse.

            const Transform& X_IF0 = isReversed ? gzJoint.X_CB : gzJoint.X_PA;
            const Transform& X_OM0 = isReversed ? gzJoint.X_PA : gzJoint.X_CB;

            const MobilizedBody::Direction direction =
                isReversed ? MobilizedBody::Reverse : MobilizedBody::Forward;

            // Here we are using a Simbody Bushing joint to implement "floating"
            // so that we'll have Euler angle rotations & derivatives rather
            // than the quaternions we would get with a Simbody Free joint.
            // Of course that does mean there will be a singularity somewhere.
            // This makes it easier to apply prescribed motion. 
            if (type == freeJointName) {
                MobilizedBody::Bushing freeJoint(
                    gzInb.masterMobod,  X_IF0,
                    massProps,          X_OM0, 
                    direction);
                Transform defX_FM = isReversed ? Transform(~gzJoint.defX_AB)
                                               : gzJoint.defX_AB;
                freeJoint.setDefaultTransform(defX_FM);
                mobod = freeJoint;

            // Pin joint with no limits
            } else if (type == "continuous") {
                Xml::Element axisElt = gzJoint.element.getRequiredElement("axis");
                UnitVec3 axis = 
                    UnitVec3(axisElt.getRequiredAttributeValueAs<Vec3>("xyz")); 
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
            // Pin joint with limits
            } else if (type == "revolute") {
                Xml::Element axisElt = gzJoint.element.getRequiredElement("axis");
                UnitVec3 axis = 
                    UnitVec3(axisElt.getRequiredAttributeValueAs<Vec3>("xyz")); 
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
                    UnitVec3(axisElt.getRequiredAttributeValueAs<Vec3>("xyz")); 
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
            } else if (type == weldJointName) {
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
                Transform X_LV = URDF::getOrigin(visuals[i]);
                Xml::Element geo = visuals[i].getRequiredElement("geometry");
                Xml::Element box = geo.getOptionalElement("box");
                if (box.isValid()) {
                    Vec3 sz = box.getRequiredAttributeValueAs<Vec3>("size");
                    mobod.addBodyDecoration(
                        X_LV, DecorativeBrick(sz/2).setOpacity(.5)
                                .setColor(color).setScale(scale));
                }
                Xml::Element sphere = geo.getOptionalElement("sphere");
                if (sphere.isValid()) {
                    Real r = sphere.getRequiredAttributeValueAs<Real>("radius");
                    mobod.addBodyDecoration(
                        X_LV, DecorativeSphere(r).setOpacity(.5)
                                .setColor(color).setScale(scale));
                }
                Xml::Element cyl = geo.getOptionalElement("cylinder");
                if (cyl.isValid()) {
                    Real r = cyl.getRequiredAttributeValueAs<Real>("radius");
                    Real l = cyl.getRequiredAttributeValueAs<Real>("length");
                    // Cylinder is along Z in URDF, Y in Simbody
                    Rotation YtoZ(Pi/2, XAxis);
                    mobod.addBodyDecoration(
                        Transform(X_LV.R()*YtoZ, X_LV.p()),
                        DecorativeCylinder(r, l/2).setOpacity(.5)
                                .setColor(color).setScale(scale));
                }
                Xml::Element meshFile = geo.getOptionalElement("mesh");
                if (meshFile.isValid()) {
                    std::string pathname = 
                        meshFile.getRequiredAttributeValue("filename");
                    const std::string::size_type spos = pathname.rfind('/');
                    if (spos != std::string::npos)
                        pathname = pathname.substr(spos+1);
                    const Vec3 scale = 
                        meshFile.getOptionalAttributeValueAs<Vec3>
                                                            ("scale", Vec3(1));
                    bool isAbsolutePath; std::string dir, fn, ext;
                    Pathname::deconstructPathname(pathname, isAbsolutePath,
                                                  dir, fn, ext);
                    PolygonalMesh polyMesh;
                    polyMesh.loadStlFile(dir + "geometry/" + fn + ".stl");
                    DecorativeMesh decMesh(polyMesh);
                    decMesh.setColor(Cyan).setScaleFactors(scale);
                    mobod.addBodyDecoration(X_LV, decMesh);
                }
            } 

            // COLLISION
            Array_<Xml::Element> coll = master.getAllElements("collision");
            const Vec3 collColor = robot.isStatic ? Green : Gray;
            //TODO: disabled
            for (unsigned i=0; i < 0*coll.size(); ++i) {
                Transform X_LC = URDF::getOrigin(coll[i]);
                Xml::Element geo = coll[i].getRequiredElement("geometry");

                // Model sphere collision surface.
                Xml::Element sphere = geo.getOptionalElement("sphere");
                if (sphere.isValid()) {
                    Real r = sphere.getRequiredAttributeValueAs<Real>("radius");
                    mobod.addBodyDecoration(
                        X_LC, DecorativeSphere(r)
                                .setRepresentation(DecorativeGeometry::DrawWireframe)
                                .setColor(collColor));
                    ContactSurface surface(ContactGeometry::Sphere(r),
                                           material);
                    if (!gzOutb.selfCollide)
                        surface.joinClique(robot.robotClique);
                    mobod.updBody().addContactSurface(X_LC, surface);
                }

                // Model cylinder collision surface (must fake with ellipsoid).
                Xml::Element cyl = geo.getOptionalElement("cylinder");
                if (cyl.isValid()) {
                    Real r   = cyl.getRequiredAttributeValueAs<Real>("radius");
                    Real len = cyl.getRequiredAttributeValueAs<Real>("length");
                    // Cylinder is along Z in URDF
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
                        surface.joinClique(robot.robotClique);
                    mobod.updBody().addContactSurface(X_LC, surface);
                }


                // Model box collision surface (must fake with ellipsoid).
                Xml::Element box = geo.getOptionalElement("box");
                if (box.isValid()) {
                    Vec3 hsz = box.getRequiredAttributeValueAs<Vec3>("size")/2;
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
                        surface.joinClique(robot.robotClique);
                    mobod.updBody().addContactSurface(X_LC, surface);
                }
            }
        }
    }

    // Weld the slaves to their masters.
    for (int lx=0; lx < robot.links.size(); ++lx) {
        URDFLinkInfo& link = robot.links.updLink(lx);
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

        URDFJointInfo& joint  = *(URDFJointInfo*)loop.getJointRef();
        URDFLinkInfo&  parent = *(URDFLinkInfo*) loop.getParentBodyRef();
        URDFLinkInfo&  child  = *(URDFLinkInfo*) loop.getChildBodyRef();

        if (joint.type == weldJointName) {
            Constraint::Weld weld(parent.masterMobod, joint.X_PA, 
                                  child.masterMobod,  joint.X_CB);
            joint.constraint = weld;
        } else if (joint.type == "ball") {
            Constraint::Ball ball(parent.masterMobod, joint.X_PA.p(), 
                                  child.masterMobod,  joint.X_CB.p());
            joint.constraint = ball;
        } else if (joint.type == freeJointName) {
            // A "free" loop constraint is no constraint at all so we can
            // just ignore it. It might be more convenient if there were
            // a 0-constraint Constraint::Free, just as there is a 0-mobility
            // MobilizedBody::Weld.
        } else
            throw std::runtime_error(
                "Unrecognized loop constraint type '" + joint.type + "'.");
    }
}







