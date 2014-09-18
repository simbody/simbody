/* -------------------------------------------------------------------------- *
 *             Simbody(tm) Example: Universal Robotics UR10 Arm               *
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
#include "UR10.h"

using namespace SimTK;

UR10::UR10()
    : m_matter(*this), m_forces(*this)
{
    setUpDirection(ZAxis);
    m_matter.setShowDefaultGeometry(false);

    //--------------------------------------------------------------------------
    //                          Constants, etc.
    //--------------------------------------------------------------------------


    // Miscellaneous.
    // --------------



    //--------------------------------------------------------------------------
    //                          Gravity
    //--------------------------------------------------------------------------

    // Gravity.
    m_gravity = Force::Gravity(m_forces, m_matter, -SimTK::ZAxis, 9.8066);

    //--------------------------------------------------------------------------
    //                          Body information
    //--------------------------------------------------------------------------
    // Mass properties.
    Body baseInfo(MassProperties(4., Vec3(0),
                                 Inertia(0.0061063308908,
                                         0.0061063308908,
                                         0.01125)));

    Body shoulderInfo(MassProperties(7.778, Vec3(0),
                                     Inertia(0.0314743125769,
                                             0.0314743125769,
                                             0.021875625)));

    const Real    upperArmMass = 12.93;
    const Vec3    upperArmCOM(0,0,.306);
    const Inertia upperArmCentral(0.421753803798,
                                  0.421753803798,
                                  0.036365625);
    Body upperArmInfo(MassProperties(upperArmMass, upperArmCOM,
        upperArmCentral.shiftFromMassCenter(-upperArmCOM, upperArmMass)));

    const Real    forearmMass = 3.87;
    const Vec3    forearmCOM(0,0,0.28615);
    const Inertia forearmCentral(0.111069694097,
                                 0.111069694097,
                                 0.010884375);
    Body forearmInfo(MassProperties(forearmMass, forearmCOM,
        forearmCentral.shiftFromMassCenter(-forearmCOM, forearmMass)));

    Body wrist1Info(MassProperties(1.96, Vec3(0),
                                   Inertia(0.0051082479567,
                                           0.0051082479567,
                                           0.0055125)));
    Body wrist2Info = wrist1Info;

    Body wrist3Info(MassProperties(0.202, Vec3(0),
                                   Inertia(0.000526462289415,
                                           0.000526462289415,
                                           0.000568125)));

    const Vec3 eeHdims(.02,.02,.02); // cube
    Body endEffectorInfo(MassProperties(.1, Vec3(0),
                            UnitInertia::brick(eeHdims)));

    // Geometry
    PolygonalMesh baseMesh, shoulderMesh, upperArmMesh, forearmMesh,
                  wrist1Mesh, wrist2Mesh, wrist3Mesh;

    baseMesh.loadObjFile("geometry/Base.obj");
    shoulderMesh.loadObjFile("geometry/Shoulder.obj");
    upperArmMesh.loadObjFile("geometry/UpperArm.obj");
    forearmMesh.loadObjFile("geometry/Forearm.obj");
    wrist1Mesh.loadObjFile("geometry/Wrist1.obj");
    wrist2Mesh.loadObjFile("geometry/Wrist2.obj");
    wrist3Mesh.loadObjFile("geometry/Wrist3.obj");

    m_matter.updGround().addBodyDecoration(Vec3(0), 
               DecorativeFrame(0.5));

    baseInfo.addDecoration(DecorativeMesh(baseMesh).setColor(Gray));
    shoulderInfo.addDecoration(DecorativeMesh(shoulderMesh).setColor(Cyan));
    upperArmInfo.addDecoration(DecorativeMesh(upperArmMesh).setColor(Gray));
    forearmInfo.addDecoration(DecorativeMesh(forearmMesh)
                              .setColor(Gray).setOpacity(.3));
    wrist1Info.addDecoration(DecorativeMesh(wrist1Mesh).setColor(Cyan));
    wrist2Info.addDecoration(DecorativeMesh(wrist2Mesh).setColor(Gray));
    wrist3Info.addDecoration(DecorativeMesh(wrist3Mesh).setColor(Cyan));

    endEffectorInfo.addDecoration(DecorativeBrick(eeHdims)
                                  .setColor(Purple).setOpacity(.5));


    //--------------------------------------------------------------------------
    //                         Mobilized Bodies
    //--------------------------------------------------------------------------
    const Rotation ZtoY(-Pi/2, XAxis);
    const Rotation ZtoY90(BodyRotationSequence, -Pi/2, XAxis, -Pi/2, ZAxis);
    m_bodies[Ground] = m_matter.updGround();

    m_bodies[Base] = MobilizedBody::Weld(
        m_matter.updGround(), Vec3(0),
        baseInfo,             Vec3(0));

    m_bodies[Shoulder] = MobilizedBody::Pin( // shoulder_pan about Z
        m_bodies[Base], Vec3(0, 0, .1273),
        shoulderInfo,   Vec3(0));

    m_bodies[UpperArm] = MobilizedBody::Pin( // shoulder_lift about Y
        m_bodies[Shoulder], Transform(ZtoY90, Vec3(0, 0.220941, 0)),
        upperArmInfo,       ZtoY);

    m_bodies[Forearm] = MobilizedBody::Pin( // elbow about Y
        m_bodies[UpperArm], Transform(ZtoY, Vec3(0, -0.1719, 0.612)),
        forearmInfo,        ZtoY);

    m_bodies[Wrist1] = MobilizedBody::Pin( // wrist1 about Y
        m_bodies[Forearm], Transform(ZtoY90, Vec3(0, 0, 0.5723)),
        wrist1Info,        ZtoY);

    m_bodies[Wrist2] = MobilizedBody::Pin( // wrist2 about Z
        m_bodies[Wrist1], Vec3(0, 0.1149, 0),
        wrist2Info,       Vec3(0));

    m_bodies[Wrist3] = MobilizedBody::Pin( // wrist3 about Y
        m_bodies[Wrist2], Transform(ZtoY, Vec3(0, 0, 0.1157)),
        wrist3Info,       ZtoY);

    m_bodies[EndEffector] = MobilizedBody::Weld(
        m_bodies[Wrist3],     Vec3(0, 0.1149, 0),
        endEffectorInfo,      Vec3(0));

    //TODO: joint stops

}




