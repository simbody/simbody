#ifndef SIMTK_SUBSYSTEM_H_
#define SIMTK_SUBSYSTEM_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/** @file
 * User-visible, client side declaration of Subsystem.
 */

#include "SimbodyCommon.h"
#include "Placement.h"
#include "simtk/SimTK.h"
#include <iostream>

namespace simtk {

class Placement;
class PlacementValue;

// Declared below and elsewhere. Indenting shows inheritance structure.
class Subsystem;
class   Feature;
class       Station;
class       Direction;
class       Orientation;
class       FrameFeature;
class           Body;
class               RigidBody;
class               DeformableBody;
class       RealMeasure;
class           RealParameter;
class       Vec3Measure;
class           Vec3Parameter;
class       StationMeasure;
class           StationParameter;
class       DirectionMeasure;
class       OrientationMeasure;
class       MassElement;
class           PointMassElement;
class           SphereMassElement;
class           CylinderMassElement;
class           BrickMassElement;
class   Joint;
class       PinJoint;
class       BallJoint;
class   Multibody;

/**
 * Subsystems form a tree because many of them have sub-Subsystems (children).
 * Parent subsystems own their children; i.e., destructing the parent destructs
 * all the children.
 *
 * Low-level "elemental" subsystems are called Features. Features require
 * a Placement in order to be useful (e.g., a Station has to have a location).
 * A Feature's Placement may have a constant value, or can involve sibling,
 * parent, or ancestor Features but not its own children. 
 *
 * Every Subsystem has a name and an index by which it is known
 * to its parent. Parentless subsystems can exist; they still have a
 * name but they do not have an index. We call these "root" or "prototype"
 * subsystems depending on intended use.
 */
class Subsystem {
public:
    Subsystem() : rep(0) { }
    Subsystem(const Subsystem&);    // external placements are not copied or assigned
    Subsystem& operator=(const Subsystem&);
    ~Subsystem();

    // calculate values for all fully-resolved placements
    void realize(Stage) const;

    const Placement&      getPlacement()  const; // fails if Subsystem is not a Feature
    const PlacementValue& getValue()      const; //   "
    void place(const Placement&);                //   "
    void replace(const Placement&);              //   "

    // Read-only access to child subsystems.
    const Subsystem&        getSubsystem       (const String&) const; // generic
    const Feature&          getFeature         (const String&) const; // slightly less generic
    const RealParameter&    getRealParameter   (const String&) const; // type checked
    const Vec3Parameter&    getVec3Parameter   (const String&) const; //   "
    const StationParameter& getStationParameter(const String&) const;
    const RealMeasure&      getRealMeasure     (const String&) const;
    const Vec3Measure&      getVec3Measure     (const String&) const;
    const StationMeasure&   getStationMeasure  (const String&) const;
    const Station&          getStation         (const String&) const;
    const Direction&        getDirection       (const String&) const;
    const Orientation&      getOrientation     (const String&) const;
    const FrameFeature&     getFrame           (const String&) const;

    // Writable access to subsystems, e.g. allowing feature placement.
    Subsystem&              updSubsystem       (const String&);   // generic
    Feature&                updFeature         (const String&);   // less generic
    RealParameter&          updRealParameter   (const String&);   // type checked
    Vec3Parameter&          updVec3Parameter   (const String&);   //   "
    StationParameter&       updStationParameter(const String&);
    RealMeasure&            updRealMeasure     (const String&);
    Vec3Measure&            updVec3Measure     (const String&);
    StationMeasure&         updStationMeasure  (const String&);
    Station&                updStation         (const String&);
    Direction&              updDirection       (const String&);
    Orientation&            updOrientation     (const String&);
    FrameFeature&           updFrame           (const String&);

    // Create a new feature on this subsystem with a given name and type, and
    // optionally create a placement for it using the prototype placement supplied.
    RealParameter&    addRealParameter
                        (const String&, const Placement& = Placement());
    Vec3Parameter&    addVec3Parameter
                        (const String&, const Placement& = Placement());
    RealMeasure&      addRealMeasure
                        (const String&, const Placement& = Placement());
    Vec3Measure&      addVec3Measure
                        (const String&, const Placement& = Placement());
    StationParameter& addStationParameter
                        (const String&, const Placement& = Placement());
    StationMeasure&   addStationMeasure
                        (const String&, const Placement& = Placement());
    Station&          addStation
                        (const String&, const Placement& = Placement());
    Direction&        addDirection
                        (const String&, const Placement& = Placement());
    Orientation&      addOrientation
                        (const String&, const Placement& = Placement());
    FrameFeature&     addFrame
                        (const String&, const Placement& = Placement());

    // This is similar to the "add" routines above, except that the newly created
    // Subsystem is modeled on the prototype Subsystem supplied here.
    Subsystem&        addSubsystemLike(const Subsystem&, const String&);

    // This generic routine will create a Feature of the specific type supplied
    // in the Subsystem argument, which must turn out to be a Feature. If a 
    // Placement is provided, it must be appropriate for that kind of Feature.
    Feature&          addFeatureLike
                        (const Subsystem&,
                         const String&, const Placement& = Placement());
   
    bool              hasParentSubsystem() const;
    int               getIndexInParent()   const; // -1 if no parent
    const Subsystem&  getParentSubsystem() const;

    bool isSameSubsystem(const Subsystem& s) const;

    String getName()              const;
    String getFullName()          const; // "ancestors/parent/name"

    // Subsystems which are children of this subsystem
    int              getNSubsystems()  const;
    const Subsystem& getSubsystem(int) const;
    Subsystem&       updSubsystem(int);

    const Subsystem& operator[](int i) const           {return getSubsystem(i);}
    Subsystem&       operator[](int i)                 {return updSubsystem(i);}
    const Subsystem& operator[](const String& s) const {return getSubsystem(s);}
    Subsystem&       operator[](const String& s)       {return updSubsystem(s);}

    String toString(const String& linePrefix="") const;

    // For internal use only.
    bool                      hasRep() const            {return rep != 0;}
    const class SubsystemRep& getRep() const            {assert(rep); return *rep;}
    class SubsystemRep&       updRep()                  {assert(rep); return *rep;}
    void                      setRep(SubsystemRep* fp)  {assert(!rep); rep=fp;}
    void checkSubsystemConsistency(const Subsystem* expParent,
                                   int              expIndexInParent,
                                   const Subsystem& expRoot) const;
protected:
    class SubsystemRep* rep;
    friend class SubsystemRep;
};

std::ostream& operator<<(std::ostream& o, const Subsystem&);


} // namespace simtk

#endif // SIMTK_SUBSYSTEM_H_
