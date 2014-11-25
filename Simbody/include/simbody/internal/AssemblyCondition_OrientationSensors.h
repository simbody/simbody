#ifndef SimTK_SIMBODY_ASSEMBLY_CONDITION_ORIENTATION_SENSORS_H_
#define SimTK_SIMBODY_ASSEMBLY_CONDITION_ORIENTATION_SENSORS_H_

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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Assembler.h"
#include "simbody/internal/AssemblyCondition.h"

#include <map>

namespace SimTK {

//------------------------------------------------------------------------------
//                           ORIENTATION SENSORS
//------------------------------------------------------------------------------
/** This AssemblyCondition specifies a correspondence between orientation
sensors fixed on mobilized bodies ("osensors") and Ground-relative orientation 
sensor readings ("observations"). Each osensor represents the three orthogonal
axis directions of a frame fixed to a mobilized body; only orientation and not
location is represented. Orientation sensing is commonly performed by
Inertial Measurements Units (IMUs). The behavior here is designed to be as
similar as possible to position sensors handled by the Markers 
AssemblyCondition. You can combine OrientationSensors and Markers on the same
body.

The idea is to adjust the q's so that each osensor is oriented similarly to its 
corresponding observation. This is normally used as a goal since we don't expect
to be able to obtain a perfect match, but you can use these as a set of 
assembly error conditions if there are enough degrees of freedom to achieve a 
near-perfect solution. 

Osensors are defined one at a time and assigned sequential index values
of type OrientationSensors::OSensorIx. They may optionally be given unique, 
case-sensitive names, and we will keep a map from name to OSensorIx. A default 
name will be assigned if none is given. A weight is assigned to every osensor, 
with default weight=1. We do not expect that all the osensors will be used; 
osensors with weights of zero will not be included in the study, nor will 
osensors for which no observation is given.

Once specified, the osensor definitions do not change during a series of inverse
kinematic (tracking) steps. The observations, on the other hand, are expected 
to come from a time series of experimental measurements of osensor orientations
and will be different at every step. They typically come from a file organized 
by "frame", meaning an observation time and a set of observed orientations, one
per sensor, corresponding to that time. During initial setup, the number of 
observations per frame and their correspondence to the defined osensors is 
specified. They can be in any order, may skip some osensors, and may include 
data for osensors that are not defined. However, once initialized each frame 
must supply the same information in the same order. Data for an unobserved 
osensor can be provided as NaN in which case it will be ignored in that frame. 
The frame time is supplied to the track() method which initiates assembly for 
a frame.

Observation-osensor correspondence maps an ObservationIx to a unique OSensorIx. 
By default, we'll expect to get an observation for each osensor and that the 
observation order and the osensor order are the same, i.e. 
ObservationIx==OSensorIx for every osensor. However, you can instead define 
observation/osensor correspondence yourself, (\e after all osensors have been 
defined), via one of the defineObservationOrder() methods. This is done by 
supplying an array of OSensorIx values, or an array of osensor names, with the 
array elements ordered by ObservationIx. Any invalid osensor index or 
unrecognized osensor name means we will ignore values provide for that 
observation; similarly, any osensors whose index or name is not specified at 
all will be ignored. 
**/
class SimTK_SIMBODY_EXPORT OrientationSensors : public AssemblyCondition {

// This is a private class used in the implementation below but not
// accessible through the API.
struct OSensor {
    OSensor(const String& name, MobilizedBodyIndex bodyB, 
            const Rotation& orientationInB, Real weight = 1)
    :   name(name), bodyB(bodyB), orientationInB(orientationInB), weight(weight) 
    { assert(weight >= 0); }

    OSensor(MobilizedBodyIndex bodyB, const Rotation& orientationInB, 
            Real weight=1)
    :   name(""), bodyB(bodyB), orientationInB(orientationInB), weight(weight) 
    { assert(weight >= 0); }

    String              name;
    MobilizedBodyIndex  bodyB;
    Rotation            orientationInB;
    Real                weight; 
};

public:

/** Define the OSensorIx type which is just a uniquely-typed int. **/
SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(OrientationSensors,OSensorIx);
/** Define the ObservationIx type which is just a uniquely-typed int. **/
SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(OrientationSensors,ObservationIx);



//------------------------------------------------------------------------------
/** @name                Construction and setup
These methods are used as an extended construction phase for OrientationSensors
objects, defining the osensors and observations that will be used in the
subsequent tracking steps. **/
/*@{*/

/** The default constructor creates an empty OrientationSensors 
AssemblyCondition object that should be filled in with calls to addOSensor() 
and optionally defineObservationOrder(). **/
OrientationSensors() : AssemblyCondition("OrientationSensors") {}

/** Define a new orientation sensor (osensor) attached to a particular 
MobilizedBody. Note that an osensor will be ignored unless an observation is 
provided for it.
@param[in]      name
    A unique name to be used to identify this osensor. If the name is
    empty or blank, a default name will be supplied.
@param[in]      bodyB
    The MobilizedBody to which this osensor is fixed. OSensors on Ground
    are allowed but will be ignored.
@param[in]      orientationInB
    This is the orientation of the osensor in \a bodyB's local frame.
@param[in]      weight
    An optional weight for use in defining the objective function, which
    combines errors in this osensor's orientation with errors in other osensors'
    orientation. If the weight is zero this osensor is ignored.
@return The unique osensor index number assigned to this osensor. These are
assigned sequentially as the osensors are added. 
@note Adding an osensor invalidates any observation/osensor correspondence; be
sure to call defineObservationOrder() \e after defining all your osensors. **/
OSensorIx addOSensor(const String& name, MobilizedBodyIndex bodyB, 
                     const Rotation& orientationInB, Real weight=1)
{   SimTK_ERRCHK1_ALWAYS(isFinite(weight) && weight >= 0, 
        "OrientationSensors::addOSensor()", 
        "Illegal orientation sensor weight %g.", weight);
    uninitializeAssembler();
    // Forget any previously-established observation/osensor correspondence.
    observation2osensor.clear(); osensor2observation.clear(); 
    observations.clear();
    const OSensorIx ix(osensors.size());
    String nm = String::trimWhiteSpace(name);
    if (nm.empty())
        nm = String("_UNNAMED_") + String(ix);

    std::pair< std::map<String,OSensorIx>::iterator, bool >
        found = osensorsByName.insert(std::make_pair(nm,ix));
    SimTK_ERRCHK2_ALWAYS(found.second, // true if insertion was done
        "OSensors::addOSensor()",
        "OSensor name '%s' was already use for OSensor %d.",
        nm.c_str(), (int)found.first->second); 

    osensors.push_back(OSensor(nm,bodyB,orientationInB,weight));
    return ix; 
}

/** Define an unnamed osensor. A default name will be assigned; that name 
will be "_UNNAMED_XX" where XX is the OSensorIx assigned to that osensor 
(don't use names of that form yourself).  
@see addOSensor(name,...) for more information. **/
OSensorIx addOSensor(MobilizedBodyIndex bodyB, const Rotation& orientationInB,
                     Real weight=1)
{   return addOSensor("", bodyB, orientationInB, weight); }


/** Define the meaning of the observation data by giving the OSensorIx 
associated with each observation. The length of the array of osensor indices 
defines the expected number of observations to be provided for each observation
frame. Any osensor index that is supplied with an invalid value means that the
corresponding observation will be present in the supplied data but should be 
ignored.
@param[in]          observationOrder
    This is an array of osensor index values, one per observation, that defines
    both the number of expected observations and the osensor corresponding
    to each observation. OSensors can be in any order; an invalid osensor index
    means that observation will be provided but should be ignored; 
    osensors whose indices are never listed are ignored. If \a observationOrder
    is supplied as a zero-length array, then we'll assume there are as
    many observations as osensors and that their indices match.

@note If you don't call this method at all, a default correspondence will
be defined as described for a zero-length \a observationOrder array (that is,
same number of observations and osensors with matching indices). Whenever you 
add a new osensor, any previously defined observation order is forgotten so the 
default correspondence will be used unless you call this again. **/
void defineObservationOrder(const Array_<OSensorIx>& observationOrder) {
    uninitializeAssembler();
    if (observationOrder.empty()) {
        observation2osensor.resize(osensors.size());
        for (OSensorIx mx(0); mx < osensors.size(); ++mx)
            observation2osensor[ObservationIx(mx)] = mx;
    } else 
        observation2osensor = observationOrder;
    osensor2observation.clear(); 
    // We might need to grow this more, but this is an OK starting guess.
    osensor2observation.resize(observation2osensor.size()); // all invalid
    for (ObservationIx ox(0); ox < observation2osensor.size(); ++ox) {
        const OSensorIx mx = observation2osensor[ox];
        if (!mx.isValid()) continue;

        if (osensor2observation.size() <= mx)
            osensor2observation.resize(mx+1);
        SimTK_ERRCHK4_ALWAYS(!osensor2observation[mx].isValid(),
            "OSensors::defineObservationOrder()", 
            "An attempt was made to associate OSensor %d (%s) with" 
            " Observations %d and %d; only one Observation per OSensor"
            " is permitted.",
            (int)mx, getOSensorName(mx).c_str(), 
            (int)osensor2observation[mx], (int)ox);

        osensor2observation[mx] = ox;
    }
    // Make room for osensor observations.
    observations.clear();
    observations.resize(observation2osensor.size(),
                        Rotation().setRotationToNaN());
}

/** Define the meaning of the observations by giving the osensor name 
corresponding to each observation, as a SimTK::Array_<String>. The length of 
the array of osensor indices defines the expected number of observations. Any
osensor name that is unrecognized or empty means that the corresponding 
observation will be present in the supplied data but should be ignored. **/
void defineObservationOrder(const Array_<String>& observationOrder) 
{   Array_<OSensorIx> osensorIxs(observationOrder.size());
    for (ObservationIx ox(0); ox < observationOrder.size(); ++ox)
        osensorIxs[ox] = getOSensorIx(observationOrder[ox]);
    defineObservationOrder(osensorIxs); }

/** Define observation order using an std::vector of SimTK::String. **/
// no copy required
void defineObservationOrder(const std::vector<String>& observationOrder)
{   defineObservationOrder(ArrayViewConst_<String>(observationOrder)); }


/** Define observation order using an Array_ of std::string. **/
// must copy
void defineObservationOrder(const Array_<std::string>& observationOrder) 
{   const Array_<String> observations(observationOrder); // copy
    defineObservationOrder(observations); }

/** Define observation order using an std::vector of std::string. **/
// must copy
void defineObservationOrder(const std::vector<std::string>& observationOrder) 
{   const Array_<String> observations(observationOrder); // copy
    defineObservationOrder(observations); }

/** Define observation order using a C array of const char* names. **/
void defineObservationOrder(int n, const char* const observationOrder[]) 
{   Array_<OSensorIx> osensorIxs(n);
    for (ObservationIx ox(0); ox < n; ++ox)
        osensorIxs[ox] = getOSensorIx(String(observationOrder[ox]));
    defineObservationOrder(osensorIxs); }
/*@}*/



//------------------------------------------------------------------------------
/** @name               Retrieve setup information
These methods are used to query information associated with the construction
and setup of this OrientationSensors object. This information does not normally
change during an osensor-tracking study, although osensor weights may be changed
by some inverse kinematics methods. **/
/*@{*/

/** Return a count n of the number of currently-defined osensors. Valid
osensor index values (of type OrientationSensors::OSensorIx) are 0..n-1. **/
int getNumOSensors() const {return osensors.size();}

/** Return the unique osensor name assigned to the osensor whose index
is provided. If the osensor was defined without a name, this will return
the default name that was assigned to it. **/
const String& getOSensorName(OSensorIx ix) 
{   return osensors[ix].name; }

/** Return the osensor index associated with the given osensor name. If the
name is not recognized the returned index will be invalid (test with
index.isValid()). **/
const OSensorIx getOSensorIx(const String& name) 
{   std::map<String,OSensorIx>::const_iterator p = osensorsByName.find(name);
    return p == osensorsByName.end() ? OSensorIx() : p->second; }

/** Get the weight currently in use for the specified osensor; this can
be changed dynamically via changeOSensorWeight(). **/
Real getOSensorWeight(OSensorIx mx)
{   return osensors[mx].weight; }

/** Get the MobilizedBodyIndex of the body associated with this osensor. **/
MobilizedBodyIndex getOSensorBody(OSensorIx mx) const
{   return osensors[mx].bodyB; }

/** Get the orientation (coordinate axes fixed in its body frame) of the given
osensor. **/
const Rotation& getOSensorStation(OSensorIx mx) const
{   return osensors[mx].orientationInB; }

/** Return the number of observations that were defined via the last call to
defineObservationOrder(). These are not necessarily all being used. If 
defineObservationOrder() was never called, we'll expect the same number of
observations as osensors although that won't be set up until the Assembler has
been initialized. **/
int getNumObservations() const {return observation2osensor.size();}

/** Return the ObservationIx of the observation that is currently associated
with the given osensor, or an invalid index if the osensor doesn't have any
corresponding observation (in which case it is being ignored). An exception 
will be thrown if the given OSensorIx is not in the range 
0..getNumOSensors()-1. **/
ObservationIx getObservationIxForOSensor(OSensorIx mx) const 
{   return osensor2observation[mx]; }

/** Return true if the supplied osensor is currently associated with an 
observation. @see getObservationIxForOSensor() **/
bool hasObservation(OSensorIx mx) const 
{   return getObservationIxForOSensor(mx).isValid(); }

/** Return the OSensorIx of the osensor that is associated with the 
given observation, or an invalid index if the observation doesn't correspond
to any osensor (in which case it is being ignored). An exception will be
thrown if the given ObservationIx is not in the range 
0..getNumObservations()-1. **/
OSensorIx getOSensorIxForObservation(ObservationIx ox) const 
{   return observation2osensor[ox]; }

/** Return true if the supplied observation is currently associated with a 
osensor. @see getOSensorIxForObservation() **/
bool hasOSensor(ObservationIx ox) const 
{   return getOSensorIxForObservation(ox).isValid();}

/** The OSensors assembly condition organizes the osensors by body after
initialization; call this to get the list of osensors on any particular body.
If necessary the Assembler will be initialized. It is an error if this 
assembly condition has not yet been adopted by an Assembler. **/
const Array_<OSensorIx>& getOSensorsOnBody(MobilizedBodyIndex mbx) {
    static const Array_<OSensorIx> empty;
    SimTK_ERRCHK_ALWAYS(isInAssembler(), "OSensors::getOSensorsOnBody()",
        "This method can't be called until the OSensors object has been"
        " adopted by an Assembler.");
    initializeAssembler();
    PerBodyOSensors::const_iterator bodyp = bodiesWithOSensors.find(mbx);
    return bodyp == bodiesWithOSensors.end() ? empty : bodyp->second;
}
/*@}*/



//------------------------------------------------------------------------------
/** @name                Execution methods
These methods can be called between tracking steps to make step-to-step
changes without reinitialization, and to access the current values of
step-to-step data including the resulting osensor errors. **/
/*@{*/

/** Move a single osensor's observed orientation without moving any of the 
others. If the value contains a NaN, this osensor/observation pair will be 
ignored the next time the assembly goal cost function is calculated. **/
void moveOneObservation(ObservationIx ox, const Rotation& observation) {
    SimTK_ERRCHK_ALWAYS(!observations.empty(), "Assembler::moveOneObservation()",
        "There are currently no observations defined. Either the Assembler"
        " needs to be initialized to get the default observation order, or you"
        " should call defineObservationOrder() explicitly.");
    SimTK_ERRCHK2_ALWAYS(ox.isValid() && ox < observations.size(),
        "Assembler::moveOneObservation()", "ObservationIx %d is invalid or"
        " out of range; there are %d observations currently defined. Use"
        " defineObservationOrder() to specify the set of observations and how"
        " they correspond to osensors.", 
        (int)ox, (int)observations.size()); 
    observations[ox] = observation; 
}

/** Set the observed osensor orientations for a new observation frame. These are
the orientations to which we will next attempt to rotate all the corresponding 
osensors. Note that not all observations necessarily have corresponding osensors
defined; orientations of those osensors must still be provided here but they 
will be ignored. The length of the \a allObservations array must be the same as
the number of defined observations; you can obtain that using 
getNumObservations(). Any observations that contain a NaN will be ignored; that
osensor/observation pair will not be used in the next calculation of the 
assembly goal cost function. **/
void moveAllObservations(const Array_<Rotation>& observations) {
    SimTK_ERRCHK2_ALWAYS(   (int)observations.size() 
                         == (int)observation2osensor.size(),
        "OSensors::moveAllObservations()",
        "Number of observations provided (%d) differs from the number of"
        " observations (%d) last defined with defineObservationOrder().",
        observations.size(), observation2osensor.size());
    this->observations = observations;
}

/** Change the weight associated with a particular osensor. If this is just
a quantitative change (e.g., weight was 0.3 now it is 0.4) then this does
not require any reinitialization and will affect the goal calculation next
time it is done. If the weight changes to or from zero (a qualitative change)
then this will uninitialize the Assembler and all the internal data structures
will be changed to remove or add this osensor from the list of active osensors.
If you want to temporarily ignore an osensor without reinitializing, you can
set its corresponding observation to NaN in which case it will simply be
skipped when the goal value is calculated. **/
void changeOSensorWeight(OSensorIx mx, Real weight) {
   SimTK_ERRCHK1_ALWAYS(isFinite(weight) && weight >= 0, 
        "OSensors::changeOSensorWeight()", 
        "Illegal osensor weight %g.", weight);

    OSensor& osensor = osensors[mx];
    if (osensor.weight == weight)
        return;

    if (osensor.weight == 0 || weight == 0)
        uninitializeAssembler(); // qualitative change

    osensor.weight = weight;
}

/** Return the current value of the orientation for this observation. This
is the orientation to which we will try to rotate the corresponding osensor if 
there is one. The result might be NaN if there is no current value for this 
observation; you can check using Rotation's isFinite() method. **/
const Rotation& getObservation(ObservationIx ox) const 
{   return observations[ox]; }

/** Return the current values of all the observed orientations. These are the
orientations to which we will try to rotate the corresponding osensors, for 
those observations that have corresponding osensors defined. Some of the values
may be NaN if there is currently no corresponding observation. Note that these
are indexed by ObservationIx; use getObservationIxForOSensor() to map a 
OSensorIx to its corresponding ObservationIx. **/
const Array_<Rotation,ObservationIx>& getAllObservations() const
{   return observations; }

/** Using the current value of the internal state, calculate the ground
frame orientation of a particular osensor. The difference between this 
orientation and the corresponding observation is the current error for this 
osensor. **/
Rotation findCurrentOSensorOrientation(OSensorIx mx) const;

/** Using the current value of the internal state, calculate the error 
between the given osensor's current orientation and its corresponding observed
orientation (unweighted), as a nonnegative angle in radians between 0 and Pi. If
the osensor is not associated with an observation, 
or if the observed location is missing (indicated by a NaN value), then the 
error is reported as zero. **/
Real findCurrentOSensorError(OSensorIx mx) const {
    const ObservationIx ox = getObservationIxForOSensor(mx);
    if (!ox.isValid()) return 0; // no observation for this osensor
    const Rotation& R_GO = getObservation(ox);
    if (!R_GO.isFinite()) return 0; // NaN in observation; error is ignored
    const Rotation R_GS = findCurrentOSensorOrientation(mx);
    const Rotation R_SO = ~R_GS*R_GO; // orientation error, in S
    const Vec4 aa_SO = R_SO.convertRotationToAngleAxis();
    return std::abs(aa_SO[0]);
}
/*@}*/



//------------------------------------------------------------------------------
/** @name              AssemblyCondition virtuals
These methods are the implementations of the AssemblyCondition virtuals. **/
/*@{*/
int initializeCondition() const OVERRIDE_11;
void uninitializeCondition() const OVERRIDE_11;
int calcErrors(const State& state, Vector& err) const OVERRIDE_11;
int calcErrorJacobian(const State& state, Matrix& jacobian) const OVERRIDE_11;
int getNumErrors(const State& state) const OVERRIDE_11;
int calcGoal(const State& state, Real& goal) const OVERRIDE_11;
int calcGoalGradient(const State& state, Vector& grad) const OVERRIDE_11;
/*@}*/

//------------------------------------------------------------------------------
                                    private:
//------------------------------------------------------------------------------
const OSensor& getOSensor(OSensorIx i) const {return osensors[i];}
OSensor& updOSensor(OSensorIx i) {uninitializeAssembler(); return osensors[i];}

                                // data members                               
                               
// OSensor definition. Any change here except a quantitative change to the
// osensor's weight uninitializes the Assembler.
Array_<OSensor,OSensorIx>       osensors;
std::map<String,OSensorIx>      osensorsByName;

// Observation-osensor corresondence specification. Any change here 
// uninitializes the Assembler.
Array_<OSensorIx,ObservationIx> observation2osensor;

// For convience in mapping from an osensor to its corresponding observation.
// ObservationIx will be invalid if a particular osensor has no associated
// observation.
Array_<ObservationIx,OSensorIx> osensor2observation;

// This is the current set of osensor orientation observations, one per entry in 
// the observation2osensor array. Changing the values here does not uninitialize
// the Assembler.            
Array_<Rotation,ObservationIx>  observations;

// After initialize, this groups the osensors by body and weeds out
// any zero-weighted osensors. TODO: skip low-weighted osensors, at
// least at the start of the assembly.
typedef std::map<MobilizedBodyIndex,Array_<OSensorIx> > PerBodyOSensors;
mutable PerBodyOSensors         bodiesWithOSensors;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_ASSEMBLY_CONDITION_ORIENTATION_SENSORS_H_
