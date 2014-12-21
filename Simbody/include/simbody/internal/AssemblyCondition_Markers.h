#ifndef SimTK_SIMBODY_ASSEMBLY_CONDITION_MARKERS_H_
#define SimTK_SIMBODY_ASSEMBLY_CONDITION_MARKERS_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-14 Stanford University and the Authors.        *
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
//                                  MARKERS
//------------------------------------------------------------------------------
/** This AssemblyCondition specifies a correspondence between stations on
mobilized bodies ("markers") and fixed ground-frame locations ("observations").
The idea is to adjust the q's so that each marker is located close to its
corresponding observation. This is normally used as a goal since we don't
expect a perfect fit, but you can use these as a set of assembly error 
conditions if there are enough degrees of freedom to achieve a near-perfect 
solution. 

Markers are defined one at a time and assigned sequential marker index values
of type Markers::MarkerIx. They may optionally be given unique, case-sensitive
names, and we will keep a map from name to MarkerIx. A default name will be 
assigned if none is given. A weight is assigned to every marker, with default 
weight=1. We do not expect that all the markers will be used; markers with 
weights of zero will not be included in the study, nor will markers for which 
no observation is given.

Once specified, the marker definitions do not change during a series of inverse
kinematic (tracking) steps. The observations, on the other hand, are expected 
to come from a time series of experimental measurements of marker locations and
will be different at every step. They typically come from a file organized by 
"frame", meaning an observation time and a set of observed locations, one per 
marker, corresponding to that time. During initial setup, the number of 
observations per frame and their correspondence to the defined markers is 
specified. They can be in any order, may skip some markers, and may include 
data for markers that are not defined. However, once initialized each frame 
must supply the same information in the same order. Data for an unobserved 
marker can be provided as NaN in which case it will be ignored in that frame. 
The frame time is supplied to the track() method which initiates assembly for 
a frame.

Observation-marker correspondence maps a ObservationIx to a unique MarkerIx. 
By default, we'll expect to get an observation for each marker and that the 
observation order and the marker order are the same, i.e. 
ObservationIx==MarkerIx for every marker. However, you can instead define 
observation/marker correspondence yourself, (\e after all markers have been 
defined), via one of the defineObservationOrder() methods. This is done by 
supplying an array of MarkerIx values, or an array of Marker names, with the 
array elements ordered by ObservationIx. Any invalid marker index or 
unrecognized marker name means we will ignore values provide for that 
observation; similarly, any markers whose index or name is not specified at 
all will be ignored. 
**/
class SimTK_SIMBODY_EXPORT Markers : public AssemblyCondition {

// This is a private class used in the implementation below but not
// accessible through the API.
struct Marker {
    Marker(const String& name, MobilizedBodyIndex bodyB, 
           const Vec3& markerInB, Real weight = 1)
    :   name(name), bodyB(bodyB), markerInB(markerInB), weight(weight) 
    { assert(weight >= 0); }

    Marker(MobilizedBodyIndex bodyB, const Vec3& markerInB, Real weight=1)
    :   name(""), bodyB(bodyB), markerInB(markerInB), weight(weight) 
    { assert(weight >= 0); }

    String              name;
    MobilizedBodyIndex  bodyB;
    Vec3                markerInB;
    Real                weight; 
};

public:

/** Define the MarkerIx type which is just a uniquely-typed int. **/
SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Markers,MarkerIx);
/** Define the ObservationIx type which is just a uniquely-typed int. **/
SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Markers,ObservationIx);



//------------------------------------------------------------------------------
/** @name                Construction and setup
These methods are used as an extended construction phase for Markers
objects, defining the markers and observations that will be used in the
subsequent tracking steps. **/
/*@{*/

/** The default constructor creates an empty Markers AssemblyCondition
object that should be filled in with calls to addMarker() and optionally
defineObservationOrder(). **/
Markers() : AssemblyCondition("Markers") {}

/** Define a new marker attached to a particular MobilizedBody. Note that
a marker will be ignored unless an observation is provided for it.
@param[in]      name
    A unique name to be used to identify this marker. If the name is
    empty or blank, a default name will be supplied.
@param[in]      bodyB
    The MobilizedBody to which this marker is fixed. Markers on Ground
    are allowed but will be ignored.
@param[in]      markerInB
    This is the position vector of the marker in \a bodyB's local frame,
    also known as the marker's "station" on \a bodyB.
@param[in]      weight
    An optional weight for use in defining the objective function, which
    combines errors in this marker's position with errors in other markers'
    positions. If the weight is zero this marker is ignored.
@return The unique marker index number assigned to this marker. These are
assigned sequentially as the marker are added. 
@note Adding a marker invalidates any observation/marker correspondence; be
sure to call defineObservationOrder() \e after defining all your markers. **/
MarkerIx addMarker(const String& name, MobilizedBodyIndex bodyB, 
                   const Vec3& markerInB, Real weight=1)
{   SimTK_ERRCHK1_ALWAYS(isFinite(weight) && weight >= 0, 
        "Markers::addMarker()", "Illegal marker weight %g.", weight);
    uninitializeAssembler();
    // Forget any previously-established observation/marker correspondence.
    observation2marker.clear(); marker2observation.clear(); 
    observations.clear();
    const MarkerIx ix(markers.size());
    String nm = String::trimWhiteSpace(name);
    if (nm.empty())
        nm = String("_UNNAMED_") + String(ix);

    std::pair< std::map<String,MarkerIx>::iterator, bool >
        found = markersByName.insert(std::make_pair(nm,ix));
    SimTK_ERRCHK2_ALWAYS(found.second, // true if insertion was done
        "Markers::addMarker()",
        "Marker name '%s' was already use for Marker %d.",
        nm.c_str(), (int)found.first->second); 

    markers.push_back(Marker(nm,bodyB,markerInB,weight));
    return ix; 
}

/** Define an unnamed marker. A default name will be assigned; that name 
will be "_UNNAMED_XX" where XX is the MarkerIx assigned to that marker 
(don't use names of that form yourself).  
@see addMarker(name,...) for more information. **/
MarkerIx addMarker(MobilizedBodyIndex bodyB, const Vec3& markerInB,
                   Real weight=1)
{   return addMarker("", bodyB, markerInB, weight); }


/** Define the meaning of the observation data by giving the MarkerIx 
associated with each observation. The length of the array of marker indices 
defines the expected number of observations to be provided for each observation
frame. Any marker index that is supplied with an invalid value means that the
corresponding observation will be present in the supplied data but should be 
ignored.
@param[in]          observationOrder
    This is an array of marker index values, one per observation, that defines
    both the number of expected observations and the marker corresponding
    to each observation. Markers can be in any order; an invalid marker index
    means that observation will be provided but should be ignored; 
    markers whose indices are never listed are ignored. If \a observationOrder
    is supplied as a zero-length array, then we'll assume there are as
    many observations as markers and that their indices match.

@note If you don't call this method at all, a default correspondence will
be defined as described for a zero-length \a observationOrder array (that is,
same number of observations and markers with matching indices). Whenever you 
add a new marker, any previously defined observation order is forgotten so the 
default correspondence will be used unless you call this again. **/
void defineObservationOrder(const Array_<MarkerIx>& observationOrder) {
    uninitializeAssembler();
    if (observationOrder.empty()) {
        observation2marker.resize(markers.size());
        for (MarkerIx mx(0); mx < markers.size(); ++mx)
            observation2marker[ObservationIx(mx)] = mx;
    } else 
        observation2marker = observationOrder;
    marker2observation.clear(); 
    // We might need to grow this more, but this is an OK starting guess.
    marker2observation.resize(observation2marker.size()); // all invalid
    for (ObservationIx ox(0); ox < observation2marker.size(); ++ox) {
        const MarkerIx mx = observation2marker[ox];
        if (!mx.isValid()) continue;

        if (marker2observation.size() <= mx)
            marker2observation.resize(mx+1);
        SimTK_ERRCHK4_ALWAYS(!marker2observation[mx].isValid(),
            "Markers::defineObservationOrder()", 
            "An attempt was made to associate Marker %d (%s) with" 
            " Observations %d and %d; only one Observation per Marker"
            " is permitted.",
            (int)mx, getMarkerName(mx).c_str(), 
            (int)marker2observation[mx], (int)ox);

        marker2observation[mx] = ox;
    }
    // Make room for marker observations.
    observations.clear();
    observations.resize(observation2marker.size(),Vec3(NaN));
}

/** Define the meaning of the observations by giving the marker name 
corresponding to each observation, as a SimTK::Array_<String>. The length of 
the array of marker indices defines the expected number of observations. Any
marker name that is unrecognized or empty means that the corresponding 
observation will be present in the supplied data but should be ignored. **/
void defineObservationOrder(const Array_<String>& observationOrder) 
{   Array_<MarkerIx> markerIxs(observationOrder.size());
    for (ObservationIx ox(0); ox < observationOrder.size(); ++ox)
        markerIxs[ox] = getMarkerIx(observationOrder[ox]);
    defineObservationOrder(markerIxs); }

/** Define observation order using an std::vector of SimTK::String. */
// no copy required
void defineObservationOrder(const std::vector<String>& observationOrder)
{   defineObservationOrder(ArrayViewConst_<String>(observationOrder)); }


/** Define observation order using an Array_ of std::string. */
// must copy
void defineObservationOrder(const Array_<std::string>& observationOrder) 
{   const Array_<String> observations(observationOrder); // copy
    defineObservationOrder(observations); }

/** Define observation order using an std::vector of std::string. */
// must copy
void defineObservationOrder(const std::vector<std::string>& observationOrder) 
{   const Array_<String> observations(observationOrder); // copy
    defineObservationOrder(observations); }

/** Define observation order using a C array of const char* names. */
void defineObservationOrder(int n, const char* const observationOrder[]) 
{   Array_<MarkerIx> markerIxs(n);
    for (ObservationIx ox(0); ox < n; ++ox)
        markerIxs[ox] = getMarkerIx(String(observationOrder[ox]));
    defineObservationOrder(markerIxs); }
/*@}*/



//------------------------------------------------------------------------------
/** @name               Retrieve setup information
These methods are used to query information associated with the construction
and setup of this Markers object. This information does not normally change
during a marker-tracking study, although marker weights may be changed by 
some inverse kinematics methods. **/
/*@{*/

/** Return a count n of the number of currently-defined markers. Valid
marker index values (of type Markers::MarkerIx) are 0..n-1. **/
int getNumMarkers() const {return markers.size();}

/** Return the unique marker name assigned to the marker whose index
is provided. If the marker was defined without a name, this will return
the default name that was assigned to it. **/
const String& getMarkerName(MarkerIx ix) 
{   return markers[ix].name; }
/** Return the marker index associated with the given marker name. If the
name is not recognized the returned index will be invalid (test with
index.isValid()). **/
const MarkerIx getMarkerIx(const String& name) 
{   std::map<String,MarkerIx>::const_iterator p = markersByName.find(name);
    return p == markersByName.end() ? MarkerIx() : p->second; }

/** Get the weight currently in use for the specified marker; this can
be changed dynamically via changeMarkerWeight(). **/
Real getMarkerWeight(MarkerIx mx)
{   return markers[mx].weight; }

/** Get the MobilizedBodyIndex of the body associated with this marker. **/
MobilizedBodyIndex getMarkerBody(MarkerIx mx) const
{   return markers[mx].bodyB; }

/** Get the station (fixed location in its body frame) of the given marker. **/
const Vec3& getMarkerStation(MarkerIx mx) const
{   return markers[mx].markerInB; }

/** Return the number of observations that were defined via the last call to
defineObservationOrder(). These are not necessarily all being used. If 
defineObservationOrder() was never called, we'll expect the same number of
observations as markers although that won't be set up until the Assembler has
been initialized. **/
int getNumObservations() const {return observation2marker.size();}

/** Return the ObservationIx of the observation that is currently associated
with the given marker, or an invalid index if the marker doesn't have any
corresponding observation (in which case it is being ignored). An exception 
will be thrown if the given MarkerIx is not in the range 
0..getNumMarkers()-1. **/
ObservationIx getObservationIxForMarker(MarkerIx mx) const 
{ return marker2observation[mx]; }

/** Return true if the supplied marker is currently associated with an 
observation. @see getObservationIxForMarker() **/
bool hasObservation(MarkerIx mx) const 
{ return getObservationIxForMarker(mx).isValid(); }

/** Return the MarkerIx of the marker that is associated with the 
given observation, or an invalid index if the observation doesn't correspond
to any marker (in which case it is being ignored). An exception will be
thrown if the given ObservationIx is not in the range 
0..getNumObservations()-1. **/
MarkerIx getMarkerIxForObservation(ObservationIx ox) const 
{ return observation2marker[ox]; }

/** Return true if the supplied observation is currently associated with a 
marker. @see getMarkerIxForObservation() **/
bool hasMarker(ObservationIx ox) const 
{ return getMarkerIxForObservation(ox).isValid();}

/** The Markers assembly condition organizes the markers by body after
initialization; call this to get the list of markers on any particular body.
If necessary the Assembler will be initialized. It is an error if this 
assembly condition has not yet been adopted by an Assembler. **/
const Array_<MarkerIx>& getMarkersOnBody(MobilizedBodyIndex mbx) {
    static const Array_<MarkerIx> empty;
    SimTK_ERRCHK_ALWAYS(isInAssembler(), "Markers::getMarkersOnBody()",
        "This method can't be called until the Markers object has been"
        " adopted by an Assembler.");
    initializeAssembler();
    PerBodyMarkers::const_iterator bodyp = bodiesWithMarkers.find(mbx);
    return bodyp == bodiesWithMarkers.end() ? empty : bodyp->second;
}
/*@}*/



//------------------------------------------------------------------------------
/** @name                Execution methods
These methods can be called between tracking steps to make step-to-step
changes without reinitialization, and to access the current values of
step-to-step data including the resulting marker errors. **/
/*@{*/

/** Move a single marker's observed location without moving any of the others.
If the value contains a NaN, this marker/observation pair will be ignored the
next time the assembly goal cost function is calculated. **/
void moveOneObservation(ObservationIx ox, const Vec3& observation) 
{   SimTK_ERRCHK_ALWAYS(!observations.empty(), "Assembler::moveOneObservation()",
        "There are currently no observations defined. Either the Assembler"
        " needs to be initialized to get the default observation order, or you"
        " should call defineObservationOrder() explicitly.");
    SimTK_ERRCHK2_ALWAYS(ox.isValid() && ox < observations.size(),
        "Assembler::moveOneObservation()", "ObservationIx %d is invalid or"
        " out of range; there are %d observations currently defined. Use"
        " defineObservationOrder() to specify the set of observations and how"
        " they correspond to markers.", 
        (int)ox, (int)observations.size()); 
    observations[ox] = observation; 
}

/** Set the observed marker locations for a new observation frame. These are
the locations to which we will next attempt to move all the corresponding 
markers. Note that not all observations necessarily have corresponding markers
defined; locations of those markers must still be provided here but they will 
be ignored. The length of the \a allObservations array must be the same as the 
number of defined observations; you can obtain that using getNumObservations().
Any observations that contain a NaN will be ignored; that marker/observation 
pair will not be used in the next calculation of the assembly goal cost 
function. **/
void moveAllObservations(const Array_<Vec3>& observations) 
{   SimTK_ERRCHK2_ALWAYS((int)observations.size() == (int)observation2marker.size(),
        "Markers::moveAllObservations()",
        "Number of observations provided (%d) differs from the number of"
        " observations (%d) last defined with defineObservationOrder().",
        observations.size(), observation2marker.size());
    this->observations = observations; }

/** Change the weight associated with a particular marker. If this is just
a quantitative change (e.g., weight was 0.3 now it is 0.4) then this does
not require any reinitialization and will affect the goal calculation next
time it is done. If the weight changes to or from zero (a qualitative change)
then this will uninitialize the Assembler and all the internal data structures
will be changed to remove or add this marker from the list of active markers.
If you want to temporarily ignore a marker without reinitializing, you can
set its corresponding observation to NaN in which case it will simply be
skipped when the goal value is calculated. **/
void changeMarkerWeight(MarkerIx mx, Real weight) {
   SimTK_ERRCHK1_ALWAYS(isFinite(weight) && weight >= 0, 
        "Markers::changeMarkerWeight()", "Illegal marker weight %g.", weight);

    Marker& marker = markers[mx];
    if (marker.weight == weight)
        return;

    if (marker.weight == 0 || weight == 0)
        uninitializeAssembler(); // qualitative change

    marker.weight = weight;
}

/** Return the current value of the location for this observation. This
is where we will try to move the corresponding marker if there is one. 
The result might be NaN if there is no current value for this observation;
you can check using Vec3's isFinite() method. **/
const Vec3& getObservation(ObservationIx ox) const {return observations[ox];}
/** Return the current values of all the observed locations. This is where we 
will try to move the corresponding markers, for those observations that have 
corresponding markers defined. Some of the values may be NaN if there is
currently no corresponding observation. Note that these are indexed by
ObservationIx; use getObservationIxForMarker() to map a MarkerIx to its
corresponding ObservationIx. **/
const Array_<Vec3,ObservationIx>& getAllObservations() const
{   return observations; }

/** Using the current value of the internal state, calculate the ground
frame location of a particular marker. The difference between this location
and the corresponding observation is the current error for this marker. **/
Vec3 findCurrentMarkerLocation(MarkerIx mx) const;

/** Using the current value of the internal state, calculate the distance 
between the given marker's current location and its corresponding observed
location (unweighted). If the marker is not associated with an observation, 
or if the observed location is missing (indicated by a NaN value), then the 
error is reported as zero. 
@note If you actually want the square of the distance, you can save some
computation time by using findCurrentMarkerErrorSquared() which avoids the
square root needed to find the actual distance.
@see findCurrentMarkerErrorSquared() **/
Real findCurrentMarkerError(MarkerIx mx) const
{   return std::sqrt(findCurrentMarkerErrorSquared(mx)); }

/** Using the current value of the internal state, calculate the (unweighted)
square of the distance between the given marker's current location and its 
corresponding observed location (the squared distance is less expensive to 
compute than the distance). If the marker is not associated with an 
observation, or if the observed location is missing (indicated by a NaN 
value), then the error is reported as zero. 
@see findCurrentMarkerError() **/
Real findCurrentMarkerErrorSquared(MarkerIx mx) const {
    const ObservationIx ox = getObservationIxForMarker(mx);
    if (!ox.isValid()) return 0; // no observation for this marker
    const Vec3& loc = getObservation(ox);
    if (!loc.isFinite()) return 0; // NaN in observation; error is ignored
    return (findCurrentMarkerLocation(mx) - loc).normSqr();
}
/*@}*/



//------------------------------------------------------------------------------
/** @name              AssemblyCondition virtuals
These methods are the implementations of the AssemblyCondition virtuals. **/
/*@{*/
int initializeCondition() const override;
void uninitializeCondition() const override;
int calcErrors(const State& state, Vector& err) const override;
int calcErrorJacobian(const State& state, Matrix& jacobian) const override;
int getNumErrors(const State& state) const override;
int calcGoal(const State& state, Real& goal) const override;
int calcGoalGradient(const State& state, Vector& grad) const override;
/*@}*/

//------------------------------------------------------------------------------
                                    private:
//------------------------------------------------------------------------------
const Marker& getMarker(MarkerIx i) const {return markers[i];}
Marker& updMarker(MarkerIx i) {uninitializeAssembler(); return markers[i];}

                                // data members                               
                               
// Marker definition. Any change here except a quantitative change to the
// marker's weight uninitializes the Assembler.
Array_<Marker,MarkerIx>         markers;
std::map<String,MarkerIx>       markersByName;

// Observation-marker corresondence specification. Any change here 
// uninitializes the Assembler.
Array_<MarkerIx,ObservationIx>  observation2marker;

// For convience in mapping from a marker to its corresponding observation.
// ObservationIx will be invalid if a particular marker has no associated
// observation.
Array_<ObservationIx,MarkerIx>  marker2observation;

// This is the current set of marker location observations, one per entry in 
// the observation2marker array. Changing the values here does not uninitialize
// the Assembler.            
Array_<Vec3,ObservationIx>      observations;

// After initialize, this groups the markers by body and weeds out
// any zero-weighted markers. TODO: skip low-weighted markers, at
// least at the start of the assembly.
typedef std::map<MobilizedBodyIndex,Array_<MarkerIx> > PerBodyMarkers;
mutable PerBodyMarkers          bodiesWithMarkers;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_ASSEMBLY_CONDITION_MARKERS_H_
