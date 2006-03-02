/**@file
 * This file contains the code used to build the various constraint types.
 */

#include "RigidBodyTree.h"
#include "ConstraintNode.h"

#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::setprecision;

///////////////////////////////////////////////
// Implementation of ConstraintNode methods. //
///////////////////////////////////////////////



/////////////////////////////////////////////////
// Define classes derived from ConstraintNode. //
/////////////////////////////////////////////////

/**
 * This class represents a single constraint equation, maintaining a fixed distance
 * between stations on two different bodies (or a body and ground).
 */
class ConstantDistanceConstraint : public ConstraintNode {
public:
    ConstantDistanceConstraint()
      : ConstraintNode() {}
    ~ConstantDistanceConstraint() {}

    /*virtual*/ int getMaxNMult() const {return 1;}
    /*virtual*/ int getNMult(const SBStateRep&) const {return 1;}
};

/**
 * This class represents three constraint equations, together conspiring to hold
 * a station on one body coincident with one on another (or on ground).
 * The current implementation uses three distance constraints between the
 * pair of bodies.
 */
class CoincidentStationsConstraint : public ConstraintNode {
public:
    CoincidentStationsConstraint()
      : ConstraintNode() {}
    ~CoincidentStationsConstraint() {}

    /*virtual*/ int getMaxNMult() const {return 3;}
    /*virtual*/ int getNMult(const SBStateRep&) const {return 3;}
};

/**
 * This class represents six constraint equations, together serving to weld
 * one body to another (or to ground).
 */
class WeldConstraint : public ConstraintNode {
public:
    WeldConstraint()
      : ConstraintNode() {}
    ~WeldConstraint() {}

    /*virtual*/ int getMaxNMult() const {return 6;}
    /*virtual*/ int getNMult(const SBStateRep&) const {return 6;}
};

