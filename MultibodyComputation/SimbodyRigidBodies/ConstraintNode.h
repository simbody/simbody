#ifndef SIMTK_SIMBODY_CONSTRAINT_NODE_H_
#define SIMTK_SIMBODY_CONSTRAINT_NODE_H_

#include "simbody/Simbody.h"
#include "SimbodyTree.h"
#include "SimbodyTreeState.h"

#include <cassert>
#include <vector>

using namespace simtk;

/**
 * This class represents a "constraint", which is in general a set of related
 * constraint equations. 
 */
class ConstraintNode {
public:
    virtual ~ConstraintNode() {}

    ConstraintNode& operator=(const ConstraintNode&);

    virtual const char* type()     const {return "unknown";}

    void setConstraintNum(int n) {constraintNum=n;}

        // TOPOLOGICAL INFO: no State needed

    int getConstraintNum()   const {return constraintNum;}
    int getMultIndex() const {return multIndex;}
    virtual int getMaxNMult() const=0; // # allocated lambda slots

        // MODELING INFO
    virtual int getNMult(const SBStateRep&) const=0; // # lambda slots in use

protected:
    // This is the constructor for the abstract base type for use by the derived
    // concrete types in their constructors.
    ConstraintNode() : multIndex(-1), constraintNum(-1) { }

    typedef std::vector<RigidBodyNode*>   RigidBodyNodeList;

    int               multIndex;      // index into lambda array
    int               constraintNum;  // unique ID number in RigidBodyTree

    friend std::ostream& operator<<(std::ostream& s, const ConstraintNode&);
};

#endif // SIMTK_SIMBODY_CONSTRAINT_NODE_H_
