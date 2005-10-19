/**@file
 * Implementations of high level multibody modeling objects for Simbody.
 */

#include "SimTK_Frame.h"

#include <string>
#include <vector>

using namespace simtk;


// Declare the 'rep' classes which actually represent the objects pointed
// to by the wrapper classes:
//    Frame,Station,Direction,MassElement,Body,Joint,RigidBody,DeformableBody,
//    Multibody,MultibodySystem

// TODO should be parameterizable.
class StationPlacement {
public:
    StationPlacement() : parentFrame(0) { }
    virtual ~StationPlacement() { }

    bool hasPlacement() const { return parentFrame != 0; }
    const Frame& getParentFrame() const {assert(parentFrame); return *parentFrame;}

    // These should allow for state to be passed in.
    virtual Vec3  getMeasureNumbers()     const = 0;

private:
    Frame*  parentFrame;    // reference to an existing Frame, don't delete!
};

class FramePlacement {
public:
    FramePlacement() : parentFrame(0) { }
    virtual ~FramePlacement() { }

    bool hasPlacement() const { return parentFrame != 0; }
    const Frame& getParentFrame() const {assert(parentFrame); return *parentFrame;}

    // These should allow for state to be passed in.
    virtual Vec3  getStation()     const = 0;   // origin measure numbers
    virtual Mat33 getOrientation() const = 0;   // columns are x,y,z measure numbers

private:
    Frame*  parentFrame;    // reference to an existing Frame, don't delete!
};

// A concrete Placement in which there are no variables.
class FixedFramePlacement : public FramePlacement {
public:
    FixedPlacement(const Vec3& sta, const Mat33& ori)
      : station(sta), orientation(ori)
    { }

    // Implementations of pure virtuals.
    Vec3  getStation()     const { return station; }
    Mat33 getOrientation() const { return orientation; }

private:
    Vec3    station;
    Mat33   orientation;
};

class FrameRep {
public:
    FrameRep(const char* nm) {
        name = nm;
        childStations.push_back  (Station("O"));
        childDirections.push_back(Direction("x"));
        childDirections.push_back(Direction("y"));
        childDirections.push_back(Direction("z"));

        childStations[0].setPlacement  (*this, Vec3(0,0,0));
        childDirections[0].setPlacement(*this, Vec3(1,0,0));
        childDirections[1].setPlacement(*this, Vec3(0,1,0));
        childDirections[2].setPlacement(*this, Vec3(0,0,1));
    }

    const string getFullName() const { 
        const string s;
        if (place.hasPlacement())
            s = place.getParent().getName() + "/";
        return s + getLocalName(); 
    }

    const string getLocalName() const { return name; }

    const Station& getOrigin() const { return childStations[0]; }
    const Direction& getAxis(int i) const 
      { assert(0<=i && i<=2); return childDirections[i]; }
    const Direction& x() const { return childDirections[0]; }
    const Direction& y() const { return childDirections[1]; }
    const Direction& z() const { return childDirections[2]; }

    const Station&   addStation(const StationPlacement&, const char* nm);
    const Direction& addDirection(const DirectionPlacement&, const char* nm);
    const Frame&     addFrame(const FramePlacement&, const char* nm);

private:
    string       name;

    // Objects wholly owned by this Frame and placed on it. The first
    // station is this Frame's origin, and the first three Directions
    // are this Frame's x,y,z axes, respectively. The rest are in no
    // particular order.
    vector<Station>   childStations;
    vector<Direction> childDirections;
    vector<Frame>     childFrames;

    // If this Frame has been placed with respect to some other
    // Frame, this is the placement information.
    FramePlacement place;
};

class StationRep {
public:
private:
    string           name;
    StationPlacement station;
};

class DirectionRep {
public:
private:
    string name;
    DirectionPlacement direction;
};

class MassElementRep {
public:
private:
    string name;
};

class JointRep {
public:
private:
    string name;
};

class MultibodySystemRep {
public:
private:
    string name;
};

// Body extends Frame.
class BodyRep {
public:
    BodyRep();
    void addMassElementLike(const MassElement&, const Placement&);
private:
    vector<MassElement> childMassElements;
};

// RigidBody extends Body.
class RigidBodyRep {
public:
private:
};

// DeformableBody extends Body.
class DeformableBodyRep {
public:
private:
    string name;
};

// Multibody extends Body.
class MultibodyRep {
public:
private:
    string name;
};

