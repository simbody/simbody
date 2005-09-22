#ifndef SIMTK_FRAME_H_
#define SIMTK_FRAME_H_

class Vec3;
class Mat33;

typedef FixedVector<Vec3,2> SpatialVector;
typedef float_type Real;

/**
 * This class represents an orthogonal, right-handed coordinate frame F, 
 * measured from and expressed in a reference coordinate frame R. F consists of
 * 3 perpendicular unit vectors defining its axes as viewed from R, 
 * and a vector from R's origin point OR to F's origin point OF. Note that
 * the meaning of "R" comes from the context in which frame F is used.
 * We use the phrase "frame F is in frame R" to describe the above relationship,
 * that is, "in" means both measured in and expressed in. 
 *
 * The axis vectors are ordered 1-2-3 or x-y-z as you prefer, with
 * z = x X y, making a right-handed set. These axes are arranged as
 * columns of a 3x3 matrix Rot_RF = [ x y z ] which is a direction cosine
 * (rotation) matrix useful for conversions between frame R and F. For
 * example, given a vector vF expressed in the F frame, that same vector
 * re-expressed in R is given by vR = Rot_RF*vF. The origin point Loc_RF is 
 * stored as the vector Loc_RF=(OF-OR) and expressed in R.
 *
 * This is a "POD" (plain old data) class with a well-defined memory
 * layout on which a client of this class may depend: There are 
 * exactly 4 consecutive, packed 3-vectors in the order x,y,z,O.
 * That is, this class is equivalent to an array of 12 Reals with 
 * the order x1,x2,x3,y1,y2,y3,z1,z2,z3,O1,O2,O3. It is expressly allowed
 * to reinterpret Frame objects in any appropriate manner that depends
 * on this memory layout.
 */
class Frame {
public:
    Frame() { }
    Frame(const Mat33& axesInR, const Vec3& originInR) {setFrame(axesInR,originInR);}
    // default copy, assignment, destructor

    setFrame(const Mat33& axesInR, const Vec3& originInR) {Rot_RF=axesInR; Loc_RF=originInR;}

    // Transform various items measured and expressed in F to those same
    // items measured and expressed in F's reference frame R.
    Vec3  xformVector2Ref  (const Vec3& vF)      const { return Rot_RF*vF; }
    Vec3  xformStation2Ref (const Vec3& sF)      const { return Loc_RF + xformVector2Ref(sF); }
    Mat33 xformRotation2Ref(const Mat33& Rot_FX) const { return Rot_RF*Rot_FX; }
    Frame xformFrame2Ref(const Frame& fF) const 
      { return Frame(xformRotation2Ref(fF.Rot_RF), xformStation2Ref(fF.Loc_RF)); }

    const Mat33& getRot_RF() const { return Rot_RF; }
    Mat33&       updRot_RF()       { return Rot_RF; }

    const Vec3&  getLoc_RF() const { return Loc_RF; }
    Vec3&        updLoc_RF()       { return Loc_RF; }

    // Computation-free conversions
    const Real*  getFrameAsArray(const Frame& f) const {return reinterpret_cast<const Real*>(&f);}
    Real*        updFrameAsArray(Frame& f)             {return reinterpret_cast<Real*>(&f);}
    const Frame& getArrayAsFrame(const Real* r)  const {return *reinterpret_cast<const Frame*>(r);}
    Frame&       updArrayAsFrame(Real* r)              {return *reinterpret_cast<Frame*>(r);}

private:
    Mat33 Rot_RF;   // rotation matrix that expresses F's axes in R
    Vec3  Loc_RF;   // location of F's origin measured from R's origin, expressed in R 
};

class MassProperties {
public:
    MassProperties() { }
    MassProperties(const Real& m, const Vec3& com, const Mat33& inertia)
      { setMassProperties(m,com,inertia); }

    void setMassProperties(const Real& m, const Vec3& com, const Mat33& inertia)
      { mass=m; comInB=com; inertiaInB=inertia); }

    bool isMassless()   const { return mass==0.; }
    bool isCentroidal() const { return comInB.isZero(); }

private:
    Real  mass;
    Vec3  comInB;
    Mat33 inertiaInB;   // about B origin (symmetric matrix)
};

#endif /* SIMTK_FRAME_H_ */
