#ifndef SimTK_COORDINATEAXIS_H 
#define SimTK_COORDINATEAXIS_H 

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-10 Stanford University and the Authors.        *
 * Authors: Michael Sherman and Paul Mitiguy                                  *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/** @file
Defines the CoordinateAxis and CoordinateDirection classes. **/

#include <cassert>

namespace SimTK {

/** This class, along with its sister class CoordinateDirection, provides 
convenient manipulation of the three coordinate axes via the definition of 
three constants XAxis, YAxis, and ZAxis each with a unique subtype and implicit 
conversion to the integers 0, 1, and 2 whenever necessary.\ Methods are 
provided to allow code to be written once that can be used to work with the
axes in any order.

There are also three CoordinateDirection constants NegXAxis, NegYAxis, and
NegZAxis, also with unique types permitting efficient compile time 
manipulation. These do not correspond to integers, however. Instead, they are
objects containing one of the CoordinateAxis objects combined with an integer
that is 1 or -1 to indicate the direction along that axis. The unary negation
operator is overloaded so that -XAxis is NegXAxis and -NegZAxis is ZAxis.
There are implicit conversions to UnitVec3 for any CoordinateAxis or 
CoordinateDirection object, yielding the equivalent (normalized) unit vector
corresponding to any of the six directions, without doing any computation
(and in particular, without normalizing).
@see CoordinateDirection **/
class CoordinateAxis {
public:
    /** Explicit construction of a CoordinateAxis from a calculated integer
    that must be 0, 1, or 2 representing XAxis, YAxis, or ZAxis. **/
    explicit CoordinateAxis( int i ) : m_myAxisId(i) {assertIndexIsInRange(i);}

    /** Implicit conversion of a CoordinateAxis to int 0, 1, or 2. **/
    operator int() const {return m_myAxisId;}

    /** Return the "next" coordinate axis after this one:
        - XAxis.getNextAxis()  returns YAxis 
        - YAxis.getNextAxis()  returns ZAxis 
        - ZAxis.getNextAxis()  returns XAxis **/
    CoordinateAxis getNextAxis() const 
    {   return CoordinateAxis((m_myAxisId+1) % 3); }

    /** Return the "previous" coordinate axis before this one:
        - XAxis.getPreviousAxis()  returns ZAxis 
        - YAxis.getPreviousAxis()  returns XAxis
        - ZAxis.getPreviousAxis()  returns YAxis **/
    CoordinateAxis getPreviousAxis() const 
    {   return CoordinateAxis((m_myAxisId+2) % 3); }

    /** Given this coordinate axis and one other, return the missing one:
        - XAxis.getThirdAxis(YAxis) returns ZAxis (and vice versa)
        - XAxis.getThirdAxis(ZAxis) returns YAxis (and vice versa)
        - YAxis.getThirdAxis(ZAxis) returns XAxis (and vice versa)
    @param[in] axis2    A coordinate axis that must be distinct from the
        current one; it is a fatal error to provide the same axis.
    @return The unmentioned third axis. **/
    CoordinateAxis getThirdAxis( const CoordinateAxis& axis2 ) const { 
       assert( isDifferentAxis(axis2) );
       CoordinateAxis nextAxis = getNextAxis();
       return nextAxis.isDifferentAxis(axis2) ? nextAxis : axis2.getNextAxis();
    }

    /** Return true if this is the X axis. **/
    bool isXAxis() const {return m_myAxisId == 0;}
    /** Return true if this is the Y axis. **/
    bool isYAxis() const {return m_myAxisId == 1;}
    /** Return true if this is the Z axis. **/
    bool isZAxis() const {return m_myAxisId == 2;}
    /** Return true if the given \a axis2 is the one following this one as
    would be reported by getNextAxis(). **/
    bool isNextAxis( const CoordinateAxis& axis2 ) const
    {   return int(getNextAxis()) == int(axis2); }
    /** Return true if the given \a axis2 is the one preceding this one as
    would be reported by getPreviousAxis(). **/
    bool isPreviousAxis( const CoordinateAxis& axis2 ) const
    {   return int(getPreviousAxis()) == int(axis2); }
    /** Return true if the given \a axis2 is the same as this one.\ You
    can use operator==() to perform the same comparison. **/
    bool isSameAxis( const CoordinateAxis& axis2 ) const
    {   return m_myAxisId == int(axis2); }
    /** Return true if both \a axis2 and \a axis3 are the same as this one. **/
    bool areAllSameAxes( const CoordinateAxis& axis2, const CoordinateAxis &axis3 ) const       
    {   return isSameAxis(axis2) && isSameAxis(axis3); }
    /** Return true if the given \a axis2 is not the same one as this 
    one.\ You can use operator!=() to perform the same comparison.  **/
    bool isDifferentAxis( const CoordinateAxis& axis2 ) const                                   
    {   return m_myAxisId != int(axis2); }
    /** Return true if neither \a axis2 nor \a axis3 is the same as this
    axis nor each other; that is, (this,axis2,axis3) together cover all three
    axes. **/
    bool areAllDifferentAxes( const CoordinateAxis& axis2, const CoordinateAxis& axis3 ) const  
    {   return isDifferentAxis(axis2) && isDifferentAxis(axis3) && axis2.isDifferentAxis(axis3); }
    /** Return true if the given \a axis2 is the one following this one in a
    forward cyclical direction, that is, if \a axis2 is the one that would be
    reported by getNextAxis(). **/
    bool isForwardCyclical( const CoordinateAxis& axis2 ) const                                
    {   return isNextAxis(axis2); }
    /** Return true if the given \a axis2 is the one following this one in a
    reverse cyclical direction, that is, if \a axis2 is the one that would be
    reported by getPreviousAxis(). **/
    bool isReverseCyclical( const CoordinateAxis& axis2 ) const                                 
    {   return isPreviousAxis(axis2); }

    /** Perform a specialized dot product between this axis and \a axis2; 
    returning one if they are the same axis and zero otherwise, without
    performing any floating point operations. **/
    int dotProduct(  const CoordinateAxis& axis2 ) const
    {   return isSameAxis(axis2) ? 1 : 0; }
    /** Return the sign that would result from a cross product between this 
    axis and \a axis2: zero if \a axis2 is the same as this axis; one if the 
    result would be in the positive direction along the third axis; -1 if it 
    would be in the negative direction. No floating point computations are
    performed. @see crossProductAxis() **/
    int crossProductSign( const CoordinateAxis& axis2 ) const
    {   return isSameAxis(axis2) ? 0 : (isNextAxis(axis2) ? 1 : -1); }
    /** Return the coordinate axis along which the cross product of this axis 
    and \a axis2 would lie: same as this if \a axis2 is the same as this axis
    (doesn't matter because the sign would be zero); otherwise, the third
    axis that is neither this one nor \a axis2. But note that the actual
    result may be along that axis or in the negative direction along that
    axis.  No floating point computations are performed. 
    @see crossProductSign(). **/
    CoordinateAxis crossProductAxis( const CoordinateAxis& axis2 ) const
    {   return isSameAxis(axis2) ? CoordinateAxis(m_myAxisId) 
                                 : getThirdAxis(axis2); }
    /** Return the axis and sign along that axis that would result from a
    cross product between this axis and \a axis2; this combines the functions
    of both crossProductAxis() and crossProductSign(). Note that if \a axis2 is the same as this
    axis we'll just return this as the axis but the sign is zero since the
    magnitude of the result would be zero. No floating point calculations are 
    performed. @see crossProductSign(), crossProductAxis() **/
    CoordinateAxis crossProduct( const CoordinateAxis& axis2, int& sign ) const
    {   sign = crossProductSign(axis2); return crossProductAxis(axis2); }

    /** Return a reference to the CoordinateAxis constant XAxis, YAxis, or
    ZAxis corresponding to the given integer index which must be 0, 1, or 2. **/
    static const CoordinateAxis& getCoordinateAxis( int i );

    /** Return true if the given integer is suitable as a coordinate axis, 
    meaning it is one of 0, 1, or 2 designating XAxis, YAxis, or ZAxis, 
    respectively. **/
    static bool  isIndexInRange( int i )        { return 0<=i && i<=2; }
    /** When in Debug mode, throw an assertion if the given integer is not
    suited as a coordinate axis, as defined by isIndexInRange(). **/
    static void  assertIndexIsInRange( int i )  { assert( isIndexInRange(i) ); } 

    // Forward declarations for subsequent helper classes
    class XCoordinateAxis; class YCoordinateAxis; class ZCoordinateAxis;
protected:
    /** @cond **/ // turn off doxygen here; these aren't for users
    class XTypeAxis{};
    class YTypeAxis{};
    class ZTypeAxis{};

    CoordinateAxis( const XTypeAxis& ) : m_myAxisId(0) {}
    CoordinateAxis( const YTypeAxis& ) : m_myAxisId(1) {}
    CoordinateAxis( const ZTypeAxis& ) : m_myAxisId(2) {}
    /** @endcond **/
private:            

    int m_myAxisId;
};


// Helper classes that allow compile time recognition of axis directions.
class CoordinateAxis::XCoordinateAxis : public CoordinateAxis {
  public: XCoordinateAxis() : CoordinateAxis(XTypeAxis()) {}
};
class CoordinateAxis::YCoordinateAxis : public CoordinateAxis {
  public: YCoordinateAxis() : CoordinateAxis(YTypeAxis()) {}
};
class CoordinateAxis::ZCoordinateAxis : public CoordinateAxis {
  public: ZCoordinateAxis() : CoordinateAxis(ZTypeAxis()) {}
};

/** Constant representing the X coordinate axis; will implicitly convert to
the integer 0 when used in a context requiring an integer. **/
static const CoordinateAxis::XCoordinateAxis  XAxis;
/** Constant representing the Y coordinate axis; will implicitly convert to
the integer 1 when used in a context requiring an integer. **/
static const CoordinateAxis::YCoordinateAxis  YAxis;
/** Constant representing the Z coordinate axis; will implicitly convert to
the integer 2 when used in a context requiring an integer. **/
static const CoordinateAxis::ZCoordinateAxis  ZAxis;

inline const CoordinateAxis& CoordinateAxis::getCoordinateAxis(int i) {
    assertIndexIsInRange(i);
    return (i==0 ? static_cast<const CoordinateAxis&>(XAxis) 
         : (i==1 ? static_cast<const CoordinateAxis&>(YAxis) 
                 : static_cast<const CoordinateAxis&>(ZAxis)));
}

/// Compare two CoordinateAxis objects. @relates CoordinateAxis 
inline bool operator==(const CoordinateAxis& a1, const CoordinateAxis& a2)
{   return a1.isSameAxis(a2); }

/// Compare two CoordinateAxis objects. @relates CoordinateAxis 
inline bool operator!=(const CoordinateAxis& a1, const CoordinateAxis& a2)
{   return a1.isDifferentAxis(a2); }


/** A CoordinateDirection is a CoordinateAxis plus a direction indicating the 
positive or negative direction along that axis. There are only six possible 
values for a CoordinateDirection, and there are predefined constants available 
covering all of them:
  - XAxis, YAxis, ZAxis are the CoordinateAxis types; they will implicitly 
    convert to positive axis directions.
  - NegXAxis, NegYAxis, NegZAxis are the negative directions.
  - The unary negation operator is overloaded so that -XAxis produces
    NegXAxis and -NegYAxis produces YAxis.
You can also produce CoordinateDirections at compile time or run time from
calculated axes and directions. 
@see CoordinateAxis **/
class CoordinateDirection {
public:
    /** Use for compile-time construction of a negative CoordinateDirection
    along one of the coordinate axes. **/
    class Negative {};

    /** Implicit conversion of a CoordinateAxis to a positive 
    CoordinateDirection along that axis. **/
    CoordinateDirection(const CoordinateAxis& axis)
    :   m_axis(axis), m_direction(1) {}

    /** Explicit creation of a negative CoordinateDirection from a 
    CoordinateAxis. **/
    CoordinateDirection(const CoordinateAxis& axis, Negative)
    :   m_axis(axis), m_direction(-1) {}

    /** Explicit creation of a CoordinateDirection from a CoordinateAxis
    and a direction calculated at run time.
    @param[in] axis         XAxis, YAxis, or ZAxis
    @param[in] direction    Must be -1 or 1.
    @note Zero is not allowed for \a direction, meaning that
    you must not try to produce one of these from the "sign" result of one of
    the cross product methods, because there the sign can be -1, 0, or 1. **/
    CoordinateDirection(const CoordinateAxis& axis, int direction)
    :   m_axis(axis), m_direction(direction) 
    {   assert(direction==1 || direction==-1); }

    /** This is the coordinate axis XAxis, YAxis, or ZAxis contained in this
    CoordinateDirection.\ Use getDirection() to determine whether this is the
    positive or negative direction. **/
    CoordinateAxis getAxis() const {return m_axis;}
    /** Returns 1 or -1 to indicate the direction along the coordinate
    axis returned by getAxis(). **/
    int getDirection() const {return m_direction;}

    /** Return true if this direction and \a dir2 are along the same axis,
    even if the direction along that axis is not the same. **/
    bool hasSameAxis(const CoordinateDirection& dir2) const
    {   return m_axis.isSameAxis(dir2.getAxis()); }

    /** Return true if this direction and \a dir2 are along the same axis,
    and in the same direction along that axis.\ You can also
    use operator==() for this comparison. **/
    bool isSameAxisAndDirection(const CoordinateDirection& dir2) const
    {   return m_axis==dir2.getAxis() && m_direction==dir2.getDirection(); }

    /** Perform a specialized dot product between this coordinate direction
    and \a dir2; returning 1 or -1 if they contain the same axis and 0 
    otherwise, without performing any floating point operations. **/
    int dotProduct(  const CoordinateDirection& dir2 ) const
    {   if (m_axis != dir2.getAxis()) return 0;
        return m_direction == dir2.getDirection() ? 1 : -1; }

    /** Return the sign that would result from a cross product between this 
    coordinate direction and \a dir2: 0 if they are along the same axis; 
    1 if the result would be in the positive direction along the third axis;
    -1 if it would be in the negative direction. No floating point 
    computations are performed. @see crossProductAxis() **/
    int crossProductSign( const CoordinateDirection& dir2 ) const
    {   if (m_axis == dir2.getAxis()) return 0;
        return m_axis.crossProductSign(dir2.getAxis())
               * m_direction * dir2.getDirection(); }

    /** Return the coordinate axis along which the cross product of this 
    coordinate direction and \a dir2 would lie: same as this if both contain
    the same axis (doesn't matter because the sign would be zero); otherwise,
    the third axis that neither this one nor \a dir2 contains. But note that
    the actual result may be along that axis or in the negative direction 
    along that axis.  No floating point computations are performed. 
    @see crossProductSign(). **/
    CoordinateAxis crossProductAxis( const CoordinateDirection& dir2 ) const
    {   return m_axis.crossProductAxis(dir2.getAxis()); }

    /** Return the axis and sign along that axis that would result from a
    cross product between this coordinate direction and \a dir2; this 
    combines the functions of both crossProductAxis() and crossProductSign().
    Note that if \a dir2 is along the same axis as this one, we'll just
    return this as the axis but the sign is zero since the magnitude of the 
    result would be zero. No floating point calculations are 
    performed. @see crossProductSign(), crossProductAxis() **/
    CoordinateAxis crossProduct( const CoordinateDirection& dir2, int& sign ) const
    {   sign = crossProductSign(dir2); return crossProductAxis(dir2); }

    // Local class declarations for helper classes.
    class NegXDirection; class NegYDirection; class NegZDirection;
private:
    CoordinateAxis m_axis;      // XAxis, YAxis, or ZAxis
    int            m_direction; // 1 or -1
};


// Helper classes that allow compile time recognition of negative axis
// directions.
class CoordinateDirection::NegXDirection : public CoordinateDirection {
  public: NegXDirection() : CoordinateDirection(XAxis,Negative()) {}
};
class CoordinateDirection::NegYDirection : public CoordinateDirection {
  public: NegYDirection() : CoordinateDirection(YAxis,Negative()) {}
};
class CoordinateDirection::NegZDirection : public CoordinateDirection {
  public: NegZDirection() : CoordinateDirection(ZAxis,Negative()) {}
};

// Predefine constants for the negative X,Y,Z directions.
static const CoordinateDirection::NegXDirection  NegXAxis;
static const CoordinateDirection::NegYDirection  NegYAxis;
static const CoordinateDirection::NegZDirection  NegZAxis;

/// Compare two CoordinateDirection objects. @relates CoordinateDirection 
inline bool operator==(const CoordinateDirection& d1, 
                       const CoordinateDirection& d2)
{   return d1.isSameAxisAndDirection(d2); }

/// Compare two CoordinateDirection objects. @relates CoordinateDirection 
inline bool operator!=(const CoordinateDirection& d1, 
                       const CoordinateDirection& d2)
{   return !d1.isSameAxisAndDirection(d2); }

/// Create the NegXAxis direction by negating XAxis. No computation
/// is necessary. @relates CoordinateAxis
inline const CoordinateDirection::NegXDirection&
operator-(const CoordinateAxis::XCoordinateAxis&){return NegXAxis;}
/// Create the NegYAxis direction by negating YAxis. No computation
/// is necessary.  @relates CoordinateAxis
inline const CoordinateDirection::NegYDirection&
operator-(const CoordinateAxis::YCoordinateAxis&){return NegYAxis;}
/// Create the NegZAxis direction by negating ZAxis. No computation
/// is necessary.  @relates CoordinateAxis
inline const CoordinateDirection::NegZDirection&
operator-(const CoordinateAxis::ZCoordinateAxis&){return NegZAxis;}

/// Create the negative direction along the given axis. No computation
/// is necessary.  @relates CoordinateAxis
inline CoordinateDirection
operator-(const CoordinateAxis& axis)
{   return CoordinateDirection(axis,CoordinateDirection::Negative()); }

/// Create the XAxis direction by negating NegXAxis. No computation
/// is necessary. @relates CoordinateDirection
inline const CoordinateAxis::XCoordinateAxis&
operator-(const CoordinateDirection::NegXDirection&){return XAxis;}
/// Create the YAxis direction by negating NegYAxis. No computation
/// is necessary.  @relates CoordinateDirection
inline const CoordinateAxis::YCoordinateAxis&
operator-(const CoordinateDirection::NegYDirection&){return YAxis;}
/// Create the ZAxis direction by negating NegZAxis. No computation
/// is necessary.  @relates CoordinateDirection
inline const CoordinateAxis::ZCoordinateAxis&
operator-(const CoordinateDirection::NegZDirection&){return ZAxis;}

/// Create the opposite direction from the given direction. No computation
/// is necessary.  @relates CoordinateDirection
inline CoordinateDirection
operator-(const CoordinateDirection& dir)
{   return CoordinateDirection(dir.getAxis(), -dir.getDirection()); }

}  // End of namespace

#endif // SimTK_COORDINATEAXIS_H



