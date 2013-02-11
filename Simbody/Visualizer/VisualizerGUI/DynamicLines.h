#ifndef _DYNAMIC_LINES_H_
#define _DYNAMIC_LINES_H_
 
#include "DynamicRenderable.h"
#include <vector>
 
class DynamicLines : public DynamicRenderable
{
  typedef Ogre::Vector3 Vector3;
  typedef Ogre::Quaternion Quaternion;
  typedef Ogre::Camera Camera;
  typedef Ogre::Real Real;
  typedef Ogre::RenderOperation::OperationType OperationType;
 
public:
  /// Constructor - see setOperationType() for description of argument.
  DynamicLines(OperationType opType=Ogre::RenderOperation::OT_LINE_STRIP);
  virtual ~DynamicLines();
 
  /// Add a point to the point list
  void addPoint(const Ogre::Vector3 &p);
  /// Add a point to the point list
  void addPoint(Real x, Real y, Real z);
 
  /// Change the location of an existing point in the point list
  void setPoint(unsigned short index, const Vector3 &value);
 
  /// Return the location of an existing point in the point list
  const Vector3& getPoint(unsigned short index) const;
 
  /// Return the total number of points in the point list
  unsigned short getNumPoints(void) const;
 
  /// Remove all points from the point list
  void clear();
 
  /// Call this to update the hardware buffer after making changes.  
  void update();
 
  /** Set the type of operation to draw with.
   * @param opType Can be one of 
   *    - RenderOperation::OT_LINE_STRIP
   *    - RenderOperation::OT_LINE_LIST
   *    - RenderOperation::OT_POINT_LIST
   *    - RenderOperation::OT_TRIANGLE_LIST
   *    - RenderOperation::OT_TRIANGLE_STRIP
   *    - RenderOperation::OT_TRIANGLE_FAN
   *    The default is OT_LINE_STRIP.
   */
  void setOperationType(OperationType opType);
  OperationType getOperationType() const;
 
protected:
  /// Implementation DynamicRenderable, creates a simple vertex-only decl
  virtual void createVertexDeclaration();
  /// Implementation DynamicRenderable, pushes point list out to hardware memory
  virtual void fillHardwareBuffers();
 
private:
  std::vector<Vector3> mPoints;
  bool mDirty;
};
 
#endif
