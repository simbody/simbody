#ifndef DYNAMIC_RENDERABLE_H
#define DYNAMIC_RENDERABLE_H
 
#include <OgreSimpleRenderable.h>
 
/// Abstract base class providing mechanisms for dynamically growing hardware buffers.
class DynamicRenderable : public Ogre::SimpleRenderable
{
public:
  /// Constructor
  DynamicRenderable();
  /// Virtual destructor
  virtual ~DynamicRenderable();
 
  /** Initializes the dynamic renderable.
   @remarks
      This function should only be called once. It initializes the
      render operation, and calls the abstract function
      createVertexDeclaration().
   @param operationType The type of render operation to perform.
   @param useIndices Specifies whether to use indices to determine the
          vertices to use as input. */
  void initialize(Ogre::RenderOperation::OperationType operationType,
                  bool useIndices);
 
  /// Implementation of Ogre::SimpleRenderable
  virtual Ogre::Real getBoundingRadius(void) const;
  /// Implementation of Ogre::SimpleRenderable
  virtual Ogre::Real getSquaredViewDepth(const Ogre::Camera* cam) const;
 
protected:
  /// Maximum capacity of the currently allocated vertex buffer.
  size_t mVertexBufferCapacity;
  /// Maximum capacity of the currently allocated index buffer.
  size_t mIndexBufferCapacity;
 
  /** Creates the vertex declaration.
   @remarks
      Override and set mRenderOp.vertexData->vertexDeclaration here.
      mRenderOp.vertexData will be created for you before this method
      is called. */
  virtual void createVertexDeclaration() = 0;
 
  /** Prepares the hardware buffers for the requested vertex and index counts.
   @remarks
      This function must be called before locking the buffers in
      fillHardwareBuffers(). It guarantees that the hardware buffers
      are large enough to hold at least the requested number of
      vertices and indices (if using indices). The buffers are
      possibly reallocated to achieve this.
   @par
      The vertex and index count in the render operation are set to
      the values of vertexCount and indexCount respectively.
   @param vertexCount The number of vertices the buffer must hold.
 
   @param indexCount The number of indices the buffer must hold. This
          parameter is ignored if not using indices. */
  void prepareHardwareBuffers(size_t vertexCount, size_t indexCount);
 
  /** Fills the hardware vertex and index buffers with data.
   @remarks
      This function must call prepareHardwareBuffers() before locking
      the buffers to ensure the they are large enough for the data to
      be written. Afterwards the vertex and index buffers (if using
      indices) can be locked, and data can be written to them. */
  virtual void fillHardwareBuffers() = 0;
};
 
#endif // DYNAMIC_RENDERABLE_H
