#ifndef _GAZEBO_VISUALIZER_H_
#define _GAZEBO_VISUALIZER_H_

#include "VisualizerProtocol.h"
#include "server.h"
#include "publisher.h"

// Gazebo visualizer. This class does not perform any drawing, except it'll be in charge of sending messages
namespace SimTK {

class GazeboVisualizer : public VisualizerProtocol
{
public:
	GazeboVisualizer();
	virtual ~GazeboVisualizer();

public:

    virtual void drawBox(const Transform& transform, const Vec3& scale, 
                 const Vec4& color, int representation);
    virtual void drawEllipsoid(const Transform& transform, const Vec3& scale, 
                       const Vec4& color, int representation, 
                       unsigned short resolution);
    virtual void drawCylinder(const Transform& transform, const Vec3& scale, 
                      const Vec4& color, int representation, 
                      unsigned short resolution);
    virtual void drawCircle(const Transform& transform, const Vec3& scale, 
                    const Vec4& color, int representation, 
                    unsigned short resolution);
    virtual void drawPolygonalMesh(const PolygonalMesh& mesh, 
                           const Transform& transform, const Vec3& scale, 
                           const Vec4& color, int representation);
    virtual void drawLine(const Vec3& end1, const Vec3& end2, const 
                  Vec4& color, Real thickness);
    virtual void drawText(const Vec3& position, const Vec3& scale, const Vec4& color, 
                  const std::string& string, bool faceCamera, bool isScreenText);
    virtual void drawCoords(const Transform& transform, const Vec3& axisLengths, 
                    const Vec4& color);

	virtual void drawMesh();

private:
	void initServer();
	void initPub();
private:
	FakeServer * server;
	Publisher * publisher;
};

}
#endif //_GAZEBO_VISUALIZER_H_
