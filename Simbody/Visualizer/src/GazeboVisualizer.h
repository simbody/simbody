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

	virtual void shutdownGUI() {};
	virtual void shakeHandsWithGUI(int inPipe, int outPipe) {}
	virtual void beginScene(Real simTime);
	virtual void finishScene();
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

	virtual void setMaxFrameRate(Real rateInFPS) const {}
    virtual void setBackgroundColor(const Vec3& color) const {}
    virtual void setBackgroundType(Visualizer::BackgroundType type) const {}
    virtual void setSystemUpDirection(const CoordinateDirection& upDir) {}
    virtual void addMenu(const String& title, int id, 
                 const Array_<std::pair<String, int> >& items) {}

public:
	void initServer();
	void initPub();
private:
	FakeServer * server;
	Publisher * publisher;

	CustomMesh* makeSphere(unsigned short resolution);
	CustomMesh* makeBox();
	CustomMesh* makeCylinder(unsigned short resolution);
	CustomMesh* makeCircle(unsigned short resolution);
	CustomMesh* makePolygonalMesh(const PolygonalMesh& mesh);

	void addVec(std::vector<float>& data, float x, float y, float z);
	void addVec(std::vector<unsigned short>& data, float x, float y, float z);

};

}
#endif //_GAZEBO_VISUALIZER_H_
