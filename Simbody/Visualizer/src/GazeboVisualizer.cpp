#include "GazeboVisualizer.h"
#include "gazebo/gazebo.hh"

using namespace SimTK;
using namespace gazebo;

// Create gazebo server
GazeboVisualizer::GazeboVisualizer()
{
//	gazebo::load();

//    pthread_t serverThread;
//    pthread_create(&serverThread, NULL, initP, this);

//	initServer();
//	std::cout << "server init complete" << std::endl;

	initPub();
	std::cout << "publisher init complete" << std::endl;
}

GazeboVisualizer::~GazeboVisualizer()
{
	delete server;
	delete publisher;
}

void GazeboVisualizer::beginScene(Real simTime)
{
}

void GazeboVisualizer::finishScene()
{
}

void GazeboVisualizer::initServer()
{
	server = new FakeServer();

	gazebo::common::Console::Instance()->Init("simbodyserver.log");
	if(!server->ParseArgs()) std::cout << "Parse args error\n";

	server->Run();
}
 
void GazeboVisualizer::initPub()
{
	std::cout << "initializing publisher" << std::endl;
	publisher = new Publisher();
//	publisher->makeBox();
}

void GazeboVisualizer::drawBox(const Transform& transform, const Vec3& scale, 
                 const Vec4& color, int representation)
{
	publisher->makeBox(transform, scale, color, representation);
}
void GazeboVisualizer::drawEllipsoid(const Transform& transform, const Vec3& scale, 
		const Vec4& color, int representation, unsigned short resolution)
{
	publisher->makeEllipsoid();
}
void GazeboVisualizer::drawCylinder(const Transform& transform, const Vec3& scale, 
                 const Vec4& color, int representation, unsigned short resolution)
{
	publisher->makeCylinder(transform, scale, color, representation, resolution);
}
void GazeboVisualizer::drawCircle(const Transform& transform, const Vec3& scale, 
                 const Vec4& color, int representation, unsigned short resolution)
{
	publisher->makeCircle();
}
void GazeboVisualizer::drawPolygonalMesh(const PolygonalMesh& mesh, const Transform& transform, const Vec3& scale, 
                 const Vec4& color, int representation)
{
	publisher->makePolygonalMesh();
}
void GazeboVisualizer::drawLine(const Vec3& end1, const Vec3& end2, 
                 const Vec4& color, Real thickness)
{
}
void GazeboVisualizer::drawText(const Vec3& position, const Vec3& scale, 
                 const Vec4& color, const std::string& string, bool faceCamera, bool isScreenText)
{

}
void GazeboVisualizer::drawCoords(const Transform& transform, const Vec3& axisLengths, 
                 const Vec4& color)
{

}
void GazeboVisualizer::drawMesh()
{
}

