#include "GazeboVisualizer.h"
#include "gazebo/gazebo.hh"

#include "SimTKcommon.h"

using namespace SimTK;
using namespace gazebo;
using namespace std;

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
}

void GazeboVisualizer::drawBox(const Transform& transform, const Vec3& scale, 
                 const Vec4& color, int representation)
{
	publisher->makeBox(transform, scale, color, representation);
}
void GazeboVisualizer::drawEllipsoid(const Transform& transform, const Vec3& scale, 
		const Vec4& color, int representation, unsigned short resolution)
{
	CustomMesh * mesh = makeSphere(resolution);

	publisher->makeEllipsoid(transform, scale, color, representation, mesh);

	delete mesh;
}
void GazeboVisualizer::drawCylinder(const Transform& transform, const Vec3& scale, 
                 const Vec4& color, int representation, unsigned short resolution)
{
	publisher->makeCylinder(transform, scale, color, representation, resolution);
}
void GazeboVisualizer::drawCircle(const Transform& transform, const Vec3& scale, 
                 const Vec4& color, int representation, unsigned short resolution)
{
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

void GazeboVisualizer::addVec(vector<float>& data, float x, float y, float z) {
	data.push_back(x);
	data.push_back(y);
	data.push_back(z);
}

void GazeboVisualizer::addVec(vector<unsigned short>& data, float x, float y, float z) {
	data.push_back((unsigned short) x);
	data.push_back((unsigned short) y);
	data.push_back((unsigned short) z);
}

CustomMesh* GazeboVisualizer::makeSphere(unsigned short resolution) {
    const int numLatitude = 4*resolution;
    const int numLongitude = 6*resolution;
    const float radius = 1.0f;
    vector<float> vertices;
    vector<float> normals;
    vector<unsigned short> faces;
    addVec(vertices, 0, radius, 0);
    addVec(normals, 0, 1, 0);
    for (int i = 0; i < numLatitude; i++) {
        float phi = (float) (((i+1)*SimTK_PI)/(numLatitude+1));
        float sphi = std::sin(phi);
        float cphi = std::cos(phi);
        float y = radius*cphi;
        float r = radius*sphi;
        for (int j = 0; j < numLongitude; j++) {
            float theta = (float) ((j*2*SimTK_PI)/numLongitude);
            float stheta = std::sin(theta);
            float ctheta = std::cos(theta);
            addVec(vertices, r*ctheta, y, r*stheta);
            addVec(normals, sphi*ctheta, cphi, sphi*stheta);
        }
    }
    addVec(vertices, 0, -radius, 0);
    addVec(normals, 0, -1, 0);
    for (int i = 1; i < numLongitude; i++)
        addVec(faces, 0, i+1, i);
    addVec(faces, 0, 1, numLongitude);
    for (int i = 1; i < numLatitude; i++) {
        int base = (i-1)*numLongitude+1;
        for (int j = 0; j < numLongitude; j++) {
            int v1 = base+j;
            int v2 = (j == numLongitude-1 ? base : v1+1);
            int v3 = v1+numLongitude;
            int v4 = v2+numLongitude;
            addVec(faces, v1, v4, v3);
            addVec(faces, v1, v2, v4);
        }
    }
    int first = (numLatitude-1)*numLongitude+1;
    int last = numLatitude*numLongitude+1;
    for (int i = first; i < last-1; i++)
        addVec(faces, i, i+1, last);
    addVec(faces, last-1, first, last);
    return new CustomMesh(vertices, normals, faces);
}


