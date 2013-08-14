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
                 const Vec4& colour, int representation)
{
//	CustomMesh * mesh = this->makeBox();
//	publisher->makeBox(transform, scale, colour, representation, mesh);
	publisher->makeBox(transform, scale, colour, representation, "Box", 0);
}
void GazeboVisualizer::drawEllipsoid(const Transform& transform, const Vec3& scale, 
		const Vec4& colour, int representation, unsigned short resolution)
{
//	CustomMesh * mesh = this->makeSphere(resolution);
//	publisher->makeEllipsoid(transform, scale, colour, representation, mesh);
	publisher->makeEllipsoid(transform, scale, colour, representation, "Ellipsoid", resolution);
}
void GazeboVisualizer::drawCylinder(const Transform& transform, const Vec3& scale, 
                 const Vec4& colour, int representation, unsigned short resolution)
{
//	CustomMesh * mesh = this->makeCylinder(resolution);
//	publisher->makeCylinder(transform, scale, colour, representation, resolution, mesh);
	publisher->makeCylinder(transform, scale, colour, representation, "Cylinder", resolution);
}
void GazeboVisualizer::drawCircle(const Transform& transform, const Vec3& scale, 
                 const Vec4& colour, int representation, unsigned short resolution)
{
//	CustomMesh * mesh = makeCircle(resolution);
//	publisher->makeCircle(transform, scale, colour, representation, resolution, mesh);
	publisher->makeCircle(transform, scale, colour, representation, "Circle", resolution);
}
void GazeboVisualizer::drawPolygonalMesh(const PolygonalMesh& polyMesh, const Transform& transform, const Vec3& scale, 
                 const Vec4& colour, int representation)
{
	CustomMesh * mesh = this->makePolygonalMesh(polyMesh);
//	publisher->makePolygonalMesh(transform, scale, colour, representation, mesh);
	publisher->makePolygonalMesh(transform, scale, colour, representation, mesh);
}
void GazeboVisualizer::drawLine(const Vec3& end1, const Vec3& end2, 
                 const Vec4& colour, Real thickness)
{
}
void GazeboVisualizer::drawText(const Vec3& position, const Vec3& scale, 
                 const Vec4& colour, const std::string& string, bool faceCamera, bool isScreenText)
{
}
void GazeboVisualizer::drawCoords(const Transform& transform, const Vec3& axisLengths, 
                 const Vec4& colour)
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

CustomMesh* GazeboVisualizer::makeBox()  {
    const float halfx = 1;
    const float halfy = 1;
    const float halfz = 1;
    vector<float> vertices;
    vector<float> normals;
    vector<unsigned short> faces;

    // lower x face
    addVec(vertices, -halfx, -halfy, -halfz);
    addVec(vertices, -halfx, -halfy, halfz);
    addVec(vertices, -halfx, halfy, halfz);
    addVec(vertices, -halfx, halfy, -halfz);
    addVec(normals, -1, 0, 0);
    addVec(normals, -1, 0, 0);
    addVec(normals, -1, 0, 0);
    addVec(normals, -1, 0, 0);
    addVec(faces, 0, 1, 2);
    addVec(faces, 2, 3, 0);
    // upper x face
    addVec(vertices, halfx, halfy, halfz);
    addVec(vertices, halfx, -halfy, halfz);
    addVec(vertices, halfx, -halfy, -halfz);
    addVec(vertices, halfx, halfy, -halfz);
    addVec(normals, 1, 0, 0);
    addVec(normals, 1, 0, 0);
    addVec(normals, 1, 0, 0);
    addVec(normals, 1, 0, 0);
    addVec(faces, 4, 5, 6);
    addVec(faces, 6, 7, 4);
    // lower y face
    addVec(vertices, -halfx, -halfy, -halfz);
    addVec(vertices, halfx, -halfy, -halfz);
    addVec(vertices, halfx, -halfy, halfz);
    addVec(vertices, -halfx, -halfy, halfz);
    addVec(normals, 0, -1, 0);
    addVec(normals, 0, -1, 0);
    addVec(normals, 0, -1, 0);
    addVec(normals, 0, -1, 0);
    addVec(faces, 8, 9, 10);
    addVec(faces, 10, 11, 8);
    // upper y face
    addVec(vertices, halfx, halfy, halfz);
    addVec(vertices, halfx, halfy, -halfz);
    addVec(vertices, -halfx, halfy, -halfz);
    addVec(vertices, -halfx, halfy, halfz);
    addVec(normals, 0, 1, 0);
    addVec(normals, 0, 1, 0);
    addVec(normals, 0, 1, 0);
    addVec(normals, 0, 1, 0);
    addVec(faces, 12, 13, 14);
    addVec(faces, 14, 15, 12);
    // lower z face
    addVec(vertices, -halfx, -halfy, -halfz);
    addVec(vertices, -halfx, halfy, -halfz);
    addVec(vertices, halfx, halfy, -halfz);
    addVec(vertices, halfx, -halfy, -halfz);
    addVec(normals, 0, 0, -1);
    addVec(normals, 0, 0, -1);
    addVec(normals, 0, 0, -1);
    addVec(normals, 0, 0, -1);
    addVec(faces, 16, 17, 18);
    addVec(faces, 18, 19, 16);
    // upper z face
    addVec(vertices, halfx, halfy, halfz);
    addVec(vertices, -halfx, halfy, halfz);
    addVec(vertices, -halfx, -halfy, halfz);
    addVec(vertices, halfx, -halfy, halfz);
    addVec(normals, 0, 0, 1);
    addVec(normals, 0, 0, 1);
    addVec(normals, 0, 0, 1);
    addVec(normals, 0, 0, 1);
    addVec(faces, 20, 21, 22);
    addVec(faces, 22, 23, 20);
    return new CustomMesh(vertices, normals, faces);
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

CustomMesh* GazeboVisualizer::makeCylinder(unsigned short resolution) {
	const int numSides = 6*resolution;
	const float halfHeight = 1;
	const float radius = 1;
	vector<float> vertices;
	vector<float> normals;
	vector<unsigned short> faces;

	// Create the top face.

	addVec(vertices, 0, halfHeight, 0);
	addVec(normals, 0, 1.0, 0);
	for (int i = 0; i < numSides; i++) {
		float theta = (float) ((i*2*SimTK_PI)/numSides);
		float stheta = std::sin(theta);
		float ctheta = std::cos(theta);
		addVec(vertices, radius*ctheta, halfHeight, radius*stheta);
		addVec(normals, 0, 1.0, 0);
	}
	for (int i = 1; i < numSides; i++)
		addVec(faces, i, 0, i+1);
	addVec(faces, numSides, 0, 1);

	// Create the bottom face.

	int bottomStart = numSides+1;
	addVec(vertices, 0, -halfHeight, 0);
	addVec(normals, 0, -1.0, 0);
	for (int i = 0; i < numSides; i++) {
		float theta = (float) ((i*2*SimTK_PI)/numSides);
		float stheta = std::sin(theta);
		float ctheta = std::cos(theta);
		addVec(vertices, radius*ctheta, -halfHeight, radius*stheta);
		addVec(normals, 0, -1.0, 0);
	}
	for (int i = 1; i < numSides; i++)
		addVec(faces, bottomStart+i, bottomStart+i+1, bottomStart);
	addVec(faces, bottomStart+numSides, bottomStart+1, bottomStart);

	// Create the sides.

	for (int i = 0; i < numSides; i++) {
		float theta = (float) ((i*2*SimTK_PI)/numSides);
		float stheta = std::sin(theta);
		float ctheta = std::cos(theta);
		float x = radius*ctheta;
		float z = radius*stheta;
		addVec(vertices, x, halfHeight, z);
		addVec(normals, ctheta, 0, stheta);
		addVec(vertices, x, -halfHeight, z);
		addVec(normals, ctheta, 0, stheta);
	}
	int sideStart = 2*numSides+2;
	for (int i = 0; i < numSides-1; i++) {
		int base = sideStart+2*i;
		addVec(faces, base, base+2, base+1);
		addVec(faces, base+1, base+2, base+3);
	}
	addVec(faces, sideStart+2*numSides-2, sideStart, sideStart+2*numSides-1);
	addVec(faces, sideStart+2*numSides-1, sideStart, sideStart+1);
	return new CustomMesh(vertices, normals, faces);
}

CustomMesh* GazeboVisualizer::makeCircle(unsigned short resolution) 
{
	const int numSides = 6*resolution;
	const float radius = 1;
	vector<float> vertices;
	vector<float> normals;
	vector<unsigned short> faces;

	// Create the front face.

	addVec(vertices, 0, 0, 0);
	addVec(normals, 0, 0, -1.0);
	for (int i = 0; i < numSides; i++) {
		float theta = (float) ((i*2*SimTK_PI)/numSides);
		float stheta = std::sin(theta);
		float ctheta = std::cos(theta);
		addVec(vertices, radius*ctheta, radius*stheta, 0);
		addVec(normals, 0, 0, -1.0);
	}
	for (int i = 1; i < numSides; i++)
		addVec(faces, i, 0, i+1);
	addVec(faces, numSides, 0, 1);

	// Create the back face.

	addVec(vertices, 0, 0, 0);
	addVec(normals, 0, 0, 1.0);
	for (int i = 0; i < numSides; i++) {
		float theta = (float) ((i*2*SimTK_PI)/numSides);
		float stheta = std::sin(theta);
		float ctheta = std::cos(theta);
		addVec(vertices, radius*ctheta, radius*stheta, 0);
		addVec(normals, 0, 0, 1.0);
	}
	int backStart = numSides+1;
	for (int i = 1; i < numSides; i++)
		addVec(faces, backStart, backStart+i, backStart+i+1);
	addVec(faces, backStart, backStart+numSides, backStart+1);
	return new CustomMesh(vertices, normals, faces);
}

CustomMesh * GazeboVisualizer::makePolygonalMesh(const PolygonalMesh& mesh)
{
	vector<float> vertices;
    vector<unsigned short> faces;
    for (int i = 0; i < mesh.getNumVertices(); i++) {
        Vec3 pos = mesh.getVertexPosition(i);
        vertices.push_back((float) pos[0]);
        vertices.push_back((float) pos[1]);
        vertices.push_back((float) pos[2]);
    }
    for (int i = 0; i < mesh.getNumFaces(); i++) {
        int numVert = mesh.getNumVerticesForFace(i);
        if (numVert < 3)
            continue; // Ignore it.
        if (numVert == 3) {
            faces.push_back((unsigned short) mesh.getFaceVertex(i, 0));
            faces.push_back((unsigned short) mesh.getFaceVertex(i, 1));
            faces.push_back((unsigned short) mesh.getFaceVertex(i, 2));
        }
        else if (numVert == 4) {
            // Split it into two triangles.

            faces.push_back((unsigned short) mesh.getFaceVertex(i, 0));
            faces.push_back((unsigned short) mesh.getFaceVertex(i, 1));
            faces.push_back((unsigned short) mesh.getFaceVertex(i, 2));
            faces.push_back((unsigned short) mesh.getFaceVertex(i, 2));
            faces.push_back((unsigned short) mesh.getFaceVertex(i, 3));
            faces.push_back((unsigned short) mesh.getFaceVertex(i, 0));
        }
        else {
            // Add a vertex at the center, then split it into triangles.

            Vec3 center(0);
            for (int j = 0; j < numVert; j++) {
                Vec3 pos = mesh.getVertexPosition(mesh.getFaceVertex(i,j));
                center += pos;
            }
            center /= numVert;
            vertices.push_back((float) center[0]);
            vertices.push_back((float) center[1]);
            vertices.push_back((float) center[2]);
            const unsigned newIndex = vertices.size()/3-1;
            for (int j = 0; j < numVert-1; j++) {
                faces.push_back((unsigned short) mesh.getFaceVertex(i, j));
                faces.push_back((unsigned short) mesh.getFaceVertex(i, j+1));
                faces.push_back((unsigned short) newIndex);
            }
        }
    }
}
