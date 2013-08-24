#ifndef _PUBLISHER_H_
#define _PUBLISHER_H_

#include <gazebo/transport/transport.hh>
#include <gazebo/msgs/msgs.hh>
#include <math/gzmath.hh>

#include "simbody/internal/common.h"
#include "SimbodyScene.h"

namespace SimTK 
{

class Publisher 
{

public:
	Publisher();

//	These functions are used for custom drawing messages
	void makeBox(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, CustomMesh* mesh);
	void makeEllipsoid(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, CustomMesh * mesh);
	void makeCylinder(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, unsigned short resolution, CustomMesh * mesh);
	void makeCircle(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, unsigned short resolution, CustomMesh * mesh);
	void makePolygonalMesh(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, CustomMesh * mesh);

	void makeLine(const Vec3& end1, const Vec3& end2, const Vec4 colour, Real thickness);

public:
//	These functions are used for drawing using visual messages
	void makeBox(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, const std::string& geometry, int resolution);
	void makeEllipsoid(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, const std::string& geometry, int resolution);
	void makeCylinder(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, const std::string& geometry, unsigned short resolution);
	void makeCircle(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, const std::string& geometry, unsigned short resolution);

private:
	void makeMesh(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, const CustomMesh* mesh, gazebo::msgs::Drawing* drawingMsg);	

	void prepareMessage(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, const std::string& geometry, unsigned short resolution, gazebo::msgs::Visual* visualMsg); 

private:
	gazebo::transport::NodePtr node;

	gazebo::transport::PublisherPtr drawingPub;
	gazebo::transport::PublisherPtr visPub;
	gazebo::transport::PublisherPtr requestPub;
};

}

#endif //_PUBLISHER_H_
