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

	void makeBox(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation);
	void makeEllipsoid(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, CustomMesh * mesh);
	void makeCylinder(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, unsigned short resolution);
	void makeCircle();
	void makePolygonalMesh();

private:
	void makeMesh(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, const CustomMesh* mesh, gazebo::msgs::Drawing* drawingMsg);	

private:
	gazebo::transport::NodePtr node;

	gazebo::transport::PublisherPtr drawingPub;
	gazebo::transport::PublisherPtr visPub;
	gazebo::transport::PublisherPtr makerPub;
	gazebo::transport::PublisherPtr requestPub;
};

}

#endif //_PUBLISHER_H_
