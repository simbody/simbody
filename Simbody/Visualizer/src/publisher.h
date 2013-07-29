#ifndef _PUBLISHER_H_
#define _PUBLISHER_H_

#include <gazebo/transport/transport.hh>
#include <gazebo/msgs/msgs.hh>
#include <math/gzmath.hh>

namespace SimTK 
{

class Publisher 
{

public:
	Publisher();

	void makeBox();
	void makeEllipsoid();
	void makeCylinder();
	void makeCircle();
	void makePolygonalMesh();
};

}

#endif //_PUBLISHER_H_
