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

	void Run();
private:
	std::string GetBoxSDF();

private:
	gazebo::transport::PublisherPtr pub;
	gazebo::transport::PublisherPtr visPub;
	gazebo::transport::PublisherPtr makerPub;

};

}

#endif //_PUBLISHER_H_
