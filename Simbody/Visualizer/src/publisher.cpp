#include "publisher.h"
#include <gazebo/gazebo.hh>
#include <iostream>

#include "gazebo/gui/BoxMaker.hh"
#include "gazebo/gui/CylinderMaker.hh"
	   	
using namespace SimTK;

Publisher::Publisher()
{
	gazebo::load();
		// Create our node for communication
		gazebo::transport::NodePtr node(new gazebo::transport::Node());
		node->Init();
	
		// Start transport
		std::cout << "Starting transport" << std::endl;
		gazebo::transport::run();
		std::cout << "Finishing transport" << std::endl;

}

void Publisher::makeBox()
{
	gazebo::gui::BoxMaker maker;
	maker.CreateTheEntity();

}

void Publisher::makeEllipsoid() 
{

}

void Publisher::makeCylinder()
{
	gazebo::gui::CylinderMaker maker;
	maker.CreateTheEntity();	
}

void Publisher::makeCircle()
{

}
void Publisher::makePolygonalMesh()
{

}

/*
void Publisher::Run()
{

	std::cout << "Running \n";
	  // Publisher loop...replace with your own code.
	  while (true)
	  {
		// Throttle Publication
		gazebo::common::Time::MSleep(1000);

		// Generate a pose
		gazebo::math::Pose pose(1, 2, 3, 4, 5, 6);

		// Convert to a pose message
		gazebo::msgs::Pose msg;
		gazebo::msgs::Set(&msg, pose);

		pub->Publish(msg);

		makeBox();
	  }

	  // Make sure to shut everything down.
	  gazebo::transport::fini();
}
*/
