#include "publisher.h"
#include <gazebo/gazebo.hh>
#include <iostream>

#include "gazebo/gui/BoxMaker.hh"
#include "gazebo/gui/CylinderMaker.hh"
	   	
using namespace SimTK;

Publisher::Publisher()
{
		// Create our node for communication
		gazebo::transport::NodePtr node(new gazebo::transport::Node());
		node->Init();
	
		// Start transport
//		gazebo::transport::run();

		// Publish to a Gazebo topic
		pub = node->Advertise<gazebo::msgs::Pose>("~/pose_example");

		visPub = node->Advertise<gazebo::msgs::Visual>("~/visual");

		makerPub = node->Advertise<gazebo::msgs::Factory>("~/factory");
}

void Publisher::makeBox()
{
	gazebo::gui::BoxMaker maker;
	maker.CreateTheEntity();

	/*
		std::cout << "Drawing boxes \n";
		gazebo::msgs::Visual* visualMsg = new gazebo::msgs::Visual();
		visualMsg->mutable_geometry()->set_type(gazebo::msgs::Geometry::BOX);
		visualMsg->mutable_material()->mutable_script()->set_name("Gazebo/TurquoiseGlowOutline");

		gazebo::msgs::Factory msg;
		gazebo::math::Vector3 p = gazebo::msgs::Convert(visualMsg->pose().position());
		gazebo::math::Vector3 size = gazebo::msgs::Convert(visualMsg->geometry().box().size());

		msg.set_sdf(this->GetBoxSDF());

		this->makerPub->Publish(msg);
		*/
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
	/////////////////////////////////////////////////
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
