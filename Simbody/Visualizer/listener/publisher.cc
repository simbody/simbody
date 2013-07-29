/*
 * Copyright 2012 Open Source Robotics Foundation
 * Copyright 2013 Dereck Wonnacott
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
*/

#include <gazebo/gazebo.hh>
#include <gazebo/transport/transport.hh>
#include <gazebo/msgs/msgs.hh>
#include <math/gzmath.hh>

#include <iostream>

class FakePublisher
{

public:
	   	
	FakePublisher() 
	{
	  // Load gazebo
	  gazebo::load();

		// Create our node for communication
		gazebo::transport::NodePtr node(new gazebo::transport::Node());
		node->Init();
	
		// Start transport
		gazebo::transport::run();

		// Publish to a Gazebo topic
		pub = node->Advertise<gazebo::msgs::Pose>("~/pose_example");

		visPub = node->Advertise<gazebo::msgs::Visual>("~/visual");

		makerPub = node->Advertise<gazebo::msgs::Factory>("~/factory");
	  // Wait for a subscriber to connect
//	  pub->WaitForConnection();

//	  visPub->WaitForConnection();
//	  makerPub->WaitForConnection();
	}

	void makeBox()
	{
		std::cout << "Drawing boxes \n";
		gazebo::msgs::Visual* visualMsg = new gazebo::msgs::Visual();
		visualMsg->mutable_geometry()->set_type(gazebo::msgs::Geometry::BOX);
		visualMsg->mutable_material()->mutable_script()->set_name("Gazebo/TurquoiseGlowOutline");

		gazebo::msgs::Factory msg;
		gazebo::math::Vector3 p = gazebo::msgs::Convert(visualMsg->pose().position());
		gazebo::math::Vector3 size = gazebo::msgs::Convert(visualMsg->geometry().box().size());

		msg.set_sdf(this->GetSDFString());

		this->makerPub->Publish(msg);

	}

	std::string GetSDFString()
	{
		int counter = 0;
		std::ostringstream newModelStr;
		newModelStr << "<sdf version ='" << "1.4" << "'>"
			<< "<model name='unit_box_" << counter << "'>"
			<< "<pose>0 0 0.5 0 0 0</pose>"
			<< "<link name ='link'>"
//			<<   "<inertial><mass>1.0</mass></inertial>"
//			<<   "<collision name ='collision'>"
//			<<     "<geometry>"
//			<<       "<box>"
//			<<         "<size>1.0 1.0 1.0</size>"
//			<<       "</box>"
//			<<     "</geometry>"
//			<< "</collision>"
			<< "<visual name ='visual'>"
			<<     "<geometry>"
			<<       "<box>"
			<<         "<size>1.0 1.0 1.0</size>"
			<<       "</box>"
			<<     "</geometry>"
			<<     "<material>"
			<<       "<script>"
			<<         "<uri>file://media/materials/scripts/gazebo.material</uri>"
			<<         "<name>Gazebo/Grey</name>"
			<<       "</script>"
			<<     "</material>"
			<<   "</visual>"
			<< "</link>"
			<< "<static>true</static>"
			<< "</model>"
			<< "</sdf>";

  		return newModelStr.str();
	}

	/////////////////////////////////////////////////
	void Run()
	{

		std::cout << "Running \n";
	  // Publisher loop...replace with your own code.
	  while (true)
	  {
		// Throttle Publication
		gazebo::common::Time::MSleep(100);

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
private:
	
	gazebo::transport::PublisherPtr pub;
	gazebo::transport::PublisherPtr visPub;
	gazebo::transport::PublisherPtr makerPub;
};

int main(int argc, char ** argv)
{
	FakePublisher pub;
	pub.Run();	
}

