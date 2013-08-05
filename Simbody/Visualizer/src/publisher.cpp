#include "publisher.h"
#include <gazebo/gazebo.hh>
#include <iostream>

#include "gazebo/gui/MeshMaker.hh"
#include "gazebo/math/Quaternion.hh"

using namespace SimTK;
using namespace gazebo;

Publisher::Publisher()
{
		gazebo::load();
		// Create our node for communication
		this->node = transport::NodePtr(new gazebo::transport::Node());
		node->Init();
	
		gazebo::transport::run();

		this->drawPub = this->node->Advertise<msgs::Drawing>("~/draw");
		this->visPub = this->node->Advertise<msgs::Visual>("~/visual");
		this->makerPub = this->node->Advertise<msgs::Factory>("~/factory");
		this->requestPub = this->node->Advertise<msgs::Request>("~/request");

}

void Publisher::makeBox(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation)
{
//	gazebo::gui::BoxMaker maker;
//	maker.CreateTheEntity();	
	
	// Simbody geometry to gazebo 
	// Position
//	Vec3 rot = transform.R().convertRotationToBodyFixedXYZ();

//	position.updR().setRotationToBodyFixedXYZ(fVec3(rot[0], rot[1], rot[2]));
	
	msgs::Visual* visualMsg = new msgs::Visual();
	visualMsg->set_name("visuals");
	visualMsg->set_parent_name("default");

	// Position is set by Pose-Vector3d
	visualMsg->mutable_geometry()->set_type(msgs::Geometry::BOX);
	math::Vector3 position(transform.p()[0], transform.p()[1], transform.p()[2]);

	msgs::Set(visualMsg->mutable_pose()->mutable_position(), position);
	msgs::Set(visualMsg->mutable_pose()->mutable_orientation(), math::Quaternion());
	
	// Setting scale
	math::Vector3 scale_g(scale[0], scale[1], scale[2]);
	msgs::Set(visualMsg->mutable_geometry()->mutable_box()->mutable_size(), scale_g);

	// Setting colour
	common::Color colour_g(colour[0], colour[1], colour[2]);
	msgs::Set(visualMsg->mutable_material()->mutable_ambient(), colour_g);

	std::cout << "Making box" << std::endl;
	visPub->Publish(*visualMsg);
}

void Publisher::makeEllipsoid() 
{

}

void Publisher::makeCylinder(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, unsigned short resolution)
{
	msgs::Visual* visualMsg = new msgs::Visual();
	visualMsg->set_name("visuals");
	visualMsg->set_parent_name("default");

	// Position is set by Pose-Vector3d
	visualMsg->mutable_geometry()->set_type(msgs::Geometry::CYLINDER);
	math::Vector3 position(transform.p()[0], transform.p()[1], transform.p()[2]);

	msgs::Set(visualMsg->mutable_pose()->mutable_position(), position);
	msgs::Set(visualMsg->mutable_pose()->mutable_orientation(), math::Quaternion());
	
	// Setting scale
	//math::Vector3 scale_g(scale[0], scale[1], scale[2]);
	visualMsg->mutable_geometry()->mutable_cylinder()->set_radius(scale[0]);
	visualMsg->mutable_geometry()->mutable_cylinder()->set_length(scale[1]);

	// Setting colour
	common::Color colour_g(colour[0], colour[1], colour[2]);
	msgs::Set(visualMsg->mutable_material()->mutable_ambient(), colour_g);

	std::cout << "Making cylinder" << std::endl;
	visPub->Publish(*visualMsg);
}

void Publisher::makeCircle()
{

}
void Publisher::makePolygonalMesh()
{
//	std::cout << "Making mesh" << std::endl;
//	gazebo::gui::MeshMaker maker;
//	maker.CreateTheEntity();
}
