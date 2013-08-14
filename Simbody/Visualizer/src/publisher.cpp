#include "publisher.h"
#include <gazebo/gazebo.hh>
#include <iostream>

#include "gazebo/gui/MeshMaker.hh"
#include "gazebo/math/Quaternion.hh"

using namespace SimTK;
using namespace gazebo;
using namespace std;

string getRandomName() {
    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";

	std::string s;

    for (int i = 0; i < 5; ++i) {
        s.append(1, alphanum[rand() % (sizeof(alphanum) - 1)]);
    }

    return s;
}

Publisher::Publisher()
{
	gazebo::load();
	// Create our node for communication
	this->node = transport::NodePtr(new gazebo::transport::Node());
	node->Init();

	gazebo::transport::run();

	this->drawingPub = this->node->Advertise<msgs::Drawing>("~/draw");
	this->visPub = this->node->Advertise<msgs::Visual>("~/visual");
	this->makerPub = this->node->Advertise<msgs::Factory>("~/factory");
	this->requestPub = this->node->Advertise<msgs::Request>("~/request");

}

void Publisher::makeBox(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, const std::string& geometry, int resolution)
{
	msgs::Visual* visualMsg = new msgs::Visual();
	prepareMessage(transform, scale, colour, representation, geometry, resolution, visualMsg);

	visualMsg->mutable_geometry()->set_type(msgs::Geometry::BOX);
	// Setting scale
	math::Vector3 scale_g(scale[0], scale[1], scale[2]);
	msgs::Set(visualMsg->mutable_geometry()->mutable_box()->mutable_size(), scale_g);

	visPub->Publish(*visualMsg);
}

void Publisher::makeEllipsoid(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, const std::string& geometry, int resolution)
{
	msgs::Visual* visualMsg = new msgs::Visual();
	prepareMessage(transform, scale, colour, representation, geometry, resolution, visualMsg);

	// Position is set by Pose-Vector3d
	visualMsg->mutable_geometry()->set_type(msgs::Geometry::SPHERE);
	// Setting scale
	math::Vector3 scale_g(scale[0], scale[1], scale[2]);
	visualMsg->mutable_geometry()->mutable_sphere()->set_radius(scale[0]);

	visPub->Publish(*visualMsg);

}
void Publisher::makeCylinder(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, const std::string& geometry, unsigned short resolution)
{
	msgs::Visual* visualMsg = new msgs::Visual();
	prepareMessage(transform, scale, colour, representation, geometry, resolution, visualMsg);

	// Position is set by Pose-Vector3d
	visualMsg->mutable_geometry()->set_type(msgs::Geometry::CYLINDER);
	visualMsg->mutable_geometry()->mutable_cylinder()->set_radius(scale[0]);
	visualMsg->mutable_geometry()->mutable_cylinder()->set_length(2*scale[1]); //half height?

	visPub->Publish(*visualMsg);

}
void Publisher::makeCircle(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, const std::string& geometry, unsigned short resolution)
{
	msgs::Visual* visualMsg = new msgs::Visual();
	prepareMessage(transform, scale, colour, representation, geometry, resolution, visualMsg);

	// Position is set by Pose-Vector3d
	visualMsg->mutable_geometry()->set_type(msgs::Geometry::PLANE);
	msgs::Set(visualMsg->mutable_geometry()->mutable_plane()->mutable_normal(), math::Vector3(0,0,1));
	msgs::Set(visualMsg->mutable_geometry()->mutable_plane()->mutable_size(), math::Vector2d(scale[0], scale[1]));

	visPub->Publish(*visualMsg);

}

void Publisher::prepareMessage(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, const std::string& geometry, unsigned short resolution, msgs::Visual* visualMsg)
{	

	std::string name = "visuals" + geometry;
	ostringstream convert;
	convert << name << resolution;
	name = convert.str();

	std::cout << "Sending visual: " + name << std::endl;

	visualMsg->set_name(name);
	visualMsg->set_parent_name("default");

	// Simbody geometry to gazebo 
	// Position is set by Pose-Vector3d
	math::Vector3 position(transform.p()[0], transform.p()[1], transform.p()[2]);

	Vec3 rot = transform.R().convertRotationToBodyFixedXYZ();
	math::Vector3 rotation(rot[0], rot[1], rot[2]);

	msgs::Set(visualMsg->mutable_pose()->mutable_position(), position);
	msgs::Set(visualMsg->mutable_pose()->mutable_orientation(), rotation);

	// Setting colour
	std::cout << "transparency " << colour[3] << std::endl;
//	common::Color colour_g(colour[0], colour[1], colour[2], colour[3]);
//	msgs::Set(visualMsg->mutable_material()->mutable_ambient(), colour_g);
//	visualMsg->set_transparency(colour[3]);

}

void Publisher::makeBox(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, CustomMesh* mesh)
{

	/*
	sleep(1);
	msgs::Request * requestMsg = msgs::CreateRequest("entity_delete", name);
	requestPub->Publish(*requestMsg);
*/

/*	
	msgs::Drawing* drawingMsg = new msgs::Drawing();
	drawingMsg->set_name("box" + getRandomName());
	drawingMsg->set_visible(true);
	this->makeMesh(transform, scale, colour, representation, mesh, drawingMsg);
	drawingPub->Publish(*drawingMsg);
*/	
	
}

void Publisher::makeEllipsoid(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, CustomMesh * mesh) 
{

	/*
	sleep(1);
	msgs::Request * requestMsg = msgs::CreateRequest("entity_delete", name);
	requestPub->Publish(*requestMsg);
*/
	/*	
	msgs::Drawing* drawingMsg = new msgs::Drawing();	
	
	drawingMsg->set_name("ellipsoid" + getRandomName());
	drawingMsg->set_visible(true);
	this->makeMesh(transform, scale, colour, representation, mesh, drawingMsg);
	drawingPub->Publish(*drawingMsg);
*/	

}

void Publisher::makeCylinder(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, unsigned short resolution, CustomMesh * mesh)
{
	/*
	sleep(1);
	msgs::Request * requestMsg = msgs::CreateRequest("entity_delete", name);
	requestPub->Publish(*requestMsg);
*/	/*
	msgs::Drawing* drawingMsg = new msgs::Drawing();
	drawingMsg->set_name("cylinder" + getRandomName());
	drawingMsg->set_visible(true);
	this->makeMesh(transform, scale, colour, representation, mesh, drawingMsg);
	drawingPub->Publish(*drawingMsg);
	*/
	
}

void Publisher::makeCircle(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, unsigned short resolution, CustomMesh * mesh)
{

	/*
	sleep(1);
	msgs::Request * requestMsg = msgs::CreateRequest("entity_delete", name);
	requestPub->Publish(*requestMsg);
*/
	/*
	msgs::Drawing* drawingMsg = new msgs::Drawing();
	drawingMsg->set_name("circle" + getRandomName());
	drawingMsg->set_visible(true);
	this->makeMesh(transform, scale, colour, representation, mesh, drawingMsg);
	drawingPub->Publish(*drawingMsg);
*/
}
void Publisher::makePolygonalMesh(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, CustomMesh* mesh)
{
//	std::cout << "Making mesh" << std::endl;
//	gazebo::gui::MeshMaker maker;
//	maker.CreateTheEntity();
	msgs::Drawing* drawingMsg = new msgs::Drawing();
	drawingMsg->set_name("polygonMesh" + getRandomName());
	drawingMsg->set_visible(true);
	std::cout << "Making polygonal mesh" << std::endl;
	this->makeMesh(transform, scale, colour, representation, mesh, drawingMsg);
	drawingPub->Publish(*drawingMsg);
}

void Publisher::makeMesh(const Transform& transform, const Vec3& scale, const Vec4& colour, int representation, const CustomMesh* mesh, msgs::Drawing* drawingMsg)
{
	
	Vec3 rot = transform.R().convertRotationToBodyFixedXYZ();

	msgs::Set(drawingMsg->mutable_transform()->mutable_translate(), math::Vector3(transform.p()[0],transform.p()[1],transform.p()[2]));
	msgs::Set(drawingMsg->mutable_transform()->mutable_rotate(), math::Vector3(rot[0], rot[1], rot[2]));
	msgs::Set(drawingMsg->mutable_transform()->mutable_scale(), math::Vector3(scale[0], scale[1], scale[2]));

	//msgs::Set(drawingMsg->mutable_pose()->mutable_position(), math::Vector3(transform.p()[0], transform.p()[1], transform.p()[2]));

	//msgs::Set(drawingMsg->mutable_pose()->mutable_orientation(), math::Quaternion());

	common::Color colour_g(colour[0], colour[1], colour[2], colour[3]);

	if(representation == DecorativeGeometry::DrawSurface)
	{
		drawingMsg->set_mode(msgs::Drawing::TRIANGLE_LIST);

		const vector<unsigned short> faces = mesh->getFaces();
		const vector<float> vertices = mesh->getVertices();
		const vector<float> normals = mesh->getNormals();

		for(int i=0; i < vertices.size(); i+=3)
		{	
			msgs::Drawing::Point* p = drawingMsg->add_point();
			msgs::Set(p->mutable_position(), math::Vector3(vertices.at(i), vertices.at(i+1), vertices.at(i+2)));
			msgs::Set(p->mutable_color(), colour_g);
		}

		for(int i=0; i < normals.size(); i+=3)
		{
			msgs::Drawing::Normal* n = drawingMsg->add_normal();
			msgs::Set(n->mutable_point(), math::Vector3(normals.at(i), normals.at(i+1), normals.at(i+2)));
		}

		for(int i=0; i < faces.size(); i+=3)
		{
			msgs::Drawing::Face* f = drawingMsg->add_edge();
			msgs::Set(f->mutable_edge(), math::Vector3(faces.at(i), faces.at(i+1), faces.at(i+2)));
		}

	}
	else if(representation == DecorativeGeometry::DrawPoints)
	{
		drawingMsg->set_mode(msgs::Drawing::POINT_LIST);

		vector<float> vertices = mesh->getVertices();
		for(int i=0; i < vertices.size(); i+=3)
		{	
			msgs::Drawing::Point* p = drawingMsg->add_point();
			msgs::Set(p->mutable_position(), math::Vector3(vertices.at(i), vertices.at(i+1), vertices.at(i+2)));
			msgs::Set(p->mutable_color(), colour_g);
		}

	}
	else if(representation == DecorativeGeometry::DrawWireframe)
	{
		drawingMsg->set_mode(msgs::Drawing::LINE_LIST);

		vector<float> vertices = mesh->getVertices();
		vector<float> normals = mesh->getNormals();

		for(int i=0; i < vertices.size(); i+=3)
		{	
			msgs::Drawing::Point* p = drawingMsg->add_point();
			msgs::Set(p->mutable_position(), math::Vector3(vertices.at(i), vertices.at(i+1), vertices.at(i+2)));
			msgs::Set(p->mutable_color(), colour_g);
		}

		for(int i=0; i < normals.size(); i+=3)
		{
			msgs::Drawing::Normal* n = drawingMsg->add_normal();
			msgs::Set(n->mutable_point(), math::Vector3(normals.at(i), normals.at(i+1), normals.at(i+2)));
		}

	}
	else
	{
		assert(false);
		drawingMsg->set_visible(false);
	}

}

