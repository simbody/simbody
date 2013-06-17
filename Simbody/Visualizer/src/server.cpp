/*
 * Copyright 2012 Open Source Robotics Foundation
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
#include <stdio.h>
#include <signal.h>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "gazebo/gazebo.hh"
#include "gazebo/transport/transport.hh"

#include "gazebo/common/Timer.hh"
#include "gazebo/common/Exception.hh"
#include "gazebo/common/Plugin.hh"
#include "gazebo/common/Common.hh"

#include "gazebo/sdf/sdf.hh"

#include "gazebo/physics/PhysicsFactory.hh"
#include "gazebo/physics/Physics.hh"
#include "gazebo/physics/World.hh"
#include "gazebo/physics/Base.hh"

#include "gazebo/Master.hh"
#include "server.h"


using namespace gazebo;

using namespace SimTK;

bool FakeServer::stop = true;

/////////////////////////////////////////////////
FakeServer::FakeServer()
{
  this->receiveMutex = new boost::mutex();
  gazebo::print_version();

  if (signal(SIGINT, FakeServer::SigInt) == SIG_ERR)
    std::cerr << "signal(2) failed while setting up for SIGINT" << std::endl;
  
//  _publisher = new FakePublisher();
}

/////////////////////////////////////////////////
FakeServer::~FakeServer()
{
  fflush(stdout);
  delete this->receiveMutex;
  delete this->master;
}

/////////////////////////////////////////////////
void FakeServer::PrintUsage()
{
  std::cerr << "Run the OpenSim server.\n\n"
    << "Usage: opensimserver [options] <world_file>\n\n";
}

/////////////////////////////////////////////////
bool FakeServer::ParseArgs(int argc, char **argv)
{
	// Get the world file name from the command line, or use "empty.world"
	// if no world file is specified.
	std::string configFilename = "worlds/empty.world";

	// Get the physics engine name specified from the command line, or use ""
	// if no physics engine is specified.
	std::string physics;

	// Load the server
	if (!this->LoadFile(configFilename, physics))
	  return false;

  this->ProcessParams();
  this->Init();

  return true;
}

/////////////////////////////////////////////////
bool FakeServer::GetInitialized() const
{
  return !this->stop && !transport::is_stopped();
}

/////////////////////////////////////////////////
bool FakeServer::LoadFile(const std::string &_filename,
                      const std::string &_physics)
{
  // Quick test for a valid file
  FILE *test = fopen(common::find_file(_filename).c_str(), "r");
  if (!test)
  {
    gzerr << "Could not open file[" << _filename << "]\n";
    return false;
  }
  fclose(test);

  // Load the world file
  sdf::SDFPtr sdf(new sdf::SDF);
  if (!sdf::init(sdf))
  {
    gzerr << "Unable to initialize sdf\n";
    return false;
  }

  if (!sdf::readFile(_filename, sdf))
  {
    gzerr << "Unable to read sdf file[" << _filename << "]\n";
    return false;
  }

  return this->LoadImpl(sdf->root, _physics);
}

/////////////////////////////////////////////////
bool FakeServer::LoadString(const std::string &_sdfString)
{
  // Load the world file
  sdf::SDFPtr sdf(new sdf::SDF);
  if (!sdf::init(sdf))
  {
    gzerr << "Unable to initialize sdf\n";
    return false;
  }

  if (!sdf::readString(_sdfString, sdf))
  {
    gzerr << "Unable to read SDF string[" << _sdfString << "]\n";
    return false;
  }

  return this->LoadImpl(sdf->root);
}

/////////////////////////////////////////////////
bool FakeServer::LoadImpl(sdf::ElementPtr _elem,
                      const std::string &_physics)
{
  std::string host = "";
  unsigned int port = 0;

  gazebo::transport::get_master_uri(host, port);

  this->master = new gazebo::Master();
  this->master->Init(port);
  this->master->RunThread();

  // Load gazebo
  gazebo::load();

  physics::load();
  
  // If a physics engine is specified,
  if (_physics.length())
  {
    // Check if physics engine name is valid
    // This must be done after physics::load();
    if (!physics::PhysicsFactory::IsRegistered(_physics))
    {
      gzerr << "Unregistered physics engine [" << _physics
            << "], the default will be used instead.\n";
    }
    // Try inserting physics engine name if one is given
    else if (_elem->HasElement("world") &&
             _elem->GetElement("world")->HasElement("physics"))
    {
      _elem->GetElement("world")->GetElement("physics")
           ->GetAttribute("type")->Set(_physics);
    }
    else
    {
      gzerr << "Cannot set physics engine: <world> does not have <physics>\n";
    }
  }

  sdf::ElementPtr worldElem = _elem->GetElement("world");
  if (worldElem)
  {
    physics::WorldPtr world = physics::create_world();

    // Create the world
    try
    {
      physics::load_world(world, worldElem);
    }
    catch(common::Exception &e)
    {
      gzthrow("Failed to load the World\n"  << e);
    }
  }

  this->node = transport::NodePtr(new transport::Node());
  this->node->Init("/gazebo");
  this->serverSub = this->node->Subscribe("/gazebo/server/control",
                                          &FakeServer::OnControl, this);

  this->worldModPub =
    this->node->Advertise<msgs::WorldModify>("/gazebo/world/modify");

  // Run the gazebo, starts a new thread
  gazebo::run();

  return true;
}

/////////////////////////////////////////////////
void FakeServer::Init()
{

  gazebo::init();

  this->stop = false;
}

/////////////////////////////////////////////////
void FakeServer::SigInt(int)
{
  stop = true;
}

/////////////////////////////////////////////////
void FakeServer::Stop()
{
  this->stop = true;
}

/////////////////////////////////////////////////
void FakeServer::Fini()
{
  this->Stop();

  gazebo::fini();

  if (this->master)
    this->master->Fini();
  delete this->master;
  this->master = NULL;
}

/////////////////////////////////////////////////
void FakeServer::Run()
{
  if (this->stop)
    return;

  unsigned int iterations = 0;
  common::StrStr_M::iterator piter = this->params.find("iterations");
  if (piter != this->params.end())
  {
    try
    {
      iterations = boost::lexical_cast<unsigned int>(piter->second);
    }
    catch(...)
    {
      iterations = 0;
      gzerr << "Unable to cast iterations[" << piter->second << "] "
        << "to unsigned integer\n";
    }
  }

  // Run each world. Each world starts a new thread
  physics::run_worlds(iterations);

  
 // _publisher = new FakePublisher();
  // Update the sensors.
  while (!this->stop && physics::worlds_running())
  {
//	  _publisher->makeBox();
    this->ProcessControlMsgs();
    common::Time::MSleep(1);
  }

  // Stop all the worlds
  physics::stop_worlds();

  // Stop gazebo
  gazebo::stop();

  // Stop the master
  this->master->Stop();
}

/////////////////////////////////////////////////
void FakeServer::ProcessParams()
{
  common::StrStr_M::const_iterator iter;
  for (iter = this->params.begin(); iter != this->params.end(); ++iter)
  {
    if (iter->first == "pause")
    {
      bool p = false;
      try
      {
        p = boost::lexical_cast<bool>(iter->second);
      }
      catch(...)
      {
        // Unable to convert via lexical_cast, so try "true/false" string
        std::string str = iter->second;
        boost::to_lower(str);

        if (str == "true")
          p = true;
        else if (str == "false")
          p = false;
        else
          gzerr << "Invalid param value[" << iter->first << ":"
                << iter->second << "]\n";
      }

      physics::pause_worlds(p);
    }
    else if (iter->first == "record")
    {
    }
  }
}

/////////////////////////////////////////////////
void FakeServer::SetParams(const common::StrStr_M &_params)
{
  common::StrStr_M::const_iterator iter;
  for (iter = _params.begin(); iter != _params.end(); ++iter)
    this->params[iter->first] = iter->second;
}

/////////////////////////////////////////////////
void FakeServer::OnControl(ConstServerControlPtr &_msg)
{
  boost::mutex::scoped_lock lock(*this->receiveMutex);
  this->controlMsgs.push_back(*_msg);
}

/////////////////////////////////////////////////
void FakeServer::ProcessControlMsgs()
{
  std::list<msgs::ServerControl>::iterator iter;
  for (iter = this->controlMsgs.begin();
       iter != this->controlMsgs.end(); ++iter)
  {
    if ((*iter).has_save_world_name())
    {
      physics::WorldPtr world = physics::get_world((*iter).save_world_name());
      if ((*iter).has_save_filename())
        world->Save((*iter).save_filename());
      else
        gzerr << "No filename specified.\n";
    }
    else if ((*iter).has_new_world() && (*iter).new_world())
    {
      this->OpenWorld("worlds/empty.world");
    }
    else if ((*iter).has_open_filename())
    {
      this->OpenWorld((*iter).open_filename());
    }
  }
  this->controlMsgs.clear();
}

/////////////////////////////////////////////////
bool FakeServer::OpenWorld(const std::string & /*_filename*/)
{
  gzerr << "Open World is not implemented\n";
  return false;
}

