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
#ifndef _FAKESERVER_HH_
#define _FAKESERVER_HH_

#include <string>
#include <vector>
#include <list>

#include <boost/program_options.hpp>
#include <boost/thread.hpp>

#include "transport/TransportTypes.hh"
#include "common/CommonTypes.hh"
#include "physics/PhysicsTypes.hh"
#include "physics/World.hh"
#include "gazebo/Master.hh"
#include "sdf/sdf.hh"

namespace boost
{
  class mutex;
}

  /// \class Master Master.hh gazebo_core.hh
  /// \brief Base class for simulation server that handles commandline options,
  /// starts a Master, runs World update and sensor generation loops.
namespace SimTK { 
  class FakeServer
  {
    public: FakeServer();
    public: virtual ~FakeServer();

    public: void PrintUsage();
	public: bool ParseArgs(int argc = 0, char **argv = 0);

    /// \brief Load a world file and optionally override physics engine type.
    /// \param[in] _filename Name of the world file to load.
    /// \param[in] _physics Type of physics engine to use (ode|bullet).
    public: bool LoadFile(const std::string &_filename="worlds/empty.world",
                          const std::string &_physics="");

    public: bool LoadString(const std::string &_sdfString);
    public: void Init();
    public: void Run();
    public: void Stop();
    public: void Fini();

	public: void SetParams(const gazebo::common::StrStr_M &params);

    public: bool GetInitialized() const;

    /// \brief Load implementation.
    /// \param[in] _elem Description of the world to load.
    /// \param[in] _physics Type of physics engine to use (ode|bullet).
	private: bool LoadImpl(sdf::ElementPtr _elem,
                           const std::string &_physics="");

    private: static void SigInt(int _v);

    private: void ProcessParams();

	private: void OnControl(ConstServerControlPtr &_msg);

    private: bool OpenWorld(const std::string &_filename);

    private: void ProcessControlMsgs();

    public: static bool stop;

	private: gazebo::Master *master;
    private: boost::thread *masterThread;
	private: gazebo::transport::NodePtr node;
	private: gazebo::transport::SubscriberPtr serverSub;
	private: gazebo::transport::PublisherPtr worldModPub;

    private: boost::mutex *receiveMutex;
	private: std::list<gazebo::msgs::ServerControl> controlMsgs;

    private: gazebo::common::StrStr_M params;

    // save argc and argv for access by system plugins
    public: int systemPluginsArgc;
    public: char** systemPluginsArgv;
  };
}
#endif
