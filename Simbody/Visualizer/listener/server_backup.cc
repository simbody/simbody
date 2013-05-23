
#include <gazebo/transport/transport.hh>
#include <gazebo/msgs/msgs.hh>
#include <gazebo/Master.hh>

#include <iostream>

class FakeServer
{
public:
	FakeServer() : simTime(0),
				pauseTime(0),
				realTime(0),
				paused(0),
				iterations(1111)
	{

		std::string host = "";
		unsigned int port = 0;

		gazebo::transport::get_master_uri(host, port);
		this->master = new gazebo::Master();
		this->master->Init(port);
		this->master->RunThread();

    gazebo::common::Console::Instance()->Init("fakeserver.log");
		gazebo::transport::init();
		node = gazebo::transport::NodePtr(new gazebo::transport::Node());
		node->Init("/opensim");
		// Advertise
		statPub = node->Advertise<gazebo::msgs::WorldStatistics>("~/world_stats");

		gazebo::transport::run();
	}

	void publishWorldStats()
	{
		// Set stuff in the protobuf
		gazebo::msgs::Set(worldStatsMsg.mutable_sim_time(), simTime);
		gazebo::msgs::Set(worldStatsMsg.mutable_real_time(), realTime);
		gazebo::msgs::Set(worldStatsMsg.mutable_pause_time(), pauseTime);

		worldStatsMsg.set_iterations(iterations);
		worldStatsMsg.set_paused(paused);

		// Publish
		statPub->Publish(this->worldStatsMsg);
	}

private:
	// Create transport node
	gazebo::transport::NodePtr node;
	// Create publisher for world_stats
	gazebo::transport::PublisherPtr statPub;
	// Create protobuf message
	gazebo::msgs::WorldStatistics worldStatsMsg;
	gazebo::Master * master;

	// Create fake world stat protobuf
	gazebo::common::Time simTime;
	gazebo::common::Time pauseTime;
	gazebo::common::Time realTime;
	bool paused;
	int iterations;
};

int main(int argc, char** argv)
{
	FakeServer s;

	while(true)
	{

	std::cout << "Server publishing messages" << std::endl;
		s.publishWorldStats();
	std::cout << "sleeping" << std::endl;
		usleep(10);
	}

}
