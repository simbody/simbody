
#include "simulationWorld.h"

#include <xplorWrap.h>

SimulationWorld* SimulationWorld::world_=0;

SimulationWorld::SimulationWorld()
{
}

class RandomClassSeeder : public RandomSeeder {
  RandomNum* randomp;
public:
  RandomClassSeeder(RandomNum& r) : randomp(&r) {}
  virtual ~RandomClassSeeder() {}
  virtual void setSeed(int seed) {randomp->setSeed(seed);}
};

void
SimulationWorld::init(const float_type& timeFac,
		      const float_type& kBoltzmann,
		      const LogLevel&   logLevel)
{
 if ( world() )
   return;

 world_ = new SimulationWorld();
 world_->timeFactor_ = timeFac;
 world_->kBoltzmann_ = kBoltzmann;
 world_->logLevel_   = logLevel;

 world()->registerRandomSeeder( new RandomClassSeeder(world()->random) );

} /* init */

void
SimulationWorld::setRandomSeed(int seed)
{
 for (int i=0 ; i<randomSeederList.size() ; i++)
   randomSeederList[i]->setSeed(seed);
} /* setRandomSeed */

void
SimulationWorld::registerRandomSeeder(RandomSeeder* seeder)
{
 randomSeederList.append( seeder );
} /* registerRandomSeeder */

namespace EnumNamespace {
  EnumNameMap LogLevel[] = {{"none",SimulationWorld::LOG_NONE},
			    {"normal",SimulationWorld::LOG_NORMAL},
			    {"all",SimulationWorld::LOG_ALL},
			    {"debug",SimulationWorld::LOG_DEBUG},
			    {"last enum",0}};
};
