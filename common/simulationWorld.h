
#ifndef __simulationWorld_hh__
#define __simulationWorld_hh__

#include <randomNum.h>
#include <cdsList.h>
#include <enumNameMap.h>

// singleton class containing common infrastructure items

class RandomSeeder;

class SimulationWorld {

public:
  float_type timeFactor() const { return timeFactor_; }
  float_type kBoltzmann() const { return kBoltzmann_; }

  // global shared random number generator
  RandomNum random;
  void setRandomSeed(int seed);
  void registerRandomSeeder(RandomSeeder*);

  // log level
  enum  LogLevel { LOG_NONE=   0, 
		   LOG_NORMAL= 5,
		   LOG_ALL=    10,
		   LOG_DEBUG=  15 };
  
  LogLevel logLevel() const { return logLevel_;}
  void setLogLevel(const LogLevel& logLevel) {logLevel_ = logLevel;}
  
  

  static SimulationWorld* world() { return world_; }

  static void init(const float_type& timeFac   ,
		   const float_type& kBoltzmann,
		   const LogLevel&   logLevel);

private:
  SimulationWorld();

  float_type timeFactor_;
  float_type kBoltzmann_;
  LogLevel   logLevel_;

  CDSList<RandomSeeder*> randomSeederList;
  
  static SimulationWorld* world_;
};  

//
// callback class used so that a random seed is updated by
// SimulationWorld's setRandomSeed
//
class RandomSeeder {
public:
  RandomSeeder() {}
  virtual ~RandomSeeder() {}
  virtual void setSeed(int) = 0;
};

namespace EnumNamespace {
  extern EnumNameMap LogLevel[];
};
#endif /* __simulationWorld_hh__ */
