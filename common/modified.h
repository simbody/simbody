
#ifndef __modified_hh__
#define __modified_hh__



#include <cdsList.h>

class ModifiedBase;
class Simulation;

// Modified: class to help track the state of one object relative to another.
//
//
//
// classes which use a Modified member are generally derived from
// ModifiedBase
//
class Modified {
  bool value_;            //true if modified
  ModifiedBase* instance; 
public:
  Modified() : value_(0), instance(0) { }
  Modified(ModifiedBase* instance) : value_(0), instance(instance) {}
  void setInstance(ModifiedBase* instp) { instance=instp; }

  //set state to modified
  void set()   { value_=1; }    

  //set state to unmodified
  void clear() { value_=0; }    

  //call ModifiedBase::updateValues and clear modified.
  inline void update();

  // return the modified state
  bool value() const { return value_; }

  //a synonym for value()
  bool operator()() const { return value(); }
};

//
// ModifiedBase: virtual base class used to help track whether quantities 
//  are modified. 
//
// the updateValues method must be overridden in the derived class.
class ModifiedBase {
public:
  mutable Modified           modified;
  CDSList<const Simulation*> registeredSimulations;

  ModifiedBase() { modified.setInstance(this); }

  virtual ~ModifiedBase();

  // register the current object to a Simulation
  void registerTo(const Simulation*);
  // unregister the current object from a Simulation
  void unRegister(const Simulation*);

  // method to be called to get current object into a consistent state.
  virtual void updateValues() = 0;
};

void
Modified::update()
{ 
 if ( !value() ) 
   return;

 clear(); //here to avoid recursion
 instance->updateValues(); 
} /* update */

#endif /* __modified_hh__ */
