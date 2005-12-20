
#ifndef __atomSelAction_hh__
#define __atomSelAction_hh__

class Simulation;

#include "cdsVec3.h"
#include "cdsMat33.h"
#include <cdsVector.h>
#include <cdsMap.h>
#include <atomSel.h>


namespace AtomSelAction {

//
// base class: all atomSelActions inherit from this. Overide the
// run class to perform a useful task on an atom
//

  class Base { 
  public:
    virtual ~Base() {}

    virtual void run(Simulation* sim  ,
		     int         index) {}

    virtual void init(const AtomSel& sel) {}

    virtual void finish() {}
    
  };

  class SetVelAction : public Base {
    CDSVec3 vel;
  public:
    SetVelAction(const CDSVec3& vel) : vel(vel) {}
    virtual void run(Simulation* sim,
		     int         index);
  };

  class Translate : public Base {
    CDSVec3 trans;
  public:
    Translate(const CDSVec3& trans) : trans(trans) {}
    virtual void run(Simulation* sim,
		     int         index);
  };

  class Rotate : public Base {
    CDSMat33 rot;
    CDSVec3  center;
  public:
    Rotate(const CDSMat33& rot,
	   const CDSVec3 center=CDSVec3(0,0,0)) : 
      rot(rot), center(center) {}
    virtual void run(Simulation* sim,
		     int         index);
  };

  class Fit : public Base {
    CDSMat33 rot;
    CDSVec3  trans;
    
    CDSVector<CDSVec3> fitTo;
    AtomSel       fitBy;
    CDSVector<CDSVec3>* altCoords;
  public:
    Fit(const CDSVector<CDSVec3>& fitTo,
	const AtomSel&       fitBy,
	      CDSVector<CDSVec3>* altCoords=0) : 
      rot(1,0,0,
	  0,1,0,
	  0,0,1), 
      trans(0,0,0), 
      fitTo(fitTo), 
      fitBy(fitBy),
      altCoords(altCoords) {}

    virtual void init(const AtomSel& sel);

    virtual void run(Simulation* sim,
		     int         index);

    CDSMat33 rotation()    { return rot; }
    CDSVec3  translation() { return trans; }
  };

  class RMSD : public Base {

    int        size;
    float_type sum;
    CDSMap<int, int> count;
    CDSMap<int, float_type> byResidue_;

    CDSVector<CDSVec3> compare;
  public:
    RMSD(const CDSVector<CDSVec3>& compare) : 
      size(0), sum(0.), compare(compare) {}

    virtual void init(const AtomSel& sel);
    virtual void run(Simulation* sim,
		     int         index);

    virtual void finish();
    
    float_type rmsd() const { return sum; }

    //rmsd for each residue
    CDSMap<int, float_type> byResidue() const { return byResidue_; }

  };
    

}; // end of namespace AtomSelAction declarations


#endif /* __atomSelAction_hh__ */
