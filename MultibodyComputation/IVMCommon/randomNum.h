#ifndef __randomNum_hh__
#define __randomNum_hh__

//
// a random number generator class,
// based on the IMSL call GGUBFS
//

class RandomNum {

 public:

  // constructor

  RandomNum();

  // accessors

  void setSeed(int newVal);
  int  seed() const;

  // uniform random number in the range 0..1
  double uniform();

 private:

  double seed_;

};

#endif /* __randomNum_hh__ */
