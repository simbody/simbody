

#include "SimTKcommon.h"

#include <iostream>
using std::cout;
using std::endl;

using namespace SimTK;

int main(int argc, char **argv) {
 
  cout << "Initializing Vector_ with braced-init-list {1.1,2.2,3.3,4.4,5.5,4.3,33,45.4}" << endl;
  Vector_<Real> testvec = {1.1,2.2,3.3,4.4,5.5,4.3,33,45.4};
  
  int i=0;
  cout << "idx : val" << endl;
  for(auto const& n : testvec){
    cout << i << " : " << n << endl;
    i++;
  }
}
