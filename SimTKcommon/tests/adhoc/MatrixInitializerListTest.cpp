#include "SimTKcommon.h"

#include <iostream>
using std::cout;
using std::endl;

using namespace SimTK;

int main(int argc, char **argv) {
 
  cout << "Initializing Matrix_ with brace-init-list {{1.5,24.0,5.5,6.5},{3,4,7,8},{4.1,5.2 ,6.3}}" << endl;
  Matrix_<Real> matrix = {{1.5,24.0,5.5,6.5},
                          {3  ,4   ,7  ,8  },
		          {4.1,5.2 ,6.3    }};
  
  cout << "Matrix *** " << endl;
  cout << matrix.toString() << endl;
  cout << "{{***}}" << endl;
}