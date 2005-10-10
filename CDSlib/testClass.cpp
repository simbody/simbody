#define TESTING

#include "cdsString.h"
#include "cdsSStream.h"
#include "cdsList.h"
#include "fixedVector.h"
#include "fixedMatrix.h"
#include "fixedSymMatrix.h"
#include "subMatrix.h"
#include "cdsVector.h"
#include "cdsMatrix.h"
#include "symMatrix.h"
#include "cdsRegex.h"
#include "cdsString.h"
#include "cdsMap.h"
#include "cg.h"
#include "cdsAuto_arr.h"
#include "cdsAuto_ptr.h"
#include "cdsRc_ptr.h"
#include "cdsFdstream.h"
#include "matrixTools.h"
#include "circDLinkedList.h"
#include "chull.h"

int
main(int argc, char* argv[])
{
    int exit=0;

    try {
        exit |= CDSList_test();
        exit |= CDSString_test();
        exit |= CDSStringStream_test();
        exit |= FixedVector_test();
        exit |= FixedMatrix_test();
        exit |= FixedSymMatrix_test();
        exit |= CDSVector_test();
        exit |= CDSMatrix_test();
        exit |= SymMatrix_test();
        exit |= SubMatrix_test();
        //exit |= CDSRegex::test();
        //exit |= CDSMap<String,int>::test();
        //exit |= CDS::exception::test();
        exit |= auto_arr_test();
        //exit |= CDS::rc_ptr<int>::test();
        //exit |= CDS::auto_ptr<int>::test();
        //exit |= CDS::test_fdstreams();
        exit |= MatrixTools::matrixToolsTest();
        //exit |= CDS::CircDLinkedList<int>::test();
        //exit |= testCHull();
        //   CG::test();
    }

    // catch( FixedVectorRangeError ) {
    //   cout << "caught range error. Fatal error.\n";
    // }   
    // catch( IOError ) {
    //   cout << "caught io error. Fatal error.\n";
    // }   
    // catch ( CDS::exception e ) {
    //   cout << "caught exception: " << e.mess << '\n';
    // }

    catch (...) {
        cout << "caught exception. Fatal error.\n";
        exit = 1;
    }

    return exit;
}

