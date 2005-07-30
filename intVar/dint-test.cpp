/**@file
 *
 * Test internal dynamics routines.
 *
 * CDS - 3/30/00
 */

#include "dinternal.h"
#include "dint-step.h"
#include "dint-loop.h"
#include "dint-atom.h"

#include <cdsAuto_arr.h>
#include <cdsAuto_ptr.h>

#include <sthead.h>
#include <cdsList.h>
#include <cdsVector.h>
#include <cdsString.h>
#include <fixedVector.h>
#include <cdsIomanip.h>
#include <math.h>

static int fcnt = 0;

/**
 * A concrete class which implements an "IVM" which is the internal
 * variable mechanical model for the molecule system. A concrete IVM
 * must implement the calcEnergy() method which calculates both energy
 * and forces. This testing class just hangs a spring between each pair
 * of atoms.
 */
class TestIVM : public IVM {
public:
    TestIVM() : IVM() {}
    virtual void calcEnergy();

    void init();
    void initCycle();

    int test();
};

void
TestIVM::calcEnergy()
{
    double k = 1.0;
    fcnt++;

    Epotential_ = 0.0;
    //FIX: check derivatives
    for (int i=1 ; i<atoms.size() ; i++) 
        atoms[i]->deriv.set( 0.0 );
    for (l_int i=1 ; i<atoms.size() ; i++) {
        for (int j=i+1 ; j<atoms.size() ; j++) {
            double dist2 = abs2( atoms[i]->pos - atoms[j]->pos );
            Epotential_ += k * dist2;
            atoms[i]->deriv += k * (atoms[i]->pos - atoms[j]->pos);
            atoms[j]->deriv += k * (atoms[j]->pos - atoms[i]->pos);
            //     atoms[i]->deriv *= -1;
        }
    }
    Epotential_ *= 0.5;
    calcTemperature(); //calcs Ekinetic
    Etotal_ = Ekinetic_ + Epotential_;
}


//
// initialize InternalDynamics quantities from xplor values
//
// 12 -- 34
// 5
// 67
//
void 
TestIVM::init()
{
    int natom = 7;
    double masses[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
    double kBoltzman = 1.0;

    initAtoms(natom,
              masses,
              kBoltzman); //create atoms array

    atoms[1]->bonds.append( atoms[2] );
    atoms[2]->bonds.append( atoms[1] );
    atoms[2]->bonds.append( atoms[3] );
    atoms[3]->bonds.append( atoms[2] );
    atoms[3]->bonds.append( atoms[4] );
    atoms[4]->bonds.append( atoms[3] );

    atoms[6]->bonds.append( atoms[7] );
    atoms[7]->bonds.append( atoms[6] );

    atoms[1]->pos = Vec3(22.304,  25.009,  14.716);
    atoms[2]->pos = Vec3(22.861,  24.772,  13.871);
    atoms[3]->pos = Vec3(22.882,  25.590,  13.228);
    atoms[4]->pos = Vec3(22.240,  23.610,  13.175);
    atoms[5]->pos = Vec3(23.831,  24.530,  14.156);
    atoms[6]->pos = Vec3(22.441,  22.352,  14.024);
    atoms[7]->pos = Vec3(23.241,  22.330,  14.937);

    atoms[1]->vel = Vec3(1.0 , 0, 0);
    atoms[2]->vel = Vec3(0.0 , 0, 0);
    atoms[3]->vel = Vec3(0.0 , 0, 0);
    atoms[4]->vel = Vec3(0.0 , 0, 0);
    atoms[5]->vel = Vec3(0.0 , 0, 0);
    atoms[6]->vel = Vec3(0.0 , 0, 0);
    atoms[7]->vel = Vec3(0.0 , 0, 0);

    atoms[1]->deriv = Vec3(0.0 , 0, 0);
    atoms[2]->deriv = Vec3(0.0 , 0, 0);
    atoms[3]->deriv = Vec3(0.0 , 0, 0);
    atoms[4]->deriv = Vec3(0.0 , 0, 0);
    atoms[5]->deriv = Vec3(0.0 , 0, 0);
    atoms[6]->deriv = Vec3(0.0 , 0, 0);
    atoms[7]->deriv = Vec3(0.0 , 0, 0);

    groupList.resize(0);
    hingeList.resize(0);
    oldBaseAtoms.resize(0);
    // InternalDynamics::initTree();
    initDynamics(0);

    groupTorsion();
    InternalDynamics::HingeSpec hingeSpec("torsion");
    for (int i=0 ; i<=natom ; i++)
        hingeSpec.aList.append( i );
    hingeList.append( hingeSpec );
    // verbose |= printNodeDef;
    // InternalDynamics::initTree();
    // InternalDynamics::initDynamics();
}

//
// initialize InternalDynamics quantities for ring
//
//       2---15
//    /     
// 1-9   3--6   14-13-12
//   |   |             |
//   5---4---8--7--10-11
//
// constraints 6<->15   10<->13 7<->14
//
void 
TestIVM::initCycle()
{
    int natom = 15;
    double masses[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                        10.,11.,12.,13.,14.,15.};
    double kBoltzman = 1.0;

    initAtoms(natom,
    masses,
    kBoltzman); //create atoms array

    atoms[1]->bonds.append( atoms[9] );
    atoms[9]->bonds.append( atoms[1] );
    atoms[9]->bonds.append( atoms[2] );
    atoms[2]->bonds.append( atoms[9] );
    atoms[3]->bonds.append( atoms[4] );
    atoms[4]->bonds.append( atoms[3] );
    atoms[4]->bonds.append( atoms[5] );
    atoms[5]->bonds.append( atoms[4] );
    atoms[5]->bonds.append( atoms[9] );
    atoms[9]->bonds.append( atoms[5] );
    atoms[3]->bonds.append( atoms[6] );
    atoms[6]->bonds.append( atoms[3] );
    atoms[7]->bonds.append( atoms[8] );
    atoms[8]->bonds.append( atoms[7] );
    atoms[8]->bonds.append( atoms[4] );
    atoms[4]->bonds.append( atoms[8] );
    atoms[7]->bonds.append( atoms[10] );
    atoms[10]->bonds.append( atoms[7] );
    atoms[10]->bonds.append( atoms[11] );
    atoms[11]->bonds.append( atoms[10] );
    atoms[11]->bonds.append( atoms[12] );
    atoms[12]->bonds.append( atoms[11] );
    atoms[12]->bonds.append( atoms[13] );
    atoms[13]->bonds.append( atoms[12] );
    atoms[2]->bonds.append( atoms[15] );
    atoms[15]->bonds.append( atoms[2] );
    atoms[13]->bonds.append( atoms[14] );
    atoms[14]->bonds.append( atoms[13] );

    atoms[1 ]->pos = Vec3(51.407,    48.981,     12.253);         
    atoms[2 ]->pos = Vec3(51.486,    47.935,     13.271);         
    atoms[3 ]->pos = Vec3(50.940,    47.075,     12.913);         
    atoms[4 ]->pos = Vec3(52.935,    47.479,     13.621);         
    atoms[5 ]->pos = Vec3(53.614,    47.873,     12.880);         
    atoms[6 ]->pos = Vec3(53.055,    46.062,     13.646);         
    atoms[7 ]->pos = Vec3(53.397,    45.808,     14.505);         
    atoms[8 ]->pos = Vec3(53.354,    48.022,     14.977);         
    atoms[9 ]->pos = Vec3(53.150,    49.082,     15.018);         
    atoms[10]->pos = Vec3(54.410,    47.853,     15.122);         
    atoms[11]->pos = Vec3(52.799,    47.517,     15.753);         
    atoms[12]->pos = Vec3(50.707,    48.281,     14.531);         
    atoms[13]->pos = Vec3(50.914,    49.332,     15.148);         
    atoms[14]->pos = Vec3(52.914,    50.332,     15.148);         
    atoms[15]->pos = Vec3(50.914,    48.332,     14.148);         

    atoms[1 ]->vel = Vec3(20.0, 0, 0);
    atoms[2 ]->vel = Vec3(0.0 , 0, 0);
    atoms[3 ]->vel = Vec3(0.0 , 0, 0);
    atoms[4 ]->vel = Vec3(0.0 , 0, 0);
    atoms[5 ]->vel = Vec3(0.0 , 0, 0);
    atoms[6 ]->vel = Vec3(0.0 , 0, 0);
    atoms[7 ]->vel = Vec3(0.0 , 0, 0);
    atoms[8 ]->vel = Vec3(0.0 , 0, 0);
    atoms[9 ]->vel = Vec3(0.0 , 0, 0);
    atoms[10]->vel = Vec3(0.0 , 0, 0);
    atoms[11]->vel = Vec3(0.0 , 0, 0);
    atoms[12]->vel = Vec3(0.0 , 0, 0);
    atoms[13]->vel = Vec3(0.0 , 0, 0);
    atoms[14]->vel = Vec3(0.0 , 0, 0);
    atoms[15]->vel = Vec3(0.0 , 0, 0);

    atoms[1 ]->deriv = Vec3(0.0 , 0, 0);
    atoms[2 ]->deriv = Vec3(0.0 , 0, 0);
    atoms[3 ]->deriv = Vec3(0.0 , 0, 0);
    atoms[4 ]->deriv = Vec3(0.0 , 0, 0);
    atoms[5 ]->deriv = Vec3(0.0 , 0, 0);
    atoms[6 ]->deriv = Vec3(0.0 , 0, 0);
    atoms[7 ]->deriv = Vec3(0.0 , 0, 0);
    atoms[8 ]->deriv = Vec3(0.0 , 0, 0);
    atoms[9 ]->deriv = Vec3(0.0 , 0, 0);
    atoms[10]->deriv = Vec3(0.0 , 0, 0);
    atoms[11]->deriv = Vec3(0.0 , 0, 0);
    atoms[12]->deriv = Vec3(0.0 , 0, 0);
    atoms[13]->deriv = Vec3(0.0 , 0, 0);
    atoms[14]->deriv = Vec3(0.0 , 0, 0);
    atoms[15]->deriv = Vec3(0.0 , 0, 0);

    // break cycles and add constraints
    constraintList.resize(0);
    constraintList.append( Pair(6,15) );
    constraintList.append( Pair(10,13) );
    constraintList.append( Pair(7,14) );

    groupList.resize(0);
    hingeList.resize(0);
    oldBaseAtoms.resize(0);
    // InternalDynamics::initTree();
    initDynamics(0);

    groupTorsion();
    InternalDynamics::HingeSpec hingeSpec("torsion");
    for (int i=0 ; i<=natom ; i++)
        hingeSpec.aList.append( i );
    hingeList.append( hingeSpec );

    // verbose |= printNodeDef;
    // InternalDynamics::initTree();
    // InternalDynamics::initDynamics();
}


int 
TestIVM::test()
{
    int exit = 0;
    double tol = 1e-14;
    // using namespace InternalDynamics;
    setVerbose(InternalDynamics::printNodeDef|InternalDynamics::printLoopInfo); // XXX
    init();
    initDynamics(0);
    const AtomList& atoms = getAtoms();

    {
    cout << "initial velocities...";
    int lexit = 0;
    CDSList<Vec3> t;
    t.append(Vec3(0,0,0));
    t.append(Vec3( .7907862337485985,  .09013290236325137,  .3165658097354623));
    t.append(Vec3( .1357399466529780, -.04150486007265110, -.0783012106372617));
    t.append(Vec3(-.0115291096445745,  .00003488899585902, -.0302656514525142));
    t.append(Vec3(-.0069196995302082, -.0018069623013814 , -.0172916085258486));
    t.append(Vec3(0,0,0));
    t.append(Vec3(0,0,0));
    t.append(Vec3(0,0,0));
    t.append(Vec3(-.1947872267530600, .09268407916055810, .1082815629728186));

    for (int i=0 ; i<atoms.size() ; i++)
        if ( abs2(atoms[i]->vel - t[i]) > tol*tol ) {
            cout << "FAILED";
            cout << atoms[i]->vel << " != " << t[i] << '\n';
            lexit=1;
        }
    if ( lexit==0 ) cout << "ok";
    else exit = 1;
    }
    cout << endl;

    {
    cout << "resetCM...";  cout.flush();
    int lexit = 0;
    //verbose |= printResetCM;
    resetCMflag=1;
    resetCM();
    CDSList<Vec3> t;
    t.append(Vec3(0,0,0));
    t.append(Vec3( .6523185106045116 , .06017487086985825, .2844882418348770 ));
    t.append(Vec3( .0508587296756875 ,-.04505921778863320,-.08246160309446333));
    t.append(Vec3(-.1036573357528213 ,-.00280589578215078,-.03375494675339982));
    t.append(Vec3(-.00115774397098383,-.03555440171077225,-.05191902572074229));
    t.append(Vec3(-.08746869633676819, .04323713586703393, .04437784061572216));
    t.append(Vec3( .02388992715731684,-.02368709773860604,-.02401235532135897));
    t.append(Vec3(-.02063272401050720, .01521653799906133, .01593726233525430));
    t.append(Vec3(-.2509439280717113 , .07462619327086268, .08934397222309618));
      
    for (int i=0 ; i<atoms.size() ; i++)
        if ( abs2(atoms[i]->vel - t[i]) > tol*tol ) {
            cout << "FAILED";
            cout << i << ": " << atoms[i]->vel << " != " << t[i] << '\n';
            lexit=1;
        }
    if ( lexit==0 ) cout << "ok";
    else exit = 1;
    } 
    cout << endl;

    using namespace InternalDynamics;
    {
    cout << "step...";cout.flush();
    bool debugging=0;
    if (debugging) 
        setVerbose( verbose() | 
                    printCoords | printEnergy | printNodePos | 
                    printNodeTheta | printNodeDef | printStepDebug | 
                    printStepInfo | printLoopInfo | printLoopDebug );
    else 
        setVerbose( 0 );

    int lexit = 0;
    double stepsize=1e-6;
    initDynamics(0);
    double E0 = Etotal();

    step(stepsize);

    CDSList<Vec3> tq;
    tq.append(Vec3(0,0,0));

    tq.append(Vec3(22.30400065231726, 25.00900006017106, 14.71600028448686));
    tq.append(Vec3(22.86100005086041, 24.77199995493966, 13.87099991753850));
    tq.append(Vec3(22.88199989634295, 25.58999999719345, 13.22799996624573));
    tq.append(Vec3(22.23999999884240, 23.60999996444499, 13.17499994808156));
    tq.append(Vec3(23.83099991253047, 24.53000004323679, 14.15600004437767));
    tq.append(Vec3(22.44100002388992, 22.35199997631407, 14.02399997598730));
    tq.append(Vec3(23.24099997936724, 22.33000001521754, 14.93700001593694));
    tq.append(Vec3(22.56319974905955, 24.57630007462422, 13.48420008934329));

    CDSList<Vec3> tv;
    tv.append(Vec3(0,0,0));
    tv.append(Vec3(    .6523192444544815 , .06016935435394602, .2844877178141606 ));
    tv.append(Vec3(    .05086118154323828,-.04506153298146813,-.08246129196430015));
    tv.append(Vec3(-.1036569920762576 ,-.00280736791994937,-.03375358837610575));
    tv.append(Vec3(-.00115761337467347,-.03555578463202729,-.05191819456361482));
    tv.append(Vec3(-.08747009973676708, .04323643246703309, .04437764361572128));
    tv.append(Vec3(    .02388993712409297,-.02368513074477553,-.02401287014865009));
    tv.append(Vec3(-.02063275741060132, .01521822514720666, .01593677318721820));
    tv.append(Vec3(-.2509401063298388 , .07462330037371790, .08934386979265429));
      
    double tol = 1e-12;
    for (int i=0 ; i<atoms.size() ; i++) {
        if ( abs2(atoms[i]->pos - tq[i])
             > 0.5*tol*tol*abs2(atoms[i]->pos + tq[i]) )
        {
            cout << "FAILED\n";
            cout << "pos(" << i << ") : " 
                 << atoms[i]->pos << " != " << tq[i] << '\n';
            lexit=1;
        }
        tol = 1e-8; //velocity determined to lower order
        if ( abs2(atoms[i]->vel - tv[i])
             > 0.5*tol*tol*abs2(atoms[i]->vel + tv[i]) )
        {
            cout << "FAILED\n";
            cout << "vel(" << i << ") : " 
                 << atoms[i]->vel << " != " << tv[i] << '\n';
            lexit=1;
        }
    }
    step(stepsize); step(stepsize); step(stepsize);

    stepsize=1e-4;
    for (l_int i=0 ; i<100 ; i++) {
        if (debugging)
            cout << i << ' ' << Etotal() << '\n';
        if (debugging)
            for (int i=1 ; i<atoms.size() ; i++)
                cout << '\t' << atoms[i]->pos << '\n';
        try {
            step(stepsize);
            printCM();
        }
        catch ( Solver::Finished ) {
            break;
        }
    }
    if (fabs(E0 - Etotal()) > 1e-3) {
        cout << "initial/final energy: " 
             << E0 << '/' << Etotal() << '\n';
        lexit=1;
    }
    if ( lexit==0 )
        cout << "ok";
    else {
        exit = 1;
        cout << "FAILED";
    }
    }
    cout << endl;

    // {
    //   cout << "minimizcg...";
    //   bool debugging=0;
    //   if (debugging) 
    //     verbose |= printCoords | printEnergy | printNodePos | printNodeTheta |
    //            printNodeDef | printStepInfo | printStepDebug;
    //   else
    //     verbose = 0;
    //   fcnt = 0;
    //   init();
    //   for (int ecnt=0 ; ecnt<1 ; ecnt++) {
    //     InternalDynamics::dEpred=5;
    //     //     InternalDynamics::Etolerance=1e-2;
    //     InternalDynamics::Gtolerance=0.01;
    //     integrateType = "MinimizeCG";
    ////     integrateType = "Minimize";
    ////     integrateType = "ConMin";
    //     InternalDynamics::initDynamics(0);
    //     if (debugging)
    //       InternalDynamics::printCM();
    //     for (int i=0 ; i<100 ; i++) {
    //       double gg = 0;
    //       for (int j=1 ; j<=7 ; j++)
    //         gg += abs2(atoms[j]->deriv);
    //       if (debugging)
    //         cout << i << ' ' << InternalDynamics::Epotential << ' ' << gg << '\n';
    //       if (debugging)
    //         for (int i=1 ; i<=7 ; i++)
    //           cout << '\t' << atoms[i]->pos << '\n';
    //       try {
    //         double s = 1.0;
    //         step(s);
    //         InternalDynamics::printCM();
    //       }
    //       catch ( Solver::Finished ) {
    //         break;
    //       }
    //     }
    //   }
    //   double gg = 0;
    //   for (int i=1 ; i<=7 ; i++)
    //     gg += abs2(atoms[i]->deriv);
    //   if ( fabs(Epotential-12.91) < 0.01 ) {
    //     cout << "ok";
    //   } else {
    //     exit = 1;
    //     cout << "FAILED";
    //   }
    //   
    //   if (debugging) {
    //     cout << "final: " << ' ' << Epotential << ' ' << gg << '\n';
    //     for (l_int i=1 ; i<=7 ; i++)
    //       cout << '\t' << atoms[i]->pos << '\n';
    //   }
    //   cout << "  [function evaluations: " << fcnt << "]\n";
    // }

    {
    cout << "powell...";cout.flush();
    bool debugging=0;
    if (debugging) 
        setVerbose( verbose() | printCoords | printEnergy | 
                    printNodePos | printNodeTheta | printNodeDef | 
                    printStepInfo | printStepDebug );
    else
        setVerbose( 0 );
    fcnt = 0;
    init();
    for (int ecnt=0 ; ecnt<1 ; ecnt++) {
        dEpred_=5;
        //     InternalDynamics::Etolerance=1e-2;
        Gtolerance_=0.01;
        solverType = "Powell";
        //     integrateType = "Minimize";
        //     integrateType = "ConMin";
        initDynamics(0);
        if (debugging)
            printCM();
        for (int i=0 ; i<100 ; i++) {
            double gg = 0;
            for (int j=1 ; j<=7 ; j++)
                gg += abs2(atoms[j]->deriv);
                if (debugging)
                    cout << i << ' ' << Epotential() << ' ' << gg << '\n';
                if (debugging)
                    for (int i=1 ; i<=7 ; i++)
                        cout << '\t' << atoms[i]->pos << '\n';
                try {
                    double s = 1.0;
                    step(s);
                    printCM();
                }
                catch ( Solver::Finished ) {
                    break;
                }
        }
    }
    double gg = 0;
    for (int i=1 ; i<=7 ; i++)
        gg += abs2(atoms[i]->deriv);
    if ( fabs(Epotential()-12.91) < 0.01 ) {
        cout << "ok";
    } else {
        exit = 1;
        cout << "FAILED";
    }
      
    if (debugging) {
        cout << "final: " << ' ' << Epotential() << ' ' << gg << '\n';
        for (l_int i=1 ; i<=7 ; i++)
            cout << '\t' << atoms[i]->pos << '\n';
    }
    cout << "  [function evaluations: " << fcnt << "]";
    }
    cout << endl;

    {
    cout << "conmin...";
    bool debugging=0;
    if (debugging) 
        setVerbose( verbose() | printCoords | printEnergy | 
                    printNodePos | printNodeTheta );
    else 
        setVerbose( 0 );
    fcnt = 0;
    init();
    for (int ecnt=0 ; ecnt<1 ; ecnt++) {
        dEpred_=5;
        //     InternalDynamics::Etolerance=1e-2;
        Gtolerance_=0.01;
        solverType = "Powell";
        solverType = "Minimize";
        solverType = "ConMin";
        initDynamics(0);
        if (debugging)
            printCM();
        for (int i=0 ; i<100 ; i++) {
            double gg = 0;
            for (int j=1 ; j<=7 ; j++)
                gg += abs2(atoms[j]->deriv);
            if (debugging)
                cout << i << ' ' << Epotential() << ' ' << gg << '\n';
            if (debugging)
                for (int i=1 ; i<=7 ; i++)
                    cout << '\t' << atoms[i]->pos << '\n';
            try {
                double s = 1.0;
                step(s);
                printCM();
            }
            catch ( Solver::Finished ) {
                break;
            }
        }
    }
    double gg = 0;
    for (int i=1 ; i<=7 ; i++)
        gg += abs2(atoms[i]->deriv);
    if ( fabs(Epotential()-12.91) < 0.01 ) {
        cout << "ok";
    } else {
        exit = 1;
        cout << "FAILED";
    }
      
    if (debugging) {
        cout << "final: " << ' ' << Epotential() << ' ' << gg << '\n';
        for (l_int i=1 ; i<=7 ; i++)
            cout << '\t' << atoms[i]->pos << '\n';
    }
    cout << "  [function evaluations: " << fcnt << "]\n";
    }

    {
    cout << "simplex...";
    bool debugging=0;
    if (debugging) 
        setVerbose( verbose() | printCoords | printEnergy | 
                    printNodePos | printNodeTheta );
    else 
        setVerbose( 0 );
    fcnt = 0;
    init();
    solverType = "Simplex";
    for (int ecnt=0 ; ecnt<5 ; ecnt++) {
        initDynamics(0);
        dEpred_=1;
        //double stepsize=0.1;
        Etolerance_=1e-8;
        if (debugging)
            printCM();
        for (int i=0 ; i<100 ; i++) {
            double gg = 0;
            for (int j=1 ; j<=7 ; j++)
                gg += abs2(atoms[j]->deriv);
            if (debugging)
                cout << i << ' ' << Epotential() << ' ' << gg << '\n';
            if (debugging)
                for (int i=1 ; i<=7 ; i++)
                    cout << '\t' << atoms[i]->pos << '\n';
            try {
                double s = 1.0;
                step(s);
                printCM();
            }
            catch ( Solver::Finished ) {
                break;
            }
        }
    }
    double gg = 0;
    for (int i=1 ; i<=7 ; i++)
        gg += abs2(atoms[i]->deriv);
    if ( fabs(Epotential()-12.91) < 0.01 ) {
        cout << "ok";
    } else {
        exit = 1;
        cout << "FAILED";
    }
      
    if (debugging) {
        cout << "final: " << ' ' << Epotential() << ' ' << gg << '\n';
        for (l_int i=1 ; i<=7 ; i++)
            cout << '\t' << atoms[i]->pos << '\n';
    }
    cout << "  [function evaluations: " << fcnt << "]\n";
    }

    {
    cout << "length constraints...";
    bool debugging=0;
    if (debugging) 
        verbose_ |= printCoords | printEnergy | printNodePos | printNodeTheta |
        printNodeDef | printStepDebug | printStepInfo |
        printLoopInfo | printLoopDebug;
    else 
        verbose_ = 0;
    fcnt = 0;
    //double stepsize=0.1;
    Etolerance_=1e-8;
    Gtolerance_=0.01;
    //FIX:   LengthConstraints::maxIters = 40;
    solverType = "Powell";
    initCycle();
    CDSList<double> bondLengths;
    for (int i=0 ; i<constraintList.size() ; i++) {
        IVMAtom* a1 = atoms[ constraintList[i].a ];
        IVMAtom* a2 = atoms[ constraintList[i].b ];
        bondLengths.append( norm(a1->pos - a2->pos) );
    }
    bool ok=1;
    for (int ecnt=0 ; ecnt<1 ; ecnt++) {
        initDynamics(0);
        if (debugging)
            printCM();
        for (int i=0 ; i<100 ; i++) {
            double gg = 0;
            for (int j=1 ; j<atoms.size() ; j++)
                gg += abs2(atoms[j]->deriv);
            if (debugging)
                cout << i << ' ' << Epotential() << ' ' << gg << '\n';
            if (debugging)
                for (int i=1 ; i<atoms.size() ; i++)
                    cout << '\t' << atoms[i]->pos << '\n';
            try {
                double s = 1.0;
                step(s);
                printCM();
            }
            catch ( Solver::Finished f ) {
                ok = f.ok;
                if ( !ok ) 
                    cout << "minimizer exited with error status: " << ok << '\n';
                break;
            }
        }
    }
    for (l_int i=0 ; i<constraintList.size() ; i++) {
        IVMAtom* a1 = atoms[ constraintList[i].a ];
        IVMAtom* a2 = atoms[ constraintList[i].b ];
        if ( fabs(bondLengths[i] - norm(a1->pos - a2->pos))
             > 1e-7*bondLengths[i] )
        {
            cout << "bond between atoms " << a1->index << " and "
                 << a2->index << " is: "
                 << setprecision(8)
                 << norm(a1->pos - a2->pos) 
                 << "  was: " << bondLengths[i] << '\n';
            ok=0;
        }
    }
    if (ok)
        cout << "ok";
    else {
        exit = 1;
        cout << "FAILED";
    }
      
    if (debugging) {
        double gg = 0;
        for (int j=1 ; j<atoms.size() ; j++)
            gg += abs2(atoms[j]->deriv);
        cout << "final: " << ' ' << Epotential() << ' ' << gg << '\n';
        for (l_int i=1 ; i<atoms.size() ; i++)
            cout << '\t' << atoms[i]->pos << '\n';
    }
    cout << "  [function evaluations: " << fcnt << "]\n";
    }

    {
    cout.flush();
    cout << "length constraints- integrated...";
    bool debugging=0;
    if (debugging) 
        verbose_ = printCoords | printEnergy | printNodePos | printNodeTheta |
                   printNodeDef | printStepDebug | printStepInfo |
                   printLoopInfo | printLoopDebug;
    else 
        verbose_ = 0;
    fcnt = 0;
    double stepsize=1e-3;
    Etolerance_=1e-8;
    Gtolerance_=0.01;
    //FIX: LengthConstraints::maxIters = 40;
    solverType = "PC6";
    initCycle();
    CDSList<double> bondLengths;
    for (int i=0 ; i<constraintList.size() ; i++) {
        IVMAtom* a1 = atoms[ constraintList[i].a ];
        IVMAtom* a2 = atoms[ constraintList[i].b ];
        bondLengths.append( norm(a1->pos - a2->pos) );
    }
    double E0;
    for (int ecnt=0 ; ecnt<1 ; ecnt++) {
        initDynamics(0);
        calcEnergy();
        E0 = Etotal();
        if (debugging)
            printCM();
        for (int i=0 ; i<100 ; i++) {
            if (debugging)
                cout << i << ' ' << Etotal() << '\n';
            if (debugging)
                for (int i=1 ; i<atoms.size() ; i++)
                    cout << '\t' << atoms[i]->pos << '\n';
            try {
                step(stepsize);
                printCM();
            }
            catch ( Solver::Finished ) {
                break;
            }
        }
    }
    bool ok=1;
    for (l_int i=0 ; i<constraintList.size() ; i++) {
        IVMAtom* a1 = atoms[ constraintList[i].a ];
        IVMAtom* a2 = atoms[ constraintList[i].b ];
        if ( fabs(bondLengths[i] - norm(a1->pos - a2->pos)) > 1e-7 ) {
            cout << "bond between atoms " << a1->index << " and "
                 << a2->index << " is: "
                 << norm(a1->pos - a2->pos) 
                 << "  was: " << bondLengths[i] << '\n';
            ok=0;
        }
    }
    if (fabs(E0 - Etotal()) > 1e-2) {
        cout << "initial/final energy: " 
             << E0 << '/' << Etotal() << '\n';
        ok=0;
    }
    if (ok)
        cout << "ok\n";
    else {
        exit = 1;
        cout << "FAILED\n";
    }
      
      
    if (debugging) {
        calcEnergy();
        cout << "final: " << ' ' << Etotal() << '\n';
        for (l_int i=1 ; i<atoms.size() ; i++)
            cout << '\t' << atoms[i]->pos << '\n';
    }
    }

    return exit;
}


int
main(int argc, char* argv[])
{
    CDS::auto_ptr<TestIVM> ivm(new TestIVM);
    return  ivm->test();
}

// define null xplor symbols
extern "C" {
void FORTRAN(energy,ENERGY)(){}
void FORTRAN(assfil,ASSFIL)(const char* name,
                  int&  coordUnit,
            const char* mode,
            const char* form,
                  int&  err,
            const int   lname,
            const int   lmode,
            const int   lform){}
void FORTRAN(writtc,WRITTC)(const char* str,
            const int& natom,
            const double* x,
            const double* y,
            const double* z,
            const int& crdnum,
            const int* crdind,
            const int& stepnum,
            const int& qfirsc,
            const double& timeStep,
            const int& nsavc,
            const int& coordUnit,
            const int &formattedOutput,
            const char* format,
            const double& scale,
            const double& offset,
            const int     lstr,
            const int     lformat){}
void FORTRAN(declar,DECLAR)(const char*   name,
            const char*   type,
            const char*   stinit,
            const char*   dont_use_me,
            const double& dpinit,
            const int     lname,
            const int     ltype,
            const int     lstinit,
            const int     lcdinit){}
void FORTRAN(pusend,PUSEND)(const char* prompt,
            const int   len){}
void FORTRAN(nextwd,NEXTWD)(const char* prompt,
            const int   len){}
void FORTRAN(nextst,NEXTST)(const char* prompt,
                  char* ret,
            const int   lprompt,
            const int   lret){}
void FORTRAN(miscom,MISCOM)(const char *prompt,
                  int&  used,
            const int   len){}
void FORTRAN(nextlo,NEXTLO)(const char *prompt,
                  int&  used,
            const int   len){}
void FORTRAN(nexti,NEXTI)(const char *prompt,
                 int&  used,
               const int   len){}
void FORTRAN(nextf,NEXTF)(const char*   prompt,
                 double& used,
               const int     len){}
void FORTRAN(nextfi,NEXTFI)(const char* prompt,
                  char* ret,
            const int lprompt,
            const int lret){}
void FORTRAN(chkend,CHKEND)(const char *prompt,
                  int&  used,
            const int   len){}
void FORTRAN(selcta,SELCTA)(      int* flags,
                  int& nflags,
            const double* x,
            const double* y,
            const double* z,
            const int& blah){}
void FORTRAN(makind,MAKIND)(      int* flags,
            const int& natom,
            const int& nflags){} 
void FORTRAN(vclose,VCLOSE)(const int&  unit,
            const char* str,
                  int&  error,
            const int   lstr){}
void FORTRAN(printe,PRINTE)(double* renrl){}
void FORTRAN(tdynkin,TDYNKIN)(const int& dof,
                 double* xv,
                 double* yv,
                 double* zv){}
}

