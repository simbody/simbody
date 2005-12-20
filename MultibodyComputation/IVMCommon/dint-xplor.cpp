/**@file
 *
 * Dynamics in internal coordinates - xplor interface.
 * 
 * CDS - 9/15/99
 */
//FIX:  TODO: use the XplorWrap::xplorVars() singleton

#include "dint-xplor.h"
#include "dint-step.h"
#include "dint-atom.h"

#include <cdsAuto_arr.h>

#include <sthead.h>
#include <cdsList.h>
#include <cdsString.h>
#include <cdsSStream.h>
#include <fixedVector.h>
#include <string.h>
#include <cdsIomanip.h>

#include "prototypes.h"


void
XplorIVM::syncPos() {
    for (int i=1 ; i<=xplorVars.natom ; i++) {
        xplorVars.x[i-1] = atoms[i]->pos.x();
        xplorVars.y[i-1] = atoms[i]->pos.y();
        xplorVars.z[i-1] = atoms[i]->pos.z();
    }
}

void
XplorIVM::syncVel() {
    for (int i=1 ; i<=xplorVars.natom ; i++) {
        xplorVars.xv[i-1] = atoms[i]->vel.x();
        xplorVars.yv[i-1] = atoms[i]->vel.y();
        xplorVars.zv[i-1] = atoms[i]->vel.z();
    }
}

void
XplorIVM::syncDeriv() {
    atoms[0]->force.set( 0.0 );
    for (int i=1 ; i<=xplorVars.natom ; i++) {
        atoms[i]->force(0) = xplorVars.dx[i-1];
        atoms[i]->force(1) = xplorVars.dy[i-1];
        atoms[i]->force(2) = xplorVars.dz[i-1];
    }
}

void
XplorIVM::calcEnergy() {
    eCount++;
    syncPos();
    syncVel();
    FORTRAN(energy,ENERGY)(xplorVars.dx,xplorVars.dy,xplorVars.dz);
    //FIX: kinetic energy, total energy and temperature are wrong.
    FORTRAN(tdynkin,TDYNKIN)(dof(),xplorVars.xv,xplorVars.yv,xplorVars.zv,
            xplorVars.mass,xplorVars.imove);
    calcTemperature(); //calcs Ekinetic
    Epotential_ = xplorVars.renr[2]; //FIX: pretty unsafe
    Etotal_ = Ekinetic() + Epotential();
    syncDeriv();
}

//
// initialize InternalDynamics quantities from xplor values
//
void 
XplorIVM::initXplor()
{
    int natom = xplorVars.natom;

    for (int i=0 ; i<atoms.size() ; i++) //clean up any old mess
        delete atoms[i];
    atoms.resize(0);

    initAtoms(  natom,
                xplorVars.mass,
                xplorVars.kBoltzman); //create atoms array
    for (int bnum=0 ; bnum<xplorVars.nbond ; bnum++)
        if ( !bondExclude.contains( Pair(xplorVars.ib[bnum],
                                    xplorVars.jb[bnum]) ) )
        {
            atoms[ xplorVars.ib[bnum] ]->bonds.append( atoms[xplorVars.jb[bnum]] );
            atoms[ xplorVars.jb[bnum] ]->bonds.append( atoms[xplorVars.ib[bnum]] );
        }

    for (l_int i=1 ; i<=natom ; i++) {
        IVMAtom* a = atoms[i];
        a->pos = CDSVec3(xplorVars.x[i-1]  , xplorVars.y[i-1]  , xplorVars.z[i-1] );
        a->vel = CDSVec3(xplorVars.xv[i-1] , xplorVars.yv[i-1] , xplorVars.zv[i-1]);
        a->force = CDSVec3(xplorVars.dx[i-1] , 
                        xplorVars.dy[i-1] , 
                        xplorVars.dz[i-1]);
        a->fric = frictionCoeff() * xplorVars.fbeta[i-1] 
                    * a->mass * xplorVars.timeFac;
    }

    initTopology();
}

CDSString
XplorIVM::idAtom(int id) const {
    StringStream ret;
    ret << id << ' ';
    if (id<0)
        ret << "(unknown)";
    else if (id==0)
        ret << "(base atom)";
    else if (id>xplorVars.natom)
        ret << "[IVM atom]";
    else {
        int lid = id-1;
        ret << CDSString(xplorVars.resName+lid*4,4) << ' ' 
            << CDSString(xplorVars.resNum+lid*4,4)  << ' ' 
            << CDSString(xplorVars.atomName+lid*4,4);
    }
    ret << ends;
    return ret.str();
}

class VerboseList {
public:
    struct ValPair {
        int val ; const char* name; 
        ValPair() {}
        ValPair(int val, const char* name) : val(val), name(name) {}
    };
    CDSList<ValPair> list;
    VerboseList() {
    using namespace InternalDynamics;
    list.append( ValPair(   printCoords          , "coords"         ));
    list.append( ValPair(   printResetCM         , "resetcm"        ));
    list.append( ValPair(   printVelFromCartCost , "velfromcartcost"));
    list.append( ValPair(   printTemperature     , "temperature"    ));
    list.append( ValPair(   printEnergy          , "energy"         ));
    list.append( ValPair(   printCMVel           , "cmvel"          ));
    list.append( ValPair(   printNodeForce       , "nodeforce"      ));
    list.append( ValPair(   printNodePos         , "nodepos"        ));
    list.append( ValPair(   printNodeTheta       , "nodetheta"      ));
    list.append( ValPair(   printStepDebug       , "stepdebug"      ));
    list.append( ValPair(   printStepInfo        , "stepinfo"       ));
    list.append( ValPair(   printNodeDef         , "nodedef"        ));
    list.append( ValPair(   printLoopDebug       , "loopdebug"      ));
    list.append( ValPair(   printLoopInfo        , "loopinfo"       ));
    }
    int size() { return list.size(); }
    const ValPair& operator[](const int i) { return list[i]; }
};
static VerboseList verboseList;


static
ostream& help(ostream& ostr) 
{
    class Helper {
        ostream& ostr;
    public:
        Helper(ostream& ostr) : ostr(ostr) {}
        void line(const char* command, const char* text) {
            ostr.setf(ios::left);
            ostr << '\t' << setw(25) << command << " ! " << text << '\n';
        }
    } h(ostr);
                
    ostr << "<Internal-dynamics-statement>:== repeat [" << "\n";
    h.line("HELP"                 , "this help info");
    h.line("CSELect=<select>"     , "coords to write to traj file");
    h.line("DEPRed=<real>"        , "predicted drop in energy in minimization");
    h.line("BREAk=<sel1> <sel2>"  , "bonds to break");
    h.line("CBREak=<sel1> <sel2>" , "bonds to break - and constrain fixed");
    h.line("ENDTime=<real>"       , "time duration for integration.");
    h.line("ETOLerance=<real>"    , "energy tol for auto-stepsize adjust and min.");
    h.line("GTOLerance=<real>"    , "gradient tolerance for minimization");
    h.line("CTOLerance=<real>"    , "tolerance for bond contraints");
    h.line("CLOOp=<logical>"      , "use loop constraints");
    h.line("MAXEnergy=<real>"     , "max energy error before stepsize halving");
    h.line("MAXTimestep=<real>"   , "max factor for timestep increase (1.025)");
    h.line("MAXCalls=<int>"       , "max number of energy evals in minimization");
    h.line("RESPonse time=<real>" , "response time in multiples of timestep");
    h.line("NPRInt=<int>"         , "how often to print energy info");
    h.line("NSAVC=<int>"          , "how often to save coords to traj file");
    h.line("NSAVV=<int>"          , "how often to save velocities to traj file");
    h.line("NSTEp=<int>"          , "number of dynamics steps to take");
    h.line("NTRFrq=<int>"         , "how often to reset CM velocity/AM to zero");
    h.line("RESEt"                , "reset topology info");
    h.line("REUSe=<logical>"      , "reuse topology and traj. info");
    h.line("TRAJectory=<file>"    , "name of coord trajectory file");
    h.line("VELocity-file=<file>" , "name of velocity trajectory file");
    h.line("SCALe vel=<logical>"  , "use velocity scaling");
    h.line("ADJUst TS=<logical>"  , "use auto-timestep adjustment");
    h.line("TIMEstep=<real>"      , "size of dynamics timestep");
    h.line("TBATh=<real>"         , "bath temperature");
    h.line("FRICtion=<real>"      , "frict coeff. don't use with vel rescaling");
    h.line("AUTO TORSion"         , "auto hinge/node setup for torsion angles");
    h.line("FIX=<sel>"            , "fix atoms in space");
    h.line("GROUp=<sel>"          , "group atoms into node");
    h.line("HINGe=<type> <sel>"   , "selected atoms given hinge type");
    h.line("ITYPe=<type>"         , "solver selection (defaults to PC6)");
    h.line("     \tPC6"           , "6th order predictor-corrector");
    h.line("     \tGear"          , "Gear predictor-corrector");
    h.line("     \tMilne"         , "Milne predictor-corrector");
    h.line("     \tRungeKutta"    , "Runge-Kutta");
    h.line("     \tVerlet"        , "Velocity Verlet");
    h.line("     \tMinimizeCG"    , "conjugate-gradient with gradient test");
    h.line("     \tPowell"        , "Powell's method conjugate-gradient");
    h.line("     \tConMin"        , "simple conjugate-gradient");
    h.line("     \tSimplex"       , "Simplex method");
    h.line("PRINt=GROUP"          , "print atoms in each specified group");
    h.line("PRINt=HINGE"          , "print user-defined hinge types");
    h.line("     =NODE"           , "print contents and type of each node");
    h.line("     =verbose="       , " print verbose data:");
    h.line("     \treset"          , "");
    h.line("     \tcoords"         , "");
    h.line("     \tresetcm"        , "");
    h.line("     \tvelfromcartcost", "");
    h.line("     \tenergy"            , "");
    h.line("     \ttemperature"    , "");
    h.line("     \tcmvel"            , "");
    h.line("     \tnodeforce"      , "");
    h.line("     \tnodepos"    , "");
    h.line("     \tnodetheta"      , "");
    h.line("     \tnodedef"        , "");
    h.line("     \tstepinfo"       , "");
    h.line("     \tstepdebug"      , "");
    h.line("     \tloopinfo"       , "");
    h.line("     \tloopdebug"      , "");
    ostr << "\n\tafter the dynamics internal statement has been processed,\n"
         << "\tvariables $DINT_TIME and $DINT_STEPS are set, as appropriate.\n";
    return ostr;
}

CDSString
formatHingeSpec(const IVM* ivm, const InternalDynamics::HingeSpec& spec)
{
    StringStream ret;
    ret.setf( ios::left );
    // s.setf( ios::left );
    ret << " " << setw(9) << spec.type << "   { ";
    for (int i=0 ; i<spec.aList.size() ; i++)
        ret << ivm->idAtom(spec.aList[i])  
            << ((i<spec.aList.size()-1)?", ":" }");
    if ( spec.atom0>=0 )
        ret << " ( " << ivm->idAtom(spec.atom0) << " ) ";
    if ( spec.atom1>=0 )
        ret << "( " << ivm->idAtom(spec.atom1) << " )";
    ret << ends;
    // s << '\n';
    // s.flags( oflags );
    return ret.str();
}

ostream& 
XplorIVM::info( ostream&    ostr,
                double&     finalTime,
                const double& stepsize,
                const int&  nprint,
                const int&  nsavc,
                const int&  nsavv,
                const int&  reuseIntVar,
                const int&  steps,
                const int&  nCMresets,
                const char* trajFilename,
                const char* velFilename) 
{
    // class Info {
    //   ostream& ostr;
    // public:
    //   Info(ostream& ostr) : ostr(ostr) {}
    //   void line(const char* command,
    //         const char* text   ) {
    //    ostr.setf(ios::left);
    //    ostr << '\t' << setw(25) << command << " ! " << text << '\n';
    //   }
    // } h(ostr);
                
    ostr << "Current Settings for Internal Variable Dynamics:\n";
    ostr << "DEPRed= " << dEpred() << '\n';
    ostr << "Bonds to break:\n";
    for (int i=0 ; i<bondExclude.size() ; i++) 
        ostr << '\t' << idAtom(bondExclude[i].a) 
             << ' ' << idAtom(bondExclude[i].b) << '\n';
    ostr << "Bonds to constrain fixed:\n";
    for (l_int i=0 ; i<constraintList.size() ; i++) 
        ostr << '\t' << idAtom(constraintList[i].a) << ' ' 
             << idAtom(constraintList[i].b) << '\n';
    ostr << "ENDTime= " << finalTime * xplorVars.timeFac<< '\n';
    ostr << "ETOLerance= " << Etolerance() << '\n';
    ostr << "GTOLerance= " << Gtolerance() << '\n';
    ostr << "CTOLerance= " << Ctolerance() << '\n';
    ostr << "CLOOp= "  << (useLengthConstraints()?"True":"False") << '\n';
    ostr << "MAXEnergy= "  << maxDeltaE() << '\n';
    ostr << "MAXTimestep= " << maxTSFactor() << '\n';
    ostr << "MAXCalls= " << maxCalls() << '\n';
    ostr << "RESPonse time= " << responseTime() << '\n';
    ostr << "NPRInt= " << nprint << '\n';
    ostr << "NSAVC= " << nsavc << '\n';
    ostr << "NSAVV= " << nsavv << '\n';
    ostr << "NSTEp= " << steps << '\n';
    ostr << "NTRFrq= " << nCMresets << '\n';
    ostr << "REUSe= " << (reuseIntVar?"True":"False") << '\n';
    ostr << "SCALe vel= " << (scaleVel()?"True":"False") << '\n';
    ostr << "ADJUst Timestep= " << (adjustTS()?"True":"False") << '\n';
    ostr << "TRAJectory= " << trajFilename << '\n';
    ostr << "VELocity-file= " << velFilename << '\n';
    ostr << "TIMEstep= " << stepsize * xplorVars.timeFac << '\n';
    ostr << "TBATh= " << bathTemp() << '\n';
    ostr << "FRICtion= " << frictionCoeff() / xplorVars.timeFac << '\n';
    ostr << "Fixed atoms:";
    for (l_int i=0 ; i<groupList.size() ; i++) 
        if ( groupList[i].contains(0) ) 
            for (int j=0 ; j<groupList[i].size() ; j++)
                ostr << ' ' << idAtom(groupList[i][j]);
    ostr << '\n';
    ostr << "Nonfixed atoms held rigid wrt each other:\n";
    for (l_int i=0 ; i<groupList.size() ; i++)
        if ( !groupList[i].contains(0) ) {
            ostr << '\t';
            for (int j=0 ; j<groupList[i].size() ; j++)
                ostr << idAtom(groupList[i][j]) << ' ';
            ostr << '\n';
        }
    ostr << "hinge specifications:\n";
    for (l_int i=0 ; i<hingeList.size() ; i++) {
        ostr << "\t\t" << formatHingeSpec(this,hingeList[i]) << '\n';
    }
    ostr << "ITYPe= " << solverType << '\n';
    ostr << "Verbose info to print:\n\t";
    for (l_int i=0 ; i<verboseList.size() ; i++)
        if ( verbose()&verboseList[i].val )
            ostr << ' ' << verboseList[i].name;
    ostr << '\n';

    return ostr;
}

void
XplorIVM::outputStepInfo(const int     step,
                         const double& stepsize,
                         const CDSString& type)
{
    FMTFLAGS_TYPE oflags = cout.flags();
    // CDSString type = integrateType + " ";
    cout.setf(ios::left,ios::adjustfield);
    cout << " -- " << setw(11) << setfill('-') << type+" " 
         << setfill(' ') << setw(0);
    cout.setf(ios::right,ios::adjustfield);
    cout << "-- step=" << setw(7) << step;
    if ( minimization() )
        cout << " --- stepsize=" 
             << setw(10) << setprecision(5) << stepsize;
    else
        cout << " --- timestep=" 
             << setw(10) << setprecision(5) << stepsize*xplorVars.timeFac;
    cout << " --- energy evals="
         << setw(5) << eCount
         << " --\n";
    if ( !minimization() ) {
        cout.setf(ios::left,ios::adjustfield);
        cout << " | E(kin)+E(poten)=" << setprecision(7) << setw(14) << Etotal()
             << " E(kin)=" << setprecision(7) << setw(14) << Ekinetic()
             << " temperature=" << setprecision(7) << setw(11) << currentTemp()
             << "|\n";
        cout.setf(ios::right,ios::adjustfield);
    }

    FORTRAN(printe,PRINTE)(xplorVars.renr);
    // cout.setf(ios::right);
    // cout << " -------------------------- itype=" 
    //      << setw(10) << integrateType
    //      << " -----------------------------------\n";
    cout << " ---------------------------------" 
         << "----------"
         << "------------------------------------\n";
    printCM();
    cout.flags(oflags);
    eCount = 0; //reset
}
 
// add fixed atoms from imove array
void
XplorIVM::fixAtomsFromExternal(const int   natom,
                               const long* imove)
{
    CDSList<int> aList; aList.append( 0 );
    for (int i=0 ; i<natom ; i++)
        if ( imove[i] == 1 )
            aList.append( i+1 );
    if (aList.size() > 1)
        groupList.append( aList );
}


void
XplorIVM::entry()
{
    // cout << "---- energy terms before conversion to internal coords ----\n";
    // FORTRAN(energy)();
    // FORTRAN(tdynkin)(3*natom,xv,yv,zv);
    // FORTRAN(printd)(renr);

    int natom = xplorVars.natom;

    static double stepsize=0.;
    if (stepsize==0.0)                      //FIX: not appropriate for minimization
        stepsize = 1e-3/xplorVars.timeFac;  //units: A/(AKMA time unit)
    double    finalTime=0.0;
    int       steps=0;
    long      done=0;
    CDSString prompt = "Internal Dynamics>";
    CDSString trajFilename;
    CDSString velFilename;
    int    formattedOutput=0;
    CDSString format="12Z6      ";
    double scale=10000.0;
    double offset=800;
    int    nCMresets=0;  //how often to remove CM translation/rotation
    int    nprint=1;
    int    nsavc=1;//??
    int    nsavv=1;//??
    int    reuseIntVar=0;
    long   crdnum = natom;
    auto_arr<long> crdind(new long[crdnum]);
    fixAtomsFromExternal(natom,xplorVars.imove);
    for (int i=0 ; i<natom ; i++) crdind[i] = 1;
        FORTRAN(pusend,PUSEND)(prompt,prompt.length());
    while (! done) {
        FORTRAN(nextwd,NEXTWD)(prompt,prompt.length());
        long used=0;
        FORTRAN(miscom,MISCOM)(prompt,used,prompt.length());
        const CDSString wd(xplorVars.wd,4);
        if (! used) {
            if ( xplorVars.wd[0] == '?' || wd == "INFO" ) {
                info(cout,finalTime,stepsize,nprint,nsavc,nsavv,reuseIntVar,steps,
                nCMresets,trajFilename,velFilename);
            } else if ( wd == "HELP" ) {
                cout << help << '\n';
                cout.flush();
                //       cout << endl;
            } else if ( wd == "CSEL" ) {
                FORTRAN(selcta,SELCTA)(crdind,crdnum,xplorVars.x,xplorVars.y,xplorVars.z,1);
            } else if ( wd == "BASE" ) {
                auto_arr<long> batom(new long[natom]);
                long nbatoms=natom;
                FORTRAN(selcta,SELCTA)(batom,nbatoms,xplorVars.x,xplorVars.y,xplorVars.z,1);
                FORTRAN(makind,MAKIND)(batom,natom,nbatoms); 
                for (int i=0 ; i<nbatoms ; i++) 
                oldBaseAtoms.append( batom[i] );
            //     if ( nbatoms>1 )
            //       cout << "Internal Dynamics> BASE: (sel1):\n"
            //        << "\tselection must be for a single atom\n";
            //     else 
                //     oldBaseAtoms.append( batom[0] );
            } else if ( wd == "BREA" ) {
                auto_arr<long> batom1(new long[natom]);
                long nbatoms1=natom;
                FORTRAN(selcta,SELCTA)(batom1,nbatoms1,xplorVars.x,xplorVars.y,xplorVars.z,1);
                FORTRAN(makind,MAKIND)(batom1,natom,nbatoms1); 
                auto_arr<long> batom2(new long[natom]);
                long nbatoms2=natom;
                FORTRAN(selcta,SELCTA)(batom2,nbatoms2,xplorVars.x,xplorVars.y,xplorVars.z,1);
                FORTRAN(makind,MAKIND)(batom2,natom,nbatoms2); 
                if ( nbatoms1>1 || nbatoms2>1 )
                    cout << "Internal Dynamics> BREAk: (sel1) (sel2):\n"
                         << "\teach selection must be for a single atom\n";
                else 
                    bondExclude.append( Pair(batom1[0],batom2[0]) );
            } else if ( wd == "CBRE" ) {
                auto_arr<long> batom1(new long[natom]);
                long nbatoms1=natom;
                FORTRAN(selcta,SELCTA)(batom1,nbatoms1,xplorVars.x,xplorVars.y,xplorVars.z,1);
                FORTRAN(makind,MAKIND)(batom1,natom,nbatoms1); 
                auto_arr<long> batom2(new long[natom]);
                long nbatoms2=natom;
                FORTRAN(selcta,SELCTA)(batom2,nbatoms2,xplorVars.x,xplorVars.y,xplorVars.z,1);
                FORTRAN(makind,MAKIND)(batom2,natom,nbatoms2); 
                if ( nbatoms1>1 || nbatoms2>1 )
                    cout << "Internal Dynamics> BREAk: (sel1) (sel2):\n"
                         << "\teach selection must be for a single atom\n";
                else {
                    bondExclude.append( Pair(batom1[0],batom2[0]) );
                    constraintList.append( Pair(batom1[0],batom2[0]) );
                }
            } else if ( wd == "DEPR" ) {
                FORTRAN(nextf,NEXTF)("DEPRedict=",dEpred_,10);
            } else if ( wd == "ETOL" ) {
                FORTRAN(nextf,NEXTF)("ETOLerance=",Etolerance_,11);
            } else if ( wd == "GTOL" ) {
                FORTRAN(nextf,NEXTF)("GTOLerance=",Gtolerance_,11);
            } else if ( wd == "CTOL" ) {
                FORTRAN(nextf,NEXTF)("CTOLerance=",Ctolerance_,11);
            } else if ( wd == "CLOO" ) {
                long tmp;
                FORTRAN(nextlo,NEXTLO)("CLOOp=",tmp,6);
                useLengthConstraints_= (tmp != 0);
            } else if ( wd == "MAXE" ) {
                FORTRAN(nextf,NEXTF)("MAXEnergy error=",maxDeltaE_,16);
            } else if ( wd == "MAXT" ) {
                FORTRAN(nextf,NEXTF)("MAXTimestep factor=",maxTSFactor_,19);
            } else if ( wd == "RESP" ) {
                FORTRAN(nextf,NEXTF)("RESPonse time=",responseTime_,14);
            } else if ( wd == "MAXC" ) {
                long tmp;
                FORTRAN(nexti,NEXTI)("MAXCalls=",tmp,9); 
                maxCalls_ = tmp;
            } else if ( wd == "NPRI" ) {
                long tmp;
                FORTRAN(nexti,NEXTI)("NPRInt=",tmp,7); 
                nprint = (tmp<1) ? 1 : tmp;
            } else if ( CDSString(xplorVars.wd,5) == "NSAVC" ) {
                long tmp;
                FORTRAN(nexti,NEXTI)("NSAVC=",tmp,6);
                nsavc = tmp;
            } else if ( CDSString(xplorVars.wd,5) == "NSAVV" ) {
                long tmp;
                FORTRAN(nexti,NEXTI)("NSAVV=",tmp,6);
                nsavv = tmp;
            } else if ( wd == "NSTE" ) {
                long tmp;
                FORTRAN(nexti,NEXTI)("NSTEp=",tmp,6);
                steps = tmp;
            } else if ( wd == "NTRF" ) {
                long tmp;
                FORTRAN(nexti,NEXTI)("NTRFrq=",tmp,7);
                nCMresets = tmp;
            } else if ( wd == "RESE" ) { //reset
                //       stepsize = 0.0;
                bondExclude.resize(0);
                groupList.resize(0);
                hingeList.resize(0);
                constraintList.resize(0);
                oldBaseAtoms.resize(0);
                // read atoms fixed externally
                fixAtomsFromExternal(natom,xplorVars.imove);
            } else if ( wd == "REUS" ) {
                long tmp;
                FORTRAN(nextlo,NEXTLO)("REUSe=",tmp,6);
                reuseIntVar = tmp;
            } else if ( wd == "SCAL" ) {
                long tmp;
                FORTRAN(nextlo,NEXTLO)("SCALe velocity=",tmp,15);
                scaleVel_ = (tmp != 0);
            } else if ( wd == "ADJU" ) {
                long tmp;
                FORTRAN(nextlo,NEXTLO)("ADJUst timestep=",tmp,16);
                adjustTS_ = (tmp != 0);
            } else if ( wd == "ASCI" ) {
                long tmp;
                FORTRAN(nextlo,NEXTLO)("ASCIi=",tmp,strlen("ASCIi="));
                formattedOutput = tmp;
            } else if ( wd == "TRAJ") {
                const int lfstr=80;
                char fstr[lfstr]; 
                FORTRAN(nextfi,NEXTFI)("TRAJectory-file=",fstr,16,lfstr);
                for (int i=0 ; i<lfstr ; i++) 
                    if ( fstr[i]==' ' ) fstr[i] = '\0';       
                        trajFilename = fstr;
            } else if ( wd == "VELO") {
                const int lfstr=80;
                char fstr[lfstr]; 
                FORTRAN(nextfi,NEXTFI)("VELocity-file=",fstr,14,lfstr);
                for (int i=0 ; i<lfstr ; i++) 
                    if ( fstr[i]==' ' ) fstr[i] = '\0';       
                        velFilename = fstr;
            } else if ( wd == "TIME" ) {
                FORTRAN(nextf,NEXTF)("TIMEstep=",stepsize,9);
                stepsize /= xplorVars.timeFac;
            } else if ( wd == "ENDT" ) {
                FORTRAN(nextf,NEXTF)("ENDTime=",finalTime,8);
                finalTime /= xplorVars.timeFac;
            } else if ( wd == "STEP" ) {
                FORTRAN(nextf,NEXTF)("STEPsize=",stepsize,9);
            } else if ( wd == "TBAT" ) {
                FORTRAN(nextf,NEXTF)("TBATh=",bathTemp_,6);
            } else if ( wd == "FRIC" ) {
                FORTRAN(nextf,NEXTF)("FRICtion coefficient=",frictionCoeff_,21);
            //FIX: need to be able to clear groupings
            } else if ( wd == "AUTO" ) { 
                // for ``automate'': auto-setup torsion dynamics
                const int ltypestr=15;
                char typestr[ltypestr]; 
                FORTRAN(nextst,NEXTST)("TYPE=",typestr,5,ltypestr);
                for (int i=0 ; i<ltypestr ; i++) 
                    if ( typestr[i]==' ' ) typestr[i] = '\0';       
                        initXplor();
                try {
                    groupTorsion();
                    InternalDynamics::HingeSpec hingeSpec;
                    hingeSpec.type = "torsion";
                    for (l_int i=1 ; i<=natom ; i++) {
                        bool inHinge=0;
                        for (int j=0 ; j<hingeList.size() ; j++)
                            if ( hingeList[j].aList.contains(i) ) {
                                inHinge=1;
                                break;
                            }
                        if ( !inHinge )
                            hingeSpec.aList.append( i );
                    }
                    hingeList.append( hingeSpec );
                }
                catch ( InternalDynamics::Exception e ) {
                    cout << "Caught exception: " << e.mess << '\n';
                    cout << "Failed to automatically determine torsion angles.\n";
                }
            } else if ( CDSString(xplorVars.wd,3) == "FIX" ) {
                auto_arr<long> flags(new long[xplorVars.natom]);
                long  nflags;
                FORTRAN(selcta,SELCTA)(flags,nflags,
                        xplorVars.x,xplorVars.y,xplorVars.z,1);
                FORTRAN(makind,MAKIND)(flags,xplorVars.natom,nflags); 
                CDSList<int> aList;
                aList.append( 0 );
                for (int i=0 ; i<nflags ; i++)
                    aList.append( flags[i] );
                if (aList.size() > 1)
                    groupList.append( aList );
            } else if ( wd == "GROU" ) {
                auto_arr<long> flags(new long[xplorVars.natom]);
                long nflags;
                FORTRAN(selcta,SELCTA)(flags,nflags,
                        xplorVars.x,xplorVars.y,xplorVars.z,1);
                FORTRAN(makind,MAKIND)(flags,xplorVars.natom,nflags); 
                CDSList<int> aList;
                for (int i=0 ; i<nflags ; i++)
                    aList.append( flags[i] );
                if (aList.size() > 1)
                    groupList.append( aList );
            } else if ( wd == "HING" ) {
                const int ltypestr=80;
                char typestr[ltypestr]; 
                FORTRAN(nextst,NEXTST)("TYPE=",typestr,5,ltypestr);
                for (int i=0 ; i<ltypestr ; i++) 
                    if ( typestr[i]==' ' ) typestr[i] = '\0';       
                auto_arr<long> flags(new long[xplorVars.natom]);
                long nflags;
                FORTRAN(selcta,SELCTA)(flags,nflags,
                        xplorVars.x,xplorVars.y,xplorVars.z,1);
                FORTRAN(makind,MAKIND)(flags,xplorVars.natom,nflags); 
                InternalDynamics::HingeSpec hingeSpec;
                hingeSpec.type = typestr;
                hingeSpec.type.downcase();
                for (l_int i=0 ; i<nflags ; i++)
                    hingeSpec.aList.append( flags[i] );
                bool ok=0;
                if ( hingeSpec.type == "bendtorsion" || hingeSpec.type == "bend") {
                    auto_arr<long> flags(new long[xplorVars.natom]);
                    long nflags;
                    FORTRAN(selcta,SELCTA)(flags,nflags,
                            xplorVars.x,xplorVars.y,xplorVars.z,1);
                    FORTRAN(makind,MAKIND)(flags,xplorVars.natom,nflags); 
                    if ( nflags != 1 ) 
                        cout << "Internal Dynamics> HINGe> second selection should "
                             << "contain a single atom.\n";
                    else {
                        hingeSpec.atom0 = flags[0];
                        ok = 1;
                    }
                    if ( ok ) {
                        ok = 0;
                        FORTRAN(selcta,SELCTA)(flags,nflags,
                                xplorVars.x,xplorVars.y,xplorVars.z,1);
                        FORTRAN(makind,MAKIND)(flags,xplorVars.natom,nflags); 
                        if ( nflags != 1 ) 
                            cout << "Internal Dynamics> HINGe> third selection should "
                                 << "contain a single atom.\n";
                        else {
                            hingeSpec.atom1 = flags[0];
                            ok = 1;
                        }
                    }
                } else
                ok = 1;
                if (ok)
                    hingeList.append( hingeSpec );
            } else if ( wd == "ITYP" ) {
                const int largstr=80;
                char argstr[largstr]; 
                FORTRAN(nextst,NEXTST)("Solver type=",argstr,17,largstr);
                for (int i=0 ; i<largstr ; i++) 
                    if ( argstr[i]==' ' ) argstr[i] = '\0';       
                        solverType = argstr;
            } else if ( wd == "PRIN" ) {
                const int largstr=80;
                char argstr[largstr]; 
                FORTRAN(nextst,NEXTST)("ARG=",argstr,5,largstr);
                for (int i=0 ; i<largstr ; i++) 
                    if ( argstr[i]==' ' ) argstr[i] = '\0';       
                if ( CDSString(argstr) == "GROUP" ) {
                    cout << "groups: \n";
                    for (int i=0 ; i<groupList.size() ; i++) {
                        for (int j=0 ; j<groupList[i].size() ; j++)
                            cout << idAtom(groupList[i][j]) << ' ';
                        cout << '\n';
                    }
                }
                //     } else if ( CDSString(argstr) == "NODE" ) {
                //       FORTRAN(nextlo,NEXTLO)("PRINt NODE=",printNodeInfo,11);
                if ( CDSString(argstr) == "HINGE" ) {
                    cout << "Hinge Type   Hinges with the following atoms: \n";
                    FMTFLAGS_TYPE oflags = cout.flags();
                    cout.setf( ios::left );
                    for (int i=0 ; i<hingeList.size() ; i++) 
                        cout << "\t" 
                             << formatHingeSpec(this,hingeList[i]) << '\n';
                //         cout << " " << setw(9) << hingeList[i].type << "   "
                //          << hingeList[i].aList << '\n';
                    cout.flags( oflags );
                } else if ( CDSString(argstr) == "VERBOSE" ) {
                    const int largstr=80;
                    char argstr[largstr]; 
                    FORTRAN(nextst,NEXTST)("PRINt VERBOSE=",argstr,14,largstr);
                    for (int i=0 ; i<largstr ; i++) 
                        if ( argstr[i]==' ' ) argstr[i] = '\0';    
                    CDSString s(argstr); s.downcase();
                    if ( s == "reset"   )    verbose_ =0;
                    else if ( s == "all"   ) verbose_ =0xffff;
                    else {
                        int ok=0;
                        for (int i=0 ; i<verboseList.size() ; i++)
                            if ( s == verboseList[i].name ) {
                                verbose_ |= verboseList[i].val;
                                ok=1;
                            }
                        if ( !ok )
                            cout << "PRINt VERBOSE: invalid input\n";
                    }
                }
            } else 
                FORTRAN(chkend,CHKEND)(prompt,done,prompt.length());
        }
    }

    if ( steps < 1 && finalTime <= 0. )
        return;
    if ( stepsize <= 0. ) {
        cout << "stepsize=0. No steps taken.\n";
        return;
    }

    try {
        if (nCMresets)
            resetCMflag=1;

        if ( !reuseIntVar ) initXplor();
        initDynamics(1);

        // print out energy info
        outputStepInfo(0,stepsize,solverType);

        long coordUnit=0;
        long velUnit=0;
        long error=0;
        long qfirsc=1;//is this first pass through loop? - reset by writtc
        long vfirsc=1;//is this first pass through loop? - reset by writtc
        if (trajFilename != "") {
            if (formattedOutput) 
                FORTRAN(assfil,ASSFIL)(trajFilename,coordUnit,"WRITE","FORMATTED",error,
                        trajFilename.length(),5,9);
            else
                FORTRAN(assfil,ASSFIL)(trajFilename,coordUnit,"WRITE","UNFORMATTED",error,
                        trajFilename.length(),5,11);
            FORTRAN(writtc,WRITTC)("CORD",natom,xplorVars.x,xplorVars.y,xplorVars.z,crdnum,crdind,0,qfirsc,
                    stepsize,nsavc,coordUnit,formattedOutput,
                    format,scale,offset,4,format.length());
        }
        if (velFilename != "") {
            if (formattedOutput) 
                FORTRAN(assfil,ASSFIL)(velFilename,velUnit,"WRITE","FORMATTED",error,
                        velFilename.length(),5,9);
            else
                FORTRAN(assfil,ASSFIL)(velFilename,velUnit,"WRITE","UNFORMATTED",error,
                        velFilename.length(),5,11);
            FORTRAN(writtc,WRITTC)("VELD",natom,xplorVars.xv,xplorVars.yv,xplorVars.zv,crdnum,crdind,0,vfirsc,
                    stepsize,nsavv,velUnit,formattedOutput,
                    format,scale,offset,4,format.length());
        }
        FORTRAN(makind,MAKIND)(crdind,natom,crdnum);//change sel setup
            
        bool   done=0;
        int    iter=0;
        double time=0.0;
        while ( !done ) {
            iter++;
            //   for (l_int i=1 ; i<=steps && !done ; i++) {
            if ( nCMresets && iter%nCMresets==0 )
                resetCMflag=1;
            
            time += stepsize;
            try {
                step(stepsize);
            }
            catch ( Solver::Finished ) {
                done=1;
            }
            if ( finalTime>0.0 && time >= finalTime )
                done=1;
            if ( steps>0 && iter>=steps)
                done=1;
            
            
            syncPos();
            syncVel();

            if ( nsavc>0 && coordUnit>0 && iter%nsavc==0 ) 
                // print out coord info
                FORTRAN(writtc,WRITTC)("CORD",natom,xplorVars.x,xplorVars.y,xplorVars.z,crdnum,crdind,iter,qfirsc,
                        stepsize,nsavc,coordUnit,formattedOutput,
                        format,scale,offset,4,format.length());
            if ( nsavv>0 && velUnit>0 && iter%nsavv==0 )
                FORTRAN(writtc,WRITTC)("VELD",natom,xplorVars.xv,xplorVars.yv,xplorVars.zv,crdnum,crdind,iter,vfirsc,
                        stepsize,nsavv,velUnit,formattedOutput,
                        format,scale,offset,4,format.length());
            
            if (iter%nprint==0)
                outputStepInfo(iter,stepsize,solverType);
        }
        const char* dum = "";
        FORTRAN(declar,DECLAR)("DINT_TIME" ,"DP",
                dum,dum,time*xplorVars.timeFac,9,2,0);
        FORTRAN(declar,DECLAR)("DINT_STEPS","DP",dum,dum,(double)iter,10,2,0);
        if ( coordUnit>0 ) FORTRAN(vclose,VCLOSE)(coordUnit,"KEEP",error,4);
    }
    catch ( InternalDynamics::Exception e) {
        cout << "InternalDynamics> uncaught exception: " << e.mess << endl;
        throw;
    }
    catch ( ... ) {
        cout << "InternalDynamics> unknown uncaught exception" << endl;
        throw;
    }
} 

extern "C" void
FORTRAN(internal_dynamics,INTERNAL_DYNAMICS)
       (const int&          natom,
        double *const       x,
        double *const       y,
        double *const       z,
        double *const       xv,
        double *const       yv,
        double *const       zv,
        double *const       dx,
        double *const       dy,
        double *const       dz,
        const int&          nenert,
        long*               qener,    
        const char* const   aner,    
        const double*       renr,    
        const int&          nbond,
        const long*         ib,
        const long*         jb,
        double*             amass,
        double*             fbeta,
        char* const         resName,
        char* const         resNum,
        char* const         atomName,
        char * const        wd,
        const long*         imove,
        const double&       timeFac,
        const double&       kBoltzman,
        const int           lwd)
{
    //sgi implementation of sync_with_stdio multiple times results in memory leak
    static bool first=1; if (first) { ios::sync_with_stdio(); first=0; }

    double dummy=0.0;
    long dummyl=0;

    static XplorIVM* ivm=0;
    if ( ivm==0 ) 
        ivm = new XplorIVM( XplorVars(natom,x,y,z,xv,yv,zv,dx,dy,dz,0,
                            resName,resNum,atomName,0,
                            nenert,qener,aner,renr,
                            nbond,ib,jb,amass,fbeta,0,imove,wd,timeFac,
                            kBoltzman,dummy,dummyl,dummyl) );

    ivm->entry();
}

