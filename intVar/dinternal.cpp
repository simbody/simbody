/**@file
 *
 *  Dynamics in internal coordinates.
 * 
 *  ****add references****
 * 
 *  8/25/99 - C. Schwieters
 * 
 *  outline:
 *  
 *  -need simple constraint mechanism (for these coordinates)
 *  -add more general constraint mechanism (for arbitrary constraints)
 *  -convert from internal to cartesian coordinates
 */

#include "dinternal.h"
#include "AtomTree.h"
#include "AtomClusterNode.h"

#include "dint-step.h"
#include "LengthConstraints.h"

#include "linemin.h"
#include "dint-atom.h"
#include "vec3.h"
#include "vec4.h"

#include "sthead.h"
#include "cdsMath.h"
#include "cdsString.h"
#include "cdsSStream.h"
#include "cdsVector.h"
#include "fixedVector.h"
#include "subVector.h"
#include "fixedMatrix.h"
#include "matrixTools.h"
#include "cdsIomanip.h"
#include "cdsFstream.h"
#include "cdsAuto_arr.h"

using namespace InternalDynamics;
using MatrixTools::inverse;

typedef FixedVector<double,6> Vec6;
typedef FixedMatrix<double,6> Mat66;

IVM::IVM() 
  : tree_(0),                  solver_(0),
    dof_(0),                   dim_(0),            verbose_( printStepInfo ), 
    useLengthConstraints_(0),  adjustTS_(1),       scaleVel_(1),
    resetCMflag(0),
    Etotal_(0.0),              Epotential_(0.0),   Ekinetic_(0.0),
    Etolerance_(0),            Gtolerance_(1e-6),  Ctolerance_(1e-8),
    minStepSize_(1e-10),       maxTSFactor_(1.025),
    maxDeltaE_(100),           maxCalls_(5000),
    kBoltzmann_(0.0),          bathTemp_(298.0),   frictionCoeff_(0.0),
    responseTime_(0.0),        currentTemp_(-1.0),
    dEpred_(0.001),            solverType( "PC6" ),
    rvecSize_(defaultRVecSize), rvecProd_(defaultRVecProd),
    vecVec3Size_(defaultVecVec3Size), vecVec3Prod_(defaultVecVec3Prod)
{
    tree_ = new AtomTree(this);
}

IVM::~IVM()
{
    delete solver_;
    delete tree_;
    for (int i=0 ; i<atoms.size() ; i++)
        delete atoms[i];
}

//
// Y is used for length constraints.
//
void
IVM::calcY() {
    tree()->calcY();
}

bool
IVM::minimization() const { 
    if ( !solver_ ) {
        cout << "IVM::minimization: solver not yet defined.\n" << flush;
        throw Exception("IVM::minimization: solver not yet defined.");
    }
    return getSolver()->minimization(); 
}

int
IVM::defaultRVecSize(const RVec& v) {
    return v.size();
}

double
IVM::defaultRVecProd(const RVec& v1, const RVec& v2) {
    return dot(v1,v2);
}

int
IVM::defaultVecVec3Size(const CDSVector<Vec3>& v) {
    return v.size(); // FIX: *3???
}

double
IVM::defaultVecVec3Prod(const CDSVector<Vec3>& v1, const CDSVector<Vec3>& v2) {
    // return dot(v1,v2);// didn't work??
    double sum=0.;
    for (int i=0 ; i<v1.size() ; i++)
        sum += v1[i]*v2[i];
    return sum;
}

//
// recursively group all atoms in a node with neighboring nodes 
// appropriate for torsion angle dynamics. 
// does not create a group if there's only a single atom.
//
void
IVM::groupTorsion(const AtomClusterNode* n) {
    CDSList<int> aList;
    const IVMAtom* a;

    //fixed atoms
    if ( n->isGroundNode() ) {
        if ( n->getAtom(1) ) 
            for (int i=0 ; ( a=n->getAtom(i) ) ; i++)
                if ( a->index>=0 )
                    aList.append( a->index );

        const AtomClusterNode* c;
        for (int i=0 ; ( c=n->getChild(i) ) ; i++) { 
            if ( c->getChild(0) || c->getAtom(1) )     
                //has children or contains more than one atom
                groupTorsion( c );       
            else
                if (c->getParentAtom()->index!=0 && c->getAtom(0)->index>=0)
                    //one atom, no children and not free atom : group
                    aList.append( c->getAtom(0)->index );
        }
    } else { //not ground node
        //base node w/ 1 or 2 atoms- w/ no bond to base and one or fewer children
        if ( n->isBaseNode() && !n->getAtom(2) &&
             n->getParentAtom()->index==0 && !n->getChild(1) )
        { 
            if ( n->getAtom(1) || n->getChild(0) )
                for (int i=0 ; ( a=n->getAtom(i) ) ; i++)
                    if ( a->index>=0 )
                        aList.append( a->index );
            n = n->getChild(0);
        }

        if (n) {
            const AtomClusterNode* c;
            for (int i=0 ; ( c=n->getChild(i) ) ; i++) {
                if ( c->getChild(0) || c->getAtom(1) ) {    
                    //has children or contains more than one atom
                    groupTorsion( c );       
                } else {
                    //no children: group
                    for (int j=0 ; ( a=c->getAtom(j) ) ; j++)
                        if ( a->index>=0 ) {
                            bool add = 1;
                            for (int k=0 ; k<hingeList.size() ; k++)
                                if ( hingeList[k].aList.contains(a->index) &&
                                     hingeList[k].type != "torsion")
                                {
                                    add = 0;
                                    break;
                                }
                            if ( add )
                                aList.append( a->index );
                        }
                }
            }

            if ( aList.size() || n->getAtom(1) )
                for (int i=0 ; ( a=n->getAtom(i) ) ; i++)
                    if ( a->index>=0 )
                        aList.append( a->index );
        }
    }

    if ( aList.size() ) groupList.append( aList );
}

//
// GroupList is overwritten, but pre-existing groupings are honored.
//
void
IVM::groupTorsion()
{
    groupList.resize(0);   //clobber existing groupList
    groupTorsion( tree()->nodeTree[0][0] );
}

void
IVM::step(double& timeStep) {
    bool ok=0;
    do {
        try {
            getSolver()->step(timeStep);
            ok=1;
        }
        catch ( Exception a ) {
            if ( String(a.mess).contains("large time-step") ) {
                cout << "InternalDynamics::step: large timestep detected. Halving.\n";
                timeStep *= 0.5;
                if ( timeStep < minStepSize() ) {
                    cout << "InternalDynamics::step: stepsize too small." << endl;
                    throw InternalDynamics::Exception("stepsize too small");
                }
            } else {
                cout << "InternalDynamics::step: caught assertion failure:\n"
                << "\t" << a.mess << '\n';
                cout.flush();
                throw;
            }
        }
    } while ( !ok );
}
 
//void
//InternalDynamics::approxTemperature()
//{
// approxTemp = 0.0;
// for (int i=0 ; i<tree->nodeTree.size() ; i++)
//   for (int j=0 ; j<tree->nodeTree[i].size() ; j++) 
//     approxTemp += tree->nodeTree[i][j]->approxKE();
//     
// approxTemp /= (dof * kBoltzman);
//}
 
void
IVM::calcTemperature() {
    Ekinetic_ = 0.0;
    int enumerateKE=0;
    if (enumerateKE) 
        cout << "thermal energy contributions (by node):\n";

    for (int i=0 ; i<tree()->nodeTree.size() ; i++)
        for (int j=0 ; j<tree()->nodeTree[i].size() ; j++) {
            if (enumerateKE) 
                cout << i << ' ' << j << ' ';
            double ke = tree()->calcClusterKineticEnergy(i,j);
            if (enumerateKE) 
                cout << ke << '\n';
            Ekinetic_ += ke;
        }


    if ( dof() )
        currentTemp_ = 2.0 * Ekinetic() / (dof() * kBoltzmann());
    else
        currentTemp_ = 0;
    // approxTemperature(); //FIX: temporary. remove after testing.
}
 
typedef SubVector<Vec6> RSubVec6;
typedef SubVector<const Vec6> ConstRSubVec6;

//
// Using current atom positions/velocities.
// This resets the CM vel and angular momentum. This operation is performed 
// in Cartesian space, and then the internal coordinates are recalculated.
// Thus, the linear vel. and AM will not generally be exactly zero after
// converting to internal coords and reconverting to Cartesians.
//
// do this in spatial coords. L = L_l + L_a
//
// correct base velocity and then recalc spatial and cartesian vels.
//
// Remove only the linear velocity component from each body. 
void
IVM::resetCM()
{
    if ( !resetCMflag ) return;
    resetCMflag=0;

    // exit if there are any fixed atoms: in this case, 
    //  removing CM velocity will change relative (internal) velocities
    if ( tree()->nodeTree[0][0]->getAtom(1) ) return;


    // Here we are going to misuse the IVMAtom class to turn each movable
    // AtomClusterNode into an "atom", with a single mass located
    // at the node's COM, having the node COM's (linear) velocity.
    static AtomList nodes;
    nodes.resize(0);
    for (int i=1 ; i<tree()->nodeTree.size() ; i++)
        for (int j=0 ; j<tree()->nodeTree[i].size() ; j++) {
            const AtomClusterNode* n  = tree()->nodeTree[i][j];
            const double acMass       = tree()->getClusterMass(i,j);
            const Vec3&  acCOM_G      = tree()->getClusterCOM_G(i,j);
            const Vec6&  acSpatialVel = tree()->getClusterSpatialVel(i,j);

            IVMAtom *a = new IVMAtom( nodes.size() , acMass );
            a->pos = acCOM_G;
            a->vel = Vec3( ConstRSubVec6(acSpatialVel,3,3).vector() )
                     - cross( ConstRSubVec6(acSpatialVel,0,3).vector() , 
                              n->getAtom(0)->pos - acCOM_G ); 
            nodes.append( a );
        }

    double mass=0.0;             //overall mass, posCM, velCM
    Vec3 posCM(0.0), velCM(0.0);
    for (l_int i=0 ; i<nodes.size() ; i++) {
        mass += nodes[i]->mass;
        posCM += nodes[i]->mass * nodes[i]->pos;
        velCM += nodes[i]->mass * nodes[i]->vel;
    }
    posCM /= mass;
    velCM /= mass;

    Vec3 L(0.0); // angular momentum excluding internal rotations
    for (int i=0 ; i<nodes.size() ; i++)
        L += nodes[i]->mass * cross( nodes[i]->pos - posCM , 
                                     nodes[i]->vel - velCM );

    if (verbose_&printResetCM) 
        cout << "resetCM: velCM: " << velCM << "overall AM: " << L << '\n';
    Vec3 omega(0.0); // the overall rotational velocity
    if ( nodes.size() >2 ) {
        InertiaTensor inertia;
        inertia.calc( posCM , nodes );
        try {
            omega = inverse(inertia) * L;
        } 
        catch ( CDS_NAMESPACE(SingularError) ) {
            cout << "InternalDynamics::resetCM: singular inertia Tensor.\n";
        }
    } else if ( nodes.size() == 2) {
        Vec3 r = nodes[1]->pos - nodes[0]->pos;
        omega = 1.0/abs2(r) * cross( r , nodes[1]->vel - nodes[0]->vel );
    }   

    // FIX: if a node lies on a principle axis
    //   add the appropriate component(s) to the total AM;

    // getSolver()->setPos( tree()->getPos() ); //FIX: ??
    // getSolver()->setVel( tree()->getVel() );

    for (int i=0; i<tree()->nodeTree[1].size(); i++) {
        AtomClusterNode* n = tree()->nodeTree[1][i];
        const Vec6& acSpatialVel = tree()->getClusterSpatialVel(1,i);

        Vec6 sVel = acSpatialVel 
                    - blockVec(omega, velCM + cross(omega, n->getAtom(0)->pos-posCM));

        tree()->setClusterVelFromSVel(1,i,sVel);
    }
    getSolver()->setPos( tree()->getPos() );
    getSolver()->setVel( tree()->getVel() );

    tree()->enforceConstraints(getSolver()->getPos(),getSolver()->getVel());
    tree()->setPosVel(getSolver()->getPos(),getSolver()->getVel());

    if (verbose_&printResetCM) { //REMOVE: for testing only
        Vec3 velCM(0.0);
        for (l_int i=1 ; i<tree()->nodeTree.size() ; i++)
            for (int j=0 ; j<tree()->nodeTree[i].size() ; j++) {
                const AtomClusterNode* n  = tree()->nodeTree[i][j];
                const double acMass       = tree()->getClusterMass(i,j);
                const Vec3&  acCOM_G      = tree()->getClusterCOM_G(i,j);
                const Vec6&  acSpatialVel = tree()->getClusterSpatialVel(i,j);

                Vec3 vCM = Vec3( ConstRSubVec6(acSpatialVel,3,3).vector() )
                                 - cross( ConstRSubVec6(acSpatialVel,0,3).vector() , 
                                          n->getAtom(0)->pos - acCOM_G ); 
                velCM += acMass * vCM;
            }
        velCM /= mass;

        Vec3 L(0.0); // angular momentum excluding internal rotations
        for (l_int i=1 ; i<tree()->nodeTree.size() ; i++)
            for (int j=0 ; j<tree()->nodeTree[i].size() ; j++) {
                const AtomClusterNode* n  = tree()->nodeTree[i][j];
                const double acMass       = tree()->getClusterMass(i,j);
                const Vec3&  acCOM_G      = tree()->getClusterCOM_G(i,j);
                const Vec6&  acSpatialVel = tree()->getClusterSpatialVel(i,j);

                Vec3 vCM = Vec3( ConstRSubVec6(acSpatialVel,3,3).vector() ) - 
                                cross( ConstRSubVec6(acSpatialVel,0,3).vector() , 
                                        n->getAtom(0)->pos - acCOM_G ); 
                L += acMass * cross( acCOM_G - posCM , vCM - velCM );
            }
        cout << "resetCM: after reset: velCM: " << velCM << "overall AM: " << L << '\n';
    }

    for (l_int i=0 ; i<nodes.size() ; i++) delete nodes[i]; //clean up.

    // --these next two lines are necessary if the pot. depends on absolute
    //   position or orientation
    // calcEnergy();        //recalc these- so eq. of motion are correct
    // calcAccel(*trees);
}

void
IVM::printCM()
{
    if ( verbose_&printCMVel ) {
        Vec3 posCM(0.0);
        double mass=0.0;;
        for (int i=1 ; i<atoms.size() ; i++) {
            IVMAtom* a = atoms[i];
            posCM += a->mass * a->pos;
            mass += a->mass;
        }
        posCM /= mass;

        Vec3 velCM(0.0);
        for (l_int i=1 ; i<atoms.size() ; i++)        //remove CM velocity
            velCM += atoms[i]->mass * atoms[i]->vel;
        velCM /= mass;

        cout.setf( ios::fixed | ios::right);

        cout << "COM vel: " << setw(6) << setprecision(3) << velCM << ' ';

        Vec3 L(0.0);
        for (l_int i=1 ; i<atoms.size() ; i++)
            L += atoms[i]->mass * cross( atoms[i]->pos-posCM , atoms[i]->vel );

        cout << "ang mom: " << setw(6) << L << '\n';
    }

    //FIX: this should be removed
    calcTemperature();

    if ( verbose_&printTemperature )
        cout << "temperature: " << currentTemp() << '\n';
    if ( verbose_&printEnergy )
        cout << "Energy: " << setprecision(14) << Etotal() 
             << " (" << Epotential() << ' ' << Ekinetic() << ')'
             << '\n';
    if ( verbose_&printCoords )
        cout << "pos: " << setprecision(14) << getSolver()->getPos() << '\n';
    if ( verbose_&(printNodeForce | printNodePos | printNodeTheta) )
        for (int i=1 ; i<tree()->nodeTree.size() ; i++)
            for (int j=0 ; j<tree()->nodeTree[i].size() ; j++)
                tree()->nodeTree[i][j]->print(verbose_);
}

// FIX: this is really xplor-specific
void
IVM::initAtoms(const int     natom,
               const double* massA,
               const double& kB)
{
    MALLOC_DEBUG;
    INSTALL_FPE_HANDLER;

    kBoltzmann_ = kB;
    for (int i=0 ; i<atoms.size() ; i++)
        delete atoms[i];
    atoms.resize(natom+1);
    atoms[0] = new IVMAtom(0,0.);   // Ground's origin 'atom' (station)
    for (l_int i=1 ; i<=natom ; i++)
        atoms[i] = new IVMAtom(i,massA[i-1]);
    cout.precision(8);
}

void
IVM::updateAccel() {
    tree()->updateAccel();
}

void
IVM::initTopology()
{
    // solver must be defined to determine hinge setup
    //  (whether to use euler angles or quaternions)
    if (getSolver()) delete solver_;
    solver_ = Solver::create(solverType,this);

    if (tree())
        delete tree_;
    tree_ = new AtomTree(this);

    if ( verbose_&printNodeDef ) {
        cout << "Node info:\n";
        cout << *tree() << '\n';
        cout << "total number of degrees of freedom: " << dof() << '\n';
    }

    dim_ = tree()->getDim();
}

void
IVM::initDynamics(bool reuseTopology)
{
    RVec pos, vel, acc;
    if ( !reuseTopology ) initTopology();

    pos = tree()->getPos();

    if ( minimization() ) {
        vel.resize( pos.size() );
        vel.set(0.0);
        tree()->enforceConstraints(pos,vel);    //FIX: probably redundant
        tree()->setPosVel(pos,vel);
    } else { // integration
        vel = tree()->getVel();
        tree()->velFromCartesian(pos,vel);
cout << "initDynamics: before enforceConstraints, pos=" << pos << endl;
        tree()->enforceConstraints(pos,vel);    //FIX: probably redundant
cout << "initDynamics: after enforceConstraints, pos=" << pos << endl;
        resetCM();
        tree()->setPosVel(pos,vel);
    }

    calcEnergy();
    if ( !minimization() )
        acc = tree()->getAccel();
    getSolver()->init(pos,vel,acc);

cout << "initDynamics: end vel=" << vel << endl;
}

String
IVM::idAtom(int id) const {
    StringStream ret;
    ret << id << ends;
    return ret.str();
}
