
#include "publicIVM.h"

#include <cdsIostream.h>
#include <cdsIomanip.h>

#include <sthead.h>
#include <cdsExcept.h>
#include <cdsSStream.h>
#include <dinternal.h>
#include <dint-step.h>
#include <dint-atom.h>
#include <dint-node.h>
#include <dint-loop.h>

#include <simulationWorld.h>
#include <simulation.h>
#include <derivList.h>
#include <math.h>

// Called when the simulation is modified outside of the IVM
// note that the modified flag is cleared when the simulation sets
// atom properties.
void 
PublicIVM::updateValues()
{
    if ( SimulationWorld::world()->logLevel() >= SimulationWorld::LOG_ALL ) 
        cout << "PublicIVM::updateValues: resyncing pos/vel from simulation\n";

    init();
}


//
// initialize InternalDynamics quantities from simulation values
//
void 
PublicIVM::init()
{
    sim->registerCallbacks(this);

    shadowPos = sim->atomPosArr(); //necessary only for space allocation
    shadowVel = sim->atomVelArr();

    for (int i=0 ; i<atoms.size() ; i++) //clean up any old mess
        delete atoms[i];
    atoms.resize(0);

    kBoltzmann_ = SimulationWorld::world()->kBoltzmann();

    atoms.resize(sim->numAtoms()+1);
    atoms[0] = new IVMAtom(0,0.);   //origin atom
    for (l_int i=1 ; i<=sim->numAtoms() ; i++) {
        atoms[i] = new IVMAtom(i,sim->atomByID(i-1));
    //   *atoms[i] = sim->atomByID(i-1); //FIX: remove repetition
        atoms[i]->fric *= frictionCoeff() * atoms[i]->mass 
                          * SimulationWorld::world()->timeFactor();
    }

    for (int bnum=0 ; bnum<sim->numBonds() ; bnum++)
        if ( !bondExclude().contains( sim->bondPairByID(bnum) ) ) {
            int a = sim->bondPairByID(bnum).a+1;
            int b = sim->bondPairByID(bnum).b+1;
            atoms[ a ]->bonds.append( atoms[ b ] );
            atoms[ b ]->bonds.append( atoms[ a ] );
        }

    loops.resize(0);
    initTopology();
    modified.clear();
}

void 
PublicIVM::initDynamics(bool reuseTopology)  { 
    modified.update();
    // modified.clear();
    IVM::initDynamics(reuseTopology); 
}

void
PublicIVM::syncPos() {
    for (int i=0 ; i<sim->numAtoms() ; i++)
        shadowPos[i] = atoms[i+1]->pos;
    sim->setAtomPosArr( shadowPos );
    modified.clear();
}

void
PublicIVM::syncVel() {
    for (int i=0 ; i<sim->numAtoms() ; i++)
        shadowVel[i] = atoms[i+1]->vel;
    sim->setAtomVelArr( shadowVel );
    modified.clear();
} 

void
PublicIVM::calcEnergy() {
    // modified not consulted!

    syncPos();  //synchronize pos/vel to simulation
    syncVel();

    static DerivList derivList;

    derivList.clear();
    Epotential_ = potList_.calcEnergyAndDerivs(derivList).energy;

    calcTemperature(); //calcs Ekinetic
    Etotal_ = Ekinetic() + Epotential();

    for (int i=0 ; i<sim->numAtoms() ; i++)
        atoms[i+1]->deriv = derivList[ sim ](i);
    gradMagnitude_ = sqrt( vecVec3Abs2(derivList[sim]) );
    if ( sim->numAtoms() )
        gradMagnitude_ /= sqrt( (double)(3*sim->numAtoms()) );

    eCount_++; 
}

void
PublicIVM::calcTemperature() {
    //modified not consulted!

    syncVel();
    Ekinetic_ = sim->kineticEnergy();

    // double tmp=0.;
    // for (int i=0 ; i<tree()->nodeTree.size() ; i++)
    //   for (int j=0 ; j<tree()->nodeTree[i].size() ; j++)
    //     tmp += tree()->nodeTree[i][j]->kineticE();

    if ( dof() )
        currentTemp_ = 2.0 * Ekinetic() / (dof() * kBoltzmann());
    else
        currentTemp_ = 0;
}



void
PublicIVM::groupTorsion() {
    init();
    IVM::groupTorsion();
}

bool
PublicIVM::step(double& stepsize)
{
    modified.update();

    bool done=0;
    if ( !minimization() ) 
        stepsize /= SimulationWorld::world()->timeFactor();
    try {
        IVM::step(stepsize);
    }
    catch ( Integrator::Finished ) {
        done=1;
    }
    if ( !minimization() ) 
        stepsize *= SimulationWorld::world()->timeFactor();

    //update simulation values
    syncPos();
    syncVel();

    return done;
}

//
// accessors: correct for offset-1 atoms
//

CDSList< CDSList<int> > 
PublicIVM::groupList() {
    CDSList< CDSList<int> > ret = IVM::groupList;
    for (int i=0 ; i<ret.size() ; i++)
        for (int j=0 ; j<ret[i].size() ; j++)
            ret[i][j]--;
    return ret;
}

CDSList<Pair> 
PublicIVM::constraintList() {
    CDSList<Pair> ret = IVM::constraintList; 
    for (int i=0 ; i<ret.size() ; i++) {
        ret[i].a--;
        ret[i].b--;
    }
    return ret;   
}

CDSList<int> 
PublicIVM::oldBaseAtoms() {
    CDSList<int> ret = IVM::oldBaseAtoms; 
    for (int i=0 ; i<ret.size() ; i++)
        ret[i]--;
    return ret;
}

CDSList<HingeSpec> 
PublicIVM::hingeList() { 
    CDSList<HingeSpec> ret = IVM::hingeList;
    for (int i=0 ; i<ret.size() ; i++)
        for (int j=0 ; j<ret[i].aList.size() ; j++) {
            ret[i].aList[j]--;
            ret[i].atom0--; // XXX (sherm) this looks wrong -- shouldn't be in inner loop?
            ret[i].atom1--;
        }
    return ret;
}

void 
PublicIVM::setGroupList(const CDSList< CDSList<int> >& list) {
    IVM::groupList = list;
    for (int i=0 ; i<list.size() ; i++)
        for (int j=0 ; j<list[i].size() ; j++)
            IVM::groupList[i][j]++;
}

void 
PublicIVM::setConstraintList(const CDSList<Pair>& list ) { 
    IVM::constraintList = list; 
    for (int i=0 ; i<list.size() ; i++) {
        IVM::constraintList[i].a++;
        IVM::constraintList[i].b++;
    }
}

//void
//PublicIVM::setBondExclude(const CDSList<Pair,1>& list ) 
//{ 
// IVM::bondExclude = list;
// for (int i=0 ; i<list.size() ; i++) {
//   IVM::bondExclude[i].a++;
//   IVM::bondExclude[i].b++;
// }
//} 

void 
PublicIVM::setOldBaseAtoms(const CDSList<int>& list ) { 
    IVM::oldBaseAtoms = list; 
    for (int i=0 ; i<list.size() ; i++)
        IVM::oldBaseAtoms[i]++;
}
 
void 
PublicIVM::setHingeList(const CDSList<HingeSpec>& list ) { 
    IVM::hingeList = list; 
    for (int i=0 ; i<list.size() ; i++)
        for (int j=0 ; j<list[i].aList.size() ; j++) {
            IVM::hingeList[i].aList[j]++;
            IVM::hingeList[i].atom0++;
            IVM::hingeList[i].atom1++;
        }
}

String
PublicIVM::idAtom(int t_id) const
{
    StringStream ret;
    int id = t_id-1;
    ret << id << ' ';
    if (id<-1)
        ret << "(unknown)";
    else if (id==-1)
        ret << "(base atom)";
    else if (id>=sim->numAtoms())
        ret << "(unknown)";
    else
        ret << sim->residueName(id) << ' ' << sim->residueNum(id)
            << ' ' << sim->atomName(id);
    ret << ends;
    return ret.str();
}

CDSVector<double>
PublicIVM::pos() const {
    modified.update();
    return tree()->getPos();
}


void
PublicIVM::setPos(const CDSVector<double>& newPos) {
    CDSVector<double> oldPos = tree()->getPos();

    if ( oldPos.size() != newPos.size() ) {
        cerr << "PublicIVM::setPos: size mismatch" << endl;
        throw InternalDynamics::Exception("PulbicIVM::setPos: size mismatch");
    }

    RVec pos = newPos;
    RVec vel;
    if ( minimization() ) {
        vel = pos;
        vel.set(0.);
    } else
        vel = tree()->getVel();
    tree()->enforceConstraints(pos,vel);
    tree()->setPosVel(pos,vel);
    syncPos();
}

void
PublicIVM::velFromCartesian() {
    modified.update();
    RVec vel = tree()->getVel();
    tree()->velFromCartesian(tree()->getPos(),vel);
    tree()->setVel(vel);
    // modified.clear(); //???
}

static int
toAtomIndex(const Simulation* sim,
            const IVMAtom*    atom) {
    if ( atom->index>0 && atom->index<=sim->numAtoms() )
        return atom->index-1;
    else
        return -1;
}

CDSList<PublicNode>
PublicIVM::nodeList()
{
    CDSList<PublicNode> ret;
    for (int i=0 ; i<tree()->nodeTree.size() ; i++)
        for (int j=0 ; j<tree()->nodeTree[i].size() ; j++) {
            HingeNode* n = tree()->nodeTree[i][j];
            PublicNode pn;
            pn.setDim( n->getDim() );
            pn.setStartIndex( n->offset()-1 ); //offsets differ
            pn.setType( n->type() );
            CDSList<int> tmpAtoms;
            for (int i=0 ; n->getAtom(i) ; i++)
                tmpAtoms.append( toAtomIndex(sim, n->getAtom(i)));
            pn.setAtoms( tmpAtoms );
            if ( n->parentAtom )
                pn.setParentAtom( toAtomIndex(sim,  n->parentAtom ) );
            else
                pn.setParentAtom( -1 );
            ret.append( pn );
        }
    return ret;
}

namespace EnumNamespace {

  EnumNameMap VerboseFlags[] = 
    {{ "printCoords"         , PublicIVM::printCoords          },
     { "printResetCM"        , PublicIVM::printResetCM         },
     { "printVelFromCartCost", PublicIVM::printVelFromCartCost },
     { "printTemperature"    , PublicIVM::printTemperature     },
     { "printEnergy"         , PublicIVM::printEnergy          },
     { "printCMVel"          , PublicIVM::printCMVel           },
     { "printNodeForce"      , PublicIVM::printNodeForce       },
     { "printNodePos"        , PublicIVM::printNodePos         },
     { "printNodeTheta"      , PublicIVM::printNodeTheta       },
     { "printStepDebug"      , PublicIVM::printStepDebug       },
     { "printStepInfo"       , PublicIVM::printStepInfo        },
     { "printNodeDef"        , PublicIVM::printNodeDef         },
     { "printLoopDebug"      , PublicIVM::printLoopDebug       },
     { "printLoopInfo"       , PublicIVM::printLoopInfo        },
     { "last enum"           , 0                               }};
};
