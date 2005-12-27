#ifdef NOTDEF

#include "ensembleSimulation.h"
#include <cdsSStream.h>
#include <cdsList.h>
#include <cdsVector.h>
#include <cdsIomanip.h>
#include <math.h>
#include <errno.h>
#include <sys/types.h>
//#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
//#include <unistd.h>
#include <signal.h>
#include <string.h>
#include <stdio.h>
#include <mutex.h>
#include <simulationWorld.h>
#include <cdsMap.h>

#include <publicIVM.h>

static bool debugThreads=0;


//id of singleton in use
EnsembleSimulation* EnsembleSimulation::currentSimulation_=0;

static CDSMap<int,EnsembleSimulation*> sizeOneSimulations;

typedef int    (*RVecSizeType)(const RVec&);
typedef double (*RVecProdType)(const RVec&, const RVec&);
typedef int    (*VecVec3SizeType)(const CDSVector<CDSVec3>&);
typedef double (*VecVec3ProdType)(const CDSVector<CDSVec3>&, 
				  const CDSVector<CDSVec3>&);

static RVecSizeType    rVecSize=0;
static RVecProdType    rVecProd=0;
static VecVec3SizeType vecVec3Size=0;
static VecVec3ProdType vecVec3Prod=0;

// the following singleton is to assure that the current simulation is
// shut down (and processes reaped).
// ideally, the EnsembleSimulation object would be reaped by the appropriate
// interpreter (Python, TCL, etc) upon shutdown. However, objects are not 
// always destroyed (due to cyclic dependencies, etc).
// CDS 5/12/04
class EnsembleSimulationShutdown {
public:
  EnsembleSimulationShutdown() {}
  ~EnsembleSimulationShutdown() {
   if ( EnsembleSimulation::currentSimulation() )
     EnsembleSimulation::currentSimulation()->shutdown();
  }
};
static EnsembleSimulationShutdown ensembleSimulationShutdown;

//
// memory allocation routines
//

struct AllocData {
  // note the definition of new/delete: this is allocated with MMap, and hence 
  // this must happen before a fork in order for this data structure to be 
  // shared.
  CDS::Mutex mutex;
  int visitCnt;
  int allocCnt;
  
  AllocData() : visitCnt(0), allocCnt(0) {}

  void* operator new(size_t s) { return MMapAlloc::alloc(s); }
  void  operator delete(void* p) { MMapAlloc::free(p); }
};

void*
EnsembleSimulation::SharedAlloc::alloc(size_t size)
{
 if ( !currentSimulation() ||
      currentSimulation()->size()==1 )
   return MMapAlloc::alloc(size);

 currentSimulation_->allocData->mutex.hold();

 // first to enter increments the allocCnt
 if ( !currentSimulation_->allocData->visitCnt )
   currentSimulation_->allocData->allocCnt++;

 currentSimulation_->allocData->visitCnt++;
 
// if (debugThreads)
//   cerr << "EnsembleSimulation::SharedAlloc::alloc: size: " << size
//	<< " visitCnt: " << currentSimulation_->allocData->visitCnt
//	<< " allocCnt: " << currentSimulation_->allocData->allocCnt
//	<< endl;

 bool deleteFile=0;
 if ( currentSimulation_->allocData->visitCnt==currentSimulation_->size() ) {
   currentSimulation_->allocData->visitCnt=0;
   deleteFile=1;
 }

 void* ret = MMapAlloc::alloc(size,deleteFile,
			      currentSimulation_->allocData->allocCnt);
 
 currentSimulation_->allocData->mutex.release();

 currentSimulation_->barrier(__FILE__,__LINE__);
 return ret;
} /* SharedAlloc::alloc */

void
EnsembleSimulation::SharedAlloc::free(void* p)
{
 MMapAlloc::free(p);
} /* SharedAlloc::free */


//
// distributed inner product
//

struct DistributedProduct {
  EnsembleSimulation* sim;
  int                 sharedSize;
  int                 resultSize;
  CDSVector<double,0,MMapAlloc> sharedSum;
  float_type          result;
  int                 sharedCnt;
  CDS::Mutex          mutex;

  void* operator new(size_t s) { return MMapAlloc::alloc(s); }
  void  operator delete(void* p) { MMapAlloc::free(p); }

  DistributedProduct(EnsembleSimulation* ensembleSim,int size) : 
    sim(ensembleSim), sharedSize(0), sharedCnt(0) 
  {
   sharedSum.resize( size );
  }

  template<class VEC>
  int vecSize(const VEC& v);

  template<class VEC>
  float_type vecProd(const VEC& v1, const VEC& v2);
  
};

template<class VEC>
inline int
DistributedProduct::vecSize(const VEC& v)
{
 int size = v.size();
 
 { //protected region
   mutex.hold();

   sharedSize += size;
   sharedCnt++;
   
   if ( debugThreads )
     cerr << "DistributedProduct::vecSize: " << sharedSize << ' ' 
	  << sharedCnt << endl;
   
   if (sharedCnt==sim->size() ) {
     resultSize = sharedSize;
     sharedSize = 0;
     sharedCnt = 0;
   }
   mutex.release();
 }


 sim->barrier(__FILE__,__LINE__);

 //note: result can't change here, because it only changes at end of next
 //      product calc.
 //
 //

 return resultSize;
} /* DistributedProduct::vecSize */

template<class VEC>
inline float_type
DistributedProduct::vecProd(const VEC& v1,
			    const VEC& v2)
{
 int index = sim->member()->memberIndex();
 sharedSum[index] = dot(v1,v2); 
// sharedSum[index] = 0.;
// for (int i=v1.offset() ; i<v1.size()+v1.offset() ; i++)
//   sharedSum[index] += v1[i] * v2[i];

 sim->barrier(__FILE__,__LINE__);
 
 float_type result = 0.;
 for (int i=0 ; i<sim->size() ; i++)
   result += sharedSum[i];

 sim->barrier(__FILE__,__LINE__);

 return result;
} /* DistributedProduct::vecProd */



template<class VEC>
static int
vectorSize(const VEC& v)
{
 EnsembleSimulation* sim = EnsembleSimulation::currentSimulation();
 return sim->distributedProduct->vecSize(v);
} /* vectorSize */

template<class VEC>
static double
vectorProduct(const VEC& v1,
	      const VEC& v2)
{
 EnsembleSimulation* sim = EnsembleSimulation::currentSimulation();
 return sim->distributedProduct->vecProd(v1,v2);
} /* vectorProduct */

static void
signal_childDied(int sig)
{
 cerr << "child exited unexpectedly.\n";
 abort();
} /* signal_childDied */

static void (*oldTermSignal)(int)=0;
static void
signalTerm_childKill(int sig)
{
 cerr << "killing children...\n";
 EnsembleSimulation::currentSimulation()->shutdown();
 (*oldTermSignal)(sig);
} /* signalTerm_childKill */

static void (*oldIntSignal)(int)=0;
static void
signalInt_childKill(int sig)
{
 cerr << "killing children...\n";
 EnsembleSimulation::currentSimulation()->shutdown();
 (*oldIntSignal)(sig);
} /* signalInt_childKill */

EnsembleSimulation::EnsembleSimulation(const char*       name,
					     int         size,
					     int         numThreads,
					     Simulation* sim) :
  Simulation(name), subSim_(sim), size_(size), numThreads_(numThreads),
  deleted(0),
  distributedProduct( (size>1)?new DistributedProduct(this,size):0 ),
  allocData( (size>1)?new AllocData():0 )
{

 if ( size>1 ) {
   if ( currentSimulation_ )
     throw CDS::exception(String("EnsembleSimulation: \n\t") +
			  "ensemble already created. Cannot create " +
			  "ensemble of ensembles");
   Simulation::makeCurrent(this);
 } else {
   if ( sizeOneSimulations.keys().contains( sim->id() ) )
     throw CDS::exception(String("EnsembleSimulation: \n\t") +
			  "attempting to create a second simulation based on "+
			  sim->name());
   id_ = sim->id();
   sizeOneSimulations[ id() ] = this;
 }


 if (numThreads_<0) {
   char* envVal = getenv("NUM_THREADS");
   if (envVal) {
     IStringStream str(envVal);
     str >> numThreads_;
   }
   if (numThreads_<0)
     numThreads_=1;
 }
   
   
 if ( numThreads_>size ) 
   throw CDS::exception(String("EnsembleSimulation: numThreads cannot be ") +
			"larger than size");

 members_.resize(size);
 weights_.resize(size);
 for (int i=0 ; i<size ; i++)
   weights_[i] = 1.0 / size;

 //signal(SIGALRM,signal_continue);
 // create size() shared regions, w/ semaphores
 for (int i=0 ; i<size ; i++) {
   OStringStream tmpName; 
   tmpName << name << '-' << i << ends;
   members_[i].reset( new EnsembleMemberSimulation(tmpName.str(),this,i) );
 }

 //fork
 members(0)->memberData_->pid = getpid();
 member_ = members(0);
 // create new processes which start off asleep

 // output redirection:
 //  child process output can go to one of three places
 //  1) stdout 
 //     set THREAD_FILENAME to /dev/tty
 //  2) /dev/null
 //     delete THREAD_FILENAME
 //  3) a separate file for each process
 //     specify THREAD_FILENAME: owner, parent process number, ensemble
 //     member numebr appended.

 const char* outputNamePrefix = getenv("THREAD_FILENAME");
 OStringStream outputName;
 if ( String(outputNamePrefix) == "/dev/tty" ) {
   ; // nothing
 } else if ( outputNamePrefix ) {
   char* logName = getlogin();
   outputName << outputNamePrefix << '-' << (logName?logName:"unknown")
	      << '-' << getpid() << "-";
 } else
   outputName << "/dev/null" << ends;


 for (int i=1 ; i<size ; i++) {
   pid_t pid = fork();
   if (pid==0) { //child
     member_ = members(i);
     if ( outputName.str()!= "" ) {
       if (outputName.str() != "/dev/null")
	 outputName << i << ends;
       freopen(outputName.str(),"aw",stdout);
       freopen(outputName.str(),"aw",stderr);
     }
     //freopen("/dev/tty","r",stdin);
     //fdopen(dup(0), "r");
     //setbuf(stdin,0);
     member()->sleep();
     //if (debugThreads)
     //  sleep(30);
     break;
   } else { // parent
     members(i)->memberData_->pid = pid;
     if (debugThreads)
       cerr << "construct: spawned " << pid << '\n';
   }
 }

 SimulationWorld* world = SimulationWorld::world();
 // each child get uniq random number seed
 world->setRandomSeed( world->random.seed()+member()->memberIndex() );

 if (debugThreads)
   cerr << "construct: " << member()->memberIndex()<< ' ' 
	<< getpid() << endl;

 // start up non-master threads
 if ( member()->memberIndex()==0 ) 
   for (int i=1 ; i<this->numThreads() ; i++)
     members(i)->wake();

 barrier(__FILE__,__LINE__);
 if ( size>1 ) {
   currentSimulation_=this;
   rVecSize    = vectorSize<RVec>;
   rVecProd    = vectorProduct<RVec>;
   vecVec3Size = vectorSize< CDSVector<CDSVec3> >;
   vecVec3Prod = vectorProduct< CDSVector<CDSVec3> >;
   // enable SIGCHLD
   signal(SIGCHLD,signal_childDied);
   signal(SIGABRT,SIG_DFL);
   oldTermSignal = signal(SIGTERM,signalTerm_childKill);
   oldIntSignal  = signal(SIGINT,signalInt_childKill);
 } 
 resize(); // this relies on currentSimulation being set.
 barrier(__FILE__,__LINE__);
} /* constructor */

EnsembleSimulation*
EnsembleSimulation::sizeOneSimulation(Simulation* sim)
{
 if ( ! sim )
   sim = Simulation::currentSimulation();

 if ( sim->type() == String("EnsembleSimulation") )
   sim = ((EnsembleSimulation*)sim)->subSim();

 if ( ! sizeOneSimulations.keys().contains( sim->id() ) ) {
   //temporary variable to work around bug in old versions of gcc (2.96)
   EnsembleSimulation* esim =
     new EnsembleSimulation(String("SizeOneEnsemble_")+sim->name(),1,1,sim);
   sizeOneSimulations[ sim->id() ] = esim;
 }

 return sizeOneSimulations[ sim->id() ];
} /* sizeOneEnsembleSimulation */

void
EnsembleSimulation::shutdown()
{
 if (debugThreads)
   cerr << "EnsembleSimulation::shutdown: entry\n";

 if (deleted) return;
 deleted=1;

 if ( size()==1 ) {
   //cleanup 
   CDSList<int> keys = sizeOneSimulations.keys();
   for (int i=0 ; i<keys.size() ; i++)
     if ( sizeOneSimulations[keys[i]] == this )
       sizeOneSimulations.remove(keys[i]);
 }

 // disable SIGCHLD
 if ( size()>1 )
   signal(SIGCHLD,SIG_DFL);


 if (debugThreads)
   cerr << "EnsembleSimulation::shutdown: sigchld disabled.\n";
 barrier(__FILE__,__LINE__);

 // pause non-master shepards
 if ( member()->memberIndex()>0           &&
      member()->memberIndex()<numThreads() )
   member()->sleep();

 Simulation::makeCurrent( subSim() );

 member()->memberData_->barrierID=0;
 
 if (debugThreads)
   cerr << "in destructor: " << member()->memberIndex()<< ' ' 
	<< getpid() << endl;

 if ( member()->memberIndex() )
   exit(0);
 
 for (int i=1 ; i<size() ; i++) {
   if (debugThreads)
     cerr << "stopping process " << members(i)->pid() << endl;
   members(i)->wake();
   int status=1;
   while ( 0==waitpid(members(i)->pid(),&status,WNOHANG) ) {
     if (debugThreads)
       cerr << "waiting for " << members(i)->pid() << "..." << endl;
     sleep(1);
   }
   //if ( ! WIFEXITED(status) )
   //  cerr << "EnsembleSimulation::shutdown"
   //	  << ": waitpid() returned abnormal status: " << status <<"\n";
 }
 if ( size()>1 ) {
   currentSimulation_=0;
   rVecSize    = 0;
   rVecProd    = 0;
   vecVec3Size = 0;
   vecVec3Prod = 0;
 }
} /* shutdown */

EnsembleSimulation::~EnsembleSimulation()
{
 shutdown();
} /* destructor */


void
EnsembleSimulation::barrier(const char* module,
				  int   id)
  // barrier for ensemble member processes
  // ensemble members run consecutively, in order of member->number
  // run one process group at a time.
  // Process group: group of numThreads processes
{
 member()->barrierIncr();
 member()->memberData_->barrierID = id;

 if ( size()==1 )
   return;

 for ( int i=0 ; i<size() ; i++) 
   if ( members(i)->barrierID()==-1 ) {
     OStringStream msg;
     msg << "EnsembleSimulation::barrier: \n\t"
	 << "process " << member()->memberIndex()
	 << "(" << member()->pid() << ")  exited unexpectedly\n"
	 << "\tat " << module << ':' << id << ends;
     
     cerr << msg;
     throw CDS::exception(msg.str());
   }     

 // sanity check:
 //  make sure that all processes assigned to this thread
 //  have already reached this barrier.
 for ( int i=member()->memberIndex() % numThreads() ;
       i<member()->memberIndex()                    ; 
       i+=numThreads()                              )
   if ( members(i)->barrierCnt() != member()->barrierCnt() ||
	members(i)->barrierID()  != member()->barrierID()    ) {
     OStringStream msg;
     msg << "EnsembleSimulation::barrier: \n\t"
	 << "process " << member()->memberIndex()
	 << "(" << member()->pid() << ")  lost synch with process "
	 << i << '(' << members(i)->pid() << ")\n" 
	 << "\tat " << module << ':' << id << ends;
     cerr << msg.str();
     cerr << "process count  ID:\n";
     for (int j=0 ; j<size() ; j++)
       cerr << j << ' ' << members(j)->barrierCnt() 
	    << ' ' << members(j)->barrierID() << '\n';
     abort();
     //throw CDS::exception(msg.str());
   }

 int nextMember = member()->memberIndex()+numThreads();

 if (debugThreads)
   cerr << "barrier: " << member()->memberIndex()<< ' ' 
	<< getpid() << ' ' << nextMember << ' ' 
	<< member()->barrierCnt() 
	<< "  [at " << module << ':' << id << ']' << endl;

 // this part is tricky:
 // all threads must be finished before returning from this barrier
 if ( nextMember >= size() ) {
   bool wait=0;
   do {      // expensive spin
     if ( ( member()->memberIndex()==0 ) &&
	  waitpid((pid_t)-1,0,WNOHANG) ) {
       OStringStream msg;
       msg << "EnsembleSimulation::barrier: \n\t"
	   << "a child process has terminated unexpectedly\n"
	   << "\tat " << module << ':' << id << ends;
       cerr << msg.str();
       throw CDS::exception(msg.str());
     }
     // FIX: child should die if parent is dead.
     //if ( ( member()->memberIndex()>0 ) &&
     //	    waitpid(members(0)->pid(),0,WNOHANG) ) {
     //	 OStringStream msg;
     //	 msg << "EnsembleSimulation::barrier: \n\t"
     //	     << "parent process has terminated unexpectedly\n";
     //	 cerr << msg;
     //	 throw CDS::exception(msg.str());
     //}
     wait=0;
     for (int i=size()-numThreads() ; i<size() ; i++)
       if ( members(i)->barrierCnt() < member()->barrierCnt()  )
	 wait=1;
   } while ( wait );
   nextMember %= numThreads();
 }

 members( nextMember )->wake();
 member()->sleep();

} /* barrier */




void
EnsembleSimulation::registerCallbacks(PublicIVM* ivm)
{
 if ( size()==1 )
   return;

 ivm->setRVecSize( rVecSize );
 ivm->setRVecProd( rVecProd );
 ivm->setVecVec3Size( vecVec3Size );
 ivm->setVecVec3Prod( vecVec3Prod );
} /* registerCallbacks */


float_type
EnsembleSimulation::kineticEnergy()
{
 
 //FIX: Is it possible to get by with only a single barrier?
 //     NOTE THE FUNKY DEFITION of KE!!!
 float_type ret = 0.;
 member()->calcKE();

 barrier(__FILE__,__LINE__);

 for (int i=0 ; i<size() ; i++)
   ret += weight(i) * members(i)->getKE();

 barrier(__FILE__,__LINE__);
 return ret;
} /* kineticEnergy */


CDSVector<CDSVec3> 
EnsembleSimulation::meanAtomPosArr() const
{ 
 CDSVector<CDSVec3> ret = members(0)->atomPosArr();
 for (int i=0 ; i<numAtoms() ; i++)
   ret[i] *= weight(0);

 for (int index=1 ; index<size() ; index++)
   for (int i=0 ; i<numAtoms() ; i++)
     ret[i] += weight(index) * members(index)->atomPos(i);

 return ret;
} /* meanAtomPosArr */


void
EnsembleSimulation::setAtomPosArr(const CDSVector<CDSVec3>& list)
{ 
 member()->setAtomPosArr(list); 
 markAsModified();
 barrier(__FILE__,__LINE__);
} /* setAtomPosArr */

void
EnsembleSimulation::setAtomVelArr(const CDSVector<CDSVec3>& list)
{ 
 member()->setAtomVelArr(list); 
 markAsModified();
 barrier(__FILE__,__LINE__);
} /* setAtomVelArr */

void
EnsembleSimulation::sync()
{ 
 barrier(__FILE__,__LINE__);
 if ( subSim()->numAtoms() != numAtoms() ) {
   resize();
   barrier(__FILE__,__LINE__);
 } else {
   setAtomPosArr( subSim()->atomPosArr() );
   setAtomVelArr( subSim()->atomVelArr() );
 }
} /* sync */

void
EnsembleSimulation::resize()
{

 numAtoms_ = subSim()->numAtoms();
 numBonds_ = subSim()->numBonds();

 for (int i=0 ; i<size() ; i++)
   members(i)->resize();

} /* resize */





void*
EnsembleMemberSimulation::MemberData::operator new(size_t s)
{
 return MMapAlloc::alloc(s);
}

void
EnsembleMemberSimulation::MemberData::operator delete(void* p)
{
 MMapAlloc::free(p);
}

EnsembleMemberSimulation::EnsembleMemberSimulation(const char*               n,
							 EnsembleSimulation* s,
						   const int                 i)
  : Simulation(n), ensembleSim_(s), memberIndex_(i), 
    sem( (s->size()>1)?new CDS::Semaphore : 0 )
{
 using CDS::exception;
 String routineName = "EnsembleMemberSimulation::EnsembleMemberSimulation";

 memberData_ = new MemberData;

 if ( s->size()==1 )
   return;



 // now, initialize the elements of memberData
 memberData_->pid = -2;
 memberData_->barrierCnt = 0;
 memberData_->barrierID  = 1;

} /* member constructor */

EnsembleMemberSimulation::~EnsembleMemberSimulation() {}

void 
EnsembleMemberSimulation::sleep()
{

 if (debugThreads)
   cerr << "sleep: entry: " << memberIndex()<< ' ' << getpid() 
	<< endl;
 //while ( memberData_->busy) {}

 sem->wait();
} /* sleep */

void
EnsembleMemberSimulation::wake()
{
 sem->wake();
} /* wake */


void
EnsembleMemberSimulation::resize()
{

 numAtoms_ = subSim()->numAtoms();
 numBonds_ = subSim()->numBonds();


 atomPosVec_.resize(numAtoms());
 atomVelVec_.resize(numAtoms());
 if ( ensembleSim()->member() == this ) {
   for (int i=0 ; i< numAtoms() ; i++)
     atomPosVec_(i) = subSim()->atomPos(i);
   for (int i=0 ; i< numAtoms() ; i++)
     atomVelVec_(i) = subSim()->atomVel(i);
 }

} /* resize */

CDSVector<CDSVec3> 
EnsembleMemberSimulation::atomPosArr() const
{ 
 if ( ensembleSim()->size()==1 )
   return subSim()->atomPosArr();

 CDSVector<CDSVec3> ret( numAtoms() );

 for (int i=0 ; i<numAtoms() ; i++)
     ret[i] = atomPos(i);  //FIX: this is lazy/slow

 return ret;
} /* atomPosArr */

CDSVector<CDSVec3> 
EnsembleMemberSimulation::atomVelArr() const
{ 
 if ( ensembleSim()->size()==1 )
   return subSim()->atomVelArr();

 CDSVector<CDSVec3> ret( numAtoms() );

 for (int i=0 ; i<numAtoms() ; i++)
     ret[i] = atomVel(i);  //FIX: this is lazy/slow

 return ret;
} /* atomVelArr */

void 
EnsembleMemberSimulation::setAtomPosArr(const CDSVector<CDSVec3>& v)
{
 if ( ensembleSim()->size()==1 ) {
   subSim()->setAtomPosArr(v);
   return;
 }

 for (int i=0 ; i<numAtoms() ; i++) // very lazy
   setAtomPos(i,v[i]);
 markAsModified();
} /* setAtomPosArr */

void 
EnsembleMemberSimulation::setAtomVelArr(const CDSVector<CDSVec3>& v)
{
 if ( ensembleSim()->size()==1 ) {
   subSim()->setAtomVelArr(v);
   return;
 }

 for (int i=0 ; i<numAtoms() ; i++) // very lazy
   setAtomVel(i,v[i]);
 markAsModified();
} /* setAtomVelArr */

void 
EnsembleMemberSimulation::setAtomPos(      int i, 
				     const CDSVec3& v)
{
 if ( ensembleSim()->size()==1 ) {
   subSim()->setAtomPos(i,v);
   return;
 }

 atomPosVec_(i) = v;
 if ( ensembleSim()->member() == this )
   subSim()->setAtomPos(i,v);
 markAsModified();
}

void 
EnsembleMemberSimulation::setAtomVel(      int i, 
				     const CDSVec3& v)
{
 if ( ensembleSim()->size()==1 ) {
   subSim()->setAtomVel(i,v);
   return;
 }

 atomVelVec_(i) = v;
 if ( ensembleSim()->member() == this )
   subSim()->setAtomVel( i, v );
 markAsModified();
}

void
EnsembleMemberSimulation::setAtomMass(      int         i, 
				      const float_type& v)
{
 subSim()->setAtomMass(i,v);
}

void
EnsembleMemberSimulation::setAtomFric(      int         i, 
				      const float_type& v)
{
 subSim()->setAtomFric(i,v);
}

void
EnsembleMemberSimulation::setAtomCharge(      int         i, 
					const float_type& v)
{
 subSim()->setAtomCharge(i,v);
}

void
EnsembleMemberSimulation::setSegmentName(      int         i, 
					 const char* v)
{
 subSim()->setSegmentName(i,v);
}

void
EnsembleMemberSimulation::setResidueName(      int         i, 
					 const char* v)
{
 subSim()->setResidueName(i,v);
}

void
EnsembleMemberSimulation::setResidueNum(      int i, 
					const int v)
{
 subSim()->setResidueNum(i,v);
}

void
EnsembleMemberSimulation::setAtomName(      int   i, 
				      const char* v)
{
 subSim()->setAtomName(i,v);
}

void
EnsembleMemberSimulation::setChemType(      int   i, 
				      const char* v)
{
 subSim()->setChemType(i,v);
}


float_type
EnsembleMemberSimulation::calcKE()
{
 if ( ensembleSim()->size() != 1     &&
      ensembleSim()->member() != this )
   subSim()->setAtomVelArr( atomVelArr() );

 memberData_->ke = subSim()->kineticEnergy();
 
 if ( ensembleSim()->size() != 1     &&
      ensembleSim()->member() != this )
   setAtomVelArr( subSim()->atomVelArr() );
 return getKE();
} /* kineticEnergy */



EnsembleRMSD::EnsembleRMSD(Simulation* sim) :
  rmsd_(0.), esim(0)
{
 if ( !sim )
   sim = Simulation::currentSimulation();

 if ( sim->type() != "EnsembleSimulation" )
   return;

 esim = (EnsembleSimulation*)sim;
 
 if ( esim->size() < 2 )
   esim = 0;

} /* constructor */

void
EnsembleRMSD::run(Simulation* sim,
		  int         index)
{
 if ( !esim )
   return;

 // FIX: should check that sim is esim
 Atom atom(esim,index);

 int resid = atom.residueNum();
   
 if ( !count.exists(resid) ) {
   count[resid] = 0;
   byResidue_[resid] = 0.;
 }

 count[resid]++;
   
 
 for (int index1=0 ; index1<esim->size() ; index1++) {
   Simulation* member1 = esim->members( index1 );
   for (int index2=index1+1 ; index2<esim->size() ; index2++) {
     Simulation* member2 = esim->members( index2 );
     CDSVec3 vec1 = member1->atomPos( atom.index() );
     CDSVec3 vec2 = member2->atomPos( atom.index() );
     double val = abs2( vec1-vec2 );
     rmsd_ += val;
     byResidue_[resid] += val;
   }
 }

} /* EnsembleRMSD::run */


float_type 
EnsembleRMSD::rmsd() const 
{
 double ret = 0.;

 if ( !esim )
   return ret;

 int totalNum = 0;
 for (int i=0 ; i<count.keys().size() ; i++) {
   int resid = count.keys()[i];
   totalNum += count[resid];
   ret += byResidue_[resid];
 }

 ret = sqrt( 2*ret / (esim->size()*(esim->size()-1)*totalNum) );
 
 return ret;
} /* EnsembleRMSD::rmsd */

CDSMap<int, float_type>
EnsembleRMSD::byResidue() const 
{
 CDSMap<int, float_type> ret;

 if ( !esim )
   return ret;

 for (int i=0 ; i<count.keys().size() ; i++) {
   int resid = count.keys()[i];
   ret[resid] = sqrt( 2*byResidue_[resid] / 
		      (esim->size()*(esim->size()-1)*count[resid]) );
 }

 return ret;
} /* EnsembleRMSD::byResidue */

#endif
