#ifdef NOTDEF


//
// System V semaphore implementation
//

#include "semaphore.h"

#include <cdsExcept.h>
#include <cdsString.h>
#include <cdsSStream.h>

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <errno.h>
#include <string.h>

#if _SEM_SEMUN_UNDEFINED || OSF || SUN
/* according to X/OPEN we have to define it ourselves */
union semun {
  int val;                    /* value for SETVAL */
  struct semid_ds *buf;       /* buffer for IPC_STAT, IPC_SET */
  unsigned short int *array;  /* array for GETALL, SETALL */
  struct seminfo *__buf;      /* buffer for IPC_INFO */
};
#endif

CDS_BEGIN_NAMESPACE

Semaphore::Semaphore(bool awake)
{
 using CDS::exception;
 String routineName = "Semaphore::Semaphore";

 sem = new int;

 *(int*)sem = semget(IPC_PRIVATE, 1, 0600 | IPC_CREAT);
 if ( *(int*)sem<0 )
   throw exception(routineName + ": error in semget():\n\t" +
		   strerror(errno));

 union semun un;
 un.val = awake;
 semctl(*(int*)sem,0,SETVAL,un);
 
// cout << "Semaphore::Semaphore: sem value: " 
//      << sem << ' ' << *(int*)sem << '\n';
} /* constructor */

Semaphore::~Semaphore() 
{
// cout << "Semaphore::~Semaphore: sem value: " 
//      << sem << ' ' << *(int*)sem << '\n';
 union semun un; un.val=0; // may give an undefined warning
 if ( semctl(*(int*)sem,0,IPC_RMID,un) == -1 ) {
   OStringStream mess;
   mess << "Semaphore::~Semaphore: error in semctl():\n\t"
	<<  strerror(errno) << " [sem value: " << sem 
	<< ' ' << *(int*)sem<< "]";
   cout << mess.str() << '\n';
   throw exception(mess.str());
 }
 delete (int*)sem;
} /* destructor */

void
Semaphore::wait()
{
// cout << "Semaphore::init: sem value: " 
//      << sem << ' ' << *(int*)sem << '\n';
 struct sembuf    op_op[1] = {
   { 0, -1, 0 }
 }; 

 // if semval==0, this blocks until semval is incremented, then
 // subtracts one.
 while ( semop(*(int*)sem, &op_op[0], 1) <0 ) 
   if ( errno != EINTR )
     throw exception(String("Semaphore::wait: error in semop():\n\t") +
		     strerror(errno));
} /* wait */

void
Semaphore::wake()
{
// cout << "Semaphore::wake: sem value: " 
//      << sem << ' ' << *(int*)sem << '\n';
 struct sembuf    op_op[1] = {
   { 0, 1, 0 } 
 }; 

 // adds one to semval and does not block
 //
 if ( semop(*(int*)sem, &op_op[0], 1) <0 )
   throw exception(String("Semaphore::wake: error in semop():\n\t") +
		   strerror(errno));
} /* wake */

CDS_END_NAMESPACE

#endif
