
#ifdef NOTDEF

//
// POSIX semaphore implementation
//

#include "semaphore.h"
#include <unistd.h>
#ifdef OSF                    /* ugly OSF hack */
#  define SEM_NAME_MAX 106
#endif
#include <semaphore.h>
#include <cdsExcept.h>
#include <cdsString.h>
#include <errno.h>
#include <string.h>
#include <mmapAlloc.h>

CDS_BEGIN_NAMESPACE

Semaphore::Semaphore(bool awake)
{
 using CDS::exception;
 String routineName = "Semaphore::Semaphore";

 sem = MMapAlloc::alloc(sizeof(sem_t));
 new (sem) sem_t;

 if ( sem_init((sem_t*)sem,1,awake) != 0 )
   throw exception(routineName + ": error in sem_init():\n\t" +
		   strerror(errno));
} /* constructor */

Semaphore::~Semaphore()
{
 ((sem_t*)sem)->~sem_t();
 MMapAlloc::free(sem);
}

void
Semaphore::wait()
{
 if ( sem_wait((sem_t*)sem) != 0 )
   throw exception(String("Semaphore::wait: error in sem_wait():\n\t") +
		   strerror(errno));
} /* wait */

void
Semaphore::wake()
{
 if ( sem_post((sem_t*)sem) != 0 )
   throw exception(String("Semaphore::wake: error in sem_wake():\n\t") +
		   strerror(errno));
} /* wake */

int
Semaphore::query()
{
 int ret;

 if ( sem_getvalue((sem_t*)sem,&ret) != 0 )
   throw exception(String("Semaphore::query: error in sem_getvalue():\n\t") +
		   strerror(errno));
 
 return ret;
} /* wake */

CDS_END_NAMESPACE

#endif
