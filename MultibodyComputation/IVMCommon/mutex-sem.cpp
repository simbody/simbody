
/*
 * mutex-sem.cc
 *   a mutex class implemented using Semaphores
 *
 */


#include "mutex.h"

#include <cdsIostream.h>
#include <semaphore.h>
#include <errno.h>

CDS_BEGIN_NAMESPACE

Mutex::Mutex() {
 mx = new Semaphore(1);
} /* Mutex::Mutex */

Mutex::~Mutex() {
 delete (Semaphore*)mx;
} /* Mutex::~Mutex */


void 
Mutex::hold() 
{
 ((Semaphore*)mx)->wait();
} /* Mutex::hold */

void
Mutex::release() 
{
 ((Semaphore*)mx)->wake();
} /* Mutex::release */

CDS_END_NAMESPACE
