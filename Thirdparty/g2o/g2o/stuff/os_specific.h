
#ifndef G2O_OS_SPECIFIC_HH_
#define G2O_OS_SPECIFIC_HH_

#ifdef WINDOWS
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#ifndef _WINDOWS
#include <sys/time.h>
#endif
#define drand48() ((double) rand()/(double)RAND_MAX)

#ifdef __cplusplus
extern "C" {
#endif

int vasprintf(char** strp, const char* fmt, va_list ap);

#ifdef __cplusplus
}
#endif

#endif

#ifdef UNIX
#include <sys/time.h>
#endif

#endif
