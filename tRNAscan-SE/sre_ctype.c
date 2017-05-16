/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* sre_ctype.c
 * 
 * For portability. Some systems have functions tolower, toupper
 * as macros (for instance, MIPS M-2000 RISC/os!)
 */

#include <ctype.h>
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

int
sre_tolower(int c)
{
  if (isupper(c)) return tolower(c);
  else return c;
}

int
sre_toupper(int c)
{
  if (islower(c)) return toupper(c);
  else return c;
}

