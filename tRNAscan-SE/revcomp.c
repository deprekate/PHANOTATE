/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* revcomp.c
 * 
 * Reverse complement of a IUPAC character string
 * 
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "squid.h"


#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


char *
revcomp(char *comp, char *seq)
{
  long  bases;
  char *bckp, *fwdp;
  int   idx;
  long  pos;
  int   c;

  if (comp == NULL) return NULL;
  if (seq == NULL)  return NULL;
  bases = strlen(seq);

  fwdp = comp;
  bckp = seq + bases -1;
  for (pos = 0; pos < bases; pos++)
    {
      c = *bckp;
      c = sre_toupper(c);
      for (idx = 0; c != iupac[idx].sym && idx < IUPACSYMNUM; idx++);
      if (idx == IUPACSYMNUM)
	{
	  Warn("Can't reverse complement an %c, pal. Using N.", c);
	  *fwdp = 'N';
	}
      else
	*fwdp = iupac[idx].symcomp;
      if (islower((int) *bckp)) *fwdp = sre_tolower((int) *fwdp);
      fwdp++;
      bckp--;
    }
  *fwdp = '\0';
  return comp;
}
  
