/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* file: types.c
 * 
 * Finicky type checkers for strings. Return 1 (TRUE) if ok, 0 elsewise.
 * 
 */

#include <string.h>
#include <ctype.h>
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Function: IsInt()
 * 
 * Returns TRUE if s points to something that atoi() will parse
 * completely and convert to an integer.
 */
int
IsInt(char *s)
{
  int hex = 0;

  if (s == NULL) {squid_errno = SQERR_PARAMETER; return 0; }

				/* skip whitespace */
  while (isspace(*s)) s++;      
				/* skip leading sign */
  if (*s == '-' || *s == '+') s++;
				/* skip leading conversion signals */
  if ((strncmp(s, "0x", 2) == 0 && (int) strlen(s) > 2) ||
      (strncmp(s, "0X", 2) == 0 && (int) strlen(s) > 2))
    {
      s += 2;
      hex = 1;
    }
  else if (*s == '0' && (int) strlen(s) > 1)
    s++;
				/* examine remainder for garbage chars */
  if (!hex)
    while (*s != '\0')
      {
	if (!isdigit(*s)) return 0;
	s++;
      }
  else
    while (*s != '\0')
      {
	if (!isxdigit(*s)) return 0;
	s++;
      }

  return 1;
}


/* Function: IsReal()
 * 
 * Purpose:  Returns TRUE if s is a string representation
 *           of a valid floating point number.
 */
int
IsReal(char *s)
{
  int gotdecimal = 0;
  int gotexp     = 0;
  int gotreal    = 0;

  if (s == NULL) return 0;

  while (isspace(*s)) s++;         /* skip leading whitespace */
  if (*s == '-' || *s == '+') s++; /* skip leading sign */

  /* Examine remainder for garbage. Allowed one '.' and
   * one 'e' or 'E'; if both '.' and e/E occur, '.'
   * must be first.
   */
  while (*s != '\0')
    {
      if (isdigit(*s)) 
	gotreal++;
      else if (*s == '.')
	{
	  if (gotdecimal) return 0; /* can't have two */
	  if (gotexp) return 0;	/* e/E preceded . */
	  else gotdecimal++;
	}
      else if (*s == 'e' || *s == 'E')
	{
	  if (gotexp) return 0;	/* can't have two */
	  else gotexp++;
	}
      else if (isspace(*s))
	break;

      s++;
    }

  while (isspace(*s)) s++;         /* skip trailing whitespace */
  if (*s == '\0' && gotreal) return 1;
  else return 0;
}
