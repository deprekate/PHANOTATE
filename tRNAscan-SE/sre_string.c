/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* sre_string.c
 * 
 * my library of extra string functions. Some for portability
 * across UNIXes
 */

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "gnuregex.h"
#include "squid.h"


#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Obsolete. Anyone who doesn't have strstr() is
 * not ANSI-compliant and must die.
 */
#ifdef NOSTR
char  *strstr(char *s, char *subs)
{
  int  i;

  for ( ; *s != 0; s++) {
    if (*s == *subs) {
      for (i = 1; subs[i] != 0 && subs[i] == s[i]; i++) ;
      if (subs[i] == 0) return(s);
      }
    }
  return (NULL);
}
#endif /* NOSTR */


char *
Strdup(char *s)
{
  char *new;
  if ((new = (char *) malloc (strlen(s) +1)) == NULL) return NULL;
  strcpy(new, s);
  return new;
}

int
Strinsert(char  *s1,            /* string to insert a char into  */
	  char   c,		/* char to insert                */
	  int    pos)		/* position in s1 to insert c at */
{
  char    oldc;
  char   *s;

  for (s = s1 + pos; c; s++)
    {
				/* swap current char for inserted one */
      oldc = *s;		/* pick up current */
      *s   = c;   		/* put down inserted one    */
      c    = oldc;		/* old becomes next to insert */
    }
  *s = '\0';

  return 1;
}


int
Strdelete(char *s1,             /* string to delete a char from       */
	  int   pos)		/* position of char to delete 0..n-1  */
{
  char *s;                      

  for (s = s1 + pos; *s; s++)
    *s = *(s + 1);

  return 1;
}

void
s2lower(char *s)
{
  for (; *s != '\0'; s++)
    *s = sre_tolower((int) *s);
}

void
s2upper(char *s)
{
  for (; *s != '\0'; s++)
    *s = sre_toupper((int) *s);
}


void *
MallocOrDie(size_t size)
{
  void *ptr;
  if ((ptr = malloc (size)) == NULL)
    Die("malloc failed");
  return ptr;
}

void *
ReallocOrDie(void *p, size_t size)
{
  void *ptr;
  if ((ptr = realloc(p, size)) == NULL)
    Die("realloc failed");
  return ptr;
}


/* Function: Strparse()
 * 
 * Purpose:  Match a regexp to a string.
 *           Return 0 if it matches, REG_NOMATCH if it doesn't.
 *           The caller may request a copy of the text that matched by 
 *           passing a non-NULL pointer to a string pointer.
 *           The called may also request copies of the text that matched 
 *           sub-regexps (parenthesized expressions) by passing
 *           ntok > 0, and pointers to ntok different string pointers.
 * 
 *           Uses the GNU regexp library in extended POSIX compatibility
 *           mode.
 *           
 *           I built this for ease of use, not speed nor efficiency.
 *
 * Example:  Strparse("foo-...-baz", "foo-bar-baz", NULL, 0)  returns 0
 *           Strparse("foo-\(...\)-baz", "foo-bar-baz", &buf1, 1, &buf2) 
 *              returns 0, copies "foo-bar-baz" into buf1, and "bar" 
 *              into buf2.
 *              
 * Args:     rexp  - regular expression, extended POSIX form
 *           s     - string to match against
 *           buf   - if non-NULL, where to put copy of matching text
 *           ntok  - number of sub-regexps returned.
 *           ...   - variable number of pointers to strings to keep
 *                   copies of matched subtexts
 *                   
 * Return:   0 on match, REG_NOMATCH on failure to match
 *           buf and the ptrs in varargs list are malloc'ed here,
 *           must be free'd by caller.
 */
int
Strparse(char *rexp, char *s, char **buf, int ntok, ...)
{
  va_list     ap;
  regex_t     pat;
  int         code;
  regmatch_t *pmatch;
  char      **sp;
  int         len;
  int         i;

  code = regcomp(&pat, rexp, REG_EXTENDED);
  if (code > 0) {
    fprintf(stderr, "regular expression compilation failed\n");
    exit(1);
  }

  if ((pmatch = (regmatch_t *) malloc (sizeof(regmatch_t) * (ntok+1))) == NULL) {
    fprintf(stderr, "malloc failed\n");
    exit(1);
  }

  code = regexec(&pat, s, ntok+1, pmatch, 0);
  if (code == 0) {
				/* make copy of full matched text */
    if (buf != NULL) {
      len = pmatch[0].rm_eo - pmatch[0].rm_so;
      if ((*buf = (char *) malloc(sizeof(char) * (len+1))) == NULL) {
	fprintf(stderr, "malloc failed\n");
	exit(1);
      }
      strncpy(*buf, s+pmatch[0].rm_so, len);
      (*buf)[len] = '\0';
    }
				/* make copies of subtexts */
    if (ntok > 0) {
      va_start(ap, ntok);
      for (i = 1; i <= ntok; i++) {
	sp = va_arg(ap, char **);
	len = pmatch[i].rm_eo - pmatch[i].rm_so;
	if ((*sp = (char *) malloc(sizeof(char) * (len+1))) == NULL) {
	  fprintf(stderr, "malloc failed\n");
	  exit(1);
	}
	strncpy(*sp, s+pmatch[i].rm_so, len);
	(*sp)[len] = '\0';
      }
    }
  }

  va_end(ap);
  free(pmatch);
  regfree(&pat);
  return code;
}


/* Function: StrShuffle()
 * 
 * Purpose:  Returns a shuffled version of s2, in s1.
 *  
 * Args:     s1 - allocated space for shuffled string.
 *           s2 - string to shuffle.
 *           
 * Return:   void
 */
void
StrShuffle(char *s1, char *s2)
{
  int  len;
  int  pos;
  char c;
  
  strcpy(s1, s2);
  for (len = strlen(s1); len > 1; len--)
    {				
      pos       = CHOOSE(len);
      c         = s1[pos];
      s1[pos]   = s1[len-1];
      s1[len-1] = c;
    }
}
  
