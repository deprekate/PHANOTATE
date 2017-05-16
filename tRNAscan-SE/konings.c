/* konings.c
 * 1.0: SRE, Wed Jul  7 08:53:37 1993
 * adapted for 2.0: SRE, Thu Sep  9 13:38:13 1993
 * 
 * Representation of secondary structure and secondary structural 
 * alignments using Danielle Konings' string notation, and Hogeweg's
 * mountain notation.
 * 
 * See: Konings and Hogeweg, J. Mol. Biol. 207:597-614 1989
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "structs.h"
#include "funcs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Function: Align2kh()
 * 
 * Purpose:  Convert an alignment (align_s linked list) into
 *           a secondary structure string. The symbols used
 *           are > and < for the two sides of a MATP; '.' 
 *           for other symbols.
 * 
 *           Also may return the "aligned" sequence; this 
 *           sequence is only aligned in the sense that deleted
 *           consensus positions are represented as '.'
 *           and insert positions are lower case.
 *
 *           Either ret_aseq or ret_khseq may be passed as NULL
 *           if you don't want one of them.
 */
int
Align2kh(struct align_s  *ali,
	 char           **ret_aseq,
	 char           **ret_khseq)
{
  struct align_s *curr;         /* ptr into alignment list     */
  char           *aseq;         /* RETURN: aligned seq         */
  char           *khseq;        /* RETURN: structure string    */
  int             len;		/* length of alignment, khseq  */
  int             pos;		/* position in khseq           */

  
  /* Count the length of the list and malloc for khseq.
   */
  len = 0;
  for (curr = ali->nxt; curr->nxt != NULL; curr = curr->nxt)
    len++;
  if ((khseq = (char *) malloc ((len+1) * sizeof(char))) == NULL ||
      (aseq  = (char *) malloc ((len+1) * sizeof(char))) == NULL)
    Die("Memory allocation failed, line %d of %s", __LINE__, __FILE__);

  /* This used to be harder. Now align_s already has a field for
   * secondary structure annotation, and we just have to copy it.
   */
  pos  = 0;
  for (curr = ali->nxt; curr->nxt != NULL; curr = curr->nxt)
    {
      switch (curr->type) 
	{
	case uBEGIN_ST:
	case uBIFURC_ST:
	  break;		/* neither should appear in an align_s! */

	case uDEL_ST:  
	  khseq[pos] = ' ';       
	  aseq[pos] = '-';                    
	  break;

	case uINSL_ST: 
        case uINSR_ST: 
	  khseq[pos] = curr->ss;  
	  aseq[pos] = sre_tolower(curr->sym);
	  break;

        case uMATL_ST: 
        case uMATR_ST: 
	case uMATP_ST: 
	  khseq[pos] = curr->ss;  
	  aseq[pos] = sre_toupper(curr->sym); 
	  break;

	default:
	  Die("unrecognized state type %d in Align2kh()", curr->type);
	}
      pos++;
    }
      
  khseq[pos] = '\0';
  aseq[pos]  = '\0';

  if (ret_khseq == NULL) free(khseq); else *ret_khseq = khseq;
  if (ret_aseq  == NULL) free(aseq);  else *ret_aseq  = aseq;
  return 1;
}




/* Function: PrintAliLandscape()
 * 
 * Purpose:  Print an alignment of sequence to model in the form 
 *           of a Konings/Hogeweg "mountain".
 */
int
PrintAliLandscape(FILE           *fp, 
		  struct cm_s    *cm,
		  struct align_s *ali)
{
  int   altitude;		/* current height on the mountain */
  struct align_s *curr;         /* ptr to current ali element     */
  struct align_s *prev;         /* ptr to previous ali element    */
  int   i;

  altitude = 0;
  prev     = NULL;
  for (curr = ali->nxt; curr->nxt != NULL; curr = curr->nxt)
    {
      if (curr->pos >= 0)
	fprintf(fp, "%4d  %c ", curr->pos+1, curr->sym);
      else
	fprintf(fp, "      %c ", curr->sym);

      for (i = 0; i < altitude; i++)
	fputc(' ', fp);
      
      switch (curr->type) 
	{
	case uBEGIN_ST:
	case uBIFURC_ST: break;
	case uDEL_ST:    fputs("   DEL ", fp);  break;
	case uINSL_ST:   fputs("`  INSL", fp);  break;
        case uINSR_ST:   fputs("\'  INSR", fp); break;
	case uMATL_ST:   fputs("\\  MATL", fp); break;
	case uMATR_ST:   fputs("/  MATR", fp);  break;

	case uMATP_ST: 
	  if (prev == NULL || curr->nodeidx > prev->nodeidx)
	    {
	      fputs(" v  MATP", fp);
	      altitude++;
	    }
	  else
	    {
	      fputs("^  MATP", fp);
	      altitude--;
	    }
	  break;

	  
	default:
	  Die("unrecognized state type %d in PrintAliLandscape()", curr->type);
	}

      printf("  %d\n", curr->nodeidx);
      prev = curr;
    }

  return 1;
}
    
  

/* Function: Trace2KHS()
 * 
 * Purpose:  From a traceback tree of seq, produce a
 *           secondary structure string. ">" and "<" are
 *           used for pairwise emissions; "." for single-stranded stuff.
 *           Note that structure is defined by pairwise emissions,
 *           not by Watson-Crick-isms and stacking rules.
 *           
 * Args:     tr          - traceback structure
 *           seq         - sequence, 0..rlen-1
 *           rlen        - length of seq and returned ss string
 *           watsoncrick - TRUE to annotate canonical pairs only
 *           ret_ss      - RETURN: alloc'ed secondary structure string
 *
 * Return:   ret_ss contains a string 0..rlen-1 containing the
 *           secondary structure. Must be free'd by caller.
 */
void
Trace2KHS(struct trace_s *tr, char *seq, int rlen, int watsoncrick, 
          char **ret_ss)  
{
  struct tracestack_s *dolist;
  struct trace_s      *curr;
  char                *ss;

  if ((ss = (char *) malloc (sizeof(char) * rlen+1)) == NULL)
    Die("malloc failed");
  memset(ss, '.', rlen);
  ss[rlen] = '\0';

  dolist = InitTracestack();
  PushTracestack(dolist, tr->nxtl);

  while ((curr = PopTracestack(dolist)) != NULL)
    {
      if ( curr->type      == uMATP_ST )
	{
	  if (! watsoncrick  ||
	      IsRNAComplement(seq[curr->emitl], seq[curr->emitr], TRUE))
	    {
	      ss[curr->emitl] = '>';
	      ss[curr->emitr] = '<';
	    }
	}

      if (curr->nxtr) PushTracestack(dolist, curr->nxtr);
      if (curr->nxtl) PushTracestack(dolist, curr->nxtl);
    }
  FreeTracestack(dolist);
  *ret_ss = ss;
}

/* Function: IsRNAComplement()
 * 
 * Purpose:  Returns TRUE if sym1, sym2 are Watson-Crick complementary.
 *           If allow_gu is TRUE, GU pairs also return TRUE.
 */
int
IsRNAComplement(char sym1, char sym2, int allow_gu)
{
  sym1 = toupper(sym1);
  sym2 = toupper(sym2);
  if (sym1 == 'T') sym1 = 'U';
  if (sym2 == 'T') sym2 = 'U';

  if ((sym1 == 'A' && sym2 == 'U') ||
      (sym1 == 'C' && sym2 == 'G') ||
      (sym1 == 'G' && sym2 == 'C') ||
      (sym1 == 'U' && sym2 == 'A') ||
      (allow_gu && sym1 == 'G' && sym2 == 'U') ||
      (allow_gu && sym1 == 'U' && sym2 == 'G'))
    return TRUE;
  else
    return FALSE;
}


/* Function: KHS2ct()
 * 
 * Purpose:  Convert a secondary structure string to an array of integers
 *           representing what position each position is base-paired 
 *           to (0..len-1), or -1 if none. This is off-by-one from a
 *           Zuker .ct file representation.
 *           
 *           The .ct representation can accomodate pseudoknots but the 
 *           secondary structure string cannot easily; the string contains
 *           "Aa", "Bb", etc. pairs as a limited representation of
 *           pseudoknots. The string contains "><" for base pairs.
 *           Other symbols are ignored. If allow_pseudoknots is FALSE,
 *           the pseudoknot symbols will be ignored and these positions
 *           will be treated as single stranded.
 *           
 * Return:   ret_ct is allocated here and must be free'd by caller.
 *           Returns 1 on success, 0 if ss is somehow inconsistent.
 */
int 
KHS2ct(char *ss, int len, int allow_pseudoknots, int **ret_ct)
{
  struct intstack_s *dolist[27];
  int *ct;
  int  i;
  int  pos, pair;
  int  status = 1;		/* success or failure return status */

  for (i = 0; i < 27; i++)
    dolist[i] = InitIntStack();

  if ((ct = (int *) malloc (len * sizeof(int))) == NULL)
    Die("malloc failed");
  for (pos = 0; pos < len; pos++)
    ct[pos] = -1;

  for (pos = 0; ss[pos] != '\0'; pos++)
    {
      if (ss[pos] > 127) status = 0; /* bulletproof against SGI buggy ctype.h */

      else if (ss[pos] == '>')	/* left side of a pair: push onto stack 0 */
	PushIntStack(dolist[0], pos);
      else if (ss[pos] == '<')	/* right side of a pair; resolve pair */
	{
	  if (! PopIntStack(dolist[0], &pair))
	    { status = 0; }
	  else
	    {
	      ct[pos]  = pair;
	      ct[pair] = pos;
	    }
	}
				/* same stuff for pseudoknots */
      else if (allow_pseudoknots && isupper((int) ss[pos]))
	PushIntStack(dolist[ss[pos] - 'A' + 1], pos);
      else if (allow_pseudoknots && islower((int) ss[pos]))
	{
	  if (! PopIntStack(dolist[ss[pos] - 'a' + 1], &pair))
	    { status = 0; }
	  else
	    {
	      ct[pos]  = pair;
	      ct[pair] = pos;
	    }
	}
      else if (allow_pseudoknots && !isgap(ss[pos])) status = 0; /* bad character */
    }

  for (i = 0; i < 27; i++)
    if ( FreeIntStack(dolist[i]) > 0)
      status = 0;

  *ret_ct = ct;
  return status;
}


