/* align.c
 * SRE, Tue Jun 29 14:05:48 1993
 * 2.0: Thu Sep 30 14:43:05 1993
 * 
 * Code for producing a multiple sequence alignment from tracebacks. 
 * 
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "squid.h"
#include "structs.h"
#include "funcs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif



/* Function: create_master()
 * 
 * Purpose:  Produce a "master copy" of the linear order in
 *           which a model's states produce a sequence. Used
 *           for reference during alignment construction.
 *           
 *           This linked list is created with Init_align()
 *           and must be free'd by Free_align(). Each structure
 *           in the list represents a possible emission from
 *           a match state, and contains the index of the responsible
 *           state. The sym field is used as an internal flag
 *           (temporary dummy structures are used for bifurc's, and
 *           flagged in the sym field for later removal). The
 *           pos and substate fields are not used or meaningful. ret_len
 *           is the maximum number of symbols that can be
 *           emitted by non-insert states.
 */
static int
create_master(struct cm_s     *cm, 
	      struct align_s **ret_master,
	      int             *ret_len)
{
  struct m2ali_s *stack;	/* stack used to traverse the model cv */
  struct align_s *master;       /* RETURN: maximal match emission alignment */
  struct align_s *curr;
  int             oldidx;
  int             oldtype;
  struct align_s *oldafter;
  int             newidx;
  int             newtype;
  struct align_s *newafter;
  int             len;

  /* Initialize the linked list of state order, master
   */
  master = Init_align();

  /* Initialize pushdown stack for traversing the model.
   */
  stack   = Init_m2ali();
  newidx  = (cm->nodes > 1) ? 1 : -1;
  newtype = uMATP_ST;		/* assume 'worst' case, all nodes generate pairs */
  Push_m2ali(stack, newidx, newtype, master);

  /* While there are still active nodes on the stack, pop one off
   * and deal with it.
   */
  len = 0;
  while (Pop_m2ali(stack, &oldidx, &oldtype, &oldafter))
    {
      if (oldidx == -1) continue; /* END */
      
      if (cm->nd[oldidx].type == BIFURC_NODE)
	{
				/* deal with right branch.
				 * Gotta insert a dummy temporarily */
	  newafter = Insafter_align(-1, '-', ' ', oldidx, uBIFURC_ST, oldafter);
	  Push_m2ali(stack, cm->nd[oldidx].nxt2, uBEGIN_ST, newafter);
	  len++;
				/* deal with left branch */
	  Push_m2ali(stack, cm->nd[oldidx].nxt,  uBEGIN_ST, oldafter);
	}

      else if (cm->nd[oldidx].type == BEGINL_NODE ||
	       cm->nd[oldidx].type == BEGINL_NODE ||
	       cm->nd[oldidx].type == ROOT_NODE)
	{			/* BEGIN's aren't reponsible for any MAT states */
	  Push_m2ali(stack, cm->nd[oldidx].nxt,  uBEGIN_ST, oldafter);
	}

      else
	{
	  (void) Insafter_align(-1, '.', ' ', oldidx, uMATP_ST, oldafter);
	  newafter = Insafter_align(-1, '.', ' ', oldidx, uMATP_ST, oldafter);
	  len += 2;
	  Push_m2ali(stack, cm->nd[oldidx].nxt, uBEGIN_ST, newafter);
	}
    }
	  
  /* Remove the temporary dummies used to BIFURC
   */
  for (curr = master->nxt; curr->nxt != NULL; curr = curr->nxt)
    if (curr->nxt->sym == '-')
      {
	Delafter_align(curr);
	len--;
      }

  Free_m2ali(stack);
  *ret_len    = len;
  *ret_master = master;
#ifdef DEBUG
/*  print_align(master); */
#endif 

  return 1;
}


/* Function: Trace2ali()
 * 
 * Purpose:  Given a traceback (tree-structure alignment of a model
 *           to a sequence), construct a linear linked list representation
 *           (align_s) of the sequence alignment to the model.
 *           
 * Args:     seq         - 0..len-1 sequence to align
 *           tr          - traceback containing tree-wise alignment
 *           watsoncrick - if TRUE, only canonical pairs get structure annotation
 *           ret_ali     - RETURN: linear list alignment               
 */
int
Trace2ali(char *seq, struct trace_s *tr, int watsoncrick, struct align_s **ret_ali)   
{
  struct align_s      *ali;     /* RETURN: linear list of alignment        */
  struct t2ali_s      *stack;	/* stack used to traverse the traceback tr */
  struct trace_s      *currtr;
  struct align_s      *newafter;
  struct align_s      *oldafter;
  struct align_s      *curr;
  char                 ssl, ssr;/* symbols <.> for secondary structure rep. */

  /* Initialize the linked list for the alignment of sequence to model
   */
  ali = Init_align();

  /* Initialize the pushdown stack for traversal of the traceback
   */
  stack = Init_t2ali();
  Push_t2ali(stack, tr, ali);

  while (Pop_t2ali(stack, &currtr, &oldafter))
    {
      switch (currtr->type) {
      case END_ST: 
	break;	/* ignore END states */

      case uBIFURC_ST:
				/* deal with right branch; insert a dummy */
	newafter = Insafter_align(-1, '*', ' ', currtr->nodeidx, currtr->type, oldafter);
	Push_t2ali(stack, currtr->nxtr, newafter);
				/* deal with left branch */
	Push_t2ali(stack, currtr->nxtl, oldafter);
	break;

      case uBEGIN_ST:
	Push_t2ali(stack, currtr->nxtl, oldafter);
	break;

      case uDEL_ST:		
	Insafter_align(-1, '-', '.', currtr->nodeidx, currtr->type, oldafter);
	Push_t2ali(stack, currtr->nxtl, oldafter);
	break;

      case uMATP_ST:
	if (! watsoncrick ||
	    IsRNAComplement(seq[currtr->emitr], seq[currtr->emitl], TRUE))
	  { ssr = '<'; ssl = '>'; }
	else
	  { ssr = '.'; ssl = '.'; }
	
	(void) Insafter_align(currtr->emitr, seq[currtr->emitr], ssr,
			      currtr->nodeidx, currtr->type, oldafter);
	newafter = Insafter_align(currtr->emitl, seq[currtr->emitl], ssl,
				  currtr->nodeidx, currtr->type, oldafter);
	Push_t2ali(stack, currtr->nxtl, newafter);
	break;

      case uINSL_ST:
      case uMATL_ST:
	newafter = Insafter_align(currtr->emitl, seq[currtr->emitl], '.',
				      currtr->nodeidx, currtr->type, oldafter);
	Push_t2ali(stack, currtr->nxtl, newafter);
	break;

      case uINSR_ST:
      case uMATR_ST:
	(void) Insafter_align(currtr->emitr, seq[currtr->emitr], '.',
			      currtr->nodeidx, currtr->type, oldafter);
	Push_t2ali(stack, currtr->nxtl, oldafter);
	break;
      }
    }
	 
  /* Remove the temporary dummies used to BIFURC
   */
  for (curr = ali->nxt; curr->nxt != NULL; curr = curr->nxt)
    if (curr->nxt->sym == '*')
      Delafter_align(curr);

  Free_t2ali(stack);
  *ret_ali = ali;

#ifdef DEBUG
/*   print_align(ali); */
#endif
  return 1;
}


/* Function: Traces2Alignment()
 * 
 * Purpose:  Given a set of tracebacks for alignments of multiple sequences to 
 *           a model, construct a multiple sequence alignment.
 *           
 *           The tricky bit, which involves some precalculations, is allowing 
 *           the proper amount of space for insertions.
 */
int
Traces2Alignment(char           **rseqs,
		 SQINFO          *sqinfo,
		 struct trace_s **tr,
		 int              nseq,
		 struct cm_s     *cm,
		 int              watsoncrick, /* TRUE to annotate only canonical pairs */
		 char          ***ret_aseqs,
		 AINFO           *ainfo)
{
  struct align_s  *master;      /* representation of how MAT columns map onto the model  */
  struct align_s **ali;         /* array of align_s alignments of the sequences to model */
  int             *matuse;      /* 0 if MAT column is never used, 1 otherwise, 1..len */
  int             *insuse;      /* per sequence, count use of inserts between MAT columns; 0..len */
  int             *max_insuse;  /* overall maxima, keep track of ins use between columns */
  int             *matpos;      /* array of MAT column positions in the alignment, 1..len */
  int              len;		/* length of master - maximum number of MAT columns */
  int              idx;		/* counter for sequences */
  struct align_s  *currmaster;
  struct align_s  *currali;
  int              aseqlen;	/* length of multiple sequence alignment */
  char           **aseqs;       /* RETURN: multiple sequence alignment   */
  char           **ss;          /* secondary structures                  */
  int              apos;	/* position in absolute alignment columns */
  int              matcol;	/* position in MAT column coord arrays (matuse, matpos) */

  /* First we use the model to calculate "master", which will define
   * the columns of the multiple sequence alignment (MAT-produced) and
   * represents how they map onto the model states. Only the stateidx field of 
   * the align_s structures is meaningful. 
   */
  if (! create_master(cm, &master, &len)) return 0;

  /* Next we "invert" each traceback (from seq->model tree alignments 
   * to model->seq linear alignments) and create an array "ali" of individual
   * alignments of the model to the sequences.
   */
  if ((ali = (struct align_s **) malloc (nseq * sizeof(struct align_s *))) == NULL)
    { Die("Memory failure, line %d of %s", __LINE__, __FILE__); return 0; }
  for (idx = 0; idx < nseq; idx++)
    if (! Trace2ali(rseqs[idx], tr[idx], watsoncrick, &ali[idx])) return 0;
	
  /* Now we're ready to start counting MAT and INS use.
   * 
   * For MAT use, all we're doing is determining whether a given 
   * column in the master alignment is used or not (matuse[1..len] is 
   * 1 if yes, 0 if no)
   * 
   * For INS use, we are counting the maximum number of occurrences
   * of insert emissions between each MAT column of the master. max_insuse[0..len]
   * keeps these numbers. max_insuse[5] is the maximum number of inserted
   * symbols between columns 5 and 6, for example. Because there is 
   * some redundancy in the model -- different INS states may emit in 
   * the same place -- we have to first increment a counter array, insuse,
   * for each individual sequence.
   */
  if (((matuse     = (int *) calloc (len+1, sizeof(int))) == NULL) ||
      ((insuse     = (int *) calloc (len+1, sizeof(int))) == NULL) ||
      ((max_insuse = (int *) calloc (len+1, sizeof(int))) == NULL) )
    { Die("Memory failure, line %d of %s", __LINE__, __FILE__); return 0; }

  for (idx = 0; idx < nseq; idx++)
    {
      for (matcol = 0; matcol <= len; matcol++)
	insuse[matcol] = 0;

      matcol = 0;
      currmaster = master->nxt;
      for (currali = ali[idx]->nxt; currali != NULL; currali = currali->nxt)
	{
	  switch (currali->type) {
	  case uMATP_ST:
	  case uMATR_ST:
	  case uMATL_ST:
	  case uDEL_ST:
	    matcol++; 
	    while (currmaster->nodeidx != currali->nodeidx)
	      {	currmaster = currmaster->nxt; matcol++; }
	    matuse[matcol] = 1;
	    currmaster = currmaster->nxt; 
	    break;

	  case uINSR_ST:
	  case uINSL_ST:
	    insuse[matcol]++;
	    break;
	  }
	}
				/* update max_insuse with new maxima, if any*/
      for (matcol = 0; matcol <= len; matcol++)
	if (insuse[matcol] > max_insuse[matcol])
	  max_insuse[matcol] = insuse[matcol];
    }

				/* calculate length of mult seq alignment, and alloc */
  aseqlen = 0;
  for (matcol = 0; matcol <= len; matcol++)
    {
      if (matuse[matcol] == 1) aseqlen++;
      aseqlen += max_insuse[matcol];
    }
  if ((aseqs = (char **) malloc (nseq * sizeof(char *))) == NULL ||
      (ss    = (char **) malloc (nseq * sizeof(char *))) == NULL)
    { Die("Memory failure, line %d of %s", __LINE__, __FILE__); return 0; }
  for (idx = 0; idx < nseq; idx++)
    if ((aseqs[idx] = (char *) malloc ((aseqlen+1) * sizeof(char))) == NULL ||
	(ss[idx]    = (char *) malloc ((aseqlen+1) * sizeof(char))) == NULL)
      { Die("Memory failure, line %d of %s", __LINE__, __FILE__); return 0; }


  /* Now we use matuse and max_insuse to calculate an array for the 
   * coordinates of the MAT columns in the multiple alignment.
   */
  if ((matpos = (int *) calloc (len+1 , sizeof(int))) == NULL)
    { Die("Memory failure, line %d of %s", __LINE__, __FILE__); return 0; }
  for (matcol = 1; matcol <= len; matcol++)
    matpos[matcol] = matpos[matcol-1] + max_insuse[matcol-1] + matuse[matcol-1];

  /* And finally, we're ready to actually construct the multiple sequence
   * alignment. The resulting alignment is flushed-right with gaps.
   */
  for (idx = 0; idx < nseq; idx++)
    {
      matcol = 0;
      apos = 0;
      currmaster = master;
      for (currali = ali[idx]->nxt; currali != NULL; currali = currali->nxt)
	{
	  switch (currali->type) {
	  case uMATP_ST:
	  case uMATR_ST:
	  case uMATL_ST:
				/* goes in a MAT column */
	    while (currmaster->nodeidx != currali->nodeidx)
	      {	currmaster = currmaster->nxt; matcol++; }
	    for (; apos < matpos[matcol]; apos++)
	      {
		aseqs[idx][apos] = '.';
		ss[idx][apos]    = ' ';
	      }

	    aseqs[idx][apos] = toupper((int)currali->sym);
	    ss[idx][apos]    = currali->ss;

	    apos++;
	    currmaster = currmaster->nxt; 
	    matcol++;
	    break;

	  case uINSR_ST:
	  case uINSL_ST:
	    aseqs[idx][apos] = tolower((int)currali->sym);
	    ss[idx][apos]    = currali->ss;
	    apos++;
	    break;

	  case uDEL_ST:
	    aseqs[idx][apos] = '.';
	    ss[idx][apos]    = ' ';
	    apos++;
	    currmaster = currmaster->nxt; 
	    matcol++;
	    break;
	  }
	}
				/* flush right */
      for (; apos < aseqlen; apos++)
	{
	  aseqs[idx][apos] = '.';
	  ss[idx][apos]    = ' ';
	}
      aseqs[idx][apos] = '\0';
      ss[idx][apos]    = '\0';
    }

  for (idx = 0; idx < nseq; idx++)
    Free_align(ali[idx]);
  Free_align(master);
  free(matpos);
  free(max_insuse);
  free(insuse);
  free(matuse);
  free(ali);

  if (ainfo != NULL)
    {
      int   leftcount;
      int   rightcount;

      ainfo->flags = 0;

      strcpy(ainfo->au, "CM RNA automatic alignment");
      ainfo->flags |= AINFO_AUTH;

      ainfo->alen = aseqlen;
      ainfo->flags |= AINFO_ALEN;

				/* copy sqinfo structure array */
      if ((ainfo->sqinfo = (SQINFO *) malloc (sizeof(SQINFO) * nseq)) == NULL)
	Die("malloc failed");
      for (idx = 0; idx < nseq; idx++)
	SeqinfoCopy(&(ainfo->sqinfo[idx]), &(sqinfo[idx]));
      
      /* Construct a consensus structure string.
       * Secondary structure strings, ss, are currently aligned to
       * the aseqs. Calculate an aligned consensus structure from
       * them. 
       */
      if ((ainfo->cs = (char *) malloc (sizeof(char) * (aseqlen+1))) == NULL)
	Die("malloc failed");
      for (apos = 0; apos < aseqlen; apos ++)
	{
	  leftcount = rightcount = 0;
	  for (idx = 0; idx < nseq; idx++)
	    if      (ss[idx][apos] == '<') rightcount++;
	    else if (ss[idx][apos] == '>') leftcount++;
	  
	  if      (rightcount > nseq / 2) ainfo->cs[apos] = '<';
	  else if (leftcount  > nseq / 2) ainfo->cs[apos] = '>';
	  else    ainfo->cs[apos] = '.';
	}
      ainfo->cs[aseqlen] = '\0';
      ainfo->flags |= AINFO_CS;

      /* Construct individual secondary structure strings by de-aligning
       * the individuals.
       */
      for (idx = 0; idx < nseq; idx++)
	{
	  MakeDealignedString(aseqs[idx], aseqlen, ss[idx], &(ainfo->sqinfo[idx].ss));
	  ainfo->sqinfo[idx].flags |= SQINFO_SS;
	}
    }

  Free2DArray(ss, nseq);
  *ret_aseqs = aseqs;
  return 1;
}


