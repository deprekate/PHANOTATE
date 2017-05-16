/* emit.c
 * 1.0: Fri Jun 11 12:59:33 1993
 * 2.0: SRE, Thu Sep  9 13:44:18 1993
 * 
 * generate sequences randomly from a model.
 * 
 * The growing sequence is kept as a linked list (align_s).
 * The model tree is traversed by pushing nodes onto a stack
 * (m2ali_s). The information kept for each active node
 * is nodeidx, state type, and a pointer into the growing linked
 * list where the next emissions should go.
 * 
 * The recursion is to pop an active node off;
 * then, switch (statetype)
 *    MATP:       pick symbol pair. insert right symbol. insert
 *                left symbol. new insertion pointer on left symbol.
 *    INSL, MATL: pick symbol. insert symbol. new insertion pointer
 *                on new symbol.
 *    INSR, MATR: pick symbol. insert symbol. new insertion pointer
 *                stays where it was.
 *    BIFURC:     insert dummy symbol. one new insertion pointer
 *                stays where it was (BIFL), other points to the
 *                new dummy (BIFR).
 *    DELETE:     no symbol. new insertion pointer stays where it
 *                was.                                                   
 */

#include <stdio.h>
#include <stdlib.h>
#include "squid.h"

#include "structs.h"
#include "funcs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static void pick_transit(struct cm_s *cm, int oldidx, int oldtype,
			 int *ret_newidx, int *ret_newtype);
static void pick_double_emit(struct cm_s *cm, int nodeidx, char *ret_syml, char *ret_symr);
static void pick_single_emit(struct cm_s *cm, int nodeidx, int type, char *ret_sym);

static void pick_best_transit(struct cm_s *cm, int oldidx, int oldtype,
			      int *ret_newidx, int *ret_newtype);
static void pick_best_double(struct cm_s *cm, int nodeidx, char *ret_syml, char *ret_symr);
static void pick_best_single(struct cm_s *cm, int nodeidx, int type, char *ret_sym);


/* Function:  EmitSequence()
 * 
 * Purpose:   Generate a sequence probabilistically from a model.
 *            The returned sequence contains upper case letters
 *            for MATCH-generated positions, lower case letters
 *            for INSERT-generated positions.
 *            
 *            Caller is reponsible for free'ing memory allocated
 *            to the sequence.
 *            
 */
int
EmitSequence(struct cm_s     *cm,        /* model                                  */
	     int             watsoncrick,/* if TRUE, annotate only canonical pairs */
	     struct align_s **ret_ali,   /* RETURN: generated "alignment"          */
	     char           **ret_khseq, /* RETURN: generated structure            */
	     char           **ret_seq)   /* RETURN: generated sequence             */ 
{
  struct m2ali_s *emstack;
  struct align_s *emlist;
  int             oldidx;
  int             oldtype;
  struct align_s *oldafter;
  int             newidx;
  int             newtype;
  struct align_s *newafter;
  char            syml;
  char            symr;
  int             pos;
  char           *seq;
  char           *khseq;
  char            ssl, ssr;	/* secondary structure annotation */

  /* Initialize the linked list of emitted sequence, emlist.
   */
  emlist = Init_align();

  /* Initialize the pushdown stack for traversing the model,
   * emstack.
   */
  emstack = Init_m2ali();
  pick_transit(cm, 0, uBEGIN_ST, &newidx, &newtype);
  Push_m2ali(emstack, newidx, newtype, emlist); 


  /* While there's still active model nodes in the stack,
   * pop one off and deal with it.
   */
  while (Pop_m2ali(emstack, &oldidx, &oldtype, &oldafter))
    {
				/* look out for end */
      if (oldidx == -1)	continue;

      /* check for BIFURC, which makes automatic transits to BEGIN states
       * of next two segments. 
       */
      if (cm->nd[oldidx].type == BIFURC_NODE)
	{
				/* deal with right branch */
	  newafter = Insafter_align(0, '-', ' ', oldidx, uBIFURC_ST, oldafter); /* insert a dummy */
	  Push_m2ali(emstack, cm->nd[oldidx].nxt2, uBEGIN_ST, newafter);
				/* deal with left branch */
	  Push_m2ali(emstack, cm->nd[oldidx].nxt, uBEGIN_ST, oldafter);
	}

      else 
	{
	  switch (oldtype) {
	  case uDEL_ST:
	    pick_transit(cm, oldidx, oldtype, &newidx, &newtype);
	    Push_m2ali(emstack, newidx, newtype, oldafter);
	    break;

	  case uMATP_ST:
	    pick_double_emit(cm, oldidx, &syml, &symr);
	    if (! watsoncrick ||
		IsRNAComplement(syml, symr, TRUE))
	      { ssl = '>'; ssr = '<'; }
	    else
	      { ssl = ssr = '.'; }
	    (void) Insafter_align(0, symr, ssr, oldidx, oldtype, oldafter);
	    newafter = Insafter_align(0, syml, ssl, oldidx, oldtype, oldafter);

	    pick_transit(cm, oldidx, oldtype, &newidx, &newtype);
	    Push_m2ali(emstack, newidx, newtype, newafter);
	    break;
	
	  case uINSL_ST:
	  case uMATL_ST:
	    pick_single_emit(cm, oldidx, oldtype, &syml);
	    newafter = Insafter_align(0, syml, '.', oldidx, oldtype, oldafter);

	    pick_transit(cm, oldidx, oldtype, &newidx, &newtype);
	    Push_m2ali(emstack, newidx, newtype, newafter);
	    break;

	  case uINSR_ST:
	  case uMATR_ST:
	    pick_single_emit(cm, oldidx, oldtype, &symr);
	    (void) Insafter_align(0, symr, '.', oldidx, oldtype, oldafter);

	    pick_transit(cm, oldidx, oldtype, &newidx, &newtype);
	    Push_m2ali(emstack, newidx, newtype, oldafter);
	    break;
	
	  default:
	    Die("Unrecognized state type %d in model.", oldtype);
	  }
	}
    }
  Free_m2ali(emstack);
  
  /* Now go through and write the correct 'pos' fields in *emlist,
   * because the caller might expect them for some reason. 
   */
  pos = 0;
  for (newafter = emlist->nxt; newafter->nxt != NULL; newafter = newafter->nxt)
    {
      if (newafter->type == uDEL_ST)
	newafter->pos = -1;
      else
	{
	  newafter->pos = pos;
	  pos++;
	}
    }

  /* Now we extract the sequence from the linked list.
   * For now, we leave the dummy characters in.
   */
  if (! Align2kh(emlist, &seq, &khseq))
    Warn("Align2kh() failed");
  
  *ret_ali   = emlist;
  *ret_khseq = khseq;
  *ret_seq   = seq;
  return 1;
}



/* Function: EmitBestSequence()
 * 
 * Purpose:  Generate the most probable sequence from a model by picking 
 *           the most probable transitions and emissions.
 *           Very similar to EmitSequence(), above.
 *           
 */
int
EmitBestSequence(struct cm_s     *cm,         /* model                         */
		 int             watsoncrick, /* if TRUE, annotate only canonical pairs */
		 struct align_s **ret_ali,    /* RETURN: generated "alignment" */
		 char           **ret_khseq,  /* RETURN: generated structure   */
		 char           **ret_seq)    /* RETURN: generated sequence    */
{
  struct align_s *ali;          /* generated "alignment" linked list       */
  struct align_s *curr;         /* ptr to current insertion pt in ali      */
  struct align_s *new;          /* ptr to newly inserted pt in ali         */
  struct m2ali_s *stack;        /* pushdown stack for traversing model     */
  int             oldidx;	/* stateidx of current insertion pt in ali */
  int             oldtype;	/* subtype of current insertion pt in ali  */
  int             newidx;	/* stateidx of newly inserted pt in ali    */
  int             newtype;	/* subtype of newly inserted pt in ali     */
  char            syml;		/* emitted symbol to the left              */
  char            symr;		/* emitted symbol to the right             */
  int             pos;		/* position in seq                         */
  char           *seq;          /* RETURN: generated most probable sequence*/
  char           *khseq;        /* RETURN: structure rep. of seq           */
  char            ssl, ssr;	/* secondary structure annotation */

  /* Initialize the linked list of emitted sequence, ali
   */
  ali = Init_align();

  /* Initialize the pushdown stack for traversing the model
   */
  stack = Init_m2ali();
  pick_best_transit(cm, 0, uBEGIN_ST, &newidx, &newtype);
  Push_m2ali(stack, newidx, newtype, ali); 


  /* While there's still active model nodes in the stack,
   * pop one off and deal with it.
   */
  while (Pop_m2ali(stack, &oldidx, &oldtype, &curr))
    {
				/* look out for end */
      if (oldidx == -1)	continue;

      /* check for BIFURC, which makes automatic transits to BEGIN states
       * of next two segments. 
       */
      if (cm->nd[oldidx].type == BIFURC_NODE)
	{
				/* deal with right branch */
	  new = Insafter_align(0, '-', ' ', oldidx, uBIFURC_ST, curr); /* insert a dummy */
	  Push_m2ali(stack, cm->nd[oldidx].nxt2, uBEGIN_ST, new);
				/* deal with left branch */
	  Push_m2ali(stack, cm->nd[oldidx].nxt, uBEGIN_ST, curr);
	}

      else 
	{
	  switch (oldtype) {
	  case uDEL_ST:
	    pick_best_transit(cm, oldidx, oldtype, &newidx, &newtype);
	    Push_m2ali(stack, newidx, newtype, curr);
	    break;

	  case uMATP_ST:
	    pick_best_double(cm, oldidx, &syml, &symr);
	    if (! watsoncrick ||
		IsRNAComplement(syml, symr, TRUE))
	      { ssl = '>'; ssr = '<'; }
	    else
	      { ssl = '.'; ssr = '.'; }

	    (void) Insafter_align(0, symr, ssr, oldidx, oldtype, curr);
	    new = Insafter_align(0, syml, ssl, oldidx, oldtype, curr);

	    pick_best_transit(cm, oldidx, oldtype, &newidx, &newtype);
	    Push_m2ali(stack, newidx, newtype, new);
	    break;
	
	  case uINSL_ST:
	  case uMATL_ST:
	    pick_best_single(cm, oldidx, oldtype, &syml);
	    new = Insafter_align(0, syml, '.', oldidx, oldtype, curr);

	    pick_best_transit(cm, oldidx, oldtype, &newidx, &newtype);
	    Push_m2ali(stack, newidx, newtype, new);
	    break;

	  case uINSR_ST:
	  case uMATR_ST:
	    pick_best_single(cm, oldidx, oldtype, &symr);
	    (void) Insafter_align(0, symr, '.', oldidx, oldtype, curr);

	    pick_best_transit(cm, oldidx, oldtype, &newidx, &newtype);
	    Push_m2ali(stack, newidx, newtype, curr);
	    break;
	
	  default:
	    Die("Unrecognized state type %d in model.", oldtype);
	  }
	}
    }
  Free_m2ali(stack);
  
  /* Now go through and write the correct 'pos' fields in *ali,
   * because the caller might expect them for some reason. 
   */
  pos = 0;
  for (curr = ali->nxt; curr->nxt != NULL; curr = curr->nxt)
    {
      if (curr->type == uDEL_ST)
	curr->pos = -1;
      else
	{
	  curr->pos = pos;
	  pos++;
	}
    }

  /* Now we extract the sequence from the linked list.
   * For now, we leave the dummy characters in.
   */
  if (! Align2kh(ali, &seq, &khseq))
    Warn("Align2kh() failed");

  *ret_ali   = ali;
  *ret_khseq = khseq;
  *ret_seq   = seq;
  return 1;
}
		 



/* Function: pick_transit()
 * 
 * Purpose:  Pick a random state transition, given a current state
 *           (specified by a stateidx and a subtype). Pass back
 *           the new state (newidx, newtype).
 */
static void
pick_transit(struct cm_s *cm,
	     int          oldidx,
	     int          oldtype,
	     int         *ret_newidx,
	     int         *ret_newtype)
{
  int    newidx;
  int    newtype;
  double sum;
  double roll;

  /* Picking a new subtype involves rolling a random
   * fraction and examining the appropriate row of the
   * 7x7 state transition matrix.
   */
  sum  = 0.0;
  roll = sre_random();
  for (newtype = 0; newtype < STATETYPES; newtype++)
    {
      sum += cm->nd[oldidx].tmx[oldtype][newtype];
      if (roll <= sum) break;
    }
  if (newtype == STATETYPES)
    Die("Failed to transit from stateidx %d subtype %d, roll %.2f",
	oldidx, oldtype, roll);
  
  /* Picking a new nodeidx is a function of the current
   * state type. This function should never be called for BIFURCs.
   */
  if (newtype == INSL_ST || newtype == INSR_ST)
    newidx = oldidx;
  else
    newidx = cm->nd[oldidx].nxt;
  
  *ret_newidx  = newidx;
  *ret_newtype = UniqueStatetype(cm->nd[newidx].type, newtype);
}



/* Function: pick_double_emit()
 * 
 * Purpose:  Given a model and a current state (stateidx, subtype),
 *           which must be an INSC or MATC, pick a pairwise emission
 *           (syml, symr) according to the probabilities in the 
 *           appropriate emission matrix.
 */
static void
pick_double_emit(struct cm_s *cm,
		 int          nodeidx,
		 char        *ret_syml,
		 char        *ret_symr)
{
  double   sum;
  double   roll;
  int      i, j;

  sum  = 0.0;
  roll = sre_random();
  for (i = 0; i < ALPHASIZE; i++)
    for (j = 0; j < ALPHASIZE; j++)
      {
	sum += cm->nd[nodeidx].mp_emit[i][j];
	if (roll <= sum) goto breakout;
      }
 breakout:
  *ret_syml = ALPHABET[i];
  *ret_symr = ALPHABET[j];
}




/* Function: pick_single_emit()
 * 
 * Purpose:  Given a model and a current state (nodeidx, type),
 *           which must be an INS(L/R) or MAT(L/R), pick an emission
 *           (sym) according to the probabilities in the 
 *           appropriate emission vector.
 */
static void
pick_single_emit(struct cm_s *cm,
		 int          nodeidx,
		 int          type,
		 char        *ret_sym)
{
  double   sum;
  double   roll;
  int      i;
  double  *emit;

				/* find correct emission vector */
  switch (type) {
  case uINSL_ST: emit = cm->nd[nodeidx].il_emit; break;
  case uINSR_ST: emit = cm->nd[nodeidx].ir_emit; break;
  case uMATL_ST: emit = cm->nd[nodeidx].ml_emit; break;
  case uMATR_ST: emit = cm->nd[nodeidx].mr_emit; break;
  default:
    Die("can't single emit from state type %d", type);
  }

  sum  = 0.0;
  roll = sre_random();
  for (i = 0; i < ALPHASIZE; i++)
    {
      sum += emit[i];
      if (roll <= sum) break;
    }
  *ret_sym = ALPHABET[i];
}




/* Function: pick_best_transit()
 * 
 * Purpose:  Pick most probable state transition, given a current state
 *           (specified by a nodeidx and a type). Pass back
 *           the new state (newidx, newtype).
 */
static void
pick_best_transit(struct cm_s *cm,
		  int          oldidx,
		  int          oldtype,
		  int         *ret_newidx,
		  int         *ret_newtype)
{
  int    y;
  int    newidx;
  int    newtype;
  double best;
  
				/* find maximum probability */
  best = 0.0;
  for (y = 0; y < STATETYPES; y++)
    if (cm->nd[oldidx].tmx[oldtype][y] > best)
      {
	best    = cm->nd[oldidx].tmx[oldtype][y];
	newtype = y;
      }
  
  /* Picking a new nodeidx is a function of the current
   * type. This function should never be called for BIFURCs.
   */
  if (newtype == INSL_ST || newtype == INSR_ST)
    newidx = oldidx;
  else
    newidx = cm->nd[oldidx].nxt;
  
  *ret_newidx  = newidx;
  *ret_newtype = UniqueStatetype(cm->nd[newidx].type, newtype);
}



/* Function: pick_best_double()
 * 
 * Purpose:  Given a model and a current state (nodeidx, type),
 *           (which must a MATP), pick most probable pairwise emission
 *           (syml, symr) according to the probabilities in the 
 *           appropriate emission matrix.
 *           
 */
static void
pick_best_double(struct cm_s *cm,
		 int          nodeidx,
		 char        *ret_syml,
		 char        *ret_symr)
{
  double   best;
  int      i, j;
  int      besti, bestj;

  best = 0.0;
  for (i = 0; i < ALPHASIZE; i++)
    for (j = 0; j < ALPHASIZE; j++)
      if (cm->nd[nodeidx].mp_emit[i][j] > best)
	{
	  best  = cm->nd[nodeidx].mp_emit[i][j];
	  besti = i;
	  bestj = j;
	}
  *ret_syml = ALPHABET[besti];
  *ret_symr = ALPHABET[bestj];
}




/* Function: pick_best_single()
 * 
 * Purpose:  Given a model and a current state (nodeidx, type),
 *           which must be an INS(L/R) or MAT(L/R), pick most probable emission
 *           (sym) according to the probabilities in the 
 *           appropriate emission vector.
 */
static void
pick_best_single(struct cm_s *cm,
		 int          nodeidx,
		 int          type,
		 char        *ret_sym)
{
  double   best;
  int      besti;
  int      i;
  double  *emit;

				/* find correct emission vector */
  switch (type) {
  case uINSL_ST: emit = cm->nd[nodeidx].il_emit; break;
  case uINSR_ST: emit = cm->nd[nodeidx].ir_emit; break;
  case uMATL_ST: emit = cm->nd[nodeidx].ml_emit; break;
  case uMATR_ST: emit = cm->nd[nodeidx].mr_emit; break;
  default:   Die("can't single emit from type %d", type);
  }

  best = 0.0;
  for (i = 0; i < ALPHASIZE; i++)
    if (emit[i] > best)
      {
	best  = emit[i];
	besti = i;
      }
  *ret_sym = ALPHABET[besti];
}








