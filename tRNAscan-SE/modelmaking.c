/* modelmaking.c
 * Tue Oct  4 15:33:21 1994
 * 
 * Bring together common elements of the model construction process.
 * Also, provides EasyModelmaker() for making a model given a structure.
 * 
 * All model makers have in common that they construct a "master" traceback
 * for the alignment, specifying which columns are match vs. insert and
 * how the model tree branches. This traceback is assigned a numbering
 * system by NumberMasterTrace(), which returns the number of nodes;
 * the caller then allocates a new CM. This new model is numbered (assigned
 * a branching structure) by TopofyNewCM(). Then individual tracebacks 
 * are constructed from individual aligned sequences by Transmogrify(). 
 * The individual tracebacks are counted into a new model with TraceCount()
 * and the counts converted to probabilities with ProbifyCM().
 * 
 * The master tree is a (slightly misused) trace_s structure with the following 
 * properties:
 *     insert columns are not represented at all. Transmogrify() must deal.
 *         
 *     emitl, emitr == 0..alen-1 coords of assigned columns. Set and valid for all
 *                     nodes, even non-emitters. END values will be
 *                     on the off-diagonal, emitl = emitr+1. (If this is
 *                     not true, Transmogrify() breaks.) The trace
 *                     construction function is responsible for this.
 *                     
 *     nodeidx == this is numbered by preorder traversal by NumberMasterTrace(),
 *                which also numbers a model. ENDs are not explicitly 
 *                represented in a CM, so they get nodeidx = -1.
 *                
 *     type == a number 0..6 for *node* type. (Usually this is a unique
 *             state type identifier.) ENDs are -1. The trace construction
 *             function is reponsible for this.
 */

#include <stdio.h>
#include <stdlib.h>

#include "version.h"
#include "squid.h"
#include "structs.h"
#include "funcs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Function: NumberMasterTrace()
 * 
 * Purpose:  Given a master trace for an alignment, number the trace tree
 *           (in the nodeidx field) in preorder traversal. END nodes
 *           will not be represented explicitly in the final CM. They
 *           get numbered -1.
 */
void
NumberMasterTrace(struct trace_s *mtr, int *ret_nodes)
{
  struct trace_s      *curr;
  struct tracestack_s *dolist;
  int nodes = 0;  

  dolist = InitTracestack();
  PushTracestack(dolist, mtr->nxtl); /* push root onto stack */

  while ((curr = PopTracestack(dolist)) != NULL)
    {
      if (curr->nxtl == NULL)   /* END node */
	curr->nodeidx = -1; 
      else		            /* other nodes */
	curr->nodeidx = nodes++;

      if (curr->nxtr != NULL)  PushTracestack(dolist, curr->nxtr);
      if (curr->nxtl != NULL)  PushTracestack(dolist, curr->nxtl);
    }      

  FreeTracestack(dolist);
  *ret_nodes = nodes;
}


/* Function: TopofyNewCM()
 * 
 * Purpose:  Given the mtr master traceback tree, which defines the
 *           topology of the model, write the nxt and nxt2 connections
 *           into the model. For the most part, these are already
 *           contained in mtr thanks to NumberMasterTrace(); the
 *           only tricky bit is converting END states from multiple 
 *           real states (in mtr) to -1 nxt flags in the cm.
 *           
 * Return:   1 on success, 0 on failure.         
 */
void
TopofyNewCM(struct cm_s *cm, struct trace_s *mtr)
{
  struct tracestack_s *dolist;
  struct trace_s      *curr;

  dolist = InitTracestack();
  PushTracestack(dolist, mtr->nxtl); /* push ROOT onto stack */

  while ((curr = PopTracestack(dolist)) != NULL)
    {
      if (curr->nxtl == NULL) continue; /* ignore ENDs */
      if (curr->nxtr != NULL)		/* deal with BIFURC states */
	{
	  cm->nd[curr->nodeidx].nxt2 = curr->nxtr->nodeidx;
	  PushTracestack(dolist, curr->nxtr);
	}
      else
	cm->nd[curr->nodeidx].nxt2 = -1;

				/* watch out for curr pointing to END states. */
      cm->nd[curr->nodeidx].type = curr->type;
      cm->nd[curr->nodeidx].nxt = (curr->nxtl->nxtl == NULL) ? -1 : curr->nxtl->nodeidx;
      PushTracestack(dolist, curr->nxtl);
    }      
  FreeTracestack(dolist);
}


/* Function: Transmogrify()
 * 
 * Purpose:  Given a master consensus traceback, create an individual
 *           "fake" traceback. The fake traceback contains inserts
 *           and converts the type field of mtr (which contains NODE
 *           type indices) into _ST type indices, including proper
 *           classification of nodes into DEL_ST or the various match
 *           states depending on what aseq[idx] looks like.
 *           
 * Args:     mtr      - master consensus traceback tree
 *           aseq     - 0..alen-1 aligned sequence.
 *           ret_tr   - RETURN: individual traceback
 *           ret_pool - RETURN: memory pool for traceback
 *           
 * Return:   (void). *ret_tr must be free'd by the caller
 */
void
Transmogrify(struct trace_s *mtr, char *aseq, struct trace_s **ret_tr, struct trmem_s **ret_pool)
{
  struct trace_s *tr;
  struct trmem_s *pool;
  struct tracestack_s *mtr_stack;
  struct tracestack_s *tr_stack;
  struct trace_s      *curr_mtr;
  struct trace_s      *curr_tr;
  int                 i2,j2;

  mtr_stack = InitTracestack();
  tr_stack  = InitTracestack();
  InitTrace(&tr, &pool);

  /* Push ROOT onto both stacks
   */
  PushTracestack(mtr_stack, mtr->nxtl);
  PushTracestack(tr_stack, AttachTrace(tr, pool, -1, -1, 0, uBEGIN_ST));

  while ((curr_mtr = PopTracestack(mtr_stack)) != NULL)
    {
      curr_tr = PopTracestack(tr_stack);

      switch (curr_mtr->type) {
      case uEND_ST: 
	DeleteTracenode(curr_tr, pool);
	break;

      case MATP_NODE:
	if (isgap(aseq[curr_mtr->emitl]))
	  {
	    if (isgap(aseq[curr_mtr->emitr])) curr_tr->type = uDEL_ST;
	    else                              curr_tr->type = uMATR_ST;
	  }
	else
	  {
	    if (isgap(aseq[curr_mtr->emitr])) curr_tr->type = uMATL_ST;
	    else                              curr_tr->type = uMATP_ST;
	  }
		
	/* May have to deal with INSL and INSR; INSL precedes INSR
	 */
	for (i2 = curr_mtr->emitl+1; i2 < curr_mtr->nxtl->emitl; i2++)
	  if (!isgap(aseq[i2]))
	    curr_tr = AttachTrace(curr_tr, pool, i2, curr_mtr->emitr, curr_mtr->nodeidx, uINSL_ST);

	/* May have to deal with INSR
	 */
	for (j2 = curr_mtr->emitr-1; j2 > curr_mtr->nxtl->emitr; j2--)
	  if (! isgap(aseq[j2]))
	    curr_tr = AttachTrace(curr_tr, pool, curr_mtr->nxtl->emitl, j2, curr_mtr->nodeidx, uINSR_ST);
	break;
	
      case MATL_NODE:
	if (isgap(aseq[curr_mtr->emitl])) curr_tr->type = uDEL_ST;
	else                              curr_tr->type = uMATL_ST;
	    
	/* May have to deal with INSL
	 */
	for (i2 = curr_mtr->emitl+1; i2 < curr_mtr->nxtl->emitl; i2++)
	  if (!isgap(aseq[i2]))
	    curr_tr = AttachTrace(curr_tr, pool, i2, curr_mtr->emitr, curr_mtr->nodeidx, uINSL_ST);
	break;

      case MATR_NODE:
	if (isgap(aseq[curr_mtr->emitr])) curr_tr->type = uDEL_ST;
	else                              curr_tr->type = uMATR_ST;

	/* May have to deal with INSR */
	for (j2 = curr_mtr->emitr-1; j2 > curr_mtr->nxtl->emitr; j2--)
	  if (! isgap(aseq[j2]))
	    curr_tr = AttachTrace(curr_tr, pool, curr_mtr->nxtl->emitl, j2, curr_mtr->nodeidx, uINSR_ST);
	break;

      case BIFURC_NODE: curr_tr->type  = uBIFURC_ST;  break;
      case BEGINL_NODE: curr_tr->type  = uBEGIN_ST;   break;

      case BEGINR_NODE:
	curr_tr->type  = uBEGIN_ST; 

	/* May have to deal with INSL.
	 * Inserts from BEGINR are *inclusive* of i
	 */
	for (i2 = curr_mtr->emitl; i2 < curr_mtr->nxtl->emitl; i2++)
	  if (!isgap(aseq[i2]))
	    curr_tr = AttachTrace(curr_tr, pool, i2, curr_mtr->emitr, curr_mtr->nodeidx, uINSL_ST);

	break;
	    
      case ROOT_NODE:
	curr_tr->type  = uBEGIN_ST; 

	/* May have to deal with INSL and INSR; note INSL precedes INSR
	 * inserts from root are inclusive of i  
	 */
	for (i2 = curr_mtr->emitl; i2 < curr_mtr->nxtl->emitl; i2++)
	  if (!isgap(aseq[i2]))
	    curr_tr = AttachTrace(curr_tr, pool, i2, curr_mtr->emitr, curr_mtr->nodeidx, uINSL_ST);

	/* May have to deal with INSR
	 */
	for (j2 = curr_mtr->emitr; j2 > curr_mtr->nxtl->emitr; j2--)
	  if (! isgap(aseq[j2]))
	    curr_tr = AttachTrace(curr_tr, pool, curr_mtr->nxtl->emitl, j2, curr_mtr->nodeidx, uINSR_ST);
	break;

      default: Die("Invalid node type %d", curr_mtr->type);
      }

      /* Push the children onto stacks, if they're not END nodes
       */
      if (curr_mtr->nxtr != NULL) 
	{
	  PushTracestack(mtr_stack, curr_mtr->nxtr);
	  PushTracestack(tr_stack, AttachTrace(curr_tr, pool, curr_mtr->nxtr->emitl, curr_mtr->nxtr->emitr,
					       curr_mtr->nxtr->nodeidx, curr_mtr->nxtr->type));
	}
      if (curr_mtr->nxtl != NULL) 
	{
	  PushTracestack(mtr_stack, curr_mtr->nxtl);
	  PushTracestack(tr_stack, AttachTrace(curr_tr, pool, curr_mtr->nxtl->emitl, curr_mtr->nxtl->emitr,
					       curr_mtr->nxtl->nodeidx, curr_mtr->nxtl->type));
	}
    }
  FreeTracestack(mtr_stack);
  FreeTracestack(tr_stack);
  *ret_pool = pool;
  *ret_tr   = tr;
}



/* Function: EasyModelmaker()
 * 
 * Purpose:  The customer always knows best.
 * 
 *           Construct a model given a stated structure. The structure
 *           is provided via a "cs" (consensus sequence) line, as would
 *           occur in an annotated SELEX file. Only > and < characters
 *           in this line are interpreted (as base pairs).
 *           
 *           Match vs. insert can be determined one of two ways. By default,
 *           the assignment is made by "gapthresh"; for columns with
 *           fractional occurence of gaps greater than this, the column
 *           is assigned to insert. If "use_rf" is TRUE, the rf (reference)
 *           line is interpreted as the assignment -- columns with non-space
 *           characters in the rf line are assigned to MATCH.
 *           
 *           Both rf and cs are provided in the ainfo structure.
 *           
 * Args:     aseq    - aligned sequences. [0..nseq-1] by [0..alen-1]
 *           ainfo   - info about the alignment, including alen, cs, 
 *                     and rf
 *           nseq    - number of seqs in aseq
 *           prior   - prior distributions for CM construction
 *           gapthresh - over this fraction of gaps, assign column as INS
 *           use_rf  - if TRUE, use rf field of ainfo for MAT/INS assignment
 *           ret_cm  - RETURN: new model                      (maybe NULL)
 *           ret_mtr - RETURN: master traceback for alignment (maybe NULL)
 *           
 * Return:   void
 *           cm is allocated here. FreeCM(*ret_cm).
 *           tr is allocated here. FreeTrace() on each one, then free(*ret_tr).
 */
void
EasyModelmaker(char **aseq, AINFO *ainfo, int nseq, struct prior_s *prior, 
	       double gapthresh, int use_rf, struct cm_s **ret_cm, struct trace_s **ret_mtr)
{
  struct cm_s         *cm;	/* new covariance model                */
  struct trace_s      *mtr;	/* master traceback tree for alignment */
  struct trace_s      *tr;      /* individual sequence traceback tree  */
  struct trmem_s      *pool;	/* memory pool for traceback tree      */
  struct tracestack_s *dolist;
  struct trace_s      *cur;
  int   *matassign;
  int    nodes;
  int    idx, apos;
  int   *ct;
  int    i,j, nxti, nxtj;

  if (! (ainfo->flags & AINFO_CS)) Die("No cs (consensus structure) line available for that alignment.");

  /* Determine match/insert assignments
   * matassign is 0..alen-1. Values are 1 if MAT, 0 if INS.
   */
  matassign = (int *) MallocOrDie(sizeof(int) * ainfo->alen);
  if (use_rf)
    {
      if (! (ainfo->flags & AINFO_RF)) Die("No rf (reference coord) line available for that alignment.");
      for (apos = 0; apos < ainfo->alen; apos++)
	matassign[apos] = (ainfo->rf[apos] == ' ') ? 0 : 1;
    }
  else
    {
      int gaps;
      for (apos = 0; apos < ainfo->alen; apos++)
	{
	  for (gaps = 0, idx = 0; idx < nseq; idx++)
	    if (isgap(aseq[idx][apos])) gaps++;
	  matassign[apos] = ((double) gaps / (double) nseq > gapthresh) ? 0 : 1;
	}
    }

  /* Determine a "ct" array, base-pairing partners for each position
   */
  if (! KHS2ct(ainfo->cs, ainfo->alen, FALSE, &ct))  Die("Consensus structure string is inconsistent"); 

  /* Make sure the consensus structure "ct" is consistent with the match assignments.
   * Wipe out all structure under INS; including the base-paired 
   * partner of INS-assigned positions
   */
  for (apos = 0; apos < ainfo->alen; apos++)
    if (! matassign[apos])
      { 
	if (ct[apos] != -1)  ct[ct[apos]] = -1;
	ct[apos] = -1;
      }

  /* Construct a master traceback tree.
   * This code is borrowed from yarn's KHS2Trace().
   * mtr's emitl, emitr, and type are properly set by this section.
   */
  InitTrace(&mtr, NULL);
  dolist = InitTracestack();
  cur = AttachTrace(mtr, NULL, 0, ainfo->alen-1, -1, ROOT_NODE);    /* attach the root       */
  PushTracestack(dolist, cur);

  while ((cur = PopTracestack(dolist)) != NULL)
    {
      i = cur->emitl;
      j = cur->emitr;
      
      /* This node accounts for i..j, but we don't know how yet.
       * Six possibilities:
       *    i > j; this is an END state; do nothing.
       *    this is already assigned as a BEGIN; push i,j
       *    i is unpaired; this is a MATL state; push i+1, j
       *    j is unpaired; this is a MATR state; push i,j-1
       *    i,j pair to each other; this is a MATP state; push i+1,j-1
       *    i,j pair but not to each other; this is a BIFURC state;
       *        pick mid ip <= mid < jp; push BEGIN i,mid and working i,mid,
       *        and push BEGIN mid+1,j and working mid+1,j
       */
      if (i > j) cur->type = uEND_ST;

      else if (cur->type == ROOT_NODE)
	{ /* try to push i,j; but deal with INSL and INSR */
	  for (nxti = i; nxti <= j; nxti++)    if (matassign[nxti]) break;
	  for (nxtj = j; nxtj >= nxti; nxtj--) if (matassign[nxtj]) break;
	  if (nxti <= nxtj) PushTracestack(dolist, AttachTrace(cur, NULL, nxti, nxtj, -1, uEND_ST));
	  else { cur->nxtl->emitl = nxti; cur->nxtl->emitr = nxtj; } /* deal with END_ST */             
	}

      else if (cur->type == BEGINL_NODE) /* no inserts */
	{
	  if (i <= j)
	    PushTracestack(dolist, AttachTrace(cur, NULL, i, j, -1, uEND_ST));
	  else
	    { cur->nxtl->emitl = nxti; cur->nxtl->emitr = j; }
	}

      else if (cur->type == BEGINR_NODE) /* INSL */
	{
	  for (nxti = i; nxti <= j; nxti++)    if (matassign[nxti]) break;
	  if (nxti <= j) PushTracestack(dolist, AttachTrace(cur, NULL, nxti, j, -1, uEND_ST));
	  else { cur->nxtl->emitl = nxti; cur->nxtl->emitr = j; } /* deal with END_ST */   
	}

      else if (ct[i] == -1) 	/* i unpaired. This is a MATL node; allow INSL */
	{
	  cur->type = MATL_NODE;
	  for (nxti = i+1; nxti <= j; nxti++)  if (matassign[nxti]) break;
	  if (nxti <= j) PushTracestack(dolist, AttachTrace(cur, NULL, nxti, j, -1, uEND_ST));
	  else { cur->nxtl->emitl = nxti; cur->nxtl->emitr = j; } /* deal with END_ST */
	}

      else if (ct[j] == -1) 	/* j unpaired. MATR node. Deal with INSR */
	{
	  cur->type = MATR_NODE;
	  for (nxtj = j-1; nxtj >= i; nxtj--) if (matassign[nxtj]) break;
	  if (i <= nxtj) PushTracestack(dolist, AttachTrace(cur, NULL, i, nxtj, -1, uEND_ST));
	  else { cur->nxtl->emitl = i; cur->nxtl->emitr = nxtj; } /* deal with END_ST */
	}

      else if (ct[i] == j) 	/* i,j paired to each other. MATP. deal with INSL, INSR */
	{
	  cur->type = MATP_NODE;
	  for (nxti = i+1; nxti <= j; nxti++)    if (matassign[nxti]) break;
	  for (nxtj = j-1; nxtj >= nxti; nxtj--) if (matassign[nxtj]) break;
	  if (nxti <= nxtj) PushTracestack(dolist, AttachTrace(cur, NULL, nxti, nxtj, -1, uEND_ST));
	  else { cur->nxtl->emitl = nxti; cur->nxtl->emitr = nxtj; } /* deal with END_ST */
	}

      else /* i,j paired but not to each other. BIFURC. no INS. */
	{  /* by convention, right side of bifurc deals with insert in middle */
	  cur->type = BIFURC_NODE;
	  PushTracestack(dolist, AttachTrace(cur, NULL, ct[i]+1, j, -1, BEGINR_NODE));
	  PushTracestack(dolist, AttachTrace(cur, NULL, i, ct[i], -1, BEGINL_NODE));
	}
    }	/* while something's on dolist stack */
  FreeTracestack(dolist);
  free(ct);

  /* Now, do the drill for constructing a model using this master trace.
   */
  NumberMasterTrace(mtr, &nodes);
  if ((cm = AllocCM(nodes)) == NULL)
    Die("failed to allocate for new model of %d nodes\n", nodes);
  TopofyNewCM(cm, mtr);
  
  for (idx = 0; idx < nseq; idx++)
    {
      Transmogrify(mtr, aseq[idx], &tr, &pool);
      if (! TraceCount(cm, aseq[idx], 
		       (ainfo->sqinfo[idx].flags & SQINFO_WGT) ? (double) ainfo->sqinfo[idx].weight : 1.0,
		       tr))
	Die("TraceCount() failed");
      FreeTrace(tr, pool);
    }
  ProbifyCM(cm, prior);
  
  free(matassign);
  if (ret_cm  != NULL) *ret_cm  = cm;  else FreeCM(cm);
  if (ret_mtr != NULL) *ret_mtr = mtr; else FreeTrace(mtr, NULL);
}
