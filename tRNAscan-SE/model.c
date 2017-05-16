/* model.c
 * Allocation, initialization, free'ing of models.
 *
 * Includes code supporting both original node-based CM structure, as well
 * as the modified, state-based CM structure used by the newer alignment
 * implementations.
 *
 * SRE, Tue Sep  7 09:22:03 1993
 * 
 */


#include <stdlib.h>
#include <string.h>
#include "structs.h"
#include "funcs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static void fill_state(struct istate_s *st, int nodeidx, int statetype, int offset);
static void copy_singlet_emissions(struct istate_s *st, double *emvec, double *rfreq);
static void copy_pairwise_emissions(struct istate_s *st, double *em, double *rfreq);
static void copy_state_transitions(struct istate_s *st, double *tvec, int tflags);
static void copy_pstate_transitions(struct pstate_s *st, double *tvec, int tflags);
static void fill_pstate(struct pstate_s *st, int nodeidx, int statetype, int offset);

/* Function: AllocCM()
 * 
 * Purpose:  Allocate for a model containing some number of nodes,
 *           inclusive of root but exclusive of ENDs. Blank the model.
 *           
 * Args:     nodes - number of nodes to allocate for
 *                   
 * Return:   pointer to the new model. Caller must free, with FreeCM()
 */
struct cm_s *
AllocCM(int nodes)
{
  struct cm_s *cm;
  int          k;

  if ((cm = (struct cm_s *) malloc (sizeof(struct cm_s))) == NULL)
    Die("malloc failed");

  cm->nodes  = nodes;

  if ((cm->nd = (struct node_s *) malloc (nodes * sizeof(struct node_s))) == NULL)
    Die("malloc failed");

  /* fast way to wipe everything to zero  */
  memset(cm->nd, 0, (nodes * sizeof(struct node_s)));
  
  /* set all the topology connections to -1 */
  for (k = 0; k < nodes; k++)
    cm->nd[k].nxt = cm->nd[k].nxt2 = -1;
      
  return cm;
}
      

/* Function: FreeCM()
 * 
 * Purpose:  Free memory allocated to a covariance model.
 * 
 * Return:   (void)
 */
void
FreeCM(struct cm_s *cm)
{
 free(cm->nd);
 free(cm);
}



/* Function: NormalizeCM()
 * 
 * Purpose:  Normalize all the probability distributions in a model.
 *           Only normalizes the meaningful ones: i.e., matr_emit
 *           emission statistics are ignored for MATL_NODEs, etc.
 *           
 * Return:  (void)
 */
void
NormalizeCM(struct cm_s *cm)
{
  int    k;			/* counter over nodes            */
  int    fy;        		/* from statetype, to statetype  */


  for (k = 0; k < cm->nodes; k++)
    {
      for (fy = 0; fy < STATETYPES; fy++)
	DNorm(cm->nd[k].tmx[fy], STATETYPES);       /* state transitions */
      
      DNorm((double *) cm->nd[k].mp_emit, ALPHASIZE * ALPHASIZE); /* MATP emissions */
      DNorm(cm->nd[k].ml_emit, ALPHASIZE);             /* MATL emissions */
      DNorm(cm->nd[k].mr_emit, ALPHASIZE);             /* MATR emissions */
      DNorm(cm->nd[k].il_emit, ALPHASIZE);             /* INSL emissions */
      DNorm(cm->nd[k].ir_emit, ALPHASIZE);             /* INSR emissions */
    }
}



/* Function: VerifyCM()
 * 
 * Purpose:  Have a look at a CM and make sure nothing stupid
 *           is wrong with it. Returns 0 if something's wrong
 *           and writes diagnostics to stderr. Returns 1 if
 *           everything looks OK.
 */
int
VerifyCM(struct cm_s *cm)
{
  int status = 1;
  int    k;			/* counter over nodes            */

  for (k = 0; k < cm->nodes; k++)
    {
      if (cm->nd[k].type < 0 || cm->nd[k].type >= NODETYPES)
	{
	  status = 0;
	  fprintf(stderr, "Node %d has invalid type %d\n", k, cm->nd[k].type);
	}

      if ((cm->nd[k].nxt <= k && cm->nd[k].nxt != -1) ||
	  (cm->nd[k].nxt >= cm->nodes))
	{
	  status = 0;
	  fprintf(stderr, "Node %d points to invalid left child %d\n", k, cm->nd[k].nxt);
	}

      if ((cm->nd[k].nxt2 <= k && cm->nd[k].nxt2 != -1) ||
	  (cm->nd[k].nxt2 >= cm->nodes))
	{
	  status = 0;
	  fprintf(stderr, "Node %d points to invalid right child %d\n", k, cm->nd[k].nxt2);
	}
    }
  return status;
}


/* Function: RearrangeCM()
 * 
 * Purpose:  Convert a cm into an "integer cm", a specialized structure
 *           used only in the alignment algorithms.
 *
 *           The integer CM is an array of istate_s structures; i.e., rather 
 *           than a node-oriented form, a state-oriented form. Rearrange
 *           transition tables to optimize the recursion in recurse_mx().
 * 
 *           The node is expanded into states in proper order (uDEL_ST,
 *           uMATP_ST, uMATL_ST, uMATR_ST, uINSL_ST, uINSR_ST). However, the 
 *           state transition vectors are rearranged such that INSL, INSR
 *           are the first elements.
 *
 *           Insert states are explicitly assumed to have a zero emission
 *           score!
 *           
 * Args:     cm -      a covariance model, probability form
 *           rfreq  -  frequencies to use as a random model (expected background)
 *           ret_icm - RETURN: array of istate_s structures for states 
 *                     in model. Contains a state 0 for the root;
 *                     does not contain anything for the end
 *           ret_statenum - number of states in ret_icm, 0..statenum-1            
 *           
 * Return:   1 on success, 0 on failure.
 *           ret_icm is malloc'ed here and must be free'ed by caller;
 *           use free(*ret_icm).
 */          
int  
RearrangeCM(struct cm_s      *cm, 
	    double           *rfreq,
	    struct istate_s **ret_icm, 
	    int              *ret_statenum)
{
  struct istate_s *icm;         /* new int-lod states-based model */
  struct istate_s *smallicm;    /* streamlined (realloc'ed) icm   */
  struct m2ali_s  *bifstack;    /* pda for deferring bifurc connection assignment */
  int        bifidx;		/* state index of a BEGINR's parent bifurc        */
  int        k;			/* counter for nodes              */
  int        y;			/* counter for states             */
  int        tflags;		/* flags for which to state transitions are used   */
  int        fflags;		/* flags for which from state transitions are used */
  int        offset;		/* offset to next connected state */

  /* We know we can fit the new model into cm->nodes * STATETYPES 
   * states. We'll give back the excess memory later.
   */
  if ((icm = (struct istate_s *) calloc ((cm->nodes * STATETYPES), sizeof(struct istate_s))) == NULL)
    return 0;
  
  bifstack = Init_m2ali();
  
  y = 0;
  for (k = 0; k < cm->nodes; k++)
    {
				/* figure out what we're connected to */
      if (cm->nd[k].nxt == -1)
	tflags = uEND_ST;
      else
	{
	  switch (cm->nd[cm->nd[k].nxt].type) {
	  case BIFURC_NODE: 
	    tflags = uBIFURC_ST; 
	    break;

	  case MATP_NODE:  
	    tflags = uDEL_ST | uMATP_ST | uMATR_ST | uMATL_ST;
	    break; 

	  case MATL_NODE:   
	    tflags = uDEL_ST | uMATL_ST;
	    break;

	  case MATR_NODE:  
	    tflags = uDEL_ST | uMATR_ST;
	    break;

	  case BEGINL_NODE:
	  case BEGINR_NODE: 
	    tflags = uBEGIN_ST;
	    break;

	  default: Die("no such node type %d", cm->nd[cm->nd[k].nxt].type);
	  }
	}

				/* figure out what we're coming from */
      switch (cm->nd[k].type) {
      case BIFURC_NODE: 
	fflags = uBIFURC_ST;
	offset = 1;
	break;

      case MATP_NODE:
	fflags = uDEL_ST | uMATP_ST | uMATL_ST | uMATR_ST | uINSL_ST | uINSR_ST;
	tflags |= uINSL_ST | uINSR_ST;
	offset = 4;
	break;

      case MATL_NODE:
	fflags = uDEL_ST | uMATL_ST | uINSL_ST;
	tflags |= uINSL_ST;
	offset = 2;
	break;
	
      case MATR_NODE:
	fflags = uDEL_ST | uMATR_ST | uINSR_ST;
	tflags |=  uINSR_ST;
	offset = 2;
	break;

      case BEGINL_NODE:
	fflags = uBEGIN_ST;
	offset = 1;
	break;
	
      case BEGINR_NODE:
	fflags  = uBEGIN_ST | uINSL_ST;
	offset  = 1;
	tflags |= uINSL_ST;
	break;

      case ROOT_NODE:
	fflags = uBEGIN_ST | uINSL_ST | uINSR_ST;
	offset = 1;
	tflags |= uINSL_ST | uINSR_ST;
	break;

      default: Die("No such node type %d\n", cm->nd[k].type);
	
      }

      if (fflags & uDEL_ST)
	{
	  fill_state(&icm[y],   k, uDEL_ST,  offset);
	  copy_state_transitions(&icm[y], cm->nd[k].tmx[DEL_ST], tflags);
	  offset--;
	  y++;
	}
      else if (fflags & uBIFURC_ST)
	{
	  fill_state(&icm[y],   k, uBIFURC_ST,  offset);
	  /* A hack. tmx[0] gets the state index of the left connected BEGIN
	   * child; tmx[1] gets the right connected BEGIN child. The left
           * child is guaranteed to be the next state, but the assignment
           * of the right state must be deferred: we push the bifurc state index
	   * into a pda
	   */
	  icm[y].tmx[0] = y+1;
	  Push_m2ali(bifstack, y, 0, NULL);
	  y++;
	}
      else if (fflags & uBEGIN_ST)
	{
	  fill_state(&icm[y],   k, uBEGIN_ST,  offset);
	  copy_state_transitions(&icm[y], cm->nd[k].tmx[DEL_ST], tflags);

	  /* continuation of the above commentary. If we're a right BEGIN_ST,
	   * then we pop the state index of our parent bifurc off the pda
	   */
	  if (cm->nd[k].type == BEGINR_NODE)
	    {
	      Pop_m2ali(bifstack, &bifidx, (int *) NULL, (struct align_s **) NULL);
	      icm[bifidx].tmx[1] = y;
	    }

	  offset--;
	  y++;
	}

      if (fflags & uMATP_ST)
	{
	  fill_state(&icm[y], k, uMATP_ST, offset);
	  copy_pairwise_emissions(&icm[y], (double *) cm->nd[k].mp_emit, rfreq);
	  copy_state_transitions(&icm[y], cm->nd[k].tmx[MATP_ST], tflags);
	  offset--;
	  y++;
	}

      if (fflags & uMATL_ST)
	{
	  fill_state(&icm[y], k, uMATL_ST, offset);
	  copy_singlet_emissions(&icm[y], cm->nd[k].ml_emit, rfreq);
	  copy_state_transitions(&icm[y], cm->nd[k].tmx[MATL_ST], tflags);
	  offset--;
	  y++;
	}

      if (fflags & uMATR_ST)
	{
	  fill_state(&icm[y], k, uMATR_ST, offset);
	  copy_singlet_emissions(&icm[y], cm->nd[k].mr_emit, rfreq);
	  copy_state_transitions(&icm[y], cm->nd[k].tmx[MATR_ST], tflags);
	  offset--;
	  y++;
	}

      if (fflags & uINSL_ST)
	{
	  fill_state(&icm[y], k, uINSL_ST, 0);
	  copy_singlet_emissions(&icm[y], cm->nd[k].il_emit, rfreq);
	  copy_state_transitions(&icm[y], cm->nd[k].tmx[INSL_ST], tflags);
	  y++;
	}

      if (fflags & uINSR_ST)
	{
				/* beware an asymmetry: INSR->INSL transits are disallowed */
	  fill_state(&icm[y], k, uINSR_ST, 0);
	  copy_singlet_emissions(&icm[y], cm->nd[k].ir_emit, rfreq);
	  copy_state_transitions(&icm[y], cm->nd[k].tmx[INSR_ST], tflags & ~uINSL_ST);
	  y++;
	}
      

      /* End states must be added
       */
      if  (cm->nd[k].nxt == -1)
	{
	  fill_state(&icm[y], -1, uEND_ST, 0);
	  y++;
	}
    } /* end loop over nodes */


  Free_m2ali(bifstack);
  
  /* Return some of the alloc'ed memory
   */
  smallicm = (struct istate_s *) realloc (icm, y * sizeof(struct istate_s));
  *ret_icm = (smallicm != NULL) ? smallicm : icm;
  *ret_statenum = y;
  return 1;
}



/* Function: fill_state()
 * 
 * Purpose:  fill in values in a state: node index, state type unique
 *           identifier, offset to first of the next states, and number 
 *           of ynext connections.
 *           
 *           transition and emission probabilities are dealt with
 *           elsewhere.
 *           
 */
static void
fill_state(struct istate_s *st,
	   int              nodeidx,
	   int              statetype,
	   int              offset)
{
  st->nodeidx    = nodeidx;
  st->statetype  = statetype;
  st->offset     = offset;
}



/* Function: copy_singlet_emissions()
 * 
 * Purpose:  Copy a singlet emission vector into a state structure,
 *           converting the probabilities into integer log odds.
 */
static void
copy_singlet_emissions(struct istate_s *st, double *emvec, double *rfreq)
{
  int x;

  for (x = 0; x < ALPHASIZE; x++)
    st->emit[x] = ILOG2(emvec[x] / rfreq[x]);
}




/* Function: copy_pairwise_emissions()
 * 
 * Purpose:  Copy a pairwise emission table into a state structure,
 *           converting the probabilities into integer log odds.
 *           Beware the funny business with the pairwise emission
 *           array; it was mp_emit[4][4], now cast to a pointer,
 *           and accessed like a vector.
 */
static void
copy_pairwise_emissions(struct istate_s *st, double *em, double *rfreq)
{
  int x; 

  for (x = 0; x < ALPHASIZE * ALPHASIZE; x++)
    st->emit[x] = ILOG2(em[x] / (rfreq[x % ALPHASIZE] * rfreq[x / ALPHASIZE]));
}




/* Function: copy_state_transitions()
 * 
 * Purpose:  Copy a state transition vector from a CM into a state
 *           structure, copying only the used state transitions
 *           as given by tflags. The state transition vector is
 *           rearranged for an optimization: transits to INSL, INSR 
 *           are placed first.
 */
static void
copy_state_transitions(struct istate_s *st,
		       double          *tvec,
		       int              tflags)
{
  int stx;		/* counter for state vector */

  stx = 0; 
  if (tflags & uINSL_ST)
    { st->tmx[stx] = ILOG2(tvec[INSL_ST]);  stx++; }
  
  if (tflags & uINSR_ST)
    { st->tmx[stx] = ILOG2(tvec[INSR_ST]);  stx++; }

  if (tflags & uDEL_ST || tflags & uBIFURC_ST ||
      tflags & uBEGIN_ST || tflags & uEND_ST)
    { st->tmx[stx] = ILOG2(tvec[DEL_ST]);   stx++; }
      
  if (tflags & uMATP_ST)
    { st->tmx[stx] = ILOG2(tvec[MATP_ST]);  stx++; }

  if (tflags & uMATL_ST)
    { st->tmx[stx] = ILOG2(tvec[MATL_ST]);  stx++; }

  if (tflags & uMATR_ST)
    { st->tmx[stx] = ILOG2(tvec[MATR_ST]);  stx++; }

  st->connectnum = stx;
}    


/* Function: MakePCM()
 * 
 * Purpose:  Like RearrangeCM(), but leaving the model
 *           in floating-point probabilities in struct pstate_s
 *           structures.
 *           
 * Args:     cm -      a covariance model, probability form
 *           ret_pcm - RETURN: array of pstate_s structures for states 
 *                     in model. Contains a state 0 for the root.
 *                     end states are explicit.
 *           ret_statenum - number of states in ret_pcm, 0..statenum-1            
 *           
 * Return:   1 on success, 0 on failure.
 *           ret_pcm is malloc'ed here and must be free'ed by caller;
 *           use free(*ret_pcm).
 */          
int  
MakePCM(struct cm_s      *cm, 
	struct pstate_s **ret_pcm, 
	int              *ret_statenum)
{
  struct pstate_s *pcm;         /* new states-based model */
  struct pstate_s *smallpcm;    /* streamlined (realloc'ed) pcm   */
  struct intstack_s *bifstack;  /* pda for deferring bifurc connection assignment */
  int        bifidx;		/* state index of a BEGINR's parent bifurc        */
  int        k;			/* counter for nodes              */
  int        y;			/* counter for states             */
  int        tflags;		/* flags for which to state transitions are used   */
  int        fflags;		/* flags for which from state transitions are used */
  int        offset;		/* offset to next connected state */

  /* We know we can fit the new model into cm->nodes * STATETYPES 
   * states. We'll give back the excess memory later.
   */
  if ((pcm = (struct pstate_s *) calloc ((cm->nodes * STATETYPES), sizeof(struct pstate_s))) == NULL)
    return 0;

  bifstack = InitIntStack();
  
  y = 0;
  for (k = 0; k < cm->nodes; k++)
    {
				/* figure out what we're connected to */
      if (cm->nd[k].nxt == -1)
	tflags = uEND_ST;
      else
	{
	  switch (cm->nd[cm->nd[k].nxt].type) {
	  case BIFURC_NODE: tflags = uBIFURC_ST; 	                       break;
	  case MATP_NODE:   tflags = uDEL_ST | uMATP_ST | uMATR_ST | uMATL_ST; break;
	  case MATL_NODE:   tflags = uDEL_ST | uMATL_ST;                       break;
	  case MATR_NODE:   tflags = uDEL_ST | uMATR_ST;                       break;
	  case BEGINL_NODE: tflags = uBEGIN_ST;                                break;
	  case BEGINR_NODE: tflags = uBEGIN_ST;                                break;
	  default: Die("no such node type %d", cm->nd[cm->nd[k].nxt].type);
	  }
	}

				/* figure out what we're coming from */
      switch (cm->nd[k].type) {
      case BIFURC_NODE: 
	fflags = uBIFURC_ST;
	offset = 1;
	break;

      case MATP_NODE:
	fflags = uDEL_ST | uMATP_ST | uMATL_ST | uMATR_ST | uINSL_ST | uINSR_ST;
	tflags |= uINSL_ST | uINSR_ST;
	offset = 4;
	break;

      case MATL_NODE:
	fflags = uDEL_ST | uMATL_ST | uINSL_ST;
	tflags |= uINSL_ST;
	offset = 2;
	break;
	
      case MATR_NODE:
	fflags = uDEL_ST | uMATR_ST | uINSR_ST;
	tflags |=  uINSR_ST;
	offset = 2;
	break;

      case BEGINL_NODE:
	fflags = uBEGIN_ST;
	offset = 1;
	break;
	
      case BEGINR_NODE:
	fflags  = uBEGIN_ST | uINSL_ST;
	offset  = 1;
	tflags |= uINSL_ST;
	break;

      case ROOT_NODE:
	fflags = uBEGIN_ST | uINSL_ST | uINSR_ST;
	offset = 1;
	tflags |= uINSL_ST | uINSR_ST;
	break;

      default: Die("No such node type %d\n", cm->nd[k].type);
	
      }

      if (fflags & uDEL_ST)
	{
	  fill_pstate(&pcm[y],   k, uDEL_ST,  offset);
	  copy_pstate_transitions(&pcm[y], cm->nd[k].tmx[DEL_ST], tflags);
	  offset--;
	  y++;
	}
      else if (fflags & uBIFURC_ST)
	{
	  fill_pstate(&pcm[y],   k, uBIFURC_ST,  offset);
				/* We defer the assignment of bifr */
	  PushIntStack(bifstack, y);
	  y++;
	}
      else if (fflags & uBEGIN_ST)
	{
	  fill_pstate(&pcm[y],   k, uBEGIN_ST,  offset);
	  copy_pstate_transitions(&pcm[y], cm->nd[k].tmx[DEL_ST], tflags);

	  /* continuation of the above commentary. If we're a right BEGIN_ST,
	   * then we pop the state index of our parent bifurc off the pda
	   */
	  if (cm->nd[k].type == BEGINR_NODE)
	    {
	      PopIntStack(bifstack, &bifidx);
	      pcm[bifidx].bifr = y;
	    }
	  offset--;
	  y++;
	}

      if (fflags & uMATP_ST)
	{
	  fill_pstate(&pcm[y], k, uMATP_ST, offset);
	  memcpy(pcm[y].emit, cm->nd[k].mp_emit, sizeof(double) * ALPHASIZE * ALPHASIZE);
	  copy_pstate_transitions(&pcm[y], cm->nd[k].tmx[MATP_ST], tflags);
	  offset--;
	  y++;
	}

      if (fflags & uMATL_ST)
	{
	  fill_pstate(&pcm[y], k, uMATL_ST, offset);
	  memcpy(pcm[y].emit, cm->nd[k].ml_emit, sizeof(double) * ALPHASIZE);
	  copy_pstate_transitions(&pcm[y], cm->nd[k].tmx[MATL_ST], tflags);
	  offset--;
	  y++;
	}

      if (fflags & uMATR_ST)
	{
	  fill_pstate(&pcm[y], k, uMATR_ST, offset);
	  memcpy(pcm[y].emit, cm->nd[k].mr_emit, sizeof(double) * ALPHASIZE);
	  copy_pstate_transitions(&pcm[y], cm->nd[k].tmx[MATR_ST], tflags);
	  offset--;
	  y++;
	}

      if (fflags & uINSL_ST)
	{
	  fill_pstate(&pcm[y], k, uINSL_ST, 0);
	  memcpy(pcm[y].emit, cm->nd[k].il_emit, sizeof(double) * ALPHASIZE);
	  copy_pstate_transitions(&pcm[y], cm->nd[k].tmx[INSL_ST], tflags);
	  y++;
	}

      if (fflags & uINSR_ST)
	{
				/* beware an asymmetry: INSR->INSL transits are disallowed */
	  fill_pstate(&pcm[y], k, uINSR_ST, 0);
	  memcpy(pcm[y].emit, cm->nd[k].ir_emit, sizeof(double) * ALPHASIZE);
	  copy_pstate_transitions(&pcm[y], cm->nd[k].tmx[INSR_ST], tflags & ~uINSL_ST);
	  y++;
	}
      

      /* End states must be added
       */
      if  (cm->nd[k].nxt == -1)
	{
	  fill_pstate(&pcm[y], -1, uEND_ST, 0);
	  y++;
	}
    } /* end loop over nodes */


  FreeIntStack(bifstack);
  
  /* Return some of the alloc'ed memory
   */
  smallpcm = (struct pstate_s *) realloc (pcm, y * sizeof(struct pstate_s));
  *ret_pcm = (smallpcm != NULL) ? smallpcm : pcm;
  *ret_statenum = y;
  return 1;
}




/* Function: copy_pstate_transitions()
 * 
 * Purpose:  Copy a state transition vector from a CM into a state
 *           structure, copying only the used state transitions
 *           as given by tflags. The state transition vector is
 *           rearranged for an optimization: transits to INSL, INSR 
 *           are placed first.
 */
static void
copy_pstate_transitions(struct pstate_s *st,
		       double          *tvec,
		       int              tflags)
{
  int stx;		/* counter for state vector */

  stx = 0; 
  if (tflags & uINSL_ST)
    { st->tmx[stx] = tvec[INSL_ST];  stx++; }
  
  if (tflags & uINSR_ST)
    { st->tmx[stx] = tvec[INSR_ST];  stx++; }

  if (tflags & uDEL_ST || tflags & uBIFURC_ST ||
      tflags & uBEGIN_ST || tflags & uEND_ST)
    { st->tmx[stx] = tvec[DEL_ST];   stx++; }
      
  if (tflags & uMATP_ST)
    { st->tmx[stx] = tvec[MATP_ST];  stx++; }

  if (tflags & uMATL_ST)
    { st->tmx[stx] = tvec[MATL_ST];  stx++; }

  if (tflags & uMATR_ST)
    { st->tmx[stx] = tvec[MATR_ST];  stx++; }

  st->connectnum = stx;
}    


/* Function: fill_pstate()
 * 
 * Purpose:  fill in values in a state: node index, state type unique
 *           identifier, offset to first of the next states.
 *           
 *           transition and emission probabilities are dealt with
 *           elsewhere.
 *           
 */
static void
fill_pstate(struct pstate_s *st,
	    int              nodeidx,
	    int              statetype,
	    int              offset)
{
  st->nodeidx    = nodeidx;
  st->statetype  = statetype;
  st->offset     = offset;
}


/* Function: NormalizePCM()
 * 
 * Purpose:  Make damn sure a probability-form, states-based CM is
 *           properly normalized. Workaround for a bug!
 */
void
NormalizePCM(struct pstate_s *pcm, int M)
{
  int    y;

  for (y = 0; y < M; y++)
    {
				/* emission distributions */
      switch (pcm[y].statetype) {
      case uMATP_ST: DNorm(pcm[y].emit, ALPHASIZE * ALPHASIZE); break;
      case uMATL_ST:
      case uMATR_ST:
      case uINSL_ST:
      case uINSR_ST: DNorm(pcm[y].emit, ALPHASIZE); break;
      }

				/* transition distributions */
      if (pcm[y].statetype != uBIFURC_ST)
	DNorm(pcm[y].tmx, pcm[y].connectnum);
    }
}
