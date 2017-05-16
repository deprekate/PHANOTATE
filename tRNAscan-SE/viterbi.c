/* viterbi.c
 * for 2.0:  SRE, Tue Sep 28 09:15:07 1993
 * from 1.0: SRE, Wed Jun 30 17:42:44 1993
 *           as revised Tue Aug 24 12:06:12 1993: integer two-matrix version
 * 
 * Implementation of the three-dimensional dynamic programming
 * algorithm for aligning a covariance model to a sequence.
 *
 * To optimize memory access patterns, the score storage is implemented
 * as a two-matrix version. amx is the
 * main storage. bmx is a smaller auxiliary matrix with a different
 * access pattern, holding scores of BEGIN state alignments; it
 * is used when calculating BIFURC scores.
 *
 * amx is [j = 0..N] [diff = 0..j] [y = 0..statenum]
 *  diff == 0 is for off-diagonal boundary conditions (this is why diff is shifted +1)
 *  diff == 1 is for the diagonal, i==j
 *
 * bmx is [y = 0..statenum] [j = 0..N] [ diff = 0..j]
 *   a j,diff matrix exists only where y is a BEGIN state
 *
 * The 2.0 implementation allows variable storage per node rather
 * than storing and calculating a fixed max number of states per node, 
 * which should save up to 2x in both time and space.
 *
 * An optimization is made which requires END states to be explicitly
 * added, so statenum (the number of states in the integer model)
 * is *inclusive* of ENDs.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "squid.h"
#include "structs.h"
#include "funcs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif
#ifdef DEBUG
#include <assert.h>
#endif

static int  allocate_mx(struct istate_s *icm,int statenum, int seqlen,	
			int ****ret_amx, int ****ret_bmx);
static int  init_mx(struct istate_s *icm, int statenum, int N, int ***amx, int ***bmx);
static int  recurse_mx(struct istate_s *icm, int statenum, char *seq, int N, int ***amx, int ***bmx);
static int  trace_mx(struct istate_s *icm, char *seq, int N, 
			 int ***amx, int ***bmx, struct trace_s **ret_tr);
static void free_mx(int ***amx, int ***bmx, int statenum, int seqlen);


/* Function: ViterbiAlign()
 * 
 * Purpose:  Align a sequence to a model, using the alignment
 *           algorithm. Return the score of the alignment and
 *           the traceback.
 *
 * Args:     icm       - the model to align sequence to
 *           statenum  - # of states in the model
 *           seq       - sequence to align model to
 *           ret_score - RETURN: global alignment score
 *           ret_trace - RETURN: traceback tree
 *                       
 * Return:   1 on success, 0 on failure.                      
 */       
int
ViterbiAlign(struct istate_s *icm,       
	     int              statenum,
	     char            *seq,      
	     double          *ret_score,
	     struct trace_s **ret_trace)
{
  int ***amx;			/* the main score matrix    */
  int ***bmx;                   /* the BEGIN score matrix   */
  int    N;			/* length of sequence       */

  N = strlen(seq);
  seq--;			/* convert to 1..N. Ugh! */

  if (! allocate_mx(icm, statenum, N, &amx, &bmx)) return 0;
#ifdef DEBUG
  printf("allocated matrices\n");
#endif

  if (! init_mx(icm, statenum, N, amx, bmx)) return 0;
#ifdef DEBUG
  printf("matrices initialized\n");
#endif
  
  if (! recurse_mx(icm, statenum, seq, N, amx, bmx)) return 0;
#ifdef DEBUG
  printf("recursion finished\n");
  PrintViterbiAMX(stdout, icm, statenum, seq, N, amx);
#endif

  *ret_score = ((double) bmx[0][N][N] / INTPRECISION);
#ifdef DEBUG
  printf("have a score of %.2f, starting traceback\n", *ret_score);
#endif

  if (! trace_mx(icm, seq, N, amx, bmx, ret_trace)) return 0;
#ifdef DEBUG
  printf("trace complete\n");
  PrintTrace(stdout, *ret_trace);
#endif
  free_mx(amx, bmx, statenum, N);

#ifdef SRE_REMOVED
  /* OK, boys, crank up the MasPar DPU.
   * 24 bytes = 6 addresses passed * 4 bytes/address (always right?)
   * MPViterbiAlign() copies ret_score out, but doesn't do anything
   * with the traceback yet
   */
  callRequest(MPViterbiAlign, 24, icm, &statenum, seq, &N, ret_score, ret_trace);
#endif /* MASPAR */

  return 1;
}

/* Function: allocate_cvmx()
 * 
 * Purpose:  Malloc space for the score matrices.
 *           amx is indexed as j, i, y.
 *           bmx is indexed as k, j, i.      
 *           In the two sequence dimensions j, i they are
 *           diagonal (+1 off diagonal) matrices with
 *           rows j = 0..N, i = 1..j+1.
 *           In the node dimension k bmx is k = 0..M.
 *           In the state dimension y amx is y = 0..numstates.
 *           
 * Args:     icm      - the int, log-odds, state-based model
 *           statenum - number of states in model
 *           seqlen   - length of sequence
 *           ret_amx  - RETURN: main score matrix
 *           ret_bmx  - RETURN: BEGIN score matrix
 * 
 * Return:   Ptr to allocated scoring matrix, or
 *           dies and exits.
 */
static int 
allocate_mx(struct istate_s *icm,
	    int     statenum,
	    int     seqlen,	
	    int ****ret_amx,
	    int ****ret_bmx)	
{
  int  ***amx;
  int  ***bmx;
  int     diag, j, y;


  /* Main matrix, amx: fastest varying index is y (j,i,y)
   */
				/* malloc for j = 0..seqlen rows */
  if ((amx = (int ***) malloc ((seqlen + 1)* sizeof(int **))) == NULL)
    Die("Memory allocation error in %s line %d", __FILE__, __LINE__);

  for (j = 0; j <= seqlen; j++)	/* loop over rows j = 0..N */
    {
				/* malloc for diag = 0..j cols */
      if ((amx[j] = (int **) malloc ((j + 1) * sizeof(int *))) == NULL)
	Die("Memory allocation error in %s line %d", __FILE__, __LINE__);
      
				/* loop over cols diag = 0..j */
      for (diag = 0; diag <= j; diag++)
				/* malloc for y = 0..statenum-1 decks */
	  if ((amx[j][diag] = (int *) malloc ((statenum) * sizeof (int))) == NULL)
	    Die("Memory allocation error in %s line %d", __FILE__, __LINE__);
    }

  
  /* BEGIN auxiliary matrix: fastest varying index is diag (y,j,diag)
   */
				/* 0..statenum-1 decks */
  if ((bmx = (int ***) malloc (statenum * sizeof(int **))) == NULL)
    Die("Memory allocation error in %s line %d", __FILE__, __LINE__);

  for (y = 0; y < statenum; y++)
    {
      bmx[y] = NULL;
      if (icm[y].statetype == uBEGIN_ST)
	{
	  if ((bmx[y] = (int **) malloc ((seqlen+1) * sizeof(int *))) == NULL)
	    Die("Memory allocation error in %s line %d", __FILE__, __LINE__);

	  for (j = 0; j <= seqlen; j++)
	    if ((bmx[y][j] = (int *) malloc ((j+1) * sizeof(int))) == NULL)
	      Die("Memory allocation error in %s line %d", __FILE__, __LINE__);
	}
    }

  *ret_amx = amx;
  *ret_bmx = bmx;
  return 1;
}



/* Function: free_mx()
 * 
 * Purpose:  Free the space allocated to the two scoring matrices.
 *           Precisely mirrors the allocations above in allocate_cvmx().
 * 
 * Return:   (void)
 */
static void
free_mx(int  ***amx,
	int  ***bmx,
	int     statenum,
	int     seqlen)
{
  int diag, j, y;

  /* Free the main matrix, amx
   */
  for (j = 0; j <= seqlen; j++)
    {
      for (diag = 0; diag <= j; diag++)
	free(amx[j][diag]);
      free(amx[j]);
    }
  free(amx);

  /* Free the auxiliary BEGIN matrix, bmx
   */
  for (y = 0; y < statenum; y++)
    if (bmx[y] != NULL)
      {
	for (j = 0; j <= seqlen; j++)
	  free(bmx[y][j]);
	free(bmx[y]);
      }
  free(bmx);
}



/* Function: init_mx()
 * 
 * Purpose:  Initialization of the scoring matrices. We initialize the off-diagonal,
 *           the diagonal, and the "floor" (end states) of the cube.
 * 
 * Return:   1 on success, 0 on failure.
 */
static int
init_mx(struct istate_s *icm,          /* integer model  */
	int              statenum,     /* number of states in icm */
	int              N,	       /* length of seq  */
	int           ***amx,
	int           ***bmx)
{
  int diag, j, y;		/* counters for indices over the cvmx            */
  int ynext;			/* index of next state k+1                       */
  int *beam;                    /* z-axis vector of numbers in amx               */

  /* Init the whole amx to -Infinity. We do this with memcpy, trying
   * to be fast. We fill in j=0,diag=0 by hand, then memcpy() the other
   * columns.
   */
  for (y = 0; y < statenum; y++)
    amx[0][0][y] = NEGINFINITY;
  for (j = 1; j <= N; j++)
    for (diag = 0; diag <= j; diag++)
      memcpy(amx[j][diag], amx[0][0], statenum * sizeof(int));

  /* Init the whole bmx to -Inf. We know state 0 is a begin (it's ROOT), so we
   * start there, and memcpy rows as needed.
   */
  for (diag = 0; diag <= N; diag++)  bmx[0][N][diag] = NEGINFINITY;
  for (j = 0; j < N; j++)
    memcpy(bmx[0][j], bmx[0][N], (j+1) * sizeof(int));

  for (y = 1; y < statenum; y++)
    if (bmx[y] != NULL)
      for (j = 0; j <= N; j++)
	memcpy(bmx[y][j], bmx[0][N], (j+1) * sizeof(int));
  
  /* Init the off-diagonal (j = 0..N; diag == 0) with -log P scores.
   * End state = 0;
   * del, bifurc states are calc'ed
   * begin states same as del's
   */
  for (j = 0; j <= N; j++)
    for (y = statenum-1; y >= 0; y--)
      {
	if (icm[y].statetype == uEND_ST)
	  amx[j][0][y] = 0; 

	else if (icm[y].statetype == uBIFURC_ST)
	  amx[j][0][y] = bmx[icm[y].tmx[0]][j][0] + bmx[icm[y].tmx[1]][j][0];

	else if (icm[y].statetype == uDEL_ST || icm[y].statetype == uBEGIN_ST)
	  {
				/* only calc DEL-DEL and BEGIN-DEL transitions. Since
				 * we optimized the state transition tables, removing
				 * the unused ones, we don't know where the number
				 * for "to DEL" is! But we can find it, because it'll
				 * be the connection to a non-infinite score */
	    beam = amx[j][0] + y + icm[y].offset;
	    for (ynext = 0; ynext < icm[y].connectnum; ynext++)
	      {
		if (*beam != NEGINFINITY)
		  amx[j][0][y] = *beam + icm[y].tmx[ynext];
		beam++;
	      }

				/* make a copy into bmx if y is a BEGIN */
	    if (icm[y].statetype == uBEGIN_ST)
	      bmx[y][j][0] = amx[j][0][y];
	  }
      }

  return 1;
}



/* Function: recurse_mx()
 * 
 * Purpose:  Carry out the fill stage of the dynamic programming
 *           algorithm.
 *           
 * Returns:  1 on success, 0 on failure.
 */
static int
recurse_mx(struct istate_s *icm,      /* integer, state-form model */
	   int              statenum, /* number of states in icm   */
	   char            *seq,      /* sequence, 1..N            */
	   int              N,	      /* length of seq             */
	   int           ***amx,      /* main scoring matrix       */
	   int           ***bmx)      /* bifurc scoring matrix     */
{
  int i, j, y;		        /* indices for 4 dimensions                    */
  int diff;			/* loop counter for difference: diff = j-i + 1 */
  int symi, symj;		/* symbol indices for seq[i], seq[j]           */
  int sc;			/* tmp for a score                             */
  int ynext;			/* index of next state y                       */

  int *beam;                    /* ptr to a beam (z-axis vector)               */
  int  leftdiff;		/* diff coord of BEGIN_L of a bifurc     */
  int  leftj;			/* j coord of BEGIN_L of a bifurc        */
  int **left_p;			/* pointer into whole 2D deck of BEGINL's of a bifurc */
  int *right_p;                 /* ptr into row of BEGIN_R's of a bifurc */
  int *scp;			/* score pointer: ptr into beam of scores being calc'ed */
  struct istate_s *st;		/* state pointer: ptr at current state in icm */
  int *tmx;
  int  emitsc;

  for (j = 1; j <= N; j++)
    {
	symj = SymbolIndex(seq[j]);      
	for (diff = 1; diff <= j; diff++)
	  {
	    i = j - diff + 1;
	    if (i < 1) break;

	    symi = SymbolIndex(seq[i]);

#ifdef DEBUG
	    assert(symi >= 0 && symi < ALPHASIZE);
	    assert(symj >= 0 && symj < ALPHASIZE);
#endif

	    scp = &amx[j][diff][statenum-1];
	    st  = &icm[statenum-1];
	    for (y = statenum-1; y >= 0; y--, scp--, st--)
	      { /* loop over states */

		if (st->statetype != uBIFURC_ST)	/* a normal (non-BIFURC) state */
		  {
		    /* Connect the "beam" pointer to the appropriate
		     * starting place in the ynext scores we're connecting
		     * y to
		     */
		    switch (st->statetype) {
		    case uBEGIN_ST:
		    case uDEL_ST:
		      beam   = amx[j][diff];    
		      emitsc = 0;    
		      break;
		    case uMATP_ST:  
		      if (diff == 1) continue;
		      beam   = amx[j-1][diff-2]; 
		      emitsc = st->emit[symi * ALPHASIZE + symj];
		      break; 
		    case uMATR_ST:
		    case uINSR_ST:
		      beam   = amx[j-1][diff-1];
		      emitsc = st->emit[symj];
		      break;
		    case uMATL_ST:
		    case uINSL_ST:  
		      beam   = amx[j][diff-1];   
		      emitsc = st->emit[symi];
		      break;
		    case uEND_ST:   continue;
		    default: Die("no such state type %d", st->statetype);
		    }
		    beam += y + st->offset;
		    tmx  = st->tmx;


		    /* Init for ynext == 0 case 
		     */		
		    *scp = *beam + *tmx;
		    
		    /* Calculate remaining cases
		     */
		    for (ynext = 1; ynext < st->connectnum; ynext++)
		      {
			beam++;
			tmx++;
			if (*beam > *scp)
			  {
			    sc = *beam + *tmx;
			    if (sc > *scp)
			      *scp = sc;
			  }
		      }
		    
		    /* Add emission scores now
		     */
		    *scp += emitsc;
		    
		    /* Make a copy into bmx if necessary
		     */
		    if (st->statetype == uBEGIN_ST)
		      bmx[y][j][diff] = *scp;
		  } /* end block of normal state stuff */
		
		else		/* a BIFURC state */
		  {
		    leftdiff = diff;
		    leftj    = j;
		    right_p  = bmx[st->tmx[1]][j];
		    left_p   = bmx[st->tmx[0]];
		    
		    /* init w/ case that left branch emits it all */
		    *scp = left_p[leftj][leftdiff] + *right_p;
		    while (leftdiff > 0)
		      {
			leftdiff--;
			leftj--;
			right_p++;
			
			sc = left_p[leftj][leftdiff] + *right_p;
			if (sc > *scp)
			  *scp = sc;
		      }
		  }
		
	      } /* end loop over states */
	  } /* end loop over diff */
      } /* end loop over j */
  return 1;
}


/* Function: trace_cvmx()
 * 
 * Purpose:  Trace stage of the dynamic programming: starting
 *           at j=N, i=1, k=0/BEGIN, trace back the optimal
 *           path. Returns a binary tree, ret_trace.
 *           Caller is reponsible for free'ing ret_trace.
 */
static int
trace_mx(struct istate_s *icm,       /* the model to align               */   
	 char            *seq,       /* sequence to align it to  1..N    */
	 int              N,
	 int           ***amx,      
	 int           ***bmx,
	 struct trace_s **ret_trace) /* RETURN: the traceback tree       */
{
  struct trace_s *tr;           /* the traceback tree under construction */
  struct trace_s *curr_tr;      /* ptr to node of tr we're working on    */
  struct tracestack_s *dolist;  /* pushdown stack of active tr nodes     */
  int diff,i, j;		/* coords in mx (0..N)                   */
  int y;			/* counter for states (0..statenum-1)    */
  int ynext;			/* holds "k+1" value                     */
  int symi, symj;		/* array indices for left, right symbols */
  int leftdiff;
  int leftj;
  int *right_p;
  int *beam;
  int  conni, connj;
  int  sc;

  /* Initialize.
   * Start at i = 1, j = N and work towards diagonal
   */
  InitTrace(&tr, NULL);         /* start a trace tree */
  dolist = InitTracestack();	/* start a stack for traversing the trace tree */

  curr_tr = AttachTrace(tr, NULL, 0, N-1, 0, uBEGIN_ST);
  PushTracestack(dolist, curr_tr);

  /* Recursion. While there's active nodes in the stack, trace from them.
   * 
   * This is cribbed from recurse_cvmx(); it's almost the exact reverse.
   * We know the best score, we just have to figure out where it came from.
   */
  while ((curr_tr = PopTracestack(dolist)) != NULL)
    {
				/* get some useful numbers, mostly for clarity */
				/* which is important, since we're sort of misusing
				 * fields in the trace structures! */
      i    = curr_tr->emitl+1;
      j    = curr_tr->emitr+1;
      y    = curr_tr->nodeidx;
      diff = j - i + 1;

      /* During use here, nodeidx field is used to hold a *state* index,
       * when we leave, everything must look like the rest of the package
       * expects, so we clean up here.
       */
      curr_tr->nodeidx   = icm[y].nodeidx;

				/* We used an END state here.
				 * (We'd better be near the diagonal!)
				 * We're done here. */
      if (icm[y].statetype == uEND_ST)
	{
	  if (i <= j) Warn("trace: didn't reach off-diag, stop at i=%d j=%d y=%d", i,j,y);
	  curr_tr->nodeidx = -1;
	  continue;
	}

      else if (icm[y].statetype == uBIFURC_ST) /* bifurc state */
	{
				/* We used a BIFURC state here. 
				 * It came from two branches. Redo the recurse_cvmx()
				 * calculation to find them. */
	  if (i > j)
	    {
	      PushTracestack(dolist, AttachTrace(curr_tr, NULL, i-1, j-1, icm[y].tmx[1], uBEGIN_ST));
	      PushTracestack(dolist, AttachTrace(curr_tr, NULL, i-1, j-1, icm[y].tmx[0], uBEGIN_ST));
	    }

	  else
	    {
	      leftdiff = diff;
	      leftj    = j;
	      right_p  = bmx[icm[y].tmx[1]] [j];
	      
	      while (leftdiff >= 0)
		{
		  if (amx[j][diff][y] == bmx[icm[y].tmx[0]][leftj][leftdiff] + *right_p)
		    {
		      PushTracestack(dolist, AttachTrace(curr_tr, NULL, i + leftdiff-1, j-1, 
							 icm[y].tmx[1], uBEGIN_ST));
		      PushTracestack(dolist, AttachTrace(curr_tr, NULL, i -1, i+leftdiff-2, 
							 icm[y].tmx[0], uBEGIN_ST));
		      break;
		    }
		  leftdiff--;
		  leftj--;
		  right_p++;
		}
	      if (leftdiff < 0)
		Die("bifurc reconstruction failed at ijy %d,%d,%d", i,j,y);
	    }
	}
      
      else			/* a normal (non-BIFURC) state */
	{
	  if (i > 0 && i <= N) symi  = SymbolIndex(seq[i]);
	  if (j > 0 && j <= N) symj  = SymbolIndex(seq[j]);

	  switch (icm[y].statetype) {
	  case uBEGIN_ST:
	  case uDEL_ST:   beam = amx[j][diff];     conni = i;   connj = j;   break;
	  case uMATP_ST:  beam = amx[j-1][diff-2]; conni = i+1; connj = j-1; break;
	  case uMATR_ST:
	  case uINSR_ST:  beam = amx[j-1][diff-1]; conni = i;   connj = j-1; break;
	  case uMATL_ST:
	  case uINSL_ST:  beam = amx[j][diff-1];   conni = i+1; connj = j;   break;
	  default: Die("no such state type %d", icm[y].statetype);
	  }
	  beam += y + icm[y].offset;

				/* Calculate the score we'll try to match,
				   by subtracting emission score as needed */
	  sc = amx[j][diff][y];
	  switch (icm[y].statetype) {
	  case uBEGIN_ST:
	  case uDEL_ST:   break;
	  case uMATP_ST:  sc -= icm[y].emit[symi * ALPHASIZE + symj]; break;
	  case uMATR_ST: 
	  case uINSR_ST:  sc -= icm[y].emit[symj]; break;
	  case uMATL_ST:
	  case uINSL_ST:  sc -= icm[y].emit[symi]; break; 
	  default: Die("no such state type %d", icm[y].statetype);
	  }

				/* find the right connection */
	  for (ynext = 0; ynext < icm[y].connectnum; ynext++, beam++)
	    if (sc == *beam + icm[y].tmx[ynext])
	      {
		PushTracestack(dolist, AttachTrace(curr_tr, NULL, conni-1, connj-1, 
						   ynext + y + icm[y].offset, 
						   icm[ynext + y + icm[y].offset].statetype));
		break;
	      }
	  if (ynext == icm[y].connectnum)
	    { Warn("can't continue traceback"); return 0; }
	      
	} /* (a normal statetype) */

    } /* (while something is in the tracestack) */

  FreeTracestack(dolist);

  *ret_trace = tr;
  return 1;
}




