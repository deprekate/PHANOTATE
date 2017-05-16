/* smallviterbi.c
 * Mon Jan 24 11:38:21 1994
 * 
 * Small-memory version of viterbi.c
 * 
 * In the two-matrix version of the alignment algorithm, we keep
 * a matrix for the BEGIN states (B matrix) and a full matrix for all
 * the scores (A matrix). The database scanning version of the algorithm takes
 * advantage of the fact that, for scoring, we don't need to keep
 * a full cube around for matrix A; we only need the current row and the
 * last row. We only need a full cube of information for the BEGIN
 * state scores, which we get from the B matrix. 
 * 
 * This trick saves a large amount of memory, depending on the structure
 * of the model. (The fewer bifurcations and begins in the model, the more
 * memory this trick can save.) Unfortunately, it gives up the ability
 * to trace back and recover an alignment.
 * 
 * In this module, we add back just enough information to the scanning
 * algorithm to enable a traceback, at the expense of re-doing some 
 * calculation. The goal is to be able to fit 200-400 nt RNA sequences
 * into memory.
 * 
 * I will refer to a model "segment". A segment is a linear (unbranched)
 * chunk of the model, starting at a BEGIN/ROOT state, ending at a 
 * BIFURC/END state.
 * 
 * Matrix cells now carry traceback information. This number, tback, indicates
 * the i,j coords that this segment's BIFURC/END aligns to. 
 * tback is determined recursively. A BIFURC/END's tback points to itself.
 * All other tbacks are copied from the state that is connected to by
 * the maximally likely path.
 * 
 * A traceback is done by using bmx as a framework, and recalculating bits
 * of the alignment in between BEGINs and BIFURC/ENDs.  Thus, if
 * we know that a BEGIN aligns to i1, j1, we can get two numbers i2,j2 from
 * the BEGIN's tback, and know that the segment aligns to the sequence
 * (i1..i2-1)(j2+1..j1). Then we can recalculate the alignment of this
 * model segment to that subsequence, and reconstruct a full traceback.
 * 
 * tback is a single unsigned integer. The two numbers i,j are restricted
 * to 16 bits and are packed into tback by bit-shifting, i<<16.
 * 
 * 
 */

/*
 * amx and atr are [j = 0..N] [diff = 0..j] [y = 0..statenum]
 *  diff == 0 is for off-diagonal boundary conditions (this is why diff is shifted +1)
 *  diff == 1 is for the diagonal, i==j
 *
 * bmx and btr are [y = 0..statenum] [j = 0..N] [ diff = 0..j]
 *   a j,diff matrix exists only where y is a BEGIN state
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

/* This is how we pack tracebacks into a single machine word.
 * We assume a 32 bit int. 
 * If you're porting this code, all you have to do is make sure
 * pack_tb() puts two 16-bit ints in one data type TBACK, and unpack_tb()
 * gets them back.
 */
typedef unsigned int TBACK;
#define pack_tb(i,j)           ((i)<<16 | (j))
#define PACKED_I               0xFFFF0000U
#define PACKED_J               0x0000FFFFU
static void unpack_tb(TBACK tback, int *ret_i, int *ret_j)
{
  *ret_j =  tback & PACKED_J;
  *ret_i = (tback & PACKED_I) >> 16;
}



static int  allocate_mx(struct istate_s *icm,int statenum, int seqlen,	
			int ****ret_amx, TBACK ****ret_atr,
			int ****ret_bmx, TBACK ****ret_btr);
static int  init_mx    (struct istate_s *icm, int statenum, int N, 
			int ***amx, TBACK ***atr,
			int ***bmx, TBACK ***btr);
static int  recurse_mx (struct istate_s *icm, int statenum, char *seq, int N, 
			int ***amx, TBACK ***atr,
			int ***bmx, TBACK ***btr);
static int  trace_mx   (struct istate_s *icm, char *seq, int N, 
			int ***bmx, TBACK ***btr, struct trace_s **ret_tr);
static void free_mx    (int ***amx, TBACK ***atr, 
			int ***bmx, TBACK ***btr,
			int statenum, int seqlen);

#ifdef DEBUG
static void  print_tb(int tb);
static void  print_small_mx(FILE *fp, struct istate_s *icm, int statenum,
			    char *seq, int N, int ***bmx, TBACK ***btr);
#endif /* DEBUG */


/* Function: SmallViterbiAlign()
 * 
 * Purpose:  Align a sequence to a model, using the small-memory
 *           variant of the alignment algorithm. Return the score 
 *           of the alignment and the traceback.
 *
 * Args:     icm       - the model to align sequence to
 *           statenum  = number of states in icm
 *           seq       - sequence to align model to
 *           ret_score - RETURN: global alignment score
 *           ret_trace - RETURN: traceback tree
 *                       
 * Return:   1 on success, 0 on failure.                      
 */       
int
SmallViterbiAlign(struct istate_s *icm,       
		  int              statenum,
		  char            *seq,      
		  double          *ret_score,
		  struct trace_s **ret_trace)
{
  int   ***amx;			/* the main score matrix     */
  TBACK ***atr;                 /* amx's traceback pointers  */
  int   ***bmx;                 /* the BEGIN score matrix    */
  TBACK ***btr;                 /* bmx's traceback pointers  */
  int      N;			/* length of sequence        */

  N = strlen(seq);
  seq--;			/* convert to 1..N. Ugh! */

  if (! allocate_mx(icm, statenum, N, &amx, &atr, &bmx, &btr)) return 0;
#ifdef DEBUG
  printf("allocated matrices\n");
#endif

  if (! init_mx(icm, statenum, N, amx, atr, bmx, btr)) return 0;
#ifdef DEBUG
  printf("matrices initialized\n");
  print_small_mx(stdout, icm, statenum, seq, N, bmx, btr);
#endif
  
  if (! recurse_mx(icm, statenum, seq, N, amx, atr, bmx, btr)) return 0;
#ifdef DEBUG
  printf("recursion finished\n");
  print_small_mx(stdout, icm, statenum, seq, N, bmx, btr);
#endif

  *ret_score = ((double) bmx[0][N][N] / INTPRECISION);
#ifdef DEBUG
  printf("have a score of %.2f, starting traceback\n", *ret_score);
#endif

  if (! trace_mx(icm, seq, N, bmx, btr, ret_trace)) return 0;
#ifdef DEBUG
  printf("trace complete\n");
#endif

  free_mx(amx, atr, bmx, btr, statenum, N);
  return 1;
}



/* Function: allocate_mx()
 * 
 * Purpose:  Malloc space for the score matrices.
 *           amx and atr are indexed as j, i, y.
 *           bmx and btr are indexed as k, j, i.      
 *           In the two sequence dimensions j, i they are
 *           diagonal (+1 off diagonal) matrices with
 *           rows j = 0..N, i = 1..j+1.
 *           In the node dimension k bmx and btr are k = 0..M.
 *           In the state dimension y amx and atr are y = 0..numstates.
 *           
 * Args:     icm      - the int, log-odds, state-based model
 *           statenum - number of states in model
 *           seqlen   - length of sequence
 *           ret_amx  - RETURN: main score matrix
 *           ret_atr  - RETURN: amx's traceback pointers
 *           ret_bmx  - RETURN: BEGIN/BIFURC/END score matrix
 *           ret_btr  - RETURN: bmx's traceback pointers
 * 
 * Return:   Ptr to allocated scoring matrix, or
 *           dies and exits.
 */
static int 
allocate_mx(struct istate_s *icm,
	    int       statenum,
	    int       seqlen,	
	    int   ****ret_amx,
	    TBACK ****ret_atr,
	    int   ****ret_bmx,
	    TBACK ****ret_btr)	
{
  int    ***amx;              
  TBACK  ***atr;
  int    ***bmx;
  TBACK  ***btr;
  int     diag, j, y;		

  /* Main matrix, amx: fastest varying index is y (j,i,y)
   * we only keep two rows for j, 0 and 1.
   */
				/* malloc for j = 0..1 rows */
  if ((amx = (int   ***) malloc (2 * sizeof(int   **))) == NULL ||
      (atr = (TBACK ***) malloc (2 * sizeof(TBACK **))) == NULL)
    Die("Memory allocation error in %s line %d", __FILE__, __LINE__);


  for (j = 0; j <= 1; j++)	/* loop over rows j = 0..1 */
    {
				/* malloc for diag = 0..j (0..seqlen) cols */
      if ((amx[j] = (int **)   malloc ((seqlen + 1) * sizeof(int   *))) == NULL ||
	  (atr[j] = (TBACK **) malloc ((seqlen + 1) * sizeof(TBACK *))) == NULL)
	Die("Memory allocation error in %s line %d", __FILE__, __LINE__);
      
				/* loop over cols diag = 0..seqlen */
      for (diag = 0; diag <= seqlen; diag++)
				/* malloc for y = 0..statenum-1 decks */
	  if ((amx[j][diag] = (int *)   malloc ((statenum) * sizeof (int  ))) == NULL ||
	      (atr[j][diag] = (TBACK *) malloc ((statenum) * sizeof (TBACK))) == NULL)
	    Die("Memory allocation error in %s line %d", __FILE__, __LINE__);
    }

  
  /* B auxiliary matrices: fastest varying index is diag (y,j,diag)
   * bmx, btr keeps score, traceback decks for BEGIN states
   */
				/* 0..statenum-1 decks */
  if ((bmx = (int ***)   malloc (statenum * sizeof(int   **))) == NULL ||
      (btr = (TBACK ***) malloc (statenum * sizeof(TBACK **))) == NULL)
    Die("Memory allocation error in %s line %d", __FILE__, __LINE__);

  for (y = 0; y < statenum; y++)
    {
      bmx[y] = NULL;
      btr[y] = NULL;

				/* we keep score info for BEGIN and BIFURC states */
      if (icm[y].statetype == uBEGIN_ST || icm[y].statetype == uBIFURC_ST)
	{
				/* j= 0..seqlen rows  */
	  if ((bmx[y] = (int   **) malloc ((seqlen+1) * sizeof(int   *))) == NULL)
	    Die("Memory allocation error in %s line %d", __FILE__, __LINE__);
				/* i = 0..j columns */
	  for (j = 0; j <= seqlen; j++)
	    if ((bmx[y][j] = (int   *) malloc ((j+1) * sizeof(int  ))) == NULL)
	      Die("Memory allocation error in %s line %d", __FILE__, __LINE__);
	}
				/* We keep traceback info only for BEGIN states */
      if (icm[y].statetype == uBEGIN_ST)
	{
	  			/* j= 0..seqlen rows  */
	  if ((btr[y] = (TBACK **) malloc ((seqlen+1) * sizeof(TBACK *))) == NULL)
	    Die("Memory allocation error in %s line %d", __FILE__, __LINE__);
				/* i = 0..j columns */
	  for (j = 0; j <= seqlen; j++)
	    if ((btr[y][j] = (TBACK *) malloc ((j+1) * sizeof(TBACK))) == NULL)
	      Die("Memory allocation error in %s line %d", __FILE__, __LINE__);
	}

    }

  *ret_amx = amx;
  *ret_atr = atr;
  *ret_bmx = bmx;
  *ret_btr = btr;
  return 1;
}



/* Function: free_mx()
 * 
 * Purpose:  Free the space allocated to the scoring and traceback matrices.
 *           Precisely mirrors the allocations above in allocate_cvmx().
 * 
 * Return:   (void)
 */
static void
free_mx(int    ***amx,
	TBACK  ***atr,
	int    ***bmx,
	TBACK  ***btr,
	int     statenum,
	int     seqlen)
{
  int diag, j, y;

  /* Free the main matrix, amx: 
   * amx[j][i][y] = [0..1] [0..seqlen] [0..statenum-1]
   */
  for (j = 0; j <= 1; j++)
    {
      for (diag = 0; diag <= seqlen; diag++)
	{
	  free(amx[j][diag]);
	  free(atr[j][diag]);
	}
      free(amx[j]);
      free(atr[j]);
    }
  free(amx);
  free(atr);

  /* Free the auxiliary matrices, bmx and btr
   * bmx[y][j][i] = [0..statenum-1] [0..seqlen] [0..seqlen]
   */
  for (y = 0; y < statenum; y++)
    {
      if (bmx[y] != NULL)
	{
	  for (j = 0; j <= seqlen; j++)
	    free(bmx[y][j]);
	  free(bmx[y]);
	}
      if (btr[y] != NULL)
	{
	  for (j = 0; j <= seqlen; j++)
	    free(btr[y][j]);
	  free(btr[y]);
	}
    }
  free(bmx);
  free(btr);
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
	TBACK         ***atr,
	int           ***bmx,
	TBACK         ***btr)
{
  int diag, j, y;		/* counters for indices over the cvmx            */
  int ynext;			/* index of next state k+1                       */
  int *beam;                    /* z-axis vector of numbers in amx               */

  /* Init the whole amx to -Infinity. We do this with memcpy, trying
   * to be fast. We fill in j=0,diag=0 by hand, then memcpy() the other
   * columns.
   */
  for (y = 0; y < statenum; y++)
    amx[0][0][y] = amx[1][0][y] = NEGINFINITY;
  for (diag = 1; diag <= N; diag++)
    {
      memcpy(amx[0][diag], amx[0][0], statenum * sizeof(int));
      memcpy(amx[1][diag], amx[0][0], statenum * sizeof(int));
    }

  /* atr END and BIFURC traceback pointers point to themselves.
   * just set everything to point at itself.
   */
  for (j = 0; j <= 1; j++)
    for (diag = 0; diag <= N; diag++)
      for (y = 0; y < statenum; y++)
	atr[j][diag][y] = pack_tb(diag, j);

  /* Init the whole bmx to -Inf. We know state 0 is a begin (it's ROOT), so we
   * start there, and memcpy rows as needed.
   */
  for (diag = 0; diag <= N; diag++)  
    bmx[0][N][diag] = NEGINFINITY;
  for (j = 0; j < N; j++)
    memcpy(bmx[0][j], bmx[0][N], (j+1) * sizeof(int));

  for (y = 1; y < statenum; y++)
    if (bmx[y] != NULL)
      for (j = 0; j <= N; j++)
	memcpy(bmx[y][j], bmx[0][N], (j+1) * sizeof(int));
  
  /* Set all btr traceback ptrs to point at themselves
   */
  for (y = 0; y < statenum; y++)
    if (btr[y] != NULL)
      for (j = 0; j <= N; j++)
	for (diag = 0; diag <= j; diag++)
	  btr[y][j][diag] = pack_tb(diag,j);

  /* Init the off-diagonal (j = 0..N; diag == 0) with -log P scores.
   * End state = 0;
   * del, bifurc states are calc'ed
   * begin states same as del's
   * THIS IS WASTEFUL AND SHOULD BE CHANGED.
   */
  for (j = 0; j <= N; j++)
    for (y = statenum-1; y >= 0; y--)
      {
	if (icm[y].statetype == uEND_ST)
	  amx[j%2][0][y] = 0; 

      else if (icm[y].statetype == uBIFURC_ST)
	amx[j%2][0][y] = bmx[icm[y].tmx[0]][j][0] + bmx[icm[y].tmx[1]][j][0];

      else if (icm[y].statetype == uDEL_ST || icm[y].statetype == uBEGIN_ST)
	{
				/* only calc DEL-DEL and BEGIN-DEL transitions. Since
				 * we optimized the state transition tables, removing
				 * the unused ones, we don't know where the number
				 * for "to DEL" is! But we can find it, because it'll
				 * be the connection to a non-infinite score */
	  beam = amx[j%2][0] + y + icm[y].offset;
	  for (ynext = 0; ynext < icm[y].connectnum; ynext++)
	    {
	      if (*beam != NEGINFINITY)
		amx[j%2][0][y] = *beam + icm[y].tmx[ynext];
	      beam++;
	    }
	}
				/* make a copy into bmx if y is a BEGIN or BIFURC */
      if (icm[y].statetype == uBEGIN_ST ||
	  icm[y].statetype == uBIFURC_ST )
	bmx[y][j][0] = amx[j%2][0][y];
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
	   TBACK         ***atr,      /* tracebacks for amx        */
	   int           ***bmx,      /* bifurc scoring matrix     */
	   TBACK         ***btr)      /* tracebacks for btr        */
{
  int i, j, y;		        /* indices for 3 dimensions                    */
  int aj;			/* 0 or 1, index for j in A matrices           */
  int diff;			/* loop counter for difference: diff = j-i + 1 */
  int symi, symj;		/* symbol indices for seq[i], seq[j]           */
  int sc;			/* tmp for a score                             */
  int ynext;			/* index of next state y                       */

  int *beam;                    /* ptr to a beam (z-axis vector)               */
  TBACK *beamt;                 /* ptr into connected traceback beam           */
  int  leftdiff;		/* diff coord of BEGIN_L of a bifurc     */
  int  leftj;			/* j coord of BEGIN_L of a bifurc        */
  int **left_p;			/* pointer into whole 2D deck of BEGINL's of a bifurc */
  int *right_p;                 /* ptr into row of BEGIN_R's of a bifurc */
  int   *scp;			/* score pointer: ptr into beam of scores being calc'ed */
  TBACK *sct;                   /* tback beam being calc'ed */
  struct istate_s *st;		/* state pointer: ptr at current state in icm */
  int *tmx;
  int  emitsc;

  for (j = 1; j <= N; j++)
    {
      aj = j % 2;
      symj = SymbolIndex(seq[j]);      

				/* we have to init END and BIF states to point
				   at themselves in this row of atr */
      for (diff = 0; diff <= j; diff++)
	for (y = 0; y < statenum; y++)
	  if (icm[y].statetype == uBIFURC_ST ||
	      icm[y].statetype == uEND_ST)
	    atr[aj][diff][y] = pack_tb(diff, j);

      for (diff = 1; diff <= j; diff++)
	{
	  i = j - diff + 1;
	  symi = SymbolIndex(seq[i]);


	  scp = &amx[aj][diff][statenum-1];
	  sct = &atr[aj][diff][statenum-1];
	  st  = &icm[statenum-1];
	  for (y = statenum-1; y >= 0; y--, scp--, sct--, st--)
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
		    beam   = amx[aj][diff];    
		    beamt  = atr[aj][diff];
		    emitsc = 0;    
		    break;
		  case uMATP_ST: /* !aj toggles from 0 to 1 and vice versa */
		    if (diff == 1) continue;
		    beam   = amx[!aj][diff-2]; 
		    beamt  = atr[!aj][diff-2]; 
		    emitsc = st->emit[symi * ALPHASIZE + symj];
		    break; 
		  case uMATR_ST:
		  case uINSR_ST:
		    beam   = amx[!aj][diff-1];
		    beamt  = atr[!aj][diff-1];
		    emitsc = st->emit[symj];
		    break;
		  case uMATL_ST:
		  case uINSL_ST:  
		    beam   = amx[aj][diff-1];   
		    beamt  = atr[aj][diff-1];
		    emitsc = st->emit[symi];
		    break;
		  case uEND_ST:   
		    continue;
		  default: Die("no such state type %d", st->statetype);
		  }
		  beam  += y + st->offset;
		  beamt += y + st->offset;
		  tmx  = st->tmx;


		  /* Init for ynext == 0 case 
		   */		
		  *scp = *beam + *tmx;
		  *sct = *beamt;
		    
		  /* Calculate remaining cases
		   */
		  for (ynext = 1; ynext < st->connectnum; ynext++)
		    {
		      beam++;
		      beamt++;
		      tmx++;
		      if (*beam > *scp)
			{
			  sc = *beam + *tmx;
			  if (sc > *scp)
			    {
			      *scp = sc;
			      *sct = *beamt;
			    }
			}
		    }
		    
		  /* Add emission scores now
		   */
		  *scp += emitsc;
		    
		  /* Make a copy into bmx, btr if necessary
		   */
		  if (st->statetype == uBEGIN_ST)
		    {
		      bmx[y][j][diff] = *scp;
		      btr[y][j][diff] = *sct;
		    }
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
				/* keep copy of score in bmx, for tracing */
		    bmx[y][j][diff] = *scp;
		  }
		
	      } /* end loop over states */
	  } /* end loop over diff */
      } /* end loop over j */
  return 1;
}


/* Function: trace_mx()
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
	 int           ***bmx,       /* matrix of BEGIN scores           */
	 TBACK         ***btr,       /* matrix of BIFURC/END tbacks      */
	 struct trace_s **ret_trace) /* RETURN: the traceback tree       */
{
  struct trace_s *tr;           /* the traceback tree under construction */
  struct trace_s *curr_tr;      /* ptr to node of tr we're working on    */
  struct tracestack_s *dolist;  /* pushdown stack of active tr nodes     */
  int diff,i, j;		/* coords in mx (0..N)                   */
  int y;			/* counter for states (0..statenum-1)    */
  int leftdiff;
  int leftj;
  int *right_p;
  int  i2, j2;			/* what's left unaccounted for at segment end */
  int  diff2;
  int  end_y;			/* index of state that ends segment           */

  /* Initialize.
   * Start at i = 1, j = N and work towards diagonal
   */
  InitTrace(&tr, NULL);         /* start a trace tree */
  dolist = InitTracestack();	/* start a stack for traversing the trace tree */

  curr_tr = AttachTrace(tr, NULL, 0, N-1, 0, BEGIN_ST);
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
      
				/* find out which bifurc/end state terminates this
				   segment */
      end_y = y+1;
      while (icm[end_y].statetype != uBIFURC_ST &&
	     icm[end_y].statetype != uEND_ST)
	end_y++;

				/* find out i2,j2 that the terminal bifurc/end
				   aligns to */
      unpack_tb(btr[y][j][diff], &diff2, &j2);
      i2 = j2 - diff2 + 1;

      /* For now, just write out what the traceback looks like.
       */
      printf("Segment from state %d to %d: aligns to %d..%d/%d..%d\n",
	     y, end_y, i, i2-1, j2+1, j);


				/* push next BEGINs onto stack; they are connected
				   to BIFURC end_y */
      if (icm[end_y].statetype == uBIFURC_ST)
	{
	  if (i2 > j2)
	    {
	      PushTracestack(dolist, AttachTrace(curr_tr, NULL, i2-1, j2-1, icm[end_y].tmx[1], BEGIN_ST));
	      PushTracestack(dolist, AttachTrace(curr_tr, NULL, i2-1, j2-1, icm[end_y].tmx[0], BEGIN_ST));
	    }
	  else
	    {
	      leftdiff = diff2;
	      leftj    = j2;
	      right_p  = bmx[icm[end_y].tmx[1]][j2];
	      
	      while (leftdiff >= 0)
		{
		  if (bmx[end_y][j2][diff] == bmx[icm[end_y].tmx[0]][leftj][leftdiff] + *right_p)
		    {
		      printf("found the bifurc: it is %d-%d and %d-%d\n",
			     i2 -1, i2+leftdiff-2, i2 + leftdiff-1, j2-1);
		      PushTracestack(dolist, AttachTrace(curr_tr, NULL, i2 + leftdiff-1, j2-1, 
							 icm[end_y].tmx[1], BEGIN_ST));
		      PushTracestack(dolist, AttachTrace(curr_tr, NULL, i2 -1, i2+leftdiff-2, 
							 icm[end_y].tmx[0], BEGIN_ST));
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
    } /* (while something is in the tracestack) */

  FreeTracestack(dolist);

  *ret_trace = tr;
  return 1;
}



#ifdef DEBUG

/* Function: print_tb()
 * 
 * Purpose:  Print the two numbers of a packed traceback.
 */
static void
print_tb(int tb)
{
  int i, j;

  unpack_tb(tb, &i, &j);
  printf("%d %d\n", i, j);
}

/* Function: PrintSmallMX()
 * 
 * Purpose:  Debugging output; print out the three-dimensional
 *           auxiliary alignment matrix produced by the 
 *           small-memory version.
 */
static void
print_small_mx(FILE            *fp,    /* open file or just stdout/stderr */
	       struct istate_s *icm,   /* the model to align           */   
	       int              statenum, /*  number of states in icm  */
	       char            *seq,   /* sequence, 1..N               */   
	       int              N,     /* length of seq                */
	       int           ***bmx,   /* auxiliary matrix             */
	       TBACK         ***btr)   /* traceback ptrs for bmx       */
{
  int j, diff, y;		/* indices in 3D matrix */
  int tbdiff, tbj;              /* traceback pointers to a j, diff position */

  for (y = 0; y < statenum; y++)
    if (bmx[y] != NULL)
      {
	fprintf(fp, "### B Matrix for state %d, type %d (%s), from node %d\n",
		y, icm[y].statetype, UstatetypeName(icm[y].statetype), icm[y].nodeidx);
	fprintf(fp, "     ");
	for (diff = 0; diff <= N; diff++)
	  fprintf(fp, "%6d  ", diff);
	fprintf(fp, "\n");
	
	for (j = 0; j <= N; j++)
	  {
	    fprintf(fp, "%c %3d ", (j > 0) ? seq[j] : (char) '*', j);
	    for (diff = 0; diff <= j; diff++)
	      fprintf(fp, "%6d  ", bmx[y][j][diff]);
	    fprintf(fp, "\n      ");
	    if (icm[y].statetype == uBEGIN_ST)
	      for (diff = 0; diff <= j; diff++)
		{
		  unpack_tb(btr[y][j][diff], &tbdiff, &tbj);
		  fprintf(fp, "%3d/%3d ", tbdiff, tbj);
		}
	    fprintf(fp, "\n");
	  }
	fprintf(fp, "\n\n");
      }
}
#endif /* DEBUG */


