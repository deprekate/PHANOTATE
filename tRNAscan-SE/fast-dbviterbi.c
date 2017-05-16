/* fast-dbviterbi.c
 * SRE, Fri Sep 30 15:09:06 1994
 * 
 * Search with a covariance model
 * Fast "banded search" version of dbviterbi.c; subsequence
 * lengths are bounded by probabilistically determined bounds.
 *
 **********
 *
 * To optimize memory access patterns, the score storage is implemented
 * as a two-matrix version. amx is the
 * main storage. bmx is a smaller auxiliary matrix with a different
 * access pattern, holding scores of BEGIN state alignments; it
 * is used when calculating BIFURC scores.
 *
 * amx is [j = 0..1] [y = 0..statenum] [diff = 0..j] 
 *  diff == 0 is for off-diagonal boundary conditions (this is why diff is shifted +1)
 *  diff == 1 is for the diagonal, i==j
 *  We only need to keep two j rows in memory (current and previous).
 * Note that this is yet *another* memory access pattern and it's different
 * from dbviterbi.c!!!
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

static int  allocate_mx(struct istate_s *icm, int statenum, int window,	
			int ****ret_amx, int ****ret_bmx);
static int  init_mx    (struct istate_s *icm, int statenum, int N, 
			int ***amx, int ***bmx);
static int  recurse_mx (struct istate_s *icm, int statenum, int *minb, int *maxb,
			char *seq, int seqlen, 
			int window, int ***amx, int ***bmx, int ithresh, 
			int (*gotone_f)(int, int, double));
static void free_mx    (int ***amx, int ***bmx, int statenum, int window);


/* Function: FastViterbiScan()
 * 
 * Purpose:  Scanning version of the Viterbi alignment algorithm,
 *           for finding matches in a long sequence. 
 *
 * Args:     icm       - the model to align sequence to (int log-odds)
 *           statenum  - length of model in states (inclusive of END)
 *           minb      - minimum length bounds for states
 *           maxb      - maximum length bounds for states            
 *           seq       - sequence to align model to
 *           window    - scanning window size (nucleotides)
 *           thresh    - scores above this are reported through gotone_f()
 *           gotone_f  - function which gets told about a match
 *                       
 * Return:   1 on success, 0 on failure.                      
 */       
int
FastViterbiScan(struct istate_s *icm, int statenum, int *minb, int *maxb,
		char *seq, int window, double thresh,
		int (*gotone_f)(int, int, double))
{
  int   ***amx;			/* the main score matrix     */
  int   ***bmx;                 /* the BEGIN score matrix    */
  int      N;			/* length of sequence        */
  int      ithresh;             /* thresh, converted and scaled to int */

  N = strlen(seq);
  seq--;			/* convert to 1..N. Ugh! */
  ithresh = (int) (thresh * INTPRECISION);

  if (! allocate_mx(icm, statenum, window, &amx, &bmx)) return 0;
#ifdef DEBUG
  printf("allocated matrices\n");
#endif

  if (! init_mx(icm, statenum, window, amx, bmx)) return 0;
#ifdef DEBUG
  printf("matrices initialized\n");
#endif
  
  if (! recurse_mx(icm, statenum, minb, maxb, seq, N, window, amx, bmx, ithresh, gotone_f)) 
    return 0;
#ifdef DEBUG
  printf("recursion finished\n");
#endif
				/* terminate scanning hit reporting */
  ReportScanHit(-1,-1, 0.0, gotone_f);
  free_mx(amx, bmx, statenum, window);
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
 *           window   - length of scanning window
 *           ret_amx  - RETURN: main score matrix
 *           ret_bmx  - RETURN: BEGIN score matrix
 * 
 * Return:   Ptr to allocated scoring matrix, or
 *           dies and exits.
 */
static int 
allocate_mx(struct istate_s *icm,
	    int              statenum,
	    int              window,	
	    int          ****ret_amx,
	    int          ****ret_bmx)

{
  int    ***amx;              
  int    ***bmx;
  int     diag, j, y;		

  /* Main matrix, amx: fastest varying index is y (j,i,y)
   * we only keep two rows for j, 0 and 1.
   */
				/* malloc for j = 0..1 rows */
  if ((amx = (int   ***) malloc (2 * sizeof(int   **))) == NULL)
    Die("Memory allocation error in %s line %d", __FILE__, __LINE__);

  for (j = 0; j <= 1; j++)	/* loop over rows j = 0..1 */
    {
				/* malloc for diag = 0..window cols */
      if ((amx[j] = (int **)   malloc ((window + 1) * sizeof(int   *))) == NULL)
	Die("Memory allocation error in %s line %d", __FILE__, __LINE__);
      
				/* loop over cols diag = 0..window */
      for (diag = 0; diag <= window; diag++)
				/* malloc for y = 0..statenum-1 decks */
	  if ((amx[j][diag] = (int *)   malloc ((statenum) * sizeof (int  ))) == NULL)
	    Die("Memory allocation error in %s line %d", __FILE__, __LINE__);
    }

  
  /* B auxiliary matrix: fastest varying index is diag (y,j,diag)
   * bmx keeps score decks for BEGIN states
   */
				/* 0..statenum-1 decks */
  if ((bmx = (int ***)   malloc (statenum * sizeof(int   **))) == NULL)
    Die("Memory allocation error in %s line %d", __FILE__, __LINE__);

  for (y = 0; y < statenum; y++)
    {
      bmx[y] = NULL;
				/* we keep score info for BEGIN states */
      if (icm[y].statetype == uBEGIN_ST)
	{
				/* j= 0..window-1 rows  */
	  if ((bmx[y] = (int   **) malloc ((window) * sizeof(int   *))) == NULL)
	    Die("Memory allocation error in %s line %d", __FILE__, __LINE__);
				/* diff = 0..window columns */
	  for (j = 0; j < window; j++)
	    if ((bmx[y][j] = (int   *) malloc ((window+1) * sizeof(int  ))) == NULL)
	      Die("Memory allocation error in %s line %d", __FILE__, __LINE__);
	}
    }
  *ret_amx = amx;
  *ret_bmx = bmx;
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
	int    ***bmx,
	int     statenum,
	int     window)
{
  int diag, j, y;

  /* Free the main matrix, amx: 
   * amx[j][i][y] = [0..1] [0..window] [0..statenum-1]
   */
  for (j = 0; j <= 1; j++)
    {
      for (diag = 0; diag <= window; diag++)
	free(amx[j][diag]);
      free(amx[j]);
    }
  free(amx);

  /* Free the auxiliary matrix, bmx 
   * bmx[y][j][i] = [0..statenum-1] [0..window] [0..window]
   */
  for (y = 0; y < statenum; y++)
    {
      if (bmx[y] != NULL)
	{
	  for (j = 0; j < window; j++)
	    free(bmx[y][j]);
	  free(bmx[y]);
	}
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
	int              window,       /* size of scanning window on sequence */
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
    amx[0][0][y] = amx[1][0][y] = NEGINFINITY;
  for (diag = 1; diag <= window; diag++)
    {
      memcpy(amx[0][diag], amx[0][0], statenum * sizeof(int));
      memcpy(amx[1][diag], amx[0][0], statenum * sizeof(int));
    }

  /* Init the whole bmx to -Inf. We know state 0 is a begin (it's ROOT), so we
   * start there, and memcpy rows as needed.
   */
  for (diag = 0; diag <= window; diag++)  
    bmx[0][0][diag] = NEGINFINITY;
  for (j = 1; j < window; j++)
    memcpy(bmx[0][j], bmx[0][0], (window+1) * sizeof(int));

  for (y = 1; y < statenum; y++)
    if (bmx[y] != NULL)
      for (j = 0; j < window; j++)
	memcpy(bmx[y][j], bmx[0][0], (window+1) * sizeof(int));
  
  /* Init the off-diagonal (j = 0..window-1; diag == 0) with -log P scores.
   * End state = 0;
   * del, bifurc states are calc'ed
   * begin states same as del's
   * THIS IS WASTEFUL AND SHOULD BE CHANGED.
   */
  for (j = 0; j < window; j++)
    for (y = statenum-1; y >= 0; y--)
      {
	/* Set the alignment of END states to the off-diagonal (diag = 0)
	 * to be zero, and never touch them again.
	 */
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
				/* make a copy into bmx if y is a BEGIN */
	if (icm[y].statetype == uBEGIN_ST)
	  bmx[y][j][0] = amx[j%2][0][y];
      }

  return 1;
}



/* Function: recurse_mx()
 * 
 * Purpose:  Carry out the fill stage of the dynamic programming
 *           algorithm. After each j row is filled in, check the score
 *           of best full alignment ending at this row; if greater
 *           than threshold (ithresh), report it.
 *           
 * Returns:  1 on success, 0 on failure.
 */
static int
recurse_mx(struct istate_s *icm,      /* integer, state-form model */
	   int              statenum, /* number of states in icm   */
	   int             *minb,
	   int             *maxb,
	   char            *seq,      /* sequence, 1..seqlen       */
	   int              seqlen,   /* length of seq             */
	   int              window,   /* length of scanning window on seq */
	   int           ***amx,      /* main scoring matrix       */
	   int           ***bmx,      /* bifurc scoring matrix     */
	   int              ithresh,  /* reporting threshold       */
	   int            (*gotone_f)(int, int, double))
{
  int i, j, y;		        /* indices for 3 dimensions                    */
  int aj;			/* 0 or 1, index for j in A matrix             */
  int bj;			/* 0..window-1, index for j in B matrix        */
  int diff;			/* loop counter for difference: diff = j-i + 1 */
  int symi, symj;		/* symbol indices for seq[i], seq[j]           */
  int sc;			/* tmp for a score                             */
  int ynext;			/* index of next state y                       */
  int bestdiff, bestscore;

  int *beam;                    /* ptr to a beam (z-axis vector)               */
  int  leftdiff;		/* diff coord of BEGIN_L of a bifurc     */
  int  leftj;			/* j coord of BEGIN_L of a bifurc        */
  int **left_p;			/* pointer into whole 2D deck of BEGINL's of a bifurc */
  int *right_p;                 /* ptr into row of BEGIN_R's of a bifurc */
  int   *scp;			/* score pointer: ptr into beam of scores being calc'ed */
  struct istate_s *st;		/* state pointer: ptr at current state in icm */
  int *tmx;
  int  emitsc;

  for (j = 1; j <= seqlen; j++)
    {
      aj = j % 2;		/* 0 or 1 index in amx      */
      bj = j % window;          /* 0..window-1 index in bmx */
      symj = SymbolIndex(seq[j]);      

      for (y = statenum-1; y >= 0; y--)
	{
	  st  = &icm[y];
	  for (diff = minb[y]; diff <= maxb[y] && diff <= j; diff++)
	    {
	      i = j - diff + 1;
	      symi = SymbolIndex(seq[i]);
	      scp = &amx[aj][diff][y];

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
		    emitsc = 0;    
		    break;
		  case uMATP_ST: /* !aj toggles from 0 to 1 and vice versa */
		    if (diff == 1) continue;
		    beam   = amx[!aj][diff-2]; 
		    emitsc = st->emit[symi * ALPHASIZE + symj];
		    break; 
		  case uMATR_ST:
		  case uINSR_ST:
		    beam   = amx[!aj][diff-1];
		    emitsc = st->emit[symj];
		    break;
		  case uMATL_ST:
		  case uINSL_ST:  
		    beam   = amx[aj][diff-1];   
		    emitsc = st->emit[symi];
		    break;
		  case uEND_ST:   
		    continue;
		  default: Die("no such state type %d", st->statetype);
		  }
		  beam  += y + st->offset;
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
		    
		  /* Make a copy into bmx, btr if necessary
		   */
		  if (st->statetype == uBEGIN_ST)
		    bmx[y][bj][diff] = *scp;
		} /* end block of normal state stuff */
		
	      else		/* a BIFURC state */
		{
		  leftdiff = diff;
		  leftj    = bj;
		  right_p  = bmx[st->tmx[1]][leftj];
		  left_p   = bmx[st->tmx[0]];

				/* init w/ case that left branch emits it all */
		  *scp = left_p[leftj][leftdiff] + *right_p;
		  while (leftdiff > 0)
		    {
		      leftdiff--;
		      leftj = leftj ? leftj-1 : window-1; /* scan window wraparound */
		      right_p++;
			
		      sc = left_p[leftj][leftdiff] + *right_p;
		      if (sc > *scp)
			*scp = sc;
		    }
		}
		
	    } /* end loop over states */
	} /* end loop over diff */

      /* We've completed a row. Now we can examine the scores in diff, 
       * aj, ROOT_ST to decide whether to report this row. If we do, 
       * we report the 1..seqlen i, j coords of the matching subsequence
       * in seq, as well as the score converted to double-precision bits.
       */
      bestdiff  = 1;
      bestscore = bmx[0][bj][1];
      for (diff = 2; diff <= window; diff++)
	if (bmx[0][bj][diff] > bestscore)
	  {
	    bestscore = bmx[0][bj][diff];
	    bestdiff  = diff;
	  }
      if (bestscore > ithresh)
	if (! ReportScanHit(j - bestdiff + 1, j, (double)(bestscore / INTPRECISION), gotone_f))
	  Warn("caller ignored report of a match!");
      } /* end loop over j */
  
  return 1;
}

