/* mp-dbviterbi.m
 * MasPar MPL code for parallelized covariance model alignment algorithm
 * Database scanning version
 * Prototype started Fri Aug 19 08:44:43 1994
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ppeio.h>    /* parallel i/o */

#include "squid.h"
#include "structs.h"
#include "funcs.h"
#include "maspar.h"

#define PROHIBIT  -99999999	/* "negative infinity" */
#define MAXDEBUGWIDTH  10

/* mpcovels_main.c: 
 * functions that the DPU calls back to
 */
visible extern int DPUReportsHit(int *dpu_i, int *dpu_j, int *dpu_score);
visible extern int NextSequenceBlock(void);

/* Debugging functions.
 */
static void  mp_print_icm(FILE *fp, struct istate_s *icm, int nstates);
static char *mp_UstatetypeName(int ustatetype);
static char *mp_StatetypeName(int statetype);
static char *mp_NodetypeName(int nodetype);
static void  mp_print_mx(plural int *curr[VPENUM],
			struct istate_s *icm, int statenum,
			plural char symi[VPENUM], plural char symj);
static char NucSymbol(int symidx);

/* Function: MPViterbiScan()
 * 
 * Purpose:  MasPar MPL parallelized version of ViterbiScan().
 * 
 *           Scan a database of sequences with a covariance model.
 *           See Notes, dbviterbi-parallel, for more detailed
 *           commentary.
 *           
 * Args:     fe_icm:       ptr to the model (on front end)
 *           fe_statenum:  ptr to length of model (on front end)
 *           fe_threshold: report scores above this (scaled integer)
 *           fe_buffer:    where we'll fetch sequence from on FE
 *                     
 * Return:   nothing (I think this is proper for an MPL subroutine)
 */
visible int
MPViterbiScan(struct istate_s *fe_icm, int *fe_statenum,
	      int  *fe_threshold, char *fe_buffer)
{
  struct istate_s *icm;         /* covariance model, integer log odds form  */
  int         statenum;         /* length of icm                            */
  int         N;                /* maximum match size                       */
  int         vpe;		/* index of VPE: 0..VPENUM-1                */
  int         y;		/* index of state: 0..statenum-1            */
  int         j;		/* index of database position               */
  int         bpos;		/* index in bif arm: 0..N                   */
  int         sbuff_idx;        /* index in sequence buffers                */
  int         xdist;		/* how far we're xnetting                   */
  int         vpeconn;		/* index of VPE to connect to               */
  int         yl, yr;		/* connected state indices for a bifurc     */
  int         threshold;        /* hit reporting threshold                  */
  int         done;		/* flag for if there's any database left    */
  int         yidx, ynext;	/* indices for states                       */
  int         extend;		/* gap-extend cost for an insert state      */
  double      dtime;
  int 	      v;
  int	      sbifx;

  plural int *atmp; 	/* temporary storage while shoving bif over  */
  plural int *btmp;
  plural int  temp[VPENUM];   
  plural int  *prev[VPENUM];    /* beam of scores for previous row */
  plural int  *curr[VPENUM];    /* beam of scores for current row  */
  plural int  *tmp[VPENUM];     /* tmp pointers used for swapping prev, curr */
  plural int ***bif;            /* bifurcation "arms", remembering BEGINL diagonals */
  plural char  seqbuffer[BLOCKSIZE]; /* a block of sequence characters */
  plural char  symi[VPENUM];	/* index of symbol i, 0..3 for ACGU  */
  plural char  symj;		/* index of symbol j, 0..3 for ACGU  */
  plural int   score;           /* score                             */
  plural int   d;		/* d -- virtual column position, ixproc * VPENUM + vpe */
  plural int   bestconn;	/* best insert-left score to connect to on row */
  plural int   bifx;		/* plural circularly permuted bifurc arm position  */
  plural int   lasti;		/* i position of last hit */
  plural int   lastj;		/* j position of last hit */
  plural int   lastscore;	/* score of last hit      */
  plural int   myj;		/* what sequence position this row is working on */
  plural char  tmpi;		/* temporary storage while shoving symi over */
  plural int   seeme;		/* TRUE if we're reporting a hit in this row */

  dpuTimerStart();

  /* Transfer model and threshold from the front end to the ACU
   * fe_buffer address will be used for blockIn() calls later.
   */
  copyIn(fe_statenum,  &statenum,  sizeof(int));
  copyIn(fe_threshold, &threshold, sizeof(int));

  if ((icm     = (struct istate_s *) malloc (statenum * sizeof(struct istate_s))) == NULL)
    { fprintf(stderr, "ACU malloc failed"); exit(1); }
  copyIn(fe_icm,     icm,     statenum * sizeof(struct istate_s));


  /* Allocation of the scoring matrix
   * Each PE keeps a number of virtual PE's, #defined as VPENUM
   * each VPE keeps two scoring beams, prev and curr, and one
   * BEGINL diagonal, bif.
   */
  N = nxproc * VPENUM - 1;

  if ((atmp = (plural int *)  p_malloc ((N+1) * sizeof(int))) == NULL)
    { fprintf(stderr, "PE malloc failed, atmp\n"); exit(1); }

  if ((bif  = (plural int ***) malloc (VPENUM * sizeof(int **))) == NULL)
    { fprintf(stderr, "PE malloc failed, bifurc arm, level 1\n"); exit(1); }
  for (vpe = 0; vpe < VPENUM; vpe++)
    {
      if ((prev[vpe] = (plural int *)  p_malloc (statenum * sizeof(int))) == NULL ||
	  (curr[vpe] = (plural int *)  p_malloc (statenum * sizeof(int))) == NULL)
	{ fprintf(stderr, "PE malloc failed, score beam; vpe %d\n", vpe); exit(1); }
      if ((bif[vpe]  = (plural int **) malloc (statenum * sizeof(int *))) == NULL)
	{ fprintf(stderr, "PE malloc failed, bifurc arm, level 2; vpe %d\n", vpe); exit(1); }

      /* we don't have BEGINL states labeled directly, but we can
       * find BEGINL indices via their parent BIFURC state
       */
      for (y = statenum-1; y >= 0; y--)
	{
	  bif[vpe][y] = NULL;
	  if (icm[y].statetype == uBIFURC_ST)
	    {
	      yl = icm[y].tmx[0]; /* yl guaranteed > y */
	      if ((bif[vpe][yl] = (plural int *) p_malloc ((N+1) * sizeof(int))) == NULL)
		{ fprintf(stderr, "PE malloc failed, bifurc arm, level 3; vpe %d state %d\n", vpe, yl);
		  exit(1); }
	    }
	}
    }


  /* Initialize the scoring matrix.
   *   1) set the whole thing to -Infinity
   *   2) set the off-diagonal (d = 0): end state gets set to zero,
   *                                    begin, delete and bifurc's are calculated
   *      the d=0 column (VPE0 in ixproc 0) will not be modified ever again.
   */
				/* everything to negative infinity */
  for (y = 0; y < statenum; y++)
    for (vpe = 0; vpe < VPENUM; vpe++)
      curr[vpe][y] = prev[vpe][y] = PROHIBIT;
  for (vpe = 0; vpe < VPENUM; vpe++)
    for (y = statenum-1; y >= 0; y--)
      if (icm[y].statetype == uBIFURC_ST)
	{ 
	  yl = icm[y].tmx[0];
	  for (bpos = 0; bpos < N+1; bpos++)
	    bif[vpe][yl][bpos] = PROHIBIT;
	}
				/* d=0 column initialized in both prev and curr */
  for (y = statenum-1; y >= 0; y--)
    if (ixproc == 0)
      {
	switch (icm[y].statetype) {
	case uEND_ST:    
	  curr[0][y] = prev[0][y] = 0;
	  break;

	case uBIFURC_ST: 
	  curr[0][y] = prev[0][y] = curr[0][icm[y].tmx[0]] + curr[0][icm[y].tmx[1]]; 
	  break;

	case uDEL_ST:
	case uBEGIN_ST:  
	  for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; ynext++, yidx++)
	    if (curr[0][ynext] != PROHIBIT)
	      curr[0][y] = prev[0][y] = curr[0][ynext] + icm[y].tmx[yidx];

	  if (bif[0][y] != NULL) /* is this a BEGINL state? */
	    bif[0][y][0] = curr[0][y];
	  break;

	case uMATP_ST:
	case uMATL_ST:
	case uMATR_ST:
	case uINSL_ST:
	case uINSR_ST:
	  break;

	default:
	  fprintf(stderr, "init stage: unrecognized state type %d (%s), state %d\n", 
		  icm[y].statetype, mp_UstatetypeName(icm[y].statetype), y);
	  exit(1);
	}
      }

  /* Initialize coords for each of [nyproc] target sequences
   */
  myj = 0;

  /* Initialize scoring
   */
  if (ixproc == 0)
    lastj  = -1;

  /* Ready to sweep across the database.
   * j treats the database as (nyproc) streams of continuous sequence
   * rather than individual sequences.
   */
  for (j = 1, done = FALSE; !done; j++)
    {
      seeme = FALSE;		/* init: no rows reporting hits yet for this j */

      /* Check if we need to load new blocks of sequence from the
       * database. If we do, we read a new block into the first
       * PE of each row. If the read fails, we're done.
       */
      if ((sbuff_idx = (j-1) % BLOCKSIZE) == 0)
	{
	  if (callRequest(NextSequenceBlock, 0))
	    blockIn(fe_buffer, seqbuffer, 0, 0, 1, NYPROC, BLOCKSIZE);
	  else
	    break;		/* no more sequence. search terminates */
	}


      /* set up symbol indices symi, symj
       *   - shift symi one to the right
       *      (don't worry about wraparound, we'll overwrite the start)
       *   - xnet copy the new symj all the way thru 
       *   - set symi = symj in d = 1 (vpe 1, ixproc 0)
       */
      xnetE[1].tmpi = symi[VPENUM-1];      /* inter-PE shift to temp storage */
      for (vpe = VPENUM-1; vpe > 0; vpe--) /* intra-PE shifts */
	symi[vpe] = symi[vpe-1];
      symi[0] = tmpi;		           /* get vpe=0 from temp storage */
      if (ixproc == 0)		           /* take special care with ixproc=0 */
	{
	  xnetcE[N].symj = seqbuffer[sbuff_idx];
	  symi[1] = symj;
	}

      /* Swap curr, prev scoring rows
       * and shove the bif arrays eastwards by one.
       */
      for (vpe = 0; vpe < VPENUM; vpe++)
	{ tmp[vpe] = prev[vpe]; prev[vpe] = curr[vpe]; curr[vpe] = tmp[vpe]; }

      for (y = 0; y < statenum; y++)
	if (icm[y].statetype == uBIFURC_ST)
	  {
	    yl = icm[y].tmx[0];

            if (ixproc != nxproc-1)
	       for (sbifx = 0; sbifx <= N; sbifx++) 
		    xnetE[1].atmp[sbifx] = bif[VPENUM-1][yl][sbifx];

	    btmp = bif[VPENUM-1][yl];
	    for (vpe = VPENUM-1; vpe > 0; vpe--)      /* intra-PE shifts */
		  bif[vpe][yl] = bif[vpe-1][yl]; 
	    bif[0][yl] = btmp;

            if (ixproc != 0) {
		  btmp = bif[0][yl];
		  bif[0][yl] = atmp;	         /* get vpe=0 from temp storage */
	 	  atmp = btmp;
            }
	  }

      /* Reinitialize the current row -- blank it out to negative infinity.
       */
      if (ixproc > 0)	/* don't reinit vpe 0 in PE's ixproc == 0 */
	for (y = 0; y < statenum; y++)
	  curr[0][y] = PROHIBIT;
      for (y = 0; y < statenum; y++)
	for (vpe = 1; vpe < VPENUM; vpe++)
	  curr[vpe][y] = PROHIBIT;
      

      /* Check if we need to completely re-initialize a row (new sequence coming)
       * Reinitialization might go faster by recursive doubling with memcpy's.
       * This may not be terrific -- anyone who's asking for reinitialization forces
       * the other processors to idle.
       */
      if (symj == 4)		/* 4 is special flag, saying 'end of sequence' */
	{
				/* do we have a score to report? */
	  if (lastj != -1 && ixproc == 0)
	    {
	      seeme = TRUE;	/* sets seeme TRUE for all valid rows */
	      callRequest(DPUReportsHit, 16, &seeme, &lasti, &lastj, &lastscore);
	      lastj = -1;
	    }

	  if (ixproc > 0)
	    {			/* don't reinit vpe 0 in PE's ixproc == 0 */
	      for (y = 0; y < statenum; y++)
		curr[0][y] = prev[0][y] = PROHIBIT;
	      for (y = statenum-1; y >= 0; y--)
		if (icm[y].statetype == uBIFURC_ST)
		  {
		    yl = icm[y].tmx[0];
		    for (bpos = 0; bpos < N+1; bpos++)
		      bif[0][yl][bpos] = PROHIBIT;
		  }
	    }
	  
	  for (y = 0; y < statenum; y++)
	    for (vpe = 1; vpe < VPENUM; vpe++)
	      curr[vpe][y] = prev[vpe][y] = PROHIBIT;
	  for (vpe = 1; vpe < VPENUM; vpe++)
	    for (y = statenum-1; y >= 0; y--)
	      if (icm[y].statetype == uBIFURC_ST)
		{ 
		  yl = icm[y].tmx[0];
		  for (bpos = 0; bpos < N+1; bpos++)
		    bif[vpe][yl][bpos] = PROHIBIT;
		}
				/* reset individual sequence coord counters */
	  myj = 0;
	}


      /* Otherwise, go on to do the main recursion, stepping
       * through the states from end to root.
       */
      else
	{
	  for (y = statenum-1; y >= 0; y--)
	    switch (icm[y].statetype) {

	    case uBEGIN_ST:
	    case uDEL_ST:
	      /* OK to do VPE0 ixproc 0 init row */
	      for (vpe = 0; vpe < VPENUM; vpe++)
		{
		  if (ixproc * VPENUM + vpe <= myj+1) /* inactivate off-diag PE's */
		    {
		      for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; 
			   ynext++, yidx++)
			{
			  score = curr[vpe][ynext] + icm[y].tmx[yidx];
			  if (score > curr[vpe][y]) curr[vpe][y] = score;
			}
		      
		      if (bif[vpe][y] != NULL) /* is this a BEGINL state? save score if it is */
			bif[vpe][y][j % (N+1)] = curr[vpe][y];
		    }
		}
	      break;
	      
	    case uMATP_ST:
	      for (vpe = 0; vpe < VPENUM; vpe++)
		{
		  if (ixproc * VPENUM + vpe <= myj+1) /* inactivate off-diag PE's */
		    {
		      d = ixproc * VPENUM + vpe;
		      if      (VPENUM == 1 && vpe == 0) { xdist = 2; vpeconn = 0; }
		      else if (VPENUM > 1  && vpe < 2)  { xdist = 1; vpeconn = vpe + VPENUM - 2; }
		      else { xdist = 0; vpeconn = vpe-2; }
		      
		      if (d > 1)
			{
			  if (xdist > 0)
			    for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; 
				 ynext++, yidx++)
			      {
				score = xnetW[xdist].prev[vpeconn][ynext] + icm[y].tmx[yidx];
				if (score > curr[vpe][y]) curr[vpe][y] = score;
			      }
			  else	/* else intra PE communication */
			    for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; 
				 ynext++, yidx++)
			      {
				score = prev[vpeconn][ynext] + icm[y].tmx[yidx];
				if (score > curr[vpe][y]) curr[vpe][y] = score;
			      }
			  curr[vpe][y] += icm[y].emit[symi[vpe] * ALPHASIZE + symj];   
			}
		    }
		}
	      break;
	      
	      
	    case uMATR_ST:
	    case uINSR_ST:
	      /* VPE 0 case is interprocessor communication */
	      if (ixproc > 0 && ixproc * VPENUM <= myj+1) /* inactivate off-diag PE's */
		{
		  for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; ynext++, yidx++)
		    {
		      score = xnetW[1].prev[VPENUM-1][ynext] + icm[y].tmx[yidx];
		      if (score > curr[0][y]) curr[0][y] = score;
		    }
		  curr[0][y] += icm[y].emit[symj]; 
		}
	      
	      for (vpe = 1; vpe < VPENUM; vpe++) /* others are intraprocessor */
		{
		  if (ixproc * VPENUM + vpe <= myj+1) /* inactivate off-diag PE's */
		    {
		      for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; ynext++, yidx++)
			{
			  score = prev[vpe-1][ynext] + icm[y].tmx[yidx];
			  if (score > curr[vpe][y]) curr[vpe][y] = score;
			}
		      curr[vpe][y] += icm[y].emit[symj]; 
		    }
		}
	      break;
	      
	      
	    case uMATL_ST:
	      /* VPE 0 case is interprocessor communication */
	      if (ixproc > 0 && ixproc * VPENUM <= myj+1)
		{
		  for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; ynext++, yidx++)
		    {
		      score = xnetW[1].curr[VPENUM-1][ynext] + icm[y].tmx[yidx];
		      if (score > curr[0][y]) curr[0][y] = score;
		    }
		  curr[0][y] += icm[y].emit[symi[0]]; 
		}
	      
	      for (vpe = 1; vpe < VPENUM; vpe++) /* others are intraprocessor */
		{
		  if (ixproc * VPENUM + vpe <= myj+1) /* inactivate off-diag PE's */
		    {
		      for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; ynext++, yidx++)
			{
			  score = curr[vpe-1][ynext] + icm[y].tmx[yidx];
			  if (score > curr[vpe][y]) curr[vpe][y] = score;
			}
		      curr[vpe][y] += icm[y].emit[symi[vpe]]; 
		    }
		}
	      break;
	      
	      
	    case uINSL_ST:
	      /* INSL states are tricky because the recursion isn't 
	       * valid for y == ynext. We calculate for the connections
	       * that are valid first (yidx = 1..connectnum-1;
	       * then we do the y->y connections using a recursive doubling 
	       * trick.
	       */
	      /* VPE 0 is interprocessor */
	      if (ixproc > 0 && ixproc * VPENUM <= myj+1)
		{
		  for (ynext = y + 1, yidx = 1; yidx < icm[y].connectnum; ynext++, yidx++)
		    {
		      score = xnetW[1].curr[VPENUM-1][ynext] + icm[y].tmx[yidx];
		      if (score > curr[0][y]) curr[0][y] = score;
		    }
		  curr[0][y] += icm[y].emit[symi[0]]; 
		}
	      /* remaining VPEs intraprocessor */
	      for (vpe = 1; vpe < VPENUM; vpe++) /* others are intraprocessor */
		{
		  if (ixproc * VPENUM + vpe <= myj+1) /* inactivate off-diag PE's */
		    {
		      for (ynext = y + 1, yidx = 1; yidx < icm[y].connectnum; ynext++, yidx++)
			{
			  score = curr[vpe-1][ynext] + icm[y].tmx[yidx];
			  if (score > curr[vpe][y]) curr[vpe][y] = score;
			}
		      curr[vpe][y] += icm[y].emit[symi[vpe]]; 
		    }
		}
	      
	      /* Now the recursive doubling trick:
	       */
	      /* step 1. scan VPE's to get best score in a PE */
	      extend = icm[y].tmx[0];
	      for (vpe = 1; vpe < VPENUM; vpe++)
		if ((score = curr[vpe-1][y] + extend) > curr[vpe][y])
		  curr[vpe][y] = score;
	      
	      /* step 2. recursive double to find best connect */
	      bestconn = PROHIBIT;
	      for (xdist = 1; xdist < nxproc; xdist += xdist)
		if (ixproc > xdist) /* active PEs */
		  if ((score = xnetW[xdist].curr[VPENUM-1][y] + extend * VPENUM * (xdist-1)) > bestconn)
		    bestconn = score;
	      
	      /* step 3. scan bestconn back onto vpe's */
	      for (vpe = 0; vpe < VPENUM; vpe++)
		if (ixproc * VPENUM + vpe <= myj+1 &&
		    (score = bestconn + extend * (vpe+1)) > curr[vpe][y])
		  curr[vpe][y] = score;
	      
	      break;
	      
	    case uBIFURC_ST:
	      yl = icm[y].tmx[0]; /* deconvolute the hack to get connected indices */
	      yr = icm[y].tmx[1];

#ifdef SRE_REMOVED 	/* debug */
	      if (myj == 5)
		printf("bifurc connects to states %d and %d\n", yl, yr);
#endif
	      
	      for (vpe = 0; vpe < VPENUM; vpe++)
		{
		  for (v = 0; v < VPENUM; v++) temp[v]=curr[v][yr];
		  if (ixproc * VPENUM + vpe <= myj+1) /* inactivate off-diag PE's */
		    {
		      vpeconn = vpe;
		      xdist   = 0;
		      bifx    = (j - (ixproc * VPENUM + vpe)) % (N+1);
		      while (xdist <= nxproc)
			{
			  if (xdist == 0)
			    {
			      if ((score = bif[vpe][yl][bifx] + curr[vpeconn][yr]) > curr[vpe][y])
				curr[vpe][y] = score;
			    }
			  else if (xdist <= ixproc)
			    {
				all temp[vpeconn]=xnetW[1].temp[vpeconn];
				score = bif[vpe][yl][bifx] + temp[vpeconn];
				if (score > curr[vpe][y]) curr[vpe][y] = score;
			    }
			  vpeconn--;
			  if (vpeconn < 0) 
			    { vpeconn = VPENUM-1;
			      xdist++;
			    }
			  bifx++;
			  if (bifx > N) bifx = 0;
			}
		    }
		}
	      break;

	    case uEND_ST:
	      break;
	      
	    default:
	      fprintf(stderr, "Bogus statetype %d\n", icm[y].statetype);
	      exit(1);

	    } /* end switch over states */

	  
	  
	  
	  /* Now we've finished a row j.
	   * Collect best score, and report it if it's over the threshold
	   * Might be faster to check if there is such a score first?
	   */
	  /* step 1. scan vpe's to get best score in a PE */
	  bestconn = curr[0][0];
	  d        = ixproc * VPENUM;
	  for (vpe = 1; vpe < VPENUM; vpe++)
	    if (curr[vpe][0] > bestconn)
	      { 
		bestconn = curr[vpe][0];
		d = ixproc * VPENUM + vpe;
	      }
	  /* step 2. recursive doubling to get best score */
	  for (xdist = 1; xdist < nxproc; xdist += xdist)
	    {
	      if (ixproc + xdist <= nxproc)
		if (xnetW[xdist].bestconn < bestconn)
		  {
		    xnetW[xdist].bestconn = bestconn;
		    xnetW[xdist].d        = d;
		  }
	    }
	  
	  /* step 3. best scores now in ixproc 0.
	     check value, report if necessary */
	  if (ixproc == 0 && bestconn > threshold)
	    {
	      /* check for overlap with last hit: */
	      if ((myj - d + 1) > lastj) 
		{			/* no overlap. Report old, store new. */
		  if (lastj != -1)
		    {
		      seeme = TRUE;
		      callRequest(DPUReportsHit, 16, &seeme, &lasti, &lastj, &lastscore);
		    }
		  lasti      = myj - d + 1;
		  lastj      = myj;
		  lastscore  = bestconn;
		}
	      else if (bestconn > lastscore)
		{	/* overlap. compare 'em. keep new if it's better */
		  lasti = myj - d + 1;
		  lastj      = myj;
		  lastscore  = bestconn;
		}
	    }
	  

#ifdef SRE_REMOVED
	  mp_print_mx(curr, icm, statenum, symi, symj);  /* debugging */
	  printf("myj on row 0 is %d (%c)\n\n", proc[0][0].myj, NucSymbol(proc[0][0].symj)); 
#endif


	  myj++;		/* bump the individual sequence coord counters */

	} /* end "else we're not reinitializing, do main recursion stuff" */
    } /* end loop over database position j */


/*  mp_print_icm(stderr, icm, statenum); */  /* debugging */

  /* Cleanup and return
   */
  for (vpe = 0; vpe < VPENUM; vpe++)
    {
      p_free(prev[vpe]);
      p_free(curr[vpe]);
      for (y = statenum-1; y >= 0; y--)
	if (icm[y].statetype == uBIFURC_ST)
	  p_free(bif[vpe][y]);
      p_free(bif[vpe]);
    }
  free(icm);
  dtime=dpuTimerElapsed();
  printf("time for MPViterbiScan on DPU is %g \n",dtime);
  return 1;
}
  

      
/* Function: mp_print_icm()
 * 
 * Purpose:  Print an integer-version CM, as used by the alignment
 *           algorithms. MPL version of a function in debug.c.
 */
static void
mp_print_icm(FILE *fp, struct istate_s *icm, int nstates)
{
  int y;
  int x;

  for (y = 0; y < nstates; y++)
    {
      fprintf(fp, "node %d state %d (%s)\n", icm[y].nodeidx,
	      y, mp_UstatetypeName(icm[y].statetype));
      fprintf(fp, "     connectnum %d offset %d (connections start at %d)\n",
	      icm[y].connectnum, icm[y].offset, y + icm[y].offset);
      fprintf(fp, "     Transitions:  ");
      for (x = 0; x < icm[y].connectnum; x++)
	fprintf(fp, "%d ", icm[y].tmx[x]);
      fputs("\n", fp);
      
      fprintf(fp, "     Emissions:  ");
      switch (icm[y].statetype)
	{
	case uMATP_ST:
	  for (x = 0; x < ALPHASIZE * ALPHASIZE; x++)
	    fprintf(fp, "%d ", icm[y].emit[x]);
	  fputs("\n", fp);
	  break;

	case uMATR_ST:
	case uMATL_ST:
	case uINSR_ST:
	case uINSL_ST:
	  for (x = 0; x < ALPHASIZE; x++)
	    fprintf(fp, "%d ", icm[y].emit[x]);
	  fputs("\n", fp);
	  break;

	default:
	  fputs("NONE\n", fp);
	  break;
	}
    }
}

static char *
mp_UstatetypeName(int ustatetype)
{
  switch (ustatetype) {
  case uDEL_ST:    return "uDEL_ST"; 
  case uMATP_ST:   return "uMATP_ST";
  case uMATL_ST:   return "uMATL_ST";
  case uMATR_ST:   return "uMATR_ST";
  case uINSL_ST:   return "uINSL_ST";
  case uINSR_ST:   return "uINSR_ST";
  case uBEGIN_ST:  return "uBEGIN_ST";
  case uEND_ST:    return "uEND_ST";
  case uBIFURC_ST: return "uBIFURC_ST";
  default:         return "Unknown state type";
  }
}

static char *
mp_StatetypeName(int statetype)
{
  switch (statetype) {
  case DEL_ST:   return "DEL/BEG/BIF/END";
  case MATP_ST:  return "MATP_ST";
  case MATL_ST:  return "MATL_ST";
  case MATR_ST:  return "MATR_ST";
  case INSL_ST:  return "INSL_ST";
  case INSR_ST:  return "INSR_ST";
  default:       return "Unknown State";
  }
}

static char *
mp_NodetypeName(int nodetype)
{
  switch (nodetype) {
  case BIFURC_NODE:  return "BIF/END NODE";
  case MATP_NODE:    return "MATP_NODE";
  case MATL_NODE:    return "MATL_NODE";
  case MATR_NODE:    return "MATR_NODE";
  case BEGINL_NODE:  return "BEGINL_NODE";
  case BEGINR_NODE:  return "BEGINR_NODE";
  case ROOT_NODE:    return "ROOT_NODE";
  default:           return "Unknown Node";
  }
}


static char
NucSymbol(int symidx)
{
  switch (symidx) {
  case 0: return 'A'; 
  case 1: return 'C';
  case 2: return 'G';
  case 3: return 'U';
  case 4: return '*';
  default: return '?';
  }
}

static void
mp_print_mx(plural int *curr[VPENUM],
	    struct istate_s *icm, int statenum,
	    plural char symi[VPENUM], plural char symj)
{
  int d, vpe, xpr, y;

  /* Print out header
   */
  printf("%13s ", "");
  for (d = 0; d < MAXDEBUGWIDTH; d++)
    printf("    %2d     ", d);
  printf("\n%13s ", "");
  for (d = 0; d < MAXDEBUGWIDTH; d++)
    {
      xpr = d / VPENUM;
      vpe = d % VPENUM;
      printf("     %c     ", NucSymbol(proc[0][xpr].symi[vpe]));
    }
  printf("\n");

  /* Print out curr for all states
   */
  printf("\n");
  for (y = 0; y < statenum; y++)
    {
      printf("%2d %10s ", y, mp_UstatetypeName(icm[y].statetype));
      for (d = 0; d < MAXDEBUGWIDTH; d++)
	{
	  xpr = d / VPENUM;
	  vpe = d % VPENUM;
	  printf("%10d ", proc[0][xpr].curr[vpe][y]);
	}
      printf("\n");
    }

}
  

