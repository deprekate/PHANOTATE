/* lengthdist.c
 * SRE, Fri Sep 30 09:51:14 1994
 * 
 * Calculate length distributions expected at each state.
 */


#include <stdio.h>
#include <stdlib.h>

#include "squid.h"
#include "structs.h"
#include "funcs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif
#ifdef DEBUG
#include <assert.h>
#endif


/* Function: LengthDistribution()
 * 
 * Purpose:  Given a covariance model, calculate the
 *           length distribution that we expect each state
 *           to be aligned to. Uses a full forward 
 *           (summed probabilities) calculation.
 *           
 * Args:     cm      - covariance model (probability form)
 *           N       - maximum length to look at.                 
 *           ret_lmx - RETURN: (0..M-1) by (0..N) matrix of probabilities
 *                     for states 0..M-1 emitting lengths of (0..N).
 *                     
 * Return:   (void)
 *           ret_lmx is alloc'ed here. Free2DArray(*ret_lmx, M)
 */
void
LengthDistribution(struct pstate_s *pcm, int statenum, int N, double ***ret_lmx)
{
  double **lmx;
  int      y, len;
  int      ynext;
  int      mid;

  /* Allocate the matrix for storing probability
   * distributions.
   */
  lmx = (double **) MallocOrDie(statenum * sizeof(double *));
  for (y = 0; y < statenum; y++)
    lmx[y] = (double *) MallocOrDie ( (N+1) * sizeof(double));

  /* Set whole thing to zero.
   */
  for (y = 0; y < statenum; y++)
    for (len = 0; len <= N; len++)
      lmx[y][len] = 0.0;

  /* Initialize at length = 0.
   */
  for (y = statenum-1; y >= 0; y--)
    switch (pcm[y].statetype) {
    case uEND_ST:
      lmx[y][0] = 1.0;
      break;
      
    case uBIFURC_ST:
      lmx[y][0] = lmx[y+1][0] * lmx[pcm[y].bifr][0];
      break;
      
    case uDEL_ST:
    case uBEGIN_ST:
      for (ynext = 0; ynext < pcm[y].connectnum; ynext++)
	lmx[y][0] += pcm[y].tmx[ynext] * lmx[y + pcm[y].offset + ynext][0];
      break;
    }


  /* Recurse for lengths 1..N.
   */
  for (len = 1; len <= N; len++)
    for (y = statenum-1; y >= 0; y--)
      switch (pcm[y].statetype) {
      case uEND_ST: 
	break;

      case uBIFURC_ST:
	for (mid = 0; mid <= len; mid++)
	  lmx[y][len] += lmx[y+1][mid] * lmx[pcm[y].bifr][len-mid];
	break;
      
      case uDEL_ST:
      case uBEGIN_ST:
	for (ynext = 0; ynext < pcm[y].connectnum; ynext++)
	  lmx[y][len] += pcm[y].tmx[ynext] * lmx[y + pcm[y].offset + ynext][len];
	break;
      
      case uMATP_ST:
	if (len > 1)
	  for (ynext = 0; ynext < pcm[y].connectnum; ynext++)
	    lmx[y][len] += pcm[y].tmx[ynext] * lmx[y + pcm[y].offset + ynext][len-2];
	break;

      case uMATR_ST:
      case uMATL_ST:
      case uINSL_ST:
      case uINSR_ST:
	for (ynext = 0; ynext < pcm[y].connectnum; ynext++)
	  lmx[y][len] += pcm[y].tmx[ynext] * lmx[y + pcm[y].offset + ynext][len-1];
	break;

      default: Die("unrecognized state type %d", pcm[y].statetype);
      }

  *ret_lmx = lmx;
  return;
}



/* Function: LengthBounds()
 * 
 * Purpose:  Takes the probability distributions produced by 
 *           LengthDistribution() and produces a set of
 *           minimum and maximum length bounds, within which
 *           lies a specified amount of the probability.
 *           
 *           The algorithm is simple. Find the peak of the 
 *           probability distribution and include it. Then
 *           look left and right; choose whichever one is
 *           higher P, and include it. Continue until the
 *           included probability exceeds the target.
 *           
 * Args:     lmx:      probability distributions from LengthDistribution()
 *           statenum: # of rows in lmx
 *           N:        max length (lmx is [0..statenum-1[0..N])
 *           epsilon:  target probability is 1.0 - epsilon.
 *           ret_min:  RETURN: [0..statenum-1] array of minimum 
 *                     lengths for each state
 *           ret_max:  RETURN: [0..statenum-1] array of maximum 
 *                     lengths for each state        
 *                     
 * Return:   (void)
 *           ret_min, ret_max alloced here. free().
 */                    
void
LengthBounds(double **lmx, int statenum, int N, double epsilon, 
	     int **ret_min, int **ret_max)
{
  int *min, *max;
  int y, len;
  int cmin, cmax;		/* current min and max */
  double best;
  double p_remain;

  min = (int *) MallocOrDie (sizeof(int) * statenum);
  max = (int *) MallocOrDie (sizeof(int) * statenum);

  for (y = 0; y < statenum; y++)
    {
				/* danger! assuming that N was large enough! */
      DNorm(lmx[y], (N+1));

				/* step 1. find peak */
      best = lmx[y][0];
      cmin = 0;
      for (len = 1; len <= N; len++)
	if (lmx[y][len] > best) { best = lmx[y][len]; cmin = len; }

				/* that's where we start */
      cmax     = cmin;
      p_remain = 1.0 - lmx[y][cmin];

				/* extend to find bounds */
      while (p_remain > epsilon)
	{
	  if (cmin == 0 && cmax == N)
	    break;
/*
	    Die("state %d distribution not within set bound of %d; %f remains\n",
		 y, N, p_remain);
*/
	  else if (cmin == 0)	/* must look right */
	    p_remain -= lmx[y][++cmax];
	  else if (cmax == N)   /* must look left */
	    p_remain -= lmx[y][--cmin];
	  else			/* look left and right */
	    {
	      if (lmx[y][cmin-1] > lmx[y][cmax+1])
		p_remain -= lmx[y][--cmin];
	      else
		p_remain -= lmx[y][++cmax];
	    }
	}
      min[y] = cmin;
      max[y] = cmax;
    }
  *ret_min = min;
  *ret_max = max;
}
