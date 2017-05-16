/* mpviterbi.m
 * MasPar MPL code for parallelized covariance model alignment algorithm
 * Prototype started Thu Jul  7 12:04:04 1994, SRE
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpl.h>

#include "squid.h"
#include "structs.h"
#include "funcs.h"

/* Debugging functions.
 */
static void  mp_print_icm(FILE *fp, struct istate_s *icm, int nstates);
static void  mp_print_mx(FILE *fp, struct istate_s *icm, int statenum,
			 char *seq, int N, plural int *mx);
static char *mp_UstatetypeName(int ustatetype);
static char *mp_StatetypeName(int statetype);
static char *mp_NodetypeName(int nodetype);


/* Function: MPViterbiAlign()
 * 
 * Purpose:  MasPar MPL parallelized version of ViterbiAlign().
 * 
 *           Align a covariance model to a sequence. Return the alignment
 *           score and a traceback.
 *           
 * Args:     fe_icm:      ptr to the model (on front end)
 *           fe_statenum: ptr to length of model (on front end)
 *           fe_seq:      ptr to sequence (on front end)
 *           fe_N:        ptr to length of sequence (on front end)
 *           fe_score:    RETURN: ptr to score of alignment
 *           fe_trace:    RETURN: ptr to traceback
 *                     
 * Return:   nothing (I think this is proper for an MPL subroutine)
 */
visible int
MPViterbiAlign(struct istate_s *fe_icm, int *fe_statenum, char *fe_seq, 
	       int *fe_N, double *fe_score, struct trace_s *fe_trace)
{
  struct istate_s *icm;         /* covariance model, integer log odds form */
  int         statenum;         /* length of icm */
  char       *seq;              /* sequence to align to */
  int         N;                /* length of seq, 1..N */
  plural int  score;		/* temp variable for holding calc'ed scores */
  plural int *mx;               /* 3D scoring matrix, distributed across PE's */
  plural int *bscores;          /* tmp array used for bifurc calculations, 0..N */
  plural int  i,j;              /* each PE has its own base pair i,j */
  plural int  symi, symj;       /* each PE keeps an index for the base symbol at i,j */
  int         y, ynext, yidx;	/* indices for states */
  int         yl, yr;		/* state indices for bifurc children */
  int         hdist, vdist;     /* distances from i,j for bifurc calculation */
  int         extend;		/* length of insert extension */

  double      alignscore;       /* calculated alignment score */

  /* Transfer data from the front end to the ACU
   */
  copyIn(fe_statenum, &statenum, sizeof(int));
  copyIn(fe_N,        &N,        sizeof(int));
  if ((icm = (struct istate_s *) malloc (statenum * sizeof(struct istate_s))) == NULL ||
      (seq = (char *)            malloc ((N+2)    * sizeof(char)))            == NULL)
    { fprintf(stderr, "malloc failed"); exit(1); }
  copyIn(fe_icm, icm, statenum * sizeof(struct istate_s));
  copyIn(fe_seq, seq, (N+2) * sizeof(char));     /* (N+1) picks up the trailing '\0' */
  
  /* Allocation of the scoring matrix
   * Each PE keeps a "beam" of 0..statenum-1 scores. 
   */
  i  = ixproc;
  j  = iyproc;
  if (i <= N+1 && j <= N)
    {
      mx      = (plural int *) p_malloc (statenum * sizeof(int));
      bscores = (plural int *) p_malloc ((N+1)    * sizeof(int));
    }

  /* Set up the sequence in the PE's
   */
  if (i > 0 && j > 0 && i <= N+1 && j <= N)
    {
      switch (seq[i]) {
      case 'A': symi = 0; break;
      case 'C': symi = 1; break;
      case 'G': symi = 2; break;
      case 'T': symi = 3; break;
      case 'U': symi = 3; break;
      default:  symi = 0; break; /* sloppy garbage! */
      }
      switch (seq[j]) {
      case 'A': symj = 0; break;
      case 'C': symj = 1; break;
      case 'G': symj = 2; break;
      case 'T': symj = 3; break;
      case 'U': symj = 3; break;
      default:  symj = 0; break; /* sloppy garbage! */
      }
    }

  /* Initialize the scoring matrix.
   *   1) set the whole thing to -Infinity
   *   2) set the off-diagonal (i=j+1): end state gets set to zero,
   *                                    begin, delete and bifurc's are calculated
   */
  if (i <= N+1 && j <= N)
    for (y = 0; y < statenum; y++)
      mx[y] = -99999999;

  for (y = statenum-1; y >= 0; y--)
    if (i == j+1)
      {
	switch (icm[y].statetype) {
	case uEND_ST:    
	  mx[y] = 0;
	  break;

	case uBIFURC_ST: 
	  mx[y] = mx[icm[y].tmx[0]] + mx[icm[y].tmx[1]]; 
	  break;

	case uDEL_ST:
	case uBEGIN_ST:  
	  for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; ynext++, yidx++)
	    if (mx[ynext] != -99999999)
	      mx[y] = mx[ynext] + icm[y].tmx[yidx];
	  break;

	case uMATP_ST:
	case uMATL_ST:
	case uMATR_ST:
	case uINSL_ST:
	case uINSR_ST:
	  break;

	default:
	  fprintf(stderr, "init stage: unrecognized state type %d\n", icm[y].statetype);
	  exit(1);
	}
      }

  /* Recursion
   */
  if (j > 0 && j <= N && i <= j) 
    for (y = statenum-1; y >= 0; y--)
      {
	switch (icm[y].statetype) {
	case uBEGIN_ST:
	case uDEL_ST:
	  for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; ynext++, yidx++)
	    {
	      score = mx[ynext] + icm[y].tmx[yidx];
	      if (score > mx[y]) mx[y] = score;
	    }
	  break;
	  
	case uMATP_ST:
	  if (j > i)
	    {
	      for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; ynext++, yidx++)
		{
		  score = xnetNE[1].mx[ynext] + icm[y].tmx[yidx];
		  if (score > mx[y]) mx[y] = score;
		}
	      mx[y] += icm[y].emit[symi * ALPHASIZE + symj]; 
	    }
	  break;
	  
	case uMATR_ST:
	  for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; ynext++, yidx++)
	    {
	      score = xnetN[1].mx[ynext] + icm[y].tmx[yidx];
	      if (score > mx[y]) mx[y] = score;
	    }
	  mx[y] += icm[y].emit[symj]; 
	  break;
	  
	case uMATL_ST:
	  for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; ynext++, yidx++)
	    {
	      score = xnetE[1].mx[ynext] + icm[y].tmx[yidx];
	      if (score > mx[y]) mx[y] = score;
	    }
	  mx[y] += icm[y].emit[symi]; 
	  break;
	  
	case uINSR_ST:
	  for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; ynext++, yidx++)
	    {
	      score = xnetN[1].mx[ynext] + icm[y].tmx[yidx];
	      if (score > mx[y]) mx[y] = score;
	    }
	  /* that only gave us insert-opens. Now check all
	     lengths of insert-extend */
	  for (extend = 1; extend < N; extend++)

/*	    if (j - extend > 0) This is probably a bug. Thu Aug  4 09:52:15 1994 */
	    if (i + extend <= j) /* This is probably the right fix. */
	      {
		score = xnetN[extend].mx[y] + extend * icm[y].tmx[0];
		if (score > mx[y]) mx[y] = score;
	      }
	  mx[y] += icm[y].emit[symi]; 
	  
	  
	  break;
	  
	case uINSL_ST:
	  for (ynext = y + icm[y].offset, yidx = 0; yidx < icm[y].connectnum; ynext++, yidx++)
	    {
	      score = xnetE[1].mx[ynext] + icm[y].tmx[yidx];
	      if (score > mx[y]) mx[y] = score;
	    }
	  /* that only gave us insert-opens. Now check all
	     lengths of insert-extend */
	  for (extend = 1; extend < N; extend++)
	    if (i + extend <= j)
	      {
		score = xnetE[extend].mx[y] + extend * icm[y].tmx[0];
		if (score > mx[y]) mx[y] = score;
	      }
	  mx[y] += icm[y].emit[symi]; 
	  break;
	  
	case uBIFURC_ST:
	  yl = icm[y].tmx[0];
	  yr = icm[y].tmx[1];
				/* First: suck in the horizontal row */
	  for (hdist = 0; hdist <= N; hdist++)
	    if (i + hdist <= j+1) /* choose active PE's */
	      bscores[hdist] = xnetE[hdist].mx[yr];

				/* Second: suck in vertical row, backwards;
				   indices are more convoluted */
	  for (vdist = N; vdist >= 0; vdist--)
	    if (j-vdist >= i-1)	/* choose active PE's */
	      {
		score = bscores[j - vdist - i + 1] + xnetN[vdist].mx[yl];
		if (score > mx[y]) mx[y] = score;
	      }
	  break;
	  
	case uEND_ST:
	  break;

	default:
	  fprintf(stderr, "Bogus statetype %d\n", icm[y].statetype);
	  exit(1);
	}
      }

  /* Uncomment these next two lines to get debugging output.
   */
/*  mp_print_icm(stderr, icm, statenum);  */
/*  mp_print_mx(stderr, icm, statenum, seq, N, mx);  */


  /* Cleanup and return
   */
  alignscore = (double) proc[N][1].mx[0] / INTPRECISION;
  copyOut(&alignscore, fe_score, sizeof(double));
  p_free(mx);
  p_free(bscores);
  free(icm);
  free(seq);
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

/* Function: mp_print_mx()
 * 
 * Purpose:  Debugging: print out a (small) 3D scoring matrix.
 * 
 * Args:     fp  - usually stdout or stderr
 *           seq - 1..N sequence
 *           mx  - scoring matrix on the PE's     
 */
static void
mp_print_mx(FILE *fp, struct istate_s *icm, int statenum,
	    char *seq, int N, plural int *mx)
{

  int i, j, y;		/* indices for three dimensions */

  for (y = 0; y < statenum; y++)
    {
      fprintf(fp, "### Matrix for state %d, type %d (%s), from node %d\n",
	      y, icm[y].statetype, mp_UstatetypeName(icm[y].statetype), icm[y].nodeidx);
      fprintf(fp, "     ");
      for (i = 1; i <= N+1; i++)
	fprintf(fp, "%6d  ", i);
      fprintf(fp, "\n");

      for (j = 0; j <= N; j++)
	{
	  fprintf(fp, "%c %3d ", (j > 0) ? seq[j] : '*', j);
	  for (i = 1; i <= j+1; i++)
	    {
	      if (proc[j][i].mx[y] == -99999999)
		fprintf(fp, "%6s  ", "-INF");
	      else
		fprintf(fp, "%6d  ", proc[j][i].mx[y]);
	    }
	  fprintf(fp, "\n");
	}
      fprintf(fp, "\n\n");
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
