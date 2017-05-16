/* fastmodelmaker.c
 *
 * Construct a covariance model from an alignment, using an approximate 
 * but fast algorithm. Complexity is N^2 in memory, N^3 in time. Get
 * back a covariance model, and information contents in bits for
 * both primary and secondary structure.
 *
 *   - use Stormo/Gutell MIXY algorithm to calculate pairwise covariances
 *     for all aligned column pairs i,j. This produces the matrix **mxy;
 *     it is indexed 0..alen-1 by 0..alen-1, as a flipped diagonal matrix
 *     mxy[j][i], j > i. The diagonal (i == j) is unused (and unallocated 
 *     for).
 *
 *   - use a quick and dirty Zuker-like algorithm to reduce this N^2 
 *     matrix to the optimal model of non-overlapping chords. This 
 *     produces the matrix **zmat. zmat is a flipped diagonal matrix, 
 *     zmat[j][i], j >= i. The diagonal (i==j) is initialized to zero. 
 *     All other cells take on positive values, summing covariances.
 *     
 *   - traceback thru the zmat matrix and construct a dynamic binary 
 *     tree *ztr, excluding inserted columns. *ztr 
 *     is constructed much the same as a normal traceback tree is constructed.
 *     emitl and emitr store the 0..alen-1 indices of the column this node is
 *     responsible for, even for BEGIN and BIFURC nodes. type is
 *     a node type (MATP_NODE), not a state type.
 *     Note that these are exceptions to the indexing scheme used
 *     by traceback structures everywhere else!
 *
 *   - Then, for each *individual* sequence, construct a fake traceback based
 *     on *ztr. Use these tracebacks to collect statistics from the alignment
 *     and initialize probabilities in a new model.
 */    

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "version.h"
#include "squid.h"
#include "structs.h"
#include "funcs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


static void mixy(char **aseq, int nseq, int alen, int ***ret_mixy);
static void free_mixy(int **mxy, int acol);
static void zfill(int **mxy, int acol, int ***ret_zmat, double *ret_sum);
static void free_zmat(int **zmat, int acol);
static void ztrace(int **zmat, double *gapfq, double gapthresh, int **mxy, 
		   int acol, struct trace_s **ret_ztr);

#ifdef DEBUG
static void dump_mixy(int **mxy, int alen);
static void dump_zmat(int **zmat, int alen);
static void dump_ztr(struct trace_s *ztr);
#endif /* DEBUG */

/* Function: Fastmodelmaker()
 * 
 * Purpose:  Given a sequence alignment, construct a reasonable
 *           starting model. Returns the model, as well as numbers
 *           for the primary and secondary structure information
 *           content of the alignment.
 *
 * Args:     aseqs       - sequence alignment. 
 *           ainfo       - info about the alignment
 *           nseq        - number of sequences
 *           prior       - prior distributions for CM construction
 *           gapthresh   - over this fraction of gaps, assign to INS
 *           ret_secinfo - RETURN: sec structure info content (bits) (maybe NULL)
 *           ret_cm      - RETURN: new model                         (maybe NULL)
 *           ret_mtr     - RETURN: master traceback for aseqs        (maybe NULL)
 *           
 * Return:   1 on success, 0 on failure.
 */
int
Fastmodelmaker(char **aseqs, AINFO *ainfo, int nseq, struct prior_s *prior,
	       double gapthresh, double *ret_secinfo, 
	       struct cm_s **ret_cm, struct trace_s **ret_mtr)
{
  struct cm_s *cm;
  int        **mxy;
  int        **zmat;
  int          M;
  double       secinfo;
  double      *gapfq;           /* frequencies of gap occurrence per column */
  struct trace_s  *ztr;
  struct trace_s  *tr;          /* individual traceback      */
  struct trmem_s  *pool;	/* memory pool for traceback */
  int          idx;
  int          gaps;
  int          apos;

  /* Precalculate gap frequencies seen in each column of alignment
   */
  if ((gapfq = (double *) malloc (sizeof(double) * ainfo->alen)) == NULL)
    Die("malloc failed, line %d file %s", __LINE__, __FILE__);
  for (apos = 0; apos < ainfo->alen; apos++)
    {
      gaps = 0;
      for (idx = 0; idx < nseq; idx++)
	if (isgap(aseqs[idx][apos])) gaps++;
      gapfq[apos] = (double) gaps / (double) nseq;
    }

				/* build Mxy matrix */
  mixy(aseqs, nseq, ainfo->alen, &mxy); 
				/* dynamic programming fill of zmat, using mxy */
  zfill(mxy, ainfo->alen, &zmat, &secinfo);
				/* traceback of zmat to make ztr tree */
  ztrace(zmat, gapfq, gapthresh, mxy, ainfo->alen, &ztr);
/*  PrintTrace(stdout, ztr); */

  NumberMasterTrace(ztr, &M);
				/* allocate for a model */
  if ((cm = AllocCM(M)) == NULL) Die("AllocCM() failed");
  TopofyNewCM(cm, ztr);

  /* For each sequence: convert consensus tree ztr to individual fake traceback.
   * Count traceback into new model. 
   */
  for (idx = 0; idx < nseq; idx++)
    {
      Transmogrify(ztr, aseqs[idx], &tr, &pool);

      if (! TraceCount(cm, aseqs[idx], 
		       (ainfo->sqinfo[idx].flags & SQINFO_WGT) ? 
		       ainfo->sqinfo[idx].weight : 1.0, 
		       tr))  
	Die("TraceCount() failed");

      FreeTrace(tr, pool);
    }

  if (! VerifyCM(cm))
    Die("Bad cm after trace counts.");

				/* convert CM to probabilities */
  ProbifyCM(cm, prior);
				/* garbage collect & return */
  free_mixy(mxy, ainfo->alen);
  free_zmat(zmat, ainfo->alen);
  free(gapfq);

  if (ret_mtr != NULL)     *ret_mtr     = ztr; else FreeTrace(ztr, NULL);
  if (ret_cm != NULL)      *ret_cm      = cm;  else FreeCM(cm);
  if (ret_secinfo != NULL) *ret_secinfo = secinfo / (double) nseq;
  return 1;
}




/* Function: mixy()
 * 
 * Purpose:  given a set of N aligned sequences aseq, calculate
 *           pairwise covariances (mutual information). ret_mixy
 *           is allocated, filled, and returned, as a diagonal 2D 
 *           (NxN) matrix of values. It must be freed by 
 *           the caller. It is a lower diagonal matrix mxy[j][i],
 *           j > i, 0..alen-1 by 0..j-1.
 *
 *           The values in mxy are integers. They are the average
 *           secondary structure information content (i.e. weighted for
 *           the number of pairs actually occurring in columns i,j)
 *           in bits, to two decimal places (i.e. info*100).
 *           
 * Returns:  mxy, which must be free'd by caller with free_mixy().
 */          
static void
mixy(char    **aseq,            /* array of aligned sequences, flushed right  */
     int       nseq,		/* number of aligned sequences */
     int       alen,		/* length of each sequence (all the same)  */
     int    ***ret_mxy)         /* RETURN: mxy array           */
{
  int    **mxy;                 /* RETURN: diagonal covariance matrix  */
  float   fx[ALPHASIZE];        /* singlet frequency vector            */
  float   fy[ALPHASIZE];	/* another singlet frequency vector    */
  float   fxy[ALPHASIZE][ALPHASIZE]; /* pairwise frequency 2D array    */
  int     idx;			/* counter for sequences               */
  int     i, j;			/* counters for columns x,y            */
  int     symi, symj;		/* counters for symbols                */
  int     pairs;		/* counter for pairs in which there are no gaps */
  long    test;

  /* Allocate for mxy
   */
  if ((mxy = (int **) malloc (alen * sizeof(int *))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  for (j = 1; j < alen; j++)
    if ((mxy[j] = (int *) malloc (j * sizeof(int))) == NULL)
      Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);

  test = (long) mxy[20];

  /* calculate mxy
   */
  for (j = 1; j < alen; j++)
    for (i = 0; i < j; i++)
      {
				/* zero counter array */
	for (symj = 0; symj < ALPHASIZE; symj++)
	  {
	    fx[symj] = fy[symj] = 0.0;
	    for (symi = 0; symi < ALPHASIZE; symi++)
	      fxy[symj][symi] = 0.0;
	  }
				/* count symbols in a column */
	pairs = 0;
	for (idx = 0; idx < nseq; idx++)
	  {
	    if (isgap(aseq[idx][i]) || isgap(aseq[idx][j]))
	      continue;

	    symi = SymbolIndex(aseq[idx][i]);
	    symj = SymbolIndex(aseq[idx][j]);
	    
	    fx[symi] += 1.0;
	    fy[symj] += 1.0;
	    fxy[symi][symj] += 1.0;
	    pairs++;
	  }

				/* convert to frequencies */
	if (pairs > 0)
	  for (symi = 0; symi < ALPHASIZE; symi++)
	    {
	      fx[symi] /= (float) pairs;
	      fy[symi] /= (float) pairs;
	      for (symj = 0; symj < ALPHASIZE; symj++)
		fxy[symi][symj] /= (float) pairs;
	    }

	/* calculate mxy. 144.269504 is a conversion of ln's into
         * bits * 100: i.e. 100 * (1/log(2)) 
	 */
	mxy[j][i] = 0;
	for (symi = 0; symi < ALPHASIZE; symi++)
	  for (symj = 0; symj < ALPHASIZE; symj++)
	    {
	      if (fxy[symi][symj] > 0.0)
		mxy[j][i] += (int) (144.269504 * fxy[symi][symj] *
				    log((fxy[symi][symj] / (fx[symi] * fy[symj]))));
	    }

	/* Sat Jul 17 22:17:17 1993:  We weight by pairs to get an expected score
	 * over all the sequences. Fixes a problem that columns with few symbols
         * could dominate the calculation just because of noise. 
	 */
	mxy[j][i] =  (mxy[j][i] * pairs) / nseq;
      }

  /* dump debugging info
   */

#ifdef DEBUG
  dump_mixy(mxy, alen);
#endif
  *ret_mxy = mxy;
}



/* Function: free_mixy()
 * 
 * Purpose:  free the space allocated for a flipped diagonal
 *           covariance matrix. To do this we also need to
 *           know alen, the number of columns in the starting
 *           sequence alignment.
 *
 * Returns:  (void)          
 */              
static void
free_mixy(int    **mxy,
	  int      alen)
{
  int j;

  for (j = 1; j < alen; j++)
    free(mxy[j]);
  free(mxy);
}




/* Function: zfill()
 * 
 * Purpose:  Calculate the optimal structure for a covariance matrix
 *           produced by mixy(). Uses a way-simplified form of the 
 *           Zuker/Nussinov dynamic programming RNA folding algorithm
 *           to find the structure which a) the emitted pairs sum
 *           to a maximum number of bits of covariance and b)
 *           has no overlapping chords (no pseudoknots). The dynamic
 *           programming matrix is allocated, filled, and returned.
 *           
 * Returns:  ret_zmat is returned thru a passed pointer; it must be 
 *           free'd by the caller using free_zmat().
 */
static void
zfill(int     **mxy,             /* diagonal covariance matrix from mixy()    */
      int       acol,		 /* size of mxy; number of aligned columns    */
      int    ***ret_zmat,        /* RETURN: filled dynamic programming matrix */
      double   *ret_sum)	 /* RETURN: total sum of Mxy for tree (bits)  */
{
  int    **zmat;
  int      i, j;
  int      diff;
  int      mid;

  /* Allocations.
   * zmat is a flipped diagonal matrix, inclusive of the diagonal;
   *  positions on both axes are numbered 0..acol-1.
   */  
  if ((zmat = (int **) malloc (acol * sizeof (int *))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  for (j = 0; j < acol; j++)
    if ((zmat[j] = (int *) malloc ((j+1) * sizeof(int))) == NULL)
      Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  
  /* Initialization.
   * We initialize the diagonal to 0.
   */
  for (j = 0; j < acol; j++)
    zmat[j][j] = 0;

  /* Dynamic programming stage.
   * Our recursion is:
   *    Sij = max { Si+1,j  (emit left, no covariance)
   *                Si,j-1  (emit right, no covariance)
   *                Si+1,j-1 + mxy[j][i].
   *                max over mid: Si,mid + Smid+1,j (bifurcation)
   *                }
   */
  for (diff = 1; diff < acol; diff++)
    for (i = 0; (j = i+diff) < acol; i++)
      {
	if (j-1 >= i+1)
	  zmat[j][i] = zmat[j-1][i+1] + mxy[j][i];
	else
	  zmat[j][i] = mxy[j][i];
	if (zmat[j][i+1] > zmat[j][i])
	  zmat[j][i] = zmat[j][i+1];
	if (zmat[j-1][i] > zmat[j][i])
	  zmat[j][i] = zmat[j-1][i];
	for (mid = i+1; mid < j-1; mid++)
	  if (zmat[mid][i] + zmat[j][mid+1] > zmat[j][i])
	    zmat[j][i] = zmat[mid][i] + zmat[j][mid+1];
      }

  /* Return
   */
#ifdef DEBUG
  dump_zmat(zmat, acol);
#endif

  *ret_sum  = (double) zmat[acol-1][0] / 100.0;
  *ret_zmat = zmat;
}


static void
free_zmat(int    **zmat,
	  int      acol)
{
  int j;
  
  for (j = 0; j < acol; j++)
    free(zmat[j]);
  free(zmat);
}



/* Function: ztrace()
 * 
 * Purpose:  Traceback through the dynamic programming matrix constructed
 *           by zfill(). Constructs a dynamic binary tree (ztr) of ztrace_s 
 *           structures, which keep track of both topology and the
 *           order in which various aligned columns are emitted.
 *           
 *           ztr ends up being the "shell" or template upon which the
 *           final model is built. ztr captures the branching structure
 *           of the model tree.
 *           
 *           Inserts are dealt with at this point. Columns with a gap frequency
 *           exceeding gapthresh are excluded from ztr. The final tree ztr contains
 *           MATR, MATL, MATP nodes only (with BEGIN and BIFURC of course).
 *         
 * Data:     ztr: a traceback tree.
 *                emitl = index of column emitted left (0..acol-1) 
 *                        or -1          
 *                emitr = index of column emitted right (0..acol-1)
 *                        or -1
 *                nodeidx = index of node in new model
 *                type   = type of node (MATP_NODE, etc.)
 *                        
 * Return:   ret_ztr is allocated here and must be free'd by the caller.
 */
static void
ztrace(int    **zmat,           /* dynamic programming matrix from zfill()  */
       double  *gapfq,          /* frequencies of gaps in columns 0..acol-1 */
       double   gapthresh,      /* above this, column is INS-generated      */
       int    **mxy,            /* the covariance matrix from mixy()        */
       int      acol,		/* number of aligned columns (size of zmat) */
       struct trace_s **ret_ztr)/* RETURN: binary tree of best structure    */
{
  struct trace_s *ztr;          /* binary tree of best structure       */
  struct trace_s *curr_ztr;     /* ptr to node of ztr we're working on */
  struct tracestack_s *dolist;  /* pushdown stack of active ztr nodes  */
  int i,j;			/* coords in zmat (0..acol-1)          */
  int mid;                      /* midpoint of a bifurcation           */

  /* Initialize.
   * Start at i = 0, j = acol-1 and work towards diagonal.
   */
  InitTrace(&ztr, NULL);        /* start a trace tree */
  dolist  = InitTracestack();	/* start a stack for traversing the trace tree */
				/* start with root aligned to 0..acol-1 */
  curr_ztr = AttachTrace(ztr, NULL, 0, acol-1, -1, ROOT_NODE);
  curr_ztr = AttachTrace(curr_ztr, NULL, 0, acol-1, -1, -1);
  PushTracestack(dolist, curr_ztr);
  
  while ((curr_ztr = PopTracestack(dolist)) != NULL)
    {
				/* where we are now in the traceback. */
      i = curr_ztr->emitl;
      j = curr_ztr->emitr;

				/* dummy END state on trace tree leaves. */
      if (i > j) curr_ztr->type = uEND_ST; /* never executes. */

				/* watch out for diagonal, where j-1,i+1 is nonsense */
      else if (i == j)		/* default to push-left; i is explained */
	{
	  curr_ztr->type  = MATL_NODE;
	  curr_ztr->nxtl->emitl = i+1;
	  curr_ztr->nxtl->emitr = j;
	}

      else if (j-1 == i && zmat[j][i] == mxy[j][i])
	{
	  curr_ztr->type = MATP_NODE;
	  curr_ztr->nxtl->emitl = i+1;
	  curr_ztr->nxtl->emitr = j-1;
	}

      else if (zmat[j][i] == zmat[j-1][i+1] + mxy[j][i])
	{		
	  curr_ztr->type = MATP_NODE;
	  if (j-1 >= i+1)
	    PushTracestack(dolist, AttachTrace(curr_ztr, NULL, i+1, j-1, -1,-1));
	  else
	    {
	      curr_ztr->nxtl->emitl = i+1;
	      curr_ztr->nxtl->emitr = j-1;
	    }
	}      
      else if (zmat[j][i] == zmat[j][i+1])
	{		
	  curr_ztr->type  = MATL_NODE;
	  if ( j >= i+1)
	    PushTracestack(dolist, AttachTrace(curr_ztr, NULL, i+1, j, -1,-1));
	  else
	    {
	      curr_ztr->nxtl->emitl = i+1;
	      curr_ztr->nxtl->emitr = j;
	    }
	}

      else if (zmat[j][i] == zmat[j-1][i])
	{
	  curr_ztr->type = MATR_NODE;
	  if ( j-1 >= i)
	    PushTracestack(dolist, AttachTrace(curr_ztr, NULL, i, j-1, -1,-1));
	  else
	    {
	      curr_ztr->nxtl->emitl = i;
	      curr_ztr->nxtl->emitr = j-1;
	    }
	}

      else
	{
	  for (mid = i+1; mid < j-1; mid++)
	    if (zmat[j][i] == zmat[mid][i] + zmat[j][mid+1])
	      {
		struct trace_s *branch;

				/* the current node is a bifurc node */
		curr_ztr->type  = BIFURC_NODE;
		
				/* it will connect to BEGIN-SEGMENT 
				 * nodes on either side, which are followed
				 * by normal segments again */
				/* right branch */
		branch = AttachTrace(curr_ztr, NULL, mid+1, j, -1, BEGINR_NODE);
		if (mid+1 <= j)
		  PushTracestack(dolist, AttachTrace(branch, NULL, mid+1, j, -1,-1));
		else
		  { branch->nxtl->emitl = mid+1; branch->nxtl->emitr = j; }
				/* left branch */
		branch = AttachTrace(curr_ztr, NULL, i, mid, -1, BEGINL_NODE);
		if (i <= mid)
		  PushTracestack(dolist, AttachTrace(branch, NULL, i, mid, -1,-1));
		else
		  { branch->nxtl->emitl = i; branch->nxtl->emitr = mid; }
		break;
	      }
	}

				/* clean up current node: deal with insertions */
      if (curr_ztr->type == MATP_NODE)
	{
	  if (gapfq[i] > gapthresh && gapfq[j] > gapthresh)
	    DeleteTracenode(curr_ztr, NULL);

	  else if (gapfq[i] > gapthresh)
	    curr_ztr->type = MATR_NODE;

	  else if (gapfq[j] > gapthresh)
	    curr_ztr->type = MATL_NODE;
	}
      else if ( (curr_ztr->type == MATR_NODE && gapfq[j] > gapthresh) ||
	        (curr_ztr->type == MATL_NODE && gapfq[i] > gapthresh))
	DeleteTracenode(curr_ztr, NULL);
    }

  FreeTracestack(dolist);

  *ret_ztr = ztr;
}
  



#ifdef DEBUG
static void
dump_mixy(int    **mxy,
	  int      alen)
{
  int i, j;

  for (j = 1; j < alen; j++)
    {
      for (i = 0; i < j; i++)
	printf("%6d", mxy[j][i]);
      puts("\n");
    }
}

static void
dump_zmat(int    **zmat,
	  int      alen)
{
  int i, j;

  for (j = 0; j < alen; j++)
    {
      for (i = 0; i <= j; i++)
	printf("%6d", zmat[j][i]);
      puts("\n");
    }
}

static void 
dump_ztr(struct trace_s *ztr)
{
  struct tracestack_s *dolist;
  struct trace_s      *curr;

  dolist = InitTracestack();
  PushTracestack(dolist, ztr->nxtl);

  while ((curr = PopTracestack(dolist)) != NULL)
    {
      printf("## ZTR STATE %#x\n", curr);
      printf("emitl    : %d\n", curr->emitl);
      printf("emitr    : %d\n", curr->emitr);
      printf("nodeidx  : %d\n", curr->nodeidx);
      printf("type     : %d\n", curr->type);
      if (curr->nxtl != NULL)
	printf("nxtl     : %#x\n", curr->nxtl);
      else
	printf("nxtl     : NULL\n");
      if (curr->nxtr != NULL)
	printf("nxtr     : %#x\n", curr->nxtr);
      else
	printf("nxtr     : NULL\n");
      
      if (curr->nxtr != NULL)
	PushTracestack(dolist, curr->nxtr);
      if (curr->nxtl != NULL)
	PushTracestack(dolist, curr->nxtl);
    }
  FreeTracestack(dolist);
}
#endif


