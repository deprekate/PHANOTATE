/* maxmodelmaker.c - maximum likelihood construction of a covariance model
 * Tue Aug 31 14:54:35 1993
 * 
 * Given a multiple sequence alignment, construct the model
 * which generates that alignment with maximal likelihood.
 * Uses a dynamic programming algorithm to assign a score
 * to the optimal subtree with a "root" node MATP, MATR, MATL,
 * BEGINL, BEGINR, or BIFURC aligned to matrix cell i,j,
 * where i and j are column and row positions in a multiple
 * sequence alignment. The alignment of ROOT to i=0,j=N-1 is
 * the score of the best model; a traceback from this point
 * creates the optimal tree structure. The tree structure is
 * explicitly aligned to the multiple sequence alignment,
 * so it is trivial to estimate the parameters of the model.
 * 
 * 
 * 
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "version.h"
#include "structs.h"
#include "funcs.h"


#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#define MAXINSERT 6		/* maximum INSL or INSR path length between MATP nodes */


/* Structure maxmx_s
 * 
 * One per cell of the 2D diagonal matrix of the alignment
 * against itself.
 */
struct maxmx_s {
  int   sc[NODETYPES-1];	/* scores of assigning each possible node type */

				/*** Traceback info:  */
  short matp_i2;		/* matp assignment connects to node ftype at i2,j2 */
  short matp_j2;
  char  matp_ftype;
  short matl_i2;		/* matl assignment connects to node ftype at i2,j  */
  char  matl_ftype;
  short matr_j2;		/* matr assignment connects to node ftype at i,j2  */
  char  matr_ftype;
  char  begl_ftype;		/* begl assignment connects to node ftype at i,j   */
  short begr_i2;		/* begr assignment connects to node ftype at i2, j */
  char  begr_ftype;
  short bifurc_mid;		/* best bifurcation is into i,mid, mid+1,j */
};
  
static struct maxmx_s **alloc_maxmx(int alen);
static void init_maxmx(struct maxmx_s **mmx, int nseq, int alen,
		       struct prior_s *prior, int *mscore, double *gapcount);
static void recurse_maxmx(int **aseqsT, float *weights, int alen, int nseq, 
			  struct prior_s *prior, int *mscore, double *gapcount,
			  double gapthresh, struct maxmx_s **mmx);
static void trace_maxmx(struct maxmx_s **mmx, int alen, struct trace_s **ret_mtr);
static void transpose_alignment(char **aseqs, int alen, int nseq, int ***ret_aseqsT);
static void singlet_emissions(int **aseqsT, float *weights, int alen,int nseq,
			      struct prior_s *prior, int **ret_mscore, double **ret_gapcount);
static int  pair_emissioncost(int *coli, int *colj,float *weights, int nseq, struct prior_s *prior);
static void frommatp_transtable(int **aseqsT, float *weights, int nseq,
				int i, int j, int i2, int j2, int *accum_insl,
				int *accum_insr, 
				double  trans[STATETYPES][STATETYPES]);
static void frommatl_transtable(int **aseqsT, float *weights, int  nseq,
				int i, int  j, int i2, int    *accum_insl,
				double  trans[STATETYPES][STATETYPES]);
static void frommatr_transtable(int **aseqsT, float *weights, int nseq, int i, 
				int j, int j2, int *accum_insr, 
				double  trans[STATETYPES][STATETYPES]);
static void frombeginr_transtable(int **aseqsT, float *weights, int nseq,
				  int j, int i2, int *accum_insl,
				  double  trans[STATETYPES][STATETYPES]);
static void frombeginl_transtable(int **aseqsT, float *weights, int nseq,
				  int i, int j, double  trans[STATETYPES][STATETYPES]);
static void fromroot_transtable(int **aseqsT, float *weights, int nseq, int i2,int j2,
				int *accum_insl, int *accum_insr,
				double  trans[STATETYPES][STATETYPES]);
static int assign_cell(int i, int j, int symi, int symj);
static void to_matp_transtable(double master_table[STATETYPES][STATETYPES],
			       double trans[STATETYPES][STATETYPES]);
static void to_matr_transtable(double master_table[STATETYPES][STATETYPES],
			       double trans[STATETYPES][STATETYPES]);
static void to_matl_transtable(double master_table[STATETYPES][STATETYPES],
			       double trans[STATETYPES][STATETYPES]);
static void to_bifurc_transtable(double master_table[STATETYPES][STATETYPES],
				 double trans[STATETYPES][STATETYPES]);
static int dot_score(double *cvec, double *pvec, int veclen);

#ifdef DEBUG
static void print_mmx(struct maxmx_s **mmx, int alen);
static void print_assignments(struct trace_s *mtr, int nseq, int alen, double *gapcount);
#endif

/* MACROS: copy_transtable()
 *         copy_singlet()
 *         copy_pairwise()
 *         
 * For speed, we use memcpy() to do the operations, relying
 * on the fact that C stores 2D arrays in 1D.
 */
#define copy_transtable(tomx, frommx) memcpy((void *) tomx, (void *) frommx, sizeof(double) * STATETYPES * STATETYPES)
#define copy_singlet(tovec, fromvec)  memcpy((void *) tovec, (void *) fromvec, sizeof(double) * ALPHASIZE)
#define copy_pairwise(tomx, frommx)   memcpy((void *) tomx, (void *) frommx, sizeof(double) * ALPHASIZE * ALPHASIZE)

/* MACROS: zero_transtable()
 *         zero_singlet()
 *         zero_pairwise()
 *         
 * For speed, we use memset() to initialize, relying
 * on the fact that C stores 2D arrays in 1D. Are doubles
 * really eight 0 bytes when they're == 0.0?
 */
#define zero_transtable(mx) memset(mx,  (char) 0, sizeof(double) * STATETYPES * STATETYPES)
#define zero_singlet(vec)   memset(vec, (char) 0, sizeof(double) * ALPHASIZE)
#define zero_pairwise(mx)   memset(mx,  (char) 0, sizeof(double) * ALPHASIZE * ALPHASIZE)

/* Function: Maxmodelmaker()
 * 
 * Purpose:  Create a maximally likely model structure from a multiple
 *           sequence alignment.
 *           
 * Args:     aseqs       - flushed sequence alignment; each seq is 0..alen-1
 *           ainfo       - information about the alignment
 *           nseq        - number of sequences
 *           gapthresh   - heuristic: fractional occupancy <= a column must be a MAT of some kind
 *           prior       - prior probability distributions, and alphabet info
 *           ret_ssinfo  - RETURN: total info content of the alignment (may pass NULL)
 *           ret_cm      - RETURN: the new covariance model            (may pass NULL)
 *           ret_mtr     - RETURN: master traceback for aseqs          (may pass NULL)
 *           
 * Return:   1 on success, 0 on failure. 
 *           ret_cm is alloc'ed here and must be free'd by the caller.
 */
int
Maxmodelmaker(char **aseqs, AINFO  *ainfo, int nseq, double gapthresh, struct prior_s *prior, 
	      double *ret_ssinfo, struct cm_s **ret_cm, struct trace_s **ret_mtr)       
{
  int            **aseqsT;      /* transpose of alignment, [1..alen+1][0..nseq-1]    */
  int             *mscore;      /* emission scores of single columns as MATR or MATL */
  double          *gapcount;    /* weighted count of gap occurrence in each column   */
  struct maxmx_s **mmx;         /* saved score and traceback pointer matrix          */
  struct cm_s     *cm;          /* RETURN: new model                                 */
  int              nodes;	/* size of new model in nodes                        */
  struct trace_s  *mtr;         /* consensus master traceback                        */
  struct trace_s  *tr;          /* individual "fake" traceback                       */
  struct trmem_s  *pool;	/* memory pool for traceback                         */
  int              idx;		/* counter for sequences                             */
  float           *wt;		/* array of 0..nseq-1 weights                        */
  double           ssinfo;      /* information content                               */

  /* Set up an array of sequence weights
   */
  if ((wt = (float *) malloc (sizeof(float) * nseq)) == NULL)
    Die("malloc failed");
  for (idx = 0; idx < nseq; idx++)
    wt[idx] = (ainfo->sqinfo[idx].flags & SQINFO_WGT) ? 
      ainfo->sqinfo[idx].weight : 1.0;

  /* Transpose (and shift by 1 column index) aseqs[0..nseq-1][0..alen-1] to 
   * aseqsT[1..alen+1][0..nseq-1].
   */
  transpose_alignment(aseqs, ainfo->alen, nseq, &aseqsT);

  /* Pre-calculate expected singlet emission scores for each column;
   * also, pick up a weighted count of gap occurrences in each column.
   */
  singlet_emissions(aseqsT, wt, ainfo->alen, nseq, prior, &mscore, &gapcount);
  

  /* Allocate the scoring matrix. It is j=0..alen rows by i=1..j+1 columns;
   * i.e., a lower left diagonal matrix for 1..alen, with an extra off-diagonal
   * (j+1,j) for boundary conditions.
   */
  mmx = alloc_maxmx(ainfo->alen);

  /* Initialize the off-diagonal (j+1,j) and the diagonal (j,j), to
   * set our boundary conditions.
   */
  init_maxmx(mmx, nseq, ainfo->alen, prior, mscore, gapcount);

  /* The heart of the calculation: recursively calculate scores for all
   * subsequences i..j, and save traceback pointers.
   */
  recurse_maxmx(aseqsT, wt, ainfo->alen, nseq, prior, mscore, gapcount, gapthresh, mmx);
#ifdef DEBUG
  print_mmx(mmx,ainfo->alen);
#endif

  /* Now we know the info content
   */
  ssinfo = (double) (mmx[ainfo->alen][0].sc[MATP_NODE] / INTPRECISION);

  /* Traceback, constructing a consensus tree.
   */
  trace_maxmx(mmx, ainfo->alen, &mtr);
/*  PrintTrace(stdout, mtr); */

  /* Count nodes in the consensus tree and number them
   */
  NumberMasterTrace(mtr, &nodes);
  
  /* Create a new model
   */
  if ((cm = AllocCM(nodes)) == NULL)
    Die("failed to allocate for new model of %d nodes\n", nodes);
  TopofyNewCM(cm, mtr);
  
  /* For each sequence, construct an individual "fake" traceback
   * using the master, and count it into the new model.
   */
  for (idx = 0; idx < nseq; idx++)
    {
      Transmogrify(mtr, aseqs[idx], &tr, &pool);

      if (! TraceCount(cm, aseqs[idx], wt[idx], tr))
	Die("TraceCount() failed");

      FreeTrace(tr, pool);
    }
#ifdef DEBUG
  printf("Checking model after trace counts\n");
  if (! VerifyCM(cm))
    Die("Verification fails -- baaaad model\n");
#endif

  /* And, finally, convert the new model to probabilities. 
   * There. Wasn't that simple?
   */
  ProbifyCM(cm, prior);

#ifdef DEBUG
  printf("Checking model after probify\n");
  if (! VerifyCM(cm))
    Die("Verification fails -- baaaad model\n");
#endif


  free(mscore);
  free(gapcount);
  Free2DArray(mmx, ainfo->alen+1);
  Free2DArray(aseqsT, ainfo->alen+2);
  free(wt);

  if (ret_mtr != NULL)    *ret_mtr    = mtr; else FreeTrace(mtr, NULL);
  if (ret_ssinfo != NULL) *ret_ssinfo = ssinfo;
  if (ret_cm != NULL)     *ret_cm     = cm;  else FreeCM(cm);
  return 1;
}
 



/* Function: alloc_maxmx()
 * 
 * Purpose:  allocate the storage matrix. It is a lower left diagonal
 *           matrix with inverted indexing, mmx[j][i], i <= j+1.
 *           
 * Args:     alen - width of alignment
 *                  
 * Return:   mmx: allocated storage matrix. Can be free'd by Free2DArray(mmx, alen+1).
 */
static struct maxmx_s **
alloc_maxmx(int alen)
{
  struct maxmx_s **mmx;
  int            i, j, y;

  if ((mmx = (struct maxmx_s **) malloc ((alen+1) * sizeof(struct maxmx_s))) == NULL)
    Die("malloc failed");
  for (j = 0; j <= alen; j++)
    if ((mmx[j] = (struct maxmx_s *) malloc ((j+2) * sizeof(struct maxmx_s))) == NULL)
      Die("malloc failed");

  /* Set the whole matrix so that scores are NEGINFINITY and
   * each cell's traceback pointers point at itself. 
   */
  for (j = 0; j <= alen; j++)
    for (i = 0; i <= j+1; i++)
      {
	for (y = 0; y < NODETYPES-1; y++)
	  mmx[j][i].sc[y] = NEGINFINITY;

	mmx[j][i].matp_i2    = (short) i;
	mmx[j][i].matp_j2    = (short) j;
	mmx[j][i].matp_ftype = (char)  MATP_NODE;
	mmx[j][i].matl_i2    = (short) i;
	mmx[j][i].matl_ftype = (char)  MATL_NODE;
	mmx[j][i].matr_j2    = (short) j;
	mmx[j][i].matr_ftype = (char)  MATR_NODE;
	mmx[j][i].begl_ftype = (char)  BEGINL_NODE;
	mmx[j][i].begr_i2    = (short) i;
	mmx[j][i].begr_ftype = (char)  BEGINR_NODE;
	mmx[j][i].bifurc_mid = (short) i;
      }
  return mmx;
}


/* Function: init_maxmx()
 * 
 * Purpose:  Initialize the scoring matrix. The offdiagonal j+1,j and the
 *           diagonal j,j are initialized.
 *           
 *           In the offdiagonal, we use BIFURC to represent END, and set
 *           its score to zero.
 *           
 *           In the diagonal, it doesn't matter whether we use MATL or
 *           MATR to generate single symbols; we use MATL for an implementation-
 *           specific reason (if we use MATR, we need an extra row for i=0,j=-1)
 *           MATL's get calculated scores and traceback pointers to 
 *           (j+1,j,END). BEGINL, BEGINR are also calculated here.
 *           MATR, MATP, BIFURC are left at NEGINFINITY.
 *           
 * Args:     mmx:      saved score/traceback pointer matrix
 *           nseq:     number of sequences in alignment
 *           alen:     number of columns in alignment
 *           prior:    prior probability distributions
 *           mscore:   singlet emission costs for each column
 *           gapcount: weighted counts of gaps in each column
 *           
 * Return:   (void). maxmx() is initialized.          
 */
static void 
init_maxmx(struct maxmx_s **mmx,
	   int              nseq,
	   int              alen,
	   struct prior_s  *prior,
	   int             *mscore,
	   double           *gapcount)
{
  int    j;
  int    fromstate, tostate;
  double trans[STATETYPES][STATETYPES];	/* state transition matrix */

  /* Do the offdiagonal (j+1,j).
   * set BIFURC/END alignment costs to zero.
   * Everything else is left at -oo.
   */
  for (j = 0; j <= alen; j++)
    mmx[j][j+1].sc[END_NODE]    = 0;

  /* Do the diagonal (j,j).
   * MATL is calculated; MATR is the same; then BEGINL, BEGINR are calculated.
   */
  for (j = 1; j <= alen; j++)
    {
      /* Make a transition matrix for MATL -> END for this (j,j)
       */
      for (fromstate = 0; fromstate < STATETYPES; fromstate++)
	for (tostate = 0; tostate < STATETYPES; tostate++)
	  trans[fromstate][tostate] = 0.0;
      trans[MATL_ST][END_ST] = (double) nseq - gapcount[j];
      trans[DEL_ST][END_ST]  = gapcount[j];
      ProbifyTransitionMatrix(trans, MATL_NODE, END_NODE, prior);

      /* Score = sum P(j | MATL) + sum T(END | j,j,(DEL|MATL))
       * Set traceback pointers.
       */
      mmx[j][j].sc[MATL_NODE] = mscore[j] +
	(int) (INTPRECISION * 
	       ((log(trans[MATL_ST][END_ST]) * ((double) nseq - gapcount[j])) +
		(log(trans[DEL_ST][END_ST]) * gapcount[j]))
	      );
      mmx[j][j].matl_i2    = j+1;
      mmx[j][j].matl_ftype = END_NODE;

      /* MATR_NODE scores are exactly the same as MATL_NODE on
       * the diagonal
       */
      mmx[j][j].sc[MATR_NODE] = mmx[j][j].sc[MATL_NODE];
      mmx[j][j].matr_j2       = j-1;
      mmx[j][j].matr_ftype    = END_NODE;

      /* Calculate BEGINL -> MATL;
       * set traceback pointer.
       */
      for (fromstate = 0; fromstate < STATETYPES; fromstate++)
	for (tostate = 0; tostate < STATETYPES; tostate++)
	  trans[fromstate][tostate] = 0.0;
      trans[BEGIN_ST][MATL_ST] = (double) nseq - gapcount[j];
      trans[BEGIN_ST][DEL_ST]  = gapcount[j];
      ProbifyTransitionMatrix(trans, BEGINL_NODE, MATL_NODE, prior);
      mmx[j][j].sc[BEGINL_NODE] = mmx[j][j].sc[MATL_NODE] 
	+ (int) (INTPRECISION * 
		 ((log(trans[BEGIN_ST][MATL_ST]) * ((double) nseq - gapcount[j])) +
		  (log(trans[BEGIN_ST][DEL_ST]) * gapcount[j]))
	      );
      mmx[j][j].begl_ftype      = MATL_NODE;

      /* Make a transition matrix for BEGINR -> MATL for this (j,j)
       */
      for (fromstate = 0; fromstate < STATETYPES; fromstate++)
	for (tostate = 0; tostate < STATETYPES; tostate++)
	  trans[fromstate][tostate] = 0.0;
      trans[BEGIN_ST][DEL_ST]  = gapcount[j];
      trans[BEGIN_ST][MATL_ST] = (double) nseq - gapcount[j];
      ProbifyTransitionMatrix(trans, BEGINR_NODE, MATL_NODE, prior);
      
      /* Score for (j,j,BEGINR_NODE), and set traceback pointers
       */
      mmx[j][j].sc[BEGINR_NODE] = mmx[j][j].sc[MATL_NODE] + 
	(int) (INTPRECISION * 
	       ((log(trans[BEGIN_ST][MATL_ST]) * ((double) nseq - gapcount[j])) +
		(log(trans[BEGIN_ST][DEL_ST]) * gapcount[j]))
	      );
      mmx[j][j].begr_i2    = j;
      mmx[j][j].begr_ftype = MATL_NODE;
    }
}


/* Function: recurse_maxmx()
 * 
 * Purpose:  Recursion calculations of the maximum likelihood CM construction
 *           algorithm.
 *           
 * Args:     aseqsT  - transposed alignment; [1..alen+1][0..nseq-1]
 *           weights - weights assigned to each sequence, usually 1.0
 *           alen    - number of columns in alignment
 *           nseq    - number of aseqs
 *           prior   - structure containing prior probability distributions
 *           mscore  - singlet match emission costs
 *           mmx     - scoring matrix we fill in; mmx[j][i], lower diagonal
 *           
 * Return:   (void). mmx is filled with scores and traceback pointers.
 */
static void
recurse_maxmx(int    **aseqsT, 
	      float   *weights, 
	      int      alen, 
	      int      nseq, 
	      struct prior_s *prior, 
	      int     *mscore, 
	      double  *gapcount,
	      double   gapthresh,
	      struct maxmx_s **mmx)
{
  int **insl_accum;     /* accumulated INSL counts for all starting positions, all seqs */
  int **insr_accum;     /* accumulated INSR counts for all starting positions, all seqs */
  double tmaster[STATETYPES][STATETYPES]; /* master copy of state transition table, counts */
  double tcounts[STATETYPES][STATETYPES]; /* state transition table, counts */
  double tmx[STATETYPES][STATETYPES]; /* state transition probabilities post-regularization */
  int    tonode;		/* node type to connect to       */
  int    i,j;                   /* current cell column, row      */
  int    i2, j2;		/* i', j': cell to connect to    */
  int    idx;			/* counter for sequences         */
  int    sc;			/* temp variable holding a score */
  int    mid;			/* midpoint for a bifurcation    */

  if ((insr_accum = (int **) malloc ((alen+1) * sizeof(int *))) == NULL ||
      (insl_accum = (int **) malloc ((alen+2) * sizeof(int *))) == NULL)
    Die("malloc failed");
  for (i = 0; i <= alen; i++)
    if ((insr_accum[i] = (int *) malloc (nseq * sizeof(int))) == NULL)
      Die("malloc failed");
  for (i = 0; i <= alen+1; i++)
    if ((insl_accum[i] = (int *) malloc (nseq * sizeof(int))) == NULL)
      Die("malloc failed");

  gapthresh *= (double) nseq; /* scale gapthresh to be comparable to counts in gapcount array */

  /* Initialize insr_accum. (Vertical, j insertion accumulator)
   * insr_accum contains, for each sequence, a count of how many symbols
   * must be inserted to get from any row fromj to the current row (exclusive
   * of fromj and the current column).
   * insr_accum is therefore [0..alen-1][0..nseq-1]
   */
  for (j2 = 1; j2 <= alen; j2++)
    for (idx = 0; idx < nseq; idx++)
      insr_accum[j2][idx] = 0;
  for (idx = 0; idx < nseq; idx++)
    if (aseqsT[1][idx] >= 0)
      insr_accum[0][idx] = 1;

  for (j = 2; j <= alen; j++)
    {
      /* Initialize insl_accum (horizontal, i insertion) array each time we start
       * a new row. insl_accum contains, for each sequence, how many symbols
       * must be inserted to get from any column i2 to the current column,
       * exclusive. insl_accum is therefore [1..alen+1][0..nseq]
       */
      for (i2 = 1; i2 <= j+1; i2++)
	for (idx = 0; idx < nseq; idx++)
	  insl_accum[i2][idx] = 0;
      for (idx = 0; idx < nseq; idx++)
	if (aseqsT[j][idx] >= 0)
	  insl_accum[j+1][idx]++;
      
      for (i = j-1; i > 0; i--)
	{

	  /* BIFURC: explain i,j as sum of i,mid,BEGINL + mid+1,j,BEGINR
	   * i <= mid <= j
	   */
	  for (mid = i; mid <= j; mid++)
	    {
	      sc = mmx[mid][i].sc[BEGINL_NODE] + mmx[j][mid+1].sc[BEGINR_NODE];
	      if (sc > mmx[j][i].sc[BIFURC_NODE])
		{
		  mmx[j][i].sc[BIFURC_NODE] = sc;
		  mmx[j][i].bifurc_mid      = mid;
		}
	    }

	  /* MATP: Score subsequence i,j, given that i,j are emitted by MATP
           * Look at all possible connections i'j': i < i' < j, i < j' < j,
           * i' <= j'+1
           * Could be i < i' <= j' < j: the extra i' <= j'+1 condition
           * allows for checking all the ways of generating i'..j'
           * as entirely insertion, and we may just as well define a
           * default...
	   */
	  for (j2 = j-1; j2 >= i && j - j2 - 1 <= MAXINSERT ; j2--)
/*	  for (j2 = j-1; j2 >= i; j2--) */
	    {
	      for (i2 = i+1; i2 <= j2+1 && i2 - i - 1 <= MAXINSERT; i2++) 
/*	      for (i2 = i+1; i2 <= j2+1; i2++) */
		{
		  if (mmx[j2][i2].sc[MATP_NODE]   < mmx[j][i].sc[MATP_NODE] &&
		      mmx[j2][i2].sc[MATL_NODE]   < mmx[j][i].sc[MATP_NODE] &&
		      mmx[j2][i2].sc[MATR_NODE]   < mmx[j][i].sc[MATP_NODE] &&
		      mmx[j2][i2].sc[BIFURC_NODE] < mmx[j][i].sc[MATP_NODE])
		    continue;

		  frommatp_transtable(aseqsT, weights, nseq, i,j, i2, j2, insl_accum[i2],
				      insr_accum[j2], tmaster);
		  
		  for (tonode = 0; tonode < 4; tonode++)
		    {
		      if (i2 >  j2 && tonode != BIFURC_NODE) continue;
		      if (i2 == j2 && tonode == MATP_NODE)   continue;

		      if (mmx[j2][i2].sc[tonode] < mmx[j][i].sc[MATP_NODE]) continue;

		      switch (tonode) {
		      case MATP_NODE:   to_matp_transtable(tmaster, tcounts);   break;
		      case MATL_NODE:   to_matl_transtable(tmaster, tcounts);   break;
		      case MATR_NODE:   to_matr_transtable(tmaster, tcounts);   break;
		      case BIFURC_NODE: to_bifurc_transtable(tmaster, tcounts); break;
		      default: Die("Gotcha. MATP, MATL, MATR, BIFURC nodes must be numbered 0..3");
		      }
		      
		      copy_transtable(tmx, tcounts);
		      ProbifyTransitionMatrix(tmx, MATP_NODE, tonode, prior);
		      
		      sc = dot_score((double *) tcounts, (double *) tmx, STATETYPES*STATETYPES)
			+ mmx[j2][i2].sc[tonode];
		      
		      if (sc > mmx[j][i].sc[MATP_NODE])
			{
			  mmx[j][i].sc[MATP_NODE] = sc;
			  mmx[j][i].matp_i2       = i2;
			  mmx[j][i].matp_j2       = j2;
			  mmx[j][i].matp_ftype    = tonode;
			}
		    }		
		  if (gapcount[i2] <= gapthresh) break;
		}
	      if (gapcount[j2] <= gapthresh) break;
	    }
	  mmx[j][i].sc[MATP_NODE] += pair_emissioncost(aseqsT[i], aseqsT[j], weights, nseq, prior);

	  /* MATR: i,j is accounted for by emitting j and connecting 
	   * to some i,j2.
	   */
	  for (j2 = j-1; j2 >= i-1; j2--)
	    {
	      if (mmx[j2][i].sc[MATP_NODE]   < mmx[j][i].sc[MATR_NODE] &&
		  mmx[j2][i].sc[MATL_NODE]   < mmx[j][i].sc[MATR_NODE] &&
		  mmx[j2][i].sc[MATR_NODE]   < mmx[j][i].sc[MATR_NODE] &&
		  mmx[j2][i].sc[BIFURC_NODE] < mmx[j][i].sc[MATR_NODE])
		continue;
	      frommatr_transtable(aseqsT, weights, nseq, i, j, j2, insr_accum[j2], tmaster);
	      
	      for (tonode = 0; tonode < 4; tonode++)
		{
		  if (i >  j2 && tonode != BIFURC_NODE) continue;
		  if (i == j2 && tonode == MATP_NODE)   continue;

		  if (mmx[j2][i].sc[tonode] < mmx[j][i].sc[MATR_NODE]) continue;

		  switch (tonode) {
		  case MATP_NODE:   to_matp_transtable(tmaster, tcounts);   break;
		  case MATL_NODE:   to_matl_transtable(tmaster, tcounts);   break;
		  case MATR_NODE:   to_matr_transtable(tmaster, tcounts);   break;
		  case BIFURC_NODE: to_bifurc_transtable(tmaster, tcounts); break;
		  default: Die("Gotcha. MATP, MATL, MATR, BIFURC nodes must be numbered 0..3");
		  }
		  
		  copy_transtable(tmx, tcounts);
		  ProbifyTransitionMatrix(tmx, MATR_NODE, tonode, prior);

		  sc = dot_score((double *) tcounts, (double *) tmx, STATETYPES*STATETYPES)
		    + mmx[j2][i].sc[tonode];
		  if (sc > mmx[j][i].sc[MATR_NODE])
		    {
		      mmx[j][i].sc[MATR_NODE] = sc;
		      mmx[j][i].matr_j2       = j2;
		      mmx[j][i].matr_ftype    = tonode;
		    }
		}		
	      if (gapcount[j2] <= gapthresh) break;
	    }
	  mmx[j][i].sc[MATR_NODE] += mscore[j];


	  /* MATL: account for i,j by emitting i and connecting to some (i2,j)
	   */
	  for (i2 = i+1; i2 <= j+1; i2++)
	    {
	      if (mmx[j][i2].sc[MATP_NODE]   < mmx[j][i].sc[MATL_NODE] &&
		  mmx[j][i2].sc[MATL_NODE]   < mmx[j][i].sc[MATL_NODE] &&
		  mmx[j][i2].sc[MATR_NODE]   < mmx[j][i].sc[MATL_NODE] &&
		  mmx[j][i2].sc[BIFURC_NODE] < mmx[j][i].sc[MATL_NODE])
		continue;

	      frommatl_transtable(aseqsT, weights, nseq, i, j, i2, insl_accum[i2], tmaster);

	      for (tonode = 0; tonode < 4; tonode++)
		{
		  if (i2 >  j && tonode != BIFURC_NODE) continue;
		  if (i2 == j && tonode == MATP_NODE)   continue;
		  
		  if (mmx[j][i2].sc[tonode] < mmx[j][i].sc[MATL_NODE]) continue;

		  switch (tonode) {
		  case MATP_NODE:   to_matp_transtable(tmaster, tcounts);   break;
		  case MATL_NODE:   to_matl_transtable(tmaster, tcounts);   break;
		  case MATR_NODE:   to_matr_transtable(tmaster, tcounts);   break; 
		  case BIFURC_NODE: to_bifurc_transtable(tmaster, tcounts); break;
		  default: Die("Gotcha. MATP, MATL, MATR, BIFURC nodes must be numbered 0..3");
		  }
		  
		  copy_transtable(tmx, tcounts);
		  ProbifyTransitionMatrix(tmx, MATL_NODE, tonode, prior);

		  sc = dot_score((double *) tcounts, (double *) tmx, STATETYPES*STATETYPES) + 
		    mmx[j][i2].sc[tonode];

		  if (sc > mmx[j][i].sc[MATL_NODE])
		    {
		      mmx[j][i].sc[MATL_NODE] = sc;
		      mmx[j][i].matl_i2       = i2;
		      mmx[j][i].matl_ftype    = tonode;
		    }
		}		
	      if (gapcount[i2] <= gapthresh) break;
	    }
	  mmx[j][i].sc[MATL_NODE] += mscore[i];


	  /* bump insl_accum: add column i to horizontal accumulator as insertion.
	   */
	  for (i2 = i+1; i2 <= j+1; i2++)
	    for (idx = 0; idx < nseq; idx++)
	      if (aseqsT[i][idx] >= 0)
		insl_accum[i2][idx]++;


	  /* BEGINR: has an INSL state, so it can connect to any
	   * (i2,j) *inclusive* of (i,j) -- that's why we just bumped
	   * the insl_accum insert counters
	   */
	  for (i2 = i; i2 <= j+1; i2++)
	    {
	      if (mmx[j][i2].sc[MATP_NODE]   < mmx[j][i].sc[BEGINR_NODE] &&
		  mmx[j][i2].sc[MATL_NODE]   < mmx[j][i].sc[BEGINR_NODE] &&
		  mmx[j][i2].sc[MATR_NODE]   < mmx[j][i].sc[BEGINR_NODE] &&
		  mmx[j][i2].sc[BIFURC_NODE] < mmx[j][i].sc[BEGINR_NODE])
		continue;

	      frombeginr_transtable(aseqsT, weights, nseq, j, i2, insl_accum[i2], tmaster);

	      for (tonode = 0; tonode < 4; tonode++)
		{
		  if (i2 >  j && tonode != BIFURC_NODE) continue;
		  if (i2 == j && tonode == MATP_NODE)   continue;

		  if (mmx[j][i2].sc[tonode] < mmx[j][i].sc[BEGINR_NODE]) continue;

		  switch (tonode) {
		  case MATP_NODE:   to_matp_transtable(tmaster, tcounts);   break;
		  case MATL_NODE:   to_matl_transtable(tmaster, tcounts);   break;
		  case MATR_NODE:   to_matr_transtable(tmaster, tcounts);   break;
		  case BIFURC_NODE: to_bifurc_transtable(tmaster, tcounts); break;
		  default: Die("Gotcha. MATP, MATL, MATR, BIFURC nodes must be numbered 0..3");
		  }
		  
		  copy_transtable(tmx, tcounts);
		  ProbifyTransitionMatrix(tmx, BEGINR_NODE, tonode, prior);

		  sc = dot_score((double *) tcounts, (double *) tmx, STATETYPES*STATETYPES)
		    + mmx[j][i2].sc[tonode];

		  if (sc > mmx[j][i].sc[BEGINR_NODE])
		    {
		      mmx[j][i].sc[BEGINR_NODE] = sc;
		      mmx[j][i].begr_i2         = i2;
		      mmx[j][i].begr_ftype      = tonode;
		    }
		}		
	      if (gapcount[i2] <= gapthresh) break;
	    }

	      
	  /* BEGINL: has no inserts, so must connect to (i,j)
	   */
	  frombeginl_transtable(aseqsT, weights, nseq, i, j, tmaster);
	  for (tonode = 0; tonode < 4; tonode++)
	    {
	      if (mmx[j][i].sc[MATP_NODE]   < mmx[j][i].sc[BEGINL_NODE] &&
		  mmx[j][i].sc[MATL_NODE]   < mmx[j][i].sc[BEGINL_NODE] &&
		  mmx[j][i].sc[MATR_NODE]   < mmx[j][i].sc[BEGINL_NODE] &&
		  mmx[j][i].sc[BIFURC_NODE] < mmx[j][i].sc[BEGINL_NODE])
		continue;

	      switch (tonode) {
	      case MATP_NODE:   to_matp_transtable(tmaster, tcounts);   break;
	      case MATL_NODE:   to_matl_transtable(tmaster, tcounts);   break;
	      case MATR_NODE:   to_matr_transtable(tmaster, tcounts);   break;
	      case BIFURC_NODE: to_bifurc_transtable(tmaster, tcounts); break;
	      default: Die("Gotcha. MATP, MATL, MATR, BIFURC nodes must be numbered 0..3");
	      }
		  
	      copy_transtable(tmx, tcounts);
	      ProbifyTransitionMatrix(tmx, BEGINL_NODE, tonode, prior);

	      sc = dot_score((double *) tcounts, (double *) tmx, STATETYPES*STATETYPES)
		+ mmx[j][i].sc[tonode];

	      if (sc > mmx[j][i].sc[BEGINL_NODE])
		{
		  mmx[j][i].sc[BEGINL_NODE] = sc;
		  mmx[j][i].begl_ftype      = tonode;
		}
	    }


	}			/* end loop over columns i */

      /* bump insr_accum: add row j to vertical accumulator as insertion.
       */
      for (j2 = 0; j2 < j; j2++)
	for (idx = 0; idx < nseq; idx++)
	  if (aseqsT[j][idx] >= 0)
	    insr_accum[j2][idx]++;

    }				/* end loop over rows j */
  

  /* Termination. ROOT can connect anywhere.
   * We hack here: ROOT alignment info is stored in the cell mmx[alen][0]
   * (otherwise, the i==0 column is unused), and score/traceback info is kept as if
   * for MATP_NODE (because ROOT and MATP have similar traceback requirements,
   * since they permit inserts on both sides.)
   * 
   */
  for (j2 = alen; j2 >= 0; j2--)
    {
      for (i2 = 1; i2 <= j2+1; i2++)
	{
	  if (mmx[j2][i2].sc[MATP_NODE]   < mmx[alen][0].sc[MATP_NODE] &&
	      mmx[j2][i2].sc[MATL_NODE]   < mmx[alen][0].sc[MATP_NODE] &&
	      mmx[j2][i2].sc[MATR_NODE]   < mmx[alen][0].sc[MATP_NODE] &&
	      mmx[j2][i2].sc[BIFURC_NODE] < mmx[alen][0].sc[MATP_NODE])
	    continue;
	  
	  fromroot_transtable(aseqsT, weights, nseq, i2, j2, insl_accum[i2],
			      insr_accum[j2], tmaster);
	  
	  for (tonode = 0; tonode < 4; tonode++)
	    {
	      if (i2 >  j2 && tonode != BIFURC_NODE) continue;
	      if (i2 == j2 && tonode == MATP_NODE)   continue;
	      if (mmx[j2][i2].sc[tonode] < mmx[alen][0].sc[MATP_NODE]) continue;
	      
	      switch (tonode) {
	      case MATP_NODE:   to_matp_transtable(tmaster, tcounts);   break;
	      case MATL_NODE:   to_matl_transtable(tmaster, tcounts);   break;
	      case MATR_NODE:   to_matr_transtable(tmaster, tcounts);   break;
	      case BIFURC_NODE: to_bifurc_transtable(tmaster, tcounts); break;
	      default: Die("Gotcha. MATP, MATL, MATR, BIFURC nodes must be numbered 0..3");
	      }
	      
	      copy_transtable(tmx, tcounts);
	      ProbifyTransitionMatrix(tmx, ROOT_NODE, tonode, prior);
	      
	      sc = dot_score((double *) tcounts, (double *) tmx, STATETYPES*STATETYPES)
		+ mmx[j2][i2].sc[tonode];
	      if (sc > mmx[alen][0].sc[MATP_NODE])
		{
		  mmx[alen][0].sc[MATP_NODE] = sc;
		  mmx[alen][0].matp_i2       = i2;
		  mmx[alen][0].matp_j2       = j2;
		  mmx[alen][0].matp_ftype    = tonode;
		}
	    }		
	  if (gapcount[i2] <= gapthresh) break;
	}
      if (gapcount[j2] <= gapthresh) break;
    }

  Free2DArray(insr_accum, alen+1);
  Free2DArray(insl_accum, alen+2);
}



/* Function: trace_maxmx()
 * 
 * Purpose:  Traceback of the filled matrix. Constructs a consensus 
 *           tree structure. The trace_s structures are only partially
 *           used: emitl and emitr hold 0..alen-1 column coords of emitted
 *           columns (even if not responsible for the emission); nodeidx 
 *           is unused; type holds a *node* type, not a state type.
 *           
 *           Note that the alignment of ROOT_NODE has been stored in
 *           mmx[alen][0] like a MATP_NODE, due to some convenient
 *           hacking by the recursion routine.
 *      
 *           The mmx[] scoring matrix traceback pointers are 1..alen coords.
 *           They have to be converted to 0..alen-1 for the traceback tree.
 *           
 * Args:     mmx     - filled scoring matrix: mmx[j][i], lower diagonal
 *           alen    - width of alignment
 *           ret_mtr - RETURN: master (consensus) traceback          
 * 
 * Return:   (void). ret_mtr must be free'd by the caller.
 */
static void
trace_maxmx(struct maxmx_s **mmx,
	    int              alen,
	    struct trace_s **ret_mtr)
{
  struct trace_s *mtr;	        /* master (consensus) traceback               */
  struct trace_s *curr_mtr;     /* current node on traceback tree             */
  struct tracestack_s *dolist;  /* pushdown stack for traversing mtr          */
  int i,j;			/* coords to connect to, next trace ssegment  */
  int nxti, nxtj;

  InitTrace(&mtr, NULL);
  dolist = InitTracestack();
  
  /* Initialization. First attach a root node. Then, trace first segment,
   * attach it, and start the tracestack to-do list with it.
   */
  curr_mtr = AttachTrace(mtr, NULL, 0, alen-1, 0, ROOT_NODE); 
  PushTracestack(dolist, AttachTrace(curr_mtr, NULL, mmx[alen][0].matp_i2 -1,
				     mmx[alen][0].matp_j2 -1, 
				     0, mmx[alen][0].matp_ftype));
		 

  while ((curr_mtr = PopTracestack(dolist)) != NULL)
    {
      i = curr_mtr->emitl + 1;	/* i,j are 1..alen */
      j = curr_mtr->emitr + 1;

      /* avoid dummy END node on trace tree leaves, and
       * avoid tracing back from off-diagonal */
      if (curr_mtr->nxtl == NULL || i > j)
	continue;
  
      switch (curr_mtr->type) {
      case MATP_NODE:
	nxti = mmx[j][i].matp_i2 - 1;
	nxtj = mmx[j][i].matp_j2 - 1;
	if (nxti <= nxtj)
	  PushTracestack(dolist, AttachTrace(curr_mtr, NULL, nxti, nxtj, 0, mmx[j][i].matp_ftype));
	else
	  { curr_mtr->nxtl->emitl = nxti; curr_mtr->nxtl->emitr = nxtj; }
	break;

      case MATL_NODE:
	nxti = mmx[j][i].matl_i2 - 1;
	nxtj = j-1;
	if (nxti <= nxtj)
	  PushTracestack(dolist, AttachTrace(curr_mtr, NULL, nxti, nxtj, 0, mmx[j][i].matl_ftype));
	else
	  { curr_mtr->nxtl->emitl = nxti; curr_mtr->nxtl->emitr = nxtj; }
	break;

      case MATR_NODE:
	nxti = i-1;
	nxtj = mmx[j][i].matr_j2 -1;
	if (nxti <= nxtj)
	  PushTracestack(dolist, AttachTrace(curr_mtr, NULL, nxti, nxtj, 0, mmx[j][i].matr_ftype));
	else
	  { curr_mtr->nxtl->emitl = nxti; curr_mtr->nxtl->emitr = nxtj; }
	break;

      case BIFURC_NODE:		/* BIFURC must attach right side first */
	PushTracestack(dolist, AttachTrace(curr_mtr, NULL, mmx[j][i].bifurc_mid, j -1,
					   0, BEGINR_NODE));
	PushTracestack(dolist, AttachTrace(curr_mtr, NULL, i -1, mmx[j][i].bifurc_mid -1,
					   0, BEGINL_NODE));
	break;

      case BEGINL_NODE:
	nxti = i-1;
	nxtj = j-1;
	if (nxti <= nxtj)
	  PushTracestack(dolist, AttachTrace(curr_mtr, NULL, nxti, nxtj, 0, mmx[j][i].begl_ftype));
	else
	  { curr_mtr->nxtl->emitl = nxti; curr_mtr->nxtl->emitr = nxtj; }
	break;

      case BEGINR_NODE:
	nxti = mmx[j][i].begr_i2 - 1;
	nxtj = j-1;
	if (nxti <= nxtj)
	  PushTracestack(dolist, AttachTrace(curr_mtr, NULL, nxti, nxtj, 0, mmx[j][i].begr_ftype));
	else
	  { curr_mtr->nxtl->emitl = nxti; curr_mtr->nxtl->emitr = nxtj; }
	break;

      default: Die("Invalid node type %d", curr_mtr->type);
      }
    }

  FreeTracestack(dolist);
  *ret_mtr = mtr;
}


/* Function: transpose_alignment()
 * 
 * Purpose:  Alignments are indexed [seqidx][position]; it turns out to
 *           be more convenient here to index them as [position][seqidx],
 *           because of memory access patterns. This transpose also
 *           lets us switch to a 1..alen indexing scheme for the alignment
 *           columns, from 0..alen-1; this is important for implementing
 *           boundary conditions properly in the scoring matrix. And
 *           finally, we store the symbols as indices to save time in
 *           further lookups: -1 for gaps, 0..3 for ACGU (or 0..19 for aminos)
 *           
 * Args:     aseqs      - flushed sequence alignment, each seq indexed 0..alen-1
 *           alen       - number of columns in alignment
 *           nseq       - number of sequences
 *           prior      - contains alphabet info
 *           ret_aseqsT - RETURN: transposed alignment, [1..alen][0..nseq-1]
 *           
 * Return:   (void). ret_aseqsT is malloc'ed here and must be free'd 
 *           by caller.          
 */
static void
transpose_alignment(char  **aseqs,
		    int     alen,
		    int     nseq,
		    int  ***ret_aseqsT)
{
  int  **aseqsT;
  int    acol;
  int    seqidx;

  if ((aseqsT = (int **) malloc ((alen+2) * sizeof(int *))) == NULL)
    Die("malloc failed");
  for (acol = 0; acol <= alen+1; acol++)
    if ((aseqsT[acol] = (int *) malloc (nseq * sizeof(int))) == NULL)
      Die("malloc failed");

				/* "guard" columns 0 and alen+1 */
  for (seqidx = 0; seqidx < nseq; seqidx++)
    {
      aseqsT[0][seqidx]      = -1;
      aseqsT[alen+1][seqidx] = -1;
    }

				/* gaps are assigned value -1 in aseqsT */
  for (seqidx = 0; seqidx < nseq; seqidx++)
    for (acol = 0; acol < alen; acol++)
      aseqsT[acol+1][seqidx] = isgap(aseqs[seqidx][acol]) ? 
	-1 : SymbolIndex(aseqs[seqidx][acol]);

  *ret_aseqsT = aseqsT;
}





/* Function: singlet_emissions()
 * 
 * Purpose:  Count emission statistics for all the columns of 
 *           a multiple alignment; calculate and return an
 *           array of emission scores for each column emitted
 *           by MATL. (The caller can assume that MATR is scored
 *           the same way.)
 *           
 * Args:     aseqsT        - [1..alen][0..nseqs-1] transpose of sequence alignment          
 *           weights       - weights on sequences (usually just 1.0 for each)
 *           alen          - number of columns in alignment
 *           nseq          - number of sequences    
 *           prior         - prior probability distributions, and alphabet info
 *           ret_mscore    - RETURN: array of singlet emission costs
 *           ret_gapcount  - RETURN: weighted counts of gaps in each column
 *           
 * Return:   (void). ret_mscore is passed back; it is malloc'ed here
 *           and must be free'd by the caller.
 */          
static void
singlet_emissions(int    **aseqsT,
		  float   *weights,
		  int      alen,
		  int      nseq,
		  struct prior_s *prior,
		  int    **ret_mscore,
		  double **ret_gapcount)
{
  double **emcounts;
  double  *gapcount;
  int     *mscore;
  double   emvec[ALPHASIZE];
  int i;
  int idx;
  int sym;

  /* Allocations
   */
  if ((emcounts = (double **) malloc ((alen+1) *             sizeof(double *))) == NULL ||
      (mscore   = (int *)     malloc ((alen+1) *             sizeof(int)))      == NULL ||
      (gapcount = (double *)  malloc ((alen+2) *             sizeof(double)))   == NULL)
    Die("malloc failed");
  for (i = 0; i <= alen; i++)
    if ((emcounts[i] = (double *) malloc (ALPHASIZE * sizeof(double))) == NULL)
      Die("malloc failed");

  /* Count symbol occurrences in each column
   */
  for (i = 1; i <= alen; i++)
    {
      gapcount[i] = 0.0;
      for (sym = 0; sym < ALPHASIZE; sym++)
	emcounts[i][sym] = 0.0;

      for (idx = 0; idx < nseq; idx++)
	{
	  if (aseqsT[i][idx] >= 0)
	    emcounts[i][aseqsT[i][idx]] += weights[idx];
	  else
	    gapcount[i] += weights[idx];
	}
    }
  gapcount[0]      = 0.0; 
  gapcount[alen+1] = 0.0;

  /* For each column, create an emission probability vector, 
   * and calculate score using emcounts.
   */
  for (i = 1; i <= alen; i++)
    {
      for (sym = 0; sym < ALPHASIZE; sym++)
	emvec[sym] = emcounts[i][sym];

      ProbifySingletEmission(emvec, uMATL_ST, prior);
      
      mscore[i] = dot_score(emcounts[i], emvec, ALPHASIZE);
      mscore[i] += (nseq - gapcount[i]) * (int) (INTPRECISION * log((float) ALPHASIZE));
    }

  Free2DArray(emcounts, alen+1);
  *ret_gapcount = gapcount;
  *ret_mscore = mscore;
}




/* Function: pair_emissioncost()
 * 
 * Purpose:  Count emission statistics for a given MATP_NODE-assigned
 *           i,j, and assign a cost.
 *           
 * Args:     coli     - ptr to column i from aseqsT transposed alignment, [0..nseq-1]
 *           colj     - ptr to column j from aseqsT transposed alignment, [0..nseq-1]
 *           weights  - weights on sequences (usually just 1.0 for each)
 *           nseq     - number of sequences
 *           prior    - prior probability distributions  
 * 
 * Return:   The cost of emitting the (i,j) column pair (integer).
 */
static int
pair_emissioncost(int            *coli,
		  int            *colj,
		  float          *weights,
		  int             nseq,
		  struct prior_s *prior)
{
  double matp_count[ALPHASIZE][ALPHASIZE];
  double matl_count[ALPHASIZE];
  double matr_count[ALPHASIZE];
  double matp_emit[ALPHASIZE][ALPHASIZE];
  double matl_emit[ALPHASIZE];
  double matr_emit[ALPHASIZE];
  int sc;
  int symi, symj;
  int idx;

  /* Zero the counter arrays
   */
  for (symi = 0; symi < ALPHASIZE; symi++)
    {
      matr_count[symi] = 0.0;
      matl_count[symi] = 0.0;
      for (symj = 0; symj < ALPHASIZE; symj++)
	matp_count[symi][symj] = 0.0;
    }
  
  /* Count pairs, singlets
   */
  for (idx = 0; idx < nseq; idx++)
    {
      if (coli[idx] == -1)
	{
	  if (colj[idx] == -1)  continue;
	  else                  matl_count[colj[idx]] += weights[idx];
	}
      else if (colj[idx] == -1)	matr_count[coli[idx]] += weights[idx];
      else                      matp_count[coli[idx]][colj[idx]] += weights[idx];
    }

  /* Create probability matrices
   */
  copy_singlet(matl_emit, matl_count);
  copy_singlet(matr_emit, matr_count);
  copy_pairwise(matp_emit, matp_count);
  ProbifySingletEmission(matl_emit, uMATL_ST, prior);
  ProbifySingletEmission(matr_emit, uMATR_ST, prior);
  ProbifyPairEmission(matp_emit, prior);

  /* convert probs to log odds 
   */
  for (symi = 0; symi < ALPHASIZE; symi++)
    {
      matl_emit[symi] *= (double) ALPHASIZE;
      matr_emit[symi] *= (double) ALPHASIZE;
      for (symj = 0; symj < ALPHASIZE; symj ++)
	matp_emit[symi][symj] *= (double) (ALPHASIZE*ALPHASIZE);
    }
  
  /* Score is the sum of dot products of counts and probability matrices
   */
  sc = dot_score((double *) matl_count, (double *) matl_emit, ALPHASIZE) + 
       dot_score((double *) matr_count, (double *) matr_emit, ALPHASIZE) +
       dot_score((double *) matp_count, (double *) matp_emit, ALPHASIZE*ALPHASIZE);
  return sc;
}





/* Function: frommatp_transtable()
 * 
 * Purpose:  Given an starting cell i,j and an ending cell
 *           i',j' (i < i' < j' < j), calculate the
 *           6x6 transition table between these pairs,
 *           assuming both (i,j) and (i',j') are assigned
 *           to MATP. Transitions to other types of nodes
 *           at i',j' (MATR, MATL, BIFURC) are easily derived from this one.
 *           
 * Args:     aseqsT     - transpose of alignment, [1..alen][0..nseq-1]
 *           weights    - weights on sequences (usually just 1.0 for each)
 *           nseq       - number of sequences                   
 *           i,j        - i,j
 *           i2,j2      - i',j'
 *           accum_insl - (0..nseq-1) array of number of inserted symbols
 *                        between i,i' for each sequence
 *           accum_insr - (0..nseq-1) array of number of inserted symbols
 *                        between j',j for each sequence
 *           trans[][]  - filled in: state transition matrix (counts)
 *           
 * Return:   (void)
 *           trans[][] is filled in.          
 */
static void
frommatp_transtable(int   **aseqsT,
		    float  *weights,
		    int     nseq,
		    int     i, 
		    int     j,
		    int     i2, 
		    int     j2,
		    int    *accum_insl,
		    int    *accum_insr,
		    double  trans[STATETYPES][STATETYPES])
{
  int fy, ty;			/* from state index, to state index */
  int idx;			/* counter for sequences */

  /* Zero the counter array
   */
  for (fy = 0; fy < STATETYPES; fy++)
    for (ty = 0; ty < STATETYPES; ty++)
      trans[fy][ty] = 0.0;

  /* For each sequence:
   *    assign fy, based on symbols vs. gaps at i,j
   *    assign ty, based on symbols vs. gaps at i',j'
   *    use accum_insr and accum_insl and bump counters appropriately.
   */
  for (idx = 0; idx < nseq; idx++)
    {
      fy = assign_cell(i,j,aseqsT[i][idx], aseqsT[j][idx]);
      ty = assign_cell(i2,j2,aseqsT[i2][idx], aseqsT[j2][idx]);

      if (accum_insl[idx] == 0 && accum_insr[idx] == 0)
	{
	  trans[fy][ty] += weights[idx];
	}
      else if (accum_insl[idx] > 0)
	{
	  trans[fy][INSL_ST]      += weights[idx];
	  trans[INSL_ST][INSL_ST] += (accum_insl[idx]-1) * weights[idx];
	  if (accum_insr[idx] > 0)
	    {
	      trans[INSL_ST][INSR_ST] += weights[idx];
	      trans[INSR_ST][INSR_ST] += (accum_insr[idx]-1) * weights[idx];
	      trans[INSR_ST][ty]      += weights[idx];
	    }
	  else
	    {
	      trans[INSL_ST][ty] += weights[idx];
	    }
	}
      else if (accum_insr[idx] > 0)
	{
	  trans[fy][INSR_ST]      += weights[idx];
	  trans[INSR_ST][INSR_ST] += (accum_insr[idx]-1) * weights[idx];
	  trans[INSR_ST][ty]      += weights[idx];
	}
    } /* end loop over all sequences */
}




/* Function: frommatl_transtable()
 * 
 * Purpose:  Given an starting cell i,j assigned to MATL_NODE 
 *           and an ending cell i',j (i < i' < j), calculate the
 *           6x6 transition table between these pairs,
 *           assuming (i',j) is assigned to MATP. 
 *           Transitions to other types of nodes
 *           at i',j (MATR, MATL, BIFURC) are easily derived from this one.
 *           
 * Args:     aseqsT     - transpose of alignment, [1..alen][0..nseq-1]
 *           weights    - weights on sequences (usually just 1.0 for each)
 *           nseq       - number of sequences                   
 *           i,j        - i,j
 *           i2         - i'
 *           accum_insl - (0..nseq-1) array of number of inserted symbols
 *                        between i,i' for each sequence
 *           trans[][]  - filled in: state transition matrix (counts)
 *           
 * Return:   (void)
 *           trans[][] is filled in.          
 */
static void
frommatl_transtable(int   **aseqsT,
		    float  *weights,
		    int     nseq,
		    int     i, 
		    int     j,
		    int     i2, 
		    int    *accum_insl,
		    double  trans[STATETYPES][STATETYPES])
{
  int fy, ty;			/* from state index, to state index */
  int idx;			/* counter for sequences */

  /* Zero the counter array
   */
  for (fy = 0; fy < STATETYPES; fy++)
    for (ty = 0; ty < STATETYPES; ty++)
      trans[fy][ty] = 0.0;


  /* For each sequence:
   *    assign fy, based on symbols vs. gaps at i
   *    assign ty, based on symbols vs. gaps at i',j
   *    use accum_insl and bump counters appropriately.
   */
  for (idx = 0; idx < nseq; idx++)
    {
      fy = (aseqsT[i][idx] == -1) ? DEL_ST : MATL_ST;
      ty = assign_cell(i2, j, aseqsT[i2][idx], aseqsT[j][idx]);

      if (accum_insl[idx] == 0)
	{
	  trans[fy][ty] += weights[idx];
	}
      else 
	{
	  trans[fy][INSL_ST]      += weights[idx];
	  trans[INSL_ST][INSL_ST] += (accum_insl[idx]-1) * weights[idx];
	  trans[INSL_ST][ty] += weights[idx];
	}
    } /* end loop over all sequences */
}




/* Function: frommatr_transtable()
 * 
 * Purpose:  Given an starting cell i,j assigned to MATR_NODE
 *           and an ending cell i,j' (i < j' < j), calculate the
 *           6x6 transition table between these pairs,
 *           assuming (i,j') is assigned to MATP. 
 *           Transitions to other types of nodes
 *           at i,j' (MATR, MATL, BIFURC) are easily derived from this one.
 *           
 * Args:     aseqsT     - transpose of alignment, [1..alen][0..nseq-1]
 *           weights    - weights on sequences (usually just 1.0 for each)
 *           nseq       - number of sequences                   
 *           i,j        - i,j
 *           j2         - j'
 *           accum_insr - (0..nseq-1) array of number of inserted symbols
 *                        between j',j for each sequence
 *           trans[][]  - filled in: state transition matrix (counts)
 *           
 * Return:   (void)
 *           trans[][] is filled in.          
 */
static void
frommatr_transtable(int   **aseqsT,
		    float  *weights,
		    int     nseq,
		    int     i, 
		    int     j,
		    int     j2,
		    int    *accum_insr,
		    double  trans[STATETYPES][STATETYPES])
{
  int fy, ty;			/* from state index, to state index */
  int idx;			/* counter for sequences */

  /* Zero the counter array
   */
  for (fy = 0; fy < STATETYPES; fy++)
    for (ty = 0; ty < STATETYPES; ty++)
      trans[fy][ty] = 0.0;

  /* For each sequence:
   *    assign fy, based on symbols vs. gaps at j
   *    assign ty, based on symbols vs. gaps at i,j'
   *    use accum_insr and bump counters appropriately.
   */
  for (idx = 0; idx < nseq; idx++)
    {
      fy = (aseqsT[j][idx] == -1) ? DEL_ST : MATR_ST;
      ty = assign_cell(i, j2, aseqsT[i][idx], aseqsT[j2][idx]);

      if (accum_insr[idx] == 0)
	{
	  trans[fy][ty] += weights[idx];
	}
      else 
	{
	  trans[fy][INSR_ST]      += weights[idx];
	  trans[INSR_ST][INSR_ST] += (accum_insr[idx]-1) * weights[idx];
	  trans[INSR_ST][ty]      += weights[idx];
	}
    } /* end loop over all sequences */
}




/* Function: frombeginr_transtable()
 * 
 * Purpose:  Given an starting cell i,j assigned to BEGINR_NODE
 *           and an ending cell i',j (i < i' < j), calculate the
 *           6x6 transition table between these pairs,
 *           assuming (i',j) is assigned to MATP. 
 *           Transitions to other types of nodes
 *           at i',j (MATR, MATL, BIFURC) are easily derived from this one.
 *           
 * Args:     aseqsT     - transpose of alignment, [1..alen][0..nseq-1]
 *           weights    - weights on sequences (usually just 1.0 for each)
 *           nseq       - number of sequences                   
 *           i,j        - i,j
 *           i2         - i'
 *           accum_insl - (0..nseq-1) array of number of inserted symbols
 *                        between j',j for each sequence
 *           trans[][]  - filled in: state transition matrix (counts)
 *           
 * Return:   (void)
 *           trans[][] is filled in.          
 */
static void
frombeginr_transtable(int   **aseqsT,
		      float  *weights,
		      int     nseq,
		      int     j,
		      int     i2,
		      int    *accum_insl,
		      double  trans[STATETYPES][STATETYPES])
{
  int fy, ty;			/* from state index, to state index */
  int idx;			/* counter for sequences */

  /* Zero the counter array
   */
  for (fy = 0; fy < STATETYPES; fy++)
    for (ty = 0; ty < STATETYPES; ty++)
      trans[fy][ty] = 0.0;

  /* For each sequence:
   *    assign ty, based on symbols vs. gaps at i',j
   *    use accum_insr and bump counters appropriately.
   */
  fy = BEGIN_ST;
  for (idx = 0; idx < nseq; idx++)
    {
      ty = assign_cell(i2, j, aseqsT[i2][idx], aseqsT[j][idx]);

      if (accum_insl[idx] == 0)
	{
	  trans[fy][ty] += weights[idx];
	}
      else 
	{
	  trans[fy][INSL_ST]      += weights[idx];
	  trans[INSL_ST][INSL_ST] += (accum_insl[idx]-1) * weights[idx];
	  trans[INSL_ST][ty]      += weights[idx];
	}
    } /* end loop over all sequences */
}


/* Function: frombeginl_transtable()
 * 
 * Purpose:  Given an starting cell i,j assigned to BEGINR_NODE (which must
 *           connect to i,j itself),
 *           calculate the 6x6 transition table between the BEGINR
 *           and the other states in the cell, assuming i,j is
 *           assigned to MATP.
 *           Transitions to other types of nodes
 *           at i,j (MATR, MATL, BIFURC) are easily derived from this one.
 *           
 * Args:     aseqsT     - transpose of alignment, [1..alen][0..nseq-1]
 *           weights    - weights on sequences (usually just 1.0 for each)
 *           nseq       - number of sequences                   
 *           i,j        - i,j
 *           trans[][]  - filled in: state transition matrix (counts)
 *           
 * Return:   (void)
 *           trans[][] is filled in.          
 */
static void
frombeginl_transtable(int   **aseqsT,
		      float  *weights,
		      int     nseq,
		      int     i, 
		      int     j,
		      double  trans[STATETYPES][STATETYPES])
{
  int fy, ty;			/* from state index, to state index */
  int idx;			/* counter for sequences */

  /* Zero the counter array
   */
  for (fy = 0; fy < STATETYPES; fy++)
    for (ty = 0; ty < STATETYPES; ty++)
      trans[fy][ty] = 0.0;

  /* For each sequence:
   *    assign ty, based on symbols vs. gaps at i',j
   *    use accum_insr and bump counters appropriately.
   */
  fy = BEGIN_ST;
  for (idx = 0; idx < nseq; idx++)
    {
      ty = assign_cell(i, j, aseqsT[i][idx], aseqsT[j][idx]);
      trans[fy][ty] += weights[idx];
    }
}



/* Function: fromroot_transtable()
 * 
 * Purpose:  Given an ending cell i2,j2, calculate the
 *           6x6 transition table between 1,alen,ROOT and this pair,
 *           assuming (i',j') is assigned
 *           to MATP. Transitions to other types of nodes
 *           at i',j' (MATR, MATL, BIFURC) are easily derived from this one.
 *           
 * Args:     aseqsT     - transpose of alignment, [1..alen][0..nseq-1]
 *           weights    - weights on sequences (usually just 1.0 for each)
 *           nseq       - number of sequences                   
 *           i2,j2      - i',j'
 *           accum_insl - (0..nseq-1) array of number of inserted symbols
 *                        between i,i' for each sequence
 *           accum_insr - (0..nseq-1) array of number of inserted symbols
 *                        between j',j for each sequence
 *           trans[][]  - filled in: state transition matrix (counts)
 *           
 * Return:   (void)
 *           trans[][] is filled in.          
 */
static void
fromroot_transtable(int   **aseqsT,
		    float  *weights,
		    int     nseq,
		    int     i2, 
		    int     j2,
		    int    *accum_insl,
		    int    *accum_insr,
		    double  trans[STATETYPES][STATETYPES])
{
  int fy, ty;			/* from state index, to state index */
  int idx;			/* counter for sequences */

  /* Zero the counter array
   */
  for (fy = 0; fy < STATETYPES; fy++)
    for (ty = 0; ty < STATETYPES; ty++)
      trans[fy][ty] = 0.0;

  /* For each sequence:
   *    assign ty, based on symbols vs. gaps at i',j'
   *    use accum_insr and accum_insl and bump counters appropriately.
   */
  for (idx = 0; idx < nseq; idx++)
    {
      ty = assign_cell(i2,j2,aseqsT[i2][idx], aseqsT[j2][idx]);

      if (accum_insl[idx] == 0 && accum_insr[idx] == 0)
	{
	  trans[BEGIN_ST][ty] += weights[idx];
	}
      else if (accum_insl[idx] > 0)
	{
	  trans[BEGIN_ST][INSL_ST] += weights[idx];
	  trans[INSL_ST][INSL_ST]  += (accum_insl[idx]-1) * weights[idx];
	  if (accum_insr[idx] > 0)
	    {
	      trans[INSL_ST][INSR_ST] += weights[idx];
	      trans[INSR_ST][INSR_ST] += (accum_insr[idx]-1) * weights[idx];
	      trans[INSR_ST][ty]      += weights[idx];
	    }
	  else
	    {
	      trans[INSL_ST][ty] += weights[idx];
	    }
	}
      else if (accum_insr[idx] > 0)
	{
	  trans[BEGIN_ST][INSR_ST] += weights[idx];
	  trans[INSR_ST][INSR_ST]  += (accum_insr[idx]-1) * weights[idx];
	  trans[INSR_ST][ty]       += weights[idx];
	}
    } /* end loop over all sequences */
}



/* Function: assign_cell()
 * 
 * Purpose:  Given that we assign a cell i,j to MATP for the
 *           whole alignment, return the actual assignment
 *           for the given sequence seq. This will be MATP
 *           if seq has symbols in both columns i,j; 
 *           MATL if j is a gap; MATR if i is a gap;
 *           DEL if both i,j are gaps.
 *
 *           If i == j, assign_cell returns MATP or DEL even though
 *           MATP is not a valid assignment; the to_matr or to_matl functions
 *           straighten this out later.
 *           
 * Args:     i,j : coordinates (columns in alignment)
 *           symi: symbol in position i
 *           symj: symbol in position j
 *           
 * Return:   DEL_ST, MATP_ST, MATR_ST, or MATL_ST
 */
static int
assign_cell(int i, int j, int symi, int symj)
{
  if (i > j)          return END_ST;
  else if (symi >= 0) 
    {
      if (symj >= 0)  return MATP_ST;
      else            return MATL_ST;
    }
  else if (symj >= 0) return MATR_ST;
  else                return DEL_ST;
}




/* Function: to_matp_transtable()
 * 
 * Purpose:  Given a master state transition tables containing
 *           counts, create a table specific for transitions
 *           to an (i',j') assigned to MATP_NODE.
 *           
 * Args:     master_table[][]  - the master table containing counts
 *           trans[][]         - the new table, returned containing counts
 * 
 * Return:   (void)
 *           trans[][] is filled.
 */                              
static void
to_matp_transtable(double master_table[STATETYPES][STATETYPES],
		   double trans[STATETYPES][STATETYPES])
{
  int fy,ty;			/* indices for from state, to state */

  /* Just copy the table.
   */
  for (fy = 0; fy < STATETYPES; fy++)
    for (ty = 0; ty < STATETYPES; ty++)
      trans[fy][ty] = master_table[fy][ty];
}


/* Function: to_matr_transtable()
 * 
 * Purpose:  Given a master state transition tables containing
 *           counts, create a table specific for transitions
 *           to an (i',j') assigned to MATR_NODE.
 *           
 * Args:     master_table[][]  - the master table containing counts
 *           trans[][]         - the new table, returned containing counts
 * 
 * Return:   (void)
 *           trans[][] is filled.
 */                              
static void
to_matr_transtable(double master_table[STATETYPES][STATETYPES],
		   double trans[STATETYPES][STATETYPES])
{
  int fy;			/* indices for from state, to state */

  for (fy = 0; fy < STATETYPES; fy++)
    {
      trans[fy][DEL_ST]  = master_table[fy][DEL_ST] + master_table[fy][MATL_ST];
      trans[fy][MATP_ST] = 0.0;
      trans[fy][MATL_ST] = 0.0;
      trans[fy][MATR_ST] = master_table[fy][MATR_ST] + master_table[fy][MATP_ST];
      trans[fy][INSL_ST] = master_table[fy][INSL_ST];
      trans[fy][INSR_ST] = master_table[fy][INSR_ST];
    }
}


/* Function: to_matl_transtable()
 * 
 * Purpose:  Given a master state transition tables containing
 *           counts, create a table specific for transitions
 *           to an (i',j') assigned to MATL_NODE.
 *           
 * Args:     master_table[][]  - the master table containing counts
 *           trans[][]         - the new table, returned containing counts
 * 
 * Return:   (void)
 *           trans[][] is filled.
 */                              
static void
to_matl_transtable(double master_table[STATETYPES][STATETYPES],
		   double trans[STATETYPES][STATETYPES])
{
  int fy;			/* indices for from state, to state */

  for (fy = 0; fy < STATETYPES; fy++)
    {
      trans[fy][DEL_ST]  = master_table[fy][DEL_ST] + master_table[fy][MATR_ST];
      trans[fy][MATP_ST] = 0.0;
      trans[fy][MATL_ST] = master_table[fy][MATL_ST] + master_table[fy][MATP_ST];
      trans[fy][MATR_ST] = 0.0;
      trans[fy][INSL_ST] = master_table[fy][INSL_ST];
      trans[fy][INSR_ST] = master_table[fy][INSR_ST];
    }
}


/* Function: to_bifurc_transtable()
 * 
 * Purpose:  Given a master state transition tables containing
 *           counts, create a table specific for transitions
 *           to an (i',j') assigned to BIFURC_NODE.
 *           
 * Args:     master_table[][]  - the master table containing counts
 *           trans[][]         - the new table, returned containing counts
 * 
 * Return:   (void)
 *           trans[][] is filled.
 */                              
static void
to_bifurc_transtable(double master_table[STATETYPES][STATETYPES],
		     double trans[STATETYPES][STATETYPES])
{
  int fy;			/* indices for from state, to state */

  for (fy = 0; fy < STATETYPES; fy++)
    {
      trans[fy][BIFURC_ST] = master_table[fy][DEL_ST] + master_table[fy][MATR_ST] +
	                     master_table[fy][MATL_ST] + master_table[fy][MATP_ST];
      trans[fy][MATP_ST] = 0.0;
      trans[fy][MATL_ST] = 0.0;
      trans[fy][MATR_ST] = 0.0;
      trans[fy][INSL_ST] = master_table[fy][INSL_ST];
      trans[fy][INSR_ST] = master_table[fy][INSR_ST];
    }
}





/* Function: dot_score()
 * 
 * Purpose:  Calculate the dot product of counts by probabilities
 *           for s set of state transitions: i.e. a total score
 *           
 * Args:     cmx:  count vector
 *           pmx:  probability vector
 *           
 * Return:   score, the dot product, as an integer
 */          
static int
dot_score(double *cvec, 
	  double *pvec,
	  int     veclen)
{
  int i;
  double score = 0.0;

  for (i = 0; i < veclen; i++, cvec++, pvec++)
    if (*pvec > 0.0) score +=  *cvec * log(*pvec);
  return ((int) INTPRECISION * score);
}
      



#ifdef DEBUG

static void
print_mmx(struct maxmx_s **mmx, int alen)
{
  int i,j,y;

  for (j = 0; j <= alen; j++)
    {
      for (y = 0; y < NODETYPES-1; y++)
	{
	  for (i = 0; i <= j+1; i++)
	    printf("%5d ", mmx[j][i].sc[y] * 10 / (int) INTPRECISION);
	  printf("\n");
	}
      puts("");
    }
}

/* Function: print_assignments()
 * 
 * Purpose:  Using a master max likelihood tree constructed from an
 *           alignment, print a line for each column of the alignment:
 *           (column)  (fractional occupancy)  (assignment)
 *           
 *           The master tree, mtr, is a traceback tree with some of
 *           the fields misused. emitl, emitr contain 0..alen-1 column
 *           coordinates (even if only one of them is emitted by
 *           this tree node). nodeidx is unused. type contains
 *           a node type, not a state type. 
 *           
 *           Code partly borrowed from Trace2ali().
 */
static void
print_assignments(struct trace_s *mtr, int nseq, int alen, double *gapcount)
{
  struct align_s      *ali;     /* linear list of alignment        */
  struct t2ali_s      *stack;	/* stack used to traverse the traceback tr */
  struct trace_s      *currtr;
  struct align_s      *currali;
  struct align_s      *newafter;
  struct align_s      *oldafter; 
  int                  currpos;

  /* First, we generate a linear linked align_s list from the tree.
   * The align_s structure fields are used as follows:
   *   pos     = 0..alen-1 position in alignment
   *   sym     = unused
   *   nodeidx = unused
   *   type    = node type (e.g. MATP_NODE)
   */

  /* Initialize the linked list for the alignment of sequence to model
   */
  ali = Init_align();

  /* Initialize the pushdown stack for traversal of the traceback
   */
  stack = Init_t2ali();
  Push_t2ali(stack, mtr->nxtl, ali);

  while (Pop_t2ali(stack, &currtr, &oldafter))
    {
      if (currtr->nxtl == NULL) continue;	/* ignore END nodes */

      if (currtr->nxtr != NULL)	/* BIFURC node */
	{
				/* deal with right branch; insert a dummy */
	  newafter = Insafter_align(-1, '-', ' ', 0, currtr->type, oldafter);
	  Push_t2ali(stack, currtr->nxtr, newafter);
				/* deal with left branch */
	  Push_t2ali(stack, currtr->nxtl, oldafter);
	}
      else{
	  switch (currtr->type) {
	  case BEGINL_NODE:
	  case BEGINR_NODE:
	  case ROOT_NODE:
	    Push_t2ali(stack, currtr->nxtl, oldafter);
	    break;

	  case MATP_NODE:
	    (void) Insafter_align(currtr->emitr, ' ', ' ', 0, currtr->type, oldafter);
	    newafter = Insafter_align(currtr->emitl, ' ', ' ', 0, currtr->type, oldafter);
	    Push_t2ali(stack, currtr->nxtl, newafter);
	    break;

	  case MATL_NODE:
	    newafter = Insafter_align(currtr->emitl, ' ', ' ', 0, currtr->type, oldafter);
	    Push_t2ali(stack, currtr->nxtl, newafter);
	    break;

	  case MATR_NODE:
	    (void) Insafter_align(currtr->emitr, ' ', ' ', 0, currtr->type, oldafter);
	    Push_t2ali(stack, currtr->nxtl, oldafter);
	    break;
	    
	  case END_NODE:
	    break;

	  default:
	    Die("no such node type %d", currtr->type);
	  }
	}
    }
  Free_t2ali(stack);

  currpos = 0;
  for (currali = ali->nxt; currali != NULL; currali = currali->nxt)
    {
      currpos++;
      while (currpos < currali->pos)
	{
	  printf("%4d %.3f INS\n",
		 currpos,
		 ((double) nseq - gapcount[currpos]) / (double) nseq);
	  currpos++;
	}


      switch (currali->type) {
      case BEGINL_NODE:
      case BEGINR_NODE:
      case ROOT_NODE:
	break;
	
      case MATL_NODE:
      case MATR_NODE:
	printf("%4d %.3f MAT singlet\n",
	       currali->pos,
	       ((double) nseq - gapcount[currali->pos]) / (double) nseq);
	break;

      case MATP_NODE:
	printf("%4d %.3f MAT pairwise\n",
	       currali->pos,
	       ((double) nseq - gapcount[currali->pos]) / (double) nseq);
	break;
      }
    }

  while (currpos <= alen)
    {
      printf("%4d %.3f INS\n",
	     currpos,
	     ((double) nseq - gapcount[currpos]) / (double) nseq);
      currpos++;
    }

  Free_align(ali);
}

#endif /* DEBUG */



