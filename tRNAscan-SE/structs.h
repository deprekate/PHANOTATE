#ifndef STRUCTSH_INCLUDED
#define STRUCTSH_INCLUDED
/* structs.h - declarations of data structures
 * SRE, Tue Aug 31 15:17:12 1993
 */
#include "squid.h"

/* Alphabet information.
 * The package is designed to be configurable for protein analysis
 * just by changing these define's. Dunno if it would be *useful*
 * to apply it to protein work -- but the possibility's there.
 */
#define ALPHATYPE   kRNA
#define ALPHASIZE   4
extern char *ALPHABET;		/* defined at top of misc.c */

/* 
 * Node types.
 * These are used for clarity in the code, not to make it easy
 * to change them: the program makes assumptions about the order
 * they come in, so DON'T CHANGE THESE.
 *
 * Specifically, the following assumptions are made:
 *  - that they come in *exactly* this order (static state transition array in prior.h)
 *  - that BIFURC, MATP, MATL, MATR are all less than 4 
 *       (second index, prior.h; also some loops in maxmodelmaker.c)
 *  - ROOT is last (indexing of scores of mmx model construction matrix, maxmodelmaker.c)
 */
#define BIFURC_NODE 0
#define MATP_NODE   1
#define MATL_NODE   2
#define MATR_NODE   3
#define BEGINL_NODE 4		
#define BEGINR_NODE 5
#define ROOT_NODE   6		
#define NODETYPES   7		/* number of different node types  */

#define END_NODE    BIFURC_NODE

/* 
 * State types.
 * These are used for clarity in the code, not to make it easy
 * to change them: the program makes assumptions about the order
 * they come in, so DON'T CHANGE THESE.
 */
#define DEL_ST         0
#define MATP_ST        1
#define MATL_ST        2
#define MATR_ST        3
#define INSL_ST        4
#define INSR_ST        5
#define STATETYPES     6	/* MATP nodes contain 6 states */

#define BEGIN_ST       DEL_ST
#define BIFURC_ST      DEL_ST
#define END_ST         DEL_ST

/* Unique identifiers for state types, used as flags not indexes
 * in the alignment algorithms.
 */
#define uDEL_ST     (1<<0)
#define uMATP_ST    (1<<1)
#define uMATL_ST    (1<<2)
#define uMATR_ST    (1<<3)
#define uINSL_ST    (1<<4) 
#define uINSR_ST    (1<<5)
#define uBEGIN_ST   (1<<6)
#define uEND_ST     (1<<7)
#define uBIFURC_ST  (1<<8)


/* Structure: node_s
 * 
 * Purpose:   Contains all the information necessary to describe a node.
 */
struct node_s {
  int  type;

  double tmx[STATETYPES][STATETYPES];	 /* up to 49 transition probs */
  double mp_emit[ALPHASIZE][ALPHASIZE];  /* 4x4 MATP emission probs   */
  double il_emit[ALPHASIZE];             /* 4 INSL emission probs     */
  double ir_emit[ALPHASIZE];             /* 4 INSR emission probs     */
  double ml_emit[ALPHASIZE];             /* 4 MATL emission probs     */
  double mr_emit[ALPHASIZE];             /* 4 MATR emission probs     */

  int nxt;			/* connection to left child */
  int nxt2;                     /* connection to right child */
};



/* Structure: cm_s
 * 
 * Purpose:   A covariance model.
 */
struct cm_s {
  int nodes;			/* number of nodes  */
  struct node_s *nd;            /* array of nodes 0..nodes-1 */
};


/* Structure: istate_s
 * 
 * In the alignment algorithms, a CM is converted to an array of states, 
 * each represented by one of these structures. Each state contains 
 * probability info as integers instead of floating point.
 * 
 * The order of the state transition vector is different than in
 * the CM. INSL and INSR are first: INSL, INSR, DEL, MATP, MATL, MATR.
 */ 
struct istate_s {
  int   nodeidx;		/* index of node this state belongs to        */
  int   statetype;		/* unique id for type of this state (uMATP_ST, etc.) */
  int   offset;			/* offset in state array to first INS state   */
  int   connectnum;		/* number of elements in tmx                  */
  int   tmx[STATETYPES];	/* rearranged transition vector, int log-odds */
  int   emit[ALPHASIZE*ALPHASIZE]; /* int lod emission vector (4 or 16) or NULL  */
};

struct pstate_s {
  int   nodeidx;		/* index of node this state belongs to        */
  int   statetype;		/* unique id for type of this state (uMATP_ST, etc.) */
  int   offset;			/* offset in state array to first INS state   */
  int   connectnum;		/* number of elements in tmx                  */
  int   bifr;			/* (uBIF_ST only) index of right connection   */
  double tmx[STATETYPES];	/* rearranged transition vector               */
  double emit[ALPHASIZE*ALPHASIZE]; /* emission vector (4 or 16) or NULL       */
};



/* Structure: prior_s
 * 
 * Purpose:   Contains the prior probability distributions for
 *            state transitions and symbol emissions, as well
 *            as the alpha "confidence" values applied during
 *            regularization, and alphabet information.
 */
struct prior_s {
  double tprior[7][4][STATETYPES][STATETYPES];	/* state transitions  */
  double matp_prior[ALPHASIZE][ALPHASIZE];      /* MATP_ST emissions  */
  double matl_prior[ALPHASIZE];                 /* MATL_ST emissions  */
  double matr_prior[ALPHASIZE];                 /* MATR_ST emissions  */
  double insl_prior[ALPHASIZE];                 /* INSL_ST emissions  */
  double insr_prior[ALPHASIZE];                 /* INSR_ST emissions  */

  double talpha[STATETYPES];	/* alpha's for state transitions      */
  double emalpha[STATETYPES];   /* alpha's for symbol emissions       */

  double rfreq[ALPHASIZE];	/* background symbol freqs for random model */
};

  

/* Structure: trace_s
 * 
 * Binary tree structure for storing a traceback of an alignment;
 * also used for tracebacks of model constructions.
 */
struct trace_s {
  int emitl;			/* i position (1..N) or 0 if nothing */
  int emitr;			/* j position (1..N) or 0 if nothing */
  
  int nodeidx;			/* index of node responsible for this alignment */
  int type;			/* type of substate (uMATP_ST, etc.) used (unique) */
  
  struct trace_s *nxtl;		/* ptr to left (or only) branch, or NULL for end */
  struct trace_s *nxtr;		/* ptr to right branch, BIFURC only, else NULL */
  struct trace_s *prv;          /* ptr to parent                               */
};


/* Structure: trmem_s
 * 
 * It's expensive in malloc()'s to build trace trees. This structure
 * allows trace.c to cut down malloc overhead, by keeping a pool
 * of trace_s structures.
 */
struct trmem_s {
  int next;			/* index of next trace_s to use in pool        */
  int num;			/* how many trace_s total in pool              */
  struct trace_s *pool;         /* alloced array of trace_s structs            */
  struct tracestack_s *used;    /* old (fully used) pools, waiting to be freed */
};
#define TMEM_BLOCK 256		/* how many trace_s to alloc per malloc() call */


/* Structure: tracestack_s
 *
 * Formerly a pushdown stack used for traversing a binary tree of trace_s structures.
 * Reimplemented as an array for malloc efficiency.
 */
struct tracestack_s {
  int next;			/* index of next trace_s pointer to use */
  int num;			/* number of trace_s pointers alloc'ed  */
  struct trace_s **list;        /* array of trace_s pointers            */
};
#define TSTACK_BLOCK 64

/* A struct align_s implements a linked list describing the alignment
 * of a model to a sequence. Note that this is the inverse of what
 * trace_s trees are for; align_s is a linear representation of
 * the alignment (from the sequence's point of view, if you will)
 */
struct align_s {
  int              pos;		/* pos in seq emitted (0..N-1; -1 if none)       */
  char             sym;		/* symbol emitted (ACGU, . if none)              */
  char             ss;          /* secondary structure character, <>.            */
  int              nodeidx;     /* index of model state aligned to this position */
  int              type;	/* type of substate reponsible for this emission (unique) */
  struct align_s  *nxt;
};

/* A struct m2ali_s implements a pushdown stack used for traversing
 * a model and producing an align_s alignment list. 
 */
struct m2ali_s {
  int              nodeidx;	/* index of position in model (0..M) */
  int              type;	/* subtype of position in model      */
  struct align_s  *after;       /* position in align_s list          */
  struct m2ali_s  *nxt;         
};


/* A struct t2ali_s implements a pushdown stack for traversing a
 * traceback tree and producing an align_s alignment list.
 */
struct t2ali_s {
  struct trace_s *tracenode;
  struct align_s *after;
  struct t2ali_s *nxt;
};

/* some stuff used when we store sums of log scores as integers,
 * for speed and precision
 */
#define INTPRECISION 1000.0	/* pick up three decimal places in our ints */
#define NEGINFINITY -999999     /* -999.999 is small enough for -Inf        */
#define POSINFINITY  999999     /* +999.999 is large enough for +Inf        */
#define ILOG2(a) (((a) > 0.0) ? (log(a) / 0.69314718 * INTPRECISION) : NEGINFINITY)

#endif /* STRUCTSH_INCLUDED */
