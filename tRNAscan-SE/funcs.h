/* funcs.h
 * Declarations and prototypes of functions
 * SRE, Fri Sep  3 09:19:50 1993
 * 
 */

#include <stdio.h>
#include "structs.h"

/* 
 * from align.c
 */
extern int Trace2ali(char *seq, struct trace_s *tr, int watsoncrick, struct align_s **ret_ali);
extern int Traces2Alignment(char **rseqs, SQINFO *sqinfo, struct trace_s **tr, int nseq,
			    struct cm_s *cm, int watsoncrick, char ***ret_aseqs,
			    AINFO *ainfo);


/* from dbviterbi.c
 */
extern int ViterbiScan(struct istate_s *icm, int statenum, char *seq, int window, 
		       double thresh, int (*gotone_f)(int, int, double));


/* from debug.c
 */
extern char *UstatetypeName(int ustatetype);
extern char *StatetypeName(int statetype);
extern char *NodetypeName(int nodetype);
extern void  PrintViterbiAMX(FILE *fp, struct istate_s *icm, int statenum, 
			     char *seq, int N, int ***amx);
extern void  PrintTrace(FILE *fp, struct trace_s *tr);
extern void  PrintAli(FILE *fp, struct align_s *ali);
extern void  PrintICM(FILE *fp, struct cm_s *cm, struct istate_s *icm, int nstates);



/* 
 * from emit.c
 * 
 * Generate sequences from a CVHMM
 */
extern int EmitSequence(struct cm_s *cm, int watsoncrick, struct align_s **ret_ali, 
			char **ret_khseq, char **ret_seq);
extern int EmitBestSequence(struct cm_s *cm, int watsoncrick, struct align_s **ret_ali, 
			    char **ret_khseq, char **ret_seq);

/* 
 * from fastmodelmaker.c
 */
extern int Fastmodelmaker(char **aseqs, AINFO *ainfo, int nseq, struct prior_s *prior, 
			  double gapthresh, double *ret_secinfo, 
			  struct cm_s **ret_cm, struct trace_s **ret_mtr);


/* from fast-dbviterbi.c
 */
extern int FastViterbiScan(struct istate_s *icm, int statenum, int *minb, int *maxb,
			   char *seq, int window, double thresh, 
			   int (*gotone_f)(int, int, double));


/* 
 * from konings.c
 * 
 * ASCII representation of structures and structural alignments
 */
extern int  Align2kh(struct align_s  *ali, char **ret_aseq, char **ret_khseq);
extern int  PrintAliLandscape(FILE *fp, struct cm_s *cm, struct align_s *ali);
extern void Trace2KHS(struct trace_s *tr, char *seq, int rlen, int watsoncrick, 
		      char **ret_ss);
extern int  KHS2ct(char *ss, int len, int allow_pknots, int **ret_ct);
extern int  IsRNAComplement(char sym1, char sym2, int allow_gu);

/* 
 * from lengthdist.c
 */
extern void LengthDistribution(struct pstate_s *pcm, int statenum, int N, double ***ret_lmx);
extern void LengthBounds(double **lmx, int statenum, int N, double epsilon, 
			 int **ret_min, int **ret_max);


/* 
 * from maxmodelmaker.c
 */
extern int Maxmodelmaker(char **aseqs, AINFO *ainfo, int nseq, double gapthresh, 
			 struct prior_s *prior, double *ret_ssinfo,
			 struct cm_s **ret_cm, struct trace_s **ret_mtr);

/* 
 * from misc.c
 */
extern int SymbolIndex(char  sym);
extern int PrepareSequence(char *seq);

/* 
 * from model.c
 */
extern struct cm_s *AllocCM(int nodes);
extern void FreeCM(struct cm_s *cm);
extern void NormalizeCM(struct cm_s *cm);
extern int  VerifyCM(struct cm_s *cm);
extern int  RearrangeCM(struct cm_s *cm, double *rfreq, struct istate_s **ret_icm, int *ret_statenum);
extern int  MakePCM(struct cm_s *cm, struct pstate_s **ret_pcm, int *ret_statenum);
extern void NormalizePCM(struct pstate_s *pcm, int M);

/* from modelmaking.c
 */
extern void NumberMasterTrace(struct trace_s *mtr, int *ret_nodes);
extern void TopofyNewCM(struct cm_s *cm, struct trace_s *mtr);
extern void Transmogrify(struct trace_s *mtr, char *aseq, struct trace_s **ret_tr, struct trmem_s **ret_pool);
extern void EasyModelmaker(char **aseq, AINFO *ainfo, int nseq, struct prior_s *prior, 
			   double gapthresh, int use_rf, 
			   struct cm_s **ret_cm, struct trace_s **ret_mtr);

/*
 * from prior.c
 */
extern int  DefaultPrior(struct prior_s **ret_prior);
extern int  WritePrior(FILE *fp, struct prior_s *prior);
extern int  ReadPrior(FILE *fp, struct prior_s **ret_prior);
extern void NormalizePrior(struct prior_s *prior);

/* 
 * from probify.c
 */
extern void ProbifyCM(struct cm_s *cm, struct prior_s *prior);
extern void ProbifyTransitionMatrix(double tmx[STATETYPES][STATETYPES],int from_node, 
				    int to_node, struct prior_s *prior);
extern void ProbifySingletEmission(double emvec[ALPHASIZE], int statetype, struct prior_s *prior);
extern void ProbifyPairEmission(double emx[ALPHASIZE][ALPHASIZE], struct prior_s *prior);

/* 
 * from save.c
 */
extern int WriteCM(FILE *fp, struct cm_s *cm);
extern int WriteBinaryCM(FILE  *fp, struct cm_s *cm);
extern int ReadCM(char *filename, struct cm_s **ret_cm);

/* from scorestack.c
 */
extern int ReportScanHit(int left, int right, double score, int (*print_hit)(int,int,double));

/* from smallviterbi.c
 */
extern int SmallViterbiAlign(struct istate_s *icm, int statenum, char *seq, 
			     double *ret_score, struct trace_s **ret_trace);

/* 
 * from structs.c
 * 
 * Implementation of data structures: trees, stacks, and linked lists
 */
extern int StatetypeIndex(int type);
extern int UniqueStatetype(int nodetype, int stidx);
extern struct m2ali_s *Init_m2ali(void);
extern void Push_m2ali(struct m2ali_s *stack, int nodeidx, int type, struct align_s *after);
extern int  Pop_m2ali(struct m2ali_s *stack, int *ret_nodeidx, int *ret_type, 
		      struct align_s **ret_after);
extern void Free_m2ali( struct m2ali_s *stack );

extern struct t2ali_s *Init_t2ali(void);
extern void Push_t2ali(struct t2ali_s *stack, struct trace_s *tracenode, struct align_s *after);
extern int  Pop_t2ali(struct t2ali_s  *stack, struct trace_s **ret_tracenode, 
		      struct align_s **ret_after);
extern void Free_t2ali( struct t2ali_s *stack );

extern struct align_s *Init_align(void);
extern struct align_s *Insafter_align(int pos, char sym, char ss, int nodeidx, int type, 
				      struct align_s *after);
extern void Delafter_align(struct align_s *after);
extern void Free_align(struct align_s *head);
extern struct intstack_s *InitIntStack(void);
extern void PushIntStack(struct intstack_s *stack, int data);
extern int  PopIntStack(struct intstack_s  *stack, int *ret_data);
extern int  FreeIntStack( struct intstack_s *stack );


/* 
 * from trace.c
 */
extern void            InitTrace(struct trace_s **ret_new, struct trmem_s **ret_pool);
extern struct trace_s *AttachTrace(struct trace_s *parent, struct trmem_s *pool, 
				   int emitl, int emitr,
				   int nodeidx, int type);
extern void            FreeTrace(struct trace_s *tr, struct trmem_s *pool);
extern void            DeleteTracenode(struct trace_s *oldtr, struct trmem_s *pool);
extern void            InitTracepool(struct trmem_s **ret_pool);
extern struct trace_s *PopTracepool(struct trmem_s *pool);
extern void            FreeTracepool(struct trmem_s *pool);
extern struct tracestack_s *InitTracestack(void);
extern void                 PushTracestack(struct tracestack_s *stack, struct trace_s *node);
extern struct trace_s      *PopTracestack(struct tracestack_s *stack);
extern void                 FreeTracestack(struct tracestack_s *stack);
extern int                  TraceCount(struct cm_s *cm, char *seq, double weight, struct trace_s *tr);
extern int TraceCountPrior(struct cm_s *cm, struct prior_s *prior, char *seq,
			   double weight, struct trace_s   *tr);

extern double TraceScore(struct cm_s *cm, char *seq, struct trace_s *tr);

/* 
 * from viterbi.c
 */
extern int ViterbiAlign(struct istate_s *cm, int statenum, char *seq, 
			double *ret_score, struct trace_s **ret_trace);




