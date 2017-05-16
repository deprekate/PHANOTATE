/* trace.c
 * cove 1.0: Mon May 17 09:38:14 1993
 * moved to cove 2.0, Mon Sep  6 13:34:55 1993
 * 
 * Unlike a traceback of a normal HMM alignment, which is linear,
 * the traceback of a covariance HMM is a tree structure. Here
 * we provide support for the traceback data structures: the
 * tree itself, and a pushdown stack used for traversing the
 * tree.
 * 
 * The trace tree structure has a dummy node at its beginning,
 * and dummy end nodes at the termination of each branch. Non-BIFURC
 * states have a NULL right branch. 
 * 
 * The pushdown stack structure has a dummy begin node, and the
 * end is signified by a final NULL ptr.
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "structs.h"		/* struct trace_s and struct tracestack_s */
#include "funcs.h"
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif
#ifdef DEBUG
#include <assert.h>
#endif

/* Function: InitTrace()
 * 
 * Purpose:  Initialize a traceback tree structure.
 *           ret_tmem may be passed as NULL for default behavior;
 *           if ret_tmem is passed, enables optimized memory
 *           behavior for traces.
 *
 * Return:   ptr to the new tree.
 */          
void
InitTrace(struct trace_s **ret_new, struct trmem_s **ret_tmem)
{
  struct trace_s *new;
  struct trace_s *end;
  struct trmem_s *pool;

  if (ret_tmem != NULL)
    {
      InitTracepool(&pool);
      new = PopTracepool(pool);
    }
  else if ((new = (struct trace_s *) malloc (sizeof(struct trace_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  
  new->emitl = new->emitr = -1;
  new->nodeidx = 0;
  new->type    = uBEGIN_ST;
  new->nxtr    = NULL;
  new->prv     = NULL;
  
  if (ret_tmem != NULL) 
    end = PopTracepool(pool);
  else if ((end = (struct trace_s *) malloc (sizeof(struct trace_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  end->type = uEND_ST;
  end->emitl = end->emitr = end->nodeidx = -1;
  end->nxtr = end->nxtl = NULL;

  end->prv  = new;
  new->nxtl = end;

  *ret_new = new;
  if (ret_tmem != NULL) *ret_tmem = pool;
}

/* Function: AttachTrace()
 * 
 * Purpose:  attach a new node to a tracetree node.
 *           There are dummy END nodes. 
 *           
 *           Because of the mechanics of tracebacks through a Viterbi matrix,
 *           we have to make sure that BIFURC children are attached
 *           right first, left second.
 *
 *           trmem_s may be NULL (default behavior) or an active 
 *           trace pool (optimized memory behavior)
 *           
 * Returns:  ptr to the new node, or NULL on failure.
 */          
struct trace_s *
AttachTrace(struct trace_s *parent,
	    struct trmem_s *pool,
	    int             emitl,
	    int             emitr,
	    int             nodeidx,
	    int             type)
{
  struct trace_s *new;
  struct trace_s *end;

  if (parent->nxtr != NULL)
    Die("That trace node is already full, fool.");

  /* If left branch is already connected to something, swap it over to the
   * right (thus enforcing the necessary rule that BIFURCS attach to the right
   * branch first), and attach a new dummy end to the left branch. 
   */
  if (parent->nxtl->nxtl != NULL)
    {
      parent->nxtr = parent->nxtl;

      if (pool != NULL)
	end = PopTracepool(pool);
      else if ((end = (struct trace_s *) malloc (sizeof(struct trace_s))) == NULL)
	Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
      end->type = uEND_ST;
      end->emitl = end->emitr = end->nodeidx = -1;
      end->nxtl = end->nxtr = NULL;

      end->prv     = parent;
      parent->nxtl = end;
    }

  if (pool != NULL)
    new = PopTracepool(pool);
  else if ((new = (struct trace_s *) malloc (sizeof(struct trace_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  new->nxtr = NULL;

  new->nxtl         = parent->nxtl;
  new->prv          = parent;
  parent->nxtl->prv = new;	/* end state also points back, to new */
  parent->nxtl      = new;

  new->emitl   = emitl;
  new->emitr   = emitr;
  new->nodeidx = nodeidx;
  new->type    = type;

  return new;
}

void
FreeTrace(struct trace_s *tr, struct trmem_s *pool)
{
  if (pool == NULL)
    {
      struct tracestack_s *stack;
      struct trace_s      *currtr;

      stack = InitTracestack();
      PushTracestack(stack, tr);

      while ((currtr = PopTracestack(stack)) != NULL)
	{
	  if (currtr->nxtr != NULL)
	    PushTracestack(stack, currtr->nxtr);
	  if (currtr->nxtl != NULL)
	    PushTracestack(stack, currtr->nxtl);
	  free(currtr);
	}
      FreeTracestack(stack);
    }
  else
    FreeTracepool(pool);
}

void
DeleteTracenode(struct trace_s *oldtr, struct trmem_s *pool)
{
  struct trace_s *parent;

  parent = oldtr->prv;

  parent->nxtl = oldtr->nxtl;
  parent->nxtr = oldtr->nxtr;
  oldtr->nxtl->prv = parent;
  if (oldtr->nxtr) oldtr->nxtr->prv = parent;
  if (pool == NULL) free(oldtr);
}

/* Functions: InitTracepool(), PopTracepool(), FreeTracepool()
 * 
 * Purpose:   Malloc() optimizations for building lots of
 *            trace trees. A "trace pool" just lets me malloc
 *            a lot of trace_s structures at once (InitTracepool).
 *            They are retrieved one at a time using PopTracepool().
 *            When done (free'ing the trace), one would call
 *            FreeTracepool().
 *            
 *            Make one trace pool per trace, if using this optimization.
 */
void
InitTracepool(struct trmem_s **ret_tmem)
{
  struct trmem_s *tmem;

  tmem = (struct trmem_s *) MallocOrDie (sizeof(struct trmem_s));
  tmem->next = 0;
  tmem->num  = TMEM_BLOCK;
  tmem->pool = (struct trace_s *) MallocOrDie (TMEM_BLOCK * sizeof(struct trace_s));
  tmem->used = InitTracestack();
  *ret_tmem = tmem;
}
struct trace_s *
PopTracepool(struct trmem_s *tmem)
{
  struct trace_s *tr;
  if (tmem->next == tmem->num)
    {				/* need a new pool */
      PushTracestack(tmem->used, tmem->pool);
      tmem->next = 0;
      tmem->num  = TMEM_BLOCK;
      tmem->pool = (struct trace_s *) MallocOrDie (TMEM_BLOCK * sizeof(struct trace_s));
    }
  tr = tmem->pool + tmem->next;
  tmem->next++;
  return tr;
}
void
FreeTracepool(struct trmem_s *tmem)
{
  struct trace_s *pop;

  while ((pop = PopTracestack(tmem->used)) != NULL)
    free(pop);
  FreeTracestack(tmem->used);
  free(tmem->pool);
  free(tmem);
}


/* Functions: InitTracestack()
 *            PushTracestack()
 *            PopTracestack()
 *            FreeTracestack()
 *            
 * Purpose:   Implementation of the pushdown stack for
 *            traversing traceback trees.
 */            
struct tracestack_s *
InitTracestack(void)
{
  struct tracestack_s *stack;

  stack = (struct tracestack_s *) MallocOrDie (sizeof(struct tracestack_s));
  stack->next = 0;
  stack->num  = TSTACK_BLOCK;
  stack->list = (struct trace_s **) MallocOrDie (sizeof(struct trace_s *) * TSTACK_BLOCK);
  return stack;
}

void 
PushTracestack(struct tracestack_s *stack, struct trace_s *tracenode)
{
  if (stack->next == stack->num)
    {
      stack->num += TSTACK_BLOCK;
      stack->list = (struct trace_s **) ReallocOrDie (stack->list, sizeof(struct trace_s *) * stack->num);
    }
  stack->list[stack->next] = tracenode;
  stack->next++;
}
struct trace_s *
PopTracestack(struct tracestack_s *stack)
{
  struct trace_s *pop;

  if (stack->next == 0) return NULL;
  stack->next--;
  pop = stack->list[stack->next];
  return pop;
}
void
FreeTracestack(struct tracestack_s *stack)
{
  free(stack->list);
  free(stack);
}



/* Function: TraceCount()
 * 
 * Purpose:  Given a trace structure, and the sequence it traces across,
 *           and a nascent model (counts form), bump the appropriate
 *           emission and transition counters in the model.
 *           
 * Return:   1 on success, 0 on failure.
 */          
int
TraceCount(struct cm_s   *cm,        /* model               */
	   char          *seq,       /* sequence, 0..len-1  */
	   double         weight,    /* weight on sequence  */
	   struct trace_s *tr)       /* traceback           */
{
  struct tracestack_s *dolist;  /* stack for traversal of traceback tree */
  struct trace_s      *curr;    /* current node in the tree              */
  int   symr, syml;
#ifdef DEBUG
  int   len;
  len = strlen(seq);
#endif

  dolist = InitTracestack();
  PushTracestack(dolist, tr->nxtl);
  
  while ((curr = PopTracestack(dolist)) != NULL)
    {
				/* ignore END states */
      if (curr->nodeidx == -1 || curr->nxtl == NULL)
	continue;			
      
				/* BIFURC states: no transits, no emission */
      if (curr->nxtr != NULL)
	{
#ifdef DEBUG
	  assert(curr->nxtr != NULL && curr->nxtl != NULL);
#endif
	  PushTracestack(dolist, curr->nxtr);
	  PushTracestack(dolist, curr->nxtl);
	}

      else if (curr->type == uINSL_ST)
	{
#ifdef DEBUG
	  assert(curr->emitl >= 0 && curr->emitl < len);
#endif
	  syml = SymbolIndex(seq[curr->emitl]);
#ifdef DEBUG
	  assert(syml >= 0 && syml < 4);
	  assert(curr->nodeidx >= 0 && curr->nodeidx < cm->nodes);
	  assert(curr->nxtl != NULL);
#endif
	  cm->nd[curr->nodeidx].tmx[INSL_ST][StatetypeIndex(curr->nxtl->type)] += weight;
	  cm->nd[curr->nodeidx].il_emit[syml] += weight;
	  PushTracestack(dolist, curr->nxtl);
	}

      else if (curr->type == uINSR_ST)
	{
#ifdef DEBUG
	  assert(curr->emitr >= 0 && curr->emitr < len);
#endif
	  symr = SymbolIndex(seq[curr->emitr]);
#ifdef DEBUG
	  assert(symr >= 0 && symr < 4);
	  assert(curr->nodeidx >= 0 && curr->nodeidx < cm->nodes);
	  assert(curr->nxtl != NULL);
#endif
	  cm->nd[curr->nodeidx].tmx[INSR_ST][StatetypeIndex(curr->nxtl->type)] += weight;
	  cm->nd[curr->nodeidx].ir_emit[symr] += weight;
	  PushTracestack(dolist, curr->nxtl);
	}

      else if (curr->type == uMATP_ST)
	{
#ifdef DEBUG
	  assert(curr->emitr >= 0 && curr->emitr < len);
	  assert(curr->emitl >= 0 && curr->emitl < len);
#endif
	  syml = SymbolIndex(seq[curr->emitl]);
	  symr = SymbolIndex(seq[curr->emitr]);
#ifdef DEBUG
	  assert(syml >= 0 && syml < 4);
	  assert(symr >= 0 && symr < 4);
	  assert(curr->nodeidx > 0 && curr->nodeidx < cm->nodes);
	  assert(curr->nxtl != NULL);
#endif
	  cm->nd[curr->nodeidx].tmx[MATP_ST][StatetypeIndex(curr->nxtl->type)] += weight;
	  cm->nd[curr->nodeidx].mp_emit[syml][symr] += weight;
	  PushTracestack(dolist, curr->nxtl);
	}
      
      else if (curr->type == uMATL_ST)
	{
#ifdef DEBUG
	  assert(curr->emitl >= 0 && curr->emitl < len);
#endif
	  syml = SymbolIndex(seq[curr->emitl]);
#ifdef DEBUG
	  assert(syml >= 0 && syml < 4);
	  assert(curr->nodeidx > 0 && curr->nodeidx < cm->nodes);
	  assert(curr->nxtl != NULL);
#endif
	  cm->nd[curr->nodeidx].tmx[MATL_ST][StatetypeIndex(curr->nxtl->type)] += weight;
	  cm->nd[curr->nodeidx].ml_emit[syml] += weight;
	  PushTracestack(dolist, curr->nxtl);
	}

      else if (curr->type == uMATR_ST)
	{
#ifdef DEBUG
	  assert(curr->emitr >= 0 && curr->emitr < len);
#endif
	  symr = SymbolIndex(seq[curr->emitr]);
#ifdef DEBUG
	  assert(symr >= 0 && symr < 4);
	  assert(curr->nodeidx > 0 && curr->nodeidx < cm->nodes);
	  assert(curr->nxtl != NULL);
#endif
	  cm->nd[curr->nodeidx].tmx[MATR_ST][StatetypeIndex(curr->nxtl->type)] += weight;
	  cm->nd[curr->nodeidx].mr_emit[symr] += weight;
	  PushTracestack(dolist, curr->nxtl);
	}

      else			/* DEL or BEGIN state */
	{
#ifdef DEBUG
	  assert(curr->nodeidx >= 0 && curr->nodeidx < cm->nodes);
	  assert(curr->nxtl->type >= 0 && curr->nxtl->type < STATETYPES);
	  assert(curr->nxtl != NULL);
#endif
	  cm->nd[curr->nodeidx].tmx[DEL_ST][StatetypeIndex(curr->nxtl->type)] += weight;
	  PushTracestack(dolist, curr->nxtl);
	}
    }

  FreeTracestack(dolist);
  return 1;
}



/* Function: TraceCountPrior()
 * 
 * Purpose:  Same as above, except that we register the counts
 *           in a prior instead of a model. Used for "training"
 *           new priors.
 *
 * Return:   1 on success, 0 on failure.
 */
int
TraceCountPrior(struct cm_s      *cm,       /* covariance model    */
		struct prior_s   *prior,    /* prior to count into */
		char             *seq,      /* sequence, 0..len-1  */
		double            weight,   /* weight on sequence  */
		struct trace_s   *tr)       /* traceback           */
{
  struct tracestack_s *dolist;  /* stack for traversal of traceback tree */
  struct trace_s      *curr;    /* current node in the tree              */
  int   symr, syml;
  int   fnode, tnode;
  int   fs, ts;

  dolist = InitTracestack();
  PushTracestack(dolist, tr->nxtl);

  while ((curr = PopTracestack(dolist)) != NULL)
    {
				/* ignore END states */
      if (curr->nodeidx == -1 || curr->nxtl == NULL)
	continue;	

			/* BIFURC states: no transits, no emission */
      if (curr->nxtr != NULL)
	{
	  PushTracestack(dolist, curr->nxtr);
	  PushTracestack(dolist, curr->nxtl);
	  continue;
	}

      syml = symr = 0;
      if (curr->emitl != -1 && !isgap(seq[curr->emitl]))
	syml = SymbolIndex(seq[curr->emitl]);
      if (curr->emitr != -1 && !isgap(seq[curr->emitr])) 
	symr = SymbolIndex(seq[curr->emitr]);
      fnode = cm->nd[curr->nodeidx].type;
      tnode = (cm->nd[curr->nodeidx].nxt != -1) ? cm->nd[cm->nd[curr->nodeidx].nxt].type : END_NODE;
      fs    = StatetypeIndex(curr->type);
      ts    = (cm->nd[curr->nodeidx].nxt != -1) ? StatetypeIndex(curr->nxtl->type) : END_ST;

      /* Verify where we're writing in memory. Had some problems here!
       */
      if (fnode < 0 || fnode > 6) Die("fnode is %d", fnode);
      if (tnode < 0 || tnode > 3) Die("tnode is %d", tnode);
      if (fs < 0 || fs >= STATETYPES) Die("fs is %d", fs);
      if (ts < 0 || ts >= STATETYPES) Die("ts is %d", ts);
      if (syml < 0 || syml >= ALPHASIZE) Die("syml is %d", syml);
      if (symr < 0 || symr >= ALPHASIZE) Die("symr is %d", symr);


      prior->tprior[fnode][tnode][fs][ts] += weight;

      switch (curr->type) {
      case uMATP_ST: prior->matp_prior[syml][symr] += weight; break;
      case uMATL_ST: prior->matl_prior[syml] += weight;       break;
      case uMATR_ST: prior->matr_prior[symr] += weight;       break;
      case uINSL_ST: prior->insl_prior[syml] += weight;       break;
      case uINSR_ST: prior->insr_prior[symr] += weight;       break;
      case uDEL_ST:  break;
      default: Die("no such state type %d", curr->type);
      }
	
      PushTracestack(dolist, curr->nxtl);
    }
  FreeTracestack(dolist);
  return 1;
}
      





/* Function: TraceScore()
 * 
 * Purpose:  Given a trace structure, and the sequence it traces across,
 *           and a model (probability form), calculate the log-odds
 *           probability score.
 * 
 *           
 * Return:   1 on success, 0 on failure.
 */          
double
TraceScore(struct cm_s   *cm,        /* model               */
	   char          *seq,       /* sequence, 0..len-1  */
	   struct trace_s *tr)       /* traceback           */
{
  struct tracestack_s *dolist;  /* stack for traversal of traceback tree */
  struct trace_s      *curr;    /* current node in the tree              */
  int   symr, syml;
  double  score;

  score  = 0;
  dolist = InitTracestack();
  PushTracestack(dolist, tr->nxtl);
  while ((curr = PopTracestack(dolist)) != NULL)
    {
				/* ignore END states */
      if (curr->nodeidx == -1 || curr->nxtl == NULL)
	continue;			
      
				/* BIFURC states: no transits, no emission */
      if (curr->nxtr != NULL)
	{
	  PushTracestack(dolist, curr->nxtr);
	  PushTracestack(dolist, curr->nxtl);
	}

      else if (curr->type == uINSL_ST)
	{
	  syml = SymbolIndex(seq[curr->emitl]);
	  score += log(cm->nd[curr->nodeidx].tmx[INSL_ST][StatetypeIndex(curr->nxtl->type)]);
	  score += log(cm->nd[curr->nodeidx].il_emit[syml]);
	  score += log(4.0);		/* for log-odds */
	  PushTracestack(dolist, curr->nxtl);
	}

      else if (curr->type == uINSR_ST)
	{
	  symr = SymbolIndex(seq[curr->emitr]);
	  score += log(cm->nd[curr->nodeidx].tmx[INSR_ST][StatetypeIndex(curr->nxtl->type)]);
	  score += log(cm->nd[curr->nodeidx].ir_emit[symr]);
	  score += log(4.0);		/* for log-odds */
	  PushTracestack(dolist, curr->nxtl);
	}

      else if (curr->type == uMATP_ST)
	{
	  syml = SymbolIndex(seq[curr->emitl]);
	  symr = SymbolIndex(seq[curr->emitr]);
	  score += log(cm->nd[curr->nodeidx].tmx[MATP_ST][StatetypeIndex(curr->nxtl->type)]);
	  score += log(cm->nd[curr->nodeidx].mp_emit[syml][symr]);
	  score += log(16.0);		/* for log-odds */
	  PushTracestack(dolist, curr->nxtl);
	}
      
      else if (curr->type == uMATL_ST)
	{
	  syml = SymbolIndex(seq[curr->emitl]);
	  score += log(cm->nd[curr->nodeidx].tmx[MATL_ST][StatetypeIndex(curr->nxtl->type)]);
	  score += log(cm->nd[curr->nodeidx].ml_emit[syml]);
	  score += log(4.0);		/* for log-odds */
	  PushTracestack(dolist, curr->nxtl);
	}

      else if (curr->type == uMATR_ST)
	{

	  symr = SymbolIndex(seq[curr->emitr]);
	  score += log(cm->nd[curr->nodeidx].tmx[MATR_ST][StatetypeIndex(curr->nxtl->type)]);
	  score += log(cm->nd[curr->nodeidx].mr_emit[symr]);
	  score += log(4.0);		/* for log-odds */
	  PushTracestack(dolist, curr->nxtl);
	}

      else			/* DEL or BEGIN state */
	{
	  score += log(cm->nd[curr->nodeidx].tmx[DEL_ST][StatetypeIndex(curr->nxtl->type)]);
	  PushTracestack(dolist, curr->nxtl);
	}
    }

  FreeTracestack(dolist);
  score = score / log(2.0);	/* convert to bits */
  return score;
}

