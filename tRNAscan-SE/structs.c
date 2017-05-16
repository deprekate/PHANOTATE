/* structs.c
 * 1.0: SRE, Tue Jul  6 18:52:34 1993
 * 2.0: SRE, Thu Sep  9 14:19:19 1993
 *
 * Boring stuff which had better be flawless. 
 * 
 * Implementation of data structures. Includes
 * various pushdown stacks used for traversing model trees and
 * traceback trees, and linked lists used for collapsing trees
 * into linear alignments/strings.
 * 
 * Stacks have dummy start states. The end is just a NULL
 * pointer off the last state in the stack.
 * 
 * Linked lists have dummy start and end states, to facilitate
 * insertion and deletion.
 *
 * Pop functions only return values when the passed pointers
 * are non-NULL, so you can ask for whatever fields you want.
 *
 * For implementation of traceback tree structures, see trace.c.
 * For implementation of model structures, see model.c.
 */

#include <stdio.h>
#include <stdlib.h>

#include "structs.h"
#include "funcs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Function: StatetypeIndex()
 * 
 * Purpose:  Convert a unique statetype identifier to a valid
 *           array index
 */
int
StatetypeIndex(int type)
{
  switch (type) {
  case uBEGIN_ST:  return BEGIN_ST;
  case uBIFURC_ST: return BIFURC_ST;
  case uDEL_ST:    return DEL_ST;
  case uEND_ST:    return END_ST;
  case uMATP_ST:   return MATP_ST;
  case uMATL_ST:   return MATL_ST;
  case uMATR_ST:   return MATR_ST;
  case uINSR_ST:   return INSR_ST;
  case uINSL_ST:   return INSL_ST;

  default: Die("no such state, %d", type);
  }
 /*NOTREACHED*/
  return 0;
}

/* Function: UniqueStatetype()
 * 
 * Purpose: Convert an array index statetype into a unique statetype,
 *          using the additional information of what kind of node
 *          the state is from.
 */
int 
UniqueStatetype(int nodetype, int stidx)
{
  switch (stidx) {
  case DEL_ST:
    switch (nodetype) {
    case -1:          return uEND_ST;
    case BIFURC_NODE: return uBIFURC_ST;
    case BEGINL_NODE:
    case BEGINR_NODE: return uBEGIN_ST;
    default:          return uDEL_ST;
    }
  case MATP_ST:   return uMATP_ST;
  case MATL_ST:   return uMATL_ST;
  case MATR_ST:   return uMATR_ST;
  case INSR_ST:   return uINSR_ST;
  case INSL_ST:   return uINSL_ST;
  default: Die("no such state index %d", stidx);
  }
  /*NOTREACHED*/
  return 0;
}



/************************************************************
 * m2ali_s implementation.
 * 
 * Functions: Init_m2ali()
 *            Push_m2ali()
 *            Pop_m2ali()
 *            Free_m2ali()
 *            
 * Implementation of the pushdown stack for traversing a model
 * and producing an alignment as a linked list of align_s
 * structures. Must keep track of a current node in the model
 * tree (stateidx, subtype) and a current insertion point in
 * the alignment (insafter)
 *************************************************************/           
struct m2ali_s *
Init_m2ali(void)
{
  struct m2ali_s *stack;

  if ((stack = (struct m2ali_s *) malloc (sizeof(struct m2ali_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  stack->nxt = NULL;
  return stack;
}

void 
Push_m2ali(struct m2ali_s      *stack,
	   int                  nodeidx,
	   int                  type,
	   struct align_s      *after)
{
  struct m2ali_s *new;

  if ((new = (struct m2ali_s *) malloc (sizeof(struct m2ali_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  new->nodeidx= nodeidx;
  new->type   = type;
  new->after   = after;

  new->nxt     = stack->nxt;
  stack->nxt   = new;
}
int
Pop_m2ali(struct m2ali_s   *stack,
	  int              *ret_nodeidx,
	  int              *ret_type,
	  struct align_s  **ret_after)
{
  struct m2ali_s *old;

  if (stack->nxt == NULL)
    return 0;

  old = stack->nxt;
  stack->nxt = old->nxt;

  if (ret_nodeidx != NULL) *ret_nodeidx = old->nodeidx;
  if (ret_type    != NULL) *ret_type    = old->type;
  if (ret_after   != NULL) *ret_after   = old->after;
  free(old);
  return 1;
}

void
Free_m2ali( struct m2ali_s *stack )
{
  while (Pop_m2ali(stack, (int *) NULL, (int *) NULL, (struct align_s **) NULL))
    ;
  free(stack);
}



/***************************************************************
 * t2ali_s implementation.
 * 
 * Functions: Init_t2ali()
 *            Push_t2ali()
 *            Pop_t2ali()
 *            Free_t2ali()
 *            
 * Implementation of the pushdown stack for traversing a traceback 
 * and producing a linked list of align_s structures.          
 ****************************************************************/

struct t2ali_s *
Init_t2ali(void)
{
  struct t2ali_s *stack;

  if ((stack = (struct t2ali_s *) malloc (sizeof(struct t2ali_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  stack->nxt = NULL;
  return stack;
}

void 
Push_t2ali(struct t2ali_s      *stack,
	   struct trace_s      *tracenode,
	   struct align_s      *after)
{
  struct t2ali_s *new;

  if ((new = (struct t2ali_s *) malloc (sizeof(struct t2ali_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  new->tracenode = tracenode;
  new->after     = after;

  new->nxt       = stack->nxt;
  stack->nxt     = new;
}

int
Pop_t2ali(struct t2ali_s   *stack,
	  struct trace_s  **ret_tracenode,
	  struct align_s  **ret_after)
{
  struct t2ali_s *old;

  if (stack->nxt == NULL)
    return 0;

  old = stack->nxt;
  stack->nxt = old->nxt;

  if (ret_tracenode != NULL) *ret_tracenode = old->tracenode;
  if (ret_after     != NULL) *ret_after     = old->after;
  free(old);
  return 1;
}

void
Free_t2ali( struct t2ali_s *stack )
{
  while (Pop_t2ali(stack, (struct trace_s **) NULL, (struct align_s **) NULL))
    ;
  free(stack);
}




/************************************************************
 * align_s implementation
 * 
 * Functions:  Init_align()
 *             Insafter_align()
 *             Free_align()
 *             Print_align()
 * 
 * Implementation of a forward-linked list for alignment of 
 * a model to a sequence.
 ************************************************************/

struct align_s *
Init_align(void)
{
  struct align_s *head;
  struct align_s *tail;

  if ((head = (struct align_s *) malloc (sizeof(struct align_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  if ((tail = (struct align_s *) malloc (sizeof(struct align_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);

  head->sym = tail->sym = ' ';
  head->ss  = tail->ss  = ' ';
  head->pos = tail->pos = -1;
  head->nodeidx = tail->nodeidx = -1;
  head->type  = tail->type  = -1;
  head->nxt = tail;
  tail->nxt = NULL;
  return head;
}

struct align_s * 
Insafter_align(int              pos,
	       char             sym,     /* ACGU base character               */
	       char             ss,      /* <.> secondary structure character */
	       int              nodeidx,
	       int              type,
	       struct align_s  *after)
{
  struct align_s *new;

  if ((new = (struct align_s *) malloc (sizeof(struct align_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
    
  new->pos     = pos;
  new->sym     = sym;
  new->ss      = ss;
  new->nodeidx = nodeidx;
  new->type    = type;

  new->nxt     = after->nxt;
  after->nxt   = new;

  return new;
}

void
Delafter_align(struct align_s *after)
{
  struct align_s *old;

  old = after->nxt;
  after->nxt = old->nxt;
  free(old);
}

void
Free_align(struct align_s *head)
{
  struct align_s *old;

  while (head != NULL)
    {
      old  = head;
      head = head->nxt;
      free(old);
    }
}

#ifdef DEBUG
void
Print_align(struct align_s *head)
{
  struct align_s *curr;

  for (curr = head->nxt; curr->nxt != NULL; curr = curr->nxt)
    fprintf(stderr, "%2d %c %c %2d %2d\n", 
	    curr->pos,
	    curr->sym,
	    curr->ss,
	    curr->nodeidx,
	    curr->type);
}
#endif /* DEBUG */



