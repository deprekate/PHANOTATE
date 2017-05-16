/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* stack.c
 * SRE, Thu Mar  3 10:08:48 1994
 * 
 * Implementation of generic stack structures.
 */

#include <stdlib.h>
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/************************************************************
 * intstack_s implementation.
 * 
 * Functions: InitIntStack() - returns ptr to new stack
 *            PushIntStack() - (void)
 *            PopIntStack()  - returns 1 on success, 0 if stack empty
 *            FreeIntStack() - returns number of elements free'd, or 0 if 
 *                             stack was empty.
 *            
 * Implementation of the pushdown stack for storing single
 * integers.
 *************************************************************/  
struct intstack_s *
InitIntStack(void)
{
  struct intstack_s *stack;

  if ((stack = (struct intstack_s *) malloc (sizeof(struct intstack_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  stack->nxt = NULL;
  return stack;
}
void 
PushIntStack(struct intstack_s *stack, int data)
{
  struct intstack_s *new;

  if ((new = (struct intstack_s *) malloc (sizeof(struct intstack_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  new->data = data;

  new->nxt     = stack->nxt;
  stack->nxt   = new;
}

int
PopIntStack(struct intstack_s  *stack, int *ret_data)
{
  struct intstack_s *old;

  if (stack->nxt == NULL) return 0;

  old = stack->nxt;
  stack->nxt = old->nxt;

  *ret_data = old->data;
  free(old); 
  return 1;
}

void
ReverseIntStack(struct intstack_s *stack)
{
  struct intstack_s *old;
  struct intstack_s *new;

  old        = stack->nxt;
  stack->nxt = NULL;
  while (old != NULL)
    {
      new        = old;		/* remove one from top of old stack */
      old        = old->nxt;
      new->nxt   = stack->nxt;  /* push it onto new stack */
      stack->nxt = new;
    }
}

int
FreeIntStack( struct intstack_s *stack )
{
  int data;
  int count = 0;

  while (PopIntStack(stack, &data))
    count++;
  free(stack);
  return count;
}
