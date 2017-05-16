/* scorestack.c
 * SRE, Tue Aug 17 09:50:39 1993
 * 
 * For unidirectional scanning search procedures, implement
 * a score reporting system that filters out hits that overlap
 * with higher-scoring printed scores.
 * 
 * The simplest rule to use would be to keep a record of the last hit;
 * on receiving a new hit: 1) if new hit overlaps with last hit and
 * new hit is higher, replace last with new; 2) if new hit overlaps
 * with last hit and new hit is lower, ignore new; 3) if new hit
 * does not overlap with last, report last and assign new to last.
 * At end, report last. This is essentially the rule used by the
 * original hmm and cove scanning procedures.
 * 
 * There is a small weakness in this strategy, in that for three
 * hits A > B > C which all overlap, only A will be reported; but
 * although A overlaps B and B overlaps C, A may not overlap C.
 * (Will this ever happen in reality? I dunno. I don't want to
 * be surprised.)
 * 
 * Thus, this more complicated strategy.
 *    Keep a stack of last hits.
 *    On receiving a new hit:
 *      1) if new overlaps last and new >  last, push new onto stack
 *      2) if new overlaps last and new <= last, ignore new
 *      3) if new doesn't overlap, resolve stack; start new stack and
 *         push new onto it.
 *    At end: resolve stack.
 *    
 *    Stack resolution:
 *      set "previously reported hit" to -1,-1 so it won't overlap w/ anything
 *      while something is in the stack:
 *         pop top hit off stack
 *         if it overlaps with previously reported hit, continue;
 *         if it doesn't overlap, report it 
 * 
 *    Testing overlap:
 *      Given two subsequences with endpoints al,ar and bl, br,
 *      with no other knowledge, we would need to test whether any of 
 *      the four endpoints are within the opposing subsequence. 
 *      However, because we're scanning unidirectionally, we know
 *      that the new right end is greater than the old right end,
 *      so we only need to test whether the old right end  >= new left 
 *      end.
 *      
 * External function:
 * 
 * ReportScanHit()   - report a hit
 *                     or, if reported coords are -1,-1, resolve old
 *                     stack, cleanup and exit.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#include "squid.h"

/* Data structure for the stack of previous hits;
 * declarations of the functions to manipulate it.
 */
struct hitstack_s {
  int    left;			/* left coord of matched segment  */
  int    right;			/* right coord of matched segment */
  double score;			/* score of match                 */
  struct hitstack_s *nxt;       /* pointer to next elem in stack  */
};
static struct hitstack_s *init_hitstack(void);
static void push_hitstack(struct hitstack_s *hstack,int left,int right, double score);
static int  pop_hitstack(struct hitstack_s *hstack, int *ret_left, int *ret_right, double *ret_score);
static void free_hitstack(struct hitstack_s *hstack);




/* Function: ReportScanHit()
 * 
 * Purpose:  Caller reports a hit during a search scan, and
 *           provides a pointer to a function we can call to
 *           print non-overlapping hits. Caller reports
 *           -1,-1 for coords to request cleanup and end.
 *           
 *           Two special cases must be dealt with:
 *             INIT: If the hit stack hasn't been started yet,
 *                   we need to initialize it before doing
 *                   anything else
 *             END:  If coords are -1,-1, we resolve the stack
 *                   and cleanup; caller is finished with us
 *                   for now.
 *                   
 * Args:     left       - left coord of hit segment
 *           right      - right coord of hit segment
 *           score      - score of the hit
 *           print_hit  - pointer to a function to print 
 *                        nonoverlapping hits
 *                        
 * Return:   1 on success, 0 on failure.
 */          
int
ReportScanHit(int      left,
	      int      right,
	      double   score,
	      int    (*print_hit)(int,int,double))
{
  static struct hitstack_s *hstack = NULL;   /* static local handle; set to NULL on 1st entry   */
  static int oldright = -1;	             /* -1 is guaranteed not to overlap w/ first report */
  int    oldleft;
  int    newleft, newright;
  double newscore;

  /* Check whether this is first entry; 
   * init hit stack if so.
   */
  if (hstack == NULL) hstack = init_hitstack();

  
  /* Check whether we have to resolve the old stack:
   * if caller is reporting it's done (-1,-1 coords),
   * or if new hit doesn't overlap last stacked hit.
   */
  if (left > oldright || (left == -1 && right == -1))
    {
      /* Stack resolution.
       */
      oldleft = INT_MAX;
      while (pop_hitstack(hstack, &newleft, &newright, &newscore))
	{
				/* does this hit not overlap w/ previous printed one? */
	  if (newright < oldleft)
	    {
	      (*print_hit)(newleft, newright, newscore);
	      oldleft = newleft;
	    }
	}
      free_hitstack(hstack);
      hstack   = NULL;
      oldright = -1;

      /* Start new stack, if not done.
       */
      if (left != -1 || right != -1)
	{
	  hstack = init_hitstack();
	  push_hitstack(hstack, left, right, score);
	  oldright = right;
	}
    }


  /* else, they overlap; if new reported score is better than last one,
   * push new one. We're guaranteed to have something in
   * the stack, so we can use the score in hstack->nxt->score.
   * Reset oldright to be the new right edge of the stack, if we add something.
   */
  else if (score > hstack->nxt->score)
    {
      push_hitstack(hstack, left, right, score);
      oldright = right;
    }
  
  /* else, they overlap and the newly reported score
   * isn't better, so we ignore it.
   */
  return 1;
}
  







/* Functions: init_hitstack()
 *            push_hitstack()
 *            pop_hitstack()
 *            free_hitstack()
 *            
 * Purpose:   Implementation of the pushdown stack for
 *            keeping old hit positions and scores.
 *
 *            The hitstack has a dummy begin element,
 *            so the first legitimate element is
 *            hstack->nxt. The last legitimate element
 *            has a NULL nxt pointer.
 */            
static struct hitstack_s *
init_hitstack(void)
{
  struct hitstack_s *hstack;

  if ((hstack = (struct hitstack_s *) malloc (sizeof(struct hitstack_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  hstack->nxt = NULL;
  return hstack;
}
static void
push_hitstack(struct hitstack_s *hstack,
	      int                left,
	      int                right,
	      double             score)
{
  struct hitstack_s *new;

  if ((new = (struct hitstack_s *) malloc (sizeof(struct hitstack_s))) == NULL)
    Die("Memory allocation failure at %s line %d", __FILE__, __LINE__);

  new->left   = left;
  new->right  = right;
  new->score  = score;

  new->nxt    = hstack->nxt;
  hstack->nxt = new;
}
static int
pop_hitstack(struct hitstack_s *hstack,
	     int               *ret_left,
	     int               *ret_right,
	     double            *ret_score)
{
  struct hitstack_s *old;

  if (hstack->nxt == NULL) return 0;

  old         = hstack->nxt;
  hstack->nxt = old->nxt;

  *ret_left  = old->left;
  *ret_right = old->right;
  *ret_score = old->score;

  free(old);
  return 1;
}
static void
free_hitstack(struct hitstack_s *hstack)
{
  int left, right;
  double score;

  while (pop_hitstack(hstack, &left, &right, &score) != 0)
    ; /* do nothing */
  free(hstack);
}
