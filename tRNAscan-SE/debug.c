/* debug.c
 * Fri Jan 28 14:10:59 1994
 * 
 * Code specifically for debugging the package.
 */


#include <stdio.h>

#include "structs.h"
#include "funcs.h"

/* Function: UstatetypeName()
 * 
 * Purpose:  "Ustatetypes" -- unique state types -- are used in the 
 *           integer-style models of the alignment algorithms.
 *           Given such a flag, return a string representation of
 *           the unique statetype.
 */
char *
UstatetypeName(int ustatetype)
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


/* Function: StatetypeName()
 * 
 * Purpose:  Given a statetype integer, return a string representation
 *           for that statetype.
 */
char *
StatetypeName(int statetype)
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


/* Function: NodetypeName()
 * 
 * Purpose:  Given a node type integer, return a printable name
 *           for the node type.
 */
char *
NodetypeName(int nodetype)
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


/* Function: PrintViterbiAMX()
 * 
 * Purpose:  Print out a normal main matrix from the original
 *           Viterbi alignment algorithm.
 * 
 */          
void
PrintViterbiAMX(FILE            *fp,       /* usually stderr/stdout     */
		struct istate_s *icm,      /* integer model             */
		int              statenum, /* length of model in states */
		char            *seq,      /* sequence, 1..N */
		int              N,        /* length of seq  */
		int           ***amx)      /* 'A' matrix */
{
  int diff, j, y;		/* indices for three dimensions */

  for (y = 0; y < statenum; y++)
    {
      fprintf(fp, "### A Matrix for state %d, type %d (%s), from node %d\n",
	      y, icm[y].statetype, UstatetypeName(icm[y].statetype), icm[y].nodeidx);
      fprintf(fp, "     ");
      for (diff = 0; diff <= N; diff++)
	fprintf(fp, "%6d  ", diff);
      fprintf(fp, "\n");

      for (j = 0; j <= N; j++)
	{
	  fprintf(fp, "%c %3d ", ((j > 0) ? seq[j] : '*'), j);
	  for (diff = 0; diff <= j; diff++)
	    fprintf(fp, "%6d  ", amx[j][diff][y]);
	  fprintf(fp, "\n");
	}
      fprintf(fp, "\n\n");
    }
}



/* Function: PrintTrace()
 * 
 * Purpose:  Debugging tool. Print a traceback tree.
 */
void
PrintTrace(FILE *fp, struct trace_s *tr)
{
  struct tracestack_s *stack;
  struct trace_s      *currtr;

  stack = InitTracestack();
  PushTracestack(stack, tr->nxtl);
  
  fprintf(fp, "  address     emitl emitr nodeidx type     nxtl          nxtr          prv\n");

  while ((currtr = PopTracestack(stack)) != NULL)
    {
      fprintf(fp, "(%p) %3d  %3d    %3d    %3d    %p  %p  %p  %s\n", 
	      currtr,
	      currtr->emitl, currtr->emitr,
	      currtr->nodeidx, currtr->type,
	      currtr->nxtl, currtr->nxtr, currtr->prv,
	      UstatetypeName(currtr->type));
      if (currtr->nxtr != NULL)
	PushTracestack(stack, currtr->nxtr);
      if (currtr->nxtl != NULL)
	PushTracestack(stack, currtr->nxtl);
    }
  FreeTracestack(stack);
}


/* Function: PrintAli()
 * 
 * Purpose:  Debugging tool. Print out an align_s structure
 */
void
PrintAli(FILE *fp, struct align_s *ali)
{
  struct align_s *curr;

  for (curr = ali; curr != NULL; curr = curr->nxt)
    fprintf(fp, "%4d %c %c %4d %s\n",
	    curr->pos, curr->sym, curr->ss, curr->nodeidx, 
	    UstatetypeName(curr->type));
}



/* Function: PrintICM()
 * 
 * Purpose:  Print an integer-version CM, as used by the alignment
 *           algorithms.
 */
void
PrintICM(FILE *fp, struct cm_s *cm, struct istate_s *icm, int nstates)
{
  int y;
  int x;

  for (y = 0; y < nstates; y++)
    {
      fprintf(fp, "node %d (%s)   state %d (%s)\n",
	      icm[y].nodeidx, NodetypeName(cm->nd[icm[y].nodeidx].type),
	      y, UstatetypeName(icm[y].statetype));
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



