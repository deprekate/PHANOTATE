/* probify.c
 * Convert counts to probabilities, using regularization.
 * SRE, Fri Sep  3 08:02:49 1993
 *
 * Because our training sequence set is finite (and often small)
 * and our number of free parameters is large, we have a small
 * sample statistics problem in estimating each parameter.
 *
 * The simplest way to deal with this problem is to use the
 * so-called Laplace law of succession. If I have N sequences
 * with a base symbol assigned to some state, and see n A's,
 * I calculate the emission probability of A as (n+1)/(N+4),
 * i.e., adding 1 to the numerator and the number of possible
 * outcomes (4) to the denominator. A summary of the proof of
 * this appears in Berg & von Hippel (J Mol Biol 193:723-750, 1987).
 * It is referred to as the "plus-one prior" by David Haussler
 * and by us.
 *
 * The plus-one prior implies that we have no knowledge at all about
 * the prior probabilities; in absence of any data, the probabilities
 * default to 1/4. What if we do have prior information? (For instance,
 * for state transitions, we know that deletes and inserts are
 * relatively rare.) We use a generalization of the Laplace law
 * of succession:
 *                n(x) + alpha * R(x)
 *         P(x) = -------------------
 *                ---
 * 		  \   n(i) + alpha * R(i)
 * 		  /__
 * 		   i
 * 
 * Here, R(x) is a "regularizer" and alpha is a weight applied
 * to the regularizer. Both were 1.0 in the plus-one prior.
 * Now, we can bias R(x) to reflect our prior expectations
 * about the probability distribution P(x). (In practice,
 * we calculate R(x) by specifying a prior probability distribution
 * and multiplying each term by the number of possible outcomes.)
 * alpha is a "confidence" term; the higher alpha is, the more
 * data it takes to outweigh the prior. We usually set alpha to
 * 1.0, but sometimes -- such as for insert state emissions, 
 * where we may assert that the emission probabilities are
 * the same as random regardless of the data -- we might use
 * arbitrarily high alpha's to freeze certain probability distributions
 * at their priors.
 * 
 * All this follows the description in Krogh et. al's HMM paper
 * (in press, JMB, 1993). 
 * 
 */


#include "structs.h"
#include "funcs.h"


/* Function: ProbifyCM()
 * 
 * Purpose:  Convert all the state transitions and symbol emissions in
 *           a covariance model from counts to probabilities.
 *           
 * Args:     cm - the model to convert
 *                
 * Return:   (void). Counts in cm become probabilities.
 */
void
ProbifyCM(struct cm_s *cm, struct prior_s *prior)
{
  int k;

  for (k = 0; k < cm->nodes; k++)
    {
      if (cm->nd[k].type != BIFURC_NODE)
	{
	  if (cm->nd[k].nxt == -1)
	    ProbifyTransitionMatrix(cm->nd[k].tmx, cm->nd[k].type, END_NODE, prior);
	  else
	    ProbifyTransitionMatrix(cm->nd[k].tmx, cm->nd[k].type, cm->nd[k+1].type, prior);
	}
      
      switch (cm->nd[k].type) {
      case MATP_NODE:
	ProbifySingletEmission(cm->nd[k].il_emit, uINSL_ST, prior);
	ProbifySingletEmission(cm->nd[k].ir_emit, uINSR_ST, prior);
	ProbifySingletEmission(cm->nd[k].ml_emit, uMATL_ST, prior);
	ProbifySingletEmission(cm->nd[k].mr_emit, uMATR_ST, prior);
	ProbifyPairEmission(cm->nd[k].mp_emit, prior);
	break;

      case MATL_NODE:
	ProbifySingletEmission(cm->nd[k].il_emit, uINSL_ST, prior);
	ProbifySingletEmission(cm->nd[k].ml_emit, uMATL_ST, prior);
	break;
	
      case MATR_NODE:
	ProbifySingletEmission(cm->nd[k].ir_emit, uINSR_ST, prior);
	ProbifySingletEmission(cm->nd[k].mr_emit, uMATR_ST, prior);
	break;

      case BEGINR_NODE:
	ProbifySingletEmission(cm->nd[k].il_emit, uINSL_ST, prior);
	break;

      case ROOT_NODE:
	ProbifySingletEmission(cm->nd[k].il_emit, uINSL_ST, prior);
	ProbifySingletEmission(cm->nd[k].ir_emit, uINSR_ST, prior);
	break;
	
      case BIFURC_NODE: break;
      case BEGINL_NODE: break;
      default: Die("Unrecognized node type %d at node %d", cm->nd[k].type, k);
      }
    }
}
	
	

/* Function: ProbifyTransitionMatrix()
 * 
 * Purpose:  Convert the state transition matrix between two nodes
 *           from counts to probabilities.
 *           
 * Args:     tmx:       6x6 state transition matrix of counts
 *           from_node: e.g. MATP_NODE, type of node we transit from
 *           to_node:   type of node we transit to
 *           prior:     prior probability distributions
 *           
 * Return:   (void). Values in tmx become probabilities.
 */
void
ProbifyTransitionMatrix(double          tmx[STATETYPES][STATETYPES],
			int             from_node,
			int             to_node,
			struct prior_s *prior)
{
  int    i,j;
  double denom;
  
  for (i = 0; i < STATETYPES; i++)
    {
      /* if no transitions to DEL in prior, this must be an unused vector */
      if (prior->tprior[from_node][to_node][i][0] > 0.0)
	{
	  denom = 0.0;
	  for (j = 0; j < STATETYPES; j++)
	    { 
	      tmx[i][j] = tmx[i][j] + prior->talpha[i] * prior->tprior[from_node][to_node][i][j];
	      denom += tmx[i][j];
	    }
	  for (j = 0; j < STATETYPES; j++)
	    tmx[i][j] /= denom;
	}
    }
}



/* Function: ProbifySingletEmission()
 * 
 * Purpose:  Convert a singlet emission vector from counts to probabilities.
 * 
 * Args:     emvec:     the emission vector
 *           statetype: type of state: uMATL_ST, uMATR_ST, uINSL_ST, uINSR_ST
 *           prior:     prior probability distributions 
 *           
 * Return:   (void). Values in emvec become probabilities.
 */          
void
ProbifySingletEmission(double          emvec[ALPHASIZE],
		       int             statetype,
		       struct prior_s *prior)
{
  int x;
  double denom;
  double *em_prior;

  /* Choose the correct prior probability distribution to use.
   */
  switch (statetype) {
  case uMATL_ST: em_prior = prior->matl_prior; break;
  case uMATR_ST: em_prior = prior->matr_prior; break;
  case uINSL_ST: em_prior = prior->insl_prior; break;
  case uINSR_ST: em_prior = prior->insr_prior; break;
  default: Die("statetype %d is not a singlet emitting state\n", statetype);
  }

  denom = 0.0;
  for (x = 0; x < ALPHASIZE; x++)
    {
      emvec[x] = emvec[x] + prior->emalpha[StatetypeIndex(statetype)] * em_prior[x];
      denom += emvec[x];
    }
  if (denom > 0.0)
    for (x = 0; x < ALPHASIZE; x++)
      emvec[x] /= denom;
}


/* Function: ProbifyPairEmission()
 * 
 * Purpose:  Convert a MATP pairwise emission matrix from counts to probabilities.
 * 
 * Args:     emx:       the emission matrix
 *           prior:     prior probability distributions 
 *           
 * Return:   (void). Values in emx become probabilities.
 */          
void
ProbifyPairEmission(double          emx[ALPHASIZE][ALPHASIZE],
		    struct prior_s *prior)
{
  int x,y;
  double denom;

  denom = 0.0;
  for (x = 0; x < ALPHASIZE; x++)
    for (y = 0; y < ALPHASIZE; y++)
      {
	emx[x][y] = emx[x][y] + prior->emalpha[MATP_ST] * prior->matp_prior[x][y];
	denom += emx[x][y];
      }
  if (denom > 0.0)
    for (x = 0; x < ALPHASIZE; x++)
      for (y = 0; y < ALPHASIZE; y++)
	emx[x][y] /= denom;      
}
