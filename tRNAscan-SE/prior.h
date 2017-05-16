#ifndef PRIORH_INCLUDED
#define PRIORH_INCLUDED
/* prior.h -- prior probability distributions, for regularization
 * SRE, Wed Sep  1 15:47:42 1993
 * 
 * There are 24 different node-node transitions; 6 kinds of node
 * we can come from (BIFURC_NODES are the exception, always with p=1.0
 * transitions), and four kinds of node we can go to (MATP, MATL,
 * MATR, BIFURC; END is treated the same as BIFURC).
 * 
 * There are 256 different kinds of state transitions, total.
 * For ease of indexing, we keep them in a [7][4][6][6] sparse
 * array, indexed [from node][to node][from statetype][to statetype].
 * Thus there are 1008 numbers, 752 are 0's. Any given row
 * tprior[fn][tn][fs] sums to 1.0 * the number of possible
 * transitions.
 * 
 * The tables are currently set up to do Laplace corrections; i.e.,
 * a "plus-one" prior.
 *
 * Christ. This really _is_ the best way to do it, I think.
 */

/* DO NOT CHANGE THE ORDER OF THE DEFINITIONS IN THIS FILE.
 * The prior.h header files are sometimes parsed at run-time,
 * and they're expected to be in this format.
 */

#include "structs.h"

                                      /* A     C     G     T */
static double def_rfreq[ALPHASIZE] = { 0.25, 0.25, 0.25, 0.25 };

static double def_talpha[STATETYPES] =
/* DEL_ST  MATP_ST MATL_ST MATR_ST INSL_ST INSR_ST */
{   1.00,   1.00,   1.00,   1.00,   1.00,   1.00 };

static double def_emalpha[STATETYPES] =
/* DEL_ST  MATP_ST MATL_ST MATR_ST INSL_ST INSR_ST */
{   1.00,   1.00,   1.00,   1.00,   1.00,   1.00 };

static double def_matp_prior[ALPHASIZE][ALPHASIZE] = 
{
  { 1.00, 1.00, 1.00, 1.00 },
  { 1.00, 1.00, 1.00, 1.00 },
  { 1.00, 1.00, 1.00, 1.00 },
  { 1.00, 1.00, 1.00, 1.00 },
};
				          /* A      C    G      T */
static double def_matl_prior[ALPHASIZE] = { 1.00, 1.00, 1.00, 1.00 };
static double def_matr_prior[ALPHASIZE] = { 1.00, 1.00, 1.00, 1.00 };
static double def_insl_prior[ALPHASIZE] = { 1.00, 1.00, 1.00, 1.00 };
static double def_insr_prior[ALPHASIZE] = { 1.00, 1.00, 1.00, 1.00 };

static double def_tprior[7][4][STATETYPES][STATETYPES] =
{
  {
  /* BIFURC_NODE --> BIFURC_NODE or END: never happens
   */
/* fs:     */ /* ts:  BIFURC_ST MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATP_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSR_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
		},

  /* BIFURC_NODE --> MATP_NODE: never happens.
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATP_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSR_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
		},
	      
  /* BIFURC_NODE --> MATL_NODE: never happens.
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATP_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSR_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
		},      

  /* BIFURC_NODE --> MATR_NODE: never happens
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATP_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSR_ST */   {     0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
		},
	    },


  {
  /* MATP_NODE --> BIFURC_NODE or END
   */
/* fs:     */ /* ts:  BIFURC_ST MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,     0.00,    0.00,    0.00,    1.00,    1.00 },
/* MATP_ST */   {      1.00,     0.00,    0.00,    0.00,    1.00,    1.00 },
/* MATL_ST */   {      1.00,     0.00,    0.00,    0.00,    1.00,    1.00 },
/* MATR_ST */   {      1.00,     0.00,    0.00,    0.00,    1.00,    1.00 },
/* INSL_ST */   {      1.00,     0.00,    0.00,    0.00,    1.00,    1.00 },
/* INSR_ST */   {      1.00,     0.00,    0.00,    0.00,    0.00,    1.00 },
		},

  /* MATP_NODE --> MATP_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,   1.00,    1.00,    1.00,    1.00,    1.00 },
/* MATP_ST */   {      1.00,   1.00,    1.00,    1.00,    1.00,    1.00 },
/* MATL_ST */   {      1.00,   1.00,    1.00,    1.00,    1.00,    1.00 },
/* MATR_ST */   {      1.00,   1.00,    1.00,    1.00,    1.00,    1.00 },
/* INSL_ST */   {      1.00,   1.00,    1.00,    1.00,    1.00,    1.00 },
/* INSR_ST */   {      1.00,   1.00,    1.00,    1.00,    0.00,    1.00 },
		},
	      
  /* MATP_NODE --> MATL_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,   0.00,    1.00,    0.00,    1.00,    1.00 },
/* MATP_ST */   {      1.00,   0.00,    1.00,    0.00,    1.00,    1.00 },   
/* MATL_ST */   {      1.00,   0.00,    1.00,    0.00,    1.00,    1.00 },     
/* MATR_ST */   {      1.00,   0.00,    1.00,    0.00,    1.00,    1.00 },     
/* INSL_ST */   {      1.00,   0.00,    1.00,    0.00,    1.00,    1.00 },     
/* INSR_ST */   {      1.00,   0.00,    1.00,    0.00,    0.00,    1.00 },     
		},      

  /* MATP_NODE --> MATR_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,   0.00,    0.00,    1.00,    1.00,    1.00 },
/* MATP_ST */   {      1.00,   0.00,    0.00,    1.00,    1.00,    1.00 },
/* MATL_ST */   {      1.00,   0.00,    0.00,    1.00,    1.00,    1.00 },
/* MATR_ST */   {      1.00,   0.00,    0.00,    1.00,    1.00,    1.00 },
/* INSL_ST */   {      1.00,   0.00,    0.00,    1.00,    1.00,    1.00 },
/* INSR_ST */   {      1.00,   0.00,    0.00,    1.00,    0.00,    1.00 },
		},
	    },

  {
  /* MATL_NODE --> BIFURC_NODE or END
   */
/* fs:     */ /* ts:  BIFURC_ST MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,     0.00,    0.00,    0.00,    1.00,    0.00 },
/* MATP_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      1.00,     0.00,    0.00,    0.00,    1.00,    0.00 },
/* MATR_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      1.00,     0.00,    0.00,    0.00,    1.00,    0.00 },
/* INSR_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
		},

 /* MATL_NODE --> MATP_NODE
   */
/* fs:     */ /* ts:  DEL_ST MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,   1.00,    1.00,    1.00,    1.00,    0.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      1.00,   1.00,    1.00,    1.00,    1.00,    0.00 },
/* MATR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      1.00,   1.00,    1.00,    1.00,    1.00,    0.00 },
/* INSR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
		},

  /* MATL_NODE --> MATL_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,   0.00,    1.00,    0.00,    1.00,    0.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      1.00,   0.00,    1.00,    0.00,    1.00,    0.00 },
/* MATR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      1.00,   0.00,    1.00,    0.00,    1.00,    0.00 },
/* INSR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
		},

  /* MATL_NODE --> MATR_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,   0.00,    0.00,    1.00,    1.00,    0.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      1.00,   0.00,    0.00,    1.00,    1.00,    0.00 },
/* MATR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      1.00,   0.00,    0.00,    1.00,    1.00,    0.00 },
/* INSR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
		},
            },


  {
  /* MATR_NODE --> BIFURC_NODE or END
   */
/* fs:     */ /* ts:  BIFURC_ST MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,     0.00,    0.00,    0.00,    0.00,    1.00 },
/* MATP_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      1.00,     0.00,    0.00,    0.00,    0.00,    1.00 },
/* INSL_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSR_ST */   {      1.00,     0.00,    0.00,    0.00,    0.00,    1.00 },
		},

  /* MATR_NODE --> MATP_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,   1.00,    1.00,    1.00,    0.00,    1.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },  
/* MATR_ST */   {      1.00,   1.00,    1.00,    1.00,    0.00,    1.00 },
/* INSL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSR_ST */   {      1.00,   1.00,    1.00,    1.00,    0.00,    1.00 },
		},
	      
  /* MATR_NODE --> MATL_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,   0.00,    1.00,    0.00,    0.00,    1.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      1.00,   0.00,    1.00,    0.00,    0.00,    1.00 },
/* INSL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSR_ST */   {      1.00,   0.00,    1.00,    0.00,    0.00,    1.00 },
		},      

  /* MATR_NODE --> MATR_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,   0.00,    0.00,    1.00,    0.00,    1.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      1.00,   0.00,    0.00,    1.00,    0.00,    1.00 },
/* INSL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSR_ST */   {      1.00,   0.00,    0.00,    1.00,    0.00,    1.00 },
		},
	    },


  {
  /* BEGINL_NODE --> BIFURC_NODE or END
   */
/* fs:     */ /* ts:  BIFURC_ST MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* BEGIN   */ { {      1.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATP_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSR_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
		},

  /* BEGINL_NODE --> MATP_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* BEGIN   */ { {      1.00,   1.00,    1.00,    1.00,    0.00,    0.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
		},
	      
  /* BEGINL_NODE --> MATL_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* BEGIN   */ { {      1.00,   0.00,    1.00,    0.00,    0.00,    0.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
		},      

  /* BEGINL_NODE --> MATR_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* BEGIN   */ { {      1.00,   0.00,    0.00,    1.00,    0.00,    0.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
		},
	    },

  {
  /* BEGINR_NODE --> BIFURC_NODE or END
   */
/* fs:     */ /* ts:  BIFURC_ST MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,     0.00,    0.00,    0.00,    1.00,    0.00 },
/* MATP_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      1.00,     0.00,    0.00,    0.00,    1.00,    0.00 },
/* INSR_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
		},

  /* BEGINR_NODE --> MATP_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,   1.00,    1.00,    1.00,    1.00,    0.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      1.00,   1.00,    1.00,    1.00,    1.00,    0.00 },
/* INSR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
		},
	      
  /* BEGINR_NODE --> MATL_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,   0.00,    1.00,    0.00,    1.00,    0.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      1.00,   1.00,    1.00,    1.00,    1.00,    1.00 },
/* INSR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
		},      

  /* BEGINR_NODE --> MATR_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* DEL_ST  */ { {      1.00,   0.00,    0.00,    1.00,    1.00,    0.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      1.00,   0.00,    0.00,    1.00,    1.00,    0.00 },
/* INSR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
		},
	    },


  {
  /* ROOT_NODE --> BIFURC_NODE or END
   */
/* fs:     */ /* ts:  BIFURC_ST MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* BEGIN   */ { {      1.00,     0.00,    0.00,    0.00,    1.00,    1.00 },
/* MATP_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      0.00,     0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      1.00,     0.00,    0.00,    0.00,    1.00,    1.00 },
/* INSR_ST */   {      1.00,     0.00,    0.00,    0.00,    0.00,    1.00 },
		},

  /* ROOT_NODE --> MATP_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* BEGIN   */ { {      1.00,   1.00,    1.00,    1.00,    1.00,    1.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      1.00,   1.00,    1.00,    1.00,    1.00,    1.00 },
/* INSR_ST */   {      1.00,   1.00,    1.00,    1.00,    0.00,    1.00 },
		},
	      
  /* ROOT_NODE --> MATL_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* BEGIN   */ { {      1.00,   0.00,    1.00,    0.00,    1.00,    1.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      1.00,   0.00,    1.00,    0.00,    1.00,    1.00 },     
/* INSR_ST */   {      1.00,   0.00,    1.00,    0.00,    0.00,    1.00 },     
		},      

  /* ROOT_NODE --> MATR_NODE
   */
/* fs:     */ /* ts:  DEL_ST  MATP_ST  MATL_ST  MATR_ST  INSL_ST  INSR_ST */
/* BEGIN   */ { {      1.00,   0.00,    0.00,    1.00,    1.00,    1.00 },
/* MATP_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATL_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* MATR_ST */   {      0.00,   0.00,    0.00,    0.00,    0.00,    0.00 },
/* INSL_ST */   {      1.00,   0.00,    0.00,    1.00,    1.00,    1.00 },
/* INSR_ST */   {      1.00,   0.00,    0.00,    1.00,    0.00,    1.00 },
		},
	    },
};


#endif /* PRIORH_INCLUDED */
