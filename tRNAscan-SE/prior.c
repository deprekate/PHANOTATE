/* prior.c
 * Configure prior probability distributions
 * Mon Sep  6 09:42:16 1993
 * 
 * Designed to make it fairly easy to replace the prior with
 * different numbers, or ever completely different structures,
 * but for now we just copy a default prior out of prior.h
 */

#include <stdio.h>
#include <string.h>
#include "structs.h"
#include "prior.h"

static int get_real(FILE *fp, double *ret_real);



/* Function: DefaultPrior()
 * 
 * Purpose:  Copy the default prior from prior.h into a structure.
 * 
 * Args:     ret_prior: RETURN: new struct containing prior prob. distributions
 * 
 * Return:   1 on success, 0 on failure. ret_prior is malloc'ed here and
 *           must be free'd by caller by a free(*ret_prior).
 */
int
DefaultPrior(struct prior_s **ret_prior)
{
  struct prior_s *prior;
  int    i,j,fy,ty;		/* counters */

  if ((prior = (struct prior_s *) malloc (sizeof(struct prior_s))) == NULL)
    Die("malloc failed");
  
  for (i = 0; i < 7; i++)
    for (j = 0; j < 4; j++)
      for (fy = 0; fy < STATETYPES; fy++)
	for (ty = 0; ty < STATETYPES; ty++)
	  prior->tprior[i][j][fy][ty] = def_tprior[i][j][fy][ty];
  
  for (i = 0; i < ALPHASIZE; i++)
    {
      prior->rfreq[i]      = def_rfreq[i];
      prior->matl_prior[i] = def_matl_prior[i];
      prior->matr_prior[i] = def_matr_prior[i];
      prior->insl_prior[i] = def_insl_prior[i];
      prior->insr_prior[i] = def_insr_prior[i];
      for (j = 0; j < ALPHASIZE; j++)
	prior->matp_prior[i][j] = def_matp_prior[i][j];
    }

  for (ty = 0; ty < STATETYPES; ty++)
    {
      prior->talpha[ty]  = def_talpha[ty];
      prior->emalpha[ty] = def_emalpha[ty];
    }
  
  *ret_prior = prior;
  return 1;
}



/* Function: ReadPrior()
 * 
 * Purpose:  Get a prior from a file.
 * 
 * Return:   1 on success, 0 on failure.
 *           ret_prior is alloced here, must be free'd by caller
 */
int
ReadPrior(FILE *fp, struct prior_s **ret_prior)
{
  struct prior_s *prior;
  double         param;
  int    fnode,tnode,fs,ts;		/* counters */
  int    i,j;    

  if ((prior = (struct prior_s *) malloc (sizeof(struct prior_s))) == NULL)
    Die("malloc failed");

  /* Read rfreq's: expected symbol emission probability distribution for 
   * unrelated background sequence (the "random model")
   */
  if (!get_real(fp, &param)) return 0;
  prior->rfreq[0] = param;
  for (i = 1; i < ALPHASIZE; i++)
    {
      if (!get_real(NULL, &param)) return 0;
      prior->rfreq[i] = param;
    }

  /* Read talpha's: weights attached to state transition priors.
   * Often all 1.0
   */
  for (i = 0; i < STATETYPES; i++)
    {
      if (!get_real(NULL, &param)) return 0;
      prior->talpha[i] = param;
    }
  
  /* Read emalpha's: weights attached to symbol emission priors.
   * Often all 1.0. INSR_ST, INSL_ST sometimes set very high to
   * fix insert states as unlearnable
   */
  for (i = 0; i < STATETYPES; i++)
    {
      if (!get_real(NULL, &param)) return 0;
      prior->emalpha[i] = param;
    }

  /* Read matp_prior: symbol emission priors for matp
   */
  for (i = 0; i < ALPHASIZE; i++)
    for (j = 0; j < ALPHASIZE; j++)
      {
	if (!get_real(NULL, &param)) return 0;
	prior->matp_prior[i][j] = param;
      }

  /* Read singlet emission priors, in order MATL, MATR, INSL, INSR
   */
  for (i = 0; i < ALPHASIZE; i++)
    {
      if (!get_real(NULL, &param)) return 0;
      prior->matl_prior[i] = param;
    }
  for (i = 0; i < ALPHASIZE; i++)
    {
      if (!get_real(NULL, &param)) return 0;
      prior->matr_prior[i] = param;
    }
  for (i = 0; i < ALPHASIZE; i++)
    {
      if (!get_real(NULL, &param)) return 0;
      prior->insl_prior[i] = param;
    }
  for (i = 0; i < ALPHASIZE; i++)
    {
      if (!get_real(NULL, &param)) return 0;
      prior->insr_prior[i] = param;
    }

  /* Read tprior: state transition priors.
   * In order [7 from nodes][4 to nodes][from state][to state]
   */
  for (fnode = 0; fnode < 7; fnode++)
    for (tnode = 0; tnode < 4; tnode++)
      for (fs = 0; fs < STATETYPES; fs++)
	for (ts = 0; ts < STATETYPES; ts++)
	  {
	    if (!get_real(NULL, &param)) return 0;
	    prior->tprior[fnode][tnode][fs][ts] = param;
	  }

  *ret_prior = prior;
  return 1;
}



/* Function: WritePrior()
 * 
 * Purpose:  Write a prior to an open file pointer.
 *           The file is usable as a prior.h header file.
 * 
 * Return:   1 on success, 0 on failure.
 */
int
WritePrior(FILE *fp, struct prior_s *prior)
{
  int    fnode,tnode,fs,ts;		/* counters */
  int    i,j;

  if (fp == NULL) return 0;


  fprintf(fp, "#ifdef PRIORH_INCLUDED\n\n");
  fprintf(fp, "#include \"structs.h\"\n\n");
  
  /* Write rfreq
   */
  fprintf(fp, "static double def_rfreq[ALPHASIZE] = { ");
  for (i = 0; i < ALPHASIZE; i++)
    fprintf(fp, "%f%s", prior->rfreq[i],
	    (i == ALPHASIZE-1) ? " };\n\n" : ", ");

  /* Write talpha
   */
  fprintf(fp, "static double def_talpha[STATETYPES] =\n{ ");
  for (i = 0; i < STATETYPES; i++)
    fprintf(fp, "%f%s", prior->talpha[i],
	    (i == STATETYPES-1) ? " };\n\n" : ", ");

  /* Write emalpha
   */
  fprintf(fp, "static double def_emalpha[STATETYPES] =\n{ ");
  for (i = 0; i < STATETYPES; i++)
    fprintf(fp, "%f%s", prior->emalpha[i],
	    (i == STATETYPES-1) ? " };\n\n" : ", ");

  /* Write matp_prior
   */
  fprintf(fp, "static double def_matp_prior[ALPHASIZE][ALPHASIZE] =\n{\n");
  for (i = 0; i < ALPHASIZE; i++)
    {
      fprintf(fp, "  { ");
      for (j = 0; j < ALPHASIZE; j++)
	fprintf(fp, "%f%s", prior->matp_prior[i][j],
		(j == ALPHASIZE-1) ? " },\n" : ", ");
    }
  fprintf(fp, "};\n\n");

  /* Write the singlet emission priors, in order MATL, MATR, INSL, INSR
   */
  fprintf(fp, "static double def_matl_prior[ALPHASIZE] = { ");
  for (i = 0; i < ALPHASIZE; i++)
    fprintf(fp, "%f%s", prior->matl_prior[i],
	    (i == ALPHASIZE-1) ? " };\n\n" : ", ");
  fprintf(fp, "static double def_matr_prior[ALPHASIZE] = { ");
  for (i = 0; i < ALPHASIZE; i++)
    fprintf(fp, "%f%s", prior->matr_prior[i],
	    (i == ALPHASIZE-1) ? " };\n\n" : ", ");
  fprintf(fp, "static double def_insl_prior[ALPHASIZE] = { ");
  for (i = 0; i < ALPHASIZE; i++)
    fprintf(fp, "%f%s", prior->insl_prior[i],
	    (i == ALPHASIZE-1) ? " };\n\n" : ", ");
  fprintf(fp, "static double def_insr_prior[ALPHASIZE] = { ");
  for (i = 0; i < ALPHASIZE; i++)
    fprintf(fp, "%f%s", prior->insr_prior[i],
	    (i == ALPHASIZE-1) ? " };\n\n" : ", ");

  /* Write the state transition priors
   */
  fprintf(fp, "static double def_tprior[7][4][STATETYPES][STATETYPES] =\n{\n");
  for (fnode = 0; fnode < 7; fnode++)
    {
      fprintf(fp, "  {\n");	/* open block of 4 tonodes */
      for (tnode = 0; tnode < 4; tnode++)
	{
	  fprintf(fp, "   { "); /* open block of 6 from states */
	  for (fs = 0; fs < STATETYPES; fs++)
	    {
	      fprintf(fp, "%s", fs == 0? " { " : "      { ");
	      for (ts = 0; ts < STATETYPES; ts++)
		fprintf(fp, "%f%s", prior->tprior[fnode][tnode][fs][ts],
			(ts == STATETYPES-1) ? " },\n" : ", ");
	    }
	  fprintf(fp, "   },\n\n"); /* end block of 6 from states */
	}
      fprintf(fp, "  },\n\n"); /* end block of 4 tonodes */
    }
  fprintf(fp, "};\n");
  fprintf(fp, "#endif /* PRIORH_INCLUDED */\n");
  return 1;
}



/* Function: NormalizePrior()
 * 
 * Purpose:  convert a prior containing counts to one suitable for saving;
 *           the sum of each vector is the number of possibilities.
 *           i.e., an emission vector sums to ALPHASIZE, a state transition
 *           vector sums to the number of downstream states, and the MATP
 *           emission table sums to ALPHASIZE*ALPHASIZE.
 */
void
NormalizePrior(struct prior_s *prior)
{
  int fnode, tnode, fs, ts;     /* counters */
  double sum;
  double conn;			/* count of downstream connected states */


  /* Normalize the vectors
   */
  DNorm(prior->rfreq,      ALPHASIZE);             /* rfreq, random model */
  DNorm((double *) prior->matp_prior, ALPHASIZE * ALPHASIZE); /* MATP emission prior */
  DNorm(prior->matl_prior, ALPHASIZE);             /* MATL emission prior */
  DNorm(prior->matr_prior, ALPHASIZE);             /* MATR emission prior */
  DNorm(prior->insl_prior, ALPHASIZE);             /* INSL emission prior */
  DNorm(prior->insr_prior, ALPHASIZE);             /* INSR emission prior */

  /* Scale them to ALPHASIZE or ALPHASIZE*ALPHASIZE
   */
  DScale(prior->rfreq,      ALPHASIZE, (double) ALPHASIZE);         /* rfreq, random model */
  DScale((double *) prior->matp_prior, ALPHASIZE*ALPHASIZE, (double)(ALPHASIZE * ALPHASIZE)); /* MATP */
  DScale(prior->matl_prior, ALPHASIZE, (double) ALPHASIZE);         /* MATL emission prior */
  DScale(prior->matr_prior, ALPHASIZE, (double) ALPHASIZE);         /* MATR emission prior */
  DScale(prior->insl_prior, ALPHASIZE, (double) ALPHASIZE);         /* INSL emission prior */
  DScale(prior->insr_prior, ALPHASIZE, (double) ALPHASIZE);         /* INSR emission prior */

  /* state transition priors
   * deal with specially
   */
  for (fnode = 0; fnode < 7; fnode++)
    for (tnode = 0; tnode < 4; tnode++)
      for (fs = 0; fs < STATETYPES; fs++)
	{
	  sum = conn = 0.0;
	  for (ts = 0; ts < STATETYPES; ts++)
	    if (prior->tprior[fnode][tnode][fs][ts] > 0.0)
	      {
		conn += 1.0;
		sum  += prior->tprior[fnode][tnode][fs][ts];
	      }

	  if (sum > 0.0)
	    for (ts = 0; ts < STATETYPES; ts++)
	      prior->tprior[fnode][tnode][fs][ts] = prior->tprior[fnode][tnode][fs][ts] * conn / sum;
	}
}


/* Function: get_real()
 * 
 * Purpose:  Read next parameter from fp.
 *           This is a very general reading function for reading
 *           parameters from files with C-style comments.
 *           As long as the caller
 *           knows the order of what he's reading, he can format
 *           the file any way he wants.
 *
 *           Works somewhat like strtok. Call it with fp
 *           on first invocation; call with NULL on subsequent
 *           invocations. Can't work on multiple files
 *           simultaneously, and can't do anything else
 *           to them (lest we lose track of begin/end comments)
 *
 *           Also note that we use strtok internally, so
 *           the caller can't call strtok() between get_real calls
 *           on the same file.
 *
 * Return:   1 on success; 0 on failure (such as end of file)
 */
static int
get_real(FILE *fp, double *ret_real)
{
  static int in_comment   = 0;
  static FILE *internalfp = NULL;
  static char *lineptr    = NULL;
  static char stripbuffer[512];	/* buffer with comments removed */
  char buffer[512];
  char *sptr;
  char *stripptr;

  if (fp != NULL) { internalfp = fp; lineptr = NULL; }

  while (1)
    {
      if (lineptr == NULL)
	{
	  /* Get next line.
	   */
	  if (fgets(buffer, 512, internalfp) == NULL) return 0;
  
	  /* Preprocess the line, stripping out comments.
	   * The stripped copy goes into stripbuffer.
	   */
	  stripptr = stripbuffer;
	  for (sptr = buffer; *sptr; sptr++)
	    {
	      /* If we're in a comment, we're ignoring stuff
	       * until we see end-comment. Else, we're saving
	       * stuff until we see start-comment.
	       */
	      if (in_comment)
		{
		  if (*sptr == '*' && *(sptr+1) == '/')
		    { in_comment = 0; sptr++; }
		}
	      else
		{
		  if (*sptr == '/' && *(sptr+1) == '*')
		    { in_comment = 1; sptr++; }
		  else
		    { *stripptr = *sptr; stripptr++; }
		}
	    }
	  *stripptr = '\0';
	  lineptr = strtok(stripbuffer, WHITESPACE);
	}

      /* Now, get the first real and return it
       */
      while (lineptr != NULL)
	{
	  if (IsReal(lineptr))
	    {
	      *ret_real = atof(lineptr);
	      lineptr = strtok(NULL, WHITESPACE);
	      return 1;
	    }
	  else
	    lineptr = strtok(NULL, WHITESPACE);
	}
    }

  /*NOTREACHED*/
  return 0;
}
