/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* alignio.c
 * SRE, Mon Jul 12 11:57:37 1993
 * 
 * Input/output of sequence alignments.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Function: FreeAlignment()
 * 
 * Purpose:  Free the space allocated to alignment, names, and optional
 *           information. 
 *           
 * Args:     aseqs - sequence alignment
 *           nseq  - number of sequences
 *           ainfo - optional extra data. May be NULL.
 */                  
void
FreeAlignment(char **aseqs, int nseq, struct aliinfo_s *ainfo)
{
  int i;

  for (i = 0; i < nseq; i++)
    {
      if (ainfo->sqinfo[i].flags & SQINFO_SS)   free(ainfo->sqinfo[i].ss);
      if (ainfo->sqinfo[i].flags & SQINFO_SA)   free(ainfo->sqinfo[i].sa);
    }
  if (ainfo->flags & AINFO_CS) free(ainfo->cs);
  if (ainfo->flags & AINFO_RF) free(ainfo->rf);
  free(ainfo->sqinfo);
  Free2DArray(aseqs, nseq);
}

/* Function: MakeAlignedString()
 * 
 * Purpose:  Given a raw string of some type (secondary structure, say),
 *           align it to a given aseq by putting gaps wherever the
 *           aseq has gaps. 
 *           
 * Args:     aseq:  template for alignment
 *           alen:  length of aseq
 *           ss:    raw string to align to aseq
 *           ret_s: RETURN: aligned ss
 *           
 * Return:   1 on success, 0 on failure (and squid_errno is set.)
 *           ret_ss is malloc'ed here and must be free'd by caller.
 */
int
MakeAlignedString(char *aseq, int alen, char *ss, char **ret_s)
{
  char *new; 
  int   apos, rpos;

  if ((new = (char *) malloc ((alen+1) * sizeof(char))) == NULL)
    { squid_errno = SQERR_MEM; return 0; }
  for (apos = rpos = 0; apos < alen; apos++)
    if (! isgap(aseq[apos]))
      {
	new[apos] = ss[rpos];
	rpos++;
      }
    else
      new[apos] = '.';
  new[apos] = '\0';

  if (rpos != strlen(ss))
    { squid_errno = SQERR_PARAMETER; free(new); return 0; }
  *ret_s = new;
  return 1;
}


/* Function: MakeDealignedString()
 * 
 * Purpose:  Given an aligned string of some type (either sequence or 
 *           secondary structure, for instance), dealign it relative
 *           to a given aseq. Return a ptr to the new string.
 *           
 * Args:     aseq  : template alignment 
 *           alen  : length of aseq
 *           ss:   : string to make dealigned copy of; same length as aseq
 *           ret_s : RETURN: dealigned copy of ss
 *           
 * Return:   1 on success, 0 on failure (and squid_errno is set)
 *           ret_s is alloc'ed here and must be freed by caller
 */
int
MakeDealignedString(char *aseq, int alen, char *ss, char **ret_s)
{
  char *new; 
  int   apos, rpos;

  if ((new = (char *) malloc ((alen+1) * sizeof(char))) == NULL)
    { squid_errno = SQERR_MEM; return 0; }
  for (apos = rpos = 0; apos < alen; apos++)
    if (! isgap(aseq[apos]))
      {
	new[rpos] = ss[apos];
	rpos++;
      }
  new[rpos] = '\0';
  if (alen != strlen(ss))
    { squid_errno = SQERR_PARAMETER; free(new); return 0; }
  *ret_s = new;
  return 1;
}



/* Function: WritePairwiseAlignment()
 * 
 * Purpose:  Write a nice formatted pairwise alignment out,
 *           with a BLAST-style middle line showing identities
 *           as themselves (single letter) and conservative
 *           changes as '+'.
 *           
 * Args:     ofp          - open fp to write to (stdout, perhaps)
 *           aseq1, aseq2 - alignments to write (not necessarily 
 *                          flushed right with gaps)
 *           name1, name2 - names of sequences
 *           spos1, spos2 - starting position in each (raw) sequence
 *           pam          - PAM matrix; positive values define
 *                          conservative changes
 *           indent       - how many extra spaces to print on left
 *           
 * Return:  1 on success, 0 on failure
 */
int
WritePairwiseAlignment(FILE *ofp,
		       char *aseq1, char *name1, int spos1,
		       char *aseq2, char *name2, int spos2,
		       int **pam, int indent)
{
  char sname1[11];              /* shortened name               */
  char sname2[11];             
  int  still_going;		/* True if writing another block */
  char buf1[61];		/* buffer for writing seq1; CPL+1*/
  char bufmid[61];              /* buffer for writing consensus  */
  char buf2[61];
  char *s1, *s2;                /* ptrs into each sequence          */
  int  count1, count2;		/* number of symbols we're writing  */
  int  rpos1, rpos2;		/* position in raw seqs             */
  int  rawcount1, rawcount2;	/* number of nongap symbols written */
  int  apos;

  strncpy(sname1, name1, 10);
  sname1[10] = '\0';
  strtok(sname1, WHITESPACE);

  strncpy(sname2, name2, 10);
  sname2[10] = '\0';
  strtok(sname2, WHITESPACE);

  s1 = aseq1; 
  s2 = aseq2;
  rpos1 = spos1;
  rpos2 = spos2;

  still_going = True;
  while (still_going)
    {
      still_going = False;
      
				/* get next line's worth from both */
      strncpy(buf1, s1, 60); buf1[60] = '\0';
      strncpy(buf2, s2, 60); buf2[60] = '\0';
      count1 = strlen(buf1);
      count2 = strlen(buf2);

				/* is there still more to go? */
      if ((count1 == 60 && s1[60] != '\0') ||
	  (count2 == 60 && s2[60] != '\0'))
	still_going = True;

				/* shift seq ptrs by a line */
      s1 += count1;
      s2 += count2;

				/* assemble the consensus line */
      for (apos = 0; apos < count1 && apos < count2; apos++)
	{
	  if (!isgap(buf1[apos]) && !isgap(buf2[apos]))
	    {
	      if (buf1[apos] == buf2[apos])
		bufmid[apos] = buf1[apos];
	      else if (pam[buf1[apos] - 'A'][buf2[apos] - 'A'] > 0)
		bufmid[apos] = '+';
	      else
		bufmid[apos] = ' ';
	    }
	  else
	    bufmid[apos] = ' ';
	}
      bufmid[apos] = '\0';

      rawcount1 = 0;
      for (apos = 0; apos < count1; apos++)
	if (!isgap(buf1[apos])) rawcount1++;
      
      rawcount2 = 0;
      for (apos = 0; apos < count2; apos++)
	if (!isgap(buf2[apos])) rawcount2++;

      (void) fprintf(ofp, "%*s%-10.10s %5d %s %5d\n", indent, "", 
		     sname1, rpos1, buf1, rpos1 + rawcount1 -1);
      (void) fprintf(ofp, "%*s                 %s\n", indent, "",
		     bufmid);
      (void) fprintf(ofp, "%*s%-10.10s %5d %s %5d\n", indent, "", 
		     sname2, rpos2, buf2, rpos2 + rawcount2 -1);
      (void) fprintf(ofp, "\n");

      rpos1 += rawcount1;
      rpos2 += rawcount2;
    }

  return 1;
}


/* Function: MingapAlignment()
 * 
 * Purpose:  Remove all-gap columns from a multiple sequence alignment
 *           and its associated data. The alignment is assumed to be
 *           flushed (all aseqs the same length).
 */
int
MingapAlignment(char **aseqs, int num, struct aliinfo_s *ainfo)
{
  int apos;			/* position in original alignment */
  int mpos;			/* position in new alignment      */
  int idx;

  /* We overwrite aseqs, using its allocated memory.
   */
  for (apos = 0, mpos = 0; aseqs[0][apos] != '\0'; apos++)
    {
				/* check for all-gap in column */
      for (idx = 0; idx < num; idx++)
	if (! isgap(aseqs[idx][apos]))
	  break;
      if (idx == num) continue;
	
				/* shift alignment and ainfo */
      if (mpos != apos)
	{
	  for (idx = 0; idx < num; idx++)
	    aseqs[idx][mpos] = aseqs[idx][apos];
	  
	  if (ainfo->flags & AINFO_CS) ainfo->cs[mpos] = ainfo->cs[apos];
	  if (ainfo->flags & AINFO_RF) ainfo->rf[mpos] = ainfo->rf[apos];
	}
      mpos++;
    }
				/* null terminate everything */
  for (idx = 0; idx < num; idx++)
    aseqs[idx][mpos] = '\0';
  ainfo->alen = mpos;	/* set new length */
  if (ainfo->flags & AINFO_CS) ainfo->cs[mpos] = '\0';
  if (ainfo->flags & AINFO_RF) ainfo->rf[mpos] = '\0';
  return 1;
}



/* Function: RandomAlignment()
 * 
 * Purpose:  Create a random alignment from raw sequences.
 * 
 *           Ideally, we would like to sample an alignment from the
 *           space of possible alignments according to its probability,
 *           given a prior probability distribution for alignments.
 *           I don't see how to describe such a distribution, let alone
 *           sample it.
 *           
 *           This is a rough approximation that tries to capture some
 *           desired properties. We assume the alignment is generated
 *           by a simple HMM composed of match and insert states.
 *           Given parameters (pop, pex) for the probability of opening
 *           and extending an insertion, we can find the expected number
 *           of match states, M, in the underlying model for each sequence.
 *           We use an average M taken over all the sequences (this is
 *           an approximation. The expectation of M given all the sequence
 *           lengths is a nasty-looking summation.)
 *           
 *           M = len / ( 1 + pop ( 1 + 1/ (1-pex) ) ) 
 *           
 *           Then, we assign positions in each raw sequence onto the M match
 *           states and M+1 insert states of this "HMM", by rolling random
 *           numbers and inserting the (rlen-M) inserted positions randomly
 *           into the insert slots, taking into account the relative probability
 *           of open vs. extend.
 *           
 *           The resulting alignment has two desired properties: insertions
 *           tend to follow the HMM-like exponential distribution, and
 *           the "sparseness" of the alignment is controllable through
 *           pop and pex.
 *           
 * Args:     rseqs     - raw sequences to "align", 0..nseq-1
 *           sqinfo    - array of 0..nseq-1 info structures for the sequences
 *           nseq      - number of sequences   
 *           pop       - probability to open insertion (0<pop<1)
 *           pex       - probability to extend insertion (0<pex<1)
 *           ret_aseqs - RETURN: alignment (flushed)
 *           ainfo     - fill in: alignment info
 * 
 * Return:   1 on success, 0 on failure. Sets squid_errno to indicate cause
 *           of failure.
 */                      
int
RandomAlignment(char **rseqs, SQINFO *sqinfo, int nseq, float pop, float pex,
		char ***ret_aseqs, AINFO *ainfo)
{
  char **aseqs;                 /* RETURN: alignment   */
  int    alen;			/* length of alignment */
  int   *rlen;                  /* lengths of each raw sequence */
  int    M;			/* length of "model"   */
  int  **ins;                   /* insertion counts, 0..nseq-1 by 0..M */
  int   *master_ins;            /* max insertion counts, 0..M */
  int    apos, rpos, idx;
  int    statepos;
  int    count;
  int    minlen;

  /* calculate expected length of model, M
   */
  if ((rlen = (int *) malloc (sizeof(int) * nseq)) == NULL)
    { squid_errno = SQERR_MEM; return 0; }
  M = 0;
  minlen = 9999999;
  for (idx = 0; idx < nseq; idx++)
    {
      rlen[idx] = strlen(rseqs[idx]);
      M += rlen[idx];
      minlen = (rlen[idx] < minlen) ? rlen[idx] : minlen;
    }
  M = (int) ((float) M / (1.0 + pop * (1.0 + 1.0 / (1.0 - pex))));
  M /= nseq;
  if (M > minlen) M = minlen;

  /* make arrays that count insertions in M+1 possible insert states
   */
  if ((ins = (int **) malloc (sizeof(int *) * nseq)) == NULL ||
      (master_ins = (int *) malloc (sizeof(int) * (M+1))) == NULL)
    { squid_errno = SQERR_MEM; return 0; }
  for (idx = 0; idx < nseq; idx++)
    {
      if ((ins[idx] = (int *) malloc (sizeof(int) * (M+1))) == NULL)
	{ squid_errno = SQERR_MEM; return 0; }
      for (rpos = 0; rpos <= M; rpos++)
	ins[idx][rpos] = 0;
    }
				/* normalize */
  pop = pop / (pop+pex);
  pex = 1.0 - pop;
				/* make insertions for individual sequences */
  for (idx = 0; idx < nseq; idx++)
    {
      apos = -1;
      for (rpos = 0; rpos < rlen[idx]-M; rpos++)
	{
	  if (sre_random() < pop || apos == -1)	/* open insertion */
	    apos = CHOOSE(M+1);        /* choose 0..M */
	  ins[idx][apos]++;
	}
    }
				/* calculate master_ins, max inserts */
  alen = M;
  for (apos = 0; apos <= M; apos++)
    {
      master_ins[apos] = 0;
      for (idx = 0; idx < nseq; idx++)
	if (ins[idx][apos] > master_ins[apos])
	  master_ins[apos] = ins[idx][apos];
      alen += master_ins[apos];
    }


  /* Now, construct alignment
   */
  if ((aseqs = (char **) malloc (sizeof (char *) * nseq)) == NULL)
    { squid_errno = SQERR_MEM; return 0; }
  for (idx = 0; idx < nseq; idx++)
    if ((aseqs[idx] = (char *) malloc (sizeof(char) * (alen+1))) == NULL)
      { squid_errno = SQERR_MEM; return 0; }
  
  for (idx = 0; idx < nseq; idx++)
    {
      apos = rpos = 0;

      for (statepos = 0; statepos <= M; statepos++)
	{
	  for (count = 0; count < ins[idx][statepos]; count++)
	    aseqs[idx][apos++] = rseqs[idx][rpos++];
	  for (; count < master_ins[statepos]; count++)
	    aseqs[idx][apos++] = ' ';

	  if (statepos != M)
	    aseqs[idx][apos++] = rseqs[idx][rpos++];
	}
      aseqs[idx][alen] = '\0';
    }
  ainfo->flags = 0;
  ainfo->alen  = alen; ainfo->flags |= AINFO_ALEN;

  if ((ainfo->sqinfo = (SQINFO *) malloc (sizeof(SQINFO) * nseq)) == NULL)
    Die("malloc failed");
  for (idx = 0; idx < nseq; idx++)
    SeqinfoCopy(&(ainfo->sqinfo[idx]), &(sqinfo[idx]));

  free(rlen);
  free(master_ins);
  Free2DArray(ins, nseq);
  *ret_aseqs = aseqs;
  return 1;
}
