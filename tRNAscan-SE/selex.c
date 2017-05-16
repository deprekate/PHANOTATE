/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* selex.c 
 * 
 * Fri Dec  4 17:43:24 1992, SRE:
 * Reading and writing aligned sequences to/from disk files.
 * Implements a new, broader specification of SELEX format
 * and supercedes alignio.c.
 *
 * SRE, Tue Nov  9 17:40:50 1993: 
 * major revision. #= special comments and aliinfo_s optional
 * alignment info support added. Support for #=CS (consensus
 * secondary structure), #=SS (individual secondary structure),
 * #=RF (reference coordinate system), #=SQ (per-sequence header info),
 * and #=AU ("author") added.
 *
 * SRE, Mon Jan 30 14:41:49 1995:
 * #=SA side chain % surface accessibility annotation supported
 * 
 * SELEX format is documented in Docs/formats.tex.
 ****************************************************************************
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <memory.h>
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static void homogenize_gapsym(char *s, char gapsym);
static int  copy_alignment_line(char *aseq, int apos, int name_rcol, 
				char *buffer, int lcol, int rcol, char gapsym);

static char commentsyms[] = "%#";

/* Function: ReadSELEX()
 * 
 * Read multiple aligned sequences from the file seqfile.
 * Store aligned sequences in aseqs, names in names, and
 * the number of sequences in num. 
 * 
 * Memory is allocated for aseqs and names, and they must be
 * free'd by the caller.
 * 
 * If optional information is desired, a non-NULL ainfo
 * pointer is passed.
 *
 * Returns 1 on success. Returns 0 on failure and sets
 * squid_errno to indicate the cause of the failure.
 */
int
ReadSELEX(char *seqfile, char ***ret_aseqs, int *ret_num, AINFO   *ainfo)
{
  FILE    *fp;                  /* ptr to opened seqfile        */
  char   **aseqs;               /* aligned seqs                 */
  int      num;			/* number of seqs read          */
  char     buffer[LINEBUFLEN];	/* input buffer for lines       */
  char     bufcpy[LINEBUFLEN];	/* strtok'able copy of buffer   */
  struct block_struc {          /** alignment data for a block: */
    int lcol;			/* furthest left aligned sym    */
    int rcol;			/* furthest right aligned sym   */
  } *blocks;
  int      blocknum;		/* number of blocks in file     */
  char    *nptr;                /* ptr to start of name on line */
  char    *sptr;                /* ptr into sequence on line    */
  int      currnum;		/* num. seqs in given block     */
  int      currblock;		/* index for blocks             */
  int      i;			/* loop counter                 */
  int      seqidx;		/* counter for seqs             */
  int      alen;                /* length of alignment          */
  int      warn_names;          /* becomes TRUE if names don't match between blocks */
  int      headnum;		/* seqidx in per-sequence header info */
  int      currlen;
  int      count;

  /***************************************************
   * First pass across file. 
   * Count seqs, get names, determine column info
   * Determine what sorts of info are active in this file.
   ***************************************************/
  ainfo->flags = 0;
				/* open the file for reading */
  fp = fopen(seqfile, "r");
  if (fp == NULL) { squid_errno = SQERR_NOFILE; return 0; }

				/* get first line of the block 
				 * (non-comment, non-blank) */
  do
    {
      if (fgets(buffer, LINEBUFLEN, fp) == NULL)
	{ squid_errno = SQERR_NODATA; return 0; }
      strcpy(bufcpy, buffer);
      if (*buffer == '#')
	{
	  if      (strncmp(buffer, "#=CS", 4) == 0) ainfo->flags |= AINFO_CS;
	  else if (strncmp(buffer, "#=RF", 4) == 0) ainfo->flags |= AINFO_RF;
	}
    }
  while ((nptr = strtok(bufcpy, WHITESPACE)) == NULL || 
	 (strchr(commentsyms, *nptr) != NULL));

  blocknum   = 0;
  warn_names = FALSE;
  while (!feof(fp))
    {
				/* allocate for info about this block. */
      if (blocknum == 0)
	blocks = (struct block_struc *) malloc  (sizeof(struct block_struc));
      else 
	blocks = (struct block_struc *) realloc (blocks, (blocknum+1) * sizeof(struct block_struc));
      if (blocks == NULL) { squid_errno = SQERR_MEM; return 0; }
      blocks[blocknum].lcol = LINEBUFLEN+1;
      blocks[blocknum].rcol = -1;
	
      currnum = 0;
      while (nptr != NULL)	/* becomes NULL when this block ends. */
      {
				/* First block only: save names */
	if (blocknum == 0)
	  {
	    if (currnum == 0)
	      ainfo->sqinfo = (SQINFO *) malloc (sizeof(SQINFO));
	    else 
	      ainfo->sqinfo = (SQINFO *) realloc (ainfo->sqinfo, (currnum + 1) * sizeof(SQINFO));
	    if (ainfo->sqinfo == NULL)
	      { squid_errno = SQERR_MEM; return 0; }

	    ainfo->sqinfo[currnum].flags = 0;
	    SetSeqinfoString(&(ainfo->sqinfo[currnum]), nptr, SQINFO_NAME);
	  }
	else			/* in each additional block: check names */
	  {
	    if (strcmp(ainfo->sqinfo[currnum].name, nptr) != 0)
	      warn_names = TRUE;
	  }
	currnum++;

				/* check rcol, lcol */
	if ((sptr = strtok(NULL, WHITESPACE)) != NULL)
	  {
				/* is this the furthest left we've
				   seen word 2 in this block? */
	    if (sptr - bufcpy < blocks[blocknum].lcol) 
	      blocks[blocknum].lcol = sptr - bufcpy;
				/* look for right side in buffer */
	    for (sptr = buffer + strlen(buffer) - 1;  
		 strchr(WHITESPACE, *sptr) != NULL;
		 sptr --)
	      /* do nothing */ ;
	    if (sptr - buffer > blocks[blocknum].rcol)
	      blocks[blocknum].rcol = sptr - buffer;
	  }

				/* get the next line; blank line means end of block */
	do
	  {
	    if (fgets(buffer, LINEBUFLEN, fp) == NULL) 
	      { nptr = NULL; break; }
	    strcpy(bufcpy, buffer);

	    if (strncmp(buffer, "#=SS", 4) == 0) ainfo->sqinfo[currnum-1].flags |= SQINFO_SS;
	    else if (strncmp(buffer, "#=SA", 4) == 0) ainfo->sqinfo[currnum-1].flags |= SQINFO_SA;
	    else if (strncmp(buffer, "#=CS", 4) == 0) ainfo->flags |= AINFO_CS;
	    else if (strncmp(buffer, "#=RF", 4) == 0) ainfo->flags |= AINFO_RF;

	    if ((nptr = strtok(bufcpy, WHITESPACE)) == NULL) 
	      break;
	  } while (strchr(commentsyms, *nptr) != NULL);
      }


				/* check that number of sequences matches expected */
      if (blocknum == 0)
	num = currnum;
      else if (currnum != num)
	{ squid_errno = SQERR_FORMAT; return 0; }
      blocknum++;

				/* get first line of next block 
				 * (non-comment, non-blank) */
      do
	{
	  if (fgets(buffer, LINEBUFLEN, fp) == NULL) { nptr = NULL; break; }
	  strcpy(bufcpy, buffer);
	}
      while ((nptr = strtok(bufcpy, WHITESPACE)) == NULL || 
	     (strchr(commentsyms, *nptr) != NULL));
    }

  
  /***************************************************
   * Get ready for second pass:
   *   figure out the length of the alignment
   *   malloc space
   *   rewind the file
   ***************************************************/

  alen = 0;
  for (currblock = 0; currblock < blocknum; currblock++)
    alen += blocks[currblock].rcol - blocks[currblock].lcol + 1;

  rewind(fp);

  /* allocations
   */
  if ((aseqs     = (char **) malloc (num * sizeof(char *))) == NULL)
    { squid_errno = SQERR_MEM; return 0; }
  if ((ainfo->flags & AINFO_CS) &&
      (ainfo->cs = (char *) malloc ((alen+1) * sizeof(char))) == NULL)
    { squid_errno = SQERR_MEM; return 0; }
  if ((ainfo->flags & AINFO_RF) &&
      (ainfo->rf = (char *) malloc ((alen+1) * sizeof(char))) == NULL)
    { squid_errno = SQERR_MEM; return 0; }
  
  for (i = 0; i < num; i++)
    {
      if ((aseqs[i]     = (char *) malloc ((alen+1) * sizeof(char))) == NULL)
	{ squid_errno = SQERR_MEM; return 0; }
      if ((ainfo->sqinfo[i].flags & SQINFO_SS) &&
	  (ainfo->sqinfo[i].ss = (char *) malloc ((alen+1) * sizeof(char))) == NULL)
	{ squid_errno = SQERR_MEM; return 0; }
      if ((ainfo->sqinfo[i].flags & SQINFO_SA) &&
	  (ainfo->sqinfo[i].sa = (char *) malloc ((alen+1) * sizeof(char))) == NULL)
	{ squid_errno = SQERR_MEM; return 0; }
    }
  
  ainfo->alen   = alen;
  ainfo->flags |= AINFO_ALEN;

  /***************************************************
   * Second pass across file. Parse header; assemble sequences
   ***************************************************/
  /* We've now made a complete first pass over the file. We know how
   * many blocks it contains, we know the number of seqs in the first
   * block, and we know every block has the same number of blocks;
   * so we can be a bit more cavalier about error-checking as we
   * make the second pass.
   */

  /* Look for header
   */
  headnum = 0;
  for (;;)
    {
      if (fgets(buffer, LINEBUFLEN, fp) == NULL)
	{ squid_errno = SQERR_NODATA; return 0; }
      strcpy(bufcpy, buffer);
      if ((nptr = strtok(bufcpy, WHITESPACE)) == NULL) continue; /* skip blank lines */

      if (strcmp(nptr, "#=AU") == 0 && /* "author" info */
	  (sptr = strtok(NULL, "\n")) != NULL)
	{
          strncpy(ainfo->au,sptr,63);
          ainfo->au[63] = '\0';
          ainfo->flags |= AINFO_AUTH;
        }

      else if (strcmp(nptr, "#=SQ") == 0)      /* per-sequence header info */
	{
				/* first field is the name */
	  if ((sptr = strtok(NULL, WHITESPACE)) == NULL)
	    { squid_errno = SQERR_FORMAT; return 0; }
	  if (strcmp(sptr, ainfo->sqinfo[headnum].name) != 0) warn_names = TRUE;

				/* second field is the weight */
	  if ((sptr = strtok(NULL, WHITESPACE)) == NULL)
	    { squid_errno = SQERR_FORMAT; return 0; }
	  SetSeqinfoString(&(ainfo->sqinfo[headnum]), sptr, SQINFO_WGT);

				/* third field is database source id */
	  if ((sptr = strtok(NULL, WHITESPACE)) == NULL)
	    { squid_errno = SQERR_FORMAT; return 0; }
	  SetSeqinfoString(&(ainfo->sqinfo[headnum]), sptr, SQINFO_ID);

				/* fourth field is database accession number */
	  if ((sptr = strtok(NULL, WHITESPACE)) == NULL)
	    { squid_errno = SQERR_FORMAT; return 0; }
	  SetSeqinfoString(&(ainfo->sqinfo[headnum]), sptr, SQINFO_ACC);

				/* fifth field is start..stop::olen */
	  if ((sptr = strtok(NULL, ".:")) == NULL)
	    { squid_errno = SQERR_FORMAT; return 0; }
	  SetSeqinfoString(&(ainfo->sqinfo[headnum]), sptr, SQINFO_START);

	  if ((sptr = strtok(NULL, ".:")) == NULL)
	    { squid_errno = SQERR_FORMAT; return 0; }
	  SetSeqinfoString(&(ainfo->sqinfo[headnum]), sptr, SQINFO_STOP);
	  
	  if ((sptr = strtok(NULL, ":\t ")) == NULL)
	    { squid_errno = SQERR_FORMAT; return 0; }
	  SetSeqinfoString(&(ainfo->sqinfo[headnum]), sptr, SQINFO_OLEN);

				/* rest of line is optional description */
	  if ((sptr = strtok(NULL, "\n")) != NULL)
	    SetSeqinfoString(&(ainfo->sqinfo[headnum]), sptr, SQINFO_DESC);
	  
	  headnum++;
	}
      else if (strcmp(nptr, "#=CS") == 0) break;
      else if (strcmp(nptr, "#=RF") == 0) break;
      else if (strchr(commentsyms, *nptr) == NULL) break; /* non-comment, non-header */
    }
  

  currlen = 0;
  for (currblock = 0 ; currblock < blocknum; currblock++)
    {
				/* parse the block */
      seqidx = 0;
      while (nptr != NULL)
	{
				/* Consensus structure */
	  if (strcmp(nptr, "#=CS") == 0)
	    {
	      if (! copy_alignment_line(ainfo->cs, currlen, strlen(nptr)-1, 
					buffer, blocks[currblock].lcol, blocks[currblock].rcol, (char) '.'))
		{ squid_errno = SQERR_FORMAT; return 0; }
	    }

				/* Reference coordinates */
	  else if (strcmp(nptr, "#=RF") == 0)
	    {
	      if (! copy_alignment_line(ainfo->rf, currlen, strlen(nptr)-1, 
					buffer, blocks[currblock].lcol, blocks[currblock].rcol, (char) ' '))
		{ squid_errno = SQERR_FORMAT; return 0; }
	    }

				/* Individual secondary structure */
	  else if (strcmp(nptr, "#=SS") == 0)
	    {
	      if (! copy_alignment_line(ainfo->sqinfo[seqidx-1].ss, currlen, strlen(nptr)-1,
					buffer, blocks[currblock].lcol, 
					blocks[currblock].rcol, (char) '.'))
		{ squid_errno = SQERR_FORMAT; return 0; }
	    }

				/* Side chain % surface accessibility code */
	  else if (strcmp(nptr, "#=SA") == 0)
	    {
	      if (! copy_alignment_line(ainfo->sqinfo[seqidx-1].sa, currlen, strlen(nptr)-1,
					buffer, blocks[currblock].lcol, 
					blocks[currblock].rcol, (char) '.'))
		{ squid_errno = SQERR_FORMAT; return 0; }
	    }
				/* Aligned sequence; avoid unparsed machine comments */
	  else if (strncmp(nptr, "#=", 2) != 0)
	    {
	      if (! copy_alignment_line(aseqs[seqidx], currlen, strlen(nptr)-1, 
					buffer, blocks[currblock].lcol, blocks[currblock].rcol, (char) '.'))
		{ squid_errno = SQERR_FORMAT; return 0; }
	      seqidx++;
	    }

				/* get next line */
	  for (;;)
	    {
	      nptr = NULL;
	      if (fgets(buffer, LINEBUFLEN, fp) == NULL) break;	/* EOF */
	      strcpy(bufcpy, buffer);
	      if ((nptr = strtok(bufcpy, WHITESPACE)) == NULL) break; /* blank */
	      if (strncmp(buffer, "#=", 2) == 0) break;      /* machine comment */
	      if (strchr(commentsyms, *nptr) == NULL) break; /* data */
	    }
	} /* end of a block */

      currlen += blocks[currblock].rcol - blocks[currblock].lcol + 1;

				/* get line 1 of next block */
      for (;;)
	{
	  if (fgets(buffer, LINEBUFLEN, fp) == NULL) break; /* no data */
	  strcpy(bufcpy, buffer);
	  if ((nptr = strtok(bufcpy, WHITESPACE)) == NULL) continue; /* blank */
	  if (strncmp(buffer, "#=", 2) == 0)       break; /* machine comment */
	  if (strchr(commentsyms, *nptr) == NULL) break; /* non-comment */
	}
    } /* end of the file */

  /* Lengths in sqinfo are for raw sequence (ungapped),
   * and SS, SA are 0..rlen-1 not 0..alen-1.
   * Only the seqs with structures come out of here with lengths set.
   */
  for (seqidx = 0; seqidx < num; seqidx++)
    {
      int apos, rpos;
				/* secondary structures */
      if (ainfo->sqinfo[seqidx].flags & SQINFO_SS)
	{
	  for (apos = rpos = 0; apos < alen; apos++)
	    if (! isgap(aseqs[seqidx][apos]))
	      {
		ainfo->sqinfo[seqidx].ss[rpos] = ainfo->sqinfo[seqidx].ss[apos];
		rpos++;
	      }
	  ainfo->sqinfo[seqidx].ss[rpos] = '\0';
	  homogenize_gapsym(ainfo->sqinfo[seqidx].ss, (char) '.');
	}
				/* Surface accessibility */
      if (ainfo->sqinfo[seqidx].flags & SQINFO_SA)
	{
	  for (apos = rpos = 0; apos < alen; apos++)
	    if (! isgap(aseqs[seqidx][apos]))
	      {
		ainfo->sqinfo[seqidx].sa[rpos] = ainfo->sqinfo[seqidx].sa[apos];
		rpos++;
	      }
	  ainfo->sqinfo[seqidx].sa[rpos] = '\0';
	  homogenize_gapsym(ainfo->sqinfo[seqidx].sa, (char) '.');
	}
    }

				/* NULL-terminate all the strings */
  if (ainfo->flags & AINFO_RF) ainfo->rf[alen] = '\0';
  if (ainfo->flags & AINFO_CS) {
    ainfo->cs[alen] = '\0';
    homogenize_gapsym(ainfo->cs, (char) '.');
  }
  for (seqidx = 0; seqidx < num; seqidx++)
    {
      aseqs[seqidx][alen]            = '\0';
      homogenize_gapsym(aseqs[seqidx], (char) '.');
    }
  
				/* find raw sequence lengths for sqinfo */
  for (seqidx = 0; seqidx < num; seqidx++)
    {
      count = 0;
      for (sptr = aseqs[seqidx]; *sptr != '\0'; sptr++)
	if (!isgap(*sptr)) count++;
      ainfo->sqinfo[seqidx].len    = count;
      ainfo->sqinfo[seqidx].flags |= SQINFO_LEN;
    }


  /***************************************************
   * Garbage collection and return
   ***************************************************/
  fclose(fp);
  free(blocks);

  if (warn_names) 
    Warn("sequences may be in different orders in blocks of %s?", seqfile);

  *ret_num   = num;
  *ret_aseqs = aseqs;
  return 1;
}




/* Function: WriteSELEX()
 * 
 * Write aligned sequences to an open file pointer,
 * breaking into multiple blocks if the sequences are
 * long. Number of symbols written per line is set by cpl.
 * The alignment must be flushed (all aseqs the same length).
 *
 * cpl cannot exceed 255.
 *
 * May also write optional information from ainfo;
 * ainfo may be NULL.
 * 
 * Returns 1 on success. Returns 0 on failure, and sets
 * squid_errno to indicate the cause.
 */
int
WriteSELEX(FILE *fp, char **aseqs, int num, struct aliinfo_s *ainfo, int cpl)
{
  int    idx;			/* counter for sequences         */
  int    namelen;		/* maximum name length used      */
  int    len;			/* tmp variable for name lengths */
  char   buffer[256];     	/* buffer for writing seq        */
  int    alen;
  int    currpos;
  char **ss;                    /* aligned secondary structure strings */
  char **sa;			/* aligned accessibility strings       */

  alen = (ainfo->flags & AINFO_ALEN) ? ainfo->alen : strlen(aseqs[0]);

			/* calculate max namelen used */
  namelen = 0;
  for (idx = 0; idx < num; idx++)
    if ((len = strlen(ainfo->sqinfo[idx].name)) > namelen) 
      namelen = len;
  if (namelen < 6) namelen = 6;


  /* Make aligned secondary structure strings
   */
  ss = (char **) MallocOrDie(sizeof(char *) * num);
  sa = (char **) MallocOrDie(sizeof(char *) * num);
  for (idx = 0; idx < num; idx++)
    {
      if (ainfo->sqinfo[idx].flags & SQINFO_SS)
	MakeAlignedString(aseqs[idx], alen, ainfo->sqinfo[idx].ss, &(ss[idx]));
      if (ainfo->sqinfo[idx].flags & SQINFO_SA)
	MakeAlignedString(aseqs[idx], alen, ainfo->sqinfo[idx].sa, &(sa[idx]));
    }

  /* Write header info
   */
  if (ainfo->flags & AINFO_AUTH)
    fprintf(fp, "#=AU %s\n", ainfo->au);

  if ((ainfo->sqinfo[0].flags & SQINFO_WGT)   ||
      (ainfo->sqinfo[0].flags & SQINFO_ID)    ||
      (ainfo->sqinfo[0].flags & SQINFO_ACC)   ||
      (ainfo->sqinfo[0].flags & SQINFO_START) ||
      (ainfo->sqinfo[0].flags & SQINFO_STOP)  ||
      (ainfo->sqinfo[0].flags & SQINFO_OLEN)  ||
      (ainfo->sqinfo[0].flags & SQINFO_DESC))
    for (idx = 0; idx < num; idx++)
      fprintf(fp, "#=SQ %-*.*s %6.4f %s %s %d..%d::%d %s\n", 
	      namelen, namelen, ainfo->sqinfo[idx].name,
	      (ainfo->sqinfo[idx].flags & SQINFO_WGT)   ? ainfo->sqinfo[idx].weight : 1.0,
	      (ainfo->sqinfo[idx].flags & SQINFO_ID)    ? ainfo->sqinfo[idx].id     : "-",
	      (ainfo->sqinfo[idx].flags & SQINFO_ACC)   ? ainfo->sqinfo[idx].id     : "-",
	      (ainfo->sqinfo[idx].flags & SQINFO_START) ? ainfo->sqinfo[idx].start  : 0,
	      (ainfo->sqinfo[idx].flags & SQINFO_STOP)  ? ainfo->sqinfo[idx].stop   : 0,
	      (ainfo->sqinfo[idx].flags & SQINFO_OLEN)  ? ainfo->sqinfo[idx].olen   : 0,
	      (ainfo->sqinfo[idx].flags & SQINFO_DESC)  ? ainfo->sqinfo[idx].desc   : "-");
  fprintf(fp, "\n");

				/* main loop: write seqs in blocks. */
  for (currpos = 0; currpos < alen; currpos += cpl)
    {
				/* Reference coord system */
      if (ainfo->flags & AINFO_RF)
	{
	  strncpy(buffer, ainfo->rf + currpos, cpl);
	  buffer[cpl] = '\0';
	  fprintf(fp, "%-*.*s  %s\n", namelen, namelen, "#=RF", buffer);
	}

				/* Consensus secondary structure */
      if (ainfo->flags & AINFO_CS)
	{
	  strncpy(buffer, ainfo->cs + currpos, cpl);
	  buffer[cpl] = '\0';
	  fprintf(fp, "%-*.*s  %s\n", namelen, namelen, "#=CS", buffer);
	}      
      
      for (idx = 0; idx < num; idx++)
	{
				/* Aligned sequence */
	  strncpy(buffer, aseqs[idx] + currpos, cpl);
	  buffer[cpl] = '\0';
	  fprintf(fp, "%-*.*s  %s\n", namelen, namelen, 
		  ainfo->sqinfo[idx].name, buffer);

				/* Individual secondary structure */
	  if (ainfo->sqinfo[idx].flags & SQINFO_SS)
	    {
	      strncpy(buffer, ss[idx] + currpos, cpl);
	      buffer[cpl] = '\0';
	      fprintf(fp, "%-*.*s  %s\n", namelen, namelen, "#=SS", buffer);
	    }

				/* Surface accessibility */
	  if (ainfo->sqinfo[idx].flags & SQINFO_SA)
	    {
	      strncpy(buffer, sa[idx] + currpos, cpl);
	      buffer[cpl] = '\0';
	      fprintf(fp, "%-*.*s  %s\n", namelen, namelen, "#=SA", buffer);
	    }
	}
				/* put blank line between blocks */
      fprintf(fp, "\n");
    }

  /* Garbage collection
   */
  for (idx = 0; idx < num; idx++)
    if (ainfo->sqinfo[idx].flags & SQINFO_SS)
      free(ss[idx]);
  free(ss);

  return 1;
}




/* Function: homogenize_gapsym()
 * 
 * Purpose:  Make gap symbols homogeneous.
 */
static void 
homogenize_gapsym(char *s, char gapsym)
{
  for (; *s != '\0'; s++)
    if (isgap(*s)) *s = gapsym; 
}
      

/* Function: copy_alignment_line()
 * 
 * Purpose:  Given a line from an alignment file, and bounds lcol,rcol
 *           on what part of it may be sequence, save the alignment into
 *           aseq starting at position apos.
 *           
 *           name_rcol is set to the rightmost column this aseqs's name
 *           occupies; if name_rcol >= lcol, we have a special case in
 *           which the name intrudes into the sequence zone.
 */
static int
copy_alignment_line(char *aseq, int apos, int name_rcol, 
		    char *buffer, int lcol, int rcol, char gapsym)
{
  char *s1, *s2;
  int   i;
  
  s1 = aseq + apos;
  s2 = buffer;			/* be careful that buffer doesn't end before lcol! */
  for (i = 0; i < lcol; i++)
    if (*s2) s2++;

  for (i = lcol; i <= rcol; i++)
    {
      if (*s2 == '\t') {
	Warn("TAB characters will corrupt a SELEX alignment! Please remove them first.");
	return 0;
      }
      if (name_rcol >= i)	/* name intrusion special case: pad left w/ gaps */
	*s1 = gapsym;
				/* short buffer special case: pad right w/ gaps  */
      else if (*s2 == '\0' || *s2 == '\n')
	*s1 = gapsym;

      else			/* normal case: copy buffer into aseq */
	*s1 = *s2;

      s1++;
      if (*s2) s2++;
    }
  return 1;
}

  
      


/* Function: DealignAseqs()
 * 
 * Given an array of (num) aligned sequences aseqs,
 * strip the gaps, represented by ' ' space characters.
 * Store the raw sequences in a new allocated array.
 * 
 * Caller is responsible for free'ing the memory allocated to
 * rseqs.
 * 
 * Returns 1 on success. Returns 0 and sets squid_errno on
 * failure.
 */
int
DealignAseqs(char **aseqs, int num, char ***ret_rseqs)
{
  char **rseqs;                 /* de-aligned sequence array   */
  int    idx;			/* counter for sequences       */
  int    depos; 		/* position counter for dealigned seq*/
  int    apos;			/* position counter for aligned seq */
  int    seqlen;		/* length of aligned seq */

				/* alloc space */
  if ((rseqs = (char **) malloc (num * sizeof(char *))) == NULL) 
    { squid_errno = SQERR_MEM; return 0; }

				/* main loop */
  for (idx = 0; idx < num; idx++)
    {
      seqlen = strlen(aseqs[idx]);
				/* alloc space */
      if ((rseqs[idx] = (char *) malloc ((seqlen + 1) * sizeof(char))) == NULL) 
	{ squid_errno = SQERR_MEM; return 0; }

				/* strip gaps */
      depos = 0;
      for (apos = 0; aseqs[idx][apos] != '\0'; apos++)
	if (!isgap(aseqs[idx][apos]))
	  {
	    rseqs[idx][depos] = aseqs[idx][apos];
	    depos++;
	  }
      rseqs[idx][depos] = '\0';
    }
  *ret_rseqs = rseqs;
  return 1;
}


/* Function: IsSELEXFormat()
 * 
 * Return TRUE if filename may be in SELEX format.
 * 
 * Accuracy is sacrificed for speed; a TRUE return does
 * *not* guarantee that the file will pass the stricter
 * error-checking of ReadSELEX(). All it checks is that
 * the first 500 non-comment lines of a file are 
 * blank, or if there's a second "word" on the line
 * it looks like sequence (i.e., it's not kOtherSeq).
 * 
 * Returns TRUE or FALSE.
 */
int
IsSELEXFormat(char *filename)
{
  FILE *fp;                     /* ptr to open sequence file */
  char  buffer[LINEBUFLEN];
  char *sptr;                   /* ptr to first word          */
  int   linenum;


  if ((fp = fopen(filename, "r")) == NULL)
    { squid_errno = SQERR_NOFILE; return 0; }

  linenum = 0;
  while (linenum < 500 && 
	 fgets(buffer, LINEBUFLEN, fp) != NULL)
    {
      linenum++;
				/* dead giveaways for extended SELEX */
      if      (strncmp(buffer, "#=AU", 4) == 0) goto DONE;
      else if (strncmp(buffer, "#=SQ", 4) == 0) goto DONE;
      else if (strncmp(buffer, "#=SS", 4) == 0) goto DONE;
      else if (strncmp(buffer, "#=CS", 4) == 0) goto DONE;
      else if (strncmp(buffer, "#=RF", 4) == 0) goto DONE;

				/* a comment? */
      if (strchr(commentsyms, *buffer) != NULL) continue;

				/* a blank line? */
      if ((sptr = strtok(buffer, WHITESPACE)) == NULL) continue;

				/* a one-word line (name only)
				   is possible, though rare */
      if ((sptr = strtok(NULL, "\n")) == NULL) continue;
      
      if (Seqtype(sptr) == kOtherSeq) {fclose(fp); return 0;}
    }

 DONE:
  fclose(fp);
  return 1;
}


/* Function: TruncateNames()
 * 
 * Make sure all names are a single word. 
 *   - if they are blank, make a name up (use the number of the sequence)
 *   - if it's already one word, leave it alone
 *   - if it's more than one word, put a terminator '\0' after the 
 *         first word 
 *         
 * Used to check an array of names before writing a SELEX-format file.
 * 
 * Returns 1 on success. Returns 0 on failure and sets squid_errno
 * to indicate the cause.
 */
int
TruncateNames(char **names, int N)
{
  int  idx;
  char newname[32];

  for (idx = 0; idx < N; idx++)
    if (names[idx] == NULL || strtok(names[idx], " \t\n") == NULL)
      {
	(void) sprintf(newname, "%d", idx);
	if (names[idx] != NULL) free(names[idx]);
	names[idx] = Strdup(newname);
      }
  return 1;
}






