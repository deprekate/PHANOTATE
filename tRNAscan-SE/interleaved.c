/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* interleaved.c 
 * I/O of interleaved format multiple alignments. 
 * Modified from selex.c
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
 * SRE, Mon Sep 11 09:20:08 1995
 * selex.c generalized and simplified to make interleaved.c
 * 
 * SELEX format is documented in Docs/formats.tex.
 ****************************************************************************
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <memory.h>
#include <unistd.h>		/* SunOS 4.x isn't fully ANSI-compliant. */
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


static void homogenize_gapsym(char *s, char gapsym);
static int  copy_alignment_line(char *aseq, int apos, int name_rcol, 
				char *buffer, int lcol, int rcol);
static char commentsyms[] = "#%";



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
		    char *buffer, int lcol, int rcol)
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

      if (name_rcol >= i)                  *s1 = '.'; /* name intrusion: pad left w/ gaps */
      else if (*s2 == '\0' || *s2 == '\n') *s1 = '.'; /* short buffer: pad right w/ gaps  */
      else			           *s1 = *s2; /* normal: copy buffer into aseq */

      s1++;
      if (*s2) s2++;
    }
  return 1;
}



/* Function: is_blankline()
 * 
 * Return TRUE if line is all whitespace.
 */
static int
is_blankline(char *buffer)
{
  for (; *buffer != '\0'; buffer++) if (! isspace(*buffer)) return 0;
  return 1;
}


/* CLUSTALV and CLUSTALW support
 * 
 * parse_header() and is_dataline() functions for ClustalV and ClustalW
 * interleaved multiple alignment format
 */
/*ARGSUSED*/
static int
parse_clustal(FILE *fp, AINFO *ainfo, int *got_sqinfo)
{
  char buffer[LINEBUFLEN];

  do {
    if (fgets(buffer, LINEBUFLEN, fp) == NULL) { squid_errno = SQERR_FORMAT; return 0; }
  } while (strncmp(buffer, "CLUSTAL ", 8) != 0 ||
	   strstr(buffer, "multiple sequence alignment") == NULL);

  *got_sqinfo = FALSE;
  return 1;
}
static int 
dataline_clustal(char *buf, char *expected_name) 
{
  while (*buf && isspace(*buf)) buf++;
  if (*buf == '\0' || strchr(commentsyms, *buf) != NULL) 
    return 0;			/* blank or comment */
  if (expected_name != NULL && strncmp(buf, expected_name, strlen(expected_name) == 0))
    return 1;			/* matches expected seq name */
  for (; *buf != '\0'; buf++)
    {				/* Clustal has no coord lines to worry about */
      if (*buf == '*' || *buf == '.') continue;   /* possible consensus line */
      if (isalnum(*buf))              return 1;   /* name or seq character   */
      if (*buf != ' ' && isgap(*buf)) return 1;   /* possible all-gap line   */
    }
  return 0;
}


/* GCG MSF support
 * 
 * parse_header() and is_dataline() routines for GCG MSF alignments
 */
static int 
parse_MSF(FILE *fp, AINFO *ainfo, int *got_sqinfo)
{
  char buffer[LINEBUFLEN];
  char *sptr;
  int  nseq;
  
  /* Get first dividing line. MSF format specifies ints after MSF: and Check:
   * but we don't make sure of this 
   */
  do {
    if (fgets(buffer, LINEBUFLEN, fp) == NULL) { squid_errno = SQERR_FORMAT; return 0; }
  } while (strstr(buffer, " MSF: ") == NULL ||
	   strstr(buffer, " Check: ") == NULL ||
	   strstr(buffer, " ..") == NULL);

  /* Get names, weight from header
   */
  nseq = 0;
  /*CONSTCOND*/
  while (1) 
    {
      if (fgets(buffer, LINEBUFLEN, fp) == NULL) { squid_errno = SQERR_FORMAT; return 0; }
      if (is_blankline(buffer)) continue; 
      if (strncmp(buffer, "//", 2) == 0) break;

      sptr = strtok(buffer, WHITESPACE);
      if (sptr == NULL || strcmp(sptr, "Name:") != 0 || strstr(sptr+5, "Weight:") != 0)
	{ squid_errno = SQERR_FORMAT; return 0; }
      
      if (nseq == 0)  ainfo->sqinfo = (SQINFO *) malloc (sizeof(SQINFO));
      else  ainfo->sqinfo = (SQINFO *) realloc (ainfo->sqinfo, (nseq + 1) * sizeof(SQINFO));
      if (ainfo->sqinfo == NULL) { squid_errno = SQERR_MEM; return 0; }
      ainfo->sqinfo[nseq].flags = 0;
      
      if ((sptr = strtok(NULL, WHITESPACE)) == NULL) {squid_errno=SQERR_FORMAT; return 0; }
      SetSeqinfoString(&(ainfo->sqinfo[nseq]), sptr, SQINFO_NAME);

      while (sptr != NULL && strcmp(sptr, "Weight:") != 0) sptr = strtok(NULL, WHITESPACE);
      if ((sptr = strtok(NULL, WHITESPACE)) == NULL) {squid_errno=SQERR_FORMAT; return 0; }
      SetSeqinfoString(&(ainfo->sqinfo[nseq]), sptr, SQINFO_WGT);
      
      nseq++;
    }

  *got_sqinfo = TRUE;
  return 1;
}
static int
dataline_MSF(char *buf, char *expected_name)
{
  while (*buf && isspace(*buf)) buf++;
  if (*buf == '\0' || strchr(commentsyms, *buf) != NULL) 
    return 0;			/* blank or comment */
  if (expected_name != NULL && strncmp(buf, expected_name, strlen(expected_name) == 0))
    return 1;			/* matches expected seq name */
  for (; *buf != '\0'; buf++)
    {				/* MSF has coordinate lines to worry about */
      if (isspace(*buf))              continue;   /* no info from spaces     */
      if (isalpha(*buf)||isgap(*buf)) return 1;   /* has data on it          */
    }
  return 0;
}



/* Function: ReadInterleaved()
 * 
 * Purpose:  Read multiple aligned sequences from the file seqfile.
 *           Store aligned sequences in aseqs, names in names, and
 *           the number of sequences in num. 
 * 
 *           Memory is allocated for aseqs and names, and they must be
 *           free'd by the caller.
 * 
 *           If optional information is desired, a non-NULL ainfo
 *           pointer is passed. 
 * 
 * Args:     seqfile:        name of alignment file to read.
 *           parse_header(): routine to parse the header of the file
 *           is_dataline():  routine to determine if a line contains data
 *           ret_aseqs:      RETURN: 2D array of aligned sequences
 *           ret_anum:       RETURN: number of aligned sequences
 *           ainfo:          RETURN: optional alignment information
 *
 * Returns 1 on success. Returns 0 on failure and sets
 * squid_errno to indicate the cause of the failure.
 */
int
ReadInterleaved(char *seqfile, 
		int (*parse_header)(FILE *, AINFO *, int *),
		int (*is_dataline)(char *, char *),
		char ***ret_aseqs, int *ret_num, AINFO *ainfo)
{
  FILE    *fp;                  /* ptr to opened seqfile        */
  char   **aseqs;               /* aligned seqs                 */
  int      nseq;		/* number of seqs read          */
  char     buffer[LINEBUFLEN];	/* input buffer for lines       */
  char     bufcpy[LINEBUFLEN];	/* copy of buffer for strtok    */
  struct block_struc {          /** alignment data for a block: */
    int lcol;			/* furthest left aligned sym    */
    int rcol;			/* furthest right aligned sym   */
  } *blocks;
  int      blocknum;		/* number of blocks in file     */
  char    *sptr;                /* ptr into line during parsing */
  int      currblock;		/* index for blocks             */
  int      seqidx;		/* counter for seqs             */
  int      alen;                /* length of alignment          */
  int      warn_names;          /* becomes TRUE if names don't match between blocks */
  int      currlen;
  int      count;
  int      inblock;		/* TRUE if in a block of data   */
  long     offset;		/* used for skipping header     */
  int      got_sqinfo;		/* TRUE if header gave us sqinfo */

  /***************************************************
   * Parse the header of the file, according to caller-supplied function.
   * The parser is responsible for making sure there are
   * no non-comment lines before the first block. Comment
   * or blank lines are OK. The parser may also fill fields
   * into ainfo.
   ***************************************************/
				/* open the file for reading */
  fp = fopen(seqfile, "r");
  if (fp == NULL) { squid_errno = SQERR_NOFILE; return 0; }

  if (! (*parse_header) (fp, ainfo, &got_sqinfo)) return 0;
  offset = ftell(fp);		/* where are we in the file? */

  /***************************************************
   * First pass across file. 
   * Count seqs, get names, determine column info
   ***************************************************/
  ainfo->flags = 0;

  blocknum   = 0;
  nseq       = 0;
  inblock    = FALSE;
  blocks     = NULL;
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
	
      seqidx = 0;
      /*CONSTCOND*/
      while (1)			/* breaks out when blank line or EOF is hit, see below */
      {
				/* get a data line */
	do
	  {
	    if (fgets(buffer, LINEBUFLEN, fp) == NULL) goto BREAKOUT;
	    if (inblock && is_blankline(buffer))       goto BREAKOUT;
	  } while (! (*is_dataline)(buffer, 
				    (got_sqinfo || (seqidx < nseq && blocknum > 0)) ?
				    ainfo->sqinfo[seqidx].name : NULL));

				/* copy line for strtok'ing. set sptr to first word */
	inblock = TRUE;
	strcpy(bufcpy, buffer);
	sptr = strtok(bufcpy, WHITESPACE);

				/* First block only: save names */
	if (blocknum == 0 && !got_sqinfo)
	  {
	    if (seqidx == 0)  ainfo->sqinfo = (SQINFO *) malloc (sizeof(SQINFO));
	    else  ainfo->sqinfo = (SQINFO *) realloc (ainfo->sqinfo, (seqidx + 1) * sizeof(SQINFO));
	    if (ainfo->sqinfo == NULL) { squid_errno = SQERR_MEM; return 0; }

	    ainfo->sqinfo[seqidx].flags = 0;
	    SetSeqinfoString(&(ainfo->sqinfo[seqidx]), sptr, SQINFO_NAME);
	  }
	else			/* in each additional block: check names */
	  {
	    if (strcmp(ainfo->sqinfo[seqidx].name, sptr) != 0)
	      warn_names = TRUE;
	  }
	seqidx++;		/* bump sequence counter */


				/* check rcol, lcol */
	if ((sptr = strtok(NULL, WHITESPACE)) != NULL)
	  {
				/* is this the furthest left we've
				   seen word 2 in this block? */
	    if (sptr - bufcpy < blocks[blocknum].lcol) 
	      blocks[blocknum].lcol = sptr - bufcpy;
				/* look for right side in buffer */
	    for (sptr = buffer + strlen(buffer) - 1; isspace(*sptr); sptr --)
	      /* do nothing */ ;
	    if (sptr - buffer > blocks[blocknum].rcol)
	      blocks[blocknum].rcol = sptr - buffer;
	  }
      }

				/* check that number of sequences matches expected */
    BREAKOUT:
      if (inblock) {
	if (blocknum != 0 && seqidx != nseq) 
	  { squid_errno = SQERR_FORMAT; return 0; }
	nseq = seqidx;
	blocknum++;
	inblock = FALSE;
      }
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

  fseek(fp, offset, SEEK_SET);	/* rewind to first data block */

  /* allocations
   */
  if ((aseqs     = (char **) malloc (nseq * sizeof(char *))) == NULL)
    { squid_errno = SQERR_MEM; return 0; }
  for (seqidx = 0; seqidx < nseq; seqidx++)
    if ((aseqs[seqidx] = (char *) malloc ((alen+1) * sizeof(char))) == NULL)
      { squid_errno = SQERR_MEM; return 0; }
  ainfo->alen   = alen;
  ainfo->flags |= AINFO_ALEN;

  /***************************************************
   * Second pass across file. Parse alignment
   ***************************************************/
  /* We've now made a complete first pass over the file. We know how
   * many blocks it contains, we know the number of seqs in the first
   * block, and we know every block has the same number of blocks;
   * so we can be a bit more cavalier about error-checking as we
   * make the second pass.
   */

  currlen = 0;
  for (currblock = 0 ; currblock < blocknum; currblock++)
    {
      for (seqidx = 0; seqidx < nseq; seqidx++)
	{
				/* get next line */
	  do {
	    if (fgets(buffer, LINEBUFLEN, fp) == NULL) 
	      { squid_errno = SQERR_FORMAT; return 0; }
	  } while (! (*is_dataline)(buffer, ainfo->sqinfo[seqidx].name));

				/* find right boundary of name */
	  sptr = buffer;
	  while (*sptr && isspace(*sptr))  sptr++;
	  while (*sptr && !isspace(*sptr)) sptr++; 

				/* parse line */
	  if (! copy_alignment_line(aseqs[seqidx], currlen, sptr - buffer,
				    buffer, blocks[currblock].lcol, 
				    blocks[currblock].rcol))
	    { squid_errno = SQERR_FORMAT; return 0; }
      }
      currlen += blocks[currblock].rcol - blocks[currblock].lcol + 1;
    } 

				/* NULL-terminate all the strings */
  for (seqidx = 0; seqidx < nseq; seqidx++)
    {
      aseqs[seqidx][alen]            = '\0';
      homogenize_gapsym(aseqs[seqidx], (char) '.');
    }

				/* find raw sequence lengths for sqinfo */
  for (seqidx = 0; seqidx < nseq; seqidx++)
    {
      count = 0;
      for (sptr = aseqs[seqidx]; *sptr != '\0'; sptr++)
	if (!isgap(*sptr)) count++;
      ainfo->sqinfo[seqidx].len    = count;
      ainfo->sqinfo[seqidx].flags |= SQINFO_LEN;
    }

				/* tidy up the alignment */
  MingapAlignment(aseqs, nseq, ainfo);

  /***************************************************
   * Garbage collection and return
   ***************************************************/
  fclose(fp);
  free(blocks);

  if (warn_names) 
    Warn("sequences may be in different orders in blocks of %s?", seqfile);

  *ret_num   = nseq;
  *ret_aseqs = aseqs;
  return 1;
}



/* Function: ReadAlignment()
 * 
 * Purpose:  Given a seqfile name and format, hand it off to appropriate
 *           parser.
 *           
 *           Currently, squid can parse alignments from the following
 *           multiple sequence alignment formats:
 *               MSF     (U. of Wisconsin GCG package MSF format)
 *               SELEX   (NeXagen/CU Boulder SELEX format)
 *               CLUSTAL (Des Higgins' CLUSTALV and CLUSTALW programs)
 *               
 * Return:   1 on success; 0 on failure.
 *           Returned data should be freed by caller with FreeAlignment()
 */
int
ReadAlignment(char             *seqfile, 
	      int               format,
	      char           ***ret_aseqs,
	      int              *ret_num,
	      struct aliinfo_s *ret_ainfo)
{
  switch (format) {
  case kMSF: 
    if (! ReadInterleaved(seqfile, parse_MSF, dataline_MSF, ret_aseqs, ret_num, ret_ainfo)) 
      return 0;
    break;
  case kSelex:
    if (! ReadSELEX(seqfile, ret_aseqs, ret_num, ret_ainfo)) 
      return 0;
    break;
  case kClustal: 
    if (! ReadInterleaved(seqfile, parse_clustal, dataline_clustal, 
			  ret_aseqs, ret_num, ret_ainfo)) 
      return 0;
    break;
  default: squid_errno = SQERR_FORMAT; return 0;
  }
  return 1;
}
