/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* File: sqio.c
 * From: ureadseq.c in Don Gilbert's sequence i/o package
 *
 * Reads and writes nucleic/protein sequence in various
 * formats. Data files may have multiple sequences.
 *
 * Heavily modified from READSEQ package
 * Copyright (C) 1990 by D.G. Gilbert
 * Biology Dept., Indiana University, Bloomington, IN 47405
 * email: gilbertd@bio.indiana.edu
 * Thanks Don!
 *
 * SRE: Modifications as noted. Fri Jul  3 09:44:54 1992
 *      Packaged for squid, Thu Oct  1 10:07:11 1992
 *      ANSI conversion in full swing, Mon Jul 12 12:22:21 1993
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifndef SEEK_SET
#include <unistd.h>	/* may Sun Microsystems rot in hell */
#endif

#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Our hack for "sequential" access of sequences in 
 * interleaved-format alignment files (MSF, SELEX):
 * static pointers at data for the entire alignment.
 */
static char  **ali_aseqs;
static char  **ali_rseqs;
static struct aliinfo_s ali_ainfo;
static int     ali_num;
static int     curridx;

#define kStartLength  500

static char *aminos      = "ABCDEFGHIKLMNPQRSTVWXYZ*";
static char *primenuc    = "ACGTUN";
static char *protonly    = "EFIPQZ";
/* static char  stdsymbols[6]   = "_.-*?"; */
static char  allsymbols[32]  = "_.-*?<>{}[]()!@#$%^&=+;:'|`~\"\\";
static char *seqsymbols  = allsymbols;
/*
    use general form of isseqchar -- all chars + symbols.
    no formats except nbrf (?) use symbols in data area as
    anything other than sequence chars.
    (wrong. PIR-CODATA does. Remove /)
*/



void
FreeSequence(char *seq, SQINFO *sqinfo)
{
  if (seq != NULL) free(seq);
  if (sqinfo->flags & SQINFO_SS)   free(sqinfo->ss);
  if (sqinfo->flags & SQINFO_SA)   free(sqinfo->sa);
}

int
SetSeqinfoString(SQINFO *sqinfo, char *sptr, int flag)
{
  int len;
  int pos;

  while (*sptr == ' ') sptr++; /* ignore leading whitespace */
  for (pos = strlen(sptr)-1; pos >= 0; pos--)
    if (! isspace(sptr[pos])) break;
  sptr[pos+1] = '\0';	       /* ignore trailing whitespace */

  switch (flag) {
  case SQINFO_NAME:
    if (*sptr != '-')
      { 
	strncpy(sqinfo->name, sptr, SQINFO_NAMELEN-1);
	sqinfo->name[SQINFO_NAMELEN-1] = '\0';
	sqinfo->flags   |= SQINFO_NAME;
      }
    break;

  case SQINFO_ID:
    if (*sptr != '-')
      { 
	strncpy(sqinfo->id, sptr, SQINFO_NAMELEN-1);
	sqinfo->id[SQINFO_NAMELEN-1] = '\0';
	sqinfo->flags |= SQINFO_ID;
      }
    break;

  case SQINFO_ACC:
    if (*sptr != '-')
      { 
	strncpy(sqinfo->acc, sptr, SQINFO_NAMELEN-1);
	sqinfo->acc[SQINFO_NAMELEN-1] = '\0';
	sqinfo->flags   |= SQINFO_ACC;
      }
    break;

  case SQINFO_DESC:
    if (*sptr != '-')
      { 
	if (sqinfo->flags & SQINFO_DESC) /* append? */
	  {
	    len = strlen(sqinfo->desc);
	    if (len < SQINFO_DESCLEN-2)	/* is there room? */
	      {
		strncat(sqinfo->desc, " ", SQINFO_DESCLEN-1-len); len++;
		strncat(sqinfo->desc, sptr, SQINFO_DESCLEN-1-len);
	      }
	  }
	else			/* else copy */
	  strncpy(sqinfo->desc, sptr, SQINFO_DESCLEN-1);
	sqinfo->desc[SQINFO_DESCLEN-1] = '\0';
	sqinfo->flags   |= SQINFO_DESC;
      }
    break;

  case SQINFO_START:
    if (!IsInt(sptr)) { squid_errno = SQERR_FORMAT; return 0; }
    sqinfo->start = atoi(sptr);
    if (sqinfo->start != 0) sqinfo->flags |= SQINFO_START;
    break;

  case SQINFO_STOP:
    if (!IsInt(sptr)) { squid_errno = SQERR_FORMAT; return 0; }
    sqinfo->stop = atoi(sptr);
    if (sqinfo->stop != 0) sqinfo->flags |= SQINFO_STOP;
    break;

  case SQINFO_OLEN:
    if (!IsInt(sptr)) { squid_errno = SQERR_FORMAT; return 0; }
    sqinfo->olen = atoi(sptr);
    if (sqinfo->olen != 0) sqinfo->flags |= SQINFO_OLEN;
    break;

  case SQINFO_WGT:
    if (*sptr != '-')
      {
	if (!IsReal(sptr)) { squid_errno = SQERR_FORMAT; return 0; }
	sqinfo->weight = (float) atof(sptr);
	sqinfo->flags |= SQINFO_WGT;
      }
    break;

  default:
    Die("Invalid flag %d to SetSeqinfoString()");
  }
  return 1;
}

void
SeqinfoCopy(SQINFO *sq1, SQINFO *sq2)
{
  sq1->flags = sq2->flags;
  if (sq2->flags & SQINFO_NAME)  strcpy(sq1->name, sq2->name);
  if (sq2->flags & SQINFO_ID)    strcpy(sq1->id,   sq2->id);
  if (sq2->flags & SQINFO_ACC)   strcpy(sq1->acc,  sq2->acc);
  if (sq2->flags & SQINFO_DESC)  strcpy(sq1->desc, sq2->desc);
  if (sq2->flags & SQINFO_LEN)   sq1->len    = sq2->len;
  if (sq2->flags & SQINFO_START) sq1->start  = sq2->start;
  if (sq2->flags & SQINFO_STOP)  sq1->stop   = sq2->stop;
  if (sq2->flags & SQINFO_OLEN)  sq1->olen   = sq2->olen;
  if (sq2->flags & SQINFO_TYPE)  sq1->type   = sq2->type;
  if (sq2->flags & SQINFO_WGT)   sq1->weight = sq2->weight;
  if (sq2->flags & SQINFO_SS)    sq1->ss     = Strdup(sq2->ss);
  if (sq2->flags & SQINFO_SA)    sq1->sa     = Strdup(sq2->sa);
}

/* Function: ToDNA()
 * 
 * Purpose:  Convert a sequence to DNA.
 *           U --> T
 */
void
ToDNA(char *seq)
{
  for (; *seq != '\0'; seq++)
    {
      if      (*seq == 'U') *seq = 'T';
      else if (*seq == 'u') *seq = 't';
    }
}

/* Function: ToRNA()
 * 
 * Purpose:  Convert a sequence to RNA.
 *           T --> U
 */
void
ToRNA(char *seq)
{
  for (; *seq != '\0'; seq++)
    {
      if      (*seq == 'T') *seq = 'U';
      else if (*seq == 't') *seq = 'u';
    }
}


static int 
isSeqChar(int c)
{
  if (c > 127) return 0;	/* IRIX 4.0 bug! isascii(255) returns TRUE */
  return (isalpha(c) || strchr(seqsymbols,c));
}

static void 
readline(FILE *f, char *s)
{
  char  *cp;

  if (NULL == fgets(s, LINEBUFLEN, f))
    *s = 0;
  else {
    cp = strchr(s, '\n');
    if (cp != NULL) *cp = 0;
    }
}

static void 
GetLine(struct ReadSeqVars *V)
{
  readline(V->f, V->sbuffer);
}


/* Function: addseq()
 * 
 * Purpose:  Add a line of sequence to the growing string in V.
 */
static void 
addseq(char *s, struct ReadSeqVars *V)
{
  char  *ptr;

  while (*s != 0) {
    if (isSeqChar((int) *s)) {
      if (*s == '-' && V->dash_equals_n) *s = 'N';
      if (V->seqlen >= V->maxseq) {
        V->maxseq += kStartLength;
        ptr = (char*) realloc(V->seq, V->maxseq+1);
        if (ptr==NULL) {
          squid_errno = SQERR_MEM;
          return;
	}
        else V->seq = ptr;
        }
      V->seq[(V->seqlen)++] = *s;
      }
    s++;
  }
}

static void
addstruc(char *s, struct ReadSeqVars *V)
{
  char *sptr;

  if (! (V->sqinfo->flags & SQINFO_SS))
    {
      if ((V->sqinfo->ss = (char *) malloc ((V->maxseq+1) * sizeof(char))) == NULL)
	{ squid_errno = SQERR_MEM; return; }
      V->sqinfo->flags |= SQINFO_SS;
      sptr = V->sqinfo->ss;
    }      
  else
    { 
      if ((V->sqinfo->ss = (char *) realloc(V->sqinfo->ss, V->maxseq+1)) == NULL)
	{ squid_errno = SQERR_MEM; return; }
      sptr = V->sqinfo->ss;
      while (*sptr != '\0') sptr++;
    }

  while (*s != 0)
    {
      if (isSeqChar((int)*s)) { *sptr = *s; sptr++; }
      s++;
    }
  *sptr = '\0';
}


static void 
readLoop(int addfirst, int (*endTest)(char *,int *), struct ReadSeqVars *V)
{
  int addend = 0;
  int done   = 0;

  V->seqlen = 0;
  if (addfirst) addseq(V->sbuffer, V);
  do {
    GetLine(V);
    done = feof(V->f);
    done |= (*endTest)(V->sbuffer, &addend);
    if (addend || !done)
      addseq(V->sbuffer, V);
  } while (!done);
}


static int
endPIR(char *s, int  *addend)
{
  *addend = 0;
  if ((strncmp(s, "///", 3) == 0) || 
      (strncmp(s, "ENTRY", 5) == 0))
    return 1;
  else
    return 0;
}

static void
readPIR(struct ReadSeqVars *V)
{
  char *sptr;
				/* load first line of entry  */
  while (!feof(V->f) && strncmp(V->sbuffer, "ENTRY", 5) != 0)
    GetLine(V);
  if (feof(V->f)) return;

  if ((sptr = strtok(V->sbuffer + 15, "\n\t ")) != NULL)
    {
      SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);
      SetSeqinfoString(V->sqinfo, sptr, SQINFO_ID);
    }
  do {
    GetLine(V);
    if (!feof(V->f) && strncmp(V->sbuffer, "TITLE", 5) == 0)
      SetSeqinfoString(V->sqinfo, V->sbuffer+15, SQINFO_DESC);
    else if (!feof(V->f) && strncmp(V->sbuffer, "ACCESSION", 9) == 0)
      {
	if ((sptr = strtok(V->sbuffer+15, " \t\n")) != NULL)
	  SetSeqinfoString(V->sqinfo, sptr, SQINFO_ACC);
      }
  } while (! feof(V->f) && (strncmp(V->sbuffer,"SEQUENCE", 8) != 0));
  GetLine(V);			/* skip next line, coords */

  readLoop(0, endPIR, V);

  /* reading a real PIR-CODATA database file, we keep the source coords
   */
  V->sqinfo->start = 1;
  V->sqinfo->stop  = V->seqlen;
  V->sqinfo->olen  = V->seqlen;
  V->sqinfo->flags |= SQINFO_START | SQINFO_STOP | SQINFO_OLEN;

  /* get next line
   */
  while (!feof(V->f) && strncmp(V->sbuffer, "ENTRY", 5) != 0)
    GetLine(V);
}



static int 
endIG(char *s, int  *addend)
{
  *addend = 1; /* 1 or 2 occur in line w/ bases */
  return((strchr(s,'1')!=NULL) || (strchr(s,'2')!=NULL));
}

static void 
readIG(struct ReadSeqVars *V)
{
  char *nm;
				/* position past ';' comments */
  do {
    GetLine(V);
  } while (! (feof(V->f) || ((*V->sbuffer != 0) && (*V->sbuffer != ';')) ));

  if (!feof(V->f))
    {
      if ((nm = strtok(V->sbuffer, "\n\t ")) != NULL)
	SetSeqinfoString(V->sqinfo, nm, SQINFO_NAME);

      readLoop(0, endIG, V);
    }
  
  while (!(feof(V->f) || ((*V->sbuffer != '\0') && (*V->sbuffer == ';'))))
    GetLine(V);
}

static int 
endStrider(char *s, int *addend)
{
  *addend = 0;
  return (strstr( s, "//") != NULL);
}

static void 
readStrider(struct ReadSeqVars *V)
{ 
  char *nm;
  
  while ((!feof(V->f)) && (*V->sbuffer == ';')) 
    {
      if (strncmp(V->sbuffer,"; DNA sequence", 14) == 0)
	{
	  if ((nm = strtok(V->sbuffer+16, ",\n\t ")) != NULL)
	    SetSeqinfoString(V->sqinfo, nm, SQINFO_NAME);
	}
      GetLine(V);
    }

  if (! feof(V->f))
    readLoop(1, endStrider, V);

  /* load next line
   */
  while ((!feof(V->f)) && (*V->sbuffer != ';')) 
    GetLine(V);
}


static int 
endGB(char *s, int *addend)
{
  *addend = 0;
  return ((strstr(s,"//") != NULL) || (strstr(s,"LOCUS") == s));
}

static void 
readGenBank(struct ReadSeqVars *V)
{
  char *sptr;
  int   in_definition;

  while (strncmp(V->sbuffer, "LOCUS", 5) != 0)
    GetLine(V);

  if ((sptr = strtok(V->sbuffer+12, "\n\t ")) != NULL)
    {
      SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);
      SetSeqinfoString(V->sqinfo, sptr, SQINFO_ID);
    }

  in_definition = FALSE;
  while (! feof(V->f))
    {
      GetLine(V);
      if (! feof(V->f) && strstr(V->sbuffer, "DEFINITION") == V->sbuffer)
	{
	  if ((sptr = strtok(V->sbuffer+12, "\n")) != NULL)
	    SetSeqinfoString(V->sqinfo, sptr, SQINFO_DESC);
	  in_definition = TRUE;
	}
      else if (! feof(V->f) && strstr(V->sbuffer, "ACCESSION") == V->sbuffer)
	{
	  if ((sptr = strtok(V->sbuffer+12, "\n\t ")) != NULL)
	    SetSeqinfoString(V->sqinfo, sptr, SQINFO_ACC);
	  in_definition = FALSE;
	}
      else if (strncmp(V->sbuffer,"ORIGIN", 6) != 0)
	{
	  if (in_definition)
	    SetSeqinfoString(V->sqinfo, V->sbuffer, SQINFO_DESC);
	}
      else
	break;
    }

  readLoop(0, endGB, V);

  /* reading a real GenBank database file, we keep the source coords
   */
  V->sqinfo->start = 1;
  V->sqinfo->stop  = V->seqlen;
  V->sqinfo->olen  = V->seqlen;
  V->sqinfo->flags |= SQINFO_START | SQINFO_STOP | SQINFO_OLEN;


  while (!(feof(V->f) || ((*V->sbuffer!=0) && (strstr(V->sbuffer,"LOCUS") == V->sbuffer))))
    GetLine(V);
				/* SRE: V->s now holds "//", so sequential
				   reads are wedged: fixed Tue Jul 13 1993 */
  while (!feof(V->f) && strstr(V->sbuffer, "LOCUS  ") != V->sbuffer)
    GetLine(V);
}

static int
endNBRF(char *s, int *addend)
{
  char  *a;

  if ((a = strchr(s, '*')) != NULL) { /* end of 1st seq */
    /* "*" can be valid base symbol, drop it here */
    *a = 0;
    *addend = 1;
    return(1);
    }
  else if (*s == '>') { /* start of next seq */
    *addend = 0;
    return(1);
    }
  else
    return(0);
}

static void 
readNBRF(struct ReadSeqVars *V)
{
  char *sptr;

  if ((sptr = strtok(V->sbuffer+4, "\n\t ")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);

  GetLine(V);   /*skip title-junk line*/

  readLoop(0, endNBRF, V);

  while (!(feof(V->f) || (*V->sbuffer != 0 && *V->sbuffer == '>')))
    GetLine(V);
}


static int
endGCGdata(char *s, int *addend) 
{
  *addend = 0;
  return (*s == '>');
}

static void
readGCGdata(struct ReadSeqVars *V)
{
  char *sptr, *lptr;
  int   binary = FALSE;		/* whether data are binary or not */
  int   blen;			/* length of binary sequence */
  
				/* first line contains ">>>>" followed by name */
  if (Strparse(">>>>([^ ]+) .+2BIT +Len: ([0-9]+)", 
	       V->sbuffer, NULL, 2, &sptr, &lptr) == 0) {
    binary = TRUE;
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);
    blen = atoi(lptr);
    free(sptr);
    free(lptr);
  } else if (Strparse(">>>>([^ ]+) .+ASCII +Len: [0-9]+",
		      V->sbuffer, NULL, 1, &sptr) == 0) {  
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);
    free(sptr);
  } else Die("bogus GCGdata format? %s", V->sbuffer);

				/* second line contains free text description */
  GetLine(V);
  SetSeqinfoString(V->sqinfo, V->sbuffer, SQINFO_DESC);

  if (binary) {
    /* allocate for blen characters +3... (allow for 3 bytes of slop) */
    if (blen >= V->maxseq) {
      V->maxseq = blen;
      if ((V->seq = (char *) realloc (V->seq, sizeof(char)*(V->maxseq+4)))==NULL)
	Die("malloc failed");
    }
				/* read (blen+3)/4 bytes from file */
    if (fread(V->seq, sizeof(char), (blen+3)/4, V->f) < ((blen+3)/4))
      Die("fread failed");
    V->seqlen = blen;
				/* convert binary code to seq */
    GCGBinaryToSequence(V->seq, blen);
  }
  else readLoop(0, endGCGdata, V);
  
  while (!(feof(V->f) || ((*V->sbuffer != 0) && (*V->sbuffer == '>'))))
    GetLine(V);
}

static int
endPearson(char *s, int *addend)
{
  *addend = 0;
  return(*s == '>');
}

static void 
readPearson(struct ReadSeqVars *V)
{
  char *sptr;

				/* check for my special FASTA format */
  if (strstr(V->sbuffer, "..") != NULL && strstr(V->sbuffer, "::") != NULL)
    {
      if ((sptr = strtok(V->sbuffer+1, "\n\t ")) != NULL)
	SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);
      if ((sptr = strtok(NULL, "\n\t ")) != NULL)
	SetSeqinfoString(V->sqinfo, sptr, SQINFO_WGT);
      if ((sptr = strtok(NULL, "\n\t ")) != NULL)
	SetSeqinfoString(V->sqinfo, sptr, SQINFO_ID);
      if ((sptr = strtok(NULL, "\n\t ")) != NULL)
	SetSeqinfoString(V->sqinfo, sptr, SQINFO_ACC);
      if ((sptr = strtok(NULL, ".")) != NULL)
	SetSeqinfoString(V->sqinfo, sptr, SQINFO_START);
      if ((sptr = strtok(NULL, ".:")) != NULL)
	SetSeqinfoString(V->sqinfo, sptr, SQINFO_STOP);
      if ((sptr = strtok(NULL, ":\t\n ")) != NULL)
	SetSeqinfoString(V->sqinfo, sptr, SQINFO_OLEN);
      if ((sptr = strtok(NULL, "\n")) != NULL)
	SetSeqinfoString(V->sqinfo, sptr, SQINFO_DESC);
    }
  else			/* else we're normal FASTA format */
    {
      if ((sptr = strtok(V->sbuffer+1, "\n\t ")) != NULL)
	SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);
      if ((sptr = strtok(NULL, "\n")) != NULL)
	SetSeqinfoString(V->sqinfo, sptr, SQINFO_DESC);
    }

    readLoop(0, endPearson, V);

  while (!(feof(V->f) || ((*V->sbuffer != 0) && (*V->sbuffer == '>'))))
    GetLine(V);
}


static int
endEMBL(char *s, int *addend)
{
  *addend = 0;
  /* Some people (Berlin 5S rRNA database, f'r instance) use
   * an extended EMBL format that attaches extra data after
   * the sequence -- watch out for that. We use the fact that
   * real EMBL sequence lines begin with five spaces.
   * 
   * We can use this as the sole end test because readEMBL() will
   * advance to the next ID line before starting to read again.
   */
  return (strncmp(s,"     ",5) != 0);
/*  return ((strstr(s,"//") != NULL) || (strstr(s,"ID   ") == s)); */
}

static void 
readEMBL(struct ReadSeqVars *V)
{
  char *sptr;

				/* make sure we have first line */
  while (!feof(V->f) && strncmp(V->sbuffer, "ID  ", 4) != 0)
    GetLine(V);

  if ((sptr = strtok(V->sbuffer+5, "\n\t ")) != NULL)
    {
      SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);
      SetSeqinfoString(V->sqinfo, sptr, SQINFO_ID);
    }

  do {
    GetLine(V);
    if (!feof(V->f) && strstr(V->sbuffer, "AC  ") == V->sbuffer)
      {
	if ((sptr = strtok(V->sbuffer+5, ";  \t\n")) != NULL)
	  SetSeqinfoString(V->sqinfo, sptr, SQINFO_ACC);
      }
    else if (!feof(V->f) && strstr(V->sbuffer, "DE  ") == V->sbuffer)
      {
	if ((sptr = strtok(V->sbuffer+5, "\n")) != NULL)
	  SetSeqinfoString(V->sqinfo, sptr, SQINFO_DESC);
      }
  } while (! feof(V->f) && strncmp(V->sbuffer,"SQ",2) != 0);
  
  readLoop(0, endEMBL, V);

  /* reading a real EMBL database file, we keep the source coords
   */
  V->sqinfo->start = 1;
  V->sqinfo->stop  = V->seqlen;
  V->sqinfo->olen  = V->seqlen;
  V->sqinfo->flags |= SQINFO_START | SQINFO_STOP | SQINFO_OLEN;

				/* load next record's ID line */
  while (!feof(V->f) && strncmp(V->sbuffer, "ID  ", 4) != 0)
    GetLine(V);
}


static int
endZuker(char *s, int *addend)
{
  *addend = 0;
  return( *s == '(' );
}

static void
readZuker(struct ReadSeqVars *V)
{
  char *sptr;

  GetLine(V);  /*s == "seqLen seqid string..."*/

  if ((sptr = strtok(V->sbuffer+6, " \t\n")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);

  if ((sptr = strtok(NULL, "\n")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_DESC);

  readLoop(0, endZuker, V);

  while (!(feof(V->f) | ((*V->sbuffer != '\0') & (*V->sbuffer == '('))))
    GetLine(V);
}

static void 
readUWGCG(struct ReadSeqVars *V)
{
  char  *si;
  char  *sptr;
  int    done;

  V->seqlen = 0;

  /*writeseq: "    %s  Length: %d  (today)  Check: %d  ..\n" */
  /*drop above or ".." from id*/
  if ((si = strstr(V->sbuffer,"  Length: ")) != NULL) *si = 0;
  else if ((si = strstr(V->sbuffer,"..")) != NULL)    *si = 0;

  if ((sptr = strtok(V->sbuffer, "\n\t ")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);

  do {
    done = feof(V->f);
    GetLine(V);
    if (! done) addseq(V->sbuffer, V);
  } while (!done);
}

    
static void
readSquid(struct ReadSeqVars *V)
{
  char *sptr;
  int   dostruc = FALSE;

  while (strncmp(V->sbuffer, "NAM ", 4) != 0) GetLine(V);

  if ((sptr = strtok(V->sbuffer+4, "\n\t ")) != NULL)
    SetSeqinfoString(V->sqinfo, sptr, SQINFO_NAME);

  /*CONSTCOND*/
  while (1)
    {
      GetLine(V);
      if (feof(V->f)) {squid_errno = SQERR_FORMAT; return; }

      if (strncmp(V->sbuffer, "SRC ", 4) == 0)
	{
	  if ((sptr = strtok(V->sbuffer+4, " \t\n")) != NULL)
	    SetSeqinfoString(V->sqinfo, sptr, SQINFO_ID);
	  if ((sptr = strtok(NULL, " \t\n")) != NULL)
	    SetSeqinfoString(V->sqinfo, sptr, SQINFO_ACC);
	  if ((sptr = strtok(NULL, ".")) != NULL)
	    SetSeqinfoString(V->sqinfo, sptr, SQINFO_START);
	  if ((sptr = strtok(NULL, ".:")) != NULL)
	    SetSeqinfoString(V->sqinfo, sptr, SQINFO_STOP);
	  if ((sptr = strtok(NULL, ": \t\n")) != NULL)
	    SetSeqinfoString(V->sqinfo, sptr, SQINFO_OLEN);
	}
      else if (strncmp(V->sbuffer, "DES ", 4) == 0)
	{
	  if ((sptr = strtok(V->sbuffer+4, "\n")) != NULL)
	    SetSeqinfoString(V->sqinfo, sptr, SQINFO_DESC);
	}
      else if (strncmp(V->sbuffer,"SEQ", 3) == 0)
	break;
  }
  
  if (strstr(V->sbuffer, "+SS") != NULL) dostruc = TRUE;

  V->seqlen = 0;
  /*CONSTCOND*/
  while (1)
    {
				/* sequence line */
      GetLine(V);
      if (feof(V->f) || strncmp(V->sbuffer, "++", 2) == 0) 
	break;
      addseq(V->sbuffer, V);
				/* structure line */
      if (dostruc)
	{
	  GetLine(V);
	  if (feof(V->f)) { squid_errno = SQERR_FORMAT; return; }
	  addstruc(V->sbuffer, V);
	}
    }


  while (!feof(V->f) && strncmp(V->sbuffer, "NAM ", 4) != 0)
    GetLine(V);
}


/* Function: SeqfileOpen()
 * 
 * Purpose : Open a sequence database file and prepare for reading
 *           sequentially.
 *           
 * Args:     filename - name of file to open
 *           format   - format of file
 *           env      - environment variable for path (e.g. BLASTDB)                     
 *
 *           Returns opened SQFILE ptr, or NULL on failure.
 */
SQFILE *
SeqfileOpen(char *filename, int format, char *env)
{
  SQFILE *dbfp;

  if ((dbfp = (SQFILE *) malloc (sizeof(SQFILE))) == NULL)
    { squid_errno = SQERR_MEM; return NULL; }

  if ((dbfp->f = fopen(filename, "r")) == NULL &&
      (dbfp->f = EnvFileOpen(filename, env)) == NULL)
    {  squid_errno = SQERR_NOFILE; return NULL; }

  /* If it's a SELEX- or MSF-formatted alignment file, we
   * hack a way to fake sequential access: read the whole
   * alignment into a static at once. Relies on the prediction
   * that these files will always be relatively small, unlike
   * GenBank or EMBL databases.
   */
  /* Note (Tue Aug 30 11:57:19 1994): This cause a bug:
   * ReadSeq() is not reentrant for alignment files.
   */
  if (format == kSelex || format == kMSF || format == kClustal)
    {
      if (! ReadAlignment(filename, format, &ali_aseqs, &ali_num, &ali_ainfo))
	return NULL;
      if (! DealignAseqs(ali_aseqs, ali_num, &ali_rseqs))
	return NULL;
      curridx = 0;
    }

  /* Load the first line.
   */
  GetLine(dbfp);

  return dbfp;
}

/* Function: SeqfilePosition()
 * 
 * Purpose:  Move to a particular offset in a seqfile.
 *           Will not work on interleaved files (SELEX, MSF).
 */
void
SeqfilePosition(SQFILE *sqfp, long offset)
{
  fseek(sqfp->f, offset, SEEK_SET);
  GetLine(sqfp);
}


void
SeqfileClose(SQFILE *sqfp)
{
  fclose(sqfp->f);
  free(sqfp);
}


/* Function: ReadSeq()
 * 
 * Purpose:  Read next sequence from an open database file.
 *           Return the sequence and associated info.
 *           
 * Args:     fp      - open sequence database file pointer          
 *           format  - format of the file (previously determined
 *                      by call to SeqfileFormat())
 *           ret_seq - RETURN: sequence
 *           sqinfo  - RETURN: filled in w/ other information  
 *           
 * Return:   1 on success, 0 on failure.
 *           ret_seq and some field of sqinfo are allocated here,
 *           The preferred call mechanism to properly free the memory is:
 *           
 *           SQINFO sqinfo;
 *           char  *seq;
 *           
 *           ReadSeq(fp, format, &seq, &sqinfo);
 *           ... do something...
 *           FreeSequence(seq, &sqinfo);
 */
int
ReadSeq(SQFILE *V, int format, char **ret_seq, SQINFO *sqinfo)
{
  int    gotuw;

  squid_errno = SQERR_OK;
  if (format < kMinFormat || format > kMaxFormat) 
    {
      squid_errno = SQERR_FORMAT;
      *ret_seq    = NULL;
      return 0;
    }

  /* Here's the hack for accessing sequences from
   * the multiple sequence alignment formats
   */
  if (format == kMSF || format == kSelex || format == kClustal)
    {

      if (curridx >= ali_num)
	{			/* none left; free static ptrs */
	  FreeAlignment(ali_aseqs, ali_num, &ali_ainfo);
	  free(ali_rseqs);      /* only free rseqs pointer array */
	  return 0;
	}
      SeqinfoCopy(sqinfo, &(ali_ainfo.sqinfo[curridx]));
      *ret_seq = ali_rseqs[curridx];
      curridx++;
      return 1;
    }
  else
    {
      if (feof(V->f)) return 0;

      V->seq           = (char*) calloc (kStartLength+1, sizeof(char));
      V->maxseq        = kStartLength;
      V->seqlen        = 0;
      V->sqinfo        = sqinfo;
      V->sqinfo->flags = 0;
      V->dash_equals_n = (format == kEMBL) ? TRUE : FALSE;

      switch (format) {
      case kIG      : readIG(V);      break;
      case kStrider : readStrider(V); break;
      case kGenBank : readGenBank(V); break;
      case kNBRF    : readNBRF(V);    break;
      case kPearson : readPearson(V); break;
      case kXPearson: readPearson(V); break;
      case kEMBL    : readEMBL(V);    break;
      case kZuker   : readZuker(V);   break;
      case kPIR     : readPIR(V);     break;
      case kSquid   : readSquid(V);   break;
      case kGCGdata : readGCGdata(V); break; 
	
      case kGCG:
	do {			/* skip leading comments on GCG file */
	  gotuw = (strstr(V->sbuffer,"..") != NULL);
	  if (gotuw) readUWGCG(V);
	  GetLine(V);
	} while (! feof(V->f));
	break;

      case kIdraw:   /* SRE: no attempt to read idraw postscript */
      default:
	squid_errno = SQERR_FORMAT;
	free(V->seq);
	return 0;
      }
      V->seq[V->seqlen] = 0; /* stick a string terminator on it */
    }

  /* Cleanup
   */
  sqinfo->len    = V->seqlen; 
  sqinfo->flags |= SQINFO_LEN;
  *ret_seq = V->seq;
  if (squid_errno == SQERR_OK) return 1; else return 0;
}  

/* Function: SeqfileFormat()
 * 
 * Purpose:  Determine format of seqfile, and return it
 *           through ret_format. From Gilbert's seqFileFormat().
 *
 * Args:     filename   - name of sequence file      
 *           ret_format - RETURN: format code for file, see squid.h 
 *                        for codes.
 *           env        - name of environment variable containing
 *                        a directory path that filename might also be
 *                        found in. "BLASTDB", for example. Can be NULL.
 *           
 * Return:   1 on success, 0 on failure.
 */          
int
SeqfileFormat(char *filename, int  *ret_format, char *env)
{
  int   foundIG      = 0;
  int   foundStrider = 0;
  int   foundGB      = 0; 
  int   foundEMBL    = 0; 
  int   foundNBRF    = 0;
  int   foundPearson = 0;
  int   foundXPearson= 0;
  int   foundZuker   = 0;
  int   gotGCGdata   = 0;
  int   gotPIR       = 0;
  int   gotSquid     = 0;
  int   gotuw        = 0;
  int   gotMSF       = 0;
  int   gotClustal   = 0;
  int   done         = 0;
  int   format       = kUnknown;
  int   nlines= 0, dnalines= 0;
  int   splen = 0;
  char  sp[LINEBUFLEN];
  FILE *fseq;

#define ReadOneLine(sp)   \
  { done |= (feof(fseq)); \
    readline( fseq, sp);  \
    if (!done) { splen = (int) strlen(sp); ++nlines; } }

  if ((fseq = fopen(filename, "r")) == NULL &&
      (fseq = EnvFileOpen(filename, env)) == NULL)
    { squid_errno = SQERR_NOFILE;  return 0; }

  /* Look at a line at a time
   */
  while ( !done ) {
    ReadOneLine(sp);

    if (sp==NULL || *sp=='\0')
      /*EMPTY*/ ; 

    /* high probability identities: */
    
    else if (strstr(sp, " MSF:")   != NULL &&
	     strstr(sp, " Type:")  != NULL &&
	     strstr(sp, " Check:") != NULL)
      gotMSF = 1;

    else if (strncmp(sp, "CLUSTAL ", 8) == 0 && 
	     strstr( sp, "multiple sequence alignment"))
      gotClustal = 1;

    else if (strstr(sp," Check: ") != NULL)
      gotuw= 1;

    else if (strncmp(sp, "///", 3) == 0 || strncmp(sp, "ENTRY ", 6) == 0)
      gotPIR = 1;

    else if (strncmp(sp, "++", 2) == 0 || strncmp(sp, "NAM ", 4) == 0)
      gotSquid = 1;

    else if (strncmp(sp, ">>>>", 4) == 0 && strstr(sp, "Len: "))
      gotGCGdata = 1;

    /* uncertain identities: */

    else if (*sp ==';') {
      if (strstr(sp,"Strider") !=NULL) foundStrider= 1;
      else foundIG= 1;
    }
    else if (strncmp(sp,"LOCUS",5) == 0 || strncmp(sp,"ORIGIN",5) == 0)
      foundGB= 1;

    else if (*sp == '>') {
      if (sp[3] == ';')                              foundNBRF     = 1;
      else if (strstr(sp, "::") && strstr(sp, "..")) foundXPearson = 1;
      else                                           foundPearson  = 1;
    }

    else if (strstr(sp,"ID   ") == sp || strstr(sp,"SQ   ") == sp)
      foundEMBL= 1;

    else if (*sp == '(')
      foundZuker= 1;

    else {
      switch (Seqtype( sp )) {
      case kDNA:
      case kRNA: if (splen>20) dnalines++; break;
      default:   break;
      }
    }

    if      (gotMSF)     {format = kMSF;     done = 1; }
    else if (gotClustal) {format = kClustal; done = 1; }
    else if (gotSquid)   {format = kSquid;   done = 1; }
    else if (gotPIR)     {format = kPIR;     done = 1; }
    else if (gotGCGdata) {format = kGCGdata; done = 1; }
    else if (gotuw)  
      {
	if (foundIG) format= kIG;  /* a TOIG file from GCG for certain */
	else format= kGCG;
	done= 1;
      }
    else if ((dnalines > 1) || done || (nlines > 500)) {
      /* decide on most likely format */
      /* multichar idents: */
      if (foundStrider)      format= kStrider;
      else if (foundGB)      format= kGenBank;
      else if (foundEMBL)    format= kEMBL;
      else if (foundNBRF)    format= kNBRF;
      /* single char idents: */
      else if (foundIG)      format= kIG;
      else if (foundPearson) format= kPearson;
      else if (foundXPearson)format= kXPearson;
      else if (foundZuker)   format= kZuker;
      /* spacing ident: */
      else if (IsSELEXFormat(filename)) format= kSelex;
      /* no format chars: */
      else 
	{
	  squid_errno = SQERR_FORMAT;
	  return 0;
	}

      done= 1;
    }
  }

  if (fseq!=NULL) fclose(fseq);

  *ret_format = format;
  return 1;
#undef  ReadOneLine
}

/* Function: GCGBinaryToSequence()
 * 
 * Purpose:  Convert a GCG 2BIT binary string to DNA sequence.
 *           0 = C  1 = T  2 = A  3 = G
 *           4 nts/byte
 *           
 * Args:     seq  - binary sequence. Converted in place to DNA.
 *           len  - length of DNA. binary is (len+3)/4 bytes
 */
int
GCGBinaryToSequence(char *seq, int len)
{
  int   bpos;			/* position in binary   */
  int   spos;			/* position in sequence */
  char  twobit;
  int   i;

  for (bpos = (len-1)/4; bpos >= 0; bpos--) 
    {
      twobit = seq[bpos];
      spos   = bpos*4;

      for (i = 3; i >= 0; i--) 
	{
	  switch (twobit & 0x3) {
	  case 0: seq[spos+i] = 'C'; break;
	  case 1: seq[spos+i] = 'T'; break;
	  case 2: seq[spos+i] = 'A'; break;
	  case 3: seq[spos+i] = 'G'; break;
	  }
	  twobit = twobit >> 2;
	}
    }
  seq[len] = '\0';
  return 1;
}


int
GCGchecksum(char *seq, int   seqlen)
{
  int  check = 0, count = 0, i;

  for (i = 0; i < seqlen; i++) {
    count++;
    check += count * sre_toupper((int) seq[i]);
    if (count == 57) count = 0;
    }
  return (check % 10000);
}


/* Function: GCGMultchecksum()
 * 
 * Purpose:  Simple modification of GCGchecksum(),
 *           to create a checksum for multiple sequences.
 *           Gaps count.
 *           
 * Args:     seqs - sequences to be checksummed
 *           nseq - number of sequences
 *           
 * Return:   the checksum, a number between 0 and 9999
 */                      
int
GCGMultchecksum(char **seqs, int nseq)
{
  int check = 0;
  int count = 0;
  int idx;
  char *sptr;

  for (idx = 0; idx < nseq; idx++)
    for (sptr = seqs[idx]; *sptr; sptr++)
      {
	count++;
	check += count * sre_toupper((int) *sptr);
	if (count == 57) count = 0;
      }
  return (check % 10000);
}




/* Function: Seqtype()
 * 
 * Purpose:  Returns a (very good) guess about type of sequence:
 *           kDNA, kRNA, kAmino, or kOtherSeq.
 *           
 *           Modified from, and replaces, Gilbert getseqtype().
 */
int
Seqtype(char *seq)
{
  int  saw;			/* how many non-gap characters I saw */
  char c;
  int  po = 0;			/* count of protein-only */
  int  nt = 0;			/* count of t's */
  int  nu = 0;			/* count of u's */
  int  na = 0;			/* count of nucleotides */
  int  aa = 0;			/* count of amino acids */
  int  no = 0;			/* count of others */
  
  /* Look at the first 300 non-gap characters
   */
  for (saw = 0; *seq != '\0' && saw < 300; seq++)
    {
      c = sre_toupper((int) *seq);
      if (! isgap(c)) 
	{
	  if (strchr(protonly, c)) po++;
	  else if (strchr(primenuc,c)) {
	    na++;
	    if (c == 'T') nt++;
	    else if (c == 'U') nu++;
	  }
	  else if (strchr(aminos,c)) aa++;
	  else if (isalpha(c)) no++;
	  saw++;
	}
    }

  if (no > 0) return kOtherSeq;
  else if (po > 0) return kAmino;
  else if (na > aa) {
    if (nu > nt) return kRNA;
    else return kDNA;
    }
  else return kAmino;
}

int
WriteSeq(FILE *outf, int outform, char *seq, SQINFO *sqinfo)
{
  int   numline = 0;
  int   lines = 0, spacer = 0, width  = 50, tab = 0;
  int   i, j, l, l1, ibase;
  char  endstr[10];
  char  s[100];			/* buffer for sequence  */
  char  ss[100];		/* buffer for structure */
  int   checksum = 0;
  int   seqlen;   
  int   which_case;    /* 0 = do nothing. 1 = upper case. 2 = lower case */
  int   dostruc;		/* TRUE to print structure lines*/

  which_case = 0;
  dostruc    = FALSE;		
  seqlen     = (sqinfo->flags & SQINFO_LEN) ? sqinfo->len : strlen(seq);

				/* intercept Selex-format requests - SRE */
  if (outform == kSelex) {
    fprintf(outf, "%10s %s\n", sqinfo->name, seq);
    return 1;
  }

  if (outform == kClustal || outform == kMSF) {
    Warn("Tried to write Clustal or MSF with WriteSeq() -- bad, bad.");
    return 1;
  }

  strcpy( endstr,"");
  l1 = 0;

  /* 10Nov91: write this out in all possible formats: */
  checksum = GCGchecksum(seq, seqlen);

  switch (outform) {

    case kUnknown:    /* no header, just sequence */
      strcpy(endstr,"\n"); /* end w/ extra blank line */
      break;

    case kGenBank:
      fprintf(outf,"LOCUS       %s       %d bp\n", 
	      (sqinfo->flags & SQINFO_ID) ? sqinfo->id : sqinfo->name,
	      seqlen);
      fprintf(outf,"DEFINITION  %s\n", 
	      (sqinfo->flags & SQINFO_DESC) ? sqinfo->desc : "-");
      fprintf(outf,"ACCESSION   %s\n", 
	      (sqinfo->flags & SQINFO_ACC) ? sqinfo->acc : "-");
      fprintf(outf,"ORIGIN      \n");
      spacer = 11;
      numline = 1;
      strcpy(endstr, "\n//");
      break;

    case kGCGdata:
      fprintf(outf, ">>>>%s  9/95  ASCII  Len: %d\n", sqinfo->name, seqlen);
      fprintf(outf, "%s\n", (sqinfo->flags & SQINFO_DESC) ? sqinfo->desc : "-");
      break;

    case kNBRF:
      if (Seqtype(seq) == kAmino)
        fprintf(outf,">P1;%s\n", sqinfo->name);
      else
        fprintf(outf,">DL;%s\n", sqinfo->name);

      fprintf(outf,"%s %s\n", 
	      sqinfo->name,
	      (sqinfo->flags & SQINFO_DESC) ? sqinfo->desc : "");
      spacer = 11;
      strcpy(endstr,"*");
      break;

    case kPIR:
      fprintf(outf, "ENTRY          %s\n", 
	      (sqinfo->flags & SQINFO_ID) ? sqinfo->id : sqinfo->name);
      fprintf(outf, "TITLE          %s\n", 
	      (sqinfo->flags & SQINFO_DESC) ? sqinfo->desc : "-");
      fprintf(outf, "ACCESSION      %s\n",
	      (sqinfo->flags & SQINFO_ACC) ? sqinfo->acc : "-");
      fprintf(outf, "SUMMARY                                #Length %d  #Checksum  %d\n",
	      sqinfo->len, checksum);
      fprintf(outf, "SEQUENCE\n");
      fprintf(outf, "                  5        10        15        20        25        30\n");
      spacer  = 2;		/* spaces after every residue */
      numline = 1;              /* number lines w/ coords     */
      width   = 30;             /* 30 aa per line             */
      strcpy(endstr, "\n///");
      break;

    case kSquid:
      fprintf(outf, "NAM  %s\n", sqinfo->name);
      if (sqinfo->flags & (SQINFO_ID | SQINFO_ACC | SQINFO_START | SQINFO_STOP | SQINFO_OLEN))
	fprintf(outf, "SRC  %s %s %d..%d::%d\n",
		(sqinfo->flags & SQINFO_ID)    ? sqinfo->id     : "-",
		(sqinfo->flags & SQINFO_ACC)   ? sqinfo->acc    : "-",
		(sqinfo->flags & SQINFO_START) ? sqinfo->start  : 0,
		(sqinfo->flags & SQINFO_STOP)  ? sqinfo->stop   : 0,
		(sqinfo->flags & SQINFO_OLEN)  ? sqinfo->olen   : 0);
      if (sqinfo->flags & SQINFO_DESC)
	fprintf(outf, "DES  %s\n", sqinfo->desc);
      if (sqinfo->flags & SQINFO_SS)
	{
	  fprintf(outf, "SEQ  +SS\n");
	  dostruc = TRUE;	/* print structure lines too */
	}
      else
	fprintf(outf, "SEQ\n");
      numline = 1;                /* number seq lines w/ coords  */
      strcpy(endstr, "\n++");
      break;

    case kEMBL:
      fprintf(outf,"ID   %s\n",
	      (sqinfo->flags & SQINFO_ID) ? sqinfo->id : sqinfo->name);
      fprintf(outf,"AC   %s\n",
	      (sqinfo->flags & SQINFO_ACC) ? sqinfo->acc : "-");
      fprintf(outf,"DE   %s\n", 
	      (sqinfo->flags & SQINFO_DESC) ? sqinfo->desc : "-");
      fprintf(outf,"SQ             %d BP\n", seqlen);
      strcpy(endstr, "\n//"); /* 11Oct90: bug fix*/
      tab = 5;     /** added 31jan91 */
      spacer = 11; /** added 31jan91 */
      break;

    case kGCG:
      fprintf(outf,"%s\n", sqinfo->name);
      if (sqinfo->flags & SQINFO_ACC)
	fprintf(outf,"ACCESSION   %s\n", sqinfo->acc); 
      if (sqinfo->flags & SQINFO_DESC)
	fprintf(outf,"DEFINITION  %s\n", sqinfo->desc);
      fprintf(outf,"    %s  Length: %d  (today)  Check: %d  ..\n", 
	      sqinfo->name, seqlen, checksum);
      spacer = 11;
      numline = 1;
      strcpy(endstr, "\n");  /* this is insurance to help prevent misreads at eof */
      break;

    case kStrider: /* ?? map ?*/
      fprintf(outf,"; ### from DNA Strider ;-)\n");
      fprintf(outf,"; DNA sequence  %s, %d bases, %d checksum.\n;\n", 
	      sqinfo->name, seqlen, checksum);
      strcpy(endstr, "\n//");
      break;

			/* SRE: Don had Zuker default to Pearson, which is not
			   intuitive or helpful, since Zuker's MFOLD can't read
			   Pearson format. More useful to use kIG */
    case kZuker:
      which_case = 1;			/* MFOLD requires upper case. */
      /*FALLTHRU*/
    case kIG:
      fprintf(outf,";%s %s\n", 
	      sqinfo->name,
	      (sqinfo->flags & SQINFO_DESC) ? sqinfo->desc : "");
      fprintf(outf,"%s\n", sqinfo->name);
      strcpy(endstr,"1"); /* == linear dna */
      break;

    case kRaw:			/* Raw: just print the whole sequence. */
      fprintf(outf, "%s\n", seq);
      return 1;

    case kXPearson:
      if (sqinfo->flags & SQINFO_WGT)
	fprintf(outf, "> %s %f ", sqinfo->name, sqinfo->weight);
      else
	fprintf(outf, "> %s - ", sqinfo->name);

      fprintf(outf, "%s %s %d..%d::%d %s\n",
	      (sqinfo->flags & SQINFO_ID)    ? sqinfo->id     : "-",
	      (sqinfo->flags & SQINFO_ACC)   ? sqinfo->acc    : "-",
	      (sqinfo->flags & SQINFO_START) ? sqinfo->start  : 0,
	      (sqinfo->flags & SQINFO_STOP)  ? sqinfo->stop   : 0,
	      (sqinfo->flags & SQINFO_OLEN)  ? sqinfo->olen   : 0,
	      (sqinfo->flags & SQINFO_DESC)  ? sqinfo->desc   : "-");
      break;

    default :
    case kPearson:
      fprintf(outf,">%s  %s\n", sqinfo->name,
	      (sqinfo->flags & SQINFO_DESC)  ? sqinfo->desc   : "");
      break;
    }

  if (which_case == 1) s2upper(seq);
  if (which_case == 2) s2lower(seq);


  width = MIN(width,100);
  for (i=0, l=0, ibase = 1, lines = 0; i < seqlen; ) {
    if (l1 < 0) l1 = 0;
    else if (l1 == 0) {
      if (numline) fprintf(outf,"%8d ",ibase);
      for (j=0; j<tab; j++) fputc(' ',outf);
      }
    if ((spacer != 0) && ((l+1) % spacer == 1)) 
      { s[l] = ' '; ss[l] = ' '; l++; }
    s[l]  = seq[i];
    ss[l] = (sqinfo->flags & SQINFO_SS) ? sqinfo->ss[i] : '.';
    l++; i++;
    l1++;                 /* don't count spaces for width*/
    if (l1 == width || i == seqlen) {
      s[l] = ss[l] = '\0';
      l = 0; l1 = 0;
      if (dostruc)
	{
	  fprintf(outf, "%s\n", s);
	  if (numline) fprintf(outf,"         ");
	  for (j=0; j<tab; j++) fputc(' ',outf);
	  if (i == seqlen) fprintf(outf,"%s%s\n",ss,endstr);
	  else fprintf(outf,"%s\n",ss);
	}
      else
	{
	  if (i == seqlen) fprintf(outf,"%s%s\n",s,endstr);
	  else fprintf(outf,"%s\n",s);
	}
      lines++;
      ibase = i+1;
      }
    }
  return lines;
} 


/* Function: ReadMultipleRseqs()
 * 
 * Purpose:  Open a data file and
 *           parse it into an array of rseqs (raw, unaligned
 *           sequences).
 * 
 *           Caller is responsible for free'ing memory allocated
 *           to ret_rseqs, ret_weights, and ret_names.
 *           
 *           Weights are currently only supported for MSF format.
 *           Sequences read from all other formats will be assigned
 *           weights of 1.0. If the caller isn't interested in
 *           weights, it passes NULL as ret_weights.
 * 
 * Returns 1 on success. Returns 0 on failure and sets
 * squid_errno to indicate the cause.
 */
int
ReadMultipleRseqs(char              *seqfile,
		  int                fformat,
		  char            ***ret_rseqs,
		  SQINFO **ret_sqinfo,
		  int               *ret_num)
{
  SQINFO *sqinfo;               /* array of sequence optional info         */
  SQFILE *dbfp;                 /* open ptr for sequential access of file  */
  char  **rseqs;                /* sequence array                          */
  char  **aseqs;                /* aligned sequences, if file is aligned   */
  struct  aliinfo_s ainfo;      /* alignment-associated information        */
  int     num;			/* number of seqs                          */
  int     numalloced;           /* num of seqs currently alloced for       */
  int     idx;

  if (fformat == kSelex || fformat == kMSF || fformat == kClustal)
    {
      if (! ReadAlignment(seqfile, fformat, &aseqs, &num, &ainfo)) return 0;
      if (! DealignAseqs(aseqs, num, &rseqs))             return 0;

      /* copy the sqinfo array
       */
      if ((sqinfo=(SQINFO *)malloc(sizeof(SQINFO)*num)) == NULL) 
	{ squid_errno = SQERR_MEM; return 0; }
      for (idx = 0; idx < num; idx++)
	SeqinfoCopy(&(sqinfo[idx]), &(ainfo.sqinfo[idx]));
      FreeAlignment(aseqs, num, &ainfo);
    }
  else
    {
				/* initial alloc */
      num        = 0;
      numalloced = 16;
      if ((rseqs  = (char **) malloc (numalloced * sizeof(char *))) == NULL ||
	  (sqinfo = (SQINFO *) malloc 
	   (numalloced * sizeof(SQINFO))) == NULL)
	{ squid_errno = SQERR_MEM; return 0; }
      if ((dbfp = SeqfileOpen(seqfile, fformat, NULL)) == NULL) return 0;      

      while (ReadSeq(dbfp, fformat, &rseqs[num], &(sqinfo[num])))
	{
	  num++;
	  if (num == numalloced) /* more seqs coming, alloc more room */
	    {
	      numalloced += 16;
	      if ((rseqs  = (char **) realloc 
		   (rseqs, numalloced * sizeof(char *))) == NULL ||
		  (sqinfo = (SQINFO *) realloc 
		   (sqinfo, numalloced * sizeof(SQINFO))) == NULL)
		{ squid_errno = SQERR_MEM; return 0; }
	    }
	}
      SeqfileClose(dbfp);
    }

  *ret_rseqs  = rseqs;
  *ret_sqinfo = sqinfo;
  *ret_num    = num;
  return 1;
}




char *
SeqFormatString(int code)
{
  switch (code) {
  case kUnknown: return "Unknown";
  case kIG:      return "Intelligenetics";
  case kGenBank: return "GenBank";
  case kNBRF:    return "PIR-NBRF";
  case kEMBL:    return "EMBL";
  case kGCG:     return "GCG";
  case kGCGdata: return "GCG datalibrary";
  case kStrider: return "Strider";
  case kXPearson:return "Extended FASTA";
  case kPearson: return "FASTA";
  case kZuker:   return "Zuker";
  case kIdraw:   return "Idraw";
  case kSelex:   return "SELEX";
  case kMSF:     return "MSF";
  case kClustal: return "Clustal";
  case kPIR:     return "PIR-CODATA";
  case kRaw:     return "raw";
  case kSquid:   return "Squid";
  default:       return "(bad code)";
  }
}



/* Function: GSIOpen()
 * 
 * Purpose:  Open a GSI file. Returns the number of records in
 *           the file and a file pointer. Returns NULL on failure.
 *           The file pointer should be fclose()'d normally.
 */
GSIFILE *
GSIOpen(char *gsifile)
{
  GSIFILE *gsi;
  char     magic[GSI_KEYSIZE];


  gsi = (GSIFILE *) MallocOrDie (sizeof(GSIFILE));
  if ((gsi->gsifp = fopen(gsifile, "r")) == NULL)
    { squid_errno = SQERR_NOFILE; return NULL; }
  if (! fread(magic, GSI_KEYSIZE, 1, gsi->gsifp))
    { squid_errno = SQERR_NODATA; return NULL; }
  if (strcmp(magic, "GSI") != 0)
    { squid_errno = SQERR_FORMAT; return NULL; }
  if (! fread(&(gsi->nfiles), sizeof(short), 1, gsi->gsifp))
    { squid_errno = SQERR_NODATA; return NULL; }
  if (! fread(&(gsi->recnum), sizeof(long), 1, gsi->gsifp))
    { squid_errno = SQERR_NODATA; return NULL; }
  return gsi;
}

/* Function: GSIGetOffset()
 * 
 * Purpose:  From a key (sequence name), find a disk offset
 *           in an open general sequence index file by binary
 *           search. Presumably GSI indexing could be even faster
 *           if we used hashing.
 *   
 * Args:     gsi         - GSI index file, opened by GSIOpen()
 *           key         - name of key to retrieve indices for
 *           ret_seqfile - pre-alloced char[32] array for seqfile name
 *           ret_offset  - return: disk offset in seqfile.         
 */
int
GSIGetOffset(GSIFILE *gsi, char *key, char *ret_seqfile, long *ret_offset)
{
  long left, right, mid;
  int cmp;
  char name[GSI_KEYSIZE + 1];
  short filenum;
  name[GSI_KEYSIZE] = '\0';

  left  = gsi->nfiles + 1;
  right = gsi->nfiles + gsi->recnum;
  mid   = (left + right) / 2;
  fseek(gsi->gsifp, mid * GSI_RECSIZE, SEEK_SET);

  while (fread(name, GSI_KEYSIZE, 1, gsi->gsifp))
    {
      cmp = strcmp(name, key);
      if (cmp == 0) break;	/* name same as key: done! */
      else if (left >= right)
	return 0;		/* failed to find key: fail! */
      else if (cmp < 0)         
	left = mid + 1;		/* name greater than key */
      else if (cmp > 0)		
	right = mid - 1;	/* name less than key */
      mid = (left + right) / 2;
      fseek(gsi->gsifp, mid * GSI_RECSIZE, SEEK_SET);
    }

  /* Read the remainder of the record
   */
  if (! fread(&filenum, sizeof(short), 1, gsi->gsifp))
    { squid_errno = SQERR_NODATA; return 0; }
  if (! fread(ret_offset, sizeof(long), 1, gsi->gsifp))
    { squid_errno = SQERR_NODATA; return 0; }

  /* Look up the sequence name
   */
  fseek(gsi->gsifp, filenum * GSI_RECSIZE, SEEK_SET);
  if (! fread(ret_seqfile, GSI_KEYSIZE, 1, gsi->gsifp))
    { squid_errno = SQERR_NODATA; return 0; }
  
  return 1;
}
    
/* Function: GSIClose()
 * 
 * Purpose:  Close an open GSI sequence index file.
 */
void
GSIClose(GSIFILE *gsi)
{
  fclose(gsi->gsifp);
  free(gsi);
}


/* Function: EnvFileOpen()
 * Date:     Sun Feb 12 10:55:29 1995            
 * 
 * Purpose:  Open a file, given a file name and an environment
 *           variable that contains a directory path. Files
 *           are opened read-only. Uses UNIX getenv() and
 *           thus is UNIX-specific.
 *           
 *           For instance: 
 *             fp = EnvFileOpen("BLOSUM45", "BLASTMAT");
 *           or:
 *             fp = EnvFileOpen("swiss", "BLASTDB");  
 *             
 * Args:     fname - name of file to open
 *           env   - name of environment variable containing path
 *  
 * Return:   FILE * to open file, or NULL on failure -- same as fopen()
 */
FILE *
EnvFileOpen(char *fname, char *env)
{
  FILE *fp;
  char *path;
  char  full[512];

  if (env == NULL) return NULL;
  if ((path = getenv(env)) == NULL) return NULL;
  if (((int) strlen(fname) + (int) strlen(env) + 2) > 512) return NULL;
  sprintf(full, "%s/%s", path, fname);
  fp = fopen(full, "r");
  return fp;
}
