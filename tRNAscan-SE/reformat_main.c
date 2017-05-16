/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* reformat_main.c
 * Mon Sep 13 13:06:51 1993
 * 
 * reformat - reformat sequence files.
 */


#include <stdio.h>
#include <string.h>
#include "squid.h"

#define OPTIONS "dhlru"

char usage[]  = "Usage: reformat [-options] <format> <seqfile>\n\
  Convert between sequence file formats.\n\
  Available formats are:\n\
    embl\n\
    fasta\n\
    xfasta\n\
    genbank\n\
    gcg\n\
    gcgdata\n\
    msf\n\
    strider\n\
    zuker\n\
    ig\n\
    nbrf\n\
    pir\n\
    selex\n\
    squid\n\
    raw\n\n\
  Available options are:\n\
    -d  : force DNA alphabet for nucleic acid sequence\n\
    -r  : force RNA alphabet for nucleic acid sequence\n\
    -l  : force lower case\n\
    -u  : force upper case\n\
    -h  : print short help and usage info\n";

struct seqfmt_s {  char *formatname; int fmt; } seqfmt[] =
{
  { "embl",    kEMBL    },
  { "fasta",   kPearson },
  { "xfasta",  kXPearson},
  { "genbank", kGenBank },
  { "gcg",     kGCG     },
  { "gcgdata", kGCGdata },
  { "msf",     kMSF     },
  { "strider", kStrider },
  { "zuker",   kZuker   },
  { "ig",      kIG      },
  { "nbrf",    kNBRF    },
  { "pir",     kPIR     },
  { "selex",   kSelex   },
  { "squid",   kSquid   },
  { "raw",     kRaw     },
};
#define NUMFORMATS  (sizeof(seqfmt) / sizeof(struct seqfmt_s))


int
main(int argc, char **argv)
{
  char     *seqfile;            /* name of sequence file */
  char     *format;
  SQFILE   *dbfp;		/* open sequence file */
  int       fmt;		/* format of seqfile  */
  int       outfmt;		/* output format */
  char     *seq;		/* sequence */
  SQINFO    sqinfo;
  int       i;

  int    force_rna;		/* TRUE to force RNA alphabet */
  int    force_dna;		/* TRUE to force DNA alphabet */
  int    force_lower;		/* TRUE to force lower case   */
  int    force_upper;		/* TRUE to force upper case   */

  int    optchar;		/* option character, command line */
  extern int   optind;
  
  /***********************************************
   * Parse command line
   ***********************************************/

  force_rna   = FALSE;
  force_dna   = FALSE;
  force_upper = FALSE;
  force_lower = FALSE;

  while ((optchar = getopt(argc, argv, OPTIONS)) != -1)
    switch (optchar) {

    case 'd': force_dna   = TRUE; break;
    case 'l': force_lower = TRUE; break;
    case 'r': force_rna   = TRUE; break;
    case 'u': force_upper = TRUE; break;

    case 'h': 
      printf("reformat %s, %s\n%s\n", squid_version, squid_date, usage);
      exit(EXIT_SUCCESS);
    default:
      Die("%s\n", usage);
    }

  if (argc - optind != 2)
    Die("%s\n", usage); 
  if (force_lower && force_upper)
    Die("Can't force both upper case and lower case. Stop trying to confuse me.\n%s", 
	usage);
  if (force_rna && force_dna)
    Die("Can't force both RNA and DNA. Stop trying to find bugs, you'll be sorry.\n%s", 
	usage);

  format  = argv[optind]; optind++;
  seqfile = argv[optind]; optind++;
  
  /***********************************************
   * Figure out what format we're supposed to write
   ***********************************************/

  outfmt = kUnknown;
  for (i = 0; i < NUMFORMATS; i++)
    if (strcasecmp(format, seqfmt[i].formatname) == 0)
      outfmt = seqfmt[i].fmt;
  
  if (outfmt == kUnknown)
    Die("Unknown output format %s\n%s", format, usage);

  /***********************************************
   * Reformat the file, printing to stdout.
   ***********************************************/

  if (! SeqfileFormat(seqfile, &fmt, NULL))
    Die("Can't determine format of file %s\n", seqfile);

  if ((fmt == kMSF || fmt == kSelex || fmt == kClustal) &&
      (outfmt == kMSF || outfmt == kSelex))
    {
      char **aseqs;
      int    num;
      AINFO  ainfo;

      ReadAlignment(seqfile, fmt, &aseqs, &num, &ainfo);

      for (i = 0; i < num; i++)
	{
	  if (force_dna)   ToDNA(aseqs[i]);
	  if (force_rna)   ToRNA(aseqs[i]);
	  if (force_lower) s2lower(aseqs[i]);
	  if (force_upper) s2upper(aseqs[i]);
	}
      
      switch (outfmt) {
      case kMSF:   WriteMSF(stdout, aseqs, num, &ainfo);       break;
      case kSelex: WriteSELEX(stdout, aseqs, num, &ainfo, 50); break;
      }
      FreeAlignment(aseqs, num, &ainfo);
    }
  else if (outfmt == kMSF || outfmt == kSelex)
    {
      Die("Sorry, you can't make alignment files from unaligned files");
    }
  else
    {
      if ((dbfp = SeqfileOpen(seqfile, fmt, NULL)) == NULL)
	Die("Failed to open sequence file %s for reading", seqfile);
  
      while (ReadSeq(dbfp, fmt, &seq, &sqinfo))
	{
	  if (force_dna)   ToDNA(seq);
	  if (force_rna)   ToRNA(seq);
	  if (force_lower) s2lower(seq);
	  if (force_upper) s2upper(seq);

	  WriteSeq(stdout, outfmt, seq, &sqinfo);
	  FreeSequence(seq, &sqinfo);
	}
      SeqfileClose(dbfp);
    }
  return 0;
}
