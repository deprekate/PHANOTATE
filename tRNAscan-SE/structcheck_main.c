/* structcheck_main.c
 * SRE, Mon Dec 20 07:48:07 1993
 * 
 * Check a set of individual RNA structures for non-Watson-Crick
 * or GU base pairs. Keep statistics on the number of such "errors" found 
 * overall and per sequence. Convert the offending base pairs to "*"
 * characters in the structure string and print out the structures
 * and sequences.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#ifdef NEED_GETOPTH
#include <getopt.h>
#endif

#include "structs.h"
#include "funcs.h"
#include "squid.h" 
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#define OPTIONS "hmo:psu"

static char usage[]  = "\
Usage: structcheck [-options] <SELEX RNA alignment file>\n\
where options are:\n\
   -h        : print short help and usage info\n\
   -m        : only check for more possible base pairs\n\
   -o <file> : save flagged structure annotation alignment to <file>\n\
   -p        : only check for non-Watson-Crick pair assignments\n\
   -s        : only check for agreement with consensus structure\n\
   -u        : only check for upper-case structured positions\n";

static char banner[] = "structcheck: check RNA secondary structures, flag questionables";

extern int VerifyKHS(char *ss);

int
main(int argc, char **argv)
{ 
  char     **aseqs;             /* RNA sequences           */
  AINFO      ainfo;             /* misc. associated alignment info */
  int        nseq;		/* number of seqs               */ 
  char      *seqfile;           /* sequence file                */
  int        idx;		/* index for sequences          */
  int        pos;		/* position index in a seq      */
  int       *ct;                /* CT0 representation of a structure */
  int        badseq;		/* number of bad structures     */
  int        badpairs;		/* total bad base pairs         */
  int        pairs;		/* total base pairs             */
  int        is_bad_seq;
  int        structure_agrees;  /* TRUE if consensus and secondary structure are same */
  int        nonconsensus;      /* count seqs w/ differing cons/indiv ss assignments */
  int        noncons_positions; /* count of pos w/ differing cons/indiv ss assignments */
  int        npos;		/* number of non-gap positions  */
  int        has_morepairs;     /* TRUE if more base pairs are possible than indicated in ss */
  int        morepairs;		/* how many more pairs should've been made in alignment */
  int        morepair_seqs;	/* how many seqs should've had more pairs made in them */
  int        badupper_seqs;
  int        badupper_bases;
  int        has_badupper;
  char      *ss;                /* aligned secondary structure string */
  
  char      *outfile;	        /* save flagged annotated alignment to */
  FILE      *ofp;		/* open outfile for writing            */
  int        check_pairs;	/* if TRUE, check that pairs are complementary */
  int        check_consensus;	/* if TRUE, compare indiv structs against consensus */
  int        check_morepairs;   /* if TRUE, check for more obvious base-pairing in ss */    
  int        check_isupper;     /* if TRUE, check structured pos's are upper, single-s is lower */

  int          optc;		/* for getopt() */
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

#ifdef MEMDEBUG			/* for Cahill's dbmalloc */
  unsigned long histid1, histid2, orig_size, current_size;
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/

  outfile         = NULL;
  check_pairs     = FALSE;
  check_consensus = FALSE;
  check_morepairs = FALSE;
  check_isupper   = FALSE;

  while ((optc = getopt(argc, argv, OPTIONS)) != -1)
    switch (optc) {


    case 'o': outfile         = optarg; break;
    case 'm': check_morepairs = TRUE;   break;
    case 'p': check_pairs     = TRUE;   break;
    case 's': check_consensus = TRUE;   break;
    case 'u': check_isupper   = TRUE;   break;
    case 'h': 
      printf("%s\n  version %s (%s)\n%s", banner, RELEASE, RELEASEDATE, usage);
      exit(0);

    default:
      Die("%s", usage);
    }

				/* Default behaviour if we're not doing a single thing */
  if (! check_morepairs && 
      ! check_pairs     &&
      ! check_consensus &&
      ! check_isupper)
    { 
      check_morepairs = TRUE;
      check_pairs     = TRUE;
      check_consensus = TRUE;
      check_isupper   = FALSE;	/* did this only for SRP-RNA alignment */
    }

  if (argc - optind != 1)
    Die("Wrong number of command line arguments.\n%s\n", usage);
  
  seqfile = argv[argc-1]; 

#ifdef MEMDEBUG
  orig_size = malloc_size(&histid1);
#endif

  /*********************************************** 
   * Get sequence data 
   ***********************************************/

				/* read the training seqs from file */
  if (! ReadSELEX(seqfile, &aseqs, &nseq, &ainfo))
    Die("Failed to read aligned sequence file %s", seqfile);
  
  for (idx = 0; idx < nseq; idx++)
    if (ainfo.sqinfo[idx].flags & SQINFO_SS) break;
  if (idx == nseq)
    Die("No secondary structure info in sequence file %s", seqfile);

  if ( check_consensus == TRUE && !(ainfo.flags & AINFO_CS))
    {
      check_consensus = FALSE;
      Warn("No consensus structure in %s; can't check against it", seqfile);
    }

  /*********************************************** 
   * Print banner
   ***********************************************/

  puts(banner);
  printf("     release %s, %s\n", RELEASE, RELEASEDATE);
  printf("---------------------------------------------------\n");
  puts("");

  
  badpairs = pairs = badseq = 0;
  nonconsensus = npos = noncons_positions = 0;
  morepairs = morepair_seqs = 0;
  badupper_seqs = badupper_bases = 0;
  for (idx = 0; idx < nseq; idx++)
    {
      if (! (ainfo.sqinfo[idx].flags & SQINFO_SS)) continue;

      /* Make an aligned secondary structure string for our tests
       */
      MakeAlignedString(aseqs[idx], ainfo.alen, ainfo.sqinfo[idx].ss, &ss);

      /* Get a ct structure to do the other tests
       */

      if (! KHS2ct(ss, ainfo.sqinfo[idx].len, FALSE, &ct))
	{
	  printf("sequence %-10.10s (#%d) has an improper secondary structure\n",
		 ainfo.sqinfo[idx].name, idx);
	  free(ct);
	  ct = NULL;
	  VerifyKHS(ss);
	}  
		 
      /* Check if structured positions are upper-case and unstructured 
       * positions are not. I used this for checking the SRP-RNA alignment
       * of Larsen and Zwieb, and it may come in handy later too.
       */
      if (check_isupper)
	{
	  has_badupper = FALSE;
	  for (pos = 0; pos < ainfo.alen; pos++)
	    {
	      if (isgap(aseqs[idx][pos])) continue;
	      if ( (ss[pos] == '.'  && ! islower(aseqs[idx][pos])) ||
		   (ss[pos] != '.'  && ! isupper(aseqs[idx][pos])) )
		{
		  has_badupper = TRUE;
		  badupper_bases++;
		  ss[pos] = '*';
		}
	    }
	  if (has_badupper)
	    badupper_seqs++;
	}

      /* Test for non-Watson-Crick base pairs
       */
      if (check_pairs && ct != NULL)
	{
	  is_bad_seq = False;
	  for (pos = 0; pos < ainfo.alen; pos++)
				/* second test makes sure we only look at
				 * each bp once */
	    if (ct[pos] != -1 && pos < ct[pos])
	      {
		pairs++;
		if (! IsRNAComplement(aseqs[idx][pos], aseqs[idx][ct[pos]], TRUE))
		  {
		    ss[pos] = ss[ct[pos]] = '*';
		    is_bad_seq = True;
		    badpairs++;
		  }
	      }
	  if (is_bad_seq) badseq++;
	}

      /* Test for disagreement with consensus structure
       */
      if (check_consensus)
	{
	  structure_agrees = True;
	  for (pos = 0; pos < ainfo.alen; pos++)
	    {
	      if (isgap(aseqs[idx][pos])) continue;
	      npos++;
	      if (ss[pos] != ainfo.cs[pos] )
		{
		  structure_agrees = False;
		  ss[pos] = '*';
		  noncons_positions++;
		}
	    }
	  if (! structure_agrees) nonconsensus++;
	}


      /* Test for whether more pairs are obviously possible in the structure.
       * *Very* crude. For each base pair i,j, if (i-1,j+1) or 
       * (i+1,j-1) are unpaired but complementary, flag them. This is *not*
       * a full-blown structure optimization algorithm (such a thing 
       * is possible, but would require dynamic programming), but it should
       * flag most suspicious spots.
       */
      if (check_morepairs && ct != NULL)
	{
	  has_morepairs = FALSE;
	  for (pos = 0; pos < ainfo.alen; pos++)
				/* second test makes sure we only look at
				 * each bp once */
	    if (ct[pos] != -1 && pos < ct[pos])
	      {
				/* check i-1,j+1 pair; careful of ends */
		if (pos > 0 && ct[pos] < ainfo.alen &&
		    ct[pos-1] == -1 && ct[ct[pos]+1] == -1 &&
		    IsRNAComplement(aseqs[idx][pos-1], aseqs[idx][ct[pos]+1], TRUE))
		  {
		    has_morepairs = TRUE;
		    morepairs++;
		    ss[pos-1] = ss[ct[pos]+1] = '*';
		  }
		    
				/* check i+1,j-1 pair; don't need to worry about ends */
		if (ct[pos+1] == -1 && ct[ct[pos]-1] == -1 &&
		    IsRNAComplement(aseqs[idx][pos+1], aseqs[idx][ct[pos]-1], TRUE))
		  {
		    has_morepairs = TRUE;
		    morepairs++;
		    ss[pos+1] = ss[ct[pos]-1] = '*';
		  }
	      }
	  if (has_morepairs) morepair_seqs++;
	}
      free(ct);

      /* Convert aligned ss back to dealigned ss
       */
      free(ainfo.sqinfo[idx].ss);
      MakeDealignedString(aseqs[idx], ainfo.alen, ss, &(ainfo.sqinfo[idx].ss));
      free(ss);
    }

  if (outfile == NULL)
    {
      if (! WriteSELEX(stdout, aseqs, nseq, &ainfo, 60))
	Die("Failed to write alignment to stdout");
    }
  else
    {
      if ((ofp = fopen(outfile, "w")) == NULL)
	Die("Failed to open flagged annotated alignment file %s", outfile);
      if (! WriteSELEX(ofp, aseqs, nseq, &ainfo, 60))
	Die("Failed to write alignment to %s", outfile);
      fclose(ofp);
      printf("Wrote flagged annotated alignment file to %s\n", outfile);
    }

  if (check_pairs)
    {
      printf("\nComplementarity check:\n");
      printf("%d/%d structures contain non-Watson-Crick, non-GU pairs\n", badseq, nseq);
      printf("%d/%d base pairs are questionable\n", badpairs, pairs);
    }
  if (check_consensus)
    {
      printf("\nConsensus structure check:\n");
      printf("%d/%d structures disagree with consensus\n", nonconsensus, nseq);
      printf("%d/%d non-gap sequence positions disagree\n", noncons_positions, npos);
    }
  if (check_morepairs)
    {
      printf("\nAdditional structure check:\n");
      printf("%d/%d structures have obvious additional pairings\n", morepair_seqs, nseq);
      printf("%d additional base pairs are predicted\n", morepairs);
    }
  if (check_isupper)
    {
      printf("\nCheck that structured positions are upper case:\n");
      printf("%d/%d structures have conflicts\n", badupper_seqs, nseq);
      printf("%d conflicts are detected\n", badupper_bases);
    }

  FreeAlignment(aseqs, nseq, &ainfo);

#ifdef MEMDEBUG
  current_size = malloc_size(&histid2);
  if (current_size != orig_size) malloc_list(2, histid1, histid2);
  else                           fprintf(stderr, "No memory leaks, sir.\n");
#endif

  return 0;
}


/* Function: VerifyKHS()
 * 
 * Purpose:  Examine a bad structure string and print out diagnostics 
 *           about it.
 *
 * Return:   1 if string is OK, 0 if string is bad.
 */
int
VerifyKHS(char *ss)
{
  int symcount[27];		/* 0 is normal pairs. 1-26 for pseudoknots */
  int i;
  int pos;
  int status = 1;

  for (i = 0; i < 27; i++)
    symcount[i] = 0;

  for (pos = 0; ss[pos] != '\0'; pos++)
    {
      if (ss[pos] > 127)	/* evade SGI ctype.h islower(), isupper() bug */
	{
	  status = 0;
	  fprintf(stderr, "  structure has garbage symbol (val %d) at position %d\n",
		  (int) ss[pos], pos);
	}
      else if (ss[pos] == '>')
	symcount[0] ++;
      else if (ss[pos] == '<')
	symcount[0] --;
      else if (isupper((int) ss[pos]))
	symcount[ss[pos] - 'A' + 1] ++;
      else if (islower((int) ss[pos]))
	symcount[ss[pos] - 'a' + 1] --;
      else if (ss[pos] != '.')
	{
	  status = 0;
	  fprintf(stderr, "  structure has invalid symbol %c at position %d\n",
		  ss[pos], pos);
	}
	  
    }
      
  if (symcount[0] != 0)
    {
      status = 0;
      fprintf(stderr, "  structure has extra paired bases: %d on %s\n",
	      abs(symcount[0]), 
	      (symcount[0] > 0) ? "left" : "right");
    }

  for (i = 1; i < 27; i++)
    if (symcount[i] != 0)
      {
	status = 0;
	fprintf(stderr, "  structure has extra paired bases for pseudoknot %c: %d on %s\n",
		(char) (i + 'A' - 1),
		abs(symcount[i]), 
		(symcount[i] > 0) ? "left" : "right");
      }
  return status;
}
  
