/* misc.c
 * Stuff with no obvious other place to go;
 * mostly alphabet-related functions.
 * SRE, Mon Sep  6 10:50:46 1993
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "squid.h"
#include "structs.h"
#include "funcs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

char *ALPHABET = RNA_ALPHABET;


/* Function: SymbolIndex()
 * 
 * Purpose:  Given a sequence symbol, return an array index. 
 *           Used for retrieving emission statistics from the
 *           model. This would be trivial, except that we must
 *           expect to see degenerate codes for both RNA and
 *           protein sequence. The rule we follow for degenerate
 *           codes is like BLAST -- we choose one possibility
 *           at random. Note that this is unlike the HMM package,
 *           which goes to the trouble of calculating a weighted
 *           average of the possibilities.
 *           
 *           If the symbol is not recognized, print a warning but
 *           treat it as a fully ambiguous position (N or X).
 *           
 * Return:   the index, usually 0..3 or 0..19.
 */
int
SymbolIndex(char  sym)
{
  char *sptr;
				/* trivial case: sym is in alphabet */
  if ((sptr = strchr(ALPHABET, sym)) != NULL)
    return (sptr - ALPHABET);

				/* non-trivial case: possible degenerate symbol */
  if (ALPHATYPE == kAmino)
    {
      switch (sym) {
      case 'B': return (strchr(ALPHABET, "ND"[CHOOSE(2)]) - ALPHABET);
      case 'Z': return (strchr(ALPHABET, "QE"[CHOOSE(2)]) - ALPHABET);
      default:
	Warn("Warning: unrecognized character %c in sequence\n", sym);
				/* break thru to case 'X' */
      case 'X':
	return(CHOOSE(20));
      }
    }

  else if (ALPHATYPE == kDNA)
    {
      switch (sym) {
      case 'B': return (strchr(ALPHABET, "CGT"[CHOOSE(3)]) - ALPHABET);
      case 'D': return (strchr(ALPHABET, "AGT"[CHOOSE(3)]) - ALPHABET);
      case 'H': return (strchr(ALPHABET, "ACT"[CHOOSE(3)]) - ALPHABET);
      case 'K': return (strchr(ALPHABET, "GT" [CHOOSE(2)]) - ALPHABET);
      case 'M': return (strchr(ALPHABET, "AC" [CHOOSE(2)]) - ALPHABET);
      case 'R': return (strchr(ALPHABET, "AG" [CHOOSE(2)]) - ALPHABET);
      case 'S': return (strchr(ALPHABET, "CG" [CHOOSE(2)]) - ALPHABET);
      case 'U': return (strchr(ALPHABET, 'T')                 - ALPHABET);
      case 'V': return (strchr(ALPHABET, "ACG"[CHOOSE(3)]) - ALPHABET);
      case 'W': return (strchr(ALPHABET, "AT" [CHOOSE(2)]) - ALPHABET);
      case 'Y': return (strchr(ALPHABET, "CT" [CHOOSE(2)]) - ALPHABET);

      case 'N': return (CHOOSE(4));
      case 'X': return (CHOOSE(4)); /* X is not IUPAC, but that doesn't stop 
					  biologists from using it. */
      default:	
	Warn("Warning: unrecognized character %c in sequence\n", sym);
	return (CHOOSE(4));
      }
    }

  else if (ALPHATYPE == kRNA)
    {
      switch (sym) {
      case 'B': return (strchr(ALPHABET, "CGU"[CHOOSE(3)]) - ALPHABET);
      case 'D': return (strchr(ALPHABET, "AGU"[CHOOSE(3)]) - ALPHABET);
      case 'H': return (strchr(ALPHABET, "ACU"[CHOOSE(3)]) - ALPHABET);
      case 'K': return (strchr(ALPHABET, "GU" [CHOOSE(2)]) - ALPHABET);
      case 'M': return (strchr(ALPHABET, "AC" [CHOOSE(2)]) - ALPHABET);
      case 'R': return (strchr(ALPHABET, "AG" [CHOOSE(2)]) - ALPHABET);
      case 'S': return (strchr(ALPHABET, "CG" [CHOOSE(2)]) - ALPHABET);
      case 'T': return (strchr(ALPHABET, 'U')                 - ALPHABET);
      case 'V': return (strchr(ALPHABET, "ACG"[CHOOSE(3)]) - ALPHABET);
      case 'W': return (strchr(ALPHABET, "AU" [CHOOSE(2)]) - ALPHABET);
      case 'Y': return (strchr(ALPHABET, "CU" [CHOOSE(2)]) - ALPHABET);

      case 'N': return (CHOOSE(4));
      case 'X': return (CHOOSE(4));  /* X is not IUPAC, but that doesn't stop 
					  biologists from using it. */
      default:	
	Warn("Warning: unrecognized character %c in sequence\n", sym);
	return (CHOOSE(4));
      }
    }
  return 0;			/* not reached */
}


/* Function: PrepareSequence()
 * 
 * Purpose:  Ran into a severe bug caused by degenerate symbols. Original
 *           strategy was to randomly assign a single symbol as we do
 *           Viterbi calculations, but since we don't keep traceback pointers
 *           when the trace tries to recalculate, it doesn't know the
 *           random choices made by VitFill(). 
 *           
 *           This is a fix, and it's a bit more extreme. Go through a
 *           sequence and *replace* degenerate symbols once
 *           and for all with single randomly chosen ones. Also,
 *           we convert to upper case ALPHATYPE alphabet.
 *           
 * Args:     seq - sequence to prepare.
 *           
 * Return:   1 on success, 0 on failure.
 */
int
PrepareSequence(char *seq)
{
  char *sym;

  for (sym = seq; *sym != '\0'; sym++)
    {
      *sym = toupper((int)*sym);
				/* sym is in alphabet, or a gap? ok, go to next one */
      if (strchr(ALPHABET, *sym) != NULL ||
	  isgap(*sym)) 
	continue;

      /* then it's a degenerate symbol.
       * According to alphabet, choose a single symbol to represent it.
       * watch out for too-clever scheme for random choice: "ABC"[random() % 3]
       */
      if (ALPHATYPE == kRNA)
	{
	  switch (*sym) {
	  case 'B': *sym = "CGU"[CHOOSE(3)]; break;
	  case 'D': *sym = "AGU"[CHOOSE(3)]; break;
	  case 'H': *sym = "ACU"[CHOOSE(3)]; break;
	  case 'K': *sym = "GU" [CHOOSE(2)]; break;
	  case 'M': *sym = "AC" [CHOOSE(2)]; break;
	  case 'R': *sym = "AG" [CHOOSE(2)]; break;
	  case 'S': *sym = "CG" [CHOOSE(2)]; break;
	  case 'T': *sym = 'U';                 break;
	  case 'V': *sym = "ACG"[CHOOSE(3)]; break;
	  case 'W': *sym = "AU" [CHOOSE(2)]; break;
	  case 'Y': *sym = "CU" [CHOOSE(2)]; break;
	  default:	Warn("Warning: unrecognized character %c in sequence\n", *sym);
	    /* break through to case 'N' */
	  case 'N': *sym = ALPHABET[CHOOSE(4)]; break;
	  }
	}
      else if (ALPHATYPE == kDNA)
	{
	  switch (*sym) {
	  case 'B': *sym = "CGT"[CHOOSE(3)]; break;
	  case 'D': *sym = "AGT"[CHOOSE(3)]; break;
	  case 'H': *sym = "ACT"[CHOOSE(3)]; break;
	  case 'K': *sym = "GT" [CHOOSE(2)]; break;
	  case 'M': *sym = "AC" [CHOOSE(2)]; break;
	  case 'R': *sym = "AG" [CHOOSE(2)]; break;
	  case 'S': *sym = "CG" [CHOOSE(2)]; break;
	  case 'U': *sym = 'T';                 break;
	  case 'V': *sym = "ACG"[CHOOSE(3)]; break;
	  case 'W': *sym = "AT" [CHOOSE(2)]; break;
	  case 'Y': *sym = "CT" [CHOOSE(3)]; break;
	  default:	Warn("Warning: unrecognized character %c in sequence\n", *sym);
	    /* break through to case 'N' */
	  case 'N': *sym = ALPHABET[CHOOSE(4)]; break;
	  }
	}
      else
	{
	  Warn("Warning: non-nucleic acid alphabet, unrecognized character %c in sequence\n", *sym);
	  *sym = ALPHABET[CHOOSE(ALPHASIZE)];
	}
    }
  return 1;
}



