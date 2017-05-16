/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* iupac.c
 * 
 * Globally defines the IUPAC symbols for nucleic acid sequence
 * Slowly evolving into a repository of globals. Tue Apr 20 1993
 */
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

char squid_version[] = "1.5g";
char squid_date[]    = "March 1996";

/* Dayhoff f(i) amino acid occurrence frequencies. From BLAST 1.0.5, NCBI.
 * In alphabetic order by single-letter code.
 */
float aafq[20] = {
  .08713, .03347, .04687, .04953, .03977, .08861, .03362, .03689, .08048, .08536, 
  .01475, .04043, .05068, .03826, .04090, .06958, .05854, .06472, .01049, .02992
};

char aa_alphabet[] = AMINO_ALPHABET;
				/* aa_index converts to pam's 27x27 scheme */
int  aa_index[20]  = { 0,  2,  3,  4,  5,  6,  7,  8, 10, 11, 
	              12, 13, 15, 16, 17, 18, 19, 21, 22, 24 };

				/* IUPAC code translations */
				/* note: sequence chars are UPPER CASE */
struct iupactype iupac[] = {
  { 'A', 'T', NTA, NTT, },
  { 'C', 'G', NTC, NTG, },
  { 'G', 'C', NTG, NTC, },
  { 'T', 'A', NTT, NTA, },
  { 'U', 'A', NTU, NTA, },
  { 'N', 'N', NTN, NTN, },
  { ' ', ' ', NTGAP, NTGAP, },
  { 'R', 'Y', NTR, NTY, },
  { 'Y', 'R', NTY, NTR, },
  { 'M', 'K', NTM, NTK, },
  { 'K', 'M', NTK, NTM, },
  { 'S', 'S', NTS, NTS, },
  { 'W', 'W', NTW, NTW, },
  { 'H', 'D', NTH, NTD, },
  { 'B', 'V', NTB, NTV, },
  { 'V', 'B', NTV, NTB, },
  { 'D', 'H', NTD, NTH, },
  };


char *stdcode1[65] = {
  "K",				/* AAA */
  "N",				/* AAC */
  "K",				/* AAG */
  "N",				/* AAU */
  "T",				/* ACA */
  "T",				/* ACC */
  "T",				/* ACG */
  "T",				/* ACU */
  "R",				/* AGA */
  "S",				/* AGC */
  "R",				/* AGG */
  "S",				/* AGU */
  "I",				/* AUA */
  "I",				/* AUC */
  "M",				/* AUG */
  "I",				/* AUU */
  "Q",				/* CAA */
  "H",				/* CAC */
  "Q",				/* CAG */
  "H",				/* CAU */
  "P",				/* CCA */
  "P",				/* CCC */
  "P",				/* CCG */
  "P",				/* CCU */
  "R",				/* CGA */
  "R",				/* CGC */
  "R",				/* CGG */
  "R",				/* CGU */
  "L",				/* CUA */
  "L",				/* CUC */
  "L",				/* CUG */
  "L",				/* CUU */
  "E",				/* GAA */
  "D",				/* GAC */
  "E",				/* GAG */
  "D",				/* GAU */
  "A",				/* GCA */
  "A",				/* GCC */
  "A",				/* GCG */
  "A",				/* GCU */
  "G",				/* GGA */
  "G",				/* GGC */
  "G",				/* GGG */
  "G",				/* GGU */
  "V",				/* GUA */
  "V",				/* GUC */
  "V",				/* GUG */
  "V",				/* GUU */
  "*",				/* UAA */
  "Y",				/* UAC */
  "*",				/* UAG */
  "Y",				/* UAU */
  "S",				/* UCA */
  "S",				/* UCC */
  "S",				/* UCG */
  "S",				/* UCU */
  "*",				/* UGA */
  "C",				/* UGC */
  "W",				/* UGG */
  "C",				/* UGU */
  "L",				/* UUA */
  "F",				/* UUC */
  "L",				/* UUG */
  "F",				/* UUU */
  "X",				/* unknown */
};




char *stdcode3[65] = {
  "Lys",			/* AAA */
  "Asn",			/* AAC */
  "Lys",			/* AAG */
  "Asn",			/* AAU */
  "Thr",			/* ACA */
  "Thr",			/* ACC */
  "Thr",			/* ACG */
  "Thr",			/* ACU */
  "Arg",			/* AGA */
  "Ser",			/* AGC */
  "Arg",			/* AGG */
  "Ser",			/* AGU */
  "Ile",			/* AUA */
  "Ile",			/* AUC */
  "Met",			/* AUG */
  "Ile",			/* AUU */
  "Gln",			/* CAA */
  "His",			/* CAC */
  "Gln",			/* CAG */
  "His",			/* CAU */
  "Pro",			/* CCA */
  "Pro",			/* CCC */
  "Pro",			/* CCG */
  "Pro",			/* CCU */
  "Arg",			/* CGA */
  "Arg",			/* CGC */
  "Arg",			/* CGG */
  "Arg",			/* CGU */
  "Leu",			/* CUA */
  "Leu",			/* CUC */
  "Leu",			/* CUG */
  "Leu",			/* CUU */
  "Glu",			/* GAA */
  "Asp",			/* GAC */
  "Glu",			/* GAG */
  "Asp",			/* GAU */
  "Ala",			/* GCA */
  "Ala",			/* GCC */
  "Ala",			/* GCG */
  "Ala",			/* GCU */
  "Gly",			/* GGA */
  "Gly",			/* GGC */
  "Gly",			/* GGG */
  "Gly",			/* GGU */
  "Val",			/* GUA */
  "Val",			/* GUC */
  "Val",			/* GUG */
  "Val",			/* GUU */
  "***",			/* UAA */
  "Tyr",			/* UAC */
  "***",			/* UAG */
  "Tyr",			/* UAU */
  "Ser",			/* UCA */
  "Ser",			/* UCC */
  "Ser",			/* UCG */
  "Ser",			/* UCU */
  "***",			/* UGA */
  "Cys",			/* UGC */
  "Trp",			/* UGG */
  "Cys",			/* UGU */
  "Leu",			/* UUA */
  "Phe",			/* UUC */
  "Leu",			/* UUG */
  "Trp",			/* UUU */
  "XXX",			/* unknown */
};
