/* version.h
 * SRE, Mon Sep  6 09:32:36 1993
 * 
 * During a real build, these won't be used at all -- they'll
 * be overridden from the Makefile. This is only used during
 * manual compilations and to keep lint quiet. It also might
 * make a good record of the revisions the package goes through.
 * 
 */

/* #define DEBUG*/			/* turn all debugging output on */


#ifndef RELEASE
#define RELEASE      "2.4.4"
#define RELEASEDATE  "January 1996"
#endif
