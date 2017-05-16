# Makefile for tRNAscan-SE 
# TML, Feb 11, 1996
#
 
COV_RELEASE    = "2.4.4"
EUFIND_RELEASE = "1.1"
SE_RELEASE = "1.3.1"
RELEASEDATE= "January 2012"
RFLAGS     = -DRELEASE=$(COV_RELEASE) -DRELEASEDATE=$(RELEASEDATE)

## Note: if you want to use the -i option, use "make no-ambig"
## This compiles versions of the pre-scanners that skip over N's,
## greatly speeding scans of unfinished sequence

## where Perl v.5 (or greater) is installed
PERLDIR = /usr/bin

## what you call the Perl v.5 executable on your system
## (this could be 'perl5' or otherwise on some systems)
PERLBIN = perl

## where you want things installed
BINDIR  = $(HOME)/bin
LIBDIR  = $(HOME)/lib/tRNAscan-SE
MANDIR  = $(HOME)/man

## NOTE !!  If you later manually move the location of 
## binaries or data files in the BINDIR or LIBDIR directories,
## you need to delete the tRNAscan-SE executable, update those
## locations in this Makefile, and re-make tRNAscan-SE 

## where you want temporary files stored
TEMPDIR = /tmp

## your compiler
CC = gcc	  # GNU cc (if available) otherwise use vendor's cc

## any special compiler flags you want
CFLAGS = -O    # ok for most machines (remove -O for DEC OSF/1 cc compiler)

## machine specific definitions
#  You shouldn't need any. The specific #define's in squid are historical.
#  However, if you do run into compilation trouble, you might need:
# -DNO_STRDUP    no strdup() function (lots of BSDish Unix)
# -DNOSTR        no strstr() function (rarer)
# 
MDEFS =

## how to install the man pages 
## either cp -- to copy unformatted man page
## or a script with identical syntax to cp, to format & install formatted page
INSTMAN   = cp
# INSTMAN   = instman
MANSUFFIX = 1

#######
## should not need to modify below this line
#######

SHELL  = /bin/sh
LIBS   = -lm
.SUFFIXES : .c .o

DOCS = README MANUAL INSTALL COPYING GNULICENSE FILES Release.history

MANSRC= tRNAscan-SE.man

PROGS = trnascan-1.4 covels-SE coves-SE eufindtRNA
UTILS = reformat revcomp seqstat shuffle
SCRIPTS = fasta2gsi.pl sstofa.pl checkversion.pl instman.pl

TSCAN13 = trnascan.c TPCsignal Dsignal

EUFIND = eufind_const.h eufind_main.c pavesi.c

SQUID-UTILS = seqstat_main.c reformat_main.c\
	      shuffle_main.c revcomp_main.c getopt.c

COVE_SUITE = covea coveb covee covet covels coves structcheck


MAIN =  align_main.c build_main.c emit_main.c learnprior_main.c train_main.c \
	scan_main.c score_main.c structcheck_main.c 

SQUIDSRC = alignio.c sqerror.c sqio.c iupac.c msf.c revcomp.c\
	selex.c sre_ctype.c sre_math.c sre_string.c stack.c types.c\
	gnuregex.c interleaved.c getopt.c

SQUIDOBJ = alignio.o sqerror.o sqio.o iupac.o msf.o revcomp.o\
	selex.o sre_ctype.o sre_math.o sre_string.o stack.o types.o\
	gnuregex.o interleaved.o getopt.o

SQUIDHDRS = squid.h sqfuncs.h gnuregex.h

SRC =   align.c dbviterbi.c debug.c emit.c fast-dbviterbi.c fastmodelmaker.c\
	konings.c lengthdist.c maxmodelmaker.c\
	misc.c model.c modelmaking.c prior.c probify.c save.c smallviterbi.c\
	scorestack.c structs.c trace.c viterbi.c $(SQUIDSRC)

MPSRC = mpviterbi.m mp-dbviterbi.m

OBJ =   align.o dbviterbi.o debug.o emit.o fast-dbviterbi.o fastmodelmaker.o\
	konings.o lengthdist.o maxmodelmaker.o\
	misc.o model.o modelmaking.o prior.o probify.o save.o smallviterbi.o\
	scorestack.o structs.o trace.o viterbi.o $(SQUIDOBJ)

MPOBJ = mpviterbi.o mp-dbviterbi.o

all: 	$(PROGS) tRNAscan-SE setpaths

covels-SE:  $(OBJ) scan_main.o
	$(CC) $(CFLAGS) $(RFLAGS) -o covels-SE scan_main.o $(OBJ) $(LIBS)

coves-SE: $(OBJ) score_main.o
	$(CC) $(CFLAGS) $(RFLAGS) -o coves-SE score_main.o $(OBJ) $(LIBS)

eufindtRNA: $(SQUIDOBJ) pavesi.o eufind_main.o
	$(CC) $(CFLAGS)  -o eufindtRNA eufind_main.o \
	pavesi.o $(SQUIDOBJ) $(LIBS)

trnascan-1.4: trnascan.o
	$(CC) $(CFLAGS) -DTSCANDIR=\"$(LIBDIR)\" -o trnascan-1.4 trnascan.c

tRNAscan-SE:
	$(PERLDIR)/$(PERLBIN) checkversion.pl
	sed 's#/tmp#$(TEMPDIR)#g' tRNAscan-SE.src | \
	sed 's#bindir = ""#bindir =\"$(BINDIR)/"#g' | \
	sed 's#/usr/local/lib/tRNAscanSE#$(LIBDIR)#g' | \
	sed 's#version = ""#version = $(SE_RELEASE)#' | \
	sed 's#release_date = ""#release_date = $(RELEASEDATE)#' | \
	sed 's#/usr/bin/perl#$(PERLDIR)/$(PERLBIN)#' > tRNAscan-SE
	sed 's#/usr/bin/perl#$(PERLDIR)/$(PERLBIN)#' instman.pl > instman 
	sed 's#/usr/bin/perl#$(PERLDIR)/$(PERLBIN)#' sstofa.pl > sstofa 
	@chmod +x tRNAscan-SE
	@chmod +x instman
	@echo ""
	@echo "tRNAscan-SE built"
	@echo ""
	@echo ""

setpaths:
	@echo 'setenv PATH "$$PATH"":""$(BINDIR)"' > setup.tRNAscan-SE
	@echo 'setenv PERL5LIB "$$PERL5LIB"":""$(BINDIR)"' >> setup.tRNAscan-SE
	@echo 'setenv MANPATH "$$MANPATH"":""$(MANDIR)"' >> setup.tRNAscan-SE
	@echo ""
	@echo "The file \"setup.tRNAscan-SE\" has been created."
	@echo "If you type \"source setup.tRNAscan-SE; rehash\", the paths"
	@echo "in your environment will be updated for this session."
	@echo ""
	@echo "To update your environment upon every login, you should"
	@echo "add the line:"
	@echo ""
	@echo "source $(PWD)/setup.tRNAscan-SE"
	@echo ""
	@echo "to your \".cshrc\" file."
	@echo ""
	@echo "If you prefer to manually update your environment variables,"
	@echo "be sure to make the following changes:"
	@echo ""
	@echo "1) Add $(BINDIR) to your PATH variable"
	@echo "2) Add $(BINDIR) to your PERl5LIB variable"
	@echo "3) Add $(MANDIR) to your MANPATH variable"
	@echo ""

install: $(PROGS) tRNAscanSE
	@if test -d $(BINDIR); then echo .; else mkdir -p $(BINDIR); fi
	@if test -d $(LIBDIR); then echo .; else mkdir -p $(LIBDIR); fi
	@if test -d $(TEMPDIR); then echo .; else mkdir -p $(TEMPDIR); fi
	@if test -d $(MANDIR)/man$(MANSUFFIX); then echo .; else mkdir -p $(MANDIR)/man$(MANSUFFIX); fi
	cp $(PROGS) tRNAscan-SE $(BINDIR)/.
	cp -R tRNAscanSE $(BINDIR)/
	cp TPCsignal Dsignal *.cm gcode.* $(LIBDIR)/.
	@if test -r trnascan-1.4-NA; then cp trnascan-1.4-NA $(BINDIR)/.; fi
	@if test -r eufindtRNA-NA; then cp eufindtRNA-NA $(BINDIR)/.; fi
	@if test -r mpcovels; then cp mpcovels $(BINDIR)/.; fi
	$(INSTMAN) $(MANSRC) $(MANDIR)/man$(MANSUFFIX)/tRNAscan-SE.$(MANSUFFIX) 

testrun:
	@echo ""
	@echo "Running tRNAscan-SE on sample sequence F22B2.fa"
	@echo ""
	@rm -f testrun.out
	@echo "tRNAscan-SE -d -y -o testrun.out Demo/F22B7.fa"
	@if tRNAscan-SE -d -y -o testrun.out Demo/F22B7.fa; \
	then echo ""; \
	else echo "tRNAscan-SE did NOT complete properly."; fi
	@if diff testrun.out testrun.ref; \
	then echo ""; cat testrun.out; \
	echo ""; echo ""; \
	echo "tRNAscan-SE up and running correctly!"; 	echo "";  \
	else echo "Test run result file differs from reference result file!";\
	echo "tRNAscan-SE may not have executed correctly."; echo ""; \
	echo "Please double check for correct compilation & installation"; \
	echo " and try 'make testrun' again."; 	echo ""; fi

# compiles versions of pre-scanners that do not conservatively
#  call tRNAs with ambiguous bases
#  Not really needed since Oct. 2000 modification to tRNAscan-SE (v 1.21)

noambig: trnascan-1.4-NA eufindtRNA-NA

trnascan-1.4-NA:  trnascan.o
	$(CC) $(CFLAGS) -DNO_AMBIG -DTSCANDIR=\"$(LIBDIR)\" -o trnascan-1.4-NA trnascan.c

eufindtRNA-NA: $(SQUIDOBJ) eufind_main.o
	$(CC) $(CFLAGS) -DNO_AMBIG -c -o pavesi-NA.o pavesi.c 
	$(CC) $(CFLAGS) -o eufindtRNA-NA eufind_main.o \
	pavesi-NA.o $(SQUIDOBJ) $(LIBS)	


clean:
	-rm -f *.o *~ Makefile.bak core testrun.out \
	$(PROGS) $(UTILS) $(COVE_SUITE) mpcovels tRNAscan-SE instman

uninstall:
	(cd $(BINDIR); rm -f $(PROGS) $(UTILS))
	(cd $(LIBDIR); rm -f TPCsignal Dsignal *.cm gcode.*)
	rmdir -ps $(LIBDIR)
	rm -f $(MANDIR)/man$(MANSUFFIX)/tRNAscan-SE.$(MANSUFFIX)
	rmdir -ps $(MANDIR)

.c.o:
	$(CC) $(CFLAGS) $(MDEFS) -c $<		

## programs from Sean Eddy's sequence i/o function library not
## needed for tRNAscan-SE but included for their utility

utils: $(UTILS)

install-utils:  $(UTILS)
	@if test -d $(BINDIR); then echo .; else mkdir -p $(BINDIR); fi
	sed 's#/usr/bin/perl#$(PERLDIR)/$(PERLBIN)#' fasta2gsi.pl > fasta2gsi
	cp $(UTILS) $(BINDIR)/.

reformat: $(SQUIDOBJ) reformat_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o reformat $(SQUIDOBJ) reformat_main.o $(LIBS)

revcomp: $(SQUIDOBJ) revcomp_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o revcomp $(SQUIDOBJ) revcomp_main.o $(LIBS)

seqstat: $(SQUIDOBJ) seqstat_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o seqstat $(SQUIDOBJ) seqstat_main.o $(LIBS)

shuffle: $(SQUIDOBJ) shuffle_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o shuffle $(SQUIDOBJ) shuffle_main.o $(LIBS)

## other programs in Cove package (below) not needed for
## tRNAscan-SE, but are included for users who wish to apply
## them for other RNA prediction work

all-cove: $(COVE_SUITE)

install-cove:  $(COVE_SUITE)
	@if test -d $(BINDIR); then echo .; else mkdir -p $(BINDIR); fi
	cp $(COVE_SUITE) $(BINDIR)/.

covea: $(OBJ) align_main.o
	$(CC) $(CFLAGS) $(RFLAGS) -o covea align_main.o $(OBJ) $(LIBS)

coveb: $(OBJ) build_main.o
	$(CC) $(CFLAGS) $(RFLAGS) -o coveb build_main.o $(OBJ) $(LIBS)

covee: $(OBJ) emit_main.o
	$(CC) $(CFLAGS) $(RFLAGS) -o covee emit_main.o $(OBJ) $(LIBS)

covet: $(OBJ) train_main.o
	$(CC) $(CFLAGS) $(RFLAGS) -o covet train_main.o $(OBJ) $(LIBS)

covels:  $(OBJ) scan_main.o
	$(CC) $(CFLAGS) $(RFLAGS) -o covels scan_main.o $(OBJ) $(LIBS)

coves: $(OBJ) score_main.o
	$(CC) $(CFLAGS) $(RFLAGS) -o coves score_main.o $(OBJ) $(LIBS)

structcheck: $(OBJ) structcheck_main.o
	$(CC) $(CFLAGS) $(RFLAGS) -o structcheck structcheck_main.o $(OBJ) $(LIBS)


## Maspar memory limits
#  do "mpctl" to see the system defaults for CMEM (ACU) and
#  PMEM (DPU) memory. The make will do an "mplimit" to reset
#  resources according to what you give for CMEM and PMEM needs.
#  
CMEM = 480K
PMEM = 63K

## Maspar version of covels
## need to set constants in "maspar.h" for your
##  maspar machine size (i.e. 128x128, 64x64, ...)

mpcovels: $(OBJ) $(MPOBJ) mpcovels_main.o
	mpl_cc $(RFLAGS) -o mpcovels mpcovels_main.o $(MPOBJ) $(OBJ) $(LIBS)
	mplimit mpcovels cmem $(CMEM) pmem $(PMEM)

