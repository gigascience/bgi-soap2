SHELL         = /bin/sh
PROG          = soap
DEBUG         = NO
PROFILE       = NO
PTHREADS      = YES
CC            = gcc
DEBUG_FLAGS   = -g3 -Wall -O2
PROFILE_FLAGS = -fprofile-arcs -ftest-coverage -pg 
RELEASE_FLAGS = -msse3 -O3 -static -funroll-loops -maccumulate-outgoing-args -fomit-frame-pointer 
STATIC_FLAGS  = -static
DFLAGS        = -DMAKE_TIME=\""`date`"\"
LIBS          = -lm
#TARBALL_EXCLUDE = "*.(o,gz,zip)"
#ZIP_EXCLUDE     = *.o *.gz *.zip

ifeq (YES, $(DEBUG))
        CFLAGS   = $(DEBUG_FLAGS) $(STATIC_FLAGS)
        DFLAGS  += -DDEBUG
#        PTHREADS = NO
else
        CFLAGS   = $(RELEASE_FLAGS) $(STATIC_FLAGS)
endif

ifeq (YES, $(PTHREADS))
        LIBS   +=  -lpthread
        DFLAGS +=  -DPTHREADS
endif

ifeq (YES, $(PROFILE))
       DFLAGS += $(PROFILE_FLAGS)
endif

OBJ = SeqIO.o MiscUtilities.o MemManager.o TextConverter.o r250.o DNACount.o HSP.o Timing.o BWT.o extratools.o soapio.o BWTAln.o Match.o PairMatch.o stdaln.o kstring.o
.SUFFIX:
.SUFFIX: .c .o

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@

all: $(PROG)

$(PROG): $(OBJ) soap.o
	$(CC) $(CFLAGS) $(DFLAGS) $(OBJ) soap.o -o $@ $(LIBS)

SeqIO.o:SeqIO.h
r250.o: r250.h
DNACount.o:DNACount.h
HSP.o:HSP.h
MiscUtilities.o:MiscUtilities.h
MemManager.o:MemManager.h
TextConverter.o:TextConverter.h
extratools.o:extratools.h BWT.h MiscUtilities.h MemManager.h TextConverter.h Timing.h HSP.h kstring.h
soapio.o:soapio.h SeqIO.h
BWT.o:BWT.h 
BWTAln.o:BWTAln.h BWT.h
Match.o:Match.h BWTAln.h soapio.h
PairMatch.o:Match.h BWTAln.h stdaln.h
MiscUtilities.o:MiscUtilities.h
MemManager.o:MemManager.h
TextConverter.o:TextConverter.h
stdaln.o:stdaln.h
kstring.o:kstring.h

clean:
	rm -f *.o $(PROG)
