# Remote Access Memory Channels (RAMC) 0.5
#   Copyright (2026) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
#   Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
#   certain rights in this software. */

BINDIR=./bin/
OBJDIR=./obj/
RAMC=../src/lib/

OBJSUFFIX=.o
BINSUFFIX=.x

CC=cc
CPPC=CC+
CFLAGS=-I</path/to/libfabric/headers> -I</path/to/pmi/headers> -D__CXI__
LDFLAGS=-L</path/to/libfabric/library> -L</path/to/pmi/library> -L</path/to/ramc/library> -lfabric -lm -lpmi -lramc

SRCS=$(wildcard *.c)
OBJTARGETS=$(addsuffix $(OBJSUFFIX), $(basename $(SRCS)))
BINTARGETS=$(addsuffix $(BINSUFFIX), $(basename $(SRCS)))

all: $(BINTARGETS)

$(OBJTARGETS): %.o: %.c
	mkdir -p $(OBJDIR)
	$(CC) -c $< -o $(OBJDIR)$@ $(CFLAGS)

$(BINTARGETS): %.x: %.o
	mkdir -p $(BINDIR)
	$(CC) -o ${BINDIR}$@ -L${OBJDIR} -L$(RAMC) $(LDFLAGS) $(addprefix ${OBJDIR}, $^) ${LIBFLAGS}


clean:
	rm -f *.o
	rm -f *.x
	rm -f $(OBJDIR)*.o
	rm -f $(BINDIR)*.x

.PHONY: all clean

