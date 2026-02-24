# Remote Access Memory Channels (RAMC) 0.5
# Copyright (2026) National Technology  Engineering Solutions of Sandia, LLC (NTESS). 
# Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains 
# certain rights in this software. */

OBJDIR=./obj/
LIBDIR=./lib/

OBJSUFFIX = .o

#
# Substitute correct compilers and paths
# -D__RAMC_DEBUG__ may report some additional info
#
CC=cc
CFLAGS=-I</path/to/libfabric/headers> -I</path/to/pmi/headers> -fPIC -D__CXI__
LDFLAGS=-L</path/to/libfabric/library> -L</path/to/pmi/library> -lfabric -lm -lpmi -shared

SRCS=$(wildcard *.c)
LIBTARGET=libramc.so
OBJTARGETS=$(addsuffix $(OBJSUFFIX), $(basename $(SRCS)))

all: $(LIBTARGET)

$(LIBTARGET): $(OBJTARGETS)
	mkdir -p $(LIBDIR)
	$(CC) -shared -o $(LIBDIR)$@ $(OBJDIR)*.o

%.o: %.c
	mkdir -p $(OBJDIR)
	$(CC) -c $< -o $(OBJDIR)$@ $(CFLAGS)
	#$(CC) -c $< -o $@ $(CFLAGS)

clean:
	rm -f *.o
	rm -f *.x
	rm -f $(OBJDIR)*.o
	rm -f $(LIBDIR)*.so

.PHONY: all clean

