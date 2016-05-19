install: scan_coreblocks scan_orphanerrs scan_seqerrs scan_blockerrs realign_seqerrs \
realign_orphanerrs realign_nonblock realign_blockerrs tomsf totfa

clean:
	rm *.o

ROBJECTS = readseq.o gcgcheck.o 

OBJECTS = readseq.o writeseq.o init.o util.o gcgcheck.o readtree.o malign.o \
        seqweight.o readmat.o quicknj.o prfalign.o calcprf1.o calcprf2.o calcgapcoeff.o \
	rascal_util.o

HEADERS = rascal.h

CC	= cc
CFLAGS  = -c -O 
LFLAGS	= -O -lm 


tomsf : $(OBJECTS) tomsf.o
	$(CC) -o $@ $(OBJECTS) tomsf.o $(LFLAGS)

totfa : $(OBJECTS) totfa.o
	$(CC) -o $@ $(OBJECTS) totfa.o $(LFLAGS)

scan_coreblocks : $(OBJECTS) scan_coreblocks.o
	$(CC) -o $@ $(OBJECTS) scan_coreblocks.o $(LFLAGS)

scan_orphanerrs : $(OBJECTS) scan_orphanerrs.o
	$(CC) -o $@ $(OBJECTS) scan_orphanerrs.o $(LFLAGS)

scan_seqerrs : $(OBJECTS) scan_seqerrs.o
	$(CC) -o $@ $(OBJECTS) scan_seqerrs.o $(LFLAGS)

scan_blockerrs : $(OBJECTS) scan_blockerrs.o
	$(CC) -o $@ $(OBJECTS) scan_blockerrs.o $(LFLAGS)

realign_seqerrs : $(OBJECTS) realign_seqerrs.o
	$(CC) -o $@ $(OBJECTS) realign_seqerrs.o $(LFLAGS)

realign_orphanerrs : $(OBJECTS) realign_orphanerrs.o
	$(CC) -o $@ $(OBJECTS) realign_orphanerrs.o $(LFLAGS)

realign_nonblock : $(OBJECTS) realign_nonblock.o
	$(CC) -o $@ $(OBJECTS) realign_nonblock.o $(LFLAGS)

realign_blockerrs : $(OBJECTS) realign_blockerrs.o
	$(CC) -o $@ $(OBJECTS) realign_blockerrs.o $(LFLAGS)

pairalign.o : pairalign.c pairalign.h $(HEADERS)
	$(CC) $(COPT) $(CFLAGS) $*.c

.c.o :	$(HEADERS)
	$(CC) $(CFLAGS) $?

