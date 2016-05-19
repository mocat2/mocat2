all: normd normd_sw normd_range normd_aln normd_subaln normd_aln1

CC	= cc
CFLAGS  = -c -g 
LFLAGS	= -g -lm 

ROBJECTS = init.o readseq.o gcgcheck.o writeseq.o

normd : normd.c $(ROBJECTS) readmatrix.o showpair.o utils.o
	$(CC) normd.c $(ROBJECTS) readmatrix.o showpair.o utils.o -o normd $(LFLAGS)

normd_aln : normd_aln.c $(ROBJECTS) readmatrix.o showpair.o utils.o
	$(CC) normd_aln.c $(ROBJECTS) readmatrix.o showpair.o utils.o -o normd_aln $(LFLAGS)

normd_aln1 : normd_aln1.c $(ROBJECTS) readmatrix.o showpair.o utils.o
	$(CC) normd_aln1.c $(ROBJECTS) readmatrix.o showpair.o utils.o -o normd_aln1 $(LFLAGS)

normd_subaln : normd_subaln.c $(ROBJECTS) readmatrix.o utils.o
	$(CC) normd_subaln.c $(ROBJECTS) readmatrix.o utils.o -o normd_subaln $(LFLAGS)

normd_sw : normd_sw.c $(ROBJECTS) readmatrix.o showpair.o utils.o
	$(CC) normd_sw.c $(ROBJECTS) readmatrix.o showpair.o utils.o -o normd_sw $(LFLAGS)

normd_range : normd_range.c $(ROBJECTS) readmatrix.o showpair.o utils.o
	$(CC) normd_range.c $(ROBJECTS) readmatrix.o showpair.o utils.o -o normd_range $(LFLAGS)

.c.o :
	$(CC) $(CFLAGS) $?

