all: cluspack

OPT = -O3

cluspack: principal.o lecfi.o filtering.o kmeans.o ecrfi.o dpc.o space.o divide.o \
	2means.o tools.o graphpc.o alignment.o treebionj.o bionj.o dissimilarities_clustering.o \
	secatorbionj.o weighting.o wardbionj.o treeward.o treecut.o treedpc.o normalized_2cut.o \
	symeig.o mixturemodels.o matrice.o
	cc -o cluspack principal.o lecfi.o filtering.o kmeans.o ecrfi.o dpc.o \
	space.o divide.o 2means.o tools.o graphpc.o alignment.o treebionj.o bionj.o \
	dissimilarities_clustering.o secatorbionj.o weighting.o wardbionj.o treeward.o \
	treecut.o treedpc.o normalized_2cut.o symeig.o mixturemodels.o matrice.o -lm $(OPT)

clean:
	rm -f *.o

principal.o: principal.c principal.h lecfi.h filtering.h ecrfi.h divide.h alignment.h treebionj.h \
	secatorbionj.h treeward.h treecut.h treedpc.h
	cc -o principal.o -c principal.c $(OPT)

lecfi.o: lecfi.c principal.h lecfi.h
	cc -o lecfi.o -c lecfi.c $(OPT)

filtering.o: filtering.c principal.h
	cc -o filtering.o -c filtering.c $(OPT)

kmeans.o: kmeans.c principal.h kmeans.h tools.h
	cc -o kmeans.o -c kmeans.c $(OPT)

ecrfi.o: ecrfi.c principal.h
	cc -o ecrfi.o -c ecrfi.c $(OPT)

divide.o: divide.c principal.h divide.h 2means.h dpc.h space.h tools.h graphpc.h \
	normalized_2cut.h
	cc -o divide.o -c divide.c $(OPT)

dpc.o: dpc.c principal.h divide.h space.h tools.h
	cc -o dpc.o -c dpc.c $(OPT)

space.o: space.c principal.h tools.h
	cc -o space.o -c space.c $(OPT)

tools.o: tools.c principal.h tools.h
	cc -o tools.o -c tools.c $(OPT)

2means.o: 2means.c principal.h divide.h 2means.h tools.h dpc.h
	cc -o 2means.o -c 2means.c $(OPT)

graphpc.o: graphpc.c principal.h divide.h tools.h graphpc.h
	cc -o graphpc.o -c graphpc.c $(OPT)

alignment.o: alignment.c principal.h
	cc -o alignment.o -c alignment.c $(OPT)

bionj.o: bionj.c bionj.h
	cc -o bionj.o -c bionj.c $(OPT)

treebionj.o: treebionj.c principal.h bionj.h treebionj.h
	cc -o treebionj.o -c treebionj.c $(OPT)

secatorbionj.o: secatorbionj.c principal.h secatorbionj.h dissimilarities_clustering.h \
	weighting.h wardbionj.h treecut.h
	cc -o secatorbionj.o -c secatorbionj.c $(OPT)

dissimilarities_clustering.o: dissimilarities_clustering.c principal.h tools.h \
	dissimilarities_clustering.h
	cc -o dissimilarities_clustering.o -c dissimilarities_clustering.c $(OPT)

weighting.o: weighting.c principal.h
	cc -o weighting.o -c weighting.c $(OPT)

wardbionj.o: wardbionj.c principal.h wardbionj.h
	cc -o wardbionj.o -c wardbionj.c $(OPT)

treeward.o: treeward.c principal.h treeward.h
	cc -o treeward.o -c treeward.c $(OPT)

treecut.o: treecut.c principal.h treecut.h dissimilarities_clustering.h
	cc -o treecut.o -c treecut.c $(OPT)

treedpc.o: treedpc.c principal.h treecut.h graphpc.h
	cc -o treedpc.o -c treedpc.c $(OPT)

normalized_2cut.o: normalized_2cut.c principal.h symeig.h tools.h normalized_2cut.h
	cc -o normalized_2cut.o -c normalized_2cut.c $(OPT)

symeig.o: symeig.c principal.h symeig.h
	cc -o symeig.o -c symeig.c $(OPT)

matrice.o: matrice.c principal.h matrice.h
	cc -o matrice.o -c matrice.c $(OPT)

mixturemodels.o: mixturemodels.c principal.h mixturemodels.h matrice.h kmeans.h
	cc -o mixturemodels.o -c mixturemodels.c $(OPT)









