LDFLAGS        = -lm
CCFLAGS        = -fPIC
CC             = gcc

all: RNAentropy

RNAentropy: convert_Vienna.o energy_par.o misc.o pair_mat.o params.o utils.o McCaskillSimple.o fold_vars.o fold.o structuralEntropy.o
	$(CC) convert_Vienna.o energy_par.o McCaskillSimple.o misc.o pair_mat.o params.o utils.o  fold_vars.o fold.o structuralEntropy.o $(LDFLAGS) -o RNAentropy 

convert_Vienna.o: convert_Vienna.c convert_Vienna.h
	$(CC) $(CCFLAGS) -c convert_Vienna.c

energy_par.o: energy_par.c energy_par.h
	$(CC) $(CCFLAGS) -c energy_par.c

McCaskillSimple.o: McCaskillSimple.c 
	$(CC) $(CCFLAGS) -c McCaskillSimple.c

misc.o: misc.c misc.h
	$(CC) $(CCFLAGS) -c misc.c

pair_mat.o: pair_mat.c pair_mat.h
	$(CC) $(CCFLAGS) -c pair_mat.c

params.o: params.c params.h
	$(CC) $(CCFLAGS) -c params.c

utils.o: utils.c utils.h
	$(CC) $(CCFLAGS) -c utils.c

fold_vars.o: fold_vars.c fold_vars.h
	$(CC) $(CCFLAGS) -c fold_vars.c

fold.o: fold.c fold.h
	$(CC) $(CCFLAGS) -c fold.c

structuralEntropy.o: structuralEntropy.c 
	$(CC) $(CCFLAGS) -c structuralEntropy.c

clean:
	rm -fr *.o RNAentropy


