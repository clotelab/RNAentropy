LDFLAGS        = -lm
CCFLAGS        = -fPIC
CC             = gcc

all: RNAentropy

RNAentropy: convert_Vienna.o IL.o HP.o energy_par.o misc.o params.o utils.o McCaskillSimple.o fold_vars.o fold.o read_epars.o gquad.o aln_util.o structuralEntropy.o
	$(CC) convert_Vienna.o IL.o HP.o energy_par.o misc.o params.o utils.o McCaskillSimple.o fold_vars.o fold.o read_epars.o gquad.o aln_util.o structuralEntropy.o $(LDFLAGS) -o RNAentropy 

McCaskillSimple.o: McCaskillSimple.c 
	$(CC) $(CCFLAGS) -c McCaskillSimple.c

misc.o: misc.c misc.h
	$(CC) $(CCFLAGS) -c misc.c

convert_Vienna.o: convert_Vienna.c convert_Vienna.h fold.h energy_const.h  \
	fold_vars.h pair_mat.h params.h 
	$(CC) $(CFLAGS) -c convert_Vienna.c

HP.o: HP.c HP.h convert_Vienna.h pair_mat.h
	$(CC) $(CFLAGS) -c HP.c

IL.o: IL.c IL.h
	$(CC) $(CFLAGS) -c IL.c

fold_vars.o: fold_vars.c fold_vars.h
	$(CC) $(CFLAGS) -c fold_vars.c

fold.o: fold.c utils.h energy_par.h fold_vars.h pair_mat.h params.h
	$(CC) $(CFLAGS) -c fold.c

utils.o: utils.c config.h
	$(CC) $(CFLAGS) -c utils.c

energy_par.o: energy_par.c energy_const.h 
	$(CC) $(CFLAGS) -c energy_par.c

params.o: params.c config.h energy_par.h fold_vars.h utils.h params.h
	$(CC) $(CFLAGS) -c params.c

read_epars.o: read_epars.c
	$(CC) $(CFLAGS) -c read_epars.c

gquad.o: gquad.c
	$(CC) $(CFLAGS) -c gquad.c

aln_util.o: aln_util.c
	$(CC) $(CFLAGS) -c aln_util.c

structuralEntropy.o: structuralEntropy.c 
	$(CC) $(CCFLAGS) -c structuralEntropy.c

clean:
	rm -fr *.o RNAentropy


