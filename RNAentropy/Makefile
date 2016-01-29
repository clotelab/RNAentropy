LDFLAGS        = -lm
CCFLAGS        = -fPIC
CC             = gcc

all: RNAentropy 

RNAentropy: convert_Vienna.o loop_energies.o IL.o HP.o energy_par.o misc.o params.o utils.o McCaskillSimple.o fold_vars.o fold.o read_epars.o gquad.o aln_util.o structuralEntropy.o
	$(CC) convert_Vienna.o loop_energies.o IL.o HP.o energy_par.o misc.o params.o utils.o McCaskillSimple.o fold_vars.o fold.o read_epars.o gquad.o aln_util.o structuralEntropy.o $(LDFLAGS) -o RNAentropy 

McCaskillSimple.o: McCaskillSimple.c 
	$(CC) $(CCFLAGS) -c McCaskillSimple.c

misc.o: misc.c misc.h
	$(CC) $(CCFLAGS) -c misc.c

convert_Vienna.o: convert_Vienna.c convert_Vienna.h fold.h energy_const.h  \
	fold_vars.h pair_mat.h params.h 
	$(CC) $(CCFLAGS) -c convert_Vienna.c

loop_energies.o: loop_energies.c utils.h energy_par.h fold_vars.h pair_mat.h params.h
	$(CC) $(CCFLAGS) -c loop_energies.c

HP.o: HP.c HP.h convert_Vienna.h pair_mat.h
	$(CC) $(CCFLAGS) -c HP.c

IL.o: IL.c IL.h
	$(CC) $(CCFLAGS) -c IL.c

fold_vars.o: fold_vars.c fold_vars.h
	$(CC) $(CCFLAGS) -c fold_vars.c

fold.o: fold.c utils.h energy_par.h fold_vars.h pair_mat.h params.h
	$(CC) $(CCFLAGS) -c fold.c

utils.o: utils.c config.h
	$(CC) $(CCFLAGS) -c utils.c

energy_par.o: energy_par.c energy_const.h 
	$(CC) $(CCFLAGS) -c energy_par.c

params.o: params.c config.h energy_par.h fold_vars.h utils.h params.h
	$(CC) $(CCFLAGS) -c params.c

read_epars.o: read_epars.c
	$(CC) $(CCFLAGS) -c read_epars.c

gquad.o: gquad.c
	$(CC) $(CCFLAGS) -c gquad.c

aln_util.o: aln_util.c
	$(CC) $(CCFLAGS) -c aln_util.c

structuralEntropy.o: structuralEntropy.c 
	$(CC) $(CCFLAGS) -c structuralEntropy.c

clean:
	rm -fr *.o RNAentropy


