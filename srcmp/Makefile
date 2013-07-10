include ../Makefile.in

#VPATH = ../srccom

PROGS   = pmgsol mockmp

INSTALLCMD = install -s

OBJ     =  mproclib.o 

OBJ2    =  main_mp.o 

OBJ3    =  mock_mp.o

FLIB    = -lmge -lvisc -lroe -ltyp -lmpp $(MYLIB) $(MYMPP)

all:    $(PROGS)

pmgsol: $(OBJ) $(OBJ2)
	$(LINK) $(LFLAGS) -o $@ $^ $(FLIB) ;\
	mv $@ $(BINDIR)/$@

mockmp: $(OBJ) $(OBJ3)
	$(LINK) $(LFLAGS) -o $@ $^ $(FLIB) ;\
	mv $@ $(BINDIR)/$@

mgsol:  pmgsol

tst:
	echo TEST_MPI;\
	printenv PWD

.PHONY : clean
clean:
	rm -f *.o  *.mod

