include ../Makefile.in

INSTALLCMD = install -s

OBJ   = bclib.o  codepak.o

LMGE  = libmge.a

LIBS  = -lmge -lvisc -lroe -ltyp $(MYLIB) 

$(LMGE): $(OBJ)
	$(AR) $@ $^ ;\
	$(RANLIB) $@

install: $(LMGE) 
	mv $^ $(LIBDIR) ;\
	$(RANLIB)  $(LIBDIR)/$^ ;\
	mv *.mod $(MODDIR)

mgsol:  main.o
	$(LINK) $(LFLAGS) -o $@ $^ $(LIBS) ;\
	mv $@ $(BINDIR)/$@

tst:    test_iforce.o
	$(LINK) $(LFLAGS) -o $@ $^ $(LIBS) ;\
	mv $@ $(BINDIR)/$@

tst_sing:test_rowechelon.o
	$(LINK) $(LFLAGS) -o $@ $^ $(LIBS) ;\
	mv $@ $(BINDIR)/$@


.PHONY : clean
clean:
	rm -f *.o *.mod

