##
## NOTE: EDIT STR TO THE PATH TO SOURCE CODES 
## 
STR    := $(HOME)/code_mge
MY_DIR := $(shell printenv MY_DIR)
ifeq ($(MY_DIR), $(NULL))
export MY_DIR:=$(STR)
$(info $(MY_DIR) ) 
endif


DIR1   = $(MY_DIR)/srcf90
DIRROE = $(MY_DIR)/plugin_roe
DIRVIS = $(MY_DIR)/plugin_visc
DIRTYP = $(MY_DIR)/srctype

ifdef MPI
  DIR2   = $(MY_DIR)/srcmp
else
  DIR2   = $(MY_DIR)/srcf90
endif

$(info MY_DIR=$(MY_DIR) ) 

INC    = $(MY_DIR)/inc
MOD    = $(MY_DIR)/mod

ltyp:
	cd $(DIRTYP); pwd;\
	make -k -f Makefile clean;\
	make -k -f Makefile install

lroe:
	cd $(DIRROE); pwd;\
	make -k -f Makefile clean;\
	make -k -f Makefile install

lvisc:
	cd $(DIRVIS); pwd;\
	make -k -f Makefile clean;\
	make -k -f Makefile install

lmge:
	cd $(DIR1); pwd;\
	make -k -f Makefile clean;\
	make -k -f Makefile install

mgsol:
	cd $(DIR2); pwd;\
	make -k -f Makefile clean;\
	make -k -f Makefile $@;\

mockmp:
	cd $(DIR2); pwd;\
	make -k -f Makefile clean;\
	make -k -f Makefile $@;\

tst:
	cd $(DIR2); pwd;\
	make -k -f Makefile clean;\
	make -k -f Makefile $@;\

tst_sing:
	cd $(DIR2); pwd;\
	make -k -f Makefile clean;\
	make -k -f Makefile $@;\


#all: cleanall ltyp lroe lvisc lmge mgsol
fresh: clean lmge mgsol

.PHONY : clean

clean:
	rm -f $(INC)/*.mod ;\
	echo "rm *.mod in $(INC)"

cleanall:
	rm -f $(MOD)/*.mod ;\
        echo "rm *.mod in $(INC)"

