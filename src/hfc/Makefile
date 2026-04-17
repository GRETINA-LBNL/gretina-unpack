CC = g++

GRROOT_DIR := $(shell if [ -z $(GRROOT_DIR) ]; then echo "../.."; else echo $(GRROOT_DIR); fi)
LIB_DIR = $(shell if [ -z $(GRROOT_LIBDIR) ]; then echo "$(GRROOT_DIR)/lib/lib32"; else echo $(GRROOT_LIBDIR); fi)
BIN_DIR = $(shell if [ -z $(GRROOT_BINDIR) ]; then echo "$(GRROOT_DIR)/bin/bin32"; else echo $(GRROOT_BINDIR); fi)

FLAG =  -O3 -g
MFILE_LIB = 
RUNFILE = GEB_HFC

INCLUDE =
LIBS = 

OBJFILES = GEB_HFC.o HFC.o

../../$(RUNFILE): GEB_HFC.cpp $(OBJFILES) 
	$(CC) $(FLAG) $(OBJFILES) -o ../../$(RUNFILE) $(LIBS) 

GEB_HFC.o: GEB_HFC.cpp HFC.h
	$(CC) $(FLAG) -c $<

HFC.o: HFC.cpp HFC.h
	$(CC) $(FLAG) -c $<


clean:
	rm -f ../../$(RUNFILE) $(OBJFILES) 
