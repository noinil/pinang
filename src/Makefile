# variables
INSTALL_DIR = ../bin/
LIB_DIR = ../lib
CC = cc
CXX = c++
LINK = $(CXX)
INCPATH = -I/usr/local/include -I../include
CXXFLAGS = -std=c++11 -O3
SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
TARGETS = $(patsubst %.cpp,p_%,$(SOURCES))
PINANGLIBS = $(patsubst %.cpp,%.o,$(wildcard plib/*.cpp))
CHK_DIR_EXISTS= test -d
MKDIR = mkdir -p
MOVE = mv -f
RM = rm -fr

# searching path
vpath %.h ../indluce
vpath %.a plib
vpath %.o plib

# rules
.PHONY : all install clean uninstall purge

all : $(PINANGLIBS) $(TARGETS)

install:
	@$(CHK_DIR_EXISTS) $(INSTALL_DIR) || $(MKDIR) $(INSTALL_DIR)
	cp -f p_* $(INSTALL_DIR)

$(PINANGLIBS) :
	@cd plib; make
$(TARGETS) : p_% : %.o
	$(LINK) $< $(LIB_DIR)/libpinang.a -o $@

$(OBJECTS) : %.o : %.cpp
	@echo " ------------------------------------------------------------ "
	@echo " Compiling $< ..."
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $< -o $@

clean:
	@echo " Cleaning ..."
	@cd plib; make clean
	@$(RM) $(OBJECTS) $(TARGETS)
uninstall:
	@$(RM) ../bin/*
purge: uninstall clean
	@echo " Cleaning everything ..."
