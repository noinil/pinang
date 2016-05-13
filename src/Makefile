# variables
INSTALL_DIR = ../bin/
CC = cc
CXX = c++
LINK = $(CXX)
INCPATH = -I/usr/local/include -I../include
CXXFLAGS = -std=c++11 -O3
OBJECTS = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
TARGETS = $(patsubst %.cpp,p_%,$(wildcard *.cpp))
CHK_DIR_EXISTS= test -d
MKDIR = mkdir -p
MOVE = mv -f
RM = rm -fr

# searching path
# vpath %.h ../indluce
# vpath %.h /usr/local/include

# rules
.PHONY : all install clean uninstall purge

all : $(TARGETS)

install:
	@$(CHK_DIR_EXISTS) $(INSTALL_DIR) || $(MKDIR) $(INSTALL_DIR)
	cp p_* $(INSTALL_DIR)

$(TARGETS) : p_% : %.o
	$(LINK) $< plib/libpinang.a -o $@

$(OBJECTS) : %.o : %.cpp
	@echo " ------------------------------------------------------------ "
	@echo " Compiling $< ..."
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $< -o $@

clean:
	@echo " Cleaning ..."
	@$(RM) $(OBJECTS) $(TARGETS)
uninstall:
	@$(RM) ../bin/*
purge: uninstall clean
	@echo " Cleaning everything ..."
