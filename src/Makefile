# variables
INSTALL_DIR = ../bin/
CC = clang
CXX = clang++
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
.PHONY : all install clean

all : $(TARGETS)

install :
	@$(CHK_DIR_EXISTS) $(INSTALL_DIR) || $(MKDIR) $(INSTALL_DIR)
	mv p_* $(INSTALL_DIR)

$(TARGETS) : p_% : %.o
	$(LINK) $< -o $@
	@$(RM) $<

$(OBJECTS) : %.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $< -o $@

clean:
	@echo " Cleaning ..."
	@$(RM) ../bin/* $(OBJECTS) $(TARGETS)

