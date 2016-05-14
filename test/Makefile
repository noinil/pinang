# variables
CC = cc
CXX = c++
LINK = $(CXX)
INCPATH = -I/usr/local/include/eigen3 -I../include
CXXFLAGS = -std=c++11 -O3
LIB_DIR = ../lib
OBJECTS = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
TARGETS = $(patsubst %.cpp,t_%,$(wildcard *.cpp))
MOVE = mv -f
RM = rm -fr

# rules
.PHONY : all install clean

all : $(TARGETS)

$(TARGETS) : t_% : %.o
	$(LINK) $< $(LIB_DIR)/libpinang.a -o $@

$(OBJECTS) : %.o : %.cpp
	@echo " ------------------------------------------------------------ "
	@echo " Compiling $< ..."
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $< -o $@

clean:
	@echo " Cleaning ..."
	@$(RM) ../bin/* $(OBJECTS) $(TARGETS)

