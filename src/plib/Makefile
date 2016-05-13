# variables
# INSTALL_DIR = ../bin/
CC = cc
CXX = c++
LINK = $(CXX)
INCPATH = -I/usr/local/include/eigen3 -I../../include
CXXFLAGS = -std=c++11 -O3 $(INCPATH)
OBJECTS = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
RM = rm -fr

# rules
.PHONY : all clean

all : $(OBJECTS)
	ar r libpinang.a $?

# $(OBJECTS) : %.o : %.cpp
# 	@echo " ------------------------------------------------------------ "
# 	@echo " Compiling $< ..."
# 	$(CXX) -c $(CXXFLAGS) $(INCPATH) $< -o $@

clean:
	@echo " Cleaning ..."
	@$(RM) $(OBJECTS) libpinang.a
