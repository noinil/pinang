# variables
CC = cc
CXX = c++
LINK = $(CXX)
INCPATH = -I/usr/local/include -I../include
CXXFLAGS = -std=c++11 -O3
OBJECTS = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
TARGETS = $(patsubst %.cpp,t_%,$(wildcard *.cpp))
MOVE = mv -f
RM = rm -fr

# rules
.PHONY : all install clean

all : $(TARGETS)

$(TARGETS) : t_% : %.o
	$(LINK) $< -o $@
	@$(RM) $<

$(OBJECTS) : %.o : %.cpp
	@echo " ------------------------------------------------------------ "
	@echo " Compiling $< ..."
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $< -o $@

clean:
	@echo " Cleaning ..."
	@$(RM) ../bin/* $(OBJECTS) $(TARGETS)

