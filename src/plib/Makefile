# variables
CC = cc
CXX = c++
LINK = $(CXX)
HEADFILES = $(wildcard ../../include/*.hpp)
INCPATH = -I/usr/local/include/eigen3 -I../../include
CXXFLAGS = -std=c++11 -O3
OBJECTS = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
LIB_DIR = ../../lib
CHK_DIR_EXISTS= test -d
MKDIR = mkdir -p
RM = rm -fr

# rules
.PHONY : all clean

all : $(OBJECTS)
	ar r libpinang.a $?
	@$(CHK_DIR_EXISTS) $(LIB_DIR) || $(MKDIR) $(LIB_DIR)
	cp libpinang.a $(LIB_DIR)

$(OBJECTS) : %.o : %.cpp $(HEADFILES)
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $< -o $@

clean:
	@echo " Cleaning pinang lib ..."
	@$(RM) $(OBJECTS) libpinang.a
	@$(RM) $(LIB_DIR)/libpinang.a
