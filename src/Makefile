default:

all: 	seq pdb_output pdb_cg_top matrix_transpose

seq:
	mkdir -p ../bin
	icpc -Wall -I../include/ -o p_seq get_sequence.cpp
	mv p_seq ../bin/
pdb_output:
	mkdir -p ../bin
	icpc -Wall -I../include/ -o p_pdb_extract pdb_output.cpp
	mv p_pdb_extract ../bin/
pdb_cg_top:
	mkdir -p ../bin
	icpc -Wall -I../include/ -o p_pdb2pin pdb_cg_top.cpp
	mv p_pdb2pin ../bin/
matrix_transpose:
	mkdir -p ../bin
	icpc -Wall -o p_transpose matrix_transpose.cpp
	mv p_transpose ../bin/
clean:
	rm -f ../bin/*
