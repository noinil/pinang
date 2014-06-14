default:

all: 	seq pdb_output pdb_cg_top

seq:
	icpc -Wall -I../include/ -o p_seq get_sequence.cpp
	mv p_seq ../bin/
pdb_output:
	icpc -Wall -I../include/ -o p_pdb_extract pdb_output.cpp
	mv p_pdb_extract ../bin/
pdb_cg_top:
	icpc -Wall -I../include/ -o p_pdb2pin pdb_cg_top.cpp
	mv p_pdb2pin ../bin/
clean:
	rm -f ../bin/*
