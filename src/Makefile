default:

all: 		seq pdb_output pdb_cg_pos

seq:
	icpc -Wall -I../include/ -o p_seq get_sequence.cpp
	mv p_seq ../bin/
pdb_output:
	icpc -Wall -I../include/ -o p_pdb_extract pdb_output.cpp
	mv p_pdb_extract ../bin/
pdb_cg_pos:
	icpc -Wall -I../include/ -o p_pdb2pos pdb_cg_pos.cpp
	mv p_pdb2pos ../bin/
clean:
	rm -f ../bin/p_seq
