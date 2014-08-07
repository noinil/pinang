default:

all: 	seq pdb_output pdb_cg_top matrix_transpose dcd_distance_0 ts_read pdb_dna_curve dcd_distance_1

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
dcd_distance_0:
	mkdir -p ../bin
	icpc -Wall -I../include/ -o p_dist_0 dcd_distances_0.cpp
	mv p_dist_0 ../bin/
dcd_distance_1:
	mkdir -p ../bin
	icpc -Wall -I../include/ -o p_dist_com dcd_distances_com.cpp
	mv p_dist_com ../bin/
ts_read:
	mkdir -p ../bin
	icpc -Wall -o p_ts_read cafemol_ts_read.cpp
	mv p_ts_read ../bin/
pdb_dna_curve:
	mkdir -p ../bin
	icpc -Wall -I../include/ -o p_dna_curve pdb_dna_curve.cpp
	mv p_dna_curve ../bin/
clean:
	rm -f ../bin/*
cct:
	rm -f ../bin/p_dist_com
