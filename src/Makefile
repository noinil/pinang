default:

all: 	install_dir pdb_related dcd_related

pdb_related: seq pdb_output pdb_cg_top matrix_transpose  ts_read pdb_dna_curve

dcd_related:  dcd_contact_num dcd_distance_0 dcd_distance_1 persistence_length dcd_dna_bending

install_dir:
	mkdir -p ../bin
seq:
	icpc -Wall -I../include/ -o p_seq get_sequence.cpp
	mv p_seq ../bin/
pdb_output:
	icpc -Wall -I../include/ -o p_pdb_extract pdb_output.cpp
	mv p_pdb_extract ../bin/
pdb_cg_top:
	icpc -Wall -I../include/ -o p_pdb2pin pdb_cg_top.cpp
	mv p_pdb2pin ../bin/
matrix_transpose:
	icpc -Wall -o p_transpose matrix_transpose.cpp
	mv p_transpose ../bin/
dcd_distance_0:
	icpc -Wall -I../include/ -o p_dist_0 dcd_distances_0.cpp
	mv p_dist_0 ../bin/
dcd_distance_1:
	icpc -Wall -I../include/ -o p_dist_com dcd_distances_com.cpp
	mv p_dist_com ../bin/
ts_read:
	icpc -Wall -o p_ts_read cafemol_ts_read.cpp
	mv p_ts_read ../bin/
pdb_dna_curve:
	icpc -Wall -I../include/ -o p_dna_curve pdb_dna_curve.cpp
	mv p_dna_curve ../bin/
dcd_contact_num:
	icpc -Wall -I../include/ -o p_contact_num dcd_contact_number.cpp
	mv p_contact_num ../bin/
persistence_length:
	icpc -Wall -I../include/ -o p_persistence_length dcd_persistence_length.cpp
	mv p_persistence_length ../bin/
dcd_dna_bending:
	icpc -Wall -I../include/ -o p_dcd_dna_bending dcd_dna_bending_curve.cpp
	mv p_dcd_dna_bending ../bin/
clean:
	rm -f ../bin/*
cct:
	rm -f ../bin/p_dna_curve
