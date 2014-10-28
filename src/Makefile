CPPINC=-I/usr/local/include

default:

all: 	install_dir pdb_related dcd_related

pdb_related: seq pdb_output pdb_cg_top matrix_transpose  ts_read pdb_dna_curve pdb_pdi_st pdb_wat_st

dcd_related:  dcd_contact_num dcd_distance_0 dcd_distance_1 persistence_length dcd_dna_bpn

install_dir:
	mkdir -p ../bin
seq:
	c++ -Wall -I../include/ -o p_seq get_sequence.cpp
	mv p_seq ../bin/
pdb_output:
	c++ -Wall -I../include/ -o p_pdb_extract pdb_output.cpp
	mv p_pdb_extract ../bin/
pdb_cg_top:
	c++ -Wall -I../include/ -o p_pdb2pin pdb_cg_top.cpp
	mv p_pdb2pin ../bin/
matrix_transpose:
	c++ -Wall -o p_transpose matrix_transpose.cpp
	mv p_transpose ../bin/
dcd_distance_0:
	c++ -Wall -I../include/ $(CPPINC) -o p_dist_0 dcd_distances_0.cpp
	mv p_dist_0 ../bin/
dcd_distance_1:
	c++ -Wall -I../include/ $(CPPINC) -o p_dist_com dcd_distances_com.cpp
	mv p_dist_com ../bin/
ts_read:
	c++ -Wall -o p_ts_read cafemol_ts_read.cpp
	mv p_ts_read ../bin/
pdb_dna_curve:
	c++ -Wall -I../include/ -o p_dna_curve pdb_dna_curve.cpp
	mv p_dna_curve ../bin/
dcd_contact_num:
	c++ -Wall -I../include/ $(CPPINC) -o p_contact_num dcd_contact_number.cpp
	mv p_contact_num ../bin/
persistence_length:
	c++ -Wall -I../include/ $(CPPINC) -o p_persistence_length dcd_persistence_length.cpp
	mv p_persistence_length ../bin/
dcd_dna_bpn:
	c++ -Wall -I../include/ $(CPPINC) -o p_dcd_dna_bp_number dcd_base_pairing_percentage.cpp
	mv p_dcd_dna_bp_number ../bin/
pdb_pdi_st:
	c++ -Wall -I../include/ -o p_pdi_st pdb_PDI_st.cpp
	mv p_pdi_st ../bin/
pdb_wat_st:
	c++ -Wall -I../include/ -o p_pwd_st pdb_PWD_st.cpp
	c++ -Wall -I../include/ -o p_pwp_st pdb_PWP_st.cpp
	mv p_pwp_st ../bin/
	mv p_pwd_st ../bin/

clean:
	rm -f ../bin/*
cct:
	rm -f ../bin/p_wat_st
