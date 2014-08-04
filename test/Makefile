default:

all: test_vec3d test_atom test_resid test_pdb_read test_dcd_read dna_curve_pdb

test_vec3d:
	icpc -Wall -I../include/ -o t1 vec3d_test.cpp
	./t1
test_atom:
	icpc -Wall -I../include/ -o t2 atom_test.cpp
	./t2
test_resid:
	icpc -Wall -I../include/ -o t3 resid_test.cpp
	./t3
test_pdb_read:
	icpc -Wall -I../include/ -o t5 pdb_read_test2.cpp
	./t5 ../data/1J7O.pdb tt.sop
test_dcd_read:
	icpc -Wall -I../include/ -o t6 read_dcd_test.cpp
	./t6 ../data/test.dcd
dna_curve_pdb:
	icpc -Wall -I../include/ -o t7 pdb_dna_curve.cpp
	./t7
clean:
	rm -f t1 t2 t3 t4 t5 tt.sop t6 test*.coor t7 _curve.pdb
