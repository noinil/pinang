default:

all: test1 test2 test3 test4

test1:
	icpc -Wall -I../include/ -o t1 vec3d_test.cpp
	./t1
test2:
	icpc -Wall -I../include/ -o t2 atom_test.cpp
	./t2
test3:
	icpc -Wall -I../include/ -o t3 resid_test.cpp
	./t3
test4:
	icpc -Wall -I../include/ -o t4 pdb_read_test.cpp
	./t4 ../data/aa.pdb
test5:
	icpc -Wall -I../include/ -o t5 pdb_read_test2.cpp
	./t5 ../data/aa.pdb
test6:
	icpc -Wall -I../include/ -o t6 constants_test.cpp
	./t6
clean:
	rm -f t1 t2 t3 t4 t5
