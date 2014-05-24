default:

all: test1 test2 test3

test1:
	icpc -Wall -I../include/ -o t1 vec3d_test.cpp
	./t1
test2:
	icpc -Wall -I../include/ -o t2 atom_test.cpp
	./t2
test3:
	icpc -Wall -I../include/ -o t3 resid_test.cpp
	./t3
clean:
	rm -f t1 t2 t3
