default:

all: test1 test2

test1:
	g++ -Wall -I../include/ -o t1 vec3d_test.cpp
	./t1
test2:
	g++ -Wall -I../include/ -o t2 atom_test.cpp
	./t2
clean:
	rm -f t1 t2
