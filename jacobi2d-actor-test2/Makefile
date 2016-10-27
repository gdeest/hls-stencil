

actor-test: actor-test.cpp jacobi2d.h common-util.cpp
	g++ -O3 -o $@ $< -g -fpermissive -std=c++11

ma-test: mem-alloc-test.cpp jacobi2d.h common-util.cpp
	g++ -O3 -o $@ $< -g -fpermissive


clean:
	rm ma-test actor-test
