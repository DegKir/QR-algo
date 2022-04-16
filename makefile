all: tests
tests: tests.o matrix.o
	g++ -O3 tests.o matrix.o -o tests
tests.o: tests.cpp
	g++ -c tests.cpp -o tests.o
matrix.o: matrix.cpp
	g++ -c matrix.cpp -o matrix.o
clean:
	rm -rf tests. *.0
