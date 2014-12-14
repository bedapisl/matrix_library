all:
	clr
	g++ main.cpp test.cpp -ggdb -std=c++11 -Wall
	if ./a.out; then echo "tests OK"; else echo "tests FAILED"; fi
