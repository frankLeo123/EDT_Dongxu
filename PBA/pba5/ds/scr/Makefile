edtbin:main
	./main ${CASE} && gprof -b -p main
main:*.c
	gcc *.c  -lm -g -o main -std=c99 -pg
