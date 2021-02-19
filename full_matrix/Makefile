CFLAGS = -g -Wall

eom: eom.o diag.o mmult.o init_4d_array.o block_matrix.o sort.o F.o W.o \
init_6d_array.o W3.o
	cc $^ -llapack -o $@

clean:
	rm -f *.o eom
