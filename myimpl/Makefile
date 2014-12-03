outseq: general_io.c seq_impl.c
	gcc -lm general_io.c seq_impl.c -o outseq
outmpi: mpi_io.c general_io.c mpi_impl.c
	mpicc -lm mpi_io.c general_io.c mpi_impl.c -o outmpi