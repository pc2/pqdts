#adapt the following to your system
CC = g++ # C++ compiler
FCC = gfortran #fortran compiler
PFCC = mpif90 #fortran compile with MPI
CFLAGS = -O3 -march=native -fopenmp -g -flto -fno-strict-aliasing
CFFLAGS = -ffree-line-length-none -Werror=aliasing -Werror=ampersand -Werror=c-binding-type -Werror=intrinsic-shadow -Werror=intrinsics-std -Werror=line-truncation -Werror=tabs -Werror=target-lifetime -Werror=underflow -Werror=unused-but-set-variable -Werror=unused-variable -Werror=unused-parameter -Werror=unused-label -Werror=conversion -Werror=zerotrip -Wno-maybe-uninitialized -Wuninitialized -Wuse-without-only -fno-strict-aliasing

#don't change the lines below
NAME = pqdts_omp.x
PNAME = pqdts_omp_mpi.x
OBJECTS = condat_simplexproj.o pqdts.o pqdts_mpi.o

pqdts: condat_simplexproj.o pqdts.o
	$(FCC) -o $(NAME) condat_simplexproj.o pqdts.o $(CFLAGS)

pqdts_mpi: condat_simplexproj.o pqdts_mpi.o
	$(PFCC) -o $(PNAME) -DMPI condat_simplexproj.o pqdts_mpi.o $(CFLAGS)

condat_simplexproj.o: condat_simplexproj.c
	$(CC) -c $(CFLAGS) $<

pqdts.o: pqdts.f90
	$(FCC) -c $(CFLAGS) $(CFFLAGS) -cpp pqdts.f90 -o pqdts.o

pqdts_mpi.o: pqdts.f90
	$(PFCC) -c $(CFLAGS) -DMPI $(CFFLAGS) -cpp pqdts.f90 -o pqdts_mpi.o

pretty:
	fprettify -i 2 pqdts.f90

clean:
	rm -f condat_simplexproj.o pqdts.o pqdts_mpi.o basics.mod $(NAME) $(PNAME)
