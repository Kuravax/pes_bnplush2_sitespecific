FC = ifx
CFLAGS = -qmkl  -i8  -I"${MKLROOT}/include" -L${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lm -ldl  -qopenmp -liomp5 -lpthread -O3
sampling: noveltysampling.f90 pipnn.f90
	$(FC) -o pessample pipnn.f90 noveltysampling.f90 $(CFLAGS) 