FC = mpif90
FFLAGS = -O3 -Wall -Wextra

OBJS = parameters_mod.o domain_mod.o solver_mod.o output_mod.o heat.o

all: heat results_dir viz_dir

heat: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

results_dir:
	mkdir -p results

viz_dir:
	mkdir -p viz

parameters_mod.o: parameters_mod.f90
	$(FC) $(FFLAGS) -c $<

domain_mod.o: domain_mod.f90 parameters_mod.o
	$(FC) $(FFLAGS) -c $<

solver_mod.o: solver_mod.f90 parameters_mod.o
	$(FC) $(FFLAGS) -c $<

output_mod.o: output_mod.f90 parameters_mod.o
	$(FC) $(FFLAGS) -c $<

heat.o: heat.f90 parameters_mod.o domain_mod.o solver_mod.o output_mod.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f heat *.o *.mod *.dat
	rm -rf results viz

run:
	mpirun -np 4 ./heat

.PHONY: all clean run results_dir viz_dir
