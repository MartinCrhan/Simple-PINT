all: PINT_simm

PINT_simm_omp: simmulate_PINT.f90
	@echo "Compiling PINT_simm"
	@gfortran -fopenmp -fbacktrace -fopt-info-all-omp -I/home/martin/bin/fftw-3.3.10/api/ -L/usr/local/lib/ -lfftw3 simmulate_PINT.f90 -lfftw3 -o PINT_simm_omp

PINT_simm_no_omp: simmulate_PINT.f90
	@echo "Compiling PINT_simm"
	@gfortran -I/home/martin/bin/fftw-3.3.10/api/ -L/usr/local/lib/ -lfftw3 simmulate_PINT.f90 -lfftw3 -o PINT_simm_no_omp

PINT_simm_acc: simmulate_PINT.f90
	@echo "Compiling PINT_simm"
	@nvfortran -acc -Minfo=acc -I/home/martin/bin/fftw-3.3.10/api/ -L/usr/local/lib/ -lfftw3 simmulate_PINT.f90 -lfftw3 -o PINT_simm_acc


clean:
	@echo "Cleaning up"
	@rm -f PINT_simm*
