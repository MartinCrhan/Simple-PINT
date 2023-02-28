all: PINT_simm

PINT_simm: simmulate_PINT.f90
	@echo "Compiling PINT_simm"
	@gfortran -I/home/crhan/bin/fftw-3.3.10/api/ -L/usr/local/lib/ -lfftw3 simmulate_PINT.f90 -lfftw3 -o PINT_simm

clean:
	@echo "Cleaning up"
	@rm -f PINT_simm
