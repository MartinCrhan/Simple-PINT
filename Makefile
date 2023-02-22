all: PINT_simm

PINT_simm: simmulate_PINT.f90
	@echo "Compiling PINT_simm"
	@gfortran -I/usr/include -L/usr/lib/x86_64-linux-gnu -lfftw3 simmulate_PINT.f90 -lfftw3 -o PINT_simm

clean:
	@echo "Cleaning up"
	@rm -f PINT_simm