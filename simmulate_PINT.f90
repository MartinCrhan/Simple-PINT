program simmulate_PINT
implicit none
! A simple currently unspecified dimensional program for running PINT simmulations in simple potentials
integer :: N 
integer :: t
integer :: stride
real :: tau
real :: m
real :: temperature
real :: alpha
real :: gamma
real :: box
real :: target_freq
integer :: nbead
integer :: n_dir
integer :: parameter_number
character(len = 8) :: int_type
character(len = 8) :: inp_type
character(len = 9) :: buffer
character(len = 4) :: freq_type
character(len = 3) :: thermostating
character(len = 3) :: centroid_constraint
character(len = 3) :: out_bead
character(len = 3) :: out_pos
character(len = 3) :: out_mom
character(len = 3) :: out_force
! probly does not have to be declared ?
open(1, file = 'input.inp', status = 'old')
read(1,*) buffer, n_dir
read(1,*) buffer, nbead
read(1,*) buffer, temperature
read(1,*) buffer, gamma ! centroid friction coefficient
read(1,*) buffer, alpha ! friction coefficient scaling factor for non-centroid DOFs
read(1,*) buffer, int_type
read(1,*) buffer, inp_type
read(1,*) buffer, freq_type
read(1,*) buffer, thermostating
read(1,*) buffer, centroid_constraint
read(1,*) buffer, t
read(1,*) buffer, tau
read(1,*) buffer, m
read(1,*) buffer, box
read(1,*) buffer, target_freq
read(1,*) buffer, stride
read(1,*) buffer, out_bead
read(1,*) buffer, out_pos
read(1,*) buffer, out_mom
read(1,*) buffer, out_force


close(1)

if (inp_type == 'fromgrid') then

    N = n_dir ** 3

elseif (inp_type == 'frominpt') then

    N = n_dir ! there is no reason for it to be a number cubed, or to really create a new variable just for this ...

end if

if (int_type == 'harmonic') then

    parameter_number = 1

elseif (int_type == '1Ddouble') then

    parameter_number = 2

elseif (int_type == 'McKenzie') then

    parameter_number = 6

elseif (int_type == '2D_Morse') then

    parameter_number = 3

end if

call run_simmulation(n_dir, N, nbead, thermostating, temperature, gamma, alpha, int_type, inp_type, freq_type, target_freq, &
 parameter_number, centroid_constraint, t, tau, m, box, stride, out_bead, out_pos, out_mom, out_force)

end program

subroutine run_simmulation(n_dir, N, nbead, thermostating, temperature, gamma, alpha, interaction, inp_type, &
 freq_type, target_freq, parameter_number, centroid_constraint, t, tau, m, box, stride, out_bead, out_pos, out_mom, &
  out_force)
implicit none
integer :: N
integer :: t
integer n_dir
integer :: stride 
real :: tau
real :: box
real :: temperature
real :: m
real :: gamma
real :: alpha
real :: target_freq
double precision, dimension(nbead) :: nm_masses
real, dimension(parameter_number) :: parameters
double precision, dimension(nbead, 3, N) :: p, p_nm, q_nm, F
double precision, dimension(nbead, nbead) :: nm_matrix
double precision, dimension(nbead) :: frequencies
integer :: i
integer :: k
integer :: l
integer :: Nsteps
integer :: nbead
integer :: parameter_number
character(len = 8) :: fmt
character(len = 8) :: x
character(len = 8) :: interaction
character(len = 8) :: inp_type
character(len = 4) :: freq_type
character(len = 3) :: thermostating
character(len = 3) :: centroid_constraint
character(len = 3) :: out_bead
character(len = 3) :: out_pos
character(len = 3) :: out_mom
character(len = 3) :: out_force
real(8), parameter :: PI = 4 * atan (1.0_8)
fmt = '(I2.2)'

if (out_pos == 'yay') then
    open(2*nbead + 3, file = 'positions_centroid.xyz', status = 'new')
    if (out_bead == 'yay') then
        do i = 1,nbead
            write (x,fmt) (i-1) 
            open(2*i + 1, file = 'positions_'//trim(x)//'.xyz', status = 'new')
        end do
    end if
end if
if (out_mom == 'yay') then
    open(2*nbead + 4, file = 'momenta_centroid.xyz', status = 'new')
    if (out_bead == 'yay') then
        do i = 1,nbead
            write (x,fmt) (i-1)
            open(2*i + 2, file = 'momenta_'//trim(x)//'.xyz', status = 'new')
        end do
    end if
end if
if (out_force == 'yay') then
    open(4*nbead + 5, file = 'forces_centroid.xyz', status = 'new')
    if (out_bead == 'yay') then
        do i = 1,nbead
            write (x,fmt) (i-1) 
            open(i + 3*nbead + 4, file = 'forces_'//trim(x)//'.xyz', status = 'new')
        end do
    end if
end if

open(2, file = 'energies.dat', status = 'new')

open(4*nbead + 6, file = 'energies_classical.dat', status = 'new')

! Initialize nm_matrix and frequencies

call init_nm(nm_matrix, frequencies, nbead, temperature, freq_type, target_freq, m , nm_masses)

! Initialize positions

if (inp_type == 'fromgrid') then
    call grid(q_nm, N, nbead, n_dir, box, nm_matrix)
elseif ( inp_type == 'frominpt' ) then
    call init_coordinates(q_nm, nbead, N, nm_matrix)
end if

! Initialize forces

call init_interactions(interaction, parameters, parameter_number)

call calc_forces(q_nm, F, N, nbead, m, interaction, parameters, parameter_number, nm_matrix)

! Initialize momenta

call init_momenta(p_nm, nm_masses, nm_matrix, temperature, nbead, N)

Nsteps = t/tau

! Obabo numerical evolution of the system

do i=1,Nsteps
    write(1,*) interaction
    if (modulo((i - 1), stride) == 0) then
        call do_output(q_nm,p_nm,F,nbead,N,nm_matrix,frequencies,nm_masses,temperature,interaction, &
         parameters, parameter_number, out_bead, out_pos, out_mom, out_force)
    end if
    if (thermostating == 'yay') then
        call thermostat(p_nm, temperature, gamma, alpha, N, nbead, nm_masses, frequencies, tau/2)
    end if
    call momentum_step(p_nm,F,N,nbead,k,l,nm_matrix,tau/2)
    call replica_step(q_nm,p_nm,N,nbead,frequencies,nm_masses,tau,centroid_constraint)
    call calc_forces(q_nm,F,N,nbead,m,interaction,parameters,parameter_number, nm_matrix)
    call momentum_step(p_nm,F,N,nbead,k,l,nm_matrix,tau/2)
    if (thermostating == 'yay') then
        call thermostat(p_nm, temperature, gamma, alpha, N, nbead, nm_masses, frequencies, tau/2)
    end if
end do

end subroutine

!
!  Introduced for a perhaps better modularity
!

subroutine to_nm(x, x_nm , nm_matrix, nbead, N)
double precision, dimension(nbead, 3, N) :: x, x_nm
double precision, dimension(nbead, nbead) :: nm_matrix
integer :: nbead
integer :: N
do k = 1,3
    do j = 1,N
        do i = 1,nbead
            x_nm(i,k,j) = 0
            do l = 1, nbead
                x_nm(i,k,j) = x_nm(i,k,j) + nm_matrix(i, l) * x(l,k,j)
            end do
        end do 
    end do
end do

end subroutine

subroutine from_nm(x, x_nm , nm_matrix, nbead, N)
double precision, dimension(nbead, 3, N) :: x, x_nm
double precision, dimension(nbead, nbead) :: nm_matrix
integer :: nbead
integer :: N

do k = 1,3
    do j = 1,N
        do i = 1,nbead
            x(i,k,j) = 0
            do l = 1, nbead
                x(i,k,j) = x(i,k,j) + nm_matrix(l, i) * x_nm(l,k,j)
            end do
        end do 
    end do
end do

end subroutine

subroutine to_nm_fftw(x, x_nm, nbead, N)
! Kinda awkward at the moment
use, intrinsic :: iso_c_binding 
implicit none
include 'fftw3.f03'
double precision, dimension(nbead, 3, N) :: x, x_nm
complex (kind = 8), dimension(nbead, 3, N) :: fft_array
type(c_ptr) :: plan
integer :: nbead
integer :: N
integer :: i, j, k, l, halfpoint_1, halfpoint_2

halfpoint_1 = nbead/2
halfpoint_2 = (nbead/2 + 1)

do k = 1,3
    do j = 1,N
        do i = 1,nbead
            fft_array(i,k,j) = complex(x(i,k,j), 0)
        end do 
    end do
end do

plan = fftw_plan_dft_1d(nbead, fft_array(1:nbead, 1, 1), fft_array(1:nbead, 1, 1) , FFTW_FORWARD,FFTW_ESTIMATE)

do k = 1,3
    do j = 1,N
        call fftw_execute_dft(plan, fft_array(1:nbead, k, j), fft_array(1:nbead, k, j))
    end do
end do

call fftw_destroy_plan(plan)

do k = 1,3
    do j = 1,N
        if (modulo(nbead, 2) == 0) then
            x_nm(1, k, j) = realpart(fft_array(1, k, j)) / (nbead ** 0.5)
            do l = 2,halfpoint_1
                x_nm(l, k, j) = realpart(fft_array(l, k, j)) * (( 2 ** 0.5) / (nbead ** 0.5))
                x_nm(nbead - l + 2, k, j) = imagpart(fft_array(l, k, j)) * (( 2 ** 0.5) / (nbead ** 0.5))
            end do
            x_nm(halfpoint_2, k , j) = realpart(fft_array(halfpoint_2, k, j)) / (nbead ** 0.5)
        else
            x_nm(1, k, j) = realpart(fft_array(1, k, j)) / (nbead ** 0.5)
            do l = 2,halfpoint_2
                x_nm(l, k, j) = realpart(fft_array(l, k, j)) * (( 2 ** 0.5) / (nbead ** 0.5))
                x_nm(nbead - l + 2, k, j) = imagpart(fft_array(l, k, j)) * (( 2 ** 0.5) / (nbead ** 0.5))
            end do
        end if
    end do
end do

end subroutine

subroutine from_nm_fftw(x, x_nm, nbead, N)
use, intrinsic :: iso_c_binding 
implicit none
include 'fftw3.f03'
double precision, dimension(nbead, 3, N) :: x, x_nm
complex (kind = 8), dimension(nbead, 3, N) :: fft_array
type(c_ptr) :: plan
integer :: nbead
integer :: N
integer :: i, j, k, l, halfpoint_1, halfpoint_2

halfpoint_1 = nbead/2
halfpoint_2 = (nbead/2 + 1)

do k = 1,3
    do j = 1,N
        if (modulo(nbead, 2) == 0) then
            fft_array(1, k, j) = complex(x_nm(1, k, j) * (nbead ** 0.5), 0)
            do l = 2,halfpoint_1
                fft_array(l, k, j) = complex(x_nm(l, k, j) * ((nbead ** 0.5) / (2 ** 0.5)), &
                 x_nm(nbead - l + 2, k, j) * ((nbead ** 0.5) / (2 ** 0.5)))
                fft_array(nbead - l + 2, k, j) = complex(x_nm(l, k, j) * ((nbead ** 0.5) / (2 ** 0.5)), &
                 - x_nm(nbead - l + 2, k, j) * ((nbead ** 0.5) / (2 ** 0.5)))
            end do
            fft_array(halfpoint_2, k, j) = complex(x_nm(halfpoint_2, k, j) * (nbead ** 0.5), 0)
        else
            fft_array(1, k, j) = complex(x_nm(1, k, j) * (nbead ** 0.5), 0)
            do l = 2,halfpoint_2
                fft_array(l, k, j) = complex(x_nm(l, k, j) * ((nbead ** 0.5) / (2 ** 0.5)), &
                 x_nm(nbead - l + 2, k, j) * ((nbead ** 0.5) / (2 ** 0.5)))
                fft_array(nbead - l + 2, k, j) = complex(x_nm(l, k, j) * ((nbead ** 0.5) / (2 ** 0.5)), &
                 - x_nm(nbead - l + 2, k, j) * ((nbead ** 0.5) / (2 ** 0.5)))
            end do
        end if
    end do
end do

plan = fftw_plan_dft_1d(nbead, fft_array(1:nbead, 1, 1), fft_array(1:nbead, 1, 1) , FFTW_BACKWARD, FFTW_ESTIMATE)

do k = 1,3
    do j = 1,N
        call fftw_execute_dft(plan, fft_array(1:nbead, k, j), fft_array(1:nbead, k, j))
    end do
end do

call fftw_destroy_plan(plan)

do k = 1,3
    do j = 1,N
        do i = 1,nbead
            ! This differs from th pyhton case, as numpy actually
            ! includes the division by nbead in it's definition of the fft
            x(i,k,j) = realpart(fft_array(i, k, j)) / nbead
        end do 
    end do
end do

end subroutine

subroutine init_momenta(p_nm, nm_masses, nm_matrix, temperature, nbead, N)
double precision, dimension(nbead, 3, N) :: p_nm
double precision, dimension(nbead, 3, 2 * N) :: init
double precision, dimension(nbead, nbead) :: nm_matrix
double precision, dimension(nbead) :: nm_masses
real :: temperature
integer :: nbead
integer :: N
integer :: i,l,k,j

call random_number(init)
do l = 1,nbead
    do k = 1,N
        do j = 1,3
            p_nm(l,j,k) = ((- 2 * log(init(l,j,k))) ** 0.5) * cos(2 * PI * init(l,j, k + N)) &
             * (nbead * temperature * nm_masses(l))**0.5
        end do
    end do
end do

end subroutine

subroutine grid(q_nm, N, nbead, n_dir, box, nm_matrix)
double precision, dimension(nbead, 3, N) :: q, q_nm
double precision, dimension(nbead, nbead) :: nm_matrix
integer :: N
integer :: nbead
integer n_dir
real :: box
double precision :: increment
integer :: i
integer :: j
integer :: k
integer :: l
integer :: loc
increment = box / (n_dir + 1)
do l = 1,nbead
    do i = 0,(n_dir - 1)
        do j = 0,(n_dir - 1)
            do k = 0,(n_dir - 1)
            loc = i + n_dir * j + (n_dir ** 2) * k + 1
            q(l,1,loc) = increment * (i + 1)
            q(l,2,loc) = increment * (j + 1)
            q(l,3,loc) = increment * (k + 1)
            end do
        end do
    end do
end do

call to_nm_fftw(q, q_nm, nbead, N)

end subroutine

subroutine init_coordinates(q_nm, nbead, N, nm_matrix)
double precision, dimension(nbead, 3, N) :: q, q_nm
double precision, dimension(nbead, nbead) :: nm_matrix
integer :: N
integer :: nbead
character(len = 8) :: fmt
character(len = 8) :: x
character(len = 2) :: buffer
fmt = '(I2.2)'

do i = 1,nbead
    write (x,fmt) (i-1) 
    open(i + 2*nbead + 4, file = 'input_'//trim(x)//'.xyz', status = 'old')
    read(i + 2*nbead + 4,*)
    read(i + 2*nbead + 4,*)
    do j = 1,N
        read(i + 2*nbead + 4,*) buffer, q(i,1,j), q(i,2,j), q(i,3,j)
    end do
    close(i + 2*nbead + 4)
end do

call to_nm_fftw(q, q_nm, nbead, N)

end subroutine


subroutine init_interactions(interaction, parameters, parameter_number)
integer :: parameter_number
real, dimension(parameter_number) :: parameters
character(len = 8) :: interaction
open(1, file = 'interactions.inp', status = 'old')
if (interaction == 'harmonic') then
    ! omega
    read(1,*)
    read(1,*) parameters(1)
elseif (interaction == '1Ddouble') then
    ! D, a
    read(1,*)
    read(1,*) parameters(1), parameters(2)
elseif (interaction == 'McKenzie') then
    ! D1, a1, D2, a2, G, g
    read(1,*)
    read(1,*) parameters(1), parameters(2), parameters(3), parameters(4), parameters(5), parameters(6)
elseif (interaction == '2D_Morse') then
    ! D, a, r0
    read(1,*)
    read(1,*) parameters(1), parameters(2), parameters(3)
end if
close(1)
end subroutine


subroutine init_nm(nm_matrix, frequencies, nbead, temperature, freq_type, target_freq, m, nm_masses)
double precision, dimension(nbead, nbead) :: nm_matrix
double precision, dimension(nbead) :: frequencies
integer :: nbead
integer :: i, j
real :: m 
double precision, dimension(nbead) :: nm_masses
real :: temperature
real :: target_freq
character(len = 4) :: freq_type
real(8), parameter :: PI = 4 * atan (1.0_8)

do i = 1,nbead
    nm_matrix(1,i) = 1 / (nbead ** 0.5)
end do 

do j = 1,nbead
    do i = 2,(floor(real(nbead/2)) + 1)
        nm_matrix(i,j) = (2**0.5) * cos(2 * PI * (i - 1) * (j - 1) / nbead) / (nbead ** 0.5)
    end do
    do i = (floor(real(nbead/2)) + 2), nbead
        nm_matrix(i,j) = (2**0.5) * sin(2 * PI * (i - 1) * (j - 1) / nbead) / (nbead ** 0.5)
    end do
end do
if (modulo(nbead, 2) == 0) then
    do i = 1,nbead
        nm_matrix(floor(real(nbead/2)) + 1, i) = (2*modulo(i,2) - 1) / (nbead ** 0.5)
    end do
end if

! NM frequencies, hbar = 1, k_b = 1 are used

if (freq_type == 'rpmd') then
    ! Set to physical frequencies in rpmd case
    do i = 1,nbead
        frequencies(i) = 2*(temperature * nbead) * sin((i - 1) * PI / nbead)
        nm_masses(i) = m
    end do

    do i = 1,nbead
        write(1,*) frequencies(i)
    end do
elseif (freq_type == 'pcmd') then
    ! Set all besides the centroid to specified frequency in the pa-cmd case
    frequencies(1) = 0

    nm_masses(1) = m

    do i = 2,nbead
        frequencies(i) = target_freq
        nm_masses(i) = m * (((2*(temperature * nbead) * sin((i - 1) * PI / nbead)) ** 2) /(target_freq ** 2))
    end do
end if
end subroutine

subroutine calc_forces(q_nm, F, N, nbead, m, interaction, parameters, parameter_number, nm_matrix)
implicit none
double precision, dimension(nbead,3,N) :: q, F, q_nm
double precision, dimension(nbead, nbead) :: nm_matrix
real, dimension(parameter_number) :: parameters
character(len = 8) :: interaction
integer :: N
integer :: nbead
integer :: parameter_number
real :: m

call from_nm_fftw(q, q_nm , nbead, N)

if (interaction == 'harmonic') then
    call calc_forces_harmonic(q,F,N,nbead,m,parameters(1))
elseif (interaction == '1Ddouble') then
    call calc_forces_1Ddouble(q,F,N,nbead,m,parameters(1),parameters(2))
elseif (interaction ==  'McKenzie') then
    call calc_forces_McKenzie(q,F,N,nbead,m,parameters(1),parameters(2),parameters(3),parameters(4),parameters(5),parameters(6))
elseif (interaction ==  '2D_Morse') then
    call calc_forces_2D_Morse(q,F,N,nbead,m,parameters(1),parameters(2),parameters(3))
end if

end subroutine

subroutine calc_PE(q, PE, N, nbead, m, interaction, parameters, parameter_number)
implicit none
double precision, dimension(nbead,3,N) :: q
double precision :: PE
real, dimension(parameter_number) :: parameters
character(len = 8) :: interaction
integer :: N
integer :: nbead
integer :: parameter_number
double precision :: m

if (interaction == 'harmonic') then
    call calc_PE_harmonic(q,PE,N,nbead,m,parameters(1))
elseif (interaction == '1Ddouble') then
    call calc_PE_1Ddouble(q,PE,N,nbead,m,parameters(1),parameters(2))
elseif (interaction == 'McKenzie') then
    call calc_PE_McKenzie(q,PE,N,nbead,m,parameters(1),parameters(2),parameters(3),parameters(4),parameters(5),parameters(6))
elseif (interaction == '2D_Morse') then
    call calc_PE_2D_Morse(q,PE,N,nbead,m,parameters(1),parameters(2),parameters(3))
end if

end subroutine


subroutine do_output(q_nm,p_nm,F,nbead,N,nm_matrix,frequencies,nm_masses,temperature,interaction, &
parameters,parameter_number, out_bead, out_pos, out_mom, out_force)
double precision, dimension(nbead, 3, N) :: q, p, q_nm, p_nm, F, F_nm 
double precision, dimension(nbead, nbead) :: nm_matrix
double precision, dimension(nbead) :: frequencies
double precision :: KE, PE, temp, KE_class, PE_class_aux
real :: temperature
double precision, dimension(nbead) :: nm_masses
integer :: parameter_number
real, dimension(parameter_number) :: parameters
character(len = 8) :: interaction
integer :: nbead
integer :: N 
integer :: l, y, z, j, k
character(len = 3) :: out_bead
character(len = 3) :: out_pos
character(len = 3) :: out_mom
character(len = 3) :: out_force

if (out_pos == 'yay') then
    call from_nm_fftw(q, q_nm , nbead, N)
    ! Write the file headers
    write(2*nbead + 3,*) N
    write(2*nbead + 3,*)
    if (out_bead == 'yay') then
        do l = 1,nbead
            write(2*l + 1,*) N
            write(2*l + 1,*)
        end do
    end if
    do k = 1,N
        if (out_bead == 'yay') then
            do l = 1,nbead
                ! Write the bead positions to the output files
                write(2*l + 1,*) 'X', q(l,1,k), q(l,2,k), q(l,3,k)
            end do
        end if
        ! Write the centroid positions to the output file
        write(2*nbead + 3,*) 'X', q_nm(1,1,k) / (nbead ** 0.5), q_nm(1,2,k) / (nbead ** 0.5), q_nm(1,3,k) / (nbead ** 0.5)
    end do
end if
if (out_mom == 'yay') then
    call from_nm_fftw(p, p_nm , nbead, N)
    write(2*nbead + 4,*) N
    write(2*nbead + 4,*)
    if (out_bead == 'yay') then
        do l = 1,nbead
            write(2*l + 2,*) N 
            write(2*l + 2,*)
        end do
    end if
    do k = 1,N
        if (out_bead == 'yay') then
            do l = 1,nbead
                write(2*l + 2,*) 'X', p(l,1,k), p(l,2,k), p(l,3,k)
            end do
        end if
        write(2*nbead + 4,*) 'X', p_nm(1,1,k) / (nbead ** 0.5), p_nm(1,2,k) / (nbead ** 0.5), p_nm(1,3,k) / (nbead ** 0.5)
    end do
end if
if (out_force == 'yay') then
    call to_nm_fftw(F, F_nm , nbead, N)
    write(4*nbead + 5,*) N 
    write(4*nbead + 5,*)
    if (out_bead == 'yay') then
        do l = 1,nbead
            write(l + 3*nbead + 4,*) N 
            write(l + 3*nbead + 4,*)
        end do
    end if
    do k = 1,N
        if (out_bead == 'yay') then
            do l = 1,nbead
                write(l + 3*nbead + 4,*) 'X', F(l,1,k), F(l,2,k), F(l,3,k)
            end do
        end if
        write(4*nbead + 5,*) 'X', F_nm(1,1,k) / (nbead ** 0.5), F_nm(1,2,k) / (nbead ** 0.5), F_nm(1,3,k) / (nbead ** 0.5)
    end do
end if

call calc_KE_classical(KE_class,p_nm,nm_masses,N,nbead)
call calc_PE_classical_auxiliary(PE_class_aux,q_nm,frequencies,nm_masses,N,nbead)
call calc_KE(q,F,KE,N,nbead,temperature)
call calc_PE(q,PE,N,nbead,nm_masses(1),interaction,parameters,parameter_number)
temp = (2 * KE_class / (3 * N * nbead))
write(2,*) KE, PE, KE + PE, temp
write(4*nbead+6,*) KE_class, PE_class_aux + PE*nbead, KE_class + PE_class_aux + PE*nbead
end subroutine

subroutine calc_forces_harmonic(q,F,N,nbead,m,omega)
implicit none
double precision, dimension(nbead,3,N) :: q, F
integer :: N
integer :: nbead
integer :: i 
integer :: k
integer :: j
real :: omega
real :: m
do i = 1,nbead
    do j = 1,3
        do k = 1,N
            F(i,j,k) = - m * (omega**2) * q(i,j,k)
        end do
    end do
end do
end subroutine

subroutine calc_PE_harmonic(q,PE,N,nbead,m,omega)
implicit none
double precision, dimension(nbead,3,N) :: q
double precision :: PE
integer :: N
integer :: nbead
integer :: i 
integer :: k
real :: omega
double precision :: m

PE = 0

do i = 1,nbead
    do k = 1,N
        PE = PE + 0.5 * m * (omega**2) * (q(i,1,k)**2 + q(i,2,k)**2 + q(i,3,k)**2)
    end do
end do

PE = PE / nbead

end subroutine

subroutine calc_forces_1Ddouble(q,F,N,nbead,m,D,a)
implicit none
double precision, dimension(nbead,3,N) :: q, F
integer :: N
integer :: nbead
integer :: i, l, j
real :: D, a, m

!
! Each "particle" represents in reality three decoupled DOFs each moving in
! a 1D double well potential.
!

do i = 1,nbead
    do j = 1,3
        do l = 1,N
            F(i,j,l) = - D * q(i,j,l) * ((q(i,j,l) ** 2) - a ** 2)
        end do
    end do
end do


end subroutine

subroutine calc_PE_1Ddouble(q,PE,N,nbead,m,D,a)
implicit none
double precision, dimension(nbead,3,N) :: q
double precision :: PE
integer :: N 
integer :: nbead
integer :: i, j, k
double precision :: m
real :: D, a

!
! Each "particle" represents in reality three decoupled DOFs each moving in
! a 1D double well potential.
!

PE = 0

do i = 1,nbead
    do j = 1,N
        do k = 1,3
            PE = PE + 0.25 * D * (q(i,k,j) ** 2 - a ** 2) ** 2
        end do
    end do
end do

PE = PE / nbead

end subroutine

subroutine calc_forces_McKenzie(q,F,N,nbead,m,D1,a1,D2,a2,g1,g2)
    implicit none
    double precision, dimension(nbead,3,N) :: q, F
    integer :: N
    integer :: nbead
    integer :: i, l
    real :: D1, a1, D2, a2, g1, g2, m
    
    !
    ! XY of particles move in a potential from https://doi.org/10.1063/1.2785186
    ! the z component is just a confining harmonic potential
    !
    
    do i = 1,nbead
        do l = 1,N
            F(i,1,l) = - 2 * D1 * q(i,1,l) * (q(i,1,l) ** 2 - a1 ** 2) + 4 * g1 * q(i,2,l) - &
             3 * g2 * q(i,2,l) * q(i,1,l) ** 2 - g2 * q(i,2,l) ** 3
            F(i,2,l) = - 2 * D2 * q(i,2,l) * (q(i,2,l) ** 2 - a2 ** 2) + 4 * g1 * q(i,1,l) - &
            3 * g2 * q(i,1,l) * q(i,2,l) ** 2 - g2 * q(i,1,l) ** 3
            F(i,3,l) = - q(i,3,l)
        end do
    end do
    
    
end subroutine

subroutine calc_PE_McKenzie(q,PE,N,nbead,m,D1,a1,D2,a2,g1,g2)
    implicit none
    double precision, dimension(nbead,3,N) :: q
    double precision :: PE
    integer :: N 
    integer :: nbead
    integer :: i, j, k
    double precision :: m
    real :: D1, a1, d2, a2, g1, g2
    
    !
    ! XY of particles move in a potential from https://doi.org/10.1063/1.2785186
    ! the z component is just a confining harmonic potential
    !
    
    PE = 0
    
    do i = 1,nbead
        do j = 1,N
            PE = PE + 0.5 * (D1 * (q(i,1,j) ** 2 - a1 ** 2) ** 2 + D2 * (q(i,2,j) ** 2 - a2 ** 2) ** 2) - &
            q(i,1,j) * q(i,2,j) * (4 * g1 - g2 * (q(i,1,j) ** 2 + q(i,2,j) ** 2)) + 0.5 * q(i,3,j) ** 2
        end do
    end do
    
    PE = PE / nbead
    
end subroutine

subroutine calc_forces_2D_Morse(q,F,N,nbead,m,D,a,r0)
    implicit none
    double precision, dimension(nbead,3,N) :: q, F
    integer :: N
    integer :: nbead
    integer :: i, l
    real ::D, a, r0, m
    
    !
    ! A simple Morse potential in the usual form in 2D dimensions
    ! the z component is just a confining harmonic potential
    !
    
    do i = 1,nbead
        do l = 1,N
            F(i,1,l) = - D * a * (q(i,1,l) / ((q(i,1,l) ** 2 + q(i,2,l) ** 2) ** 0.5)) * &
             (1 - exp(- a * ((q(i,1,l) ** 2 + q(i,2,l) ** 2) ** 0.5 - r0))) * &
             exp(- a * ((q(i,1,l) ** 2 + q(i,2,l) ** 2) ** 0.5 - r0))
            F(i,2,l) = - D * a * (q(i,2,l) / ((q(i,1,l) ** 2 + q(i,2,l) ** 2) ** 0.5)) * &
             (1 - exp(- a * ((q(i,1,l) ** 2 + q(i,2,l) ** 2) ** 0.5 - r0))) * &
             exp(- a * ((q(i,1,l) ** 2 + q(i,2,l) ** 2) ** 0.5 - r0))
            F(i,3,l) = - q(i,3,l)
        end do
    end do
    
    
end subroutine

subroutine calc_PE_2D_Morse(q,PE,N,nbead,m,D,a,r0)
    implicit none
    double precision, dimension(nbead,3,N) :: q
    double precision :: PE
    integer :: N 
    integer :: nbead
    integer :: i, j, k
    double precision :: m
    real :: D, a, r0
    
    !
    ! A simple Morse potential in the usual form in 2D dimensions
    ! the z component is just a confining harmonic potential
    !
    
    PE = 0
    
    do i = 1,nbead
        do j = 1,N
            PE = PE + D * (1 - exp(- a * ((q(i,1,j) ** 2 + q(i,2,j) ** 2) ** 0.5 - r0))) ** 2 + &
            0.5 * q(i,3,j) ** 2
        end do
    end do
    
    PE = PE / nbead
    
end subroutine

subroutine calc_KE(q,F,KE,N,nbead,temperature)
implicit none
double precision, dimension(nbead,3,N) :: q, F
double precision, dimension(3,N) :: q_centr
double precision :: KE
real :: temperature
integer :: N
integer :: nbead
integer :: i 
integer :: k
integer :: j

do i = 1,N
    do j = 1,3
        q_centr(j,i) = 0
        do k = 1,nbead
            q_centr(j,i) = q_centr(j,i) + q(k,j,i)
        end do
        q_centr = q_centr / nbead
    end do
end do

KE = N * temperature / 2

do i = 1,N
    do j = 1,3
        do k = 1,nbead
            KE = KE - (0.5 / nbead) * (q(k,j,i) - q_centr(j,i)) * F(k,j,i)
        end do
    end do
end do

end subroutine

subroutine calc_KE_classical(KE,p,nm_masses,N,nbead)
double precision, dimension(nbead,3,N) :: p
double precision :: KE
double precision, dimension(nbead) :: nm_masses
integer :: N 
integer :: k
integer :: j
KE = 0
do k = 1,N
    do j = 1,nbead
        KE = KE + (p(j,1,k)**2 + p(j,2,k)**2 + p(j,3,k)**2)/(2*nm_masses(j))
    end do
end do
end subroutine

subroutine calc_PE_classical_auxiliary(PE,q,frequencies,nm_masses,N,nbead)
!Takes normal mode coordinates as an input
double precision, dimension(nbead,3,N) :: q
double precision, dimension(nbead) :: frequencies
double precision :: PE
double precision, dimension(nbead) :: nm_masses 
integer :: N 
integer :: nbead
integer :: i 
integer :: j 

PE = 0

do i = 1,N
    do j = 1,nbead
        PE = PE + 0.5 * ((frequencies(j)) ** 2) * nm_masses(j) * (q(j,1,i) ** 2 + q(j,2,i) ** 2 + q(j,3,i) ** 2)
    end do
end do

end subroutine

! Possible support for LJ interactions

!subroutine calc_forces_LJ()
!end subroutine

subroutine momentum_step(p_nm,F,N,nbead,k,l,nm_matrix,tau)
implicit none
double precision, dimension(nbead,3,N) :: p, p_nm , F
double precision, dimension(nbead, nbead) :: nm_matrix
integer :: N
integer :: nbead
integer :: k
integer :: j
integer :: l
real :: tau

call from_nm_fftw(p, p_nm , nbead, N)

do k = 1,N
    do j = 1,3
        do l = 1,nbead
            p(l,j,k) = p(l,j,k) + F(l,j,k)*tau
        end do 
    end do
end do

call to_nm_fftw(p, p_nm , nbead, N)

end subroutine

subroutine replica_step(q_nm,p_nm,N,nbead,frequencies,nm_masses,tau,centroid_constraint)
implicit none
double precision, dimension(nbead,3,N) ::q_nm, p_nm
double precision, dimension(nbead) :: frequencies
double precision :: q_nm_temp
double precision :: p_nm_temp
integer :: N
integer :: nbead
integer :: i
integer :: k
integer :: j
integer :: l
real :: tau
double precision, dimension(nbead) :: nm_masses
character(len = 3) :: centroid_constraint
do k = 1,3
    do j = 1,N
        ! Centroid update
        if (centroid_constraint == 'nay') then
            q_nm_temp = q_nm(1,k,j)
            p_nm_temp = p_nm(1,k,j)
            q_nm(1,k,j) = (1 / nm_masses(1)) * p_nm_temp * tau + q_nm_temp
        end if
        do i = 2,nbead
            q_nm_temp = q_nm(i,k,j)
            p_nm_temp = p_nm(i,k,j)
            p_nm(i,k,j) = cos(frequencies(i) * tau) * p_nm_temp - nm_masses(i) * frequencies(i) &
             * sin(frequencies(i) * tau) * q_nm_temp
            q_nm(i,k,j) = (1 / (nm_masses(i) * frequencies(i))) * sin(frequencies(i) * tau) * p_nm_temp
            q_nm(i,k,j) = q_nm(i,k,j) + cos(frequencies(i) * tau) * q_nm_temp
        end do 
    end do
end do
end subroutine

subroutine thermostat(p_nm, temperature, gamma, scale, N, nbead, nm_masses, frequencies, tau)
! Langevin
double precision, dimension(nbead,3,N) :: p_nm
double precision, dimension(nbead,3,2*N) :: init
double precision, dimension(nbead) :: frequencies
double precision :: KE
real :: gamma
real :: scale
double precision, dimension(nbead) :: nm_masses
real :: tau
real :: temperature
integer :: N
integer :: nbead
double precision :: A 
double precision :: B
double precision :: alpha
double precision :: beta
logical :: local
real(8), parameter :: PI = 4 * atan (1.0_8)

local = .TRUE.

call random_number(init)

do i = 2,nbead
    A = exp(-2 * scale * frequencies(i) * tau)
    B = (1 - A ** 2) ** 0.5
    do k = 1,N
        do j = 1,3
            p_nm(i,j,k) = p_nm(i,j,k) * A + ((nm_masses(i) * temperature * nbead) ** 0.5) &
             * B * ((- 2 * log(init(i,j,k))) ** 0.5) * cos(2 * PI * init(i,j, k + N))
        end do
    end do
end do

if (local .eqv. .TRUE.) then
    A = exp(- gamma * tau)
    B = (1 - A ** 2) ** 0.5
    do k = 1,N
        do j = 1,3
            p_nm(1,j,k) = p_nm(1,j,k) * A + ((nm_masses(1) * temperature * nbead) ** 0.5) &
             * B * ((- 2 * log(init(1,j,k))) ** 0.5) * cos(2 * PI * init(1,j, k + N))
        end do
    end do
else

    ! A global CSVR thermostat for the centroid degs of freedom. Included for completeness.

    KE = 0

    do k = 1,N
        do j = 1,3
            KE = KE + (p_nm(1,j,k) ** 2) / (2 * m)
        end do
    end do

    A = exp(- 2 * gamma * tau)

    B = 0

    do k = 1,3
        do j = 1,N
            B = B + (((- 2 * log(init(1,j,k))) ** 0.5) * cos(2 * PI * init(1,j, k + N)) ** 2)
        end do 
    end do

    alpha = (A + ((1 - A) * (temperature * nbead) * B / (2 * KE)) + 2 * ((- 2 * log(init(1,1,1))) ** 0.5) &
     * cos(2 * PI * init(1,1, 1 + N)) * ((A * (nbead * temperature) * (1 - A) / (2 * KE)) ** 0.5)) ** 0.5

    beta = ((- 2 * log(init(1,1,1))) ** 0.5) * cos(2 * PI * init(1,1, 1 + N)) &
     + ((2 * KE * A / ((1 - A) * temperature * nbead)) ** 0.5)

    alpha = sign(alpha, beta)

    do k = 1,3
        do j = 1,N
            p_nm(1,k,j) = alpha * p_nm(1, k, j)
        end do 
    end do

end if


end subroutine







