!---------------------------------------------------------Initial state module-----------------------------------------------------!
module Initial_state_module

    implicit none
    
    contains

    subroutine sc_lattice(npart,density,positions,length,space,celdim)

        !--------------------------------------------------------------------------------------------------------------------------!
        ! Author: Oliver Loveday
        ! About:
        ! The sc_lattice subroutine generates a crystal system (disposition-like) with a simple cubic unit cell, centered at the 
        ! origin of coordinates.
        ! Input variables:
        !   npart: total number of particles in the system
        !   density: density
        ! Output variables:
        !   positions(npart,3): matrix containing the positions of all the particles (npart)
        !   length: length of each side of the box
        !   space: distance between two nodes in the same dimension
        !   celdim: number of nodes in each dimension
        !--------------------------------------------------------------------------------------------------------------------------!
        
        implicit none

        integer, intent(in) :: npart
        double precision, intent(in) :: density
        integer, intent(out) :: celdim
        double precision, intent(out) :: positions(npart,3), length
        double precision :: space
        integer :: i, j, k, numpart
        
        ! Parameters are determined
        length = (dble(npart)/density)**(1.d0/3.d0)
        celdim = int(npart**(1.d0/3.d0)+0.5d0)
        space = length/dble(celdim)

        ! Particles are introduced in each position
        numpart = 0
        do i = 1, celdim
            do j = 1, celdim
                do k = 1, celdim
                    numpart = numpart + 1
                    ! Counter is updated each time a new particle is introduced
                    positions(numpart,1) = -0.5d0*length + dble(i-1)*space
                    positions(numpart,2) = -0.5d0*length + dble(j-1)*space
                    positions(numpart,3) = -0.5d0*length + dble(k-1)*space
    
                end do
            end do
        end do
    
        return
    
    end subroutine sc_lattice

    subroutine bimodal_vel(nparts,temp,vel)

        !--------------------------------------------------------------------------------------------------------------------------!
        ! Author: Oliver Loveday
        ! Information
        ! The bimodal_vel initializes the velocities of a particle system following the bimodal distribution, which is compatible
        ! for a given temperature
        ! Input variables:
        !   nparts: total number of particles in the system
        !   temp: temperature of the system
        ! Input and output variables (modified):
        !   vel(nparts,3): matrix which contains the velocities of the particles for each direction (x, y and z)
        !--------------------------------------------------------------------------------------------------------------------------!
    
        implicit none
        ! Input and output variables
        integer :: nparts
        double precision :: temp, vel(nparts,3)
        ! Other variables
        integer :: ipos, ineg, i, j
        double precision :: vpos, vneg, xi
    
        vpos = dsqrt(temp)
        vneg = -vpos
        ipos = 0
        ineg = 0
        do i = 1, nparts
            call random_number(xi)
            if ((xi < 0.5d0).and.(ipos < nparts/2)) then
                ipos = ipos + 1
                do j = 1, 3
                    vel(i,j) = vpos
                end do
            else if (ineg < nparts/2) then
                ineg = ineg + 1
                do j = 1, 3
                    vel(i,j) = vneg
                end do
            else
                ipos = ipos + 1
                do j = 1, 3
                    vel(i,j) = vpos
                end do
            end if
        end do
    
        return
    
    end subroutine bimodal_vel

    ! read input file
    subroutine Read_parameters(params_file, nparts, density, mass, sigma, epsilon, cutoff, temp_room, temp_init, bimodal,&
        disorder_system, thermostat, integrator, time_step, init_steps, steps, measure_steps, traj_steps, boxes, observables_file, &
        MSD_file, positions_file, gofr_file)

        implicit none
        ! input
        character(len=50), intent(in) :: params_file
        ! output
        integer, intent(out) :: nparts
        double precision, intent(out) :: density ! g/cm^3
        double precision, intent(out) :: mass ! g/mol
        double precision, intent(out) :: sigma ! A
        double precision, intent(out) :: epsilon ! KJ/mol
        double precision, intent(out) :: cutoff ! relative to cell length
        double precision, intent(out) :: temp_room ! K
        double precision, intent(out) :: temp_init ! K
        character(len=3), intent(out) :: bimodal ! Si/No
        character(len=3), intent(out) :: disorder_system ! Yes/No
        character(len=3), intent(out) :: thermostat ! Yes/No
        character(len=6), intent(out) :: integrator ! Euler/Verlet
        double precision, intent(out) :: time_step ! ps
        integer, intent(out) :: init_steps ! initial steps in case of disordered system
        integer, intent(out) :: steps ! steps of simulation
        integer, intent(out) :: measure_steps ! steps between two measures
        integer, intent(out) :: traj_steps ! steps between two trajectories
        integer, intent(out) :: boxes ! number of boxes histogram (g(r))
        character(len=80), intent(out) :: observables_file ! energy, temperature, etc evolution
        character(len=80), intent(out) :: MSD_file ! mean square displacement
        character(len=80), intent(out) :: positions_file ! system evolution 
        character(len=80), intent(out) :: gofr_file ! gofr.f90 output file
        ! others
        integer :: my_unit

        open(newunit=my_unit, file=params_file)
        read(my_unit,*)
        read(my_unit,*)
        read(my_unit,*)
        read(my_unit,*)
        read(my_unit,*)
        read(my_unit,*)
        read(my_unit,*)
        read(my_unit,*)
        read(my_unit,*) 
        read(my_unit,*) nparts
        read(my_unit,*)
        read(my_unit,*) density
        read(my_unit,*) 
        read(my_unit,*) mass
        read(my_unit,*) 
        read(my_unit,*) sigma
        read(my_unit,*) 
        read(my_unit,*) epsilon
        read(my_unit,*)
        read(my_unit,*)
        read(my_unit,*)
        read(my_unit,*) 
        read(my_unit,*) cutoff
        read(my_unit,*)
        read(my_unit,*) temp_room
        read(my_unit,*)
        read(my_unit,*) temp_init
        read(my_unit,*)
        read(my_unit,*) bimodal
        read(my_unit,*) 
        read(my_unit,*) disorder_system
        read(my_unit,*) 
        read(my_unit,*) thermostat
        read(my_unit,*) 
        read(my_unit,*) integrator
        read(my_unit,*) 
        read(my_unit,*) time_step
        read(my_unit,*)
        read(my_unit,*) init_steps
        read(my_unit,*) 
        read(my_unit,*) steps
        read(my_unit,*) 
        read(my_unit,*) measure_steps
        read(my_unit,*)
        read(my_unit,*) traj_steps
        read(my_unit,*)
        read(my_unit,*) boxes
        read(my_unit,*)
        read(my_unit,*)
        read(my_unit,*)
        read(my_unit,*) 
        read(my_unit,*) observables_file
        read(my_unit,*) 
        read(my_unit,*) MSD_file
        read(my_unit,*) 
        read(my_unit,*) positions_file
        read(my_unit,*)
        read(my_unit,*) gofr_file
        close(my_unit)

        return

    end subroutine Read_parameters

end module Initial_state_module
