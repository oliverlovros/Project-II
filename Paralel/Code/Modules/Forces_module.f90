module forces_module

use Constants_module
use pbc_module
use Parallel_module
!--------------------------!
implicit none

contains

subroutine LJ_potential(npart,positions,cutoff,length,pbc_on,Upot,force)

    implicit none
    integer, intent(in) :: npart, pbc_on
    double precision, intent(in) :: positions(npart,3), cutoff, length
    double precision :: Upot, force(npart,3)
    integer :: i, j, k
    integer :: part_1, part_2
    double precision :: dr(3), dr2, diff
    double precision :: fx, fy, fz, ftot(npart,3)

    force = 0.d0
    Upot  = 0.d0
    
    do i=interactions(rank,1),interactions(rank,2)

        part_1=particles_interaction(i,1)
        part_2=particles_interaction(i,2)

        ! difference of positions
        dr(1) = positions(part_1,1) - positions(part_2,1)
        dr(2) = positions(part_1,2) - positions(part_2,2)
        dr(3) = positions(part_1,3) - positions(part_2,3)

        ! PBC are applied
        call pbc(dr,length)

        ! Distance calculation
        dr2 = dr(1)**2 + dr(2)**2 + dr(3)**2

        ! Cutoff correction
        if(dr2 <= cutoff**2) then

            diff = (48.d0/dr2**7 - 24.d0/dr2**4)
            ! force
            fx = diff*dr(1)
            fy = diff*dr(2)
            fz = diff*dr(3)
            ! Force part_1
            force(part_1,1) = force(part_1,1) + fx
            force(part_1,2) = force(part_1,2) + fy
            force(part_1,3) = force(part_1,3) + fz
            ! Force part_2
            force(part_2,1) = force(part_2,1) - fx
            force(part_2,2) = force(part_2,2) - fy
            force(part_2,3) = force(part_2,3) - fz
            ! Potential energy
            Upot = Upot+4.d0*(1.d0/dr2**6 - 1.d0/dr2**3) - 4.d0*( 1.d0/cutoff**12 - 1.d0/cutoff**6)

        end if
    end do
    
    call MPI_allreduce(force,ftot,size(force),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
    call MPI_allreduce(Upot,Upot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)

    force = ftot

    return 
end subroutine LJ_potential

subroutine Kinetic_Energy(nparts,velocity,kin_E)

    !------------------------------------------------------------------------------------------------------------------------------!
    ! Author: Daniel Conde
    ! About:
    ! The Kinetic_Energy subroutine calculates the kinetic energy for a particle system.
    ! Input variables:
    !   nparts: total number of particles in the system
    !   velocity(nparts,3): matrix containing the velocities in each direction for each particle 
    ! Output variables:
    !   kin_E: kinetic energy of the system
    !------------------------------------------------------------------------------------------------------------------------------!

    implicit none
    ! Input and output variables
    integer, intent(in) :: nparts
    double precision, intent(in) :: velocity(nparts,3)
    double precision :: kin_E
    ! Internal variables (subroutine)
    integer :: n, i
    double precision :: vel2
   
    kin_E = 0.d0
    do n = particles(rank,1), particles(rank,2)
        vel2 = 0.d0
        do i = 1, 3
            vel2 = vel2 + velocity(n,i)**2
        end do
        kin_E = kin_E + 0.5d0*vel2
    end do
    call MPI_allreduce(kin_E,kin_E,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)

    return
   
end subroutine Kinetic_Energy



subroutine pressure(n,rho,temp,l,cutoff,positions,pres)
  
    implicit none
    ! Variables
    integer :: n
    double precision :: rho, temp, l, cutoff, positions(n,3)
    double precision :: pres,UPot
    integer :: i, part_1, part_2
    double precision :: dr(3), dr2, force(3)

    pres = 0.d0
    
    do i=interactions(rank,1),interactions(rank,2)

        part_1=particles_interaction(i,1)
        part_2=particles_interaction(i,2)
  
        !dr = positions(part_1,:) - position(part_2,:)
        dr(1) = positions(part_1,1) - positions(part_2,1)
        dr(2) = positions(part_1,2) - positions(part_2,2)
        dr(3) = positions(part_1,3) - positions(part_2,3)
  
        ! PBC are applied
          call pbc(dr,l)
  
        ! Distance calculation
          dr2 = dr(1)**2 + dr(2)**2 + dr(3)**2

        ! Cutoff correction
          if(dr2 <= cutoff**2) then
            force = (48.d0/dr2**7 - 24.d0/dr2**4)*dr
            pres = pres + sum(dr*force)
          end if
    end do

    call MPI_allreduce(pres,pres,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)

    pres=pres+16/3.*pi*rho**2*(2/3.*(1./cutoff)**9)-(1/cutoff)**3/(1/(3*l**3))! Added for cut-off corrections
    pres = rho*temp + pres/(3.d0*l**3)
    
    return

end subroutine pressure


end module forces_module

