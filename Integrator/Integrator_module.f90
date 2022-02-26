module Integrator_module
    
    ! parameters
    implicit none

    ! functions and subroutines
    contains
    
    subroutine velocity_verlet(nparts,box_L,cutoff,dt,positions,velocities,pot_E,forces)

     !-----------------------------------------------------------------------------------------------------------------------------!
     ! Codi escrit per Alex Teruel
     ! Informació
     ! La subrutina aplica el algorisme de velocity Verlet per integrar un pas de temps per un sistema de partícules interaccionant
     ! sota un potencial de Lennard-Jones.
     ! Variables d'entrada:
     !   nparts: número de partícules del sistema
     !   box_L: longitud de cada costat de la caixa
     !   cutoff: radi a partir del qual es negligeixen les interaccions 
     !   dt: el salt temporal
     ! Variables d'entrada i sortida (modificades):
     !   positions(npart,3): matriu que conte les posicions de les npart partícules
     !   velocities(nparts,3): matriu que conte la velocitat en cada direcció de cada partícula
     !   pot_E: energia potencial del sistema
     !   forces(nparts,3): matriu que conte la força en cada dimensió que rep cada partícula
     !-----------------------------------------------------------------------------------------------------------------------------!
    
        implicit none
        ! variables d'entrada i sortida
        integer :: nparts
        double precision :: box_L, cutoff, dt
        double precision :: positions(nparts,3), velocities(nparts,3)
        double precision :: pot_E, forces(nparts,3)
        ! variables internes subrutina
        double precision :: newforces(nparts,3)
    
        ! calculem noves posicions
        positions = positions + velocities*dt + 0.5d0*forces*dt**2
        ! apliquem les condicions de contorn periodiques
        call pbc_nparts(nparts,positions,box_L)
        ! calculem noves forces i energia potencial
        call LJ_potential(nparts,positions,cutoff,box_L,1,pot_E,newforces)
        ! calculem noves velocitats
        velocities = velocities + 0.5d0*dt*(forces + newforces)
        ! reasignem valors de les forces
        forces = newforces
    
        return
    
    end subroutine velocity_verlet

    subroutine velocity_verlet_andersen(nparts,box_L,cutoff,nu,sigma,dt,positions,velocities,pot_E,forces)

     !-----------------------------------------------------------------------------------------------------------------------------!
     ! Codi escrit per Alex Teruel
     ! Informació
     ! La subrutina aplica el algorisme de velocity Verlet per integrar un pas de temps per un sistema de partícules interaccionant
     ! sota un potencial de Lennard-Jones i un termostat Andersen.
     ! Variables d'entrada:
     !   nparts: número de partícules del sistema
     !   box_L: longitud de cada costat de la caixa
     !   cutoff: radi a partir del qual es negligeixen les interaccions
     !   nu: probabilitat d'acceptar un canvi de velocitat
     !   sigma: desviació de la distribució
     !   dt: el salt temporal
     ! Variables d'entrada i sortida (modificades):
     !   positions(npart,3): matriu que conte les posicions de les npart partícules
     !   velocities(nparts,3): matriu que conte la velocitat en cada direcció de cada partícula
     !   pot_E: energia potencial del sistema
     !   forces(nparts,3): matriu que conte la força en cada dimensió que rep cada partícula
     !-----------------------------------------------------------------------------------------------------------------------------!
    
        implicit none
        ! variables d'entrada i sortida
        integer :: nparts
        double precision :: box_L, cutoff, nu, sigma, dt
        double precision :: positions(nparts,3), velocities(nparts,3)
        double precision :: pot_E, forces(nparts,3)
        ! variables internes subrutina
        double precision :: newforces(nparts,3)
    
        ! calculem noves posicions
        positions = positions + velocities*dt + 0.5d0*forces*dt**2
        ! apliquem les condicions de contorn periodiques
        call pbc_nparts(nparts,positions,box_L)
        ! calculem noves forces i energia potencial
        call LJ_potential(nparts,positions,cutoff,box_L,1,pot_E,newforces)
        ! calculem noves velocitats
        velocities = velocities + 0.5d0*dt*(forces + newforces)
        ! reasignem valors de les forces
        forces = newforces
        ! apliquem el termostat
        call therm_Andersen(nparts,velocities,nu,sigma)
    
        return
    
    end subroutine velocity_verlet_andersen

    subroutine euler(nparts,box_L,cutoff,dt,positions,velocities,pot_E,forces)

        implicit none
        integer :: nparts
        double precision :: positions(nparts,3), velocities(nparts,3), forces(nparts,3), dt, box_L, cutoff, pot_E
    
        positions = positions + velocities*dt + 0.5d0*forces*dt**2
        velocities = velocities + forces*dt
        call pbc_nparts(nparts,positions,box_L)
        call LJ_potential(nparts,positions,cutoff,box_L,1,pot_E,forces)

    
    
    end subroutine euler

    subroutine euler_andersen(nparts,box_L,cutoff,nu,sigma,dt,positions,velocities,pot_E,forces)

        implicit none
        integer :: nparts
        double precision :: positions(nparts,3), velocities(nparts,3), forces(nparts,3), dt, box_L, cutoff, pot_E, nu, sigma
    
        positions = positions + velocities*dt + 0.5d0*forces*dt**2
        velocities = velocities + forces*dt
        call pbc_nparts(nparts,positions,box_L)
        call LJ_potential(nparts,positions,cutoff,box_L,1,pot_E,forces)
        call therm_Andersen(nparts,velocities,nu,sigma)
    
    end subroutine euler_andersen

end module Integrator_module