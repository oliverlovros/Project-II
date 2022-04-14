! This module contains 4 subroutines:
!
!   1. Initialize MPI
!   2. Finalize MPI
!   3. Create the array that assigns each processor its particles
!   4. Create the array that assigns each processor its interactions
!
!   And the global variables needed to run the program in parallel (see below).
!   This module must be "used" in the rest of the code files (main program and modules).

module parallel_module

    use MPI
    implicit none
    ! global variables
    integer :: nproc                                    ! number of processors
    integer :: rank                                     ! processor number (from 0 up to nproc-1)
    integer :: ierror                                   ! control variables
    integer, allocatable :: particles(:,:)              ! array that contains the particles range of each processor
    integer, allocatable :: interactions(:,:)           ! array that contains the interactions range of each processor
    integer, allocatable :: particles_interaction(:,:)  ! array that contains the two particles of each interaction
    integer, allocatable :: displac(:)                  ! displacement of each processor (in particles array)
    integer, allocatable :: num_send(:)                 ! particles of each processor

    contains

    ! 1. Subroutine that initializes MPI.
    subroutine initialize_MPI()

        ! initialize MPI
        call MPI_INIT(ierror)
        ! Number of processors available
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierror)
        ! Identify processors
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

        return

    end subroutine initialize_MPI

    ! 2. Subroutine that finalizes MPI.
    subroutine finalize_MPI()

        ! Finalize MPI
        call MPI_FINALIZE(ierror)

        return

    end subroutine finalize_MPI


    ! 3. Subroutine that assigns each processor its particles.
    subroutine assign_particles(npart)

        implicit none
        integer, intent(in)  :: npart
        integer              :: np, res, i

        ! allocate "particles"
        allocate(particles(0:nproc-1,2), displac(0:nproc-1), num_send(0:nproc-1))

        ! number of particles/processor
        np  = int(npart/nproc)
        ! remaining particles
        res = npart - np*nproc

        ! first processor
        particles(0,1) = 1
        if (res == 0) then
            particles(0,2) = np
        else
            particles(0,2) = np + 1
            res = res - 1
        end if
        ! the rest of processors
        do i = 1, nproc-1
            particles(i,1) = particles(i-1,2) + 1
            if (res == 0) then
                particles(i,2) = particles(i,1) + np - 1
            else
                particles(i,2) = particles(i,1) + np
                res = res - 1
            end if
        end do

        ! displacement and number of particles

        ! first processor
        displac(0) = 0
        num_send(0) = particles(0,2)-particles(0,1)+1
        ! the rest of processors
        do i = 1, nproc-1
            displac(i) = particles(i-1,2)
            num_send(i) = particles(i,2)-particles(i,1)+1
        end do

        return

    end subroutine assign_particles

    ! 4. Subroutine that assigns each processor its interactions and generates the array that assigns two particles to each interaction.
    subroutine assign_interactions(npart)

        implicit none
        integer, intent(in)  :: npart
        integer              :: nint, np, res, i, j, k

        ! allocate "interactions" and "particles_interaction"
        allocate(interactions(0:nproc-1,2), particles_interaction(npart*(npart-1)/2,2))

        ! number of interactions
        nint = npart*(npart-1)/2
        ! number of interactions/processor
        np  = int(nint/nproc)
        ! remaining particles
        res = nint - np*nproc

        interactions(0,1) = 1
        if (res == 0) then
            interactions(0,2) = np
        else
            interactions(0,2) = np + 1
            res = res - 1
        end if

        do i = 1, nproc-1
            interactions(i,1) = interactions(i-1,2) + 1
            if (res == 0) then
                interactions(i,2) = interactions(i,1) + np - 1
            else
                interactions(i,2) = interactions(i,1) + np
                res = res - 1
            end if
        end do

        k = 0
        do i = 1, npart-1
            do j = i+1, npart
                k = k + 1
                particles_interaction(k,1) = i
                particles_interaction(k,2) = j
            end do
        end do

        return

    end subroutine assign_interactions


end module parallel_module