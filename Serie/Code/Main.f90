program Main

    ! modules
    use Constants_module
    use Initial_state_module
    use Pbc_module
    use Forces_module
    use Integrator_module
    use pair_distribution_function_module

    implicit none
    ! files
    character(len=80) :: infile, outfile1, outfile2, outfile3, outfile4
    ! parameters
    integer :: nparts
    double precision :: in_rho            ! density of the system in g/cm^3
    double precision :: mass              ! molecular mass of the system in g/mol
    double precision :: LJ_sig            ! Lennard-Jones length parameter in A
    double precision :: LJ_eps            ! Lennard-Jones energy parameter in KJ/mol
    double precision :: cutoff            ! radius of cutoff divided by the length of the cell
    double precision :: in_temp           ! room temperature in K
    double precision :: in_temp_init      ! initial temperature of the system in K
    character(len=3) :: bimodal           ! if "Yes" the system initialize velocities using a bimodal distribution associated with 
    !                                       the initial temperature 
    character(len=3) :: melting           ! if "Yes" the system realizes init_steps steps at initial temperature
    character(len=3) :: thermostat        ! if "Yes" the system uses a thermostat at room temperature
    character(len=6) :: integrator        ! the integrators can be Euler or Velocity Verlet
    double precision :: in_dt             ! time between two steps in ps
    integer          :: init_steps        ! initial steps in case of disordered system
    integer          :: steps             ! steps of the simulation
    integer          :: measure_steps     ! steps between two measures
    integer          :: traj_steps        ! steps between two trajectories
    integer          :: boxes             ! number of boxes of the histogram
    ! parameters in reduced units
    double precision :: rho, temp, temp_init, dt
    ! observables
    double precision :: time              ! time evolution in reduced units
    double precision :: kin               ! kinetic energy per particle in reduced units
    double precision :: pot               ! potential energy per particle in reduced units
    double precision :: tot               ! total energy per particle in reduced units
    double precision :: tempi             ! temperature of the system in reduced units
    double precision :: pres              ! pressure in reduced units
    double precision :: dr2               ! square displacement in reduced units
    ! conversion factors
    double precision :: lconv, econv, rhoconv, tempconv, pconv, timeconv
    ! simulation box
    integer :: box_m
    double precision :: box_l, box_a
    ! dynamics
    double precision, allocatable :: newpos(:,:), oldpos(:,:), vel(:,:), force(:,:)
    ! pair distribution function
    integer :: measures 
    double precision :: dr
    double precision, allocatable :: gofr(:)
    ! others
    double precision :: dnparts
    double precision :: sigma, nu
    integer :: iter, l, therm_on, integrator_num

    ! read parameters
    infile = "parameters.txt"
    call Read_parameters(infile, nparts, in_rho, mass, LJ_sig, LJ_eps, cutoff, in_temp, in_temp_init, bimodal, &
    melting, thermostat, integrator, in_dt, init_steps, steps, measure_steps, traj_steps, boxes,  outfile1, outfile2, &
    outfile3, outfile4)

    ! print information on the screen
    print*, " "
    print*, "******************************************************************************"
    print*, "                     Molecular Dynamics Simulation"
    print*, "          System of Particles with Lennard-Jones Interaction"
    print*, "******************************************************************************"
    print*, " "
    print"(a28,x,i18)", " Number of particles:          ", nparts
    print"(a28,x,f18.12)", " Density (g/cm^3):          ", in_rho
    print"(a28,x,f18.12)", " Well depth (KJ/mol):       ", LJ_eps
    print"(a28,x,f18.12)", " Characteristic length (A): ", LJ_sig
    if (thermostat == "Yes") print"(a28,x,f18.12)", " Room Temperature (K):       ", in_temp
    if (bimodal == "Yes" .or. melting == "Yes") print"(a28,x,f18.12)", " Initial Temperature (K):    ", in_temp_init
    print"(a28,x,a18)", " Integrator:                ", integrator
    print"(a28,x,f18.12)", " Time step (ps):            ", in_dt
    print"(a28,x,i18)", " Steps:                     ", steps
    print*, " "
    print*, "The simulation starts:"
    print*, " "
    

    ! open files
    open(12,file=outfile1)  ! observables evolution
    open(13,file=outfile2)  ! mean square displacement
    open(36,file=outfile3)  ! system evolution

    ! determine vectors dimensions
    allocate(newpos(nparts,3), oldpos(nparts,3), vel(nparts,3), force(nparts,3), gofr(boxes))

    ! calculate conversion fators
    dnparts = dble(nparts)
    lconv = LJ_sig
    econv = LJ_eps/dnparts
    rhoconv = 1d24*mass/(avogadro*LJ_sig**3)
    tempconv = LJ_eps/(avogadro*boltzmann)
    pconv = 1d33*LJ_eps/(avogadro*LJ_sig**3d0)
    timeconv = 0.1d0*(mass*LJ_sig**2d0/LJ_eps)**0.5d0

    ! calculate parameters in reduced units
    rho = in_rho/rhoconv
    temp = in_temp/tempconv
    temp_init = in_temp_init/tempconv
    dt = in_dt/timeconv

    ! -------------------------------------------------------- prepare the system -------------------------------------------------!
    ! initialize SC geometry
    call sc_lattice(nparts,rho,newpos,box_l,box_a,box_m)

    if (thermostat == "Yes") then 
        therm_on = 1
    else
        therm_on = 0
    endif
    
    ! initialize velocities
    if (bimodal == "Yes") then
        call bimodal_vel(nparts,temp_init,vel)
    else
        vel = 0d0
    endif

    ! initialize the histogram (pair distribution function)
    call pair_distribution_function(1,nparts,measures,rho,newpos,boxes,gofr,box_l,cutoff,dr,outfile4)

    ! determine the integrator
    if (integrator == "Verlet") then
        integrator_num = 1
    else ! Euler
        integrator_num = 0
    endif

    ! Melting
    if (melting == "Yes") then
        print*, "Melting the system at initial Temperature..."
        call LJ_potential(nparts,newpos,cutoff,box_l,1, pot, force)
        nu = 5.d0*dt
        sigma = dsqrt(temp_init)
        do iter  = 1, init_steps
            if (integrator_num == 1) then 
                call velocity_verlet_andersen(nparts,box_l,cutoff,nu,sigma,dt,newpos,vel,pot,force,1)
            else
                call euler_andersen(nparts,box_l,cutoff,nu,sigma,dt,newpos,vel,pot,force,1)
            endif
        enddo
    endif
    print*, "Computing the dynamics..."
    !------------------------------------------------------------- the system is ready --------------------------------------------!
    ! Initial state
    write(12,*) '# N = ',nparts,', rho (g/cm^3) = ',in_rho,', L (A) = ',box_L*lconv,', a (A) = ',box_a*lconv,', M = ',box_m
    write(12,*) '# time (ps), kin/N (kJ/mol), pot/N (kJ/mol), tot/N (kJ/mol), temp (K), pres (Pa)'
    time = 0.d0
    call LJ_potential(nparts,newpos,cutoff,box_l,1, pot, force)
    call Kinetic_Energy(nparts, vel, kin)
    tot = kin + pot
    tempi = kin*2.d0/(3.d0*dnparts-3.d0)
    call pressure(nparts,rho,tempi,box_l,cutoff,newpos,pres)
    write(12,*) time*timeconv, kin*econv, pot*econv, tot*econv, tempi*tempconv, pres*pconv
    ! time evolution (dynamics)
    oldpos = newpos
    dr2 = 0.d0
    nu = 5.d0*dt
    sigma = dsqrt(temp)
    do iter  = 1, steps
        ! next point
        time = time + dt
        if (integrator_num == 1) then
            call velocity_verlet_andersen(nparts,box_l,cutoff,nu,sigma,dt,newpos,vel,pot,force,therm_on)
        else
            call euler_andersen(nparts,box_l,cutoff,nu,sigma,dt,newpos,vel,pot,force,therm_on)
        endif
        ! mean square displacement
        dr2 = dr2 + sum((newpos - oldpos)**2)
        if (mod(iter,100) == 0) then
            write(13,*) time*timeconv, dr2*lconv**2/dble(nparts)/100.d0
            oldpos = newpos
        end if
        if (mod(iter,measure_steps) == 0) then

            ! sample distances (pair distribution function)
            call pair_distribution_function(2,nparts,measures,rho,newpos,boxes,gofr,box_l,cutoff,dr,outfile4)

            ! calculate energy, temperature and pressure
            call Kinetic_Energy(nparts, vel, kin)
            tot = kin + pot
            tempi = kin*2.d0/(3.d0*dnparts-3.d0)
            call pressure(nparts,rho,tempi,box_l,cutoff,newpos,pres)
            write(12,*) time*timeconv, kin*econv, pot*econv, tot*econv, tempi*tempconv, pres*pconv
        end if

        ! trajectories
        if (mod(iter,traj_steps) == 0) then

            ! Write trajectory in xyz format to be read by VMD
            write(36,*) nparts
            write(36,*)
            do l=1,nparts
                write(36,*) 'A', newpos(l,1)*lconv, newpos(l,2)*lconv, newpos(l,3)*lconv
            end do

        end if
    end do

    ! normalize histogram (pair distribution function)
    call pair_distribution_function(3,nparts,measures,rho,newpos,boxes,gofr,box_l,cutoff,dr,outfile4)
    ! write results in real units
    dr = dr*lconv
    call pair_distribution_function(4,nparts,measures,rho,newpos,boxes,gofr,box_l,cutoff,dr,outfile4)

    ! close files

    close(12)
    close(13)
    close(36)

    print*, " "
    print*, "******************************************************************************"
    print*, "                          End of Simulation"
    print*, "******************************************************************************"

end program Main