program Main

    !------------------------------------------------------------------------------------------------------------------------------!
    ! Codi escrit per Alex Teruel
    ! Informació
    ! El programa simula un sistema de N partícules en contacte amb un bany tèrmic que interaccionen mitjançant un potencial de 
    ! Lennard-Jones.
    ! Les variables d'entrada es llegeixen en un fitxer que el programa demana al executar-se.
    ! El format del fitxer d'entrada ha de ser:
    !   nparts rho temp
    !   dt t_steps k_steps
    !   mass epsilon sigma
    !   outfile
    !   opt1 outfile1
    !   opt2 dr outfile2
    ! Les variables d'entrada són:
    !   nparts: número de partícules
    !   rho: la densitat en unitats reduides
    !   temp: la temperatura ambient en unitats reduides
    !   dt: el salt temporal en unitats reduides
    !   t_steps: els salts temporals totals d'integració
    !   k_steps: els salts temporals entre dos mesures
    !   mass: massa de les partícules en g/mol
    !   LJ_eps: paràmetre energia LJ en kj/mol
    !   LJ_sig: paràmetre longitud LJ en ang
    !   outfile: nom del fitxer on volem guardar els resultats
    !   opt1: 1 si volem calcular el desplaçament quadràtic mitjà, 0 en cas contrari
    !   outfile1: nom del fitxer on volem guardar la evolució temporal del desplaçament quadràtic mitjà
    !   opt2: 1 si volem calcular la funció de distribució radial, 0 en cas contrari
    !   outfile2: nom del fitxer on volem guardar la funció de distribució radial
    ! Les variables de sortida són:
    !   time: el instant de temps en picosegons
    !   kin/nparts: energia cinètica per partícula en kj/mol
    !   pot/nparts: energia potencial per partícula en kj/mol
    !   tot/nparts: energia total per partícula en kj/mol
    !   tempi: temperatura instantània en k
    !   pres: pressió en Pa
    !   dr2: desplaçament quadràtic mitjà en ang^2
    !   dr: distància entre partícules ang
    !   gofr: funcio de distribució radial normalitzada
    !------------------------------------------------------------------------------------------------------------------------------!

    use Constants_module
    use Initial_state_module
    use Pbc_module
    use Forces_module
    use Integrator_module
    implicit none
    ! fitxers
    integer :: opt1, opt2
    character(len=80) :: infile, outfile, outfile1, outfile2
    ! variables d'entrada
    integer :: nparts
    double precision :: rho, temp, temp_init
    integer :: t_steps, k_steps
    double precision :: dt
    double precision :: mass, LJ_eps, LJ_sig
    character(len=3) :: geometry
    character(len=2) :: bimodal, disorder_system, thermostat
    character(len=6) :: integrator
    ! variables de sortida
    double precision :: time, kin, pot, tot, tempi, pres
    double precision:: dr2
    ! altres variables
    integer :: box_m
    double precision :: box_l, box_a
    double precision :: dnparts, lconv, econv, rhoconv, tempconv, pconv, timeconv
    double precision, allocatable :: positions(:,:), velocities(:,:), forces(:,:)
    double precision, allocatable :: oldpos(:,:)
    double precision :: cutoff, sigma, nu
    integer :: i,l

    ! leer parametros del fichero de entrada
    infile = "parametros.dat"
    call Read_parameters(infile, nparts, geometry, rho, mass, LJ_sig, LJ_eps, cutoff, temp, temp_init, bimodal, &
    disorder_system, thermostat, integrator, dt, t_steps, k_steps, outfile, outfile1, outfile2)

    ! ficheros de salida
    open(12,file=outfile)
    open(13,file=outfile1)
    open(15,file=outfile2)
    !Escribe tray para poder leerse en vmd
    open(36,file='vmd.xyz')

    ! vectores
    allocate(positions(nparts,3), oldpos(nparts,3), velocities(nparts,3), forces(nparts,3))

    ! conversion de unidades reducidas a reales
    dnparts = dble(nparts)
    lconv = LJ_sig
    econv = LJ_eps/dnparts
    rhoconv = 1d24*mass/(navo*LJ_sig**3)
    tempconv = LJ_eps/(navo*kbol)
    pconv = 1d33*LJ_eps/(navo*LJ_sig**3d0)
    timeconv = 0.1d0*(mass*LJ_sig**2d0/LJ_eps)**0.5d0

    ! radio de cutoff
    cutoff = cutoff*box_l

    ! geometria
    if (geometry == "SC ") call sc_lattice(nparts,rho,positions,box_l,box_a,box_m)
    !if (geometry == "FCC") call fcc_lattice(nparts,rho,positions,box_l,box_a,box_m)

    ! iniciamos velocidades
    if (bimodal == "Si") then
        call bimodal_vel(nparts,temp_init,velocities)
    else
        velocities = 0d0
    endif

    ! sistema desordenado
    if (disorder_system == "Si") then
        call LJ_potential(nparts,positions,cutoff,box_l,1, pot, forces)
        nu = 5.d0*dt
        sigma = dsqrt(temp_init)
        do i  = 1, 1000000
            call velocity_verlet_andersen(nparts,box_l,cutoff,nu,sigma,dt,positions,velocities,pot,forces)
        enddo
    endif
    
    ! fuerza, energia, presion y temperatura inicial
    write(12,*) '# N = ',nparts,', rho (g/cm^3) = ',rho*rhoconv,', L (A) = ',box_L*lconv,', a (A) = ',box_a*lconv,', M = ',box_m
    write(12,*) '# time (ps), kin/N (kJ/mol), pot/N (kJ/mol), tot/N (kJ/mol), temp (K), pres (Pa)'
    time = 0.d0
    call LJ_potential(nparts,positions,cutoff,box_l,1, pot, forces)
    call Kinetic_Energy(nparts, velocities, kin)
    tot = kin + pot
    tempi = kin*2.d0/(3.d0*dnparts-3.d0)
    call pressure(nparts,rho,tempi,box_l,cutoff,positions,pres)
    write(12,*) time*timeconv, kin*econv, pot*econv, tot*econv, tempi*tempconv, pres*pconv
    ! bucle temporal
    nu = 5.d0*dt
    sigma = dsqrt(temp)
    dr2 = 0.d0
    oldpos = positions
    do i  = 1, t_steps
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Escribir trayectoria en xyz para leerse en VMD!
    write(36,*) nparts
    write(36,*)
    do l=1,nparts
          write(36,*) 'A', positions(l,1), positions(l,2), positions(l,3)
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! siguiente punto
        time = time + dt
        call velocity_verlet_andersen(nparts,box_l,cutoff,nu,sigma,dt,positions,velocities,pot,forces)
        dr2 = dr2 + sum((positions - oldpos)**2)
            if (mod(i,100) == 0) then
                write(13,*) time*timeconv, dr2*lconv**2/dble(nparts)/100.d0
                oldpos = positions
            end if
        if (mod(i,k_steps) == 0) then
            ! energia cinetica, temperatura y energia total
            call Kinetic_Energy(nparts, velocities, kin)
            tot = kin + pot
            tempi = kin*2.d0/(3.d0*dnparts-3.d0)
            call pressure(nparts,rho,tempi,box_l,cutoff,positions,pres)
            write(12,*) time*timeconv, kin*econv, pot*econv, tot*econv, tempi*tempconv, pres*pconv
        end if
    end do
    positions = positions*lconv
    write(15,*) nparts, rho*rhoconv
    do i = 1, nparts
        write(15,*) positions(i,1), positions(i,2), positions(i,3)
    end do

    close(15)
    close(12)
    close(13)
    close(36)!cerrar vmd.xyz
end program Main

!**********************************************************************************************************************************!
! read input file
subroutine Read_parameters(params_file, nparts, geometry, density, mass, sigma, epsilon, cutoff, temp_room, temp_init, bimodal, &
    disorder_system, thermostat, integrator, time_step, steps, measure_steps, observables_file, MSD_file, positions_file)

    implicit none
    ! input
    character(len=50), intent(in) :: params_file
    ! output
    integer, intent(out) :: nparts
    character(len=3), intent(out) :: geometry ! SC/FCC
    double precision, intent(out) :: density ! g/cm^3
    double precision, intent(out) :: mass ! g/mol
    double precision, intent(out) :: sigma ! A
    double precision, intent(out) :: epsilon ! KJ/mol
    double precision, intent(out) :: cutoff ! relative to cell length
    double precision, intent(out) :: temp_room ! K
    double precision, intent(out) :: temp_init ! K
    character(len=2), intent(out) :: bimodal ! Si/No
    character(len=2), intent(out) :: disorder_system ! Si/No
    character(len=2), intent(out) :: thermostat ! Si/No
    character(len=6), intent(out) :: integrator ! Euler/Verlet
    double precision, intent(out) :: time_step ! ps
    integer, intent(out) :: steps ! steps of simulation
    integer, intent(out) :: measure_steps ! steps between two measures
    character(len=80), intent(out) :: observables_file ! energy, temperature, etc evolution
    character(len=80), intent(out) :: MSD_file ! mean square displacement
    character(len=80), intent(out) :: positions_file ! histogram of particles positions
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
    read(my_unit,*) geometry
    read(my_unit,*)
    read(my_unit,*) density
    read(my_unit,*) 
    read(my_unit,*) mass
    read(my_unit,*) 
    read(my_unit,*) sigma
    read(my_unit,*) 
    read(my_unit,*) epsilon
    read(my_unit,*) 
    read(my_unit,*) cutoff
    read(my_unit,*)
    read(my_unit,*) temp_room
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*)
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
    read(my_unit,*) steps
    read(my_unit,*) 
    read(my_unit,*) measure_steps
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*)
    read(my_unit,*) 
    read(my_unit,*) observables_file
    read(my_unit,*) 
    read(my_unit,*) MSD_file
    read(my_unit,*) 
    read(my_unit,*) positions_file
    close(my_unit)

    return

end subroutine Read_parameters
!**********************************************************************************************************************************!
