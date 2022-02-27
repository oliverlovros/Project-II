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

    use constantes
    implicit none
    ! fitxers
    integer :: opt1, opt2
    character(len=80) :: infile, outfile, outfile1, outfile2
    ! variables d'entrada
    integer :: nparts
    double precision :: rho, temp
    integer :: t_steps, k_steps
    double precision :: dt
    double precision :: mass, LJ_eps, LJ_sig
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
    integer :: i

    ! llegim fitxers i iniciem variables
    print*, 'Fitxer de entrada:'
    read(*,*) infile
    open(12,file=infile)
    read(12,*) nparts, rho, temp
    read(12,*) dt, t_steps, k_steps
    read(12,*) mass, LJ_eps, LJ_sig
    read(12,*) outfile
    read(12,*) opt1, outfile1
    read(12,*) opt2, outfile2
    close(12)

    ! creem fitxers de sortida
    open(12,file=outfile)
    if (opt1 == 1) open(13,file=outfile1)
    if (opt2 == 1) open(15,file=outfile2)
    ! dimensions 
    allocate(positions(nparts,3), velocities(nparts,3), forces(nparts,3))
    ! iniciem la geometria
    call sc_lattice(nparts,rho,positions,box_l,box_a,box_m)
    if (opt1 == 1) then
        allocate(oldpos(nparts,3))
        oldpos = positions
    end if
    ! iniciem velocitats
    velocities = 0d0
    ! radi de cutoff
    cutoff = 0.5d0*box_l
    ! conversio unitats reduides a reals
    dnparts = dble(nparts)
    lconv = LJ_sig
    econv = LJ_eps/dnparts
    rhoconv = 1d24*mass/(navo*LJ_sig**3)
    tempconv = LJ_eps/(navo*kbol)
    pconv = 1d33*LJ_eps/(navo*LJ_sig**3d0)
    timeconv = 0.1d0*(mass*LJ_sig**2d0/LJ_eps)**0.5d0
    ! força, energia, pressió i temperatura inicial
    write(12,*) '# N = ', nparts, ', rho = ', rho*rhoconv, ', L = ', box_L*lconv, ', a = ', box_a*lconv, ', M = ', box_m
    write(12,*) '# time (ps), kin/N (kJ/mol), pot/N (kJ/mol), tot/N (kJ/mol), temp (K), pres (Pa)'
    time = 0.d0
    call LJ_potential(nparts,positions,cutoff,box_l,1, pot, forces)
    call Kinetic_Energy(nparts, velocities, kin)
    tot = kin + pot
    tempi = kin*2.d0/(3.d0*dnparts-3.d0)
    call pressure(nparts,rho,tempi,box_l,cutoff,positions,pres)
    write(12,*) time*timeconv, kin*econv, pot*econv, tot*econv, tempi*tempconv, pres*pconv
    ! comença el bucle temporal
    nu = 5.d0*dt
    sigma = dsqrt(temp)
    dr2 = 0.d0
    do i  = 1, t_steps
        ! següent punt
        time = time + dt
        call velocity_verlet_andersen(nparts,box_l,cutoff,nu,sigma,dt,positions,velocities,pot,forces)
        if (opt1 == 1) then
            dr2 = dr2 + sum((positions - oldpos)**2)
            if (mod(i,100) == 0) then
                write(13,*) time*timeconv, dr2*lconv**2/dble(nparts)/100.d0
                oldpos = positions
            end if
        end if
        if (mod(i,k_steps) == 0) then
            ! calculem energia cinètica, temperatura i energia total
            call Kinetic_Energy(nparts, velocities, kin)
            tot = kin + pot
            tempi = kin*2.d0/(3.d0*dnparts-3.d0)
            call pressure(nparts,rho,tempi,box_l,cutoff,positions,pres)
            write(12,*) time*timeconv, kin*econv, pot*econv, tot*econv, tempi*tempconv, pres*pconv
        end if
    end do
    if (opt2 == 1) then
        positions = positions*lconv
        write(15,*) nparts, rho*rhoconv
        do i = 1, nparts
            write(15,*) positions(i,1), positions(i,2), positions(i,3)
        end do
        close(15)
    end if

    close(12)
    if (opt1 == 1) close(13)

end program Main