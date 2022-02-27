program radial_distribution_function

    implicit none
    character(len=80) :: infile, outfile
    integer :: i, j, part, maxi, nparts, indice
    double precision :: rmax, dr, rho, integral
    double precision, allocatable :: positions(:,:), distances(:), gofr(:), vol(:)
    double precision, parameter :: pi = 4.d0*datan(1.d0)

    print*, 'Introduce el fichero de entrada.'
    read(*,*) infile
    open(12,file=infile)
    print*, 'Introduce el fichero de salida.'
    read(*,*) outfile
    print*,'Introduce el numero de cajas.'
    read(*,*) maxi
    open(12,file=infile)
    read(12,*) nparts, rho
    allocate(positions(nparts,3), distances(nparts*(nparts-1)/2), vol(maxi), gofr(maxi))
    do i = 1, nparts
        read(12,*) positions(i,1), positions(i,2), positions(i,3)
    end do
    close(12)
    part = 0
    do i = 1,nparts-1
        do j = i+1, nparts
            part = part + 1
            distances(part) = dsqrt(sum((positions(i,:)-positions(j,:))**2))
        end do
    end do
    rmax = maxval(distances)*1.03d0
    dr = rmax/dble(maxi)
    do i = 1, maxi
        vol(i) = 4.d0*pi/3.d0*dr**3*dble(i**3-(i-1)**3)
    end do
    gofr = 0.d0
    do i = 1, nparts*(nparts-1)/2
        indice = int(distances(i)/dr) + 1
        gofr(indice) = gofr(indice) + 1.d0
    end do
    gofr = 2.d0*gofr/(dble(nparts)*rho*vol)
    integral = 0.d0
    do i = 1, maxi
        integral = integral + gofr(i)*(i*dr-0.5d0*dr)**2*dr
    end do
    integral = integral/dble(nparts-1)*4.d0*pi*rho
    print*, integral
    open(15,file=outfile)
    do i = 1, maxi
        write(15,*) dr*i-0.5d0*dr, gofr(i), gofr(i)/integral
    end do
    close(15)

end program radial_distribution_function