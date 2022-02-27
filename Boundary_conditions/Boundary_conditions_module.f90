module pbc_module

    !parameters
    implicit none

    !functions and subroutines
    contains

    subroutine pbc(position,length)

        !------------------------------------------------------------------------------------------------------------------------------!
        ! Informació
        ! La subrutina aplica les condicions periòdiques de contorn (pbc) sobre una partícula.
        ! Variables d'entrada:
        !   length: longitud de cada costat de la caixa
        ! Variables d'entrada i sortida (modificades):
        !   position(3): vector que conte la posició on es troba la partícula
        !------------------------------------------------------------------------------------------------------------------------------!
    
        implicit none
        double precision, intent(in) :: length
        double precision :: position(3)
        integer :: i
       
        do i = 1, 3
            if(position(i) >  0.5d0*length) position(i) = position(i) - length
            if(position(i) < -0.5d0*length) position(i) = position(i) + length
        end do
    
        return
       
    end subroutine pbc
    
    subroutine pbc_nparts(nparts,positions,length)
    
        !------------------------------------------------------------------------------------------------------------------------------!
        ! Informació
        ! La subrutina aplica les condicions periòdiques de contorn (pbc) sobre totes les partícules del sistema.
        ! Variables d'entrada:
        !   nparts: número de partícules del sistema
        !   length: longitud de cada costat de la caixa
        ! Variables d'entrada i sortida (modificades):
        !   position(nparts,3): matriu que conte les posicions de les nparts partícules del sistema
        !------------------------------------------------------------------------------------------------------------------------------!
    
        implicit none
        integer, intent(in) :: nparts
        double precision, intent(in) :: length
        double precision :: positions(nparts,3)
        integer :: i, n
    
        do n = 1, nparts
            do i = 1, 3
                if(positions(n,i) >  0.5d0*length) positions(n,i) = positions(n,i) - length
                if(positions(n,i) < -0.5d0*length) positions(n,i) = positions(n,i) + length
            end do
        end do
    
        return
    
    end subroutine pbc_nparts
end module pbc_module
