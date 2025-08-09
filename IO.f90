! ============================================================================================ !
! IO MODULE
! ============================================================================================ !
! CONTAINS:     input and output of data
! LAST EDIT:    08/08/2025 by JZ
! ============================================================================================ !
    
module io_mod
    use parameters_mod
    implicit none
    
    contains
    
    
    ! == Console output ==
    subroutine console(message)
        character(len=*), intent(in) :: message
        character(len=10) :: time_str
        
        call date_and_time(time = time_str)
        print*, '[', time_str(1:2), ':', time_str(3:4), ':', time_str(5:6), ']: ', message
    end subroutine console
    
    
    ! == Open output files ==
    subroutine open_output_files(wavefunction_unit, flux_unit)
        integer, intent(out) :: wavefunction_unit, flux_unit
        
        open(newunit = wavefunction_unit, file = 'wavefunction.bin', form = 'unformatted', access = 'stream', status = 'replace')
        open(newunit = flux_unit, file = 'flux.bin', form = 'unformatted', access = 'stream', status = 'replace')
        
        call console('Output files have been succesfully opened.')
        
    end subroutine open_output_files
    
    
    ! == Wavefunction output ==
    subroutine wavefunction_o(psi, wfc_unit)
        complex*16, dimension(:,:), intent(in) :: psi
        integer, intent(in) :: wfc_unit
        write(wfc_unit) psi
    end subroutine wavefunction_o
    
    
    ! == Flux output ==
    subroutine flux_o(psi, flux_unit)
        complex*16, dimension(:,:), intent(in) :: psi
        integer, intent(in) :: flux_unit
        
        real*8, dimension(Ny)   :: jx_vec
        complex*16              :: dpsi_dx
        integer                 :: iy
        
        do iy = 1, Ny
            dpsi_dx = ( psi(Nx_flux+1,iy) - psi(Nx_flux-1,iy) ) / (2.0 * dx)
            jx_vec(iy) = (hbar / mass) * aimag( conjg(psi(Nx_flux,iy)) * dpsi_dx ) 
        end do
        
        write(flux_unit) jx_vec
        
    end subroutine flux_o
    
    
    ! == Close output files ==
    subroutine close_output_files(wavefunction_unit, flux_unit)
        integer, intent(in) :: wavefunction_unit, flux_unit
        
        close(wavefunction_unit)
        close(flux_unit)
        
        call console('Output files have been succesfully closed.')
    end subroutine close_output_files
    
    
    ! == Parameters export ==
    subroutine parameters_write(unit)
        integer, intent(in) :: unit
        
        write(unit,*) "=============================================================="
        write(unit,*) "                 Simulation Parameters                        "
        write(unit,*) "=============================================================="
        
        write(unit,*) ""
        write(unit,*) "[Grid]"
        write(unit,'(A25,1X,I12)') "Nx:", Nx
        write(unit,'(A25,1X,I12)') "Ny:", Ny
        write(unit,'(A25,1X,I12)') "Nt:", Nt
        write(unit,'(A25,1X,F12.6)') "xmax:", xmax
        write(unit,'(A25,1X,F12.6)') "ymax:", ymax
        write(unit,'(A25,1X,F12.6)') "dx:", dx
        write(unit,'(A25,1X,F12.6)') "dy:", dy
        write(unit,'(A25,1X,F12.6)') "dt:", dt

        ! Band
        write(unit,*) ""
        write(unit,*) "[Band matrix]"
        write(unit,'(A25,1X,I12)') "N:", N
        write(unit,'(A25,1X,I12)') "kl:", kl
        write(unit,'(A25,1X,I12)') "ku:", ku

        ! Flux
        write(unit,*) ""
        write(unit,*) "[Flux]"
        write(unit,'(A25,1X,F12.6)') "x_flux:", x_flux
        write(unit,'(A25,1X,I12)') "Nx_flux:", Nx_flux

        ! Wave packet
        write(unit,*) ""
        write(unit,*) "[Wave packet]"
        write(unit,'(A25,1X,F12.6)') "mass:", mass
        write(unit,'(A25,1X,F12.6)') "x0:", x0
        write(unit,'(A25,1X,F12.6)') "y0:", y0
        write(unit,'(A25,1X,F12.6)') "sigmax:", sigmax
        write(unit,'(A25,1X,F12.6)') "sigmay:", sigmay
        write(unit,'(A25,1X,F12.6)') "px0:", px0
        write(unit,'(A25,1X,F12.6)') "py0:", py0

        ! Barrier
        write(unit,*) ""
        write(unit,*) "[Barrier]"
        write(unit,'(A25,1X,F12.6)') "V_A:", V_A
        write(unit,'(A25,1X,F12.6)') "V_B:", V_B
        write(unit,'(A25,1X,F12.6)') "V_C:", V_C
        write(unit,'(A25,1X,F12.6)') "V_d:", V_d

        ! CAP
        write(unit,*) ""
        write(unit,*) "[CAP]"
        write(unit,'(A25,1X,F12.6)') "cap_alpha:", cap_alpha
        write(unit,'(A25,1X,F12.6)') "cap_x:", cap_x
        write(unit,'(A25,1X,F12.6)') "cap_y:", cap_y

        write(unit,*) "=============================================================="
        
    end subroutine parameters_write
    
    subroutine parameters_o()
        integer :: param_u
        integer :: console_u = 6
        
        call parameters_write(console_u)
        
        open(newunit = param_u, file = 'parameters.txt', status = 'replace', form = 'formatted')
        call parameters_write(param_u) 
        close(param_u)
        call console('Parameters have been succesfully written to parameters.txt.')
        
    end subroutine parameters_o
    
    
    ! == Potential output ==
    subroutine potential_o(potential)
    implicit none
    complex*16, intent(in) :: potential(Nx,Ny)
    integer :: ix, iy, pot_unit

    open(newunit=pot_unit, file='potential.csv', status='replace', action='write')

    write(pot_unit,*) "x,y,Re(V),Im(V)"

    do iy = 1, Ny
        do ix = 1, Nx
            write(pot_unit,'(4(1X,F15.8))')  &
                -xmax + (ix-1)*dx,   &
                -ymax + (iy-1)*dy,   &
                real(potential(ix,iy)),  &
                aimag(potential(ix,iy))
        end do
    end do

    close(pot_unit)
    call console('Potential have been succesfully written to potential.csv.')
end subroutine potential_o

    
    

    
end module io_mod