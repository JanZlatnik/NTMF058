! ============================================================================================ !
! INITIALIZATION MODULE
! ============================================================================================ !
! CONTAINS:     grid, potential and wave-packet initialization
! LAST EDIT:    08/08/2025 by JZ
! ============================================================================================ !
    
module initialization_mod
    use parameters_mod
    use io_mod


    contains
    
    ! == Grid initialization ==
    subroutine grid_init(x_coords,y_coords)
        real*8, dimension(:), intent(out) :: x_coords, y_coords
        integer :: i 
        do i = 1, Nx
            x_coords(i) = -xmax + (i-1) * dx
        end do
        do i = 1, Ny
            y_coords(i) = -ymax + (i-1) * dy
        end do
        call console('Grid has been succesfully initialized.')
    end subroutine grid_init
    
    
    ! == Potential initialization ==
    subroutine potential_init(potential, x_coords, y_coords)
        complex*16, dimension(:,:), intent(out) :: potential
        real*8, dimension(:), intent(in)        :: x_coords, y_coords
        
        integer :: ix, iy
        real*8  :: x, y, v_barrier, v_capx, v_capy
        
        do iy = 1, Ny
            y = y_coords(iy)
            
            v_capy = 0.0
            if (abs(y) > cap_y) then
                v_capy = cap_alpha * (abs(y) - cap_y)**3
            end if
            
            do ix = 1, Nx
                x = x_coords(ix)
                
                v_barrier = V_A * exp(-(x/V_B)**2) * (1 - exp(-((y-V_d)/V_C)**2 ) ) * (1 - exp(-((y+V_d)/V_C)**2 ) )
                
                v_capx = 0.0
                if (abs(x) > cap_x) then
                    v_capx = cap_alpha * (abs(x) - cap_x)**3
                end if
                
                potential(ix,iy) = v_barrier - iu * (v_capx + v_capy)
            end do
        end do
        call console('Potential has been succesfully initialized.')
    end subroutine potential_init
    
    
    ! == Wavepacket initialization ==
    subroutine wavepacket_init(psi, x_coords, y_coords)
        complex*16, dimension(:,:), intent(out) :: psi
        real*8, dimension(:), intent(in)        :: x_coords, y_coords
        
        integer :: ix, iy
        real*8  :: x, y, norm
        
        do iy = 1, Ny
            y = y_coords(iy)
            do ix = 1, Nx
                x = x_coords(ix)
                psi(ix,iy) = exp(-((x-x0)/sigmax)**2 - ((y-y0)/sigmay)**2 ) * exp(iu * (px0 * x + py0 * y) / hbar)
            end do
        end do
        
        norm = sum( abs(psi)**2 ) * dx * dy
        psi = psi / sqrt(norm)

        call console('Wave-packet has been succesfully initialized.')
    end subroutine wavepacket_init


end module initialization_mod