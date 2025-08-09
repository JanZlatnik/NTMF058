! ============================================================================================ !
! 2D DOUBLE-SLIT EXPERIMENT SIMULATION PROGRAM
! ============================================================================================ !
! CONTAINS:     solver of 2D TDSE using Crank-Nicolson method on a rectangular grid with 
!               a potential barrier and CAP, computes pobability flux on a given line
! LAST EDIT:    08/08/2025 by JZ
! ============================================================================================ !
    
program double_slit
    use parameters_mod
    use initialization_mod
    use io_mod
    use crank_nicolson_mod
    
    implicit none
    
    real*8, dimension(Nx)               :: x_coords
    real*8, dimension(Ny)               :: y_coords
    
    complex*16, dimension(Nx, Ny)       :: potential
    complex*16, dimension(Nx, Ny)       :: psi
    
    complex*16, dimension(:,:), allocatable :: A
    complex*16, dimension(:,:), allocatable :: B
    
    integer :: wfn_u, flux_u, i, j
    character(len=20) :: val_str
    real*8 :: re, im
    
    
    print*, '=============================================================='
    print*, '         2D DOUBLE-SLIT EXPERIMENT SIMULATION PROGRAM         '
    print*, '=============================================================='
    
    call parameters_o()
    
    call grid_init(x_coords, y_coords)
    call potential_init(potential, x_coords, y_coords)
    call potential_o(potential)
    call wavepacket_init(psi, x_coords, y_coords)
    
    call open_output_files(wfn_u, flux_u)
    
    allocate(A(2*kl+ku+1, N), B(kl+ku+1, N))
    call build_cn_matrices(A, B, potential)
    
    
    
    !open(unit=11111, file='B_matrix.csv', status='replace', action='write')
    ! do i = 1, kl+ku+1
    ! do j = 1, N
    !    re = real(B(i,j))
    !    im = aimag(B(i,j))
    !    write(val_str, '(F0.2, A, F0.2, A)') re, ' + ', im, 'i'
    !    if (j < N) then
    !       write(11111, '(A)', advance='no') trim(val_str) // ','
    !    else
    !       write(11111, '(A)') trim(val_str)
    !    end if
    ! end do
    !end do


    call CN_solver(psi, A, B, wfn_u, flux_u)
    call close_output_files(wfn_u, flux_u)
    deallocate(A, B)
    
    print*, '=============================================================='
    print*, '               SIMULATION SUCCESFULLY FINISHED!               '
    print*, '=============================================================='
    print*, 'Press ENTER to exit...'
    read(*,*)
    

end program double_slit