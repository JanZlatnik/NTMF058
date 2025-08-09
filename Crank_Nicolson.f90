! ============================================================================================ !
! CRANK NICOLSON MODULE
! ============================================================================================ !
! CONTAINS:     solver of Schrödinger equation using CN-method
! LAST EDIT:    08/08/2025 by JZ
! ============================================================================================ !
    
module crank_nicolson_mod
    use parameters_mod
    use io_mod
    implicit none
    
    contains
    
    
    ! == Build CN band matrices ==
    subroutine build_cn_matrices(A, B, potential)
        complex*16, dimension(:,:), intent(out) :: A, B
        complex*16, dimension(:,:), intent(in)  :: potential
        
        integer     :: ix, iy, k
        complex*16  :: alphax, alphay, beta, pot
        
        call console('Building band matrices...')
        
        alphax = iu * dt * hbar / (4.0 * mass * dx**2)
        alphay = iu * dt * hbar / (4.0 * mass * dy**2)
        beta = iu * dt / (2.0 * hbar)
        
        A = (0.0, 0.0)
        B = (0.0, 0.0)
        
        do iy = 1, Ny
            do ix = 1, Nx
                k = ix + (iy - 1) * Nx
                
                pot = beta * potential(ix, iy)
                
                A(ku + kl + 1, k) = 1.0 + 2.0 * alphax + 2.0 * alphay + pot
                B(ku + 1, k) = 1.0 - 2.0 * alphax - 2.0 * alphay - pot
                
                if (ix > 1) then
                    A(ku+kl, k) = -alphax
                    B(ku, k) = alphax
                end if
                
                if (ix < Nx) then
                    A(ku+kl+2, k) = -alphax
                    B(ku+2, k) = alphax
                end if
                
                if (iy > 1) then
                    A(ku+kl+1-Nx, k) = -alphay
                    B(ku+1-Nx, k) = alphay
                end if
                
                if (iy < Ny) then
                    A(ku+kl+1+Nx, k) = -alphay
                    B(ku+1+Nx, k) = alphay
                end if
                
            end do
        end do
        
    end subroutine build_cn_matrices
    
    
    ! == CN solver ==
    subroutine CN_solver(psi, A, B, wfn_u, flux_u)
        complex*16, dimension(:,:), intent(inout)   :: psi
        complex*16, dimension(:,:), intent(in)      :: B
        complex*16, dimension(:,:), intent(inout)   :: A
        integer, intent(in)                         :: wfn_u, flux_u
        character(len=50)                           :: message
    
        external zgbmv, zgbtrf, zgbtrs
        
        complex*16, dimension(N)    :: psi_vec, rhs_vec
        integer, dimension(N)       :: ipiv
        integer                     :: info, it
        
        call console('Starting Crank-Nicolson solver...')
        
        ! == LU decomposition ==
        call console('Starting LU decomposition...')
        call zgbtrf(N, N, kl, ku, A, 2*kl+ku+1, ipiv, info)
        if (info/=0) then
            write(message, '(A,I0)') '[ERROR]: LU decomposition failed with info: ', info
            call console(message)
            stop
        else
            call console('LU decomposition completed successfully.')
        end if
        
        call wavefunction_o(psi, wfn_u)
        call flux_o(psi, flux_u)      
        
        ! == Time loop ==
        call console('Starting evolutuion computation...')
        do it = 1, Nt
            psi_vec = reshape(psi, (/N/))
            
            call zgbmv('N', N, N, kl, ku, (1.0d0,0.0d0), B, kl+ku+1, psi_vec, 1, (0.0d0,0.0d0), rhs_vec, 1) 

            call zgbtrs('N', N, kl, ku, 1, A, 2*kl+ku+1, ipiv, rhs_vec, N, info)
            if (info/=0) then
                write(message, '(A,I0)') '[ERROR]: LU solver failed with info: ', info
                call console(message)
                stop
            end if
            
            psi_vec = rhs_vec
            psi = reshape(psi_vec, (/Nx, Ny/))
            
            call wavefunction_o(psi, wfn_u)
            call flux_o(psi, flux_u)
            
            write(message, '(A,I0,A,I0,A)') 'Successfully computed step ', it, ' out of ', Nt,'.'
            call console(message)
            
        end do
        
        
    end subroutine CN_solver
    
        

end module crank_nicolson_mod