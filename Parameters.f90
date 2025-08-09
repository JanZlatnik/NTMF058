! ============================================================================================ !
! PARAMETERS MODULE
! ============================================================================================ !
! CONTAINS:     physical and numerical parameters
! LAST EDIT:    07/08/2025 by JZ
! ============================================================================================ !
    
module parameters_mod
    implicit none
    
    
    ! == Physical constants ==
    real*8, parameter       :: hbar = 1.0
    complex*16, parameter   :: iu = (0.0,1.0)
    
    ! == Grid parameters ==
    integer, parameter  :: Nx = 501                                     
    integer, parameter  :: Nt = 2000
    real*8, parameter   :: xmax = 20.0
    real*8, parameter   :: ymax = 15.0
    real*8, parameter   :: dx = 2.0 * xmax / (Nx - 1)
    real*8, parameter   :: dy = dx
    integer, parameter  :: Ny = NINT(2.0 * ymax / dy) + 1
    real*8, parameter   :: dt = dx**2 / 4.0
    
    
    ! == Band matrices parameters ==
    integer, parameter  :: N = Nx * Ny
    integer, parameter  :: kl = Nx
    integer, parameter  :: ku = Nx
    
    
    ! == Flux parameters ==
    real*8, parameter   :: x_flux = 0.7 * xmax
    integer, parameter  :: Nx_flux = NINT(x_flux + xmax) / dx + 1
    
    
    ! == Wave-packet parameters ==
    real*8, parameter   :: mass = 1.0
    real*8, parameter   :: x0 = -0.5 * xmax
    real*8, parameter   :: y0 = 0.0
    real*8, parameter   :: sigmax = 1.5
    real*8, parameter   :: sigmay = 1.5
    real*8, parameter   :: px0 = 15.0
    real*8, parameter   :: py0 = 0.0 
    
    
    ! == Barrier parameters ==
    real*8, parameter   :: V_A = (px0**2+py0**2)/mass       ! height - twice the kinetic energy
    real*8, parameter   :: V_B = 0.75                       ! x-width of the barrier
    real*8, parameter   :: V_C = 2.0                        ! y-width of the slits
    real*8, parameter   :: V_d = 3.0                        ! half-y-distance between the slits
    
    
    ! == CAP prarameters ==
    real*8, parameter   :: cap_alpha = V_A/20.0                  ! CAP strenght coefficient
    real*8, parameter   :: cap_x = 0.9 * xmax
    real*8, parameter   :: cap_y = ymax - (xmax-cap_x)
    
    

end module parameters_mod