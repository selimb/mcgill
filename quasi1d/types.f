C       ===============================================================
C       Types
C       ===============================================================
        module types
            implicit none
C           Use double precision
            integer, parameter :: dp=kind(0.d0) 
        end module
        module input
            use types, only: dp
            use constants, only: ptot_in
            implicit none
C           User input
            integer, public, parameter :: flx_scheme = 1
            integer, public, parameter :: timestep_scheme = 1
            integer, public, parameter :: nx = 40
            real(dp), private, parameter :: p_exit_ratio = 0.8_dp
            real(dp), public, parameter :: eps = 0.8_dp
            real(dp), public, parameter :: tol = 1e-15
            real(dp), public, parameter :: cfl = 0.8
C           Calculations
            real(dp), public, parameter :: dx = 1.0_dp/nx
            real(dp), public, parameter :: p_exit = p_exit_ratio*ptot_in
        end module
