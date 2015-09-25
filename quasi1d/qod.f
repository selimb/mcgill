C       ===============================================================
C       Note on implementation.
C           We pass around `prim` between functions. This corresponds
C           to (rho, u, p e, c), the primitive variables. Generating W and
C           F is easy, and the primitive variables are readily 
C           available when necessary, such as when plotting or
C           calculating the maximum eigenvalues. 
C       Picture of grid indices:
C             node        edge
C           |---o---|---o---|---o ... ---|---o---|
C               1   1   2   2   3       J-1  J
C       ===============================================================
C       ===============================================================
C       Global variables: types, constants and inputs
C       ===============================================================
        module types
            implicit none
C           Use double precision
            integer, parameter :: dp=kind(0.d0) 
        end module
        module constants
            use types, only: dp
            implicit none
C           Math
            real(dp), public, parameter :: pi = 4.0_dp*atan(1.0_dp)
C           Grid constants
            real(dp), public, parameter :: h = 0.025_dp
            real(dp), public, parameter :: t1 = 0.8_dp
            real(dp), public, parameter :: t2 = 3
C           Fluid constants
            real(dp), public, parameter :: g = 1.4_dp
            real(dp), public, parameter :: R = 1716
C           BCs
            real(dp), public, parameter :: ttot_in = 530.2_dp
            real(dp), public, parameter :: ptot_in = 2117
            real(dp), public, parameter :: M_in = 1.2_dp
        end module
        module input
            use types, only: dp
            use constants, only: ptot_in
            implicit none
C           User input
            integer, parameter :: flx_scheme = 1
            integer, public, parameter :: nx = 40
            real(dp), private, parameter :: exit_p_ratio = 0.8_dp
            real(dp), public, parameter :: eps = 0.8_dp
            real(dp), public, parameter :: tol = 1e-15
C           Calculations
            real(dp), public, parameter :: dx = 1.0_dp/nx
            real(dp), public, parameter :: exit_p = exit_p_ratio*ptot_in
        end module
C       ===============================================================
C       Initialization
C       ===============================================================
        module init
        use types, only: dp
        implicit integer (i, n)
        public mkgrid, init_state
        contains 

C       Make initial grid.
        subroutine mkgrid(x, s)
            use input, only: dx
            use constants, only: h, pi, t1, t2
            real(dp), dimension(:), intent(out) :: x, s
            integer :: n, i
            n = size(x)
            do i=1, n
                x(i) = i*dx
                s(i) = 1.0_dp - h*(sin(pi*x(i)**t1))**t2
            end do
        end
C       Initialize field        
        subroutine init_state(prim)
            use constants, only: g, m_in, ptot_in, R
            real(dp), dimension(:, :), intent(out) :: prim
            real(dp) :: rh0, p0, t0, u0, rho_e, M_term
            n = size(prim, 2)
            M_term = ( 1.0_dp + 0.5_dp*(g - 1.0_dp)*M_in**2 )
            t0 = ttot_in/M_term
            p0 = ptot_in/( M_term**(g/(g - 1.0_dp)) )
            rho0 = p0/(R*t0)
            u0 = 0
            rho_e = p0/(g - 1.0_dp)
            do i = 1, n
                prim(1, i) = rho0
                prim(2, i) = u0
                prim(3, i) = p0
                prim(4, i) = rho_e/rho0
                prim(5, i) = sqrt(g*p0/rho0)
            end do
        end
        end module init
C       ===============================================================
C       Helper functions to calculate commonly required quantities
C       ===============================================================
        module common_calcs
        use types, only: dp
        implicit integer (i, n)
        contains 

C       Calculate prim from W
        pure function calc_prim(w) result(prim)
            use constants, only: g
            real(dp), dimension(:, :), intent(in) :: w
            real(dp), dimension(5, size(w, 2)) :: prim
            n = size(w, 2)
            do i = 1, n
                rho = w(1, i)
                u = w(2, i)/rho
                e = w(3, i)
                E = e/rho
                p = (g - 1.0_dp)*rho*(E - 0.5_dp*u**2)
                prim(1, i) = rho
                prim(2, i) = u
                prim(3, i) = p
                prim(4, i) = e
                prim(5, i) = sqrt(g*p/rho)
            end do
        end function
C       Calculate W from prim
        pure function calc_w(prim) result(w)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)) :: w
            do i=1, n
                w(1, i) = prim(1, i)
                w(2, i) = prim(1, i)*prim(2, i)
                w(3, i) = prim(4, i)
            end do
        end function
C       Calculate F
        pure function calc_f(prim) result(f)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)) :: f
C           f(1) = rho*u
C           f(2) = rho*u**2 + p
C           f(3) = (e + p)*u
            n = size(prim, 2)
            do i=1, n
                f(1, i) = prim(1, i)*prim(2, i)
                f(2, i) = prim(1, i)*prim(2, i)**2 + prim(3, i)
                f(3, i) = (prim(4, i) + prim(3, i))*prim(2, i)
            end do
        end function
C       Calculate Q
        pure function calc_q(prim, s) result(q)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: s
            real(dp), dimension(3, size(prim, 2)) :: q
            n = size(prim, 2)
            do i=2, n
                q(1, i) = 0
                q(2, i) = prim(3, i)*(s(i) - s(i-1))
                q(3, i) = 0
            end do
        end function
C       Calculate maximum eigenvalue
        pure function calc_lambda(u, c) result(lambda)
            real(dp), intent(in) :: u, c
            real(dp) :: lambda
            lambda = max(u, u + c, u - c)
        end function
C       Calculate residuals
        pure function calc_r(f_edge, s_edge, q) result(r)
            real(dp), dimension(:, :), intent(in) :: f_edge
            real(dp), dimension(size(f_edge, 2)), intent(in) :: s_edge
            real(dp), dimension(3, size(f_edge, 2) + 1), intent(in) :: q
            real(dp), dimension(3, size(f_edge, 2) + 1) :: r
            n = size(q, 2)
            do i = 2, n - 1
            do k = 1, 3
                r(k, i) = f_edge(k, i)*s_edge(i)
     &                    - f_edge(k, i - 1)*s_edge(i -1)
     &                    + q(k, i)
            end do
            end do
        end function
C       Calculate error over all residuals
        pure function calc_err(r) result(err)
            real(dp), dimension(:, :), intent(in) :: r
            real(dp) :: err
            err = maxval(r)
        end function
        end module common_calcs
C       ===============================================================
C       Schemes for flux evaluation.
C       ===============================================================
        module flx_schemes
        use types, only: dp
        use input, only: flx_scheme, eps
        use common_calcs
        implicit integer (i, n)
        contains 

C       Scalar Dissipation
        pure function flx_scalar(prim, w, f) result (f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp) :: u_avg, c_avg, lambda
            n = size(prim, 2)
            do i = 1, n - 1
                u_avg = 0.5_dp*(prim(2, i) + prim(2, i+1))
                c_avg = 0.5_dp*(prim(5, i) + prim(5, i+1))
                lambda = calc_lambda(u_avg, c_avg)
                do k = 1, 3
                    f_edge(k, i) = 0.5_dp*(
     &                  f(k, i) + f(k, i+1)
     &                  - eps*lambda*(w(k, i+1) - w(k, i)))
                end do
            end do
        end function
C       Choose a flux evaluation scheme based on input
C       1 : Scalar dissipation
        pure function flx_eval(prim, w, f) result (f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3,size(prim,2)-1) :: f_edge
            select case (flx_scheme)
                case (1)
                    f_edge = flx_scalar(prim, w, f)
                case default
                    f_edge = flx_scalar(prim, w, f)
            end select
        end function
        end module flx_schemes
C       ===============================================================
C       Time stepping
C       ===============================================================
C       Calculate residuals
C       First order euler
        module timestepping
        use types, only: dp
        use input, only: dx
        use common_calcs
        use flx_schemes, only: flx_eval
        implicit integer (i, n)
        contains 

        subroutine euler_xp(prim, s, s_edge, w_n, err)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: s 
            real(dp), dimension(size(prim,2)-1), intent(out) :: s_edge
            real(dp), dimension(3, size(prim, 2)), intent(out) :: w_n
            real(dp), intent(out) :: err
            real(dp), dimension(3, size(prim, 2)) :: w, f, q, r
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp), dimension(size(prim, 2)) :: dt
            real(dp) :: dt_v
C           TODO calculate DT

            w = calc_w(prim)
            f = calc_f(prim)
            q = calc_q(prim, s)
            f_edge = flx_eval(prim, w, f)
            r = calc_r(f_edge, s_edge, q)
            n = size(prim, 2)
            do i = 2, n-1
                dt_v = dt(i)/(s(i)*dx)
                do k = 1, 3
                    w_n(k, i) = w(k, i) - dt_v*r(k, i)
                end do
            end do
            err = calc_err(r)
        end subroutine
        end module
C       ===============================================================
C       Boundary Conditions
C       ===============================================================
        module bc
        use types, only: dp
        implicit integer (i, n)
        contains 

        subroutine update_inlet(prim)
            real(dp), dimension(:, :), intent(out) :: prim
            real(dp) :: m
C       TODO
        end subroutine
        subroutine update_outlet(prim)
            real(dp), dimension(:, :), intent(out) :: prim
            real(dp) :: m
C       TODO
        end subroutine
        subroutine update_bc(prim)
            real(dp), dimension(:, :), intent(out) :: prim
            call update_inlet(prim)
            call update_outlet(prim)
        end subroutine
        end module
C       ===============================================================
C       Main
C       ===============================================================
        program main
        use types, only: dp
        use constants
        use input, only: nx, tol
        use init, only: init_state, mkgrid
        implicit integer (i, n)
        real(dp) :: er
        real(dp), dimension(nx) :: x, s
        real(dp), dimension(5, nx) :: prim
        call mkgrid(x, s)
        call init_state(prim)
C       err = 1
C       do while (err > tol)
C           call time_step(prim, s, s_edge, w_n, err)
C           call calc_prim(w_n, prim)
C           call update_bc(prim)
C           err = calc_err
C       end do
C       TODO post process
        end program
