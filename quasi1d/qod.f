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
            real(dp), public, parameter :: pi = 4*atan(1.0_8)
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
            integer, public, parameter :: nx = 40
            real(dp), public, parameter :: dx = 1.0_dp/nx
            real(dp), private, parameter :: exit_p_ratio = 0.8_dp
            real(dp), public, parameter :: exit_p = exit_p_ratio*ptot_in
            real(dp), public, parameter :: eps = 0.8_dp
            integer, parameter :: flx_scheme = 1
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
            implicit none
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
            use input
            use constants
            implicit none
            real(dp), dimension(:, :), intent(out) :: prim
            real(dp) :: rh0, p0, t0, u0
            integer :: n, i
            n = size(prim, 2)
            t0 = ttot_in/(1 + 0.5*(g - 1.0)*M_in**2
            p0 = ptot_in/(1 + 0.5*(g - 1.0)*M_in**2)**(g/(g-1))
            rho0 = p0/(R*t0)
            u0 = 0
            rho_e = p/(g - 1)
            do i = 1, n
                prim(1, i) = rho0
                prim(2, i) = u0
                prim(3, i) = p0
                prim(4, i) = rho_e/rho
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
        subroutine calc_prim(w, prim)
            use constants, only: g
            real(dp), dimension(:, :), intent(in) :: w
            real(dp), dimension(5, size(w, 2)), intent(out) :: prim
            n = size(w, 2)
            do i = 1, n
                rho = w(1, i)
                u = w(2, i)/rho
                e = w(3, i)
                E = e/rho
                p = (g - 1)*rho*(E - 0.5*u**2)
                prim(1, i) = rho
                prim(2, i) = u
                prim(3, i) = p
                prim(4, i) = e
                prim(5, i) = sqrt(g*p/rho)
            end do
C       Calculate W from prim
        subroutine calc_w(prim, w)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(out) :: w
            do i=1, n
                w(1, i) = prim(1, i)
                w(2, i) = prim(1, i)*prim(2, i)
                w(3, i) = prim(4, i)
            end do
        end subroutine
C       Calculate F
        subroutine calc_f(prim, f)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(out) :: f
C           f(1) = rho*u
C           f(2) = rho*u**2 + p
C           f(3) = (e + p)*u
            n = size(prim, 2)
            do i=1, n
                f(1) = prim(1, i)*prim(2, i)
                f(2) = prim(1)*prim(2)**2 + prim(3)
                f(3) = (prim(4) + prim(3))*prim(2)
            end do
        end subroutine
C       Calculate Q
        subroutine calc_q(prim, s, q)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(:), intent(in), :: s
            real(dp), dimension(3, size(s)), intent(out) :: q
            n = size(s)
            do i=2, n
                q(1, i) = 0
                q(2, i) = 0
                q(3, i) = prim(3, i)(s(i) - s(i-1))
            end do
        end subroutine
C       Calculate maximum eigenvalue
        subroutine calc_lambda(u, c, lambda)
            real(dp), intent(in) :: u, c
            real(dp), intent(out) :: lambda
            lambda = max(u, u + c, u - c)
        end subroutine
C       Calculate residuals
        subroutine calc_r(f_edge, s_edge, q, r)
            real(dp), dimension(:, :), intent(in) :: f_edge, s_edge
            real(dp), dimension(:, :), intent(in) :: q
            real(dp), dimension(size(q,1), size(q,2)), intent(out) :: r
            n = size(q, 2)
            do i = 2, n - 1
            do k = 1, 3
                r(k, i) = f_edge(k, i)*s_edge(k, i)
     &                    - f_edge(k, i - 1)*s_edge(k, i -1)
     &                    + q[k, i]
            end do
            end do
        end subroutine
C       Calculate error over all residuals
        function calc_err(r) result(err)
            real(dp), dimension(:, :), intent(in) :: r
            real(dp) :: err
            err = maxval(r)
        end function
        end module
C       ===============================================================
C       Schemes for flux evaluation.
C       ===============================================================
        module flx_schemes
        use types, only: dp
        use input, only: flx_scheme, eps
        implicit integer (i, n)
        contains 

C       Scalar Dissipation
        subroutine flx_scalar(prim, w, f, f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3,size(prim,2)-1), intent(out) :: f_edge
            real(dp) :: u_avg, c_avg, lambda
            n = size(prim, 2)
            do i = 1, n - 1
                u_avg = 0.5*(prim(2, i) + prim(2, i+1))
                c_avg = 0.5*(c(i) + c(i+1))
                call calc_lambda(u_avg, c_avg, lambda)
                do k = 1, 3
                    f_edge(k, i) = 0.5*(
     &                  f(k, i) + f(k, i+1)
     &                  - eps*lambda*(w(k, i+1) - w(i)))
                end do
            end do
        end subroutine
C       Choose a flux evaluation scheme based on input
C       1 : Scalar dissipation
        subroutine flx_eval(prim, w, f, f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3,size(prim,2)-1), intent(out) :: f_edge
            select case (flx_scheme)
                case (1)
                    call flx_scalar(prim, f_edge)
                case default
                    call flx_scalar(prim, f_edge)
            end select
        end subroutine
        end module
C       ===============================================================
C       Time stepping
C       ===============================================================
C       Calculate residuals
C       First order euler
        subroutine euler_xp(prim, s, s_edge, w_n, err)
            use input
            implicit none
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: s 
            real(dp), dimension(size(prim,2)-1), intent(out) :: s_edge
            real(dp), dimension(3, size(prim, 2)), intent(out) :: w_n
            real(dp), intent(out) :: err
            real(dp), dimension(3, size(prim, 2)) :: w, f, q, r
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp) :: dt_v
C           TODO calculate DT
            call calc_w(prim, w)
            call calc_f(prim, f)
            call calc_q(prim, s, q)
            call flx_eval(prim, w, f, c, f_edge)
            call calc_r(f_edge, s_edge, q, r)
            n = size(prim, 2)
            do i = 2, n-1
                dt_v = dt(i)/(s(i)*dx)
                do k = 1, 3
                    w_n(k, i) = w(k, i) - dt_v*r(k, i)
                end do
            end do
            call calc_err(r, err)
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
        implicit integer (i, n)
        use constants
        use input, only: nx, tol
        real(dp) :: er
        real(dp), dimension(nx) :: x, s
        real(dp), dimension(5, nx) :: prim
        call mkgrid(x, s)
        call init(prim)
        err = 1
        do while (err > tol)
            call time_step(prim, s, s_edge, w_n, err)
            call calc_prim(w_n, prim)
            call update_bc(prim)
            err = calc_err
        end do
C       TODO post process
        end program
