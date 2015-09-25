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
C       Global variables: constants and input
C       ===============================================================
        module constants
            implicit none
C           Math
            real*8, public, parameter :: pi = 4*atan(1.0_8)
C           Grid constants
            real*8, public, parameter :: h = 0.025
            real*8, public, parameter :: t1 = 0.8
            real*8, public, parameter :: t2 = 3
C           Fluid constants
            real*8, public, parameter :: g = 1.4
            real*8, public, parameter :: R = 1716
C           BCs
            real*8, public, parameter :: ttot_in = 530.2
            real*8, public, parameter :: ptot_in = 2117.0
            real*8, public, parameter :: M_in = 1.2
        end module
        module input
            use constants
            implicit none
            integer, public, parameter :: nx = 40
            real*8, public, parameter :: dx = 1.0/nx
            real*8, private, parameter :: exit_p_ratio = 0.8
            real*8, public, parameter :: exit_p = 0.8*ptot_in
            real*8, public, parameter :: eps = 0.8
            integer, parameter :: flx_scheme = 1
        end module
C       ===============================================================
C       Initialization
C       ===============================================================
C       Make initial grid.
        subroutine mkgrid(x, s)
            use input
            use constants
            implicit none
            real*8, dimension(nx), intent(out) :: x, s
            integer :: i
            do i=1, nx
                x(i) = i*dx
                s(i) = 1 - h*(sin(pi*x(i)**t1))**t2
            end do
        end
C       Initialize field        
        subroutine init(prim)
            use input
            use constants
            implicit none
            real*8, dimension(5, nx), intent(out) :: prim
            real*8 :: rh0, p0, t0, u0
            t0 = ttot_in/(1 + 0.5*(g - 1.0)*M_in**2
            p0 = ptot_in/(1 + 0.5*(g - 1.0)*M_in**2)**(g/(g-1))
            rho0 = p0/(R*t0)
            u0 = 0.0
            rho_e = p/(g - 1)
            do i=1, nx
                prim(1, i) = rho0
                prim(2, i) = u0
                prim(3, i) = p0
                prim(4, i) = rho_e/rho
                prim(5, i) = sqrt(g*p0/rho0)
            end do
        end
C       ===============================================================
C       Helper functions to calculate commonly required quantities
C       ===============================================================
C       Calculate prim from W
        subroutine calc_prim(w, prim)
            use input
            use constants
            implicit none
            real*8, dimension(3, nx), intent(in) :: w
            real*8, dimension(5, nx), intent(out) :: prim
            do i=1, nx
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
C       Calculate W
        subroutine calc_w(prim, w)
            use input
            implicit none
            real*8, dimension(5, nx), intent(in) :: prim
            real*8, dimension(3, nx), intent(out) :: w
            do i=1, nx
                w(1, i) = prim(1, i)
                w(2, i) = prim(1, i)*prim(2, i)
                w(3, i) = prim(4, i)
            end do
        end subroutine
C       Calculate F
        subroutine calc_f(prim, f)
            use input
            implicit none
            real*8, dimension(5, nx), intent(in) :: prim
            real*8, dimension(3, nx), intent(out) :: f
C           f(1) = rho*u
C           f(2) = rho*u**2 + p
C           f(3) = (e + p)*u
            do i=1, nx
                f(1) = prim(1, i)*prim(2, i)
                f(2) = prim(1)*prim(2)**2 + prim(3)
                f(3) = (prim(4) + prim(3))*prim(2)
            end do
        end subroutine
C       Calculate Q
        subroutine calc_q(prim, s, q)
            use input
            implicit none
            real*8, dimension(5, nx), intent(in) :: prim
            real*8, intent(in, nx) :: s
            real*8, dimension(3, nx), intent(out) :: q
            do i=2, nx
                q(1, i) = 0
                q(2, i) = 0
                q(3, i) = prim(3, i)(s(i) - s(i-1))
            end do
        end subroutine
C       Calculate maximum eigenvalue
        subroutine calc_lambda(u, c, lambda)
            use input
            implicit none
            real*8, intent(in) :: u, c
            real*8, intent(out) :: lambda
            lambda = max(u, u + c, u - c)
        end subroutine
C       Calculate residuals
        subroutine calc_r(f_edge, s_edge, q, r)
            use input
            implicit none
            real*8, dimension(3, nx-1), intent(in) :: f_edge, s_edge
            real*8, dimension(3, nx), intent(in) :: q
            real*8, dimension(3, nx), intent(out) :: r
            do i = 2, nx - 1
            do k = 1, 3
                r(k, i) = f_edge(k, i)*s_edge(k, i)
     &                    - f_edge(k, i - 1)*s_edge(k, i -1)
     &                    + q[k, i]
            end do
            end do
        end subroutine
C       Calculate error over all residuals
        subroutine calc_err(r, err)
            use input
            implicit none
            real*8, dimension(3, nx), intent(in) :: r
            real*8, intent(out) :: err
            err = maxval(r)
        end subroutine
C       ===============================================================
C       Schemes for flux evaluation.
C       ===============================================================
C       Scalar Dissipation
        subroutine flx_scalar(prim, w, f, c, f_edge)
            use input
            implicit none
            real*8, dimension(5, nx), intent(in) :: prim
            real*8, dimension(3, nx), intent(in) :: w, f
            real*8, dimension(nx), intent(in) :: c
            real*8, dimension(3, nx-1), intent(out) :: f_edge
            real*8 :: u_avg, c_avg, lambda
            do i = 1, nx - 1
                u_avg = 0.5*(prim(2, i) + prim(2, i+1))
                c_avg = 0.5*(c(i) + c(i+1))
                calc_lambda(u_avg, c_avg, lambda)
                do k = 1, 3
                    f_edge(k, i] = 0.5*(
     &                  f(k, i) + f(k, i+1)
     &                  - eps*lambda*(w(k, i+1) - w(i)))
                end do
            end do
        end subroutine
C       Choose a flux evaluation scheme based on input
C       1 : Scalar dissipation
        subroutine flx_eval(prim, f_edge)
            use input
            implicit none
            real*8, dimension(5, nx), intent(in) :: prim
            real*8, dimension(3, nx-1), intent(out) :: f_edge
            select case (flx_scheme)
                case (1)
                    flx_scalar(prim, f_edge)
                case default
                    flx_scalar(prim, f_edge)
            end select
        end subroutine
C       ===============================================================
C       Time stepping
C       ===============================================================
C       Calculate residuals
C       First order euler
        subroutine euler_xp(prim, s, s_edge, w_n, err)
            use input
            implicit none
            real*8, dimension(5, nx), intent(in) :: prim
            real*8, dimension(nx), intent(in) :: s 
            real*8, dimension(nx-1), intent(out) :: s_edge
            real*8, dimension(3, nx), intent(out) :: w_n
            real*8, intent(out) :: err
            real*8, dimension(3, nx) :: w, f, q, r
            real*8, dimension(3, nx-1) :: f_edge
            real*8 :: dt_v
            calc_w(prim, w)
            calc_f(prim, f)
            calc_q(prim, s, q)
            flx_eval(prim, w, f, c, f_edge)
            calc_r(f_edge, s_edge, q, r)
            do i = 2, nx-1
                dt_v = dt(i)/(s(i)*dx)
                do k = 1, 3
                    w_n(k, i) = w(k, i) - dt_v*r(k, i)
                end do
            end do
            calc_err(r, err)
        end subroutine
C       ===============================================================
C       Boundary Conditions
C       ===============================================================
        subroutine update_inlet(prim, c, m)
            use input
            real*8, dimension(5), intent(in) :: prim
            real*8 :: m
        end subroutine
        subroutine update_outlet(prim, c, m)
            use input
        end subroutine
        subroutine update_bc(prim)
            use input
            real*8, dimension(5, nx), intent(out) :: nx
            real*8 :: u, c, m
            calc_c
            u = prim(2, 1)
        end subroutine
C       ===============================================================
C       Main
C       ===============================================================
        program main
            implicit none
            use constants
            use input
C       1. INITIALIZE
C           IMPOSE EXIT STATIC PRESSURE
            integer :: i, idx
            real*8 :: tol, er
            real*8, dimension(nx) :: x, s
            real*8, dimension(5, nx) :: prim
            real*8, dimension(3, nx-1) :: f
            call mkgrid(x, s)
            call init(prim)
            er = 1
            do while (er > tol)
                time_step(prim, s, s_edge, w_n, err)
                calc_prim(w_n, prim)
                update_bc(prim)
                er = calc_err
            end do
C           post process


        end program
