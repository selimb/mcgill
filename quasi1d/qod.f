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
            real(dp), public, parameter :: gam = 1.4_dp
            real(dp), public, parameter :: Rgas = 1716
            real(dp), public, parameter :: cv = Rgas/(gam - 1)
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
C       ===============================================================
C       Initialization
C       ===============================================================
        module setup
        use types, only: dp
        implicit none
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
            use constants, only: gam, M_in, ptot_in, ttot_in, Rgas
            use input, only: p_exit
            real(dp), dimension(:, :), intent(out) :: prim
            real(dp) :: rho0, p0, t0, u0, rho_e, M_term
            integer :: i, n
            n = size(prim, 2)
            M_term = ( 1.0_dp + 0.5_dp*(gam - 1.0_dp)*M_in**2 )
            t0 = ttot_in/M_term
            p0 = ptot_in/( M_term**(gam/(gam - 1.0_dp)) )
            rho0 = p0/(Rgas*t0)
            u0 = 0
            rho_e = p0/(gam - 1.0_dp)
            do i = 1, n
                prim(1, i) = rho0
                prim(2, i) = u0
                prim(3, i) = p0
                prim(4, i) = rho_e/rho0
                prim(5, i) = sqrt(gam*p0/rho0)
            end do
C           Impose static pressure
            prim(3, n) = p_exit
        end
        end module setup
C       ===============================================================
C       Helper functions to calculate commonly required quantities
C       ===============================================================
        module common_calcs
        use types, only: dp
        implicit none
        contains 

C       Calculate prim from W
        pure function calc_prim(w) result(prim)
            use constants, only: gam
            real(dp), dimension(:, :), intent(in) :: w
            real(dp), dimension(5, size(w, 2)) :: prim
            real(dp) :: rho, u, p, e
            integer :: i, n
            n = size(w, 2)
            do i = 1, n
                rho = w(1, i)
                u = w(2, i)/rho
                e = w(3, i)
                E = e/rho
                p = (gam - 1.0_dp)*rho*(E - 0.5_dp*u**2)
                prim(1, i) = rho
                prim(2, i) = u
                prim(3, i) = p
                prim(4, i) = e
                prim(5, i) = sqrt(gam*p/rho)
            end do
        end function
C       Calculate W from prim
        pure function calc_w(prim) result(w)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)) :: w
            integer :: i, n
            n = size(prim, 2)
            do i = 1, n
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
            integer :: i, n
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
            integer :: i, n
            n = size(prim, 2)
            do i=2, n
                q(1, i) = 0
                q(2, i) = prim(3, i)*(s(i) - s(i-1))
                q(3, i) = 0
            end do
        end function
C       Calculate maximum eigenvalue
        elemental function calc_lambda_max(u, c) result(lambda_max)
            real(dp), intent(in) :: u, c
            real(dp) :: lambda_max
            lambda_max = max(u, u + c, u - c)
        end function
C       Calculate residuals
        pure function calc_r(s, f_edge, q) result(r)
            real(dp), dimension(:), intent(in) :: s
            real(dp), dimension(3, size(s)-1), intent(in) :: f_edge
            real(dp), dimension(3, size(s)), intent(in) :: q
            real(dp), dimension(3, size(s)) :: r
            real(dp), dimension(size(s)-1) :: s_edge
            integer :: i, k, n
            n = size(q, 2)
            do i = 1, n - 1
                s_edge(i) = 0.5_dp*(s(i) + s(i + 1))
            end do
            do i = 2, n-1
            do k = 1, 3
                r(k, i) = f_edge(k, i)*s_edge(i)
     &                    - f_edge(k, i - 1)*s_edge(i - 1)
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
        implicit none
        contains 

C       Scalar Dissipation
        pure function flx_scalar(prim, w, f) result (f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp) :: u_avg, c_avg, lambda_max
            integer :: i, k, n
            n = size(prim, 2)
            do i = 1, n - 1
                u_avg = 0.5_dp*(prim(2, i) + prim(2, i+1))
                c_avg = 0.5_dp*(prim(5, i) + prim(5, i+1))
                lambda_max = calc_lambda_max(u_avg, c_avg)
                do k = 1, 3
                    f_edge(k, i) = 0.5_dp*(
     &                  f(k, i) + f(k, i+1)
     &                  - eps*lambda_max*(w(k, i+1) - w(k, i)))
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
C       Boundary Conditions
C       ===============================================================
        module bc
        use types, only: dp
        use input, only: dx
        implicit none
        private 
        public update_bc
        contains 

        subroutine update_inlet(prim_in, dt)
            use constants, only: M_in
            real(dp), dimension(:), intent(inout) :: prim_in
            real(dp), intent(in) :: dt
            if (M_in > 1) then
                return
            end if
            write (*,*) "Not Implemented."
            call exit(1)
        end subroutine
        subroutine update_outlet(prim, dt)
            use constants, only: gam, Rgas, cv
            real(dp), dimension(:, :), intent(inout) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: dt
            real(dp), dimension(3) :: lambdas, R
            real(dp), parameter :: p5 = 0.5_dp
            real(dp) :: rho_m, u_m, p_m, e_m, T_m, c_m, M_m, dt_dx_m
            real(dp) :: rhodiff, udiff, pdiff, cdiff
            real(dp) :: drho, du, dp
            integer :: m, m_1
            m = size(prim, 2)
            m_1 = m - 1
            rho_m = prim(1, m)
            u_m = prim(2, m)
            p_m = prim(3, m)
            e_m = prim(4, m)
            c_m = prim(5, m)
            dt_dx_m = dt(m)/dx
            rhodiff = rho_m - prim(1, m_1)
            udiff = u_m - prim(2, m_1)
            pdiff = p_m - prim(3, m_1)
            cdiff = c_m - prim(5, m_1)
C           Compute eigenvalues
            lambdas(1) = p5*dt_dx_m*udiff
            lambdas(2) = p5*dt_dx_m*(udiff + cdiff)
            lambdas(3) = p5*dt_dx_m*(udiff - cdiff)
C           Compute characteristic relations
            R(1) = -lambdas(1)*(rhodiff - pdiff/(c_m**2))
            R(2) = -lambdas(2)*(rhodiff + rho_m*c_m*udiff)
            R(3) = -lambdas(3)*(rhodiff - rho_m*c_m*udiff)
C           Compute exit mach number
            M_m = udiff/cdiff
C           Compute dp
            if (M_m > 1) then
                dp = p5*(R(2) + R(3))
            else
                dp = 0
            end if 
C           Update drho and du
            drho = R(1) + dp/(c_m**2)
            du = (R(2) - dp)/(rho_m*c_m)
C           Update flow propeties
            rho_m = rho_m + drho
            u_m = u_m + du
            p_m = p_m + dp
            t_m = p_m/(rho_m*Rgas)
            e_m = rho_m*(cv*t_m + p5*u_m**2)
            c_m = sqrt(gam*p_m/rho_m)
        end subroutine
        subroutine update_bc(prim, dt)
            real(dp), dimension(:, :), intent(inout) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: dt
            integer :: n
            n = size(prim, 2)
            call update_inlet(prim(:, 1), dt(1))
            call update_outlet(prim, dt)
        end subroutine
        end module

C       ===============================================================
C       Time stepping
C       ===============================================================
C       Calculate residuals
C       First order euler
        module timestepping
        use types, only: dp
        use input, only: dx, cfl, timestep_scheme
        use common_calcs
        use flx_schemes, only: flx_eval
        use bc, only: update_bc
        implicit none
        private 
        public timestep
        contains 

C       Euler explicit scheme
C       ---------------------------------------------------------------
        pure function euler_xp_dt(prim) result(dt)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(size(prim, 2)) :: dt
            real(dp) :: u, c, lambda_max 
            integer :: i, n
            n = size(prim, 2)
            do i = 1, n
                u = prim(2, i)
                c = prim(5, i)
                lambda_max = calc_lambda_max(u, c)
                dt(i) = dx*cfl/lambda_max
            end do
        end function
        subroutine euler_xp(prim, s, err)
            real(dp), dimension(:, :), intent(inout) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: s 
            real(dp), intent(out) :: err
            real(dp), dimension(3, size(prim, 2)) :: w, f, q, r, w_n
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp), dimension(size(prim, 2)) :: dt, lambdas
            real(dp) :: dt_v
            integer :: i, k, n
C           Compute time step
            dt = euler_xp_dt(prim)
            w = calc_w(prim)
            f = calc_f(prim)
            q = calc_q(prim, s)
C           Compute flux across edges
            f_edge = flx_eval(prim, w, f)
C           Compute residual
            r = calc_r(s, f_edge, q)
C           Update W
            n = size(prim, 2)
            do i = 2, n-1
                dt_v = dt(i)/(s(i)*dx)
                do k = 1, 3
                    w_n(k, i) = w(k, i) - dt_v*r(k, i)
                end do
            end do
C           Update BCs
            call update_bc(prim, dt)
C           Update state vector
            prim(:, 2:n-1) = calc_prim(w_n(:, 2:n-1))
C           Calculate residual
            err = calc_err(r)
        end subroutine
        subroutine timestep(prim, s, err)
            use input, only: timestep_scheme
            real(dp), dimension(:, :), intent(inout) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: s 
            real(dp), intent(out) :: err
            select case (timestep_scheme)
                case (1)
                    call euler_xp(prim, s, err)
                case default
                    call euler_xp(prim, s, err)
            end select
        end subroutine
        end module
C       ===============================================================
C       Main
C       ===============================================================
        program main
        use types, only: dp
        use constants
        use input, only: nx, tol
        use setup, only: init_state, mkgrid
        use timestepping, only: timestep
        implicit none
        integer :: iter, i, k, n
        real(dp) :: err
        real(dp), parameter :: max_iter = 100
        real(dp), dimension(nx) :: x, s
        real(dp), dimension(5, nx) :: prim
        character(len=40), parameter :: fmt_ = 'EN20.8)'
        character(len=40) :: fmt1 = '(' // fmt_
        character(len=40) :: fmt5 = '(5' // fmt_
        call mkgrid(x, s)
        call init_state(prim)
        err = 1
        iter = 0
        do while (err > tol .and. iter < max_iter)
            call timestep(prim, s, err)
            iter = iter + 1
        end do
        n = size(x)
        open(10, file='output')
        write(10,*) 'vars=x,s,rho,u,p,e,c'
        do i = 1, n
            write(10, fmt1, advance='no') x(i)
            write(10, fmt1, advance='no') s(i)
            write(10, fmt5, advance='no') (prim(k, i), k=1,5)
            write(10, *) ''
        end do
        write (*, *) iter
        write (*, *) err
C       TODO post process
        end program
