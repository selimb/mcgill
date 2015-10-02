C       ===============================================================
C       Time stepping
C       ===============================================================
        module timestepping
        use types, only: dp
        use inputs, only: params
        use common_calcs
        use flx_schemes, only: flx_eval
        use bc, only: update_bc
        implicit none
        private
        public timestep
        contains

        pure function calc_dt(prim) result(dt)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(size(prim, 2)) :: dt
            real(dp) :: u, c, lambda_max
            integer :: i, n
            n = size(prim, 2)
            do i = 1, n
                u = prim(2, i)
                c = prim(5, i)
                lambda_max = calc_lambda_max(u, c)
                dt(i) = params%dx*params%cfl/lambda_max
            end do
        end function
C       Euler explicit scheme
C       ---------------------------------------------------------------
        subroutine euler_xp(prim, s, r)
            real(dp), dimension(:, :), intent(inout) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: s
            real(dp), dimension(3, size(prim, 2)), intent(out) :: r
            real(dp), dimension(3, size(prim, 2)) :: w, f, q, w_n
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp), dimension(size(prim, 2)) :: dt
            real(dp) :: dt_v
            integer :: i, k, n
C           Compute time step
            dt = calc_dt(prim)
            w = calc_w(prim)
            f = calc_f(prim)
            q = calc_q(prim, s)
C           Compute flux across edges
            f_edge = flx_eval(prim, w, f)
C           Compute residual
            r = calc_r(s, f_edge, q)
C           Update W
            n = size(prim, 2)
            do i = 2, n - 1
                dt_v = dt(i)/(s(i)*params%dx)
                do k = 1, 3
                    w_n(k, i) = w(k, i) - dt_v*r(k, i)
                end do
            end do
C           Update BCs
            call update_bc(prim, dt)
C           Update state vector
            prim(:, 2:n-1) = calc_prim(w_n(:, 2:n-1))
        end subroutine

C       Runge Kutta Fourth-Order
C       ---------------------------------------------------------------
        subroutine rk_step(prim, w0, s, dt_2,
     &          prim_new, w_new, r_new
     &  )
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w0
            real(dp), dimension(size(prim, 2)), intent(in) :: s, dt_2
            real(dp), dimension(size(prim, 1), size(prim, 2)),
     &          intent(out) :: prim_new
            real(dp), dimension(3, size(prim, 2)), intent(out) ::
     &          w_new, r_new
            real(dp), dimension(3, size(prim, 2)) :: w, f, q
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            integer :: i, k, n
            n = size(prim, 2)
            w = calc_w(prim)
            f = calc_f(prim)
            q = calc_q(prim, s)
            f_edge = flx_eval(prim, w, f)
            r_new = calc_r(s, f_edge, q)
            do i = 2, n - 1
                do k = 1, 3
                    w_new(k, i) = w0(k, i) - dt_2(i)*r_new(k, i)
                end do
            end do
            w_new(:, 1) = w(:, 1)
            w_new(:, n) = w(:, n)
            prim_new = calc_prim(w_new)
        end subroutine
        subroutine rk4(prim, s, r_new)
            real(dp), dimension(:, :), intent(inout) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: s
            real(dp), dimension(3, size(prim, 2)), intent(out) :: r_new
            real(dp), dimension(3, size(prim, 2)) :: w0, w1, w2, w3,
     &          w_new
            real(dp), dimension(3, size(prim, 2)) :: r1, r2, r3
            real(dp), dimension(size(prim, 1), size(prim, 2)) ::
     &          prim1, prim2, prim3
            real(dp), dimension(size(prim, 2)) :: dt, dt_2
            integer :: n
            n = size(prim, 2)
            dt = calc_dt(prim)
            dt_2 = 0.5*dt
            w0 = calc_w(prim)
            call rk_step(prim, w0, dt_2, s, prim1, w1, r1)
            call rk_step(prim1, w0, dt_2, s, prim2, w2, r2)
            call rk_step(prim2, w0, dt_2, s, prim3, w3, r3)
            w_new = (1/6)*(w0 + 2*w1 + 2*w2 + w3)
            r_new = (1/6)*(r1 + r2 + r3)
            call update_bc(prim, dt)
            prim(:, 2:n-1) = calc_prim(w_new(:, 2:n-1))
        end subroutine
C       Choose a timestepping scheme based on input
C       1 : Euler Explicit
        subroutine timestep(prim, s, r)
            real(dp), dimension(:, :), intent(inout) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: s
            real(dp), dimension(3, size(prim, 2)), intent(out) :: r
            select case (params%timestep_scheme)
                case (1)
                    call euler_xp(prim, s, r)
                case (2)
                    call rk4(prim, s, r)
                case default
                    call euler_xp(prim, s, r)
            end select
        end subroutine
        end module
