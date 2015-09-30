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
                dt(i) = params%dx*params%cfl/lambda_max
            end do
        end function
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
C           Calculate residual
        end subroutine
        subroutine timestep(prim, s, r)
            real(dp), dimension(:, :), intent(inout) :: prim
            real(dp), dimension(size(prim, 2)), intent(in) :: s
            real(dp), dimension(3, size(prim, 2)), intent(out) :: r
            select case (params%timestep_scheme)
                case (1)
                    call euler_xp(prim, s, r)
                case default
                    call euler_xp(prim, s, r)
            end select
        end subroutine
        end module
