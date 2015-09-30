C       ===============================================================
C       Picture of grid indices:
C             node        edge
C           |---o---|---o---|---o ... ---|---o---|
C               1   1   2   2   3       J-1  J
C       ===============================================================
C       ===============================================================
C       Schemes for flux evaluation.
C       ===============================================================
        module flx_schemes
        use types, only: dp
        use constants, only: gam
        use inputs, only: params
        use common_calcs
        implicit none
        contains
C       Calculate the diagonalizers S^-1*C^-1 and CS
        pure subroutine calc_diagonalizers(prim, SCinv, CS)
            real(dp), dimension(:), intent(in) :: prim
            real(dp), dimension(3, 3), intent(out) :: SCinv, CS
            real(dp), dimension(3, 3) :: S, CC, Sinv, Cinv
            real(dp) :: beta, alpha, rho, u, c, c2
            S(:, :) = 0
            Sinv(:, :) = 0
            CC(:, :) = 0
            Cinv(:, :) = 0
            beta = gam - 1
            rho = prim(1)
            u = prim(2)
            c = prim(5)
            c2 = c**2
            alpha = 0.5*u**2
            S(1, 1) = 1
            S(2, 1) = -u/rho
            S(3, 1) = alpha*beta
            S(2, 2) = 1/rho
            S(3, 2) = -u*beta
            S(3, 3) = beta
            Sinv(1, 1) = 1
            Sinv(2, 1) = u
            Sinv(3, 1) = alpha
            Sinv(2, 2) = rho
            Sinv(3, 2) = rho*u
            Sinv(3, 3) = 1/beta
            CC(1, 1) = 1
            CC(2, 2) = rho*c
            CC(3, 2) = -rho*c
            CC(1, 3) = -1/(c2)
            CC(2, 3) = 1
            CC(3, 3) = 1
            Cinv(1, 1) = 1
            Cinv(1, 2) = 1/(2*c2)
            Cinv(2, 2) = 1/(2*rho*c)
            Cinv(3, 2) = 0.5
            Cinv(1, 3) = 1/(2*c2)
            Cinv(2, 3) = -1/(2*rho*c)
            Cinv(3, 3) = 0.5
            SCinv = matmul(Sinv, Cinv)
            CS = matmul(CC, S)
        end subroutine calc_diagonalizers

C       Calculate the diagonal matrices Yp and Ym according to the
C       Steger Warming scheme
        pure subroutine calc_diag_sw(prim, Yp, Ym)
            real(dp), dimension(:), intent(in) :: prim
            real(dp), dimension(3, 3), intent(out) :: Yp, Ym
            real(dp), dimension(3) :: lambdas
            real(dp) :: u, c, sq
            integer :: i
            Yp(:, :) = 0
            Ym(:, :) = 0
            u = prim(2)
            c = prim(5)
            lambdas(1) = u
            lambdas(2) = u + c
            lambdas(3) = u - c
            do i = 1, 3
                sq = sqrt(lambdas(i)**2 + params%eps**2)
                Yp(i, i) = 0.5*(lambdas(i) + sq)
                Ym(i, i) = 0.5*(lambdas(i) - sq)
            end do
        end subroutine

C       Steger-Warming
        pure function flx_sw(prim, w, f) result (f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp), dimension(3, 3, size(prim, 2)) :: Ap, Am
            real(dp), dimension(3, 3) :: Yp, Ym, CS, SCinv
            integer :: i, n
            n = size(prim, 2)
            do i = 1, n
                call calc_diagonalizers(prim(:, i), SCinv, CS)
                call calc_diag_sw(prim(:, i), Yp, Ym)
                Ap(:, :, i) = matmul( matmul(SCinv, Yp), CS )
                Am(:, :, i) = matmul( matmul(SCinv, Ym), CS )
            end do
            do i = 1, n - 1
                f_edge(:, i) = matmul(Ap(:, :, i), w(:, i))
     &                         + matmul(Am(:, :, i + 1), w(:, i + 1))
            end do
        end function flx_sw

C       Modified Steger-Warming
        pure function flx_msw(prim, w, f) result (f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp), dimension(3, 3) :: Ap_edge, Am_edge
            real(dp), dimension(5) :: prim_edge
            real(dp), dimension(3, 3) :: Yp, Ym, CS, SCinv
            integer :: i, k, n, m
            n = size(prim, 2)
            m = size(prim, 1)
            do i = 1, n - 1
                do k = 1 , m
                    prim_edge(k) = 0.5*(prim(k, i) + prim(k, i + 1))
                end do
                call calc_diagonalizers(prim_edge, SCinv, CS)
                call calc_diag_sw(prim_edge, Yp, Ym)
                Ap_edge = matmul( matmul(SCinv, Yp), CS )
                Am_edge = matmul( matmul(SCinv, Ym), CS )
                f_edge(:, i) = matmul(Ap_edge, w(:, i))
     &                         + matmul(Am_edge, w(:, i + 1))
            end do
        end function flx_msw

C       Corrected-Modified Steger-Warming
        pure function flx_cmsw(prim, w, f) result(f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp), dimension(3, size(prim, 2) - 1) :: f_msw, f_sw
            real(dp) :: omega, p, pp, dp_dx
            integer :: i, k, n
            f_msw = flx_msw(prim, w, f)
            f_sw = flx_sw(prim, w, f)
            n = size(prim, 2)
            do i = 1, n - 1
                p = prim(3, i)
                pp = prim(3, i + 1)
                dp_dx = (pp - p)/min(p, pp)
                omega = 1/(1 + dp_dx**2)
                do k = 1, 3
                    f_edge(k, i) = omega*f_msw(k, i)
     &                             + (1 - omega)*f_sw(k, i)
                end do
            end do
        end function flx_cmsw

C       Scalar Dissipation
        pure function flx_scalar(prim, w, f) result (f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            real(dp) :: u_avg, c_avg, lambda_max
            integer :: i, k, n
            n = size(prim, 2)
            do i = 1, n - 1
                u_avg = 0.5*(prim(2, i) + prim(2, i+1))
                c_avg = 0.5*(prim(5, i) + prim(5, i+1))
                lambda_max = calc_lambda_max(u_avg, c_avg)
                do k = 1, 3
                    f_edge(k, i) = 0.5*(
     &                  f(k, i) + f(k, i+1)
     &                  - params%eps*lambda_max*(w(k, i+1) - w(k, i))
     &              )
                end do
            end do
        end function flx_scalar
C       Choose a flux evaluation scheme based on input
C       1 : Scalar dissipation
C       2 : Steger Warming
C       3 : Modified Steger Warming
        pure function flx_eval(prim, w, f) result (f_edge)
            real(dp), dimension(:, :), intent(in) :: prim
            real(dp), dimension(3, size(prim, 2)), intent(in) :: w, f
            real(dp), dimension(3, size(prim, 2) - 1) :: f_edge
            select case (params%flx_scheme)
                case (1)
                    f_edge = flx_scalar(prim, w, f)
                case (2)
                    f_edge = flx_sw(prim, w, f)
                case (3)
                    f_edge = flx_msw(prim, w, f)
                case (4)
                    f_edge = flx_cmsw(prim, w, f)
                case default
                    f_edge = flx_scalar(prim, w, f)
            end select
        end function
        end module flx_schemes
