C       ===============================================================
C       Schemes for flux evaluation.
C       ===============================================================
        module flx_schemes
        use types, only: dp
        use inputs, only: flx_scheme, eps
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
                u_avg = 0.5*(prim(2, i) + prim(2, i+1))
                c_avg = 0.5*(prim(5, i) + prim(5, i+1))
                lambda_max = calc_lambda_max(u_avg, c_avg)
                do k = 1, 3
                    f_edge(k, i) = 0.5*(
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
