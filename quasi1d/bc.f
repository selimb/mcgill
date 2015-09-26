C       ===============================================================
C       Boundary Conditions
C       ===============================================================
        module bc
        use types, only: dp
        use inputs, only: params
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
            dt_dx_m = dt(m)/params%dx
            rhodiff = rho_m - prim(1, m_1)
            udiff = u_m - prim(2, m_1)
            pdiff = p_m - prim(3, m_1)
            cdiff = c_m - prim(5, m_1)
C           Compute eigenvalues
            lambdas(1) = 0.5*dt_dx_m*udiff
            lambdas(2) = 0.5*dt_dx_m*(udiff + cdiff)
            lambdas(3) = 0.5*dt_dx_m*(udiff - cdiff)
C           Compute characteristic relations
            R(1) = -lambdas(1)*(rhodiff - pdiff/(c_m**2))
            R(2) = -lambdas(2)*(rhodiff + rho_m*c_m*udiff)
            R(3) = -lambdas(3)*(rhodiff - rho_m*c_m*udiff)
C           Compute exit mach number
            M_m = udiff/cdiff
C           Compute dp
            if (M_m > 1) then
                dp = 0.5*(R(2) + R(3))
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
            e_m = rho_m*(cv*t_m + 0.5*u_m**2)
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

