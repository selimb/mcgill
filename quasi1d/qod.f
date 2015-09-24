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
            implicit none
            integer, public, parameter :: nx = 40
            real*8, public, parameter :: dx = 1.0/nx
        end module
C       Make initial grid.
        subroutine mkgrid(x, s)
            use input
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
            use constants
            implicit none
            real*8, dimension(4, nx), intent(out) :: prim
            real*8 :: rh0, p0, t0, u0
            t0 = ttot_in/(1 + 0.5*(g - 1.0)*M_in**2
            p0 = ptot_in/(1 + 0.5*(g - 1.0)*M_in**2)**(g/(g-1))
            rho0 = p0/(R*t0)
            u0 = 0.0
            do i=1, nx
                prim(1,i) = rho0
                prim(2,i) = u0
                prim(3,i) = p0
                rho_e = p/(g - 1)
                prim(4,i) = rho_e/rho
            end do
        end
C       Calculate W
        subroutine calcw(prim, w)
            implicit none
            real*8, dimension(4), intent(in) :: prim
            real*8, dimension(3), intent(out) :: w
            w[1] = prim[1]
            w[2] = prim[1]*prim[2]
            w[3] = prim[4]
        end subroutine
C       Calculate F
        subroutine calf(prim, f)
            implicit none
            real*8, dimension94), intent(in) :: prim
            real*8, dimension(3), intent(out) :: f
            f[1] = prim[1]*prim[2]
            f[2] = prim[1]*prim[2]**2 + prim[3]
            f[3] = (prim[4] + prim[3])*prim[2]
        end subroutine
        program quasi
            implicit none
            use constants
            use input
C       1. INITIALIZE
C           IMPOSE EXIT STATIC PRESSURE
            integer :: i
            real*8, dimension(nx) :: x, s
            call mkgrid(x, s)
            real*8, dimension(4, nx) :: prim
            call init

C       2. ITERATION LOOP


C
        end program
