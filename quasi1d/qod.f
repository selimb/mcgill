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
        subroutine init(w)
            use constants
            implicit none
            real*8
        end
        program quasi
            implicit none
            use constants
            use input
C       1. INITIALIZE
C           SETUP THE GRID
C           INITIALIZE STATE VECTOR
C               - Density
C               - Momentum
C               - Energy
C           IMPOSE EXIT STATIC PRESSURE
            integer :: i
            real*8, dimension(nx) :: x, s
            call mkgrid(x, s)
            call init
            open(3, file='bump3dgrid.p3dfmt',form='formatted')
            write(3,*) (x(i),i=1,nx),
       &               (s(i),i=1,nx)
C       2. ITERATION LOOP


C
        end program
