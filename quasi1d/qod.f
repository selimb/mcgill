        subroutine mkgrid(nx, x, s)
C       Make initial grid.        
           integer, intent(in) :: nx
           real*8, dimension(nx), intent(out) :: x, s
           real*8 :: h, t1, t2, pi, dx
           pi = 4*atan(1.0_8)
           h = 0.025
           t1 = 0.8
           t2 = 3
           dx = 1.0/nx
           write(*,*) nx
           write(*,*) pi
           do i=1, nx
                x(i) = i*dx
                s(i) = 1 - h*(sin(pi*x(i)**t1))**t2
           end do
           write(*,*) x(2)
           write(*,*) s(2)
        end

        program quasi
        implicit none

C       1. INITIALIZE
C           SETUP THE GRID
C           INITIALIZE STATE VECTOR
C               - Density
C               - Momentum
C               - Energy
C           IMPOSE EXIT STATIC PRESSURE

        integer :: nx, i
        real*8 :: dx, pi, g, t0t, p0t, R, M_in
        real*8, dimension(:), allocatable :: x, s

C       FLOW CONSTANTS
        t0t = 531.2
        g = 1.4
        M_in = 1.2
        R = 1716

C       GENERATE GRID
        nx = 40
        allocate(x(nx), s(nx))
        call mkgrid(nx, x, s)
        open(3, file='bump3dgrid.p3dfmt',form='formatted')
        write(3,*) (x(i),i=1,nx),
     &             (s(i),i=1,nx)
        deallocate (x, s)

C       2. ITERATION LOOP


C
        end program
