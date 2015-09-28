C       ===============================================================
C       Note on implementation.
C           We pass around `prim` between functions. This
C           corresponds to (rho, u, p e, c), the primitive variables.
C           Generating W and F is easy, and the primitive variables
C           are readily available when necessary, such as
C           when plotting or calculating the maximum eigenvalues.
C       ===============================================================
C       Main
C       ===============================================================
        program main
        use types, only: dp
        use utils, only: write1d
        use constants
        use inputs, only: read_input_file, params
        use setup, only: init_state, mkgrid
        use common_calcs, only: calc_err
        use timestepping, only: timestep
        implicit none
        real(dp), dimension(:), allocatable :: x, s
        real(dp), dimension(:, :), allocatable :: prim, r
        real(dp) :: err
        integer :: iter, i, k, n
        character(len=40), parameter :: fmt_ = 'EN20.8)'
        character(len=40) :: fmt1 = '(' // fmt_
        character(len=40) :: fmt5 = '(5' // fmt_
        call read_input_file('flow.prm')
        allocate(prim(5, params%nx))
        allocate(r(3, params%nx))
        allocate(x(params%nx))
        allocate(s(params%nx))
        call mkgrid(x, s)
        call init_state(prim)
        err = 1
        iter = 0
        do while (err > params%tol .and. iter < params%max_iter)
            iter = iter + 1
            call timestep(prim, s, r)
            err = calc_err(r)
        end do
        write(*,*) 'Number of iterations:'
        write(*,*) iter
        write(*,*) 'Error'
        write(*,*) err
        n = size(x)
        open(10, file='output')
        write(10,*) 'x s rho u p e c'
        do i = 1, n
            write(10, fmt1, advance='no') x(i)
            write(10, fmt1, advance='no') s(i)
            write(10, fmt5, advance='no') (prim(k, i), k=1,5)
            write(10, *) ''
        end do
C       TODO post process
        end program
