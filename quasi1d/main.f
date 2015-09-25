C       ===============================================================
C       Main
C       ===============================================================
        program main
        use types, only: dp
        use constants
        use inputs, only: read_input_file, nx, tol
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
        real(dp), parameter :: max_iter = 10000
        call read_input_file('input')
        allocate(prim(5, nx))
        allocate(r(3, nx))
        allocate(x(nx))
        allocate(s(nx))
        call mkgrid(x, s)
        call init_state(prim)
        err = 1
        iter = 1
        do while (err > tol .and. iter < max_iter)
            call timestep(prim, s, r)
            err = calc_err(r)
            iter = iter + 1
        end do
        n = size(x)
        open(10, file='output')
        write(10,*) 'x s rho u p e c'
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
