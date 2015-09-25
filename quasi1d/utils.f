C       ===============================================================
C       Utility functions
C       ===============================================================
        module utils
        use types, only: dp
        implicit none
        contains

        subroutine writeall(w, f, q, f_edge, r)
            real(dp), dimension(:, :), intent(in) :: w, f, q, f_edge, r
            integer :: i, n
            write(*, *) 'W1'
            write(*, *) (w(1, i), i=1,n) 
            write(*, *) 'W2'
            write(*, *) (w(2, i), i=1,n) 
            write(*, *) 'W3'
            write(*, *) (w(3, i), i=1,n) 
            write(*, *) 'F1'
            write(*, *) (f(1, i), i=1,n) 
            write(*, *) 'F2'
            write(*, *) (f(2, i), i=1,n) 
            write(*, *) 'F3'
            write(*, *) (f(3, i), i=1,n) 
            write(*, *) 'Q2'
            write(*, *) (q(2, i), i=1,n)
            write(*, *) 'Fedge1'
            write(*, *) (f_edge(1, i), i=1,n) 
            write(*, *) 'Fedge2'
            write(*, *) (f_edge(2, i), i=1,n) 
            write(*, *) 'Fedge3'
            write(*, *) (f_edge(3, i), i=1,n) 
            write(*, *) 'r1'
            write(*, *) (r(1, i), i=1,n) 
            write(*, *) 'r2'
            write(*, *) (r(2, i), i=1,n) 
            write(*, *) 'r3'
            write(*, *) (r(3, i), i=1,n) 
            write(*, *) ''
            write(*, *) ''
            write(*, *) ''
            write (*, *) 'Iteration complete'
            write(*, *) ''
            write(*, *) ''
            write(*, *) ''
        end subroutine
        end module
