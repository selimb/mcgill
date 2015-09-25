C       ===============================================================
C       Read input file
C       ===============================================================
        module inputs
            use types, only: dp
            use constants, only: ptot_in
            implicit none
C           User input
            integer, public :: flx_scheme
            integer, public :: timestep_scheme
            integer, public :: nx
            real(dp), private :: p_exit_ratio
            real(dp), public :: eps
            real(dp), public :: tol
            real(dp), public :: cfl
C           Calculations
            real(dp), public :: dx
            real(dp), public :: p_exit
            contains
            subroutine read_input_file(filename)
                character(len=*), intent(in) :: filename
                character(len=50) :: line, var_name, val
                logical :: back
                integer :: stat, idx, l
                back = .false.
                open(10, file=filename)
                do
                    read(10, "(a)", iostat=stat) line
                    if (stat < 0) exit
                    l = len_trim(line)
                    idx = scan(line, '=', back)
                    if (idx < 1) then
                        cycle
                    end if
                    var_name = line(1 : idx - 1)
                    val = line(idx+1 : l)
                    select case (var_name)
                        case ("flx_scheme")
                            read (val, *) flx_scheme
                        case ("timestep_scheme")
                            read (val, *) timestep_scheme
                        case ("nx")
                            read (val, *) nx
                        case ("p_exit_ratio")
                            read (val, *) p_exit_ratio
                        case ("eps")
                            read (val, *) eps
                        case ("tol")
                            read (val, *) tol
                        case ("cfl")
                            read (val, *) cfl
                        case default
                    end select
                end do
                dx = 1.0_dp/nx
                p_exit = p_exit_ratio*ptot_in            
            end subroutine
        end module
