program MANIAC

    use montecarlo_module
    use parameters_parser
    use output_management
    use tabulated_utils
    use prepare_utils
    use prescan_files
    use energy_utils
    use input_parser
    use data_parser
    use cli_utils

    implicit none

    integer :: i, j

    ! Step 1 : Program initialization
    call parse_command_line_arguments() ! Handle -i, -d, -p, -r, -o options
    call setup_output_files()           ! Open log file and create output directory

    ! Step 2 : Read input files
    call prescan_inputs()               ! Prescan input files to anticipate allocating size



    write(*,*) "=========================================="
    write(*,*) "           RESIDUE TYPE SUMMARY           "
    write(*,*) "=========================================="
    write(*,*)
    do i = 1, size(res_infos)
        write(*,'(A,I4)')   " Residue ID:                   ", i
        write(*,'(A,I4)')   "   Number of atom types:      ", res_infos(i)%nb_types
        write(*,'(A,I6)')   "   Total atoms in residue:    ", res_infos(i)%nb_atoms
        write(*,'(A,I6,I6)') "   Residue counts (main/resv): ", res_infos(i)%nb_res(1), res_infos(i)%nb_res(2)

        if (res_infos(i)%nb_types > 0) then
            write(*,'(A)', advance='no') "   Atom types:                "
            do j = 1, res_infos(i)%nb_types
                write(*,'(I6)', advance='no') res_infos(i)%types(j)
            end do
            write(*,*)
        end if

        write(*,*)
    end do

    write(*,*) "=========================================="


    write(*,*) "=========================================="












!    call read_input_file()              ! Read the main MANIAC input file











    ! call ReadSystemData()               ! Read topology/data file
    ! call ReadParameters()               ! Read simulation parameters (Lennard-Jones, etc.)

    ! ! Step 3 : Simulation preparation
    ! call prepare_simulation_parameters()  ! Set up MC parameters, initial checks
    ! call PrecomputeTable()              ! Precompute tables for faster calculation
    ! call compute_system_energy(primary)   ! Compute initial total energy

    ! ! Step 4 :Monte Carlo simulation
    ! call MonteCarloLoop()               ! Main MC loop

    ! ! Step 5 :Final reporting and cleanup
    ! call FinalReport()                  ! Print energy and statistics

end program MANIAC
