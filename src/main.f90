program MANIAC

    use montecarlo_module
    use parameters_parser
    use output_management
    use tabulated_utils
    use prepare_utils
    use energy_utils
    use input_parser
    use data_parser
    use cli_utils

    implicit none

    ! Step 1 : Program initialization
    call parse_command_line_arguments() ! Handle -i, -d, -p, -r, -o options
    call setup_output_files()           ! Open log file and create output directory

    ! Step 2 : Read input files
    call read_input_file()              ! Read the main MANIAC input file
    call ReadSystemData()               ! Read topology/data file
    call ReadParameters()               ! Read simulation parameters (Lennard-Jones, etc.)

    ! Step 3 : Simulation preparation
    call prepare_simulation_parameters()  ! Set up MC parameters, initial checks
    call PrecomputeTable()              ! Precompute tables for faster calculation
    call ComputeSystemEnergy(primary)   ! Compute initial total energy

    ! Step 4 :Monte Carlo simulation
    call MonteCarloLoop()               ! Main MC loop

    ! Step 5 :Final reporting and cleanup
    call FinalReport()                  ! Print energy and statistics

end program MANIAC
