module montecarlo_module

    use molecule_translation    ! Provides routines to perform translational moves on molecules
    use monte_carlo_utils       ! Provides Monte Carlo-specific utilities (e.g., step adjustment)
    use molecule_creation       ! Provides routines to insert new molecules into the system
    use molecule_deletion       ! Provides routines to remove molecules from the system
    use molecule_rotation       ! Provides routines to perform rotational moves on molecules
    use simulation_state        ! Stores and manages the current state of the simulation
    use molecule_swap           ! Provides routines to swap molecules
    use output_utils            ! Handles updating files and printing simulation status
    use random_utils            ! Provides random number generation routines
    use write_utils             ! Utilities for writing formatted output or logs
    use constants               ! Defines physical and simulation constants
    use widom                   ! Provides routines for widom insertion

    use, intrinsic :: iso_fortran_env, only: real64 ! Ensures consistent 64-bit real precision

    implicit none

contains

    !------------------------------------------------------------------------------
    ! Performs the main Monte Carlo simulation loop, executing a series of molecular
    ! moves and exchanges over multiple blocks and steps.
    !------------------------------------------------------------------------------
    subroutine monte_carlo_loop()

        implicit none

        ! Local Variables
        real(real64) :: random_draw     ! Random number used for move selection
        integer :: residue_type         ! Index of molecule to be moved
        integer :: molecule_index       ! Index of molecule copy

        ! Initialization
        call LogStartMC()               ! Log starting message
        call update_output_files(.false.)       ! Write initial topology
        status%block = 1 ! Initialize Monte Carlo counters
        status%step  = 1 ! Initialize Monte Carlo counters

        !----------------------------------------------
        ! Main Monte Carlo Loop
        !----------------------------------------------
        do
            ! Pick a molecule type and instance
            residue_type    = PickRandomResidueType(input%is_active)
            molecule_index  = PickRandomMoleculeIndex(primary%num_residues(residue_type))

            ! Perform Monte Carlo move
            random_draw = rand_uniform()

            if (random_draw <= proba%translation) then

                ! Case 1: Small translation move
                call attempt_translation_move(residue_type, molecule_index)

            else if (random_draw <= proba%rotation + proba%translation) then

                ! Case 2: Rotation move
                call attempt_rotation_move(residue_type, molecule_index)

            else if (random_draw <= proba%rotation + proba%translation + proba%swap) then

                ! Case 3: Swap move
                call attempt_swap_move(residue_type, molecule_index)

            else

                ! Insertion/deletion move or Widom
                if (proba%insertion_deletion > 0) then
                    if (rand_uniform() <= PROB_CREATE_DELETE) then

                        ! Attempt to create a molecule
                        molecule_index = primary%num_residues(residue_type) + 1
                        call attempt_creation_move(residue_type, molecule_index)

                    else

                        ! Attempt to delete a molecule
                        call attempt_deletion_move(residue_type, molecule_index)

                    end if
                else

                    ! Widom trial
                    molecule_index = primary%num_residues(residue_type) + 1
                    call widom_trial(residue_type, molecule_index)

                end if
            end if

            !----------------------------------------------
            ! Adjust Monte Carlo step counters
            !----------------------------------------------
            status%step = status%step + 1

            ! Finished all steps in the current block?
            if (status%step > status%desired_step) then

                status%step  = 1                  ! reset step
                status%block = status%block + 1   ! move to next block

                !----------------------------------------------
                ! Adjust step sizes & output status at the end of the block
                !----------------------------------------------
                call adjust_move_step_sizes()        ! Adjust MC step sizes
                call PrintStatus()                ! Print current simulation status
                call update_output_files(.true.)  ! Update output files

            end if

            ! Finished all blocks?
            if (status%block > status%desired_block) exit

        end do

    end subroutine monte_carlo_loop

end module montecarlo_module
