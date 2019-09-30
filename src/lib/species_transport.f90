MODULE species_transport
  !*Brief Description:* This module contains all the subroutines common
  !to species transport models, this includes gas exchange, gas mixing,
  !and particle transport models
  !*LICENSE:*
  !
  !
  !
  !*Full Description:*
  !More info on what the module does if necessary
  !
  USE other_consts
  IMPLICIT NONE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  PRIVATE
  PUBLIC initialise_transport

CONTAINS
  !
  !##############################################################################
  !
  SUBROUTINE initialise_transport()
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIALISE_TRANSPORT" :: INITIALISE_TRANSPORT
    USE indices
    USE arrays, ONLY: dp
    USE gas_exchange, ONLY: initial_gasexchange
    USE diagnostics, ONLY: enter_exit

    !local variables

    CHARACTER(len=60) :: sub_name

    sub_name = 'initialise_transport'
    CALL enter_exit(sub_name,1)

    CALL allocate_memory_speciestrans

    SELECT CASE (model_type)
    CASE ('gas_exchange')
       PRINT *, 'You are solving a gas exchange model'
       !Note that as V, Q are prerequisites something needs to be added here that checks
       !these have been read in and if not sets up linear gradient based on some default parameters
    CASE ('gas_mix')
       PRINT *, 'You are solving a gas mixing model'
       !Note that as V is prerequisites something needs to be added here that checks
       !these have been read in and if not sets up linear gradient based on some default parameters
    CASE ('gas_transfer')
       PRINT *, 'You are solving a gas transfer model'
       !Note that as V, Q are prerequisites something needs to be added here that checks
       !these have been read in and if not sets up linear gradient based on some default parameters
       !note a linear q gradient should  be set up to scale for shunt fraction automatically
       CALL initial_gasexchange(149.0_dp)
       CALL solve_transport

    END SELECT
    CALL enter_exit(sub_name,2)
  END SUBROUTINE initialise_transport

  !
  !##############################################################################
  !
  SUBROUTINE solve_transport()
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIALISE_TRANSPORT" :: INITIALISE_TRANSPORT
    USE indices
    USE arrays, ONLY: dp
    USE gas_exchange, ONLY: steadystate_gasexchange
    USE diagnostics, ONLY: enter_exit

    !local variables
    REAL(dp) c_art_o2, c_ven_o2,p_art_co2,p_art_o2, p_ven_co2,p_ven_o2

    CHARACTER(len=60) :: sub_name

    sub_name = 'solve_transport'
    CALL enter_exit(sub_name,1)


    SELECT CASE (model_type)
    CASE ('gas_exchange')
       PRINT *, 'Nothing implemented'
       !Note that as V, Q are prerequisites something needs to be added here that checks
       !these have been read in and if not sets up linear gradient based on some default parameters
    CASE ('gas_mix')
       PRINT *, 'Nothing implemented'
       !Note that as V is prerequisites something needs to be added here that checks
       !these have been read in and if not sets up linear gradient based on some default parameters
    CASE ('gas_transfer')
       PRINT *, 'Calling gas transfer model '
       p_art_co2=40.0_dp
       p_ven_co2=45.0_dp
       p_art_o2=100.0_dp
       p_ven_o2=40.0_dp
       CALL steadystate_gasexchange(c_art_o2,c_ven_o2,&
            p_art_co2,p_art_o2,149.0_dp,p_ven_co2,p_ven_o2,0.03_dp,&
            0.8_dp*(260.0_dp*1.0e+3_dp/60.0_dp),260.0_dp*1.0e+3_dp/60.0_dp )

    END SELECT
    CALL enter_exit(sub_name,2)
  END SUBROUTINE solve_transport


  !
  !###########################################################################################
  !
  SUBROUTINE allocate_memory_speciestrans()
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIALISE_TRANSPORT" :: INITIALISE_TRANSPORT
    USE indices
    USE arrays, ONLY: dp,gasex_field,num_units
    USE diagnostics, ONLY: enter_exit

    CHARACTER(len=60) :: sub_name
    sub_name = 'allocate_memory_speciestrans'
    CALL enter_exit(sub_name,1)

!!! allocate memory for the gasex_field array, if not already allocated
    IF(.NOT.ALLOCATED(gasex_field)) ALLOCATE(gasex_field(num_gx,num_units))


    CALL enter_exit(sub_name,2)
  END SUBROUTINE allocate_memory_speciestrans

END MODULE species_transport
