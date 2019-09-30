!> \file
!> \author Merryn Tawhai, Alys Clark
!> \brief This module handles all geometry read/write/generation.
!>
!> \section LICENSE
!>
!>
!> Contributor(s):
!>
!>\Description
!> This module handles all all geometry read/write/generation.

MODULE indices
  IMPLICIT NONE
  !parameters
  ! indices for elem_ordrs
  INTEGER :: num_ord=3,no_gen=1,no_hord=2,no_sord=3
  ! indices for node_fields
  INTEGER :: num_nj,nj_aw_press,nj_bv_press,nj_conc1,&
       nj_conc2
  ! indices for elem_field
  INTEGER ::num_ne,ne_radius,ne_length,ne_vol,&
       ne_resist,ne_t_resist,ne_Vdot,ne_Vdot0,ne_a_A,&
       ne_dvdt,ne_radius_in,ne_radius_in0,&
       ne_radius_out,ne_radius_out0,ne_group,ne_Qdot
  ! indices for unit_field
  INTEGER :: num_nu,nu_vol,nu_comp,nu_conc2,nu_Vdot0,nu_Vdot1, &
       nu_Vdot2,nu_dpdt,nu_pe,nu_vt,nu_air_press,nu_conc1,nu_vent,&
       nu_vd,nu_perf,nu_blood_press
  !indices for gas exchange field
  ! indices for gasex_field
  INTEGER,PARAMETER :: num_gx = 12
  INTEGER,PARAMETER :: ng_p_alv_o2=1      ! index for alveolar partial pressure of O2
  INTEGER,PARAMETER :: ng_p_alv_co2=2     ! index for alveolar partial pressure of CO2
  INTEGER,PARAMETER :: ng_p_ven_o2=3      ! index for local venous partial pressure of O2
  INTEGER,PARAMETER :: ng_p_ven_co2=4     ! index for local venous partial pressure of CO2
  INTEGER,PARAMETER :: ng_p_cap_o2=5      ! index for local end capillary partial pressure of O2
  INTEGER,PARAMETER :: ng_p_cap_co2=6     ! index for local end capillary partial pressure of CO2
  INTEGER,PARAMETER :: ng_source_o2=7     ! index for source (flux) of O2
  INTEGER,PARAMETER :: ng_source_co2=8    ! index for source (flux) of CO2
  INTEGER,PARAMETER :: ng_Vc=9            ! index for unit's capillary blood volume
  INTEGER,PARAMETER :: ng_sa=10           ! index for unit's capillary surface area
  INTEGER,PARAMETER :: ng_tt=11           ! index for transit time in unit
  INTEGER,PARAMETER :: ng_time=12         ! index for time elapsed for RBC in capillaries

  !model type
  CHARACTER(len=60) :: model_type

  PUBLIC num_ord,no_gen,no_hord,no_sord

  PUBLIC num_nj,nj_aw_press,nj_bv_press,nj_conc1,nj_conc2

  PUBLIC num_ne,ne_radius,ne_length,ne_vol,&
       ne_resist,ne_t_resist,ne_Vdot,ne_Vdot0,ne_a_A,&
       ne_dvdt,ne_radius_in,ne_radius_in0,ne_radius_out,&
       ne_radius_out0,ne_group,ne_Qdot

  PUBLIC num_nu,nu_vol,nu_comp, nu_conc2,nu_Vdot0,nu_Vdot1, &
       nu_Vdot2,nu_dpdt,nu_pe,nu_vt,nu_air_press,&
       nu_conc1,nu_vent,nu_vd,&
       nu_perf,nu_blood_press

  PUBLIC num_gx, ng_p_alv_o2,ng_p_alv_co2,ng_p_ven_o2,ng_p_ven_co2, &
       ng_p_cap_o2, ng_p_cap_co2,ng_source_o2,ng_source_co2, &
       ng_Vc, ng_sa, ng_tt, ng_time


  PUBLIC model_type

  !Interfaces
  PRIVATE
  PUBLIC define_problem_type,ventilation_indices, perfusion_indices, get_ne_radius, get_nj_conc1, &
       growing_indices

CONTAINS

  !> Define problem type
  SUBROUTINE define_problem_type(PROBLEM_TYPE)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_PROBLEM_TYPE" :: DEFINE_PROBLEM_TYPE
    USE other_consts, ONLY: MAX_FILENAME_LEN, MAX_STRING_LEN
    USE diagnostics, ONLY: enter_exit

    CHARACTER(len=MAX_FILENAME_LEN),INTENT(in) :: PROBLEM_TYPE

    CHARACTER(len=60) :: sub_name

    sub_name = 'define_problem_type'
    CALL enter_exit(sub_name,1)
    SELECT CASE (PROBLEM_TYPE)
    CASE ('gas_exchange')
       PRINT *, 'You are solving a gas exchange model, setting up indices'
       CALL exchange_indices
    CASE ('gas_mix')
       PRINT *, 'You are solving a gas mixing model, setting up indices'
       CALL gasmix_indices
    CASE ('gas_transfer')
       PRINT *, 'You are solving a gas transfer model, setting up indices'
       CALL exchange_indices
    CASE ('perfusion')
       PRINT *, 'You are solving a static perfusion model, setting up indices'
       CALL perfusion_indices
    CASE ('ventilation')
       PRINT *, 'You are solving a ventilation model, setting up indices'
       CALL ventilation_indices
    CASE('grow_tree')
       PRINT *, 'You are solving a growing problem, setting up indices'
       CALL growing_indices
    END SELECT
    model_type=TRIM(PROBLEM_TYPE)
    CALL enter_exit(sub_name,2)
  END SUBROUTINE define_problem_type

  !>Gas mixing indices
  SUBROUTINE exchange_indices
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GASMIX_INDICES" :: GASMIX_INDICES

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE
    CHARACTER(len=60) :: sub_name

    sub_name = 'exchange_indices'
    CALL enter_exit(sub_name,1)
    ! indices for elem_ordrs. These dont usually change.
    ! indices for node_field
    num_nj=3
    nj_conc1=2
    nj_conc2=3

    ! indices for elem_field
    num_ne=9
    ne_radius=1
    ne_length=2
    ne_vol=3
    ne_resist=4
    ne_Vdot=5
    ne_Qdot=6
    ne_dvdt=7

    ! indices for unit_field
    num_nu=7
    nu_vol=1
    nu_comp=2
    nu_Vdot0=3
    nu_vd=4
    nu_perf=5
    nu_conc1=6
    nu_conc2=7


    CALL enter_exit(sub_name,2)
  END SUBROUTINE exchange_indices

  !>Gas mixing indices
  SUBROUTINE gasmix_indices
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GASMIX_INDICES" :: GASMIX_INDICES

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE
    CHARACTER(len=60) :: sub_name

    sub_name = 'gasmix_indices'
    CALL enter_exit(sub_name,1)
    ! indices for elem_ordrs. These dont usually change.
    ! indices for node_field
    num_nj=3
    nj_aw_press=2
    nj_conc1=3
    ! indices for elem_field
    num_ne=9
    ne_radius=1
    ne_length=2
    ne_vol=3
    ne_resist=4
    ne_t_resist=5
    ne_Vdot=6
    ne_Vdot0=7
    ne_a_A=8
    ne_dvdt=9
    ! indices for unit_field
    num_nu=11
    nu_vol=1
    nu_comp=2
    nu_Vdot0=3
    nu_Vdot1=4
    nu_Vdot2=5
    nu_dpdt=6
    nu_pe=7
    nu_vt=8
    nu_air_press=9
    nu_conc1=10
    nu_vent=11
    CALL enter_exit(sub_name,2)
  END SUBROUTINE gasmix_indices

  !> Ventilation indices
  SUBROUTINE ventilation_indices
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_VENTILATION_INDICES" :: VENTILATION_INDICES

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE
    CHARACTER(len=60) :: sub_name

    sub_name = 'ventilation_indices'
    CALL enter_exit(sub_name,1)
    ! indices for elem_ordrs. These dont usually change.
    ! indices for node_field
    num_nj=2 !number of nodal fields
    nj_aw_press=2 !air pressure
    ! indices for elem_field
    num_ne=8 !number of element fields
    ne_radius=1 !radius of airway
    ne_length=2 !length of airway
    ne_vol=3 !volume
    ne_resist=4 !resistance of airway
    ne_t_resist=5
    ne_Vdot=6 !Air flow, current time step
    ne_Vdot0=7 !air flow, last timestep
    ne_dvdt=8
    ! indices for unit_field
    num_nu=10
    nu_vol=1
    nu_comp=2
    nu_Vdot0=3
    nu_Vdot1=4
    nu_Vdot2=5
    nu_dpdt=6
    nu_pe=7
    nu_vt=8
    nu_air_press=9
    nu_vent=10
    CALL enter_exit(sub_name,2)
  END SUBROUTINE ventilation_indices
  !
  !########################################################################
  !
  !> Growing indices
  SUBROUTINE growing_indices
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GROWING_INDICES" :: GROWING_INDICES

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE
    CHARACTER(len=60) :: sub_name

    sub_name = 'growing_indices'
    CALL enter_exit(sub_name,1)
    ! indices for elem_ordrs. These dont usually change.
    ! indices for node_field
    num_nj=0 !number of nodal fields
    ! indices for elem_field
    num_ne=2 !number of element fields
    ne_radius=1 !radius of airway
    ne_length=2 !length of airway
    ! indices for unit_field
    num_nu=0
    CALL enter_exit(sub_name,2)
  END SUBROUTINE growing_indices
  !
  !######################################################################
  !
  !> Perfusion indices
  SUBROUTINE perfusion_indices
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_PERFUSION_INDICES" :: PERFUSION_INDICES

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE
    CHARACTER(len=60) :: sub_name

    sub_name = 'perfusion_indices'
    CALL enter_exit(sub_name,1)

    ! indices for node_field
    num_nj=1
    nj_bv_press=1 !pressure in blood vessel
    ! indices for elem_field
    num_ne=9
    ne_radius=1 !strained average radius over whole element
    ne_radius_in=2 !strained radius into an element
    ne_radius_out=3 !strained radius out of an element
    ne_length=4!length of an elevent
    ne_radius_in0=5!unstrained radius into an element
    ne_radius_out0=6!unstrained radius out of an element
    ne_Qdot=7 !flow in an element
    ne_resist=8 !resistance of a blood vessel
    ne_group=9!Groups vessels into arteries (field=0), capillaries (field=1) and veins(field=2)
    !indices for units
    num_nu=2
    nu_perf=1
    nu_blood_press=2

    CALL enter_exit(sub_name,2)
  END SUBROUTINE perfusion_indices

  FUNCTION get_ne_radius() RESULT(res)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GET_NE_RADIUS" :: GET_NE_RADIUS

    USE diagnostics, ONLY: enter_exit

    IMPLICIT NONE
    CHARACTER(len=60) :: sub_name
    INTEGER :: res

    sub_name = 'get_ne_radius'
    CALL enter_exit(sub_name,1)

    res=ne_radius

    CALL enter_exit(sub_name,2)
  END FUNCTION get_ne_radius

  FUNCTION get_nj_conc1() RESULT(res)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GET_NJ_CONC1" :: GET_NJ_CONC1

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE
    CHARACTER(len=60) :: sub_name
    INTEGER :: res

    sub_name = 'get_nj_conc1'
    CALL enter_exit(sub_name,1)

    res = nj_conc1

    CALL enter_exit(sub_name,2)
  END FUNCTION get_nj_conc1

END MODULE indices
