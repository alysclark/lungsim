MODULE ventilation
  !*Brief Description:* This module handles all code specific to simulating ventilation
  !
  !*LICENSE:*
  !TBC
  !
  !
  !*Full Description:*
  !
  ! This module handles all code specific to simulating ventilation

  IMPLICIT NONE
  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  PRIVATE
  PUBLIC evaluate_vent
  PUBLIC evaluate_uniform_flow
  PUBLIC two_unit_test
  PUBLIC sum_elem_field_from_periphery

CONTAINS

!!!###################################################################################
  !*evaluate_vent:* Sets up and solves venilation model
  SUBROUTINE evaluate_vent
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_VENT" :: EVALUATE_VENT
    USE arrays,ONLY: dp,elem_field,elem_units_below,num_elems,num_units,units,unit_field
    USE indices,ONLY: ne_Vdot,ne_Vdot0,ne_t_resist,nu_comp,nu_dpdt,nu_pe,nu_vt
    USE exports,ONLY: export_1d_elem_field,export_terminal_solution
    USE other_consts !PI
    USE diagnostics, ONLY: enter_exit
    USE geometry,ONLY: set_initial_volume,volume_of_mesh
    IMPLICIT NONE

    ! Local variables
    INTEGER :: Gdirn,iter_step,n,ne,num_brths,num_itns,nunit
    REAL(dp) :: ChestWallRestVol,chest_wall_compliance,constrict,COV,&
         dpmus,dt,endtime,err_est,err_tol,FRC,i_to_e_ratio,init_vol,&
         last_vol,now_vol,Pcw,p_mus,pmus_factor_in,pmus_factor_ex,pmus_step,&
         ppl_current,pptrans,press_in,prev_flow,ptrans_frc,refvol,RMaxMean,&
         RMinMean,sum_dpmus,sum_dpmus_ei,sum_expid,sum_tidal,Texpn,T_interval,&
         time,Tinsp,totalc,Tpass,ttime,undef,volume_target,volume_tree,WOBe,&
         WOBr,WOBe_insp,WOBr_insp,WOB_insp
    CHARACTER :: expiration_type*(10)
    LOGICAL :: CONTINUE,converged

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'evaluate_vent'
    CALL enter_exit(sub_name,1)

    init_vol = 1.0_dp

!!! -------------  DESCRIPTION OF IMPORTANT VARIABLES ---------------
!!! pmus_factor (_in and _ex) are used to scale the driving pressureS, to converge
!!! the tidal volume and expired volume to the target volume. Note that this does
!!! not represent passive recoil of the lung. It assumes a sinusoidal pressure
!!! acting external to the lung, approximately the opposite of the inspiration
!!! muscle pressure.

!!! T_interval !the total length of the breath

!!! Gdirn !direction of gravitational gradient for tissue volumes
!!! 'Gdirn' is 1(x), 2(y), 3(z); upright lung (for our models) is z, supine is y.

!!! press_in !constant value of pressure at the entry to the model

!!! COV !coefficient of variation of tissue compliance

!!! RMaxMean !ratio max to mean volume

!!! RMinMean !ratio min to mean volume

!!! i_to_e_ratio ! ratio of inspiration time to expiration time

!!! refvol !sets the 'zero stress' state at half of our lung volume (OK for FRC)

!!! volume_target ! the target tidal volume, in mm^3

!!! pmus_step !the change in Ppl that is driving the flow, in Pa

!!! Texpn ! time for expiration

!!! Tinsp ! time for inspiration

!!! undef !!!! undef is the volume at which our stress-strain relationship has zero
!!! stress. i.e. the tissue is 'undeformed'. Lung tissue is never stress-free,
!!! however our continuum model approach requires this to be defined (a theoretical
!!! state: should be less than our minimum simulation volume, e.g. RV).

!!! sum_tidal !tidal volume for a breath

!!! sum_expid ! sum of expired volume for a breath. Should equal the sum inspired.
!!! -----------------------------------------------------------------------------

!!! set default values for the parameters that control the breathing simulation
!!! these should be controlled by user input (showing hard-coded for now)

    CALL read_params_evaluate_flow(Gdirn, chest_wall_compliance, &
         constrict, COV, FRC, i_to_e_ratio, pmus_step, press_in,&
         refvol, RMaxMean, RMinMean, T_interval, volume_target, expiration_type)
    CALL read_params_main(num_brths, num_itns, dt, err_tol)

!!! calculate key variables from the boundary conditions/problem parameters
    Texpn = T_interval / (1.0_dp+i_to_e_ratio)
    Tinsp = T_interval - Texpn

!!! store initial branch lengths, radii, resistance etc. in array 'elem_field'
    CALL update_elem_field

    CALL volume_of_mesh(init_vol,volume_tree) ! to get deadspace volume

!!! distribute the initial tissue unit volumes along the gravitational axis.
    CALL set_initial_volume(Gdirn,COV,FRC*1.0e+6_dp,RMaxMean,RMinMean)
    undef = refvol * (FRC*1.0e+6_dp-volume_tree)/DBLE(elem_units_below(1))

!!! calculate the total model volume
    CALL volume_of_mesh(init_vol,volume_tree)

    WRITE(*,'('' Anatomical deadspace = '',F8.3,'' ml'')') volume_tree/1.0e+3_dp ! in mL
    WRITE(*,'('' Respiratory volume   = '',F8.3,'' L'')') (init_vol-volume_tree)/1.0e+6_dp !in L
    WRITE(*,'('' Total lung volume    = '',F8.3,'' L'')') init_vol/1.0e+6_dp !in L

    unit_field(nu_dpdt,1:num_units) = 0.0_dp

!!! calculate the compliance of each tissue unit
    CALL tissue_compliance(chest_wall_compliance,undef)
    totalc = SUM(unit_field(nu_comp,1:num_units)) !the total model compliance

    CALL update_pleural_pressure(ppl_current) !calculate new pleural pressure
    pptrans=SUM(unit_field(nu_pe,1:num_units))/num_units

    ChestWallRestVol = init_vol + 0.2e+6_dp/98.0665_dp * (-ppl_current)
    Pcw = (ChestWallRestVol - init_vol)/(0.2e+6_dp/98.0665_dp)

!!! write out the header information for run-time output
    WRITE(*,'(2X,''Time'',3X,''Inflow'',4X,''V_t'',5X,''Raw'',5X,&
         &''Comp'',4X,''Ppl'',5X,''Ptp'',5X,''VolL'',4X,''Pmus'',&
         &4X,''Pcw'',2X,''Pmus-Pcw'')')
    WRITE(*,'(3X,''(s)'',4X,''(mL/s)'',3X,''(mL)'',1X,''(cmH/L.s)'',&
         &1X,''(L/cmH)'',1X,''(...cmH2O...)'',&
         &4X,''(L)'',5X,''(......cmH2O.......)'')')

    WRITE(*,'(F7.3,2(F8.1),8(F8.2))') &
         0.0,0.0,0.0, &  !time, flow, tidal
         elem_field(ne_t_resist,1)*1.0e+6_dp/98.0665_dp, & !airway resistance (cmH2O/L.s)
         totalC*98.0665_dp/1.0e+6_dp, & !total model compliance
         ppl_current/98.0665_dp, & !Ppl (cmH2O)
         -ppl_current/98.0665_dp, & !mean Ptp (cmH2O)
         init_vol/1.0e+6_dp, & !total model volume (L)
         0.0, & !Pmuscle (cmH2O)
         Pcw/98.0665_dp, & !Pchest_wall (cmH2O)
         (-Pcw)/98.0665_dp !Pmuscle - Pchest_wall (cmH2O)

!!! Initialise variables:
    pmus_factor_in = 1.0_dp
    pmus_factor_ex = 1.0_dp
    time = 0.0_dp !initialise the simulation time.
    n = 0 !initialise the 'breath number'. This is incremented at the start of each breath.
    sum_tidal = 0.0_dp ! initialise the inspired and expired volumes
    sum_expid = 0.0_dp
    last_vol = init_vol

    CONTINUE=.TRUE.
    DO WHILE(CONTINUE)
       n=n+1 !increment the breath number
       ttime=0.0_dp !each breath starts with ttime=0
       endtime = T_interval * n - 0.5_dp * dt !the end time of the breath
       p_mus = 0.0_dp !initialise the muscle pressure to zero
       ptrans_frc=SUM(unit_field(nu_pe,1:num_units))/num_units !ptrans at frc
       IF(n.GT.1)THEN !write out 'end of breath' information
          WRITE(*,'('' End of breath, inspired = '',F10.2,'' L'')') sum_tidal/1.0e+6_dp
          WRITE(*,'('' End of breath, expired  = '',F10.2,'' L'')') sum_expid/1.0e+6_dp
          WRITE(*,'('' Peak muscle pressure    = '',F10.2,'' cmH2O'')') &
               pmus_step*pmus_factor_in/98.0665_dp
          WRITE(*,'('' Drift in FRC from start = '',F10.2,'' %'')') &
               100*(now_vol-init_vol)/init_vol
          WRITE(*,'('' Difference from target Vt = '',F8.2,'' %'')') &
               100*(volume_target-sum_tidal)/volume_target
          WRITE(*,'('' Total Work of Breathing ='',F7.3,''J/min'')')WOB_insp
          WRITE(*,'('' elastic WOB ='',F7.3,''J/min'')')WOBe_insp
          WRITE(*,'('' resistive WOB='',F7.3,''J/min'')')WOBr_insp

          IF(DABS(volume_target).GT.1.0e-5_dp)THEN
             ! modify driving muscle pressure by volume_target/sum_tidal
             ! this increases p_mus for volume_target>sum_tidal, and
             ! decreases p_mus for volume_target<sum_tidal
             pmus_factor_in=pmus_factor_in*DABS(volume_target/sum_tidal)
             pmus_factor_ex=pmus_factor_ex*DABS(volume_target/sum_expid)
          ENDIF
          sum_tidal=0.0_dp !reset the tidal volume
          sum_expid=0.0_dp !reset the expired volume
          unit_field(nu_vt,1:num_units) = 0.0_dp !reset acinar tidal volume
          sum_dpmus=0.0_dp
          sum_dpmus_ei=0.0_dp
       ENDIF

!!! Do the simulation of each breath
       DO WHILE (time.LT.endtime)
          !          time = time + dt
          ttime = ttime + dt
          ! set the increment in driving (muscle) pressure
          IF(expiration_type(1:6).EQ.'active')THEN
             IF(ttime.LT.Tinsp)THEN
                dpmus=pmus_step*pmus_factor_in*PI* &
                     SIN(2.0_dp*pi/(2.0_dp*Tinsp)*ttime)/(2.0_dp*Tinsp)*dt
             ELSEIF(ttime.LE.Tinsp+Texpn)THEN
                dpmus=pmus_step*pmus_factor_ex*PI* &
                     SIN(2.0_dp*pi*(0.5d0+(ttime-Tinsp)/(2.0_dp*Texpn)))/ &
                     (2.0_dp*Texpn)*dt
             ENDIF
          ELSEIF(expiration_type(1:7).EQ.'passive')THEN
             IF(ttime.LE.Tinsp+0.5d0*dt)THEN
                dpmus=pmus_step*pmus_factor_in*PI*dt* &
                     SIN(pi*ttime/Tinsp)/(2.0_dp*Tinsp)
                sum_dpmus=sum_dpmus+dpmus
                sum_dpmus_ei=sum_dpmus
             ELSE
                Tpass=0.1d0
                dpmus=MIN(-sum_dpmus_ei/(Tpass*Texpn)*dt,-sum_dpmus)
                sum_dpmus=sum_dpmus+dpmus
             ENDIF
          ENDIF

          p_mus = p_mus + dpmus !current value for muscle pressure
          prev_flow=elem_field(ne_Vdot,1)

!!! Solve for a new flow and pressure field
!!! We will estimate the flow into each terminal lumped
!!! parameter unit (assumed to be an acinus), so we can calculate flow
!!! throughout the rest of the tree simply by summation. After summing
!!! the flows we can use the resistance equation (P0-P1=R1*Q1) to update
!!! the pressures throughout the tree.

          !initialise Qinit to the previous flow
          elem_field(ne_Vdot0,1:num_elems) = elem_field(ne_Vdot,1:num_elems)
          converged = .FALSE.
          iter_step=0
          DO WHILE (.NOT.converged)
             iter_step=iter_step+1 !count the iterative steps
             CALL estimate_flow(dpmus,dt,err_est) !analytic solution for Q
             IF(iter_step.GT.1.AND.err_est.LT.err_tol)THEN
                converged=.TRUE.
             ELSE IF(iter_step.GT.num_itns)THEN
                converged=.TRUE.
                WRITE(*,'('' Warning: lower convergence '// &
                     'tolerance and time step - check values, Error='',D10.3)') &
                     err_est
             ENDIF
             CALL sum_elem_field_from_periphery(ne_Vdot) !sum the flows recursively UP the tree
             CALL update_elem_field !updates resistances
             CALL update_node_pressures(press_in) !updates the pressures at nodes
             CALL update_unit_dpdt(dt) ! update dP/dt at the terminal units
          ENDDO !converged

          CALL update_unit_volume(dt,Tinsp,Texpn) ! Update tissue unit volumes, and unit tidal volumes
          CALL volume_of_mesh(now_vol,volume_tree) !calculate the mesh volume, store in 'now_vol'
          CALL update_elem_field  !update element lengths, volumes, resistances
          CALL tissue_compliance(chest_wall_compliance,undef) !update the unit compliances, uses 'undef' as input
          totalc = SUM(unit_field(nu_comp,1:num_units)) !the total model compliance
          CALL update_pleural_pressure(ppl_current) !calculate new pleural pressure
          CALL update_proximal_pressure !updates values of pressure at proximal nodes of end branches
          CALL calculate_work(now_vol-init_vol,now_vol-last_vol,WOBe,WOBr,pptrans)!calculate work of breathing
          last_vol=now_vol
          Pcw = (now_vol - ChestWallRestVol)/(0.2d6/98.0665_dp)

          time = time + dt
          WRITE(*,'(F7.3,2(F8.1),8(F8.2))') &
               time, & !time through breath (s)
               elem_field(ne_Vdot,1)/1.0e+3_dp, & !flow at the inlet (mL/s)
               (now_vol - init_vol)/1.0e+3_dp, & !current tidal volume (mL)
               elem_field(ne_t_resist,1)*1.0e+6_dp/98.0665_dp, & !airway resistance (cmH2O/L.s)
               totalC*98.0665_dp/1.0e+6_dp, & !total model compliance
               ppl_current/98.0665_dp, & !Ppl (cmH2O)
               pptrans/98.0665_dp, & !mean Ptp (cmH2O)
               now_vol/1.0e+6_dp, & !total model volume (L)
               p_mus/98.0665_dp, & !Pmuscle (cmH2O)
               -Pcw/98.0665_dp, & !Pchest_wall (cmH2O)
               (p_mus+Pcw)/98.0665_dp !Pmuscle - Pchest_wall (cmH2O)

          ! increment the tidal volume, or the volume expired
          IF(elem_field(ne_Vdot,1).GT.0.0_dp)THEN
             sum_tidal=sum_tidal+elem_field(ne_Vdot,1)*dt
          ELSE
             sum_expid=sum_expid-elem_field(ne_Vdot,1)*dt
             IF(prev_flow.GT.0.0_dp)THEN
                WOBe_insp=(WOBe+sum_tidal*ptrans_frc*1.0e-9_dp)*(30.0_dp/Tinsp)
                WOBr_insp=WOBr*(30.0_dp/Tinsp)
                WOB_insp=WOBe_insp+WOBr_insp
                WOBe=0.0_dp
                WOBr=0.0_dp
             ENDIF
          ENDIF

       ENDDO !while ttime<endtime

       !...  CHECK WHETHER TO CONTINUE
       IF(n.GE.num_brths)THEN
          CONTINUE=.FALSE.
       ELSEIF(DABS(volume_target).GT.1.0e-3_dp)THEN
          IF(DABS(100.0_dp*(volume_target-sum_tidal) &
               /volume_target).GT.0.1_dp.OR.(n.LT.2))THEN
             CONTINUE=.TRUE.
          ELSE
             CONTINUE=.FALSE.
          ENDIF
       ENDIF

    ENDDO !...WHILE(CONTINUE)

    WRITE(*,'(''End of breath, inspired = '',F10.2,'' L'')') sum_tidal/1.0e+6_dp
    WRITE(*,'(''End of breath, expired  = '',F10.2,'' L'')') sum_expid/1.0e+6_dp
    WRITE(*,'(''Peak muscle pressure    = '',F10.2,'' cmH2O'')') &
         pmus_step*pmus_factor_in/98.0665_dp
    WRITE(*,'(''Drift in FRC from start = '',F10.2,'' %'')') &
         100.0_dp*(now_vol-init_vol)/init_vol
    WRITE(*,'(''Difference from target Vt = '',F8.2,'' %'')') &
         100*(volume_target-sum_tidal)/volume_target
    WRITE(*,'(''Total Work of Breathing ='',F7.3,''J/min'')')WOB_insp
    WRITE(*,'(''elastic WOB ='',F7.3,''J/min'')')WOBe_insp
    WRITE(*,'(''resistive WOB='',F7.3,''J/min'')')WOBr_insp

!!! Transfer the tidal volume for each elastic unit to the terminal branches, and sum up the tree.
!!! Divide by inlet flow. This gives the time-averaged and normalised flow field for the tree.
    DO nunit=1,num_units !for each terminal element only (with tissue units attached)
       ne=units(nunit) !local element number
       elem_field(ne_Vdot,ne) = unit_field(nu_vt,nunit)
    ENDDO
    CALL sum_elem_field_from_periphery(ne_Vdot)
    elem_field(ne_Vdot,1:num_elems) = elem_field(ne_Vdot,1:num_elems)/elem_field(ne_Vdot,1)

    !    call export_terminal_solution(TERMINAL_EXNODEFILE,'terminals')

    CALL enter_exit(sub_name,2)

  END SUBROUTINE evaluate_vent


  !###################################################################################

  SUBROUTINE evaluate_uniform_flow
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_UNIFORM_FLOW" :: EVALUATE_UNIFORM_FLOW

    USE arrays,ONLY: dp,elem_field,num_elems,num_units,&
         units,unit_field
    USE indices,ONLY: ne_Vdot,nu_Vdot0,nu_vol
    USE other_consts
    USE diagnostics, ONLY: enter_exit
    USE geometry,ONLY: volume_of_mesh
    IMPLICIT NONE

    ! Local variables
    INTEGER :: ne,nunit
    REAL(dp) :: init_vol,volume_tree

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'evaluate_uniform_flow'
    CALL enter_exit(sub_name,1)

!!! calculate the total model volume
    CALL volume_of_mesh(init_vol,volume_tree)

!!! initialise the flow field to zero
    elem_field(ne_Vdot,1:num_elems) = 0.0_dp

!!! For each elastic unit, calculate uniform ventilation
    DO nunit=1,num_units !for each terminal element only (with tissue units attached)
       ne=units(nunit) !local element number
       unit_field(nu_Vdot0,nunit) = unit_field(nu_vol,nunit)/(init_vol-volume_tree)
       elem_field(ne_Vdot,ne) = unit_field(nu_Vdot0,nunit)
    ENDDO

    CALL sum_elem_field_from_periphery(ne_Vdot)

    CALL enter_exit(sub_name,2)

  END SUBROUTINE evaluate_uniform_flow


!!!###################################################################################

  SUBROUTINE update_unit_dpdt(dt)
    USE arrays,ONLY: dp,elem_nodes,node_field,&
         num_units,units,unit_field
    USE indices,ONLY: nj_aw_press,nu_dpdt,nu_air_press
    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    REAL(dp), INTENT(in) :: dt
    INTEGER :: ne,np1,nunit
    REAL(dp) :: est

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_unit_dpdt'
    CALL enter_exit(sub_name,1)

    ! this is the rate of change of pressure at the proximal end of the element
    ! that supplies the tissue unit, i.e. not the rate of change of pressure within the unit.
    DO nunit=1,num_units
       ne=units(nunit)
       np1=elem_nodes(1,ne)
       ! linear estimate
       est=(node_field(nj_aw_press,np1) &
            -unit_field(nu_air_press,nunit))/dt
       ! weight new estimate with the previous dP/dt
       unit_field(nu_dpdt,nunit)=0.5d0*(est+unit_field(nu_dpdt,nunit))
    ENDDO !nunit

    CALL enter_exit(sub_name,2)

  END SUBROUTINE update_unit_dpdt


!!!###################################################################################

  SUBROUTINE update_proximal_pressure
    USE arrays,ONLY: elem_nodes,node_field,num_units,units,unit_field
    USE indices,ONLY: nj_aw_press,nu_air_press
    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER :: ne,np1,nunit

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_proximal_pressure'
    CALL enter_exit(sub_name,1)

    ! update the pressure at the proximal node of the element that feeds the elastic unit
    DO nunit=1,num_units
       ne=units(nunit)
       np1=elem_nodes(1,ne)
       unit_field(nu_air_press,nunit)=node_field(nj_aw_press,np1) !store the pressure at entry node
    ENDDO !noelem

    CALL enter_exit(sub_name,2)

  END SUBROUTINE update_proximal_pressure


!!!###################################################################################

  SUBROUTINE update_pleural_pressure(ppl_current)
!!! Update the mean pleural pressure based on current Pel (=Ptp) and Palv,
!!! i.e. Ppl(unit) = -Pel(unit)+Palv(unit)
    USE arrays,ONLY: dp,elem_nodes,node_field,num_units,units,unit_field
    USE diagnostics, ONLY: enter_exit
    USE indices,ONLY: nj_aw_press,nu_pe
    IMPLICIT NONE

    REAL(dp),INTENT(out) :: ppl_current
    INTEGER :: ne,np2,nunit

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_pleural_pressure'
    CALL enter_exit(sub_name,1)

    ppl_current = 0.0_dp
    DO nunit=1,num_units
       ne=units(nunit)
       np2=elem_nodes(2,ne)
       ppl_current = ppl_current - unit_field(nu_pe,nunit) + &
            node_field(nj_aw_press,np2)
    ENDDO !noelem
    ppl_current = ppl_current/num_units

    CALL enter_exit(sub_name,2)

  END SUBROUTINE update_pleural_pressure


!!!###################################################################################


  SUBROUTINE update_node_pressures(press_in)
    USE arrays,ONLY: dp,elem_field,elem_nodes,node_field,num_elems
    USE indices,ONLY: ne_Vdot,ne_resist,nj_aw_press
    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    REAL(dp),INTENT(in) :: press_in
    !Local parameters
    INTEGER :: ne,np1,np2

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_node_pressures'
    CALL enter_exit(sub_name,1)

!!! Use the known resistances and flows to calculate nodal pressures through whole tree

    ! set the initial node pressure to be the input pressure (usually zero)
    ne=1 !element number at top of tree, usually = 1
    np1=elem_nodes(1,ne) !first node in element
    node_field(nj_aw_press,np1)=press_in !set pressure at top of tree

    DO ne=1,num_elems !for each element
       np1=elem_nodes(1,ne) !start node number
       np2=elem_nodes(2,ne) !end node number
       !P(np2) = P(np1) - Resistance(ne)*Flow(ne)
       node_field(nj_aw_press,np2)=node_field(nj_aw_press,np1) &
            -elem_field(ne_resist,ne)*elem_field(ne_Vdot,ne)
    ENDDO !noelem

    CALL enter_exit(sub_name,2)

  END SUBROUTINE update_node_pressures


!!!###################################################################################

  SUBROUTINE tissue_compliance(chest_wall_compliance,undef)
    USE arrays,ONLY: dp,num_units,units,unit_field
    USE indices,ONLY: nu_comp,nu_pe,nu_vol
    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    !     Parameter List
    REAL(dp), INTENT(in) :: chest_wall_compliance,undef

    INTEGER :: ne,nunit
    REAL(dp),PARAMETER :: a = 0.433d0, b = -0.611d0, cc = 2500.0_dp
    REAL(dp) :: exp_term,lambda,ratio

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_tissue_compliance'
    CALL enter_exit(sub_name,1)

    !.....dV/dP=1/[(1/2h^2).c/2.(3a+b)exp().(4h(h^2-1)^2)+(h^2+1)/h^2)]

    DO nunit=1,num_units
       ne=units(nunit)
       !calculate a compliance for the tissue unit
       ratio=unit_field(nu_vol,nunit)/undef
       lambda = ratio**(1.0_dp/3.0_dp) !uniform extension ratio
       exp_term=DEXP(0.75d0*(3.0_dp*a+b)*(lambda**2-1.0_dp)**2)

       unit_field(nu_comp,nunit)=cc*exp_term/6.0_dp*(3.0_dp*(3.0_dp*a+b)**2 &
            *(lambda**2-1.0_dp)**2/lambda**2+(3.0_dp*a+b) &
            *(lambda**2+1.0_dp)/lambda**4)
       unit_field(nu_comp,nunit)=undef/unit_field(nu_comp,nunit) !in units of volume/pressure
       ! add the chest wall (proportionately) in parallel
       unit_field(nu_comp,nunit) = 1.0_dp/(1.0_dp/unit_field(nu_comp,nunit)&
            +1.0_dp/(chest_wall_compliance/DBLE(num_units)))

       !estimate an elastic recoil pressure for the unit
       unit_field(nu_pe,nunit) = cc/2.0_dp*(3.0_dp*a+b)*(lambda**2.0_dp &
            -1.0_dp)*exp_term/lambda
    ENDDO !nunit

    CALL enter_exit(sub_name,2)

  END SUBROUTINE tissue_compliance


!!!###################################################################################

  SUBROUTINE sum_elem_field_from_periphery(ne_field)
    USE arrays,ONLY: dp,elem_cnct,elem_field,elem_symmetry,num_elems
    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER,INTENT(in) :: ne_field

    !Local parameters
    REAL(dp) :: field_value
    INTEGER :: i,ne,ne2

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'sum_elem_field_from_periphery'
    CALL enter_exit(sub_name,1)

    DO ne=num_elems,1,-1
       IF(elem_cnct(1,0,ne).GT.0)THEN !not terminal
          field_value=0.d0
          DO i=1,elem_cnct(1,0,ne) !for each possible daughter branch (maximum 2)
             ne2=elem_cnct(1,i,ne) !the daughter element number
             field_value=field_value+elem_symmetry(ne2)*elem_field(ne_field,ne2) !sum daughter fields
          ENDDO !noelem2
          elem_field(ne_field,ne)=field_value
       ENDIF
    ENDDO !noelem

    CALL enter_exit(sub_name,2)

  END SUBROUTINE sum_elem_field_from_periphery

!!!###################################################################################

  SUBROUTINE update_unit_volume(dt,Tinsp, Texpn)
    USE arrays,ONLY: dp,elem_field,elem_nodes,num_units,&
         units,unit_field
    USE diagnostics, ONLY: enter_exit
    USE indices,ONLY: ne_Vdot,nu_vol,nu_vt,nu_vent
    IMPLICIT NONE

    REAL(dp),INTENT(in) :: dt,Tinsp,Texpn
    INTEGER :: ne,np,nunit

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_unit_volume'
    CALL enter_exit(sub_name,1)

    DO nunit=1,num_units
       ne=units(nunit)
       np=elem_nodes(2,ne)
       ! update the volume of the lumped tissue unit
       unit_field(nu_vol,nunit)=unit_field(nu_vol,nunit)+dt* &
            elem_field(ne_Vdot,ne) !in mm^3
       IF(elem_field(ne_Vdot,1).GT.0.0_dp)THEN  !only store inspired volume
          unit_field(nu_vt,nunit)=unit_field(nu_vt,nunit)+dt* &
               elem_field(ne_Vdot,ne)
          unit_field(nu_vent,nunit)=unit_field(nu_vt,nunit)/(TInsp+Texpn) ! volume  per secondto alvolus
       ENDIF
    ENDDO !nunit

    CALL enter_exit(sub_name,2)

  END SUBROUTINE update_unit_volume

!!!####################################################################

  SUBROUTINE update_elem_field
    USE arrays,ONLY: dp,elem_field,elem_nodes,node_xyz,num_elems
    USE diagnostics, ONLY: enter_exit
    USE indices,ONLY: ne_Vdot,ne_length, &
         ne_radius,ne_resist,ne_t_resist,ne_vol
    USE other_consts
    IMPLICIT NONE

    ! Local variables
    INTEGER :: ne,np1,np2
    REAL(dp),PARAMETER :: gas_density = 0.1146d-5 ! g.mm^-3
    REAL(dp),PARAMETER :: gas_viscosity = 0.18d-4 ! Pa.s
    REAL(dp) :: gamma,resistance,reynolds,zeta
    REAL(dp) :: rad,le

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_elem_volume'
    CALL enter_exit(sub_name,1)

    DO ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)

       ! element length
       elem_field(ne_length,ne) = DSQRT((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)

       ! element volume
       elem_field(ne_vol,ne) = PI * elem_field(ne_radius,ne)**2 * &
            elem_field(ne_length,ne)

       le=elem_field(ne_length,ne)
       rad=elem_field(ne_radius,ne)

       ! element Poiseuille (laminar) resistance in units of Pa.s.mm-3
       resistance = 8.0_dp*GAS_VISCOSITY*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_radius,ne)**4) !laminar resistance

       ! element turbulent resistance (flow in bifurcating tubes)
       gamma = 0.357_dp !inspiration
       IF(elem_field(ne_Vdot,ne).LT.0.0_dp) gamma = 0.46_dp !expiration

       reynolds=DABS(elem_field(ne_Vdot,ne)*2.0_dp*GAS_DENSITY/ &
            (PI*elem_field(ne_radius,ne)*GAS_VISCOSITY))
       zeta = MAX(1.0_dp,dsqrt(2.0_dp*elem_field(ne_radius,ne)* &
            reynolds/elem_field(ne_length,ne))*gamma)
       elem_field(ne_resist,ne) = resistance * zeta

       elem_field(ne_t_resist,ne) = elem_field(ne_resist,ne)

    ENDDO !noelem

    CALL enter_exit(sub_name,2)

  END SUBROUTINE update_elem_field

!!!####################################################################

  SUBROUTINE estimate_flow(dpmus,dt,err_est)
    USE arrays,ONLY: dp,elem_field,num_units,&
         units,unit_field
    USE diagnostics, ONLY: enter_exit
    USE indices,ONLY: ne_Vdot,ne_Vdot0,ne_resist,nu_comp,&
         nu_dpdt,nu_Vdot0,nu_Vdot1,nu_Vdot2
    IMPLICIT NONE

    REAL(dp),INTENT(in) :: dpmus,dt
    REAL(dp),INTENT(out) :: err_est

    INTEGER :: ne,nunit
    REAL(dp) :: alpha,beta,flow_diff,flow_sum,Q,Qinit

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'estimate_flow'
    CALL enter_exit(sub_name,1)

    err_est = 0.0_dp
    flow_sum = 0.0_dp
!!! For each elastic unit, calculate Qbar (equation 4.13 from Swan thesis)
    DO nunit=1,num_units !for each terminal element only (with tissue units attached)
       ne=units(nunit) !local element number
       ! Calculate the mean flow into the unit in the time step
       ! alpha is the rate of change of pressure at start node of terminal element
       alpha=unit_field(nu_dpdt,nunit) !dPaw/dt, initialised to zero and updated each iter
       Qinit=elem_field(ne_Vdot0,ne) !terminal element flow, updated each dt
       ! beta is the rate of change of pleural pressure
       beta = dpmus/dt ! == dPmus/dt (-ve for inspiration), updated each dt

       !      Q = C*(alpha-beta)+(Qinit-C*(alpha-beta))*DEXP(-dt/(C*R))
       Q = unit_field(nu_comp,nunit)*(alpha-beta)+ &
            (Qinit-unit_field(nu_comp,nunit)*(alpha-beta))* &
            DEXP(-dt/(unit_field(nu_comp,nunit)*elem_field(ne_resist,ne)))

       unit_field(nu_Vdot2,nunit)=unit_field(nu_Vdot1,nunit) !flow at iter-2
       unit_field(nu_Vdot1,nunit)=unit_field(nu_Vdot0,nunit) !flow at iter-1

       ! flow estimate for current iter includes flow estimates at previous two iters
       !       unit_field(nu_Vdot0,nunit) = 0.75d0*unit_field(nu_Vdot2,nunit)+ &
       !            0.25d0*(Q+unit_field(nu_Vdot1,nunit))*0.5d0
       ! flow estimate for current iter includes flow estimate at previous iter
       !       unit_field(nu_Vdot0,nunit) = (Q + unit_field(nu_Vdot1,nunit))*0.5d0
!!! from original code:
       unit_field(nu_Vdot0,nunit) = 0.75d0*unit_field(nu_Vdot2,nunit)+ &
            0.25d0*(Q+unit_field(nu_Vdot1,nunit))*0.5d0

       flow_diff=unit_field(nu_Vdot0,nunit) - elem_field(ne_Vdot,ne)
       err_est=err_est+DABS(flow_diff)**2 !sum up the error for all elements
       flow_sum=flow_sum+unit_field(nu_Vdot0,nunit)**2

!!! ARC: DO NOT CHANGE BELOW. THIS IS NEEDED FOR THE ITERATIVE STEP
!!! - SIMPLER OPTIONS JUST FORCE IT TO CONVERGE WHEN ITS NOT
       elem_field(ne_Vdot,ne) = (unit_field(nu_Vdot0,nunit)&
            +unit_field(nu_Vdot1,nunit))/2.0_dp
       unit_field(nu_Vdot0,nunit) = elem_field(ne_Vdot,ne)
    ENDDO !nunit

    ! the estimate of error for the iterative solution
    err_est=err_est/(flow_sum*DBLE(num_units))

    CALL enter_exit(sub_name,2)

  END SUBROUTINE estimate_flow

!!!#############################################################################

  SUBROUTINE calculate_work(breath_vol,dt_vol,WOBe,WOBr,pptrans)
    USE arrays,ONLY: dp,elem_nodes,node_field,&
         num_units,units,unit_field
    USE diagnostics, ONLY: enter_exit
    USE indices,ONLY: nj_aw_press,nu_pe
    IMPLICIT NONE
    REAL(dp) :: breath_vol,dt_vol,WOBe,WOBr,pptrans
    ! Local variables
    INTEGER :: ne,np1,nunit
    REAL(dp) :: p_resis,p_trans

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'calculate_work'
    CALL enter_exit(sub_name,1)

    p_resis=0.0d0
    !estimate elastic and resistive WOB for each dt (sum dP.V)
    p_trans=SUM(unit_field(nu_pe,1:num_units))/num_units
    DO nunit=1,num_units
       ne=units(nunit)
       np1=elem_nodes(2,ne)
       p_resis=p_resis+node_field(nj_aw_press,1)-node_field(nj_aw_press,np1)
    ENDDO
    p_resis=p_resis/num_units
    ! vol in mm3 *1e-9=m3, pressure in Pa, hence *1d-9 = P.m3 (Joules)
    WOBe=WOBe+(p_trans-pptrans)*breath_vol*1d-9
    WOBr=WOBr+p_resis*dt_vol*1d-9

    pptrans=p_trans

    CALL enter_exit(sub_name,2)

  END SUBROUTINE calculate_work

!!!###########################################################################

  SUBROUTINE read_params_main(num_brths, num_itns, dt, err_tol)
    USE arrays,ONLY: dp
    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER,INTENT(out) :: num_brths, num_itns
    REAL(dp) :: dt,err_tol

    ! Input related variables
    CHARACTER(len=100) :: buffer, label
    INTEGER :: pos
    INTEGER, PARAMETER :: fh = 15
    INTEGER :: ios
    INTEGER :: line

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'read_params_main'
    CALL enter_exit(sub_name,1)

    ios = 0
    line = 0
    OPEN(fh, file='Parameters/params_main.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    DO WHILE (ios == 0)
       READ(fh, '(A)', iostat=ios) buffer
       IF (ios == 0) THEN
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = SCAN(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          SELECT CASE (label)
          CASE ('num_brths')
             READ(buffer, *, iostat=ios) num_brths
             PRINT *, 'Read num_brths: ', num_brths
          CASE ('num_itns')
             READ(buffer, *, iostat=ios) num_itns
             PRINT *, 'Read num_itns: ', num_itns
          CASE ('dt')
             READ(buffer, *, iostat=ios) dt
             PRINT *, 'Read dt: ', dt
          CASE ('err_tol')
             READ(buffer, *, iostat=ios) err_tol
             PRINT *, 'Read err_tol: ', err_tol
          CASE default
             PRINT *, 'Skipping invalid label at line', line
          END SELECT
       END IF
    END DO

    CLOSE(fh)
    CALL enter_exit(sub_name,2)

  END SUBROUTINE read_params_main

!!!###########################################################################
  SUBROUTINE read_params_evaluate_flow (Gdirn, chest_wall_compliance, &
       constrict, COV, FRC, i_to_e_ratio, pmus_step, press_in,&
       refvol, RMaxMean, RMinMean, T_interval, volume_target, expiration_type)

    USE arrays,ONLY: dp
    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER,INTENT(out) :: Gdirn
    REAL(dp),INTENT(out) :: chest_wall_compliance, constrict, COV,&
         FRC, i_to_e_ratio, pmus_step, press_in,&
         refvol, RMaxMean, RMinMean, T_interval, volume_target
    CHARACTER,INTENT(out) :: expiration_type*(*)

    ! Input related variables
    CHARACTER(len=100) :: buffer, label
    INTEGER :: pos
    INTEGER, PARAMETER :: fh = 15
    INTEGER :: ios
    INTEGER :: line

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    ios = 0
    line = 0
    sub_name = 'read_params_evaluate_flow'
    CALL enter_exit(sub_name,1)

    ! following values are examples from control.txt
    !    T_interval = 4.0_dp !s
    !    Gdirn = 3
    !    press_in = 0.0_dp !Pa
    !    COV = 0.2d0
    !    RMaxMean = 1.29d0
    !    RMinMean = 0.78d0
    !    i_to_e_ratio = 0.5d0 !dimensionless
    !    refvol = 0.6d0 !dimensionless
    !    volume_target = 8.d5 !mm^3  800 ml
    !    pmus_step = -5.4d0 * 98.0665_dp !-5.4 cmH2O converted to Pa
    !    expiration_type = 'passive' ! or 'active'
    !    chest_wall_compliance = 0.2d6/98.0665_dp !(0.2 L/cmH2O --> mm^3/Pa)

    OPEN(fh, file='Parameters/params_evaluate_flow.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    DO WHILE (ios == 0)
       READ(fh, '(A)', iostat=ios) buffer
       IF (ios == 0) THEN
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = SCAN(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          SELECT CASE (label)
          CASE ('FRC')
             READ(buffer, *, iostat=ios) FRC
             PRINT *, 'Read FRC: ', FRC
          CASE ('constrict')
             READ(buffer, *, iostat=ios) constrict
             PRINT *, 'Read constrict: ', constrict
          CASE ('T_interval')
             READ(buffer, *, iostat=ios) T_interval
             PRINT *, 'Read T_interval: ', T_interval
          CASE ('Gdirn')
             READ(buffer, *, iostat=ios) Gdirn
             PRINT *, 'Read Gdirn: ', Gdirn
          CASE ('press_in')
             READ(buffer, *, iostat=ios) press_in
             PRINT *, 'Read press_in: ', press_in
          CASE ('COV')
             READ(buffer, *, iostat=ios) COV
             PRINT *, 'Read COV: ', COV
          CASE ('RMaxMean')
             READ(buffer, *, iostat=ios) RMaxMean
             PRINT *, 'Read RMaxMean: ', RMaxMean
          CASE ('RMinMean')
             READ(buffer, *, iostat=ios) RMinMean
             PRINT *, 'Read RMinMean: ', RMinMean
          CASE ('i_to_e_ratio')
             READ(buffer, *, iostat=ios) i_to_e_ratio
             PRINT *, 'Read i_to_e_ratio: ', i_to_e_ratio
          CASE ('refvol')
             READ(buffer, *, iostat=ios) refvol
             PRINT *, 'Read refvol: ', refvol
          CASE ('volume_target')
             READ(buffer, *, iostat=ios) volume_target
             PRINT *, 'Read volume_target: ', volume_target
          CASE ('pmus_step')
             READ(buffer, *, iostat=ios) pmus_step
             PRINT *, 'Read pmus_step_coeff: ', pmus_step
          CASE ('expiration_type')
             READ(buffer, *, iostat=ios) expiration_type
             PRINT *, 'Read expiration_type: ', expiration_type
          CASE ('chest_wall_compliance')
             READ(buffer, *, iostat=ios) chest_wall_compliance
             PRINT *, 'Read chest_wall_compliance: ', chest_wall_compliance
          CASE default
             PRINT *, 'Skipping invalid label at line', line
          END SELECT
       END IF
    END DO

    CLOSE(fh)
    CALL enter_exit(sub_name,2)

  END SUBROUTINE read_params_evaluate_flow

  !###################################################################################

  SUBROUTINE two_unit_test
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_TWO_UNIT_TEST" :: TWO_UNIT_TEST
    USE arrays
    USE diagnostics, ONLY: enter_exit
    USE geometry,ONLY: append_units
    USE indices
    IMPLICIT NONE

    INTEGER ne,noelem,nonode,np

    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'two_unit_test'
    CALL enter_exit(sub_name,1)

    ! set up a test geometry. this has only three branches and two elastic units

    num_nodes=4 !four nodes (branch junctions)
    num_elems=3 !three elements (branches)
    num_units = 2

    ALLOCATE (nodes(num_nodes))
    ALLOCATE (node_xyz(3,num_nodes))
    ALLOCATE (node_field(num_nj,num_nodes))
    ALLOCATE(elems(num_elems))
    ALLOCATE(elem_cnct(-1:1,0:2,0:num_elems))
    ALLOCATE(elem_nodes(2,num_elems))
    ALLOCATE(elem_ordrs(num_ord,num_elems))
    ALLOCATE(elem_symmetry(num_elems))
    ALLOCATE(elem_units_below(num_elems))
    ALLOCATE(elems_at_node(num_nodes,0:3))
    ALLOCATE(elem_field(num_ne,num_elems))
    ALLOCATE(elem_direction(3,num_elems))
    ALLOCATE(units(num_units))
    ALLOCATE(unit_field(num_nu,num_units))

    nodes=0
    node_xyz = 0.0_dp !initialise all values to 0
    node_field = 0.0_dp
    elem_field = 0.0_dp
    unit_field=0.0_dp
    elems=0
    units=0
    elem_cnct = 0

    DO nonode=1,num_nodes !loop over all of the nodes
       np=nonode
       nodes(nonode)=np !set the local node number to be the same as order in node list
    ENDDO !nonode

    node_xyz(3,2) = -100.0_dp !setting the z coordinate of node 2
    node_xyz(2,3) = -50.0_dp !setting the y coordinate of node 3
    node_xyz(2,4) = 50.0_dp !setting the y coordinate of node 4
    node_xyz(3,3) = -150.0_dp !setting the z coordinate of node 3
    node_xyz(3,4) = -150.0_dp !setting the z coordinate of node 4

    elem_field(ne_radius,1) = 8.0_dp
    elem_field(ne_radius,2) = 6.0_dp
    elem_field(ne_radius,3) = 5.0_dp

    ! set up elems:
    DO noelem=1,num_elems !loop over all of the elements
       ne=noelem
       elems(noelem)=ne !set the local element number to be the same as order in element list
       elem_nodes(2,noelem)=ne+1
    ENDDO !noelem
    elem_nodes(1,1)=1
    elem_nodes(1,2)=2
    elem_nodes(1,3)=2

    elem_cnct(-1,0,1:num_elems) = 1 !initialise all branches to have 1 parent
    elem_cnct(-1,0,1) = 0 !element 1 has 0 adjacent branches in the -Xi1 direction
    elem_cnct(1,0,1)=2 ! element 1 has 2 adjacent branches in the +Xi1 direction
    elem_cnct(1,1,1)=2 !element number of 1st adjacent branch
    elem_cnct(1,2,1)=3 !element number of 2nd adjacent branch
    elem_cnct(-1,1,2)=1 !element number of parent branch
    elem_cnct(-1,1,3)=1 !element number of parent branch

    elem_ordrs(no_gen,1)=1 !branch generation
    elem_ordrs(no_gen,2)=2 !branch generation
    elem_ordrs(no_gen,3)=3 !branch generation
    elem_ordrs(no_Hord,1)=2 !branch Horsfield order
    elem_ordrs(no_Hord,2)=1 !branch Horsfield order
    elem_ordrs(no_Hord,3)=1 !branch Horsfield order
    elem_ordrs(no_Sord,1)=2 !branch Strahler order
    elem_ordrs(no_Sord,2)=1 !branch Strahler order
    elem_ordrs(no_Sord,3)=1 !branch Strahler order

    CALL append_units

    unit_field(nu_vol,1) = 1.5d6 !arbitrary volume for element 2 (was BBM(1,ne)
    unit_field(nu_vol,2) = 1.5d6 !arbitrary volume for element 3

    elem_units_below(1)=2
    elem_units_below(2)=1
    elem_units_below(3)=1

    elem_symmetry(1:num_elems) = 1

    CALL enter_exit(sub_name,2)

  END SUBROUTINE two_unit_test

  !###################################################################################
END MODULE ventilation
