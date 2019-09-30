MODULE wave_transmission

  !*Brief Description:* Simulating wave propagation in a 1D tree structure
  !
  !*LICENSE:*
  !
  !
  !
  !*Full Description:*
  !Simulating wave propagation in a 1D tree structure
  !
  USE arrays, ONLY: dp
  USE other_consts, ONLY: PI
  IMPLICIT NONE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  PRIVATE
  PUBLIC evaluate_wave_transmission


CONTAINS
  !
  !##############################################################################
  !
  SUBROUTINE evaluate_wave_transmission(grav_dirn,grav_factor,&
       n_time,heartrate,a0,no_freq,a,b,n_adparams,admittance_param,n_model,model_definition,cap_model)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_WAVE_TRANSMISSION: EVALUATE_WAVE_PROPAGATION
    USE indices
    USE arrays, ONLY: dp,all_admit_param,num_elems,elem_field,fluid_properties,elasticity_param,num_units,&
         units,node_xyz,elem_cnct,elem_nodes,node_field
    USE diagnostics, ONLY: enter_exit



    INTEGER, INTENT(in) :: n_time
    REAL(dp), INTENT(in) :: heartrate
    REAL(dp), INTENT(in) :: a0
    INTEGER, INTENT(in):: no_freq
    REAL(dp), INTENT(in) :: a(no_freq)
    REAL(dp), INTENT(in) :: b(no_freq)
    INTEGER, INTENT(in) :: n_adparams
    REAL(dp), INTENT(in) :: admittance_param(n_adparams)
    INTEGER, INTENT(in) :: n_model
    REAL(dp), INTENT(in) :: model_definition(n_model)
    INTEGER, INTENT(in) :: grav_dirn
    INTEGER, INTENT(in) :: cap_model

    TYPE(all_admit_param) :: admit_param
    TYPE(fluid_properties) :: fluid
    TYPE(elasticity_param) :: elast_param

    CHARACTER(len=60) :: mesh_type
    REAL(dp) :: viscosity
    REAL(dp) :: density
    REAL(dp) :: harmonic_scale
    REAL(dp) :: steady_flow
    COMPLEX(dp), ALLOCATABLE :: eff_admit(:,:)
    COMPLEX(dp), ALLOCATABLE :: char_admit(:,:)
    COMPLEX(dp), ALLOCATABLE :: reflect(:,:)
    COMPLEX(dp), ALLOCATABLE :: prop_const(:,:)
    COMPLEX(dp), ALLOCATABLE :: p_factor(:,:)
    REAL(dp), ALLOCATABLE :: forward_pressure(:)
    REAL(dp), ALLOCATABLE :: reflected_pressure(:)
    REAL(dp), ALLOCATABLE :: forward_flow(:)
    REAL(dp), ALLOCATABLE :: reflected_flow(:)
    INTEGER :: min_art,max_art,min_ven,max_ven,min_cap,max_cap,ne,nu,nt,nf,np
    CHARACTER(len=30) :: tree_direction,mechanics_type
    REAL(dp) start_time,end_time,dt,time,omega
    REAL(dp) grav_vect(3),grav_factor,mechanics_parameters(2)
    INTEGER :: AllocateStatus,fid=10,fid2=20,fid3=30,fid4=40,fid5=50
    CHARACTER(len=60) :: sub_name
    LOGICAL :: constant_sheet_height

    sub_name = 'evalulate_wave_transmission'
    CALL enter_exit(sub_name,1)
    !!MODEL TYPE AND FLUID PROPERTIES
    !mesh_type: can be simple_tree, full_plus_ladder, full_sheet, full_tube The first can be airways, arteries, veins but no special features at the terminal level, the last one has arteries and veins connected by capillary units of some type (lung ladder acinus, lung sheet capillary bed, capillaries are just tubes represented by an element)
    IF(model_definition(1).EQ.1.0_dp)THEN
       mesh_type='simple_tree'
    ELSEIF(model_definition(1).EQ.2.0_dp)THEN
       mesh_type='full_plus_ladder'
       !elseif(model_definition(1).eq.3.0_dp)then
       !  full_sheet
       !elseif(model_definition(1).eq.4.0_dp)then
       ! full_tube
    ELSE
       PRINT *, 'ERROR: Your geometry choice has not yet been implemented'
       CALL EXIT(0)
    ENDIF
    !viscosity and density of fluid
    IF(model_definition(2).EQ.1.0_dp)THEN !BLOOD
       viscosity=fluid%blood_viscosity
       density=fluid%blood_density
    ELSEIF(model_definition(2).EQ.2.0_dp)THEN !AIR
       viscosity=fluid%air_viscosity
       density=fluid%air_density
    ELSE
       viscosity=model_definition(3)
       density=model_definition(4)
    ENDIF

    !!SET UP ADMITTANCE MODEL
    IF(admittance_param(1).EQ.1.0_dp)THEN
       admit_param%admittance_type='lachase_standard'
       elast_param%vessel_type='elastic_hooke'
       elast_param%elasticity_parameters(1)=admittance_param(2)!Pa
       elast_param%elasticity_parameters(2)=admittance_param(3)!Unitless
       elast_param%elasticity_parameters(3)=admittance_param(4)!dummy
    ELSEIF(admittance_param(1).EQ.2.0_dp)THEN
       admit_param%admittance_type='lachase_modified'
       elast_param%vessel_type='elastic_hooke'
       elast_param%elasticity_parameters(1)=admittance_param(2)!Pa
       elast_param%elasticity_parameters(2)=admittance_param(3)!Unitless
       elast_param%elasticity_parameters(3)=admittance_param(4)!dummy
    ELSEIF(admittance_param(1).EQ.3.0_dp)THEN
       admit_param%admittance_type='zhu_chesler'
       elast_param%vessel_type='elastic_hooke'
       elast_param%elasticity_parameters(1)=admittance_param(2)!Pa
       elast_param%elasticity_parameters(2)=admittance_param(3)!Unitless
       elast_param%elasticity_parameters(3)=admittance_param(4)!dummy
    ELSEIF(admittance_param(1).EQ.4.0_dp)THEN
       admit_param%admittance_type='duan_zamir'
       elast_param%vessel_type='elastic_alpha'
       elast_param%elasticity_parameters(1)=admittance_param(2)!/Pa
       elast_param%elasticity_parameters(2)=admittance_param(3)!Unitless
       elast_param%elasticity_parameters(3)=admittance_param(4)!dummy
    ELSE
       PRINT *, 'ERROR: Your admittance model choice has not yet been implemented'
       CALL EXIT(0)
    ENDIF

    !! SET UP PARAMETERS DEFINING OUTLET BOUNDARY CONDITIONS
    IF(admittance_param(5).EQ.1.0_dp)THEN !note we need to check that the right number of parameters have been input
       admit_param%bc_type='two_unit_wk'
       admit_param%two_parameter%admit_P1=admittance_param(6)
       admit_param%two_parameter%admit_P2=admittance_param(7)
    ELSEIF(admittance_param(5).EQ.2.0_dp)THEN
       admit_param%bc_type='three_unit_wk'
       admit_param%three_parameter%admit_P1=admittance_param(6)
       admit_param%three_parameter%admit_P2=admittance_param(7)
       admit_param%three_parameter%admit_P3=admittance_param(8)
    ELSEIF(admittance_param(5).EQ.4.0_dp)THEN
       admit_param%bc_type='two_wk_plus'
       admit_param%four_parameter%admit_P1=admittance_param(6)
       admit_param%four_parameter%admit_P2=admittance_param(7)
       admit_param%four_parameter%admit_P3=admittance_param(8)
       admit_param%four_parameter%admit_P4=admittance_param(9)
    ELSEIF(admittance_param(5).EQ.5.0_dp)THEN
       admit_param%bc_type='zero_reflection'
    ELSE
       PRINT *, 'ERROR: Your boundary condition choice has not yet been implemented'
       CALL EXIT(0)
    ENDIF

    !! Determining the capillary model (Constant sheet height / Non-constant sheet height)
    IF (cap_model.EQ.1) THEN
       constant_sheet_height = .TRUE.
    ELSE
       constant_sheet_height = .FALSE.
    ENDIF

    mechanics_type='linear'
    IF (mechanics_type.EQ.'linear') THEN
       mechanics_parameters(1)=5.0_dp*98.07_dp !average pleural pressure (Pa)
       mechanics_parameters(2)=0.25_dp*0.1e-2_dp !pleural density, defines gradient in pleural pressure
    ELSE
       PRINT *, 'ERROR: Only linear mechanics models have been implemented to date,assuming default parameters'
       CALL EXIT(0)
    ENDIF

    grav_vect=0.d0
    IF (grav_dirn.EQ.1) THEN
       grav_vect(1)=1.0_dp
    ELSEIF (grav_dirn.EQ.2) THEN
       grav_vect(2)=1.0_dp
    ELSEIF (grav_dirn.EQ.3) THEN
       grav_vect(3)=1.0_dp
    ELSE
       PRINT *, "ERROR: Posture not recognised (currently only x=1,y=2,z=3))"
       CALL EXIT(0)
    ENDIF
    grav_vect=grav_vect*grav_factor

    !!Determine steady component of flow
    IF(a0.EQ.0.0_dp)THEN !Using steady flow solution at inlet as a0
       steady_flow=elem_field(ne_Qdot,1)!ASSUMING FIRST ELEMENT
    ELSE!otherwise input a0 is used
       steady_flow=a0
    ENDIF
    !! SET UP PARAMETERS DEFINING COMPLIANCE MODEL
    harmonic_scale=heartrate/60.0_dp !frequency of first harmonic (Hz)
    !!ALLOCATE MEMORY
    ALLOCATE (eff_admit(1:no_freq,num_elems), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for eff_admit array ***"
    ALLOCATE (char_admit(1:no_freq,num_elems), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for char_admit array ***"
    ALLOCATE (reflect(1:no_freq,num_elems), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for reflect array ***"
    ALLOCATE (prop_const(1:no_freq,num_elems), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for prop_const array ***"
    ALLOCATE (p_factor(1:no_freq,num_elems), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for p_factor array ***"
    ALLOCATE (forward_pressure(n_time), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for forward_p array ***"
    ALLOCATE (reflected_pressure(n_time), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for reflected_p array ***"
    ALLOCATE (forward_flow(n_time), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for forward_q array ***"
    ALLOCATE (reflected_flow(n_time), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for reflected_q array ***"

    !initialise admittance
    char_admit=0.0_dp
    eff_admit=0.0_dp
    !calculate characteristic admittance of each branch
    CALL characteristic_admittance(no_freq,char_admit,prop_const,harmonic_scale, &
         density,viscosity,admit_param,elast_param,mechanics_parameters,grav_vect)

    !Apply boundary conditions to terminal units
    CALL boundary_admittance(no_freq,eff_admit,char_admit,admit_param,harmonic_scale,&
         density,viscosity,elast_param,mesh_type)


    ! calculate effective admittance through the tree
    IF(mesh_type.EQ.'full_plus_ladder')THEN
       min_art=1
       ne=1
       DO WHILE(elem_field(ne_group,ne).EQ.0.0_dp)
          max_art=ne
          ne=ne+1
       ENDDO
       min_ven=ne
       DO WHILE(elem_field(ne_group,ne).EQ.2.0_dp)
          max_ven=ne
          ne=ne+1
       ENDDO
       min_cap=ne
       max_cap=num_elems

       !vein admittance
       tree_direction='converging'
       CALL tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
            min_ven,max_ven,tree_direction)
       !        !cap admittance
       CALL capillary_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
            min_cap,max_cap,elast_param,mechanics_parameters,grav_vect,cap_model)!
       !art admittance
       tree_direction='diverging'
       CALL tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
            min_art,max_art,tree_direction)
    ELSE!Assume simple tree
       tree_direction='diverging'
       min_art=1
       max_art=num_elems
       CALL tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
            min_art,max_art,tree_direction)
    ENDIF

    !
    !    !calculate pressure drop through arterial tree (note to do veins too need to implement this concept thro' whole ladder model)
    !    !Also need to implement in reverse for veins
    CALL pressure_factor(no_freq,p_factor,reflect,prop_const,harmonic_scale,min_art,max_art)
    OPEN(fid5, file = 'inputimpedance.txt',action='write')
    WRITE(fid5,fmt=*) 'input impedance:'
    DO nf=1,no_freq
       omega=nf*harmonic_scale
       WRITE(fid5,fmt=*) omega,ABS(eff_admit(nf,1)),&
            ATAN2(imagpart(eff_admit(nf,1)),realpart(eff_admit(nf,1)))
    ENDDO
    CLOSE(fid5)

    start_time=0.0_dp
    end_time=60.0_dp/heartrate
    dt=(end_time-start_time)/n_time
    time=start_time
    !consider first pressure and flow into the vessel (at x=0)
    OPEN(fid, file = 'incident_pressure.txt', action='write')
    OPEN(fid2, file = 'incident_flow.txt',action='write')
    OPEN(fid3, file = 'total_pressure.txt',action='write')
    OPEN(fid4, file = 'total_flow.txt',action='write')
    DO nu =1,num_units
       ne=units(nu)
       forward_pressure=0.0_dp
       reflected_pressure=0.0_dp
       forward_flow=0.0_dp
       reflected_flow=0.0_dp
       DO nt=1,n_time
          DO nf=1,no_freq
             omega=2*pi*nf*harmonic_scale
             forward_pressure(nt)=forward_pressure(nt)+ABS(p_factor(nf,ne))*a(nf)*COS(omega*time+b(nf)+&
                  ATAN2(imagpart(p_factor(nf,ne)),realpart(p_factor(nf,ne))))

             reflected_pressure(nt)=reflected_pressure(nt)+ABS(p_factor(nf,ne))*a(nf)*&
                  ABS(reflect(nf,ne))*EXP((-2*elem_field(ne_length,ne))*realpart(prop_const(nf,ne)))*&
                  COS(omega*time+b(nf)+&
                  ATAN2(imagpart(p_factor(nf,ne)),realpart(p_factor(nf,ne)))+&
                  (-2*elem_field(ne_length,ne))*imagpart(prop_const(nf,ne))+&
                  ATAN2(imagpart(reflect(nf,ne)),realpart(reflect(nf,ne))))

             forward_flow(nt)=forward_flow(nt)+ABS(char_admit(nf,ne))*ABS(p_factor(nf,ne))*a(nf)*&
                  COS(omega*time+b(nf)+&
                  ATAN2(imagpart(p_factor(nf,ne)),realpart(p_factor(nf,ne)))+&
                  ATAN2(imagpart(char_admit(nf,ne)),realpart(char_admit(nf,ne))))

             reflected_flow(nt)=reflected_flow(nt)+ABS(char_admit(nf,ne))*ABS(p_factor(nf,ne))*a(nf)*&
                  ABS(reflect(nf,ne))*EXP((-2*elem_field(ne_length,ne))*realpart(prop_const(nf,ne)))*&
                  COS(omega*time+b(nf)+&
                  ATAN2(imagpart(p_factor(nf,ne)),realpart(p_factor(nf,ne)))+&
                  (-2*elem_field(ne_length,ne))*imagpart(prop_const(nf,ne))+&
                  ATAN2(imagpart(reflect(nf,ne)),realpart(reflect(nf,ne)))+&
                  ATAN2(imagpart(char_admit(nf,ne)),realpart(char_admit(nf,ne))))

          ENDDO
          time=time+dt
       ENDDO
       np=elem_nodes(2,ne)
       WRITE(fid,fmt=*) ne, forward_pressure+node_field(nj_bv_press,np)
       WRITE(fid2,fmt=*) ne, forward_flow+elem_field(ne_Qdot,ne)

       WRITE(fid3,fmt=*) ne, forward_pressure+reflected_pressure + node_field(nj_bv_press,np)
       WRITE(fid4,fmt=*) ne, forward_flow-reflected_flow + elem_field(ne_Qdot,ne)


    ENDDO
    CLOSE(fid)
    CLOSE(fid2)
    CLOSE(fid3)
    CLOSE(fid4)


    !!DEALLOCATE MEMORY
    DEALLOCATE (eff_admit, STAT = AllocateStatus)
    DEALLOCATE (char_admit, STAT = AllocateStatus)
    DEALLOCATE (reflect, STAT = AllocateStatus)
    DEALLOCATE (prop_const, STAT=AllocateStatus)
    DEALLOCATE (p_factor, STAT=AllocateStatus)
    DEALLOCATE (forward_pressure, STAT=AllocateStatus)
    DEALLOCATE (reflected_pressure, STAT=AllocateStatus)
    DEALLOCATE (forward_flow, STAT=AllocateStatus)
    DEALLOCATE (reflected_flow, STAT=AllocateStatus)
    CALL enter_exit(sub_name,2)
  END SUBROUTINE evaluate_wave_transmission
  !
  !##############################################################################
  !
  !*boundary_admittance* applies chosen admittance boundary conditions at the terminal units
  SUBROUTINE boundary_admittance(no_freq,eff_admit,char_admit,admit_param,harmonic_scale,&
       density,viscosity,elast_param,mesh_type)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_boundary_admittance: boundary_admittance
    USE arrays,ONLY: num_elems,all_admit_param,units,num_units,elasticity_param,elem_cnct
    USE diagnostics, ONLY: enter_exit

    INTEGER, INTENT(in) :: no_freq
    COMPLEX(dp), INTENT(inout) :: eff_admit(1:no_freq,num_elems)
    COMPLEX(dp), INTENT(inout) :: char_admit(1:no_freq,num_elems)
    TYPE(all_admit_param) :: admit_param
    REAL(dp), INTENT(in) :: harmonic_scale
    REAL(dp), INTENT(in) :: density
    REAL(dp), INTENT(in) :: viscosity
    CHARACTER(len=60), INTENT(in) :: mesh_type

    TYPE(elasticity_param) :: elast_param

    !local variables
    INTEGER :: nf,ne,nunit
    REAL(dp) :: omega,R1,R2,C,length,radius,C_term,E,h_bar
    REAL(dp) ::  h,L_term,R_term,vein_res
    COMPLEX(dp) :: term_admit

    CHARACTER(len=60) :: sub_name
    sub_name = 'boundary_admittance'
    CALL enter_exit(sub_name,1)
    IF(admit_param%bc_type.EQ.'two_unit_wk')THEN
       R1=admit_param%two_parameter%admit_P1
       C=admit_param%two_parameter%admit_P2
       DO nf=1,no_freq !step through frequencies
          omega=nf*2*PI*harmonic_scale
          IF(mesh_type.EQ.'simple_tree')THEN
             DO nunit=1,num_units
                ne=units(nunit)
                !temporarily store in eff_admit, to be added to the char admit
                eff_admit(nf,ne)=(1.0_dp+CMPLX(0.0_dp,1.0_dp)*omega*R1*C)/R1
             ENDDO
          ELSEIF(mesh_type.EQ.'full_plus_ladder')THEN
             DO ne=1,num_elems
                IF(elem_cnct(1,0,ne).EQ.0)THEN
                   !temporarily store in eff_admit, to be added to the char admit
                   eff_admit(nf,ne)=(1.0_dp+CMPLX(0.0_dp,1.0_dp)*omega*R1*C)/R1
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ELSEIF(admit_param%bc_type.EQ.'three_unit_wk')THEN
       R1=admit_param%three_parameter%admit_P1
       R2=admit_param%three_parameter%admit_P2
       C=admit_param%three_parameter%admit_P3
       DO nf=1,no_freq !step through frequencies
          omega=nf*2*PI*harmonic_scale
          IF(mesh_type.EQ.'simple_tree')THEN
             DO nunit=1,num_units
                ne=units(nunit)
                !temporarily store in eff_admit, to be added to the char admit
                eff_admit(nf,ne)=(1+CMPLX(0.0_dp,1.0_dp)*omega*R2*C)/(R1+R2+CMPLX(0,1)*omega*R1*R2*C)
             ENDDO
          ELSEIF(mesh_type.EQ.'full_plus_ladder')THEN
             DO ne=1,num_elems
                IF(elem_cnct(1,0,ne).EQ.0)THEN
                   !temporarily store in eff_admit, to be added to the char admit
                   eff_admit(nf,ne)=(1+CMPLX(0,1)*omega*R2*C)/(R1+R2+CMPLX(0,1)*omega*R1*R2*C)
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ELSEIF(admit_param%bc_type.EQ.'two_wk_plus')THEN
       !special case for uterine arteries which are in parallel with shunts
       E=elast_param%elasticity_parameters(1) !Pa
       h_bar=elast_param%elasticity_parameters(1)!this is a fraction of the radius so is unitless
       vein_res=0.45_dp
       R1=admit_param%four_parameter%admit_P1
       C=admit_param%four_parameter%admit_P2
       length=admit_param%four_parameter%admit_P3
       radius=admit_param%four_parameter%admit_P4
       DO nf=1,no_freq !step through frequencies
          omega=nf*2*PI*harmonic_scale
          DO nunit=1,num_units
             ne=units(nunit)
             !temporarily store in eff_admit, to be added to the char admit
             !ADMITTANCE DUE TO THE TERMINAL LOAD
             eff_admit(nf,ne)=(1.0_dp+CMPLX(0,1.0_dp)*omega*R1*C)/R1
             ! A SECOND ADMITTANCE IN PARALLEL REPRESENTING SHUNTS
             h=h_bar*radius
             C_term=3.0_dp*PI*radius**3/(2.0_dp*h*E)!
             L_term=9.0_dp*density&
                  /(4.0_dp*PI*radius**2)!per unit length
             R_term=81.0_dp*viscosity/ &
                  (8.0_dp*PI*radius**4) !laminar resistance per unit length
             !G=0.0_dp
             term_admit=SQRT(CMPLX(0.0_dp,1.0_dp,8)*omega*C_term)&
                  /SQRT(R_term+CMPLX(0.0_dp,1.0_dp,8)*omega*L_term)*50.0_dp*1.0_dp
             term_admit=term_admit/(1+term_admit*vein_res)
             eff_admit(nf,ne)=term_admit+eff_admit(nf,ne)
          ENDDO
       ENDDO
    ELSEIF(admit_param%bc_type.EQ.'zero_reflection')THEN
       DO nf=1,no_freq
          IF(mesh_type.EQ.'simple_tree')THEN
             DO nunit=1,num_units
                ne=units(nunit)
                eff_admit(nf,ne)=char_admit(nf,ne) !in effective admittance subroutine this will become a 'dummy' daughter admittance and zero the reflection coefficient
             ENDDO
          ELSEIF(mesh_type.EQ.'full_plus_ladder')THEN
             DO ne=1,num_elems
                IF(elem_cnct(1,0,ne).EQ.0)THEN
                   !temporarily store in eff_admit
                   eff_admit(nf,ne)=char_admit(nf,ne)
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ENDIF
    CALL enter_exit(sub_name,2)
  END SUBROUTINE boundary_admittance

  !
  !##############################################################################
  !
  !*characteristic_admittance* calculates the characteristic admittance of each
  SUBROUTINE characteristic_admittance(no_freq,char_admit,prop_const,harmonic_scale,&
       density,viscosity,admit_param,elast_param,mechanics_parameters,grav_vect)
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAD:"SO_characteristic_admittance: characteristic_admittance
    USE other_consts, ONLY: MAX_STRING_LEN
    USE indices
    USE arrays, ONLY: num_elems,elem_field,elasticity_param,all_admit_param,elem_nodes
    USE pressure_resistance_flow, ONLY: calculate_ppl
    USE math_utilities, ONLY: bessel_complex
    USE diagnostics, ONLY: enter_exit

    INTEGER, INTENT(in) :: no_freq
    COMPLEX(dp), INTENT(inout) :: char_admit(1:no_freq,num_elems)
    COMPLEX(dp), INTENT(inout) :: prop_const(1:no_freq,num_elems)
    REAL(dp), INTENT(in) :: harmonic_scale
    REAL(dp), INTENT(in) :: density
    REAL(dp), INTENT(in) :: viscosity
    REAL(dp),INTENT(in) :: mechanics_parameters(2),grav_vect(3)

    TYPE(elasticity_param) :: elast_param
    TYPE(all_admit_param) :: admit_param

    !local variables
    REAL(dp) :: L,C,R, G,omega,gen_factor
    REAL(dp) :: E,h_bar,h,wavespeed,wolmer !should be global - maybe express as alpha (i.e. pre multiply)
    COMPLEX(dp) :: f10,bessel0,bessel1
    INTEGER :: ne,nf,nn,np
    INTEGER :: exit_status=0
    REAL(dp) :: R0,Ppl,Ptm,Rg_in,Rg_out
    CHARACTER(len=60) :: sub_name
    sub_name = 'characteristic_admittance'
    CALL enter_exit(sub_name,1)

    DO ne=1,num_elems
       DO nn=1,2
          IF(nn.EQ.1) np=elem_nodes(1,ne)
          IF(nn.EQ.2) np=elem_nodes(2,ne)
          CALL calculate_ppl(np,grav_vect,mechanics_parameters,Ppl)
          Ptm=Ppl     ! Pa
          IF(nn.EQ.1)R0=elem_field(ne_radius_in0,ne)
          IF(nn.EQ.2)R0=elem_field(ne_radius_out0,ne)
          IF(admit_param%admittance_type.EQ.'duan_zamir')THEN!alpha controls elasticity
             IF(nn.EQ.1)Rg_in=R0*(Ptm*elast_param%elasticity_parameters(1)+1.d0)
             IF(nn.EQ.2)Rg_out=R0*(Ptm*elast_param%elasticity_parameters(1)+1.d0)
          ELSE!Hooke type elasticity
             h=elast_param%elasticity_parameters(2)*R0
             IF(nn.EQ.1) Rg_in=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elast_param%elasticity_parameters(1)*h)
             IF(nn.EQ.2) Rg_out=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elast_param%elasticity_parameters(1)*h)
          ENDIF
       ENDDO
       elem_field(ne_radius_out,ne)=(Rg_in+Rg_out)/2.0_dp

    ENDDO


    E=elast_param%elasticity_parameters(1) !Pa
    h_bar=elast_param%elasticity_parameters(2)!this is a fraction of the radius so is unitless
    DO ne=1,num_elems
       IF(admit_param%admittance_type.EQ.'lachase_standard')THEN
          h=h_bar*elem_field(ne_radius_out,ne)
          C=3.0_dp*PI*elem_field(ne_radius_out,ne)**3*elem_field(ne_length,ne)/(2.0_dp*h*E)
          L=density*elem_field(ne_length,ne)/(4*PI*elem_field(ne_radius_out,ne)**2)
          R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
               (PI*elem_field(ne_radius_out,ne)**4) !laminar resistance
          G=0.0_dp
       ELSEIF(admit_param%admittance_type.EQ.'lachase_modified')THEN
          h=h_bar*elem_field(ne_radius_out,ne)
          C=3.0_dp*PI*elem_field(ne_radius_out,ne)**3/(2.0_dp*h*E)!
          L=9.0_dp*density&
               /(4.0_dp*PI*elem_field(ne_radius_out,ne)**2)!per unit length
          R=81.0_dp*viscosity/ &
               (8.0_dp*PI*elem_field(ne_radius_out,ne)**4) !laminar resistance per unit length
          G=0.0_dp
       ELSEIF(admit_param%admittance_type.EQ.'zhu_chesler')THEN
          h=h_bar*elem_field(ne_radius_out,ne)
          C=3.0_dp*PI*elem_field(ne_radius_out,ne)**3*elem_field(ne_length,ne)/(2.0_dp*h*E)
          L=9.0_dp*density*elem_field(ne_length,ne)/(4.0_dp*PI*elem_field(ne_radius_out,ne)**2)
          R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
               (PI*elem_field(ne_radius_out,ne)**4) !laminar resistance
          G=0.0_dp
       ELSEIF(admit_param%admittance_type.EQ.'duan_zamir')THEN
          DO nf=1,no_freq !radius needs to  be multipled by 1000 to go to mm (units of rest of model)
             omega=nf*2*PI*harmonic_scale!q/s
             wolmer=(elem_field(ne_radius_out,ne))*SQRT(omega*density/viscosity)
             CALL bessel_complex(wolmer*CMPLX(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp),bessel0,bessel1)
             f10=2*bessel1/(wolmer*CMPLX(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp)*bessel0)!no units
             wavespeed=SQRT(1.0_dp/(2*density*elast_param%elasticity_parameters(1)))*SQRT(1-f10)! !mm/s
             char_admit(nf,ne)=PI*(elem_field(ne_radius_out,ne))**2/(density*wavespeed/(1-f10))*SQRT(1-f10)!mm3/Pa
             prop_const(nf,ne)=CMPLX(0.0_dp,1.0_dp,8)*omega/(wavespeed)!1/mm
          ENDDO
       ELSE !Unrecognised admittance model
          PRINT *, "EXITING"
          PRINT *, "Unrecognised admittance model, please check inputs"
          CALL EXIT(exit_status)
       ENDIF
       IF(admit_param%admittance_type.EQ.'duan_zamir')THEN
       ELSE
          DO nf=1,no_freq
             omega=nf*2*PI*harmonic_scale
             char_admit(nf,ne)=SQRT(G+CMPLX(0.0_dp,1.0_dp,8)*omega*C)/SQRT(R+CMPLX(0.0_dp,1.0_dp,8)*omega*L)!mm3/Pa.s
             prop_const(nf,ne)=SQRT((G+CMPLX(0.0_dp,1.0_dp,8)*omega*C)*(R+CMPLX(0.0_dp,1.0_dp,8)*omega*L))!1/mm
          ENDDO!nf
       ENDIF
    ENDDO!ne

    CALL enter_exit(sub_name,2)
  END SUBROUTINE characteristic_admittance

  !##################################################################
  !
  !*tree_admittance:* Calculates the total admittance of a tree
  SUBROUTINE tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
       min_elem,max_elem,tree_direction)
    USE indices
    USE arrays,ONLY: dp,num_elems,elem_cnct,elem_field
    USE diagnostics, ONLY: enter_exit
    INTEGER, INTENT(in) :: no_freq
    COMPLEX(dp), INTENT(inout) :: eff_admit(1:no_freq,num_elems)
    COMPLEX(dp), INTENT(in) :: char_admit(1:no_freq,num_elems)
    COMPLEX(dp), INTENT(inout) :: reflect(1:no_freq,num_elems)
    COMPLEX(dp), INTENT(in) :: prop_const(1:no_freq,num_elems)
    REAL(dp), INTENT(in) :: harmonic_scale
    INTEGER, INTENT(in) :: min_elem,max_elem
    CHARACTER(len=30), INTENT(in) :: tree_direction

    CHARACTER(len=60) :: sub_name
    !local variables
    REAL(dp) :: invres,elem_res(num_elems),omega
    INTEGER :: num2,ne,ne2,nf,num3,ne3,ne_sist
    COMPLEX(dp) :: daughter_admit,sister_admit,sister_current,a,b,m1,m2

    sub_name = 'tree_admittance'
    CALL enter_exit(sub_name,1)
    reflect(:,:)=CMPLX(0.0_dp,0.0_dp,8)

    IF(tree_direction.EQ.'diverging')THEN
       DO nf=1,no_freq
          omega=nf*2*PI*harmonic_scale
          DO ne=max_elem,min_elem,-1!step backward through elements
             daughter_admit=CMPLX(0.0_dp,0.0_dp,8)!

             DO num2=1,elem_cnct(1,0,ne)!will only do stuff to non-terminals will add one daughter if no branching
                ne2=elem_cnct(1,num2,ne)!for each downstream element
                daughter_admit=daughter_admit+eff_admit(nf,ne2)!sum of two child elements
             ENDDO
             IF(elem_cnct(1,0,ne).GT.0)THEN !not a terminal as it has a downstream
                reflect(nf,ne)=(char_admit(nf,ne)-daughter_admit)/&
                     (char_admit(nf,ne)+daughter_admit)!double checked
                eff_admit(nf,ne)=char_admit(nf,ne)*(1&
                     -reflect(nf,ne)*EXP(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
                     (1+reflect(nf,ne)*EXP(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))
             ELSE!a terminal
                daughter_admit=eff_admit(nf,ne) !a boundary condition is applied here
                reflect(nf,ne)=(char_admit(nf,ne)-daughter_admit)/&
                     (char_admit(nf,ne)+daughter_admit)
                !now we overwrite the effective admittance of the terminal to include reflection from the daughter.
                eff_admit(nf,ne)=char_admit(nf,ne)*(1&
                     -reflect(nf,ne)*EXP(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
                     (1&
                     +reflect(nf,ne)*EXP(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))
             ENDIF
          ENDDO!ne
       ENDDO!nf
    ELSEIF(tree_direction.EQ.'converging')THEN
       DO nf=1,no_freq
          omega=nf*2*PI*harmonic_scale
          DO ne=min_elem,max_elem!step forward through elements
             daughter_admit=CMPLX(0.0_dp,0.0_dp,8)!
             sister_admit=CMPLX(0.0_dp,0.0_dp,8)!
             DO num2=1,elem_cnct(1,0,ne)!will only do stuff to non-terminals
                ne2=elem_cnct(1,num2,ne)!for each downstream element
                daughter_admit=daughter_admit+eff_admit(nf,ne2)!sum of two child elements
                DO num3=1,elem_cnct(-1,0,ne2)!sisters
                   ne3=elem_cnct(-1,num3,ne2)!for each upstream element of the daughter
                   IF(ne3.NE.ne)THEN
                      ne_sist=ne3
                      sister_admit=sister_admit+char_admit(nf,ne3)
                      sister_current=EXP(-1.0_dp*prop_const(nf,ne3)*elem_field(ne_length,ne3))/&
                           EXP(-1.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne))
                      a=EXP(-1.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne))
                      b=EXP(-1.0_dp*prop_const(nf,ne3)*elem_field(ne_length,ne3))
                   ENDIF
                ENDDO
             ENDDO
             m1=1.0_dp+2.0_dp*(b/a)*((2*a*b-a**2-1)*char_admit(nf,ne)+&
                  (a**2-1)*sister_admit+(a**2-1)*daughter_admit)/((1-b**2)*char_admit(nf,ne)+&
                  (1+b**2-2*a*b)*sister_admit+(1-b**2)*daughter_admit)
             IF(elem_cnct(1,0,ne).GT.0)THEN !not a terminal
                reflect(nf,ne)=(char_admit(nf,ne)-m1*sister_admit-daughter_admit)/&
                     (char_admit(nf,ne)+daughter_admit+sister_admit)!ARC- to check
                eff_admit(nf,ne)=char_admit(nf,ne)*(1&
                     -reflect(nf,ne)*EXP(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
                     (1+reflect(nf,ne)*EXP(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))
             ELSE!a terminal
                daughter_admit=eff_admit(nf,ne) !a boundary condition is applied here
                reflect(nf,ne)=(char_admit(nf,ne)-daughter_admit)/&
                     (char_admit(nf,ne)+daughter_admit)
                eff_admit(nf,ne)=char_admit(nf,ne)*(1&
                     -reflect(nf,ne)*EXP(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
                     (1+reflect(nf,ne)*EXP(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))
             ENDIF
          ENDDO!ne
       ENDDO!nf
    ENDIF

    CALL enter_exit(sub_name,2)
  END SUBROUTINE tree_admittance
  !##################################################################
  !
  !*capillaryadmittance:* Calculates the total admittance of a tree
  SUBROUTINE capillary_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
       min_elem,max_elem,elast_param,mechanics_parameters,grav_vect,cap_model)
    USE indices
    USE arrays,ONLY: dp,num_elems,elem_cnct,elem_field,capillary_bf_parameters,elem_nodes,&
         node_field,node_xyz,elasticity_param
    USE pressure_resistance_flow,ONLY: calculate_ppl
    USE capillaryflow, ONLY:cap_flow_admit
    USE diagnostics, ONLY: enter_exit

    INTEGER, INTENT(in) :: no_freq
    COMPLEX(dp), INTENT(inout) :: eff_admit(1:no_freq,num_elems)
    COMPLEX(dp), INTENT(inout) :: char_admit(1:no_freq,num_elems)
    COMPLEX(dp), INTENT(inout) :: reflect(1:no_freq,num_elems)
    COMPLEX(dp), INTENT(in) :: prop_const(1:no_freq,num_elems)
    REAL(dp), INTENT(in) :: harmonic_scale
    INTEGER, INTENT(in) :: min_elem,max_elem
    REAL(dp),INTENT(in) :: mechanics_parameters(2),grav_vect(3)
    INTEGER, INTENT(in) :: cap_model

    TYPE(capillary_bf_parameters) :: cap_param
    TYPE(elasticity_param) :: elast_param

    CHARACTER(len=60) :: sub_name
    !local variables
    INTEGER :: ne, ne2,num2,nf,ne0,ne1
    REAL(dp) :: daughter_admit,omega,length,Hart,Hven,Gamma_sheet,alpha_c,Ptp
    INTEGER :: grav_dirn,i
    REAL (dp) :: Ppl,P1,P2
    REAL(dp) :: Q01,Rin,Rout,Lin,Lout,x_cap,y_cap,z_cap
    COMPLEX(dp) :: eff_admit_downstream(no_freq)

    sub_name = 'capillary_admittance'
    CALL enter_exit(sub_name,1)
    reflect(:,:)=CMPLX(0.0_dp,0.0_dp,8)

    DO ne=min_elem,max_elem
       ne0=elem_cnct(-1,1,ne)!upstream element number
       ne1=elem_cnct(1,1,ne) !downstream element number
       P1=node_field(nj_bv_press,elem_nodes(2,ne0)) !pressure at start node of capillary element
       P2=node_field(nj_bv_press,elem_nodes(1,ne1))
       Q01=elem_field(ne_Qdot,ne0) !flow in element upstream of capillary element !mm^3/s
       Rin=elem_field(ne_radius_out0,ne0)!radius of upstream element
       Rout=elem_field(ne_radius_out0,ne1) !radius of downstream element
       x_cap=node_xyz(1,elem_nodes(1,ne))
       y_cap=node_xyz(2,elem_nodes(1,ne))
       z_cap=node_xyz(3,elem_nodes(1,ne))
       CALL calculate_ppl(elem_nodes(1,ne),grav_vect,mechanics_parameters,Ppl)
       Lin=elem_field(ne_length,ne0)
       Lout=elem_field(ne_length,ne1)
       Ptp=(cap_param%Palv-(-Ppl))/98.06d0 !Pa -> cmH2O
       DO i=1,no_freq
          eff_admit_downstream(i)=eff_admit(i,ne1)
       ENDDO
       CALL cap_flow_admit(ne,eff_admit(:,ne),eff_admit_downstream,Lin,Lout,P1,P2,&
            Ppl,Q01,Rin,Rout,x_cap,y_cap,z_cap,no_freq,harmonic_scale,&
            elast_param,cap_model)
    ENDDO!ne

  END SUBROUTINE capillary_admittance
  !
  !################################################
  !
  !*pressure_factor:* Calculates change in pressure through tree
  SUBROUTINE pressure_factor(no_freq,p_factor,reflect,prop_const,harmonic_scale,ne_min,ne_max)
    USE indices
    USE arrays,ONLY: dp,num_elems,elem_cnct,elem_field
    USE diagnostics, ONLY: enter_exit
    INTEGER, INTENT(in) :: no_freq
    COMPLEX(dp), INTENT(inout) :: p_factor(1:no_freq,num_elems)
    COMPLEX(dp), INTENT(inout) :: reflect(1:no_freq,num_elems)
    COMPLEX(dp), INTENT(in) :: prop_const(1:no_freq,num_elems)
    REAL(dp), INTENT(in) :: harmonic_scale
    INTEGER, INTENT(in) :: ne_min,ne_max



    CHARACTER(len=60) :: sub_name
    !local variables
    INTEGER :: ne, nf,ne_up
    REAL(dp) :: omega

    sub_name = 'pressure_factor'
    CALL enter_exit(sub_name,1)

    p_factor=1.0_dp
    DO nf=1,no_freq
       omega=nf*2*PI*harmonic_scale
       DO ne=ne_min,ne_max
          !look for upstream element
          IF(elem_cnct(-1,0,ne).EQ.0)THEN !no upstream elements, inlet, ignore
             ne_up=ne_min
             p_factor(nf,ne)=(1.0_dp)!* &!assumes input admittance is the same as characteristic admittance for this vessel
             !exp(-1.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne))!/&
             !(1+reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))
          ELSE
             ne_up=elem_cnct(-1,1,ne)
             p_factor(nf,ne)=p_factor(nf,ne_up)*(1+reflect(nf,ne_up))* &
                  EXP(-1.0_dp*elem_field(ne_length,ne_up)*prop_const(nf,ne_up))/&
                  (1+reflect(nf,ne)*EXP(-2.0_dp*elem_field(ne_length,ne)*prop_const(nf,ne)))
          ENDIF!neup
       ENDDO!ne
    ENDDO!nf

    CALL enter_exit(sub_name,2)
  END SUBROUTINE pressure_factor
END MODULE wave_transmission
