MODULE pressure_resistance_flow
  !*Brief Description:* This module contains tools that are used to solve systems of equations representing steady pressure, resistance and flow problems in any branching geometry. The subroutines in this module are core subroutines that are used in many problem types and are applicable beyond lung modelling
  !
  !*LICENSE:*
  !TBC
  !
  !
  !*Full Description:*
  !
  !This module contains tools that are used to solve systems of equations representing steady pressure, resistance and flow problems in any branching geometry. The subroutines in this module are core subroutines that are used in many problem types and are applicable beyond lung modelling
  USE solve, ONLY: BICGSTAB_LinSolv,pmgmres_ilu_cr
  USE other_consts, ONLY: TOLERANCE
  IMPLICIT NONE
  !Module parameters

  !Module types

  !Module depvar

  !Interfaces
  PRIVATE
  PUBLIC evaluate_prq,calculate_ppl
CONTAINS
  !###################################################################################
  !
  !*evaluate_PRQ:* Solves for pressure and flow in a rigid or compliant tree structure
  SUBROUTINE evaluate_prq(mesh_type,vessel_type,grav_dirn,grav_factor,bc_type,inlet_bc,outlet_bc)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_PRQ" :: EVALUATE_PRQ
    USE indices
    USE capillaryflow,ONLY: cap_flow_ladder
    USE arrays,ONLY: dp,num_elems,num_nodes,elem_field,elem_nodes,elem_cnct,node_xyz
    USE diagnostics, ONLY: enter_exit
    !local variables
    INTEGER :: mesh_dof,depvar_types
    INTEGER, ALLOCATABLE :: mesh_from_depvar(:,:,:)
    INTEGER, ALLOCATABLE :: depvar_at_node(:,:,:)
    INTEGER, ALLOCATABLE :: depvar_at_elem(:,:,:)
    INTEGER, DIMENSION(0:2,2) :: depvar_totals
    INTEGER, ALLOCATABLE :: SparseCol(:)
    INTEGER, ALLOCATABLE :: SparseRow(:)
    INTEGER, ALLOCATABLE :: update_resistance_entries(:)
    REAL(dp), ALLOCATABLE :: SparseVal(:)
    REAL(dp), ALLOCATABLE :: RHS(:)
    INTEGER :: num_vars,NonZeros,MatrixSize
    INTEGER :: AllocateStatus

    REAL(dp), ALLOCATABLE :: prq_solution(:,:),solver_solution(:)
    REAL(dp) :: viscosity,density,inlet_bc,outlet_bc,inletbc,outletbc,grav_vect(3),gamma,total_resistance,ERR
    LOGICAL, ALLOCATABLE :: FIX(:)
    LOGICAL :: ADD=.FALSE.,CONVERGED=.FALSE.
    CHARACTER(len=60) :: sub_name,mesh_type,vessel_type,mechanics_type,bc_type
    INTEGER :: grav_dirn,no,depvar,KOUNT,nz,ne,SOLVER_FLAG,ne0,ne1,nj
    REAL(dp) :: MIN_ERR,N_MIN_ERR,elasticity_parameters(3),mechanics_parameters(2),grav_factor,P1
    REAL(dp) :: P2,Q01,Rin,Rout,x_cap,y_cap,z_cap,Ppl,LPM_R,Lin,Lout
    INTEGER :: update_flow_nzz_row

    sub_name = 'evaluate_prq'
    CALL enter_exit(sub_name,1)
    !!---------DESCRIPTION OF MODEL Types -----------
    !mesh_type: can be simple_tree, full_plus_ladder, full_sheet, full_tube The first can be airways, arteries, veins but no special features at the terminal level, the last one has arteries and veins connected by capillary units of some type (lung ladder acinus, lung sheet capillary bed, capillaries are just tubes represented by an element)

    !vessel_type:
    !rigid, no elasticity, no parameters required
    !elastic_g0_beta, R=R0*((Ptm/G0)+1.d0)^(1.d0/elasticity_parameters(2)),with an optional maximum pressure beyond which the vessel radius is constant three parameters, g0, elasticity_parameters(2), elasticity_parameters(3)
    !elastic alpha,  R=R0*(alpha*Ptm+1.d0), up to a limit elasticity_parameters(3) two parameters alpha, elasticity_parameters(3)
    !elastic_hooke, two parameters E and h,R=R0+3.0_dp*R0**2*Ptm/(4.0_dp*E*h*R0)

    !mechanics type:
    !linear two parmeters, transpulmonary pressure (average) and pleural density (gradient)
    !mechanics, two parameters, pressure and stretch fields

    !bc_type:
    !pressure (at inlet and outlets)
    !flow (flow at inlet pressure at outlet).


    mechanics_type='linear'

    IF (vessel_type.EQ.'rigid') THEN
       elasticity_parameters=0.0_dp
    ELSEIF (vessel_type.EQ.'elastic_g0_beta') THEN
       elasticity_parameters(1)=6.67e3_dp!G0 (Pa)
       elasticity_parameters(2)=1.0_dp!elasticity_parameters(2)
       elasticity_parameters(3)=32.0_dp*98.07_dp !elasticity_parameters(3) (Pa)
    ELSEIF (vessel_type.EQ.'elastic_alpha') THEN
       elasticity_parameters(1)=1.503e-4_dp!alpha (1/Pa)
       elasticity_parameters(2)=32.0_dp*98.07_dp !elasticity_parameters(3) (Pa)
       elasticity_parameters(3)=0.0_dp !Not used
    ELSEIF (vessel_type.EQ.'elastic_hooke') THEN
       elasticity_parameters(1)=1.5e6_dp !Pa
       elasticity_parameters(2)=0.1_dp!this is a fraction of the radius so is unitless
       elasticity_parameters(3)=0.0_dp !Not used
    ELSE
       PRINT *, 'WARNING: Your chosen vessel type does not seem to be implemented assuming rigid'
       vessel_type='rigid'
       elasticity_parameters=0.0_dp
    ENDIF

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

    IF(bc_type.EQ.'pressure')THEN
       inletbc=inlet_bc
       outletbc=outlet_bc
    ELSEIF(bc_type.EQ.'flow')THEN
       inletbc=inlet_bc
       outletbc=outlet_bc
    ELSEIF((bc_type.NE.'pressure').AND.(bc_type.NE.'flow'))THEN
       PRINT *,"unsupported bc_type",bc_type
       CALL EXIT(1)
    ENDIF

    !!---------PHYSICAL PARAMETERS-----------
    !viscosity: fluid viscosity
    !density:fluid density
    !gamma:Pedley correction factor
    density=0.10500e-02_dp !kg/cm3
    viscosity=0.33600e-02_dp !Pa.s
    gamma = 0.327_dp !=1.85/(4*sqrt(2))

    !! Allocate memory to depvar arrays
    mesh_dof=num_elems+num_nodes
    depvar_types=2 !pressure/flow
    ALLOCATE (mesh_from_depvar(0:2,mesh_dof,0:2), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for mesh_from_depvar array ***"
    ALLOCATE (depvar_at_elem(0:2,depvar_types,num_elems), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for depvar_at_elem array ***"
    ALLOCATE (depvar_at_node(num_nodes,0:2,depvar_types), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for depvar_at_node array ***"
    ALLOCATE (prq_solution(mesh_dof,2), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for prq_solution array ***"
    prq_solution=0.0_dp !initialise
    ALLOCATE (FIX(mesh_dof), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for FIX array ***"



    !! Setting up mappings between nodes, elements and solution depvar
    CALL calc_depvar_maps(mesh_from_depvar,depvar_at_elem,&
         depvar_totals,depvar_at_node,mesh_dof,num_vars)

    !! Define boundary conditions
    !first call to define inlet boundary conditions
    CALL boundary_conditions(ADD,FIX,bc_type,grav_vect,density,inletbc,outletbc,&
         depvar_at_node,depvar_at_elem,prq_solution,mesh_dof,mesh_type)
    !second call if simple tree need to define pressure bcs at all terminal branches
    IF(mesh_type.EQ.'simple_tree')THEN
       ADD=.TRUE.
       CALL boundary_conditions(ADD,FIX,bc_type,grav_vect,density,inletbc,outletbc,&
            depvar_at_node,depvar_at_elem,prq_solution,mesh_dof,mesh_type)
    ELSEIF(mesh_type.EQ.'full_plus_ladder')THEN
       ADD=.TRUE.
       CALL boundary_conditions(ADD,FIX,bc_type,grav_vect,density,inletbc,outletbc,&
            depvar_at_node,depvar_at_elem,prq_solution,mesh_dof,mesh_type)
    ENDIF

    KOUNT=0
    !! Calculate resistance of each element
    CALL calculate_resistance(viscosity,KOUNT)

    !! Calculate sparsity structure for solution matrices
    !Determine size of and allocate solution vectors/matrices
    CALL calc_sparse_size(mesh_dof,depvar_at_elem,depvar_at_node,FIX,NonZeros,MatrixSize)
    ALLOCATE (SparseCol(NonZeros), STAT = AllocateStatus)!Note we should be able to calculate the nonzeros and matrix size analtyically then we wont need this.
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for SparseCol array ***"
    ALLOCATE (SparseRow(MatrixSize+1), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for SparseRow array ***"
    ALLOCATE (SparseVal(NonZeros), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for SparseVal array ***"
    ALLOCATE (RHS(MatrixSize), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for RHS array ***"
    ALLOCATE (solver_solution(MatrixSize), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    ALLOCATE (update_resistance_entries(num_elems), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    update_resistance_entries = 0
    !calculate the sparsity structure
    CALL calc_sparse_1dtree(bc_type,density,FIX,grav_vect,mesh_dof,depvar_at_elem, &
         depvar_at_node,NonZeros,MatrixSize,SparseCol,SparseRow,SparseVal,RHS, &
         prq_solution,update_resistance_entries,update_flow_nzz_row)
!!! --ITERATIVE LOOP--
    MIN_ERR=1.d10
    N_MIN_ERR=0
    DO WHILE(.NOT.CONVERGED)
       KOUNT=KOUNT+1
       PRINT*, 'Outer loop iterations:',KOUNT
!!! Initialise solution vector based on bcs and rigid vessel resistance
       IF(KOUNT.EQ.1)THEN!set up boundary conditions
          IF(bc_type.EQ.'pressure')THEN
             IF(mesh_type.EQ.'full_plus_ladder')THEN
                total_resistance=1000.0_dp
             ELSE
                CALL tree_resistance(total_resistance)
             ENDIF
             CALL initialise_solution(inletbc,outletbc,(inletbc-outletbc)/total_resistance, &
                  mesh_dof,prq_solution,depvar_at_node,depvar_at_elem,FIX)
             !move initialisation to solver solution (skipping BCs).
             no=0
             DO depvar=1,mesh_dof !loop over mesh dofs
                IF(.NOT.FIX(depvar))THEN
                   no=no+1
                   solver_solution(no)=prq_solution(depvar,1)
                ENDIF
             ENDDO !mesh_dof
          ELSE!flow BCs to be implemented
          ENDIF
       ELSE!Need to update just the resistance values in the solution matrix
          DO ne=1,num_elems !update for all ne
             IF(update_resistance_entries(ne).GT.0)THEN
                nz=update_resistance_entries(ne)
                SparseVal(nz)=-elem_field(ne_resist,ne) !Just updating resistance
             ENDIF
          ENDDO
          IF(bc_type.EQ.'flow')THEN !update RHS to account for element resistance
             DO ne=1,num_elems
                depvar = depvar_at_elem(1,1,ne)
                IF(FIX(depvar))THEN
                   RHS(update_flow_nzz_row) = prq_solution(depvar,1)*elem_field(ne_resist,ne)
                ENDIF
             ENDDO
             !SparseVal(nz)=-elem_field(ne_resist,ne) !Just updating resistance
          ENDIF
       ENDIF!first or subsequent iteration
       !! ----CALL SOLVER----
       CALL pmgmres_ilu_cr(MatrixSize, NonZeros, SparseRow, SparseCol, SparseVal, &
            solver_solution, RHS, 500, 500,1.d-5,1.d-4,SOLVER_FLAG)
       IF(SOLVER_FLAG == 0)THEN
          PRINT *, 'Warning: pmgmres has reached max iterations. Solution may not be valid if this warning persists'
       ELSEIF(SOLVER_FLAG ==2)THEN
          PRINT *, 'ERROR: pmgmres has failed to converge'
          DEALLOCATE (mesh_from_depvar, STAT = AllocateStatus)
          DEALLOCATE (depvar_at_elem, STAT = AllocateStatus)
          DEALLOCATE (depvar_at_node, STAT = AllocateStatus)
          DEALLOCATE (prq_solution, STAT = AllocateStatus)
          DEALLOCATE (FIX, STAT = AllocateStatus)
          DEALLOCATE (solver_solution, STAT = AllocateStatus)
          DEALLOCATE (SparseCol, STAT = AllocateStatus)
          DEALLOCATE (SparseVal, STAT = AllocateStatus)
          DEALLOCATE (SparseRow, STAT = AllocateStatus)
          DEALLOCATE (RHS, STAT = AllocateStatus)
          DEALLOCATE (update_resistance_entries, STAT=AllocateStatus)
          EXIT
       ENDIF
       !!--TRANSFER SOLVER SOLUTIONS TO FULL SOLUTIONS
       ERR=0.0_dp
       no=0
       DO depvar=1,mesh_dof
          IF(.NOT.FIX(depvar)) THEN
             no=no+1
             prq_solution(depvar,2)=prq_solution(depvar,1) !temp storage of previous solution
             prq_solution(depvar,1)=solver_solution(no) !new pressure & flow solutions
             IF(DABS(prq_solution(depvar,1)).GT.0.d-6)THEN
                ERR=ERR+(prq_solution(depvar,2)-prq_solution(depvar,1))**2.d0/prq_solution(depvar,1)**2
             ENDIF
          ENDIF
       ENDDO !no2
       !rigid vessels no need to update - tag as converged and exit
       IF(vessel_type.EQ.'rigid')THEN
          ERR=0.0_dp
          CONVERGED=.TRUE.
       ELSE
          !Update vessel radii based on predicted pressures and then update resistance through tree
          CALL calc_press_area(grav_vect,KOUNT,depvar_at_node,prq_solution,&
               mesh_dof,vessel_type,elasticity_parameters,mechanics_parameters)
          CALL calculate_resistance(viscosity,KOUNT)

          !Put the ladder stuff here --> See solve11.f
          IF(mesh_type.EQ.'full_plus_ladder')THEN
             DO ne=1,num_elems
                IF(elem_field(ne_group,ne).EQ.1.0_dp)THEN!(elem_field(ne_group,ne)-1.0_dp).lt.TOLERANCE)then
                   ne0=elem_cnct(-1,1,ne)!upstream element number
                   ne1=elem_cnct(1,1,ne)
                   P1=prq_solution(depvar_at_node(elem_nodes(2,ne0),0,1),1) !pressure at start node of capillary element
                   P2=prq_solution(depvar_at_node(elem_nodes(1,ne1),0,1),1)!pressure at end node of capillary element
                   Q01=prq_solution(depvar_at_elem(1,1,ne0),1) !flow in element upstream of capillary element !mm^3/s
                   Rin=elem_field(ne_radius_out0,ne0)!radius of upstream element
                   Rout=elem_field(ne_radius_out0,ne1) !radius of downstream element
                   x_cap=node_xyz(1,elem_nodes(1,ne))
                   y_cap=node_xyz(2,elem_nodes(1,ne))
                   z_cap=node_xyz(3,elem_nodes(1,ne))
                   CALL calculate_ppl(elem_nodes(1,ne),grav_vect,mechanics_parameters,Ppl)
                   Lin=elem_field(ne_length,ne0)
                   Lout=elem_field(ne_length,ne1)
                   CALL cap_flow_ladder(ne,LPM_R,Lin,Lout,P1,P2,&
                        Ppl,Q01,Rin,Rout,x_cap,y_cap,z_cap,&
                        .FALSE.)
                   elem_field(ne_resist,ne)=LPM_R
                ENDIF
             ENDDO
          ENDIF

          ERR=ERR/MatrixSize !sum of error divided by no of unknown depvar
          IF(ERR.LE.1.d-6.AND.(KOUNT.NE.1))THEN
             CONVERGED=.TRUE.
             PRINT *,"Convergence achieved after",KOUNT,"iterations",ERR
          ELSE !if error not converged
             IF(ERR.GE.MIN_ERR) THEN
                N_MIN_ERR=N_MIN_ERR+1
             ELSE
                MIN_ERR=ERR
             ENDIF
             PRINT *,"Not converged, error =",ERR
          ENDIF !ERR not converged
       ENDIF!vessel type
    ENDDO !notconverged

    !need to write solution to element/nodal fields for export
    CALL map_solution_to_mesh(prq_solution,depvar_at_elem,depvar_at_node,mesh_dof)
    !NEED TO UPDATE TERMINAL SOLUTION HERE. LOOP THO' UNITS AND TAKE FLOW AND PRESSURE AT TERMINALS
    CALL map_flow_to_terminals
    !EXPORT LADDER SOLUTION
    IF(mesh_type.EQ.'full_plus_ladder')THEN
       OPEN(10, file='micro_flow_ladder.out', status='replace')
       OPEN(20, file='micro_flow_unit.out', status='replace')
       DO ne=1,num_elems
          IF(elem_field(ne_group,ne).EQ.1.0_dp)THEN!(elem_field(ne_group,ne)-1.0_dp).lt.TOLERANCE)then
             ne0=elem_cnct(-1,1,ne)!upstream element number
             ne1=elem_cnct(1,1,ne)
             P1=prq_solution(depvar_at_node(elem_nodes(2,ne0),0,1),1) !pressure at start node of capillary element
             P2=prq_solution(depvar_at_node(elem_nodes(1,ne1),0,1),1)!pressure at end node of capillary element
             Q01=prq_solution(depvar_at_elem(1,1,ne0),1) !flow in element upstream of capillary element !mm^3/s
             Rin=elem_field(ne_radius_out0,ne0)!radius of upstream element
             Rout=elem_field(ne_radius_out0,ne1) !radius of downstream element
             x_cap=node_xyz(1,elem_nodes(1,ne))
             y_cap=node_xyz(2,elem_nodes(1,ne))
             z_cap=node_xyz(3,elem_nodes(1,ne))
             CALL calculate_ppl(elem_nodes(1,ne),grav_vect,mechanics_parameters,Ppl)
             Lin=elem_field(ne_length,ne0)
             Lout=elem_field(ne_length,ne1)
             CALL cap_flow_ladder(ne,LPM_R,Lin,Lout,P1,P2,&
                  Ppl,Q01,Rin,Rout,x_cap,y_cap,z_cap,&
                  .TRUE.)
          ENDIF
       ENDDO
       CLOSE(10)
       CLOSE(20)
    ENDIF

    DEALLOCATE (mesh_from_depvar, STAT = AllocateStatus)
    DEALLOCATE (depvar_at_elem, STAT = AllocateStatus)
    DEALLOCATE (depvar_at_node, STAT = AllocateStatus)
    DEALLOCATE (prq_solution, STAT = AllocateStatus)
    DEALLOCATE (FIX, STAT = AllocateStatus)
    DEALLOCATE (solver_solution, STAT = AllocateStatus)
    DEALLOCATE (SparseCol, STAT = AllocateStatus)
    DEALLOCATE (SparseVal, STAT = AllocateStatus)
    DEALLOCATE (SparseRow, STAT = AllocateStatus)
    DEALLOCATE (RHS, STAT = AllocateStatus)
    DEALLOCATE (update_resistance_entries, STAT=AllocateStatus)
    CALL enter_exit(sub_name,2)
  END SUBROUTINE evaluate_prq
  !
  !###################################################################################
  !
  !*boundary_conditions:* Defines boundary conditions for prq problems
  SUBROUTINE boundary_conditions(ADD,FIX,bc_type,grav_vect,density,inletbc,outletbc,&
       depvar_at_node,depvar_at_elem,prq_solution,mesh_dof,mesh_type)
    USE arrays,ONLY: dp,num_elems,num_nodes,elem_nodes,elem_cnct,node_xyz,units,&
         num_units
    USE diagnostics, ONLY: enter_exit

    INTEGER :: mesh_dof
    INTEGER :: depvar_at_elem(0:2,2,num_elems)
    INTEGER :: depvar_at_node(num_nodes,0:2,2)
    REAL(dp) :: prq_solution(mesh_dof,2),inletbc,outletbc,density,grav_vect(3)
    LOGICAL:: ADD
    LOGICAL :: FIX(mesh_dof)
    CHARACTER(len=60) ::bc_type,mesh_type

    ! local variables
    INTEGER :: nonode,np,ne,ny1,nj,np_in
    REAL(dp) :: grav
    CHARACTER(len=60) :: sub_name

    sub_name = 'boundary_conditions'
    CALL enter_exit(sub_name,1)
    IF(.NOT.ADD)THEN
       ! Initial values
       FIX(1:mesh_dof)=.FALSE.
       prq_solution = 0
       ! Fixed boundary conditions
       ! These are inlet BCs, apply to all inlet BCs (there should only be one)
       DO ne=1,num_elems
          !ne=elems(noelem)
          IF (elem_cnct(-1,0,ne) == 0) THEN !Entry element
             IF(BC_TYPE == 'pressure')THEN
                np=elem_nodes(1,ne)
                ny1=depvar_at_node(np,1,1) !for fixed pressure BC
                FIX(ny1)=.TRUE. !set fixed
                prq_solution(ny1,1)=inletbc !Putting BC value into solution array
                np_in=np !inlet node set here, gravity reference to this point
             ELSE IF(BC_TYPE == 'flow')THEN
                ny1=depvar_at_elem(0,1,ne) !fixed
                FIX(ny1)=.TRUE. !set fixed
                prq_solution(ny1,1)=inletbc !Putting BC value into solution array
                np_in=elem_nodes(1,ne)
             ENDIF
          ENDIF
       ENDDO
    ELSE !Add terminal pressure BC for all terminal branches
       IF(mesh_type.EQ.'simple_tree')THEN !NEED TO SET GRAVITY IN THIS CASE
          DO nonode=1,num_units
             np=elem_nodes(2,units(nonode)) !Second node in element = terminal node
             ny1=depvar_at_node(np,1,1) !for fixed pressure BC
             FIX(ny1)=.TRUE. !set fixed
             !! NB// Add gravitational factor in here
             IF(np_in.EQ.0) THEN
                PRINT *,"Warning --> np_in is not set yet, setting to first node as default"
                np_in=1 !Setting to first node as default
             ENDIF
             grav=0.d0
             DO nj=1,3
                grav=grav+density*grav_vect(nj)*9810.d0*(node_xyz(nj,np_in)-node_xyz(nj,np))
             ENDDO
             prq_solution(ny1,1)=outletbc-grav !Putting BC value into solution array
             !print *,"BC----",ny1,outletbc,grav,prq_solution(ny1,1)
          ENDDO
       ELSE
          DO ne=1,num_elems
             !ne=elems(noelem)
             IF (elem_cnct(1,0,ne) == 0) THEN !EXIT ELEMENT
                np=elem_nodes(2,ne)
                ny1=depvar_at_node(np,1,1) !for fixed pressure BC
                FIX(ny1)=.TRUE. !set fixed
                prq_solution(ny1,1)=outletbc !Putting BC value into solution array
                np_in=np !inlet node set here, gravity reference to this point
             ENDIF
          ENDDO

       ENDIF
    ENDIF
    CALL enter_exit(sub_name,2)
  END SUBROUTINE boundary_conditions
  !
  !###################################################################################
  !
  SUBROUTINE calculate_resistance(viscosity,KOUNT)
    USE arrays,ONLY: dp,num_elems,elem_nodes,node_xyz,&
         elem_field
    USE other_consts
    USE indices
    USE diagnostics, ONLY: enter_exit
    REAL(dp):: viscosity
    INTEGER :: KOUNT
    !local variables
    INTEGER :: ne,np1,np2
    REAL(dp) :: resistance,zeta
    CHARACTER(len=60) :: sub_name

    sub_name = 'calculate_resistance'
    CALL enter_exit(sub_name,1)

    !Loop over all elements in model and define resistance for that branch.
    DO ne=1,num_elems
       !ne=elems(noelem)
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       ! element length
       elem_field(ne_length,ne) = DSQRT((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)
       elem_field(ne_radius,ne)=(elem_field(ne_radius_in,ne)+elem_field(ne_radius_out,ne))/2.0_dp
       ! element Poiseuille (laminar) resistance in units of Pa.s.mm-3
       resistance = 8.d0*viscosity*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_radius,ne)**4) !laminar resistance
       ! element turbulent resistance (flow in bifurcating tubes)
       !reynolds=DABS(elem_field(ne_Qdot,ne)*2.d0*density/ &
       !   (PI*elem_field(ne_radius,ne)*viscosity))
       zeta = 1.0_dp!MAX(1.d0,dsqrt(2.d0*elem_field(ne_radius,ne)* &
       !reynolds/elem_field(ne_length,ne))*gamma)
       IF(elem_field(ne_group,ne).EQ.1.0_dp)THEN
          elem_field(ne_resist,ne) = 1000.0_dp !initialises resistance for first iteration
       ELSE
          elem_field(ne_resist,ne) = resistance * zeta
       ENDIF
       !print *,"TESTING RESISTANCE",ne,elem_field(ne_resist,ne),elem_field(ne_radius,ne),elem_ordrs(2,ne)
       !endif
    ENDDO

    CALL enter_exit(sub_name,2)
  END SUBROUTINE calculate_resistance
  !
  !##################################################################
  !
  !*calc_depvar_maps:* calculates the mapping between the nodes and elements and
  !the problem depvar that are needed for matrix setup and solution
  SUBROUTINE calc_depvar_maps(mesh_from_depvar,depvar_at_elem,&
       depvar_totals,depvar_at_node,mesh_dof,num_vars)
    USE arrays,ONLY: num_elems,num_nodes,elem_nodes
    USE diagnostics, ONLY: enter_exit
    CHARACTER(len=60) :: sub_name

    INTEGER :: mesh_dof
    INTEGER :: mesh_from_depvar(0:2,mesh_dof,0:2)
    INTEGER :: depvar_at_elem(0:2,2,num_elems)
    INTEGER :: depvar_at_node(num_nodes,0:2,2)
    INTEGER :: depvar_totals(0:2,2)
    INTEGER :: num_vars
    !     local variables
    INTEGER :: ny_start=0
    INTEGER :: nc,ne,nn,np,nrc,ny

    sub_name = 'calc_depvar_maps'
    CALL enter_exit(sub_name,1)

    depvar_totals = 0
    mesh_from_depvar = 0
    depvar_at_elem = 0
    depvar_at_node = 0

    !nrc = loops from 0,1,2

    !  Set up mapping arrays for current region:
    !  includes nested loop creating depvar_at_node & depvar_at_elem consecutively
    DO nrc=0,2 !row or column no
       ny=0 !depvar tag
       DO nc=1,2 !no of dependent depvar types
          ! if(nrc.NE.0) ny=ny_start !--> resets to ny_start only for nrc=1,2. Maybe should also do for nrc=1??
          ny=ny_start
          DO ne=1,num_elems
             !ne=elems(noelem)
             DO nn=1,2 !Loop over number of nodes in each element: always=2 in 1D element
                np=elem_nodes(nn,ne)
                IF(depvar_at_node(np,nrc,nc).EQ.0)THEN
                   !Checks if this node already been done
                   ny=ny+1
                   IF(ny.GT.mesh_dof) THEN
                      PRINT *,"1. Need to increase mesh_dof! ny=",ny
                      CALL EXIT(1)
                   ENDIF
                   IF(nrc.NE.0) THEN
                      IF(ny.GT.depvar_totals(nrc,nc)) depvar_totals(nrc,nc)=ny
                   ENDIF
                   depvar_at_node(np,nrc,nc)=ny
                   IF(nrc.NE.1.OR.nc.EQ.1) THEN
                      ! don't set mesh_from_depvar for nrc=1(rows), nc<>1(LHS) cases
                      mesh_from_depvar(0,ny,nrc)=1 !mesh dof is nodal based
                      mesh_from_depvar(1,ny,nrc)=np
                      mesh_from_depvar(2,ny,nrc)=nc
                   ENDIF !nrc.NE.1
                ENDIF !depvar_at_node.EQ.0
             ENDDO !nn (np)
             ny=ny+1
             IF(ny.GT.mesh_dof) THEN
                PRINT *,"2. Need to increase mesh_dof! ny=",ny
                CALL EXIT(1)
             ENDIF
             IF(nrc.NE.0) THEN
                IF(ny.GT.depvar_totals(nrc,nc)) depvar_totals(nrc,nc)=ny
             ENDIF
             depvar_at_elem(nrc,nc,ne)=ny
             IF(nrc.NE.1.OR.nc.EQ.1) THEN
                !                 don't set mesh_from_depvar for nrc=1(rows), nc<>1(LHS) cases
                mesh_from_depvar(0,ny,nrc)=2 !mesh dof is element based
                mesh_from_depvar(1,ny,nrc)=nc
                mesh_from_depvar(2,ny,nrc)=ne
             ENDIF
          ENDDO !noelem (ne)
       ENDDO !nc
    ENDDO !nrc
    num_vars=ny
    !print *,"Max ny number=",ny


    CALL enter_exit(sub_name,2)
  END SUBROUTINE calc_depvar_maps
  !
  !##################################################################
  !
  !*tree_resistance:* Calculates the total resistance of a tree (arterial tree only)

  SUBROUTINE tree_resistance(resistance)
    USE indices
    USE arrays,ONLY: dp,num_elems,elem_cnct,elem_field
    USE diagnostics, ONLY: enter_exit
    CHARACTER(len=60) :: sub_name
    !local variables
    REAL(dp), INTENT(out) :: resistance
    REAL(dp) :: invres,elem_res(num_elems)
    INTEGER :: num2,ne,ne2

    sub_name = 'tree_resistance'
    CALL enter_exit(sub_name,1)

    elem_res(1:num_elems)=elem_field(ne_resist,1:num_elems)
    DO ne=num_elems,1,-1
       !ne=elems(num)
       invres=0.0_dp
       DO num2=1,elem_cnct(1,0,ne)
          ne2=elem_cnct(1,num2,ne)
          invres=invres+1.0_dp/elem_res(ne2)
       ENDDO
       IF(elem_cnct(1,0,ne).GT.0)THEN
          elem_res(ne)=elem_res(ne)+1.0_dp/invres
       ENDIF
    ENDDO
    resistance=elem_res(1)

    CALL enter_exit(sub_name,2)
  END SUBROUTINE tree_resistance
  !
  !##################################################################
  !
  !*initialise_solution:* Calculates an estimate for initial solution to a prq problem based on cardiact output and pressure BCs

  SUBROUTINE initialise_solution(pressure_in,pressure_out,cardiac_output,mesh_dof,prq_solution,&
       depvar_at_node,depvar_at_elem,FIX)

    USE indices
    USE arrays,ONLY: dp,num_elems,elem_ordrs,elem_nodes,num_nodes
    USE diagnostics, ONLY: enter_exit
    INTEGER, INTENT(in) :: mesh_dof
    INTEGER,INTENT(in) :: depvar_at_elem(0:2,2,num_elems)
    INTEGER,INTENT(in) :: depvar_at_node(num_nodes,0:2,2)
    REAL(dp), INTENT(in) :: pressure_in, pressure_out,cardiac_output
    REAL(dp) :: prq_solution(mesh_dof,2)
    LOGICAL, INTENT(in) :: FIX(mesh_dof)
    !local variables
    INTEGER :: nn,ne,np,n_depvar
    CHARACTER(len=60) :: sub_name
    sub_name = 'intialise_solution'
    CALL enter_exit(sub_name,1)
    DO ne=1,num_elems
       !ne=elems(noelem)
       DO nn=1,2 !Loop over number of nodes in each element: always=2 in 1D element
          np=elem_nodes(nn,ne) !Node number
          n_depvar=depvar_at_node(np,0,1) !--> This will be a pressure because it's the depvar at a node
          IF(.NOT.FIX(n_depvar))THEN
             prq_solution(n_depvar,1)=(pressure_in+pressure_out)/2.0
          ENDIF
          n_depvar=depvar_at_elem(0,1,ne) !--> This will be a flow because it's the depvar at the element
          IF(.NOT.FIX(n_depvar))THEN
             prq_solution(n_depvar,1)=cardiac_output/(2**(elem_ordrs(1,ne)-1)) !Here you can use the generation number to split flow
          ENDIF
       ENDDO
    ENDDO
    CALL enter_exit(sub_name,2)
  END SUBROUTINE initialise_solution
  !
  !##################################################################
  !
  !*calc_sparse_1d_tree:* Calculates sparsity structure for 1d tree problems

  SUBROUTINE calc_sparse_1dtree(bc_type,density,FIX,grav_vect,mesh_dof,depvar_at_elem,&
       depvar_at_node,NonZeros,MatrixSize,SparseCol,SparseRow,SparseVal,RHS,&
       prq_solution,update_resistance_entries,update_flow_nzz_row)


    USE indices
    USE arrays,ONLY: dp,num_elems,elem_nodes,num_nodes,elems_at_node,elem_field,&
         node_xyz,elem_cnct,elem_ordrs

    USE diagnostics, ONLY: enter_exit

    CHARACTER(len=60) :: bc_type
    REAL(dp), INTENT(in) :: density
    REAL(dp), INTENT(in) :: grav_vect(3)
    INTEGER, INTENT(in) :: mesh_dof
    LOGICAL, INTENT(in) :: FIX(mesh_dof)
    INTEGER,INTENT(in) :: depvar_at_elem(0:2,2,num_elems)
    INTEGER,INTENT(in) :: depvar_at_node(num_nodes,0:2,2)

    INTEGER, INTENT(in) :: NonZeros,MatrixSize
    INTEGER, INTENT(inout) :: SparseCol(NonZeros)
    INTEGER, INTENT(inout) :: SparseRow(MatrixSize+1)
    REAL(dp), INTENT(inout) :: SparseVal(NonZeros)
    REAL(dp), INTENT(inout) :: RHS(MatrixSize)
    REAL(dp), INTENT(inout) :: prq_solution(mesh_dof,2)
    INTEGER, INTENT(inout) :: update_resistance_entries(num_elems)
    INTEGER, INTENT(inout) :: update_flow_nzz_row
    !local variables
    INTEGER :: ne,nn,np,np1,np2,depvar,depvar1,depvar2,depvar3,flow_var,fixed_var_index,offset,nzz,&
         nzz_row,ne2,noelem2,ne3
    LOGICAL :: FlowBalancedNodes(num_nodes)
    LOGICAL :: NodePressureDone(num_nodes)
    LOGICAL :: ElementPressureEquationDone(num_elems)
    LOGICAL :: elem_found,one_node_balanced
    REAL(dp) :: flow_term
    REAL(dp) :: grav
    INTEGER :: nj

    CHARACTER(len=60) :: sub_name
    sub_name = 'calc_sparse_1dtree'
    CALL enter_exit(sub_name,1)
    !Initialise matrices and indices
    SparseCol=0
    SparseRow=1
    SparseVal=0.0_dp
    RHS=0.0_dp

    nzz=1 !position in SparseCol and SparseVal
    nzz_row=1 !position in SparseRow
    FlowBalancedNodes = .FALSE. !.TRUE. for nodes which have had a conservation of flow equation done
    NodePressureDone = .FALSE.  !.TRUE. for nodes which have been processed
    ElementPressureEquationDone = .FALSE.
    offset=0!variable position offset

    DO ne=1,num_elems
       !look at pressure variables at each node
       DO nn=1,2 !2 nodes in 1D element
          np=elem_nodes(nn,ne)
          depvar = depvar_at_node(np,1,1)
          IF((.NOT.NodePressureDone(np)).AND.(.NOT.FIX(depvar)))THEN !check if this node is not fixed and hasn't already been processed (as nodes are shared between elements)
             ne2=0
             IF(nn.EQ.1)THEN !first node of the element
                ne2=ne! use the current element
             ELSEIF(nn.EQ.2)THEN !second node of the element
                IF((bc_type.EQ.'pressure').OR.(.NOT.ElementPressureEquationDone(ne)))THEN !if bc_type is pressure or element pressure equation for the current element hasn't been used
                   ne2=ne! use the current element
                ELSE
                   !look for another element connected to this node with pressure equation that hasn't been used
                   IF (elems_at_node(np,0).GT.1)THEN
                      elem_found=.FALSE.
                      noelem2 = 1
                      DO WHILE ((.NOT.elem_found).AND.(noelem2.LE.elems_at_node(np,0)))
                         ne3=elems_at_node(np,noelem2)
                         IF((ne3.NE.ne).AND.(.NOT.ElementPressureEquationDone(ne3)))THEN
                            ne2 = ne3
                            elem_found=.TRUE.
                         ENDIF
                         noelem2 = noelem2 + 1
                      END DO
                   ENDIF
                ENDIF
             ENDIF
             IF(ne2.GT.0)THEN
                !do the pressure equation for element ne2
                !pressure for node 1 - pressure for node 2 - resistance * flow at element ne2 = 0
                np1=elem_nodes(1,ne2)
                depvar1=depvar_at_node(np1,1,1) !pressure variable for first node
                np2=elem_nodes(2,ne2) !second node
                depvar2=depvar_at_node(np2,1,1) !pressure variable for second node
                depvar3=depvar_at_elem(0,1,ne2) !flow variable for element
                grav=0.d0
                IF(elem_field(ne_group,ne2).EQ.1.0_dp)THEN
                ELSEIF(elem_ordrs(no_gen,ne2).EQ.1)THEN !gravitational head not applied in inlets
                ELSE
                   DO nj=1,3
                      grav=grav+density*grav_vect(nj)*9810.0_dp*(node_xyz(nj,elem_nodes(1,ne2))-node_xyz(nj,elem_nodes(2,ne2)))!rho g L cos theta (Pa)
                   ENDDO
                ENDIF
                IF(FIX(depvar1))THEN !checking if pressure at 1st node is fixed
                   !store known variable - inlet pressure
                   RHS(nzz_row) = -prq_solution(depvar1,1) + grav
                ELSE
                   !unknown variable -pressure for node 1
                   CALL get_variable_offset(depvar1,mesh_dof,FIX,offset)
                   SparseCol(nzz) = depvar1 - offset !variable number
                   SparseVal(nzz)=1.0_dp !variable coefficient
                   nzz=nzz+1 !next column
                   RHS(nzz_row) = grav
                ENDIF
                IF(FIX(depvar2))THEN !checking if pressure at 2nd node is fixed
                   !store known variable - outlet pressure
                   RHS(nzz_row) = prq_solution(depvar2,1) + grav
                ELSE
                   !unknown variable - pressure for node 2
                   CALL get_variable_offset(depvar2,mesh_dof,FIX,offset)
                   SparseCol(nzz) = depvar2 - offset !variable number
                   SparseVal(nzz)=-1.0_dp !variable coefficient
                   nzz=nzz+1 !next column
                ENDIF
                IF(FIX(depvar3))THEN !checking if flow at element ne2 is fixed
                   !store known variable - inlet flow * resistance for element ne
                   RHS(nzz_row) = prq_solution(depvar3,1)*elem_field(ne_resist,ne2)
                   update_flow_nzz_row = nzz_row
                ELSE
                   !unknown flow
                   CALL get_variable_offset(depvar3,mesh_dof,FIX,offset)
                   SparseCol(nzz) = depvar3-offset !variable position in the unknown variable vector
                   SparseVal(nzz)=-elem_field(ne_resist,ne2) !variable coefficient = resistance for element ne2
                   update_resistance_entries(ne2) = nzz
                   nzz=nzz+1 !next column
                ENDIF
                nzz_row=nzz_row+1 !store next row position
                SparseRow(nzz_row)=nzz
                NodePressureDone(np) = .TRUE.
                ElementPressureEquationDone(ne2) = .TRUE.
             ENDIF
          ENDIF
       ENDDO !nn

       !look at flow variable for the element
       flow_var = depvar_at_elem(0,1,ne)
       IF(.NOT.FIX(flow_var))THEN !don't do anything if flow is fixed
          one_node_balanced = .FALSE.
          !check if node 1 or node 2 are unbalanced
          DO nn=1,2 !do flow balance for each element node
             np = elem_nodes(nn,ne)
             IF((elems_at_node(np,0).GT.1).AND.(.NOT.FlowBalancedNodes(np)))THEN !if there is more than one element at a node and the node is not already flow balanced
                IF((bc_type.EQ.'pressure').OR.((bc_type.EQ.'flow').AND.(.NOT.one_node_balanced)))THEN !do just one flow balance equation for bc_type flow
                   !go through each element connected to node np and add the conservation of flow equation for the elements
                   DO noelem2=1,elems_at_node(np,0)
                      ne2=elems_at_node(np,noelem2)
                      depvar=depvar_at_elem(1,1,ne2)
                      flow_term = 0
                      IF(np.EQ.elem_nodes(2,ne2))THEN !end node
                         flow_term = 1.0_dp
                      ELSEIF(np.EQ.elem_nodes(1,ne2))THEN !start node
                         flow_term = -1.0_dp
                      ENDIF
                      IF(FIX(depvar))THEN
                         RHS(nzz_row)=-prq_solution(depvar,1)*flow_term
                      ELSE
                         !populate SparseCol and SparseVal
                         CALL get_variable_offset(depvar,mesh_dof,FIX,offset)
                         SparseCol(nzz) = depvar - offset
                         SparseVal(nzz) = flow_term
                         nzz = nzz + 1
                      ENDIF
                   ENDDO
                   FlowBalancedNodes(np) = .TRUE.
                   nzz_row=nzz_row+1 !store next row position
                   SparseRow(nzz_row)=nzz
                   one_node_balanced = .TRUE.
                ENDIF !checking bc_type
             ENDIF !flow at node np is unbalanced
          ENDDO !nn

          !if flow balancing hasn't been done for any node for element ne and pressure equation hasn't already been done, do the pressure equation for the element
          IF((.NOT.one_node_balanced).AND.(.NOT.ElementPressureEquationDone(ne)))THEN

             !do the pressure equation for element ne
             !pressure for node 1 - pressure for node 2 - resistance * flow at element ne = 0
             np1=elem_nodes(1,ne)
             depvar1=depvar_at_node(np1,1,1) !pressure variable for first node
             np2=elem_nodes(2,ne) !second node
             depvar2=depvar_at_node(np2,1,1) !pressure variable for second node

             !unknown variable -pressure for node 1
             CALL get_variable_offset(depvar1,mesh_dof,FIX,offset)
             SparseCol(nzz) = depvar1 - offset !variable number
             SparseVal(nzz)=1.0_dp !variable coefficient
             nzz=nzz+1 !next column

             IF(FIX(depvar2))THEN !checking if pressure at 2nd node is fixed
                !store known variable - outlet pressure
                RHS(nzz_row) = prq_solution(depvar2,1)
             ELSE
                !unknown variable - pressure for node 2
                CALL get_variable_offset(depvar2,mesh_dof,FIX,offset)
                SparseCol(nzz) = depvar2 - offset !variable number
                SparseVal(nzz)=-1.0_dp !variable coefficient
                nzz=nzz+1 !next column
             ENDIF

             !unknown flow
             CALL get_variable_offset(flow_var,mesh_dof,FIX,offset)
             SparseCol(nzz) = flow_var-offset !variable position in the unknown variable vector
             SparseVal(nzz)=-elem_field(ne_resist,ne) !variable coefficient = resistance for element ne
             update_resistance_entries(ne) = nzz
             nzz=nzz+1 !next column

             nzz_row=nzz_row+1 !store next row position
             SparseRow(nzz_row)=nzz

             ElementPressureEquationDone(ne) = .TRUE.
          ENDIF
       ENDIF
    ENDDO !ne

    CALL enter_exit(sub_name,2)
  END SUBROUTINE calc_sparse_1dtree


  !
  !##################################################################
  !
  !*calc_sparse_size:* Calculates sparsity sizes

  SUBROUTINE calc_sparse_size(mesh_dof,depvar_at_elem,depvar_at_node,FIX,NonZeros,MatrixSize)
    USE indices
    USE arrays,ONLY: num_elems,elem_nodes,num_nodes,elems_at_node
    USE diagnostics, ONLY: enter_exit
    INTEGER, INTENT(in) :: mesh_dof
    INTEGER,INTENT(in) :: depvar_at_elem(0:2,2,num_elems)
    INTEGER,INTENT(in) :: depvar_at_node(num_nodes,0:2,2)
    LOGICAL, INTENT(in) :: FIX(mesh_dof)
    INTEGER :: NonZeros,MatrixSize
    !local variables
    INTEGER :: i,ne,np,fixed_variables, fixed_flows, fixed_pressures
    CHARACTER(len=60) :: sub_name
    sub_name = 'calc_sparse_size'
    CALL enter_exit(sub_name,1)


    fixed_variables = 0
    !count fixed variables
    DO i=1,mesh_dof
       IF(FIX(i))THEN
          fixed_variables = fixed_variables + 1
       ENDIF
    ENDDO
    MatrixSize = mesh_dof - fixed_variables

    !get count of fixed flows
    fixed_flows = 0
    DO ne=1,num_elems
       IF(FIX(depvar_at_elem(1,1,ne)))THEN
          fixed_flows = fixed_flows + 1
       ENDIF
    ENDDO

    fixed_pressures = fixed_variables - fixed_flows

    !count of pressure equations = (number of elements * 3 variables in each equation) - fixed pressures - fixed flows
    NonZeros = num_elems*3 - fixed_pressures - fixed_flows
    !count of conservation of flow equations = sum of elements connected to nodes which have at least 2 connected elements - fixed flows
    DO np=1, num_nodes
       IF(elems_at_node(np,0).GT.1)THEN
          NonZeros = NonZeros + elems_at_node(np,0)
       ENDIF
    ENDDO
    NonZeros = NonZeros - fixed_flows
    CALL enter_exit(sub_name,2)
  END SUBROUTINE calc_sparse_size

  !##################################################################
  !
  !*calc_press_area:* Calculates new radii based on pressure area relnships

  SUBROUTINE calc_press_area(grav_vect,KOUNT,depvar_at_node,prq_solution,&
       mesh_dof,vessel_type,elasticity_parameters,mechanics_parameters)

    USE indices
    USE arrays,ONLY: dp,num_nodes,num_elems,elem_field,elem_nodes,node_xyz
    USE diagnostics, ONLY: enter_exit
    CHARACTER(len=60), INTENT(in) :: vessel_type
    REAL(dp), INTENT(in) :: grav_vect(3)
    INTEGER,INTENT(in) :: KOUNT,mesh_dof
    INTEGER,INTENT(in) :: depvar_at_node(num_nodes,0:2,2)
    REAL(dp),INTENT(in) ::  prq_solution(mesh_dof,2)
    REAL(dp),INTENT(in) :: elasticity_parameters(3),mechanics_parameters(2)

    !local variables
    INTEGER :: nj,np,ne,ny,nn
    REAL(dp) :: h,Ptm,R0,Pblood,Ppl

    CHARACTER(len=60) :: sub_name
    sub_name = 'calc_press_area'
    CALL enter_exit(sub_name,1)
    IF(KOUNT.EQ.1)THEN !store initial, unstressed radius values
       DO  ne=1,num_elems
          elem_field(ne_radius_in0,ne)=elem_field(ne_radius_in,ne)
          elem_field(ne_radius_out0,ne)=elem_field(ne_radius_out,ne)
       ENDDO !elems
    ENDIF

    DO ne=1,num_elems
       DO nn=1,2
          IF(nn.EQ.1) np=elem_nodes(1,ne)
          IF(nn.EQ.2) np=elem_nodes(2,ne)
          ny=depvar_at_node(np,0,1)
          CALL calculate_ppl(np,grav_vect,mechanics_parameters,Ppl)
          Pblood=prq_solution(ny,1) !Pa
          Ptm=Pblood+Ppl     ! Pa
          IF(nn.EQ.1)R0=elem_field(ne_radius_in0,ne)
          IF(nn.EQ.2)R0=elem_field(ne_radius_out0,ne)
          IF(vessel_type.EQ.'elastic_g0_beta')THEN
             IF(Ptm.LT.elasticity_parameters(3).AND.elasticity_parameters(1).GT.0.0_dp)THEN
                IF(nn.EQ.1) elem_field(ne_radius_in,ne)=R0*((Ptm/elasticity_parameters(1))+1.d0)**(1.d0/elasticity_parameters(2))
                IF(nn.EQ.2) elem_field(ne_radius_out,ne)=R0*((Ptm/elasticity_parameters(1))+1.d0)**(1.d0/elasticity_parameters(2))
             ELSEIF(Ptm.LT.0.0_dp.OR.elasticity_parameters(1).LT.TOLERANCE)THEN
                IF(Ptm.LT.0)WRITE(*,*) 'Transmural pressure < zero',ne,Ptm,Pblood,Ppl
                IF(nn.EQ.1) elem_field(ne_radius_in,ne)=R0
                IF(nn.EQ.2) elem_field(ne_radius_out,ne)=R0
             ELSE!ptm>ptmmax
                IF(nn.EQ.1)THEN
                   elem_field(ne_radius_in,ne)=R0*((elasticity_parameters(3)/elasticity_parameters(1))+1.d0) &
                        **(1.d0/elasticity_parameters(2))
                ENDIF
                IF(nn.EQ.2)THEN
                   elem_field(ne_radius_out,ne)=R0*((elasticity_parameters(3)/elasticity_parameters(1))+1.d0) &
                        **(1.d0/elasticity_parameters(2))
                ENDIF
             ENDIF
          ELSEIF(vessel_type.EQ.'elastic_alpha')THEN
             IF(Ptm.LT.elasticity_parameters(2))THEN
                IF(nn.EQ.1) elem_field(ne_radius_in,ne)=R0*((Ptm*elasticity_parameters(1))+1.d0)
                IF(nn.EQ.2) elem_field(ne_radius_out,ne)=R0*((Ptm*elasticity_parameters(1))+1.d0)
             ELSEIF(Ptm.LT.0.0_dp)THEN
                IF(Ptm.LT.0)WRITE(*,*) 'Transmural pressure < zero',ne,Ptm,Pblood,Ppl
                IF(nn.EQ.1) elem_field(ne_radius_in,ne)=R0
                IF(nn.EQ.2) elem_field(ne_radius_out,ne)=R0
             ELSE!ptm>ptmmax
                IF(nn.EQ.1)THEN
                   elem_field(ne_radius_in,ne)=R0*((elasticity_parameters(2)/elasticity_parameters(1))+1.d0)
                ENDIF
                IF(nn.EQ.2)THEN
                   elem_field(ne_radius_out,ne)=R0*((elasticity_parameters(2)/elasticity_parameters(1))+1.d0)
                ENDIF
             ENDIF
          ELSEIF(vessel_type.EQ.'elastic_hooke')THEN
             h=elasticity_parameters(2)*R0
             IF(nn.EQ.1) elem_field(ne_radius_in,ne)=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elasticity_parameters(1)*h)
             IF(nn.EQ.2) elem_field(ne_radius_out,ne)=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elasticity_parameters(1)*h)
          ELSE
             PRINT *, 'no vessel type defined, assuming rigid'
             IF(nn.EQ.1) elem_field(ne_radius_in,ne)=R0
             IF(nn.EQ.2) elem_field(ne_radius_out,ne)=R0
          ENDIF
       ENDDO!nn
    ENDDO!ne




    CALL enter_exit(sub_name,2)
  END SUBROUTINE calc_press_area
  !##############################################################################
  !
  !*map_solution_to_mesh* maps the solution array to appropriate nodal and element fields
  SUBROUTINE map_solution_to_mesh(prq_solution,depvar_at_elem,depvar_at_node,mesh_dof)
    USE indices
    USE arrays,ONLY: dp,num_nodes,num_elems,elem_field,node_field
    USE diagnostics, ONLY: enter_exit

    INTEGER, INTENT(in) :: mesh_dof
    REAL(dp),INTENT(in) ::  prq_solution(mesh_dof,2)
    INTEGER,INTENT(in) :: depvar_at_elem(0:2,2,num_elems)
    INTEGER,INTENT(in) :: depvar_at_node(num_nodes,0:2,2)
    !local variables
    INTEGER :: np,ne,ny


    CHARACTER(len=60) :: sub_name
    sub_name = 'map_solution_to_mesh'
    CALL enter_exit(sub_name,1)
    DO  ne=1,num_elems
       ny=depvar_at_elem(1,1,ne)
       elem_field(ne_Qdot,ne)=prq_solution(ny,1)
    ENDDO !elems
    DO np=1,num_nodes
       ny=depvar_at_node(np,0,1)
       node_field(nj_bv_press,np)=prq_solution(ny,1)
    ENDDO

    CALL enter_exit(sub_name,2)
  END SUBROUTINE map_solution_to_mesh

  !##############################################################################
  !
  !*map_flow_to_terminals* maps the solution array to appropriate nodal and element fields
  SUBROUTINE map_flow_to_terminals
    USE indices
    USE arrays,ONLY: elem_field,node_field,num_units,units,unit_field,elem_nodes
    USE diagnostics, ONLY: enter_exit
    INTEGER :: nu,ne,np
    CHARACTER(len=60) :: sub_name
    sub_name = 'map_flow_to_terminals'
    CALL enter_exit(sub_name,1)

    DO nu=1,num_units
       ne=units(nu)
       np=elem_nodes(2,ne)
       unit_field(nu_perf,nu)=elem_field(ne_Qdot,ne)
       unit_field(nu_blood_press,nu)=node_field(nj_bv_press,np)
    ENDDO

    CALL enter_exit(sub_name,2)
  END SUBROUTINE map_flow_to_terminals
  !
  !
  !
  !*calculate_ppl* calculates pleural pressure at a node
  !
  SUBROUTINE calculate_ppl(np,grav_vect,mechanics_parameters,Ppl)
    USE indices
    USE arrays,ONLY: dp,node_xyz
    USE diagnostics, ONLY: enter_exit
    INTEGER, INTENT(in) :: np
    REAL(dp), INTENT(in) :: mechanics_parameters(2)
    REAL(dp), INTENT(in) :: grav_vect(3)
    REAL(dp), INTENT(out) :: Ppl
    !Local variables
    INTEGER :: nj
    REAL(dp) :: G_PLEURAL,HEIGHT(3)
    CHARACTER(len=60) :: sub_name

    sub_name = 'calculate_ppl'

    CALL enter_exit(sub_name,1)
    G_PLEURAL=0.0_dp    !gravitational force
    DO nj=1,3
       HEIGHT(nj)=node_xyz(nj,np)-node_xyz(nj,1) !ARC - where to put grav reference height?
       G_PLEURAL=G_PLEURAL+mechanics_parameters(2)*grav_vect(nj)*9810.d0*HEIGHT(nj) !kg
    ENDDO
    Ppl=mechanics_parameters(1)-G_PLEURAL !Pa

    CALL enter_exit(sub_name,2)
  END SUBROUTINE calculate_ppl

  !
  !##################################################################
  !
  SUBROUTINE get_variable_offset(depvar,mesh_dof,FIX,offset)
    !*Description*: This subroutine returns the number of fixed variables in the depvar array that came before the input depvar
    INTEGER, INTENT(in) :: depvar, mesh_dof
    LOGICAL, INTENT(in) :: FIX(mesh_dof)
    INTEGER, INTENT(inout) :: offset

    INTEGER :: depvar1

    offset = 0
    DO depvar1 = 1,depvar
       IF(FIX(depvar1))THEN
          offset = offset + 1
       ENDIF
    ENDDO

  END SUBROUTINE get_variable_offset

END MODULE pressure_resistance_flow
