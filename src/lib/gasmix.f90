!> \file
!> \author Merryn Tawhai
!> \brief This module contains code specific to running gas mixing problems
!>
!> \section LICENSE
!>
!>
!> Contributor(s):
!>
!>\Description
!> This module contains code specific to running gas mixing problems

MODULE gasmix
  USE arrays,ONLY: dp
  IMPLICIT NONE
  PRIVATE airway_mesh_deform

  INTEGER,PUBLIC :: inlet_node = 1
  INTEGER,PRIVATE :: NonZeros_unreduced
  REAL(dp),PRIVATE :: ideal_mass,initial_mass
  REAL(dp),PRIVATE :: total_volume_change

  INTEGER,ALLOCATABLE :: sparsity_col(:),reduced_col(:)
  INTEGER,ALLOCATABLE :: sparsity_row(:),reduced_row(:)
  REAL(dp),ALLOCATABLE :: global_K(:),global_M(:),global_AA(:),global_BB(:)
  REAL(dp),ALLOCATABLE :: global_R(:)

CONTAINS

!!!#########################################################################

  SUBROUTINE assemble_gasmix(diffusion_coeff,nonzeros_unreduced)

    USE arrays,ONLY: dp,elem_nodes,num_elems
    USE geometry, ONLY: volume_of_mesh
    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER,INTENT(in) :: nonzeros_unreduced
    REAL(dp),INTENT(in) :: diffusion_coeff

    INTEGER :: i,j,ncol,ne,nentry,nrow
    REAL(dp) :: elem_K(2,2),elem_M(2,2),elem_R(2)
    LOGICAL :: found
    CHARACTER(len=60) :: sub_name

!!!................................................................

    sub_name = 'assemble_gasmix'
    CALL enter_exit(sub_name,1)

    global_K(1:nonzeros_unreduced) = 0.0_dp
    global_M(1:nonzeros_unreduced) = 0.0_dp

    DO ne=1,num_elems
       CALL element_gasmix(ne,elem_K,elem_M,elem_R,diffusion_coeff)
       DO i=1,2
          nrow = elem_nodes(i,ne)
          DO j=1,2
             ncol = elem_nodes(j,ne)
             found=.FALSE.
             nentry = sparsity_row(nrow) ! start check at start of row
             DO WHILE (.NOT.found)
                IF(ncol.EQ.sparsity_col(nentry))THEN
                   found = .TRUE.
                ELSE
                   nentry = nentry+1
                ENDIF
             ENDDO
             global_K(nentry) = global_K(nentry) + elem_K(i,j)
             global_M(nentry) = global_M(nentry) + elem_M(i,j)
          ENDDO !j
       ENDDO !i
    ENDDO !noelem

    CALL enter_exit(sub_name,2)

  END SUBROUTINE assemble_gasmix

!!!################################################################################

  SUBROUTINE calc_mass(nj,nu_field,gas_mass)
    USE arrays,ONLY: dp,elem_cnct,elem_nodes,elem_symmetry,&
         node_field,num_elems,num_nodes,num_units,elem_field,units,unit_field
    USE indices,ONLY: ne_vol,nu_vol
    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER,INTENT(in) :: nj,nu_field
    REAL(dp) :: gas_mass
    !     Local Variables
    INTEGER :: ne,ne0,np1,np2,nunit
    REAL(dp) :: average_conc
    REAL(dp),ALLOCATABLE :: tree_mass(:)
    CHARACTER(len=60):: sub_name

    !...........................................................................

    sub_name = 'calc_mass'
    CALL enter_exit(sub_name,1)

    IF(.NOT.ALLOCATED(tree_mass)) ALLOCATE(tree_mass(num_nodes))

    ! initialise to the mass in each element
    DO ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       average_conc = (node_field(nj,np1)+node_field(nj,np2))/2.0_dp
       tree_mass(ne) = average_conc*elem_field(ne_vol,ne)
    ENDDO

    ! add the mass in each elastic unit to terminal elements
    DO nunit=1,num_units
       ne=units(nunit)
       tree_mass(ne) = tree_mass(ne) + &
            unit_field(nu_vol,nunit)*unit_field(nu_field,nunit)
    ENDDO

    ! sum mass recursively up the tree
    DO ne=num_elems,2,-1 ! not for the stem branch; parent = 0
       ne0=elem_cnct(-1,1,ne)
       tree_mass(ne0) = tree_mass(ne0) + DBLE(elem_symmetry(ne))*tree_mass(ne)
    ENDDO !noelem

    gas_mass = tree_mass(1)

    DEALLOCATE(tree_mass)

    CALL enter_exit(sub_name,2)

  END SUBROUTINE calc_mass

!!!###################################################################

  SUBROUTINE element_gasmix(ne,elem_K,elem_M,elem_R,diffusion_coeff)
    USE arrays,ONLY: dp,elem_field,elem_symmetry
    USE indices,ONLY: ne_a_A,ne_length,ne_radius
    USE other_consts
    IMPLICIT NONE

    INTEGER,INTENT(in) :: ne
    REAL(dp) :: elem_K(2,2),elem_M(2,2),elem_R(2)
    REAL(dp),INTENT(in) :: diffusion_coeff

    !Local variables
    REAL(dp) :: a_A_ratio,inner_area,length,outer_area,radius

    radius = elem_field(ne_radius,ne)
    length = elem_field(ne_length,ne)
    a_A_ratio = elem_field(ne_a_A,ne)
    outer_area=PI*radius**2
    inner_area=outer_area*a_A_ratio

    elem_M(1,1) = outer_area*length/3.0_dp*DBLE(elem_symmetry(ne))
    elem_M(1,2) = outer_area*length/3.0_dp/2.0_dp*DBLE(elem_symmetry(ne))
    elem_M(2,1) = outer_area*length/3.0_dp/2.0_dp
    elem_M(2,2) = outer_area*length/3.0_dp

    elem_K(1,1) = (inner_area*diffusion_coeff/length)*DBLE(elem_symmetry(ne))
    elem_K(1,2) = (-inner_area*diffusion_coeff/length)*DBLE(elem_symmetry(ne))
    elem_K(2,1) = -inner_area*diffusion_coeff/length
    elem_K(2,2) = inner_area*diffusion_coeff/length

    elem_R = 0.0_dp

  END SUBROUTINE element_gasmix

!!!########################################################################

  SUBROUTINE initial_gasmix(initial_concentration,inlet_concentration)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIAL_GASMIX" :: INITIAL_GASMIX

    USE arrays,ONLY: dp,node_field,num_nodes
    USE indices,ONLY: nj_conc1,nu_conc1
    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    REAL(dp),INTENT(in) :: initial_concentration,inlet_concentration

    CHARACTER(len=60):: sub_name

    ! #########################################################################

    sub_name = 'initial_gasmix'
    CALL enter_exit(sub_name,1)

    node_field(nj_conc1,1:num_nodes) = initial_concentration

    ! initialise the 'ideal mass' to the mass of gas initially in model
    CALL calc_mass(nj_conc1,nu_conc1,ideal_mass)

    node_field(nj_conc1,1) = inlet_concentration
    total_volume_change = 0.0_dp ! records the volume change from FRC

    ! allocate the arrays for solving
    IF(.NOT.ALLOCATED(sparsity_col)) ALLOCATE(sparsity_col(1+3*(num_nodes-1)))
    IF(.NOT.ALLOCATED(reduced_col))  ALLOCATE(reduced_col(1+3*(num_nodes-1)))
    IF(.NOT.ALLOCATED(sparsity_row)) ALLOCATE(sparsity_row(num_nodes+1))
    IF(.NOT.ALLOCATED(reduced_row))  ALLOCATE(reduced_row(num_nodes+1))
    IF(.NOT.ALLOCATED(global_K))     ALLOCATE(global_K(1+3*(num_nodes-1)))
    IF(.NOT.ALLOCATED(global_M))     ALLOCATE(global_M(1+3*(num_nodes-1)))
    IF(.NOT.ALLOCATED(global_AA))    ALLOCATE(global_AA(1+3*(num_nodes-1)))
    IF(.NOT.ALLOCATED(global_BB))    ALLOCATE(global_BB(num_nodes))
    IF(.NOT.ALLOCATED(global_R))     ALLOCATE(global_R(num_nodes))

    ! calculate the sparsity pattern for unreduced system
    CALL sparse_gasmix

    CALL enter_exit(sub_name,2)

  END SUBROUTINE initial_gasmix

!!!##########################################################################

  SUBROUTINE reduce_gasmix(MatrixSize,NonZeros,noffset_entry,noffset_row,&
       inspiration)
    USE arrays,ONLY: num_nodes
    IMPLICIT NONE

    INTEGER :: MatrixSize,NonZeros,&
         noffset_entry,noffset_row
    LOGICAL,INTENT(in) :: inspiration

    INTEGER :: i

    IF(inspiration)THEN !remove first row and column (note: also for breath-hold)

       DO i=1,num_nodes ! one more than # of rows
          reduced_row(i) = sparsity_row(i+1)-3
       ENDDO
       NonZeros = NonZeros_unreduced - 3
       DO i=1,NonZeros
          reduced_col(i) = sparsity_col(i+3)-1
       ENDDO
       reduced_row(1)=1
       MatrixSize = num_nodes - 1
       noffset_entry = 3
       noffset_row = 1

    ELSE !expiration

       DO i=1,num_nodes+1
          reduced_row(i) = sparsity_row(i)
       ENDDO
       NonZeros = NonZeros_unreduced
       DO i=1,NonZeros
          reduced_col(i) = sparsity_col(i)
       ENDDO
       reduced_row(1)=1
       MatrixSize = num_nodes
       noffset_entry = 0
       noffset_row = 0

    ENDIF

  END SUBROUTINE reduce_gasmix

!!!#######################################################################

  SUBROUTINE solve_gasmix(fileid,inr_itr_max,out_itr_max,diffusion_coeff,&
       dt,initial_volume,inlet_concentration,inlet_flow,solve_tolerance,time_end,&
       time_start,inspiration)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SOLVE_GASMIX" :: SOLVE_GASMIX

!!! Assemble matrices for 1D inert gas mixing equation, and solve. The sparsity
!!! structure is first calculated. At each time step the mesh size is changed
!!! incrementally. The equations are solved using a Lagrange-Galerkin method which separates
!!! solution of the advective and diffusive components of the equation. The advective
!!! component is solved using method of characteristics (in this case a simple tracking
!!! of material point movement); and this is followed by solution of the diffusion equation
!!! with the advective solution as an initial condition. See Tawhai, M.H., PhD thesis for
!!! further explanation.

    USE arrays,ONLY: dp,node_field,num_nodes
    USE exports,ONLY: export_node_field
    USE indices,ONLY: nj_conc1,nu_conc1
    USE other_consts
    USE diagnostics, ONLY: enter_exit
    USE geometry, ONLY: volume_of_mesh
    USE solve,ONLY: pmgmres_ilu_cr
    IMPLICIT NONE

    INTEGER,INTENT(in) :: fileid,inr_itr_max,out_itr_max
    REAL(dp),INTENT(in) :: diffusion_coeff,dt,initial_volume,&
         inlet_concentration,inlet_flow,solve_tolerance,time_end,time_start
    LOGICAL,INTENT(in) :: inspiration

    ! Local variables
    REAL(dp),ALLOCATABLE :: solution(:)

    INTEGER :: MatrixSize,nonzeros,ncol,nentry, &
         noffset_entry,noffset_row,np,nrow,nrow_BB,SolverFlag
    REAL(dp) :: AA,BB,current_mass,current_volume, &
         mass_error,theta,time,volume_error,volume_tree,zero_tol=1.0e-8_dp,&
         mass0,mass1,mass_error_deform,mass_error_track,mass_error_solve
    LOGICAL :: carryon
    CHARACTER(len=60) :: sub_name

    ! #############################################################################

    sub_name = 'solve_gasmix'
    CALL enter_exit(sub_name,1)

    ! allocatable array to store the current solution
    IF(.NOT.ALLOCATED(solution)) ALLOCATE(solution(num_nodes))

    ! the weighting for matrices in the reduced system: A = M+K*dt*theta; B = -K*c^(n)*dt
    ! can change this such that fully implicit or fully explicit, however 2/3 generally most stable
    theta=2.0_dp/3.0_dp

    ! get the sparsity arrays for the reduced system. uses compressed row format.
    CALL reduce_gasmix(MatrixSize,nonzeros,noffset_entry,noffset_row,inspiration)

    time = time_start ! initialise the time

    carryon = .TRUE. ! logical for whether solution continues

    ! main time-stepping loop:  time-stepping continues while 'carryon' is true
    DO WHILE (carryon) !

       time = time + dt ! increment time

       IF(time_end + 0.5_dp*dt - time .GT. zero_tol)THEN !continue

          CALL calc_mass(nj_conc1,nu_conc1,mass0) ! calculate model mass, for error check
          CALL airway_mesh_deform(dt,initial_volume,inlet_flow) ! change model size by dV
          CALL calc_mass(nj_conc1,nu_conc1,mass1) ! calculate model mass, for error check
          mass_error_deform = 100 * (mass1-mass0)/mass0 ! error from deforming
          CALL update_unit_mass(dt,inlet_concentration,inlet_flow) ! track gas into/out units

          CALL track_back(dt,inlet_concentration,inlet_flow,inspiration) ! track gas along branches

          CALL calc_mass(nj_conc1,nu_conc1,mass0) ! calculate model mass, for error check
          mass1 = mass1 + inlet_flow*dt*node_field(nj_conc1,1) ! expected mass following tracking
          mass_error_track = 100 * (mass0-mass1)/mass1 ! error from tracking units and branches

          ! assemble the element matrices. Element matrix calculation can be done directly
          ! (based on assumption of interpolation functions) or using Gaussian interpolation.
          CALL assemble_gasmix(diffusion_coeff,nonzeros_unreduced)

          ! initialise the values in the solution matrices
          global_AA(1:nonzeros) = 0.0_dp ! equivalent to M in Tawhai thesis
          global_BB(1:num_nodes) = 0.0_dp ! equivalent to K in Tawhai thesis
          global_R(1:num_nodes) = 0.0_dp

          ! Assemble the reduced system of matrices
          DO np=1,num_nodes ! Loop over rows of unreduced system
             nrow = np ! conveniently true for the way we set up our models
             ! different boundary conditions are applied during inspiration and
             ! expiration: Dirichlet at model entry during inspiration (concentration
             ! = inlet_concentration), and Neumann at model entry during expiration (dcdt = 0)
             IF(.NOT.inspiration)THEN
                BB = global_R(nrow)         !get reduced R.H.S.vector
                DO nentry = sparsity_row(nrow),sparsity_row(nrow+1)-1  !each row entry
                   ncol = sparsity_col(nentry)
                   BB = BB-global_K(nentry)*node_field(nj_conc1,ncol)*dt ! -K*c^(n)*dt
                   AA = global_M(nentry) + dt*theta*global_K(nentry) ! M+K*dt*theta
                   global_AA(nentry) = global_AA(nentry) + AA
                ENDDO
                global_BB(nrow) =  global_BB(nrow) + BB
             ELSEIF(inspiration.AND.np.NE.1)THEN !not first row
                BB = global_R(nrow)         !get reduced R.H.S.vector
                DO nentry = sparsity_row(nrow),sparsity_row(nrow+1)-1  !each row entry
                   ncol = sparsity_col(nentry)
                   BB = BB-global_K(nentry)*node_field(nj_conc1,ncol)*dt
                   AA = global_M(nentry) + dt*theta*global_K(nentry) !M+K*dt*theta
                   IF(ncol.NE.1)THEN ! not first column
                      global_AA(nentry-noffset_entry) = &
                           global_AA(nentry-noffset_entry) + AA
                   ENDIF
                ENDDO
                global_BB(nrow-noffset_row) = &
                     global_BB(nrow-noffset_row) + BB
             ENDIF
          ENDDO

          solution(1:num_nodes) = 0.0_dp

          ! Call a solver to solve the system of reduced equations.
          ! Here we use an iterative solver (GMRES == Generalised Minimal
          ! RESidual method). The solver requires the solution matrices to
          ! be represented in compressed row format.
          CALL pmgmres_ilu_cr (MatrixSize,NonZeros,reduced_row,&
               reduced_col,global_AA,solution,global_BB,&
               out_itr_max,inr_itr_max,solve_tolerance,solve_tolerance,SolverFlag)

          ! transfer the solver solution (in 'Solution') to the node field array
          DO np = 1,num_nodes
             IF(.NOT.inspiration.OR.(inspiration.AND.np.GT.1))THEN
                IF(inspiration)THEN
                   nrow_BB=np-1
                ELSE
                   nrow_BB=np
                ENDIF
                node_field(nj_conc1,np) = node_field(nj_conc1,np) &
                     + Solution(nrow_BB) !c^(n+1)=c^(n)+dc
                node_field(nj_conc1,np) = MAX(0.d0,node_field(nj_conc1,np))
                node_field(nj_conc1,np) = MIN(1.d0,node_field(nj_conc1,np))
             ENDIF
          ENDDO

          !          if(.not.inspiration)then
          !             call smooth_expiration(nj_conc1)
          !          endif
          ! estimate the volume and mass errors
          CALL calc_mass(nj_conc1,nu_conc1,current_mass)
          mass_error_solve = 100 * (current_mass-mass0)/mass0

          CALL volume_of_mesh(current_volume,volume_tree)
          ideal_mass = ideal_mass + inlet_flow*dt*node_field(nj_conc1,1)
          volume_error = 1.0e+2_dp*(current_volume - (initial_volume +  &
               total_volume_change))/(initial_volume + total_volume_change)

          IF(ideal_mass.GT.0.0_dp)THEN
             mass_error = 1.0e+2_dp*(current_mass - ideal_mass)/ideal_mass
          ELSE
             mass_error = 0.0_dp
          ENDIF
          ! output current solution to screen
          WRITE(*,'(F10.3,7(F10.2),F12.3)') &
               time,inlet_flow/1.0e+6_dp,total_volume_change/1.0e+6_dp, &
               current_volume/1.0e+6_dp,volume_error,mass_error_deform,&
               mass_error_track,mass_error_solve,node_field(nj_conc1,1)
          IF(.NOT.inspiration)THEN
             WRITE(fileid,'(F10.3,3(D14.5),4(F10.2),D14.5)') &
                  time,inlet_flow,total_volume_change, &
                  current_volume,volume_error,mass_error_deform,&
                  mass_error_track,mass_error_solve,node_field(nj_conc1,1)
          ENDIF
       ELSE
          carryon = .FALSE.
       ENDIF
    ENDDO !carryon

    CALL enter_exit(sub_name,2)

  END SUBROUTINE solve_gasmix

!!!########################################################################

  SUBROUTINE sparse_gasmix
    USE arrays,ONLY: elem_cnct,elem_nodes,num_elems

    IMPLICIT NONE

    INTEGER :: n_unreduced,i,ncol,ne,ne2,np1,np2,nrow

    sparsity_row(1) = 1
    n_unreduced = 1

    DO ne=1,num_elems ! note using local numbering
       IF(elem_cnct(-1,0,ne).EQ.0)THEN !at the inlet
          np1=elem_nodes(1,ne) ! start node
          nrow=np1
          DO i=1,2
             np2=elem_nodes(i,ne)
             ncol=np2
             sparsity_col(n_unreduced)=ncol
             n_unreduced=n_unreduced+1
          ENDDO
          sparsity_row(nrow+1)=n_unreduced
       ENDIF

       np1=elem_nodes(2,ne) !end node
       nrow=np1
       DO i=1,2
          np2=elem_nodes(i,ne)
          ncol=np2
          sparsity_col(n_unreduced)=ncol
          n_unreduced=n_unreduced+1
       ENDDO
       DO i=1,elem_cnct(1,0,ne) ! for each child branch
          ne2=elem_cnct(1,i,ne)
          np2=elem_nodes(2,ne2)
          ncol=np2
          sparsity_col(n_unreduced)=ncol
          n_unreduced=n_unreduced+1
       ENDDO
       sparsity_row(nrow+1)=n_unreduced
    ENDDO !noelem

    NonZeros_unreduced = n_unreduced - 1

  END SUBROUTINE sparse_gasmix

!!!####################################################################

  SUBROUTINE update_unit_mass(dt,inlet_concentration,inlet_flow)
    USE arrays,ONLY: dp,elem_field,elem_nodes,num_units,units,unit_field
    USE indices,ONLY: ne_Vdot,nu_conc1,nu_vol
    USE geometry, ONLY: volume_of_mesh
    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    REAL(dp),INTENT(in) :: dt,inlet_concentration,inlet_flow

    ! Local parameters
    INTEGER :: ne,np,nunit
    REAL(dp) :: concentration,mass_change,&
         unit_mass,volume_change,sum_mass_change
    CHARACTER(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_unit_mass'
    CALL enter_exit(sub_name,1)

    ! calculate the change in volume of elastic units
    ! note that conducting airways are assumed rigid
    sum_mass_change=0.0_dp
    DO nunit=1,num_units
       ne = units(nunit) ! local element number
       np = elem_nodes(2,ne) ! end node
       unit_mass = unit_field(nu_vol,nunit)*unit_field(nu_conc1,nunit)
       volume_change = elem_field(ne_Vdot,ne)*inlet_flow*dt
       IF(inlet_flow.GE.0.0_dp)THEN !filling units
          CALL track_to_location(ne,concentration,mass_change,&
               dt,inlet_concentration,inlet_flow)
       ELSE ! emptying units
          mass_change = unit_field(nu_conc1,nunit)*volume_change
       ENDIF
       unit_mass = unit_mass + mass_change
       unit_field(nu_conc1,nunit) = unit_mass/unit_field(nu_vol,nunit)
       sum_mass_change = sum_mass_change+mass_change
    ENDDO

    CALL enter_exit(sub_name,2)

  END SUBROUTINE update_unit_mass

!!!####################################################################

  SUBROUTINE airway_mesh_deform(dt,initial_volume,inlet_flow)

!!! changes the size of elastic units and airways, based on pre-computed
!!! fields for flow into the units (ne_Vdot) and element volume change (ne_dvdt).
!!! The concentration in the units is adjusted to make sure that mass is conserved
!!! when the unit changes volume. If a unit changes in volume, then both the radius
!!! and length scale as the cube root of volume change. For alveolated airways (where
!!! a/A is < 1) the outer radius (indexed as nj_radius) is scaled and a/A unchanged.

    USE arrays,ONLY: dp,elem_field,elem_nodes,&
         num_elems,num_units,units,unit_field
    USE indices,ONLY: ne_dvdt,ne_Vdot,ne_length,ne_radius,&
         ne_vol,nu_conc1,nu_vol
    USE geometry, ONLY: volume_of_mesh
    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    REAL(dp),INTENT(in) :: dt,initial_volume,inlet_flow

    ! Local parameters
    INTEGER :: ne,np,nunit
    REAL(dp) :: current_volume,ratio,tree_volume,unit_mass,volume_change
    CHARACTER(len=60) :: sub_name

    ! #################################################################

    sub_name = 'airway_mesh_deform'
    CALL enter_exit(sub_name,1)

!!! calculate the change in volume of elastic units:
    ! mass(unit) = volume(unit) * concentration(unit)
    ! dv(unit) = flow(unit) * dt
    ! volume(new) = volume(old) + dv
    ! concentration(new) = mass/volume(new)
    DO nunit=1,num_units !for each of the elastic units
       ne = units(nunit) ! local element number
       np = elem_nodes(2,ne) ! end node, attaches to unit
       unit_mass = unit_field(nu_vol,nunit)*unit_field(nu_conc1,nunit) ! m = v(unit)*c
       volume_change = elem_field(ne_Vdot,ne)*inlet_flow*dt ! dv = q(unit)*dt
       unit_field(nu_vol,nunit) = unit_field(nu_vol,nunit) + volume_change ! v(new)
       ! important: adjust the concentration to maintain mass conservation
       unit_field(nu_conc1,nunit) = unit_mass/unit_field(nu_vol,nunit)
    ENDDO

!!! calculate the inflation/deflation of multi-branching acini (or other branches)
    ! dv(element) = dvdt(element) * total_volume_change
    ! ratio_of_volumes = (volume(element)+dv(element))/volume(element)
    ! volume(new) = volume(old) * ratio_of_volumes
    ! length(new) = length(old) * cuberoot(ratio_of_volumes)
    ! radius(new) = radius(old) * cuberoot(ratio_of_volumes)
    DO ne=1,num_elems
       volume_change = elem_field(ne_dvdt,ne)*inlet_flow*dt
       ratio = (volume_change+elem_field(ne_vol,ne))/elem_field(ne_vol,ne)
       elem_field(ne_vol,ne)=elem_field(ne_vol,ne)*ratio
       elem_field(ne_length,ne)=elem_field(ne_length,ne)*(ratio)**(1/3)
       elem_field(ne_radius,ne)=elem_field(ne_radius,ne)*(ratio)**(1/3)
    ENDDO

!!! calculate the total model volume
    CALL volume_of_mesh(current_volume,tree_volume)

    total_volume_change = current_volume-initial_volume

    CALL enter_exit(sub_name,2)

  END SUBROUTINE airway_mesh_deform

!!! ######################################################################

  SUBROUTINE smooth_expiration(nj_field)
    USE arrays,ONLY: elem_cnct,elem_ordrs,num_elems
    IMPLICIT NONE

    INTEGER,INTENT(in) :: nj_field

    INTEGER :: ne,ne_next,ngen,ngen_next,num_smooth
    INTEGER,ALLOCATABLE :: nsmooth_list(:)
    LOGICAL,ALLOCATABLE :: smoothed(:)

    ALLOCATE(nsmooth_list(num_elems))
    ALLOCATE(smoothed(num_elems))
    smoothed(1:num_elems)=.FALSE.

    DO ne=1,num_elems
       ngen=elem_ordrs(1,ne)
       IF(ngen.GT.1)THEN
          IF(elem_cnct(1,0,ne).EQ.1)THEN !do only for discretised branch
             IF(.NOT.smoothed(ne))THEN  !hasn't been done as part of branch
                ne_next=elem_cnct(1,1,ne) !element # of adjacent element
                ngen_next=elem_ordrs(1,ne_next) !generation of adjacent element
                IF(ngen.EQ.ngen_next)THEN !same generation, therefore same branch
                   num_smooth = 1
                   nsmooth_list(1) = ne
                   DO WHILE(ne_next.NE.0)
                      num_smooth = num_smooth+1
                      nsmooth_list(num_smooth)=ne_next
                      smoothed(ne_next)=.TRUE.
                      IF(elem_cnct(1,0,ne_next).EQ.1)THEN
                         ne_next=elem_cnct(1,1,ne_next)
                         ngen_next=elem_ordrs(1,ne_next)
                         IF(ngen.NE.ngen_next)  ne_next=0
                      ELSE
                         ne_next=0
                      ENDIF
                   ENDDO     !WHILE
                   CALL smooth_expiration_linear(nj_field,nsmooth_list,num_smooth)
                ENDIF          !gen
             ENDIF             !SMOOTHED
          ELSE                 !(NXI=0 or 2)
             smoothed(ne)=.TRUE.
          ENDIF                !NXI
       ENDIF
    ENDDO !noelem(ne)

    DEALLOCATE(nsmooth_list)
    DEALLOCATE(smoothed)

  END SUBROUTINE smooth_expiration

!!! ######################################################################

  SUBROUTINE smooth_expiration_linear(nj_field,nsmooth_list,num_smooth)
    USE arrays,ONLY: dp,elem_field,elem_nodes,node_field,num_elems
    USE indices,ONLY: ne_length
    IMPLICIT NONE

!!! Parameter List
    INTEGER,INTENT(in) :: nj_field,num_smooth,nsmooth_list(*)

!!! Local Variables
    INTEGER :: n,ne,np
    REAL(dp) :: conc_end,conc_start,intercept,slope
    REAL(dp),ALLOCATABLE :: length(:)

    ALLOCATE(length(num_elems))

    !for first node
    ne=nsmooth_list(1)
    np=elem_nodes(1,ne)
    conc_start=node_field(nj_field,np)

    !for last node
    ne=nsmooth_list(num_smooth)
    np=elem_nodes(2,ne)
    conc_end=node_field(nj_field,np)

    DO n=1,num_smooth
       ne=nsmooth_list(n)
       IF(n.EQ.1)THEN
          length(n)=elem_field(ne_length,ne)
       ELSE
          length(n)=length(n-1)+elem_field(ne_length,ne)
       ENDIF
    ENDDO
    !calculate linear equation
    intercept=conc_start
    slope=(conc_end-intercept)/length(num_smooth)

    !update solution values for the branch
    ne=nsmooth_list(1)
    np=elem_nodes(1,ne)

    intercept=node_field(nj_field,np)
    ne=nsmooth_list(num_smooth)
    np=elem_nodes(2,ne)

    slope=(node_field(nj_field,np)-intercept)/length(num_smooth)

    node_field(nj_field,np)=intercept
    DO n=1,num_smooth
       ne=nsmooth_list(n)
       np=elem_nodes(2,ne)
       node_field(nj_field,np)=intercept+slope*length(n)
    ENDDO

    DEALLOCATE(length)

  END SUBROUTINE smooth_expiration_linear

!!! ##################################################################

  SUBROUTINE transfer_flow_vol_from_units()
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_TRANSFER_FLOW_VOL_FROM_UNITS" :: TRANSFER_FLOW_VOL_FROM_UNITS

    USE arrays,ONLY: dp,elem_cnct,elem_field,elem_symmetry,&
         elem_units_below,expansile,units,num_elems,num_units,&
         unit_field
    USE indices,ONLY: ne_dvdt,ne_Vdot,ne_length,ne_radius,ne_vol,nu_vol
    IMPLICIT NONE

!!! Parameters

!!! Local variables
    INTEGER :: nchild,ne_child,ne_parent,ne_stem,&
         nparent,nunit,num_list,num_list_temp,num_list_total
    INTEGER,ALLOCATABLE :: elem_list(:),elem_list_temp(:),&
         elem_list_total(:)
    REAL(dp),ALLOCATABLE :: vol_below(:)
    REAL(dp) :: ratio


    ALLOCATE(vol_below(num_elems))
    vol_below(1:num_elems) = elem_field(ne_vol,1:num_elems)
    elem_field(ne_dvdt,1:num_elems) = 0.0_dp !rate of change of element volume
    expansile(1:num_elems) = .FALSE.

    ALLOCATE(elem_list(num_elems))
    ALLOCATE(elem_list_temp(num_elems))
    ALLOCATE(elem_list_total(num_elems))

    DO nunit = 1,num_units
       ne_stem=units(nunit) ! the element number supplying a unit
       num_list=1
       elem_list(1)=ne_stem
       elem_list(2:elem_units_below(ne_stem)) = 0

       num_list_temp=0
       elem_list_temp(1:elem_units_below(ne_stem)) = 0

       num_list_total = 0
       elem_list_total = 0 !record all elements below the unit

       DO WHILE(num_list.GT.0)
          DO nparent=1,num_list ! for each 'parent' element in the list
             ne_parent=elem_list(nparent) ! the 'parent' element number
             DO nchild=1,elem_cnct(1,0,ne_parent) ! for each child branch of 'parent'
                ne_child = elem_cnct(1,nchild,ne_parent)
                num_list_temp=num_list_temp+1
                elem_list_temp(num_list_temp)=ne_child
                num_list_total=num_list_total+1
                elem_list_total(num_list_total)=ne_child
                expansile(ne_child)=.TRUE.
             ENDDO !nchild
          ENDDO !nparent
          num_list = num_list_temp
          num_list_temp = 0
          elem_list(1:num_list) =elem_list_temp(1:num_list)
          elem_list_temp(1:elem_units_below(ne_stem)) = 0
       ENDDO !while

       ! calculate the volume below each branch
       DO nchild = num_list_total,1,-1
          ne_child=elem_list_total(nchild)
          ne_parent=elem_cnct(-1,1,ne_child)
          vol_below(ne_parent) = vol_below(ne_parent) &
               + DBLE(elem_symmetry(ne_child))*vol_below(ne_child)
       ENDDO !nchild

!!! scale the elements so same total size as original unit size, and set the
!!! rate of volume change for an element to be its volume/unit volume*stem flow
       ratio = unit_field(nu_vol,nunit)/(vol_below(ne_stem)&
            -elem_field(ne_vol,ne_stem))
       DO nchild=1,num_list_total
          ne_child=elem_list_total(nchild)
          elem_field(ne_vol,ne_child)=elem_field(ne_vol,ne_child)*ratio
          elem_field(ne_length,ne_child)=elem_field(ne_length,ne_child)*(ratio)**(1/3)
          elem_field(ne_radius,ne_child)=elem_field(ne_radius,ne_child)*(ratio)**(1/3)
          elem_field(ne_dvdt,ne_child)=elem_field(ne_Vdot,ne_stem)*&
               elem_field(ne_vol,ne_child)/&
               (unit_field(nu_vol,nunit)-elem_field(ne_vol,ne_stem))
       ENDDO
!!! calculate the flow field. flow(in) = flow(out) + DV
       DO nparent=num_list_total,1,-1
          ne_parent=elem_list_total(nparent)
          elem_field(ne_Vdot,ne_parent)=elem_field(ne_dvdt,ne_parent)
          DO nchild=1,elem_cnct(1,0,ne_parent)
             ne_child=elem_cnct(1,nchild,ne_parent)
             elem_field(ne_Vdot,ne_parent)=elem_field(ne_Vdot,ne_parent)+&
                  elem_field(ne_Vdot,ne_child)*elem_symmetry(ne_child)
          ENDDO !nchild
       ENDDO !nparent
    ENDDO !nunit

    num_units = 0

    DEALLOCATE(elem_list)
    DEALLOCATE(elem_list_temp)
    DEALLOCATE(elem_list_total)
    DEALLOCATE(vol_below)

  END SUBROUTINE transfer_flow_vol_from_units

!!! ##################################################################

  SUBROUTINE track_back(dt,inlet_concentration,inlet_flow,inspiration)
    USE arrays,ONLY: dp,elem_cnct,elem_field,elem_nodes,elem_symmetry,&
         elems_at_node,node_field,num_nodes,num_units,unit_field,units
    USE indices,ONLY: ne_dvdt,ne_Vdot,ne_vol,nj_conc1,nu_conc1
    IMPLICIT NONE

!!! Parameters
    REAL(dp),INTENT(in) :: dt,inlet_concentration,inlet_flow
    LOGICAL,INTENT(in) :: inspiration

!!! Local variables
    INTEGER :: elems_below(1024),i,j,ne,ne_parent,nextra,np,np1,np2,nunit,num_to_check
    REAL(dp) :: branch_fraction(1024),branch_time(1024),concentration,cumulative_mass,&
         flow_fraction,flow_parent,local_xi,time_through_element,total_time
    REAL(dp),ALLOCATABLE :: initial_conc(:)

    ALLOCATE(initial_conc(num_nodes))
    initial_conc = 0.0_dp

    IF(inspiration)THEN
       DO np=2,num_nodes
          ne=elems_at_node(np,1) ! first element; will be most proximal branch
          CALL track_to_location(ne,concentration,cumulative_mass,&
               dt,inlet_concentration,inlet_flow)
          initial_conc(np) = concentration
       ENDDO !np
       initial_conc(1) = inlet_concentration

    ELSE !tracking for expiration
       ! set the concentration at terminal branch to the same as in elastic unit (if attached)
       DO nunit=1,num_units
          ne=units(nunit)
          np=elem_nodes(2,ne)
          initial_conc(np) = unit_field(nu_conc1,nunit)
          node_field(nj_conc1,np) = unit_field(nu_conc1,nunit)
       ENDDO
       DO np=1,num_nodes
          IF(np.EQ.1)THEN
             elems_below(1) = elems_at_node(np,1)
             branch_time(1) = 0.0_dp
             branch_fraction(1) = 1.0_dp
             num_to_check = 1
             ne_parent = elems_at_node(np,1)
          ELSE
             IF(elems_at_node(np,0).EQ.1)THEN  !only one adjoining element, so a terminal node
                num_to_check = 0
                initial_conc(np) = node_field(nj_conc1,np)
             ELSE
                ne_parent = elems_at_node(np,1)
                DO j=2,elems_at_node(np,0)
                   elems_below(j-1) = elems_at_node(np,j) ! check each element
                   branch_time(j-1) = 0.0_dp
                   branch_fraction(j-1) = 1.0_dp
                ENDDO
                num_to_check = elems_at_node(np,0)-1
             ENDIF
          ENDIF
          DO WHILE(num_to_check.GT.0) ! while still some to check
             nextra=0            !records the addition of elements to check
             DO j=1,num_to_check
                ne=elems_below(j)
                IF(np.GT.1) ne_parent = elem_cnct(-1,1,ne)
                total_time=branch_time(j)
                time_through_element = dabs(elem_field(ne_vol,ne)/&
                     (inlet_flow*elem_field(ne_Vdot,ne)))
                IF(total_time+time_through_element.GE.dt)THEN
                   !     A material point is within this element
                   local_xi = (dt-total_time)/time_through_element
                   np1 = elem_nodes(1,ne)
                   np2 = elem_nodes(2,ne)
                   concentration = node_field(nj_conc1,np1)*(1.0_dp-local_xi)&
                        +node_field(nj_conc1,np2)*local_xi
!!! this isn't quite right: needs to take into account the rate of change of parent branch size
!!! also should do this when estimating time to pass through a branch.
                   flow_parent = elem_field(ne_Vdot,ne_parent) - elem_field(ne_dvdt,ne_parent)
                   flow_fraction = branch_fraction(j)*(elem_field(ne_Vdot,ne)*&
                        elem_symmetry(ne))/flow_parent
                   initial_conc(np) = initial_conc(np) + concentration*flow_fraction
                ELSE
                   !     need to check the next elements
                   IF(elem_cnct(1,0,ne).EQ.0)THEN !set to terminal node concentration
                      np2 = elem_nodes(2,ne)
                      concentration = node_field(nj_conc1,np2)
                      flow_parent = elem_field(ne_Vdot,ne_parent) - elem_field(ne_dvdt,ne_parent)
                      flow_fraction = branch_fraction(j)*(elem_field(ne_Vdot,ne)*&
                           elem_symmetry(ne))/flow_parent
                      initial_conc(np) = initial_conc(np) + &
                           concentration*flow_fraction
                   ELSE
                      DO i=1,elem_cnct(1,0,ne) ! for each child element
                         nextra=nextra+1
                         elems_below(num_to_check+nextra) = elem_cnct(1,i,ne)
                         branch_time(num_to_check+nextra) = &
                              total_time+time_through_element
                         flow_parent = elem_field(ne_Vdot,ne_parent) - elem_field(ne_dvdt,ne_parent)
                         branch_fraction(num_to_check+nextra) = branch_fraction(j)*&
                              (elem_field(ne_Vdot,ne)*elem_symmetry(ne))/flow_parent
                      ENDDO   ! i
                   ENDIF      ! terminal or not
                ENDIF         ! exceed the time criteria?
             ENDDO               ! j

             IF(nextra.GT.64)THEN
                WRITE(*,*) 'Need to decrease the time step'
             ENDIF
             DO j=1,nextra
                elems_below(j) = elems_below(num_to_check+j)
                branch_time(j) = branch_time(num_to_check+j)
                branch_fraction(j)= branch_fraction(num_to_check+j)
             ENDDO               !j, for nextra
             num_to_check=nextra
          ENDDO                  ! WHILE num_to_check.GT.0
       ENDDO                     !np
    ENDIF !inspiration

    node_field(nj_conc1,1:num_nodes) = initial_conc(1:num_nodes)

    DEALLOCATE(initial_conc)

  END SUBROUTINE track_back

!!! ##################################################################

  SUBROUTINE track_to_location(ne,concentration,cumulative_mass,&
       dt,inlet_concentration,inlet_flow)

    USE arrays,ONLY: dp,elem_cnct,elem_field,elem_nodes,node_field
    USE indices,ONLY: ne_Vdot,ne_vol,nj_conc1
    IMPLICIT NONE

!!! Parameters
    INTEGER,INTENT(in) :: ne
    REAL(dp),INTENT(in) :: dt,inlet_concentration,inlet_flow
    REAL(dp) :: concentration,cumulative_mass

!!! Local variables
    INTEGER :: ne0,np1,np2
    REAL(dp) :: local_xi,mean_concentration,proportion_flow,&
         time_through_element,total_time
    LOGICAL :: CONTINUE

    cumulative_mass = 0.0_dp !used only if there is a 'unit' appended

    time_through_element = elem_field(ne_vol,ne)/(inlet_flow*elem_field(ne_Vdot,ne))
    total_time = time_through_element
    np1 = elem_nodes(1,ne)
    np2 = elem_nodes(2,ne)

    IF(total_time.GE.dt)THEN !location is within this element
       local_xi = 1.0_dp-dt/time_through_element
       concentration = (1.0_dp-local_xi)*node_field(nj_conc1,np1) + &
            local_xi*node_field(nj_conc1,np2)

       mean_concentration = 0.5_dp*(concentration+node_field(nj_conc1,np2))
       cumulative_mass = cumulative_mass + mean_concentration* &
            elem_field(ne_vol,ne)*(1.0_dp-local_xi)
    ELSE

       mean_concentration = 0.5_dp*(node_field(nj_conc1,np1)+node_field(nj_conc1,np2))
       cumulative_mass = cumulative_mass + mean_concentration*elem_field(ne_vol,ne)
       CONTINUE = .TRUE.
       ne0=ne

       DO WHILE (CONTINUE)

          IF(elem_cnct(-1,1,ne0).GT.0)THEN
             ne0 = elem_cnct(-1,1,ne0) !parent element
             np1 = elem_nodes(1,ne0)
             np2 = elem_nodes(2,ne0)

             time_through_element = elem_field(ne_vol,ne0)/&
                  (inlet_flow*elem_field(ne_Vdot,ne0))
             total_time = total_time + time_through_element

             proportion_flow = elem_field(ne_Vdot,ne)/elem_field(ne_Vdot,ne0)

             IF(total_time.GE.dt)THEN !location is within this element
                local_xi = (total_time-dt)/time_through_element !logic correct
                concentration = (1.0_dp-local_xi)*node_field(nj_conc1,np1) + &
                     local_xi*node_field(nj_conc1,np2)

                mean_concentration = 0.5_dp*(concentration+node_field(nj_conc1,np2))
                cumulative_mass = cumulative_mass + mean_concentration* &
                     elem_field(ne_vol,ne0)*(1.0_dp-local_xi)*proportion_flow
                CONTINUE=.FALSE.
             ELSE
                concentration = 0.5_dp*(node_field(nj_conc1,np1)+node_field(nj_conc1,np2))
                cumulative_mass = cumulative_mass + concentration* &
                     elem_field(ne_vol,ne0)*proportion_flow
             ENDIF
          ELSE
             CONTINUE=.FALSE.
             concentration = inlet_concentration

             proportion_flow = elem_field(ne_Vdot,ne)/elem_field(ne_Vdot,1)
             cumulative_mass = cumulative_mass + concentration* &
                  inlet_flow*(dt-total_time)*proportion_flow
          ENDIF
       ENDDO
    ENDIF

  END SUBROUTINE track_to_location

!!!##############################################################################

  SUBROUTINE normalised_slope(num_breaths,dt,t_fit_0,t_fit_1,&
       time_expiration,vol_expired,vol_frc,resultsfile)

    USE arrays,ONLY: dp

    IMPLICIT NONE

!!! Parameters
    INTEGER,INTENT(in) :: num_breaths
    REAL(dp),INTENT(in) :: dt,t_fit_0,t_fit_1,time_expiration,vol_expired,vol_frc
    CHARACTER(len=*),INTENT(in):: resultsfile

!!! Local variables
    INTEGER :: fh,i,ibeg,iend,ios=0, i_ss_end,nbreath,nlabel,num_data
    REAL(dp) :: exp_vol,exp_conc,exp_time,exp_vol_last,flow,&
         norm_slope,Sacin,Scond,slope,Sn0,Sn1,Sn2,sum_gas_expired,sum_n2_expired,&
         sum_x,sum_xsq,sum_xy,sum_y,turnovers,&
         turnover_0,turnover_1,turnover_2,vol_end,vol_start
    CHARACTER(len=100) :: buffer
    LOGICAL :: first,found_1,found_2

    vol_start = t_fit_0 * vol_expired ! the expired volume where fit starts from
    vol_end = t_fit_1 * vol_expired !the expired volume where fit goes to

    fh=20
    OPEN(fh, file=resultsfile, status='old') !open existing file

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.
    ios=0
    first = .TRUE.
    found_1 = .FALSE.
    found_2 = .FALSE.

    DO nbreath = 1,num_breaths
       exp_time = 0.0_dp
       exp_vol_last = 0.0_dp
       sum_gas_expired = 0.0_dp
       sum_x = 0.0_dp
       sum_y = 0.0_dp
       sum_xsq = 0.0_dp
       sum_xy = 0.0_dp
       i = 0
       DO WHILE (exp_time < time_expiration)
          READ(fh, '(A)', iostat=ios) buffer
          IF(ios .NE. 0)THEN
             WRITE(*,'('' Premature end of file, at breath'',I3,'' time'',F8.3)') nbreath,exp_time
             STOP
          ENDIF
          ! line contains: time, flow, DV, vol, 4 x %errors, expired conc
          i_ss_end = LEN(buffer)
          DO nlabel = 1,9
             ibeg = INDEX(buffer," ") + 1 !get location of first integer beyond ws in string
             buffer = ADJUSTL(buffer(ibeg:i_ss_end)) ! get info beyond ws, remove leading ws
             iend = INDEX(buffer," ") !get location of next ws in sub-string
             SELECT CASE (nlabel)
             CASE(2)
                READ (buffer(1:iend-1), *, iostat=ios) flow
             CASE(9)
                READ (buffer(1:iend-1), *, iostat=ios) exp_conc
             END SELECT
          ENDDO
          exp_vol = exp_vol_last + dabs(flow * dt)
          IF(exp_vol.GE.vol_start.AND.exp_vol.LE.vol_end)THEN
             i=i+1
             sum_x = sum_x + exp_vol
             sum_y = sum_y + exp_conc
             sum_xsq = sum_xsq + exp_vol**2
             sum_xy = sum_xy + exp_vol * exp_conc
          ENDIF
          sum_gas_expired = sum_gas_expired + dabs(flow) * dt * exp_conc
          exp_time = exp_time + dt
          exp_vol_last = exp_vol
       ENDDO ! while time
       num_data = i

       slope = -(num_data*sum_xy-sum_x*sum_y)/(num_data*sum_xsq-sum_x**2)*1.0e+6_dp !fractional conc/L
       sum_n2_expired = 1.0_dp - sum_gas_expired/vol_expired !fractional
       norm_slope=slope/sum_n2_expired

       turnovers = nbreath*vol_expired/vol_frc

       WRITE(*,'(I3,F7.2,F12.4,D13.4,F10.4)') nbreath,turnovers,norm_slope,slope,sum_n2_expired

       IF(first)THEN
          first = .FALSE.
          turnover_0 = turnovers
          Sn0 = norm_slope
       ENDIF
       IF(.NOT.found_1.AND.turnovers.GE.1.5_dp)THEN
          found_1 = .TRUE.
          turnover_1 = turnovers
          Sn1 = norm_slope
       ELSE IF(.NOT.found_2.AND.turnovers.GE.6.0_dp)THEN
          found_2 = .TRUE.
          turnover_2 = turnovers
          Sn2 = norm_slope
       ENDIF

    ENDDO ! num_breaths

    IF(found_1.AND.found_2)THEN
       Scond=(Sn2-Sn1)/(turnover_2-turnover_1)
       Sacin=Sn0-Scond*turnover_0
       WRITE(*,'('' Sacin = '',F10.4,'' L.-1, Scond = '',F10.4)') Sacin,Scond
       WRITE(*,'('' Calculated using TO = '',F10.4,'' and '',F10.4)') turnover_1,turnover_2
    ELSE
       WRITE(*,'('' Turnovers is'',F10.4,''; too small for calculating Scond and Sacin'')') turnovers
    ENDIF

    CLOSE(fh)

  END SUBROUTINE normalised_slope

END MODULE gasmix
