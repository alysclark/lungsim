MODULE arrays
  !*Brief Description:* This module defines arrays.
  !
  !*LICENSE:*
  !
  !
  !*Contributor(s):* Merryn Tawhai, Alys Clark
  !
  !*Full Description:*
  !
  !This module defines arrays

  IMPLICIT NONE

  INTEGER :: num_elems,num_elems_2d,num_nodes,num_data,num_nodes_2d,num_units,num_lines_2d,maxgen
  INTEGER, PARAMETER :: dp=KIND(0.d0) !  for double precision
  REAL(dp),PARAMETER :: zero_tol = 1.0e-12_dp
  REAL(dp),PARAMETER :: loose_tol = 1.0e-6_dp

  INTEGER,ALLOCATABLE :: nodes(:) !allocated in define_node_geometry
  INTEGER,ALLOCATABLE :: nodes_2d(:) !allocated in define_node_geometry_2d
  INTEGER,ALLOCATABLE :: node_versn_2d(:) !allocated in define_node_geometry_2d
  INTEGER,ALLOCATABLE :: elems(:) !allocated in define_1d_elements
  INTEGER,ALLOCATABLE :: lines_2d(:)
  INTEGER,ALLOCATABLE :: parentlist(:)
  INTEGER,ALLOCATABLE :: line_versn_2d(:,:,:)
  INTEGER,ALLOCATABLE :: lines_in_elem(:,:)
  INTEGER,ALLOCATABLE :: nodes_in_line(:,:,:)
  INTEGER,ALLOCATABLE :: elems_2d(:) !allocated in define_elem_geometry_2d
  INTEGER,ALLOCATABLE :: elem_cnct(:,:,:)  !NXI(-ni:ni,1,ne)
  INTEGER,ALLOCATABLE :: elem_cnct_2d(:,:,:)
  INTEGER,ALLOCATABLE :: elem_nodes(:,:)
  INTEGER,ALLOCATABLE :: elem_nodes_2d(:,:)
  INTEGER,ALLOCATABLE :: elem_versn_2d(:,:)
  INTEGER,ALLOCATABLE :: elem_lines_2d(:,:)
  INTEGER,ALLOCATABLE :: elem_ordrs(:,:)
  INTEGER,ALLOCATABLE :: elem_symmetry(:)
  INTEGER,ALLOCATABLE :: elem_units_below(:)
  INTEGER,ALLOCATABLE :: elems_at_node(:,:)
  INTEGER,ALLOCATABLE :: elems_at_node_2d(:,:)
  INTEGER,ALLOCATABLE :: units(:)

  REAL(dp),ALLOCATABLE :: arclength(:,:)
  REAL(dp),ALLOCATABLE :: elem_field(:,:) !properties of elements
  REAL(dp),ALLOCATABLE :: elem_direction(:,:)
  REAL(dp),ALLOCATABLE :: node_xyz(:,:)
  REAL(dp),ALLOCATABLE :: data_xyz(:,:)
  REAL(dp),ALLOCATABLE :: data_weight(:,:)
  REAL(dp),ALLOCATABLE :: node_xyz_2d(:,:,:,:)
  REAL(dp),ALLOCATABLE :: gasex_field(:,:) !gasexchange specific fields
  REAL(dp),ALLOCATABLE :: unit_field(:,:) !properties of elastic units
  REAL(dp),ALLOCATABLE :: node_field(:,:)
  REAL(dp),ALLOCATABLE :: scale_factors_2d(:,:)

  LOGICAL,ALLOCATABLE :: expansile(:)

  TYPE capillary_bf_parameters
     INTEGER :: num_symm_gen=9 !no units
     REAL(dp) :: total_cap_area=0.63000e02_dp !m
     REAL(dp) :: Palv=0.0_dp!Pa
     REAL(dp) :: H0=0.35000e-05_dp !m
     REAL(dp) :: K_cap=0.12000e02_dp
     REAL(dp) :: F_cap=0.18000e01_dp
     REAL(dp) :: F_sheet=0.10400e00_dp
     REAL(dp) :: sigma_cap=0.43637e03_dp !Pa
     REAL(dp) :: mu_c=0.19200e-02_dp !Pa.s
     REAL(dp) :: alpha_a=2.33e-08_dp !m/Pa
     REAL(dp) :: alpha_v=2.33e-08_dp !m/Pa
     REAL(dp) :: F_rec=0.64630e00_dp
     REAL(dp) :: sigma_rec=0.22300e04_dp
     REAL(dp) :: L_c=0.11880e-02_dp !m
     REAL(dp) :: Plb_c=0.0_dp !Pa
     REAL(dp) :: Pub_c=3138.24_dp !Pa
     REAL(dp) :: Pub_a_v=3138.24_dp !Pa
     REAL(dp) :: L_art_terminal=0.13000e-03_dp !m
     REAL(dp) :: L_vein_terminal=0.13000e-03_dp !m
     REAL(dp) :: R_art_terminal=0.10000e-04_dp !m
     REAL(dp) :: R_vein_terminal=0.90000e-05!m
  END TYPE capillary_bf_parameters

  TYPE admittance_param
     CHARACTER (len=20) :: admittance_type
     CHARACTER (len=20) :: bc_type
  END TYPE admittance_param
  TYPE, EXTENDS (admittance_param) :: two_parameter
  REAL(dp) :: admit_P1=1.0_dp
  REAL(dp) :: admit_P2=1.0_dp
END TYPE two_parameter
TYPE, EXTENDS (two_parameter) :: three_parameter
REAL(dp) :: admit_P3=1.0_dp
END TYPE three_parameter
TYPE, EXTENDS (three_parameter) :: four_parameter
REAL(dp) :: admit_P4=1.0_dp
END TYPE four_parameter
TYPE,EXTENDS (four_parameter) :: all_admit_param
END TYPE all_admit_param

TYPE elasticity_vessels
CHARACTER(len=20) ::vessel_type
END TYPE elasticity_vessels
TYPE, EXTENDS(elasticity_vessels) :: elasticity_param
REAL(dp) :: elasticity_parameters(3)=0.0_dp
END TYPE elasticity_param

TYPE fluid_properties
REAL(dp) :: blood_viscosity=0.33600e-02_dp !Pa.s
REAL(dp) :: blood_density=0.10500e-02_dp !kg/cm3
REAL(dp) :: air_viscosity
REAL(dp) :: air_density
END TYPE fluid_properties

! temporary, for debugging:
REAL(dp) :: unit_before

PRIVATE

PUBLIC set_node_field_value, elem_field, num_elems, num_elems_2d, elem_nodes, node_xyz, &
nodes,nodes_2d, elems, num_nodes, num_nodes_2d, num_data, data_xyz, data_weight, &
node_xyz_2d, node_versn_2d, units, num_units, unit_field, node_field, dp, &
elem_cnct, elem_ordrs, elem_direction, elems_at_node, elem_symmetry, expansile, &
elem_units_below, maxgen,capillary_bf_parameters, zero_tol,loose_tol,gasex_field, &
num_lines_2d, lines_2d, line_versn_2d, lines_in_elem, nodes_in_line, elems_2d, &
elem_cnct_2d, elem_nodes_2d, elem_versn_2d, elem_lines_2d, elems_at_node_2d, arclength, &
scale_factors_2d, parentlist, fluid_properties, elasticity_vessels, admittance_param, &
elasticity_param, all_admit_param

CONTAINS
SUBROUTINE set_node_field_value(row, col, VALUE)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SET_NODE_FIELD_VALUE" :: SET_NODE_FIELD_VALUE
IMPLICIT NONE

INTEGER, INTENT(in) :: row, col
REAL(dp), INTENT(in) :: VALUE

node_field(row, col) = VALUE

END SUBROUTINE set_node_field_value


END MODULE arrays
