PROGRAM DARTRAY_TRUSTI
  use smooth_grid_routines
  use user_routines_trustI
  use rt_routines
  use io_routines
  use dartray_hub
  IMPLICIT NONE
  
  !! Initialize MPI 
  call initialize_mpi
  
  !! INPUT VARIABLE INITIALIZE  
  call input_initialize

!!! Read input file 
  call read_input_file

!!! Read wavelength grid 
  call read_lambda_list
  
!!! CHECK INPUT
  call check_input 
  
!!! GRID INITIALIZATION 
  call grid_initialize

!!! set TRUST Opacity and star luminosity
  call set_trustI
  
!!! write info file
  call write_file_info
  
!!! RT calculation
  call dartray_switch

!!! Terminate mpi 
  call terminate_mpi
  
CONTAINS
  
END PROGRAM DARTRAY_TRUSTI
