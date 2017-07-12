!> This module contains subroutines that execute the RT algorithms. The rt_type() keywork allows all the subroutines in DART-Ray to distinguish between the different parts of the RT calculation. The rt_algorithm() keyword is used to distinguish between the different RT algorithms.  
MODULE dartray_hub
  use smooth_grid_routines
  use rt_routines
  use io_routines  
  use sed_routines
  use visual_routines
  IMPLICIT NONE 
CONTAINS 

   !> This is the procedure for the main RT algorithm. Used both for the stellar emission and dust emission RT calculations. In the case of the dust emission RT, it is called multiple times until convergence criteria is reached. 
 subroutine dartray_main
  
   rt_type = rtt_start

   do   !!! 

      select case (rt_type)
      case(rtt_start)
         call rt_prepare
      case(rtt_precalc_cell)
         call rt_loop
      case(rtt_precalc_src)        
         call rt_loop
         call rt_mpi
      case(rtt_output_part1)
         call make_output          
      case(rtt_prep_part2)
         call rt_prepare         
      case(rtt_dir_cell)
         call rt_loop
      case(rtt_dir_src)               
         call rt_loop
         call rt_mpi         
      case(rtt_i_obs_dir_cell)
         call rt_prepare
         call rt_loop_iobs
      case(rtt_i_obs_dir_src)
         call rt_loop_iobs
      case(rtt_output_part2)
         call rt_prepare
         call calc_sed
         call make_maps
         call make_output       
      case(rtt_scatt)
         do 
            call rt_prepare
            if (cnflag) exit ! exit if scattering iteration convergence reached 
            call rt_loop
            call rt_mpi            
         end do         
         call make_output
      case(rtt_i_obs)
         call rt_prepare
         call rt_loop_iobs
         call calc_sed
         call make_maps
         call make_output

      end select

      call select_rt_type
      if (done) exit
      
   end do

 end subroutine dartray_main

 !> This is the procedure for the 2D algorithm. The calculation of the i_obs() arrays is not implemented because the scaspe() arrays are symmetrized only for the directions of the HEALpix sphere, not the added directions for the observer lines-of-sight. 
subroutine dartray_main_2D
   
   rt_type = rtt_start

   do   !!! 

      select case (rt_type)
      case(rtt_start)
         call rt_prepare
      case(rtt_precalc_cell) 
         call rt_loop_2D
      case(rtt_precalc_src)        
         call rt_loop
         call rt_mpi      
      case (rtt_output_part1)
         call make_output
      case(rtt_prep_part2)
         call rt_prepare
      case(rtt_dir_cell)
         call rt_loop_2D
      case(rtt_dir_src)         
         call rt_loop
         call rt_mpi
      case(rtt_i_obs_dir_src)
         !call rt_loop_iobs             
      case(rtt_i_obs_dir_cell)
         call rt_prepare
         !call rt_loop_iobs
      case(rtt_output_part2)
         call rt_prepare
         call make_output
      case(rtt_scatt)
         do 
            call rt_prepare
            if (cnflag) exit ! exit if convergence reached 
            call rt_loop_2D
            call rt_mpi           
         end do         
         call make_output
      case(rtt_i_obs)
         call rt_prepare
         !call rt_loop_iobs   
         call make_output

      end select

      call select_rt_type
      if (done) exit
      
   end do

 end subroutine dartray_main_2D 

 !> This DART-Ray procedure calculates i_obs() arrays from the input arrays and pre-calculated scaspe_tot array. When calculating i_obs() with this procedure, the scattered luminosities considered are those of the HEALPix directions in the scaspe_tot() array, not the extra directions for the observer line-of-sight. This allows to calculate maps for observer directions and positions other than those used in the "main" RT calculation. To avoid overwriting of the i_obs() file name, the output file names are slightly changed compared to those set in set_filenames(). 
!> \todo for the i_obs_dust algorithm, you should also add the option of adding the radiation field due to dust emission when calculating the dust emission luminosity. At the moment only the radiation field due to stellar emission is considered (see grid_initialize_dust). 
 subroutine dartray_i_obs

   ! set filenames and units in case of dust algorithm
   if (rt_algorithm_ID == rta_i_obs_dust) then
      call grid_initialize_dust
      call set_dust_emission
   endif
   
   ! change filenames 
   call reset_filenames
   
   ! calculate i_obs part2 
   call prepare_p_src
   call create_scaspe
   call calc_total_luminosity
   call create_i_obs    ! escaping brigthness array
   rt_type =rtt_i_obs_dir_src
   call rt_loop_iobs
   rt_type =rtt_i_obs_dir_cell
   call rt_loop_iobs
   rt_type = rtt_output_part2
   if (print_output_part2) then  
      call calc_sed
      call make_maps
      call make_output
   endif

   if (only_direct_rt) return ! 

   ! read scaspe_tot array 
   call create_scaspe_tot
   rt_type=rtt_read_scaspe_tot
   call read_output 

   ! calculate i_obs
   rt_type=rtt_i_obs
   call rt_prepare
   call rt_loop_iobs
   call calc_sed
   call make_maps
   call make_output
   
 end subroutine dartray_i_obs


!> This dartray procedure calculates and prints the SEDs for the direct and total radiation luminosity. In order to work both i_obs_part2 and i_obs files have to be printed (it outputs an error if they don't exist). If the keyword print_maps is set, then it also calculates the surface brightness maps for the external observer. 
 subroutine dartray_sed

   ! set filenames and units in case of dust algorithm
   if (rt_algorithm_ID == rta_sed_dust) then
      rt_type = rtt_grid_init_dust
      lnum = lnum_dust
      call set_filenames
      call set_units
   endif

   call prepare_p_src
   if (rt_algorithm_ID == rta_sed_dust) tot_p_src = 0 ! no point source for dust RT
   call prepare_scaspe_splitting  ! this has to appear before create_i_obs (see lnum_node)
   call create_i_obs    ! escaping brigthness array  
   if (tot_ndir > 0) call calc_sed_arrays
   if (main_prc) then 
      if (print_output_part2) call print_sed_arr_dir
      call print_sed_arr
   endif
  
   if (print_maps .or. print_maps_in) then
      if (print_output_part2) then 
         rt_type=rtt_read_i_obs_part2
         call read_output
         call make_maps
         if (main_prc .and. tot_ndir > 0) call print_map_arr_out(file_maps_part2)
         if (main_prc .and. tot_ndir_in > 0) call print_map_in_arr_out(file_maps_in_part2)
      endif

      rt_type=rtt_read_i_obs
      call read_output
      call make_maps
      if (main_prc .and. tot_ndir > 0)  call print_map_arr_out(file_maps)
      if (main_prc .and. tot_ndir_in > 0) call print_map_in_arr_out(file_maps_in)
   endif
  
 end subroutine dartray_sed

 !> Dartray procedure for the dust emission RT calculation.
 !> \todo Add isotropic scattering option and change scaspe arrays accordingly. You might consider making a data structure with two types of scaspe arrays, one depending on direction and one not. THIS IS FOR NEXT VERSION
 subroutine dartray_dust

   ! return if no_dust_rt is set
   if (no_dust_rt) then
      if (main_prc) print *, 'STOP: no_dust_rt = TRUE. No dust emission RT will be performed!'
      return
   endif
      
   ! initialize dust RT grids 
   call grid_initialize_dust

   ! dust heating iterations 
   do 

      ! set dens_stars for dust emission RT and check convergence criteria 
      call set_dust_emission
      
      if (cnflag_dust) exit
      
      if (rt_algorithm_ID == rta_dust) then 
      ! run RT main algorithm
         call dartray_main
      elseif (rt_algorithm_ID == rta_dust2d) then 
         call dartray_main_2D
      else
         if (main_prc) print *, 'STOP(dartray_dust): rt_algorithm_ID not recognized!'
         call stop_prc
      endif

      ! update i_obs_arr_tot
      call update_i_obs_arr
      
   end do

   ! print message
   if (main_prc) print *, 'Dust RT iterations completed' 

   ! make output
   call make_output_dust
      
 end subroutine dartray_dust 

!> Algorithm to project physical quantities on the observer maps, such as optical depth or intrinsic stellar luminosity (without dust extinction).  
subroutine dartray_projection 

! set rt_type for projection initialization 
rt_type = rtt_grid_init_projection

! set filenames 
call set_filenames

! create i_obs arrays 
call prepare_scaspe_splitting ! always before create_i_obs
call create_i_obs

! assign values to i_obs arrays 
call assign_i_obs_to_project

! set number of point sources = 0 if projection of optical depth
if (param_to_project == 'optical_depth') tot_p_src = 0

! calculate maps   
call make_maps

! print maps 
if (main_prc .and. tot_ndir > 0 .and. print_maps)  call print_map_arr_out(file_maps)
if (main_prc .and. tot_ndir_in > 0 .and. print_maps_in) call print_map_in_arr_out(file_maps_in)

end subroutine dartray_projection



 !> Selects Dartray algorithm depending on input rt_algorithm().
 subroutine dartray_switch

   select case(rt_algorithm)
   case('main')
      call dartray_main
      call dartray_dust
   case('2D')
      call dartray_main_2D
      call dartray_dust
   case('i_obs', 'i_obs_dust')
      call dartray_i_obs
   case('sed', 'sed_dust')
      call dartray_sed
   case('dust', 'dust2D')
      call dartray_dust
   case('projection')
      call dartray_projection
   case default
      print *, 'STOP: something wrong here: dartray_switch!'
      STOP
   end select

 end subroutine dartray_switch


 
 !> This subroutine is used to select the right rt_type(). It simply set rt_type to the next value in the list.
 !> /todo It would be good to print an appropriate message when a new part of the calculation is started. For now it prints the ID numbers for the rt calculation step. 
  subroutine select_rt_type
   integer, parameter  :: ntype=12
   integer :: i 
   integer :: rt_type_arr(ntype)
   integer :: ierr 
   ! this select the type of rt calculation to be performed
   ! 1) precalc_src   ( precalculation for stellar point sources )
   ! 2) precalc_cell  ( )
   ! 3) dir_src       ( direct light calculation for stellar point sources )
   ! 4) dir_cell      ( direct light calculation for emitting cells )
   ! 5) scatt         ( scattering iterations )

   done = .FALSE.
   ! Start program 
   rt_type_arr(1) = rtt_start
   ! Precalculation for stellar point sources. Note: the order cell - > src is not important for the 3D mode but it is for the 2D mode because of fix symmetry routines. 
   rt_type_arr(2) = rtt_precalc_cell
   ! Precalculation for emitting cells. 
   rt_type_arr(3) = rtt_precalc_src
   ! Write output part 1
   rt_type_arr(4) = rtt_output_part1
   ! Prepare for part 2.
   rt_type_arr(5) = rtt_prep_part2
   ! Direct light calculation for emitting cells.
   rt_type_arr(6) = rtt_dir_cell
   ! Direct light calculation for point sources.
   rt_type_arr(7) = rtt_dir_src
   ! Outgoing radiation specific intensity calculation for emitting cells. 
   rt_type_arr(8) = rtt_i_obs_dir_cell
   ! Outgoing radiation specific intensity calculation for point sources. 
   rt_type_arr(9) = rtt_i_obs_dir_src
   ! Write output part 2. 
   rt_type_arr(10) = rtt_output_part2
   ! Scattered radiation calculation. 
   rt_type_arr(11) = rtt_scatt
   ! Outgoing radiation specific intensity calculation for scattered radiation sources. 
   rt_type_arr(12) = rtt_i_obs

   ! in case of dust RT algorithm, skip first part if iterations_dustem > 1
   if (rt_algorithm_ID == rta_dust .or. rt_algorithm_ID == rta_dust2d) then
      if (rt_type == rtt_start .and. iterations_dustem > 1) then
         rt_type = rtt_output_part1
      endif
   endif

   ! only direct light calculation if only_direct_rt = .TRUE.
   if (only_direct_rt .and. rt_type == rtt_output_part2) then 
      rt_type = rtt_i_obs
   endif

   do i =1, ntype

      if (rt_type == rt_type_arr(i)) exit 

   enddo

   if (i < ntype) then
      rt_type = rt_type_arr(i+1)
      if (main_prc) print *, 'Start step', rt_type
      call mpi_barrier(MPI_COMM_WORLD, ierr)
   else
      done = .TRUE.
      if (main_prc) print *, 'RT calculation completed' 

   endif

 end subroutine select_rt_type


 !> This routine performs several checks before the beginning of some steps in the rt calculation (look at each case for more info). For example, it checks whether another scattering iteration is necessary or whether the intermediate output files have to be restored and the direct light calculation to be skipped.
 subroutine rt_prepare
 
  select case (rt_type)

  case(rtt_start)
     !call reduce_grid_res
     call set_en_lim
     call prepare_p_src 
     call calc_total_luminosity
     call check_files
     call check_memory
          
     if (rt_algorithm_ID == rta_2D .or. rt_algorithm_ID == rta_dust2D) then 
        call check_grid_symmetry
        call check_2d_src
     endif
      
     if (file_restore_part1) then
        call read_output
        rt_type = rtt_output_part1
     endif

     if (file_restore) then
        !call restore_grid_original_res
        call create_scaspe   ! scattered energy array
        call create_i_obs    ! escaping brigthness array
        call pre_calc_ads_arr
        call create_psel_av_arr
        call read_output
        call create_mpi_type 
        rt_type=rtt_output_part2  !! in this way you go directly to step 'scatt'
     endif

  case(rtt_prep_part2)     
     call create_scaspe   ! scattered energy array   
     call create_i_obs    ! escaping brigthness array
     call pre_calc_ads_arr
     call create_psel_av_arr
     call create_mpi_type 
     
  case(rtt_i_obs_dir_cell)
     !call restore_grid_original_res
     
  case(rtt_output_part2)
     call update_i_obs_arr
     
  case(rtt_scatt) 
   !  npix_hp=12*nside_sca*nside_sca
     u_fest_arr=u_final_arr
     lum_lost_prev=lum_lost
     iterations= iterations +1

     if (iterations == 1) then
        !call reduce_grid_res  !!! this is here because grid resolution has been restored to create output part2 
        call create_scaspe_tot
        if (rt_algorithm_ID == rta_2D .or. rt_algorithm_ID == rta_dust2d .or. sequential_scattering) call create_scaspe_prev                   
        bm_par=bm_par_sca  !!! different bm_par value during scattering iterations       
     endif

     if (rt_algorithm_ID == rta_2D  .or. rt_algorithm_ID == rta_dust2d .or. sequential_scattering) then 
        scaspe_prev_arr=scaspe_arr       
     endif

     call calc_total_luminosity_sca

     if (.not.cnflag) then
        if (main_prc) print *, 'start scattering iteration number', iterations 
        call expand_psel_av_arr
        !read(*,*)
     else
        !call restore_grid_original_res  !!! restore grid res before going to next step 
     endif

  case(rtt_i_obs)
     call calc_total_luminosity_sca
     
  case default  

     print *, 'you should not get here: rt_prepare'
     call stop_prc

  end select

end subroutine rt_prepare

!> This routine calls all the routines to initialize the grid and sources/observer parameters.
subroutine grid_initialize

   rt_type = rtt_grid_init_stars

   !iterations = 0  ! these are initialized  in read_main_grid
   !iterations_dustem = 0
   !nside_sca=2**kp_sca
   
   call read_main_grid
   if (use_lambda_grid) call read_lambda_grid
   call set_filenames   
   call set_base
   call make_csize_arr
   call read_dir_out
   call read_pos_obs
   call read_p_src
   call set_chunk_size
   call set_lambda_arr_si
   call set_units   
   if (print_sed) then
      call create_sed_arr
   endif
   call prepare_dust_model
   call scale_dens_arr
   call set_npix_arr  
   call set_walls
   
end subroutine grid_initialize

!> Initializes arrays for dust RT algorithms. Note that lnum() is still equal to lnum_stars() in the first called subroutines.
subroutine grid_initialize_dust

  rt_type = rtt_grid_init_dust
  
  ! read u_final_arr for the stellar emission if necessary. 
  if (rt_algorithm_ID == rta_dust .or. rt_algorithm_ID == rta_dust2D .or. rt_algorithm_ID == rta_i_obs_dust) then
     call read_ufield_arr
  endif

  ! reshape arrays 
  call store_reshape_arrays

  ! start of the dust RT algorithm
  
!!!--------------------!!!!
  ! Firstly, some needed initializations
  lnum = lnum_dust   ! set lnum to dust emission value
  if (rt_algorithm_ID == rta_main) rt_algorithm_ID = rta_dust
  if (rt_algorithm_ID == rta_2D) rt_algorithm_ID = rta_dust2d
  iterations = 0  !(see comment on this in read_main_grid)
  iterations_dustem = 0 !(see set_dust_emission)
  !kp_sca_max = 0  ! low resolution scattering for dust emission
  !nside_sca=2**kp_sca
  tot_p_src = 0  ! no point sources in dust emission RT calculation.
  src_cell = -1 ! no cells hosting point sources! 
  bm_par = bm_par_sca ! less ray angular density for dust RT
  cnflag_dust = .FALSE. 
  tau_cell_max = 0  ! No grid reduction for the moment 
!!! -------------------!!!!

  ! set units
  call set_units
    
  ! scale dens_arr
  call scale_dens_arr

  ! transform stellar emission U field into I field
  call convert_ufield_ifield

  ! set new output filenames and check existence
  call set_filenames

  ! initialise sed output arrays 
  if (print_sed) then
      call create_sed_arr
   endif
  
 end subroutine grid_initialize_dust

!> Calls the subroutines to print the output files for the dust RT algorithms.
 !> \todo DONE set right units for u_final_arr output files (always J/m/pc^3).
 !> \todo add SED for internal observer to calc_sed
subroutine make_output_dust

  ! part 2 (direct light)
  rt_type = rtt_output_part2
  if (tot_ndir > 0) then
     i_obs_arr = i_obs_arr_dir
  endif
  if (tot_ndir_in > 0) then
     i_obs_in_arr = i_obs_in_arr_dir
  endif
  call calc_sed
  call make_maps
  call make_output

  if (only_direct_rt) return

  rt_type = rtt_scatt
  call make_output

  rt_type = rtt_i_obs
   if (tot_ndir > 0) then
     i_obs_arr = i_obs_arr_tot
  endif
  if (tot_ndir_in > 0) then
     i_obs_in_arr = i_obs_in_arr_tot
  endif
  call calc_sed
  call make_maps
  call make_output
  
end subroutine make_output_dust
 


END MODULE dartray_hub
