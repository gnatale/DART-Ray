!> Contains the all the routines for I/O.  
MODULE io_routines
  use HDF5
  use smooth_grid_routines
  use rt_routines
  IMPLICIT NONE
  !> @param dir_grid Directory where the grid files are located 
  character(LEN=lcar) :: dir_grid
  !> @param dir_runs Directory for the output files. If it does not exist, it is created automatically in check_input().
  character(LEN=lcar) :: dir_runs
  !> @param grid_file Main grid file name
  character(LEN=lcar) :: grid_file
  !> @param grid_info_file Grid info file name 
  character(LEN=lcar) :: grid_info_file
  !> @param label_model Label model within grid_file()
  character(LEN=lcar) :: label_model
  !> @param label_model_out Label contained in all the output file names. 
  character(LEN=lcar) :: label_model_out
  !> @param label_model_out_i_obs Label contained in all the output file names for the i_obs and i_obs_dust RT algorithm. If not provided, it is set automatically to label_model_out().
  character(LEN=lcar) :: label_model_out_i_obs
  !> @param label_model_lambda_grid Label model within grid_file_lambda() 
  character(LEN=lcar) :: label_model_lambda_grid
  !> @param label_wave Wavelength label for a single wavelength.
  character(LEN=lcar) :: label_wave
  !> @param label_wave_arr Wavelength labels for all wavelengths.
  character(LEN=lcar), allocatable :: label_wave_arr(:)
  !> @param grid_file_lambda Lambda grid file name
  character(LEN=lcar) :: grid_file_lambda
  !> @param grid_file_lambda_arr Lambda grid file name for all wavelengths.
  character(LEN=lcar), allocatable :: grid_file_lambda_arr(:)
  !> @param file_dir_out File containing the angular directions \f$ \theta \f$ and \f$ \phi \f$ specifying the lines-of-sight of the far-away observers. The first line has to be a comment. The following lines have to be the list of \f$ \theta \f$ (first column) and \f$ \phi \f$ (second column) in radiants. Note that \f$ \theta \f$ is defined as the angle between the line-of-sight and the Z-axis. \f$ \phi \f$ is the line-of-sight position angle on the XY plane measured from the X-axis.    
  character(LEN=lcar) :: file_dir_out
  !> @param file_pos_obs File containing the 3D coordinates of the internal observers within the RT model. The first line has to be a comment. The following lines contain the X, Y,Z coordinates for each internal observer.
  character(LEN=lcar) :: file_pos_obs
  !> @param file_p_src File containing the 3D coordinates for the point sources within the RT model The first line has to be a comment. The following lines contain the X, Y,Z coordinates for each point source. 
  character(LEN=lcar) :: file_p_src
  !> @param file_u_fest_part1 File name of the output for part1. 
  character(LEN=lcar) :: file_u_fest_part1
  !> @param file_u_fest_part1_arr File names of the output for part1 for all wavelengths.
  character(LEN=lcar), allocatable :: file_u_fest_part1_arr(:)
  !> @param file_ufield_part2 File containing u_final() as calculated at the end of the direct light processing. 
  character(LEN=lcar) :: file_ufield_part2
  !> @param file_ufield_part2_arr File names of the output containing u_final() as calculated at the end of the direct light processing for all wavelengths. 
  character(LEN=lcar), allocatable :: file_ufield_part2_arr(:)
  !> @param file_ufield File containing u_final() as calculated at the end of the scattered light processing. 
  character(LEN=lcar) :: file_ufield
  !> @param file_ufield_arr File names of the output containing u_final() as calculated at the end of the scattered light processing for all wavelengths. 
  character(LEN=lcar), allocatable :: file_ufield_arr(:)
  !> @param file_info File where the input parameters of the RT calculation are stored. 
  character(LEN=lcar) :: file_info
  !> @param file_scaspe_part2 File containing scaspe() as calculated at the end of the direct light processing. It can be restored but, usually, it should not be printed because it can be very large. Set print_scaspe_part2 = .TRUE. if you want this file to be printed. 
  character(LEN=lcar) :: file_scaspe_part2
  !> @param file_scaspe_part2_arr File names of the output containing scaspe() as calculated at the end of the direct light processing for all wavelengths. It can be restored but, usually, it should not be printed because it can be very large. Set print_scaspe_part2 = .TRUE. if you want these files to be printed. 
  character(LEN=lcar), allocatable :: file_scaspe_part2_arr(:)
  !> @param file_i_obs File containing i_obs() as calculated at the end of the scattered light iterations.
  character(LEN=lcar) :: file_i_obs
  !> @param file_i_obs_arr File names of the output containing i_obs() as calculated at the end of the scattered light iterations for all wavelengths.
  character(LEN=lcar), allocatable :: file_i_obs_arr(:)
  !> @param file_i_obs_part2 File containing i_obs() as calculated at the end of the direct light processing.
  character(LEN=lcar) :: file_i_obs_part2
  !> @param file_i_obs_part2_arr File names of the output containing i_obs() as calculated at the end of the direct light processing for all wavelengths.
  character(LEN=lcar), allocatable :: file_i_obs_part2_arr(:)
  !> @param file_lum_lost File containing the percentage of luminosity lost during the entire RT calculation at each wavelength (the sum of the luminosity of the rays blocked because of the energy condition \f$ \delta U < f_U*U_{LL} \f$).
  character(LEN=lcar) :: file_lum_lost
  !> @param file_lum_lost_part2 File containing the percentage of luminosity lost during the direct light processing (the sum of the luminosity of the rays blocked because of the energy condition \f$ \delta U < f_U*U_{LL} \f$). 
  character(LEN=lcar) :: file_lum_lost_part2
  !> @param file_i_obs_in File containing i_obs_in() as calculated at the end of the scattering iterations. 
  character(LEN=lcar) :: file_i_obs_in
  !> @param file_i_obs_in_arr File names of the output containing i_obs_in() as calculated at the end of the scattering iterations for all wavelengths. 
  character(LEN=lcar), allocatable :: file_i_obs_in_arr(:)  
  !> @param file_i_obs_in_part2 File containing i_obs_in() as calculated at the end of the direct light processing. 
  character(LEN=lcar) :: file_i_obs_in_part2
  !> @param file_i_obs_in_part2_arr File names of the output containing i_obs_in() as calculated at the end of the direct light processing for all wavelengths. 
  character(LEN=lcar), allocatable :: file_i_obs_in_part2_arr(:)
  !> @param file_scaspe_tot File containing scaspe_tot(). 
  character(LEN=lcar) :: file_scaspe_tot
  !> @param file_scaspe_tot_arr File names of the output containing scaspe_tot() for all wavelengths.
  character(LEN=lcar), allocatable :: file_scaspe_tot_arr(:)
  !> @param file_psel_av_part2 File containing psel_av_arr() as derived at the end of the direct light processing. 
  character(LEN=lcar) :: file_psel_av_part2
  !> @param file_psel_av File containing psel_av_arr() as derived at the end of the scattered light iterations. 
  character(LEN=lcar) :: file_psel_av
  !> @param file_sed_arr_dir File containing emission SEDs, for each observer line of sight, obtained at the end of direct light processing. 
  character(LEN=lcar) :: file_sed_arr_dir
  !> @param file_sed_arr File containing emission SEDs, for each observer line of sight, obtained at the end of scattered light processing.
  character(LEN=lcar) :: file_sed_arr
  !> @param file_sed_dust_part2 File containing dust emission SEDs, for each observer line of sight, obtained at the end of direct light processing. Note that for the dust emission an additional iteration procedure is necessary to find temperature convergence.  
  character(LEN=lcar) :: file_sed_dust_part2
  !> @param file_sed_dust File containing dust emission SEDs, for each observer line of sight, obtained at the end of scattered light processing. Note that for the dust emission an additional iteration procedure is necessary to find temperature convergence.  
  character(LEN=lcar) :: file_sed_dust
  !> @param file_maps_part2 File containing the surface brightness maps for the direct light obtained by an observer located far away from the model. Units in MJy/sr. 
  character(LEN=lcar) :: file_maps_part2
  !> @param file_maps File containing the surface brightness maps for the direct+scattered light obtained by an observer located far away from the model. Units in MJy/sr. 
  character(LEN=lcar) :: file_maps
  !> @param file_maps_in_part2 File containing the surface brightness maps for the direct light obtained by an observer located inside the model. Units in MJy/sr. 
  character(LEN=lcar) :: file_maps_in_part2
  !> @param file_maps_in File containing the surface brightness maps for the direct+scattered light obtained by an observer located inside the model. Units in MJy/sr. 
  character(LEN=lcar) :: file_maps_in
  !> @param file_gra_fa File containing input grain size distribution for Graphite grains. To use input grain size distributions, set dust_model = 'user'.
  character(LEN=lcar) :: file_gra_fa
  !> @param file_sil_fa File containing input grain size distribution for Silicate grains. To use input grain size distributions, set dust_model = 'user'.
  character(LEN=lcar) :: file_sil_fa
  !> @param file_pah_neu_fa File containing input grain size distribution for neutral PAH molecules. To use input grain size distributions, set dust_model = 'user'.
  character(LEN=lcar) :: file_pah_neu_fa  
  !> @param file_pah_ion_fa File containing input grain size distribution for ionized PAH molecules. To use input grain size distributions, set dust_model = 'user'.
  character(LEN=lcar) :: file_pah_ion_fa
  !> @param file_q_gra File containing the \f$Q_\lambda\f$ opacity and \f$g_\lambda\f$ anisotropy coefficients for Graphite grains.
  character(LEN=lcar) :: file_q_gra
  !> @param file_q_sil File containing the \f$Q_\lambda\f$ opacity and \f$g_\lambda\f$ anisotropy coefficients for Silicates grains.
  character(LEN=lcar) :: file_q_sil
  !> @param file_q_pah_neu File containing the \f$Q_\lambda\f$ opacity and \f$g_\lambda\f$ anisotropy coefficients for neutral PAH molecules.
  character(LEN=lcar) :: file_q_pah_neu
  !> @param file_q_pah_ion File containing the \f$Q_\lambda\f$ opacity and \f$g_\lambda\f$ anisotropy coefficients for ionized PAH molecules.
  character(LEN=lcar) :: file_q_pah_ion 
  !> @param file_calorimetry_Gra File containing the grain density, specific enthalpy and specific heat capacity for Graphite and PAH grains 
  character(LEN=lcar) :: file_calorimetry_Gra
  !> @param file_calorimetry_Sil File containing the grain density, specific enthalpy and specific heat capacity for Silicates grains 
  character(LEN=lcar) :: file_calorimetry_Sil
  !> @param file_param_src File containing the properties used to define the spectra of the point sources. It can be arbitrarily structured since it is read in the "user_routines" modules. 
  character(LEN=lcar) :: file_param_src

  !> @param narr_grid Number of arrays in the main grid file grid_file(). 
  integer, parameter :: narr_grid =9
  !> @param dims Dimensions of the arrays to be read or printed.
  INTEGER(HSIZE_T), allocatable  :: dims(:)  ! Dataset dimensions
  !> @param print_output_part1 Equal to TRUE if array u_fest() has to be printed at the end of precalculation phase (output file file_u_fest_part1()).  
  logical :: print_output_part1
  !> @param print_output_part2 Equal to TRUE if the output of part2 has to be written on disk (output files: file_ufield_part2(), file_i_obs_part2(), file_i_obs_in_part2()). This should be done if one needs the u_final(), i_obs() and i_obs_in() arrays due only to the direct light. Note that the scaspe() array is not printed by default. To print that as well at the end of the direct light processing, you need to set print_scaspe_part2() = .TRUE. If you print all the part 2 output, including the scaspe() array, you will be able later to restart the calculation from the scattering iterations.      
  logical :: print_output_part2
  !> @param print_scaspe_part2 Equal to TRUE if scaspe() array has to be written on disk at the end of the direct light processing. This is not included by default in print_output_part2() because the scaspe() arrays can take a LOT of disk space. If you want to be able to restart the calculation from the scattering iterations, the scaspe() 
  logical :: print_scaspe_part2
  
  !> @param print_sed TRUE if file_sed() has to be printed.
  logical :: print_sed
  !> @param restore_file_mpi TRUE if previously written files have to be restored during runs with more than one MPI process. In case only one MPI process is present, choice is requested in the terminal.
  logical restore_file_mpi
  !> @param file_restore Logical variable equal to TRUE when all part 2 output files exists (that is the files printed at the end of the direct light processing). In this case the code can read them and start directly from the scattering iterations. 
  logical :: file_restore
  !> @param file_restore_part1 TRUE when part1 output file already exists and file_restore = .FALSE. In this case the code can read this file and start from the direct light processing.   
  logical :: file_restore_part1 
  !> @param input_av_opacities TRUE if integrated/averaged opacity coefficients kext(), kabs(), kext(), gsca() have to be read from input file file_av_opacities
  logical :: input_av_opacities
  !> @param file_av_opacities File containing the integrated/averaged opacity coefficients kext(), kabs(), kext(), gsca(). It is read only if keyword input_av_opacities is set TRUE in the input file. 
  character (LEN=lcar) :: file_av_opacities
  !> @param no_dust_rt TRUE if no dust emission RT has to be performed.
  logical :: no_dust_rt
  
  
  !> @param file_nbody_sph File name of the imported Nbody-SPH simulation. This file contains the following arrays:
!> 'agestar': the age of the stellar particles in Gyr;
!> 'fehgas': [Fe/H] ratio of the gas particles;
!> 'fehstar': [Fe/H] ratio of the stellar particles;
!> 'gascoord': 3D coordinates of the gas particles;
!> 'gastemp': Gas particle temperature [K];
!> 'mgas': Gas particle mass in M(sun);
!> 'mstar': Star particle mass in M(sun);
!> 'ofegas': [O/Fe] ratio for the gas particles;
!> 'starcoord': 3D coordinates of the stellar particles. 
!>This file can be created from a tipsy file using the python script tipsy2dartray.py. Note that the gas metallicity is inferred from the O abundance, the stellar metallicity from the Fe abundance. 
  character (LEN=lcar) :: file_nbody_sph

  !> @param stellar_library Name of the stellar library that has to be used to convert mass into luminosity. Choices are: 'starburst99' (standard starburst99 Single age stellar population SED, Kroupa IMF), 'maraston05_kr_rhb' (Maraston et al.2005 SSP SEDs, Kroupa IMF, red horizontal branch), 'user' (user provided table, see file_stellar_library().
  character (LEN=lcar) :: stellar_library

  !> @param file_stellar_library Name of the HDF5 file containing the stellar library specified in stellar_library(). This is an input parameter only if stellar_library() = 'user'. This file contains three arrays: lambda_lib_arr (wavelength list in [um]),  met_arr (metallicity), lum_to_mass_arr (stellar luminosity-to-mass ratio in erg/s/Hz/M(sun)). See e.g. STELLAR_LIBRARIES/starburst99/table_lum_mass_vs_age_met_starburst99.h5. 
  character (LEN=lcar) :: file_stellar_library 

  !> @param use_stellar_library TRUE if stellar library has to be loaded. 
  logical :: use_stellar_library

  !> @param file_pcell File containing the pcell_star() and pcell_gas() arrays. 
  character (LEN=lcar) :: file_pcell
  
  
  ! INPUT NAMELIST ! 
  namelist /dartray_input_strings/ label_model_lambda_grid,label_model_out,label_model_out_i_obs, file_dir_out, file_pos_obs, file_p_src, file_lambda_list, dir_runs, dir_grid,   grid_file, rt_algorithm, units_luminosity, units_csize, units_lambda,dust_model,dust_opacity_tables, file_gra_fa, file_sil_fa, file_pah_neu_fa, file_pah_ion_fa, file_av_opacities, dust_heating_type, file_q_gra, file_q_sil, file_q_pah_neu, file_q_pah_ion, file_calorimetry_Gra, file_calorimetry_Sil,file_nbody_sph,file_stellar_library,stellar_library, param_to_project, file_param_src
  namelist /dartray_input_var/ kp_sca_max, rad_lim, accuracy, conv_en_lim, bm_par,   bm_par_sca, bm_par_max, lambda_ref,max_lambda_stars, min_lambda_dust, dist_obs,ind_i_obs, ind_out_maps, n_dust_size_qabs, n_dust_wave_qabs, tau_cell_max, n_dust_temp_cal, npixel_maps, map_size_factor, kp_maps, x_wall_coord, y_wall_coord, z_wall_coord, z_sun, max_sca_iterations, n_int_rf_bins
  namelist /dartray_input_logical/ print_scaspe_tot, print_output_part1, print_output_part2, print_scaspe_part2, restore_file_mpi, use_lambda_grid, use_dir_out,  use_pos_obs, use_p_src,print_psel_av, sequential_scattering, print_sed,input_av_opacities, no_communications, no_dust_rt,only_direct_rt, test_run, print_maps, print_maps_in, x_wall_on, y_wall_on, z_wall_on, use_stellar_library, limit_scattering_iterations 


CONTAINS


 !> This subroutine prepares the dsetname and rank for all the elements in grid_file().    
  subroutine make_dsetname_main_grid(dsetname,rank)
    INTEGER     ::   rank(0:)
    character(LEN=lcar) :: dsetname(0:)

    dsetname(0)='ncell'     ; rank(0)=1 
    dsetname(1)='ccoord'    ; rank(1)=2 
    dsetname(2)='cchild'    ; rank(2)=1 
    dsetname(3)='cindex'    ; rank(3)=1 
    dsetname(4)='lvl'       ; rank(4)=1 
    dsetname(5)='csize'     ; rank(5)=1 
    dsetname(6)='dens'      ; rank(6)=1 
    dsetname(7)='dens_stars'  ; rank(7)=1
    dsetname(8)='base'      ; rank(8)=1

  end subroutine make_dsetname_main_grid
  
  !> This subroutine prints the main grid in grid_file().
  subroutine print_3d_grid_file

    USE ISO_C_BINDING
    
    INTEGER(HID_T) :: file_id       
     INTEGER(HID_T) :: dset_id       ! Dataset identifier
     INTEGER(HID_T) :: dspace_id     ! Dataspace identifier     
     INTEGER     ::   rank(0:narr_grid-1)                         ! Dataset rank
     INTEGER     ::   error  ! Error flag
     character(LEN=lcar) :: dsetname(0:narr_grid-1)
     integer :: i
     TYPE(C_PTR) :: f_ptr
     
     call make_dsetname_main_grid(dsetname, rank)    

!    Initialize FORTRAN interface.
     CALL h5open_f (error)

 ! Create a new file using default properties.
     ! 
     CALL h5fcreate_f(trim(adjustl(dir_grid))//trim(adjustl(grid_file)), H5F_ACC_TRUNC_F, file_id, error)    


     do i=0,narr_grid-1
        
      CALL  make_dims(dsetname(i),rank(i))
     
     ! Create the dataspace.     
     CALL h5screate_simple_f(rank(i), dims, dspace_id, error)
     
! Create the dataset with default properties.
     select case(dsetname(i))
     case('ncell', 'cchild', 'lvl')
        CALL h5dcreate_f(file_id, dsetname(i), H5T_NATIVE_INTEGER , dspace_id, &
             dset_id, error)
     case('cindex')
        CALL h5dcreate_f(file_id, dsetname(i), h5kind_to_type(int64,H5_INTEGER_KIND), dspace_id, &
             dset_id, error)
     case default   
        CALL h5dcreate_f(file_id, dsetname(i),H5T_NATIVE_DOUBLE , dspace_id, &
          dset_id, error)
     end select
 
! Open an existing dataset.
     !
     CALL h5dopen_f(file_id, dsetname(i), dset_id, error)  
   ! note the input arrays might actually be bigger than dims when created with the standard create_adap_program
  if (dsetname(i) == 'ncell' ) then 
     CALL h5dwrite_f(dset_id,H5T_NATIVE_INTEGER , ncell, dims, error)
  else if (dsetname(i) == 'ccoord' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , ccoord, dims, error)
  else if (dsetname(i) == 'cchild' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_INTEGER , cchild, dims, error)   
  else if (dsetname(i) == 'cindex' ) then
    ! CALL h5dwrite_f(dset_id,H5T_NATIVE_INTEGER , cindex, dims, error)    
     f_ptr=c_loc(cindex(0))
     CALL h5dwrite_f(dset_id,h5kind_to_type(int64, H5_INTEGER_KIND), f_ptr,  error)
  else if (dsetname(i) == 'lvl' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_INTEGER , lvl, dims, error)
  else if (dsetname(i) == 'csize' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , csize, dims, error)
  else if (dsetname(i) == 'dens' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , dens, dims, error)
  else if (dsetname(i) == 'dens_stars' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , dens_stars, dims, error)
  else if (dsetname(i) == 'base' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_INTEGER , base, dims, error)   
  endif

! End access to the dataset and release resources used by it.
     !
     CALL h5dclose_f(dset_id, error)

     !
     ! Terminate access to the data space.
     !
     CALL h5sclose_f(dspace_id, error)

  
       deallocate(dims) 
    end do
        
  ! Terminate access to the file.
     !
     CALL h5fclose_f(file_id, error)     

!    Close FORTRAN interface.
!
     CALL h5close_f(error)

       print *, '3D grid printed'
 
  end subroutine print_3d_grid_file

  !> This subroutine reads the main grid from grid_file().
  subroutine read_main_grid
    USE ISO_C_BINDING
    INTEGER(HID_T) :: file_id                            ! File identifier
    INTEGER(HID_T) :: dset_id       ! Dataset identifier
    INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
    CHARACTER (LEN=lcar) :: dsetname(0:narr_grid-1) ! Dataset name
    INTEGER     ::   rank(0:narr_grid-1)   
    INTEGER(HSIZE_T) :: num
    INTEGER     ::   error  ! Error flag
    integer :: i
    TYPE(C_PTR) :: f_ptr

    if (main_prc) print *, 'reading main grid....'
    
    call make_dsetname_main_grid(dsetname, rank)  
    
    CALL h5open_f (error)
    CALL h5fopen_f (trim(adjustl(dir_grid))//trim(adjustl(grid_file)), H5F_ACC_RDWR_F, file_id, error)
    
    if (error /= 0.) then 
       if (main_prc) print *,  'main grid file not found'
       call stop_prc
    endif

    do i=0,narr_grid-1
    
    CALL h5dopen_f(file_id, dsetname(i), dset_id, error)
    CALL h5dget_space_f(dset_id, dspace_id, error) 
    call h5sget_simple_extent_npoints_f(dspace_id, num, error)

    if (i == 0) then
       tot_ncell=num
       allocate(ncell(0:tot_ncell-1), cchild(0:tot_ncell-1), cindex(0:tot_ncell-1),lvl(0:tot_ncell-1), ccoord(3,0:tot_ncell-1) & 
, csize(0:tot_ncell-1),dens(0:tot_ncell-1),dens_stars(0:tot_ncell-1) &
, u_fest(0:tot_ncell-1), u_final(0:tot_ncell-1),src_cell(0:tot_ncell-1),lock_cell(0:tot_ncell-1))
       allocate(dens_ref(0:tot_ncell-1), dens_stars_ref(0:tot_ncell-1), dens_arr(0:lnum-1, 0:tot_ncell-1),dens_stars_arr(0:lnum-1,0:tot_ncell-1), u_fest_arr(0:lnum-1,0:tot_ncell-1), u_final_arr(0:lnum-1,0:tot_ncell-1))

    endif

    CALL  make_dims(dsetname(i),rank(i))

     if (dsetname(i) == 'ncell' ) then 
     CALL h5dread_f(dset_id,H5T_NATIVE_INTEGER , ncell, dims, error)
  else if (dsetname(i) == 'ccoord' ) then
     CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , ccoord, dims, error)
  else if (dsetname(i) == 'cchild' ) then
     CALL h5dread_f(dset_id,H5T_NATIVE_INTEGER , cchild, dims, error)   
  else if (dsetname(i) == 'cindex' ) then
     !CALL h5dread_f(dset_id,h5kind_to_type(int64,H5_INTEGER_KIND), cindex, dims, error)
     f_ptr=c_loc(cindex(0))
     CALL h5dread_f(dset_id,h5kind_to_type(int64, H5_INTEGER_KIND), f_ptr,  error)     
  else if (dsetname(i) == 'lvl' ) then
     CALL h5dread_f(dset_id,H5T_NATIVE_INTEGER , lvl, dims, error)
  else if (dsetname(i) == 'csize' ) then
     CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , csize, dims, error)
  else if (dsetname(i) == 'dens' ) then ! the dens_ref read here is overwritten in read_lambda_grid if use_lambda_grid is set. 
     CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , dens_ref, dims, error)
  else if (dsetname(i) == 'dens_stars' ) then
     CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , dens_stars_ref, dims, error)
  else if (dsetname(i) == 'base' ) then
     CALL h5dread_f(dset_id,H5T_NATIVE_INTEGER , base, dims, error)   
  endif

    CALL h5sclose_f(dspace_id, error)
    CALL h5dclose_f(dset_id, error)

    if (main_prc) print *, ' '
    if (main_prc) print *, trim(dsetname(i)), ' loaded'
    if (main_prc) print *, 'dimensions = ', dims 

    deallocate(dims)

    end do 
    
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

    u_fest = 0
    u_fest_arr = 0
    u_final = 0
    u_final_arr = 0
    dens = 0
    dens_arr = 0
    dens_stars = 0
    dens_stars_arr = 0
    lock_cell = 0
    src_cell = -1 ! not zero
    max_lvl=maxval(lvl)
    tot_spare_cells = 0
    iterations = 0  ! this is 0, so it becomes 1 at the first scattering iteration (the 0 is for the direct light phase)
    iterations_dustem = 0  

    if (main_prc) print *, 'max nesting level =', max_lvl

    call print_done 
    
  end subroutine read_main_grid


  !> This subroutine reads the lambda grid from grid_file_lambda(). In this code version, only the stellar emission luminosity array dens_stars_arr() is uploaded. The dust density array dens_arr() is derived later by scaling the reference dens() array in the main grid. 
  subroutine read_lambda_grid
    INTEGER(HID_T) :: file_id                            ! File identifier
    INTEGER(HID_T) :: dset_id       ! Dataset identifier
    INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
    integer,parameter  :: narr_grid_lambda=2
    CHARACTER (LEN=lcar) :: dsetname(0:narr_grid_lambda-1) ! Dataset name
    INTEGER     ::   rank(0:narr_grid_lambda-1)   
    INTEGER(HSIZE_T) :: num
    INTEGER     ::   error  ! Error flag
    integer :: i
    integer :: iw
    logical :: found_ref_grid

    if (main_prc) print *, 'Reading lambda grids...'

    found_ref_grid = .FALSE.
    
    dsetname(0)='dens'   ; rank(0)=1
    dsetname(1)='dens_stars'  ; rank(1)=1
    !dsetname(0)='dens_stars'   ; rank(0)=1
    
    do iw = 0, lnum-1

        if (main_prc) print *, '-- lambda =', lambda_arr(iw)

    CALL h5open_f (error)
    CALL h5fopen_f (trim(adjustl(dir_grid))//trim(adjustl(grid_file_lambda_arr(iw))), H5F_ACC_RDWR_F, file_id, error)
    
    if (error /= 0.) then 
       print *, 'lambda grid file not found'
       call stop_prc
    endif

    do i=0,narr_grid_lambda-1

       call make_dims(dsetname(i),rank(i))
    
       CALL h5dopen_f(file_id, dsetname(i), dset_id, error)
    
       if (dsetname(i) == 'dens' ) then
          CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , dens, dims, error)
       else if (dsetname(i) == 'dens_stars' ) then
          CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , dens_stars, dims, error)
       endif

       CALL h5dclose_f(dset_id, error)

       if (main_prc) print *, ' '
       if (main_prc) print *, trim(dsetname(i)), ' loaded'
       if (main_prc) print *, 'dimensions = ', dims 

       deallocate(dims)

    end do
    
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

    dens_arr(iw,:) = dens
    dens_stars_arr(iw,:) = dens_stars

    if (count(dens_stars < 0) > 0) then 
       print *, 'STOP(read_lambda_grid): negative value in dens_stars!'
       stop
    endif 

    if (abs(lambda_ref-lambda_arr(iw))/lambda_ref < 1E-4) then 
       dens_ref = dens
       found_ref_grid = .TRUE.
    endif
    
 enddo

 

 ! STOP if reference grid not found
 if (.not.found_ref_grid) then
    if (main_prc) then 
       print *, 'ERROR(read_lambda_grid): reference grid not found among the lambda grids! Did you inlude this wavelength in the input wavelength list ? If so, did you create the lambda grid for the reference wavelength?'
       print *, 'lambda_ref = ', lambda_ref
    endif
    call stop_prc
 endif
    

 call print_done

  end subroutine read_lambda_grid


  !> Counts the lines of an input ASCII file. IMPORTANT: it works only if there are no blank lines at the end of files. Check that before using. 
  !> @param [in] id File unit number
  !> @param [out] num Number of lines 
  subroutine count_lines(id,num)
    integer :: id, num,i,v

    i=0
    do 
       read(id,*,iostat=v)
       if (v /= 0) exit
       i=i+1
    end do
    num=i
    rewind(id)
    
    
  end subroutine count_lines

  !> Reads the lines-of-sight theta and phi angles from file_dir_out().
  subroutine read_dir_out
    integer :: i,j,id , num
    logical :: file_exists
    
    tot_ndir=0
    tot_ndir_scaspe=0

    if (.not. use_dir_out) return

    inquire(file=trim(adjustl(dir_grid))//trim(adjustl(file_dir_out)), exist=file_exists)
        
    if (file_exists) then 

       ! read direction escaping brigthness
       
       open(newunit=id,file=trim(adjustl(dir_grid))//trim(adjustl(file_dir_out)),status='old')
       if (main_prc) print *, 'reading dir_out file: ', trim(adjustl(file_dir_out))
       
       call count_lines(id,num)
       tot_ndir = num -1 ! The first line is a comment
       if (tot_ndir <= 0) then
          print *, 'something wrong while reading ', trim(adjustl(file_dir_out))
          call stop_prc
       endif
       if (main_prc) print *, 'number of external observer directions =', tot_ndir 

       allocate(dir_i_out(0:tot_ndir-1,2))

       read(id,*) 

       do i=0,tot_ndir-1 

          read(id,*) (dir_i_out(i,j), j=1,2)   !! 1=theta 2=phi
          if (main_prc) print *, dir_i_out(i,:)
       enddo
       close(id)

       if (rt_algorithm_ID /= rta_i_obs .and. rt_algorithm_ID /= rta_i_obs_dust) then 
      
          tot_ndir_scaspe=tot_ndir
          
       else 

          call find_out_tot_ndir_scaspe

       endif

    else 
       print *, 'ERROR: file_dir_out does not exist.'
       CALL STOP_PRC
       
    endif

    call print_done

  end subroutine read_dir_out


!> This subroutine reads the internal observer positions from file_pos_obs().  
  subroutine read_pos_obs
    integer :: i,j,id , num
    logical :: file_exists
    
    tot_ndir_in=0

    if (.not. use_pos_obs) return

    if (main_prc) print *, 'reading pos_obs file....' 
    
    inquire(file=trim(adjustl(dir_grid))//trim(adjustl(file_pos_obs)), exist=file_exists)
    
    if (file_exists) then 
       open(newunit=id,file=trim(adjustl(dir_grid))//trim(adjustl(file_pos_obs)),status='old')
       if (main_prc) print *, 'reading pos obs file: ', file_pos_obs
       call count_lines(id,num)
       tot_ndir_in = num -1 ! The first line is a comment
       if (tot_ndir_in <= 0) then
          print *, 'something wrong while reading ', file_pos_obs
          call stop_prc
       endif
       if (main_prc) print *, 'number of internal observer positions =', tot_ndir_in 

       allocate(ccoord_obs(3,0:tot_ndir_in-1)) ! first index is cartesian coordinates id and the second is the observer position id. This order because it is similar to ccoord array (3D cell array) 

       read(id,*)

       do i=0,tot_ndir_in-1 

          read(id,*) (ccoord_obs(j,i), j=1,3)   !! x , y, z 
          if (main_prc) print *, (ccoord_obs(j,i), j=1,3)
       enddo
       close(id)
    else 
       print *, 'ERROR(read_pos_obs): POS_OBS FILE DOES NOT EXIST'
       CALL STOP_PRC
   
    endif

    call print_done

  end subroutine read_pos_obs

 !> This subroutine reads the stellar point source positions inside the model from file_p_src(). 
  subroutine read_p_src
  
    integer :: i,j,id , num
    logical :: file_exists

    tot_p_src=0
    if (.not. use_p_src) return

    inquire(file=trim(adjustl(dir_grid))//trim(adjustl(file_p_src)), exist=file_exists)

    if (file_exists) then 

       open(newunit=id,file=trim(adjustl(dir_grid))//trim(adjustl(file_p_src)),status='old')
       if (main_prc) print *, 'reading point source file: ', trim(adjustl(file_p_src))

       call count_lines(id,num)
       tot_p_src = num -1 ! The first line is a comment
       if (tot_p_src <= 0) then
          print *, 'something wrong while reading ', file_p_src
          call stop_prc
       endif
       if (main_prc) print *, 'number stellar source positions =', tot_p_src 

       allocate(ccoord_p_src(3,0:tot_p_src-1),cell_src(0:tot_p_src-1)) ! first index is cartesian coordinates id and the second is the observer position id. This order because it is similar to ccoord array (3D cell array) 

       read(id,*)

       do i=0,tot_p_src-1 

          read(id,*) (ccoord_p_src(j,i), j=1,3)   !! x , y, z 
          if (main_prc) print *, (ccoord_p_src(j,i), j=1,3)
       enddo
       close(id)

       allocate(lum_p_src(0:tot_p_src-1), lum_p_src_ref(0:tot_p_src-1), lum_p_src_arr(0:lnum-1,0:tot_p_src-1))

    else 
       print *, 'ERROR: STELLAR POINT SOURCE POSITION FILE DOES NOT EXIST'
       call stop_prc
    endif

    call print_done

  end subroutine read_p_src  

 !> This subroutine prepares the arrays csize_arr(), carea_arr(), and cvol_arr(). It also sets modelsize() and pabs_max(). 
  subroutine make_csize_arr
    real(kind=real64) :: cellsize
    integer :: i

    if (main_prc) print *, 'Setting cell sizes...'
    
    ! create carea_arr and cvol_arr arrays
    allocate(csize_arr(0:max_lvl),carea_arr(0:max_lvl),cvol_arr(0:max_lvl))
    modelsize=csize(0)
    pabs_max=modelsize*10   ! this is a factor involved in the calculation of the average ray path length from each cell 


    do i=0,max_lvl
       if (i > 0) then 
          cellsize=modelsize*(baseinv(1)*baseinv(2)**(i-1))  
       else 
          cellsize=modelsize
       endif
       csize_arr(i)=cellsize
       carea_arr(i)=cellsize**2
       cvol_arr(i)=cellsize**3

    enddo

    call print_done


  end subroutine make_csize_arr

  !> This subroutine prepares dims() array for reading or writing HDF5 files.
  !> @param [in] dsetname Dsetname of the array in the HDF5 file
  !> @param [in] rank Rank of the array in the HDF5 file 
  subroutine make_dims(dsetname,rank)
    
    CHARACTER (LEN=lcar) :: dsetname  ! Dataset name
    INTEGER     ::   rank
    
    if ((rank == 1).and.(dsetname /= 'base')) then
       allocate(dims(1))
       dims=tot_ncell   
    elseif (dsetname =='base') then 
       allocate(dims(1))
       dims=2    
    else if (rank == 2) then
       allocate(dims(2))
       dims=(/3,tot_ncell/)
    endif

  end subroutine make_dims

  !> This subroutine sets the output file names. It uses the labels specified in label_model_out() and label_wave_arr().
  subroutine set_filenames
    integer :: i
    character(len = 20) :: chext
    integer :: i0

    if (rt_type == rtt_grid_init_dust) then
       if (main_prc) print *, 'setting output file names for dust emission RT....'
       chext = '_dust_'//trim(adjustl(dust_heating_type)) ! name extension for dust RT output files
       call deallocate_filename_arr  ! filename_arr already allocated previously for the stellar emission RT output files. 
       i0 = i_lambda_dust(0)
    elseif (rt_type == rtt_grid_init_stars) then
       if (main_prc) print *, 'setting output file names for stellar emission RT....'
       chext = ''
       i0 = 0
    elseif (rt_type == rtt_grid_init_projection) then 
       if (main_prc) print *, 'setting output file names for projection algorithm....'
       select case(param_to_project)
       case('stellar_emission')
          chext = '_stars_em'
       case('optical_depth')
          chext = '_opt_depth'
       case default 
          if (main_prc) then 
             print *, 'ERROR(set_filenames): param_to_project not recognized!'
             print *, 'param_to_project =', param_to_project
          endif
          call stop_prc
       end select
       call deallocate_filename_arr
       i0 = 0
    else 
       print *, 'ERROR(set_filenames): rt_type not recognized here!'
       print *, 'rt_type =', rt_type
       STOP
    endif


    file_info='grid_'//trim(adjustl(label_model_out))//'_info'//trim(adjustl(chext))//'.dat'
    file_psel_av_part2='grid_'//trim(adjustl(label_model_out))//'_psel_av_part2'//trim(adjustl(chext))//'.h5'
    file_psel_av='grid_'//trim(adjustl(label_model_out))//'_psel_av'//trim(adjustl(chext))//'.h5'

    file_sed_arr_dir = 'grid_'//trim(adjustl(label_model_out))//'_sed_dir'//trim(adjustl(chext))//'.h5'
    file_sed_arr = 'grid_'//trim(adjustl(label_model_out))//'_sed'//trim(adjustl(chext))//'.h5'

    file_lum_lost_part2='grid_'//trim(adjustl(label_model_out))//'_lum_lost_part2'//trim(adjustl(chext))//'.h5'
    file_lum_lost='grid_'//trim(adjustl(label_model_out))//'_lum_lost'//trim(adjustl(chext))//'.h5'
    file_maps_part2 = 'grid_'//trim(adjustl(label_model_out))//'_maps_part2'//trim(adjustl(chext))//'.h5'
    file_maps = 'grid_'//trim(adjustl(label_model_out))//'_maps'//trim(adjustl(chext))//'.h5'
    file_maps_in_part2 = 'grid_'//trim(adjustl(label_model_out))//'_maps_in_part2'//trim(adjustl(chext))//'.h5'
    file_maps_in = 'grid_'//trim(adjustl(label_model_out))//'_maps_in'//trim(adjustl(chext))//'.h5'
    

    call allocate_filename_arr

    
    do i = 0, lnum-1 
    
       file_u_fest_part1_arr(i)='grid_'//trim(adjustl(label_model_out))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_u_fest_part1'//trim(adjustl(chext))//'.h5'    
       file_ufield_part2_arr(i)='grid_'//trim(adjustl(label_model_out))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_ufield_part2'//trim(adjustl(chext))//'.h5'
       file_scaspe_part2_arr(i)='grid_'//trim(adjustl(label_model_out))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_scaspe_part2'//trim(adjustl(chext))//'.h5'
       file_i_obs_part2_arr(i)='grid_'//trim(adjustl(label_model_out))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_i_obs_part2'//trim(adjustl(chext))//'.h5'
       file_i_obs_in_part2_arr(i)='grid_'//trim(adjustl(label_model_out))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_i_obs_in_part2'//trim(adjustl(chext))//'.h5'
       file_i_obs_arr(i)='grid_'//trim(adjustl(label_model_out))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_i_obs'//trim(adjustl(chext))//'.h5'
       file_i_obs_in_arr(i)='grid_'//trim(adjustl(label_model_out))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_i_obs_in'//trim(adjustl(chext))//'.h5'
       file_ufield_arr(i)='grid_'//trim(adjustl(label_model_out))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_ufield'//trim(adjustl(chext))//'.h5'
       !file_llost_part2_arr(i)='grid_'//trim(adjustl(label_model_out))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_llost_part2'//trim(adjustl(chext))//'.dat'
       !file_llost_arr(i)='grid_'//trim(adjustl(label_model_out))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_llost'//trim(adjustl(chext))//'.dat'
       file_scaspe_tot_arr(i)='grid_'//trim(adjustl(label_model_out))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_scaspe_tot'//trim(adjustl(chext))//'.h5'
  
    end do
    

    call print_done 

  end subroutine set_filenames

!> Changes i_obs() array filename for i_obs RT algorithm. 
subroutine reset_filenames 
  integer :: i, i0
  character(len = 20) :: chext

if (main_prc) print *, 'resetting output filenames for outgoing specific intensity....'

if (rt_algorithm_ID == rta_i_obs_dust) then 
   chext = '_dust_'//trim(adjustl(dust_heating_type)) ! name extension for dust RT output files
   i0 = i_lambda_dust(0)
elseif (rt_algorithm_ID == rta_i_obs) then
   chext = ''
   i0 = 0
else 
   print *, 'ERROR(reset_filenames): rt_algorithm_ID not recognized here!'
   print *, 'rt_algorithm_ID =', rt_algorithm_ID
   STOP
endif


do i = 0, lnum -1 
   file_i_obs_part2_arr(i)='grid_'//trim(adjustl(label_model_out_i_obs))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_i_obs_part2'//trim(adjustl(chext))//'_int.h5'
   file_i_obs_arr(i)='grid_'//trim(adjustl(label_model_out_i_obs))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_i_obs'//trim(adjustl(chext))//'_int.h5'
   file_i_obs_in_part2_arr(i)='grid_'//trim(adjustl(label_model_out_i_obs))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_i_obs_in_part2'//trim(adjustl(chext))//'_int.h5'
   file_i_obs_in_arr(i)='grid_'//trim(adjustl(label_model_out_i_obs))//'_l'//trim(adjustl(label_wave_arr(i+i0)))//'um_i_obs_in'//trim(adjustl(chext))//'_int.h5'
enddo

file_sed_arr_dir = 'grid_'//trim(adjustl(label_model_out_i_obs))//'_sed_dir'//trim(adjustl(chext))//'_int.h5'
file_sed_arr = 'grid_'//trim(adjustl(label_model_out_i_obs))//'_sed'//trim(adjustl(chext))//'_int.h5'
file_maps_part2 = 'grid_'//trim(adjustl(label_model_out_i_obs))//'_maps_part2'//trim(adjustl(chext))//'_int.h5'
file_maps = 'grid_'//trim(adjustl(label_model_out_i_obs))//'_maps'//trim(adjustl(chext))//'_int.h5'
file_maps_in_part2 = 'grid_'//trim(adjustl(label_model_out_i_obs))//'_maps_in_part2'//trim(adjustl(chext))//'_int.h5'
file_maps_in = 'grid_'//trim(adjustl(label_model_out_i_obs))//'_maps_in'//trim(adjustl(chext))//'_int.h5'


call print_done

end subroutine reset_filenames


  !> Allocates arrays containing the output filenames
  subroutine allocate_filename_arr

    allocate(file_u_fest_part1_arr(0:lnum-1), file_ufield_part2_arr(0:lnum-1), file_scaspe_part2_arr(0:lnum-1), file_i_obs_part2_arr(0:lnum-1), file_i_obs_in_part2_arr(0:lnum-1), file_i_obs_arr(0:lnum-1), file_i_obs_in_arr(0:lnum-1), file_ufield_arr(0:lnum-1), file_scaspe_tot_arr(0:lnum-1))

  end subroutine allocate_filename_arr

  !> Deallocates arrays containing the output filenames
  subroutine deallocate_filename_arr

    deallocate(file_u_fest_part1_arr, file_ufield_part2_arr, file_scaspe_part2_arr, file_i_obs_part2_arr, file_i_obs_in_part2_arr, file_i_obs_arr, file_i_obs_in_arr, file_ufield_arr,  file_scaspe_tot_arr)

  end subroutine deallocate_filename_arr



  
  !> This subroutine prints the used RT parameters in file_info().
  subroutine write_file_info

    integer :: unit

    if (main_prc) then

       print *, 'Writing RT info file...' 
    
    open(newunit=unit, file=trim(adjustl(dir_runs))//file_info, status='unknown')
  
    write(unit,*) ' RADIATION TRANSFER PARAMETERS'
    write(unit,*) 'kp_sca_max=', kp_sca_max
    write(unit,*) 'rad_lim=', rad_lim
    write(unit,*) 'accuracy=', accuracy
    write(unit,*) 'conv_en_lim=',conv_en_lim
    write(unit,*) 'bm_par=',bm_par
    write(unit,*) 'bm_par_sca=',bm_par_sca
    write(unit,*) 'bm_par_max=', bm_par_max
    write(unit,*) 'tau_cell_max= ', tau_cell_max
    write(unit,*) 'N processor =', nproc
    write(unit,*) 'lambda_arr=',lambda_arr
    write(unit,*)
    write(unit,*) ' OPACITY PARAMETERS' 
    write(unit,*) 'kabs_arr=',kabs_arr
    write(unit,*) 'ksca_arr=',ksca_arr
    write(unit,*) 'gsca_arr=',gsca_arr
    write(unit,*)
    write(unit,*) ' OTHER PARAMETERS'
    write(unit,*) ' cs =', cs

    close(unit)

 endif

 call print_done

end subroutine write_file_info

!> This routine checks for output file existence. If any of the final output file exists, the program exists with an error message. If all of the intermediate files exist, the program restores them if restore_file_mpi() = .TRUE.. Then, the program stars directly from the direct light processing or the scattering iterations. Note that this is used only for the stellar emission RT. In the case of the dust emission RT, this is more complicated because of the dust heating iterations and it is not implemented for the moment. 
subroutine check_files
  
  
  logical :: file_exists
  logical, allocatable :: file_exists_arr(:,:)
  character(LEN=10) :: reply
  integer :: i
     
     file_restore = .FALSE.
     file_restore_part1 = .FALSE.

     ! if dust algorithm, check only at the first iteration
     if (rt_algorithm_ID == rta_dust .or. rt_algorithm_ID == rta_dust2D) then
        if (iterations_dustem > 1) return
     endif
     
     if (main_prc) print *, 'checking existence output files....'
     
!!! check u_field file

     do i = 0, lnum-1 
     
     inquire(file=trim(adjustl(dir_runs))//file_ufield_arr(i), exist=file_exists)  
     if (file_exists) then
        if (main_prc) then
           print *, trim(adjustl(dir_runs))//file_ufield_arr(i), ' already exists! '
           call stop_prc
        endif
     endif

!!! check i_obs file
     if (use_dir_out) then 
        inquire(file=trim(adjustl(dir_runs))//file_i_obs_arr(i), exist=file_exists)
        if (file_exists) then
           if (main_prc) then 
              print *, trim(adjustl(dir_runs))//file_i_obs_arr(i), ' already exists! '
              call stop_prc
           endif
        endif
     endif
  

!!! check i_obs_in file
     if (use_pos_obs) then 
        inquire(file=trim(adjustl(dir_runs))//file_i_obs_in_arr(i), exist=file_exists)
        if (file_exists) then
           if (main_prc) then 
              print *, trim(adjustl(dir_runs))//file_i_obs_in_arr(i), ' already exists! '
              call stop_prc
           endif
        endif
     endif

!!! check scaspe_tot file
     if (print_scaspe_tot) then 
        inquire(file=trim(adjustl(dir_runs))//file_scaspe_tot_arr(i), exist=file_exists)
        if (file_exists) then
           if (main_prc) then
              print *, trim(adjustl(dir_runs))//file_scaspe_tot_arr(i), ' already exists! '
              call stop_prc
           endif
        endif
     endif

  end do

  !! check psel_av file
     if (print_psel_av) then 
        inquire(file=trim(adjustl(dir_runs))//file_psel_av, exist=file_exists)
        if (file_exists) then
           if (main_prc) then 
              print *, trim(adjustl(dir_runs))//file_psel_av, ' already exists! '
              call stop_prc
           endif
        endif
     endif

     !! check lum_lost_file
     inquire(file=trim(adjustl(dir_runs))//file_lum_lost, exist=file_exists)
     if (file_exists) then
        if (main_prc) then 
           print *, trim(adjustl(dir_runs))//file_lum_lost, ' already exists! '
           call stop_prc
        endif
     endif

     if (main_prc) print *, 'Final output files not found'
     call print_done

     ! return in case of the dust emission RT algorithms (no check of intermediate files)
     if (rt_algorithm_ID == rta_dust .or. rt_algorithm_ID == rta_dust2d) return

     ! start checking output part2 in case of stellar emission RT 
     if (main_prc) print *, 'checking existence output part2....'
     
!!! Check intermediate file existence
     allocate(file_exists_arr(0:lnum-1, 6))

     do i = 0, lnum -1 
     file_exists_arr = .TRUE.
     inquire(file=trim(adjustl(dir_runs))//file_ufield_part2_arr(i), exist=file_exists)
     file_exists_arr(i,1)=file_exists
     inquire(file=trim(adjustl(dir_runs))//file_scaspe_part2_arr(i), exist=file_exists)
     file_exists_arr(i,2)=file_exists

     
     if (use_dir_out) then 
        inquire(file=trim(adjustl(dir_runs))//file_i_obs_part2_arr(i), exist=file_exists)
        file_exists_arr(i,3)=file_exists
     endif
     if (use_pos_obs) then 
        inquire(file=trim(adjustl(dir_runs))//file_i_obs_in_part2_arr(i), exist=file_exists)
        file_exists_arr(i,4)=file_exists
     endif

  end do

  if (print_psel_av) then 
     inquire(file=trim(adjustl(dir_runs))//file_psel_av_part2, exist=file_exists)
     file_exists_arr(:,5)=file_exists  ! this is a trick
  endif

  inquire(file=trim(adjustl(dir_runs))//file_lum_lost_part2, exist=file_exists)
  file_exists_arr(:,6)=file_exists  ! this is a trick
     
  if (all(file_exists_arr)) then
     deallocate(file_exists_arr)
!!$     if (np_mpi == 1) then 
!!$           
!!$        print *, 'All intermediate files exist. Restore them ? (yes/no)'
!!$        do 
!!$           read(*,*) reply
!!$           select case(reply)
!!$           case('yes')                   
!!$              file_restore=.TRUE.
!!$              exit
!!$           case('no')
!!$              exit
!!$           case default
!!$              if (main_prc) print *, 'yes or no ?'
!!$           end select
!!$        enddo
!!$              
!!$     else 

     if (restore_file_mpi) then
        file_restore = .TRUE.           
     endif
                
!!$     endif

     if (file_restore) then 
        if (main_prc) print *, 'Intermediate files will be restored'
        call print_done
        return
     endif
     
  endif

  if (main_prc) print *, 'Output files of part 2 do not exist or do not have to be restored.'
  call print_done
        
  ! check part 1 file existence
  if (main_prc) print *, 'Checking existence output part 1...'
            
  if (allocated(file_exists_arr)) deallocate(file_exists_arr)
  allocate(file_exists_arr(0:lnum-1,0:0))

  file_exists_arr = .TRUE.
        
  do i = 0, lnum-1 
        
     inquire(file=trim(adjustl(dir_runs))//file_u_fest_part1_arr(i), exist=file_exists)
     file_exists_arr(i,0) = file_exists
  enddo
        
  if (all(file_exists_arr)) then
     deallocate(file_exists_arr)
!!$     if (np_mpi == 1) then               
!!$        print *, 'Precalculation output files exist. Restore them ? (yes/no)'
!!$        do 
!!$           read(*,*) reply
!!$           select case(reply)
!!$           case('yes')  
!!$              file_restore_part1 = .TRUE.
!!$              exit
!!$           case('no')
!!$              exit
!!$           case default
!!$              if (main_prc) print *, 'yes or no ?'
!!$           end select
!!$        enddo
!!$           
!!$     else 

     if (restore_file_mpi) file_restore_part1 = .TRUE.
        
!!$     endif

     if (file_restore_part1) then 
        if (main_prc) print *, 'Output part 1 will be restored!'
        call print_done
        return
     endif
     
  endif
           
if (allocated(file_exists_arr)) deallocate(file_exists_arr)
  if (main_prc) print *, 'Output part 1 not found or not to be restored.'

  call print_done

end subroutine check_files

!>  This subroutine is used to print the output files.
!> @param [in] filename File name
!> @param [in] dsetname Dsetname of the array to be printed
!> @param [in] dims Array or scalar with the array dimensions
!> @param [in] rank Rank of the array to be written on the HDF5 file
subroutine print_big_array(filename,dsetname,dims,rank)
  
     CHARACTER (LEN=lcar) :: filename ! File name
     CHARACTER (LEN=lcar) :: dsetname ! Dataset name

     INTEGER(HID_T) :: file_id                            ! File identifier
     INTEGER(HID_T) :: dset_id       ! Dataset identifier
     INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
 
     INTEGER(HSIZE_T), DIMENSION(*)  :: dims  ! Dataset dimensions
    
     INTEGER     ::   rank                         ! Dataset rank
     INTEGER     ::   error  ! Error flag
     integer :: i, j 

   !
!    Initialize FORTRAN interface.
     CALL h5open_f (error)
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error) 
     CALL h5screate_simple_f(rank, dims, dspace_id, error)  
     CALL h5dcreate_f(file_id, dsetname,H5T_NATIVE_DOUBLE , dspace_id, &
       dset_id, error)
     !CALL h5dopen_f(file_id, dsetname, dset_id, error)

     select case(dsetname)
     case('u_final')
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , u_final, dims, error)
     case('i_obs') 
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , i_obs, dims, error)
     case('i_obs_in' )
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , i_obs_in, dims, error) 
     case('scaspe') 
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , scaspe, dims, error)
     case('scaspe_tot')
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , scaspe_tot, dims, error)
     case('u_fest') 
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , u_fest, dims, error)
     case('psel_av_arr') 
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , psel_av_arr, dims, error)
     case('sed_arr_dir') 
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , sed_arr_dir, dims, error)
     case('sed_arr') 
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , sed_arr, dims, error)
     case('lum_lost')
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , lum_lost, dims, error)
     case('map_arr_out')
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , map_arr_out, dims, error)
     case('map_in_arr_out')
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , map_in_arr_out, dims, error)
     case default 
        print *, 'ERROR(print_big_array): dsetname not recognized!'
        stop
     end select

     CALL h5dclose_f(dset_id, error)

     CALL h5sclose_f(dspace_id, error)

     ! flush buffer
     !CALL h5fflush_f(file_id, H5F_SCOPE_LOCAL_F, error)
     
     CALL h5fclose_f(file_id, error)
     
     CALL h5close_f(error)

     
end subroutine print_big_array

!> Adds extra info to printed big array. This works only for one dimensional arrays.   
subroutine add_info_big_array(filename,dsetname_add, dims_add, rank_add)

  INTEGER(HID_T) :: file_id                            ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
  integer :: narr
  INTEGER     ::   rank_add(0:), rank_in                   ! Dataset rank
  INTEGER     ::   error  ! Error flag
  INTEGER(HSIZE_T) ::  dims_add(0:),dims_in(0:0)
  character(LEN=lcar) :: dsetname_add(0:), filename
  integer :: i

  narr = size(dsetname_add)

  CALL h5open_f (error)

  CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)

  do i=0,narr-1
            
     dims_in = dims_add(i)
     rank_in = rank_add(i)
     CALL h5screate_simple_f(rank_in, dims_in, dspace_id, error)
     
     if (dsetname_add(i) /= 'nside_sca') then  
        CALL h5dcreate_f(file_id, dsetname_add(i),H5T_NATIVE_DOUBLE , dspace_id, &
          dset_id, error)   
     else 
        CALL h5dcreate_f(file_id, dsetname_add(i),H5T_NATIVE_INTEGER , dspace_id, &
          dset_id, error) 
    endif
        
     if (dsetname_add(i) == 'nside_sca' ) then
        CALL h5dwrite_f(dset_id,H5T_NATIVE_INTEGER , nside_sca, dims_in, error)  
     elseif (dsetname_add(i) == 'lambda_arr_maps' ) then
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , lambda_arr_maps, dims_in, error)   
     elseif (dsetname_add(i) == 'size_map' ) then
        CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , size_map, dims_in, error)      
        
     endif
     
     CALL h5dclose_f(dset_id, error)

     CALL h5sclose_f(dspace_id, error)
     
  end do
      
  CALL h5fclose_f(file_id, error)     

  ! flush buffer
  !CALL h5fflush_f(file_id, H5F_SCOPE_LOCAL_F, error)
  
  CALL h5close_f(error)

  if (error /= 0.) then 
     print *, 'HDF5 error =', error 
     STOP '(add_info_big_array) something did not work with HDF5 file writing'
  endif

end subroutine add_info_big_array


!> Writes u_fest_arr() on files file_u_fest_part1_arr(). 
subroutine print_u_fest_arr
  integer :: i,j,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname
   
  dsetname='u_fest'
  rank=1
  allocate(dims(rank))
  dims=tot_ncell  

  do i = 0, lnum -1 
     u_fest = u_fest_arr(i,:)
     filename=trim(adjustl(dir_runs))//file_u_fest_part1_arr(i)      
     call print_big_array(filename,dsetname,dims,rank)
  enddo
  deallocate(dims)

end subroutine print_u_fest_arr

!> Reads u_fest_arr from previously written files. 
subroutine read_u_fest_arr
  integer :: i,j,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname

  dsetname='u_fest'
  rank=1
  allocate(dims(rank))
  dims=tot_ncell
     
  do i = 0, lnum -1 
     filename=trim(adjustl(dir_runs))//file_u_fest_part1_arr(i)      
     call read_big_array(filename,dsetname,dims)
     u_fest_arr(i,:) = u_fest
  enddo
           
  deallocate(dims)

end subroutine read_u_fest_arr

!> Writes ufield_arr() on files file_ufield_part2_arr() and file_ufield_arr().
subroutine print_ufield_arr
  integer :: i,j,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname
  
  dsetname='u_final'
  rank=1
  allocate(dims(rank))
  dims=tot_ncell
  
  do i = 0, lnum-1
     u_final = u_final_arr(i,:)
     if (rt_type == rtt_output_part2) then 
        filename=trim(adjustl(dir_runs))//file_ufield_part2_arr(i)
     elseif (rt_type == rtt_scatt) then
        filename=trim(adjustl(dir_runs))//file_ufield_arr(i)
     endif
     call print_big_array(filename,dsetname,dims,rank)
  enddo
  deallocate(dims)
  
end subroutine print_ufield_arr

!> Reads ufield_arr from the previously written files. 
subroutine read_ufield_arr
  integer :: i
  integer :: ierr
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname, dir_tmp

  dsetname='u_final'
  rank=1
  allocate(dims(rank))
  dims=tot_ncell
  dir_tmp = ''

  do i = 0, lnum -1
     if (rt_type == rtt_start) then 
        filename=trim(adjustl(dir_runs))//file_ufield_part2_arr(i)
     elseif (rt_type == rtt_read_ufield .or. rt_type == rtt_grid_init_dust) then
        filename=trim(adjustl(dir_runs))//file_ufield_arr(i)
     endif
     call check_file_existence(dir_tmp, filename) 
     call read_big_array(filename,dsetname,dims)
     u_final_arr(i,:) = u_final
  enddo
  
  deallocate(dims)

end subroutine read_ufield_arr


!> Writes psel_av_arr() on file_psel_av_part2() and file_psel_av().
subroutine print_psel_av_arr
  integer :: i,j,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname
   
  dsetname='psel_av_arr'
  rank=3
  allocate(dims(rank))
  if (rt_type == rtt_output_part2) then 
     filename=trim(adjustl(dir_runs))//file_psel_av_part2
     dims=(/1,1,tot_ncell+tot_p_src/)  
  elseif (rt_type == rtt_scatt) then
     filename=trim(adjustl(dir_runs))//file_psel_av
     if (iterations_dustem /= 0) then 
        dims=(/iterations_dustem,iterations,tot_ncell+tot_p_src/) ! note: here the second dimension is iterations and not iterations+1 (the +1 is to include the direct light phase as well) because +1 is added in rt_prepare before stopping the scattering iterations.   
     else
        dims=(/1,iterations,tot_ncell+tot_p_src/)
     endif

  endif
  call print_big_array(filename,dsetname,dims,rank)
  deallocate(dims)

end subroutine print_psel_av_arr

!> Reads psel_av_arr() from previously written files. 
subroutine read_psel_av_arr
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname

  if (rt_type == rtt_start) then 
     filename=trim(adjustl(dir_runs))//file_psel_av_part2
  endif
  
  dsetname='psel_av_arr'
  rank=3
  allocate(dims(rank))
  dims=(/1,iterations,tot_ncell+tot_p_src/)
  call read_big_array(filename,dsetname,dims)
  deallocate(dims)

end subroutine read_psel_av_arr

!> Writes i_obs_arr() on file_i_obs_part2_arr() and file_i_obs_arr().
subroutine print_i_obs_arr
  integer :: i,j,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname
  integer :: i0,i1

  call set_i_opacity_arrays(i0,i1)

  dsetname='i_obs'
  rank=2
  allocate(dims(rank))
  dims=(/tot_ncell+tot_p_src,tot_ndir/)  

  do i = 0, lnum -1
     if (no_communications .and. .not. main_prc) exit !if no communications, only main MPI process makes output files.
     if (count(ind_i_obs == i+i0) == 0) cycle
     if (iq_sca_node(i)) then
        call set_wavelength_index(i,k)
        i_obs = i_obs_arr(k,:,:)
        if (rt_type == rtt_output_part2) then 
           filename=trim(adjustl(dir_runs))//file_i_obs_part2_arr(i)
        elseif (rt_type == rtt_i_obs) then
           filename=trim(adjustl(dir_runs))//file_i_obs_arr(i)
        endif

        call print_big_array(filename,dsetname,dims,rank)
     endif
  enddo

  deallocate(dims)

end subroutine print_i_obs_arr

!> Reads i_obs_arr() from previously written files. 
subroutine read_i_obs_arr
  integer :: i,k
  integer :: ierr
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname
 
  dsetname='i_obs'
  rank=2
  allocate(dims(rank))
  dims=(/tot_ncell+tot_p_src,tot_ndir/)
  
  do i = 0, lnum-1 
     if (iq_sca_node(i)) then
        call set_wavelength_index(i,k)
        if (rt_type == rtt_start .or. rt_type == rtt_read_i_obs_part2) then 
           filename=trim(adjustl(dir_runs))//file_i_obs_part2_arr(i)
        elseif (rt_type == rtt_read_i_obs) then
           filename=trim(adjustl(dir_runs))//file_i_obs_arr(i)
        endif
        call read_big_array(filename,dsetname,dims)
        i_obs_arr(k, :,:) = i_obs
     endif
  end do
  deallocate(dims)

end subroutine read_i_obs_arr

!> Writes i_obs_in_arr() on file_i_obs_in_part2_arr() and file_i_obs_in_arr().
subroutine print_i_obs_in_arr
integer :: i,j,k
INTEGER     ::   rank 
character (LEN=lcar)  :: filename,  dsetname
integer :: i0,i1

call set_i_opacity_arrays(i0,i1)

dsetname='i_obs_in'
rank=2
allocate(dims(rank))
dims=(/tot_ncell+tot_p_src,tot_ndir_in/) 

do i = 0, lnum -1
   if (no_communications .and. .not. main_prc) exit !if no communications, only main MPI process makes output files.   
   if (count(ind_i_obs == i+i0) == 0) cycle
   if (iq_sca_node(i)) then
      call set_wavelength_index(i,k)
      i_obs_in = i_obs_in_arr(k,:,:)
      if (rt_type == rtt_output_part2) then 
         filename=trim(adjustl(dir_runs))//file_i_obs_in_part2_arr(i)
      elseif (rt_type == rtt_i_obs) then
         filename=trim(adjustl(dir_runs))//file_i_obs_in_arr(i)   
      endif
      call print_big_array(filename,dsetname,dims,rank)
   endif
end do

deallocate(dims)

end subroutine print_i_obs_in_arr

!> Reads i_obs_in_arr() from previously written files. 
subroutine read_i_obs_in_arr
  integer :: i,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname

  dsetname='i_obs_in'
  rank=2
  allocate(dims(rank))
  dims=(/tot_ncell+tot_p_src,tot_ndir_in/)
  do i = 0, lnum-1 
     if (iq_sca_node(i)) then
        call set_wavelength_index(i,k)
        if (rt_type == rtt_start .or. rt_type == rtt_read_i_obs_part2) then 
           filename=trim(adjustl(dir_runs))//file_i_obs_in_part2_arr(i)
        elseif (rt_type == rtt_read_i_obs) then 
           filename=trim(adjustl(dir_runs))//file_i_obs_in_arr(i)
        endif
        call read_big_array(filename,dsetname,dims)
        i_obs_in_arr(k,:,:) = i_obs_in
     endif
  enddo
  deallocate(dims)
end subroutine read_i_obs_in_arr


!> Writes scaspe_arr() on file_scaspe_part2_arr().
subroutine print_scaspe_arr
  integer :: i,j,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname 
  integer :: ierr
  integer :: i0, i1 
  
  call set_i_opacity_arrays(i0,i1)

  dsetname='scaspe'
  rank=2
  allocate(dims(rank))
  
  do i = 0, lnum-1
     dims=(/npix_arr(i+i0),tot_ncell/)
    
     if (no_communications .and. .not. main_prc) exit !if no communications, only main MPI process makes output files.   
     if (iq_sca_node(i)) then
        allocate(scaspe(0:dims(1)-1, 0:dims(2)-1))
        call set_wavelength_index(i,k)
        scaspe(:, :) = scaspe_arr(k)%a(:,:)
        filename=trim(adjustl(dir_runs))//file_scaspe_part2_arr(i) ! the index i is right here
        call print_big_array(filename,dsetname,dims,rank)
        deallocate(scaspe)
     endif
     if (.not. no_communications) call mpi_barrier(MPI_COMM_WORLD, ierr)  ! this is here because HDF5 file writing sometimes can give problems when too much stuff is going on.  
     
  end do

  deallocate(dims)

end subroutine print_scaspe_arr

!> Reads scaspe_arr() from previously written files. 
subroutine read_scaspe_arr

  integer :: i,j,k
  integer :: tot_el, ierr
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname
  integer :: i0, i1 
  
  call set_i_opacity_arrays(i0,i1)

  dsetname='scaspe'
  rank=2
  allocate(dims(rank))
  
  do i = 0, lnum-1 
     dims=(/npix_arr(i+i0),tot_ncell/)     
     if (iq_sca_node(i)) then
        allocate(scaspe(0:dims(1)-1, 0:dims(2)-1))
        call set_wavelength_index(i,k)
        filename=trim(adjustl(dir_runs))//file_scaspe_part2_arr(i)
        call read_big_array(filename,dsetname,dims) 
        scaspe_arr(k)%a(:,:) = scaspe(:,:)
        deallocate(scaspe)
     endif
  enddo
           
  deallocate(dims)

end subroutine read_scaspe_arr


!> Writes scaspe_tot_arr() on file_scaspe_tot_arr().
subroutine print_scaspe_tot_arr
  integer :: i,j,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname
  integer, parameter :: n_add = 1
  character (LEN=lcar)  :: dsetname_add(0:n_add-1) 
  integer :: rank_add(0:n_add-1)
  integer :: ierr
  INTEGER(HSIZE_T), allocatable ::  dims_add(:)
  integer :: i0, i1 
  
  call set_i_opacity_arrays(i0,i1)

  dsetname='scaspe_tot'
  rank=2
  allocate(dims(rank))
  
  ! additional info 
  dsetname_add = 'nside_sca'
  rank_add = 1 
  allocate(dims_add(0:rank_add(0)-1))
  dims_add = 1 
  
  do i = 0, lnum-1
     dims=(/npix_arr(i+i0),tot_ncell/)
     nside_sca = 2**kp_sca_arr(i+i0)
     
     if (no_communications .and. .not. main_prc) exit !if no communications, only main MPI process makes output files.     
     if (iq_sca_node(i)) then
        allocate(scaspe_tot(0:dims(1)-1, 0:dims(2)-1))
        call set_wavelength_index(i,k)
        scaspe_tot(:, :) = scaspe_tot_arr(k)%a(:,:)
        filename=trim(adjustl(dir_runs))//file_scaspe_tot_arr(i) ! the index i is right here
        call print_big_array(filename,dsetname,dims,rank)
        call add_info_big_array(filename,dsetname_add, dims_add, rank_add) 
        deallocate(scaspe_tot)
     endif
     if (.not. no_communications) call mpi_barrier(MPI_COMM_WORLD, ierr)  ! this is here because HDF5 file writing sometimes can give problems when too much stuff is going on. 
     
  end do

  deallocate(dims, dims_add)

end subroutine print_scaspe_tot_arr

!> Reads scaspe_tot_arr() from previously written files. 
subroutine read_scaspe_tot_arr
  integer :: i,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname
  integer :: ierr
  integer :: i0, i1 
  
  call set_i_opacity_arrays(i0,i1)
  
  !npix=12*nside_sca**2+tot_ndir_scaspe
  dsetname='scaspe_tot'
  rank=2
  allocate(dims(rank))
 
  do i = 0, lnum-1
     dims=(/npix_arr(i+i0),tot_ncell/)    
     if (iq_sca_node(i)) then
        allocate(scaspe_tot(0:dims(1)-1, 0:dims(2)-1))
        call set_wavelength_index(i,k)
        filename=trim(adjustl(dir_runs))//file_scaspe_tot_arr(i)
        call read_big_array(filename,dsetname,dims)
        scaspe_tot_arr(k)%a(:,:) = scaspe_tot
        deallocate(scaspe_tot)
     endif
  enddo
  deallocate(dims)
  
end subroutine read_scaspe_tot_arr


!> Writes sed_arr_dir() on file_sed_arr_dir().
subroutine print_sed_arr_dir
  integer :: i,j,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname 

  dsetname='sed_arr_dir'
  rank=2
  allocate(dims(rank)) 
  dims=(/lnum_tot,tot_ndir/) 
  filename = trim(adjustl(dir_runs))//file_sed_arr_dir
  call print_big_array(filename,dsetname,dims,rank)
  deallocate(dims)

end subroutine print_sed_arr_dir

!> Writes sed_arr() on file_sed_arr().
subroutine print_sed_arr
  integer :: i,j,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname 

  dsetname='sed_arr'
  rank=2
  allocate(dims(rank)) 
  dims=(/lnum_tot,tot_ndir/) 
  filename = trim(adjustl(dir_runs))//file_sed_arr
  call print_big_array(filename,dsetname,dims,rank)
  deallocate(dims)

end subroutine print_sed_arr


!> Writes map_arr_out() on file_maps() and file_maps_part2().
subroutine print_map_arr_out(filename)
  integer :: i,j,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname 
  integer, parameter :: n_add = 1
  character (LEN=lcar)  :: dsetname_add(0:n_add-1) 
  integer :: rank_add(0:n_add-1)
  INTEGER(HSIZE_T), allocatable ::  dims_add(:)

  dsetname='map_arr_out'
  rank=4
  allocate(dims(rank)) 
  dims=(/npixel_maps, npixel_maps, lnum_maps,tot_ndir/) 
  filename = trim(adjustl(dir_runs))//filename

  call print_big_array(filename,dsetname,dims,rank)

  dsetname_add='lambda_arr_maps'
  rank_add(0) = 1 
  allocate(dims_add(0:rank_add(0)-1))
  dims_add = lnum_maps

  call add_info_big_array(filename,dsetname_add, dims_add, rank_add) 
  
  dsetname_add='size_map'
  rank_add(0) = 1 
  dims_add = 1

  call add_info_big_array(filename,dsetname_add, dims_add, rank_add) 

  deallocate(dims)

end subroutine print_map_arr_out

!> Writes map_in_arr_out() on file_maps_in() and file_maps_in_part2().
subroutine print_map_in_arr_out(filename)
  integer :: i,j,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname 
  integer, parameter :: n_add = 1
  character (LEN=lcar)  :: dsetname_add(0:n_add-1) 
  integer :: rank_add(0:n_add-1)
  INTEGER(HSIZE_T), allocatable ::  dims_add(:)

  dsetname='map_in_arr_out'
  rank=3
  allocate(dims(rank)) 
  dims=(/npix_maps, lnum_maps,tot_ndir_in/) 
  filename = trim(adjustl(dir_runs))//filename

  dsetname_add='lambda_arr_maps'
  rank_add(0) = 1 
  allocate(dims_add(0:rank_add(0)-1))
  dims_add = lnum_maps

  call print_big_array(filename,dsetname,dims,rank)
  call add_info_big_array(filename,dsetname_add, dims_add, rank_add) 
  deallocate(dims)

end subroutine print_map_in_arr_out



!> This routine decides whith output files to write depening on rt_type(). 
subroutine make_output

  integer :: i,j,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname
  integer :: ierr

  !skip if dust RT and cnflag_dust /= TRUE
  if (rt_algorithm_ID == rta_dust .or. rt_algorithm_ID == rta_dust2d) then
     if (.not. cnflag_dust) return
  endif

  select case (rt_type)

  case(rtt_output_part1)
     if (print_output_part1) then
        if (main_prc) then 
           print *, 'writing output part 1...'
           call print_u_fest_arr         
        endif
        call print_done
     endif
     
  case(rtt_output_part2)

     if (print_output_part2) then 

        if (main_prc) print *, 'writing output part 2...'
        
        if (rt_algorithm_ID == rta_main .or. rt_algorithm_ID == rta_2D) then ! the following output is not printed in the dust RT algorithms 
           if (main_prc) then 
              print *, '--writing ufield part2'
              call print_ufield_arr
              
              if (print_psel_av) then 
                 print *, '--writing psel_av_part2'
                 call print_psel_av_arr
              endif
              
              print *, '--writing lum_lost file for part2'
              call print_lum_lost
           endif
           
        endif
  
        if (tot_ndir > 0 ) then 
           if (main_prc) print *, '--writing i_obs part2'
           call print_i_obs_arr
        endif

        if (tot_ndir_in > 0) then 
           print *, '--writing i_obs_in part2'
           call print_i_obs_in_arr
        endif
  
        call mpi_barrier(MPI_COMM_WORLD, ierr)  ! this is just so files are printed together by the different processes 
     
        if (print_scaspe_part2) then
           if (rt_algorithm_ID == rta_main .or. rt_algorithm_ID == rta_2D) then
              if (main_prc) print *, '--writing scaspe part2'
              call print_scaspe_arr
           endif
        endif
     endif

     if (print_sed) then
        if (main_prc .and. tot_ndir > 0) then
           print *, '--writing emission SED for direct light'
           call print_sed_arr_dir
        endif
     endif

     if (print_maps) then 
        if (main_prc .and. tot_ndir > 0) then
           print *, '--writing map_arr_out for direct light'
           call print_map_arr_out(file_maps_part2)
        endif
     endif

     if (print_maps_in) then 
        if (main_prc .and. tot_ndir_in > 0) then
           print *, '--writing map_in_arr_out for direct light'
           call print_map_in_arr_out(file_maps_in_part2)
        endif
     endif
  
     call print_done

  case (rtt_scatt)

     if (main_prc) then 
        print *, 'writing output...'
        print *, '--writing ufield'
        call print_ufield_arr
     endif

     call mpi_barrier(MPI_COMM_WORLD, ierr)  ! this is just so files are printed together by the different processes 

     if (print_scaspe_tot) then 
        if (main_prc) print *, '--writing scaspe_tot'
        call print_scaspe_tot_arr
     endif

     if (main_prc) then 
        if (print_psel_av) then        
           print *, '--writing psel_av'
           call print_psel_av_arr
        endif

        print *, '--writing lum_lost file'
        call print_lum_lost

     endif

     call print_done

  case(rtt_i_obs)

     if (main_prc) print *, 'writing output scattering iterations...'
     
     if (tot_ndir > 0) then 
        if (main_prc) print *, '--writing i_obs...'
        call print_i_obs_arr
     endif

     if (tot_ndir_in > 0) then 
        if (main_prc) print *, '--writing i_obs_in...'
        call print_i_obs_in_arr
        call print_done
     endif

     if (print_sed) then
        if (main_prc .and. tot_ndir > 0) then
           print *, '--writing total emission SED'
           call print_sed_arr
        endif
     endif

     if (print_maps) then
        if (main_prc .and. tot_ndir > 0) then
           print *, '--writing map_arr_out for total light'
           call print_map_arr_out(file_maps)
        endif
     endif

     if (print_maps_in) then
        if (main_prc .and. tot_ndir_in > 0) then
           print *, '--writing map_in_arr_out for total light'
           call print_map_in_arr_out(file_maps_in)
        endif
     endif

  end select
  
end subroutine make_output

!> This routine decides which arrays to read depending on the value of rt_type().
subroutine read_output
 
  integer :: i,j,k
  integer :: tot_el, ierr
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname
  
  select case (rt_type)

  case (rtt_start)

     if (file_restore) then

        if (main_prc) then  ! read by main_prc
           print *, 'restoring output part2...'
           print *, '-- restoring u_final'
        endif
        call read_ufield_arr  ! read by all processes 

        !tot_el = lnum*tot_ncell
        !call mpi_Bcast(u_final_arr,tot_el, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,ierr)  

        if (main_prc) print *, '-- restoring scaspe' ! read by all processes
        call read_scaspe_arr

        if (main_prc) print *, '-- restoring lum_lost' ! read by all processes
        call read_lum_lost
        
        if (tot_ndir > 0) then 
           if (main_prc) print *, '-- restoring i_obs'
           call read_i_obs_arr
        endif
        
        if (tot_ndir_in > 0) then 
           
           if (main_prc) print *, '-- restoring i_obs_in'
           call read_i_obs_in_arr
           
        endif

        if (print_psel_av) then ! no print, no read
           if (main_prc) then ! this is read only by the main process. In this way, the result is correct when the array is reduced in rt_mpi. 
              print *, '-- restoring psel_av_arr'
              call read_psel_av_arr
           endif
        endif
     endif

     if (file_restore_part1) then
        if (main_prc) then  
           print *, 'restoring output part 1...'
           print *, '-- restoring u_fest'
        endif
        call read_u_fest_arr ! read by all MPI processes 

       ! tot_el = lnum*tot_ncell
       ! call mpi_Bcast(u_fest_arr,tot_el, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,ierr)   
     endif

  case (rtt_read_i_obs_part2)

     if (tot_ndir > 0) then 
        if (main_prc) print *, 'reading i_obs_part2....'
        call read_i_obs_arr    
     endif

     if (tot_ndir_in > 0) then 
        if (main_prc) print *, 'reading i_obs_in_part2....'
        call read_i_obs_in_arr    
     endif
 
  case (rtt_read_i_obs)

     if (tot_ndir > 0) then 
        if (main_prc) print *, 'reading i_obs....'
        call read_i_obs_arr
     endif

     if (tot_ndir_in > 0) then 
        if (main_prc) print *, 'reading i_obs_in....'
        call read_i_obs_in_arr    
     endif

  case (rtt_read_ufield)

     if (main_prc) print *, 'reading u_final.... '
     call read_ufield_arr
     
  case(rtt_read_scaspe_tot)

     if (main_prc) print *, 'reading scaspe_tot....'
     call read_scaspe_tot_arr
     
  end select

  call print_done

end subroutine read_output

!> This routine reads previously printed arrays.
!> @param [in] filename Name of the file containing the arrays to be read
!> @param [in] dsetname Dsetname of the array to be read
!> @param [in] dims Array of scalar speficying the dimensions of the arrays to be read
subroutine read_big_array(filename,dsetname,dims)


  CHARACTER (LEN=lcar) :: filename ! File name
  CHARACTER (LEN=lcar) :: dsetname ! Dataset name

  INTEGER(HID_T) :: file_id                            ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
 
  INTEGER(HSIZE_T), DIMENSION(*)  :: dims  ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(2)  :: dims_read_scaspe, maxdims   
  INTEGER     ::   error,error_arr(7)  ! Error flag
 
  !
!    Initialize FORTRAN interface.
!
     CALL h5open_f (error)

     error_arr(1)=error
     
     CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)

     error_arr(2)=error

     CALL h5dopen_f(file_id, dsetname, dset_id, error)

     error_arr(3)=error
  
     select case(dsetname)
     case('u_final' )     
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , u_final, dims, error)
     case('i_obs') 
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , i_obs, dims, error)
     case('i_obs_in')
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , i_obs_in, dims, error) 
     case( 'scaspe' ) 
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , scaspe, dims, error)
     case('scaspe_tot')  
        ! this is used only in the i_obs algorithm. the following lines are needed to check that the dimensions derived from kp_sca_max and tot_ndir_scaspe are correct
        call h5dget_space_f(dset_id,dspace_id,error)
        call h5sget_simple_extent_dims_f(dspace_id, dims_read_scaspe, maxdims, error) 
        if (dims_read_scaspe(1) /= dims(1) .or. dims_read_scaspe(2) /= dims(2)) then 
           print *, 'error(read_big_array): the dimensions used to read the scaspe_tot files are not correct. If you are using the i_obs algorithm, are you sure you are using the same kp_sca_max used in the main RT calculation ?'
           print *, 'dims_read_scaspe =', dims_read_scaspe(:2)
           print *, 'dims =', dims(:2)
           stop
        endif
 
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , scaspe_tot, dims, error)
     case('u_fest') 
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , u_fest, dims, error) 
     case('psel_av_arr') 
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , psel_av_arr, dims, error)    
     case('lum_lost')
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , lum_lost, dims, error)
     case default
        print *, 'ERROR(read_big_array): dsetname not recognized!'
        print *, 'dsetname = ', dsetname
        STOP
     end select

     error_arr(4)=error

     CALL h5dclose_f(dset_id, error)
     error_arr(5)=error
    
     CALL h5fclose_f(file_id, error)
     error_arr(6)=error
     
     CALL h5close_f(error)
     error_arr(7)=error
     
     if (any(error_arr /= 0)) then
        print *, 'error while opening =', filename
        call stop_prc
     endif

  end subroutine read_big_array

!> Reads additional info from output files. 
subroutine read_info_big_array(filename,dsetname_add, dims_add, rank_add)

  CHARACTER (LEN=lcar) :: filename ! File name
  INTEGER(HSIZE_T) ::  dims_add(0:),dims_in(0:0)
  character (LEN=lcar)  :: dsetname_add(0:) 
  integer :: rank_add(0:)

  INTEGER(HID_T) :: file_id                            ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
     
  INTEGER     ::   error,error_arr(7)  ! Error flag
  integer :: i, narr 
  integer :: nside_sca_read

  narr = size(rank_add)
  !
  !    Initialize FORTRAN interface.
!
  CALL h5open_f (error)

  error_arr(1)=error
     
  CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)

  error_arr(2)=error

  do i=0, narr-1 

     CALL h5dopen_f(file_id, dsetname_add(i), dset_id, error)

     error_arr(3)=error
  
     if (dsetname_add(i) == 'nside_sca' ) then   
        dims_in = dims_add(i)
        CALL h5dread_f(dset_id,H5T_NATIVE_INTEGER, nside_sca_read, dims_in, error)   
     endif

     error_arr(4)=error
     
     CALL h5dclose_f(dset_id, error)
     error_arr(5)=error
     
  enddo
    
  CALL h5fclose_f(file_id, error)
  error_arr(6)=error
     
  CALL h5close_f(error)
  error_arr(7)=error
     
  if (any(error_arr /= 0)) then
     print *, 'ERROR(read_info_big_array): error while opening =', filename
     call stop_prc
  endif

  print *, 'nside read =', nside_sca_read
  stop

end subroutine read_info_big_array

! This routine writes the lum_lost() files at the end of the direct light processing and after the scattered light processing.   
!!$  subroutine write_lum_lost
!!$
!!$    character(LEN=lcar) :: filename,warning='WARNING'
!!$    real(kind=real64) :: frac_lost
!!$    integer :: unit,i
!!$
!!$    select case (rt_type)
!!$
!!$    case (rtt_output_part2,rtt_scatt)
!!$
!!$       do i = 0, lnum -1
!!$          
!!$          frac_lost=lum_lost(i)/tot_rad_en_or(i)
!!$
!!$          if (rt_type == rtt_output_part2) filename=file_llost_part2_arr(i)
!!$          if (rt_type == rtt_scatt) filename=file_llost_arr(i)
!!$
!!$          if (frac_lost > 0.01) then 
!!$             filename=trim(adjustl(dir_runs))//trim(adjustl(warning))//'_'//filename
!!$          else 
!!$             filename=trim(adjustl(dir_runs))//filename
!!$          endif
!!$
!!$          open(newunit=unit, file=filename, status='unknown')
!!$
!!$          write(unit,*) 'Total luminosity processed'
!!$          write(unit,*) tot_rad_en_or(i)
!!$          write(unit,*) 'Lost luminosity'
!!$          write(unit,*) lum_lost(i)
!!$          write(unit,*) 'Fraction of lost luminosity'
!!$          write(unit,*) frac_lost
!!$          
!!$          close(unit)
!!$
!!$       end do
!!$       
!!$    case DEFAULT
!!$
!!$       print *, 'here you should not get : write_lum_lost'
!!$
!!$    end select
!!$
!!$  end subroutine write_lum_lost

!> Writes lost luminosity files file_lum_lost_part2() and file_lum_lost()
subroutine print_lum_lost
  integer :: i,j,k
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname
  real(kind=real64) :: frac_lum(0:lnum-1)

  dsetname='lum_lost'
  rank=1
  allocate(dims(rank))
  dims=lnum
  
  if (rt_type == rtt_output_part2) then 
     filename=trim(adjustl(dir_runs))//file_lum_lost_part2
  elseif (rt_type == rtt_scatt) then
     filename=trim(adjustl(dir_runs))//file_lum_lost
  endif
  call print_big_array(filename,dsetname,dims,rank)
  
  deallocate(dims)

  ! check on the ratio lum_lost/tot_rad_en_or
  do i = 0, lnum-1 

     if (tot_rad_en_or(i) > 0) then 
        frac_lum(i) = lum_lost(i)/tot_rad_en_or(i) 
     endif

  end do
  
  if (any(abs(frac_lum) > 0.01)) then 
     print *, 'WARNING: Lost luminosity fraction higher than 1% in at least some wavelengths. Check lost luminosity output files!'
  endif

end subroutine print_lum_lost
  
!> Reads lum_lost() from previously written files.
!> \todo this routine is returned immediately because it can give problems in case the wavelength grid is changed compared to the previous calculation (not as for i_obs and scaspe that are written on different files giving more flexibility). Fix this problem and remove the return.  
subroutine read_lum_lost
  integer :: i
  integer :: ierr
  INTEGER     ::   rank 
  character (LEN=lcar)  :: filename,  dsetname, dir_tmp

  return  !
  
  dsetname='lum_lost'
  rank=1
  allocate(dims(rank))
  dims=lnum
  dir_tmp = ''

  if (rt_type == rtt_start) then 
     filename=trim(adjustl(dir_runs))//file_lum_lost_part2
  else 
     print*, 'ERROR(read_lum_lost): rt_type not recognized!'
     STOP
  endif
  call check_file_existence(dir_tmp, filename) 
  call read_big_array(filename,dsetname,dims)
    
  deallocate(dims)

end subroutine read_lum_lost



  !> This routine reads file_lambda_list(). It has to be placed before check_input because the lambda grid files contain a label dependent on the wavelength. It assigns lnum_tot(), lnum_stars(), lnum_dust(). 
  subroutine read_lambda_list
    integer :: i, v
    integer :: unit
    integer :: n_ind_i_obs, n_ind_out_maps
    integer, allocatable :: temp_ind_i_obs(:), temp_ind_out_maps(:)

    ! check dir_grid

    if (dir_grid == 'not_provided') then
       print *, 'ERROR: Input dir_grid missing. Did you assign dir_grid in the input files ?'
       call stop_prc
    else
       call check_dir_existence(dir_grid, .FALSE.)  
    endif

    
    ! check file_lambda_list
   if (file_lambda_list == 'not_provided' ) then
      print *, 'ERROR: Input file_lambda_list missing'
      call stop_prc
   else
      call check_file_existence(dir_grid, file_lambda_list)     
   endif
   
   open(newunit=unit,file=trim(adjustl(dir_grid))//trim(adjustl(file_lambda_list)), status='old')
   
   i=0 ! there are no comments in lambda list file
   v=0
   do 
      read(unit,*,iostat=v)
      if (v /= 0) exit
      i=i+1
   end do
   lnum_tot=i
   rewind(unit)

   allocate(lambda_arr(0:lnum_tot-1), label_wave_arr(0:lnum_tot-1))

   if (main_prc) print *, 'reading wavelength grid...'
   
   do i=0,lnum_tot-1
      
      read(unit,*) lambda_arr(i)
      write(label_wave_arr(i),'(F9.3)') lambda_arr(i) 
      if (main_prc) print *, 'lambda = ', trim(adjustl(label_wave_arr(i)))
      
      if (i > 0) then 
         if (lambda_arr(i) < lambda_arr(i-1)) then
            print *, 'Lambda list has to be in ascending order. Correct this in file_lambda_list before trying again... '
            close(unit)
            call stop_prc
         endif
      endif 
      
   enddo
   
   close(unit)

   if (use_lambda_grid) then

      if (label_model_lambda_grid == 'not_provided') then
         print *, 'ERROR: Input label_model_lambda_grid missing!'
         call stop_prc
      endif

      allocate(grid_file_lambda_arr(0:lnum_tot-1))
      
      do i = 0, lnum_tot-1
         
         grid_file_lambda_arr(i)='grid_'//trim(adjustl(label_model_lambda_grid))//'_l'//trim(adjustl(label_wave_arr(i)))//'um.h5'
         
      end do 
      
   endif

   lnum_stars = lnum_tot  ! initialise
   lnum_dust = lnum_tot
   i_lambda_stars = [0, lnum_tot-1]
   i_lambda_dust = [0, lnum_tot-1]
   
   ! find lnum_stars
   if (max_lambda_stars /= -1) then 

      do i = 1, lnum_tot-1

         if (lambda_arr(i) > max_lambda_stars) then
            lnum_stars = i
            i_lambda_stars(1) = i-1
            exit
         endif

      enddo

   else

      if (main_prc) print *, 'WARNING: max_lambda_stars input missing. Using entire wavelength grid. This might mean too many wavelengths will be used for the RT calculation of stellar emission!'
      max_lambda_stars = lambda_arr(lnum_tot-1)
      
   endif

   ! find lnum_dust
   if (min_lambda_dust == -1) then 

      if (main_prc) print *, 'WARNING: min_lambda_dust input missing. Using default value min_lambda_dust = 1 um!'
      select case(units_lambda)
      case('um')
         if (main_prc) print *, "WARNING: Input units_lambda is 'um' (microns)! "
         min_lambda_dust = 1.
      case ('not_provided')
         if (main_prc) print *, 'ERROR: Input units_lambda missing'!'
         call stop_prc         
      case default
         if (main_prc) print *, 'ERROR: units_lambda should be um (microns)!'
         call stop_prc
      end select
   endif

   do i = lnum_tot-1, 0, -1

      if (lambda_arr(i) <= min_lambda_dust) then
         lnum_dust = lnum_tot-i
         i_lambda_dust(0) = lnum_tot-lnum_dust
         exit
      endif
      
   end do
   
   ! set lnum to stellar emission wavelength range value
   lnum = lnum_stars

   ! return if in grid creation mode
   if (grid_creation) return 

   ! -----------------------------------------------------------------
   ! set ind_i_obs. the reallocation is needed because there is no way to know the input size of ind_i_iobs in advance. 
   allocate(temp_ind_i_obs(0:max_n_input-1)) 
   temp_ind_i_obs = ind_i_obs
   deallocate(ind_i_obs)
   n_ind_i_obs = count(temp_ind_i_obs /= -1)
   select case(n_ind_i_obs)
   case(0)
      if (main_prc) print *, 'WARNING: all i_obs files will be printed! Input wavelength indeces using ind_i_obs to select only few elements to be printed. E.g. ind_i_obs = 0,13,45'
      allocate(ind_i_obs(0:lnum_tot-1))
      do i = 0, lnum_tot-1
         ind_i_obs(i) = i
      enddo
      n_ind_i_obs = lnum_tot 
   case(1:)
      allocate(ind_i_obs(0:n_ind_i_obs-1))
      ind_i_obs = temp_ind_i_obs(0:n_ind_i_obs-1)
   case default
      print *, 'ERROR: you should not get here: real_lambda_list!'
      STOP
   end select
   deallocate(temp_ind_i_obs)

   ! print out ind_i_obs
   if (main_prc .and. n_ind_i_obs /= lnum_tot) then 
      print *, 'WARNING: i_obs arrays will be printed only for the following wavelengths. To change this, input a different ind_i_obs array.'
      do i = 0, n_ind_i_obs-1
         print *, 'lambda = ', lambda_arr(ind_i_obs(i))
      end do
   endif

   ! checks for ind_i_obs
   if (any(ind_i_obs > lnum_tot-1)) then 
      if (main_prc) print *, 'ERROR(read_lambda_list): at least one index in ind_i_obs is higher than the total number of wavelengths.'
      STOP
   endif

   ! -----------------------------------------------------------------
   ! set ind_out_maps. the reallocation is needed because there is no way to know the input size of ind_out_maps in advance. 
   allocate(temp_ind_out_maps(0:max_n_input-1)) 
   temp_ind_out_maps = ind_out_maps
   deallocate(ind_out_maps)
   n_ind_out_maps = count(temp_ind_out_maps /= -1)
   select case(n_ind_out_maps)
   case(0)
      if (main_prc) print *, 'WARNING: all surface brightness maps will be printed! Input wavelength indeces using ind_out_maps to select only few elements to be printed. E.g. ind_out_maps = 0,13,45'
      allocate(ind_out_maps(0:lnum_tot-1))
      do i = 0, lnum_tot-1
         ind_out_maps(i) = i
      enddo
      n_ind_out_maps = lnum_tot
   case(1:)
      allocate(ind_out_maps(0:n_ind_out_maps-1))
      ind_out_maps = temp_ind_out_maps(0:n_ind_out_maps-1)
   case default
      print *, 'ERROR: you should not get here: real_lambda_list!'
      STOP
   end select
   deallocate(temp_ind_out_maps)  

   ! print out ind_out_maps
   if (main_prc .and. n_ind_out_maps /= lnum_tot) then 
      print *, 'WARNING: surface brightness maps will be printed only for the following wavelengths. To change this, input a different ind_out_maps array.'
      do i = 0, n_ind_out_maps-1
         print *, 'lambda = ', lambda_arr(ind_out_maps(i))
      end do
   endif

   ! checks for ind_out_maps
   if (any(ind_out_maps > lnum_tot-1)) then 
      if (main_prc) print *, 'ERROR(read_lambda_list): at least one index in ind_out_maps is higher than the total number of wavelengths.'
      STOP
   endif

   ! set lnum_maps 
   lnum_maps = n_ind_out_maps 

   ! check lnum_maps in comparison to lnum_stars and lnum_dust 
   if ((print_maps .or. print_maps_in).and.(lnum_maps > 10*lnum_stars .or. lnum_maps > 10*lnum_dust)) then 
      if (main_prc) print *, 'WARNING: Maps for many wavelengths will be printed with respect to the number of wavelengths used in the RT calculations for the stellar emission RT and/or the dust emission RT. You might add appropriate n_ind_out_maps to the input file.'
   endif


   ! set lambda_arr_maps
   allocate(lambda_arr_maps(0:lnum_maps-1))
   do i = 0, lnum_maps-1  
         lambda_arr_maps(i) = lambda_arr(ind_out_maps(i))
   end do

   ! ----------------------
   ! check lnum_dust
   if (lnum_dust < 10 .or. lnum_stars < 10) then
      if (main_prc) print *, 'WARNING: less than 10 wavelengths in the stellar emission and/or dust emission wavelength range. No dust RT will be performed!'
      no_dust_rt = .TRUE.
   endif
   
   call print_done
   
 end subroutine read_lambda_list

!> This routine prints the lambda grids in grid_file_lambda(). 
subroutine print_lambda_grid
    
     INTEGER(HID_T) :: file_id                            ! File identifier
     INTEGER(HID_T) :: dset_id       ! Dataset identifier
     INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
     integer,parameter  :: narr_grid_lambda=2
     INTEGER     ::   rank(0:narr_grid_lambda-1)                         ! Dataset rank
     INTEGER     ::   error  ! Error flag
     character(LEN=lcar) :: dsetname(0:narr_grid_lambda-1)
     integer :: i

     dsetname(0)='dens'   ; rank(0)=1
     dsetname(1)='dens_stars'  ; rank(1)=1
     

!    Initialize FORTRAN interface.
     CALL h5open_f (error)

 ! Create a new file using default properties.
     ! 
     CALL h5fcreate_f(trim(adjustl(dir_grid))//trim(adjustl(grid_file_lambda)), H5F_ACC_TRUNC_F, file_id, error)    


     do i=0,narr_grid_lambda-1
        
        call make_dims(dsetname(i),rank(i))
           
     ! Create the dataspace.     
     CALL h5screate_simple_f(rank(i), dims, dspace_id, error)
     
! Create the dataset with default properties.
        CALL h5dcreate_f(file_id, dsetname(i),H5T_NATIVE_DOUBLE , dspace_id, &
          dset_id, error)
    
 
! Open an existing dataset.
     !
    ! CALL h5dopen_f(file_id, dsetname(i), dset_id, error)  
   ! note the input arrays might actually be bigger than dims when created with the standard create_adap_program
  
  if (dsetname(i) == 'dens' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , dens, dims, error)
  else if (dsetname(i) == 'dens_stars' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , dens_stars, dims, error) 
  endif

! End access to the dataset and release resources used by it.
     !
     CALL h5dclose_f(dset_id, error)

     !
     ! Terminate access to the data space.
     !
     CALL h5sclose_f(dspace_id, error)

  
       deallocate(dims) 
    end do
        
  ! Terminate access to the file.
     !
     CALL h5fclose_f(file_id, error)     

!    Close FORTRAN interface.
!
     CALL h5close_f(error)

     end subroutine print_lambda_grid

     !> This routines finds out the value of tot_ndir_scaspe needed in the i_obs algorithm (in this case this value differ from tot_ndir. This is so because in the i_obs algorithm you might want to use a different number of output directions compared to the one you used in the main RT run.).
subroutine find_out_tot_ndir_scaspe

  CHARACTER (LEN=lcar) :: filename ! File name
  CHARACTER (LEN=lcar) :: dsetname ! Dataset name

  INTEGER(HID_T) :: file_id                            ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
  INTEGER(HSIZE_T) :: dims_in(0:0)
  INTEGER(HSIZE_T), DIMENSION(2)  :: dims,maxdims  ! Dataset dimensions
  INTEGER     ::   error,error_arr(7)  ! Error flag
  integer :: i 
  integer :: tot_ndir_prev    !!! used to check all tot_ndir_scaspe are the same
  integer :: nside_sca_read
  
  if (main_prc) print *, 'figuring out the tot_ndir_scaspe value from previously written files....'

  if (only_direct_rt) then 
     tot_ndir_scaspe = 0 ! scaspe array not important if only direct light is processed.
     return
  endif
  
  do i = 0, lnum -1

     call check_file_existence(dir_runs, file_scaspe_tot_arr(i))  
     
  filename=trim(adjustl(dir_runs))//file_scaspe_tot_arr(i)
  dsetname='scaspe_tot'
  
  CALL h5open_f (error)

  CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)

  CALL h5dopen_f(file_id, dsetname, dset_id, error)

  call h5dget_space_f(dset_id,dspace_id,error)

 ! Getting dims from dataspace
  call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)  

  CALL h5dclose_f(dset_id, error)

  dsetname='nside_sca'
  dims_in = 1

  CALL h5dopen_f(file_id, dsetname, dset_id, error)

  CALL h5dread_f(dset_id,H5T_NATIVE_integer , nside_sca_read, dims_in, error)

  CALL h5dclose_f(dset_id, error)

!!! ACHTUNG : the following is not needed anymore because nside can be varied with wavelength now. 
!!$  if (nside_sca_read /= nside_sca) then 
!!$     print *, 'STOP: the parameter nside_sca read from the scaspe_tot arrays does not coincide with that derived from the input kp_sca. Set the correct kp_sca and try again'
!!$     print *, 'nside_sca_read =', nside_sca_read
!!$     print *, 'nside_sca =', nside_sca
!!$     stop
!!$  endif
  
  if (nside_sca_read /= 0 ) then 

     tot_ndir_scaspe=dims(1)-12*nside_sca_read**2

  else ! this is the case of assumed isotropic scattering (just 1 element for the scattering source function) 

    tot_ndir_scaspe = dims(1) -1 

 endif
  
  CALL h5fclose_f(file_id, error)
  
  CALL h5close_f(error)

  if (i > 0 .and. tot_ndir_scaspe /= tot_ndir_prev) then

     print *, 'ERROR: the values of tot_ndir_scaspe are not the same in the different scaspe_tot files. Have you generated them in the same calculation ?'
     CALL STOP_PRC
  endif

  tot_ndir_prev = tot_ndir_scaspe 
  
end do

call print_done 
  
end subroutine find_out_tot_ndir_scaspe

!> This routine initialize all the input variables, so it is possible to check whether all required input has been provided in check_input().
subroutine input_initialize

  !label_model = 'not_provided'
  label_model_lambda_grid = 'not_provided'
  label_model_out = 'not_provided'   
  label_model_out_i_obs = 'not_provided'   
  file_dir_out = 'not_provided'  
  file_pos_obs = 'not_provided' 
  file_p_src = 'not_provided'
  file_lambda_list = 'not_provided' 
  dir_runs = 'not_provided'
  dir_grid = 'not_provided'
  grid_file = 'not_provided'
  grid_file_lambda = 'not_provided'
  units_luminosity = 'not_provided'
  units_csize = 'not_provided'
  units_lambda = 'not_provided'
  dust_model = 'not_provided'
  dust_opacity_tables = 'not_provided'
  file_gra_fa = 'not_provided'
  file_sil_fa = 'not_provided'
  file_pah_neu_fa = 'not_provided'
  file_pah_ion_fa = 'not_provided'
  file_q_gra = 'not_provided'
  file_q_sil = 'not_provided'
  file_q_pah_neu = 'not_provided'
  file_q_pah_ion = 'not_provided'
  file_av_opacities = 'not_provided'
  dust_heating_type = 'not_provided' 
  file_calorimetry_Gra = 'not_provided' 
  file_calorimetry_Sil =  'not_provided'
  file_stellar_library = 'not_provided'
  stellar_library = 'not_provided'
  param_to_project = 'not_provided'
  file_param_src = 'not_provided'
  
  kp_sca_max = -1          
  rad_lim = -1.0   
  accuracy = -1.0   
  conv_en_lim = -1.0 
  bm_par = -1 
  bm_par_sca = -1   
  bm_par_max = -1  
  tau_cell_max = -1    
  max_lambda_stars = -1.
  min_lambda_dust = -1.
  dist_obs = -1.
  lambda_ref = -1.
  n_dust_size_qabs = 0
  n_dust_wave_qabs = 0 
  n_dust_temp_cal = 0 
  npixel_maps = -1
  map_size_factor = -1
  kp_maps = -1
  x_wall_coord = (/0.,1./)
  y_wall_coord = (/0.,1./)
  z_wall_coord = (/0.,1./)
  z_sun = 0.018
  max_sca_iterations = -1 
  n_int_rf_bins = -1 
  print_scaspe_tot = .FALSE.
  print_output_part1 = .FALSE.
  print_output_part2 = .FALSE.
  print_scaspe_part2 = .FALSE.
  print_psel_av = .FALSE.
  print_sed = .FALSE.
  restore_file_mpi = .FALSE.
  use_lambda_grid = .FALSE.
  use_dir_out = .FALSE.
  use_pos_obs = .FALSE.
  use_p_src = .FALSE.
  sequential_scattering = .FALSE.
  input_av_opacities = .FALSE. 
  rt_algorithm = 'not_provided'      
  no_communications = .FALSE.
  no_dust_rt = .FALSE.
  only_direct_rt = .FALSE.
  grid_creation = .FALSE. 
  grid_creation_lambda = .FALSE. 
  test_run = .FALSE. 
  print_maps = .FALSE.
  print_maps_in = .FALSE.
  x_wall_on = .FALSE.
  y_wall_on = .FALSE.
  z_wall_on = .FALSE.
  use_stellar_library = .FALSE. 
  limit_scattering_iterations = .FALSE.
  
end subroutine input_initialize

!> This routine checks that the required input has been provided and that the values are in the allowed range. It also checks that the file and directories actually exist.
subroutine check_input
  !> @param nproc_max Number of CPUs on the machine where DART-Ray is running. Found through linux command getconf _NPROCESSORS_ONLN.  
  integer :: nproc_max
  !> @param inum Counter used to check directory existence.
  integer :: inum
  !> @param file_exists Logical variable used to check file existence. 
  logical :: file_exists
  !> @param file_tmp Temporary file to store output linux commands
  character(LEN=lcar) :: file_tmp
  !> @param i Do loop counter
  integer :: i
  !> @param ierr MPI error argument
  integer :: ierr
  character (LEN= lcar) :: dir_tmp
  integer(kind = int64) :: count, count_rate, count_max
  integer :: unit
  logical :: error

  error = .FALSE.
 
  call system_clock(count, count_rate, count_max)

  write(file_tmp, *) count+id_mpi*1000  ! adding id_mpi makes the difference  
  file_tmp = 'dartray.tmp.'//trim(adjustl(file_tmp))

  if (main_prc) print *, 'checking input string variables...'

  ! rt_algorithm
  select case(rt_algorithm)
  case('main')
     rt_algorithm_ID = rta_main
  case('2D')
     rt_algorithm_ID = rta_2D
  case('i_obs')
     rt_algorithm_ID = rta_i_obs
  case('i_obs_dust')
     rt_algorithm_ID = rta_i_obs_dust
  case('sed')
     rt_algorithm_ID = rta_sed
  case('sed_dust')
     rt_algorithm_ID = rta_sed_dust   
  case('dust')
     rt_algorithm_ID = rta_dust
  case('dust2D')
     rt_algorithm_ID = rta_dust2d
  case('projection')
     rt_algorithm_ID = rta_projection
  case default
     if (main_prc) then
        print *, 'ERROR: Input rt_algorithm missing or not valid'
        print *, 'rt_algorithm = ', rt_algorithm
     endif
     error = .TRUE.
  end select

  ! dust heating type
  select case(dust_heating_type)
  case('eff')
     dust_heating_type_ID = dt_eff
  case('equ')
     dust_heating_type_ID = dt_equ
  case('sto')
     dust_heating_type_ID = dt_sto
  case('sto_lib')
     dust_heating_type_ID = dt_sto_lib
  case ('not_provided')
     dust_heating_type_ID = dt_none ! in case of dust RT algorithm, the calculation will stop below. 
  case default
     if (main_prc) then
        print *, 'ERROR: Input dust_heating_type not recognized!'
        print *, 'dust_heating_type = ', dust_heating_type
     endif
     error = .TRUE.
  end select

  ! dust_heating type needed for dust RT calculations
  select case(rt_algorithm_ID)
  case(rta_dust, rta_dust2D, rta_i_obs_dust)
     if (dust_heating_type_ID == dt_none) then
        if (main_prc) print*, 'ERROR: specify dust_heating_type to run RT algorithm ', trim(rt_algorithm)
        error = .TRUE.
     endif

     if (no_dust_rt) then 
        if (main_prc) print*, 'WARNING: no_dust_rt = .FALSE. needed for dust RT algorithms. Set to FALSE now.'
        no_dust_rt = .FALSE.
     endif

  end select

  ! assign print_sed to true if rt_algorithm = 'sed'
  if (.not.(print_sed) .and. (rt_algorithm_ID == rta_sed .or. rt_algorithm_ID == rta_sed_dust)) then
     if (main_prc) print *, 'WARNING: SED algorithm requires print_sed = .TRUE.. This is set now.'
     print_sed = .TRUE.
  endif

  ! assign use_dir_out to false if rt_algorithm == '2D'
  if (use_dir_out .and. (rt_algorithm_ID == rta_2D .or. rt_algorithm_ID == rta_dust2D)) then
     if (main_prc) print *, 'WARNING: use_dir_out not allowed for rt_algorithm = 2D. Set to FALSE.'
     use_dir_out = .FALSE.
     print_sed = .FALSE.
  endif

  ! assign use_pos_obs to false if rt_algorithm == '2D'
  if (use_pos_obs .and. (rt_algorithm_ID == rta_2D .or. rt_algorithm_ID == rta_dust2D)) then
     if (main_prc) print *, 'WARNING: use_pos_obs not allowed for rt_algorithm = 2D. Set to FALSE.'
     use_pos_obs = .FALSE.
  endif


  ! set no_communications to false if rt_algorithm_ID = rta_i_obs, rta_i_obs_dust and rta_sed
  if (no_communications) then
     select case(rt_algorithm_ID)
     case(rta_i_obs, rta_i_obs_dust, rta_sed, rta_projection)
        if (main_prc) print *, 'WARNING: no_communications has to be set to .FALSE. for i_obs, sed and projection algorithms. Done now.'
        no_communications = .FALSE.
     end select
  endif

  ! in case of the projection algorithm, the quantity to project should also be assigned.  
  if (rt_algorithm_ID == rta_projection) then 
     select case(param_to_project)
     case('stellar_emission', 'optical_depth')

     case('not_provided')
        if (main_prc) print *, 'ERROR: param_to_project has to be input in projection algorithm'
        error = .TRUE.
     case default
        if (main_prc) print*, 'ERROR: param_to_project not recognized!'
        error = .TRUE.
     end select
  endif

  
  ! label_model  -- This does not need to be defined necessarily   
  !if (label_model == 'not_provided') then
  !   print *, 'STOP: Input label_model missing'
  !   stop
  !endif

  ! label_model_lambda_grid -- This does not need to be defined necessarily
  !if (label_model_lambda_grid == 'not_provided' .and. use_lambda_grid) then
  !   print *, 'STOP: Input label_model_lambda_grid missing'
  !   stop
  !endif

  ! label_model_out
  if (label_model_out == 'not_provided') then
     print *, 'ERROR: Input label_model_out missing'
     error = .TRUE.
  endif

  ! label_model_out_i_obs 
  if (rt_algorithm_ID == rta_i_obs .or. rt_algorithm_ID == rta_i_obs_dust) then 
     if (label_model_out_i_obs == 'not_provided') then 
        label_model_out_i_obs = label_model_out
     endif
  endif

  ! file_dir_out
  if (use_dir_out) then 
     if (file_dir_out == 'not_provided' ) then
        print *, 'ERROR: Input file_dir_out missing. If output i_obs not required, set use_dir_out = .FALSE.'
        ERROR = .TRUE.
     else
        call check_file_existence(dir_grid, file_dir_out)     
     endif
  endif

  ! file_pos_obs
  if (use_pos_obs) then 
     if (file_pos_obs == 'not_provided') then
        print *, 'ERROR: Input file_pos_obs missing. If output i_obs_in not required, set use_pos_obs = .FALSE.'
        ERROR = .TRUE.
     else
        call check_file_existence(dir_grid, file_pos_obs)     
     endif
  endif
     
  ! file_p_src
  if (use_p_src) then 
     if (file_p_src == 'not_provided') then
        print *, 'ERROR: Input file_p_src missing. If there are no point sources to be considered, then set use_p_src = .FALSE.'
        ERROR = .TRUE.
     else
        call check_file_existence(dir_grid, file_p_src)     
     endif
  endif

  ! file_lambda_list
  !! this is checked in read_lambda_list which is called before check_input. It is needed to assign names to lambda grids.  

  ! dir_runs
  if (dir_runs == 'not_provided') then
     print *, 'ERROR: Input dir_runs missing'
     error = .TRUE.
  else     
     call check_dir_existence(dir_runs, .TRUE.)     
  endif

  ! dir_grid
  ! this is checked in read_lambda_list. 

  ! grid_file
  if (grid_file == 'not_provided') then
     print *, 'ERROR: Input grid_file missing'
     error = .TRUE.
  else
     call check_file_existence(dir_grid, grid_file) 
  endif

  ! check lambda_ref
  if (lambda_ref <= 0) then
        print *, 'ERROR: Invalid lambda_ref value'
        print *, 'lambda_ref = ', lambda_ref 
        print *, 'Allowed range = (0,Inf]'
        error = .TRUE.
     endif

  
  ! grid_file_lambda_arr. It has to be checked that the lambda grid corresponding to the reference wavelength is present. This is done in read_lambda_grids
  if (use_lambda_grid) then     

     do i = 0, lnum -1 
        grid_file_lambda = grid_file_lambda_arr(i)
        call check_file_existence(dir_grid, grid_file_lambda)
     enddo
     
  endif

  ! units variables

  select case(units_luminosity)
  case('erg/s/Hz', 'W/Hz')
     if (main_prc) print *, 'Input luminosity units are ', trim(adjustl(units_luminosity)) 
  case('not_provided')
     if (main_prc) print *, 'ERROR: input units_luminosity missing!'
     error = .TRUE.
  case default
     if (main_prc) then 
        print *, 'ERROR: units_luminosity not recognized! Possible choices: erg/s/Hz, W/Hz.'
        print *, 'units_luminosity = ', trim(units_luminosity)
     endif
     error = .TRUE.
  end select
  
  select case(units_csize)
  case('pc')
     if (main_prc) print *, 'Input grid cell size units are', trim(adjustl(units_csize))
  case('not_provided')
     if (main_prc) print *, 'ERROR: input units_csize missing!'
     error = .TRUE.
  case default
     if (main_prc) then 
        print *, 'ERROR: units_csize not recognized! Possible choices are: pc.'
        print *, 'units_csize = ', trim(units_csize)
     endif
     error = .TRUE.        
  end select

  call error_stop(error)
  
  call check_input_dust_model 

  call print_done 
  
  if (main_prc) print *, 'checking input parameters...'

  error = .FALSE.

  ! kp_sca 
  if (kp_sca_max < 0 .or. kp_sca_max > 4) then
     print *, 'ERROR: Invalid kp_sca_max value'
     print *, 'kp_sca_max = ', kp_sca_max 
     print *, 'Allowed range = [0,4]'
     error = .TRUE.
  endif

  ! rad_lim
  if (rad_lim < 0 .or. rad_lim > 2) then
     print *, 'ERROR: Invalid rad_lim value'
     print *, 'rad_lim =', rad_lim
     print *, 'Allowed range = [0, 2]'
     error = .TRUE.
  endif

  ! accuracy
  if (accuracy < 0 .or. accuracy > 1) then
     print *, 'ERROR: Invalid accuracy value'
     print *, 'accuracy = ', accuracy
     print *, 'Allowed range = [0,1]'
     error = .TRUE.
  endif

  ! conv_en_lim
  if (conv_en_lim < 0. .or. conv_en_lim > 1) then
     print *, 'ERROR: Invalid conv_en_lim value'
     print *, 'conv_en_lim = ', conv_en_lim
     print *, 'Allowed range = [0,1]'
     error = .TRUE.
  endif

  ! bm_par
  if (bm_par < 0 .or. bm_par > 1000) then
     print *, 'ERROR: Invalid bm_par value'
     print *, 'bm_par = ', bm_par
     print *, 'Allowed range = [0,1000]'
     error = .TRUE.
  endif

  ! bm_par_sca
  if (bm_par_sca < 0 .or. bm_par_sca > 1000) then
     print *, 'ERROR: Invalid bm_par_sca value'
     print *, 'bm_par_sca = ', bm_par_sca
     print *, 'Allowed range = [0,1000]'
     error = .TRUE.
  endif
  
  ! bm_par_max
  if (bm_par_max < 10*maxval([bm_par, bm_par_sca])) then
     print *, 'ERROR: Invalid bm_par_max value'
     print *, 'bm_par_max = ', bm_par_max
     print *, 'Allowed range = [10*max([bm_par,bm_par_sca]), inf]'
     error = .TRUE.
  endif

  ! max_sca_iterations
  if (limit_scattering_iterations ) then 
     if (max_sca_iterations < 0 .or. max_sca_iterations > 10**6) then
        print *, 'ERROR: Invalid max_sca_iterations value'
        print *, 'max_sca_iterations = ', max_sca_iterations
        print *, 'Allowed range = [0,1E6]'
        error = .TRUE.
     endif
     if (main_prc) print *, 'WARNING: limit_scattering_iteration = .TRUE. The maximum number of scattering iterations will be :', max_sca_iterations
  endif
  
  
  ! tau_cell_max
  if (tau_cell_max < 0 .or. tau_cell_max > 0.5) then
     !print *, 'STOP: Invalid tau_cell_max value'
     !print *, 'tau_cell_max = ', tau_cell_max
     !print *, 'Allowed range = [0,0.5]'
     !error = .TRUE.
     if (main_prc) print *, 'WARNING: tau_cell_max not specified or outside allowed range of values (0,0.5). Set default value tau_cell_max = 0.05.'
     tau_cell_max = 0.05 
  endif

  ! nproc
  ! Getting number of CPUs on machine. Using "system" here because execute_command_line is implemented only in latest compiler versions. 
  call system("getconf _NPROCESSORS_ONLN > "//trim(adjustl(file_tmp)))
  open(newunit=unit, file = file_tmp, status ='old')
  read(unit,*) nproc_max
  close(unit)
  call system("rm "//trim(adjustl(file_tmp)))
  if (nproc < 1 .or. nproc > nproc_max) then
     print *, 'ERROR: Invalid nproc value'
     print *, 'nproc = ', nproc
     print *, 'Allowed range = [1,nproc_max] with nproc_max = ', nproc_max
     error = .TRUE.
  endif

  ! dist_obs
  if (print_sed) then
     if (dist_obs <= 0.1) then
        print *, 'ERROR: Invalid dist_obs'
        print *, 'dist_obs = ', dist_obs
        print *, 'Allowed range = (0.1, Inf]'
        print *, 'Make sure this parameter is defined in the input file.'
        error = .TRUE.
     endif
  endif

  ! check number of integrated energy bins for the SED library approach of the stochastically heated dust emission 
  if (dust_heating_type_ID == dt_sto_lib) then 
     if (n_int_rf_bins <= 10) then 
        print *, 'ERROR: Invalid n_int_rf_bins'
        print *, 'n_int_rf_bins = ', n_int_rf_bins
        print *, 'Allowed range = (10, 1000]'
        print *, 'Make sure this parameter is defined in the input file.'
        error = .TRUE.
     endif   
  endif

  call error_stop(error)
  
  call print_done 
  
  if (main_prc) print *, 'check input logical variables...'

  error = .FALSE.

  if (sequential_scattering) then 
     if (main_prc) print *, 'Using sequential scattering. One scattering order at the time!'
  endif
  
  if (print_output_part1) then
     if (main_prc) print *, 'Output part 1 will be written on disk!'
  else
     if (main_prc) print *, 'Output part 1  will NOT be written on disk! To change this, set print_output_part1 = .TRUE. '     
  endif

  if (print_output_part2) then
     if (main_prc) print *, 'Output part 2 will be written on disk! Note that by default the scaspe array is not written. To change that you need to set print_scaspe_part2 = .TRUE.'
  else
     if (main_prc) print *, 'Output part 2  will NOT be written on disk! To change this, set print_output_part2 = .TRUE. '     
  endif

  if (print_scaspe_part2) then
     if (main_prc) print *, 'WARNING: Scaspe array will be output on disk at the end of the direct light processing. This might take A LOT of disk space.'
  endif
  
  if (print_scaspe_tot) then
     if (main_prc) print *, 'WARNING: Scaspe arrays will be output on disk at the end of the scattering iterations. This might take A LOT of disk space.'
  endif

  if (np_mpi > 1) then 
  
     if (.not. restore_file_mpi .and. main_prc) print *, 'WARNING: Previously written intermediate wil not be restored! If you want that, set restore_file_mpi = .TRUE.' 

  endif
     
  if (.not. use_lambda_grid) then
     if (main_prc) print *, 'WARNING: This RT calculation is not using Lambda grids. Set use_lambda_grid = .TRUE. if they have to be used.'
  endif

  if (.not. use_dir_out) then
     if (main_prc) print *, 'WARNING: External observer directions not included. Set use_dir_out = .TRUE. if file_dir_out has to be read.'
  endif

  if (.not. use_pos_obs) then
     if (main_prc) print *, 'WARNING: Internal observer positions not included. Set use_pos_obs = .TRUE. if file_pos_obs has to be read.'
  endif

  if (.not. use_p_src) then
     if (main_prc) print *, 'WARNING: Point sources not included. Set use_p_src = .TRUE. if file_p_src has to be read.'
  endif

  if (.not. print_sed) then
     if (main_prc) print *, 'WARNING: Total emission SED will not be printed. Set print_sed = .TRUE. if this has to be done.'
  endif

  if (no_communications) then
     if (main_prc) print *, 'WARNING: no_communications = .TRUE. Scaspe and i_obs arrays will NOT be distributed (this might use lots of memory). This mode requires sequential_scattering = .TRUE. This is set automatically if needed.'
     if (.not. sequential_scattering) sequential_scattering = .TRUE.
  endif

  if (no_dust_rt) then
     if (main_prc) print *, 'WARNING: no_dust_rt = .TRUE. No dust emission RT calculation will be performed!'
  endif
  
  if (only_direct_rt) then 
     if (main_prc) print *, 'WARNING: only_direct_rt = .TRUE. Only direct light will be processed (no scattered light iterations)'
     print_output_part2 = .TRUE. ! by default output part2 is printed if only direct light is processed. 
  endif

  if (test_run) then 
     if (main_prc) print *, 'WARNING: test_run = .TRUE. RT subroutines will not be run! This mode is only to check that the program can run without problems. DO NOT USE THE OUTPUT FILES! THE RESULTS ARE MEANINGLESS! Set test_run = .FALSE. for running the normal RT calculations.'
  endif

  if (print_maps) then
     if (main_prc) print *, 'print_maps = .TRUE. Surface brightness maps for the external observers will be calculated and output!'
     if (npixel_maps < 200) then 
        if (main_prc) then 
           print *, 'ERROR: npixel_maps value not allowed!'
           print *, 'npixel_maps =', npixel_maps
           print *, 'allowed range = [200, Inf]'
        endif
        error = .TRUE.
     endif
     
     if (map_size_factor < 0.2) then
        if (main_prc) then 
           print *, 'ERROR: map_size_factor value not allowed!'
           print *, 'map_size_factor =', map_size_factor
           print *, 'allowed range = [0.2, Inf]'
        endif
        error = .TRUE.
     endif
     
  endif

  if (print_maps_in) then
     
     if (main_prc) print *, 'print_maps = .TRUE. Surface brightness maps for the internal observers will be calculated and output!'
     if (kp_maps < 1) then 
        if (main_prc) then 
           print *, 'ERROR: kp_maps value not allowed!'
           print *, 'kp_maps =', kp_maps
           print *, 'allowed range = [1, Inf]'
        endif
        error = .TRUE.
     endif

  endif

  
  if (x_wall_on) then
     if (main_prc) print *, 'x_wall_on = .TRUE. X-perpendicular wall will be applied'
     if (x_wall_coord(1) > x_wall_coord(2) .or. (x_wall_coord(1) < 0) .or. (x_wall_coord(2) > 1)) then
        if (main_prc) then
           print *, 'ERROR: X_wall_coord values not allowed!'
           print *, 'x_wall_coord = ', x_wall_coord
           print *, 'allowed range = [0,1]'
        endif
        error = .TRUE.
     endif
  endif

  if (y_wall_on) then
     if (main_prc) print *, 'y_wall_on = .TRUE. X-perpendicular wall will be applied'
     if (y_wall_coord(1) > y_wall_coord(2) .or. (y_wall_coord(1) < 0) .or. (y_wall_coord(2) > 1)) then
        if (main_prc) then
           print *, 'ERROR: y_wall_coord values not allowed!'
           print *, 'y_wall_coord = ', y_wall_coord
           print *, 'allowed range = [0,1]'
        endif
        error = .TRUE.
     endif
  endif

  if (z_wall_on) then
     if (main_prc) print *, 'z_wall_on = .TRUE. X-perpendicular wall will be applied'
     if (z_wall_coord(1) > z_wall_coord(2) .or. (z_wall_coord(1) < 0) .or. (z_wall_coord(2) > 1)) then
        if (main_prc) then
           print *, 'ERROR: z_wall_coord values not allowed!'
           print *, 'z_wall_coord = ', z_wall_coord
           print *, 'allowed range = [0,1]'
        endif
        error = .TRUE.
     endif
  endif

  if (use_stellar_library) then 
     select case(stellar_library)
     case('starburst99', 'maraston05_kr_rhb')
     
     case('user')
        if (main_prc) then 
           print *, 'User defined stellar library not yet implemented!' 
        endif
        error = .TRUE.
      
     case default
        if (main_prc) then 
           print *, 'ERROR(check_input): stellar_library not found!'
           print *, 'stellar library =', stellar_library
        endif
        error = .TRUE.
        
     end select
     
  endif

  call error_stop(error)

  call print_done 
  
 end subroutine check_input

!> Checks that the input regarding the input dust model is correct. 
subroutine check_input_dust_model
  integer :: i
  !> @param file_q_arr Set of input Q opacity files
  character(LEN=lcar) :: file_q_arr(0:max_n_dust_comp-1)
  character(LEN=lcar) :: file_calorimetry_arr(0:max_n_dust_cal_type-1)
  character(LEN= lcar) :: dir_tmp
  logical :: c1, c2 

! dust model and input size distribution files
  if (dust_model == 'not_provided') then
     if (main_prc) print *, 'ERROR: specify input dust_model!'
     call stop_prc
  else if (dust_model == 'user') then
     if (main_prc) print *, 'User provided grain size distributions will be used.'
     iq_dust_model = [0,0,0,0]
     dir_tmp = './'  ! here because check_file_existence accept only LEN=lcar characters. 
     if (file_gra_fa /= 'not_provided') then
        call check_file_existence(dir_tmp, file_gra_fa)
        iq_dust_model(0) = 1
        if (main_prc) print *, 'Using Graphite grain size distribution!'
     else
        if (main_prc) print *, 'WARNING: Graphite grain size distribution file not provided (file_gra_fa)!'
     endif
     if (file_sil_fa /= 'not_provided') then
        call check_file_existence(dir_tmp, file_sil_fa)
        iq_dust_model(1) = 1
        if (main_prc) print *, 'Using Silicate grain size distribution!'
     else
        if (main_prc) print *, 'WARNING: Silicate grain size distribution file not provided (file_sil_fa)!'
     endif
     if (file_pah_neu_fa /= 'not_provided') then
        call check_file_existence(dir_tmp, file_pah_neu_fa)
        iq_dust_model(2) =1
        if (main_prc) print *, 'Using Neutral PAH molecule size distribution!'
     else
        if (main_prc) print *, 'WARNING: Neutral grain size distribution file not provided (file_pah_neu_fa)!'
     endif
     if (file_pah_ion_fa /= 'not_provided') then
        call check_file_existence(dir_tmp, file_pah_ion_fa)
        iq_dust_model(3) =1
        if (main_prc) print *, 'Using ionized PAH molecule size distribution!'
     else
        if (main_prc) print *, 'WARNING: Ionized PAH molecule size distribution file not provided (file_pah_ion_fa)!'
     endif

     if (.not. any(iq_dust_model)) then
        if (main_prc) print *, "ERROR: To use dust_model = 'user', specify at least one grain size distribution file (file_gra_fa, file_sil_fa, file_pah_neu_fa, file_pah_ion_fa)."
        STOP
     endif

     select case(dust_opacity_tables) 
     case('TRUST')
        if (main_prc) print *, 'Dust opacity Q coefficients from the opacity tables of the TRUST benchmark project (INSERT REFERENCE HERE)!'
     case('DraineLi06')
        if (main_prc) print *, 'Dust opacity Q coefficients from the opacity tables of Draine & Li 2006 (CHECK THIS REFERENCE)!'
     case('user')
        ! check input opacity tables and n_dust_size_qabs
        if (main_prc) print *, 'Checking existence Q opacity files... ' 
        file_q_arr = (/file_q_gra, file_q_sil, file_q_pah_neu, file_q_pah_ion/)
        do i = 0, max_n_dust_comp-1 
           if (iq_dust_model(i)) then
              call check_file_existence(dir_tmp, file_q_arr(i))
              if (n_dust_size_qabs(i) <= 10) then
                 if (main_prc) then 
                    print *, 'ERROR(check_input): input n_dust_size_qabs not allowed!'
                    print *, 'n_dust_size_qabs = ', n_dust_size_qabs
                    print *, 'allowed range = [10, Inf)'
                 endif
                 call stop_prc
              endif
           endif
        end do
        !check n_dust_wave_qabs
        if (n_dust_wave_qabs <= 10) then
           if (main_prc) then 
              print *, 'ERROR(check_input): input n_dust_wave_qabs not allowed!'
              print *, 'n_dust_wave_qabs = ', n_dust_wave_qabs
              print *, 'allowed range = [10, Inf)'
           endif
           call stop_prc
        endif
        ! check input Calorimetry files for the case of stochastically heated dust emission 
        if (dust_heating_type_ID == dt_sto .or. dust_heating_type_ID == dt_sto_lib) then  
           if (main_prc) print *, 'Checking existence of specific enthalphy files... '
           file_calorimetry_arr = (/file_calorimetry_Gra, file_calorimetry_Sil/)
           if (iq_dust_model(0) .or. any(iq_dust_model(2:3))) then  ! check Graphite/PAHs 
              call check_file_existence(dir_tmp, file_calorimetry_arr(0))
           endif
           if (iq_dust_model(1)) then    ! Check Silicates 
              call check_file_existence(dir_tmp, file_calorimetry_arr(1))
           end if
           ! check n_dust_temp_cal
           c1 = ((n_dust_temp_cal(0) < 100).and. (iq_dust_model(0) .or. any(iq_dust_model(2:3))))
           c2 = ((n_dust_temp_cal(0) < 100).and. (iq_dust_model(1)))
           
           if (c1 .or. c2) then
              if (main_prc) then 
                 print *, 'ERROR(check_input): input n_dust_temp_cal not allowed!'
                 print *, 'n_dust_temp_cal = ', n_dust_temp_cal
                 print *, 'allowed range = [10, Inf)'
              endif
              call stop_prc
           endif
            
        endif
     case default
        if (main_prc) print *, 'ERROR(check_input): dust_opacity_tables not recognized!'
        call stop_prc
     end select     
  else if (dust_model == 'TRUST') then
     iq_dust_model = [1,1,1,0]
     dust_opacity_tables = 'TRUST'
     if (main_prc) print *, 'Dust model as in the TRUST benchmark project will be used. To change this specify different dust_model!'
  else if (dust_model == 'DraineLi06') then
     iq_dust_model = [1,1,1,1]
     dust_opacity_tables = 'DraineLi06'
     if (main_prc) print *, 'Dust model as in Draine & Li (2006) will be used. To change this specify different dust_model!'   
  else
     if (main_prc) then 
        print *, 'ERROR: dust_model not recognized!'
        print *, 'dust_model = ', trim(dust_model)
     endif
     STOP
  endif

  ! input integrated/average opacities coefficients
  if (input_av_opacities) then
     dir_tmp = './'
     if (file_av_opacities /= 'not_provided') then
        call check_file_existence(dir_tmp, file_av_opacities)
     else
        if (main_prc) print *, 'ERROR: input_av_opacities = .TRUE. but file_av_opacities has not been provided!'
        STOP
     endif

  else
     if (main_prc) print *, 'WARNING: The integrated/average opacities coefficients will be calculated from the tabulated grain coefficients and the grain size distribution. To use input file instead, set input_av_opacities = .TRUE. and provide file_av_opacities.'
  endif


end subroutine check_input_dust_model

 !> This subroutine checks that the input file exists.
 subroutine check_file_existence(dir, filename)
   !> @param [in] filename Input filename
   character(LEN=lcar) :: filename
   !> @param [in] dir Directory containing the input file
   character(LEN=lcar) :: dir 
   !> @param file_exists Logical equal to TRUE if the input file exists.
   logical :: file_exists

   inquire(file=trim(adjustl(dir))//trim(adjustl(filename)), exist=file_exists)

   if (.not. file_exists) then
      print *, 'ERROR : file '//trim(adjustl(dir))//trim(adjustl(filename))//' not found!'
      call stop_prc
   endif

 end subroutine check_file_existence

 !> This subroutine checks that a directory exists. If the logical create_it is set to TRUE, it creates the directory in case it not exists. Otherwise it output an error and stop the program.
 subroutine check_dir_existence(dir,create_it)
   !> @param dir Input directory.
   character(LEN=lcar) :: dir
   !> @param create_it Logical equal to TRUE if non existent directory has to be created.
   logical :: create_it
   !> @param file_tmp Temporary file name used to store output linux commands.
   character(LEN=lcar) :: file_tmp
   !> @param inum Counter of lines in file_tmp.
   integer :: inum
   !> @param is String index
   integer :: is 
   !> @param unit file unit number
   integer :: unit
   integer (kind=int64) :: count, count_rate, count_max
 
   call system_clock(count, count_rate, count_max)

   write(file_tmp, *) count +id_mpi*1000  ! adding id_mpi makes the difference 
   file_tmp = 'dartray.tmp.'//trim(adjustl(file_tmp))
      
   ! fix possible problem with slash at the end of dir
   is = len_trim(dir(:))
   if (dir(is:is) == '/') then
      dir(is:is) = ''
   else
      is = is +1  ! select next position if dir does not end with "/"
   endif
    
   ! check dir exists (no standard Fortran way to do that. Use unix command find.
   
!  call system('find . -path '//trim(adjustl(dir))//' | wc -l > '//trim(adjustl(file_tmp)))
   call system('find '//trim(adjustl(dir))//' '//'-type d -wholename '//trim(adjustl(dir))//' | wc -l > '//trim(adjustl(file_tmp))) ! this way it works for directories outside the current directory as well (not just subdirectories). 
  
   open(newunit=unit, file=file_tmp,status='old')
   read(unit,*) inum
   close(unit)
   
   call system("rm "//trim(adjustl(file_tmp)))
   if (inum == 0) then
      print *,  'directory '//trim(adjustl(dir))//' not found!'
      if (create_it) then 
          print *, 'Create it!'
         call system('mkdir -p '//trim(adjustl(dir)))
      else
        ! print *, 'STOP'
         call stop_prc
      endif
         
   endif

   ! add final slash
   dir(is:is) = '/'

 end subroutine check_dir_existence
 
!> This subroutine sets the chuck size for the OPENMP loops in rt_loop(), rt_loop_2D(), rt_loop_i_obs()
subroutine set_chunk_size

  if (main_prc) print *, 'setting chunk_size.... '
  
  chunk_size = tot_ncell/nproc**2/np_mpi
  if (main_prc) print *, 'chunk size = ', chunk_size

  if (chunk_size == 0) then
     print *, 'ERROR chunck_size = 0! Something wrong here: set_chunk_size'
     CALL STOP_PRC
  endif

  call print_done
  
end subroutine set_chunk_size

!> Checks that the RAM memory of the machine is sufficient to host the big scaspe() arrays required in DART-Ray. WARNING: It assumes that all MPI processes run on different cluster nodes.  
subroutine check_memory

!> @param file_tmp Temporary file to store output linux commands
  character(LEN=lcar) :: file_tmp
  character(LEN=lcar) :: s1
  integer(kind=int64) :: max_mem_kb, needed_mem_kb
  integer (kind=int64) :: counts, count_rate, count_max
  integer :: id 
  integer :: il 
  integer :: i0, i1
  real(kind=real32), parameter :: ext_factor = 1.2
  real(kind=real32), parameter :: offset = 2E6
  
  ! set temporary file name
  if (main_prc) print *, 'Checking required RAM memory...'
  
  call set_i_opacity_arrays(i0,i1)

  call system_clock(counts, count_rate, count_max)

  write(file_tmp, *) counts +id_mpi*1000  ! adding id_mpi makes the difference 
  file_tmp = 'dartray.tmp.'//trim(adjustl(file_tmp))

  ! find out RAM memory size

  call system("cat /proc/meminfo > "//trim(adjustl(file_tmp)))

  open(newunit=id, file=file_tmp,status='old')
  read(id,*) s1, max_mem_kb
  close(id)
  
  call system("rm "//trim(adjustl(file_tmp)))

  ! calculate approximately memory required assuming mostly taken by scaspe arrays plus a standard 2 Gbytes (main grid + precalculated ads_arr arrays)
  needed_mem_kb = 0 
  do il = 0, lnum -1 
     needed_mem_kb = needed_mem_kb + npix_arr(il+i0)*tot_ncell  ! one scaspe array (without tot_ndir)
     
  enddo
 
  if (.not. sequential_scattering .and. .not. (rt_algorithm_ID == rta_2D .or. rt_algorithm_ID == rta_dust2D)) then
     needed_mem_kb = needed_mem_kb*2     ! multiplied by 2 (to include scaspe_arr and scaspe_tot)
  else
     needed_mem_kb = needed_mem_kb*3     ! multiplied by 3 (to include scaspe_arr and scaspe_tot and scaspe_prev)
  endif
     
  needed_mem_kb = needed_mem_kb*8    ! 8 bytes per real64

  needed_mem_kb = needed_mem_kb/1000   ! bytes -> kb

  if (.not. no_communications) then 
     needed_mem_kb = needed_mem_kb/np_mpi  ! divided by number of MPI processes
  endif

  needed_mem_kb = needed_mem_kb*ext_factor + offset
     
  if (needed_mem_kb  > max_mem_kb) then
     print *, 'Too much memory required. Try to reduce kp_sca_max or number of input wavelengths'
     print *, 'Extimated needed memory [kb] = ', needed_mem_kb
     print *, 'Total memory available [kb] = ', max_mem_kb
     call stop_prc
  endif

  if (main_prc) print *, needed_mem_kb, ' [kb]'

  call print_done 

end subroutine check_memory

!> Initialises MPI environment.  
subroutine initialize_mpi

  integer :: required, provided 
  integer :: ierr, hl,stat
  character*(MPI_max_processor_name) :: hostname 
  character (LEN=255) :: nproc_num

  !required = MPI_THREAD_SERIALIZED
  !required = MPI_THREAD_MULTIPLE
  required = MPI_THREAD_FUNNELED

  call mpi_init_thread(required, provided,ierr)

  call mpi_comm_size(mpi_comm_world,np_mpi,ierr)

  call mpi_comm_rank(mpi_comm_world,id_mpi,ierr)

  call mpi_get_processor_name(hostname, hl,ierr)

  ! assign main_prc
  if (id_mpi == 0) main_prc = .TRUE.

  ! print number of MPI processes
  if (main_prc) print *, 'number of MPI processes =', np_mpi

  ! check required MPI thread environment 
    if (required /= provided) then 
       if (main_prc) print *, 'ERROR : Required thread environment (MPI_THREAD_FUNNELLED) not allowed.'
     call stop_prc
  endif

  ! get openmp thread number 
  call get_environment_variable('OMP_NUM_THREADS',nproc_num, status=stat )
  
  if (stat /= 0) then 
     print *, 'ERROR: OMP_NUM_THREADS variable not found!'
     stop
  endif

  ! set OpenMP thread number 
  read(nproc_num, *) nproc 

  ! set starting time 
  time_start = mpi_wtime()
  

end subroutine initialize_mpi

!> Exits the MPI environment 
subroutine terminate_mpi
  integer :: ierr

   call mpi_finalize(ierr)

end subroutine terminate_mpi

!> Handles all the collective MPI operations at the end of each RT step  
subroutine rt_mpi
  integer :: i 
  
  select case (rt_type)

  case(rtt_precalc_src)
     call reduce_u_fest_arr
        
  case(rtt_dir_src)
     call reduce_u_final_arr

     if (no_communications) then

        call reduce_scaspe_arr
              
     endif
        
     call reduce_lum_lost
     
     call reduce_psel_av_arr
        
  case(rtt_scatt)
     call reduce_u_final_arr

     if (no_communications) then 
        call reduce_scaspe_arr
     endif
        
     call reduce_lum_lost

     call reduce_psel_av_arr
        
  case default
        print *, 'You should never get here!'
        call stop_prc
        
     end select

end subroutine rt_mpi

!> Reduces u_fest_arr() stored in all MPI processes.
subroutine reduce_u_fest_arr
  integer :: tot_el
  integer :: ierr
  real(kind=real64), allocatable :: in_arr(:,:)

  tot_el = lnum*tot_ncell
  
  ! Reduce u_fest_arr and broadcast to all processes
        
  allocate(in_arr(0:lnum-1,0:tot_ncell-1))
  in_arr = 0
               
  call mpi_reduce(u_fest_arr, in_arr, tot_el, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr)
        
  if (main_prc) then 
     u_fest_arr = in_arr
  endif
  deallocate(in_arr)
        
  call mpi_Bcast(u_fest_arr,tot_el, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,ierr)
  
end subroutine reduce_u_fest_arr

!> Reduces u_final_arr() stored in all MPI processes.
subroutine reduce_u_final_arr
  integer :: tot_el
  integer :: ierr
  real(kind=real64), allocatable :: in_arr(:,:)

  tot_el = lnum*tot_ncell
  
  ! Reduce u_final_arr         
  allocate(in_arr(0:lnum-1,0:tot_ncell-1))
  in_arr = 0

  if (rt_type == rtt_scatt .or. iterations_dustem > 1) then
     u_final_arr = u_final_arr - u_fest_arr
     where(u_final_arr < 0) ! this is because of numerical accuracy. Sometimes the arrays can be very close and the difference is negative while it should be zero. 
        u_final_arr = 0
     end where

  endif
     
  call mpi_allreduce(u_final_arr, in_arr, tot_el, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)

  if (rt_type /= rtt_scatt .and. iterations_dustem <= 1) then 
     u_final_arr = in_arr
  else
     u_final_arr = in_arr + u_fest_arr
  endif
     
  deallocate(in_arr)

end subroutine reduce_u_final_arr

!> Reduces dens_stars_arr() stored in all MPI processes. Used only in dust RT algorithms. 
subroutine reduce_dens_stars_arr
  integer :: tot_el
  integer :: ierr
  real(kind=real64), allocatable :: in_arr(:,:)

  tot_el = lnum*tot_ncell
  
  ! Reduce dens_stars_arr         
  allocate(in_arr(0:lnum-1,0:tot_ncell-1))
  in_arr = 0
     
  call mpi_allreduce(dens_stars_arr, in_arr, tot_el, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)

     dens_stars_arr = in_arr
     
  deallocate(in_arr)

end subroutine reduce_dens_stars_arr

!> Reduces tot_dust_em_sed() stored in all MPI processes. Used only in dust RT algorithms. 
subroutine reduce_tot_dust_em_sed
  integer :: tot_el
  integer :: ierr
  real(kind=real64), allocatable :: in_arr(:,:,:)

  tot_el = lnum*n_int_rf_bins**2
  
  ! Reduce dens_stars_arr         
  allocate(in_arr(0:lnum-1,0:n_int_rf_bins-1, 0:n_int_rf_bins-1))
  in_arr = 0
     
  call mpi_allreduce(tot_dust_em_sed, in_arr, tot_el, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)

     tot_dust_em_sed = in_arr
     
  deallocate(in_arr)

end subroutine reduce_tot_dust_em_sed



!> Reduces scaspe_arr() stored in all MPI processes. 
subroutine reduce_scaspe_arr 
integer :: ierr
real(kind=real64), allocatable :: out_arr(:,:), in_arr(:,:)
integer :: i 
integer :: i0, i1 
  
call set_i_opacity_arrays(i0,i1)

! Reduce scaspe_arr if no_communications is set. Since this requires lots of memory, it is done one wavelength at the time.

! scaspe_arr          
do i = 0, lnum -1 

   allocate(out_arr(0:npix_arr(i+i0)-1,0:tot_ncell-1), in_arr(0:npix_arr(i+i0)-1,0:tot_ncell-1))

   out_arr = scaspe_arr(i)%a(:,:)
   in_arr = 0
   
   call mpi_allreduce(out_arr, in_arr, npix_arr(i+i0)*tot_ncell, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)

   scaspe_arr(i)%a(:,:) = in_arr 
     
   deallocate(out_arr, in_arr)

end do


  
end subroutine reduce_scaspe_arr

!> Reduces lum_lost() stored in all MPI processes.
subroutine reduce_lum_lost 
integer :: ierr
real(kind=real64), allocatable :: in_arr_1d(:)

! Reduce lum_lost
        
allocate(in_arr_1d(0:lnum-1))
in_arr_1d = 0

if (rt_type == rtt_scatt .or. rt_algorithm_ID == rta_dust .or. rt_algorithm_ID == rta_dust2D) then
   lum_lost = lum_lost - (np_mpi-1._real64)/np_mpi*lum_lost_prev
endif
                
call mpi_allreduce(lum_lost, in_arr_1d, lnum, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
        
lum_lost = in_arr_1d
deallocate(in_arr_1d)
        
end subroutine reduce_lum_lost


!> Reduces psel_av_arr() stored in all MPI processes.
subroutine reduce_psel_av_arr 
integer :: ierr
real(kind=real64), allocatable :: in_arr(:,:,:)
integer :: i0 

if (.not. print_psel_av) return

! Reduce psel_av_arr  but only to master process (which prints the array
        
allocate(in_arr(0:0,0:0,0:tot_ncell+tot_p_src-1))                 
in_arr = 0

i0 = iterations_dustem-1
if (iterations_dustem == 0) i0 = 0 ! this is because iterations_dustem is equal to 0 during the stellar emission RT. 
                
call mpi_reduce(psel_av_arr(i0,iterations,:), in_arr, tot_ncell+tot_p_src, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr)

if (main_prc) then 
   psel_av_arr(i0,iterations,:) = in_arr(0,0,:)
endif
deallocate(in_arr)

end subroutine reduce_psel_av_arr


!> Reduces i_obs_arr() stored in all MPI processes. NOT used for the moment. 
subroutine reduce_i_obs_arr 
integer :: ierr
real(kind=real64), allocatable :: out_arr(:,:), in_arr(:,:)
integer :: i


! i_obs_arr
if (tot_ndir > 0) then 
   allocate(out_arr(0:tot_ncell_p_src-1, 0:tot_ndir-1), in_arr(0:tot_ncell_p_src-1,0:tot_ndir-1))
           
   do i = 0, lnum -1 
                 
      out_arr = i_obs_arr(i, 0:tot_ncell_p_src-1, 0:tot_ndir-1)
      in_arr = 0
      
      call mpi_allreduce(out_arr, in_arr, tot_ncell_p_src*tot_ndir, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
              
      i_obs_arr(i, 0:tot_ncell_p_src-1, 0:tot_ndir-1) = in_arr 
              
   end do
   
   deallocate(out_arr, in_arr)

endif

end subroutine reduce_i_obs_arr


!> Reduces i_obs_in_arr() stored in all MPI processes. NOT used for the moment.
subroutine reduce_i_obs_in_arr 
integer :: ierr
real(kind=real64), allocatable :: out_arr(:,:), in_arr(:,:)
integer :: i


! i_obs_in_arr
if (tot_ndir_in > 0) then 
   allocate(out_arr(0:tot_ncell_p_src-1, 0:tot_ndir_in-1), in_arr(0:tot_ncell_p_src-1,0:tot_ndir_in-1))
           
   do i = 0, lnum -1 

      out_arr = i_obs_in_arr(i, 0:tot_ncell_p_src-1, 0:tot_ndir_in-1)
      in_arr = 0
              
      call mpi_allreduce(out_arr, in_arr, tot_ncell_p_src*tot_ndir_in, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)

      i_obs_in_arr(i, 0:tot_ncell_p_src-1, 0:tot_ndir_in-1) = in_arr 
              
   end do

   deallocate(out_arr, in_arr)
endif


end subroutine reduce_i_obs_in_arr


!> Reduces map_arr_out() stored in all MPI processes.
subroutine reduce_map_arr_out
  integer :: tot_el
  integer :: ierr
  real(kind=real64), allocatable :: in_arr(:,:,:,:)

  if (np_mpi == 1) return ! no need to reduce if just one MPI process

  tot_el = npixel_maps**2*lnum_maps*tot_ndir
  
  ! Reduce u_fest_arr and broadcast to all processes

  allocate(in_arr(0:npixel_maps-1, 0:npixel_maps-1, 0:lnum_maps-1, 0:tot_ndir-1))
  in_arr = 0
               
  call mpi_reduce(map_arr_out, in_arr, tot_el, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr)
        
  if (main_prc) then 
     map_arr_out = in_arr
  endif
  deallocate(in_arr)
  
end subroutine reduce_map_arr_out

!> Reduces map_in_arr_out() stored in all MPI processes.
subroutine reduce_map_in_arr_out
  integer :: tot_el
  integer :: ierr
  real(kind=real64), allocatable :: in_arr(:,:,:)

  tot_el = npix_maps*lnum_maps*tot_ndir_in
  
  ! Reduce u_fest_arr and broadcast to all processes
        
  allocate(in_arr(0:npix_maps-1, 0:lnum_maps-1, 0:tot_ndir_in-1))
  in_arr = 0
               
  call mpi_reduce(map_in_arr_out, in_arr, tot_el, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr)
        
  if (main_prc) then 
     map_in_arr_out = in_arr
  endif
  deallocate(in_arr)
  
end subroutine reduce_map_in_arr_out



!> Reads input file for RT calculation. The input file has follow the dartray command on the terminal prompt. 
subroutine read_input_file
integer :: ioun
character(len=lcar) :: input_filename

call get_command_argument(1,input_filename)

!check input_file argument has been given
if (trim(input_filename) == '') then 
   if (main_prc) print *, 'provide input file in command line...EXIT' 
   call stop_prc
endif

! allocate ind_i_obs to 1000 elements 
allocate(ind_i_obs(0:max_n_input-1), ind_out_maps(0:max_n_input-1))
ind_i_obs = -1 
ind_out_maps = -1 

! read input file 
if (main_prc) print *, 'reading input file...'

open(newunit=ioun, file = input_filename, status ='old', action ='read')
read(ioun, nml = dartray_input_strings)
read(ioun, nml = dartray_input_var)
read(ioun, nml = dartray_input_logical)
close(ioun)

call print_done

end subroutine read_input_file

!> Reads an Nbody-SPH simulation file ( file_nbody_sph()).
subroutine read_nbody_sph_simulation
  USE ISO_C_BINDING
  INTEGER(HID_T) :: file_id                            ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
  integer, parameter :: narr_grid_nbody_sph = 9
  CHARACTER (LEN=lcar) :: dsetname(0:narr_grid_nbody_sph-1) ! Dataset name
  INTEGER     ::   rank(0:narr_grid_nbody_sph-1)   
  INTEGER(HSIZE_T) :: num
  INTEGER     ::   error  ! Error flag
  integer :: i
  TYPE(C_PTR) :: f_ptr
  real(kind=real64) :: x0,x1, y0,y1,z0,z1

  if (main_prc) print *, 'reading N-body/SPH simulation file....'
    
  dsetname(0)='mstar'     ; rank(0)=1 
  dsetname(1)='starcoord'    ; rank(1)=2 
  dsetname(2)='agestar'    ; rank(2)=1 
  dsetname(3)='fehstar'    ; rank(3)=1 
  dsetname(4)='mgas'     ; rank(4)=1 
  dsetname(5)='gascoord'       ; rank(5)=2    
  dsetname(6)='gastemp'      ; rank(6)=1 
  dsetname(7)='ofegas'  ; rank(7)=1
  dsetname(8)='fehgas'      ; rank(8)=1
  
  CALL h5open_f (error)
  CALL h5fopen_f (trim(adjustl(dir_grid))//trim(adjustl(file_nbody_sph)), H5F_ACC_RDWR_F, file_id, error)
  
  if (error /= 0.) then 
     if (main_prc) then 
        print *,  'N-body/SPH simulation file not found!'
        print *, 'file =', trim(adjustl(dir_grid))//trim(adjustl(file_nbody_sph))
     endif
     call stop_prc
  endif

  do i=0,narr_grid_nbody_sph-1
     
     CALL h5dopen_f(file_id, dsetname(i), dset_id, error)
     CALL h5dget_space_f(dset_id, dspace_id, error) 
     
     call h5sget_simple_extent_npoints_f(dspace_id, num, error)

     if (i == 0) then  ! allocate stellar particles arrays
        tot_star_particles=num
        allocate(mstar(0:tot_star_particles-1), starcoord(3,0:tot_star_particles-1), agestar(0:tot_star_particles-1), fehstar(0:tot_star_particles-1), star_lum(0:tot_star_particles-1), pcell_star(0:tot_star_particles-1))
        star_lum = 0 ! this is the only array not assigned here.
        pcell_star = -1 ! this -1 is important. see user_routines_Nbody_SPH
     elseif (i>0 .and. i < 4) then ! check array size is = tot_star_particles or tot_star_particles*3
        if (num /= tot_star_particles .and. num /= tot_star_particles*3) then 
           if (main_prc) print *, 'Stellar particle arrays in the Nbody/SPH simulation file do not have the same size'
           STOP
        endif
     elseif (i == 4) then 
        tot_gas_particles=num
        allocate(mgas(0:tot_gas_particles-1), gascoord(3,0:tot_gas_particles-1), gastemp(0:tot_gas_particles-1), fehgas(0:tot_gas_particles-1), ofegas(0:tot_gas_particles-1), pcell_gas(0:tot_gas_particles-1))
        pcell_gas = -1  ! this -1 is important. see user_routines_Nbody_SPH
     elseif (i > 4) then 
        if (num /= tot_gas_particles .and. num /= tot_gas_particles*3) then 
           if (main_prc) print *, 'Gas particle arrays in the Nbody/SPH simulation file do not have the same size'
           STOP
        endif
     endif
     
     allocate(dims(rank(i)))
     
     if (dsetname(i) == 'mstar' ) then
        dims = (/tot_star_particles/)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , mstar, dims, error)
     else if (dsetname(i) == 'starcoord' ) then
        dims = (/3, tot_star_particles/)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , starcoord, dims, error)
     else if (dsetname(i) == 'agestar' ) then
        dims = (/tot_star_particles/)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , agestar, dims, error)
     else if (dsetname(i) == 'fehstar' ) then
        dims = (/tot_star_particles/)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , fehstar, dims, error)   
     else if (dsetname(i) == 'mgas' ) then
        dims = (/tot_gas_particles/)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , mgas, dims, error)
     else if (dsetname(i) == 'gascoord' ) then
        dims = (/3,tot_gas_particles/)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , gascoord, dims, error)   
     else if (dsetname(i) == 'gastemp' ) then
        dims = (/tot_gas_particles/)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , gastemp, dims, error)
     else if (dsetname(i) == 'ofegas' ) then
        dims = (/tot_gas_particles/)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , ofegas, dims, error)
     else if (dsetname(i) == 'fehgas' ) then
        dims = (/tot_gas_particles/)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , fehgas, dims, error)
     endif
     
     CALL h5sclose_f(dspace_id, error)
     CALL h5dclose_f(dset_id, error)
     
     if (main_prc) print *, ' '
     if (main_prc) print *, trim(dsetname(i)), ' loaded'
     if (main_prc) print *, 'dimensions = ', dims 
     
     deallocate(dims)
     
  end do
    
  CALL h5fclose_f(file_id, error)
  CALL h5close_f(error) 

  ! some transformations 
  starcoord = starcoord*1E3  ! kpc -> pc 
  gascoord = gascoord*1E3


  if (main_prc) then 
     print *, 'Stellar particle coordinate range [pc]'
     x0=minval(starcoord(1,0:tot_star_particles-1)) ; x1=maxval(starcoord(1,0:tot_star_particles-1))
     print *, 'X', x0,x1

     y0=minval(starcoord(2,0:tot_star_particles-1)) ; y1=maxval(starcoord(2,0:tot_star_particles-1))
     print *, 'Y', y0,y1

     z0=minval(starcoord(3,0:tot_star_particles-1)) ; z1=maxval(starcoord(3,0:tot_star_particles-1))
     print *, 'Z', z0,z1
     print *, '' 
     print *, 'Gas particle coordinate range [pc]'
     x0=minval(gascoord(1,0:tot_gas_particles-1)) ; x1=maxval(gascoord(1,0:tot_gas_particles-1))
     print *, 'X', x0,x1

     y0=minval(gascoord(2,0:tot_gas_particles-1)) ; y1=maxval(gascoord(2,0:tot_gas_particles-1))
     print *, 'Y', y0,y1

     z0=minval(gascoord(3,0:tot_gas_particles-1)) ; z1=maxval(gascoord(3,0:tot_gas_particles-1))
     print *, 'Z', z0,z1

  endif
 
  call print_done 

end subroutine read_nbody_sph_simulation

!> reads the input stellar library (see file_stellar_library()). 
subroutine read_stellar_library 
  USE ISO_C_BINDING
  INTEGER(HID_T) :: file_id                            ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
  integer, parameter :: narr_grid_lib = 4
  CHARACTER (LEN=lcar) :: dsetname(0:narr_grid_lib-1) ! Dataset name
  INTEGER     ::   rank(0:narr_grid_lib-1)   
  INTEGER(HSIZE_T) :: num
  INTEGER     ::   error  ! Error flag
  integer :: i
  
  if (main_prc) print *, 'reading stellar library file....'

  select case (stellar_library)
  case('starburst99')
     file_stellar_library = './STELLAR_LIBRARIES/starburst99/table_lum_mass_vs_age_met_starburst99.h5'
  case('maraston05_kr_rhb')
     file_stellar_library = './STELLAR_LIBRARIES/maraston2005/table_lum_mass_vs_age_met_maraston2005_kr_rhb.h5'
  case('user')
     file_stellar_library = file_stellar_library
  case default 
     if (main_prc) then 
        print *, 'ERROR(read_stellar_library): stellar_library not defined!'
        print *, 'stellar_library = ', stellar_library
     endif
     call stop_prc
     
  end select

  dsetname(0)='lambda_lib'     ; rank(0)=1 
  dsetname(1)='met_lib'    ; rank(1)=1
  dsetname(2)= 'age_lib'   ; rank(2)=1
  dsetname(3)='lum_to_mass_lib'    ; rank(3)=3  

  CALL h5open_f (error)
  CALL h5fopen_f (trim(adjustl(file_stellar_library)), H5F_ACC_RDWR_F, file_id, error)
  
  if (error /= 0.) then 
     if (main_prc) then 
        print *,  'ERROR(read_stellar_library): file not found!'
        print *, 'file =', trim(adjustl(file_stellar_library))
     endif
     call stop_prc
  endif

  do i=0,narr_grid_lib-1
     
     CALL h5dopen_f(file_id, dsetname(i), dset_id, error)
     CALL h5dget_space_f(dset_id, dspace_id, error) 
     
     call h5sget_simple_extent_npoints_f(dspace_id, num, error)

     if (i == 0) then  ! allocate stellar particles arrays
       nlambda_lib = num
       allocate(lambda_lib(0:nlambda_lib-1))
    elseif (i == 1) then
       nmet_lib = num
       allocate(met_lib(0:nmet_lib-1))
    elseif (i == 2) then
       nage_lib = num
       allocate(age_lib(0:nage_lib-1))
    elseif (i == 3) then 
       !nage_lib = num/nlambda_lib/nmet_lib 
       allocate(lum_to_mass_lib(0:nlambda_lib-1, 0:nage_lib-1, 0:nmet_lib-1))
    endif
     
    allocate(dims(rank(i)))
     
    if (dsetname(i) == 'lambda_lib' ) then
        dims = (/nlambda_lib/)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , lambda_lib, dims, error)
     else if (dsetname(i) == 'met_lib' ) then
        dims = (/nmet_lib/)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , met_lib, dims, error)
     else if (dsetname(i) == 'age_lib' ) then
        dims = (/nage_lib/)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , age_lib, dims, error)     
     else if (dsetname(i) == 'lum_to_mass_lib' ) then
        dims = (/nlambda_lib, nage_lib, nmet_lib/)
        CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE , lum_to_mass_lib, dims, error)
     endif
     
     CALL h5sclose_f(dspace_id, error)
     CALL h5dclose_f(dset_id, error)
     
     if (main_prc) print *, ' '
     if (main_prc) print *, trim(dsetname(i)), ' loaded'
     if (main_prc) print *, 'dimensions = ', dims 
     
     deallocate(dims)
     
  end do
    
  CALL h5fclose_f(file_id, error)
  CALL h5close_f(error) 
 
  call print_done


end subroutine read_stellar_library

!> Prints the time passed from the beginning of the RT calculation in seconds. 
subroutine print_time
real(kind=real64) :: delta_time

time =  mpi_wtime()

delta_time = time-time_start

if (main_prc) print *, 'TIME [s] = ', delta_time

end subroutine print_time



END MODULE io_routines
