!> Contains the main RT routines. In particular, the routines to calculate the ray-cell intersections, while a ray propagates through the 3D grid, and the corresponding contributions to the radiation field energy density and scattered radiation.  
MODULE rt_routines
  use smooth_grid_routines
  use healpix_routines
  use ray_list
  use omp_lib
  USE ISO_C_BINDING
  IMPLICIT NONE
  !> @param kp_sca_max HEALPix parameter used to determine the maximum allowed number of "spherical pixels" in the scaspe() arrays: npix_max = 12*2^(2*kp_sca_max). The allowed range is [0,4]. Values higher than 4 are not allowed because the scaspe arrays would become too big. However, this can be changed in check_input(). 
  INTEGER(kind=int32) :: kp_sca_max

  !> @param nside_sca HEALPix nside parameter. Used when calling some HEALPix subroutines. 
  INTEGER(kind=int32) :: nside_sca
  
  !> @param bm_par Minimum number of ray-intersections per cell during the pre-iteration and direct-light processing for the stellar emission RT. Allowed range = [0,1000].
  INTEGER(kind=int32) :: bm_par
  !> @param bm_par_sca Minimum number of ray-intersections per cell during the scattered light processing for the stellar emission RT and all phases for the dust emission RT. Allowed range = [0,1000]. 
  INTEGER(kind=int32) :: bm_par_sca
  !> @param bm_par_max Maximum number of ray-intersections per cell. Allowed range = [100*max([bm_par,bm_par_sca]), inf]'
  INTEGER(kind=int32) :: bm_par_max
  !> @param pabs_max Maximum value allowed for pabs variable in calc_psel(). Set to 10* modelsize()  
  INTEGER(kind=int32) :: pabs_max
  !> @param rad_lim Maximum distance relative to modelsize() that can be crossed by the rays in the precalculation phase. Allowed range = [0,2]. 
  real(KIND=real64) :: rad_lim
  !> @param en_lim This is the parameter \f$ f_U \f$ of the DART-Ray algorithm. It is used when evaluating \f$ \delta U_\lambda < f_U*U_{\lambda, LL} \f$
  real(KIND=real64) :: en_lim
  !> @param conv_en_lim This is the parameter \f$ f_L \f$ used to determine when to stop the scattering iterations. This happens when at the end of an iteration the remaining scattered luminosity to be processed is less than \f$ f_L \f$ times the total scattered luminosity (see tot_rad_en()). 
  real(KIND=real64) :: conv_en_lim
  !> @param tau_cell_max Maximum optical depth allowed for leaf cells in the grid creation program. Allowed range = [0,0.5].   
  real(KIND=real64) :: tau_cell_max
  !> @param accuracy Accuracy parameter determining value of en_lim() (that is, \f$ f_U \f$) through the formula en_lim() = accuracy()/( tot_sources*0.25) with  tot_sources = sum(leaf_cells) + tot_p_src() - tot_spare_cells(). Allowed range = [0,1].
  real(KIND=real64) :: accuracy
  !> @param tot_rad_en Total luminosity for all wavelengths. During precalculation and direct-light processing is the total luminosity of the stars. During the scattering iteration is the total scattered luminosity. 
  real(kind=real64), allocatable :: tot_rad_en(:)
  !> @param tot_rad_en_or Total stellar luminosity for all wavelengths. Used to determine lost luminosity fraction in print_lum_lost(). If this fraction is higher than 1%, a warning is output on the terminal. 
  real(KIND=real64), allocatable :: tot_rad_en_or(:)
  !> @param rt_type Integer variable used to select the different types of RT algorithm steps (e.g. precalculation or direct light processing). See select_rt_type().
  integer :: rt_type
  !> @param rt_algorithm Type of RT algorithm to be used. Choices are: 
!> 'main': this is the standard RT algorithm. Firstly, it proceeds with the calculation of the radiation fields, scattering source function, outgoing radiation specific intensity, surface brightness maps for the stellar emission. Then, it performs the same for the dust emission (unless no_dust_rt() = .TRUE.);
!> '2D': same as the 'main' RT algorithm but for axisymmetric models. In this case, the code checks that the input grid is indeed axisymmetric and then performs the RT calculations taking into account the symmetries. This mode performs RT calculations up to a factor 8 faster then the main mode; 
!> 'dust' and 'dust_2D': these RT algorithms perform the dust emission RT calculations. To use these modes, it is required that the stellar emission RT has already been performed and the output radiation fields are printed on disk. The two modes are for the normal and 2D algorithm respectively;
!> 'sed' and 'sed_dust': these modes can be used to calculate the integrated SEDs and the maps using the outgoing specific intensity files ( file_i_obs() and file_i_obs_in()) printed in a previous RT calculation. The two modes are for the stellar emission and the dust emission respectively;  
!> 'i_obs' and 'i_obs_dust': these modes can be used to calculate the outgoing specific intensity and the corresponding surface brightness maps for arbitrary observer positions without repeating a complete RT calculation. These modes can be used only if the scattering source function has be printed on disk (see print_scaspe_tot()). The two modes are for the stellar emission and dust emission respectively;
!>  'projection': this mode can be used to make maps of the stellar emission without any dust or maps of the dust optical depth. Useful to check that the input grid has been correctly calculated before starting a complete RT calculation.  
  character(LEN = 20) :: rt_algorithm
  !> @param rt_algorithm_ID Integer variable associated with the rt_algorithm() choices. See check_input() for definitions. 
  integer :: rt_algorithm_ID 
  !> @param done Logical variable equal to TRUE when the RT algorithm run has been completed.  
  logical :: done
  !> @param cnflag Logical variable equal to ( tot_lumcell/ tot_rad_en() < conv_en_lim()) witn tot_lumcell=sum( lumcell()). When TRUE, the scattered luminosity left for processing is too small and the scattering iterations are stopped. 
  logical :: cnflag
  !> @param cnflag_dust Logical variable used to stop the dust heating iterations. Equal to TRUE when dens_stars_arr() is a small fraction conv_en_lim() of the comulated emission stored in dens_stars_arr_old().   
  logical :: cnflag_dust  
  !> @param nside HEALPix nside parameter.
  integer :: nside
  !> @param clvl_old Subdivision level of the last intersected cell. 
  integer :: clvl_old
  !> @param nside_min HEALPix minimum nside parameter used to select the ray angular density in main_dir_loop().
  integer, parameter :: nside_min= 4 ! this cannot be changed (see create_high_ray_list)!!!
  !> @param npix_main Number of spherical pixels in the HEALPix sphere for nside()=1. 
  integer, parameter :: npix_main=12
  !> @param src_ccoord(3) Radiation source 3D position (it can be either a cell centre or a point source position). 
  real(kind=real64) :: src_ccoord(3)
  !> @param ccindd_nc(3,max_lvl) Tree coordinate position of the cell originating the ray. 
  integer, allocatable :: ccindd_nc(:,:)
  !> @param ccindd(3,max_lvl) Tree coordinate position of a cell (usually the intersected cell). 
  integer, allocatable :: ccindd(:,:)
  !> @param nproc Number of processors to be used for OPENMP loops. 
  integer :: nproc
  !> @param ipsel_av Counter of rays launched by the same cell. Used in the calculation of the average ray path psel_av_arr().
  integer :: ipsel_av
  !> @param ipsel_av_tot Counter of total number of rays 
  integer :: ipsel_av_tot
  !> @param glepsilon Very small factor used to determine the ray direction (positive or negative). IS THERE A BETTER WAY TO DO THIS ? 
  real (kind=real64), parameter :: glepsilon=1.0e-7
  
  !> @param ads_arr4(0:dim_npix_unique-1)%a(0:npix_unique-1, 0:npix_rays(0)-1). Precalculated ads_arr() factor for nside() =4. 
  type(var_arr_2d),allocatable :: ads_arr4(:)
  !> @param ads_arr8(0:dim_npix_unique-1)%a(0:npix_unique-1, 0:npix_rays(1)-1). Precalculated ads factor for nside() = 8. 
  type(var_arr_2d),allocatable :: ads_arr8(:)
  !> @param ads_arr16(0:dim_npix_unique-1)%a(0:npix_unique-1, 0:npix_rays(2)-1). Precalculated ads factor for nside() = 16. 
  type(var_arr_2d),allocatable :: ads_arr16(:)
  !> @param ads_arr32(0:dim_npix_unique-1)%a(0:npix_unique-1, 0:npix_rays(2)-1). Precalculated ads factor for nside() = 32. 
  type(var_arr_2d),allocatable :: ads_arr32(:)
  !> @param ads_arr64(0:dim_npix_unique-1)%a(0:npix_unique-1, 0:npix_rays(2)-1).  Precalculated ads factor for nside() = 64. 
  type(var_arr_2d),allocatable :: ads_arr64(:)
  !> @param ads_arr128(0:dim_npix_unique-1)%a(0:npix_unique-1, 0:npix_rays(2)-1). Precalculated ads factor for nside() = 128. 
  type(var_arr_2d),allocatable :: ads_arr128(:)
  !> @param ads_arr256(0:dim_npix_unique-1)%a(0:npix_unique-1, 0:npix_rays(2)-1). Precalculated ads factor for nside() = 256. 
  type(var_arr_2d),allocatable :: ads_arr256(:)

  !> @param theta_sca(0:dim_npix_unique-1)%a(0:npix_unique(i)-1) Array of theta angles (n^z) corresponding to the elements of the first index of the scaspe() arrays. 
  type(var_arr_1d), allocatable :: theta_sca(:)
  
  !> @param phi_sca(0:dim_npix_unique-1)%a(0:npix_unique(i)-1) Arrays of phi angles (n^xy) corresponding to the elements of the first index of the scaspe() arrays. 
  type(var_arr_1d), allocatable :: phi_sca(:)

  !> @param sin_theta_sca(0:dim_npix_unique-1)%a(0:npix_unique(i)-1) Array equal to sin(theta_sca()). Used when calculating ads_arr(). 
  type(var_arr_1d), allocatable :: sin_theta_sca(:)
  
  !> @param sin_phi_sca(0:dim_npix_unique-1)%a(0:npix_unique(i)-1) Array equal to sin( phi_sca()). Used when calculating ads_arr().
  type(var_arr_1d), allocatable :: sin_phi_sca(:)

  !> @param cos_theta_sca(0:dim_npix_unique-1)%a(0:npix_unique(i)-1) Array equal to cos( theta_sca()). Used when calculating ads_arr().
  type(var_arr_1d), allocatable :: cos_theta_sca(:)

  !> @param cos_phi_sca(0:dim_npix_unique-1)%a(0:npix_unique(i)-1) Array equal to sin( phi_sca()). Used when calculating ads_arr().
  type(var_arr_1d), allocatable :: cos_phi_sca(:)

  !> @param ffn_arr_mw(0:lnum_node-1)%a(0:npix_arr(i)-1). Scattering phase function values at the angular directions defined for the scaspe() arrays for a given ray direction. A variable length data type is used. The data datype argument contains the values for each line-of-sight direction (npix_arr(i) values for each wavelength i).
  type(var_arr_1d), allocatable :: ffn_arr_mw(:)
  
  !> @param ffn_arr(0:dim_npix_unique-1)%a(0:npix_unique(i)-1). This is an array used during the calculation of ffn_arr_mw at a single wavelength component. It has to be assigned outside the subroutine to avoid repeated allocation.  
  type(var_arr_1d), allocatable :: ffn_arr(:)

  !> @param ads_arr(0:dim_npix_unique-1)%a(0:npix_unique-1) Array of cos(theta) used when calculating the values of the Henyey-Greenstein function for a certain ray direction and the angular directions of the scaspe() arrays. 
  type(var_arr_1d), allocatable :: ads_arr(:)
  
  !> @param cs Light speed defined in the input.
  real(kind=real64) :: cs
  !> @param psel_av Average ray path for the rays cast from a certain cell.
  real(kind=real64) :: psel_av
  !> @param pabs_arr(3) Array storing some constant factors when calculating the ray-cell intersection in calc_psel(). These factors are the same as in the previous intersected cell if the subdivision level, stored in clvl_old(), does not change. 
  real(kind=real64) :: pabs_arr(3)
  !> @param i_obs_temp Temporary value of specific intensity when calculating i_obs(). It has a different value in each OPENMP thread. 
  real(kind=real64), allocatable :: i_obs_temp(:)
  !> @param lum_lost_temp Temporary value of the luminosity not processed within a HEALPix section (currently corresponding to nside = 4). These values accumulate on lum_lost().
  real(kind=real64), allocatable :: lum_lost_temp(:)
  !> @param lum_lost Total luminosity not processed in the entire RT run. SHARED
  real(kind=real64), allocatable :: lum_lost(:)
  !> @param lum_lost_prev Total luminosity not processed during the direct-light processing or the previous scattering iteration. SHARED
  real(kind=real64), allocatable :: lum_lost_prev(:)
  !> @param iterations Counter of scattering iterations. 
  integer :: iterations
  !> @param iterations_dustem Counter of dust self-heating iterations. 
  integer :: iterations_dustem
  !> @param vec_mod Distance between radiation source and the internal observer. 
  real(KIND=real64) :: vec_mod
  !> @param ix(0:dim_npix_unique-1)%a(0:npix_hp_unique(i)-1) Scaspe() array indices after X symmetry transformation. SHARED
  type(var_arr_1d), allocatable :: ix(:)
  !> @param iy(0:dim_npix_unique-1)%a(0:npix_hp_unique(i)-1) Scaspe() array indices after Y symmetry transformation. SHARED
  type(var_arr_1d), allocatable :: iy(:)
  !> @param iz(0:dim_npix_unique-1)%a(0:npix_hp_unique(i)-1) Scaspe() array indices after Z symmetry transformation. SHARED
  type(var_arr_1d), allocatable :: iz(:)
  !> @param ixy(0:dim_npix_unique-1)%a(0:npix_hp_unique(i)-1) Scaspe() array indices after XY symmetry transformation. SHARED
  type(var_arr_1d), allocatable :: ixy(:)
  !> @param ixz(0:dim_npix_unique-1)%a(0:npix_hp_unique(i)-1) Scaspe() array indices after XZ symmetry transformation. SHARED
  type(var_arr_1d), allocatable :: ixz(:)
  !> @param iyz(0:dim_npix_unique-1)%a(0:npix_hp_unique(i)-1) Scaspe() array indices after YZ symmetry transformation. SHARED
  type(var_arr_1d), allocatable :: iyz(:)
  !> @param ixyz(0:dim_npix_unique-1)%a(0:npix_hp_unique(i)-1) Scaspe() array indices after XYZ symmetry transformation. SHARED
  type(var_arr_1d), allocatable :: ixyz(:)
  !> @param chunk_size Chunk size for the OPENMP loops. Set in set_chunk_size().
  integer :: chunk_size
  !> @param npix Number of pixels in the scaspe() arrays, including user - defined directions (see file_dir_out())
 ! integer :: npix

  !> @param npix_arr Number of pixels in the scaspe() arrays at each wavelength, including user - defined directions (see file_dir_out())
  integer, allocatable :: npix_arr(:)
  !> @param npix_hp Number of pixels in the scaspe() arrays, NOT including user_defined directions
  integer :: npix_hp

  !> @param tot_npix_arr_local Total number of elements in npix_arr for the wavelength locally stored in the scaspe_arr() array. Used to calculate how many elements have to be stored in scaspe_temp_send().
  integer, allocatable :: tot_npix_arr_local(:)

  !> @param npix_hp_arr Number of pixels in the scaspe() arrays at each wavelength, NOT including user_defined directions
  integer, allocatable :: npix_hp_arr(:)
  !> @param npix_unique List of the single values present in npix_arr() without repetitions.
  integer, allocatable :: npix_unique(:)

  !> @param npix_hp_unique List of the single values present in npix_hp_arr() without repetitions.
  integer, allocatable :: npix_hp_unique(:)

  !> @param kp_unique List of HEALPix kp values corresponding to the values in npix_unique (remember that npix contains the number of HEALPix pixel plus the external observer lines of sight).
  integer, allocatable :: kp_unique(:)
  
  !> @param dim_npix_unique Number of elements in npix_unique(). 
  integer :: dim_npix_unique

  !> @param max_npix Maximum value of npix_arr() within the wavelength range used in the current RT algorithm. 
  integer :: max_npix
  
  !> @param kp_sca_arr Values of kp parameter used to determine the values of npix_arr() and npix_hp_arr() at each wavelength.
  integer, allocatable :: kp_sca_arr(:)
  
  !> @param ik_sca_arr For each wavelength, it gives the right subscript to use when operating with arrays with dimension dim_npix_unique(). For example, it is used to know which element of kp_unique() corresponds to a certain element in kp_sca_arr(). If an element has value -1, it means that  scattering is considered isotropic for the corresponding wavelength. 
  integer, allocatable :: ik_sca_arr(:)

  !> @param i2d External loop index in rt_loop_2D().
  integer :: i2d
  
  ! -- RT_TYPE ID VALUES -- ! 
  !> @param rtt_start rt_type() ID for the beginning of the RT run.
  integer, parameter :: rtt_start = 0 
  !> @param rtt_precalc_cell rt_type() ID for the precalculation of the radiation field due to the emitting cells. 
  integer, parameter :: rtt_precalc_cell = 1 
  !> @param rtt_precalc_src rt_type() ID for the precalculation of the radiation field due to the point sources.  
  integer, parameter :: rtt_precalc_src = 2 
  !> @param rtt_output_part1 rt_type() ID for the output of the results of the precalculation.
  integer, parameter :: rtt_output_part1 = 3 
  !> @param rtt_prep_part2 rt_type() ID  for the preparation of the direct light processing. 
  integer, parameter ::  rtt_prep_part2 = 4 
  !> @param rtt_dir_cell rt_type() ID for the direct light calculation for emitting cells.
  integer, parameter ::  rtt_dir_cell = 5 
  !> @param rtt_i_obs_dir_cell rt_type() ID for the outgoing radiation specific intensity calculation for emitting cells. 
  integer, parameter :: rtt_i_obs_dir_cell = 6 
  !> @param rtt_dir_src rt_type() ID for the direct light calculation for point sources.
  integer, parameter ::  rtt_dir_src = 7 
  !> @param rtt_i_obs_dir_src rt_type() ID for the outgoing radiation specific intensity calculation for point sources. 
  integer, parameter ::  rtt_i_obs_dir_src = 8 
  !> @param rtt_output_part2 rt_type() ID for the output of the direct light processing. 
  integer, parameter :: rtt_output_part2 = 9 
   !> @param rtt_scatt rt_type() ID for the scattered radiation calculation. 
   integer, parameter :: rtt_scatt = 10 
   !> @param rtt_i_obs rt_type() ID for the outgoing radiation specific intensity calculation for scattered radiation sources. 
   integer, parameter ::  rtt_i_obs = 11
   !> @param rtt_read_i_obs_part2 rt_type() ID for the reading of i_obs() array (direct light only) from file. 
   integer, parameter :: rtt_read_i_obs_part2 = 12
   !> @param rtt_read_i_obs rt_type() ID for the reading of i_obs() array (direct light and scattered light included) from file. 
   integer, parameter :: rtt_read_i_obs = 13
   !> @param rtt_read_ufield rt_type() ID for the reading of u_field() array from file. 
   integer, parameter :: rtt_read_ufield = 14 
   !> @param rtt_read_scaspe_tot rt_type() ID for the reading of the scaspe_tot() array from file.
   integer, parameter :: rtt_read_scaspe_tot = 15
   !> @param rtt_grid_init_stars rt_type() ID for stellar emission RT grid initilialization. 
   integer, parameter :: rtt_grid_init_stars = 16
   !> @param rtt_grid_init_dust rt_type() ID for dust emission RT grid initilialization. 
   integer, parameter :: rtt_grid_init_dust = 17
   !> @param rtt_grid_init_projection ID for projection algorithm grid initilialization. 
   integer, parameter :: rtt_grid_init_projection = 18

   
! -- RT_ALGORITHM ID VALUES -- ! 
   !> @param rta_main rt_algorithm_ID() ID value for the standard RT calculation.
   integer, parameter :: rta_main = 0
   !> @param rta_i_obs rt_algorithm_ID() ID value for the outgoing radiation specific intensity RT calculation for the stellar emission. 
   integer, parameter :: rta_i_obs = 1
   !> @param rta_2D rt_algorithm_ID() ID value for the 2D RT calculation (still 3D but it exploits the 2D symmetry of the input model)
   integer, parameter :: rta_2D = 2
   !> @param rta_sed rt_algorithm_ID() ID value for the SED calculation.
   integer, parameter :: rta_sed = 3
   !> @param rta_dust rt_algorithm_ID() ID value for the dust RT calculation.
   integer, parameter :: rta_dust = 4
   !> @param rta_dust_2D rt_algorithm_ID() ID value for the dust RT calculation in 2D mode.
   integer, parameter :: rta_dust2D = 5
   !> @param rta_i_obs_dust rt_algorithm_ID() ID value for the outgoing radiation specific intensity RT calculation for the dust emission. 
   integer, parameter :: rta_i_obs_dust = 6
   !> @param rta_sed_dust rt_algorithm_ID() ID value for the SED calculation for the dust emission.
   integer, parameter :: rta_sed_dust = 7
   !> @param rta_projection rt_algorithm_ID() ID value for projection of physical quantities as optical depth, intrinsic stellar emission.
   integer, parameter :: rta_projection = 8
   
   ! -- RAY_STATUS ID VALUES -- !
   !> @param ras_go_high ray_status ID value for rays that have to be split.
   integer, parameter :: ras_go_high = 0
   !> @param ras_gone ray_status ID value for rays that can be eliminated.
   integer, parameter :: ras_gone = 1
   !> @param ras_go_low ray_status ID value for rays that can be potentially merged.
   integer, parameter :: ras_go_low = 2
   !> @param ras_first_launch ray_status ID value for rays that have just been cast from the corresponding radiation source.
   integer, parameter :: ras_first_launch = 3
   !> @param ras_re_launched ray_status ID value for rays that are child/parent rays and that start their propagation where the corresponding parent/child rays have been blocked
   integer, parameter :: ras_re_launched = 4 
   
   !> @param print_scaspe_tot Logical equal to TRUE if scaspe_tot() has to be printed at the end of the scattering iterations. The scaspe_tot() array can then be used in RT algorithm "i_obs" (see rt_algorithm()). 
   logical :: print_scaspe_tot

   !> @param print_psel_av Equal to TRUE if psel_av_arr() array has to be in the output. 
  logical :: print_psel_av

   !> @param sequential_scattering True if, during each scattering iterations thescattered luminosity is propagated for the same scattering order. That is, the contribution to the scattered radiation, calculated after each ray-cell intersection, are stored in scaspe_arr but the assignment of scaspe_temp is done using scaspe_prev_arr which contains the scaspe_arr calculated in the previous iteration (before the current iteration starts, scaspe_arr is initialised).
   logical :: sequential_scattering

   !> param tol_p_src Tolerance parameter to find host cell of a point source
   real(kind=real64), parameter :: tol_p = 1E-6 

   !!! INDICES ARRAYS

   !> @param iq_a Array whose elements are equal to .TRUE. if the corresponding ray_intensity element is higher than zero. This array does not vary its size during the calculation. PRIVATE
   logical, allocatable :: iq_a(:)

   !> @param iq Subscript list of the ray_intensity element higher than zero. This array varies its size during the calculation. PRIVATE
   integer, allocatable :: iq(:) 

  !> @param lnum_a Number of ray_intensity elements higher than zero
   integer :: lnum_a

  !> @param lnum_a_old Number of ray_intensity elements higher than zero in the previous ray-cell intersection step
   integer :: lnum_a_old 

   !>@param lnum_node Number of wavelengths considered in the scaspe and i_obs arrays stored within the node memory. 
   integer :: lnum_node

   !>@param lnum_maps Number of wavelengths for which the surface brightness maps have to be calculated. 
   integer :: lnum_maps

   !>@param lnum_node_maps Number of wavelengths for which the surface brightness maps have to be calculated and that are stored in the local i_obs arrays. 
   integer :: lnum_node_maps

   !>@param lnum_node_arr Array of number of wavelengths considered in the scaspe and i_obs arrays in each node memory. SHARED 
   integer, allocatable :: lnum_node_arr(:)

   !> @param iq_sca_node Arrays whose elements are .TRUE. if corresponding wavelength in lambda_arr() is considered in the local scaspe and i_obs arrays. SHARED
   logical, allocatable :: iq_sca_node(:)

   !> @param iq_sca_id Array of the subscripts corresponding to the wavelengths stored in the local scaspe arrays (SHARED). 
   integer, allocatable :: iq_sca_id(:) 

   !> @param iq_maps_id Array of the subscripts corresponding to the wavelength index of the locally stored i_obs_arr() arrays for which the surface brightness maps have to be calculated (SHARED). 
   integer, allocatable :: iq_maps_id(:) 

   !> @param im_lambda_arr Array of the subscripts corresponding to the MPI process which hosts the corresponding wavelengths in the local scaspe arrays (SHARED). 
   integer, allocatable :: im_lambda_arr(:) 
   
   !> @param ray_intensity_dep Array containing the elements of ray_intensity higher than zero. It is input to the subroutine deposit (you cannot input ray_intensity(iq) into a subroutine. The array elements are not updated.) PRIVATE
   real(KIND=real64), allocatable :: ray_intensity_dep(:)

   !> @param no_communications TRUE if no communication mode has to be used. This means that the scaspe and i_obs arrays are not distributed among the different MPI processes. Instead, copies of the entire arrays are created for each MPI process. This mode requires much more memory for each node and thus it is suggested only when very few wavelengths are used or the 3D grid resolution is very small. The advantage of using this method is that no communication is performed during the RT calculation. At the end of each stage, the arrays are reduced as for the u_final() arrays. This mode requires sequential_scattering() = .TRUE. because otherwise there is a complex race condition to avoid in assign_scaspe_temp_arr().
  logical :: no_communications

   !------------------------------------------------
   ! EN_SCA storing arrays data type and variables
   !------------------------------------------------

   !> Contains the single values of en_sca which have to be sent to other nodes
   type, bind(C) :: list_en_sca 
      !> @param cc Intersected cell ID number
      integer :: cc
      !> @param nside HEALPix nside parameter associated with the ray 
      integer :: nside 
      !> @param ipix HEALPix direction ID number 
      integer :: ipix 
      !> @param il Wavelength ID number 
      integer :: il
      !> @param en_sca Scattered luminosity
      real(kind=real64) :: en_sca
   end type list_en_sca

   !> @param en_sca_list_received Contains the scattering parameters received by the local MPI process SHARED 
   type (list_en_sca), allocatable :: en_sca_list_received(:)

   !> @param en_sca_id_thread Contains the ID of the thread that has to process the corresponding element in en_sca_list_received SHARED 
   integer, allocatable :: en_sca_id_thread(:)

   !> @param en_sca_list List of scattered luminosities and associated variables that have to be sent to other nodes  PRIVATE
   type (list_en_sca), allocatable :: en_sca_list(:)

   !> @param en_sca_list_all Contains all the en_sca_list() arrays, so the MASTER thread can send this list to the other MPI processes. Note this array is OPEN MP SHARED not PRIVATE. Before the elements are transfered to other arrays, en_sca_list_all() is sorted in order of receiving MPI process.    
   type (list_en_sca), allocatable :: en_sca_list_all(:)

   !> @param temp_en_sca_list Copy of en_sca_list() used when sorting array  PRIVA
   type (list_en_sca), allocatable :: temp_en_sca_list(:)

   !> @param en_sca_arrtype MPI data type ID for en_sca_list_arr() elements SHARED
   integer :: en_sca_arrtype

   !> @param ind_en_sca_list Indeces of the MPI process to which en_sca_list()  elements are transfered PRIVATE
   integer, allocatable :: ind_en_sca_list(:) 

   !> @param count_en_sca Counter of elements stored in en_sca_list for the single OPENMP thread (PRIVATE). 
   integer :: count_en_sca

   !> @param count_en_sca_arr Array of counters of elements stored in en_sca_list for all OPENMP threads (SHARED). 
   integer, allocatable :: count_en_sca_arr(:)

   !> @param count_en_sca_tot Total number of elements to be input in en_sca_list_all() (SHARED).
   integer :: count_en_sca_tot

   !> @param i0_count_en_sca(0: nproc()* np_mpi()-1) Positions of the first element to be placed within en_sca_list_all by each OPENMP thread for each MPI process (SHARED). 
   integer, allocatable :: i0_count_en_sca(:)

   !> @param el_out_arr(0: nproc()-1, 0: np_mpi() -1) Array containing the number of en_sca_list element to be sent to each MPI process for each OPENMP thread  (SHARED). 
   integer, allocatable :: el_out_arr(:,:)

   !> @param i0_el_out_arr(0: nproc()-1, 0: np_mpi() -1) Array containing the first position of the blocks of en_sca_list elements to be sent to each MPI process for each OPENMP thread  (SHARED). 
   integer, allocatable :: i0_el_out_arr(:,:)

   !> @param el_out_mpi Array containing the number of en_sca_list element to be sent to each MPI process (SHARED)
   integer, allocatable :: el_out_mpi(:)

   !> @param el_in_mpi Array containing the number of en_sca_list element to be received by MPI process (SHARED)
   integer, allocatable :: el_in_mpi(:)

   !> @param i0_el_out_mpi Array containing the first position of the blocks of en_sca_list elements to be sent to each MPI process (SHARED).
   integer, allocatable :: i0_el_out_mpi(:)

   !> @param i0_el_in_mpi Array containing the first position of the blocks of en_sca_list elements to be received by each MPI process (SHARED).
   integer, allocatable :: i0_el_in_mpi(:)

   !> @param en_sca_confirm TRUE when the en_sca_list elements from the corresponding MPI process have been received (SHARED)
   logical, allocatable :: en_sca_confirm(:)
   
   !> @param size_en_sca_list Size of en_sca_list(). 
   integer :: size_en_sca_list 

   !> @param n_el_in_tot Total number of en_sca_list_received elements 
   integer :: n_el_in_tot

   !> @param ThreadID OPENMP thread ID number 
   integer :: ThreadID

   !> @param mpi_proc_completed It is .TRUE. if the corresponding MPI process does not have any other en_sca element to send out. If all its threads are in the waiting part at the end of the source loop in rt_loop or rt_loop_2d, it is used to exit the loop when all elements are .TRUE. SHARED
   logical, allocatable :: mpi_proc_completed(:)

   !> @param wait_thread_arr It is .TRUE. if the corresponding thread has entered the waiting region at the end of the source loop in rt_loop or rt_loop_2d. SHARED
   logical, allocatable :: wait_thread_arr(:)

   !> @param handling_mpi_thread_arr It is .TRUE. if the corresponding thread has called the handle_mpi_transfer routine SHARED
   logical, allocatable :: handling_mpi_thread_arr(:)

   !> @param count_transfer Counter of MPI transfers
   integer :: count_transfer

   !> @param count_transfer_start Starting value for count_transfer(). 
   integer, parameter :: count_transfer_start=100

   !> @param num_processed_cells Counter of processed cells
   integer :: num_processed_cells

   !-----------------------------------------------
   ! Variables for scaspe_temp transfer 
   !----------------------------------------------

   !> @param cc_list_all(0: np_mpi()-1,0: nproc()-1) ID numbers of the cells whose scaspe_arr() values have to be transmitted between the MPI processes (SHARED)
   integer, allocatable :: cc_list_all(:,:) 
   
   !> @param cc_list_local(0: nproc()-1) ID numbers of the cells whose scaspe_arr() values have to be transmitted to the local MPI processes (SHARED)
   integer, allocatable :: cc_list_local(:) 

   !> @param scaspe_temp_recv(0:lnum-1)%a(0: npix_hp()+ tot_ndir_scaspe()-1,0: nproc()*num_scaspe_pass()-1) Scaspe_arr() elements received by the local MPI process (SHARED)
   type(var_arr_2d), allocatable :: scaspe_temp_recv(:)

  !> @param scaspe_temp_send(0: sum(npix_arr(i))* (nproc()* np_mpi()*num_scaspe_pass)-1) with the sum over npix_arr only for i corresponding to the wavelengths locally stored in the scaspe_arr() array. Scaspe_arr() elements sent by the local MPI process (SHARED).  
   real(kind=real64), allocatable :: scaspe_temp_send(:)

   !> @param scaspe_temp_arr(0:lnum-1)%a(0: npix_hp_arr()+ tot_ndir_scaspe()-1) Extracted part of the scaspe_arr() array used by a single thread.  (PRIVATE)
   type(var_arr_1d), allocatable :: scaspe_temp_arr(:)

   !> @param scaspe_temp_arr_big(0:lnum-1)%a(0: npix_hp_arr()+ tot_ndir_scaspe()-1,0:num_scaspe_pass-1) Storage of the extracted parts of the scaspe_arr() array used by a single thread. This is used so the scaspe values are transfered less often (PRIVATE)
   type(var_arr_2d), allocatable :: scaspe_temp_arr_big(:)

   !> @param iscaspe_big Counter of used scaspe_temp_arr_big() elements 
   integer :: iscaspe_big 

   !> @param transfer_type_all(0: nproc()-1, 0: np_mpi()-1) Array containing the transfer_type value of each OPENMP thread and MPI process SHARED
   integer*1, allocatable :: transfer_type_all(:,:)

   !> @param transfer_type_local(0: nproc()-1) Array containing the transfer_type values for the local OPENMP threads. SHARED
   integer*1, allocatable :: transfer_type_local(:)

   !> @param  transfer_type_tot_local Sum of  transfer_type_local()
   integer*1 :: transfer_type_tot_local

   !> @param transfer_type_tot_all Sum of transfer_type_tot_local()
   integer :: transfer_type_tot_all

   !> @param num_scaspe_pass Number of cells whose scaspe values are transfered to other nodes. 
   integer :: num_scaspe_pass

   !---------------------------------
   ! MPI TRANSFER PARAMETERS
   !---------------------------------

   !> @param MPI_TRANSFER_SCASPE ID used in handle_mpi_transfers() to distinguish call requiring transfer of scaspe arrays
   integer*1, parameter :: MPI_TRANSFER_SCASPE = 0

   !> @param MPI_TRANSFER_EN_SCA ID used in handle_mpi_transfers() to distinguish call requiring process
   integer*1, parameter :: MPI_TRANSFER_EN_SCA = 1 

   !----------------------------------
   !SED parameters and arrays
   !----------------------------------

   !> @param dist_obs Observer distance [pc] used when calculating integrated SED.
   real(kind=real64) :: dist_obs  

   !> @param ind_i_obs List of indeces of the wavelengths of the i_obs arrays to be printed. By default ind_i_obs contains all the wavelength indeces. If indeces are specified in the input file, only the i_obs for the corresponding wavelengths are printed at the end of the calculation.  
   integer, allocatable :: ind_i_obs(:)

   !> @param ind_out_maps List of indeces of the wavelengths for which the surface brightness maps have to be printed. By default ind_out_maps contains all the wavelength indeces. If indeces are specified in the input file, only the maps for the corresponding wavelengths are printed at the end of the calculation.  
   integer, allocatable :: ind_out_maps(:)

   !> @param max_n_input Maximum number of input elements for allocatable input arrays as ind_i_obs() and ind_out_maps().
   integer, parameter :: max_n_input = 1000

   !> @param sed_arr(0: lnum_tot()-1, 0: tot_ndir() -1) Contains total emission SED including all components and for all directions (stellar or dust emission depending on the algorithm).
   real(kind=real64), allocatable :: sed_arr(:,:)
   !> @param sed_arr_dir(0: lnum_tot()-1, 0: tot_ndir() -1) Contains emission SED of the direct light (stellar or dust emission depending on the algorithm).
   real(kind=real64), allocatable :: sed_arr_dir(:,:)

   !> @param test_run If set TRUE it will not run the rt_loop routines. Useful to check that the all the subroutines are running fine and the output arrays are correctly printed. However, the results will not be correct and therefore should be deleted afterwards. 
  logical :: test_run

  ! ----------------------------------
  ! WALL Parameters
  ! ----------------------------------

  !> @param x_wall_on TRUE if wall perpendicular to X direction is set
  logical :: x_wall_on

  !> @param y_wall_on TRUE  if wall perpendicular to Y direction is set
  logical :: y_wall_on

  !> @param z_wall_on TRUE  if wall perpendicular to Z direction is set
  logical :: z_wall_on

  !> @param x_wall_coord Positions of the X-perpendicular walls in relative coordinates (0 = -modelsize/2., 1 = modelsize/2).
  real(kind=real64) :: x_wall_coord(2)

  !> @param y_wall_coord Positions of the Y-perpendicular walls in relative coordinates (0 = -modelsize/2., 1 = modelsize/2).
  real(kind=real64) :: y_wall_coord(2)

  !> @param z_wall_coord Positions of the Z-perpendicular walls in relative coordinates (0 = -modelsize/2., 1 = modelsize/2).
  real(kind=real64) :: z_wall_coord(2)

  !> @param only_direct_rt TRUE if only direct light has to be processed.
  logical :: only_direct_rt

  ! --------------------------------------
  ! projection algorithm parameters 
  ! -------------------------------------
  !> @param param_to_project Parameter to be projected in the 'projection' rt_algorithm(). Choices are: 'stellar_emission' and 'optical_depth'. 
  character(LEN=lcar) :: param_to_project

   
  !$OMP THREADPRIVATE(src_ccoord,nside,ccindd_nc, ccindd, ffn_arr_mw, ffn_arr, ads_arr,clvl_old, psel_av, ipsel_av, i_obs_temp, &
   !$OMP lum_lost_temp, pabs_arr,vec_mod, iq_a, iq, lnum_a,lnum_a_old, ray_intensity_dep, en_sca_list, count_en_sca, ThreadID, scaspe_temp_arr, scaspe_temp_arr_big, iscaspe_big, temp_en_sca_list, ind_en_sca_list, nside_sca)


CONTAINS

!> Finds the host cells for the stellar point sources.
!> It sets the values for src_cell() and cell_src().
subroutine prepare_p_src()
  
  integer :: i,j,in_cube(3)
  real(kind=real64) :: rel_vec(3),cellsize 

  if (tot_p_src == 0) return ! no point sources

  if (main_prc) print *, 'preparing src_cell and cell_src arrays....'

  do i=0, tot_p_src -1 

     do j=0, tot_ncell-1
        if (cchild(j) /= -1) cycle 
        in_cube=0
        rel_vec=ccoord_p_src(:,i)-ccoord(:,j)
        cellsize=csize(j) 
        
        where (abs(rel_vec) < (1.+tol_p)*cellsize/2.)  !!! is the factor 1E-6 OK ? 
           in_cube =1 
        end where

        if (sum(in_cube) == 3) then 
           src_cell(j)=i
           cell_src(i)=j

           call fix_ccoord_p_src(ccoord_p_src(:,i), rel_vec,cellsize)
           
           exit
        end if

     end do

     if (sum(in_cube) /= 3) then 
        print *, 'STOP: host cell not found. Something wrong here'
        call stop_prc
     endif

  end do

  call print_done
  
end subroutine prepare_p_src

!> Shifts slightly the source position within the host cell in case it is located on the border of the host cell.
subroutine fix_ccoord_p_src(pos, rel, cellsize)
  real(kind=real64) :: pos(3), rel(3), cellsize
  integer :: i

  do i = 1, 3
     
     if (abs(rel(i)) >= cellsize/2.) then

        if (rel(i) > 0) then ! source on the "right" 

           pos(i) = pos(i) - tol_p*cellsize

        else

           pos(i) = pos(i) + tol_p*cellsize

        endif

     endif

  end do
  
end subroutine fix_ccoord_p_src
  

!> Calculates the total stellar luminosity within the RT model and store it in lumcell(). It also sets the arrays tot_rad_en() and tot_rad_en_or().
subroutine calc_total_luminosity
IMPLICIT NONE
integer :: i,j,i0,i1

if (main_prc) print *, 'calculating total stellar luminosities....'

call set_i_opacity_arrays(i0,i1)

! define luminosity for each cell  
if (.not. allocated(lumcell)) allocate(lumcell(0 : lnum -1, 0:tot_ncell-1))
lumcell=0

do i=0, tot_ncell-1
 
   if (cchild(i) /= -1) cycle
   
   lumcell(:,i)=dens_stars_arr(:,i)*(csize(i)**3)

end do

! set to zero luminosity of grouped cells so they will not be processed in rt_loop (see reduce_grid_res).

!call nullify_lum_group_cells 

! calculate total luminosity 

if (.not. allocated(tot_rad_en)) allocate(tot_rad_en(0:lnum-1), tot_rad_en_or(0:lnum-1))

do i = 0, lnum-1 

tot_rad_en(i)=sum(lumcell(i,:))

if (tot_p_src > 0) then
   tot_rad_en(i)=tot_rad_en(i)+sum(lum_p_src_arr(i, :))
endif

if (main_prc) print *, 'L(',lambda_arr(i+i0),') = ',tot_rad_en(i)

end do 



!!! Initialize lost luminosity
if (.not. allocated(lum_lost)) allocate(lum_lost(0:lnum-1))
if (iterations_dustem <= 1) lum_lost=0  ! the <= 1 includes both the stellar emission case and the first iteration for the dust emission case 

if (.not. allocated(lum_lost_prev)) allocate(lum_lost_prev(0:lnum-1))
if (iterations_dustem <= 1) lum_lost_prev = 0

!!! Store original total luminosity (at the beginning of scattering iterations, tot_rad_en is redefined). In case of dust RT, it adds the new luminosities added in each dust heating iterations 
if (iterations_dustem <= 1) tot_rad_en_or = 0
tot_rad_en_or=tot_rad_en_or+tot_rad_en 

call print_done 

end subroutine calc_total_luminosity

!> Calculates tot_rad_en() and lumcell() at the beginning of each scattering iteration.
subroutine calc_total_luminosity_sca

  integer :: i,j,k,ierr,tot_el
  integer :: i0, i1 
  real(kind=real64), allocatable :: out_arr(:),out_arr_2d(:,:)
  real(kind=real64) :: tot_lumcell(0:lnum-1)

  ! set starting wavelength index 
  call set_i_opacity_arrays(i0,i1)

  select case(rt_type)
     case(rtt_scatt)
        if (iterations == 1) then 

           tot_rad_en = 0  ! calculate tot_rad_en at first iteration
           do i = 0, lnum -1
              if (no_communications .and. .not. main_prc) exit ! using this, the main MPI process is the only one calculating tot_rad_en in the no communications mode.
              if (iq_sca_node(i)) then 
                 k = (i - id_mpi)/np_mpi  
                 tot_rad_en(i) =sum(scaspe_arr(k)%a(0:npix_hp_arr(i+i0)-1,:))
              endif
           enddo
           allocate(out_arr(0:lnum-1))
           out_arr = 0
           call mpi_reduce(tot_rad_en, out_arr, lnum, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr)
           if (main_prc) then 
              tot_rad_en = out_arr
           endif
           deallocate(out_arr)
   
           call mpi_Bcast(tot_rad_en,lnum, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,ierr)  
           
        endif

        ! calculate total cells escaping luminosity 
        lumcell = 0
        do i=0,tot_ncell-1
           if (no_communications .and. .not. main_prc) exit 
           if (cchild(i) /= -1) cycle
           do j = 0, lnum-1
              if (iq_sca_node(j)) then
                 k = (j - id_mpi)/np_mpi
                 lumcell(j,i)=sum(scaspe_arr(k)%a(0:npix_hp_arr(j+i0)-1,i))
              endif
           enddo
        end do
        
        allocate(out_arr_2d(0 : lnum -1, 0:tot_ncell-1))
        out_arr_2d = 0 
        tot_el = lnum*tot_ncell
        call mpi_reduce(lumcell, out_arr_2d, tot_el, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr)
        if (main_prc) then 
           lumcell = out_arr_2d
        endif
        deallocate(out_arr_2d)
        
        call mpi_Bcast(lumcell,tot_el, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,ierr)  

        !call nullify_lum_group_cells 
        
        ! check convergence criteria 
    
        do i = 0, lnum-1
           tot_lumcell(i)=sum(lumcell(i,:))
        enddo

        ! flag convergence reached
        if (sum(tot_rad_en) > 0) then 
           if (main_prc) print *, 'Maximum escaping energy/total scatt energy=', maxval(tot_lumcell/tot_rad_en) 
           cnflag=(all(tot_lumcell/tot_rad_en < conv_en_lim))
        else if (sum(tot_rad_en) == 0) then 
           cnflag = .TRUE.
        else 
           if (main_prc) print *, 'weird value for tot_rad_en :', tot_rad_en
           call stop_prc
        endif

        ! check if maximum iteration has been reached
        if (limit_scattering_iterations .and. iterations > max_sca_iterations) then
           if (main_prc) print *, 'Maximum number of scattering iteration reached!'
           cnflag = .TRUE.
        endif

        ! set scaspe_arr to 0 in case of sequential scattering algorithm. Scaspe_prev_arr contains the values calculated in the previous iteration. Before that it adds the scaspe_arr of the previous iteration to scaspe_tot_arr
        if (sequential_scattering) then
           if (iterations > 1) then
              do i = 0, lnum_node -1
                 scaspe_tot_arr(i)%a = scaspe_tot_arr(i)%a + scaspe_arr(i)%a
              enddo
           endif
           do i = 0, lnum_node -1 
              scaspe_arr(i)%a = 0
           enddo
        endif
           
        call mpi_barrier(MPI_COMM_WORLD, ierr)

     case(rtt_i_obs)

        lumcell = 0
        do i=0,tot_ncell-1
           if (cchild(i) /= -1) cycle
           do j = 0, lnum-1
              if (iq_sca_node(j)) then
                 k = (j - id_mpi)/np_mpi
                 lumcell(iq_sca_id(k),i)=sum(scaspe_tot_arr(k)%a(:,i))  ! in the iobs loop, this is just needed to check that there is scattered luminosity in a cell 
              endif
           enddo
        end do
        
     end select

end subroutine calc_total_luminosity_sca

!!!!!> Assigns the right value of total luminosity and dust density to the parent cells. Needed to use lower resolution grids (see reduce_grid_res() ). Note this subroutine has to be used outside OPENMP RT loops.  
!!$subroutine assign_parent_dens_arr 
!!$  integer :: i, j,k, ib,l,cc 
!!$  integer :: lvl_i, lvl_cc
!!$  real(kind=real64) :: tot_lum_parent, tot_lum_child, tot_dust_parent, tot_dust_child 
!!$
!!$  ! initialize parent cells
!!$
!!$  do i = 0, tot_ncell-1 
!!$
!!$     if (cchild(i) == -1) cycle 
!!$
!!$     dens_stars_arr(:,i) = 0      
!!$     dens_arr(:,i) = 0      
!!$
!!$  end do
!!$
!!$  ! Add cchild value to parent cell
!!$
!!$  allocate(ccindd(3,max_lvl))
!!$ 
!!$  do i = 0, tot_ncell-1
!!$
!!$     if (cchild(i) /= -1) cycle
!!$
!!$     ccindd = 0
!!$     lvl_i = lvl(i)
!!$
!!$     call cindex_to_ccindd(i,lvl_i,ccindd)
!!$
!!$     l=1 
!!$     ib=1
!!$     do j=1,max_lvl
!!$        if (j > 1) ib=2 
!!$        k=(ccindd(3,j)*base(ib)+ccindd(2,j))*base(ib)+ccindd(1,j)+1
!!$        l=l+k-1
!!$  
!!$        if (lvl(l) == lvl_i-1) then          
!!$           cc=l
!!$           lvl_cc = lvl(cc)
!!$           exit
!!$        else
!!$           l=cchild(l)       
!!$        endif
!!$     end do
!!$
!!$     dens_stars_arr(:,cc) = dens_stars_arr(:,cc) + dens_stars_arr(:,i)*(csize_arr(lvl_i)/csize_arr(lvl_cc))**3
!!$     dens_arr(:,cc) = dens_arr(:,cc) + dens_arr(:,i)*(csize_arr(lvl_i)/csize_arr(lvl_cc))**3
!!$     
!!$  end do
!!$
!!$  deallocate(ccindd)
!!$
!!$  ! Check parent cells contain same total luminosity as child cells
!!$  do i = 0, lnum-1
!!$
!!$     tot_lum_parent = sum(dens_stars_arr(i,:)*csize(:)**3, mask = cchild /= -1)
!!$     tot_lum_child = sum(dens_stars_arr(i,:)*csize(:)**3, mask = cchild == -1)
!!$     if (abs(tot_lum_parent-tot_lum_child)/abs(tot_lum_parent) > 1E-8) then
!!$        if (main_prc) print *, 'STOP(assign_parent_dens_arr): Luminosity in parent cells do not match that in child cells! Something is wrong...'
!!$        STOP
!!$     endif
!!$
!!$     tot_dust_parent = sum(dens_arr(i,:)*csize(:)**3, mask = cchild /= -1)
!!$     tot_dust_child = sum(dens_arr(i,:)*csize(:)**3, mask = cchild == -1)
!!$     if (abs(tot_dust_parent-tot_dust_child)/abs(tot_dust_parent) > 1E-8) then
!!$        if (main_prc) print *, 'STOP(assign_parent_dens_arr): Total dust in parent cells do not match that in child cells! Something is wrong...'
!!$        STOP
!!$     endif
!!$
!!$  end do
!!$
!!$end subroutine assign_parent_dens_arr

!!!!!> Groups the smallest cells if their optical depth, at the shortest wavelength, is very small (that is, less than tau_cell_max()). In this case, it assigns a value of -1 to the cchild() element of their parent cell. The original values of cchild() are stored in cchild_or().
!!!!> \todo This feature is still experimental. 
!!$subroutine reduce_grid_res
!!$  
!!$  integer :: i,count 
!!$  real(kind=real64) :: tau_cell
!!$  
!!$  tot_spare_cells = 0
!!$
!!$  if (rt_algorithm_ID /= rta_dust .and. rt_algorithm_ID /= rta_dust2D) return ! grid resolution reduced only for dust RT calculation 
!!$  
!!$  if (max_lvl > 3) then
!!$     
!!$     if (main_prc) print *, 'reduce grid resolution if possible....'  
!!$  
!!$!!! create copy of cchild
!!$     if (.not.allocated(cchild_or)) then
!!$        allocate(cchild_or(0:tot_ncell-1))
!!$        cchild_or=cchild
!!$     endif
!!$
!!$!!! group optically thin cells with lvl = max_lvl -1 
!!$     count =0 
!!$     do i=0, tot_ncell -1
!!$        
!!$        if ((lvl(i) == max_lvl-1).and.(dens_arr(0,i) > 0.).and.(cchild(i) /= -1)) then   
!!$           tau_cell=real(dens_arr(0,i)*csize(i)) ! 0 has to be the index of the shortest wavelength. 
!!$           
!!$           if (tau_cell < tau_cell_max) then
!!$              count=count +1
!!$              cchild(i) = -1             
!!$           endif
!!$           
!!$        endif
!!$     enddo
!!$     
!!$     tot_spare_cells=count*base(2)**3-count
!!$  
!!$     if (main_prc) print *, 'number of grouped cells after reducing grid resolution=', tot_spare_cells
!!$     
!!$     call print_done
!!$     
!!$  endif
!!$  
!!$end subroutine reduce_grid_res

!!!!!> Sets to zero the luminosity of grouped cells so they will not be processed in rt_loop (see reduce_grid_res()).
!!$subroutine nullify_lum_group_cells 
!!$  integer :: i, j, ic, nls
!!$
!!$  if (rt_algorithm_ID /= rta_dust .and. rt_algorithm_ID /= rta_dust2D) return ! grid resolution reduced only for dust RT calculation 
!!$
!!$  if (tau_cell_max == 0) return 
!!$  print *, 'STOP(nullify_lum_group_cells) this is not allowed for the moment'
!!$  stop
!!$
!!$  nls=base(2)**3
!!$  
!!$  do i=0, tot_ncell-1
!!$
!!$     if (cchild(i) == -1 .and. cchild_or(i) /= -1) then
!!$        
!!$        do j=0, nls-1
!!$
!!$           ic=cchild_or(i)+j
!!$           
!!$           lumcell(:,ic) = 0 
!!$
!!$        end do
!!$
!!$     endif
!!$
!!$  end do
!!$
!!$end subroutine nullify_lum_group_cells


!!! !> Restores the original grid resolution (if needed).
!!$subroutine restore_grid_original_res
!!$
!!$  integer :: i, j, ic, nls
!!$  integer :: il 
!!$  logical :: allocated_scaspe, allocated_scaspe_tot 
!!$
!!$  if (rt_algorithm_ID /= rta_dust .and. rt_algorithm_ID /= rta_dust2D) return ! grid resolution reduced only for dust RT calculation 
!!$  
!!$  if (max_lvl > 3) then 
!!$  
!!$     if (main_prc) print *, 'restoring original grid resolution if needed'
!!$  
!!$     nls=base(2)**3
!!$     
!!$     allocated_scaspe = allocated(scaspe_arr)
!!$     allocated_scaspe_tot = allocated(scaspe_tot_arr)
!!$  
!!$     do i=0, tot_ncell -1
!!$
!!$        if ((lvl(i) == max_lvl -1).and.(dens_arr(0,i) > 0.).and.(cchild(i) == -1).and.(cchild_or(i) /= -1)) then
!!$        
!!$           do j=0, nls-1
!!$           
!!$              ic=cchild_or(i)+j
!!$              
!!$              u_final_arr(:,ic)=u_final_arr(:,i)
!!$              u_fest_arr(:,ic)=u_fest_arr(:,i)
!!$
!!$              
!!$              if (allocated_scaspe) then
!!$                 do il = 0, lnum_node-1 
!!$                    scaspe_arr(il)%a(:,ic)=scaspe_arr(il)%a(:,i)/nls
!!$                 enddo
!!$              endif
!!$
!!$              if (allocated_scaspe_tot) then
!!$                 do il = 0, lnum_node-1 
!!$                    scaspe_tot_arr(il)%a(:,ic)=scaspe_tot_arr(il)%a(:,i)/nls
!!$                 enddo
!!$              endif
!!$              
!!$           end do
!!$           
!!$           cchild(i)=cchild_or(i)
!!$               
!!$        endif
!!$  
!!$     enddo
!!$
!!$     call print_done
!!$  
!!$  endif
!!$
!!$end subroutine restore_grid_original_res
 
  !> This is the main RT loop over the sources of radiation. It is used in all calculation phases where the radiation field is derived.
 subroutine rt_loop
   
   integer :: i, im,src_id, sid
   real(kind=real64) :: src_lum(0:lnum-1)

   if (test_run) return
   
   call omp_set_num_threads(nproc)  

   !$OMP PARALLEL DEFAULT(NONE), &
   !$OMP PRIVATE(i,im,src_id,src_lum, sid), &
   !$OMP SHARED(tot_ncell,cchild, max_lvl, lvl,rt_type, scaspe_arr,npix_hp,tot_ndir,tot_ndir_scaspe,src_cell, chunk_size,lnum, np_mpi, id_mpi, lnum_node, handling_mpi_thread_arr, count_transfer, num_scaspe_pass, main_prc)

   ThreadID = OMP_GET_THREAD_NUM() 
   call allocate_rt_loop_arrays  
   
   !$OMP DO SCHEDULE(DYNAMIC,chunk_size)

   do im = 0, tot_ncell/np_mpi  !!! loop on cells NOTE: no -1 here    
   
      i = im*np_mpi+id_mpi
      call check_im(im)
      !call put_lock_to_cell(i)
      
      if (i > tot_ncell -1) cycle
       
      !if ( count_transfer > 150) cycle
      
      src_lum=src_lum_value(i)
      src_id=i

      if ((cchild(i) /= -1).or.(sum(src_lum) == 0)) cycle
     
      call count_processed_cells

         ipsel_av=0   !!! counter number of rays 
         psel_av=0    !!! average ray path
            
      call calc_src_ccoord(i)

      if (rt_type == rtt_dir_src ) then 
          sid = tot_ncell+src_cell(i)   ! iobs index          
       else 
          sid = i
       endif

      call cindex_to_ccindd(src_id,lvl(src_id),ccindd_nc)
      
      call assign_scaspe_temp_arr(i)

      !print *, src_id, id_mpi
      
      if (mod(i,1000) == id_mpi) write(*,*) src_id
!!$print *, 'source cell coordinates' 
!!$print *, 'x=', src_ccoord(1)
!!$print *, 'y=', src_ccoord(2)
!!$print *, 'z=', src_ccoord(3)

     call main_dir_loop(src_id,sid,src_lum)

   end do

   !$OMP END DO NOWAIT

   !print *, 'arrivato ', ThreadID, id_mpi
   call wait_end_thread_loop

   call deallocate_rt_loop_arrays

   !$OMP END PARALLEL

   call print_done

 end subroutine rt_loop

 !> Contains the main RT loop over the sources of radiation for the RT 2D mode. The loop is perfomed twice because the first instance processes only the cells with x,y,z > 0 while the second processes those with x =0 OR y=0 OR z=0. Repeating the loops is necessary because only the off-axis cells provide contributions toe the radiation field that can be easily symmetrised.
 subroutine rt_loop_2D
   
   integer :: i, im, src_id, sid
   real(kind=real64) :: src_lum(0:lnum-1)
   integer :: ierr 

   if (test_run) return

   call omp_set_num_threads(nproc) 

   do i2d=0,1 

   !$OMP PARALLEL DEFAULT(NONE), PRIVATE(i,im,src_id,src_lum, sid), &
   !$OMP SHARED(tot_ncell,cchild, max_lvl, lvl,rt_type, scaspe_arr,npix_hp,tot_ndir,tot_ndir_scaspe,src_cell,i2d,ccoord,chunk_size, lnum,np_mpi,id_mpi, handling_mpi_thread_arr, num_scaspe_pass) 

      call allocate_rt_loop_arrays
      ThreadID = OMP_GET_THREAD_NUM()

      !$OMP DO SCHEDULE(DYNAMIC,chunk_size)

      do im=0, tot_ncell/np_mpi  !!! loop on cells NOTE: no -1 here 

         !ThreadID = OMP_GET_THREAD_NUM()      
         i = im*np_mpi+id_mpi
         call check_im(im)
         !call put_lock_to_cell(i)  
         if (i > tot_ncell -1) cycle 
      
         src_lum=src_lum_value(i)
         src_id=i
      
         if ((cchild(i) /= -1).or.(sum(src_lum) == 0)) cycle

         call count_processed_cells

         if (skip_cell_2D_rt_loop(i)) cycle 

         ipsel_av=0   !!! counter number of rays 
         psel_av=0    !!! average ray path
         
         call calc_src_ccoord(i)

         if (rt_type == rtt_dir_src ) then 
            sid = tot_ncell+src_cell(i)   ! iobs index          
         else 
            sid = i
         endif

         call cindex_to_ccindd(src_id,lvl(src_id),ccindd_nc)

         call assign_scaspe_temp_arr(i)         
      
         !print *, src_id, id_mpi
         if (mod(i,1000) == id_mpi) write(*,*) src_id
!!$print *, 'source cell coordinates (first estimate rad-field)' 
!!$print *, 'x=', src_ccoord(1)
!!$print *, 'y=', src_ccoord(2)
!!$print *, 'z=', src_ccoord(3)

         call main_dir_loop(src_id,sid,src_lum)
      
      end do
      
      !$OMP END DO NOWAIT
   
      call wait_end_thread_loop

      call deallocate_rt_loop_arrays

      !$OMP END PARALLEL
      
      if (i2d == 0) then
         call fix_symmetry()     
      endif

   enddo ! end i2d loop

   call print_done 

 end subroutine rt_loop_2D

!> Assigns to scaspe_temp_arr() the values needed during the current OpenMP thread loop chunk. 
 subroutine assign_scaspe_temp_arr(i)
   integer :: i,il
   integer :: iw, ip
   
   if (rt_type == rtt_scatt ) then
      
      if (.not. no_communications) then
         
         if (iscaspe_big >= num_scaspe_pass) then 
            handling_mpi_thread_arr(ThreadID) = .TRUE.          
            call handle_mpi_transfers(MPI_TRANSFER_SCASPE,i)
         endif
         do il = 0, lnum -1 
            scaspe_temp_arr(il)%a = scaspe_temp_arr_big(il)%a(:,iscaspe_big)
         end do

      else
           
         do il = 0, lnum -1 
            scaspe_temp_arr(il)%a = scaspe_prev_arr(il)%a(:,i)  
         end do
      endif
         
   endif

  

 end subroutine assign_scaspe_temp_arr
 

!> Used to skip cell processing within rt_loop_2D. TRUE when the cell i does not have to be processed within 2D loop i2d. 
logical function skip_cell_2D_rt_loop(i)
integer :: i

skip_cell_2D_rt_loop = .FALSE.

if (i2d == 0) then
!!! only positive cells in the first loop. Then symmetrize
   if ((ccoord(1,i) <= 0.).or.(ccoord(2,i) <= 0.).or.(ccoord(3,i) <= 0.)) skip_cell_2D_rt_loop = .TRUE. 
else
!!! only cells on Cartesian axis in the second loop
   if ((ccoord(1,i) /= 0.).and.(ccoord(2,i) /= 0.).and.(ccoord(3,i) /= 0.))skip_cell_2D_rt_loop = .TRUE. 
endif

end function skip_cell_2D_rt_loop


 !> Waits until all threads and all MPI processes have completed the source loop and sent out the remaining en_sca elements.
 subroutine wait_end_thread_loop

   ! no need to wait in precalculations
   if (rt_type == rtt_precalc_cell .or. rt_type == rtt_precalc_src) return 

   wait_thread_arr(ThreadID) = .TRUE.
   do  ! loop until all processes have completed the rt_loop
     
      call handle_mpi_transfers(MPI_TRANSFER_EN_SCA,-1)
      
      if (all(mpi_proc_completed) .and. all(wait_thread_arr).and. all(transfer_type_all == MPI_TRANSFER_EN_SCA)) then
         exit 
      endif
      
   end do

 end subroutine wait_end_thread_loop

!> Allocates arrays used in the rt_loop(), rt_loop_2d() and rt_loop_iobs() subroutines. 
subroutine allocate_rt_loop_arrays
  integer :: i0, i1, il 
  integer :: i,k

  call set_i_opacity_arrays(i0,i1)

  !$OMP MASTER

  size_en_sca_list = 500000!*(lnum-lnum_node)/(lnum-maxval(lnum_node_arr))  ! the different scaling is because nodes which store less scaspe arrays take less time to fill en_sca_list arrays.   
  num_scaspe_pass = chunk_size
  
  count_transfer = count_transfer_start 
  transfer_type_tot_local = 0 
  allocate(count_en_sca_arr(0:nproc-1), i0_count_en_sca(0:nproc*np_mpi-1), mpi_proc_completed(0:np_mpi-1))
  count_en_sca_arr =0 ! initialise here. Important!
  mpi_proc_completed = .FALSE.
  allocate(cc_list_all(0:nproc-1,0:np_mpi-1), cc_list_local(0:nproc-1))
  allocate(scaspe_temp_recv(0:lnum-1))
  do il = 0, lnum -1 
     allocate(scaspe_temp_recv(il)%a(0:npix_arr(il+i0)-1,0:nproc*num_scaspe_pass-1))
  end do
  cc_list_all = -1  ! initialise here. Important! 
  cc_list_local = -1 
  allocate(transfer_type_all(0:nproc-1, 0:np_mpi-1), transfer_type_local(0:nproc-1))
  transfer_type_local = MPI_TRANSFER_EN_SCA ! initialise here. Important! 
  transfer_type_all = MPI_TRANSFER_EN_SCA
  if (lnum_node > 0) then 
     allocate(scaspe_temp_send(0:tot_npix_arr_local(id_mpi)*nproc*np_mpi*num_scaspe_pass-1))
  else 
     allocate(scaspe_temp_send(0:0))
  endif
  
  allocate(en_sca_list_all(0:0))  ! the size of this array is changed later
  allocate(en_sca_list_received(0:0), en_sca_id_thread(0:0)) ! size changed later
  allocate(wait_thread_arr(0:nproc-1))
  wait_thread_arr = .FALSE.
  allocate(handling_mpi_thread_arr(0:nproc-1))
  handling_mpi_thread_arr = .FALSE.
  
  num_processed_cells = 0

  allocate(el_out_arr(0:nproc-1, 0: np_mpi-1), i0_el_out_arr(0:nproc-1, 0: np_mpi-1))
  allocate(el_out_mpi(0: np_mpi-1), i0_el_out_mpi(0: np_mpi-1))
  allocate(el_in_mpi(0: np_mpi-1), i0_el_in_mpi(0: np_mpi-1))
  allocate(en_sca_confirm(0: np_mpi-1))
  en_sca_confirm = .FALSE.  
  
  !$OMP END MASTER 
  !$OMP BARRIER 
  allocate(ccindd_nc(3,max_lvl), ccindd(3,max_lvl))
  if (allocated(scaspe_arr)) then
     allocate(ads_arr(0:dim_npix_unique-1), ffn_arr(0:dim_npix_unique-1), ffn_arr_mw(0:lnum_node-1))
     do i = 0, dim_npix_unique-1
        allocate(ads_arr(i)%a(0:npix_unique(i)-1))
     enddo
     do i = 0, dim_npix_unique-1
        allocate(ffn_arr(i)%a(0:npix_unique(i)-1))
     end do
     do i = 0, lnum -1
        if (iq_sca_node(i)) then
           call set_wavelength_index(i,k)
           allocate(ffn_arr_mw(k)%a(0:npix_arr(i+i0)-1))
        endif
     end do
  endif
  
  allocate(i_obs_temp(0:lnum-1)) 
  allocate(lum_lost_temp(0:lnum-1))
  allocate(iq_a(0:lnum-1))
  allocate(iq(0:lnum-1)) ! note that this array varies its size during the calculation
  allocate(ray_intensity_dep(0:lnum-1))
  call create_en_sca_list
  if (rt_type == rtt_scatt .or. rt_type == rtt_i_obs) then 
     allocate(scaspe_temp_arr(0:lnum-1))
     do i = 0, lnum -1 
        allocate(scaspe_temp_arr(i)%a(0:npix_arr(i+i0)-1))
     end do     
  endif
  if (rt_type == rtt_scatt .and. .not. no_communications) then 
     allocate(scaspe_temp_arr_big(0:lnum-1))
     do i = 0, lnum -1 
        allocate(scaspe_temp_arr_big(i)%a(0:npix_arr(i+i0)-1, 0:num_scaspe_pass-1 ))
     end do
  endif 
  
  !$OMP BARRIER

end subroutine allocate_rt_loop_arrays

!> Deallocates arrays used in the rt_loop(), rt_loop_2d() and rt_loop_iobs() subroutines. 
subroutine deallocate_rt_loop_arrays

  deallocate(ccindd_nc, ccindd)
  if (allocated(scaspe_arr)) deallocate(ads_arr, ffn_arr,ffn_arr_mw)
  deallocate(i_obs_temp)
  deallocate(lum_lost_temp)
  deallocate(iq_a)
  deallocate(iq)
  deallocate(ray_intensity_dep)
  if (allocated(en_sca_list)) deallocate(en_sca_list,temp_en_sca_list, ind_en_sca_list)
  if (rt_type == rtt_scatt .or. rt_type == rtt_i_obs) then 
     deallocate(scaspe_temp_arr)
  endif
  if (rt_type == rtt_scatt .and. .not. no_communications) then 
     deallocate(scaspe_temp_arr_big)
  endif

  !$OMP BARRIER
  !$OMP MASTER

  deallocate(count_en_sca_arr, i0_count_en_sca, mpi_proc_completed)
  deallocate(cc_list_all, cc_list_local, scaspe_temp_recv)
  deallocate(transfer_type_all, transfer_type_local)
  deallocate(scaspe_temp_send)
  deallocate(en_sca_list_all)
  deallocate(en_sca_list_received, en_sca_id_thread)
  deallocate(wait_thread_arr)
  deallocate(handling_mpi_thread_arr)
  deallocate(el_out_arr, i0_el_out_arr)
  deallocate(el_out_mpi, i0_el_out_mpi)
  deallocate(el_in_mpi, i0_el_in_mpi)
  deallocate(en_sca_confirm)

  !$OMP END MASTER 
  !$OMP BARRIER 

  

end subroutine deallocate_rt_loop_arrays



 !> Contains the loop over the HEALPix directions used by rt_loop() and rt_loop_2D(). The minimum HEALPix resolution for the ray angular density corresponds to nside = 4. If necessary, the ray angular density is automatically increased or decreased along the ray path according to the input values for bm_par(), bm_par_max() and bm_par_sca().
 subroutine main_dir_loop(src_id,sid,src_lum)
   
   
   integer :: i, src_id, ipix, clvl,next,sid
   integer :: idir, idir2,nside_ray,isel_ray,i_ext,i0,i1,cc_old
   real(kind=real64) :: src_lum(0:lnum-1)
   real(kind=real64) :: src_lum_ray(0:lnum-1)
   real(kind=real64) :: prev_ray,dplane_ray
   real(kind=real64) :: theta, phi
   integer :: ray_type
   logical :: flag_k
   !real(kind=real64) :: scaspe_temp(0:,0:)
     
     main_dir: do idir=0, npix_main-1  !!! loop on sectors 
   
         dir2:   do idir2=0, 3 !!! each main sector is subdivided in 4 smaller sectors 
               
            nside=nside_min
  
             lum_lost_temp=0 ! initialize lost_temp array

            call create_high_ray_list(idir,idir2,src_lum,src_id,nside/2)
                
            do   !!! loop on resolution until requirement are fullfilled 
               
               call extract_ray_list(nside)
               
              !  print *, nside

               if (.not.no_ray) then  ! needed if no extracted ray
                  next=size(ext_ray_list)
               else 
                  next=0
               endif
                   
               ext_list:    do i_ext=0, next - 1 
                  nside_ray=ext_ray_list(i_ext)%nside
                  theta=ext_ray_list(i_ext)%theta
                  phi=ext_ray_list(i_ext)%phi
                  src_lum_ray=ext_ray_list(i_ext)%src_lum
                  dplane_ray=ext_ray_list(i_ext)%dplane 
                  prev_ray=ext_ray_list(i_ext)%prev
                  isel_ray=ext_ray_list(i_ext)%isel
                  ray_type=ext_ray_list(i_ext)%ray_type 
                  cc_old=ext_ray_list(i_ext)%cc_old   
            
                  if (ray_type == ray_type_high) then 
                     
                     call ang2pix_nest(nside_ray, theta, phi,i0)
                     
                     i0=i0*4  ! this gives the ipix at higher nested levels 
                     i1=i0+3
                     
                  endif

                  if (ray_type == ray_type_reco) then 
                     ! in this case theta and phi are already the right ones
                     i0=1
                     i1=1
                     
                  endif

                  if (ray_type == ray_type_low) then
                     call ang2pix_nest(nside, theta, phi,i0) 
                     i0=i0
                     i1=i0      
                  endif

                  if (ray_type == ray_type_gone) cycle 
      
                  sect_dir: do ipix=i0,i1  
       
                     if (ray_type == ray_type_high) then  
                        call pix2ang_nest(nside, ipix, theta, phi)       
                     endif

                      if ((rt_type == rtt_scatt).and.(cc_old == src_id)) then !!! THIS is important to select the correct scattered luminosity!!!
                         call assign_src_lum(theta,phi,src_lum)
                         src_lum_ray=src_lum
       
                      endif
      
                     clvl=lvl(src_id) !!! this should be here, very important! It restore the right nesting level for the emitting cell 

                     call ray_tracing(src_id, phi, theta,clvl, src_lum_ray,dplane_ray,prev_ray,isel_ray &
                          ,ray_type,ipix,cc_old)
                      

                  end do sect_dir

               end do ext_list
     
               call define_next_level(flag_k)  

               if (flag_k) then
        
                  exit 
               endif

            end do  !!! end loop on angular resolutions 

            select case (rt_type)
            case (rtt_dir_src, rtt_dir_cell, rtt_scatt)
               do i = 0, lnum -1
                  !$OMP ATOMIC      
                  lum_lost(i)=lum_lost(i)+lum_lost_temp(i)
               enddo
                  
            end select
         end do dir2

      end do  main_dir

      select case (rt_type)
      case (rtt_dir_src, rtt_dir_cell, rtt_scatt)
         if (print_psel_av) then 
            psel_av= psel_av/ ipsel_av            
            psel_av_arr(iterations_dustem,iterations, sid)=psel_av
            !OMP ATOMIC
            ipsel_av_tot=ipsel_av_tot+ipsel_av
         endif
         !print *, ipsel_av_tot
      end select

 end subroutine main_dir_loop

 !> Performs the ray-tracing loop to calculate the output specific intensity arrays i_obs(). 
 subroutine rt_loop_iobs

   integer :: i,j, il, src_id, ipix, clvl,sid
   integer :: i0, i1 
   integer :: isel_ray,cc_old
   real(kind=real64) :: src_lum_ray(0:lnum-1)
   real(kind=real64) :: prev_ray,dplane_ray,src_lum(0:lnum-1)
   real(kind=real64) :: theta, phi
   real(kind=real64) :: ro(3)
   integer :: ray_type
   
   if (tot_ndir == 0 .and. tot_ndir_in == 0.) return ! no observer
   if (lnum_node == 0) return   ! no locally stored i_obs 
   if (no_communications .and. .not. main_prc) return
   if (test_run) return

   call set_i_opacity_arrays(i0,i1)
   
   call omp_set_num_threads(nproc)  

   !$OMP PARALLEL DEFAULT(NONE), &
   !$OMP PRIVATE(il,theta,phi,src_lum_ray,i,j,src_id,clvl,dplane_ray & 
   !$OMP , prev_ray,isel_ray,ray_type,cc_old,ipix,sid,src_lum, ro ), &
   !$OMP SHARED(tot_ncell,cchild,lvl,i_obs_arr,dir_i_out,tot_ndir,tot_ndir_scaspe,max_lvl,src_cell,rt_type,scaspe_tot_arr,npix_hp,tot_ndir_in,i_obs_in_arr,rt_algorithm_ID, chunk_size, lnum, iq_sca_id, main_prc, ccoord_obs,i0, lnum_node, iq_sca_node, npix_hp_arr)

   call allocate_rt_loop_arrays

   !$OMP DO SCHEDULE(DYNAMIC,chunk_size)
   do i=0, tot_ncell-1  !!! loop on cells. Note that there is no MPI splitting of the work here. Every MPI process handles the wavelengths stored locally in the scaspe arrays 

      src_lum_ray=src_lum_value(i)
      src_id=i
      ray_type = ray_type_i_obs
      
      if ((cchild(i) /= -1).or.(sum(src_lum_ray) == 0)) cycle
      
      if (main_prc .and. mod(i,1000) ==0) write(*,*) src_id

      call calc_src_ccoord(i)

       if (rt_type == rtt_i_obs_dir_src ) then 
          sid = tot_ncell+src_cell(i)   ! iobs index          
       else 
          sid = i
       endif

      
      call cindex_to_ccindd(src_id,lvl(src_id),ccindd_nc)

       if (rt_type == rtt_i_obs) then 
          do il = 0, lnum -1 
             scaspe_temp_arr(il)%a = 0
          enddo
          do il = 0, lnum_node -1 
             scaspe_temp_arr(iq_sca_id(il))%a=scaspe_tot_arr(il)%a(:,i)  ! look scaspe_tot not scaspe. 
          end do
         
      endif


      if (tot_ndir > 0) then 

         do j=0,tot_ndir-1
            i_obs_temp=0
            theta=dir_i_out(j,1)
            phi=dir_i_out(j,2)

            if (rt_type == rtt_i_obs) then 
               if (rt_algorithm_ID /= rta_i_obs .and. rt_algorithm_ID /= rta_i_obs_dust) then 
                  do il = 0, lnum -1 
                     if (iq_sca_node(il)) then 
                        src_lum_ray(il)=scaspe_temp_arr(il)%a(npix_hp_arr(il+i0)+j)*npix_hp_arr(il+i0)
                     endif
                  end do
               else 
                  call assign_src_lum(theta,phi,src_lum)
                  src_lum_ray=src_lum
                  
               endif

            endif

            clvl=lvl(src_id) !!! this should be here, very important! It restore the right nesting level for the emitting cell 
            cc_old=src_id
            
            
            call ray_tracing(src_id, phi, theta,clvl, src_lum_ray,dplane_ray,prev_ray,isel_ray &
                  ,ray_type,ipix,cc_old)

             i_obs_arr(:,sid,j)=i_obs_arr(:,sid,j)+i_obs_temp(iq_sca_id)

           
          enddo
     
       endif

       ray_type = ray_type_i_obs_in

        if (tot_ndir_in > 0) then 

         do j=0,tot_ndir_in-1 
            i_obs_temp=0
            ro = ccoord_obs(:,j)
            
            call find_theta_phi_obs_in(ro,src_ccoord, theta,phi)         
            
            if (rt_type == rtt_i_obs) then
               ! NOTE the input theta here is in the HEALPix convention
               call assign_src_lum(theta,phi,src_lum)
               src_lum_ray=src_lum
              
            endif
            
            clvl=lvl(src_id) !!! this should be here, very important! It restore the right nesting level for the emitting cell 
            cc_old=src_id

            call ray_tracing(src_id, phi, theta,clvl, src_lum_ray,dplane_ray,prev_ray,isel_ray &
                  ,ray_type,ipix,cc_old)

             i_obs_in_arr(:,sid,j)=i_obs_in_arr(:,sid,j)+i_obs_temp(iq_sca_id)

          enddo
     
       endif

    end do

!$OMP END DO NOWAIT

    call deallocate_rt_loop_arrays

!$OMP END PARALLEL

 end subroutine rt_loop_iobs

 !> Performs the ray-tracing calculation for a single ray. The ray propagates through the model until it reach model edge or a blocking condition is met. The blocking condition depends on the RT calculation phase (see subroutine deposit()).
 !> @param [in] nc ID number of the cell originating the ray
 !> @param [in] al Phi angle of the ray direction
 !> @param [in] dl Theta angle of the ray direction (ray-direction^Z).
 !> @param [in] clvl Subdivision level of the cell originating the ray (but not used when the ray doesn't start from the cell nc
 !> @param [in] src_lum Luminosity associated with the ray
 !> @param [in] dplane_ray Distance from last intersection plane along axis perpendicular to that plane.
 !> @param [in] prev_ray Distance crossed by the blocked ray from which the current ray is derived.
 !> @param [in] isel_ray Direction perpendicular to the last intersection plane (0 = x, 1= y, 2=z).
 !> @param [in] ray_type Type of ray (e.g. ray_type_high(), ray_type_low(), ray_type_reco(), ray_type_gone(), ray_type_i_obs(), ray_type_i_obs_in() )
 !> @param [in] ipix HEALPix number of the ray
 !> @param [in] cc_old ID number of the last cell intersected by the blocked ray (this can be the same or not for the current ray)
 
subroutine ray_tracing(nc,al,dl, clvl,src_lum,dplane_ray,prev_ray,isel_ray &
                          ,ray_type,ipix,cc_old)
  
  
  integer :: i,j,k,l, incvec(3), cc, isel, clvl,isel_old
  integer :: isel_ray,npix_beam, ipix,cc_new
  integer :: cc_old, nc 
  REAL(KIND=real64),     INTENT(IN) ::  al,dl
  real(kind=real64) :: cproj(3)
  integer :: ray_status
  real(KIND=real64) :: psel, pabs(3),prev,b, length
  real(kind=real64) :: dplane_ray,prev_ray,src_lum_low_list(0:lnum-1), src_lum_high_list(0:lnum-1), beam_i
  integer :: ib
  real(KIND=real64) :: ray_intensity(0:lnum-1)
  real(kind=real64) :: cellsize
  real(KIND=real64), intent(in) :: src_lum(0:lnum-1)
  logical :: flag_beam,out, cell_found
  integer :: ray_type

  
  call calc_cproj(al,dl,cproj)
  
  ! this part calculates cc (the number of the starting cell) when the ray-tracing starts from the emitting cell (or point source) or from a different cell (when the ray has been blocked and put in the high or low ray list).  NOTE: ! This should be before cproj=1/cproj to work 
  if ((cc_old == nc)) then   
     cc=nc      ! current cell number
     prev=0.0   ! path already considered from source to last blocking point
     ray_status= ras_first_launch  ! this is used to select different cases when the ray is blocked 
     isel_old=0 ! index direction of last intersection plane
     flag_beam=.TRUE. ! this flag avoids recursion when ray crosses the first cell 
     clvl_old=-1 
     ccindd=ccindd_nc 
     cc_new = -1 
  else 
     
     call find_cc_new2(cc_old,cc_new,cproj,prev_ray,out,clvl)
     
     prev=prev_ray

     if (out) return 
     cc=cc_new  
     prev=prev   
     ray_status=ras_re_launched
     isel_old=isel_ray   
     flag_beam=.TRUE.
     clvl_old=-1
  endif
  
! determine sign increment vector 
  do i=1,3
     incvec(i)=max(min(int(cproj(i)/glepsilon),1),-1) ! this abstruse formula actually works. Takes into account when cproj is zero as well  
    
  end do

  cproj=1/cproj  !!! THIS CANNOT BE MOVED 
  
  ray_intensity=src_lum/(4*pi*csize(nc)*csize(nc)) ! ray specific intensity 
  
  iq_a = .FALSE. 
  call set_iq_list(ray_intensity, flag_beam)

  npix_beam=12*(nside*nside)
  beam_i=4*pi/npix_beam

  if (rt_type == rtt_dir_src .or. rt_type == rtt_dir_cell .or. rt_type == rtt_scatt) then 
      
     call calc_ffn_arr(al,dl,ipix)   !!! this calculates a factor needed multiple times later for the calculation of the scattered light contribution      

  endif

  ! print *, 'ray status=', ray_status 

   main: do   ! LOOP of ray-cell intersections (until border of the model or until the ray is blocked)

      ! lock cell
      !$OMP ATOMIC 
      lock_cell(cc) = lock_cell(cc) + 1 

      if (rt_type == rtt_dir_src .or. rt_type == rtt_dir_cell .or. rt_type == rtt_scatt) then 
         call set_iq_list(ray_intensity, flag_beam)
      endif
      
     ! Calculate length of ray-cell intersection by computing distances to cell walls and pick the closest intersection
     call calc_psel(incvec,psel,pabs,isel,clvl,cc,cproj)
        
!!$    print *, psel, cc,isel,clvl
!!$     print *, csize(cc)
!!$     print *, ccoord(:,cc)
!!$
!!$     do i=1,max_lvl
!!$        print *, ccindd(:,i)
!!$     enddo

     length=psel-prev
     
     if (abs(length) < 1E-5*csize_arr(clvl)) length=0
     
     if (length < 0) then 
        print *, 'length=', length
        print *, 'psel=', psel
        print *, 'prev=',prev
        print *, 'cell size=', csize(cc)
        print *, clvl
        print *, cc_old, cc_new
        print *, cc,nc
        print *, ccoord(:,cc)        
        print *, dplane_ray,isel_ray
        print *, 'what ? Something wrong in ray_tracing'
        call stop_prc
        
     endif

     if (cc /= nc .and. cc /= cc_new) flag_beam=.FALSE. 

     if (length /= 0.) then 
        call deposit(cc,al,dl,ipix,ray_intensity_dep,length,nc,psel,ray_status,flag_beam)
        ray_intensity(iq)= ray_intensity_dep
     endif
    
     select case (ray_status)
     case(ras_go_high)
        src_lum_high_list=ray_intensity*(4*pi*csize(nc)*csize(nc)) 
        dplane_ray=prev/cproj(isel_old)       
        ray_type=ray_type_high
        call add_to_high_ray_list(dl,al,nside,prev,isel_old,dplane_ray,src_lum_high_list,cc,ray_type)
        !$OMP ATOMIC 
        lock_cell(cc) = lock_cell(cc) - 1
        return
 
     case(ras_gone)
        ipsel_av=ipsel_av+1
        psel_av=psel_av+psel  !!! this is to estimate average psel
        !$OMP ATOMIC 
      lock_cell(cc) = lock_cell(cc) - 1
        return
  
     case(ras_go_low)
        src_lum_low_list=ray_intensity*(4*pi*csize(nc)*csize(nc))
        call add_to_low_ray_list(dl,al,nside,prev,isel_old,cproj,src_lum_low_list,cc)
        !$OMP ATOMIC 
        lock_cell(cc) = lock_cell(cc) - 1
        return
 
     end select
     
     if (ray_type == ray_type_i_obs_in) then 
        if ((vec_mod >= prev).and.(vec_mod <= psel)) then        
   ! correction to i_obs because observer is typically somewhere inside a cube
           i_obs_temp=ray_intensity*exp(dens_arr(:,cc)*(psel-vec_mod))
           !$OMP ATOMIC 
           lock_cell(cc) = lock_cell(cc) - 1
           return       
        endif
     endif

prev=psel  ! this CANNOT be moved 
isel_old=isel ! store old value of intersected plane direction
clvl_old=clvl

! Find the next cell 

do i=1,3  ! loop on xyz directions

   if (i.eq.isel) then   ! direction intersected plane 
          
      k=incvec(i) 
     
      do j=clvl,1,-1
         call increment(j,ccindd(i,j),k,k) 
         if (k == 0) exit  
              
      end do
         
      if (k.ne.0) then    ! if ray going outside model
           
         ipsel_av=ipsel_av+1
         psel_av=psel_av+psel    !!! this is to estimate average psel
              
         if (rt_type == rtt_i_obs_dir_cell .or. rt_type == rtt_i_obs_dir_src .or. rt_type == rtt_i_obs) then 
         i_obs_temp=ray_intensity
         endif
         !$OMP ATOMIC 
         lock_cell(cc) = lock_cell(cc) - 1                              
         return
      endif
              
      ! do j=clvl+1,max_nlevel   ! this take care of ccindd for  higher nested levels 
            
      ib=1      
     ! do j=clvl+1,min(clvl+1,max_lvl)  ! this take care of ccindd for  higher nested levels  ! achtung! This is correct only when there are no jump in grid levels between adjacent cells 
      do j=clvl+1,max_lvl
         if (j > 1) ib=2
         if (incvec(i).eq.1) then
            ccindd(i,j)=0
         elseif (incvec(i).eq.-1) then
            ccindd(i,j)=base(ib)-1
         else            
         endif
      end do
   else
      ! Get actual crossing coordinate relative to current cell
             
      cellsize=csize_arr(clvl)
      b=psel/cproj(i)+(src_ccoord(i)-ccoord(i,cc)+0.5_real64*cellsize)  
       
      ! Subdivide it into integer coordinates from the current level up
      !  do j=clvl+1,max_nlevel
      ib=1
      !do j=clvl+1,min(clvl+1,max_lvl)  ! achtung! This is correct only when there are no jump in grid levels between adjacent cells 
      do j=clvl+1,max_lvl
         
         if (j > 1) ib=2
         cellsize=csize_arr(j)
         ccindd(i,j)=int(b/cellsize)
                  
         if (ccindd(i,j) == base(ib)) then  ! this is here in case ray is intersecting cell edge
            ccindd(i,j)=base(ib)-1
                                
         endif
         b=mod(b,cellsize)
    
      end do
   endif

end do ! loop on xyz directions

! unlock cell
!$OMP ATOMIC 
lock_cell(cc) = lock_cell(cc) - 1 

call ccindd_to_cc(ccindd,cc,clvl, cell_found)

if (cell_found) cycle main 

print *, 'Here you should never get (RAY TRACING)'  
CALL STOP_PRC
  

end do main ! loop intersected cells 

end subroutine ray_tracing

!> Sets the subscript list iq() of the ray_intensity() elements higher than zero and sets ray_intensity_dep(). 
subroutine set_iq_list(ray_intensity, flag_beam)
integer :: i,k
real(kind=real64) :: ray_intensity(0:lnum-1)
logical :: flag_beam ! this is .TRUE. before the first ray-cell intersection

where (ray_intensity > 0) 
   iq_a = .TRUE.
end where

lnum_a = count(iq_a)

if (flag_beam) then  ! initialize lnum_a_old 
   lnum_a_old = size(iq)
endif

if (lnum_a_old /= lnum_a .or. flag_beam) then ! 
   deallocate(iq)
   allocate(iq(0:lnum_a-1))
   k=0
   do i = 0, lnum-1
      if (iq_a(i)) then 
         iq(k) = i
         k = k + 1
      endif
   end do
   
   deallocate(ray_intensity_dep)
   allocate(ray_intensity_dep(0:lnum_a-1))
   ray_intensity_dep = ray_intensity(iq)
   lnum_a_old = lnum_a
   
endif

end subroutine set_iq_list

!> Assigns the source luminosity within rt_loop() depending on the rt_type() and the cell/point source ID number.
!> @param [in] id Cell/point source ID number
function src_lum_value(id)

  integer :: id
  real(kind = real64) :: src_lum_value(0:lnum-1)

  if (rt_type == rtt_precalc_src .or. rt_type == rtt_dir_src .or. rt_type == rtt_i_obs_dir_src) then
     if (src_cell(id) /= -1 ) then 
        src_lum_value = lum_p_src_arr(:,src_cell(id))
     else
        src_lum_value =0
     endif

  elseif (rt_type == rtt_precalc_cell .or. rt_type == rtt_dir_cell .or. rt_type == rtt_i_obs_dir_cell) then
     src_lum_value=lumcell(:,id)

  elseif (rt_type == rtt_scatt .or. rt_type == rtt_i_obs) then
     src_lum_value=lumcell(:,id) ! this is used only to check that there is luminosity to be scattered in the cell 

  endif
        


end function src_lum_value

!> Assigns the source coordinates to src_ccoord() depending on rt_type().
subroutine calc_src_ccoord(id)
  
integer :: id 
  
  if (rt_type == rtt_precalc_src .or. rt_type == rtt_dir_src.or. rt_type == rtt_i_obs_dir_src) then
     if (src_cell(id) /= -1 ) then 
        src_ccoord=ccoord_p_src(:,src_cell(id))        
     else
        print *, 'you should not get here: calc_src_ccoord'
        call stop_prc
     endif
     
  else
     src_ccoord=ccoord(:,id)
  endif


end subroutine calc_src_ccoord

!> Calculates the ray-direction projections cproj() on the coordinate axis.
!> @param [in] al angle in rad between ray-projection on XY plane and X axis
!> @param [in] dl angle in rad between ray and Z axis
!> @param [out] cproj unit vector with ray direction components
subroutine calc_cproj(al,dl,cproj)
  real(kind=real64) :: cproj(3), al, dl

! projection coefficients for the given Healpix  direction
cproj(1)=sin(dl)*cos(al)
cproj(2)=sin(dl)*sin(al)
cproj(3)=cos(dl)

end subroutine calc_cproj

!> Finds the starting cell for a new ray (which is not necessarily the last cell crossed by the blocked ray from which the new ray is derived).
!> @param [in] cc_old ID number of the last cell intersected by the blocked ray
!> @param [out] cc_new ID number of the cell intersected by the new ray
!> @param [in] cproj Ray direction projections on coordinate axis
!> @param [in] prev_ray Distance crossed by the blocked ray
!> @param [out] out Flag equal to TRUE if ray is now outside the model
!> @param [out] clvl Subdivision level of the cell intersected by the new ray
subroutine find_cc_new2(cc_old,cc_new,cproj,prev_ray,out,clvl)
  
  real(kind=real64) :: cproj(3),dist,pcoord(3),rel_vec(3),cellsize
  real(kind=real64) :: psel,prev,prev_ray,norm_pcoord(3)
  integer :: clvl
  INTEGER :: cc_old,cc_new,outcube(3)
  integer :: i,ib
  logical :: out, cell_found

  out = .FALSE.  ! This is 1 if the vectors ends outside the model

  prev=prev_ray

  ! calculate ray previous path 
  pcoord=src_ccoord+prev*cproj

  ! check whether end vector is outside model
  do i=1,3
     if (abs(pcoord(i)) > modelsize/2.) then
        out = .TRUE. 
     
        return
     endif
  end do

  ! check whether vector end is within the current cell

  rel_vec=pcoord-ccoord(:,cc_old)

  cellsize=csize(cc_old)  !!! note that the new clvl has not been found yet. I cannot use csize_arr here  
  
  call calc_outcube(rel_vec,cellsize,outcube)

  if (sum(outcube) == 0) then
     clvl=lvl(cc_old)
     cc_new=cc_old
     call cindex_to_ccindd(cc_new,clvl,ccindd)
   
     return
  endif

  ! if not, find the different cell 

  norm_pcoord=(pcoord+modelsize/2.)

  ! find cell ccindd
  
  ib=1
  do i=1, max_lvl
     if (i > 1) ib = 2
 
     ccindd(1,i)= int(norm_pcoord(1)/csize_arr(i))
     ccindd(2,i)= int(norm_pcoord(2)/csize_arr(i))
     ccindd(3,i)= int(norm_pcoord(3)/csize_arr(i))

     where(ccindd(:,i) == base(ib))   !!! This is to handle the rare case when the ray intersect exactly the cell edge within the numerical accuracy 
        ccindd(:,i) = base(ib)-1
     endwhere
     
     norm_pcoord(1)=(norm_pcoord(1)-ccindd(1,i)*csize_arr(i))
     norm_pcoord(2)=(norm_pcoord(2)-ccindd(2,i)*csize_arr(i))
     norm_pcoord(3)=(norm_pcoord(3)-ccindd(3,i)*csize_arr(i))
 
  enddo

  ! find cell number

  call ccindd_to_cc(ccindd,cc_new,clvl,cell_found)

  if (.not.cell_found) then
     print *, 'cell not found in find_cc_new2'
     call stop_prc
     
  endif

  ! check whether result is correct 
  rel_vec=pcoord-ccoord(:,cc_new)

  cellsize=csize_arr(clvl)

  !print *, 'cellsize/2=', cellsize/2.
  call calc_outcube(rel_vec,cellsize,outcube)

  ! if the ray edge within the found cell, NO PROBLEM. Just return 
 
  if (sum(outcube) == 0) then 
     return
  else
     print *, 'something wrong: host cell not found (find_cc_new2)'
     print *, cc_old, cc_new
     print *, ccoord(:,cc_new)
     print *, pcoord
     print *, pcoord+modelsize/2.
     print *, modelsize
     print *, rel_vec
     print *, cellsize/2 
     print *, clvl
      do i=1, max_lvl
         print *, ccindd(:,i)
         end do 
     call stop_prc

  endif

end subroutine find_cc_new2

!> Calculates outcube, a xyz vector where an element is equal to 1 if the corresponding component of rel_vec is larger than the cell side divided by 2. This means that the selected point is outside the cell.
!> @param [in] rel_vec Vector from cell centre to considered point
!> @param [in] cellsize Size of the cell
!> @param [out] outcube Vector whose element are equal to 1 if the projection of rel_vec on the corresponding axis is larger than cellsize
 subroutine calc_outcube(rel_vec,cellsize,outcube)
   real(kind=real64) :: rel_vec(3), cellsize
   integer :: outcube(3),i

   outcube=0
   do i=1,3

      if ((abs(rel_vec(i))  > (1+1E-7)*cellsize/2.0d0)) then !.and.(i /= isel_ray)) then
         outcube(i)=1
      endif

   end do

 end subroutine calc_outcube

 !>  Precalculates the ads_arr() arrays for low values of nside() (<=256) Ads_arr() elements for larger nside values are calculated when needed during the RT loop. 
 subroutine pre_calc_ads_arr

   integer, parameter :: num_arr=7  ! number of precalculated ffn_arr arrays
   integer :: npix_rays(0:num_arr-1),i,j, iarr
   real(KIND=real64), allocatable :: theta_rays(:),phi_rays(:), sin_theta_rays(:),sin_phi_rays(:),cos_theta_rays(:),cos_phi_rays(:)
   real(KIND=real64) :: theta,phi,sin_dl,cos_dl,sin_al,cos_al

   ! skip if already run once
   if (allocated(ads_arr4)) return

   ! if not, continue
   if (main_prc) print *, 'precalculation of the Henyey-Greenstein cosine factors...'
   
   npix_rays(0)=12*(nside_min**2) ! nside = 4  
   npix_rays(1)=npix_rays(0)*4    ! nside = 8
   npix_rays(2)=npix_rays(1)*4    ! nside = 16
   npix_rays(3)=npix_rays(2)*4    ! nside = 32
   npix_rays(4)=npix_rays(3)*4    ! nside = 64
   npix_rays(5)=npix_rays(4)*4    ! nside = 128
   npix_rays(6)=npix_rays(5)*4    ! nside = 256

   allocate(ads_arr4(0:dim_npix_unique-1), ads_arr8(0:dim_npix_unique-1), ads_arr16(0:dim_npix_unique-1), ads_arr32(0:dim_npix_unique-1), ads_arr64(0:dim_npix_unique-1), ads_arr128(0:dim_npix_unique-1), ads_arr256(0:dim_npix_unique-1), ads_arr(0:dim_npix_unique-1), ffn_arr(0:dim_npix_unique-1))

   do i = 0, dim_npix_unique-1
      allocate(ads_arr4(i)%a(0:npix_unique(i)-1, 0:npix_rays(0)-1))
      allocate(ads_arr8(i)%a(0:npix_unique(i)-1, 0:npix_rays(1)-1))
      allocate(ads_arr16(i)%a(0:npix_unique(i)-1, 0:npix_rays(2)-1))
      allocate(ads_arr32(i)%a(0:npix_unique(i)-1, 0:npix_rays(3)-1))
      allocate(ads_arr64(i)%a(0:npix_unique(i)-1, 0:npix_rays(4)-1))
      allocate(ads_arr128(i)%a(0:npix_unique(i)-1, 0:npix_rays(5)-1))
      allocate(ads_arr256(i)%a(0:npix_unique(i)-1, 0:npix_rays(6)-1))
      allocate(ads_arr(i)%a(0:npix_unique(i)-1))
      allocate(ffn_arr(i)%a(0:npix_unique(i)-1))
   end do


   do iarr=0,num_arr-1 

      allocate(theta_rays(0:npix_rays(iarr)-1),phi_rays(0:npix_rays(iarr)-1),sin_theta_rays(0:npix_rays(iarr)-1), &
           cos_theta_rays(0:npix_rays(iarr)-1),sin_phi_rays(0:npix_rays(iarr)-1),cos_phi_rays(0:npix_rays(iarr)-1))

      do i=0,npix_rays(iarr)-1 
         call pix2ang_nest(nside_min*(2**iarr),i,theta,phi)
         theta_rays(i)=theta
         phi_rays(i)=phi
      enddo

      sin_theta_rays=sin(theta_rays) 
      cos_theta_rays=cos(theta_rays)
      sin_phi_rays=sin(phi_rays)
      cos_phi_rays=cos(phi_rays)

      
      do i=0,npix_rays(iarr)-1  ! loop on ray-directions 
 
         sin_dl=sin_theta_rays(i)
         cos_dl=cos_theta_rays(i)
         sin_al=sin_phi_rays(i)
         cos_al=cos_phi_rays(i)

         call calc_ads_arr(sin_dl, cos_dl, sin_al,cos_al)

         do j = 0, dim_npix_unique-1 
         
            if (iarr==0) then 
               ads_arr4(j)%a(:,i)= ads_arr(j)%a
            else if (iarr==1) then 
               ads_arr8(j)%a(:,i)= ads_arr(j)%a 
            else if (iarr==2) then 
               ads_arr16(j)%a(:,i)= ads_arr(j)%a
            else if (iarr==3) then 
               ads_arr32(j)%a(:,i)= ads_arr(j)%a
            else if (iarr==4) then 
               ads_arr64(j)%a(:,i)= ads_arr(j)%a
            else if (iarr==5) then 
               ads_arr128(j)%a(:,i)= ads_arr(j)%a
            else if (iarr==6) then
               ads_arr256(j)%a(:,i)= ads_arr(j)%a 
            endif
         
         end do
         
      enddo

      deallocate(theta_rays,phi_rays,sin_theta_rays,cos_theta_rays,sin_phi_rays,cos_phi_rays)

   enddo

   deallocate(ads_arr, ffn_arr)  ! Important! This subroutine is called outside the OPENMP loop. Inside the loop new ads_arr and ffn_arr are allocated and different for each thread.

   call print_done 
   
end subroutine  pre_calc_ads_arr

!> Calculates the array ffn_arr() needed to derive the angular distribution of the scattered light after a ray-cell intersection.
!> @param [in] al Phi angle of the ray direction
!> @param [in] dl Theta angle of the ray direction
!> @param [in] ipix HEALPix number of the ray direction
subroutine calc_ffn_arr(al,dl,ipix)
  
  real(kind=real64) :: al,dl,const, tot_hg
  real(kind=real64) :: sin_dl,cos_dl,sin_al,cos_al,ads
 ! type(var_arr_1d) :: ffn_arr(0:dim_npix_unique-1)
  integer :: i,ipix,k
  integer :: i0, i1
  integer :: ik_sca

  ! allocate arguments of ffn_arr
!!$  do i = 0, dim_npix_unique-1
!!$     allocate(ffn_arr(i)%a(0:npix_unique(i)-1))
!!$  end do
  

  if (nside <=256) then  

     do i = 0, dim_npix_unique-1 
        if (nside==4) then 
           ads_arr(i)%a=ads_arr4(i)%a(:,ipix)
        else if (nside ==8) then 
           ads_arr(i)%a=ads_arr8(i)%a(:,ipix)
        else if (nside ==16) then 
           ads_arr(i)%a=ads_arr16(i)%a(:,ipix)   
        else if (nside ==32) then 
           ads_arr(i)%a=ads_arr32(i)%a(:,ipix)  
        else if (nside ==64) then 
           ads_arr(i)%a=ads_arr64(i)%a(:,ipix)
        else if (nside ==128) then 
           ads_arr(i)%a=ads_arr128(i)%a(:,ipix) 
        else if (nside ==256) then 
           ads_arr(i)%a=ads_arr256(i)%a(:,ipix)       
        endif
     enddo

  else 
     
     if (dim_npix_unique > 0) then ! no need to do this in case of isotropic scattering. 

        sin_dl=sin(dl)
        cos_dl=cos(dl)
        sin_al=sin(al)
        cos_al=cos(al)

        call calc_ads_arr(sin_dl, cos_dl, sin_al,cos_al)

     endif

  endif

  ! calculate ffn_arr_mw at each wavelength

  call set_i_opacity_arrays(i0,i1)
  
  k = 0 
  do i = 0, lnum-1

    ! if (.not.iq_a(i) .or.(.not.iq_sca_node(i))) cycle   ! calculate only for the elements that will be used 

     if (.not.iq_sca_node(i)) cycle   ! calculate only for the elements stored within local mpi process 

     if (iq_a(i)) then   ! calculate only for elements that will be used  

        ik_sca = ik_sca_arr(i+i0)
        
        if (ik_sca /= -1) then ! non isotropic scattering 
           const=(1.0-gsca_arr(i+i0)*gsca_arr(i+i0))/npix_hp_arr(i+i0)
        
           ffn_arr(ik_sca)%a=(1.0+gsca_arr(i+i0)*gsca_arr(i+i0)-2.0*gsca_arr(i+i0)*ads_arr(ik_sca)%a)**(-1.5)
            
           ffn_arr(ik_sca)%a=const*ffn_arr(ik_sca)%a

           tot_hg=sum(ffn_arr(ik_sca)%a(0:npix_hp_arr(i+i0)-1))          
           
           ffn_arr(ik_sca)%a=ffn_arr(ik_sca)%a/tot_hg   ! normalization Henyey-Greenstein function ! note that this is correct also for the observer lines of sight direction because of formula in assign_src_lum. 

           ffn_arr_mw(k)%a = ffn_arr(ik_sca)%a

        else ! isotropic scattering

           ffn_arr_mw(k)%a = 1/npix_hp_arr(i+i0) ! Same comment as above about the normalization. 

        endif

     endif

     k = k + 1
     
  enddo

  ! deallocate arguments of ffn_arr
!!$  do i = 0, dim_npix_unique-1
!!$     deallocate(ffn_arr(i)%a)
!!$  end do

end subroutine calc_ffn_arr

!> Calculates the array ads_arr(), which contains the values of the factors cos(theta) within the Henyey-Greenstein function for a single ray direction and for all directions of the scaspe() arrays.
!> @param [in] sin_dl Sine of the theta angle of the ray direction
!> @param [in] cos_dl Cosine of the theta angle of the ray direction
!> @param [in] sin_al Sine of the phi angle of the ray direction
!> @param [in] cos_al Cosine of the phi angle of the ray direction
subroutine calc_ads_arr(sin_dl, cos_dl, sin_al,cos_al)
  real(kind=real64) :: sin_dl, cos_dl, sin_al,cos_al,const,tot_hg
  integer :: i

  do i = 0, dim_npix_unique -1 
  
     ads_arr(i)%a=(cos_theta_sca(i)%a*cos_dl+sin_theta_sca(i)%a*sin_dl*(cos_phi_sca(i)%a*cos_al+sin_phi_sca(i)%a*sin_al)) ! note that this formula takes into account our convention on theta (ray direction^Z)
  end do
 
end subroutine calc_ads_arr

!> Calculates the length of the current ray-cell intersection.
!> @param [in] incvec xyz vector equal to +1 or -1 depending if the ray is pointing towards the positive or negative direction of the corresponding axis. It is equal to 0 if the ray is along a perpendicular axis.
!> @param [out] psel Length of the ray-cell intersection
!> @param [out] pabs Array of possible lengths for the ray-cell intersection. These are the lengths in the cases the ray hits planes perpendicular to x,y or z directions when exiting from the intersected cell.
!> @param [out] isel Direction perpendicular to the plane hit by the ray when exiting the intersected cell (0 = x, 1= y, 2=z).
!> @param [in] clvl Subdivision level of the intersected cell
!> @param [in] cc ID number of the intersected cell
!> @param [in] cproj Ray direction projections on coordinate axis
subroutine calc_psel(incvec,psel,pabs,isel,clvl,cc,cproj)
  
  integer :: i,incvec(3),isel,clvl,cc,nc
  real(kind=real64) :: psel,pabs(3),cellsize,cproj(3)

  pabs=pabs_max
  
  if (clvl == clvl_old) then 
  
     where(incvec /=0) ! is =0 when the ray direction is along a perpendicular axis
        ! note: in this way you avoid using cproj when it is infinite 
        pabs=cproj*ccoord(:,cc)+pabs_arr
     endwhere
     
  else

     where(incvec /=0)
        pabs=cproj*(ccoord(:,cc)-src_ccoord+0.5d0*dble(incvec)*csize_arr(clvl)) 
     endwhere

     pabs_arr=cproj*((-src_ccoord)+0.5d0*dble(incvec)*csize_arr(clvl)) 

  endif

  
  if (pabs(1) < pabs(3)) then
     if (pabs(1) < pabs(2)) then
        psel = pabs(1)
        isel = 1
     else
        psel= pabs(2)
        isel =2
     endif
  else
     if (pabs(3) < pabs(2)) then
        psel=pabs(3)
        isel =3
     else
        psel =pabs(2)
        isel=2
     endif

  endif

end subroutine calc_psel

!> Adds the contributions to the radiation field energy density array and to the scattered luminosity array of the cell intersected by a ray. It also updates the specific intensity of the ray which is modified by the ray-cell intersection if the cell is dusty. Note that the ray_intensity() array in this subroutine is different from the one in the caller subroutine. 
!> @param [in] cc ID number of the intersected cell
!> @param [in,out] ray_intensity Specific intensity of the ray
!> @param [in] length Length of the ray-cell intersection
!> @param [in] nc ID number of the cell originating the ray
!> @param [in] psel Path crossed by the ray
!> @param [in, out] ray_status Label indicating whether the ray has to be blocked and if the angular density of the rays has to be increased or decreased. Possible values ('first_launch' ,'re-launched', 'gone', 'go_high' , 'go_low')
!> @param [in] flag_beam Flag with TRUE value when the ray has just started from the original cell or has been just relaunched from an intermediate point. It avoids infinite loops when new angular resolution is too low after being too high.!> @param [in] ray_type Ray_type ID  

subroutine deposit(cc,al,dl,ipix,ray_intensity,length,nc,psel,ray_status,flag_beam)

REAL(KIND=real64) :: csize_i,csize_f,area_i,beam_i,vol_f,psel,beam_f,i_av(0:lnum_a-1),area_f,beam
real(kind=real64) :: tau_c(0:lnum_a-1),exp_tau_c(0:lnum_a-1)
real(kind=real64) :: en_sca(0:lnum_a-1)
real(kind=real64) :: al,dl
integer, intent(IN) :: cc,nc,ipix
integer :: npix_beam,i, clvl
real(kind=real64) :: ray_intensity(0:lnum_a-1),u_beam(0:lnum_a-1)
real(KIND=real64), intent(IN) :: length
logical :: flag_ray,flag_beam
integer :: ray_status
real(kind=real64) :: dist_cell
integer :: i0, i1
logical :: wall_hit

npix_beam=12*(nside*nside)

tau_c=dens_arr(iq,cc)*length
exp_tau_c=exp(-tau_c)

clvl=lvl(nc)
csize_i=csize_arr(clvl)
area_i=carea_arr(clvl)
beam_i=4*pi/npix_beam  ! this npix is correct, it is the number of the initial cell

clvl=lvl(cc)
csize_f=csize_arr(clvl)
area_f=carea_arr(clvl)
vol_f=cvol_arr(clvl)
beam_f=carea_arr(clvl)/(psel*psel)
!dist_cell = calc_dist_cell(nc,cc)
!beam_f = carea_arr(clvl)/(dist_cell)**2

! average specific intensity (needed for the calculation of the radiation field contribution below)

if (cc /= nc) then 
   where (tau_c < 1E-008 .or. ray_intensity == 0)  
      i_av=ray_intensity     ! non dusty cells (trans=1)
   elsewhere 
      i_av=ray_intensity/tau_c*(1-exp_tau_c)  ! dusty cells !!! this formula can be derived by integrating I_av over a cellsize and then divide by the cellsize
   end where
elseif (cc == nc) then
    do i = 0, lnum_a -1 
    if (tau_c(i) < 1E-008 .or. ray_intensity(i) == 0) then   
       i_av(i)=ray_intensity(i)   
   else
      if (rt_type /= rtt_precalc_src .and. rt_type /= rtt_dir_src) then 
         ! this is the case of homogenous emission and absorption 
         tau_c(i)=dens_arr(iq(i),cc)*csize_i
         exp_tau_c(i)=exp(-tau_c(i))
         i_av(i)=ray_intensity(i)/(tau_c(i)*tau_c(i))*(exp_tau_c(i)+tau_c(i)-1)
      else
         ! this is the stellar point source case 
         i_av(i)=ray_intensity(i)/(dens_arr(iq(i),cc)*length)*(1-exp_tau_c(i))
      endif
   endif
enddo
endif

!!!---------------------------------------------------------------------
!!! pre-calculation until tau and psel limit (to estimate u_fest)
!!! -------------------------------------------------------------------

select case (rt_type) 
   
case(rtt_precalc_src, rtt_precalc_cell)
   ! case ray-density not high enough   
   if (psel > rad_lim*modelsize) then
      ray_status = ras_gone   ! arrived to limit precalculation region    
      return
   else

      if ((beam_f/bm_par <= beam_i).and.(.not.flag_beam)) then 
      
         ray_status = ras_go_high ! to be put in the high ray list  
         return 
         
      else

         if ((beam_f/bm_par_max < beam_i).or.(nside <= 2*nside_min).or.(flag_beam)) then

            if (lock_cell(cc) == 1) then
               u_fest_arr(iq,cc)=u_fest_arr(iq,cc)+i_av*length/cs/vol_f*area_i*beam_i
            else
               do i = 0, lnum_a-1
                  !$OMP ATOMIC
                  u_fest_arr(iq(i),cc)=u_fest_arr(iq(i),cc)+i_av(i)*length/cs/vol_f*area_i*beam_i
               enddo
            endif
               
         else
            
            ray_status = ras_go_low  ! for the low ray list 
            return

         endif
                  
      endif
      
   endif


!!!---------------------------------------------------------------------
!!! Final ray-tracing calculation  
!!! -------------------------------------------------------------------


case(rtt_dir_src, rtt_dir_cell, rtt_scatt)
   
   ! case ray-density not high enough
   if ((beam_f/bm_par <= beam_i).and.(.not.flag_beam)) then
     
      u_beam=i_av*csize_f/cs/vol_f*area_i*beam_f
      
      where (u_beam < en_lim*u_fest_arr(iq,cc))
         ray_intensity = 0
         lum_lost_temp(iq)=lum_lost_temp(iq)+i_av*beam_i*csize_i*csize_i
      end where

      wall_hit = .FALSE.
      call check_wall_hit(cc, wall_hit)
      
      if (any(ray_intensity > 0) .and. .not. wall_hit) then 
         ray_status = ras_go_high        
         return
      endif
          
      ray_status = ras_gone
      return
     
      ! case right ray-density 
   elseif ((beam_f/bm_par_max < beam_i).or.(nside <= 2*nside_min).or.(flag_beam)) then    
      
      ! contribution to radiation field energy density
      if (lock_cell(cc) == 1 ) then 
         
         u_final_arr(iq,cc)=u_final_arr(iq,cc)+i_av*length/cs/vol_f*area_i*beam_i
      else
         do i = 0, lnum_a -1
            !$OMP ATOMIC 
            u_final_arr(iq(i),cc)=u_final_arr(iq(i),cc)+i_av(i)*length/cs/vol_f*area_i*beam_i
         enddo
      endif
         
      ! scattering informations

      en_sca = 0

      call set_i_opacity_arrays(i0,i1)
      
      if (cc /= nc) then 
         where(tau_c >  1E-008)
            en_sca = ray_intensity*(1.0-exp_tau_c)*ksca_arr_norm(iq+i0)*area_i*beam_i ! NOTE: beam_i is different compared to the line above
         elsewhere ! use taylor series approximation for low tau_c
            en_sca = ray_intensity*tau_c*ksca_arr_norm(iq+i0)*area_i*beam_i
         end where
      else  
         where(tau_c >  1E-008)
            en_sca = ray_intensity/(tau_c)*(exp_tau_c+tau_c-1)*ksca_arr_norm(iq+i0)*area_i*beam_i
         elsewhere 
            en_sca = ray_intensity*tau_c/2*ksca_arr_norm(iq+i0)*area_i*beam_i
         end where
      endif
        
      call process_scatt_rad(cc,al,dl,ipix, en_sca)   
      
      ! case ray-density too high 
   else 
      
      ray_status = ras_go_low      
      return  
   endif

   
end select

!!!!--------------------
! update ray_intensity    
!!!--------------------

if (cc /= nc) then    
      ray_intensity=ray_intensity*exp_tau_c
else 
   
   where (tau_c < 1E-008)   !!! tau_c=0. because low numbers create problems in if statement
      ray_intensity=ray_intensity
         
   elsewhere  
  
      ray_intensity=ray_intensity/(tau_c)*(1-exp_tau_c)      
     
   end where
endif

end subroutine deposit

!> Calculates the distance between the centres of two cells. 
real function calc_dist_cell(nc,cc)
  integer :: nc, cc
  real(kind=real64) :: r1(3), r2(3)

  r1 = ccoord(:,nc)
  r2 = ccoord(:,cc)

  calc_dist_cell = sqrt(sum((r1-r2)**2))

end function  calc_dist_cell

!> Checks whather the cell cc is located beyond the walls set in the input (see e.g. x_wall_on() and x_wall_coord() ). Note that the wall is not effective during the point source ray-tracing loop.  
subroutine check_wall_hit(cc,wall_hit)
  integer :: cc
  logical :: wall_hit

  if (rt_type == rtt_dir_src) return 

  if (x_wall_on) then
     if (ccoord(1,cc) < x_wall_coord(1) .or. ccoord(1,cc) > x_wall_coord(2)) then
        wall_hit = .TRUE.
     endif
  endif

  if (y_wall_on) then
     if (ccoord(2,cc) < y_wall_coord(1) .or. ccoord(2,cc) > y_wall_coord(2)) then
        wall_hit = .TRUE.
     endif
  endif

  if (z_wall_on) then
     if (ccoord(3,cc) < z_wall_coord(1) .or. ccoord(3,cc) > z_wall_coord(2)) then
        wall_hit = .TRUE.
     endif
  endif

end subroutine check_wall_hit


!> Sets the values of the wall coordinates in model units.
subroutine set_walls

  if (main_prc) print *, 'setting walls if needed...'

  if (x_wall_on) then
     x_wall_coord = x_wall_coord*modelsize - modelsize/2.
     if (main_prc) then 
        print *, 'X-wall coordinates:'
        print *, x_wall_coord
     endif
  endif

  if (y_wall_on) then
     y_wall_coord = y_wall_coord*modelsize - modelsize/2.
     if (main_prc) then 
        print *, 'Y-wall coordinates:'
        print *, y_wall_coord
     endif
  endif

  if (z_wall_on) then
     z_wall_coord = z_wall_coord*modelsize - modelsize/2.
     if (main_prc) then 
        print *, 'Z-wall coordinates:'
        print *, z_wall_coord
     endif
  endif

  call print_done

end subroutine set_walls


 !> Adds the contribution of the scattered light luminosity to the crossed cell scaspe_arr() and scaspe_tot_arr() arrays. Note that this operation is normally performed through array vectorization since it is much faster than the corresponding do loop. However, when running the code in parallel, there is the possibility of race condition when two rays cross the same cell simultaneously. In order to avoid this, a do loop is performed with $OMP ATOMIC when lock_cell() = 1. 
!> Note also that if for the elements of scaspe_arr() not stored locally, the en_sca elements are stored in en_sca_list until transfering to the corresponding node.   
!> @param [in] cc ID number of the cell intersected by a ray
!> @param [in] en_Sca Total luminosity scattered by the intersected cell
!> @param [in] iw Wavelength index
subroutine process_scatt_rad(cc, al,dl,ipix,en_sca)
 
  real (Kind=real64) :: en_sca(0:lnum_a-1), en_sca_out,al,dl
  integer :: i,cc, iw,ipix
  integer :: i0, i1, k, il  
  integer :: i0l, i1l
  
  if (sum(en_sca) == 0.) return
  
  call set_i_opacity_arrays(i0l,i1l)
       
  if (lock_cell(cc) == 1) then 
      
     do iw = 0, lnum_a -1
        il = iq(iw)
        en_sca_out = en_sca(iw)
        if (iq_sca_node(il)) then
           call set_wavelength_index(il,k)
           scaspe_arr(k)%a(:,cc)=scaspe_arr(k)%a(:,cc) + en_sca_out*ffn_arr_mw(k)%a(:)
           
        else 

           call add_to_en_sca_list(cc,ipix,il,en_sca_out)

        endif
        
     enddo

        
     if (allocated(scaspe_tot_arr) .and. (.not. sequential_scattering)) then  ! this works only during scattering iterations

        do iw = 0, lnum_a-1
           il = iq(iw)
           if (iq_sca_node(il)) then 
              k = (il - id_mpi)/np_mpi
              if (print_scaspe_tot) then
                 i0 = 0
              else
                 i0 = npix_hp_arr(il+i0l)
              endif
              i1 = npix_arr(il+i0l)-1  
              scaspe_tot_arr(k)%a(i0:i1,cc)=scaspe_tot_arr(k)%a(i0:i1,cc)+ en_sca(iw)*ffn_arr_mw(k)%a(i0:i1)
           endif
       
        enddo

     endif

  else  ! multiple locked cell 

     do iw = 0, lnum_a -1
        il = iq(iw)
        en_sca_out = en_sca(iw)
        if (iq_sca_node(il)) then 
           call set_wavelength_index(il,k)
           do i = 0, npix_arr(il+i0l) -1
              !$OMP ATOMIC 
              scaspe_arr(k)%a(i,cc)=scaspe_arr(k)%a(i,cc) + en_sca_out*ffn_arr_mw(k)%a(i)   
           enddo
        else

           call add_to_en_sca_list(cc,ipix,il,en_sca_out)

           
        endif
     enddo
     
     if (allocated(scaspe_tot_arr) .and. (.not. sequential_scattering)) then  ! this works only during scattering iterations
        do iw = 0, lnum_a -1
           il = iq(iw)
           if (iq_sca_node(il)) then 
              k = (il - id_mpi)/np_mpi
              if (print_scaspe_tot) then
                 i0 = 0
              else
                 i0 = npix_hp_arr(il+i0l)
              endif
              i1 = npix_arr(il+i0l)-1  
              do i = i0 ,i1
                 !$OMP ATOMIC 
                 scaspe_tot_arr(k)%a(i,cc)=scaspe_tot_arr(k)%a(i,cc)+ en_sca(iw)*ffn_arr_mw(k)%a(i)
              enddo
           endif
        enddo
     endif

  endif

end subroutine process_scatt_rad

!> Sets the right wavelength subscript depending whether no_communications() is TRUE or FALSE. 
subroutine set_wavelength_index(il,k)
integer :: il, k 

if (.not. no_communications) then 
   k = (il - id_mpi)/np_mpi  ! this is always an integer (mathematically) 
else
   k = il
endif

   
end subroutine set_wavelength_index



!> Determines the next value of nside() in the next iteration for main_dir_loop().
!> @param [out] flag_k Logical variable equal to TRUE if there are no rays left in the ray lists. 
subroutine define_next_level(flag_k)
  
  integer :: nside_new
  logical :: flag_k 

  if (.not.no_ray) deallocate(ext_ray_list)

  flag_k=.FALSE.

  if (ihigh > 0) then 
     ! find minimum value of nside from high_ray_list and multiply by 2
     !nside_new=minval(ray_high_list(0:ihigh-1)%nside)*2
     ! define new nside 
     
     nside=nside*2
  elseif ((ihigh == 0).and.(ilow >0)) then 
     ! find maximum value of nside from low_ray_list and divide by 2
     !nside_new=maxval(ray_low_list(0:ilow-1)%nside)/2
     
     nside=nside/2
  else  
     ! choose proper flag so one exits the loop
     flag_k=.TRUE.  
  endif

end subroutine define_next_level


!> Finds maximum value in npix_arr within the range of wavelength used in the current RT algorithm 
subroutine calc_max_npix
  integer :: i0, i1 

  call set_i_opacity_arrays(i0,i1)
  
  max_npix = maxval(npix_arr(i0:i0+lnum-1))  ! here lnum not lnum_node

end subroutine calc_max_npix

!> Allocates and initialises the scattered luminosity array scaspe() as well as the sin/cos arrays ( theta_sca(), phi_sca(), sin_theta_sca(), cos_theta_sca(), sin_phi_sca(), cos_phi_sca()) for the fixed HEALPix directions and observer directions considered in the scaspe() arrays. Since these values depend only on kp_sca at each wavelength and on the observer lines-of-sight, the arrays are not stored for each wavelength but only for each kp_sca value present in kp_sca_arr().     
subroutine create_scaspe

  real(kind=real64) :: theta, phi
  integer :: i,j,k
  integer :: i0,i1
  
  if (main_prc) print *, 'create scaspe arrays...'

  call set_i_opacity_arrays(i0,i1)
    
  call prepare_scaspe_splitting 

  call calc_max_npix
  
  if ((rt_algorithm_ID == rta_i_obs .or. rt_algorithm_ID == rta_i_obs_dust).and. only_direct_rt) return ! no need to create scaspe arrays for making only direct light maps 

  !if (.not. allocated(scaspe))  allocate(scaspe(0:max_npix-1,0:tot_ncell-1))
  !scaspe=0

  if (.not. allocated(scaspe_arr)) then
     allocate(scaspe_arr(0:lnum_node-1))
     do i = 0, lnum-1
        if (iq_sca_node(i)) then
           call set_wavelength_index(i,k)
           allocate(scaspe_arr(k)%a(0:npix_arr(i+i0)-1, 0:tot_ncell-1))
           scaspe_arr(k)%a = 0
        endif
     end do
  endif
 
  if (rt_algorithm_ID == rta_i_obs .or. rt_algorithm_ID == rta_i_obs_dust) return  ! no need to create following arrays in the i_obs algorithms. Also, problem with tot_ndir_scaspe if the following lines are used. 

  if (.not. allocated(theta_sca))  then
     allocate(theta_sca(0:dim_npix_unique-1),phi_sca(0:dim_npix_unique-1),sin_theta_sca(0:dim_npix_unique-1), cos_theta_sca(0:dim_npix_unique-1),sin_phi_sca(0:dim_npix_unique-1),cos_phi_sca(0:dim_npix_unique-1))

     do i = 0, dim_npix_unique-1

        allocate(theta_sca(i)%a(0:npix_unique(i)-1), phi_sca(i)%a(0:npix_unique(i)-1), sin_theta_sca(i)%a(0:npix_unique(i)-1), cos_theta_sca(i)%a(0:npix_unique(i)-1), sin_phi_sca(i)%a(0:npix_unique(i)-1), cos_phi_sca(i)%a(0:npix_unique(i)-1))

     end do
     
  endif

  do i = 0, dim_npix_unique-1
     nside_sca = 2**kp_unique(i) 
     do j=0,npix_hp_unique(i)-1        
        call pix2ang_nest(nside_sca,j,theta,phi)
        theta_sca(i)%a(j)=theta
        phi_sca(i)%a(j)=phi
     enddo
  enddo

  if (tot_ndir_scaspe > 0 ) then
     do i = 0, dim_npix_unique-1
        theta_sca(i)%a(npix_unique(i)-tot_ndir_scaspe:npix_unique(i)-1)=dir_i_out(:,1)
        phi_sca(i)%a(npix_unique(i)-tot_ndir_scaspe:npix_unique(i)-1)=dir_i_out(:,2)
     enddo
  endif

  do i = 0, dim_npix_unique-1 
     sin_theta_sca(i)%a=sin(theta_sca(i)%a) 
     cos_theta_sca(i)%a=cos(theta_sca(i)%a)
     sin_phi_sca(i)%a=sin(phi_sca(i)%a)
     cos_phi_sca(i)%a=cos(phi_sca(i)%a)
  enddo
 
  call print_done

end subroutine create_scaspe

!> Prepares the splitting of the scaspe arrays within different node memory. It affects only calculations using multiple cluster nodes. It calculates the number of wavelengths lnum_node() for the scaspe array stored locally. 
subroutine prepare_scaspe_splitting 
  integer :: nw_per_node ! number of wavelengths per node
  integer :: nw_mod 
  integer :: diff 
  integer :: i, im,im1
  integer :: ierr 
  integer :: i0,i1 
  integer :: k

  call set_i_opacity_arrays(i0,i1)

if (.not. allocated(iq_sca_node))  allocate(iq_sca_node(0:lnum-1))  ! allocate and initialise
  iq_sca_node = .FALSE. 
  lnum_node = 0 

  if (np_mpi == 1 .or. no_communications) then   ! handle cases 1 MPI process or no_communications 
     lnum_node = lnum
     !iq_sca_node = .TRUE.
     
  else 

     diff = lnum - np_mpi   ! all other cases  
     nw_per_node = lnum/np_mpi 
     nw_mod = mod(lnum,np_mpi)

     select case(diff)
     case(0) 
        lnum_node = 1 
     case( : -1)
        if (id_mpi <= lnum-1) lnum_node = 1 
     case( 1 :)
        lnum_node = nw_per_node
        if (id_mpi <= nw_mod -1) then 
           lnum_node = lnum_node + 1
        endif
     case default
        print *, 'you should not get here: prepare_scaspe_splitting.'
        stop

     end select

  endif

! set iq_sca_node and iq_sca_id 
     
  im1 = lnum_node

  if (im1 > 0 .and. (.not. allocated(iq_sca_id))) allocate(iq_sca_id(0:im1-1))
  
  do im = 0, im1 -1 

     if (.not. no_communications) then 
        i = im*np_mpi + id_mpi
     else
        i = im
     endif
        
     if (i > lnum-1 ) cycle 
     
     iq_sca_node(i) = .TRUE.
     iq_sca_id(im) = i 
     
  enddo

! communicate number of wavelengths to all other MPI processes
  if (.not. allocated(lnum_node_arr)) allocate(lnum_node_arr(0:np_mpi-1))
  call mpi_allgather(lnum_node, 1, MPI_INTEGER, lnum_node_arr, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

! create im_lambda_arr(). Contains the ID of the MPI process hosting the scaspe and i_obs arrays for the corresponding wavelength. 
  if (.not. allocated(im_lambda_arr)) allocate(im_lambda_arr(0:lnum-1))

   
  do i = 0, lnum-1 

     im_lambda_arr(i) = mod(i, np_mpi)  ! not used in no_communications mode
   
  end do

! calculate lnum_node_maps 

lnum_node_maps = 0 
do i = 0, lnum_maps -1
   if (ind_out_maps(i)-i0 < 0 .or. ind_out_maps(i)-i0 > lnum-1) cycle ! this is to avoid segmentation fault in the line below
   if (iq_sca_node(ind_out_maps(i)-i0)) lnum_node_maps = lnum_node_maps +1 
end do 

! create iq_maps_id 
if (lnum_node_maps > 0 .and. .not. allocated(iq_maps_id)) then 
   allocate(iq_maps_id(0:lnum_node_maps-1))
   k = 0 
   do i = 0, lnum_node -1 
      if (any(ind_out_maps -i0 == iq_sca_id(i))) then 
         iq_maps_id(k) = i 
         k = k +1 
      endif
   end do
endif

! calculate tot_npix_arr_local 

  if (.not. allocated(tot_npix_arr_local)) allocate(tot_npix_arr_local(0:np_mpi-1))

  tot_npix_arr_local = 0
  do i = 0, lnum -1 

     im = im_lambda_arr(i)
     tot_npix_arr_local(im) = tot_npix_arr_local(im) + npix_arr(i+i0)

  end do
  
end subroutine prepare_scaspe_splitting


!> Allocates and initialises specific intensity arrays
!> i_obs() and i_obs_in(). 
subroutine create_i_obs
  
  if (main_prc) print *, 'create i_obs arrays...'
  
  tot_ncell_p_src = tot_ncell+tot_p_src
  
  if (tot_ndir > 0) then 

     if (.not. allocated(i_obs)) allocate(i_obs(0:tot_ncell_p_src-1,0:tot_ndir-1), i_obs_arr(0:lnum_node-1, 0:tot_ncell_p_src-1,0:tot_ndir-1))
     i_obs = 0
     i_obs_arr = 0

     if (rt_algorithm_ID == rta_dust .or. rt_algorithm_ID == rta_dust2d) then
        if (iterations_dustem == 1) then
           allocate(i_obs_arr_dir(0:lnum_node-1, 0:tot_ncell_p_src-1,0:tot_ndir-1), i_obs_arr_tot(0:lnum_node-1, 0:tot_ncell_p_src-1,0:tot_ndir-1 ))
           i_obs_arr_dir = 0
           i_obs_arr_tot = 0 
        endif
     endif

  endif

  if (tot_ndir_in > 0) then 
   if (.not. allocated(i_obs_in))  allocate(i_obs_in(0:tot_ncell_p_src-1,0:tot_ndir_in-1), i_obs_in_arr(0:lnum_node-1, 0:tot_ncell_p_src-1,0:tot_ndir_in-1))
     i_obs_in = 0
     i_obs_in_arr = 0

     if (rt_algorithm_ID == rta_dust .or. rt_algorithm_ID == rta_dust2d) then
        if (iterations_dustem == 1) then
           allocate(i_obs_in_arr_dir(0:lnum_node-1, 0:tot_ncell_p_src-1,0:tot_ndir_in-1), i_obs_in_arr_tot(0:lnum_node-1, 0:tot_ncell_p_src-1,0:tot_ndir_in-1))
           i_obs_in_arr_dir = 0
           i_obs_in_arr_tot = 0
        endif
     endif
     
  endif

  call print_done
 
end subroutine create_i_obs

!> Allocates and initialises the scaspe_tot() array. Unlike the scaspe() array, the scaspe_tot() includes all of the scattered light at each position. If printed on an output file, it can be later used to calculate scattered light maps at arbitrary line-of-sight directions. 
subroutine create_scaspe_tot
  integer :: i, k 
  integer :: i0, i1 
   
if (main_prc) print *, 'creating scaspe_tot array....'

call set_i_opacity_arrays(i0,i1)

call calc_max_npix

if (.not. allocated(scaspe_tot_arr)) then

   !allocate(scaspe_tot(0:max_npix-1,0:tot_ncell-1))

   allocate(scaspe_tot_arr(0:lnum_node-1))
   do i = 0, lnum-1 
      if (iq_sca_node(i)) then
         call set_wavelength_index(i,k)
         allocate(scaspe_tot_arr(k)%a(0:npix_arr(i+i0)-1,0:tot_ncell-1))
      endif
   end do
endif

!scaspe_tot=0
scaspe_tot_arr = scaspe_arr
  
call print_done

end subroutine create_scaspe_tot

!> Assigns the src_lum() values for rays carrying scattered light. This is done by simply locating the right element of the scaspe array. It is called in three cases: 1) for rays used in the radiation field calculation (rt_type = 'scatt'), 2) for rays used to derive i_obs() for an external observer in the rt_algorithm 'iobs', and 3) for rays used to derive the i_obs_in() value for internal observers. This version contains an interpolatation scheme over the sphere, but it is not used for the moment. 
subroutine assign_src_lum(theta,phi,src_lum)

real(kind=real64), INTENT(IN) :: theta,phi
!real(kind=real64) :: scaspe_temp(0:,0:)
real(kind=real64) :: theta0,phi0,src_lum(0:lnum-1),ang_dist(0:8),theta_p,phi_p
integer :: ip, ip_arr(0:dim_npix_unique-1), type, i_neigh(8),nneigh,i,i_arr(0:8),im1(0:0),in
integer :: i0, i1, il
integer :: ik_sca 

call set_i_opacity_arrays(i0,i1)

do i = 0, dim_npix_unique-1 
    nside_sca = 2**kp_unique(i)
    call ang2pix_nest(nside_sca,theta,phi,ip)
   ! call pix2ang_nest(nside_sca,i0,theta_p,phi_p)
    ip_arr(i) = ip
end do

if ((rt_type == rtt_scatt).or.(rt_type == rtt_i_obs)) then   

   do il = 0, lnum -1 
      ik_sca = ik_sca_arr(il+i0) 
      if (ik_sca /= -1) then ! non-isotropic scattering 
         ip = ip_arr(ik_sca)
      else  ! isotropic scattering 
         ip = 0
      endif
      src_lum(il)=scaspe_temp_arr(il)%a(ip)*npix_hp_arr(il+i0) !!! this would be the total cell luminosity if every pixel for the scattering would emit the same amount of energy. this is consistent with formulas in ray_tracing  
   end do
      
else  
   print *, 'this should not happen: assign_src_lum'
   print *, 'fix sphere interpolation to take into account variable npix.'
   CALL STOP_PRC
   
    ! interpolation from 4 closest spherical pixels weighted with distance 
   theta0=theta
   phi0=phi
   
call neighbours_nest(nside_sca, i0, i_neigh, nneigh)   !!! find negihbours central pixel

i_arr(0)=i0
i_arr(1:nneigh)=i_neigh(1:nneigh)

!!$print *, theta0,nneigh
!!$print *, i_arr


!! calculate angular distance 
do i=0,nneigh    ! no -1 here, one element more is zero element 

   
   call pix2ang_nest(nside_sca,i_arr(i),theta_p,phi_p)
   theta_p=theta_p
  ! print *, theta,phi,i_arr(i)
   ang_dist(i)=cos(theta0)*cos(theta_p)+sin(theta0)*sin(theta_p)*cos(phi_p-phi0)
   ang_dist(i)=acos(ang_dist(i))
   !print *, ang_dist(i)
   
enddo

!!! find the three shortest angular distances in the neighbour list 

if (nneigh == 8) then
in=3
elseif (nneigh ==7) then 
in=2
else
   print *, 'something wrong here: assign_src_lum'
   call stop_prc
   
endif 

do i=4,nneigh

   im1=maxloc(ang_dist(1:in))
   if (ang_dist(i) < ang_dist(im1(0))) then
      ang_dist(im1(0))=ang_dist(i)
      i_arr(im1(0))=i_arr(i)
   endif
enddo
 

!!! calculate weighted average for src_lum

do i =0, lnum-1 
 !  src_lum(i)=sum(scaspe_temp(i_arr(0:in),i)/ang_dist(0:in))/sum(1/ang_dist(0:in))*npix_hp  ! to be restored when you fix this part 
enddo

endif

end subroutine assign_src_lum

!> Sets the en_lim (\f$ f_U \f$) parameter, which checks if a ray carries negligible luminosity and can be blocked without significantly reducing the calculation accuracy.
subroutine set_en_lim

 !integer,allocatable  :: leaf_cells(:)
 integer :: tot_sources, tot_leaf_cells

! find out how many leaf cells (this takes into account reduce resolution here if appropriate)

 if (main_prc) print *, 'setting f_U.... '
 
 tot_leaf_cells = count(cchild == -1)
 
 tot_sources = tot_leaf_cells+tot_p_src-tot_spare_cells ! typically it does not change much if only cells with non zero emissivity and density are considered.  

 if (main_prc) print *, 'total number of emitters =', tot_sources 

 en_lim=accuracy/(tot_sources*0.25)

 if (main_prc) print *, 'f_U = ', en_lim

 call print_done 

end subroutine set_en_lim

!> Finds the theta and phi angles for the direction from a source (emitting cell or point source) to an internal observer.
!> @param [in] i ID number of the internal observer
!> @param [out] theta Source - observer direction theta angle
!> @param [out] phi Source - observer direction phi angle
subroutine find_theta_phi_obs_in(ro,rc,theta,phi)
IMPLICIT NONE 
real(kind= real64) :: vec(3),theta,phi,z(3),vec_xy,phi_ratio, theta_ratio,rc(3), ro(3)
integer :: i

z=[0,0,1]  ! zeta versor

!vec=ccoord_obs(:,i)-src_ccoord  ! vector from source to observer
!vec=ccoord_obs(:,i)-rc
vec = ro - rc 

!!$print *, 'start find theta'
!!$print *, ccoord_obs(:,i)
!!$print *, rc
!!$print*, 'vec =', vec


vec_mod=dot_product(vec,vec)**0.5   ! distance cell - observer

!print *, 'vec mod=', vec_mod

if (vec_mod > 1E-008) then   ! this condition is here in case observer is at the centre of a cell  
   theta_ratio=dot_product(z,vec)/vec_mod
  ! print *, 'theta ratio=', theta_ratio
   if (abs(theta_ratio-1.0) < 1E-8) theta_ratio=1.0  ! these lines are here to avoid NAN from acos below in case theta_ratio = 1 or -1 
   if (abs(theta_ratio+1.0) < 1E-8) theta_ratio=-1.0

   theta=acos(theta_ratio)

   !print *,'theta=',theta*180/pi
   vec_xy=sin(theta)  ! projection unit vector on XY plane

   if (vec_xy >  1E-8) then  
  ! print *, 'vec_xy=',vec_xy
   phi_ratio=vec(1)/(vec_mod*vec_xy)   
   else 
     phi_ratio=1.
  endif 

   if (abs(phi_ratio-1.0) < 1E-8) phi_ratio=1.0  ! these lines are here to avoid NAN from acos below in case phi_ratio = 1 or -1 
   if (abs(phi_ratio+1.0) < 1E-8) phi_ratio=-1.0

   phi=acos(phi_ratio)    ! vec(1) = vec_x  

   if (vec(2) < 0) phi=2*pi-phi    ! vec(2) = vec_y 

else 
   theta=0
   phi=0
   vec_mod=0

endif 

end subroutine find_theta_phi_obs_in 

!> Creates the arrays that stores the average path of the rays originating from each source. The array is expanded at the beginning of each scattering iteration/dust heating iteration in expand_psel_av_arr(). 
subroutine create_psel_av_arr
  
  ! return if psel_av_arr not going to be output 
  if (.not. print_psel_av) return

  ! expand instead of create if dust RT iterations > 1   
  if (iterations_dustem > 1) then 
     call expand_psel_av_arr
     return
  endif 
  
  if (main_prc) print *, 'creating psel_av_arr array....'

  allocate(psel_av_arr(0:0,0:0,0:tot_ncell+tot_p_src-1))
  psel_av_arr=0
  ipsel_av_tot=0  ! counter total number of rays 

  call print_done
  
end subroutine create_psel_av_arr

!> Expands the psel_av_arr() array at the beginning of each scattering iteration and of each dust heating iteration. 
subroutine expand_psel_av_arr

real(kind=real64),allocatable :: temp_arr(:,:,:)
integer :: n_it, n_it_sca

! return if psel_av_arr not going to be output 
if (.not. print_psel_av) return

! note "iterations" is already counting the next one
! the following if statement is here because iterations_dustem = 0 during stellar emission RT
if (iterations_dustem <= 1) then 
   n_it = 1   ! number of dust heating iteration
else 
   n_it = iterations_dustem
endif

n_it_sca = size(psel_av_arr,2)  ! number of scattering iterations

if (n_it == size(psel_av_arr,1) .and. iterations+1 <= n_it_sca) return ! no need to expand if the number of scattering iterations is less than the total number of scattering iterations in the previous dust heating iterations. 

if (main_prc) print *, 'expand psel_av_arr array...'

allocate(temp_arr(0:n_it-1, 0:n_it_sca-1,0:tot_ncell+tot_p_src-1))

temp_arr=psel_av_arr

deallocate(psel_av_arr)

n_it_sca = max(n_it_sca, iterations+1) ! this increases the scattering iteration dimension only if necessary during dust heating iterations 

allocate(psel_av_arr(0:n_it-1,0:n_it_sca-1,0:tot_ncell+tot_p_src-1))
psel_av_arr=0
psel_av_arr(0:n_it-1,0:iterations-1,:)=temp_arr

deallocate(temp_arr)

call print_done

end subroutine expand_psel_av_arr

!> Checks that the 3D grid is indeed axisymmetric in the rt_algorithm 2D. 
  subroutine check_grid_symmetry
  
    integer :: i,j,k,nc,link_list(0:6),cc

    if (main_prc) print *, 'checking grid symmetry....'
    
    do i=0, tot_ncell-1  !!! loop on cells 
 
       if ((cchild(i) /= -1)) cycle
       if ((ccoord(1,i) < 0).or.(ccoord(2,i) < 0).or.(ccoord(3,i) < 0)) cycle

       nc=i
    
       call find_linked_cells(nc,link_list)

       ! X symmetry
       cc=link_list(0)

       if ((ccoord(1,nc) /= -ccoord(1,cc)).or.(ccoord(2,nc) /= ccoord(2,cc)).or.(ccoord(3,nc) /= ccoord(3,cc))) then
          print *, 'STOP GRID NOT SYMMETRIC X', nc,cc
          CALL STOP_PRC
       endif

       ! Y symmetry
       cc=link_list(1)

       if ((ccoord(1,nc) /= ccoord(1,cc)).or.(ccoord(2,nc) /= -ccoord(2,cc)).or.(ccoord(3,nc) /= ccoord(3,cc))) then
          print *, 'STOP GRID NOT SYMMETRIC Y', nc,cc
          CALL STOP_PRC
       endif

       ! Z symmetry
       cc=link_list(2)

       if ((ccoord(1,nc) /= ccoord(1,cc)).or.(ccoord(2,nc) /= ccoord(2,cc)).or.(ccoord(3,nc) /= -ccoord(3,cc))) then
          print *, 'STOP GRID NOT SYMMETRIC Z', nc,cc
          CALL STOP_PRC
       endif

       ! XY symmetry
       cc=link_list(3)

       if ((ccoord(1,nc) /= -ccoord(1,cc)).or.(ccoord(2,nc) /= -ccoord(2,cc)).or.(ccoord(3,nc) /= ccoord(3,cc))) then
          print *, 'STOP GRID NOT SYMMETRIC XY', nc,cc
          CALL STOP_PRC
       endif

       ! XZ symmetry
       cc=link_list(4)

       if ((ccoord(1,nc) /= -ccoord(1,cc)).or.(ccoord(2,nc) /= ccoord(2,cc)).or.(ccoord(3,nc) /= -ccoord(3,cc))) then
          print *, 'STOP GRID NOT SYMMETRIC XZ', nc,cc
          CALL STOP_PRC
       endif

       ! YZ symmetry
       cc=link_list(5)

       if ((ccoord(1,nc) /= ccoord(1,cc)).or.(ccoord(2,nc) /= -ccoord(2,cc)).or.(ccoord(3,nc) /= -ccoord(3,cc))) then
          print *, 'STOP GRID NOT SYMMETRIC YZ', nc,cc
          CALL STOP_PRC
       endif

       ! XYZ symmetry
       cc=link_list(6)

       if ((ccoord(1,nc) /= -ccoord(1,cc)).or.(ccoord(2,nc) /= -ccoord(2,cc)).or.(ccoord(3,nc) /= -ccoord(3,cc))) then
          print *, 'STOP GRID NOT SYMMETRIC XYZ', nc,cc
          CALL STOP_PRC
       endif

    enddo

    if (main_prc) print *, 'GRID SYMMETRY CONFIRMED'
    call print_done 
    
end subroutine check_grid_symmetry

!> Finds the linked cells in the 2D mode (that is, the cells located at the symmetric points). 
!> @param [in] nc ID number of a cell
!> @param [out] link_list Array of ID number of the linked cells
subroutine find_linked_cells(nc,link_list)
  
  integer, allocatable :: ccindd(:,:),ccindd_new(:,:)
  integer :: nlevel,i,j,k,clvl,nc,cc,link_list(0:6), ib
  logical :: cell_found

  nlevel=lvl(nc) 
  
  allocate(ccindd(3,nlevel), ccindd_new(3,nlevel))
  ccindd=0
  ccindd_new=0

  call cindex_to_ccindd(nc,nlevel,ccindd)

  ! find -X copy

  ccindd_new=ccindd

  ib = 1 
  do i=1,nlevel

     if (i > 1) ib = 2

     ccindd_new(1,i) = base(ib) - 1 - ccindd(1,i) 
     !if (ccindd(1,i) == 2) ccindd_new(1,i) = 0
     !if (ccindd(1,i) == 0) ccindd_new(1,i) = 2   
     !  print *, i, (ccindd_new(k,i),k=1,3)
  enddo

  call ccindd_to_cc(ccindd_new,cc,clvl,cell_found)

  if ((clvl /= nlevel).or.(.not.cell_found)) then 
     print *, 'something wrong here: find_linked_cells'
     call stop_prc
  endif

  link_list(0)=cc

  ! find -Y copy

  !print *, 'y copy'

  ccindd_new=ccindd

  ib = 1 
  do i=1,nlevel

     if (i > 1) ib = 2
     ccindd_new(2,i) = base(ib) - 1 - ccindd(2,i) 
     !if (ccindd(2,i) == 2) ccindd_new(2,i) = 0
     !if (ccindd(2,i) == 0) ccindd_new(2,i) = 2
    
     !   print *, i, (ccindd_new(k,i),k=1,3)
  enddo

  call ccindd_to_cc(ccindd_new,cc,clvl,cell_found)

  if ((clvl /= nlevel).or.(.not.cell_found)) then 
     print *, 'something wrong here: find_linked_cells'
     call stop_prc
  endif

  link_list(1)=cc

  ! find -Z copy

  !print *, 'z copy'

  ccindd_new=ccindd

  ib = 1 
  do i=1,nlevel

     if (i > 1) ib = 2
     ccindd_new(3,i) = base(ib) - 1 - ccindd(3,i) 
     !if (ccindd(3,i) == 2) ccindd_new(3,i) = 0
     !if (ccindd(3,i) == 0) ccindd_new(3,i) = 2
 
     !   print *, i, (ccindd_new(k,i),k=1,3)
  enddo

  call ccindd_to_cc(ccindd_new,cc,clvl,cell_found)

  if ((clvl /= nlevel).or.(.not.cell_found)) then 
     print *, 'something wrong here: find_linked_cells'
     call stop_prc
  endif

  link_list(2)=cc 

  ! find -XY copy

  !print *, 'xy copy'

  ccindd_new=ccindd

  ib = 1 
  do i=1,nlevel

     if (i > 1) ib = 2
     ccindd_new(1,i) = base(ib) - 1 - ccindd(1,i)
     ccindd_new(2,i) = base(ib) - 1 - ccindd(2,i)
     
     !if (ccindd(1,i) == 2) ccindd_new(1,i) = 0
     !if (ccindd(1,i) == 0) ccindd_new(1,i) = 2
     !if (ccindd(2,i) == 2) ccindd_new(2,i) = 0
     !if (ccindd(2,i) == 0) ccindd_new(2,i) = 2
    
     !   print *, i, (ccindd_new(k,i),k=1,3)
  enddo

  call ccindd_to_cc(ccindd_new,cc,clvl,cell_found)

  if ((clvl /= nlevel).or.(.not.cell_found)) then 
     print *, 'something wrong here: find_linked_cells'
     call stop_prc
  endif

  link_list(3)=cc 
  ! find -XZ copy

  !print *, 'xz copy'
 
  ccindd_new=ccindd

  ib = 1 
  do i=1,nlevel

     if (i > 1) ib = 2
     ccindd_new(1,i) = base(ib) - 1 - ccindd(1,i)
     ccindd_new(3,i) = base(ib) - 1 - ccindd(3,i)
    ! if (ccindd(1,i) == 2) ccindd_new(1,i) = 0
    ! if (ccindd(1,i) == 0) ccindd_new(1,i) = 2
    ! if (ccindd(3,i) == 2) ccindd_new(3,i) = 0
    ! if (ccindd(3,i) == 0) ccindd_new(3,i) = 2
    
     !   print *, i, (ccindd_new(k,i),k=1,3)
  enddo

  call ccindd_to_cc(ccindd_new,cc,clvl,cell_found)

  if ((clvl /= nlevel).or.(.not.cell_found)) then 
     print *, 'something wrong here: find_linked_cells'
     call stop_prc
  endif

  link_list(4)=cc
 
  ! find -YZ copy

  !print *, 'yz copy'
 
  ccindd_new=ccindd

  ib = 1 
  do i=1,nlevel

     if (i > 1) ib = 2
     ccindd_new(2,i) = base(ib) - 1 - ccindd(2,i)
     ccindd_new(3,i) = base(ib) - 1 - ccindd(3,i)
    ! if (ccindd(2,i) == 2) ccindd_new(2,i) = 0
    ! if (ccindd(2,i) == 0) ccindd_new(2,i) = 2
    ! if (ccindd(3,i) == 2) ccindd_new(3,i) = 0
    ! if (ccindd(3,i) == 0) ccindd_new(3,i) = 2
    
     !   print *, i, (ccindd_new(k,i),k=1,3)
  enddo

  call ccindd_to_cc(ccindd_new,cc,clvl,cell_found)

  if ((clvl /= nlevel).or.(.not.cell_found)) then 
     print *, 'something wrong here: find_linked_cells'
     call stop_prc
  endif

  link_list(5)=cc
 
  ! find -XYZ copy

  !print *, 'xyz copy'
 
  ccindd_new=ccindd

  ib = 1 
  do i=1,nlevel

     if (i > 1) ib = 2
     ccindd_new(1,i) = base(ib) - 1 - ccindd(1,i)
     ccindd_new(2,i) = base(ib) - 1 - ccindd(2,i)
     ccindd_new(3,i) = base(ib) - 1 - ccindd(3,i)
     
   !  if (ccindd(1,i) == 2) ccindd_new(1,i) = 0
   !  if (ccindd(1,i) == 0) ccindd_new(1,i) = 2
   !  if (ccindd(2,i) == 2) ccindd_new(2,i) = 0
   !  if (ccindd(2,i) == 0) ccindd_new(2,i) = 2
   !  if (ccindd(3,i) == 2) ccindd_new(3,i) = 0
   !  if (ccindd(3,i) == 0) ccindd_new(3,i) = 2
    
     !   print *, i, (ccindd_new(k,i),k=1,3)
  enddo

  call ccindd_to_cc(ccindd_new,cc,clvl,cell_found)

  if ((clvl /= nlevel).or.(.not.cell_found)) then 
     print *, 'something wrong here: find_linked_cells'
     call stop_prc
  endif

  link_list(6)=cc

  ! print *, ccoord(:,nc)
  ! print *, ccoord(:,link_list)
  deallocate(ccindd, ccindd_new)
  
end subroutine find_linked_cells

!> Checks that no more than one source is allowed at the grid origin in the 2D mode. This is needed because it is the only case where symmetries can be exploited in the same way as for the emitting cells.  
subroutine check_2d_src

  if (main_prc) print *, 'checking point source position...'
  
  if (tot_p_src < 1) then
     call print_done
     return
  else if (tot_p_src == 1) then 
     if (sum(ccoord_p_src(:,0)) /= 0 ) then 
        print *, 'only ONE source positioned in the origin is allowed in 2D mode'
        call stop_prc
     endif
     return
  else 
     print *, 'only ONE source positioned in the origin is allowed in 2D mode'
     call stop_prc
  endif

  call print_done
    
end subroutine check_2d_src


!> Takes care of symmetrising the radiation fields and scaspe() arrays at the end of the RT loops in the RT 2D mode. R and z symmetry is exploited here. 
subroutine fix_symmetry
  integer :: ierr

  select case (rt_type)
  case(rtt_precalc_cell)
     call fix_symmetry_part1
  case(rtt_dir_cell)
     if (iterations_dustem <= 1) then 
        call fix_symmetry_part2
     else
        call fix_symmetry_part3  ! this takes into account that the arrays to symmetrize have corresponding non zero prev arrays
     endif
     !lum_lost=lum_lost*8  
  case(rtt_scatt) 
     call fix_symmetry_part3
     !lum_lost=(lum_lost-lum_lost_prev)*8+lum_lost_prev  
  case default
  end select

  ! add lum_lost contribution. Remember that lum_lost and lum_lost_prev are zero during precalculation. lum_lost_prev is zero during direct stellar light processing and during the direct dust emission processing (but only for the first dust heating iterations). The following formula shold be able to cope with all cases.  
  lum_lost=(lum_lost-lum_lost_prev)*8+lum_lost_prev 

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  
end subroutine 

!> Handles the array symmetrization after the RT precalculation in 2D mode.
subroutine fix_symmetry_part1
 
 
    integer :: i,j,k,nc,link_list(0:6),cc
    real(kind=real64) :: u_arr(0:lnum-1,0:6)
    
    do i=0, tot_ncell-1  !!! loop on cells 

       if ((cchild(i) /= -1)) cycle
       if ((ccoord(1,i) < 0.).or.(ccoord(2,i) < 0.).or.(ccoord(3,i) < 0.)) cycle

       nc=i

       call find_linked_cells(nc,link_list)

       do j = 0, lnum -1

          ! Add delta U to positive cube

          u_arr(j,:)=u_fest_arr(j,link_list)
          
          u_fest_arr(j,i)=u_fest_arr(j,i)+sum(u_arr(j,:))
    
          ! make it symmetric 
          u_fest_arr(j,link_list)=u_fest_arr(j,i)
       enddo
 
    enddo
     
end subroutine fix_symmetry_part1

!> Takes care of symmetrizing radiation fields and scaspe() arrays at the end of direct light processing for the RT 2D mode.
 subroutine fix_symmetry_part2
   
   integer :: i,j,k,cc,nc,link_list(0:6)
   real (kind=real64) :: u_arr(0:lnum-1,0:6)
   type(var_arr_1d), allocatable :: scaspe_x(:),scaspe_y(:),scaspe_z(:),scaspe_xy(:),scaspe_xz(:),scaspe_yz(:),scaspe_xyz(:),scaspe_fix(:) 
   integer :: i0, i1
   integer :: ik_sca
   integer :: icc

   ! calculate re-order indexing for scaspe arrays 

   if (.not. allocated(ix)) call calc_scaspe_indices

  ! npix=12*nside_sca**2
   call set_i_opacity_arrays(i0,i1)

   allocate(scaspe_x(0:dim_npix_unique-1), scaspe_y(0:dim_npix_unique-1), scaspe_z(0:dim_npix_unique-1), scaspe_xy(0:dim_npix_unique-1), scaspe_xz(0:dim_npix_unique-1), scaspe_yz(0:dim_npix_unique-1), scaspe_xyz(0:dim_npix_unique-1), scaspe_fix(0:dim_npix_unique-1))

   do i = 0, dim_npix_unique-1 

      allocate(scaspe_x(i)%a(0:npix_hp_unique(i)-1), scaspe_y(i)%a(0:npix_hp_unique(i)-1), scaspe_z(i)%a(0:npix_hp_unique(i)-1), scaspe_xy(i)%a(0:npix_hp_unique(i)-1), scaspe_xz(i)%a(0:npix_hp_unique(i)-1), scaspe_yz(i)%a(0:npix_hp_unique(i)-1), scaspe_xyz(i)%a(0:npix_hp_unique(i)-1), scaspe_fix(i)%a(0:npix_hp_unique(i)-1))

   end do
  

   do i=0, tot_ncell-1  !!! loop on cells 

      if ((cchild(i) /= -1)) cycle
      if ((ccoord(1,i) < 0).or.(ccoord(2,i) < 0).or.(ccoord(3,i) < 0)) cycle

      nc=i

      call find_linked_cells(nc,link_list)

      do j = 0, lnum-1 
      
         ! Add delta U to positive cube

         u_arr(j,:)=u_final_arr(j,link_list)

         u_final_arr(j,i)=u_final_arr(j,i)+sum(u_arr(j,:))

         ! make it symmetric 
         
         u_final_arr(j,link_list)=u_final_arr(j,i)
   
         ! make scaspe array symmetric (this is trickier)

         if (iq_sca_node(j)) then
            call set_wavelength_index(j,k)
            
            ik_sca = ik_sca_arr(j+i0) 

            if (ik_sca /= -1) then ! non isotropic scattering 

               ! x symmetry
               cc=link_list(0)
               scaspe_x(ik_sca)%a=scaspe_arr(k)%a(:,cc)
               scaspe_x(ik_sca)%a=scaspe_x(ik_sca)%a(ix(ik_sca)%int)

               ! y symmetry
               cc=link_list(1)
               scaspe_y(ik_sca)%a=scaspe_arr(k)%a(:,cc)
               scaspe_y(ik_sca)%a=scaspe_y(ik_sca)%a(iy(ik_sca)%int)
               
               ! z symmetry
               cc=link_list(2)
               scaspe_z(ik_sca)%a=scaspe_arr(k)%a(:,cc)
               scaspe_z(ik_sca)%a=scaspe_z(ik_sca)%a(iz(ik_sca)%int)
               
               ! xy symmetry
               cc=link_list(3)
               scaspe_xy(ik_sca)%a=scaspe_arr(k)%a(:,cc)
               scaspe_xy(ik_sca)%a=scaspe_xy(ik_sca)%a(ixy(ik_sca)%int)
            
               ! xz symmetry
               cc=link_list(4)
               scaspe_xz(ik_sca)%a=scaspe_arr(k)%a(:,cc)
               scaspe_xz(ik_sca)%a=scaspe_xz(ik_sca)%a(ixz(ik_sca)%int)
               
               ! yz symmetry
               cc=link_list(5)
               scaspe_yz(ik_sca)%a=scaspe_arr(k)%a(:,cc)
               scaspe_yz(ik_sca)%a=scaspe_yz(ik_sca)%a(iyz(ik_sca)%int)
               
               ! xyz symmetry
               cc=link_list(6)
               scaspe_xyz(ik_sca)%a=scaspe_arr(k)%a(:,cc)
               scaspe_xyz(ik_sca)%a=scaspe_xyz(ik_sca)%a(ixyz(ik_sca)%int)

               ! add scaspe from symmetric cells 
               scaspe_arr(k)%a(:,i)=scaspe_arr(k)%a(:,i)+scaspe_x(ik_sca)%a+scaspe_y(ik_sca)%a+scaspe_z(ik_sca)%a+scaspe_xy(ik_sca)%a+scaspe_xz(ik_sca)%a+scaspe_yz(ik_sca)%a+scaspe_xyz(ik_sca)%a
               
               ! now make scaspe symmetric 
               scaspe_fix(ik_sca)%a=scaspe_arr(k)%a(:,i)
            
               cc=link_list(0)
               scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(ix(ik_sca)%int)
               cc=link_list(1)
               scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(iy(ik_sca)%int)
               cc=link_list(2)
               scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(iz(ik_sca)%int)
               cc=link_list(3)
               scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(ixy(ik_sca)%int)
               cc=link_list(4)
               scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(ixz(ik_sca)%int)
               cc=link_list(5)
               scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(iyz(ik_sca)%int)
               cc=link_list(6)
               scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(ixyz(ik_sca)%int)

            else ! isotropic scattering 

               ! add scaspe from symmetric cells 
               do icc = 0, 6
                  scaspe_arr(k)%a(:,i)=scaspe_arr(k)%a(:,i)+scaspe_arr(k)%a(:,cc)
               end do

               ! symmetrize
               do icc = 0, 6
                  cc=link_list(icc)
                  scaspe_arr(k)%a(:,cc)= scaspe_arr(k)%a(:,i)                  
               end do
            
            endif

         endif
         
      enddo
   end do

   deallocate(scaspe_x,scaspe_y,scaspe_z,scaspe_xy,scaspe_xz,scaspe_yz,scaspe_xyz,scaspe_fix) 

   

 end subroutine fix_symmetry_part2

 !> Symmetrises the radiation fields and scaspe() arrays at the end of each scattering iteration in the RT 2D mode.
 subroutine fix_symmetry_part3
   
   integer :: i,j,k,cc,nc,link_list(0:6)
   real (kind=real64) :: u_arr(0:6)
   type(var_arr_1d), allocatable :: scaspe_x(:),scaspe_y(:),scaspe_z(:),scaspe_xy(:),scaspe_xz(:),scaspe_yz(:),scaspe_xyz(:),scaspe_fix(:) 
   integer :: i0, i1 
   integer :: ik_sca
   integer :: icc 
   real (kind=real64) :: scaspe_diff(0:0)

   ! calculate re-order indexing for scaspe arrays 

   if (.not. allocated(ix)) call calc_scaspe_indices

   call set_i_opacity_arrays(i0,i1)

   allocate(scaspe_x(0:dim_npix_unique-1), scaspe_y(0:dim_npix_unique-1), scaspe_z(0:dim_npix_unique-1), scaspe_xy(0:dim_npix_unique-1), scaspe_xz(0:dim_npix_unique-1), scaspe_yz(0:dim_npix_unique-1), scaspe_xyz(0:dim_npix_unique-1), scaspe_fix(0:dim_npix_unique-1))

   do i = 0, dim_npix_unique-1 

      allocate(scaspe_x(i)%a(0:npix_hp_unique(i)-1), scaspe_y(i)%a(0:npix_hp_unique(i)-1), scaspe_z(i)%a(0:npix_hp_unique(i)-1), scaspe_xy(i)%a(0:npix_hp_unique(i)-1), scaspe_xz(i)%a(0:npix_hp_unique(i)-1), scaspe_yz(i)%a(0:npix_hp_unique(i)-1), scaspe_xyz(i)%a(0:npix_hp_unique(i)-1), scaspe_fix(i)%a(0:npix_hp_unique(i)-1))

   end do

   do i=0, tot_ncell-1  !!! loop on cells 

      if ((cchild(i) /= -1)) cycle
      if ((ccoord(1,i) < 0).or.(ccoord(2,i) < 0).or.(ccoord(3,i) < 0)) cycle

      nc=i

      call find_linked_cells(nc,link_list)

      do j = 0, lnum -1 
      
         ! Add delta U to positive cube

         u_arr(:)=u_final_arr(j,link_list)-u_fest_arr(j,link_list)
         call remove_negative(u_arr)
         
         u_final_arr(j,i)=u_final_arr(j,i)+sum(u_arr(:))

         ! make it symmetric 
   
         u_final_arr(j,link_list)=u_final_arr(j,i)
  
         ! make scaspe array symmetric (this is trickier)

         if (iq_sca_node(j)) then
            call set_wavelength_index(j,k)

             ik_sca = ik_sca_arr(j+i0) 
             
             if (ik_sca /= -1) then 

                if (.not. sequential_scattering) then

                   ! x symmetry
                   cc=link_list(0)
                   scaspe_x(ik_sca)%a=scaspe_arr(k)%a(:,cc)-scaspe_prev_arr(k)%a(:,cc)
                   call remove_negative(scaspe_x(ik_sca)%a)
                   scaspe_x(ik_sca)%a=scaspe_x(ik_sca)%a(ix(ik_sca)%int)
               
                   ! y symmetry
                   cc=link_list(1)
                   scaspe_y(ik_sca)%a=scaspe_arr(k)%a(:,cc)-scaspe_prev_arr(k)%a(:,cc)
                   call remove_negative(scaspe_y(ik_sca)%a)
                   scaspe_y(ik_sca)%a=scaspe_y(ik_sca)%a(iy(ik_sca)%int)
               
                   ! z symmetry
                   cc=link_list(2)
                   scaspe_z(ik_sca)%a=scaspe_arr(k)%a(:,cc)-scaspe_prev_arr(k)%a(:,cc)
                   call remove_negative(scaspe_z(ik_sca)%a)
                   scaspe_z(ik_sca)%a=scaspe_z(ik_sca)%a(iz(ik_sca)%int)
               
                   ! xy symmetry
                   cc=link_list(3)
                   scaspe_xy(ik_sca)%a=scaspe_arr(k)%a(:,cc)-scaspe_prev_arr(k)%a(:,cc)
                   call remove_negative(scaspe_xy(ik_sca)%a)
                   scaspe_xy(ik_sca)%a=scaspe_xy(ik_sca)%a(ixy(ik_sca)%int)
               
                   ! xz symmetry
                   cc=link_list(4)
                   scaspe_xz(ik_sca)%a=scaspe_arr(k)%a(:,cc)-scaspe_prev_arr(k)%a(:,cc)
                   call remove_negative(scaspe_xz(ik_sca)%a)
                   scaspe_xz(ik_sca)%a=scaspe_xz(ik_sca)%a(ixz(ik_sca)%int)
               
                   ! yz symmetry
                   cc=link_list(5)
                   scaspe_yz(ik_sca)%a=scaspe_arr(k)%a(:,cc)-scaspe_prev_arr(k)%a(:,cc)
                   call remove_negative(scaspe_yz(ik_sca)%a)
                   scaspe_yz(ik_sca)%a=scaspe_yz(ik_sca)%a(iyz(ik_sca)%int)
         
                   ! xyz symmetry
                   cc=link_list(6)
                   scaspe_xyz(ik_sca)%a=scaspe_arr(k)%a(:,cc)-scaspe_prev_arr(k)%a(:,cc)
                   call remove_negative(scaspe_xyz(ik_sca)%a)
                   scaspe_xyz(ik_sca)%a=scaspe_xyz(ik_sca)%a(ixyz(ik_sca)%int)
               
                else

                   ! x symmetry
                   cc=link_list(0)
                   scaspe_x(ik_sca)%a=scaspe_arr(k)%a(:,cc)
                   scaspe_x(ik_sca)%a=scaspe_x(ik_sca)%a(ix(ik_sca)%int)
                   
                   ! y symmetry
                   cc=link_list(1)
                   scaspe_y(ik_sca)%a=scaspe_arr(k)%a(:,cc)
                   scaspe_y(ik_sca)%a=scaspe_y(ik_sca)%a(iy(ik_sca)%int)
                   
                   ! z symmetry
                   cc=link_list(2)
                   scaspe_z(ik_sca)%a=scaspe_arr(k)%a(:,cc)
                   scaspe_z(ik_sca)%a=scaspe_z(ik_sca)%a(iz(ik_sca)%int)
                   
                   ! xy symmetry
                   cc=link_list(3)
                   scaspe_xy(ik_sca)%a=scaspe_arr(k)%a(:,cc)
                   scaspe_xy(ik_sca)%a=scaspe_xy(ik_sca)%a(ixy(ik_sca)%int)
                   
                   ! xz symmetry
                   cc=link_list(4)
                   scaspe_xz(ik_sca)%a=scaspe_arr(k)%a(:,cc)
                   scaspe_xz(ik_sca)%a=scaspe_xz(ik_sca)%a(ixz(ik_sca)%int)
                   
                   ! yz symmetry
                   cc=link_list(5)
                   scaspe_yz(ik_sca)%a=scaspe_arr(k)%a(:,cc)
                   scaspe_yz(ik_sca)%a=scaspe_yz(ik_sca)%a(iyz(ik_sca)%int)
                   
                   ! xyz symmetry
                   cc=link_list(6)
                   scaspe_xyz(ik_sca)%a=scaspe_arr(k)%a(:,cc)
                   scaspe_xyz(ik_sca)%a=scaspe_xyz(ik_sca)%a(ixyz(ik_sca)%int)
                
                endif
         
                ! add scaspe from symmetric cells 
                scaspe_arr(k)%a(:,i)=scaspe_arr(k)%a(:,i)+scaspe_x(ik_sca)%a+scaspe_y(ik_sca)%a+scaspe_z(ik_sca)%a+scaspe_xy(ik_sca)%a+scaspe_xz(ik_sca)%a+scaspe_yz(ik_sca)%a+scaspe_xyz(ik_sca)%a
                ! if (print_scaspe_tot) then
                if (.not. sequential_scattering) then 
                   scaspe_tot_arr(k)%a(:,i)=scaspe_tot_arr(k)%a(:,i)+scaspe_x(ik_sca)%a+scaspe_y(ik_sca)%a+scaspe_z(ik_sca)%a+scaspe_xy(ik_sca)%a+scaspe_xz(ik_sca)%a+scaspe_yz(ik_sca)%a+scaspe_xyz(ik_sca)%a
                endif
             
                ! now make scaspe symmetric 
                scaspe_fix(ik_sca)%a=scaspe_arr(k)%a(:,i)
 
                cc=link_list(0)
                scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(ix(ik_sca)%int)
                cc=link_list(1)
                scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(iy(ik_sca)%int)
                cc=link_list(2)
                scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(iz(ik_sca)%int)
                cc=link_list(3)
                scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(ixy(ik_sca)%int)
                cc=link_list(4)
                scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(ixz(ik_sca)%int)
                cc=link_list(5)
                scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(iyz(ik_sca)%int)
                cc=link_list(6)
                scaspe_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(ixyz(ik_sca)%int)
             
                !if (print_scaspe_tot) then
                if (.not. sequential_scattering) then 
                   scaspe_fix(ik_sca)%a=scaspe_tot_arr(k)%a(:,i)
 
                   cc=link_list(0)
                   scaspe_tot_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(ix(ik_sca)%int)
                   cc=link_list(1)
                   scaspe_tot_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(iy(ik_sca)%int)
                   cc=link_list(2)
                   scaspe_tot_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(iz(ik_sca)%int)
                   cc=link_list(3)
                   scaspe_tot_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(ixy(ik_sca)%int)
                   cc=link_list(4)
                   scaspe_tot_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(ixz(ik_sca)%int)
                   cc=link_list(5)
                   scaspe_tot_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(iyz(ik_sca)%int)
                   cc=link_list(6)
                   scaspe_tot_arr(k)%a(:,cc)=scaspe_fix(ik_sca)%a(ixyz(ik_sca)%int)
             
                endif

             else 

                if (.not. sequential_scattering) then
                   
                   do icc = 0, 6 
                      cc=link_list(icc)
                      scaspe_diff = scaspe_arr(k)%a(:,cc)-scaspe_prev_arr(k)%a(:,cc)
                      call remove_negative(scaspe_diff)
                      scaspe_arr(k)%a(:,i) = scaspe_arr(k)%a(:,i) + scaspe_diff
                      scaspe_tot_arr(k)%a(:,i) = scaspe_tot_arr(k)%a(:,i) + scaspe_diff
                   enddo

                else 

                   do icc = 0, 6 
                      cc=link_list(icc)
                      scaspe_arr(k)%a(:,i) = scaspe_arr(k)%a(:,i) + scaspe_arr(k)%a(:,cc)
                   enddo

                endif

                ! symmetrize 
                do icc = 0, 6 
                   cc=link_list(icc)
                   scaspe_arr(k)%a(:,cc) = scaspe_arr(k)%a(:,i)
                   if (.not. sequential_scattering) then 
                      scaspe_tot_arr(k)%a(:,cc) = scaspe_tot_arr(k)%a(:,i)
                   endif
                enddo

             endif

          endif
         
       enddo
    enddo

   deallocate(scaspe_x,scaspe_y,scaspe_z,scaspe_xy,scaspe_xz,scaspe_yz,scaspe_xyz,scaspe_fix)    

 end subroutine fix_symmetry_part3

!> Assigns zero to all negative elements of an array. Used in fix_symmetry_part3() to remove negative due to numerical accuracy. Remember that this could give problems in case you want to work with negative intensities.  
 subroutine remove_negative(array)
   real(kind=real64) :: array(0:)

   where(array < 0)
      array =0
   end where

 end subroutine remove_negative

!> Same as remove_negative() but for 2D arrays. 
 subroutine remove_negative_2d(array)
   real(kind=real64) :: array(0:, 0:)

   where(array < 0)
      array =0
   end where

 end subroutine remove_negative_2d


!> Allocates and initialises a new scaspe array ( scaspe_prev_arr()) which is used while symmetrising the scaspe_arr() array in the RT 2D mode or in the sequential scattering algorithm (see sequential_scattering()).
subroutine create_scaspe_prev
  integer :: i,k 
  integer :: i0, i1

 ! integer :: npix

  if (main_prc) print *, 'creating scaspe_prev array...'
  
 ! npix=12*nside_sca**2
  call calc_max_npix

  call set_i_opacity_arrays(i0,i1)

!if (.not. allocated(scaspe_prev))  allocate(scaspe_prev(0:max_npix-1, 0:tot_ncell-1))

if (.not. allocated(scaspe_prev_arr)) then 
   allocate(scaspe_prev_arr(0:lnum_node-1))
   do i = 0, lnum-1 
      if (iq_sca_node(i)) then
         call set_wavelength_index(i,k)
         allocate(scaspe_prev_arr(k)%a(0:npix_arr(i+i0)-1, 0:tot_ncell-1))
         scaspe_prev_arr(k)%a = 0
      endif
   end do
endif

!  scaspe_prev = 0.
  scaspe_prev_arr = scaspe_arr

  call print_done

end subroutine create_scaspe_prev

!> Derives the re-ordering scaspe() indices for symmetrising the scaspe arrays in the RT 2D mode.  
subroutine calc_scaspe_indices
 
  integer :: i,iu, is(0:0)
  type(var_arr_1d), allocatable :: mask(:)  ! note that the logical array is used for this derived data type variable (l not a!)
  real(kind=real64) :: phi_sca0,theta_sca0
  
!  npix=12*nside_sca**2

  allocate(ix(0:dim_npix_unique-1),iy(0:dim_npix_unique-1), iz(0:dim_npix_unique-1), ixy(0:dim_npix_unique-1),ixz(0:dim_npix_unique-1), iyz(0:dim_npix_unique-1),ixyz(0:dim_npix_unique-1),mask(0:dim_npix_unique-1))

  do i = 0, dim_npix_unique -1 

    allocate(ix(i)%int(0:npix_hp_unique(i)-1),iy(i)%int(0:npix_hp_unique(i)-1), iz(i)%int(0:npix_hp_unique(i)-1), ixy(i)%int(0:npix_hp_unique(i)-1),ixz(i)%int(0:npix_hp_unique(i)-1), iyz(i)%int(0:npix_hp_unique(i)-1),ixyz(i)%int(0:npix_hp_unique(i)-1),mask(i)%l(0:npix_hp_unique(i)-1)) 

    ix(i)%int = 0 
    iy(i)%int = 0 
    iz(i)%int = 0 
    ixy(i)%int = 0 
    ixz(i)%int = 0
    iyz(i)%int = 0
    ixyz(i)%int = 0 

  end do
  
  do iu = 0, dim_npix_unique -1 

     do i=0, npix_hp_unique(iu)-1
        ! X symmetry
        if ((phi_sca(iu)%a(i) >= 0).and.(phi_sca(iu)%a(i) <= pi/2)) then
           phi_sca0=pi-phi_sca(iu)%a(i)
           mask(iu)%l = (theta_sca(iu)%a == theta_sca(iu)%a(i)) 
           is=minloc(abs(phi_sca(iu)%a-phi_sca0),mask(iu)%l)-1
           ix(iu)%int(i)=is(0)
           ix(iu)%int(is)=i
        
        endif

        if ((phi_sca(iu)%a(i) > pi).and.(phi_sca(iu)%a(i) <= 3./2.*pi)) then ! look here only ">" (not ">=") at first term. otherwise double swap because of previous if statement
           phi_sca0=2*pi-(phi_sca(iu)%a(i)-pi)
           mask(iu)%l = (theta_sca(iu)%a == theta_sca(iu)%a(i)) 
           is=minloc(abs(phi_sca(iu)%a-phi_sca0),mask(iu)%l)-1
           ix(iu)%int(i)=is(0)
           ix(iu)%int(is)=i
      
        endif
   
        
        ! Y symmetry
        if ((phi_sca(iu)%a(i) >= 0).and.(phi_sca(iu)%a(i) <= pi)) then
           phi_sca0=2*pi-phi_sca(iu)%a(i)
           if (phi_sca0 == 2*pi) phi_sca0=0
           mask(iu)%l = (theta_sca(iu)%a == theta_sca(iu)%a(i)) 
           is=minloc(abs(phi_sca(iu)%a-phi_sca0),mask(iu)%l)-1
           iy(iu)%int(i)=is(0)
           iy(iu)%int(is)=i
           
        endif

        ! Z symmetry
        if ((theta_sca(iu)%a(i) >= 0).and.(theta_sca(iu)%a(i) <= pi/2)) then
           theta_sca0=pi-theta_sca(iu)%a(i)
           mask(iu)%l = (phi_sca(iu)%a == phi_sca(iu)%a(i)) 
           is=minloc(abs(theta_sca(iu)%a-theta_sca0),mask(iu)%l)-1
           iz(iu)%int(i)=is(0)
           iz(iu)%int(is)=i
           
        endif
        
        ! XY symmetry
        if ((phi_sca(iu)%a(i) >= 0).and.(phi_sca(iu)%a(i) < pi)) then ! important < pi here not <= pi (otherwise swap/problems because there is no 2pi )
           phi_sca0=phi_sca(iu)%a(i)+pi
           mask(iu)%l = (theta_sca(iu)%a == theta_sca(iu)%a(i)) 
           is=minloc(abs(phi_sca(iu)%a-phi_sca0),mask(iu)%l)-1
           ixy(iu)%int(i)=is(0)
           ixy(iu)%int(is)=i
           
        endif

        ! XZ symmetry
        if ((phi_sca(iu)%a(i) >= 0).and.(phi_sca(iu)%a(i) <= pi/2)) then
           phi_sca0=pi-phi_sca(iu)%a(i)
           theta_sca0 = pi-theta_sca(iu)%a(i)
           is=minloc(abs(theta_sca(iu)%a-theta_sca0))-1
           theta_sca0=theta_sca(iu)%a(is(0))
           mask(iu)%l = (theta_sca(iu)%a == theta_sca0)
           is=minloc(abs(phi_sca(iu)%a-phi_sca0),mask(iu)%l)-1
           ixz(iu)%int(i)=is(0)
           ixz(iu)%int(is)=i
           
        endif

        if ((phi_sca(iu)%a(i) > pi).and.(phi_sca(iu)%a(i) <= 3./2.*pi)) then
           phi_sca0=2*pi-(phi_sca(iu)%a(i)-pi)
           theta_sca0 = pi-theta_sca(iu)%a(i)
           is=minloc(abs(theta_sca(iu)%a-theta_sca0))-1
           theta_sca0=theta_sca(iu)%a(is(0))
           mask(iu)%l = (theta_sca(iu)%a == theta_sca0)
           is=minloc(abs(phi_sca(iu)%a-phi_sca0),mask(iu)%l)-1
           ixz(iu)%int(i)=is(0)
           ixz(iu)%int(is)=i
        endif

        ! YZ symmetry
        if ((phi_sca(iu)%a(i) >= 0).and.(phi_sca(iu)%a(i) <= pi)) then 
           phi_sca0=2*pi-phi_sca(iu)%a(i)
           if (phi_sca0 == 2*pi) phi_sca0=0
           theta_sca0 = pi-theta_sca(iu)%a(i)
           is=minloc(abs(theta_sca(iu)%a-theta_sca0))-1
           theta_sca0=theta_sca(iu)%a(is(0))
           mask(iu)%l = (theta_sca(iu)%a == theta_sca0) 
           is=minloc(abs(phi_sca(iu)%a-phi_sca0),mask(iu)%l)-1
           iyz(iu)%int(i)=is(0)
           iyz(iu)%int(is)=i     
        endif
        
        ! XYZ symmetry
        if ((phi_sca(iu)%a(i) >= 0).and.(phi_sca(iu)%a(i) < pi)) then ! important < pi here not <= pi (otherwise swap/problems because there is no 2pi )
           phi_sca0=phi_sca(iu)%a(i)+pi
           theta_sca0 = pi-theta_sca(iu)%a(i)
           is=minloc(abs(theta_sca(iu)%a-theta_sca0))-1
           theta_sca0=theta_sca(iu)%a(is(0))
           mask(iu)%l = (theta_sca(iu)%a == theta_sca0)  
           is=minloc(abs(phi_sca(iu)%a-phi_sca0),mask(iu)%l)-1
           ixyz(iu)%int(i)=is(0)
           ixyz(iu)%int(is)=i
           
        endif
   
     end do

  enddo

  deallocate(mask)

end subroutine calc_scaspe_indices

!> Allocates and initializes en_sca_list() array and the counter count_en_sca().
subroutine create_en_sca_list

  if (rt_type == rtt_precalc_cell .or. rt_type == rtt_precalc_src) return 

  if (allocated(en_sca_list)) then 
     if (size(en_sca_list) /= size_en_sca_list) then 
        deallocate(en_sca_list,temp_en_sca_list,ind_en_sca_list)
     else
        return  
     endif
  endif
  
  allocate(en_sca_list(0:size_en_sca_list-1),temp_en_sca_list(0:size_en_sca_list-1), ind_en_sca_list(0:size_en_sca_list-1))

  en_sca_list(:)%cc = 0
  en_sca_list(:)%nside = 0 
  en_sca_list(:)%ipix = 0
  en_sca_list(:)%il = 0 
  en_sca_list(:)%en_sca = 0. 

  count_en_sca = 0

end subroutine create_en_sca_list

!> Adds additional element to the en_sca_list() of the OPENMP thread. If the list contains the maximum number of elements allowed, the transfering subroutine handle_mpi_transfer is called. 
subroutine add_to_en_sca_list(cc,ipix,il,en_sca_out)
  integer :: i,cc,ipix,il
  real(kind=real64) :: en_sca_out

  i = count_en_sca 
  
  en_sca_list(i)%cc = cc
  en_sca_list(i)%nside = nside 
  en_sca_list(i)%ipix = ipix
  en_sca_list(i)%il = il 
  en_sca_list(i)%en_sca = en_sca_out

  count_en_sca = count_en_sca + 1

! three possibilities. The first seems to be faster. 
  if (count_en_sca == size_en_sca_list .or. any(handling_mpi_thread_arr)) then 
!  if (count_en_sca == size_en_sca_list .or. count(handling_mpi_thread_arr) > nproc/2-1) then
! if (count_en_sca == size_en_sca_list) then
     handling_mpi_thread_arr(ThreadID) = .TRUE.
     call handle_mpi_transfers(MPI_TRANSFER_EN_SCA,-1)
     
  endif

end subroutine add_to_en_sca_list

!> Processes the en_sca_list arrays of each OPENMP thread. It joins them in a single array. Then the master thread of each MPI process handles the communication and processing. 
!> @param el_out Number of en_sca elements to be sent to each MPI process.
!> @param i0_el_out List of positions of the first elements in the en_sca_list_all() array to be sent to each MPI process. 
subroutine process_en_sca_list
  integer :: i,i0,i1,im,it
  integer :: el_out(0:np_mpi-1), i0_el_out(0:np_mpi-1) 
  integer :: ierr
    
  el_out = 0
  i0_el_out = 0

  if (count_en_sca > 0) then 
     
     call sort_en_sca_list(el_out, i0_el_out)  

  endif

  el_out_arr(ThreadID,:) = el_out
  i0_el_out_arr(ThreadID,:) = i0_el_out
  count_en_sca_arr(threadID) = count_en_sca

  !$OMP BARRIER
  !$OMP MASTER 

  count_en_sca_tot = sum(count_en_sca_arr)

  i0_count_en_sca = 0  

  im = 0 
  it = 0
  do i = 1, nproc*np_mpi-1 
     ! this is the list of indeces used below to allow each thread to assign its en_sca_list elements to en_sca_list_all at the right place. That is, such that the en_sca_list_all elements are sorted in terms of MPI processes to which they have to be sent.  
     i0_count_en_sca(i) = i0_count_en_sca(i-1) + el_out_arr(it,im)
     it = it + 1 
     if (it == nproc) then 
        it = 0 
        im = im + 1 
     endif
  enddo
    
 if (size(en_sca_list_all) /= count_en_sca_tot) then       
    deallocate(en_sca_list_all)
    allocate(en_sca_list_all(0:count_en_sca_tot-1))
 endif

  !$OMP END MASTER 
  !$OMP BARRIER 
   
 do im = 0, np_mpi -1 
    
    i0 = i0_count_en_sca(im*nproc+threadID)
    i1 = i0 + el_out(im) -1 
 
    if (el_out(im) > 0) then 
       en_sca_list_all(i0:i1) = en_sca_list(i0_el_out(im):i0_el_out(im)+el_out(im)-1)  ! join en_sca_list arrays 
    endif

 end do

 !$OMP BARRIER 
 !$OMP MASTER  

 call prepare_en_sca_transfer

 !$OMP END MASTER
 !$OMP BARRIER   

 call en_sca_transfer()

 count_en_sca = 0

 !$OMP BARRIER   
  
end subroutine process_en_sca_list

!> Sorts the elements of en_sca_list for each thread before passing it to en_sca_list_all(). This avoids to make the sorting later within the master block. 
!> @param el_out Number of elements to be sent to each MPI process in en_sca_list
!> @param i0_el_out After the sorting of en_sca_list, this array gives the position of the first element in en_sca_list of the blocks to be sent to each MPI process. 
subroutine sort_en_sca_list(el_out, i0_el_out)
integer :: el_out(0:np_mpi-1)
integer :: i0_el_out(0:np_mpi-1), k_arr(0:np_mpi-1) 
integer :: i, i0, i1,il,im, iel

! copy en_sca_list_all to another array (used for sorting) 
temp_en_sca_list = en_sca_list

el_out = 0
i0 = 0 
i1 = count_en_sca

 do i = i0, i1-1  ! this loop derives the destination of each en_sca_list element and the total number to be sent to each mpi process. Remember that this is for the thread value 

    il = en_sca_list(i)%il   
    
    im = im_lambda_arr(il)

    ind_en_sca_list(i) = im 

    el_out(im) = el_out(im) + 1 
    
 end do

! determine position first element in the en_sca_list to be sent for each processor 

i0_el_out = 0

do i = 1, np_mpi-1 

i0_el_out(i) = i0_el_out(i-1) + el_out(i-1) 

enddo

! sort en_sca_list_all elements in order of MPI process

 k_arr = 0
 
 do i = i0, i1-1

    im = ind_en_sca_list(i)

    iel = k_arr(im)+i0_el_out(im)

    en_sca_list(iel) = temp_en_sca_list(i)

    k_arr(im) = k_arr(im)+1
 
 end do

end subroutine sort_en_sca_list



!> Processes en_sca_list_received within each MPI process. This is done by all threads bu the master thread.  
subroutine process_en_sca_received(im)
  integer :: i,k,im, ip
  integer :: i0, i1, n_el_in
  integer :: cc, cc_old
  integer :: il, ipix,ipix_old
  integer :: nside_prev, nside_old
  real(kind=real64) :: al,dl,en_sca_out
  logical :: iq_a_prev(0:lnum-1)
  integer :: ierr
  type(var_arr_1d), allocatable :: ffn_arr_mw_prev(:)
  integer :: i0l, i1l

  call set_i_opacity_arrays(i0l,i1l)

  allocate(ffn_arr_mw_prev(0:lnum_node-1))  ! maybe this should go outside this routine so it is not repeated 
  do i = 0, lnum-1 
     if (iq_sca_node(i)) then
        call set_wavelength_index(i,k)
        allocate(ffn_arr_mw_prev(k)%a(0:npix_arr(i+i0l)-1))
     endif
  enddo

! store values of variables to be restored later 
  nside_prev = nside  ! this is the nside before the ray propagation was interrupted
  iq_a_prev = iq_a
  ffn_arr_mw_prev = ffn_arr_mw
  
! initialize loop variables 
  nside_old = -1 
  ipix_old = -1 
  cc_old = -1 
  n_el_in = el_in_mpi(im)
  i0 = i0_el_in_mpi(im) 
  i1 = i0 + n_el_in-1
 
! loop on en_sca_list_received elements 
  do i = i0, i1
     if (en_sca_id_thread(i) /= ThreadID) cycle
     cc = en_sca_list_received(i)%cc
     nside = en_sca_list_received(i)%nside  ! note: this is a module variable!
     ipix = en_sca_list_received(i)%ipix 
     if (cc /= cc_old) then 
        call put_lock_to_cell(cc)
        if (i /= i0) call remove_lock_from_cell(cc_old)
        cc_old = cc
     endif

     if (nside > 256) then ! al and dl not needed in calc_ffn_arr if nside <= 256
        call pix2ang_nest(nside,ipix,dl,al)
     else 
        al = 0 
        dl = 0
     endif
     il = en_sca_list_received(i)%il
     en_sca_out = en_sca_list_received(i)%en_sca    
     
     if (nside /= nside_old .or. ipix /= ipix_old) then  ! this should minimize the re-calculation of the scattering phase function. however, this subroutine is still not very efficient because most probably the other cells hit by the same ray are not going to be processed by the same thread. Another possibility is to assign the en_sca_list element depending on ipix instead. However, in this way the assignments to scaspe_arr and scaspe_tot_arr are not safe because the same cc element might be updated by multiple threads.         
        iq_a = iq_sca_node
        call calc_ffn_arr(al,dl,ipix)
     else 
        if (.not.iq_a(il)) then  !!! this is probably useless now  
           print *, 'here you should never get'
           stop
           iq_a(il) = .TRUE.    
           call calc_ffn_arr(al,dl,ipix)
        endif
     endif

     k = (il - id_mpi)/np_mpi  ! wavelength index 
     if (lock_cell(cc) == 1) then
        
        scaspe_arr(k)%a(:,cc)=scaspe_arr(k)%a(:,cc) + en_sca_out*ffn_arr_mw(k)%a
        if (allocated(scaspe_tot_arr) .and. (.not. sequential_scattering)) then 
           scaspe_tot_arr(k)%a(:,cc)=scaspe_tot_arr(k)%a(:,cc) + en_sca_out*ffn_arr_mw(k)%a
        endif

     else 
        
        do ip = 0, npix_arr(il+i0l) -1
           !$OMP ATOMIC 
           scaspe_arr(k)%a(ip,cc)=scaspe_arr(k)%a(ip,cc) + en_sca_out*ffn_arr_mw(k)%a(ip)   
        enddo
        if (allocated(scaspe_tot_arr) .and. (.not. sequential_scattering)) then 
           do ip = 0, npix_arr(il+i0l) -1
              !$OMP ATOMIC 
              scaspe_tot_arr(k)%a(ip,cc)=scaspe_tot_arr(k)%a(ip,cc)+ en_sca_out*ffn_arr_mw(k)%a(ip)
           enddo
        endif

     endif

     nside_old = nside 
     ipix_old = ipix 
     
  end do

  call remove_lock_from_cell(cc_old)

! restore previous values 
  nside = nside_prev 
  iq_a = iq_a_prev 
  ffn_arr_mw = ffn_arr_mw_prev

end subroutine process_en_sca_received


!> Prepares en_sca_list_all() and determines how many elements to send to each MPI process and starting indeces.  
subroutine prepare_en_sca_transfer
 integer :: i, i0,i1, il, im, iel    
 integer :: ierr
 integer :: ibuff

 ! number of elements to be sent out. 

 el_out_mpi = sum(el_out_arr, 1)  

 ! determine position first element in the en_sca_list to be sent to each MPI process 

 i0_el_out_mpi = 0

 do i = 1, np_mpi-1 
    
    i0_el_out_mpi(i) = i0_el_out_mpi(i-1) + el_out_mpi(i-1) 
    
 enddo
 
 ! communicate number of elements that will be received. If a process has no more elements to send out and it is in the waiting part at the end of the rt_loop, it sends out el_out_mpi = -1. Then each MPI process understands that that process in that state.  

 if (count_en_sca_tot == 0 .and. all(wait_thread_arr)) el_out_mpi = -1 

 do im = 0, np_mpi-1 
    
    call MPI_SCATTER(el_out_mpi, 1, MPI_INTEGER, ibuff, 1, MPI_INTEGER, im, MPI_COMM_WORLD, ierr)
    el_in_mpi(im) = ibuff

 end do

 if (count_en_sca_tot == 0) el_out_mpi = 0 

 mpi_proc_completed = .FALSE.

 where(el_in_mpi == -1) 
    el_in_mpi = 0
    mpi_proc_completed = .TRUE.
 end where

! determine position first element in the en_sca_list_received to be received by each MPI process. 

 i0_el_in_mpi = 0

 do i = 1, np_mpi-1 

    i0_el_in_mpi(i) = i0_el_in_mpi(i-1) + el_in_mpi(i-1) 

 enddo
 
end subroutine prepare_en_sca_transfer

!> Transfers en_sca_list elements between the MPI processes and process them. The transfering is done by the master threads while the other threads process the elements as soon as they arrive. This method improves efficiency because communication and computation are partly overlapped.   
subroutine en_sca_transfer() 
  integer :: i, j, im,k 
  integer :: ierr 
  integer :: i0, i1, n_el_out, n_el_in
  integer :: im_s,im_r
  type (list_en_sca), allocatable :: en_sca_buff(:)

  if (ThreadID == 0) then  ! Master thread
  
     n_el_in_tot = sum(el_in_mpi) ! total number of received en_sca elements 
     if (n_el_in_tot > 0) then 
        if (size(en_sca_list_received) /= n_el_in_tot) then 
           deallocate(en_sca_list_received, en_sca_id_thread)
           allocate(en_sca_list_received(0:n_el_in_tot-1), en_sca_id_thread(0:n_el_in_tot-1))
        endif
     endif
  
     do im = 0, np_mpi-1 ! ID of root process in mpi_scatterv
  
        i0 = 0 
        i1 = el_in_mpi(im)-1
        n_el_in = el_in_mpi(im)
        
        if (n_el_in > 0) then 
           allocate(en_sca_buff(i0:i1))
        else 
           allocate(en_sca_buff(0:0)) ! this is just to avoid possible error when calling mpi_scatterv  
        endif


        call mpi_scatterv(en_sca_list_all, el_out_mpi, i0_el_out_mpi, en_sca_arrtype, en_sca_buff,el_in_mpi(im), en_sca_arrtype, im,MPI_COMM_WORLD,ierr)


        if (n_el_in > 0) then 
           i0 = i0_el_in_mpi(im) 
           i1 = i0 + n_el_in-1
           en_sca_list_received(i0:i1) = en_sca_buff
           
           do j = i0, i1
              ! elements are distributed among the threads for processing. the "nproc-1" and the "+1" are there because only the non master threads are used to process the received elements, while the master processing keep communicating. 
             !en_sca_id_thread(j) = mod(en_sca_list_received(j)%cc,nproc-1)+1             
              en_sca_id_thread(j) = mod(en_sca_list_received(j)%ipix,nproc-1)+1  !! this is faster but not safe when updating scaspe_arr. Sometimes atomic statement needed (see process_en_sca_received)
                         
           end do
        
           
        endif

        deallocate(en_sca_buff)

        en_sca_confirm(im) = .TRUE.  ! this is used to confirm that the en_sca_list elements from MPI process im have been received. 
     
     enddo

  else  

     do im = 0, np_mpi-1 
        
        ! Important note about the following spin-lock. It seems that this works only if the spin-lock is in the current position. If you move it within the following if statement or the process_en_sca_received subroutine, it doesn't work. What happens is that en_sca_confirm, although a shared variable, is not updated in the local memory for each thread if the spinlock is in the other locations. A way to assure en_sca_confirm is always updated is to insert a OMP FLUSH statement. However, this is very inefficient. In other words, DO NOT MOVE THE SPIN-LOCK anywhere else!!! 
        do  ! wait that elements have been transfered 

           if (en_sca_confirm(im)) exit 
           
        end do
        
         if (el_in_mpi(im) > 0) then 

           ! process elements 
           call process_en_sca_received(im)
        endif
        
     end do

  endif

  !$OMP BARRIER
  !$OMP MASTER

  en_sca_confirm = .FALSE.

  !$OMP END MASTER
  !$OMP BARRIER 
  
        
  ! TEST TO CHECK COMMUNICATION DONE PROPERLY (IS THIS STILL WORKING AFTER THE CHANGES ? ) 

!!$     im_s = 3  ! sending process 
!!$     im_r = 0  ! receiving process 
!!$
!!$  if (id_mpi == im_s ) then 
!!$     write(25,*) 'sender proc'
!!$     i0 = i0_el_out(im_r)
!!$     i1 = i0 + el_out(im_r)-1
!!$     write(25,*) el_out(im_r)
!!$     do i = i0, i1 
!!$        write(25,*) en_sca_list_all(i)
!!$     enddo 
!!$     
!!$  endif
!!$
!!$  if (id_mpi == im_r) then 
!!$     write(26,*) 'received'
!!$     write(26,*) el_in(im_s)
!!$     i0 = i0_el_in(im_s)
!!$     i1 = i0 + el_in(im_s)-1
!!$     do i = i0, i1 
!!$     write(26,*) en_sca_list_received(i)
!!$     enddo 
!!$  endif
!!$
!!$  call mpi_barrier(MPI_COMM_WORLD,ierr)
!!$  stop

end subroutine en_sca_transfer

!> Creates the MPI data type used for transfering the en_sca_list elements between MPI processes. 
subroutine create_mpi_type
  integer, parameter :: nt = 5 ! number of list_en_sca data_type components 
  integer :: blocklen(0:nt-1), type(0:nt-1) 
  integer(KIND=MPI_ADDRESS_KIND) :: disp(0:nt-1), lb, extent,disp0  
  integer :: ierr 
  integer :: i
  type (list_en_sca) :: prova_list(0:1)
  integer :: en_sca_type

  call MPI_GET_ADDRESS(prova_list(0)%cc, disp(0), ierr) 
  call MPI_GET_ADDRESS(prova_list(0)%nside, disp(1), ierr) 
  call MPI_GET_ADDRESS(prova_list(0)%ipix, disp(2), ierr) 
  call MPI_GET_ADDRESS(prova_list(0)%il, disp(3), ierr) 
  call MPI_GET_ADDRESS(prova_list(0)%en_sca, disp(4), ierr) 

  disp0= disp(0)
  do i = 0, nt-1 

     disp(i) = disp(i) - disp0    

  end do

  blocklen = 1 

  type(0) = MPI_INTEGER 
  type(1) = MPI_INTEGER
  type(2) = MPI_INTEGER
  type(3) = MPI_INTEGER
  type(4) = MPI_DOUBLE_PRECISION

  call MPI_TYPE_CREATE_STRUCT(nt, blocklen, disp, type, en_sca_type, ierr) 
  call MPI_TYPE_COMMIT(en_sca_type, ierr) 

  call MPI_GET_ADDRESS(prova_list(0), disp(0), ierr) 
  call MPI_GET_ADDRESS(prova_list(1), disp(1), ierr) 
  extent = disp(1) - disp(0)  
  lb = 0 
 
  call MPI_TYPE_CREATE_RESIZED(en_sca_type, lb, extent, en_sca_arrtype, ierr) 
  call MPI_TYPE_COMMIT(en_sca_arrtype, ierr) 
  
  call mpi_barrier(MPI_COMM_WORLD,ierr)

end subroutine create_mpi_type

!> Assigns values to scaspe_temp_arr() array. It calls the subroutines to send and receive scaspe_arr() arrays values from the other MPI processes. 
subroutine assign_scaspe_temp(cc)
 integer :: cc
 integer :: i, it0,it1
 integer :: ierr 
 integer :: il
 
 cc_list_local(ThreadID) = cc  ! assign local values
 
 !$OMP BARRIER

 !$OMP MASTER

 ! communicate lists of cell IDs to all processes. These are the first cells in the blocks of scaspe_temp elements to be received by each MPI process/ OPENMP thread.  
 call mpi_allgather(cc_list_local, nproc, MPI_INTEGER, cc_list_all, nproc, MPI_INTEGER, MPI_COMM_WORLD, ierr)  

!!$ if (main_prc) print *, 'transfer started'
!!$ if (main_prc) then
!!$    print *, 'new list'
!!$    do i =0, np_mpi -1 
!!$       
!!$       print *, cc_list_all(:,i)
!!$    enddo
!!$
!!$ endif

 if (any(cc_list_all /= -1)) then 
    
    ! prepare sets of scaspe_temp arrays to be sent to each MPI process 
    call prepare_scaspe_temp_transfer

    ! transfer scaspe_temp arrays 
    call scaspe_temp_transfer
 endif 
 
 !$OMP END MASTER 

 !$OMP BARRIER

 if (cc /= -1) then !each thread get its own block of scaspe_temp arrays
    it0 = ThreadID*num_scaspe_pass
    it1 = (ThreadID+1)*num_scaspe_pass-1
    do il = 0, lnum-1 
       scaspe_temp_arr_big(il)%a = scaspe_temp_recv(il)%a(:,it0:it1)
    end do
    iscaspe_big = 0
 endif

   
end subroutine assign_scaspe_temp

!> Prepares the scaspe_temp_send() array whose parts will be transfered to other nodes using mpi_scatterv. To reduce the number of scaspe transfers, a group of scaspe_temp_arr arrays (corresponding to many cells not just those requested) are passed to the MPI process/OPENMP Thread that sent the request. This group corresponds to the scaspe_temp_arr elements of the cells that will be processed sequentially by that MPI process/OPENMP thread within the current OPENMP chunk (see OMP DO SCHEDULE). The order in which elements are assigned to the 1-dim scaspe_temp_send() array is angular directions (npix_arr(kl)), wavelength(kl), cell id (ic)  
subroutine prepare_scaspe_temp_transfer
integer :: ip,im,im0,im1, ith, il, k, ik0, ik1, ic, i0, i1   
integer :: cc 
integer :: ierr 
integer :: i0l, i1l,kl

call set_i_opacity_arrays(i0l,i1l)

if (lnum_node > 0) then   ! proceed only if there are scaspe arrays stored locally

   scaspe_temp_send = 0   

   k = 0  ! this is the index of scaspe blocks (of size num_scaspe_pass) to be extracted 
   do ip = 0, np_mpi-1   ! MPI process ID 

      do ith = 0, nproc -1   ! OPENMP THREAD ID 

         cc = cc_list_all(ith,ip)  ! REQUESTED CELL ID          
         if (cc == -1) cycle    ! no request in this case 
         
         im0 = (cc-ip)/np_mpi  ! this is the index im as in rt_loop
       !  print *, 'im0 =', im0, threadID,id_mpi
        
         im1 = im0 + num_scaspe_pass  
         ik0 = k*num_scaspe_pass  ! Starting index for each scaspe_temp block in scaspe_temp_send(ik0*tot_npix_arr_local(id_mpi)).

         do im = im0, im1-1 
          !  print *, 'im =', im, threadID, id_mpi
            ! the following IF is to avoid loading elements that will not be used because the Thread will terminate its OPENMP chunk iteration and will start a new iteration with a different im0 at this point.  
            ! mod(im,chunk_size) == 0 means im is at the beginning of an OPENMP thread chunk interval
            ! im /= im0 is to exclude the case in which the thread will just start the iteration at that im value.
            if (mod(im,chunk_size) == 0 .and. im /= im0 ) then              
               exit  
            endif
            ic = im*np_mpi+ip  ! this is the cell ID 
            if (ic > tot_ncell -1) exit   ! exit if no more cell left  
            if (rt_algorithm_ID == rta_2D .or. rt_algorithm_ID == rta_dust2d) then                 
               if (skip_cell_2D_rt_loop(ic)) cycle   ! for 2D algorithm, add only cells processed in the two rt_loop_2D do loops.
            endif
            if (.not. sequential_scattering) then 
               i0 = (ik0+im-im0)*tot_npix_arr_local(id_mpi) ! starting element for the assignment to scaspe_temp_send below. It takes into account the number of cells processed in the current loop (im-im0), the starting point corresponding to the first element to be passed to the corresponding MPI process/OPENMP Thread (ik0), and the total number of directions considered in each scaspe_temp_arr() array (tot_npix_arr_local(id_mpi)).
               do il = 0, lnum-1
                  if (iq_sca_node(il)) then
                     call set_wavelength_index(il,kl)                       
                     i1 = i0+npix_arr(il+i0l)                              
                     scaspe_temp_send(i0:i1-1) = scaspe_arr(kl)%a(:,ic)  
                     scaspe_arr(kl)%a(:,ic) = 0
                     i0 = i1 
                  endif
               end do
               
            else     
               i0 = (ik0+im-im0)*tot_npix_arr_local(id_mpi)
               do il = 0, lnum-1
                  if (iq_sca_node(il)) then
                     call set_wavelength_index(il,kl)
                     i1 = i0+npix_arr(il+i0l)
                     scaspe_temp_send(i0:i1-1) = scaspe_prev_arr(kl)%a(:,ic)  
                     scaspe_prev_arr(kl)%a(:,ic) = 0  
                     i0 = i1 
                  endif
               enddo
            endif
            
         enddo
         k = k + 1 
      end do
      
   end do

endif

end subroutine prepare_scaspe_temp_transfer

!> Sends chunks of scaspe_temp_send() to the other processes using mpi_scatterv. It stores the received values in scaspe_temp_recv().
subroutine scaspe_temp_transfer
integer :: im, i,il,k, ik,ik0,ik1 
integer :: el_out(0:np_mpi-1), el_in(0:np_mpi-1)
integer :: i0_el_out(0:np_mpi-1)
integer :: lnum_recv
real(kind=real64), allocatable :: scaspe_buff(:)
integer :: ierr
integer, allocatable :: iq_in(:)
integer :: nproc_in
integer :: nthreads
integer :: i0, i1 
integer :: i0l, i1l

! set number of elements to be sent/received and starting indeces
call set_i_opacity_arrays(i0l,i1l)

el_out = 0 
el_in = 0

do i = 0, lnum -1 
   if (iq_sca_node(i)) then
      call set_wavelength_index(i,k)
      el_out = el_out + npix_arr(i+i0l)*num_scaspe_pass   
      el_in = el_in + npix_arr(i+i0l)*lnum_node_arr*num_scaspe_pass  ! NOTE: lnum_node_arr is an array
   endif

end do 

do i = 0, np_mpi -1

   el_out(i) = el_out(i)*count(cc_list_all(:,i) /= -1)
   el_in(i) = el_in(i)*count(cc_list_all(:,id_mpi) /= -1)
   
end do

i0_el_out = 0

do i = 1, np_mpi-1 

i0_el_out(i) = i0_el_out(i-1) + el_out(i-1) 

enddo

! set index for scaspe_temp_recv

nproc_in = count(cc_list_local /= -1)  ! this is the number of OMP threads that required the assignment of scaspe 
allocate(iq_in(0:nproc_in*num_scaspe_pass-1))

k = 0 
do i = 0, nproc-1 ! loop on OMP Threads

   if (cc_list_local(i) == -1) cycle

   ik0 = k*num_scaspe_pass
   ik1 = (k+1)*num_scaspe_pass
   
   do ik = ik0, ik1-1  
      iq_in(ik) = i*num_scaspe_pass+ik-ik0   ! these are the positions of the scaspe_temp_recv arrays where the values in scaspe_buff will be assigned. Note that these are subscripts corresponding to the ID of the cells whose scaspe_temp_arr values have been requested. They are the same for all the different wavelengths and thus for all the different MPI processes for which scaspe_temp_arr values are received.   
   enddo
   k = k + 1 

end do

! scatter scaspe_temp_send array   
do im = 0, np_mpi-1 ! ID of root process in mpi_scatterv

   if (lnum_node_arr(im) == 0) cycle 

   allocate(scaspe_buff(0:tot_npix_arr_local(im)*nproc_in*num_scaspe_pass-1))
   
   call mpi_scatterv(scaspe_temp_send, el_out, i0_el_out, MPI_DOUBLE_PRECISION, scaspe_buff,el_in(im), MPI_DOUBLE_PRECISION, im,MPI_COMM_WORLD,ierr)

   if (nproc_in > 0) then 
      i0 = 0 
      do ik = 0, nproc_in*num_scaspe_pass-1
         do il = 0, lnum_node_arr(im) -1 
            k = il*np_mpi+im   ! Wavelength index         
            i1 = i0 + npix_arr(k+i0l)
            scaspe_temp_recv(k)%a(:,iq_in(ik)) = scaspe_buff(i0:i1-1)
            i0 = i1
         enddo
      enddo
   endif
   
   deallocate(scaspe_buff)

end do

deallocate(iq_in)

! check on scaspe_temp_recv (see above in prepare_scaspe_temp_transfer). Not sure this test is working after implementing transfers of blocks of scaspe_temp elements. 
!!$if (id_mpi == 2 ) then 
!!$   write(26,*) scaspe_temp_recv(:,1,2)
!!$endif

end subroutine scaspe_temp_transfer


!> Handles MPI transfers depending on the transfer_type input from each OPENMP thread and MPI PROCESS. The reason for this subroutine is that there are two different situations where a thread can ask for communication to the other MPI process: 1) the assignment of scaspe_temp_arr() arrays; 2) the transfer and processing of the en_sca_list() elements.
!> This routine is called when all the OPENMP threads within an MPI process have reached one of the two calling points (depending on the operation needed). When all the MPI processes have arrived to mpi_allreduce, they communicate transfer_type_tot_local(), which is equal to the number of threads only if no thread needs communication of scaspe_temp elements (that is, all threads requested communication of en_sca_list elements). Then the routines assign_scaspe_temp() and process_en_sca_list() are called to perform the communications. 
subroutine handle_mpi_transfers(transfer_type,cc)
  integer*1 :: transfer_type
  integer :: ierr
  integer :: i 
  integer :: cc  ! used only for scaspe processing 
     
  transfer_type_local(ThreadID) = transfer_type
  
  !$OMP BARRIER    
  !$OMP MASTER
  
  transfer_type_tot_local = sum(transfer_type_local)
  
   call mpi_allreduce(transfer_type_tot_local, transfer_type_tot_all, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  !$OMP END MASTER
  !$OMP BARRIER    

   if (transfer_type_tot_all == 0) then  
      call assign_scaspe_temp(cc)
   else if (transfer_type_tot_all == nproc*np_mpi) then 
      call process_en_sca_list
   else
      call assign_scaspe_temp(cc)
      call process_en_sca_list
   endif

   handling_mpi_thread_arr(ThreadID) = .FALSE.

   !$OMP BARRIER
   !$OMP MASTER 

  count_transfer = count_transfer +1 
  
  ! if (main_prc) print *, 'count transfer = ', count_transfer  
   

  !$OMP END MASTER
  !$OMP BARRIER

  ! the following is needed only when size_en_sca_list is let varying during the calculation. At the moment it is fixed. 
  if (size(en_sca_list) /= size_en_sca_list) call create_en_sca_list()
  
  !$OMP BARRIER

end subroutine handle_mpi_transfers

!> It counts the number of cells that have been processed so far. Only those with non zero luminosity are counted. This count can be used when estimating the average processing time for a single cell. 
subroutine count_processed_cells

!$OMP ATOMIC
num_processed_cells = num_processed_cells+1

end subroutine count_processed_cells

!> It checks whether the thread is processing a cell with an im value which is the start of an OPENMP loop chunk. If so, it sets iscaspe_big = num_scaspe_pass so that the query of scaspe_temp_arr_big is made in the rt_loop. Otherwise, it just increases the value of iscaspe_big. This subroutine is not important in the no communication mode (see no_communications()). 
subroutine check_im(im)
integer :: im

if (mod(im,chunk_size) == 0) then
   iscaspe_big = num_scaspe_pass
else 
   iscaspe_big = iscaspe_big + 1
endif 

end subroutine check_im

!> Assigns the values of dens_arr() by scaling the dens() array read from the main grid file. 
subroutine scale_dens_arr
  integer :: i, i0, i1,j 
  real(kind=real64), allocatable :: temp_arr(:)
  
  if (main_prc) print *, 'wavelength scaling of the dust density array...'

  dens_ref=dens_ref/(kext_ref)

  call set_i_opacity_arrays(i0,i1)
  
  if (.not. use_lambda_grid .or. (rt_type == rtt_grid_init_dust) ) then  ! assign dens_arr in this case
   
     do i = 0, lnum-1

        dens_arr(i,:)=dens_ref(:)*kext_arr(i+i0)  ! rescaling opacity with correct lambda 
     
     enddo

  else  ! check that the values in the input files are consistent 

     allocate(temp_arr(0:tot_ncell-1))
     
     do i = 0, lnum-1

        temp_arr = dens_ref(:)*kext_arr(i+i0)

        do j = 0, tot_ncell-1

           if (cchild(j) /= -1 .or. dens_ref(j) == 0) cycle

           if (abs(temp_arr(j)-dens_arr(i,j))/temp_arr(j) > 5E-3) then
              if (main_prc) then
                 print *, 'STOP(scale_dens_arr): the input lambda grid dust density and that calculated by rescaling the reference grid values do not match! Did you use the same dust model when creating the grids ?'
                 print *, 'lambda = ', lambda_arr(i+i0)
                 print *, 'cell ID = ', j
                 print *, 'lambda grid value = ', dens_arr(i,j)
                 print *, 'scaled value = ', temp_arr(j)
              endif
              call stop_prc
           endif

        end do

     end do

     deallocate(temp_arr)
     
  endif

  dens_ref = dens_ref*kext_ref

  call print_done 

end subroutine scale_dens_arr

!> Sets the subscript of the first element of the value segment within kext_arr(), gsca_arr() and ksca_norm_arr() used in the RT calculation. 
subroutine set_i_opacity_arrays(i0,i1)
  integer :: i0,i1
  
  i0 = 0
  i1 = lnum_stars-1
  if ((rt_algorithm_ID == rta_dust .or. rt_algorithm_ID == rta_dust2D .or. rt_algorithm_ID == rta_i_obs_dust .or. rt_algorithm_ID == rta_sed_dust).and. .not. (rt_type == rtt_grid_init_stars)) then  ! the not condition is to allow correct behaviour in call from grid_initialize
     i0 = i_lambda_dust(0)
     i1 = lnum_tot-1
  endif

end subroutine set_i_opacity_arrays

  


!> Increases lock_cell() by a unity. 
subroutine put_lock_to_cell(i)
  integer :: i 
  
  !$OMP ATOMIC 
  lock_cell(i) = lock_cell(i) + 1 

end subroutine put_lock_to_cell

!> Decreases lock_cell() by a unity. 
subroutine remove_lock_from_cell(i)

  integer :: i 
  
  !$OMP ATOMIC 
  lock_cell(i) = lock_cell(i) - 1 

end subroutine remove_lock_from_cell


!> Changes the size of the wavelength dimension of u_fest_arr() to match the dust emission wavelength range. Stores the radiation field energy density in the stellar emission wavelength range into u_final_uv_opt(). New u_final_arr and dens_arr arrays are created to cover the wavelength range of the dust emission. Note that in this subroutine, lnum is still set to lnum_stars.
subroutine store_reshape_arrays

  if (main_prc) print *, 'Reshaping arrays to match dust emission wavelength grid'
  
  ! reallocate u_fest_arr
  deallocate(u_fest_arr)
  allocate(u_fest_arr(0:lnum_dust-1,0:tot_ncell-1))
  u_fest_arr = 0 

  ! create u_final_uv_opt and pass values of u_final_arr
  allocate(u_final_uv_opt(0:lnum_stars-1,0:tot_ncell-1))
  u_final_uv_opt = u_final_arr
  
  ! reshape u_final_arr and reinitialize
  deallocate(u_final_arr)
  allocate(u_final_arr(0:lnum_dust-1,0:tot_ncell-1))
  u_final_arr = 0 
  
  ! reshape dens_arr and reinitialize
  deallocate(dens_arr)
  allocate(dens_arr(0:lnum_dust-1,0:tot_ncell-1))
  dens_arr = 0

  ! reshape dens_stars and reinitialize
  deallocate(dens_stars_arr)
  allocate(dens_stars_arr(0:lnum_dust-1,0:tot_ncell-1), dens_stars_arr_prev(0:lnum_dust-1,0:tot_ncell-1))
  dens_stars_arr = 0
  dens_stars_arr_prev = 0

  ! deallocate arrays if necessary
  if (rt_algorithm /= 'dust' .and. rt_algorithm /= 'dust2D' .and. rt_algorithm /= 'i_obs_dust') then
     deallocate(scaspe_arr)
     if (allocated(scaspe_tot_arr)) deallocate(scaspe_tot_arr)
     if (allocated(scaspe_prev_arr)) deallocate(scaspe_prev_arr)
     if (allocated(i_obs)) deallocate(i_obs, i_obs_arr)
     if (allocated(i_obs_in)) deallocate(i_obs_in, i_obs_in_arr)
     if (allocated(psel_av_arr)) deallocate(psel_av_arr) 
     if (allocated(lumcell)) deallocate(lumcell, tot_rad_en, tot_rad_en_or, lum_lost, lum_lost_prev)
     if (allocated(iq_sca_node)) deallocate(iq_sca_node, iq_sca_id, iq_maps_id, lnum_node_arr, im_lambda_arr)
     if (allocated(map_arr_out)) deallocate(map_arr_out)
     if (allocated(map_in_arr_out)) deallocate(map_in_arr_out)
  endif

  call print_done
     
end subroutine store_reshape_arrays

!> Sets units_i_obs(), units_ufield() and cs(). 
subroutine set_units

  ! set units_ufield
  select case (rt_type)
  case(rtt_grid_init_stars)

     select case (units_luminosity) ! here it assumes csize always in pc
     case('erg/s/Hz')
        units_ufield = 'erg/Hz/pc^3'
     case('W/Hz')
        units_ufield = 'J/Hz/pc^3'
     case default
        if (main_prc) print *, 'STOP(set_units): fix units_ufield type here'
        call stop_prc
     end select

     ! set units_i_obs
     select case (units_luminosity) ! here it assumes csize always in pc
     case('erg/s/Hz')
        units_i_obs = 'erg/s/Hz/pc^2/sr'
     case('W/Hz')
        units_i_obs = 'W/Hz/pc^2/sr'
     case default
        if (main_prc) print *, 'STOP(set_units): fix units_i_obs type here'
        call stop_prc
     end select

     ! set light speed 
     select case (units_csize)
     case('pc')
        cs = cspeed/parsec
     case default
        if (main_prc) print *, 'STOP(set_units): fix cs type here'
        call stop_prc
     end select
        
  case(rtt_grid_init_dust)
     
     units_i_obs = 'W/m/pc^2/sr' ! dust emission always in these units. units_ufield refers only to stellar emission radiation field. 

  case default
     if (main_prc) print *, 'STOP(set_units): rt_algorithm_ID not recognized!'
     call stop_prc
     
  end select
  
end subroutine set_units

!> Adds the i_obs_arr(), calculated durint the current dust heating iteration, to the arrays i_obs_arr_dir() and i_obs_arr_tot(). 
subroutine update_i_obs_arr

  if (rt_algorithm_ID /= rta_dust .and. rt_algorithm_ID /= rta_dust2d) return 
  
  select case(rt_type)
  case(rtt_output_part2)
     if (tot_ndir > 0) then
        i_obs_arr_dir = i_obs_arr_dir + i_obs_arr
     endif
     if (tot_ndir_in > 0) then
        i_obs_in_arr_dir = i_obs_in_arr_dir + i_obs_in_arr
     endif
  case default
     if (tot_ndir > 0) then
        i_obs_arr_tot = i_obs_arr_tot + i_obs_arr
     endif
     if (tot_ndir_in > 0) then
        i_obs_in_arr_tot = i_obs_in_arr_tot + i_obs_in_arr
     endif
  end select

end subroutine update_i_obs_arr

!> Sets the number of pixels for the scaspe() arrays at each wavelength (npix_arr() and npix_hp_arr()). Having a variable size for the scattering source function allows to save memory especially in the infrared. 
subroutine set_npix_arr
  real(kind=real64) :: g, g_th
  real(kind=real64) :: pmax, pmin, cos_theta_HM, fwhm, delta_min
  integer :: i,j,k_min
  integer :: nmin
  logical :: kp_included(0:kp_sca_max) 
  
  if (main_prc) print *, 'set HealPix number for scattering source function...'
  
  nmin = 5  ! minimum number of pixel sampling the 1D HG function within the FWHM

  g_th = 2*1E-3  ! gsca threshold value (below that consider scattering as isotropic)
  
  allocate(npix_arr(0:lnum_tot-1), npix_hp_arr(0:lnum_tot-1), kp_sca_arr(0:lnum_tot-1), ik_sca_arr(0:lnum_tot-1))

  do i = 0, lnum_tot-1

     g = abs(gsca_arr(i))  ! formulae are the same for positive and negative gsca if you take the absolute value         

     if (g < g_th) then  ! for low values assume isotropic scattering
        npix_hp_arr(i) = 1
        kp_sca_arr(i) = -1  ! note that in this case the formula for npix_hp is not valid.  
        cycle
     endif
     
     pmax = (1-g**2)/(1+g**2-2*g)**1.5  ! max and min HG phase function 
     pmin = (1-g**2)/(1+g**2+2*g)**1.5

     cos_theta_HM = ((2*(1-g**2)/(pmax+pmin))**(2./3.)-(1+g**2))/(-2*g) ! this formula can be derived by finding the cos(theta) such that the function p(theta)-p(180 deg) has half its maximum value. The subtraction of the p(180) term is equivalent to a background subtraction necessary to determine the FWHM. 

     fwhm = 2*acos(cos_theta_hm)
     
     delta_min = fwhm/nmin

     k_min = int(1/(2*log10(2.))*log10(4*pi/(12*delta_min**2))) ! minimum HEALPix k parameter needed to have about 5 pixel within the FWHM of the HG phase function. This formula can be derived from delta_min = (4*pi/npix_min)^0.5 with npix_min = 12*2^(2*k_min)

     if (k_min > kp_sca_max) k_min = kp_sca_max ! kp_sca_max is the maximum allowed k parameter 
     kp_sca_arr(i) = k_min
     npix_hp_arr(i) = 12*2**(2*k_min)

     !print *, lambda_arr(i), gsca_arr(i),  k_min, npix_hp_arr(i)   

  enddo  

  npix_arr = npix_hp_arr + tot_ndir_scaspe

  ! assign npix_unique
  
  dim_npix_unique = 0  ! find out how many different values 
  kp_included = .FALSE.
  
  do i = 0, kp_sca_max

     if (any(npix_hp_arr ==  12*2**(2*i))) then
        dim_npix_unique = dim_npix_unique + 1
        kp_included(i) = .TRUE.
     endif

  end do

  allocate(npix_hp_unique(0:dim_npix_unique-1),npix_unique(0:dim_npix_unique-1), kp_unique(0:dim_npix_unique-1))

  ik_sca_arr = -1 

  j = 0 
  do i = 0, kp_sca_max
     if (kp_included(i)) then
        kp_unique(j) = i 
        npix_hp_unique(j) = 12*2**(2*i)
        where (kp_sca_arr == kp_unique(j))
           ik_sca_arr = j
        end where
        j = j +1 
     endif
  end do
  
  npix_unique = npix_hp_unique  + tot_ndir_scaspe

  call print_done 

end subroutine set_npix_arr

!> Assigns i_obs_arr() values for the projection algorithm. 
subroutine assign_i_obs_to_project
  integer :: i,j 

select case(param_to_project)
case('stellar_emission')
   do i = 0, lnum-1 
      ! i_obs_arr
      if (iq_sca_node(i)) then 
         do j = 0, tot_ndir-1 ! same value for all inclinations 
            i_obs_arr(i,0:tot_ncell-1, j) = dens_stars_arr(i,:)*csize(:)/(4*pi) 
         enddo
         if (use_p_src) then ! add point sources if present  
            do j = 0, tot_ndir-1 ! same value for all inclinations 
               i_obs_arr(i,tot_ncell:tot_ncell_p_src-1, j) = lum_p_src_arr(i,:)/csize(cell_src)**2/(4*pi) 
            enddo
         endif
       ! i_obs_in_arr
         do j = 0, tot_ndir_in-1 ! same value for all inclinations 
            i_obs_in_arr(i,0:tot_ncell-1, j) = dens_stars_arr(i,:)*csize(:)/(4*pi) 
         enddo
         if (use_p_src) then ! add point sources if present  
            do j = 0, tot_ndir_in-1 ! same value for all inclinations 
               i_obs_in_arr(i,tot_ncell:tot_ncell_p_src-1, j) = lum_p_src_arr(i,:)/csize(cell_src)**2/(4*pi) 
            enddo
         endif
      endif
   end do

case('optical_depth')
   do i = 0, lnum-1 
      if (iq_sca_node(i)) then 
         do j = 0, tot_ndir-1 ! same value for all inclinations 
            i_obs_arr(i,0:tot_ncell-1, j) = dens_arr(i,:)*csize(:) 
         enddo

         do j = 0, tot_ndir_in-1 ! same value for all inclinations 
            i_obs_in_arr(i,0:tot_ncell-1, j) = dens_arr(i,:)*csize(:) 
         enddo        
      endif
   end do
   

case default 
   if (main_prc) print *, 'ERROR(assign_i_obs_project): param_to_project not recognized!'
   call stop_prc

end select


end subroutine assign_i_obs_to_project



END MODULE rt_routines
