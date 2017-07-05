!> This module contains many of the variables used to define the properties of the 3D grid as well as subroutines for creating the grid, finding neighbour cells and more
MODULE smooth_grid_routines  
  use iso_fortran_env
  IMPLICIT NONE 
  include 'mpif.h'
  !> @param modelsize  Size of the model (that is, size of the root cell)  
  real(kind=real64) :: modelsize
  ! Wavelengths
  !> @param lambda  Wavelength variable
  real(kind=real64) :: lambda
  !> @param lambda_ref Reference wavelength. It is used to scale the dens array at each wavelength. In case use_lambda_grid() is set, the lambda grid corresponding to the reference wavelength should be included in the input.  
  real(kind=real64) :: lambda_ref
  !> @param lamdba_ref_SI Reference wavelength in SI units [m].
  real(kind=real64) :: lambda_ref_SI
  ! Counters
  !> @param max_ncell Maximum number of cells allowed in the grid creation program
  integer (kind=int32) :: max_ncell
  !> @param ncurr Counter of cell number in grid creation program
  integer (kind=int32) :: ncurr
  !> @param nlevel Counter of sudivision levels in grid creation program
  integer (kind=int32) :: nlevel
  !> @param tot_ncell Total number of cells in the grid 
  integer (kind=int32) :: tot_ncell
  !> @param tot_ncell_p_src Total number of cells in the grid plus number of point sources 
  integer (kind=int32) :: tot_ncell_p_src
  !> @param tot_ndir Total number of line-of-sight direction specified in the input file_dir_out
  integer (kind=int32) :: tot_ndir
  !> @param tot_ndir_in Total number of observer position within the RT model specified in the input file_pos_obs
  integer (kind=int32) :: tot_ndir_in
  !> @param tot_p_src Total numer of point sources specified in the input file_p_src
  integer (kind=int32) :: tot_p_src
  !> @param max_lvl Maximum number of cell subdivisions
  integer (kind=int32) :: max_lvl
  !> @param min_lvl Minimum number of cell subdivisions
  integer (kind=int32) :: min_lvl
  !> @param tot_spare_cells This variable is equal to the difference between the total number of leaf cells in the original grid and the total number of leaf cells in the grid after the reduction in resolution has been applied. See reduce_grid_res().
  integer (kind=int32) :: tot_spare_cells
  !> @param tot_ndir_scaspe Total number of line-of-sight direction elements included in the scaspe array
  integer (kind=int32) :: tot_ndir_scaspe
  ! 3D Grid arrays
  !> @param ncell(0 : tot_ncell() -1) Array of cell ID numbers
  integer (kind=int32) , allocatable :: ncell(:)
  !> @param cchild(0 : tot_ncell()-1) For each element "i", this array gives the ID number of the first cell in the set of child cells derived from the cell "i". If the cell "i" is a leaf cell, cchild(i) = -1.  
  integer (kind=int32) , allocatable :: cchild(:)
  !> @param lvl(0 : tot_ncell() -1) This array contains the subdivision level for each cell.
  integer (kind=int32) , allocatable :: lvl(:)
  !> @param nstart(0 : max_lvl()+2) This array counts the number of cells for each subdivision loop in the grid creation program.
  integer (kind=int32) , allocatable :: nstart(:)
  !> @param cchild_or(0 : tot_ncell()-1) This is the cchild array of the grid before the reduction of grid resolution is applied.
  integer (kind=int32) , allocatable :: cchild_or(:)
  !> @param cindex(0 : tot_ncell()-1) This array contains a binary code specifying the position of the cell in the cell tree
  integer(kind=int64), allocatable, target  :: cindex(:)
  !> @param ccoord(3,0 : tot_ncell()-1) Two dimensional array containing the 3D position (x,y,z) of the centre for each cell.  
  real(kind=real64), allocatable  :: ccoord(:,:)
  !> @param csize(0 : tot_ncell()-1) This array contains the size of each cell.
  real(kind=real64), allocatable  :: csize(:)
  !> @param dens(0 : tot_ncell()-1) This array contains the extinction coefficient for each cell (in units of length ^ -1). 
  real(kind=real64), allocatable  :: dens(:)
  !> @param dens_ref(0 : tot_ncell()-1) This array contains the extinction coefficient for each cell (in units of length ^ -1) at a reference wavelength. WARNING: Just after reading the main grid (in read_main_grid()), dens_ref is equal to the dens() array in grid_file(). However, the values inside this array can be varied in the user defined routines if needed. 
  real(kind=real64), allocatable  :: dens_ref(:)
  !> @param dens_arr(0 : lnum() -1, 0 : tot_ncell()-1) This array contains the extinction coefficient for each cell (in units of length ^ -1) and for all wavelengths. 
  real(kind=real64), allocatable  :: dens_arr(:,:)
  !> @param dens_stars(0 : tot_ncell()-1) This array contains the radiation source volume emissivity (in units of luminosity/volume) for each grid cell at a single wavelength.
  real(kind=real64), allocatable  :: dens_stars(:)
  !> @param dens_stars_ref(0 : tot_ncell()-1) This array contains the radiation source volume emissivity (in units of luminosity/volume) for each grid cell at the reference wavelength.
  real(kind=real64), allocatable  :: dens_stars_ref(:)
  !> @param dens_stars_arr(0: lnum() -1, 0 : tot_ncell()-1) This array contains the radiation source volume emissivity (in units of luminosity/volume) for each grid cell at all wavelengths. In the case of the dust RT algorithm, it contains only the luminosity still to be propagated through the RT model. 
  real(kind=real64), allocatable  :: dens_stars_arr(:,:)
  !> @param dens_stars_arr_prev(0: lnum() -1, 0 : tot_ncell()-1) Dust emission volume emissivity. It is used to store the total dust emissivity derived in the last dust heating iteration. When the emissivity does not differ more than a factor conv_en_lim() from the newly calculated dust emissivity, the dust heating iterations are stopped.   
  real(kind=real64), allocatable  :: dens_stars_arr_prev(:,:)

  ! Grid Base arrays
  !> @param base(2) This array contains only two elements. The first is the subdivision factor for the root cell and the second is the subdivision factor for all subsequent subdivisions. 
  integer :: base(2)
  !> @param basediv(2) This array contains two elements equal to the minimum power of two which is bigger than the corresponding values in base(). It is used to transform the binary code in cindex() in the explicit tree coordinate of a cell.  
  integer (kind = int64) :: basediv(2)
  !> @param basemask(2) This is equal to basediv()-1. It is used in the same routines as basediv. 
  integer :: basemask(2)
  !> @param baseinv This is equal to base()^-1. It is used when deriving the cell sizes. 
  real(kind=real64) :: baseinv(2)
  ! Max cell parameters
  !> @param max_dtau Maximum optical depth allowed for the leaf cells in the grid creation program.
  real(kind=real64) ::  max_dtau
  !> @param max_dlum Maximum luminosity allowed for the leaf cells relative to the total grid luminosity (which is calculated in user defined routines).
  real(kind=real64) :: max_dlum
  ! RT calculation arrays
  !> @param psel_av_arr(0: iterations_dustem-1, 0 : iterations() , 0 : tot_ncell() + tot_p_src()-1) Array containing the average ray path for each cell and for each direct light/scattered light iteration. Note that the array is expanded in expand_psel_av_arr(). 
  REAL(KIND=real64), allocatable :: psel_av_arr(:,:,:)
  !> @param lum_p_src(0 : tot_p_src()-1) Array containing the luminosities of each point source at a single wavelength 
  REAL(KIND=real64), allocatable :: lum_p_src(:)
  !> @param lum_p_src_ref(0 : tot_p_src()-1) Array containing the luminosities of each point source at the reference wavelength 
  REAL(KIND=real64), allocatable :: lum_p_src_ref(:)
  !> @param lum_p_src_arr(0 : lnum() -1 , 0 : tot_p_src()-1) Array containing the luminosities of each point source for each wavelength. 
  REAL(KIND=real64), allocatable :: lum_p_src_arr(:, :)
  !> @param u_fest(0 : tot_ncell()-1) This array contains the lower limit of the radiation field energy density for each cell and for a single wavelength.
  real(KIND=real64), allocatable :: u_fest(:)
  !> @param u_fest_arr(0 : lnum() -1, 0 : tot_ncell()-1) This array contains the lower limit of the radiation field energy density for each cell and for each wavelength.
  real(KIND=real64), allocatable :: u_fest_arr(:, :)  
  !> @param u_final(0 : tot_ncell()-1) This array contains the radiation field energy density for each cell and for a single wavelength.  
  REAL(KIND=real64), allocatable :: u_final(:)
  !> @param u_final_arr(0 : lnum() -1 , 0 : tot_ncell()-1) This array contains the radiation field energy density for each cell and for each wavelength.  
  REAL(KIND=real64), allocatable :: u_final_arr(:, :)
  !> @param u_final_uv_opt(0 : lnum_stars() -1 , 0 : tot_ncell()-1) This array contains the radiation field energy density for each cell and for each wavelengthproduced by the stellar emission. It is used in the dust emission RT calculation to derived the energy absorbed by the dust. Note that this array is transformed into the radiation field average intensity in calc_conv_ufield_ifield().  
  REAL(KIND=real64), allocatable :: u_final_uv_opt(:, :)

  !> Data type used to create variable length arrays. 1D version. 
  type var_arr_1d
     real(kind=real64), allocatable :: a(:)
     logical, allocatable :: l(:)
     integer, allocatable :: int(:)
  end type var_arr_1d

  !> Data type used to create variable length arrays. 2D version. 
  type var_arr_2d
     real(kind=real64), allocatable :: a(:,:)
  end type var_arr_2d

  !> @param scaspe_arr(0:lnum()-1 ). Array containing the scattered luminosity (for each cell and wavelength) to be processed, propagating in a set of npix_arr directions. This is a variable size array (see data type var_arr_2d() ). The dimension of each block depends on the wavelength through npix_arr().   
  type(var_arr_2d), allocatable :: scaspe_arr(:)

  !> @param scaspe(0: npix()-1, 0: tot_ncell()-1 ) with npix = max (npix_arr()), with npix_arr including only the values for the wavelength used in the RT calculation (it differs for the stellar and dust RT). Array containing the monochromatic scattered luminosity for each cell to be processed, propagating in a set of npix directions. 
  real(KIND=real64), allocatable :: scaspe(:, :)
  
  !> @param scaspe_tot_arr(0:lnum()-1)%a(0:npix_arr(i)-1,0:tot_ncell-1). Array containing the total scattered radiation luminosity, for each cell and wavelength, propagating in a set of directions. This array can be used to recover the scattered luminosity source function and make maps at arbitrary directions, including those not initially considered in the main RT calculation.
  type(var_arr_2d), allocatable :: scaspe_tot_arr(:)
  
  !> @param scaspe_tot(0 : npix - 1, 0: tot_ncell()-1) with npix=12* nside_sca()^2+ tot_ndir_scaspe(). Array containing the total scattered radiation luminosity, for each cell and at a single wavelength, propagating in a set of directions. This array can be used to recover the scattered luminosity source function and make maps at arbitrary directions, including those not initially considered in the main RT calculation. 
  real(KIND=real64), allocatable :: scaspe_tot(:,:)
  
  !> @param scaspe_prev(0:npix-1,0:tot_ncell-1) with npix=npix_arr(i). It stores the values of scaspe() before each scattering iteration. Used for 2D RT algorithm and sequential scattering mode.   
  real (kind=real64), allocatable :: scaspe_prev(:,:)
  
  !> @param scaspe_prev_arr(0:lnum-1). It stores the values of scaspe_arr() before each scattering iteration for all wavelengths. Used for 2D RT algorithm and sequential scattering mode. 
  type(var_arr_2d), allocatable :: scaspe_prev_arr(:) 
  
  !> @param i_obs(0:tot_sources-1,0: tot_ndir()-1) with tot_sources = tot_ncell()+ tot_p_src(). Array containing, for each cell and point source, the monochromatic specific intensity received by the observers located outside and far-away from the model along the specified line-of-sight directions in file_dir_out().
  REAL(KIND=real64), allocatable :: i_obs(:,:)

  !> @param i_obs_arr(0: lnum()-1, 0:tot_sources-1,0: tot_ndir()-1) with tot_sources = tot_ncell()+ tot_p_src(). Array containing, for each cell, wavelength and point source, the specific intensity received by the observers located outside and far-away from the model along the specified line-of-sight directions in file_dir_out().
  REAL(KIND=real64), allocatable :: i_obs_arr(:,:,:)

  !> @param i_obs_arr_dir(0: lnum()-1, 0:tot_sources-1,0: tot_ndir()-1) with tot_sources = tot_ncell()+ tot_p_src(). Same as i_obs_arr() but containing only the direct light from the sources. Used in the dust RT algorithms to distinguish the direct and scattered light in the successive dust heating iterations.
  REAL(KIND=real64), allocatable :: i_obs_arr_dir(:,:,:)

  !> @param i_obs_arr_tot(0: lnum()-1, 0:tot_sources-1,0: tot_ndir()-1) with tot_sources = tot_ncell()+ tot_p_src(). Same as i_obs_arr() but used in the dust RT algorithms to accumulate the values of i_obs_arr() from the different dust heating iterations
  REAL(KIND=real64), allocatable :: i_obs_arr_tot(:,:,:)

  !> @param dir_i_out(0: tot_ndir()-1,2) Array containing the theta and phi angles specifying the line-of-sight directions towards the observers. 
  REAL(KIND=real64), allocatable :: dir_i_out(:,:)

  !> @param ccoord_obs(3,0: tot_ndir_in()-1) Array containing the internal observer positions. 
  REAL(KIND=real64), allocatable :: ccoord_obs(:,:)

  !> @param i_obs_in(0:tot_sources-1,0: tot_ndir_in()-1) with tot_sources = tot_ncell() + tot_p_src(). Array containing, for each cell and point source, the monochromatic specific intensity received by the internal observers. 
  REAL(KIND=real64), allocatable :: i_obs_in(:, :)

  !> @param i_obs_in_arr(0 : lnum() -1, 0:tot_sources-1,0: tot_ndir_in()-1) with tot_sources = tot_ncell() + tot_p_src(). Array containing, for each cell, wavelength and point source, the specific intensity received by the internal observers. 
  REAL(KIND=real64), allocatable :: i_obs_in_arr(:,:, :)

  !> @param i_obs_in_arr_dir(0 : lnum() -1, 0:tot_sources-1,0: tot_ndir_in()-1) with tot_sources = tot_ncell() + tot_p_src(). As i_obs_arr_dir() but for the internal observers. 
  REAL(KIND=real64), allocatable :: i_obs_in_arr_dir(:,:, :)

  !> @param i_obs_in_arr_tot(0 : lnum() -1, 0:tot_sources-1,0: tot_ndir_in()-1) with tot_sources = tot_ncell() + tot_p_src(). As i_obs_arr_tot() but for the internal observers. 
  REAL(KIND=real64), allocatable :: i_obs_in_arr_tot(:,:, :)

  !> @param ccoord_p_src(3, 0 : tot_p_src()-1) Array containing the positions of the point sources. 
  REAL(KIND=real64), allocatable :: ccoord_p_src(:,:)

  !> @param i_obs_dust(nl,0 : tot_ncell()-1) with nl=size[ ind_i_obs()]. Array containing the dust emission specific intensity reaching the observer. Dust emission is considered optically thin (for the moment), so the specific intensity is not direction dependent. 
  real(KIND=real64), allocatable :: i_obs_dust(:,:)

  !> @param csize_arr(0 : max_lvl()) Array of precalculated cell sizes depending on cell subdivision level.
  real(KIND=real64), allocatable :: csize_arr(:)

  !> @param carea_arr(0 : max_lvl()) Array of precalculated cell areas depending on cell subdivision level.
  REAL(KIND=real64), allocatable :: carea_arr(:)

  !> @param cvol_arr(0 : max_lvl()) Array of precalculated cell volumes depending on cell subdivision level.
  REAL(KIND=real64), allocatable :: cvol_arr(:)

  !> @param src_cell(0 : tot_ncell()-1) Array containing the ID number of the point source hosted in each cell. If no sources are hosted by a cell, src_cell = -1. 
  integer (kind=int32), allocatable :: src_cell(:)

  !> @param cell_src(0 : tot_p_src()-1) Array containing the ID number of the cells hosting each point source.  
  integer (kind=int32), allocatable :: cell_src(:)

  !> @param lumcell(0 : lnum() -1 , 0 : tot_ncell()-1) Array containing the total luminosity of each cell and wavelength.
  real(kind=real64), allocatable :: lumcell(:, :)

  !> @param lock_cell(0 : tot_ncell()-1) Array used to lock cells so there is no race condition happening in process_scatt_temp()
  integer (kind=int32), allocatable :: lock_cell(:)

  ! Character variable length
  !> @param lcar Parameter equal to the maximum number of characters for string variables.
  integer, parameter :: lcar =200
  ! Opacity parameters and arrays
  !> @param kabs Dust absorption coefficient at a single wavelength. 
  real (Kind=real64) :: kabs
  !> @param ksca Dust scattering coefficient at a single wavelength. 
  real (Kind=real64) :: ksca
  !> @param kext Dust extinction coefficient at a single wavelength.
  real (Kind=real64) :: kext
  !> @param gsca Parameter of the Heney-Greeinstein function at a single wavelength.
  real (Kind=real64) :: gsca
  !> @param kext_ref Dust extinction coefficient at the reference wavelength.
  real (Kind=real64) :: kext_ref
  !> @param tau_nh Optical depth per unit Hydrogen column density (typically in units of tau*H/cm^2). Used when calculating dust emission. 
  real (Kind=real64) :: tau_nh
  !> @param tau_nh_ref As tau_nh but at a reference wavelength
  real (Kind=real64) :: tau_nh_ref
  !> @param tot_n_dust Total number of dust grains in the assumed dust model. Needed to convert e.g. kext() from cross_section per grain to optical depth per hydrogen surface density (as e.g. tau_nh_ref()).
  real(kind=real64) :: tot_n_dust  
  !> @param kabs_arr(0 : lnum()-1) Array containing kabs() at all wavelengths. During the calculation this is in the form of cross section per grain ([m^2]). 
  real (kind=real64),allocatable :: kabs_arr(:)
  !> @param kext_arr(0 : lnum()-1) Array containing kext() at all wavelengths.  During the calculation this is in the form of cross section per grain ([m^2]). 
  real (kind=real64),allocatable :: kext_arr(:)
  !> @param ksca_arr(0 : lnum()-1) Array containing ksca() at all wavelengths.  During the calculation this is in the form of cross section per grain ([m^2]). 
  real (kind=real64),allocatable :: ksca_arr(:)
  !> @param gsca_arr(0 : lnum()-1) Array containing gsca() at all wavelengths.
  real (kind=real64),allocatable :: gsca_arr(:)
  !> @param ksca_arr_norm(0 : lnum()-1) Array containing normalized ksca() (ksca()/kext, that is, the albedo) at all wavelengths. 
  real (kind=real64),allocatable :: ksca_arr_norm(:)
  !> @param qabs_arr_in(0 : n_dust_comp()-1,0:n_dust_maxsize_qabs-1,0: lnum_tot()-1) Array of Q_abs values for n_dust_comp dust() chemical species, n_dust_maxsize_qabs() dust sizes and lnum_tot() wavelengths. Note that these values are not interpolated at the same grain sizes of the input grain size distribution.  
  real(kind=real64), allocatable :: qabs_arr_in(:,:,:)
  !> @param qabs_ref_in Same as qabs_arr_in() but for reference wavelength.
  real(kind=real64), allocatable :: qabs_ref_in(:,:)
  !> @param qsca_arr_in(0 : n_dust_comp()-1,0:n_dust_maxsize_qabs-1,0: lnum_tot()-1) Same as qabs_arr() but for Qsca. 
  real(kind=real64), allocatable :: qsca_arr_in(:,:,:)
  !> @param qsca_ref_in Same as qsca_arr_in() but for reference wavelength.
  real(kind=real64), allocatable :: qsca_ref_in(:,:)
  !> @param qext_arr_in(0 : n_dust_comp()-1,0:n_dust_maxsize_qabs-1,0: lnum_tot()-1) Same as qabs_arr() but for Qext. 
  real(kind=real64), allocatable :: qext_arr_in(:,:,:)
  !> @param qext_ref_in Same as qext_arr_in() but for reference wavelength.
  real(kind=real64), allocatable :: qext_ref_in(:,:)
  !> @param gsca_arr_in(0 : n_dust_comp()-1,0:n_dust_maxsize_qabs-1,0: lnum_tot()-1) Same as qabs_arr() but for gsca. 
  real(kind=real64), allocatable :: gsca_arr_in(:,:,:)
  !> @param gsca_ref_in Same as gsca_arr_in() but for reference wavelength.
  real(kind=real64), allocatable :: gsca_ref_in(:,:)
  !> @param dust_size_qabs(0 : n_dust_comp()-1,0:n_dust_maxsize_qabs-1) Array containing grain sizes for each value of qabs_arr_in().
  real(kind=real64), allocatable :: dust_size_qabs(:,:)
  !> @param dust_size_fa(0 : n_dust_comp()-1,0:n_dust_maxsize_fa-1) Array containing the grain sizes for the grain size distribution in dust_fa().
  real(kind=real64), allocatable :: dust_size_fa(:,:)
  !> @param dust_fa(0 : n_dust_comp()-1,0:n_dust_maxsize_fa-1) Array containing the grain size distribution for n_dust_comp() chemical species and up to n_dust_maxsize_fa() sizes (different grain species might have grain size sampling). 
  real(kind=real64), allocatable :: dust_fa(:,:)
  !> @param qabs_arr_fa(0 : n_dust_comp()-1,0:n_dust_maxsize_fa-1,0: lnum()-1) Array containing the Qabs values interpolated to the same grain sizes of the size distribution in dust_size_fa(). Although the array second dimension is n_dust_maxsize_fa, for each grain species "i" only n_dust_size(i) grain sizes are included in this array. 
  real(kind=real64), allocatable :: qabs_arr_fa(:,:,:)
   !> @param qabs_arr_planck(0 : n_dust_comp()-1,0:n_dust_maxsize_fa-1,0: n_temp()-1) Array containing the planck averaged Qabs values interpolated to the same grain sizes of the size distribution in dust_size_fa(). Although the array second dimension is n_dust_maxsize_fa, for each grain species "i" only n_dust_size(i) grain sizes are included in this array. 
  real(kind=real64), allocatable :: qabs_arr_planck(:,:,:)
  !> @param n_temp_planck Number of temperatures considered in the calculation of the Planck averaged Qabs (see qabs_arr_planck()) 
  integer, parameter :: n_temp_planck = 100
  !> @param T_arr_planck Temperature array used in the calculation of qabs_arr_planck()
  real(kind=real64),allocatable :: T_arr_planck(:)
  !> @param qabs_ref_fa Same as qabs_arr_fa() but for reference wavelength.
  real(kind=real64), allocatable :: qabs_ref_fa(:,:)
  !> @param qsca_arr_fa Same as qabs_arr_fa() but for Qsca. 
  real(kind=real64), allocatable :: qsca_arr_fa(:,:,:)
  !> @param qsca_ref_fa Same as qsca_arr_fa() but for reference wavelength.
  real(kind=real64), allocatable :: qsca_ref_fa(:,:)
  !> @param qext_arr_fa Same as qabs_arr_fa() but for Qext. 
  real(kind=real64), allocatable :: qext_arr_fa(:,:,:)
  !> @param qext_ref_fa Same as qext_arr_fa() but for reference wavelength.
  real(kind=real64), allocatable :: qext_ref_fa(:,:)
  !> @param gsca_arr_fa Same as qabs_arr_fa() but for gsca. 
  real(kind=real64), allocatable :: gsca_arr_fa(:,:,:)
  !> @param gsca_ref_fa Same as gsca_arr_fa() but for reference wavelength.
  real(kind=real64), allocatable :: gsca_ref_fa(:,:)
  !> @param delta_dust_size(0 : n_dust_comp()-1, 0:n_dust_maxsize_fa-1) Array containing the bin widths for integrations over the size distribution. Values for n_dust_comp chemical species and n_dust_maxsize_fa size distribution values.
  real(kind=real64), allocatable :: delta_dust_size(:,:)
  !> @param n_dust_comp Number of dust chemical species in the assumed dust model.
  integer :: n_dust_comp
  !> @param n_dust_wave_qabs Number of wavelengths in the input \f$Q_\lambda \f$, \f$g_\lambda\f$ tables for the assumed dust_model() and dust_opacity_tables().
  integer :: n_dust_wave_qabs
  !> @param n_dust_maxsize_qabs Maximum number of sizes in the input Qabs array for the assumed dust model. The actual number can vary depending on the grain chemical composition. 
  integer :: n_dust_maxsize_qabs
  !> @param n_dust_maxsize_fa Maximum number of grain size distribution values in the input size distribution array for the assumed dust model. The actual number can vary depending on the grain chemical composition. 
  integer :: n_dust_maxsize_fa
  !> @param n_dust_size(0: n_dust_comp()-1) Array containing the number of grain sizes in the input size distribution arrays. Values for n_dust_comp chemical compositions. 
  integer , allocatable :: n_dust_size(:)

  !> @param dust_model Name of dust model to be used. Choices: 'TRUST' for TRUST benchmark project, 'DraineLi06' for dust model of Draine & Li (2006), 'user' for user defined chemical composition and size distributions. For the standard choices ('TRUST', 'DraineLi06') there is no need to specify dust_opacity_model(). 
  character (len=lcar) :: dust_model

  !> @param dust_opacity_tables \f$Q_\lambda\f$ opacity and \f$g_\lambda\f$ tables to be used. If a standard dust_model() is selected ('TRUST', 'DraineLi06'), there is no need to input this variable. If dust_model() = 'user', this variable is required. Choices are: 'DraineLi06' (same grain opacities as in the Draine \& Li 2007 model); 'TRUST' (same opacities as in the TRUST benchmark project); 'user' (user-provided tables). In the latter case, specify the input tables in file_q_gra(), file_q_sil(), file_q_pah_neu() and file_q_pah_ion() (at least one table to be input). Also specify the corresponding size distributions in file_gra_fa(), file_sil_fa(), file_pah_neu_fa() and file_pah_ion_fa(). 
  character (len=lcar) :: dust_opacity_tables
  
  !> @param max_n_dust_comp Maximum number of dust species allowed in the input (corresponding to Silicates, Graphite, Neutral PAH and ionized PAH).
  integer, parameter :: max_n_dust_comp = 4
  !> @param iq_dust_model Elements equal to TRUE when the size distributions for the following grain species are used: Graphite, Silicates, neutral PAH molecules and ionized PAH molecules. Values set in check_input(). At least one element has to be TRUE for user defined grain size distribution. 
  logical :: iq_dust_model(0:max_n_dust_comp-1)
  !> @param n_dust_size_qabs(0: max_n_dust_comp()-1) Array containing the number of grain sizes in the input \f$Q_\lambda\f$, \f$g_\lambda\f$ tables. The values correspond to the following grain compositions (in order): Graphite, Silicates, neutral PAHs and ionized PAHs.
  integer :: n_dust_size_qabs(0:max_n_dust_comp-1)
  !> @param max_n_dust_cal_type Maximum number of dust types for the input specific enthalphy files. At the moment there are only two types that can be specified: Graphite dust (including PAHs) and Silicates dust. See file_calorimetry_Gra(), file_calorimetry_Sil() and load_cT_hT_tables().  
  integer, parameter :: max_n_dust_cal_type = 2
  !> @param grain_density_arr Grain density in g/cm^3 in the input files for the grain specific enthalpy/ specific capacity (see file_calorimetry_Gra() and file_calorimetry_Sil()). 
  real(kind=real64) :: grain_density_arr(0:max_n_dust_cal_type-1)
  !> @param n_dust_temp_cal Number of dust temperatures in the input files for the grain specific enthalpy/ specific capacity (see file_calorimetry_Gra() and file_calorimetry_Sil()). 
  integer :: n_dust_temp_cal(0:max_n_dust_cal_type-1)
  !> @param cal_temp Temperature read from the input tables for the specific enthalpy / heat capacity.
  real(kind=real64), allocatable :: cal_temp(:,:)
  !> @param grain_enthalpy Specific enthalpy values read in the input files (either standard TRUST Calorimetry files or user provided).
  real(kind=real64), allocatable :: grain_enthalpy(:,:)
  !> @param grain_heat_capacity Specific heat capacity values read in the input files (either standard TRUST Calorimetry files or user provided).
  real(kind=real64), allocatable :: grain_heat_capacity(:,:)
  !> @param iq_ct_table Indices associating each grain species with the corresponding specific enthalpy/ heat capacity table.
  integer, allocatable :: iq_ct_table(:)
  !> @param n_int_rf_bins Number of integrated UV (or optical) energy bins to be used in the SED adaptive approach for the calculation of the stochastically heated dust emission. Note that the total number of bins is equal to the square of this input variable.   
  integer :: n_int_rf_bins

  ! Pi parameters
  !> @param halfpi &pi;/2
  REAL(KIND=real64),PARAMETER :: halfpi=asin(1.)
  !> @param pi &pi;
  REAL(KIND=real64),PARAMETER :: pi=2*halfpi
  !> @param twopi 2&pi;
  REAL(KIND=real64),PARAMETER :: twopi= 2*pi
  !> @param twothird 2./3.
  REAL(KIND=real64),PARAMETER :: twothird=2./3.
  ! Wavelength grid variables
  !> @param file_lambda_list Name of the file containing the wavelength list
  character(LEN=lcar) :: file_lambda_list
  !> @param lambda_arr(0 : lnum()-1) Array containing the wavelength list read from file_lambda_list()
  real(Kind=real64), allocatable :: lambda_arr(:)
  !> @param lambda_arr_maps(0 : lnum()-1) Array containing the wavelengths for which the surface brightness maps are calculated. 
  real(Kind=real64), allocatable :: lambda_arr_maps(:)

  !> @param lambda_arr_SI(0 : lnum()-1) Array containing the wavelength list lambda_arr() transformed into SI units.
  real(Kind=real64), allocatable :: lambda_arr_SI(:)
  !> @param lambda_arr_SI_bin(0 : lnum()-2) Array containing the wavelength bin average value for integration over wavelength. By "wavelength bins" it is meant the bins defined by the lambda_arr_SI() array.  
  real(Kind=real64), allocatable :: lambda_arr_SI_bin(:)
  !> @param delta_lambda_bin(0: lnum()-1) Array containing the wavelength bin widths for integration over wavelength. Note that the bins are all containing the lambda_arr_SI() values, except the edges. In those cases the bins are either extending to the right of the lambda_arr_SI() value (for the lowest values of lambda_arr_SI()) or to the left (for the highest value). Look at the following ("|" are the wavelength grid points and "--", "__" denote a full bin. At the edges only the right or left part of the bin are considered): |-_|_-|-_|_-|-_|_-|. Note that the non-edge bin sizes are not necessarily the same. At the moment, they are the same in log space.  
  real(Kind=real64), allocatable ::  delta_lambda_bin(:)
  !> @param delta_lambda_bin_stars Same as delta_lambda_bin() but only for the wavelength range covered by the stellar emission (determined by lambda_arr() and max_lambda_stars()).
  real(Kind=real64), allocatable ::  delta_lambda_bin_stars(:)
  !> @param delta_lambda_bin_dust Same as delta_lambda_bin() but only for the wavelength range covered by the dust emission (determined by lambda_arr() and min_lambda_dust()).
  real(Kind=real64), allocatable ::  delta_lambda_bin_dust(:)
  !> @param lnum Number of wavelengths. Note that this value can change. For example, it is lnum_stars() during the RT calculation for the stellar emission and lnum_dust during the RT calculation for the dust emission.  
  integer :: lnum
  !> @param lambda_arr_SI_HD_dust Fine wavelength grid in the infrared range. Used to calculate wavelength integrated quantities which are independent on the input wavelength grid. 
  real(kind=real64), allocatable :: lambda_arr_SI_HD_dust(:)
  !> @param lambda_arr_SI_HD_dust_bin Central wavelength of the bins determined by lambda_arr_SI_HD_dust(). 
  real(kind=real64), allocatable :: lambda_arr_SI_HD_dust_bin(:)
  !> @param delta_lambda_HD_dust Wavelength bin sizes corresponding to lambda_arr_SI_HD_dust()
  real(kind=real64), allocatable :: delta_lambda_bin_HD_dust(:)
  !> @param lnum_stars Number of wavelengths for stellar emissivity grids.
  integer :: lnum_stars
  !> @param lnum_dust Number of wavelengths for dust emission wavelength grid.
  integer :: lnum_dust
  !> @param lnum_tot Total number of wavelengths in the input wavelength grid.
  integer :: lnum_tot
  !> @param max_lambda_stars Maximum wavelength used for the stellar emission wavelength grid. If not input, the code will consider all wavelengths in the input wavelength grid for the stellar emission RT calculation. 
  real(kind=real64) :: max_lambda_stars
  !> @param min_lambda_dust Minimum wavelength in the dust emission wavelength grid. If not input, the code sets default value min_lambda_dust = 1 micron.
  real(kind=real64) :: min_lambda_dust
  !> @param i_lambda_stars(0:1) Initial and final indeces of the stellar emission wavelength grid within lambda_arr().
  integer :: i_lambda_stars(0:1)
  !> @param i_lambda_dust(0:1) Initial and final indeces of the dust emission wavelength grid within lambda_arr().
  integer :: i_lambda_dust(0:1)
  
  ! Units
  !> @param units_i_obs Units of i_obs() arrays (used in sum_i_obs()). 
  character (len=20) :: units_i_obs
  !> @param units_ufield Units of radiation field energy density array u_final_uv_opt() needed when calculating dust emission. Note that these units are only for the radiation field due to stellar radiation (NOT for that due to dust radiation which is always in J/m/pc^3).  
  character (len=20) :: units_ufield
  !> @param units_lambda Units for the input wavelength array. The only choice is 'um' (microns). This is more like a reminder than an input parameter.  
  character (len=20) :: units_lambda
  !> @param units_luminosity Units of input stellar emission luminosity. Choices are 'erg/s/Hz', 'W/Hz'. 
  character (len=20) :: units_luminosity
  !> @param units_csize Units of grid cell sizes. The only possible choice is  'pc'. A reminder! 
  character (len=20) :: units_csize

  ! Units constants
  !> @param parsec parsec in [m].
  real(kind=real64), parameter ::  parsec=3.08567758E16 ! m
  !> @param parsec_cgs parsec in [cm].
  real(kind=real64), parameter ::  parsec_cgs=3.08567758E18 ! cm
  !> @param cspeed Light speed in [m/s].
  real(kind=real64), parameter :: cspeed=2.99792458E8 ! m/s
  
  ! Total mass parameters
  !> @param tot_gas_mass Total gas mass in the RT model.
  real(kind=real64) :: tot_gas_mass
  !> @param tot_dust_mass Total dust mass in the RT model.
  real(kind=real64) :: tot_dust_mass

  ! Dust heating variables
  !> @param dust_heating_type Type of dust heating. Choices are: 'eff' (effective single grain. Equilibrium dust temperature/emission calculated for a single grain with average opacity parameters); 'equ' (Equilibrium dust temperature/emission calculated for each grain specified by the input grain size distributions); 'sto' (Full non-equilibrium dust temperature/emission calculation for all grains specified by the input grain size distributions); 'sto_lib' (As 'sto' but using the adaptive SED library approach of Natale et al. (2015)).
  character (len=20) :: dust_heating_type

  !> @param dust_heating_type_ID Integer variable corresponding to the value of dust_heating_type(). 
  integer :: dust_heating_type_ID
  
  !> @param dt_none Integer variable corresponding to dust_heating_type() = 'not_provided'. In this case, no dust heating/emission calculation is performed.
  integer, parameter :: dt_none = -1
  
  !> @param dt_eff Integer variable corresponding to dust_heating_type() = 'eff'
  integer, parameter :: dt_eff = 0

  !> @param dt_equ Integer variable corresponding to dust_heating_type() = 'equ'
  integer, parameter :: dt_equ = 1

  !> @param dt_sto Integer variable corresponding to dust_heating_type() = 'sto'
  integer, parameter :: dt_sto = 2
    
  !> @param dt_sto_lib Integer variable corresponding to dust_heating_type() = 'sto_lib'
  integer, parameter :: dt_sto_lib = 3

  !> @param tot_dust_em_sed Dust emission sed library in [W/m/H] calculated during the SED adaptive library approach 
  real(kind=real64), allocatable :: tot_dust_em_sed(:,:,:)


  ! MPI parameters 

  !> @param np_mpi Number of MPI processes.
  integer :: np_mpi
  
  !> @param id_mpi ID number of the local MPI process.
  integer :: id_mpi
  
  !> @param main_prc Equal to TRUE if id_mpi = 0 
  logical :: main_prc

  ! USE variables
  !> @param use_lambda_grid Logical equal to TRUE if a lambda grid is used. If FALSE, input lambda_ref() is required. 
  logical :: use_lambda_grid
  !> @param use_dir_out Logical equal to TRUE if file_dir_out() has to be read.
  logical :: use_dir_out
  !> @param use_pos_obs Logical equal to TRUE if file_pos_obs() has to be read.
  logical :: use_pos_obs
  !> @param use_p_src Logical equal to TRUE if file_p_src() has to be read.
  logical :: use_p_src

  ! variables for grid creation
  !> @param grid_creation TRUE when subroutines are used by the grid creation programs. 
  logical :: grid_creation 
  !> @param grid_creation_lambda TRUE when subroutines are used by the grid creation programs for the calculation of the lambda grids.  
  logical :: grid_creation_lambda 

  ! constants 
  !> @param hplanck Planck constant in [J s]
  real(kind=real64), parameter :: hplanck = 6.62606957E-34

  !> @param kboltz Boltzmann constant in [J/K] 
  real(kind=real64), parameter :: kboltz=1.3806488E-23

  !> @param sigmaSB Stefan-Boltzmann constant in [W/m^2/K^4]
  real(kind=real64), parameter :: sigmaSB= 5.67036713E-8 

  !> @param msun Solar mass in [kg]
  real(kind=real64), parameter :: msun=1.9892000e+30

  !> @param m_h Atomic mass unit in [kg]
  real(kind=real64), parameter :: m_h = 1.6605402e-27


  ! VISUALIZATION ROUTINE VARIABLES
  !> @param npixel_maps Number of pixels for each side of the external observer maps.  
  integer :: npixel_maps
  !> @param print_maps TRUE if surface brightness maps for the external observers have to be calculated and printed.
  logical :: print_maps
  !> @param map_arr_out(x,y,lambda, direction) External observer surface brightness maps at the stellar emission or dust emission wavelengths stored in the local i_obs_arr() for all observer line-of-sight.  
  real(kind=real64), allocatable :: map_arr_out(:,:,:,:)
  !> @param map_size_factor The linear size of the external observer maps is equal to map_size_factor* modelsize(). Map_size_factor cannot be less than 0.2. However, a factor of at least 1.8 is suggested to avoid the possibility of cells not projected within the observer map. Note that if the i_obs array files are printed on disk, the maps can be recalculated using the 'sed' RT algorithm. 
  real(kind=real64) :: map_size_factor
  !> @param kp_maps Parameter determining the number of pixels npix_maps() on the internal observer surface brightness map. This is equal to 12*2**(2* kp_maps()).
  integer :: kp_maps
  !> @param print_maps_in TRUE if surface brightness maps for the internal observers have to be calculated and printed.
  logical :: print_maps_in
  !> @param npix_maps Number of pixel on the internal observer surface brightness maps. npix_maps = 12*2**(2* kp_maps())
  integer :: npix_maps
  !> @param map_in_arr_out(line-of-sight ID,lambda, observer) Internal observer surface brightness maps at the stellar emission or dust emission wavelengths stored in the local i_obs_in_arr() for all observer positions.  
  real(kind=real64), allocatable :: map_in_arr_out(:,:,:)
  !> @param size_map Size of each side of the output surface brightness maps for the external observer. This is map_size_factor() times modelsize(). Value saved in the output file_map(). 
  real(kind=real64) :: size_map


  ! N-body/SPH simulation variables

  !> @param tot_star_particles Number of stellar particles in the Nbody/SPH simulation file (file_nbody_sph() ).
  integer :: tot_star_particles 

  !> @param mstar Mass of the stellar particles in the Nbody/SPH simulation file (file_nbody_sph() ).
  real(kind=real64), allocatable :: mstar(:)

  !> @param starcoord Coordinates of the stellar particles in the Nbody/SPH simulation file (file_nbody_sph() ).
  real(kind=real64), allocatable :: starcoord(:,:)
  
  !> @param agestar Age of the stellar particles in the Nbody/SPH simulation file (file_nbody_sph() ).
  real(kind=real64), allocatable :: agestar(:)

  !> @param fehstar [Fe/H] of the stellar particles in the Nbody/SPH simulation file (file_nbody_sph() ).
  real(kind=real64), allocatable :: fehstar(:)

  !> @param tot_gas_particles Number of gas particles in the Nbody/SPH simulation file (file_nbody_sph() ).
  integer :: tot_gas_particles 

  !> @param mgas Mass of the gas particles in the Nbody/SPH simulation file (file_nbody_sph() ).
  real(kind=real64), allocatable :: mgas(:)
 
  !> @param gascoord Coordinates of the gas particles in the Nbody/SPH simulation file (file_nbody_sph() ).
  real(kind=real64), allocatable :: gascoord(:,:)

  !> @param gastemp Temperature of the gas particles in the Nbody/SPH simulation file (file_nbody_sph() ).
  real(kind=real64), allocatable :: gastemp(:)
  
  !> @param fehgas [Fe/H] of the gas particles in the Nbody/SPH simulation file (file_nbody_sph() ).
  real(kind=real64), allocatable :: fehgas(:)

  !> @param ofegas [O/Fe] of the gas particles in the Nbody/SPH simulation file (file_nbody_sph() ).
  real(kind=real64), allocatable :: ofegas(:)

  !> @param star_lum Stellar particle luminosity in erg/s/Hz
  real(kind=real64), allocatable :: star_lum(:)
  
  ! STELLAR LIBRARY PARAMETERS 

  !> @param lambda_lib Wavelength array in the input stellar library file (in microns). See file_stellar_library() and stellar_library(). 
  real(kind=real64), allocatable :: lambda_lib(:)

  !> @param nlambda_lib Number of wavelengths in the input stellar library file (see lambda_lib().
  integer :: nlambda_lib 

  !> @param met_lib Metallicity array in the input stellar library file. See file_stellar_library() and stellar_library(). 
  real(kind=real64), allocatable :: met_lib(:) 

  !> @param nmet_lib Number of metallicities in the input stellar library file. See met_lib(). 
  integer :: nmet_lib

  !> @param age_lib Ages of stellar populations in the input stellar library file.
  real(kind=real64), allocatable :: age_lib(:) 

  !> @param nage_lib Number of stellar ages in the input stellar library file. See age_lib(). 
  integer :: nage_lib

  !> @param lum_to_mass_lib(0:lambda_lib()-1, 0:nages_lib-1, 0:nmet_lib-1) Stellar luminosity-to-mass ratio values loaded from the input stellar library file. 
  real(kind=real64), allocatable :: lum_to_mass_lib(:,:,:)

  !> @param lum_to_mass_int(0: lnum_tot()-1, 0: nages_lib()-1, 0: nmet_lib()-1) Stellar luminosity-to-mass ratio values loaded from the input stellar library file and interpolated at the same wavelengths as in lambda_arr().  
  real(kind=real64), allocatable :: lum_to_mass_int(:,:,:)

  !> @param z_sun Solar metallicity value. The default value is 0.018, but it can be modified in the input if needed. Value checking is done within the routines that use this parameter.
  real(kind=real32) :: z_sun 

  !> @param pcell_star ID of the grid cell to which each star particles belong.
  integer, allocatable :: pcell_star(:)

  !> @param pcell_gas ID of the grid cell to which each gas particles belong.
  integer, allocatable :: pcell_gas(:)

  ! -------------------
  ! Scattering iteration parameters

  !> @param max_sca_iterations Maximum number of scattering iterations. It requires limit_scattering_iterations() = .TRUE. 
  integer :: max_sca_iterations

  !> @param limit_scattering_iterations TRUE if a maximum number of iterations has to be set using max_sca_iterations()
  logical :: limit_scattering_iterations 

CONTAINS

!> This subroutine allocates the 3D grid arrays and initialize the first elements. 
  subroutine create_grid_arrays()
    
    allocate(ncell(0:max_ncell-1), cchild(0:max_ncell-1),lvl(0:max_ncell-1),cindex(0:max_ncell-1), ccoord(3,0:max_ncell-1) &
         , csize(0:max_ncell-1),dens(0:max_ncell-1),dens_stars(0:max_ncell-1))

    allocate(nstart(0:max_lvl+2))
    
    !initialization 
    ncell=0
    cchild =0
    lvl=0
    cindex=0
    ccoord=0
    csize=0
    dens=0
    dens_stars=0
    nstart=0

    ! first cells initialization  
     tot_ncell=0     
      ncell(0)=0
      lvl(0)=0
      csize(0)=modelsize
      ncurr=0   
      nlevel=0
      cchild(0)=1
      nstart(0)=0
      nstart(1)=1

      ccoord(1,0)=0.
      ccoord(2,0)=0.
      ccoord(3,0)=0.

    end subroutine create_grid_arrays
    
    !> This subroutine derives the cellsize for a certain cell subdivision level.
    !> @param [in] clvl Cell subdivision level
    !> @param [out] cellsize Cell size 
  subroutine calc_cellsize(cellsize,clvl)
    
   real(kind=real64) :: cellsize 
   integer :: clvl 
  
   if (clvl > 0) then 
      cellsize=modelsize*baseinv(1)*baseinv(2)**(clvl-1)    
   else
      cellsize=modelsize
   endif

  end subroutine calc_cellsize 

!> This subroutine sets the base() derived quantities basediv() and basemask(). 
  subroutine set_base()
   
    integer(kind=int32) :: i,j, base3 

    if (main_prc) print *, 'Setting grid subdivision factors...' 
    
   ! set inverse of base (double number) 
    baseinv=1./dble(base)

   ! set basediv (factor used for division binary cell code)

    do i=1,2    
       j=0
       base3=base(i)**3
       do 
  
          if (base3 < 2**j) exit 
          j=j+1 
       end do
       basediv(i)=2**j
       
    enddo

    basemask=basediv-1 

    call print_done
    
  end subroutine set_base

  !> This routine finds the cell neighbour of the cell nc() along the main cell axis.
  !> @param [in] isel 1=x-dir 2=y-dir 3=z-dir
  !> @param [in] inc +1 or -1 depending if positive or negative direction.
  !> @param [out] cc ID number of the neighbour cell
  !> @param [out] out If = 1, the cell is at the border of RT model.
subroutine find_neighbours(isel, inc, cc,out)
  
IMPLICIT NONE
!INTEGER, PARAMETER :: dp= selected_int_kind(16)
integer :: i,j,k,l,ib, isel, inc  
integer(kind=int32) :: nc,nc_level,clvl,cc,out
integer, allocatable :: ccindd(:,:)
real(kind=real64) :: cellsize, b
logical :: cell_found

out =0 ! if out is equal one, the cell is at the border of scene

!isel=1 ! 1=x-dir 2=y-dir 3=z-dir
!inc=1 ! +1 or -1 depending on increment

nc=ncurr

nc_level=lvl(nc)
allocate(ccindd(3,max_lvl))
ccindd=0

! cell integer coodinates within parent cell

clvl=nc_level

!print *, 'Cell coordinates:'
call cindex_to_ccindd(nc,clvl,ccindd)

! neighbour cell integer coordinates

do i=1,3  ! loop on xyz directions

   if (i.eq.isel) then   ! direction intersected plane 
  
        k=inc
        
        do j=clvl,1,-1   !!! same size or bigger cells       
         call increment(j,ccindd(i,j),k,k)            
         end do 
    
        if (k.ne.0) then    ! if ray going outside model
          out=1 
         ! print *, 'cell=',nc 
         ! print *, 'it is OUT!'
         return
        endif

        ib=1
        do j=clvl+1,max_lvl   !!! smaller cells 
           if (j > 1) ib=2
         if (inc.eq.1) then
          ccindd(i,j)=0
         elseif (inc.eq.-1) then
          ccindd(i,j)=base(ib)-1
         else
         endif
       end do
    else
       ! Get actual crossing coordinate relative to current cell

       call calc_cellsize(cellsize,clvl)
    !cellsize=modelsize*baseinv(1)*baseinv(2)**(clvl-1) 
     b=0.5*cellsize              
     ! Subdivide it into integer coordinates from the current level up
     ib=1
     do j=clvl+1,max_lvl
        if (j >1) ib=2
      call calc_cellsize(cellsize,j) 
       ! cellsize=modelsize*baseinv(1)*baseinv(2)**(j-1)     
         ccindd(i,j)=floor(b/cellsize)          
        ! if (ccindd(i,j).ge.base(ib)) then  
        !  ccindd(i,j)=min(ccindd(i,j),base(ib)-1)          
        ! endif
         b=mod(b,cellsize)
!         print *, j
!         print *, i,ccindd(i,j)
!         read(*,*)
     end do
  endif

end do ! loop on xyz directions 

! Encode a new found cell and check it's existance
call ccindd_to_cc(ccindd,cc,clvl,cell_found)

deallocate(ccindd)

end subroutine find_neighbours


!> This routine is used to increment the cell coordinates in the tree representation along a specific axis. In this way, one finds the coordinates of the next cell.  
!> @param[in] ilvl Subdivision level
!> @param val Coordinate value (between 0 and base()-1)
!> @param[in] inc Equal to +1 or -1 depending if increment in the positive or negative direction. 
!> @param[out] ovr Flag equal to zero if the new cell is within the same parent cell block. Otherwise it stores the increment value.
subroutine increment(ilvl, val,inc,ovr)

IMPLICIT NONE 
integer :: val,inc,ovr,ib,ilvl

ib=1
if (ilvl > 1) ib=2

 val=val+inc     
 ovr=inc
 if (val.ge.base(ib)) then         
   val=0
 elseif (val.lt.0) then        
   val=base(ib)-1
 else        
   ovr=0
 endif

end subroutine increment

!> This subroutine is used by the grid creation program to compare the subdivision level of a cell neighbour to the cell ncurr() with the subdivision level that is or will be in the region of the cell ncurr(). Depending if the ncurr() cell will be further subdivided or not, the subroutine outputs a flag specifying if the neighbour cell has to be subdivided or not (see comment inside code to more precise explanations).   
!> @param [in] cc ID number of the neighbour cell.
!> @param [out] flag_jump Flag equal to 1 or 0 depending whether the neighbour cell has to be subdivided or not 
subroutine check_level_jump(cc,flag_jump)
IMPLICIT NONE 
integer :: cc,flag_jump

flag_jump=0

! this condition is valid for cells subdivided because of user defined conditions
if ((lvl(ncurr) == nlevel).and.(lvl(ncurr) > lvl(cc))) flag_jump=1

! this condition is valid for neighbour cells subdivided.
! Note that if lvl(ncurr)< nlevel then ncurr is not further subdivided. Only its neighbours if needed. For this reason there is the +1 in this formula compared to the one above. In this case lvl(ncurr) is the final subdivision level in the region covered by ncurr. In the condition above the region of ncurr is going to be further subdivided.   
if ((lvl(ncurr) < nlevel).and.(lvl(ncurr) >  (lvl(cc) +1))) flag_jump=1

end subroutine check_level_jump

!> This routines derives the cell tree coordinates from tree binary code in cindex().
!> @param [in] nc ID number of the cell
!> @param [in] clvl Subdivision level of the cell
!> @param [out] ccindd[3, max_lvl() ] Array containing cell tree coordinates for each subdivision level.  
subroutine cindex_to_ccindd(nc,clvl,ccindd)
integer :: ib, ccindd(:,:),i 
integer(kind=int32) :: clvl,nc,j


ib=1
do i=1,clvl
   if (i>1) then 
   ib=2    
      j=iand((cindex(nc))/(basediv(1)*basediv(2)**(i-2)),basemask(2))-1
   else
     j=iand((cindex(nc))/(basediv(1)**(i-1)),basemask(1))-1
   endif 
     ccindd(1,i)=mod(j,base(ib))
     ccindd(2,i)=mod(j/base(ib),base(ib))
     ccindd(3,i)=mod(j/base(ib)**2,base(ib))  
!print *, ccindd(1,i),ccindd(2,i), ccindd(3,i)

end do   

end subroutine cindex_to_ccindd

!> This routine finds the cell ID number from its tree coordinates.
!> @param [in] ccindd[3, max_lvl() ] Array containing cell tree coordinates for each subdivision level.
!> @param [out] cc Cell ID number
!> @param [out] clvl Cell subdivision level
!> @param [out] cell_found Logical variable equal to TRUE if the cell is found.
subroutine ccindd_to_cc(ccindd,cc,clvl,cell_found)
  
  integer :: ccindd(:,:),i,ib,clvl
  integer(kind=int32) :: cc,l,k
  logical :: cell_found

  l=1

  cell_found = .false.
   
  ib=1
  do i=1,max_lvl
     if (i > 1) ib=2 
     k=(ccindd(3,i)*base(ib)+ccindd(2,i))*base(ib)+ccindd(1,i)+1
     l=l+k-1
  
     if (cchild(l).eq.-1) then          
        cc=l
        clvl=i
        cell_found = .true.
 
        exit
     else
        l=cchild(l)       
     endif
  end do

!print *, 'find NB'
!print *, 'cell=',nc, '(nlvl=',nc_level,')'
!print *, 'neighbour =', cc, '(nlvl=',clvl,')'
!print *, ccindd(1,clvl),ccindd(2,clvl), ccindd(3,clvl)

end subroutine ccindd_to_cc

!> Prints 'DONE' on the terminal. It also sets a barrier for MPI processes.
 subroutine print_done
   integer :: ierr

   if (main_prc) then 
      print *, 'DONE'
      print *, '-----------------------------'
      print *, ''
   endif 

   call mpi_barrier(MPI_COMM_WORLD,ierr)
   
 end subroutine print_done

!> Exits process. 
subroutine stop_prc
  
  character*16 :: id_str
  character*80 :: out_msg

  write(id_str,*) id_mpi 

  out_msg = 'STOP process '//trim(adjustl(id_str))
  print *, out_msg

  STOP 

end subroutine stop_prc

!> Exits the program if the input variable error is TRUE. 
subroutine error_stop(error)
logical :: error

if (error) then 
   if (main_prc) print *, 'STOP: some errors found (see ERROR messages above)'
   call stop_prc
endif

end subroutine error_stop


!> Sort integer array "list" and returns array of original positions (order). Wrapper to quick_sort(). 
subroutine quick_sort_int(list, order)
  integer, DIMENSION (:), INTENT(IN OUT)  :: list
  integer, DIMENSION (:), INTENT(OUT)  :: order
  real(kind=real64), allocatable :: list_real(:)
  integer :: num 

  num = size(list)
  allocate(list_real(0:num-1))
  list_real = real(list)

  call quick_sort(list_real, order)

  list = int(list_real)

  deallocate(list_real)

end subroutine quick_sort_int


!> Sort array "list" using quick sort algorithm. It also returns an array "order" containing the position of each element in the original array. This version works for real64 variables. To sort other types of variable arrays, use the appropriate wrapper (e.g. quick_sort_int).  
RECURSIVE SUBROUTINE quick_sort(list, order)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.
! NOTE by GN: downloaded at http://jblevins.org/mirror/amiller/qsort.f90 and modified 

IMPLICIT NONE
real(kind=real64), DIMENSION (:), INTENT(IN OUT)  :: list
integer, DIMENSION (:), INTENT(OUT)  :: order

! Local variable
INTEGER :: i

DO i = 1, SIZE(list)
  order(i) = i
END DO

CALL quick_sort_1(1, SIZE(list))

CONTAINS

RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
real(kind=real64)   :: reference, temp
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort_1(left_end, j)
  IF (i < right_end) CALL quick_sort_1(i, right_end)
END IF

END SUBROUTINE quick_sort_1


SUBROUTINE interchange_sort(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
real(kind=real64)   :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    END IF
  END DO
END DO

END SUBROUTINE interchange_sort

END SUBROUTINE quick_sort


END MODULE smooth_grid_routines
