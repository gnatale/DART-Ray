!> Contains the subroutines to specify the stellar emission and dust density distribution for the analytical galaxy geometry.
MODULE user_routines_galaxy
use smooth_grid_routines
use iso_fortran_env
use HDF5
use io_routines
use sed_routines
IMPLICIT NONE  
 
! Galaxy parameters 
!> @param lambda_in Wavelength variable 
real(kind=real64) :: lambda_in
!> @param lambda_min Minimum wavelength used in the calculation of the lambda grids.  
real(kind=real64) :: lambda_min
!> @param lambda_max Maximum wavelength used in the calculation of the lambda grids. 
real(kind=real64) :: lambda_max
! Old stellar disk 
!> @param old Scaling factor applied to the input old stellar luminosity spectra.  
real(kind=real64) ::  old
!> @param lnu_old Old stellar luminosity spectra. 
real(kind=real64) :: lnu_old
!> @param  hs_disk_b B-band (4430 A) old stellar disk scale length. 
real(kind=real64) :: hs_disk_b
!> @param zs_disk Scale height for the old stellar disk (basic value). 
real(kind=real64) :: zs_disk
!> @param zs_disk_r1 Scale height for the old stellar disk at the inner radius hsin().
real(kind=real64) :: zs_disk_r1
!> @param zs_disk_rsun Scale height for the old stellar disk at the radius rsun(). 
real(kind=real64) :: zs_disk_rsun
!> @param hsin Inner radius of the old stellar disk. 
real(kind=real64) :: hsin
!> @param hs_disk Old disk scale length or X'-axis scale length if elliptical exponential model is used. 
real(kind=real64) :: hs_disk
!> @param hs_disk2  Y'-axis old disk  scale length for the elliptical exponential model.  
real(kind=real64) :: hs_disk2
!> @param eta_disk0 Central value of the old stellar disk volume emissivity. 
real(kind=real64) :: eta_disk0
!> @param chi_disk Parameter used in the "flared" profile to specify the steepness of the radial profile within the inner radius hsin().   
real(kind=real64) :: chi_disk
!> @param hs_disk_arr Scale lenghts values for the old stellar disk at each wavelength of the input wavelength grid. In the case of the elliptical exponential profile, this is the X'-axis scale length. By default, all values are set equal to hs_disk_b. If different values are required, hs_disk_arr can be specified in the input. The wavelength indeces of the input hs_disk values have to be specified in id_hs_disk_arr()
real(kind=real64), allocatable :: hs_disk_arr(:)
!> @param hs_disk2_arr Same as hs_disk_arr() but for the Y'-axis scale length required in the elliptical exponential disk profile type. 
real(kind=real64), allocatable :: hs_disk2_arr(:)
!> @param id_hs_disk_arr Wavelength grid indeces of the values of hs_disk specified in the input hs_disk_arr() and hs_disk2_arr().
integer, allocatable :: id_hs_disk_arr(:)
!> @param maxsize_hs_disk_arr Maximum number of input values in hs_disk_arr() and hs_disk2_arr() 
integer, parameter  :: maxsize_hs_disk_arr = 1000
!> @param theta_disk_ellipt Rotation angle around the Z-axis for the elliptical exponential profile for the old stellar component [deg] 
real(kind=real64) :: theta_disk_ellipt

! YOUNG STELLAR DISK
!> @param sfr Scaling factor applied to the young stellar disk luminosity spectra.
real(kind=real64) ::  sfr
!> @param  lnu_sf Young stellar disk luminosity. 
real(kind=real64) :: lnu_sf
!> @param hs_tdisk Scale length of the young stellar disk or, in the case of the ellipsoidal exponential profile, X'-axis scale length 
real(kind=real64) :: hs_tdisk
!> @param hs_tdisk2 Y'-axis young stellar disk scale length for the ellipsoidal esponential profile. 
real(kind=real64) :: hs_tdisk2
!> @param zs_tdisk Scale height of the young stellar disk (basic value). 
real(kind=real64) ::  zs_tdisk
!> @param hs1in Inner radius of the young stellar disk.  
real(kind=real64) ::  hs1in
!> @param eta_tdisk0 Central value of the young stellar disk volume emissivity.
real(kind=real64) :: eta_tdisk0
!> @param zs_tdisk_r1 Scale height of the young stellar disk at the inner radius hs1in().
real(kind=real64) :: zs_tdisk_r1
!> @param zs_tdisk_rsun Scale height of the young stellar disk at the radius rsun().
real(kind=real64) :: zs_tdisk_rsun
!> @param chi_tdisk Parameter used in the "flared" profile to specify the steepness of the radial profile within the inner radius hs1in(). 
real(kind=real64) :: chi_tdisk
!> @param theta_tdisk_ellipt Rotation angle around the Z-axis for the elliptical exponential profile for the young stellar component [deg].
real(kind=real64) :: theta_tdisk_ellipt
! BULGE
!> @param reff Bulge effective radius. 
real(kind=real64) ::  reff
!> @param lnu_bulge Bulge volume emissivity. 
real(kind=real64) :: lnu_bulge
!> @param acap_bulge Inner truncation radius for the bulge profile
real(kind=real64) :: acap_bulge
!> @param ellipt Bulge ellipticity. The normalized radial coordinate m for the bulge volume emissivity is defined as 
!> \f[
!> m=\sqrt{R^2 + (z/ellipt)^2}/ R_{eff}
!> \f]
real(kind=real64) :: ellipt
!> @param ellipt_xy Bulge XY ellipticity. The radial coordinate R for the bulge volume emissivity is defined as 
!> \f[
!> R=\sqrt{x^2 + (y/ellipt_{xy})^2}
!> \f]
real(kind=real64) :: ellipt_xy
!> @param mtrunc Outern truncation radius for the bulge in units of effective radii. 
real(kind=real64) :: mtrunc
!> @param bd_ratio Bulge-to-Disk luminosity ratio. Used to define luminosity scaling for the bulge, that is, bd_ratio*old*lnu_old.  
real(kind=real64) :: bd_ratio
!> @param eta_bulge0 Central value of the bulge volume emissivity.
real(kind=real64) :: eta_bulge0
!> @param nsersic Bulge sersic index (see inside subroutine av_star_bulge() for the volume emissivity formula used in each case).  
integer :: nsersic 
!> @param theta_bulge Rotation angle around the Z-axis for the bulge profile. 
real(kind=real64) :: theta_bulge

! total luminosity
!> @param lnu_tot Model total luminosity. 
real(kind=real64) :: lnu_tot
! THICK DUST DISK
!> @param hd_disk Scale length of the thick dust disk or X'-axis scale length in the case of the ellipsoidal exponential profile  
real(kind=real64) ::  hd_disk
!> @param hd_disk2 Y'-axis thick dust disk scale length for the ellipsoidal exponential profile. 
real(kind=real64) :: hd_disk2
!> @param zd_disk Scale height for the thick dust disk (basic value). 
real(kind=real64) ::  zd_disk
!> @param zd_disk_r1 Scale height for the thick dust disk at the inner radius hdin().
real(kind=real64) ::  zd_disk_r1
!> @param zd_disk_rsun Scale height for the thick dust disk at the radius rsun(). 
real(kind=real64) ::  zd_disk_rsun
!> @param hdin Inner radius for the thick dust disk.  
real(kind=real64) ::  hdin
!> @param tau1 B-band (4430 A) central face-on optical depth of the thick dust disk in the case of the pure exponential profile. In the case of the flared profiles, this is the optical depth one would have if the exponential radial profile for R>R1 would continue until R=0.  
real(kind=real64) ::  tau1
!> @param kext_disk0 Central value for the thick dust disk density profile.
real(kind=real64) ::  kext_disk0
!> @param chi_dust_disk Parameter used in the "flared" profile to specify the steepness of the radial profile within the inner radius hdin().
real(kind=real64) ::  chi_dust_disk
!> @param theta_dust_disk_ellipt Rotation angle around the Z-axis for the elliptical exponential profile for the thick dust disk component [deg].
real(kind=real64) :: theta_dust_disk_ellipt

! THIN DUST DISK 
!> @param hd_tdisk Scale length of the thin dust disk or X'-axis scale length for the ellipsoidal exponential profile. 
real(kind=real64) ::  hd_tdisk
!> @param hd_tdisk2 Y'-axis thin dust disk scale length for the ellipsoidal exponential profile. 
real(kind=real64) :: hd_tdisk2
!> @param zd_tdisk Scale height for the thin dust disk (basic value). 
real(kind=real64) :: zd_tdisk
!> @param hd1in Inner radius for the thin dust disk. 
real(kind=real64) :: hd1in
!> @param tau2 B-band (4430 A) central face-on optical depth of the thin dust disk in the case of the pure exponential profile. Other profiles still use this factor to set the dust density although its meaning is different. 
real(kind=real64) :: tau2
!> @param kext_tdisk0 Central value for the thin dust disk density profile.
real(kind=real64) :: kext_tdisk0
!> @param zd_tdisk_r1 Scale height for the thin dust disk at the inner radius hd1in().
real(kind=real64) :: zd_tdisk_r1
!> @param zd_tdisk_rsun Scale height for the thin dust disk at the radius rsun(). 
real(kind=real64) :: zd_tdisk_rsun
!> @param chi_dust_tdisk Parameter used in the "flared" profile to specify the steepness of the radial profile within the inner radius hd1in().  
real(kind=real64) :: chi_dust_tdisk 
!> @param theta_dust_tdisk_ellipt Rotation angle around the Z-axis for the elliptical exponential profile for the thin dust disk component [deg].
real(kind=real64) :: theta_dust_tdisk_ellipt

! Other geometrical parameters 
!> @param rtrun Radial truncation radius for all disks. 
real(kind=real64) ::  rtrun
!> @param rsun Third radial position used in the "flared" profile. Called rsun because it has been used to set the sun position in Milky Way model. 
real(kind=real64) ::  rsun
!> @param max_z Maximum |z| allowed for non zero stellar emissivity and dust density.  
real(kind=real64) ::  max_z
!> @param max_rad Maximum R allowed for non zero stellar emissivity and dust density. 
real(kind=real64) ::  max_rad
!> @param sha Factor determining the sharpness of the profiles close to the truncation radius rtrun(). The profile goes to zero within a radial interval equal to sha() in units of rtrun(). This factor is for the thick stellar and dust disks.   
real(kind=real64) ::  sha
!> @param sha1 Factor determining the sharpness of the profiles close to the truncation radius rtrun(). The profile goes to zero within a radial interval equal to sha1() in units of rtrun(). This factor is for the thin stellar and dust disks.     
real(kind=real64) ::  sha1
!> @param omega_max Maximum subtended solid angle of a cell seen from the position (rsun(), 0,0). To use this variable you have to make sure that the corresponding line in the function subdivision is used while compiling the program.  
real(kind=real64) ::  omega_max
! Filenames and keywords 
!> @param lcar_type String variable length parameter. 
integer, parameter :: lcar_type=30
! ! @param lcar2 String variable length parameter. 
!integer, parameter :: lcar2=80

!> @param old_disk_type Old stellar disk profile type. Possible choices:
!> - 'expR_expz' :   
!> \f[
!> \rho = \rho_o exp(-R/hs -|z|/z_s)
!> \f]
!> - 'expR_sech2z':
!> \f[
!> \rho(R,z) = \rho_o exp(-R/hs)sech^2(|z|/z_s)
!> \f]
!> - 'flared_expz': 
!> \f[
!> \gamma = log(\frac{zs(R_{sun})-zs(0)}{zs(R1)-zs(0)})/log(\frac{R_{sun}}{R1})
!> \f]
!> \f[
!> zs(R)=zs(0)+(zs(R1)-zs(0))*(R/R1)^{\gamma}
!> \f]
!> \f[
!> \rho(R,z) = \rho_o [\frac{R}{R1}*(1-\chi) + \chi]*\frac{zs(0)}{zs(R)}*exp(-\frac{R1}{hs})*exp(-\frac{|z|}{zs(R)})\qquad if ~  R < R1              
!> \f]
!> \f[
!> \rho(R,z) = \rho_o \frac{zs(0)}{zs(R)}*exp(-\frac{R}{hs})*exp(-\frac{|z|}{zs(R)}) \qquad if ~ R \ge R1              
!> \f]
!> - 'flared_sech2z':
!> as 'flared_expz' but with \f$sech^2(\frac{|z|}{zs(R)})\f$ instead of \f$exp(-\frac{|z|}{zs(R)})\f$; 
!> - 'ellipt_expR_expz':
!> \f$zs(R)\f$ and \f$\gamma\f$ defined as for 'flared_expz'
!> \f[ 
!> \rho(x,y,z) = \rho_o \frac{zs(0)}{zs(R)}*exp\left(-\sqrt{(\frac{x}{hs})^2+ (\frac{y}{hs2})^2}\right)*exp(-\frac{|z|}{zs(R)}) \qquad if ~  R < R1
!> \f]
!> \f[ 
!>  \rho(x,y,z) = 0 \qquad if ~ R \ge R1   
!> \f]
!> - 'ellipt_expR_sech2z':
!> as 'ellipt_expR_expz' but with \f$sech^2(\frac{|z|}{zs(R)})\f$ instead of \f$exp(-\frac{|z|}{zs(R)})\f$.   
character(LEN=lcar_type) :: old_disk_type
!> @param old_disk_type_ID ID number for old_disk_type().
integer :: old_disk_type_ID
!> @param young_disk_type Young stellar disk profile type. Same choices as for old_disk_type().
character(LEN=lcar_type) :: young_disk_type
!> @param young_disk_type_ID ID number for young_disk_type().
integer :: young_disk_type_ID
!> @param thick_disk_type Thick dust disk profile type. Same choices as for old_disk_type().
character(LEN=lcar_type) :: thick_disk_type
!> @param thick_disk_type_ID ID number for thick_disk_type().
integer::  thick_disk_type_ID
!> @param thin_disk_type Thin dust disk profile type. Same choices as for old_disk_type().
character(LEN=lcar_type) :: thin_disk_type
!> @param thin_disk_type_ID ID number for thin_disk_type().
integer::  thin_disk_type_ID
!> @param Specify which stellar emission components have to be included: 'all'= all components; 'disk' = only old stellar disk; 'tdisk' = only thin stellar disk; 'bulge' = only bulge component.    
character(LEN=lcar_type) :: grid_type

! GRID TYPES ID 
integer, parameter :: expR_expz_ID = 0 
integer, parameter :: expR_sech2z_ID = 1
integer, parameter :: flared_expz_ID = 2 
integer, parameter :: flared_sech2z_ID = 3 
integer, parameter :: ellipt_expR_expz_ID = 4 
integer, parameter :: ellipt_expR_sech2z_ID = 5  

! Extra grid arrays 
!> @param dens_disk Volume emissivity for the old stellar disk. 
real(kind=real64), allocatable  :: dens_disk(:)
!> @param dens_tdisk Volume emissivity for the young stellar disk. 
real(kind=real64), allocatable  :: dens_tdisk(:)
!> @param dens_bulge Volume emissivity for the bulge.
real(kind=real64), allocatable  :: dens_bulge(:)
!> @param dens_dust_disk Dust density for the thick dust disk. 
real(kind=real64), allocatable  :: dens_dust_disk(:)
!> @param dens_dust_tdisk Dust density for the thin dust disk. 
real(kind=real64), allocatable  :: dens_dust_tdisk(:) 
! parameters for grid subdivision
!> @param Solid angle subtended by a cell 
real(kind=real64) :: omega_cell
! integration step numbers
integer, parameter :: step_int=10
! Dust emission arrays amd parameters from 2D code
!!$real(KIND=real64), allocatable :: g_ccoord(:,:),g_lum(:)
!!$real(kind=real64) :: max_r_dem,max_z_dem
!!$integer :: nr_dem,nz_dem
!!$character(LEN=250) :: code_model_dem,file_grid_dem
!!$real(kind=real64) :: scaling_dust_comp_arr(3)

!> @param file_old_star_sed File containing the SED to be used by the old stellar components in units luminosities to be multiplied by the OLD parameter. 
character(LEN=lcar) :: file_old_star_sed
!> @param file_young_star_sed File containing the SED to be used by the young stellar components in units luminosities to be multiplied by the SFR parameter. 
character(LEN=lcar) :: file_young_star_sed
!> @param lambda_arr_old_sed Wavelength array in [um] for the old stellar luminosity sed
real(kind=real64), allocatable :: lambda_arr_old_sed(:)  ! um
!> @param nlambda_old_sed Number of elements in lambda_arr_old_sed()
integer :: nlambda_old_sed
!> @param lambda_arr_sf_sed Wavelength array in [um] for the young stellar luminosity sed
real(kind=real64), allocatable :: lambda_arr_sf_sed(:)  ! um
!> @param nlambda_sf_sed Number of elements in lambda_arr_sf_sed()
integer :: nlambda_sf_sed
!> @param lnu_old_unit Old stellar luminosity sed in unit luminosities for the old() parameter.  
real(kind=real64), allocatable :: lnu_old_unit(:)   ! W/Hz 
!> @param lnu_sf_unit Young stellar luminosity sed in unit luminosities for the young() parameter.  
real(kind=real64), allocatable :: lnu_sf_unit(:)  ! W/Hz  

!> @param z_subd_lim Z coordinate below which cells have to be subdivided if their R is less than R_subd_lim()
real(kind=real64) :: z_subd_lim 

!> @param z_subd_lim2 Z coordinate used in the 'milky_way' subdivision criteria to limit the subdivisions of cells with small subtended solid angle from the position (rsun(), 0,0). See subdivision_criteria(). 
real(kind=real64) :: z_subd_lim2 

!> @param R_subd_lim R coordinate below which cells have to be subdivided if their z is less than z_subd_lim()
real(kind=real64) :: R_subd_lim  

!> @param subdivision_criteria Subdivision criteria to be used. Choices are:
!> - 'standard': subdivision IF subdivision level < max_lvl() AND (cell optical depth > max_dtau() OR stellar disk luminosity > max_dlum()*total luminosity OR  subdivision level < min_lvl() OR (ABS(z) < z_subd_lim AND R < R_subd_lim()))
!> 'milky way' : subdivision IF subdivision level < max_lvl() AND (cell optical depth > max_dtau() OR stellar disk luminosity > max_dlum()*total luminosity OR  subdivision level < min_lvl() OR (ABS(z) < z_subd_lim AND R < R_subd_lim()) OR (cell subtended solid angle from position (rsun(), 0,0) > omega_max AND ABS(z)-half cell size < z_subd_lim2))  
character(LEN=lcar_type) :: subdivision_criteria


! input namelist 
namelist /galaxy_input_strings/ label_model_lambda_grid, dir_grid, grid_file, grid_info_file, file_lambda_list,units_lambda, grid_type, old_disk_type, young_disk_type , thick_disk_type, thin_disk_type, dust_model,dust_opacity_tables, file_gra_fa, file_sil_fa, file_pah_neu_fa, file_pah_ion_fa, file_av_opacities, file_q_gra, file_q_sil, file_q_pah_neu, file_q_pah_ion, file_old_star_sed, file_young_star_sed, subdivision_criteria

namelist /galaxy_input_var/ lambda_ref, lambda_min, lambda_max, rtrun, rsun,    max_z, max_rad, sha, sha1, omega_max,  modelsize, base, max_ncell, max_lvl, min_lvl,  max_dtau, max_dlum, n_dust_size_qabs, n_dust_wave_qabs, z_subd_lim, R_subd_lim, z_subd_lim2

namelist /galaxy_input_var_old_disk/ old, hs_disk_b, zs_disk, zs_disk_r1, zs_disk_rsun, chi_disk,  hsin, hs_disk_arr, hs_disk2_arr, id_hs_disk_arr, theta_disk_ellipt

namelist /galaxy_input_var_young_disk/ sfr, hs_tdisk, hs_tdisk2, zs_tdisk, zs_tdisk_r1,  zs_tdisk_rsun, chi_tdisk, hs1in, theta_tdisk_ellipt

namelist /galaxy_input_var_bulge/ reff, acap_bulge, ellipt, mtrunc, bd_ratio,  nsersic, theta_bulge, ellipt_xy

namelist /galaxy_input_var_thick_dust_disk/ tau1, hd_disk, hd_disk2, zd_disk, zd_disk_r1,  zd_disk_rsun, chi_dust_disk, hdin, theta_dust_disk_ellipt

namelist /galaxy_input_var_thin_dust_disk/ tau2, hd_tdisk, hd_tdisk2, zd_tdisk, zd_tdisk_r1,  zd_tdisk_rsun, chi_dust_tdisk, hd1in, theta_dust_tdisk_ellipt 

namelist /galaxy_input_logical/ input_av_opacities

CONTAINS

!> Reads input file for the grid creation for the galaxy 2D model. 
subroutine read_input_galaxy
integer :: ioun, i, id  
character(len=lcar) :: input_filename
integer :: n_hs
real(kind= real64), allocatable :: temp_arr(:), temp2_arr(:)

! allocate hs_disk_arr
allocate(hs_disk_arr(0:maxsize_hs_disk_arr-1), hs_disk2_arr(0:maxsize_hs_disk_arr-1), id_hs_disk_arr(0:maxsize_hs_disk_arr-1), temp_arr(0:maxsize_hs_disk_arr-1), temp2_arr(0:maxsize_hs_disk_arr-1))
id_hs_disk_arr = -1 
hs_disk_arr = 0
hs_disk2_arr = 0
temp_arr = 0 
temp2_arr = 0 

! Some variable initializations
call initialize_input_galaxy

! read input argument 
call get_command_argument(1,input_filename)

!check input_file argument has been given
if (trim(input_filename) == '') then 
   if (main_prc) print *, 'provide input file in command line...EXIT' 
   call stop_prc
endif

! read input file 
if (main_prc) print *, 'reading input file...'

open(newunit=ioun, file = input_filename, status ='old', action ='read')
read(ioun, nml = galaxy_input_strings )
read(ioun, nml = galaxy_input_var )
read(ioun, nml = galaxy_input_var_old_disk)
read(ioun, nml = galaxy_input_var_young_disk )
read(ioun, nml = galaxy_input_var_bulge )
read(ioun, nml = galaxy_input_var_thick_dust_disk )
read(ioun, nml = galaxy_input_var_thin_dust_disk )
read(ioun, nml = galaxy_input_logical)
close(ioun)

! check dust model
call check_input_dust_model

! set hs_disk_arr to proper values 
temp_arr = hs_disk_arr ! input hs_disk_arr
hs_disk_arr = hs_disk_b   ! set all values to B-band hs_disk
hs_disk2_arr = hs_disk_b
n_hs = count(id_hs_disk_arr /= -1)

do i = 0, n_hs-1
   id = id_hs_disk_arr(i)
   hs_disk_arr(id) = temp_arr(i)
   hs_disk2_arr(id) = temp2_arr(i)
end do

deallocate(temp_arr, temp2_arr,  id_hs_disk_arr)

! set grid creation true 
grid_creation = .TRUE.

! assign disk_type_ID 
call assign_disk_type_ID(old_disk_type, old_disk_type_ID)
call assign_disk_type_ID(young_disk_type, young_disk_type_ID)
call assign_disk_type_ID(thick_disk_type, thick_disk_type_ID)
call assign_disk_type_ID(thin_disk_type, thin_disk_type_ID)

! check input variables 
call check_input_galaxy

call print_done

end subroutine read_input_galaxy

!> Initialize variables specific for the galaxy grid creation program. 
subroutine initialize_input_galaxy

! galaxy input string 
label_model_lambda_grid = 'not_provided'
dir_grid = 'not_provided'
grid_file = 'not_provided'
grid_info_file = 'not_provided'
file_lambda_list = 'not_provided'
units_lambda = 'not_provided'
grid_type = 'not_provided'
old_disk_type = 'not_provided'
young_disk_type = 'not_provided'
thick_disk_type = 'not_provided'
thin_disk_type = 'not_provided'
file_old_star_sed = 'not_provided'
file_young_star_sed = 'not_provided'
subdivision_criteria = 'not_provided' 

! galaxy input var 
lambda_ref = 0 
lambda_min = 0 
lambda_max = 0 
rtrun = -1 
rsun = -1 
max_z = -1 
max_rad = -1 
sha = -1 
sha1 = -1 
omega_max = -1 
modelsize = -1 
base = [0,0] 
max_ncell = -1 
max_lvl = -1 
min_lvl = -1 
max_dtau = -1 
max_dlum = -1 
n_dust_size_qabs = 0 
n_dust_wave_qabs = 0 
z_subd_lim = -1 
R_subd_lim = -1 
z_subd_lim2 = -1

! old stellar disk 
old = -1 
hs_disk_b = -1
zs_disk = -1
zs_disk_r1 = -1
zs_disk_rsun = -1 
chi_disk = -100000
hsin = -1
!hs_disk_arr ! initialized above
!hs_disk2_arr
!id_hs_disk_arr
theta_disk_ellipt = -1 

! young stellar disk 
sfr = -1 
hs_tdisk = -1
hs_tdisk2 = -1
zs_tdisk = -1
zs_tdisk_r1 = -1
zs_tdisk_rsun = -1
chi_tdisk = -100000
hs1in = -1 
theta_tdisk_ellipt = -1

! bulge 
reff = -1 
acap_bulge = -1
ellipt = -1
mtrunc = -1
bd_ratio = -1
nsersic = -1
theta_bulge = -361
ellipt_xy =  -1

! thick dust disk
tau1 = -1
hd_disk = -1
hd_disk2 = -1
zd_disk = -1
zd_disk_r1 = -1
zd_disk_rsun = -1
chi_dust_disk = -100000
hdin = -1
theta_dust_disk_ellipt = -1

! thin dust disk
tau2 = -1
hd_tdisk = -1
hd_tdisk2 = -1 
zd_tdisk = -1
zd_tdisk_r1 = -1
zd_tdisk_rsun = -1
chi_dust_tdisk = -100000
hd1in = -1 
theta_dust_tdisk_ellipt = -1

! input logical 
input_av_opacities = .FALSE.


end subroutine initialize_input_galaxy

!> Checks that the input variables have allowed values. 
subroutine check_input_galaxy

logical :: error 

error = .FALSE.

! label_model_lambda_grid 
if (label_model_lambda_grid == 'not_provided') then 
   print *,  'ERROR: Input label_model_lambda_grid missing'
   error = .TRUE.
endif

! dir_grid 
if (dir_grid == 'not_provided') then 
   print *,  'ERROR: Input dir_grid missing'
   error = .TRUE.
else
   call check_dir_existence(dir_grid, .FALSE.)  
endif

! grid_file
if (grid_file == 'not_provided') then 
   print *,  'ERROR: Input grid_file missing'
   error = .TRUE.
endif

! grid_file_info
if (grid_info_file == 'not_provided') then 
   print *,  'ERROR: Input grid_info missing'
   error = .TRUE.
endif

! file_lambda_list
if (file_lambda_list == 'not_provided') then 
   print *,  'ERROR: Input file_lambda_list missing'
   error = .TRUE.
endif

! units_lambda
if (units_lambda == 'not_provided') then 
   print *,  'ERROR: Input units_lambda missing (it has to be always um)'
   error = .TRUE.
endif

!grid_type 
if (grid_type == 'not_provided') then 
   print *,  'ERROR: Input grid_type missing'
   error = .TRUE.
endif

!old_disk_type 
if (old_disk_type == 'not_provided') then 
   print *,  'ERROR: Input old_disk_type missing'
   error = .TRUE.
endif

!young_disk_type
if (young_disk_type == 'not_provided') then 
   print *,  'ERROR: Input young_disk_type missing'
   error = .TRUE.
endif

!thick_disk_type
if (thick_disk_type == 'not_provided') then 
   print *,  'ERROR: Input thick_disk_type missing'
   error = .TRUE.
endif

!thin_disk_type
if (thin_disk_type == 'not_provided') then 
   print *,  'ERROR: Input thin_disk_type missing'
   error = .TRUE.
endif

!file_old_star_sed
if (file_old_star_sed == 'not_provided') then 
   print *,  'ERROR: Input file_old_star_sed missing'
   error = .TRUE.
endif

!file_young_star_sed
if (file_young_star_sed == 'not_provided') then 
   print *,  'ERROR: Input file_young_star_sed missing'
   error = .TRUE.
endif

!subdivision_criteria
if (subdivision_criteria == 'not_provided') then 
   print *,  'ERROR: Input subdivision_criteria missing'
   error = .TRUE.
endif

! lambda_ref
if (lambda_ref < 0.089 .or. lambda_ref > 1000) then
     print *, 'ERROR: Invalid lambda_ref value'
     print *, 'lambda_ref = ', lambda_ref 
     print *, 'Allowed range = [0.09,1000]'
     error = .TRUE.
  endif

! lambda_min
if (lambda_min < 0.089 .or. lambda_min > 1000) then
     print *, 'ERROR: Invalid lambda_min value'
     print *, 'lambda_min = ', lambda_min 
     print *, 'Allowed range = [0.09,1000]'
     error = .TRUE.
  endif

! lambda_max
if (lambda_max < 0.089 .or. lambda_max > 1000) then
     print *, 'ERROR: Invalid lambda_max value'
     print *, 'lambda_min = ', lambda_max 
     print *, 'Allowed range = [0.09,1000]'
     error = .TRUE.
  endif

! lambda_min < lambda_max
if (lambda_min > lambda_max) then 
   print *, 'ERROR: lambda_min > lambda_max!' 
   error = .TRUE. 
endif

! rtrun
if (rtrun < 0.) then
   print *, 'ERROR: Invalid rtrun value'
   print *, 'rtrun = ', rtrun 
   print *, 'Allowed range = [0.,Inf]'
   error = .TRUE.
endif

! rsun
if (rsun < 0.) then
   print *, 'ERROR: Invalid rsun value'
   print *, 'rsun = ', rsun 
   print *, 'Allowed range = [0.,Inf]'
   error = .TRUE.
endif

! max_z
if (max_z < 0.) then
   print *, 'ERROR: Invalid max_z value'
   print *, 'max_z = ', max_z 
   print *, 'Allowed range = [0.,Inf]'
   error = .TRUE. 
endif

! max_rad
if (max_rad < 0.) then
   print *, 'ERROR: Invalid max_rad value'
   print *, 'max_rad = ', max_rad
   print *, 'Allowed range = (0.,Inf]'
   error = .TRUE.
endif

! sha
if (sha < 0.) then
   print *, 'ERROR: Invalid sha value'
   print *, 'sha = ', sha
   print *, 'Allowed range = [0.,Inf]'
   error = .TRUE.
endif

! sha1
if (sha1 < 0.) then
   print *, 'ERROR: Invalid sha1 value'
   print *, 'sha1 = ', sha1
   print *, 'Allowed range = [0.,Inf]'
   error = .TRUE.
endif

if (subdivision_criteria == 'milky_way') then 
   ! omega_max
   if (omega_max < 0.) then
      print *, 'ERROR: Invalid omega_max value'
      print *, 'omega_max = ', omega_max
      print *, 'Allowed range = [0.,Inf]'
      error = .TRUE.
   endif
endif

! modelsize
if (modelsize <= 0.) then
   print *, 'ERROR: Invalid modelsize value'
   print *, 'modelsize = ', modelsize
   print *, 'Allowed range = (0.,Inf]'
   error = .TRUE.
endif

! base
if (minval(base) < 2) then 
   print *, 'ERROR: Invalid base values'
   print *, 'base = ', base
   print *, 'Allowed range = [2.,Inf]'
   error = .TRUE.
endif

! max_ncell
if (max_ncell < 100000) then
   print *, 'ERROR: Invalid max_ncell value'
   print *, 'max_ncell = ', max_ncell
   print *, 'Allowed range = [100000,Inf]'
   error = .TRUE.
endif

! max_lvl
if (max_lvl < 1) then
   print *, 'ERROR: Invalid max_lvl value'
   print *, 'max_lvl = ', max_lvl
   print *, 'Allowed range = [1,Inf]'
   error = .TRUE.
endif

! min_lvl
if (min_lvl < 1) then
   print *, 'ERROR: Invalid min_lvl value'
   print *, 'min_lvl = ', min_lvl
   print *, 'Allowed range = [1,Inf]'
   error = .TRUE.
endif

! min_lvl < max_lvl
if (min_lvl > max_lvl) then 
   print *, 'ERROR: min_lvl > max_lvl!' 
   error = .TRUE.
endif

! max_dtau
if (max_dtau < 0) then
   print *, 'ERROR: Invalid max_dtau value'
   print *, 'max_dtau = ', max_dtau
   print *, 'Allowed range = [0,Inf]'
   error = .TRUE.
endif

! max_dlum
if (max_dlum < 0) then
   print *, 'ERROR: Invalid max_dlum value'
   print *, 'max_dlum = ', max_dlum
   print *, 'Allowed range = [0,Inf]'
   error = .TRUE.
endif

! z_subd_lim
if (z_subd_lim < 0.) then
   print *, 'ERROR: Invalid z_subd_lim value'
   print *, 'z_subd_lim =', z_subd_lim
   print *, 'Allowed range = (0.,Inf]'
   error = .TRUE.
endif

! R_subd_lim
if (R_subd_lim < 0.) then
   print *, 'ERROR: Invalid R_subd_lim value'
   print *, 'R_subd_lim =', R_subd_lim
   print *, 'Allowed range = (0.,Inf]'
   error = .TRUE.
endif

if (subdivision_criteria == 'milky_way') then 
   ! z_subd_lim2
   if (z_subd_lim2 < 0.) then
      print *, 'ERROR: Invalid z_subd_lim2 value'
      print *, 'z_subd_lim2 =', z_subd_lim2
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif
endif

if (grid_type == 'all' .or. grid_type == 'disk' .or. grid_type == 'bulge') then 
   ! old
   if (old < 0.) then
      print *, 'ERROR: Invalid old value'
      print *, 'old =', old
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif
endif

if (grid_type == 'all' .or. grid_type == 'disk') then 

   ! hs_disk_b
   if (hs_disk_b < 0.) then
      print *, 'ERROR: Invalid hs_disk_b value'
      print *, 'hs_disk_b =', hs_disk_b
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

   ! zs_disk
   if (zs_disk < 0.) then
      print *, 'ERROR: Invalid zs_disk value'
      print *, 'zs_disk =', zs_disk
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

   if (old_disk_type_ID == flared_expz_ID .or. old_disk_type_ID == flared_sech2z_ID .or. old_disk_type_ID == ellipt_expR_expz_ID .or. old_disk_type_ID == ellipt_expR_sech2z_ID) then 

      ! zs_disk_r1
      if (zs_disk_r1 < 0.) then
         print *, 'ERROR: Invalid zs_disk_r1 value'
         print *, 'zs_disk_r1 =', zs_disk_r1
         print *, 'Allowed range = (0.,Inf]'
         error = .TRUE.
      endif

      ! zs_disk_rsun
      if (zs_disk_rsun < 0.) then
         print *, 'ERROR: Invalid zs_disk_rsun value'
         print *, 'zs_disk_rsun =', zs_disk_rsun
         print *, 'Allowed range = (0.,Inf]'
         error = .TRUE.
      endif

      ! hsin
      if (hsin <= 0.) then
         print *, 'ERROR: Invalid hsin value'
         print *, 'hsin =', hsin
         print *, 'Allowed range = (0.,Inf]'
         error = .TRUE.
      endif

   endif

   if (old_disk_type_ID == flared_expz_ID .or. old_disk_type_ID == flared_sech2z_ID) then 

      ! chi_disk
      if (chi_disk == -100000) then
         print *, 'ERROR: chi_disk value not input'
         print *, 'chi_disk =', chi_disk
         print *, 'Allowed range = (-Inf,Inf]'
         error = .TRUE.
      endif

   endif

   if (old_disk_type_ID == ellipt_expR_expz_ID .or. old_disk_type_ID == ellipt_expR_sech2z_ID) then 

      ! theta_disk_ellipt 
      if (theta_disk_ellipt < 0.) then
         print *, 'ERROR: Invalid theta_disk_ellipt value'
         print *, 'theta_disk_ellipt =', theta_disk_ellipt
         print *, 'Allowed range = (0.,Inf]'
         error = .TRUE.
      endif
   endif

endif

if (grid_type == 'all' .or. grid_type == 'tdisk') then 
  
   ! sfr
   if (sfr < 0.) then
      print *, 'ERROR: Invalid sfr value'
      print *, 'sfr =', sfr
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

   ! hs_tdisk
   if (hs_tdisk < 0.) then
      print *, 'ERROR: Invalid hs_tdisk value'
      print *, 'hs_tdisk =', hs_tdisk
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

   if (young_disk_type_ID == ellipt_expR_expz_ID .or. young_disk_type_ID == ellipt_expR_sech2z_ID) then 

      ! hs_tdisk2
      if (hs_tdisk2 < 0.) then
         print *, 'ERROR: Invalid hs_tdisk2 value'
         print *, 'hs_tdisk2 =', hs_tdisk2
         print *, 'Allowed range = (0.,Inf]'
         error = .TRUE.
      endif

      ! theta_tdisk_ellipt 
      if (theta_tdisk_ellipt < 0.) then
         print *, 'ERROR: Invalid theta_tdisk_ellipt value'
         print *, 'theta_tdisk_ellipt =', theta_tdisk_ellipt
         print *, 'Allowed range = (0.,Inf]'
         error = .TRUE.
      endif

   endif

   ! zs_tdisk
   if (zs_tdisk < 0.) then
      print *, 'ERROR: Invalid zs_tdisk value'
      print *, 'zs_tdisk =', zs_tdisk
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

   if (young_disk_type_ID == flared_expz_ID .or. young_disk_type_ID == flared_sech2z_ID .or. young_disk_type_ID == ellipt_expR_expz_ID .or. young_disk_type_ID == ellipt_expR_sech2z_ID) then 


      ! zs_tdisk_r1
      if (zs_tdisk_r1 < 0.) then
         print *, 'ERROR: Invalid zs_tdisk_r1 value'
         print *, 'zs_tdisk_r1 =', zs_tdisk_r1
         print *, 'Allowed range = (0.,Inf]'
         error = .TRUE.
      endif

      ! zs_tdisk_rsun
      if (zs_tdisk_rsun < 0.) then
         print *, 'ERROR: Invalid zs_tdisk_rsun value'
         print *, 'zs_tdisk_rsun =', zs_tdisk_rsun
         print *, 'Allowed range = (0.,Inf]'
         error = .TRUE.
      endif

      ! hs1in
      if (hs1in <= 0.) then
         print *, 'ERROR: Invalid hs1in value'
         print *, 'hs1in =', hs1in
         print *, 'Allowed range = (0.,Inf]'
         error = .TRUE.
      endif

   endif

   if (young_disk_type_ID == flared_expz_ID .or. young_disk_type_ID == flared_sech2z_ID) then 
   
      ! chi_tdisk
      if (chi_tdisk == -100000) then
         print *, 'ERROR: chi_tdisk value not input'
         print *, 'chi_tdisk =', chi_tdisk
         print *, 'Allowed range = (-Inf,Inf]'
         error = .TRUE.
      endif
   
   endif
      
endif

if (grid_type == 'all' .or. grid_type == 'bulge') then 

   !reff 
   if (reff <= 0.) then
      print *, 'ERROR: Invalid reff value'
      print *, 'reff =', reff
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

   !acap_bulge 
   if (acap_bulge <= 0.) then
      print *, 'ERROR: Invalid acap_bulge value'
      print *, 'acap_bulge =', acap_bulge
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif
   
   !ellipt 
   if (ellipt <= 0) then
      print *, 'ERROR: Invalid ellipt value'
      print *, 'ellipt =', ellipt
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

   !mtrunc 
   if (mtrunc < 0) then
      print *, 'ERROR: Invalid mtrunc value'
      print *, 'mtrunc =', mtrunc
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

   !bd_ratio 
   if (bd_ratio < 0) then
      print *, 'ERROR: Invalid bd_ratio value'
      print *, 'bd_ratio =', bd_ratio
      print *, 'Allowed range = [0.,Inf]'
      error = .TRUE.
   endif

   !nsersic 
   if (nsersic < 1 .or. nsersic > 10 ) then
      print *, 'ERROR: Invalid nsersic value'
      print *, 'nsersic =', nsersic
      print *, 'Allowed range = [1.,10]'
      error = .TRUE.
   endif

   !theta_bulge 
   if (abs(theta_bulge) > 360) then
      print *, 'ERROR: Invalid theta_bulge value'
      print *, 'theta_bulge =', theta_bulge
      print *, 'Allowed range = [-360,360]'
      error = .TRUE.
   endif

   !ellipt_xy 
   if (ellipt_xy <= 0) then
      print *, 'ERROR: Invalid ellipt_xy value'
      print *, 'ellipt_xy =', ellipt_xy
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

endif


!tau1 
if (tau1 < 0) then
   print *, 'ERROR: Invalid tau1 value'
   print *, 'tau1 =', tau1
   print *, 'Allowed range = (0.,Inf]'
   error = .TRUE.
endif

!hd_disk 
if (hd_disk <= 0) then
   print *, 'ERROR: Invalid hd_disk value'
   print *, 'hd_disk =', hd_disk
   print *, 'Allowed range = (0.,Inf]'
   error = .TRUE.
endif

if (thick_disk_type_ID == ellipt_expR_expz_ID .or. thick_disk_type_ID == ellipt_expR_sech2z_ID) then

   !hd_disk2 
   if (hd_disk2 <= 0) then
      print *, 'ERROR: Invalid hd_disk2 value'
      print *, 'hd_disk2 =', hd_disk2
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

   !theta_dust_disk_ellipt 
   if (theta_dust_disk_ellipt < 0) then
      print *, 'ERROR: Invalid theta_dust_disk_ellipt value'
      print *, 'theta_dust_disk_ellipt =', theta_dust_disk_ellipt
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

endif


!zd_disk 
if (zd_disk <= 0) then
   print *, 'ERROR: Invalid zd_disk value'
   print *, 'zd_disk =', zd_disk
   print *, 'Allowed range = (0.,Inf]'
   error = .TRUE.
endif

if (thick_disk_type_ID == flared_expz_ID .or. thick_disk_type_ID == flared_sech2z_ID .or. thick_disk_type_ID == ellipt_expR_expz_ID .or. thick_disk_type_ID == ellipt_expR_sech2z_ID) then 

   !zd_disk_r1 
   if (zd_disk_r1 < 0) then
      print *, 'ERROR: Invalid zd_disk_r1 value'
      print *, 'zd_disk_r1 =', zd_disk_r1
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif
   
   !zd_disk_rsun 
   if (zd_disk_rsun < 0) then
      print *, 'ERROR: Invalid zd_disk_rsun value'
      print *, 'zd_disk_rsun =', zd_disk_rsun
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

   !hdin 
   if (hdin <= 0) then
      print *, 'ERROR: Invalid hdin value'
      print *, 'hdin =', hdin
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

endif

if (thick_disk_type_ID == flared_expz_ID .or. thick_disk_type_ID == flared_sech2z_ID) then 

   !chi_dust_disk 
   if (chi_dust_disk == -100000) then
      print *, 'ERROR: Invalid chi_dust_disk value'
      print *, 'chi_dust_disk =', chi_dust_disk
      print *, 'Allowed range = (-Inf,Inf]'
      error = .TRUE.
   endif
endif


!tau2 
if (tau2 < 0) then
   print *, 'ERROR: Invalid tau2 value'
   print *, 'tau2 =', tau2
   print *, 'Allowed range = (0.,Inf]'
   error = .TRUE.
endif

!hd_tdisk 
if (hd_tdisk <= 0) then
   print *, 'ERROR: Invalid hd_tdisk value'
   print *, 'hd_tdisk =', hd_tdisk
   print *, 'Allowed range = (0.,Inf]'
   error = .TRUE.
endif

if (thin_disk_type_ID == ellipt_expR_expz_ID .or. thin_disk_type_ID == ellipt_expR_sech2z_ID) then

   !hd_tdisk2 
   if (hd_tdisk2 <= 0) then
      print *, 'ERROR: Invalid hd_tdisk2 value'
      print *, 'hd_tdisk2 =', hd_tdisk2
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

   !theta_dust_tdisk_ellipt 
   if (theta_dust_tdisk_ellipt < 0) then
      print *, 'ERROR: Invalid theta_dust_tdisk_ellipt value'
      print *, 'theta_dust_tdisk_ellipt =', theta_dust_tdisk_ellipt
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

endif

!zd_tdisk 
if (zd_tdisk <= 0) then
   print *, 'ERROR: Invalid zd_tdisk value'
   print *, 'zd_tdisk =', zd_tdisk
   print *, 'Allowed range = (0.,Inf]'
   error = .TRUE.
endif

if (thin_disk_type_ID == flared_expz_ID .or. thin_disk_type_ID == flared_sech2z_ID .or. thin_disk_type_ID == ellipt_expR_expz_ID .or. thin_disk_type_ID == ellipt_expR_sech2z_ID) then 

   !zd_tdisk_r1 
   if (zd_tdisk_r1 < 0) then
      print *, 'ERROR: Invalid zd_tdisk_r1 value'
      print *, 'zd_tdisk_r1 =', zd_tdisk_r1
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

   !zd_tdisk_rsun 
   if (zd_tdisk_rsun < 0) then
      print *, 'ERROR: Invalid zd_tdisk_rsun value'
      print *, 'zd_tdisk_rsun =', zd_tdisk_rsun
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

   !hdin 
   if (hd1in <= 0) then
      print *, 'ERROR: Invalid hd1in value'
      print *, 'hdin =', hd1in
      print *, 'Allowed range = (0.,Inf]'
      error = .TRUE.
   endif

endif

if (thin_disk_type_ID == flared_expz_ID .or. thin_disk_type_ID == flared_sech2z_ID) then 

   !chi_dust_tdisk 
   if (chi_dust_tdisk == -100000) then
      print *, 'ERROR: Invalid chi_dust_tdisk value'
      print *, 'chi_dust_tdisk =', chi_dust_tdisk
      print *, 'Allowed range = (-Inf,Inf]'
      error = .TRUE.
   endif

endif


if (error) then 
   print *, 'STOP(check_input_galaxy): Some wrong input variables!'
   STOP
endif


end subroutine check_input_galaxy


!> Checks that lambda_arr contains at least two wavelengths of which one is equal to 0.443 um (B-band). This wavelength is necessary for calc_scaling_factors_dust() to work properly. 
subroutine check_lambda_list_galaxy

if (lnum_tot < 2) then 
   print *, 'STOP(read_input_galaxy): At least two wavelengths needed (of which one equal to 0.443um)!'
   stop
endif

if (count(abs(lambda_arr-0.443)/0.443 < 1E-3) == 0) then 
   print *, 'STOP(read_input_galaxy): Add lambda = 0.443um to wavelength list. This is necessary for proper scaling of opacity factors!'
   stop
endif

end subroutine check_lambda_list_galaxy

  
subroutine check_grid_type
  
  if ((grid_type /= 'all').and.(grid_type /= 'disk').and.(grid_type /= 'tdisk').and.(grid_type /= 'bulge')) then
     print *, 'grid_type not recognized =', grid_type
     stop
  endif

end subroutine check_grid_type

!> Calculates luminosity of the different stellar components depending on wavelength and intensity parameters old and sfr. The unit luminosities are those from Popescu et al. 2011 
subroutine calc_luminosities_pt11
  integer, parameter :: nlambda_pt11 = 14
  ! wavelength array in PT11. 4429.999 instead of 4430 to avoid wrong interpolation if lambda = 4430 for the old stellar disk. Same for 912. 
  real(kind=real64), parameter :: lambda_arr_pt11(0:nlambda_pt11-1) = (/911.999, 1350., 1500., 1650., 2000., 2200., 2500., 2800., 3650., 4429.999, 5640., 8090., 12590., 22E3/)*1E-10 ! A -> m
  real(kind=real64), parameter :: lnu_old_unit(0:nlambda_pt11-1) = (/0.,0.,0.,0.,0.,0.,0.,0.,0., 4.771E21, 9.382E21 , 19.54E21  , 72.20E21, 64.97E21/)  ! W/Hz 
  real(kind=real64), parameter :: lnu_sf_unit(0:nlambda_pt11-1) = (/0.344E21, 0.905E21,0.844E21, 0.863E21, 0.908E21, 0.926E21, 0.843E21 , 0.910E21 , 1.842E21, 2.271E21, 3.837E21, 5.734E21 , 0.931E21 , 0.728E21/) 

 
  integer :: i0, i1
  real(kind=real64) :: m, c 

  ! Initialization
  lnu_old=0.
  lnu_sf=0.
  lnu_bulge=0.

  if ((grid_type == 'all').or.(grid_type == 'disk').or.(grid_type == 'bulge')) then    
!!! OLD stellar disk 
!!! luminosity values from Popescu et al. 2011
     if (lambda < 4429.999E-10) then
        lnu_old =0.      
     elseif (lambda>= 4429.999E-10.and. lambda <= 2.2001E-6) then 
        call value_locate(lambda, lambda_arr_pt11, i0,i1, .FALSE.)
        m = (lnu_old_unit(i1)-lnu_old_unit(i0))/(lambda_arr_pt11(i1)- lambda_arr_pt11(i0))
        c = lnu_old_unit(i1) - m*lambda_arr_pt11(i1)
        lnu_old = m*lambda + c 
     elseif (lambda > 2.2E-6) then 
        lnu_old=lnu_old_unit(nlambda_pt11-1)*(2.2E-6/lambda)**2.   !!! longer wavelength (K_lamba/lambda)^2
     endif

     lnu_old=old*lnu_old

  endif

if ((grid_type == 'all').or.(grid_type == 'tdisk')) then

   if (lambda <= 911.999E-10) then 
      lnu_sf =0.
   elseif (lambda>= 911.999E-10 .and. lambda <= 2.2001E-6) then 
        call value_locate(lambda, lambda_arr_pt11, i0,i1, .FALSE.)
        m = (lnu_sf_unit(i1)-lnu_sf_unit(i0))/(lambda_arr_pt11(i1)- lambda_arr_pt11(i0))
        c = lnu_sf_unit(i1) - m*lambda_arr_pt11(i1)
        lnu_sf = m*lambda + c 
     elseif (lambda > 2.2E-6) then    
      lnu_sf=lnu_sf_unit(nlambda_pt11-1)*(2.2E-6/lambda)**2.    !!! longer wavelength (K_lamba/lambda)^2
   endif

   lnu_sf=sfr*lnu_sf

endif

if ((grid_type == 'all').or.(grid_type == 'bulge')) then 
   lnu_bulge=bd_ratio*lnu_old

endif

!print *, lnu_old,lnu_sf,lnu_bulge
lnu_tot=lnu_old+lnu_sf+lnu_bulge


end subroutine calc_luminosities_pt11

!> Calculates luminosity of the different stellar components depending on wavelength and intensity parameters old and sfr. The unit luminosities are those from the input files file_old_star_sed() and file_young_star_sed(). 
subroutine calc_luminosities
!!$  !> Number of wavelengths for the old stellar SED
!!$  integer :: nlambda_old
!!$  !> Number of wavelengths for the young stellar SED
!!$  integer :: nlambda_old
  ! wavelength array in PT11. 4429.999 instead of 4430 to avoid wrong interpolation if lambda = 4430 for the old stellar disk. Same for 912. 
!!$  real(kind=real64), allocatable :: lambda_arr_old_sed(:)  ! um
!!$  real(kind=real64), allocatable :: lambda_arr_sf_sed(:)  ! um
!!$  real(kind=real64), allocatable :: lnu_old_unit(:)   ! W/Hz 
!!$  real(kind=real64), allocatable :: lnu_sf_unit(:)  ! W/Hz  
  integer :: i0, i1
  real(kind=real64) :: m, c 

  ! Initialization
  lnu_old=0.
  lnu_sf=0.
  lnu_bulge=0.

  if ((grid_type == 'all').or.(grid_type == 'disk').or.(grid_type == 'bulge')) then    
!!! OLD stellar disk 
!!! luminosity values from file_old_star_sed 
 
     if (lambda < lambda_arr_old_sed(0)*0.999) then
        lnu_old =0.      
     elseif (lambda>= lambda_arr_old_sed(0)*0.999 .and. lambda <= lambda_arr_old_sed(nlambda_old_sed-1)) then 
        call value_locate(lambda, lambda_arr_old_sed, i0,i1, .FALSE.)
        m = (lnu_old_unit(i1)-lnu_old_unit(i0))/(lambda_arr_old_sed(i1)- lambda_arr_old_sed(i0))
        c = lnu_old_unit(i1) - m*lambda_arr_old_sed(i1)
        lnu_old = m*lambda + c 
     elseif (lambda > lambda_arr_old_sed(nlambda_old_sed-1)) then 
        lnu_old=lnu_old_unit(nlambda_old_sed-1)*(lambda_arr_old_sed(nlambda_old_sed-1)/lambda)**2.   !!! longer wavelength (last_lambba/lambda)^2
     endif

     lnu_old=old*lnu_old

  endif

if ((grid_type == 'all').or.(grid_type == 'tdisk')) then
! young stellar disk

   if (lambda <= lambda_arr_sf_sed(0)*0.999) then 
      lnu_sf =0.
   elseif (lambda>= lambda_arr_sf_sed(0)*0.999 .and. lambda <= lambda_arr_sf_sed(nlambda_sf_sed-1)) then 
        call value_locate(lambda, lambda_arr_sf_sed, i0,i1, .FALSE.)
        m = (lnu_sf_unit(i1)-lnu_sf_unit(i0))/(lambda_arr_sf_sed(i1)- lambda_arr_sf_sed(i0))
        c = lnu_sf_unit(i1) - m*lambda_arr_sf_sed(i1)
        lnu_sf = m*lambda + c 
     elseif (lambda > lambda_arr_sf_sed(nlambda_sf_sed-1)) then    
      lnu_sf=lnu_sf_unit(nlambda_sf_sed-1)*(lambda_arr_sf_sed(nlambda_sf_sed-1)/lambda)**2.    !!! longer wavelength (last_lambda/lambda)^2
   endif

   lnu_sf=sfr*lnu_sf

endif

if ((grid_type == 'all').or.(grid_type == 'bulge')) then 
   lnu_bulge=bd_ratio*lnu_old

endif

!print *, lnu_old,lnu_sf,lnu_bulge
lnu_tot=lnu_old+lnu_sf+lnu_bulge


end subroutine calc_luminosities



  subroutine calc_scaling_factors_stars
    ! This subroutines calculates the appropriate scaling factors for the stellar and dust components
    real(kind=real64) :: termr, termrin, termr1, termr1in

    !!! DISK types 
    !! 'exp' : pure double exponential
    !! 'nosdgrad' : flat inner galaxy until Rin and then exponential
    !! 'grad05'   : linear inner gradient such that f(R=0)=0.5*f(R=Rin), then exponential
    !! 'flared_grad05' : as grad05 but with flare    
    !! 'ellipt_exp'
    
    ! parameter initialization
    eta_disk0=0.
    eta_tdisk0=0.
    eta_bulge0=0. 
    
! Old stellar disk
    if ((grid_type == 'all').or.(grid_type == 'disk')) then 

    if (old_disk_type_ID == expR_expz_ID .or. old_disk_type_ID == expR_sech2z_ID.or. old_disk_type_ID == ellipt_expR_expz_ID .or. old_disk_type_ID == ellipt_expR_sech2z_ID) then  ! for the elliptical exponential profiles this is not exactly the right scaling but it should be close. the scaling is fixed at the end of the grid calculation.  

       termr = 1.-exp(-rtrun/hs_disk)-(rtrun/hs_disk)* exp(-rtrun/hs_disk)
       eta_disk0=lnu_old/(4*pi*zs_disk*hs_disk**2*termr)

!!$    else if (old_disk_type == 'nosdgrad') then 
!!$       termr = exp(-hsin/hs_disk)-exp(-rtrun/hs_disk)+(hsin/hs_disk)*exp(-hsin/hs_disk) &
!!$        -(rtrun/hs_disk)*exp(-rtrun/hs_disk)
!!$       termrin = 2*pi*hs1in**2*zs_disk*exp(-hs1in/hs_disk)  ! constant profile
!!$
!!$       eta_disk0=lnu_old/(4*pi*zs_disk*hs_disk**2*termr + termrin) 
!!$
!!$    else if ((old_disk_type == 'grad05').or.(old_disk_type == 'flared_grad05')) then 
!!$
!!$       termr = exp(-hsin/hs_disk)-exp(-rtrun/hs_disk)+(hsin/hs_disk)*exp(-hsin/hs_disk)&
!!$        -(rtrun/hs_disk)*exp(-rtrun/hs_disk)
!!$
!!$       termrin = 5./3.*pi*hsin**2*zs_disk*exp(-hsin/hs_disk)  ! grad 0.5 
!!$       eta_disk0=lnu_old/(4*pi*zs_disk*hs_disk**2*termr + termrin)

    else if (old_disk_type_ID == flared_expz_ID .or. old_disk_type_ID == flared_sech2z_ID ) then 

       termr = exp(-hsin/hs_disk)-exp(-rtrun/hs_disk)+(hsin/hs_disk)*exp(-hsin/hs_disk)&
        -(rtrun/hs_disk)*exp(-rtrun/hs_disk)

       termrin = 4./3.*(1+chi_disk/2)*pi*hsin**2*zs_disk*exp(-hsin/hs_disk)  ! chi factor
       
       eta_disk0=lnu_old/(4*pi*zs_disk*hs_disk**2*termr + termrin)

!!$    else if (old_disk_type_ID == ellipt_expR_expz_ID .or. old_disk_type_ID == ellipt_expR_sech2z_ID) then 
!!$
!!$      eta_disk0 = 1 

    else

       print *, 'old_disk_type not recognized: ', old_disk_type_ID
       stop 
       
    endif

 endif


 if ((grid_type == 'all').or.(grid_type == 'tdisk')) then 
    ! YOUNG stellar disc 
    if (young_disk_type_ID == expR_expz_ID .or. young_disk_type_ID == expR_sech2z_ID) then
       termr1 = 1.-exp(-rtrun/hs_tdisk)-(rtrun/hs_tdisk)* exp(-rtrun/hs_tdisk)
       eta_tdisk0=lnu_sf/(4*pi*zs_tdisk*hs_tdisk**2*termr1)

!!$    else if (young_disk_type == 'nosdgrad') then 
!!$   
!!$       termr1 = exp(-hs1in/hs_tdisk)-exp(-rtrun/hs_tdisk)+(hs1in/hs_tdisk)*exp(-hs1in/hs_tdisk)&
!!$            -(rtrun/hs_tdisk)*exp(-rtrun/hs_tdisk)
!!$       termr1in = 2*pi*hs1in**2*zs_tdisk*exp(-hs1in/hs_tdisk)  ! constant profile
!!$       eta_tdisk0=lnu_sf/(4*pi*zs_tdisk*hs_tdisk**2*termr1 + termr1in )
!!$
!!$    else if ((young_disk_type == 'grad05').or.(young_disk_type == 'flared_grad05')) then
!!$       termr1 = exp(-hs1in/hs_tdisk)-exp(-rtrun/hs_tdisk)+(hs1in/hs_tdisk)*exp(-hs1in/hs_tdisk)&
!!$            -(rtrun/hs_tdisk)*exp(-rtrun/hs_tdisk)
!!$
!!$       termr1in = 5./3.*pi*hs1in**2*zs_tdisk*exp(-hs1in/hs_tdisk)  ! grad 0.5 
!!$       eta_tdisk0=lnu_sf/(4*pi*zs_tdisk*hs_tdisk**2*termr1 + termr1in )  

   else if (young_disk_type_ID == flared_expz_ID .or. young_disk_type_ID == flared_sech2z_ID) then
       termr1 = exp(-hs1in/hs_tdisk)-exp(-rtrun/hs_tdisk)+(hs1in/hs_tdisk)*exp(-hs1in/hs_tdisk)&
            -(rtrun/hs_tdisk)*exp(-rtrun/hs_tdisk)

       termr1in = 4./3.*(1+chi_tdisk/2)*pi*hs1in**2*zs_tdisk*exp(-hs1in/hs_tdisk)  ! chi factor 
       eta_tdisk0=lnu_sf/(4*pi*zs_tdisk*hs_tdisk**2*termr1 + termr1in )     

    else if (young_disk_type_ID == ellipt_expR_expz_ID .or. young_disk_type_ID == ellipt_expR_sech2z_ID) then 

       eta_tdisk0 = 1 

    else
       print *, 'young_disk_type not recognized: ', young_disk_type_ID
       stop 
    endif

 endif

    if ((grid_type == 'all').or.(grid_type == 'bulge')) then 
      eta_bulge0=1  !!! this is scaled after grid creation 
   endif

    end subroutine calc_scaling_factors_stars

subroutine calc_scaling_factors_dust
! This calculates the scaling factor for the dust distributions (assuming double exponential shape and face-on total optical depth in the B-band speficied in the input )
real(kind=real64) :: taub_f, tau_ratio

taub_f=tau1+tau2

if (tau1 < 0 .or. tau2 < 0) then 
print*, 'tau cannot be negative!!!'
stop
endif 

lambda=4430.E-10 !!! This is because taus are for the B-band 
call interpolate_kext
taub_f=taub_f/kext  

lambda=lambda_in
call interpolate_kext
taub_f=taub_f*kext    !!! This is the face-on optical depth at lambda_in

if (tau2 > 0) then 
tau_ratio=tau1/tau2

kext_tdisk0=taub_f/(2*(1+tau_ratio)*zd_tdisk)
kext_disk0=tau_ratio*kext_tdisk0*zd_tdisk/zd_disk

else 

kext_tdisk0=0
kext_disk0=taub_f/(2*zd_disk)

endif  

end subroutine calc_scaling_factors_dust

!> Interpolates value of kext given wavelength lambda. 
subroutine interpolate_kext
real(kind=real64) :: m, c 
integer :: i0, i1 

call value_locate(lambda, lambda_arr_SI, i0,i1, .FALSE.)
m = (kext_arr(i1)-kext_arr(i0))/(lambda_arr_SI(i1)- lambda_arr_SI(i0))
c = kext_arr(i1) - m*lambda_arr_SI(i1)
kext = m*lambda + c 

end subroutine interpolate_kext



real(kind=real64) function av_disk(x,y,z,cellsize,disk_comp,disk_type_ID)
IMPLICIT NONE
real(Kind=real64) :: x,y,z,cellsize,dx,tot_sum,sum_temp,xi,yj,zk,rad
real(Kind=real64) :: x0,y0,z0,xexp,ds,rterm,zterm,zterm1,zc1,const, zc_flat_grad, norm_rad
real(kind=real64) :: a,hc,hc2,r1,zc,zc_r1,zc_rsun,zc_r,chi_par
real(kind=real64) :: theta_ellipt, cos_theta_el, sin_theta_el,xi_rot, yj_rot
integer :: steps,i,j,k,nt
character(LEN=lcar_type) :: disk_comp 
integer :: disk_type_ID

! load disk_parameters depending con disk_component and disk_type 

if (disk_comp == 'disk') then 
   a=eta_disk0
   hc=hs_disk   
   zc=zs_disk

   if (disk_type_ID == flared_expz_ID .or. disk_type_ID == flared_sech2z_ID ) then 
      zc_r1=zs_disk_r1
      zc_rsun=zs_disk_rsun
      chi_par = chi_disk
      r1=hsin
   endif

   if (disk_type_ID == ellipt_expR_expz_ID .or. disk_type_ID == ellipt_expR_sech2z_ID) then 
      hc2=hs_disk2
      zc_r1=zs_disk_r1
      zc_rsun=zs_disk_rsun
      r1=hsin
      theta_ellipt = theta_disk_ellipt
   endif

elseif (disk_comp == 'tdisk') then 
   a=eta_tdisk0
   hc=hs_tdisk   
   zc=zs_tdisk

   if (disk_type_ID == flared_expz_ID .or. disk_type_ID == flared_sech2z_ID) then 
      zc_r1=zs_tdisk_r1
      zc_rsun=zs_tdisk_rsun
      chi_par = chi_tdisk
      r1=hs1in
   endif

   if (disk_type_ID == ellipt_expR_expz_ID .or. disk_type_ID == ellipt_expR_sech2z_ID) then 
      hc2=hs_tdisk2
      zc_r1=zs_tdisk_r1
      zc_rsun=zs_tdisk_rsun
      r1=hs1in
      theta_ellipt = theta_tdisk_ellipt
   endif

elseif (disk_comp == 'dust_disk') then 
   a=kext_disk0
   hc=hd_disk   
   zc=zd_disk

   if (disk_type_ID == flared_expz_ID .or. disk_type_ID == flared_sech2z_ID) then 
      zc_r1=zd_disk_r1
      zc_rsun=zd_disk_rsun
      chi_par = chi_dust_disk
      r1=hdin
   endif

   if (disk_type_ID == ellipt_expR_expz_ID .or. disk_type_ID == ellipt_expR_sech2z_ID) then 
      hc2=hd_disk2
      zc_r1=zd_disk_r1
      zc_rsun=zd_disk_rsun
      r1=hdin
      theta_ellipt = theta_dust_disk_ellipt
   endif

else if (disk_comp == 'dust_tdisk') then 
   a=kext_tdisk0
   hc=hd_tdisk   
   zc=zd_tdisk

   if (disk_type_ID == flared_expz_ID .or. disk_type_ID == flared_sech2z_ID) then 
      zc_r1=zd_tdisk_r1
      zc_rsun=zd_tdisk_rsun
      chi_par = chi_dust_tdisk
      r1=hd1in
   endif

   if (disk_type_ID == ellipt_expR_expz_ID .or. disk_type_ID == ellipt_expR_sech2z_ID) then 
      hc2=hd_tdisk2
      zc_r1=zd_tdisk_r1
      zc_rsun=zd_tdisk_rsun
      r1=hd1in
      theta_ellipt = theta_dust_tdisk_ellipt
   endif

  
else 

print *, 'disk component non recognized:', disk_comp

stop
endif


if (a == 0) then 
av_disk=0
return
endif

if (disk_type_ID == flared_expz_ID .or. disk_type_ID == flared_sech2z_ID ) then 
   if (abs(zc_r1-zc) > 1E-5*zc) then
      xexp=log((zc_rsun-zc)/(zc_r1-zc))/(log(rsun/r1))
   else 
      xexp=0.
   endif
   
endif

if (disk_type_ID == ellipt_expR_expz_ID .or. disk_type_ID == ellipt_expR_sech2z_ID) then 

cos_theta_el = cos(theta_ellipt*pi/180.)  
sin_theta_el = sin(theta_ellipt*pi/180.)

endif


dx=cellsize
x0=x-dx/2.
y0=y-dx/2.
z0=z-dx/2.

steps=step_int
ds=dx/steps

tot_sum=0

nt=nproc
if (grid_creation_lambda) nt = 1 

! to do: here you should insert new formula for dust disk with flat centre but no taper (because this is the thick not the thin dust disk) 

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,j,k,sum_temp,xi,yj,zk,rad,zc_r,zc1,rterm,zterm,zterm1,const, norm_rad, xi_rot, yj_rot), &
!$OMP SHARED(steps,x0,ds,y0,z0,tot_sum,xexp,max_rad,max_z,rsun,zc_r1,r1,zc,hc,hc2,disk_type_ID, zc_flat_grad,disk_comp,chi_par,cos_theta_el, sin_theta_el), NUM_THREADS(nt)

!$OMP DO SCHEDULE(DYNAMIC,1)

do i=0,steps-1

   sum_temp=0

   do j=0,steps-1
      do k=0,steps-1

         !print *,x0+i*ds
         xi=x0+(i+0.5)*ds
         yj=y0+(j+0.5)*ds
         zk=z0+(k+0.5)*ds

         zk=abs(zk)  !!! valid only for axisymmetric models 

         if (zk > max_z) cycle

         rad=sqrt(xi**2+yj**2)

         if (rad > max_rad) cycle 

         if (disk_type_ID == expR_expz_ID) then  

            sum_temp=sum_temp+exp(-rad/hc)*exp(-abs(zk)/zc)

         else if (disk_type_ID == expR_sech2z_ID) then 

            sum_temp=sum_temp+exp(-rad/hc)*(cosh(-(zk)/zc))**(-2) ! Hyperbolic secant

!!$         else if (disk_type == 'grad05') then 
!!$
!!$            if (rad >= r1) then 
!!$      
!!$               sum_temp=sum_temp+exp(-rad/hc)*(cosh(-(zk)/zc))**(-2) ! Hyperbolic secant
!!$
!!$            else 
!!$
!!$               sum_temp=sum_temp+(0.5*rad/r1+0.5)*exp(-r1/hc)*(cosh(-(zk)/zc))**(-2)  ! gradient in the centre 0.5 + Hyperbolic secant
!!$
!!$            endif
            
!!$         elseif (disk_type == 'nosdgrad') then 
!!$
!!$            if (rad >= r1) then 
!!$
!!$               sum_temp=sum_temp+exp(-rad/hc)*(cosh(-(zk)/zc))**(-2) ! Hyperbolic secant
!!$
!!$            else 
!!$
!!$               sum_temp=sum_temp+exp(-r1/hc)*(cosh(-(zk)/zc))**(-2)
!!$
!!$            endif

!!$         else if (disk_type == 'flared_grad05') then 
!!$
!!$          !  xexp=log((zc_rsun-zc)/(zc_r1-zc))/(log(rsun/r1))
!!$            zc_r=zc+(zc_r1-zc)*(rad/r1)**xexp    ! flare (increasing scale height)
!!$           
!!$            if (rad >= r1) then
!!$               sum_temp=sum_temp+zc/zc_r*(exp(-rad/hc))*(cosh(-(zk)/zc_r))**(-2)  
!!$            else 
!!$               sum_temp=sum_temp+(0.5*rad/r1+0.5)*zc/zc_r*(exp(-r1/hc))*(cosh(-(zk)/zc_r))**(-2)
!!$               
!!$            endif

         else if (disk_type_ID == flared_expz_ID) then 

          !  xexp=log((zc_rsun-zc)/(zc_r1-zc))/(log(rsun/r1))
            zc_r=zc+(zc_r1-zc)*(rad/r1)**xexp    ! flare (increasing scale height)
           
            if (rad >= r1) then
               sum_temp=sum_temp+zc/zc_r*(exp(-rad/hc))*exp(-abs(zk)/zc_r)   
            else 
               sum_temp=sum_temp+(rad/r1*(1-chi_par) + chi_par)*zc/zc_r*(exp(-r1/hc))*exp(-abs(zk)/zc_r)              
            endif

            if (chi_par < 0. .and. sum_temp < 0) sum_temp = 0 ! in case of negative chi_par, the values could become negative at small radii 

         else if (disk_type_ID == flared_sech2z_ID) then 

            !  xexp=log((zc_rsun-zc)/(zc_r1-zc))/(log(rsun/r1))
            zc_r=zc+(zc_r1-zc)*(rad/r1)**xexp    ! flare (increasing scale height)
           
            if (rad >= r1) then
               sum_temp=sum_temp+zc/zc_r*(exp(-rad/hc))*(cosh(-(zk)/zc_r))**(-2)  
            else 
               sum_temp=sum_temp+(rad/r1*(1-chi_par) + chi_par)*zc/zc_r*(exp(-r1/hc))*(cosh(-(zk)/zc_r))**(-2)              
            endif

            if (chi_par < 0. .and. sum_temp < 0) sum_temp = 0 ! in case of negative chi_par, the values could become negative at small radii 
            

!!$         else if (disk_type == 'grad05_zflat') then 
!!$
!!$            if (rad >= r1) then 
!!$
!!$               rterm=exp(-rad/hc)
!!$
!!$            else 
!!$
!!$               rterm=(0.5*rad/r1+0.5)*exp(-r1/hc)
!!$
!!$            endif
!!$            
!!$            zc1= zc_flat_grad*rad
!!$            zterm1=exp(-zc1/zc)
!!$            const=zc/((zc1+zc)*zterm1)
!!$
!!$            if (zk >= zc1) then
!!$
!!$               zterm = const*exp(-zk/zc)
!!$
!!$            else
!!$
!!$               zterm = const*zterm1
!!$
!!$            endif
!!$            
!!$            sum_temp=sum_temp+rterm*zterm

         else if (disk_type_ID == ellipt_expR_expz_ID) then

            zc_r=zc+(zc_r1-zc)*(rad/r1)**xexp    ! flare (increasing scale height)

            if (rad >= r1) then

               sum_temp = sum_temp
               
            else

               norm_rad = sqrt((xi/hc)**2+(yj/hc2)**2)
               sum_temp = sum_temp + zc/zc_r*(exp(-norm_rad))*exp(-abs(zk)/zc_r) 

            endif

         else if (disk_type_ID == ellipt_expR_sech2z_ID) then

            zc_r=zc+(zc_r1-zc)*(rad/r1)**xexp    ! flare (increasing scale height)

            if (rad >= r1) then

               sum_temp = sum_temp
               
            else

               xi_rot = xi*cos_theta_el + yj*sin_theta_el ! coordinate transformation
               yj_rot = -xi*sin_theta_el + yj*cos_theta_el
               norm_rad = sqrt((xi_rot/hc)**2+(yj_rot/hc2)**2)
               sum_temp = sum_temp + zc/zc_r*(exp(-norm_rad))*(cosh(-(zk)/zc_r))**(-2) 

            endif

                        
         else 

            print *,' disk_type_ID not recognized:',disk_type_ID
            stop

         endif

      enddo
   enddo

   !$OMP ATOMIC
   tot_sum=tot_sum+sum_temp 

enddo

!$OMP END DO NOWAIT
!$OMP END PARALLEL

tot_sum=tot_sum*ds**3

av_disk=a*tot_sum/cellsize**3

end function av_disk


real(kind=real64) function av_star_bulge(x,y,z,cellsize)
real(Kind=real64) :: x,y,z,cellsize,a,br,bz,dx,x0,y0,z0,ds,f_exp,sum1,xi,yj,zk,sum2,z1,rad,bsersic,ellipt2,n,acap,m,sum_temp !m is the small r in Richard's formula
integer :: steps,i,j,k,nt  
real(Kind=real64) :: cos_theta_bulge, sin_theta_bulge, xi_rot, yj_rot
real(Kind=real64),parameter,dimension(10) :: bsersic_arr=(/1.67835 , 3.67206, 5.67017, &
     7.66925, 9.66872, 11.6684, 13.6681, 15.6679, 17.6678,19.6677/)

a=eta_bulge0

if (a == 0) then 
av_star_bulge=0
return
endif 


if (nsersic > 10) STOP 'nsersic too high' 

bsersic=bsersic_arr(nsersic)

ellipt2=ellipt**2

acap=acap_bulge/reff

n=real(nsersic)

cos_theta_bulge = cos(theta_bulge*pi/180)
sin_theta_bulge = sin(theta_bulge*pi/180)

dx=cellsize
x0=x-dx/2.
y0=y-dx/2.
z0=z-dx/2.

steps=step_int
ds=dx/steps

sum1=0

nt=nproc
if (grid_creation_lambda) nt = 1

! to do: here you should insert new formula for dust disk with flat centre and taper (because this is the thin dust disk) 

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP PRIVATE(i,j,k,sum_temp,xi,yj,zk,rad,m,xi_rot, yj_rot), &
!$OMP SHARED(steps,x0,ds,y0,z0,sum1,bsersic,n,acap,ellipt2,reff,mtrunc,max_rad,cos_theta_bulge, sin_theta_bulge,ellipt_xy), NUM_THREADS(nt)

!$OMP DO SCHEDULE(DYNAMIC,1)

do i=0,steps-1

sum_temp=0

do j=0,steps-1
do k=0,steps-1

!print *,x0+i*ds
xi=x0+(i+0.5)*ds
yj=y0+(j+0.5)*ds
zk=z0+(k+0.5)*ds

!if (abs(zk) > max_z) cycle

xi_rot = xi*cos_theta_bulge + yj*sin_theta_bulge ! coordinate transformation
yj_rot = -xi*sin_theta_bulge + yj*cos_theta_bulge

rad=sqrt(xi_rot**2+(yj_rot/ellipt_xy)**2)
!rad=sqrt(xi**2+yj**2)

if (rad > max_rad) cycle 

m=((rad**2 + (zk*zk/ellipt2))**0.5)/ reff

if (m < acap) m=acap

! bulge is truncated at mtrunc effective radii
if (m <= mtrunc) sum_temp=sum_temp+ (1.0/m**((2*n-1)/(2*n))) * exp(-bsersic*m**(1./n))

enddo
enddo

!$OMP ATOMIC
 sum1=sum1+sum_temp 

enddo

!$OMP END DO NOWAIT
!$OMP END PARALLEL

sum1=sum1*ds**3

av_star_bulge=a*sum1/cellsize**3

end function av_star_bulge

subroutine create_grid_arrays_galaxy 

if (tot_ncell == 0) then 
allocate(dens_disk(0:max_ncell-1), dens_tdisk(0:max_ncell-1), dens_bulge(0:max_ncell-1), dens_dust_disk(0:max_ncell-1),dens_dust_tdisk(0:max_ncell-1))
else 
allocate(dens_disk(0:tot_ncell-1), dens_tdisk(0:tot_ncell-1), dens_bulge(0:tot_ncell-1), dens_dust_disk(0:tot_ncell-1),dens_dust_tdisk(0:tot_ncell-1))
endif 

dens_disk=0
dens_tdisk=0
dens_bulge=0
dens_dust_disk=0
dens_dust_tdisk=0

end subroutine create_grid_arrays_galaxy 


subroutine fix_dens_stars_arrays
! this routine checks whether total luminosities are calculated correctly and sum up the dens_stars arrays 

real(kind=real64) :: tot_disk, tot_tdisk, tot_bulge

tot_disk=sum(dens_disk*csize**3,mask=(cchild == -1 .and. ncell <= tot_ncell ))

if (old_disk_type_ID == ellipt_expR_expz_ID .or. old_disk_type_ID == ellipt_expR_sech2z_ID) then
   if (tot_disk > 0) then 
      dens_disk= dens_disk/tot_disk*lnu_old
      tot_disk = lnu_old
   endif
endif

if (tot_disk > 0.) then 
print *, 'accuracy Old disk =', abs(tot_disk-lnu_old)/lnu_old
else 
print *, 'No Old stellar disk'
endif  

tot_tdisk=sum(dens_tdisk*csize**3,mask=(cchild == -1 .and. ncell <= tot_ncell))

if (young_disk_type_ID == ellipt_expR_expz_ID .or. young_disk_type_ID == ellipt_expR_sech2z_ID) then 
   if (tot_tdisk > 0) then 
      dens_tdisk= dens_tdisk/tot_tdisk*lnu_sf
      tot_tdisk = lnu_sf
   endif
endif

if (tot_tdisk > 0.) then 
print *, 'accuracy Young disk =', abs(tot_tdisk-lnu_sf)/lnu_sf
else 
print *, 'No Young stellar disk'
endif  

tot_bulge=sum(dens_bulge*csize**3,mask=(cchild == -1 .and. ncell <= tot_ncell))
 
if (tot_bulge > 0.) then 

dens_bulge=dens_bulge*lnu_bulge/tot_bulge

endif 

dens_stars=dens_disk+dens_tdisk+dens_bulge

print *, ''

end subroutine fix_dens_stars_arrays

subroutine print_3d_grid_file_galaxy
    
     INTEGER(HID_T) :: file_id                            ! File identifier
     INTEGER(HID_T) :: dset_id       ! Dataset identifier
     INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
     integer, parameter :: narr =5
     INTEGER     ::   rank(0:narr-1)                         ! Dataset rank
     INTEGER     ::   error  ! Error flag
     INTEGER(HSIZE_T) ::  dims(1)
     character(LEN=lcar_type) :: dsetname(0:narr-1)
     integer :: i

     

     dsetname(0)='dens_disk'     ; rank(0)=1 
     dsetname(1)='dens_tdisk'    ; rank(1)=1 
     dsetname(2)='dens_bulge'    ; rank(2)=1 
     dsetname(3)='dens_dust_disk' ; rank(3)=1
     dsetname(4)='dens_dust_tdisk' ; rank(4)=1
     
     dims=tot_ncell

!    Initialize FORTRAN interface.
     CALL h5open_f (error)

     ! Open an existing file.
  !
     CALL h5fopen_f (trim(adjustl(dir_grid))//trim(adjustl(grid_file)), H5F_ACC_RDWR_F, file_id, error)

      do i=0,narr-1
            
     ! Create the dataspace.     
     CALL h5screate_simple_f(rank(i), dims, dspace_id, error)
     
     ! Create the dataset with default properties.
     !     
        CALL h5dcreate_f(file_id, dsetname(i),H5T_NATIVE_DOUBLE , dspace_id, &
          dset_id, error)  
 
        ! Open an existing dataset.
     !
     !   CALL h5dopen_f(file_id, dsetname(i), dset_id, error)  
        ! note the arrays might actually be bigger than dims when created with the standard create_adap_program
  if (dsetname(i) == 'dens_disk' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , dens_disk, dims, error)
  elseif (dsetname(i) == 'dens_tdisk' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , dens_tdisk, dims, error)
  elseif (dsetname(i) == 'dens_bulge' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , dens_bulge, dims, error)
  elseif (dsetname(i) == 'dens_dust_disk' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , dens_dust_disk, dims, error)
  elseif (dsetname(i) == 'dens_dust_tdisk' ) then
     CALL h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE , dens_dust_tdisk, dims, error)

     
  endif
     
! End access to the dataset and release resources used by it.
     !
     CALL h5dclose_f(dset_id, error)

     !
     ! Terminate access to the data space.
     !
     CALL h5sclose_f(dspace_id, error)
     
    end do
      
  ! Terminate access to the file.
     !
     CALL h5fclose_f(file_id, error)     

!    Close FORTRAN interface.
!
     CALL h5close_f(error)

    if (error /= 0.) then 
   print *, 'H5df error =', error 
    STOP 'something did not work with H5DF grid writing'
    else 
       print *, 'extra arrays 3D grid printed'
    endif

  end subroutine print_3d_grid_file_galaxy

subroutine print_profiles_stars
! this is to print emissivity/opacity profiles on a file for comparison with other codes 
real(kind=real64) :: x,y,z,step,cellsize,prof
integer :: i,nstep, id
integer, parameter :: ncomp = 5
real(kind=real64) :: op_disk0,op_tdisk0
character(LEN=lcar_type) :: disk_comp(ncomp)
integer :: disk_type_ID(ncomp)

nstep=1000

disk_comp(1)='disk'    ; disk_type_ID(1)=old_disk_type_ID
disk_comp(2)='tdisk'   ; disk_type_ID(2)=young_disk_type_ID
disk_comp(3)='dust_disk'  ; disk_type_ID(3)=thick_disk_type_ID
disk_comp(4)='dust_tdisk' ; disk_type_ID(4)=thin_disk_type_ID
disk_comp(5)='bulge' !; disk_type_ID(5)='bulge'

do id=1,ncomp

! radial profile
step=rtrun/(nstep-1)

y=0.
z=0
cellsize=5

if (disk_comp(id) /= 'bulge') then 
open(15, file='prof_'//trim(adjustl(disk_comp(id)))//'_R.dat', status='unknown')
else 
open(15, file='prof_'//trim(adjustl(disk_comp(id)))//'_R.dat', status='unknown')
endif 
 
write(15,*)

do i=0, nstep-1

x=i*step

if (disk_comp(id) /= 'bulge') then 
prof=av_disk(x,y,z,cellsize,disk_comp(id),disk_type_ID(id))
else 
prof=av_star_bulge(x,y,z,cellsize)
endif 

write(15,*),x,prof

enddo

close(15)

! z profile 

step=max_z/(nstep-1)

y=0.
x=0.
cellsize=5

if (disk_comp(id) /= 'bulge') then 
open(15, file='prof_'//trim(adjustl(disk_comp(id)))//'_z.dat', status='unknown') 
else 
open(15, file='prof_'//trim(adjustl(disk_comp(id)))//'_z.dat', status='unknown') 
endif 

write(15,*)

do i=0, nstep-1

z=i*step

if (disk_comp(id) /= 'bulge') then 
prof=av_disk(x,y,z,cellsize,disk_comp(id),disk_type_ID(id))
else 
prof=av_star_bulge(x,y,z,cellsize)
endif 

write(15,*),z,prof

enddo

close(15)

end do

print *, 'profiles printed. Continue ?'
read(*,*)

end subroutine print_profiles_stars

!!$subroutine set_galaxy
!!$  ! This routines sets the correct opacity and stellar luminosity
!!$  real(kind=real64) :: lambda_swap, lumstar
!!$  integer :: i
!!$
!!$  if (main_prc) print *, 'Setting opacity parameters...'
!!$  
!!$  allocate(kext_arr(0:lnum-1),kabs_arr(0:lnum-1), ksca_arr(0:lnum-1), gsca_arr(0:lnum -1))
!!$  
!!$  do i =0, lnum-1
!!$     lambda = lambda_arr(i)
!!$     call read_opacity_table
!!$     kext_arr(i) = kext
!!$     kabs_arr(i) = kabs
!!$     ksca_arr(i) = ksca
!!$     gsca_arr(i) = gsca
!!$  enddo
!!$    
!!$  select case (rt_opacity)
!!$
!!$  case('sca')
!!$
!!$     kabs_arr=kabs_arr/kext_arr  !!! normalization : the dens array values in the grid already contain extinction coefficient 
!!$     ksca_arr=ksca_arr/kext_arr
!!$     kext_arr=1.
!!$
!!$  case('abs')
!!$
!!$     kabs_arr=1
!!$     ksca_arr=0
!!$     kext_arr=1
!!$
!!$  case default
!!$     print *, 'rt_opacity not recognized'
!!$     stop     
!!$  end select
!!$
!!$  call print_done
!!$ 
!!$end subroutine set_galaxy

subroutine set_galaxy_input

  call check_grid_type
  lambda = lambda_in
  !call calc_luminosities_pt11
  call read_stellar_sed
  call calc_luminosities
  call calc_scaling_factors_stars    
  call calc_scaling_factors_dust
  

end subroutine set_galaxy_input

!!$subroutine set_galaxy_input_dust
!!$! read dust emission grids that are required for 3D grid creation
!!$  character*6 :: label
!!$  integer :: tot_points,i
!!$  real(KIND=real64), allocatable :: g_ccoord_tot(:,:),g_lum_tot(:)
!!$  character(LEN=lcar_type) :: grid_type_dem_arr(3)
!!$
!!$  grid_type_dem_arr(1) = 'irr1'
!!$  grid_type_dem_arr(2) = 'irr2'
!!$  grid_type_dem_arr(3) = 'irr_HII'
!!$  
!!$  print *, lambda_in/1E4
!!$  !write(label,'(F6.1)') (lambda_in/1E4)   ! um
!!$  write(label,'(I6)') int(lambda_in/1E4)   ! um
!!$  print *, label
!!$
!!$  select case(grid_type_dem)
!!$     ! single component reading 
!!$  case('irr1','irr2','irr_HII')
!!$     file_grid_dem='/pollux/work/gn/data2/gn/MILKY_WAY_MAPS/INPUT_GRIDS/'//'grid_'//trim((adjustl(code_model_dem)))//'_l'//trim(adjustl(label))//'um_'//trim(adjustl(grid_type_dem))//'.dat'
!!$     print *, 'open dust emission input grid (1):'
!!$     print *, file_grid_dem
!!$     call read_grid_pt11(tot_points)
!!$     if (grid_type_dem == 'irr1') then 
!!$        g_lum = g_lum*scaling_dust_comp_arr(1)
!!$     else if (grid_type_dem == 'irr2') then
!!$        g_lum = g_lum*scaling_dust_comp_arr(2)
!!$     else if (grid_type_dem == 'irr_HII') then
!!$        g_lum = g_lum*scaling_dust_comp_arr(3)
!!$     endif
!!$  case('tot')
!!$     ! read the 3 components and sum them 
!!$     do i=1,3
!!$        file_grid_dem='/pollux/work/gn/data2/gn/MILKY_WAY_MAPS/INPUT_GRIDS/'//'grid_'//trim((adjustl(code_model_dem)))//'_l'//trim(adjustl(label))//'um_'//trim(adjustl(grid_type_dem_arr(i)))//'.dat'
!!$        print *, 'open dust emission input grid (3):'
!!$        print *, file_grid_dem
!!$        call read_grid_pt11(tot_points) 
!!$        if (i == 1) then
!!$           allocate(g_ccoord_tot(2,0:tot_points-1), g_lum_tot(0:tot_points-1))
!!$           g_ccoord_tot=g_ccoord
!!$           g_lum_tot=0
!!$           
!!$        endif
!!$       
!!$        if (sum(g_ccoord_tot) /= sum(g_ccoord)) then      ! check same coordinates
!!$           print *, 'maybe coordinates arrays do not match'
!!$           stop
!!$        endif
!!$        g_lum_tot=g_lum_tot+g_lum
!!$           
!!$     end do
!!$     ! set g_lum to total dust emission
!!$     g_lum=g_lum_tot
!!$     
!!$     deallocate(g_ccoord_tot,g_lum_tot)
!!$
!!$  case default
!!$        print *, 'grid_type_dem not recognized'
!!$        stop
!!$
!!$     end select 
!!$
!!$end subroutine set_galaxy_input_dust


subroutine assign_dens_to_parent
! This routine assigns the dens values to subdivided parent cells. This is necessary in the lambda grid creation because in that case only the leaf cells have the dens values already calculated 

integer :: ilvl,i,j,ic,nls

nls=base(2)**3

do ilvl=max_lvl-1,0,-1 

   if (ilvl == 0) nls=base(1)**3

   do i=0,tot_ncell-1 

      if (lvl(i) /= ilvl) cycle 
      
      if (cchild(i) /= -1) then 

         do j=0,nls-1

            ic=cchild(i)+j
            dens(i)=dens(i)+dens(ic)*(csize(ic)/csize(i))**3
            dens_stars(i)=dens_stars(i)+dens_stars(ic)*(csize(ic)/csize(i))**3
            
         end do

      end if

   end do

end do


end subroutine assign_dens_to_parent

!> Reads the input stellar emission SEDs in unit luminosities from the input files file_old_star_sed() and file_young_star_sed().
subroutine read_stellar_sed
logical :: file_exists
integer :: i, id, num 

if (allocated(lambda_arr_old_sed)) return  ! do not do this when creating lambda grids 

! read old stellar SED
inquire(file=trim(adjustl(dir_grid))//trim(adjustl(file_old_star_sed)), exist=file_exists)

if (file_exists) then 

   open(newunit=id,file=trim(adjustl(dir_grid))//trim(adjustl(file_old_star_sed)),status='old')

   call count_lines(id,num)

   nlambda_old_sed = num -1 ! The first line is a comment
   if (nlambda_old_sed <= 0) then
      print *, 'something wrong while reading ', file_old_star_sed
      call stop_prc
   endif

   allocate(lambda_arr_old_sed(0:nlambda_old_sed-1), lnu_old_unit(0:nlambda_old_sed-1))

   read(id,*)

   do i=0,nlambda_old_sed-1 
      
      read(id,*)  lambda_arr_old_sed(i), lnu_old_unit(i)
         
   enddo
   close(id)

else 

   print *, 'STOP: file_old_star_sed not found!'
   print *, 'file = ', trim(adjustl(dir_grid))//trim(adjustl(file_old_star_sed))
   call stop_prc
   
endif

! read young stellar SED
inquire(file=trim(adjustl(dir_grid))//trim(adjustl(file_young_star_sed)), exist=file_exists)

if (file_exists) then 

   open(newunit=id,file=trim(adjustl(dir_grid))//trim(adjustl(file_young_star_sed)),status='old')

   call count_lines(id,num)

   nlambda_sf_sed = num -1 ! The first line is a comment
   if (nlambda_sf_sed <= 0) then
      print *, 'something wrong while reading ', file_young_star_sed
      call stop_prc
   endif

   allocate(lambda_arr_sf_sed(0:nlambda_sf_sed-1), lnu_sf_unit(0:nlambda_sf_sed-1))

   read(id,*)

   do i=0,nlambda_sf_sed-1 
      
      read(id,*)  lambda_arr_sf_sed(i), lnu_sf_unit(i)
         
   enddo
   close(id)

else 

   print *, 'STOP: file_sf_star_sed not found!'
   print *, 'file = ', trim(adjustl(dir_grid))//trim(adjustl(file_young_star_sed))
   call stop_prc
   
endif

lambda_arr_old_sed = lambda_arr_old_sed*1E-6  ! um - > m 
lambda_arr_sf_sed = lambda_arr_sf_sed*1E-6 

end subroutine read_stellar_sed


!!$subroutine read_grid_pt11(tot_points) 
!!$! this reads the grid of emissivities produced by the 2D code 
!!$real(kind=real64) :: r0,r1,z0,z1
!!$real(kind=real64) :: a(3)
!!$integer :: i,j,v,iq(0:0),tot_points
!!$
!!$open(1,file=file_grid_dem,status='old')
!!$
!!$i=0
!!$do 
!!$read(1,*,iostat=v)
!!$if (v /= 0) exit
!!$i=i+1
!!$end do
!!$tot_points=i
!!$rewind(1)
!!$
!!$tot_points=tot_points-2 ! the first two rows are comments  
!!$
!!$!tot_points=4000
!!$if (.not.allocated(g_ccoord)) then 
!!$
!!$   allocate(g_ccoord(2,0:tot_points-1), g_lum(0:tot_points-1))
!!$
!!$endif 
!!$
!!$do i=0,1   ! first two are comments
!!$read(1,*)
!!$enddo 
!!$
!!$
!!$do i=0,tot_points-1 
!!$
!!$read(1,*) (g_ccoord(j,i),j=1,2), g_lum(i)
!!$
!!$enddo
!!$close(1)
!!$
!!$r0=minval(g_ccoord(1,0:tot_points-1))
!!$r1=maxval(g_ccoord(1,0:tot_points-1))
!!$print *, 'R', r0,r1
!!$
!!$z0=minval(g_ccoord(2,0:tot_points-1))
!!$z1=maxval(g_ccoord(2,0:tot_points-1))
!!$print *, 'z', z0,z1
!!$
!!$! count number of R points and z points 
!!$
!!$i=0    ! z points
!!$do 
!!$if (g_ccoord(1,i) /= g_ccoord(1,i+1)) exit  
!!$i=i+1
!!$enddo 
!!$ 
!!$nz_dem=i+1  ! number of z points   IMPORTANT: these formula are valid only if for each R the grid values are present for the same z values. It is also required that the same order is used (each nz lines are for the same R and z is in increasing order). 
!!$nr_dem=tot_points/nz_dem    ! number of R points 
!!$
!!$max_r_dem=r1
!!$max_z_dem=z1
!!$
!!$
!!$
!!$end subroutine read_grid_pt11

!!$real(kind=real64) function av_dustem(x,y,z,cellsize)
!!$  ! this calculates the average emissivity in the current cell
!!$real(Kind=real64) :: x,y,z,r0,xmin,xmax,ymin,ymax,zmin,zmax,rmin,rmax,cellsize,step
!!$real(kind=real64) :: tot_lum_dustem,xp,yp,zp,rp,step_int,fq(4),dr_dz(5)
!!$integer :: i,flagrmax,flagzmax,flagrmin,flagzmin,irmin,irmax,izmin,izmax,ir,iz,k,ix,iy,iz_int
!!$integer :: i11,i12,i21,i22
!!$real(kind=real64), parameter :: step_pc=20 ! integration step 
!!$!integer, parameter :: np=100
!!$integer :: np
!!$
!!$av_dustem=0
!!$ 
!!$xmin=x-cellsize/2.
!!$xmax=x+cellsize/2.
!!$ymin=y-cellsize/2.
!!$ymax=y+cellsize/2.
!!$zmin=z-cellsize/2.
!!$zmax=z+cellsize/2.
!!$
!!$!!! note: it would be better to make np variable depending on the size of the cell 
!!$!!!! better : the integration should be done on a logarithmic scale when the cell
!!$! is very large and on a linear scale when is small... but this is harder to implement in more dimensions   
!!$
!!$np=int(cellsize/step_pc)
!!$ 
!!$if (np < 2) np=2  ! min numbers of integration steps  
!!$if (np >200 ) np=200 ! max number of integration steps 
!!$
!!$step_int=cellsize/(np-1)   ! integration step
!!$
!!$!print *, step_int
!!$tot_lum_dustem=0.
!!$
!!$! create grid of points inside a cube for integration
!!$do ix=0,np-2  ! -2 because middle points are considered for integration
!!$   do iy=0, np-2 
!!$      do iz_int=0, np-2 
!!$ 
!!$         xp=xmin+(ix+0.5)*step_int    !!! these are the 3D coordinates of the points 
!!$         yp=ymin+(iy+0.5)*step_int    !!! inside the cube to calculate the emissivity for 
!!$         zp=zmin+(iz_int+0.5)*step_int
!!$         
!!$         rp=sqrt(xp**2+yp**2)
!!$                       
!!$         do ir=nr_dem-1,1,-1 ! until 1 and not zero because otherwise it can exit as -1      
!!$            if (rp > g_ccoord(1,ir*nz_dem)) exit              
!!$       enddo
!!$
!!$       if (rp > max_r_dem) cycle  ! zero emissivity for points beyond max R 
!!$
!!$       do iz=nz_dem-1,1,-1 
!!$            if (abs(zp) > g_ccoord(2,iz)) exit
!!$       enddo
!!$
!!$       if (abs(zp) > max_z_dem) cycle  ! zero emissivity for points beyond max z 
!!$ !!$print *, rp,zp
!!$ !!$print *, g_ccoord(1,ir*nz), g_ccoord(2,iz)
!!$
!!$       i11= ir*nz_dem+iz    ! indeces of points around rp,zp
!!$       i12= ir*nz_dem+iz+1
!!$       i21= (ir+1)*nz_dem+iz
!!$       i22= (ir+1)*nz_dem+iz+1
!!$ !!$       print *,' rp zp'
!!$ !!$       print *, rp, zp
!!$ !!$       print *, 'indeces'
!!$ !!$       print *, ir,iz,nr-1,nz-1
!!$       
!!$
!!$       fq(1)=g_lum(i11)   !! f11
!!$       fq(2)=g_lum(i12)   !! f12 
!!$       fq(3)=g_lum(i21)   !! f21 
!!$       fq(4)=g_lum(i22)   !! f22
!!$
!!$     dr_dz(1)=(g_ccoord(1,(ir+1)*nz_dem)-rp)*(g_ccoord(2,iz+1)-abs(zp))  !!! bilinear interpolation
!!$     dr_dz(2)=(g_ccoord(1,(ir+1)*nz_dem)-rp)*(abs(zp)-g_ccoord(2,iz)) 
!!$     dr_dz(3)=(rp-g_ccoord(1,ir*nz_dem))*(g_ccoord(2,iz+1)-abs(zp)) 
!!$     dr_dz(4)=(rp-g_ccoord(1,ir*nz_dem))*(abs(zp)-g_ccoord(2,iz)) 
!!$     dr_dz(5)=(g_ccoord(1,(ir+1)*nz_dem)-g_ccoord(1,ir*nz_dem))*(g_ccoord(2,iz+1)-g_ccoord(2,iz))
!!$       
!!$ !!$       print *, 'ccord'
!!$ !!$       print *, g_ccoord(1,i11),g_ccoord(2,i11)
!!$ !!$       print *, g_ccoord(1,i12),g_ccoord(2,i12)
!!$ !!$       print *, g_ccoord(1,i21),g_ccoord(2,i21)
!!$ !!$       print *, g_ccoord(1,i22),g_ccoord(2,i22)
!!$       
!!$      tot_lum_dustem=tot_lum_dustem+(sum(fq*dr_dz(1:4))/dr_dz(5))  !!! later you multiply by the step volume to get the total energy in the cell
!!$ !!$      print *, tot_lum_dustem
!!$ !!$      read(*,*)
!!$
!!$
!!$enddo
!!$enddo
!!$enddo
!!$
!!$
!!$tot_lum_dustem=tot_lum_dustem*step_int**3
!!$av_dustem=tot_lum_dustem/cellsize**3
!!$  
!!$end function av_dustem

!> Assigns the right value of the disk_type_ID() depending on the input disk_type(). Using numbers instead of character variables is more efficient in IF statement and less error prone while modifying the code. 
subroutine assign_disk_type_ID(disk_type, disk_type_ID)
character(LEN=lcar_type) :: disk_type
integer :: disk_type_ID

select case(disk_type)
case('expR_expz')
   disk_type_ID = expR_expz_ID
case('expR_sech2z')
   disk_type_ID = expR_sech2z_ID   
case('flared_expz')
   disk_type_ID = flared_expz_ID
case('flared_sech2z')
   disk_type_ID = flared_sech2z_ID
case('ellipt_expR_expz')
   disk_type_ID = ellipt_expR_expz_ID
case('ellipt_expR_sech2z')
   disk_type_ID = ellipt_expR_sech2z_ID
case default 
   STOP '(assign_disk_type_ID): input disk_type not recognized!'
end select

end subroutine assign_disk_type_ID

END MODULE user_routines_galaxy
