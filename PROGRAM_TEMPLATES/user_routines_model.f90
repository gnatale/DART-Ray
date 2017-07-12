MODULE user_routines_model
use smooth_grid_routines
use iso_fortran_env
use HDF5
use io_routines
use sed_routines
IMPLICIT NONE  
 
! Galaxy parameters 
!> @param lambda_min Minimum wavelength used in the calculation of the lambda grids.  
real(kind=real64) :: lambda_min
!> @param lambda_max Maximum wavelength used in the calculation of the lambda grids. 
real(kind=real64) :: lambda_max

! total luminosity
!> @param lnu_tot Model total luminosity. 
real(kind=real64) :: lnu_tot

! Filenames and keywords 
!> @param lcar_type String variable length parameter. 
integer, parameter :: lcar_type=30
! ! @param lcar2 String variable length parameter. 
!integer, parameter :: lcar2=80

! integration step numbers
integer, parameter :: step_int=10

!> @param z_subd_lim Z coordinate below which cells have to be subdivided if their R is less than R_subd_lim()
real(kind=real64) :: z_subd_lim 

!> @param R_subd_lim R coordinate below which cells have to be subdivided if their z is less than z_subd_lim()
real(kind=real64) :: R_subd_lim  

!> @param subdivision_criteria Subdivision criteria to be used. Choices are:
!> - 'standard': subdivision IF subdivision level < max_lvl() AND (cell optical depth > max_dtau() OR stellar luminosity > max_dlum()*total luminosity OR  subdivision level < min_lvl() OR (ABS(z) < z_subd_lim AND R < R_subd_lim())) 
character(LEN=lcar_type) :: subdivision_criteria

XXX add model variables here XXX 

! input namelist

XXX add model input parameters in the namelist blocks (you can create new blocks if you want) XXXX 

namelist /model_input_strings/ label_model_lambda_grid, dir_grid, grid_file, grid_info_file, file_lambda_list,units_lambda, dust_model,dust_opacity_tables, file_gra_fa, file_sil_fa, file_pah_neu_fa, file_pah_ion_fa, file_av_opacities, file_q_gra, file_q_sil, file_q_pah_neu, file_q_pah_ion, subdivision_criteria

namelist /model_input_var/ lambda_ref, lambda_min, lambda_max, modelsize, base, max_ncell, max_lvl, min_lvl,  max_dtau, max_dlum, n_dust_size_qabs, n_dust_wave_qabs, z_subd_lim, R_subd_lim

namelist /model_input_logical/ input_av_opacities

CONTAINS

!> Reads input file for the grid creation for the user-defined model. 
subroutine read_input_model
integer :: ioun, i, id  
character(len=lcar) :: input_filename
integer :: n_hs

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

XXX modify namelist block names if needed XXX 

read(ioun, nml = model_input_strings )
read(ioun, nml = model_input_var )
read(ioun, nml = model_input_logical)
close(ioun)

! check dust model
call check_input_dust_model

! set grid creation true 
grid_creation = .TRUE.

call print_done

end subroutine read_input_model

!> This subroutine is used to derive all the parameters that are necessary to obtain the dust density and the stellar luminosity density at the wavelength of the main grid or those of the lambda grids 
subroutine set_model_input

 XXX add calls to subroutines deriving all the necessary parameters used in the calculation of the average dust density and stellar luminosity density at a single wavelength lambda XXX 
  

end subroutine set_model_input



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


END MODULE user_routines_model
