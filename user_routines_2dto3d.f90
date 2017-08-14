!> Contains the routines to interpolate a 2D (R,z) grid into a 3D grid. 
MODULE user_routines_2dto3d
use smooth_grid_routines
use iso_fortran_env
use HDF5
use io_routines
use sed_routines
IMPLICIT NONE  
 
!> @param lambda_in Wavelength variable 
real(kind=real64) :: lambda_in

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

! integration step numbers
integer, parameter :: step_int=10

!> @param z_subd_lim Z coordinate below which cells have to be subdivided if their R is less than R_subd_lim()
real(kind=real64) :: z_subd_lim 

!> @param R_subd_lim R coordinate below which cells have to be subdivided if their z is less than z_subd_lim()
real(kind=real64) :: R_subd_lim  

!> @param subdivision_criteria Subdivision criteria to be used. Choices are:
!> - 'standard': subdivision IF subdivision level < max_lvl() AND (cell optical depth > max_dtau() OR stellar luminosity > max_dlum()*total luminosity OR  subdivision level < min_lvl() OR (ABS(z) < z_subd_lim AND R < R_subd_lim())) 
character(LEN=lcar_type) :: subdivision_criteria

!> @param dir_grid_2d Directory where the 2D grids are saved. 
character(LEN=lcar) :: dir_grid_2d

!> @param label_model_2d Label contained in the 2D grid files: 'grid_'+label_model_2d+'_l'+wavelength+'um.dat'
character(LEN=lcar) :: label_model_2d

!> @param g_coord 2D grid coordinates. The first time a 2D grid is read, the program simply assigns values to this array. The following times, the program checks that the values are exactly the same 
real(KIND=real64), allocatable :: g_ccoord(:,:)

!> @param g_lum Luminosity density values. Remember that they have to be either in [W/Hz/pc^3] or [erg/Hz/pc^3].
real(KIND=real64), allocatable :: g_lum(:)

!> @param g_dens Extinction coefficient values. Remember that it has to be always in [pc^-1]
real(KIND=real64), allocatable :: g_dens(:)

!> @param tot_points_2d Total number of 2D grid points. 
integer :: tot_points_2d

!> @param nr_dem Number of sampling points in the R dimension. 
integer :: nr_dem 

!> @param nz_dem Number of sampling points in the z dimension. 
integer :: nz_dem 

!> @param max_r_dem Max value of radial distance R. 
real(kind=real64) :: max_r_dem

!> @param max_z_dem Max value of vertical distance z. 
real(kind=real64) :: max_z_dem

!> @param quantity_em ID to identify volume emissivity calculations.
integer, parameter :: quantity_em = 0 

!> @param quantity_ext ID to identify dust extinction coefficient calculations. 
integer, parameter :: quantity_ext = 1 

! input namelist

!XXX add model input parameters in the namelist blocks (you can create new blocks if you want) XXXX 

namelist /m2dto3d_input_strings/ label_model_lambda_grid, dir_grid, grid_file, grid_info_file, file_lambda_list,units_lambda, dust_model,dust_opacity_tables, file_gra_fa, file_sil_fa, file_pah_neu_fa, file_pah_ion_fa, file_av_opacities, file_q_gra, file_q_sil, file_q_pah_neu, file_q_pah_ion, subdivision_criteria, dir_grid_2d, label_model_2d

namelist /m2dto3d_input_var/ lambda_ref, lambda_min, lambda_max, modelsize, base, max_ncell, max_lvl, min_lvl,  max_dtau, max_dlum, n_dust_size_qabs, n_dust_wave_qabs, z_subd_lim, R_subd_lim

namelist /m2dto3d_input_logical/ input_av_opacities

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
read(ioun, nml = m2dto3d_input_strings )
read(ioun, nml = m2dto3d_input_var )
read(ioun, nml = m2dto3d_input_logical)
close(ioun)

! check dust model
call check_input_dust_model

! check input variables 
call check_input_2dto3d

! set grid creation true 
grid_creation = .TRUE.

call print_done

end subroutine read_input_model

!> Checks input variable values 
subroutine check_input_2dto3d 

! check dir_grid_2d 
call check_dir_existence(dir_grid_2d, .FALSE.) 


end subroutine check_input_2dto3d



!> This subroutine is used to derive all the parameters that are necessary to obtain the dust density and the stellar luminosity density at the wavelength of the main grid or those of the lambda grids 
subroutine set_2dto3d_input
  integer :: ilambda

  !print *, 'setting 2dto3d input...'
  
  call find_ilambda(ilambda)
  
  ! read 2D grid at the specified wavelength   
  call read_grid_2d(ilambda)
 
end subroutine set_2dto3d_input

!> Find wavelength index corresponding to the input wavelength. If the wavelength is not found within the wavelength list, the program exits with an error.  
subroutine find_ilambda(ilambda) 
integer :: i, ilambda

  do i = 0, lnum

     if (abs(lambda_in-lambda_arr_SI(i))/lambda_arr_SI(i) < 1E-4) exit 

  end do

  if (i > lnum -1) then 
     print *, 'ERROR(find_ilambda): lambda_in not found in lambda_arr_SI array!'
     stop
  endif

  ilambda = i 


end subroutine find_ilambda


!> Reads the 2D grids of emissivities or dust extinction coefficients  
subroutine read_grid_2d(ilambda)

real(kind=real64) :: r0,r1,z0,z1
real(kind=real64) :: a(3)
integer :: i,j,v,iq(0:0)
integer :: ilambda
character(len=lcar) :: file_grid_2d
real(kind=real64), allocatable :: g_ccoord_temp(:,:), g_lum_temp(:), g_dens_temp(:)
logical :: first_read 
integer :: id 
character(len=lcar) :: grid_quantity

! set first_read variable 
first_read = .FALSE.

! set file name 
file_grid_2d = 'grid_'//trim(adjustl(label_model_2d))//'_l'//trim(adjustl(label_wave_arr(ilambda)))//'um.dat'

! check file existence 
call check_file_existence(dir_grid_2d, file_grid_2d)

! open file 
print *, 'reading file '//trim(adjustl(file_grid_2d))

open(newunit=id,file=trim(adjustl(dir_grid_2d))//trim(adjustl(file_grid_2d)),status='old')

! count lines 
call count_lines(id,tot_points_2d)

tot_points_2d=tot_points_2d-1 ! the first row is a comment  

! allocate temporary grid arrays 
allocate(g_ccoord_temp(2,0:tot_points_2d-1), g_lum_temp(0:tot_points_2d-1), g_dens_temp(0:tot_points_2d-1))

! allocate permanent arrays if needed 
if (.not.allocated(g_ccoord)) then 

   first_read = .TRUE.
   allocate(g_ccoord(2,0:tot_points_2d-1), g_lum(0:tot_points_2d-1), g_dens(0:tot_points_2d-1))

endif 

! read file 
read(id,*) ! first line is a comment 

do i=0,tot_points_2d-1 
   
   read(id,*) (g_ccoord_temp(j,i),j=1,2), g_lum_temp(i), g_dens_temp(i)
   
enddo
close(id)

! Assign to permanent grid point array and write values range

if (first_read) then 

   g_ccoord = g_ccoord_temp

   ! Figure out ranges of R and z 
   r0=minval(g_ccoord(1,0:tot_points_2d-1))
   r1=maxval(g_ccoord(1,0:tot_points_2d-1))
   print *, 'R range: ', r0,r1

   z0=minval(g_ccoord(2,0:tot_points_2d-1))
   z1=maxval(g_ccoord(2,0:tot_points_2d-1))
   print *, 'z range: ', z0,z1

   max_r_dem=r1
   max_z_dem=z1

   ! count number of R points and z points
   i=0    ! z points
   do 
      if (g_ccoord(1,i) /= g_ccoord(1,i+1)) exit  
      i=i+1
   enddo

   nz_dem=i+1  ! number of z points   
   nr_dem=tot_points_2d/nz_dem    ! number of R points 

   print *, 'N(R) = ', nr_dem
   print *, 'N(z) =', nz_dem 

   ! check ascending order for R values 
   do i = 1, tot_points_2d-1
      
      if (g_ccoord(1,i) < g_ccoord(1,i-1)) then 
         print *, 'ERROR(read_grid_2d): R values are not all in ascending order!'
         print *, 'check row number ', i+1
         stop
      endif

   enddo

   ! check repetition for R values 
   do i = 0, nr_dem -1 
       
      if (count(g_ccoord(1,i*nz_dem: (i+1)*nz_dem-1) /= g_ccoord(1,i*nz_dem)) >0) then 
         print *, 'ERROR(read_grid_2d): R values are not repeating as expected!'
         print *, 'check rows from ', i*nz_dem+2, 'to ', (i+1)*nz_dem+1
         stop
      endif
   
   enddo

   ! check ascending order for z values 
   do i = 1, nz_dem -1 

      if (g_ccoord(2,i) < g_ccoord(2,i-1)) then
         print *, 'ERROR(read_grid_2d): z values are not all in ascending order!'
         print *,'check row number ', i+1
         stop
      endif

   end do

   ! check repetition for z values 
   do i = 0, nz_dem -1 

      do j = 0, nr_dem -1 
    
         if (g_ccoord(2,i) /= g_ccoord(2,i+j*nz_dem)) then 
            print *, 'ERROR(read_grid_2d): z values are not repeating as expected' 
            print *, 'check row number ', i+j*nz_dem+2
            stop
         endif

      end do

   end do
   
else   ! or check grid points are the same as 

    if (count(g_ccoord /= g_ccoord_temp) > 0) then 
       print *, 'ERROR(read_grid_2d): The 2D grid points just read are different from those in the previous tables! The values have to be exactly the same in all input 2D grids!'
       stop
    endif

endif

! assign volume emissivity and dust extinction coefficients  
g_lum = g_lum_temp
g_dens = g_dens_temp


 
!!$nz_dem=i+1  ! number of z points   IMPORTANT: these formula are valid only if for each R the grid values are present for the same z values. It is also required that the same order is used (each nz lines are for the same R and z is in increasing order). 
!!$nr_dem=tot_points/nz_dem    ! number of R points 
!!$
!!$max_r_dem=r1
!!$max_z_dem=z1

! deallocate temporary arrays 
deallocate(g_ccoord_temp, g_lum_temp, g_dens_temp)

end subroutine read_grid_2d

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

!> Calculates the average emissivity ([W/Hz/pc^3] or [erg/Hz/pc^3]) or dust extinction coefficient ([pc^-1]) in the current cell.
real(kind=real64) function av_dens_2dto3d(x,y,z,cellsize, quantity_ID)
  
real(Kind=real64) :: x,y,z,r0,xmin,xmax,ymin,ymax,zmin,zmax,rmin,rmax,cellsize,step
real(kind=real64) :: tot_dustem,xp,yp,zp,rp,step_int,fq(4),dr_dz(5)
integer :: i,flagrmax,flagzmax,flagrmin,flagzmin,irmin,irmax,izmin,izmax,ir,iz,k,ix,iy,iz_int
integer :: i11,i12,i21,i22
real(kind=real64) :: step_pc ! integration step 
integer :: np 
integer :: quantity_ID


av_dens_2dto3d=0

select case(quantity_ID)
case(quantity_em, quantity_ext)
 
case default 
   print *, 'ERROR(av_dens_2dto3d): variable quantity_ID not recognized!'
   print *, 'quantity_ID =', quantity_ID
   stop
end select

xmin=x-cellsize/2.
xmax=x+cellsize/2.
ymin=y-cellsize/2.
ymax=y+cellsize/2.
zmin=z-cellsize/2.
zmax=z+cellsize/2.

step_pc = modelsize/1000.

np=int(cellsize/step_pc)
 
if (np < 2) np=2  ! min numbers of integration steps  
if (np >200 ) np=200 ! max number of integration steps 

step_int=cellsize/(np-1)   ! integration step

tot_dustem=0.

! create grid of points inside a cube for integration
do ix=0,np-2  ! -2 because middle points are considered for integration
   do iy=0, np-2 
      do iz_int=0, np-2 
 
         xp=xmin+(ix+0.5)*step_int    !!! these are the 3D coordinates of the points 
         yp=ymin+(iy+0.5)*step_int    !!! inside the cube to calculate the emissivity for 
         zp=zmin+(iz_int+0.5)*step_int
         
         rp=sqrt(xp**2+yp**2)
                       
         do ir=nr_dem-1,1,-1 ! until 1 and not zero because otherwise it can exit as -1      
            if (rp > g_ccoord(1,ir*nz_dem)) exit              
       enddo

       if (rp > max_r_dem) cycle  ! zero emissivity for points beyond max R 

       do iz=nz_dem-1,1,-1 
            if (abs(zp) > g_ccoord(2,iz)) exit
       enddo

       if (abs(zp) > max_z_dem) cycle  ! zero emissivity for points beyond max z 
 !!$print *, rp,zp
 !!$print *, g_ccoord(1,ir*nz), g_ccoord(2,iz)

       i11= ir*nz_dem+iz    ! indeces of points around rp,zp
       i12= ir*nz_dem+iz+1
       i21= (ir+1)*nz_dem+iz
       i22= (ir+1)*nz_dem+iz+1
 !!$       print *,' rp zp'
 !!$       print *, rp, zp
 !!$       print *, 'indeces'
 !!$       print *, ir,iz,nr-1,nz-1
       
       select case (quantity_ID)
       case(quantity_em)
          fq(1)=g_lum(i11)   !! f11
          fq(2)=g_lum(i12)   !! f12 
          fq(3)=g_lum(i21)   !! f21 
          fq(4)=g_lum(i22)   !! f22
       case(quantity_ext)
          fq(1)=g_dens(i11)   !! f11
          fq(2)=g_dens(i12)   !! f12 
          fq(3)=g_dens(i21)   !! f21 
          fq(4)=g_dens(i22)   !! f22
       case default
          print *, 'ERROR(av_dens_2dto3d): here you should never get!'
          stop
       end select

     dr_dz(1)=(g_ccoord(1,(ir+1)*nz_dem)-rp)*(g_ccoord(2,iz+1)-abs(zp))  !!! bilinear interpolation
     dr_dz(2)=(g_ccoord(1,(ir+1)*nz_dem)-rp)*(abs(zp)-g_ccoord(2,iz)) 
     dr_dz(3)=(rp-g_ccoord(1,ir*nz_dem))*(g_ccoord(2,iz+1)-abs(zp)) 
     dr_dz(4)=(rp-g_ccoord(1,ir*nz_dem))*(abs(zp)-g_ccoord(2,iz)) 
     dr_dz(5)=(g_ccoord(1,(ir+1)*nz_dem)-g_ccoord(1,ir*nz_dem))*(g_ccoord(2,iz+1)-g_ccoord(2,iz))
       
 !!$       print *, 'ccord'
 !!$       print *, g_ccoord(1,i11),g_ccoord(2,i11)
 !!$       print *, g_ccoord(1,i12),g_ccoord(2,i12)
 !!$       print *, g_ccoord(1,i21),g_ccoord(2,i21)
 !!$       print *, g_ccoord(1,i22),g_ccoord(2,i22)
       
      tot_dustem=tot_dustem+(sum(fq*dr_dz(1:4))/dr_dz(5))  !!! later you multiply by the step volume to get the total value for the cell
 !!$      print *, tot_dustem
 !!$      read(*,*)


enddo
enddo
enddo


tot_dustem=tot_dustem*step_int**3
av_dens_2dto3d=tot_dustem/cellsize**3
  
end function av_dens_2dto3d


END MODULE user_routines_2dto3d
