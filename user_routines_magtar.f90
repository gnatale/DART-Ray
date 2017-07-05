MODULE user_routines_magtar
use smooth_grid_routines
use iso_fortran_env
use io_routines
use sed_routines, only : bplanck
IMPLICIT NONE 

! Standard parameters
!> @param tau_z Optical depth per pc at the reference wavelength. 
real(kind=real64) :: tau_z
!> @param Dust number density / Extinction coefficient 
real(kind=real64) :: ndust

! Ellipsoid geometry parameters
!> @param ax Ellipsoidal semi-axis along the X-axis. 
real(kind=real32) :: ax=0
!> @param by Ellipsoidal semi-axis along the Y-axis.
real(kind=real32) :: by=0
!> @param cz Ellipsoidal semi-axis along the Z-axis.
real(kind=real32) :: cz=0 ! pc ellipsoid semiaxis
!> @param elrad_width Semi-width of the dust shell in normalized radial units. 
real(kind=real64) :: elrad_width=0 
! Ring geometry parameters 
!real(kind=real32) :: rad_ring=0, w_ring=0,h_ring=0
! Geometry label 
!> @param dust_geometry  Type of dust geometries. Choices: 'shell', 'cavity', 'wind'. See corresponding av_rho_dust subroutines for the assumed formula. 
character*8 :: dust_geometry 

!> @param R_subd_lim normalized R coordinate below which cells have to be subdivided if their optical depth is higher than max_dtau().  
real(kind=real64) :: R_subd_lim  

!> @param subdivision_criteria Subdivision criteria to be used. Choices are:
!> 'shell': subdivide IF subdivision level < max_lvl() AND (cell optical depth > max_dtau OR subdivision level < min_lvl())
!> 'cavity': subdivide IF subdivision level < max_lvl() AND (cell optical depth > max_dtau AND R < R_subd_lim)
!> 'wind' : same as for 'cavity'
character(LEN=lcar) :: subdivision_criteria

! input namelist 
namelist /magtar_input_strings/ dir_grid, grid_file, grid_info_file,  file_lambda_list, units_lambda,dust_model,dust_opacity_tables, file_gra_fa, file_sil_fa, file_pah_neu_fa, file_pah_ion_fa, file_av_opacities, file_q_gra, file_q_sil, file_q_pah_neu, file_q_pah_ion, dust_geometry, subdivision_criteria

namelist /magtar_input_var/ modelsize, tau_z, lambda_ref, base, ax, by, cz, elrad_width, max_ncell, max_lvl, min_lvl, max_dtau, max_dlum, n_dust_size_qabs, n_dust_wave_qabs, R_subd_lim  

namelist /magtar_input_logical/ input_av_opacities


CONTAINS

!> Reads input file for the grid creation for the magtar models.
!> \todo check if input_av_opacities is really needed here. 
subroutine read_input_magtar

integer :: ioun
character(len=lcar) :: input_filename

call get_command_argument(1,input_filename)

!check input_file argument has been given
if (trim(input_filename) == '') then 
   if (main_prc) print *, 'provide input file in command line...EXIT' 
   call stop_prc
endif

! read input file 
if (main_prc) print *, 'reading input file...'

open(newunit=ioun, file = input_filename, status ='old', action ='read')
read(ioun, nml = magtar_input_strings)
read(ioun, nml = magtar_input_var)
read(ioun, nml = magtar_input_logical)
close(ioun)

! check dust model
call check_input_dust_model

! set grid creation true 
grid_creation = .TRUE.

call print_done

end subroutine read_input_magtar

!> Sets dust density for the specified dust_geometry(). 
subroutine set_mag_model_density()  

  kext_ref=kext_ref/(parsec**2)   ! m^2 --> pc ^2 

!!$  if (dust_geometry == 'ring') then 
!!$     ndust=tau_z/(kext_ref*h_ring)   ! pc ^-3  !!! optical depth perpedicular to the ring (along z direction)
  if (dust_geometry == 'shell') then
     ! ndust = tau_z/(kext*2*cz*elrad_width)  !!! optical depth perpendicular to the ellipsoidal shell
     ndust = tau_z/(kext_ref)
  else if (dust_geometry == 'cavity') then  !!! optical depth per parsec 
     ndust = tau_z/(kext_ref)
  else if (dust_geometry == 'wind') then
     ndust = tau_z/(kext_ref)
  endif

  ndust=ndust*kext_ref    ! pc ^-1  this is now the extinction coefficient. This passage seems useless now but it is important when scaling at wavelengths different than lambda = 1 um
  

!  call read_planck_table_magtar(lambda,lumstar)   ! do not confuse lumstar with lum_stars !! lumstar is the luminosity of the single star in the magtar benchmark


end subroutine set_mag_model_density

!!$real(kind=real64) function av_rho_dust_ring(x,y,z,cellsize)
!!$INTEGER, PARAMETER :: nstep=10
!!$real(Kind=real64) :: x,y,z,cellsize,xmin,xmax,ymin,ymax,zmin,zmax,rad
!!$real(Kind=real64) :: xp(0:nstep-1),yp(0:nstep-1),zp(0:nstep-1),sx,sy,sz,tot_ndust,rho_dust,cvol
!!$integer :: i,j,k
!!$logical :: dust_x_lim, dust_y_lim, dust_z_lim,dust_rad_lim
!!$
!!$xmin=x-cellsize/2.
!!$xmax=x+cellsize/2.
!!$ymin=y-cellsize/2.
!!$ymax=y+cellsize/2.
!!$zmin=z-cellsize/2.
!!$zmax=z+cellsize/2.
!!$
!!$sx=abs(xmax-xmin)/(nstep-1)
!!$sy=abs(ymax-ymin)/(nstep-1)
!!$sz=abs(zmax-zmin)/(nstep-1)
!!$
!!$do i=0,nstep-2
!!$xp(i)=xmin+(i+0.5)*sx
!!$yp(i)=ymin+(i+0.5)*sy
!!$zp(i)=zmin+(i+0.5)*sz
!!$end do
!!$
!!$tot_ndust=0.
!!$cvol=sx*sy*sz
!!$
!!$
!!$do i=0,nstep-2
!!$ do j=0,nstep-2
!!$  do k=0, nstep-2
!!$
!!$   rad=sqrt(xp(i)**2+yp(j)**2)
!!$   dust_rad_lim=(abs(rad_ring-rad) < w_ring/2.) 
!!$   dust_z_lim=(abs(zp(k)) < h_ring/2.)
!!$     
!!$if ((dust_rad_lim).and.(dust_z_lim)) then 
!!$rho_dust=ndust 
!!$else
!!$rho_dust=0
!!$endif
!!$
!!$tot_ndust=tot_ndust+rho_dust
!!$
!!$end do
!!$end do
!!$end do
!!$
!!$av_rho_dust_ring=tot_ndust*cvol/cellsize**3
!!$
!!$
!!$end function av_rho_dust_ring

real(kind=real64) function av_rho_dust_shell(x,y,z,cellsize)

INTEGER, PARAMETER :: nstep=10
real(Kind=real64) :: x,y,z,cellsize,xmin,xmax,ymin,ymax,zmin,zmax,rad
real(Kind=real64) :: a,b,c,width
real(Kind=real64) :: xp(0:nstep-1),yp(0:nstep-1),zp(0:nstep-1),sx,sy,sz,tot_ndust,rho_dust,cvol
integer :: i,j,k
logical :: dust_x_lim, dust_y_lim, dust_z_lim,dust_rad_lim

xmin=x-cellsize/2.
xmax=x+cellsize/2.
ymin=y-cellsize/2.
ymax=y+cellsize/2.
zmin=z-cellsize/2.
zmax=z+cellsize/2.

sx=abs(xmax-xmin)/(nstep-1)
sy=abs(ymax-ymin)/(nstep-1)
sz=abs(zmax-zmin)/(nstep-1)

do i=0,nstep-2
xp(i)=xmin+(i+0.5)*sx
yp(i)=ymin+(i+0.5)*sy
zp(i)=zmin+(i+0.5)*sz
end do

tot_ndust=0.
cvol=sx*sy*sz

do i=0,nstep-2
 do j=0,nstep-2
  do k=0, nstep-2

   rad=sqrt((xp(i)/ax)**2+(yp(j)/by)**2+(zp(k)/cz)**2)
   dust_rad_lim=(abs(rad-1) < elrad_width) 
        
if (dust_rad_lim) then 
rho_dust=ndust 
else
rho_dust=0
endif

tot_ndust=tot_ndust+rho_dust

end do
end do
end do

av_rho_dust_shell=tot_ndust*cvol/cellsize**3


end function av_rho_dust_shell

real(kind=real64) function av_rho_dust_cavity(x,y,z,cellsize)

INTEGER, PARAMETER :: nstep=10
real(Kind=real64) :: x,y,z,cellsize,xmin,xmax,ymin,ymax,zmin,zmax,rad
real(Kind=real64) :: a,b,c,width
real(Kind=real64) :: xp(0:nstep-1),yp(0:nstep-1),zp(0:nstep-1),sx,sy,sz,tot_ndust,rho_dust,cvol
integer :: i,j,k
logical :: dust_x_lim, dust_y_lim, dust_z_lim,dust_rad_lim

xmin=x-cellsize/2.
xmax=x+cellsize/2.
ymin=y-cellsize/2.
ymax=y+cellsize/2.
zmin=z-cellsize/2.
zmax=z+cellsize/2.

sx=abs(xmax-xmin)/(nstep-1)
sy=abs(ymax-ymin)/(nstep-1)
sz=abs(zmax-zmin)/(nstep-1)

do i=0,nstep-2
xp(i)=xmin+(i+0.5)*sx
yp(i)=ymin+(i+0.5)*sy
zp(i)=zmin+(i+0.5)*sz
end do

tot_ndust=0.
cvol=sx*sy*sz

do i=0,nstep-2
 do j=0,nstep-2
  do k=0, nstep-2

   rad=sqrt((xp(i)/ax)**2+(yp(j)/by)**2+(zp(k)/cz)**2)
   dust_rad_lim=(rad > 1) 
        
if (dust_rad_lim) then 
rho_dust=ndust 
else
rho_dust=0
endif

tot_ndust=tot_ndust+rho_dust

end do
end do
end do

av_rho_dust_cavity=tot_ndust*cvol/cellsize**3


end function av_rho_dust_cavity

real(kind=real64) function av_rho_dust_wind(x,y,z,cellsize)

INTEGER, PARAMETER :: nstep=10
real(Kind=real64) :: x,y,z,cellsize,xmin,xmax,ymin,ymax,zmin,zmax,rad
real(Kind=real64) :: a,b,c,width
real(Kind=real64) :: xp(0:nstep-1),yp(0:nstep-1),zp(0:nstep-1),sx,sy,sz,tot_ndust,rho_dust,cvol
integer :: i,j,k
logical :: dust_x_lim, dust_y_lim, dust_z_lim,dust_rad_lim

xmin=x-cellsize/2.
xmax=x+cellsize/2.
ymin=y-cellsize/2.
ymax=y+cellsize/2.
zmin=z-cellsize/2.
zmax=z+cellsize/2.

sx=abs(xmax-xmin)/(nstep-1)
sy=abs(ymax-ymin)/(nstep-1)
sz=abs(zmax-zmin)/(nstep-1)

do i=0,nstep-2
xp(i)=xmin+(i+0.5)*sx
yp(i)=ymin+(i+0.5)*sy
zp(i)=zmin+(i+0.5)*sz
end do

tot_ndust=0.
cvol=sx*sy*sz

do i=0,nstep-2
 do j=0,nstep-2
  do k=0, nstep-2

   rad=sqrt((xp(i)/ax)**2+(yp(j)/by)**2+(zp(k)/cz)**2)
   !dust_rad_lim=(rad <= 1+elrad_width) 
        
if (rad <= 1 ) then 
   rho_dust=ndust*(rad**2)
   
elseif (rad > 1) then
   rho_dust=ndust*(rad**(-2))
endif

tot_ndust=tot_ndust+rho_dust

end do
end do
end do

av_rho_dust_wind=tot_ndust*cvol/cellsize**3


end function av_rho_dust_wind


!!$subroutine read_lambda_list
!!$integer :: i, v
!!$
!!$open(12,file=file_lambda_list, status='old')
!!$
!!$i=0 ! there are no comments in lambda list file
!!$do 
!!$   read(12,*,iostat=v)
!!$   if (v /= 0) exit
!!$   i=i+1
!!$end do
!!$lnum=i
!!$rewind(12)
!!$
!!$allocate(lambda_arr(0:lnum-1))
!!$
!!$do i=0,lnum-1
!!$
!!$   read(12,*) lambda_arr(i)
!!$
!!$enddo
!!$
!!$close(12)
!!$
!!$end subroutine read_lambda_list

!> Sets the correct stellar luminosity
subroutine set_magtar
  ! This routines sets the correct opacity and stellar luminosity
  real(kind=real64) :: lambda_swap, lumstar

!!$  lambda_swap=lambda
!!$  
!!$  lambda=1.  ! reference lambda 1 um 
!!$  call read_opacity_table_trustI
!!$  dens=dens/kext
!!$
!!$!! rescaling opacity and source luminosity at the correct lambda
!!$
!!$  lambda= lambda_swap
!!$  
!!$  call read_opacity_table_trustI
!!$  dens=dens*kext  ! rescaling opacity with correct lambda 

!! assign stellar luminosity   
  call read_assign_param_src()

!!$  kabs=kabs/kext  !!! normalization : the dens array values in the grid already contain extinction coefficient 
!!$  ksca=ksca/kext
!!$  kext=1.

  
end subroutine set_magtar

subroutine read_assign_param_src
! this reads stellar source parameters 

logical :: file_exists
integer :: i,j,v ,outcube(3)
integer :: il 
real(kind=real64), parameter :: sigma_SB=5.6705085E-8, lsun=3.8268000E26, cspeed=2.9979246e+08 ! [MKS] 
real(kind=real64), allocatable :: teff_arr(:), lbol_arr(:), lstar_arr(:)

! open file with source parameters  

allocate(teff_arr(0:tot_p_src-1),lbol_arr(0:tot_p_src-1), lstar_arr(0:tot_p_src-1)) 

open(1,file=trim(adjustl(dir_grid))//trim(adjustl(file_param_src)),status='old')

read(1,*)

do i=0, tot_p_src-1

read(1,*) teff_arr(i),lbol_arr(i)

enddo 

close(1)

lbol_arr=lbol_arr*lsun    ! from Lsun to W/s 

! determine stellar luminosity at reference wavelength 

do il = 0, lnum-1 

   lambda = lambda_arr(il)

   do i=0,tot_p_src-1
    
      lstar_arr(i)=bplanck(teff_arr(i),lambda*1E-6)*pi*lbol_arr(i)/(sigma_SB*teff_arr(i)**4)   ! W/s/m
   
   enddo

   lum_p_src=lstar_arr*1E7*(lambda*1E-6)**2./cspeed  ! from W/s/m to erg/s/Hz
 
   print *, 'source luminosities'
   print *, 'lambda = ', lambda
   print *, lum_p_src

   lum_p_src_arr(il,:) = lum_p_src
   
enddo

deallocate(teff_arr, lbol_arr, lstar_arr)

end subroutine read_assign_param_src

END MODULE user_routines_magtar
