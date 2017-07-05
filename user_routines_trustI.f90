MODULE user_routines_trustI
use smooth_grid_routines
use iso_fortran_env
use io_routines
!use dust_model_trustI
IMPLICIT NONE 

! TRUST I parameters
!> @param tau_z Vertical optical depth of the slab at the reference wavelength. 
real(kind=real64) :: tau_z
!> @param ndust Dust number density / Extinction coefficient 
real(kind=real64) :: ndust

real(kind=real32),parameter :: lz_slab=3.
real(kind=real32),parameter :: z0_slab=-5, z1_slab=-2 ! pc 
real(kind=real32),parameter :: x0_slab=-5, x1_slab=5 ! pc 
real(kind=real32),parameter :: y0_slab=-5, y1_slab=5 ! pc 
real(kind=real32),parameter :: x_star=0, y_star=0, z_star=4 ! pc

!> @param min_lvl_in Minimum subdivision level within the slab for the "standard" subdivision_criteria(). 
integer :: min_lvl_in

!> @param subdivision_criteria Subdivision criteria to be used. Choices are:
!> - 'standard': subdivision IF subdivision level < max_lvl() AND (cell upper border coinciding with slab upper border OR (cell optical depth > max_dtau AND subdivision level < min_lvl_in()) OR subdivision level < min_lvl())  
character(LEN=lcar) :: subdivision_criteria

! input namelist 
namelist /trustI_input_strings/ dir_grid, grid_file, grid_info_file,  file_lambda_list, units_lambda,dust_model,dust_opacity_tables, file_gra_fa, file_sil_fa, file_pah_neu_fa, file_pah_ion_fa, file_av_opacities, file_q_gra, file_q_sil, file_q_pah_neu, file_q_pah_ion, subdivision_criteria

namelist /trustI_input_var/ modelsize, tau_z, lambda_ref, base,  max_ncell, max_lvl, min_lvl, max_dtau, max_dlum, n_dust_size_qabs, n_dust_wave_qabs, min_lvl_in

 namelist /trustI_input_logical/ input_av_opacities

CONTAINS

!> Reads input file for the grid creation for the TRUST I slab benchmark.
subroutine read_input_trustI

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
read(ioun, nml = trustI_input_strings)
read(ioun, nml = trustI_input_var)
read(ioun, nml = trustI_input_logical)
close(ioun)

! check dust model
call check_input_dust_model

!dust_model = 'TRUST'
!dust_opacity_tables = 'TRUST'
!iq_dust_model = [1,1,1,0]

! set grid creation true 
grid_creation = .TRUE.

call print_done

end subroutine read_input_trustI




!> Set the slab density given the vertical optical depth at the reference wavenegth
subroutine set_slab_density()

  kext_ref=kext_ref/(parsec**2)   ! m^2 --> pc ^2 
  
  ndust=tau_z/(kext_ref*lz_slab)   ! pc ^-3

  ndust=ndust*kext_ref    ! pc ^-1  this is the extinction coefficient. This passage seems useless now but it is important when scaling at wavelengths different than lambda = 1 um 

!  call read_planck_table_trustI(lambda,lumstar)   ! do not confuse lumstar with lum_stars !! lumstar is the luminosity of the single star in the trustI benchmark


end subroutine set_slab_density

subroutine read_planck_table_trustI(lumstar)
IMPLICIT NONE
!INTEGER, PARAMETER :: dp= selected_int_kind(16), opnum=16
!INTEGER, parameter :: opnum=16
character(LEN=lcar)  :: file_lum
real (Kind=real64), allocatable :: tlambda(:), tlum(:),tlum2(:)
real (kind=real64) :: lumstar
integer :: i, j, iq(0:0),v,opnum
real (Kind=real64) :: m_t,c_t

file_lum=trim(adjustl(dir_grid))//'BB_T10000_L100000.dat' ! 1) lambda (um), 2) luminosity (erg /s/Hz) , 3) (erg /s /um)

open(12,file=file_lum, status='old')


i=-7 ! the first seven lines are comments 
do 
read(12,*,iostat=v)
if (v /= 0) exit
i=i+1
end do
opnum=i
rewind(12)

allocate(tlambda(0:opnum-1), tlum(0:opnum-1), tlum2(0:opnum-1))

do i=1,6    ! skip file header 
read (12,*)
end do

do i=0,opnum-1

! um    
read(12,*) tlambda(i), tlum(i),tlum2(i)

enddo 

close(12)


!!$iq=minloc(abs(tlambda-lambda))-1  ! arrays defined from 0
!!$
!!$kabs=tkabs(iq(0))
!!$ksca=tksca(iq(0))
!!$kext=tkext(iq(0))
!!$gsca=tgsca(iq(0))

if (lambda < tlambda(0)) then   
   print *, 'lambda smaller than tabulated lambda in table'
   stop
elseif (lambda > tlambda(opnum-1)) then
   print *, 'lambda greater than tabulated lambda in table'
   stop
else 
iq=minloc(abs(tlambda-lambda))    !!! index tabulated value for lambda    
iq=iq-1  !  this is only because of IDL type index convention 
if (lambda <= tlambda(iq(0)) ) iq=iq-1   ! this makes iq the index of the tabulated value just at the left of the real one 
j=iq(0)

tlambda=dlog10(tlambda)
tlum=dlog10(tlum)

lambda=dlog10(lambda)

m_t=(tlum(j+1)-tlum(j))/(tlambda(j+1)-tlambda(j))
c_t=tlum(j+1)-m_t*tlambda(j+1)

lumstar=m_t*lambda+c_t
lumstar=10.**lumstar

lambda=10**lambda

endif

end subroutine read_planck_table_trustI

real(kind=real64) function av_rho_dust_slab(x,y,z,cellsize)
IMPLICIT NONE 
!INTEGER, PARAMETER :: dp= selected_int_kind(16)
INTEGER, PARAMETER :: nstep=1  !!! the cell is uniform 
real(Kind=real64) :: x,y,z,cellsize,xmin,xmax,ymin,ymax,zmin,zmax
real(Kind=real64) :: xp(0:nstep-1),yp(0:nstep-1),zp(0:nstep-1),sx,sy,sz,tot_ndust,rho_dust,cvol
integer :: i,j,k
logical :: dust_x_lim, dust_y_lim, dust_z_lim

xmin=x-cellsize/2.
xmax=x+cellsize/2.
ymin=y-cellsize/2.
ymax=y+cellsize/2.
zmin=z-cellsize/2.
zmax=z+cellsize/2.

sx=cellsize/nstep
sy=sx
sz=sx


do i=0,nstep-1
xp(i)=xmin+(i+0.5)*sx
yp(i)=ymin+(i+0.5)*sy
zp(i)=zmin+(i+0.5)*sz
end do

tot_ndust=0.
cvol=sx*sy*sz


do i=0,nstep-1
 do j=0,nstep-1
  do k=0, nstep-1

     dust_z_lim=(zp(k) >= z0_slab .and. zp(k) <=z1_slab)
     dust_x_lim=(xp(i) >= x0_slab .and. xp(i) <= x1_slab)
     dust_y_lim=(yp(j) >= y0_slab .and. yp(j) <= y1_slab)
     
if ((dust_x_lim).and.(dust_y_lim).and.(dust_z_lim)) then 
   rho_dust=ndust
else
rho_dust=0
endif

tot_ndust=tot_ndust+rho_dust

end do
end do
end do

av_rho_dust_slab=tot_ndust*cvol/cellsize**3


end function av_rho_dust_slab

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

!> This routines sets the correct stellar luminosity.
subroutine set_trustI
  real(kind=real64) :: lumstar
  integer :: i

  if (main_prc) print *, 'setting opacity and source luminosity....'
  
  !lambda=lambda_ref  ! reference lambda 1 um 
  !call read_opacity_table_trustI
  !dens_ref=dens_ref/(kext_ref)  ! see load_kext_arr

!! rescaling opacity and source luminosity at the correct lambda
 
 ! call load_kext_arr
  
  !do i = 0, lnum-1

 !    dens_arr(i,:)=dens_ref(:)*kext_arr(i)  ! rescaling opacity with correct lambda 

 ! enddo 
  
  !! assign stellar luminosity

  do i = 0, lnum-1

  lambda = lambda_arr(i)   
  call read_planck_table_trustI(lumstar)

  lum_p_src_arr(i, 0)=lumstar   !! 0 is the index of the only point source in this model

  enddo
  
 ! kabs_arr=kabs_arr/kext_arr  !!! normalization : the dens array values in the grid already contain extinction coefficient 
 ! ksca_arr=ksca_arr/kext_arr
 ! kext_arr=1.

  call print_done
  
end subroutine set_trustI

END MODULE user_routines_trustI
