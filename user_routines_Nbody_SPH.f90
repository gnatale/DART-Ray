!> Contains the subroutines to specify the stellar emission and dust density distribution for N-body/SPH galaxy simulations. 
MODULE user_routines_Nbody_SPH
use smooth_grid_routines
use iso_fortran_env
use HDF5
use io_routines
use sed_routines
IMPLICIT NONE 

!> @param lambda_min Minimum wavelength [um] for which the lambda grid has to be calculated. 
real(kind=real64) :: lambda_min
!> @param lambda_max Maximum wavelength [um] for which the lambda grid has to be calculated. 
real(kind=real64) :: lambda_max

!> @param nstep_g Number of grid elements per axis used to subdivide the entire RT model when calculating the grid index arrays sec_ind_star() and sec_ind_gas().
INTEGER, PARAMETER, private :: nstep_g=10

!> @param sec_ind_star(0 :nstep_g-1, 0:nstep_g-1, 0:nstep_g-1)%(:) Contains subscripts of stellar particles located within the sector ix,iy,iz of the model. Useful for quicker search of particles while calculating the main grid.  
type(var_arr_1d), allocatable :: sec_ind_star(:,:,:)

!> @param sec_ind_gas(0 :nstep_g-1, 0:nstep_g-1, 0:nstep_g-1)%(:) Contains subscripts of gas particles located within the sector ix,iy,iz of the model. Useful for quicker search of particles while calculating the main grid.  
type(var_arr_1d), allocatable :: sec_ind_gas(:,:,:)

!> @param k_sec_ind_star(0 :nstep_g-1, 0:nstep_g-1, 0:nstep_g-1) Number of stellar particles present in each model sector 
integer, allocatable :: k_sec_ind_star(:,:,:)

!> @param k_sec_ind_gas(0 :nstep_g-1, 0:nstep_g-1, 0:nstep_g-1) Number of gas particles present in each model sector 
integer, allocatable :: k_sec_ind_gas(:,:,:)

!> @param gastemp_limit Max gas temperature for gas containing dust particles. Default value = 1E6 K. 
real(Kind=real64) :: gastemp_limit = 1E6

!> Data type used for pcell_index arrays (pcell_index_star() and pcell_index_gas()).
Type Irr_arr
   integer, allocatable :: plist(:)
end type Irr_arr

!> @param pcell_index_star Each element contains the list of stellar particles in the corresponding grid cell. 
type(Irr_arr), allocatable :: pcell_index_star(:)

!> @param pcell_index_gas Each element contains the list of gas particles in the corresponding grid cell. 
type(Irr_arr), allocatable :: pcell_index_gas(:)

!> @param subdivision_criteria Subdivision criteria to be used. Choices are:
!> - 'standard': subdivision IF subdivision level < max_lvl() AND (cell optical depth > max_dtau()*optical depth averaged over entire model OR stellar luminosity > max_dlum()*total model luminosity OR  subdivision level < min_lvl())  
character(LEN=lcar) :: subdivision_criteria


! input namelist 
namelist /Nbody_sph_input_strings/ label_model_lambda_grid,dir_grid, grid_file, grid_info_file,  file_lambda_list, units_lambda,dust_model,dust_opacity_tables, file_gra_fa, file_sil_fa, file_pah_neu_fa, file_pah_ion_fa, file_av_opacities, file_q_gra, file_q_sil, file_q_pah_neu, file_q_pah_ion,file_nbody_sph,stellar_library, file_stellar_library, subdivision_criteria

namelist /Nbody_sph_input_var/ modelsize, lambda_ref,  lambda_min, lambda_max,base,  max_ncell, max_lvl, min_lvl, max_dtau, max_dlum, n_dust_size_qabs, n_dust_wave_qabs, z_sun,gastemp_limit

namelist /Nbody_sph_input_logical/ input_av_opacities

CONTAINS

!> Reads input file for grid creation. 
subroutine read_input_Nbody_SPH

integer :: ioun
character(len=lcar) :: input_filename

! initialise lambda_min and lambda_max 
lambda_min = 0   ! [um] 
lambda_max = 3000. 

call get_command_argument(1,input_filename)

!check input_file argument has been given
if (trim(input_filename) == '') then 
   if (main_prc) print *, 'provide input file in command line...EXIT' 
   call stop_prc
endif

! read input file 
if (main_prc) print *, 'reading input file...'

open(newunit=ioun, file = input_filename, status ='old', action ='read')
read(ioun, nml = Nbody_sph_input_strings)
read(ioun, nml = Nbody_sph_input_var)
read(ioun, nml = Nbody_sph_input_logical)
close(ioun)


! check dust model
call check_input_dust_model

! set grid creation true 
grid_creation = .TRUE.

! set pcell output file 
file_pcell = 'grid_'//trim(adjustl(label_model_lambda_grid))//'_pcell.h5'



call print_done

end subroutine read_input_Nbody_SPH

!> Sets extinction coefficient kext in units of parsec^2/Msun 
subroutine set_kext_gas(lambda_in)
  real(kind=real64) :: lambda_in, conv
  integer :: il 

  call find_lambda_index(lambda_in, il)

  conv=1./(1.4*m_h*parsec**2)*msun   !!! this formula assumes n_tot=n_h+n_he=1.1*n_h
  kext=kext_arr(il)*tot_n_dust*conv  ! parsec^2/Msun  

end subroutine set_kext_gas



!!$subroutine read_input_vgal()
!!$!!! This subroutines reads the input dust opacity and luminosity for the TRUST I benchmark
!!$  
!!$  real(Kind=real64), parameter :: m_h=1.6605402e-24, msun=1.9892000e+33, parsec=3.0856780e+18 !! cgs units
!!$  real(kind=real64) :: conv
!!$  real(kind=real64) :: x0,x1,y0,y1,z0,z1
!!$  
!!$ !  call read_opacity_table_trustI  WHAT ABOUT THIS 
!!$  !kext=kext/((parsec*1.E2)**2)   ! cm^2 --> pc ^2
!!$
!!$  conv=1./(1.4*m_h*parsec**2)*msun   !!! this formula assumes n_tot=n_h+n_he=1.1*n_h
!!$  kext=tau_nh*conv   ! parsec^2/Msun  !!! look different definition of kext compared to TRUST I grid creation program
!!$
!!$  ! read input particle grid
!!$ !!$  if (.not.allocated(g_ccoord)) then   !!! this avoids repeating this when makeing lambda grids  
!!$ !!$  call read_grid_vgal(x0,x1,y0,y1,z0,z1)
!!$ !!$  endif 
!!$  
!!$  ! convert mass to stellar luminosity 
!!$ ! call stellar_mass_to_luminosity
!!$
!!$  ! read gas metallicity table. Remember that the correction for metallicity is applied to g_mass not to kext!!! )
!!$ ! call read_gas_metallicity()
!!$
!!$  ! make index to find particles in a faster way
!!$ !!$  if (.not.allocated(sec_ind_star)) then 
!!$ !!$  call create_grid_index
!!$ !!$  print *, 'index done'
!!$ !!$  endif 
!!$  
!!$end subroutine read_input_vgal


!!$subroutine read_grid_vgal(x0,x1,y0,y1,z0,z1) 
!!$IMPLICIT NONE
!!$real(kind=real64) :: x0,x1,y0,y1,z0,z1
!!$real(kind=real64) :: a(3)
!!$integer :: i,j,v,iq(0:0)
!!$
!!$open(1,file=file_input_grid,status='old')
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
!!$tot_points=tot_points-1 ! the last row contains negative numbers (victor file) 
!!$
!!$!tot_points=4000
!!$
!!$allocate(g_ccoord(3,0:tot_points-1), g_mass(0:tot_points-1),pop(0:tot_points-1),g_age_temp(0:tot_points-1),abs_loc(0:tot_points-1),z_met(0:tot_points-1), star_lum(0:tot_points-1), pcell(0:tot_points-1) )
!!$
!!$do i=0,tot_points-1 
!!$
!!$read(1,*) pop(i), (g_ccoord(j,i),j=1,3), (a(j),j=1,3), g_mass(i),g_age_temp(i)
!!$
!!$enddo
!!$close(1)
!!$
!!$g_mass=g_mass*2.33E5 !!! Msun (scaling factor used by Victor) 
!!$g_ccoord=g_ccoord*1.E3    !!! kpc --> pc 
!!$
!!$abs_loc=0.
!!$z_met=0.
!!$star_lum=0.
!!$pcell=-1
!!$
!!$x0=minval(g_ccoord(1,0:tot_points-1))
!!$x1=maxval(g_ccoord(1,0:tot_points-1))
!!$print *, 'X', x0,x1
!!$
!!$y0=minval(g_ccoord(2,0:tot_points-1))
!!$y1=maxval(g_ccoord(2,0:tot_points-1))
!!$print *, 'Y', y0,y1
!!$
!!$z0=minval(g_ccoord(3,0:tot_points-1))
!!$z1=maxval(g_ccoord(3,0:tot_points-1))
!!$print *, 'Z', z0,z1
!!$
!!$
!!$end subroutine read_grid_vgal



!!$subroutine read_gas_metallicity
!!$
!!$IMPLICIT NONE
!!$integer :: i,v,ng,j
!!$real(KIND=real64),parameter  ::  z_sun=0.018
!!$real (Kind=real64), allocatable :: fe_h(:)
!!$real(Kind=real64) :: a, z
!!$
!!$! read table gas particle metallicities 
!!$open(1,file='../RTCODE_V7/Fe_H_sq.gas.01000main',status='old')
!!$i=0
!!$do 
!!$read(1,*,iostat=v)
!!$if (v /= 0) exit
!!$i=i+1
!!$end do
!!$ng=i
!!$rewind(1)
!!$
!!$ng = ng -1  ! first line is a comment
!!$
!!$!print *, ng
!!$
!!$allocate(fe_h(0:ng-1))
!!$
!!$read(1,*)    
!!$
!!$do i=0, ng-1
!!$
!!$read(1,*), fe_h(i), a
!!$
!!$end do 
!!$
!!$close(1)
!!$
!!$! assign metallicities (relative to solar value) to z_met 
!!$
!!$j=-1
!!$
!!$do i=0,tot_points-1 
!!$
!!$if (pop(i) == 3) then 
!!$
!!$j=j+1
!!$
!!$if (fe_h(j) > -600.) then  !!! there is a flag in the table == -66
!!$!g_mass(i)=g_mass(i)*10.**fe_h(j)  !!! here z_sun is not needed 
!!$z_met(i)=10.**fe_h(j)
!!$else
!!$!g_mass(i)=0.
!!$z_met(i)=0.
!!$endif
!!$
!!$!if (g_age_temp(i) > t_limit) g_mass(i)=0.
!!$
!!$end if 
!!$
!!$end do
!!$
!!$deallocate(fe_h)
!!$
!!$
!!$end subroutine  read_gas_metallicity

!> Creates a grid index for the stellar and gas particles so it is quicker to find particles in a certain region of space. 
subroutine create_grid_index
real(KIND=real64) :: step,x0,xi(0:nstep_g-1)
integer :: i,j,k,ip,ix,iy,iz,flagx,flagy,flagz
integer, parameter :: npart_start = 1000
integer ::  nk
integer, allocatable :: temp_arr(:)

if (main_prc) print *, 'create grid index for simulation particles...'

! allocate sec_ind_star and sec_ind_gas
! The size of these variable size arrays will be varied later if necessary

allocate(sec_ind_star(0:nstep_g-1,0:nstep_g-1,0:nstep_g-1), sec_ind_gas(0:nstep_g-1,0:nstep_g-1,0:nstep_g-1))

do ix = 0, nstep_g-1
   do iy = 0, nstep_g-1
      do iz = 0, nstep_g-1
         allocate(sec_ind_star(ix,iy,iz)%int(0:npart_start-1))
         sec_ind_star(ix,iy,iz)%int = 0
         allocate(sec_ind_gas(ix,iy,iz)%int(0:npart_start-1))
         sec_ind_gas(ix,iy,iz)%int = 0
      enddo
   enddo
end do

! allocate k_sec_ind_star and k_sec_ind_gas
allocate(k_sec_ind_star(0:nstep_g-1,0:nstep_g-1,0:nstep_g-1), k_sec_ind_gas(0:nstep_g-1,0:nstep_g-1,0:nstep_g-1))
k_sec_ind_star = 0 
k_sec_ind_gas = 0 

! define model sectors 
x0=-modelsize/2.
step=modelsize/nstep_g
do i=0,nstep_g-1 
   xi(i)=x0+i*step
enddo

! loop on stellar particles 
do ip=0, tot_star_particles-1

   flagx=0
   flagy=0
   flagz=0

   do i=nstep_g-1,0,-1

      if ((starcoord(1,ip) > xi(i)).and.(starcoord(1,ip) < xi(nstep_g-1)+step).and.(flagx == 0)) then 
         ix= i
         flagx=1
      endif

      if ((starcoord(2,ip) > xi(i)).and.(starcoord(2,ip) < xi(nstep_g-1)+step).and.(flagy ==0)) then 
         iy= i
         flagy=1
      endif

      if ((starcoord(3,ip) > xi(i)).and.(starcoord(3,ip) < xi(nstep_g-1)+step).and.(flagz ==0)) then 
         iz= i
         flagz=1
      endif

   end do

   if ((flagx == 1).and.(flagy ==1).and.(flagz ==1)) then 

      k_sec_ind_star(ix,iy,iz)=k_sec_ind_star(ix,iy,iz)+1
      k=k_sec_ind_star(ix,iy,iz)

      if (k > size(sec_ind_star(ix,iy,iz)%int)-1) then ! expand sec_ind_star if necessary 
         nk = size(sec_ind_star(ix,iy,iz)%int)
         allocate(temp_arr(0:nk-1))
         temp_arr = sec_ind_star(ix,iy,iz)%int
         deallocate(sec_ind_star(ix,iy,iz)%int)
         allocate(sec_ind_star(ix,iy,iz)%int(0:2*nk-1))
         sec_ind_star(ix,iy,iz)%int(0:nk-1) = temp_arr
         deallocate(temp_arr)
      endif
 
      sec_ind_star(ix,iy,iz)%int(k)=ip

   endif

end do

! loop on gas particles 
do ip=0, tot_gas_particles-1

   flagx=0
   flagy=0
   flagz=0

   do i=nstep_g-1,0,-1

      if ((gascoord(1,ip) > xi(i)).and.(gascoord(1,ip) < xi(nstep_g-1)+step).and.(flagx == 0)) then 
         ix= i
         flagx=1
      endif

      if ((gascoord(2,ip) > xi(i)).and.(gascoord(2,ip) < xi(nstep_g-1)+step).and.(flagy ==0)) then 
         iy= i
         flagy=1
      endif

      if ((gascoord(3,ip) > xi(i)).and.(gascoord(3,ip) < xi(nstep_g-1)+step).and.(flagz ==0)) then 
         iz= i
         flagz=1
      endif

   end do

   if ((flagx == 1).and.(flagy ==1).and.(flagz ==1)) then 

      k_sec_ind_gas(ix,iy,iz)=k_sec_ind_gas(ix,iy,iz)+1
      k=k_sec_ind_gas(ix,iy,iz)

      if (k > size(sec_ind_gas(ix,iy,iz)%int)-1) then ! expand sec_ind_gas if necessary 
         nk = size(sec_ind_gas(ix,iy,iz)%int)
         allocate(temp_arr(0:nk-1))
         temp_arr = sec_ind_gas(ix,iy,iz)%int
         deallocate(sec_ind_gas(ix,iy,iz)%int)
         allocate(sec_ind_gas(ix,iy,iz)%int(0:2*nk-1))
         sec_ind_gas(ix,iy,iz)%int(0:nk-1) = temp_arr
         deallocate(temp_arr)
      endif
 
      sec_ind_gas(ix,iy,iz)%int(k)=ip

   endif

end do

call print_done

end subroutine create_grid_index

subroutine av_galaxy(x,y,z,cellsize,av_dust,av_stars,av_gas,av_st_mass,cc)
IMPLICIT NONE 

real(Kind=real64) :: x,y,z,x0,xmin,xmax,ymin,ymax,zmin,zmax,tot_mass_dust,cellsize,step,xi(0:nstep_g-1),av_dust,av_stars,av_st_mass, p_met
real(kind=real64) :: tot_lum_stars,tot_mass_gas,av_gas,tot_mass_stars,cellsize3
integer :: i,flagxmax,flagymax,flagzmax,flagxmin,flagymin,flagzmin,ixmin,ixmax,iymin,iymax,izmin,izmax,ix,iy,iz,k,cc,totflag,np

xmin=x-cellsize/2.
xmax=x+cellsize/2.
ymin=y-cellsize/2.
ymax=y+cellsize/2.
zmin=z-cellsize/2.
zmax=z+cellsize/2.

cellsize3=cellsize**3

tot_mass_dust=0.   !!! this is actually gas mass multiplied by Z/Z_sun 
tot_mass_gas=0.
tot_lum_stars=0.
tot_mass_stars=0.

if (.not.allocated(pcell_index_star)) then  !!! this is when creating main grid 

   ! define sectors 
   x0=-modelsize/2.
   step=modelsize/nstep_g
   do i=0,nstep_g-1 
      xi(i)=x0+i*step
   enddo

   flagxmax=0
   flagxmin=0
   flagymax=0
   flagymin=0
   flagzmax=0
   flagzmin=0

   ixmax=nstep_g-1
   ixmin=0
   iymax=nstep_g-1
   iymin=0
   izmax=nstep_g-1
   izmin=0

   ! find interval in index grid 
   do i=nstep_g-1,0,-1
      if ((xmax > xi(i)).and.(flagxmax == 0)) then 
         ixmax= i
         flagxmax=1
      endif

      if ((xmin > xi(i)).and.(flagxmin == 0)) then 
         ixmin= i
         flagxmin=1
      endif

      if ((ymax > xi(i)).and.(flagymax == 0)) then 
         iymax= i
         flagymax=1
      endif

      if ((ymin > xi(i)).and.(flagymin == 0)) then 
         iymin= i
         flagymin=1
      endif

      if ((zmax > xi(i)).and.(flagzmax == 0)) then 
         izmax= i
         flagzmax=1
      endif

      if ((zmin > xi(i)).and.(flagzmin == 0)) then 
         izmin= i
         flagzmin=1
      endif

      totflag=flagxmax+flagxmin+flagymax+flagymin+flagzmax+flagzmin
      if (totflag == 6) exit 

   end do

!!$print *, 'XMIN', ixmin
!!$print *, xmin, xi(ixmin),xi(ixmin)+step
!!$print *, 'YMIN', iymin
!!$print *, ymin, xi(iymin),xi(iymin)+step
!!$print *, 'ZMIN', izmin
!!$print *, zmin, xi(izmin),xi(izmin)+step
!!$print *, 'XMAX', ixmax
!!$print *, xmax, xi(ixmax),xi(ixmax)+step
!!$print *, 'YMAX', iymax
!!$print *, ymax, xi(iymax),xi(iymax)+step
!!$print *, 'ZMAX', izmax
!!$print *, zmax, xi(izmax),xi(izmax)+step

!print *, 'x=',x,'y=',y,'z=',z
   
   !----------------------------
   ! loop for stellar particles 
   do ix=ixmin,ixmax
      do iy=iymin,iymax
         do iz=izmin,izmax

            if (k_sec_ind_star(ix,iy,iz) < 1000) then
               call omp_set_num_threads(4)
            else 
               call omp_set_num_threads(16)
            endif

            !$OMP PARALLEL DEFAULT(NONE), PRIVATE(k,i ), &
            !$OMP SHARED(ix,iy,iz, xmin,xmax,ymin,ymax,zmin,zmax,tot_mass_dust,tot_mass_gas,cchild, pcell_star,tot_lum_stars, tot_mass_stars, cc, star_lum,  gastemp_limit, sec_ind_star , k_sec_ind_star, starcoord, mstar)

            !$OMP DO SCHEDULE(DYNAMIC,20)

            do k = 0, k_sec_ind_star(ix,iy,iz)-1 

               i=sec_ind_star(ix,iy,iz)%int(k)

               if ((starcoord(1,i) < xmin).or.(starcoord(1,i) >= xmax).or.(starcoord(2,i) < ymin).or.(starcoord(2,i) >= ymax) .or.(starcoord(3,i) < zmin).or.(starcoord(3,i) >= zmax)) cycle 

               !$OMP ATOMIC
               tot_lum_stars=tot_lum_stars+star_lum(i)
   
               !$OMP ATOMIC
               tot_mass_stars=tot_mass_stars+mstar(i)

               if (cchild(cc) == -1) pcell_star(i) =cc ! this is the cell where the particle is 

            enddo

            !$OMP END DO NOWAIT
            !$OMP END PARALLEL

         enddo
      enddo
   enddo

! loop for gas particles 
   do ix=ixmin,ixmax
      do iy=iymin,iymax
         do iz=izmin,izmax

            if (k_sec_ind_gas(ix,iy,iz) < 1000) then
               call omp_set_num_threads(4)
            else 
               call omp_set_num_threads(16)
            endif

            !$OMP PARALLEL DEFAULT(NONE), PRIVATE(k,i, p_met ), &
            !$OMP SHARED(ix,iy,iz, xmin,xmax,ymin,ymax,zmin,zmax,tot_mass_dust,tot_mass_gas,cchild, pcell_gas, cc, gastemp_limit, sec_ind_gas , k_sec_ind_gas, gascoord, gastemp, ofegas, fehgas, mgas)

            !$OMP DO SCHEDULE(DYNAMIC,20)

            do k = 0, k_sec_ind_gas(ix,iy,iz)-1 

               i=sec_ind_gas(ix,iy,iz)%int(k)

               if ((gascoord(1,i) < xmin).or.(gascoord(1,i) >= xmax).or.(gascoord(2,i) < ymin).or.(gascoord(2,i) >= ymax) .or.(gascoord(3,i) < zmin).or.(gascoord(3,i) >= zmax)) cycle 

               if (gastemp(i) < gastemp_limit) then   !!! temperature cut for dust presence
                  p_met = 10.**(ofegas(i)+fehgas(i))
                  !$OMP ATOMIC
                  tot_mass_dust=tot_mass_dust+p_met*mgas(i)   !!! this is  Mgas*Z/Z_sun (no need to multiply by z_sun here).
               endif
               !$OMP ATOMIC
               tot_mass_gas=tot_mass_gas+mgas(i)

               !! this is to assign the cell to which the particle belong. I don't think cchild(cc) =1 is needed here. Check before using! 
               if (cchild(cc) == -1) pcell_gas(i) =cc ! this is the cell where the particle is 

            enddo

            !$OMP END DO NOWAIT
            !$OMP END PARALLEL

         enddo
      enddo
   enddo

else   !!! This is the option when creating the the lambda grids 

   ! loop on star particles 
   np=size(pcell_index_star(cc)%plist)

   !$OMP PARALLEL DEFAULT(NONE), PRIVATE(k,i ), &
   !$OMP SHARED(cchild,tot_lum_stars, tot_mass_stars, cc, star_lum, np , pcell_index_star, mstar)

   !$OMP DO SCHEDULE(DYNAMIC,20)

   do k=0, np-1

      i= pcell_index_star(cc)%plist(k)
  
      !$OMP ATOMIC
      tot_lum_stars=tot_lum_stars+star_lum(i)
   
      !$OMP ATOMIC
      tot_mass_stars=tot_mass_stars+mstar(i)
!!$ 
!!$endif

   enddo

   !$OMP END DO NOWAIT
   !$OMP END PARALLEL

   ! loop on gas particles 

   np=size(pcell_index_gas(cc)%plist)

   !$OMP PARALLEL DEFAULT(NONE), PRIVATE(k,i, p_met ), &
   !$OMP SHARED(tot_mass_dust,tot_mass_gas,cchild, cc, np , pcell_index_gas, gastemp_limit, ofegas, fehgas, mgas, gastemp)

   !$OMP DO SCHEDULE(DYNAMIC,20)

   do k=0, np-1

      i= pcell_index_gas(cc)%plist(k)

      if (gastemp(i) < gastemp_limit) then   !!! temperature cut for dust presence
         p_met = 10.**(ofegas(i)+fehgas(i))
         !$OMP ATOMIC
         tot_mass_dust=tot_mass_dust+p_met*mgas(i)   !!! this is actually Mgas*Z/Z_sun     
      endif
      !$OMP ATOMIC
      tot_mass_gas=tot_mass_gas+mgas(i)

   enddo

   !$OMP END DO NOWAIT
   !$OMP END PARALLEL
   
endif

!print *,'GAS MASS=',tot_mass_dust
!print *, 'STAR LUM=',tot_lum_stars
av_dust=tot_mass_dust/cellsize3
av_stars=tot_lum_stars/cellsize3
av_gas=tot_mass_gas/cellsize3
av_st_mass=tot_mass_stars/cellsize3

if (av_stars < 0) then 
   print *, 'negative av_stars'
   print *, av_stars
   stop
endif

end subroutine av_galaxy

!!$subroutine set_Nbody_SPH
!!$  ! This routines sets the opacity coefficients 
!!$  
!!$!  call read_opacity_table_trustI  WHAT ABOUT THIS 
!!$  
!!$  kabs=kabs/kext  !!! normalization : the dens array values in the grid already contain extinction coefficient 
!!$  ksca=ksca/kext
!!$  kext=1.
!!$  
!!$end subroutine set_Nbody_SPH

!> Prints the pcell_star() and pcell_gas() arrays in an output file (file_pcell() ).  
subroutine print_pcell
  
  integer, parameter :: narr = 2
  CHARACTER (LEN=lcar) :: dsetname(0:narr-1) ! Dataset name
  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier 
  INTEGER(HSIZE_T) :: dims_arr(0:narr-1), dims(0:0)   ! Dataset dimensions    
  INTEGER     ::   rank_arr(0:narr-1), rank    ! Dataset rank
  INTEGER     ::   error  ! Error flag
  integer :: i, j 

  if (main_prc) print *, 'printing pcell arrays...' 

  dsetname(0) = 'pcell_star'   ; rank_arr(0)=1  ; dims_arr(0) = tot_star_particles  
  dsetname(1) = 'pcell_gas'    ; rank_arr(1)=1  ; dims_arr(1) = tot_gas_particles

  !    Initialize FORTRAN interface.
  CALL h5open_f (error)
  CALL h5fcreate_f(trim(adjustl(dir_grid))//trim(adjustl(file_pcell)), H5F_ACC_TRUNC_F, file_id, error) 

  do i = 0, narr-1 

     dims = dims_arr(i)
     rank = rank_arr(i)

     CALL h5screate_simple_f(rank, dims, dspace_id, error) 

     CALL h5dcreate_f(file_id, dsetname(i),H5T_NATIVE_INTEGER , dspace_id, &
          dset_id, error)
    
     CALL h5dopen_f(file_id, dsetname(i), dset_id, error)

     if (dsetname(i) == 'pcell_star' ) then 
        CALL h5dwrite_f(dset_id,H5T_NATIVE_INTEGER , pcell_star, dims, error)
     elseif(dsetname(i) == 'pcell_gas') then 
        CALL h5dwrite_f(dset_id,H5T_NATIVE_INTEGER , pcell_gas, dims, error)
     endif

     CALL h5dclose_f(dset_id, error)

     CALL h5sclose_f(dspace_id, error)

  enddo

  CALL h5fclose_f(file_id, error)

  CALL h5close_f(error)

  call print_done

end subroutine print_pcell

!> Reads pcell_star() and pcell_gas() from file_pcell(). 
subroutine read_pcell
  
  integer, parameter :: narr = 2
  CHARACTER (LEN=lcar) :: dsetname(0:narr-1) ! Dataset name
  INTEGER(HID_T) :: file_id                            ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
  INTEGER(HSIZE_T)  :: dims_arr(0:narr-1), dims(0:0)  ! Dataset dimensions    
  INTEGER     ::   error  ! Error flag
  integer :: i
 
  if (main_prc) print *, 'reading pcell arrays...' 
  
  dsetname(0) = 'pcell_star'    ; dims_arr(0) = tot_star_particles  
  dsetname(1) = 'pcell_gas'     ; dims_arr(1) = tot_gas_particles
  
  CALL h5open_f (error)
     
  CALL h5fopen_f (trim(adjustl(dir_grid))//trim(adjustl(file_pcell)), H5F_ACC_RDWR_F, file_id, error)
    
  do i = 0, narr -1 

     dims = dims_arr(i)

     CALL h5dopen_f(file_id, dsetname(i), dset_id, error)
     
     if (dsetname(i) == 'pcell_star' ) then    
        CALL h5dread_f(dset_id,H5T_NATIVE_INTEGER , pcell_star, dims, error)
     elseif (dsetname(i) == 'pcell_gas' ) then
        CALL h5dread_f(dset_id,H5T_NATIVE_INTEGER , pcell_gas, dims, error)
     endif
     
     CALL h5dclose_f(dset_id, error)
     
  enddo
    
  CALL h5fclose_f(file_id, error)
  
  CALL h5close_f(error)
  
  call print_done
  
end subroutine read_pcell

!> Creates an index array for each cell pointing to the particles contained in that cell. Using this index, the lambda grid creation is faster.
subroutine make_pcell_index

  integer :: i, np,ip,k,i0,i1
  integer, allocatable :: order(:),pcell_new(:)

  if (main_prc) print *, 'deriving pcell_index...'

  allocate(pcell_index_star(0:tot_ncell-1), pcell_index_gas(0:tot_ncell-1))

  ! ----------------------------------------
  ! make pcell_index for star particles 
  allocate(order(0:tot_star_particles-1), pcell_new(0:tot_star_particles-1))
  pcell_new=pcell_star 

  if (main_prc) print *, '--sorting pcell_star'
  call quick_sort_int(pcell_new, order)    

  order=order-1

  if (main_prc) print *, '--make pcell_star_index'

  !!! find first element 
  call find_first_element(pcell_new, i0)

!!! fill pcell_index_star 
  do i = 0, tot_ncell -1
   
     if (cchild(i) == -1.and.(dens_stars_ref(i)>0)) then 
 
!!$        print *, ccoord(:,i)
!!$        print *, csize(i)

        call find_i_range(pcell_new,i,i0,i1)
        
        np=i1-i0   !!! total particles in cell i 

        allocate(pcell_index_star(i)%plist(0:np-1))       

        pcell_index_star(i)%plist=order(i0:i1-1)

!!$        do k=i0,i1-1
!!$
!!$           print *, starcoord(:,order(k))
!!$
!!$        enddo      
       
        i0=i1 
       
     end if

  end do

  deallocate(order, pcell_new)

  ! ----------------------------------------
  ! make pcell_index for gas particles 
  allocate(order(0:tot_gas_particles-1), pcell_new(0:tot_gas_particles-1))
  pcell_new=pcell_gas 

  if (main_prc) print *, '--sorting pcell_gas'
  call quick_sort_int(pcell_new, order)    

  order=order-1

  if (main_prc) print *, '--make pcell_gas_index'

  !!! find first element 
  call find_first_element(pcell_new, i0)

!!! fill pcell_index_gas 
  do i = 0, tot_ncell -1
     ! print *, 'index cell =', i

     if (cchild(i) == -1.and.(dens_ref(i)>0)) then 
 
!!$        print *, ccoord(:,i)
!!$        print *, csize(i)
!!$        read(*,*)

        call find_i_range(pcell_new,i,i0,i1)
      
        np=i1-i0   !!! total particles in cell i 

        allocate(pcell_index_gas(i)%plist(0:np-1))       

        pcell_index_gas(i)%plist=order(i0:i1-1)

!!$        do k=i0,i1-1
!!$
!!$           print *, gascoord(:,order(k))
!!$
!!$        enddo      
       
        i0=i1 
       
     end if

  end do

  deallocate(order, pcell_new)

  call print_done 

end subroutine make_pcell_index

!> finds first non zero element of an integer array. 
subroutine find_first_element(arr, i0)
  integer :: arr(0:)
  integer :: i0, num 
 
  num = size(arr)

  i0=0

  do while (i0 <= num-1)
     if (arr(i0) >= 0) then        
        exit
     endif
     i0=i0+1
  end do
  

end subroutine find_first_element

!> Finds range of subscripts where an integer array has the value "i". It requires that the array is sorted and that i0 is the subscript of an element before the first with value "i". 
subroutine find_i_range(arr,i,i0,i1)
  integer :: arr(0:)
  integer :: i, i0, i1 
  integer :: num

  num = size(arr)
  
  if (arr(i0) /= i) then 
     do while (i0 < num -1)
        if (arr(i0) == i) exit 
        i0=i0+1
     end do
  endif
  
  i1=i0

  do while (i1 < num -1)
     if (arr(i1) /= i) exit
     i1=i1+1
  end do

end subroutine find_i_range


!> This routine assigns the dens values to subdivided parent cells. This is necessary in the lambda grid creation because in that case only the leaf cells have the dens values already calculated 
subroutine assign_dens_to_parent
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




END MODULE user_routines_Nbody_SPH
