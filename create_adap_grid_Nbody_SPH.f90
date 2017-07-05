PROGRAM create_adap_grid_Nbody_SPH
  ! This program creates the main 3D grid for the TRUST I benchmark
  ! The model parameters can be expressed below and user-defined parameters and function have to be expressed in user_routines module
  ! Smooth_grid_routines contains the routines to subdivide the grid 
  use iso_fortran_env
  use smooth_grid_routines
  use user_routines_Nbody_SPH
  use io_routines
  use sed_routines
  IMPLICIT NONE
  integer :: il,i0(1:1),i1(1:1)
  logical :: file_exists
  integer :: ierr 

  !! Initialize MPI 
  call initialize_mpi

  !! INPUT VARIABLE INITIALIZE  
  call input_initialize

  ! read input file 
  call read_input_Nbody_SPH

  print *, 'Are the output filenames OK? (Press RETURN to continue)'
  print *, grid_file, grid_info_file
  read(*,*)

  ! read lambda grid 
  call read_lambda_list
  call set_lambda_arr_si

  ! set opacity values 
  call prepare_dust_model

  ! read and interpolate stellar emission library 
  call prepare_stellar_emission_library

  ! read simulation file
  call read_nbody_sph_simulation
 
  ! create particle grid index 
  call create_grid_index

 ! check if main grid already exists. If not create it!  
  
  inquire(file=trim(adjustl(dir_grid))//trim(adjustl(grid_file)), exist=file_exists)

  select case (file_exists)

  case (.FALSE.)
     print *, 'Main grid does not exist. Create it !'

     ! set stellar particles luminosity at lambda = lambda_ref 
     call set_star_particle_luminosity(lambda_ref)

     ! set kext_gas
     call set_kext_gas(lambda_ref)
    
     ! --------------------------------------------------------------------! 
     ! create grid arrays     
     call create_grid_arrays()
 
     ! set base 
     call set_base

     ! --------------------------------------------------------------------!
     ! start subdivision loop    
     call subdivision_loop()

     ! multiply by kext to make it consistent with other galaxy calculations 
     dens=dens*kext  !!! NOTE: THE POSITION OF THIS LINE IS different from the other similar grid creation codes. It is here in case one wants to print gas density separately.  

     ! --------------------------------------------------------------------!
     ! print 3D grid on H5DF file      
     call print_3d_grid_file()  
 
     ! print info file 
     call print_info_file

     ! print pcell file (list of cells to which particles belong)
     call print_pcell

     print *, 'Main grid creation completed!'
     stop

     ! -----------------------------------------------
     ! if so, read it ! 
   case(.TRUE.)
      print *, 'Main grid does already exist. Read it !'
      call read_main_grid
      call read_pcell

   end select

   ! make cell-particles index so it is faster to find the right particles when creating lambda grid 
   call make_pcell_index
   
  ! ---------------------------------------------------------------------!
  ! loop on lambdas to make monochromatic emissivity/ density tables

   ! make lambda grids 
  i0=minloc(abs(lambda_arr-lambda_min)/lambda_min)-1
  i1=minloc(abs(lambda_arr-lambda_max)/lambda_max)-1

  do il=i0(1), i1(1)
   
     lambda=lambda_arr(il)
     print *, 'lambda= ', lambda
      
     write(label_wave,'(F9.3)') lambda_arr(il) ! this converts number to string 

     ! set filename lambda grid
     grid_file_lambda='grid_'//trim(adjustl(label_model_lambda_grid))//'_l'//trim(adjustl(label_wave))//'um.h5'
     
      ! set stellar particles luminosity at lambda  
     call set_star_particle_luminosity(lambda)  ! lambda NOT lambda_ref

     ! set kext_gas
     call set_kext_gas(lambda)
      
     ! create grids
     call make_lambda_grid 

     ! assign values to parent cells 
     call assign_dens_to_parent
     
     ! print them
     call print_lambda_grid

     print *, 'grid for lambda = ', lambda, ' printed' 
     call print_done

  end do

  
CONTAINS

  subroutine make_lambda_grid
    ! this subroutine makes grids and the selected lambda
    real(kind=real64) :: x,y,z,cellsize,av_rho_dust,av_rho_stars
    integer :: cc,i

    print *, 'calculating grid for lambda =', lambda 

    call OMP_SET_NESTED(.true.)
    call omp_set_num_threads(4)
    
    !$OMP PARALLEL DEFAULT(NONE), PRIVATE(i,x,y,z,cellsize,av_rho_dust,av_rho_stars,cc), &
!$OMP SHARED(tot_ncell,csize,dens,dens_stars,ccoord)

    !$OMP DO SCHEDULE(STATIC,50)

    do i=0, tot_ncell-1
       
       !print *, 'i=',i,'of ',tot_ncell-1
       
       x=ccoord(1,i)
       y=ccoord(2,i)
       z=ccoord(3,i) 
       cellsize=csize(i)
       cc=i
       
       call calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_stars,cc)
       
       dens(i)=av_rho_dust
       
       dens_stars(i)=av_rho_stars

       if (dens_stars(i) < 0) then 
          print *, 'negative value '
          print *, i
          print *, dens_stars(i)
          stop
       endif
      
    enddo

    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    dens=kext*dens  !!! important !!! 

end subroutine make_lambda_grid

!> Obtains the average dust extinction coefficient and the average stellar luminosity of the current cell from the user-defined routines.
subroutine calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_stars,cc)

  real(Kind=real64) :: x,y,z,cellsize,av_rho_dust,av_rho_stars,av_gas, av_st_mass
  integer :: cc 

  call av_galaxy(x,y,z,cellsize,av_rho_dust,av_rho_stars,av_gas, av_st_mass,cc)

end subroutine calc_dens


logical function subdivision(x,y,z,cellsize,av_rho_dust,av_rho_stars)
  ! In this routine the cell subdivision criteria are specified 
real(kind=real64) :: x,y,z,cellsize,av_rho_dust,av_rho_stars,tau,lum_stars
integer :: ib=1 

  tau=kext*dens(ncurr)*cellsize   ! tau=K_ext*rho_N*length  here rho_N is in H atoms per volume

  lum_stars=dens_stars(ncurr)*cellsize**3

select case(subdivision_criteria)
case('standard')
   subdivision= (((lvl(ncurr) < max_lvl).and.((lum_stars > max_dlum*dens_stars(0)*modelsize**3).or.(tau > max_dtau*kext*dens(0)*modelsize).or.(nlevel < min_lvl))).or.(lvl(ncurr) < nlevel)) ! do not remove term (lvl(ncurr) < nlevel). Needed for smooth cell refinement algorithm. 
case default 
   print *, 'STOP(subdivision): subdivision_criteria not recognized!'
   print *, 'subdivision_criteria =',subdivision_criteria
   stop
end select
  
end function subdivision

  subroutine print_info_file()
! this subroutine calculates maximum and minimum tau/lum cell and print info file 
integer :: k,i
real(kind=real64) :: max_dtau_eff,max_dlum_eff,mean_dtau,mean_dlum  

! calculate max and mean 
k=0
max_dtau_eff=-999.
max_dlum_eff=-999.
mean_dtau=0
mean_dlum=0

do i=0, tot_ncell-1 
if (cchild(i) == -1) then  
k=k+1
   if (dens(i)*csize(i) > max_dtau_eff) max_dtau_eff=dens(i)*csize(i)
mean_dtau=mean_dtau+dens(i)*csize(i)
  if (dens_stars(i)*csize(i)**3 > max_dlum_eff) max_dlum_eff=dens_stars(i)*csize(i)**3
mean_dlum=mean_dlum+dens_stars(i)*csize(i)**3

endif 

end do 

mean_dtau=mean_dtau/k
mean_dlum=mean_dlum/k

! print grid info on a file
open(42, file=trim(adjustl(dir_grid))//trim(adjustl(grid_info_file)), status='unknown') 
write(42,*) ' GRID PARAMETERS'
write(42,*) ' modelsize =',  modelsize
!write(42,*) ' tau_z =', tau_z
write(42,*) ' lambda=', lambda  
!write(42,*) ' dens_stars_max', dens_stars_max
write(42,*) ' base =', base  
write(42,*) ' max_ncell', max_ncell 
write(42,*) ' MAX DTAU PARAM (input)=', max_dtau
write(42,*) ' MAX DLUM PARAM (input)=', max_dlum
write(42,*) ' MIN NLVL (input)=', min_lvl
write(42,*) ' MAX NLVL (input)=', max_lvl
write(42,*) ' MAX DTAU (output)=', max_dtau_eff
write(42,*) ' MEAN DTAU (output)=', mean_dtau
write(42,*) ' MAX DLUM (output)=', max_dlum_eff
write(42,*) ' MEAN DLUM (output)=', mean_dlum
write(42,*)

close(42)

end subroutine print_info_file


!!!!-----------------------------------!!!!!
!!! DO NOT MODIFY FOLLOWING ROUTINES !!!! 
!!!------------------------------------!!!!   

  subroutine subdivision_loop()
      ! This subroutine performs the subdivision loop until all cells satisfy the input subdivision criteria 
    real(kind=real64) :: x,y,z, cellsize, av_rho_dust,av_rho_stars
    real(kind=real64) :: lum_stars, tau
    integer :: i
    integer(kind=int32) :: a,b

    ! FIRST CELL initialization
      x=ccoord(1,0)
      y=ccoord(2,0)
      z=ccoord(3,0)
      ncurr=0
    
      call calc_cellsize(cellsize,nlevel)
    
      call calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_stars,ncurr)

      dens(0) = av_rho_dust
      dens_stars(0) = av_rho_stars
      
      do 

         subdiv:  do ncurr=nstart(nlevel),nstart(nlevel+1)-1    
! check no parent cell in the list (apart from first step) 
   if ((ncurr /= 0).and.(cchild(ncurr) /= (-1))) then 
      STOP 'there is something weird going on. ' 
      endif    ! this is needed to check problems with neighbour cell subdivision 

      x=ccoord(1,ncurr)
      y=ccoord(2,ncurr)
      z=ccoord(3,ncurr)

      call calc_cellsize(cellsize,nlevel)

      !call calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_stars)
    
      ! dens(ncurr)=av_rho_dust
     
       !dens_stars(ncurr)=av_rho_stars

      av_rho_dust=dens(ncurr)
      av_rho_stars=dens_stars(ncurr)
      
 if(subdivision(x,y,z,cellsize,av_rho_dust,av_rho_stars)) then 

     ! Taking care of neighbour cells if necessary 

    if (lvl(ncurr) > 1) then

     call subdivide_neighbour_cells()

    endif
!!$
!!$! Subdivide current cell  
     
   if (lvl(ncurr) == nlevel) then  
      call subdivide_cell(ncurr,nlevel) 
   endif
!!$
  endif 

 end do subdiv

! update current level
      nlevel=nlevel+1
      nstart(nlevel+1)=tot_ncell+1
 
      if (nlevel < max_lvl) then 
      print *,'level up to',nlevel
      print *,(nstart(i),i=0,nlevel)
      
        cycle
      else          
         exit
      endif 

   end do

! ---------------------------------------------
! Check if there are subdivided neighbours in the list and check if they have neighbours which need to be subdivided
! -------------------------------------------------

   a=nstart(nlevel)
   b=nstart(nlevel+1)

print *, 'A   B'
print *, a,b
!read(*,*)

do

   do ncurr=a,b-1  
      if (lvl(ncurr) < nlevel) then 
         if (lvl(ncurr) > 1) then

            call subdivide_neighbour_cells
               
         endif
      endif
   end do

   a=b
   b=tot_ncell+1
   print *, 'A   B'
   print *, a,b
   ! read(*,*)

   if ((b-a) > 1) then      
      cycle
   else          
      exit
   endif

end do

tot_ncell=tot_ncell+1  ! this is important for io_routines

print *, 'grid subdivision completed'
   

 end subroutine subdivision_loop



subroutine subdivide_neighbour_cells
! this routine subdivide neighbour cells to cell ncurr if necessary (that is if the neighbours of ncurr do not have same nesting level of ncurr) 
  integer :: isel, inc, out, flag_jump
  integer(kind=int32) :: cc,nlevel_cc  
  
 dir_loop: do isel=1,3
  
 inc_loop: do inc=+1,-1,-2    
    
   call find_neighbours(isel,inc,cc,out)

   if (out /= 1) then

      call check_level_jump(cc,flag_jump)
      if (flag_jump == 1) then 
                 
         nlevel_cc=lvl(cc)

         call subdivide_cell(cc,nlevel_cc) 
 
      endif

   endif

   end do inc_loop
   end do dir_loop


  end subroutine subdivide_neighbour_cells

 subroutine subdivide_cell(cc,nlevel_cc) 
!This subroutine splits a cell in many child cells and creates appropriate links 
IMPLICIT NONE 

integer, pointer  :: ix,iy,iz
integer, target :: ii(3)
integer :: i,j,ib
integer(kind=int32) :: cc,nlevel_cc 
real(Kind=real64) :: cellsize,x,y,z,av_rho_dust,av_rho_stars

ix => ii(1)
iy => ii(2)
iz => ii(3)

   cchild(cc)=tot_ncell+1

   ib=1
   if (nlevel_cc > 0 ) ib=2

if (tot_ncell+base(ib)**3 > max_ncell) then
   print *, 'too many cells! Raise max_ncell and try again' 
   stop
endif

   
   do iz=0,base(ib)-1
   do iy=0,base(ib)-1
   do ix=0,base(ib)-1

     tot_ncell=tot_ncell+1
     ncell(tot_ncell)=tot_ncell
    
     !print *, cindex(cc)

     if (nlevel_cc >0) then
   cindex(tot_ncell)=ior(cindex(cc),((iz*base(ib)+iy)*base(ib)+ix+1)*basediv(1)*basediv(2)**(nlevel_cc-1))
   else 
     cindex(tot_ncell)=ior(cindex(cc),((iz*base(ib)+iy)*base(ib)+ix+1))
     endif

   if (cindex(tot_ncell) < 0) then 
      print *, 'CINDEX problem'
      print *, '< 0'
      print *, cindex(cc)
      stop
      endif 

      call calc_cellsize(cellsize, nlevel_cc+1)
    !cellsize=modelsize*dble(baseinv(1)*baseinv(2)**(nlevel_cc))
     
    csize(tot_ncell)=cellsize
    lvl(tot_ncell)=(nlevel_cc+1)
    

    do j=1,3
       if (mod(base(ib),2) == 1) then 
    ccoord(j,tot_ncell)=ccoord(j,cc)+cellsize*dble(ii(j)-base(ib)/2)
        
       elseif (mod(base(ib),2) == 0) then
    ccoord(j,tot_ncell)=ccoord(j,cc)+cellsize*dble((ii(j)-base(ib)/2+0.5))
         
       else
          STOP 'something wrong here: subdivide routine'
       endif
          
 end do

      x=ccoord(1,tot_ncell)
      y=ccoord(2,tot_ncell)
      z=ccoord(3,tot_ncell)

      cchild(tot_ncell)=-1  !!! this has to be before calc_dens for pcell to be properly assigned in av_galaxy !!! 

      call calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_stars,tot_ncell)
       
      dens(tot_ncell)=av_rho_dust    
      dens_stars(tot_ncell)=av_rho_stars      

print *, tot_ncell, (ccoord(j,tot_ncell),j=1,3)!, dens(tot_ncell),dens_stars(tot_ncell)



end do
end do
end do

end subroutine subdivide_cell  
  



END PROGRAM create_adap_grid_Nbody_SPH
