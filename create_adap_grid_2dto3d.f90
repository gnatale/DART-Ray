!> Creates a 3D grids from axisymmetric input grids of source volume emissivity and dust extinction coefficients. The 2D grid has to be defined on a grid of R and z points as in the example input files.   
PROGRAM create_adap_grid_2dto3d
  use iso_fortran_env
  use smooth_grid_routines
  use user_routines_2dto3d
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
  call read_input_model

  ! read lambda grid 
  call read_lambda_list
  call set_lambda_arr_si

  ! set opacity values 
  call prepare_dust_model
    
  ! check if main grid already exists. 
  
  inquire(file=trim(adjustl(dir_grid))//grid_file, exist=file_exists)

 ! IF not create it !!!

  select case (file_exists)

  case (.FALSE.)
     print *, 'Main grid does not exist. Create it !'
     ! --------------------------------------------------------------------! 
     
     lambda_in = lambda_ref_SI
     call set_2dto3d_input   
     
     ! --------------------------------------------------------------------! 
     ! create grid arrays     
     call create_grid_arrays()

     ! set subdivision factors 
     call set_base
     
     ! --------------------------------------------------------------------!
     ! start subdivision loop 
     call subdivision_loop() 
     
     ! --------------------------------------------------------------------!
     ! print 3D grid on H5DF file  
     call print_3d_grid_file()
     
     ! print info file 
     call print_info_file 
     
     print *, 'Main grid creation completed!'
     stop
case(.TRUE.)
   print *, 'Main grid does already exist. Read it !'
   call read_main_grid
        
end select

! ---------------------------------------------------------------------!
! loop on lambdas to make monochromatic emissivity/ density tables

! rename label_model

label_model= trim(adjustl(label_model_lambda_grid))

! make lambda grids
i0=minloc(abs(lambda_arr-lambda_min)/lambda_min)-1
i1=minloc(abs(lambda_arr-lambda_max)/lambda_max)-1

do il=i0(1), i1(1)
    
   lambda=lambda_arr_SI(il)  
   print *, 'lambda= ', lambda_arr(il)
   write(label_wave,'(F9.3)') lambda_arr(il) ! this converts number to string 

   ! set filename lambda grid
   grid_file_lambda='grid_'//trim(adjustl(label_model))//'_l'//trim(adjustl(label_wave))//'um.h5'
   print *, trim(adjustl(grid_file_lambda))
     

   ! set star and dust coefficients 
   lambda_in = lambda
   call set_2dto3d_input  
  
   ! create grids
   call make_lambda_grid    
     
   ! assign values to parent cells 
   call assign_dens_to_parent
     
   ! print them
   call print_lambda_grid

end do

CONTAINS

   subroutine make_lambda_grid
    ! this subroutine makes grids and the selected lambda
    real(kind=real64) :: x,y,z,cellsize,av_rho_dust, av_rho_stars
    integer :: cc,i
    

    grid_creation_lambda = .TRUE. 

    !print *, 'calculating grid for lambda =', lambda 

    call OMP_SET_NESTED(.true.)
    call omp_set_num_threads(nproc)
    
    !$OMP PARALLEL DEFAULT(NONE), PRIVATE(i,x,y,z,cellsize,av_rho_dust,av_rho_stars), &
!$OMP SHARED(tot_ncell,csize,dens,dens_stars,ccoord,cchild)

    !$OMP DO SCHEDULE(DYNAMIC,50)

    do i=0, tot_ncell-1
             
       if (cchild(i) /= -1 ) cycle 
       !print *, 'i=',i,'of ',tot_ncell-1
       x=ccoord(1,i)
       y=ccoord(2,i)
       z=ccoord(3,i) 
       cellsize=csize(i)

       
       call calc_dens(x,y,z,cellsize,av_rho_dust, av_rho_stars)  
     
       dens(i)=av_rho_dust              
       dens_stars(i)=av_rho_stars
       
    enddo

    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

end subroutine make_lambda_grid


!> Obtains the average dust extinction coefficient and the average luminosity density of the current cell.
subroutine calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_stars)
  
  real(Kind=real64) :: x,y,z,cellsize,av_rho_dust,av_rho_stars
     
  av_rho_dust= av_dens_2dto3d(x,y,z,cellsize, quantity_ext)    
  av_rho_stars = av_dens_2dto3d(x,y,z,cellsize, quantity_em)   
 
end subroutine calc_dens



!> Returns true if the cell with central coordinates x,y,z and size cellsize has to be subdivided. 
logical function subdivision(x,y,z,cellsize) 
real(kind=real64) :: x,y,z,cellsize,tau,lum_stars
integer :: ib=1 

  tau=dens(ncurr)*cellsize   ! tau=K_ext*rho_N*length  here rho_N is in H atoms per volume

  lum_stars=dens_stars(ncurr)*cellsize**3  
  
! subdivision for RT calculations

select case(subdivision_criteria)
case('standard') 
   subdivision=(((nlevel < max_lvl).and.((tau > max_dtau).or.& 
(lum_stars > max_dlum*lnu_tot) & 
.or.(nlevel < min_lvl).or.((abs(z) < z_subd_lim).and.(sqrt(x**2+y**2) < R_subd_lim)) )).or.(lvl(ncurr) < nlevel))  ! do not remove term (lvl(ncurr) < nlevel). Needed for smooth cell refinement algorithm. 

case default 
   print *, 'STOP(subdivision): subdivision_criteria not recognized!'
   print *, 'subdivision_criteria =',subdivision_criteria
   stop
end select
  
end function subdivision

  subroutine print_info_file()
! this subroutine calculates maximum and minimum tau/lum cell and print info file 
integer :: k,i,id
real(kind=real32) :: max_dtau_eff,max_dlum_eff,mean_dtau,mean_dlum  

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
open(newunit=id, file=trim(adjustl(dir_grid))//trim(adjustl(grid_info_file)), status='unknown') 

   write(id,*)
   write(id,*)  '! 3D grid parameters '  
   write(id,*)  'modelsize = ', modelsize
   write(id,*) ' base =', base  
   write(id,*) ' max_ncell', max_ncell 
   write(id,*) ' MAX DTAU PARAM (input)=', max_dtau
   write(id,*) ' MAX DLUM PARAM (input)=', max_dlum
   write(id,*) ' MIN NLVL (input)=', min_lvl
   write(id,*) ' MAX NLVL (input)=', max_lvl
   write(id,*) ' MAX DTAU (output)=', max_dtau_eff
   write(id,*) ' MEAN DTAU (output)=', mean_dtau
   write(id,*) ' MAX DLUM (output)=', max_dlum_eff
   write(id,*) ' MEAN DLUM (output)=', mean_dlum
   write(id,*)

close(id)

end subroutine print_info_file


!!!!-----------------------------------!!!!!
!!! DO NOT MODIFY FOLLOWING ROUTINES !!!! 
!!!------------------------------------!!!!   
!> Performs the subdivision loop until all cells satisfy the input subdivision criteria 
  subroutine subdivision_loop()
      
    real(kind=real64) :: x,y,z, cellsize, av_rho_dust,av_rho_stars
    real(kind=real64) :: lum_stars, tau
    real(kind=real64) :: cell_pos(3),rel_vec(3),dist
    integer :: i,outcube(3)
    integer(kind=int32) :: a,b
    

    ! FIRST CELL initialization
      x=ccoord(1,0)
      y=ccoord(2,0)
      z=ccoord(3,0)

      call calc_cellsize(cellsize,nlevel)

      
      call calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_stars)  

         dens(0) = av_rho_dust
         dens_stars(0) = av_rho_stars
        
         
         
      do 

         subdiv:  do ncurr=nstart(nlevel),nstart(nlevel+1)-1    
! check no parent cell in the list (apart from first step) 
   if ((ncurr /= 0).and.(cchild(ncurr) /= (-1))) then 
      print *, ncurr, cchild(ncurr)
      STOP 'there is something weird going on: subdivision loop ' 
      endif    ! this is needed to check problems with neighbour cell subdivision 

      x=ccoord(1,ncurr)
      y=ccoord(2,ncurr)
      z=ccoord(3,ncurr)

      call calc_cellsize(cellsize,nlevel)
      
      
 if(subdivision(x,y,z,cellsize)) then 

     ! Taking care of neighbour cells if necessary 

    if (lvl(ncurr) > 1) then

     call subdivide_neighbour_cells()

    endif
!!$
!!$! Subdivide current cell  
     
   if (lvl(ncurr) == nlevel) then  
 call subdivide_cell(ncurr, nlevel) 
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
         
         call calc_dens(x,y,z,cellsize,av_rho_dust,av_rho_stars)   
       
         dens(tot_ncell)=av_rho_dust    
         dens_stars(tot_ncell)=av_rho_stars
                  
         cchild(tot_ncell)=-1

         print *, tot_ncell, (ccoord(j,tot_ncell),j=1,3) !, dens(tot_ncell),dens_stars(tot_ncell)


end do
end do
end do

end subroutine subdivide_cell  
  



END PROGRAM create_adap_grid_2dto3d
