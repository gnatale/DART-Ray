!> Contains all the subroutines to handle the ray lists. A ray goes into a ray list if it is blocked because of the solid angle requirements for the ray beam not being fulfilled but still having substantial luminosity to be processed. There is a "high" list for rays to be split and a "low" list for rays that can potentially be joined together. 
MODULE ray_list
  use iso_fortran_env
  use healpix_routines
  use smooth_grid_routines
  IMPLICIT NONE 
  type list_type1
     !> @param theta Theta angle of the ray direction
     real(kind=real64) :: theta
     !> @param phi Phi angle of the ray direction
     real(kind=real64) :: phi
     !> @param src_lum Ray luminosity
     real(kind=real64), allocatable :: src_lum(:)
     !> @param nside HEALPix nside parameter associated with the ray
     integer :: nside
     !> @param dplane Distance from the ray origin to the last intersection plane along the perpendicular to that plane. 
     real(kind=real64) :: dplane
     !> @param prev Path crossed by the ray
     real(kind=real64) :: prev
     !> @param isel Index specifying the perpendicular to the last intersection plane (0=x, 1=y, 2=z)
     integer :: isel
     !> @param ray_type Integer specifying type of ray (e.g. 'ray_type_high', 'ray_type_low', 'ray_type_reco', 'ray_type_gone')
     integer :: ray_type
     !> @param cc_old ID number of the last intersected cell 
     integer :: cc_old
   
  end type list_type1
  !> @param ray_high_list High ray list. These rays are going to be split.
  type (list_type1), allocatable :: ray_high_list(:)
  !> @param ext_ray_list Extracted ray list. These rays are going to be processed in the next ray iteration.
  type (list_type1), allocatable :: ext_ray_list(:)
  !> @param ray_low_list Low ray list. These rays are going to be merged if possible.
  type (list_type1), allocatable :: ray_low_list(:)
  !> @param ihigh Counter of rays in the ray_high_list()
  integer :: ihigh
  !> @param ilow Counter of rays in the ray_low_list()
  integer :: ilow
  !> @param dim_high Dimension of ray_high_list(). The list can be expanded if necessary.
  integer :: dim_high
  !> @param di_low Dimension of ray_low_list(). The list can be expanded if necessary. 
  integer :: dim_low
  !> @param ind_ext Array containing the indeces of the rays to be extracted
  integer, allocatable :: ind_ext(:)
  !> @param inum Initial dimension of the ray lists.
  integer, parameter :: inum=2000
  !> @param no_ray Logical equal to TRUE if there are no rays to be extracted from the lists. 
  logical :: no_ray
  !> @param ray_type_high Ray type identifier for "high" rays. These rays are going to be split. 
  integer, parameter :: ray_type_high = 0 
  !> @param ray_type_low Ray type identifier for "low" rays. These rays are going to be merged if possible.  
  integer, parameter :: ray_type_low = 1 
  !> @param ray_type_reco Ray type identifier for "reco" rays. These rays could not be merged and their propagation restarted from where they were blocked.
  integer, parameter :: ray_type_reco = 2 
  !> @param ray_type_gone Ray type identifier for "gone" rays. These rays can be eliminated.  
  integer, parameter :: ray_type_gone = 3 
  !> @param ray_type_i_obs Ray type identifier for "i_obs" rays. These rays are used when calculating the values for the i_obs() array.
  integer, parameter :: ray_type_i_obs = 4 
  !> @param ray_type_i_obs_in Ray type identifier for "i_obs_in" rays. These rays are used when calculating the values for the i_obs_in() array.
  integer, parameter :: ray_type_i_obs_in = 5  

  !$OMP THREADPRIVATE(ray_high_list, ihigh,ind_ext,ext_ray_list,ray_low_list,ilow,dim_high,dim_low, no_ray)

CONTAINS  

  !> This routine creates the first ray to be processed. It is set up as an high ray type.
  !> @param [in] idir First index used for determining ray number on the HEALPix sphere. Values in the range [0, npix_main-1]. HEALPix pixel number in NESTED scheme equal to idir*nside**2+idir2*nside**2/4. 
  !> @param [in] idir2 Second index used for determining ray number on the HEALPix sphere. Values in the range [0, 3]. HEALPix pixel number in NESTED scheme equal to idir*nside**2+idir2*nside**2/4.
  !> @param [in] src_lum Source luminosity.
  !> @param [in] src_id Source ID number.
  !> @param [in] nside_low HEALPix NSIDE parameter for the resolution level immediately lower to the ray that has to be created. In this case alway equal to nside_min()/2. 
  subroutine create_high_ray_list(idir,idir2,src_lum,src_id,nside_low)
    
    integer :: idir,idir2, nside_low,isel_old,ipix,src_id,iw
    REAL(KIND=REAL64) :: src_lum(0:lnum-1),theta, phi, prev
    logical :: arr_exist 
  
    ! check if high res list already exist 
    arr_exist=allocated(ray_high_list)
    if (arr_exist) STOP 'ray_high_list should not be already allocated' 
    ihigh=0
    allocate(ray_high_list(0:inum-1))
    dim_high=inum
    
    ! ray direction 
    ipix=idir*nside_low**2+idir2*nside_low**2/4. 

    call pix2ang_nest(nside_low, ipix, theta, phi)  

    prev=0
    isel_old=0

    ! assign value to ray_high_list
    ray_high_list(ihigh)%theta=theta  
    ray_high_list(ihigh)%phi=phi
    ray_high_list(ihigh)%nside=nside_low
    ray_high_list(ihigh)%dplane=prev    !!!   these are just initialized 
    ray_high_list(ihigh)%prev=prev
    ray_high_list(ihigh)%isel=isel_old  !!!
    allocate(ray_high_list(ihigh)%src_lum(0:lnum-1))
    ray_high_list(ihigh)%src_lum=src_lum
    ray_high_list(ihigh)%ray_type=ray_type_high
    ray_high_list(ihigh)%cc_old=src_id
  
    ihigh=ihigh+1  ! counter of rays in the high ray list 

  end subroutine create_high_ray_list

  !> This routine extracts a ray_list with the input NSIDE() value for the next ray propagation loop.
  !> @param [in] nside HEALPix NSIDE parameter for the next ray loop in main_dir_loop().
  subroutine extract_ray_list(nside)
   
    integer :: i ,j,nside,isel_ray,ipix,ipix_add,isel_ray_add,icount,cc_old,cc_old_add, ni,ic_arr(4),iw,ic
    logical :: arr_exist
    REAL(KIND=REAL64) :: theta, phi, src_lum_ray(0:lnum-1), dplane_ray,dplane_ray_add,prev_ray
    integer :: ray_type
    integer, allocatable :: ipix_arr(:)

    no_ray = .FALSE.

    if (ihigh > 0) then   !!! extract rays from high list 
       ! find elements to be extracted 
       ni=0

       do i=0, ihigh-1 

          if (ray_high_list(i)%nside == nside/2) then ! note: this is divided by 2 because the next rays will have 2*nside(parent ray)
             ni=ni+1
          endif

       end do

       if (ni == 0) then
          no_ray =.TRUE.
          return
       endif

       allocate(ind_ext(0:ni-1))

       j=0
       do i=0,ihigh-1 

          if (ray_high_list(i)%nside == nside/2) then ! note: this is divided by 2

             ind_ext(j)=i  !!! this is the list of rays to be extracted !!!
             j=j+1
          endif
       end do


! check if extracted  high res list already exist 
       arr_exist=allocated(ext_ray_list)
       if (arr_exist) STOP 'ray_list should not be already allocated' 

       allocate(ext_ray_list(0:ni-1))

       do i=0,ni-1
          allocate(ext_ray_list(i)%src_lum(0:lnum-1))
          ext_ray_list(i)=ray_high_list(ind_ext(i)) ! add element to extracted list
          
       end do

! remove element from the high list 

       call compress_high_list()       
       deallocate(ind_ext)

    else if ((ihigh == 0).and.(ilow >0)) then   !!! extract rays from low list
       ! find elements to be extracted from low list

       ni = 0 
       do i=0, ilow-1 
      
          if (ray_low_list(i)%nside == nside*2) then  ! note: this is multiplied by 2 because the nside value of the rays in the next loop is nside(child)/2
             
             ni=ni+1
          endif
       end do

       if (ni == 0) then
          no_ray = .TRUE.
          return
       endif

       allocate(ind_ext(0:ni-1))
       
        j=0
       do i=0,ilow-1 

          if (ray_low_list(i)%nside == nside*2) then ! note: this is multiplied by 2
             ind_ext(j)=i  !!! this is the list of rays to be extracted !!!
             j=j+1
             
          endif
       end do

   !write(35,*) 'list low list'
   !write(35,*) ray_low_list(0:2)%dplane
  
       ! check if extracted  ray list already exist 
       arr_exist=allocated(ext_ray_list)
       if (arr_exist) STOP 'ray_list should not be already allocated' 
   
       allocate(ext_ray_list(0:ni-1))

       do i=0,size(ind_ext)-1
          allocate(ext_ray_list(i)%src_lum(0:lnum-1))
          ext_ray_list(i)=ray_low_list(ind_ext(i)) ! add element to extracted list 
       end do
     
! this part of the routine looks for the extracted rays with last intersection plane in common and falling in the same higher beam. Then assign ray_type "gone" to all these rays except one which gets the average luminosity and prev path parameter as well as the theta and phi angle of the parent ray 
       allocate(ipix_arr(0:ni-1))

       do i=0,ni-1 ! pre calculation ipix for inner loop later 
          theta=ext_ray_list(i)%theta 
          phi=ext_ray_list(i)%phi
          call ang2pix_nest(nside, theta, phi,ipix)
          ipix_arr(i)=ipix
       enddo

       
       do i=0,ni-1

          ic_arr=-1 ! initialize index array for group of merging rays 
          
          if (ext_ray_list(i)%ray_type /= ray_type_gone) then
             dplane_ray=ext_ray_list(i)%dplane             
             isel_ray=ext_ray_list(i)%isel 
             cc_old=ext_ray_list(i)%cc_old
             
             ipix=ipix_arr(i)
            
             call pix2ang_nest(nside, ipix, theta, phi) ! the output theta and phi here are different than before                   
             
             icount=1 ! numbers of merged rays
             ic_arr(icount)=i
             
             do j=i+1,ni-1
               
                dplane_ray_add=ext_ray_list(j)%dplane
                isel_ray_add=ext_ray_list(j)%isel
                cc_old_add=ext_ray_list(j)%cc_old
                
                ipix_add= ipix_arr(j)

                if ((ipix_add == ipix).and.(isel_ray == isel_ray_add).and.((abs(dplane_ray- &
                dplane_ray_add)/abs(dplane_ray) < 0.001))) then
                   icount=icount+1
                   ic_arr(icount)=j
                endif

                if (icount == 4) exit 

             enddo
             
             if (icount == 4) then ! merge rays 
                
                ext_ray_list(i)%nside=nside 
                ext_ray_list(i)%theta=theta    
                ext_ray_list(i)%phi=phi

                do iw = 0, lnum-1
                   do ic = 2, 4
                      ext_ray_list(i)%src_lum(iw) = ext_ray_list(i)%src_lum(iw) + ext_ray_list(ic_arr(ic))%src_lum(iw)
                   enddo
                   ext_ray_list(i)%src_lum(iw) = ext_ray_list(i)%src_lum(iw)/4
                 enddo
                ext_ray_list(i)%prev= sum(ext_ray_list(ic_arr)%prev)/4

                ext_ray_list(ic_arr(2:4))%ray_type = ray_type_gone
                
             else
                
               ! print *, 'here',icount
               ! read(*,*)

                do j=1,icount ! rays cannot be merged, put back into high ray list
                   theta=ext_ray_list(ic_arr(j))%theta 
                   phi=ext_ray_list(ic_arr(j))%phi
                   !nside=ext_ray_list(ic_arr(j))%nside
                   prev_ray=ext_ray_list(ic_arr(j))%prev
                   isel_ray=ext_ray_list(ic_arr(j))%isel
                   dplane_ray=ext_ray_list(ic_arr(j))%dplane
                   src_lum_ray=ext_ray_list(ic_arr(j))%src_lum
                   cc_old=ext_ray_list(ic_arr(j))%cc_old
                   ray_type=ray_type_reco

                   call add_to_high_ray_list(theta,phi,nside,prev_ray,isel_ray,dplane_ray,src_lum_ray,cc_old,ray_type)

                   ext_ray_list(ic_arr(j))%ray_type = ray_type_gone

                end do
               
             endif
             
          endif
       end do

       call compress_low_list()
       deallocate(ind_ext, ipix_arr)

    else 
       STOP ' you should never get here: EXTRACT RAY LIST'
    endif


end subroutine extract_ray_list

!> This routines expands ray_list arrays when there are too many rays to be added to the list.
!> @param [in] label ray_type_high() or ray_type_low() depending on which ray list has to be expanded.
subroutine expand_ray_list(label)

IMPLICIT NONE
type (list_type1), allocatable ::tmp_arr(:)
integer :: label
integer :: num, dim_arr   ! number of elements in the list, new array dimension
integer :: i
real (kind=real64) :: ex_par !  expansion coefficient

ex_par=1.2 ! expansion coefficient 

if (label == ray_type_high) then 
   num=size(ray_high_list)
   dim_high=int(ex_par*num)
   allocate(tmp_arr(0:num-1))
   do i = 0, num-1 
      allocate(tmp_arr(i)%src_lum(0:lnum-1))
   enddo
   tmp_arr(0:num-1)=ray_high_list
   deallocate(ray_high_list)
   allocate(ray_high_list(0:dim_high-1))
   do i = 0, num-1 
      allocate(ray_high_list(i)%src_lum(0:lnum-1))
   enddo
   ray_high_list(0:num-1)=tmp_arr

elseif (label == ray_type_low) then 
   num=size(ray_low_list)
   dim_low=int(ex_par*num)

   allocate(tmp_arr(0:num-1))
   do i = 0, num-1 
      allocate(tmp_arr(i)%src_lum(0:lnum-1))
   enddo
   tmp_arr(0:num-1)=ray_low_list
   deallocate(ray_low_list)
   allocate(ray_low_list(0:dim_low-1))
   do i = 0, num-1 
      allocate(ray_low_list(i)%src_lum(0:lnum-1))
   enddo
   ray_low_list(0:num-1)=tmp_arr

else 
print *, 'type list not recognized!!! Something wrong here'
stop
endif 

end subroutine expand_ray_list


!>  This routine removes the elements from high_list that will be processed later   
subroutine compress_high_list()

  type (list_type1), allocatable :: tmp_arr(:)
  integer :: iarr  ! size new high_list array
  integer :: i,it
  integer, allocatable :: keep_list(:)

  iarr=ihigh-size(ind_ext)
  
  if (iarr > 0 ) then  
     
     allocate(tmp_arr(0:iarr-1), keep_list(0:ihigh-1))

     keep_list = 1
     keep_list(ind_ext) = 0
     
     it=0
     do i=0, ihigh-1
   
        !   if (minval(abs(ind_ext-i)) >= 1) then  ! this selects the elements that
        if (keep_list(i) == 1) then
           allocate(tmp_arr(it)%src_lum(0:lnum-1))
           tmp_arr(it)=ray_high_list(i)        ! have to stay in the high list 
           it=it+1
        endif
     end do

     deallocate(ray_high_list, keep_list)
     allocate(ray_high_list(0:dim_high-1))
     do i = 0, iarr-1
        allocate(ray_high_list(i)%src_lum(0:lnum-1))
     end do
     ray_high_list(0:iarr-1)=tmp_arr
     ihigh=iarr

     deallocate(tmp_arr)
     
  else 
     deallocate(ray_high_list)
     ihigh=0
  endif

end subroutine compress_high_list

!> This routine removes elements from low_list that will be processed later.  
subroutine compress_low_list()

  type (list_type1), allocatable :: tmp_arr(:)
  integer :: iarr  ! size new low_list array
  integer :: i,it,nind
  integer, allocatable :: keep_list(:)
  
  nind = size(ind_ext)
  
  iarr=ilow-nind 
  
  if (iarr > 0 ) then

     allocate(tmp_arr(0:iarr-1), keep_list(0:ilow-1))

     keep_list = 1
     keep_list(ind_ext) = 0

     it=0
     do i=0, ilow-1  
        ! if (minval(abs(ind_ext-i)) >= 1) then
        if (keep_list(i) == 1) then
           allocate(tmp_arr(it)%src_lum(0:lnum-1))
           tmp_arr(it)=ray_low_list(i)
           it=it+1        
        endif
     end do
     
     deallocate(ray_low_list, keep_list)
     allocate(ray_low_list(0:dim_low-1))
     do i = 0, iarr-1
        allocate(ray_low_list(i)%src_lum(0:lnum-1))
     end do
     ray_low_list(0:iarr-1)=tmp_arr
     ilow=iarr

     deallocate(tmp_arr)

  else 
     !print *, ray_low_list(0:3*ilow-1)%nside
     deallocate(ray_low_list)
     ilow=0
  endif


end subroutine compress_low_list

!> This subroutine add the blocked ray to the high ray list.
!> @param [in] theta Theta angle of the blocked ray direction
!> @param [in] phi Phi angle of the blocked ray direction
!> @param [in] Nside HEALPix nside parameter associated with the blocked ray
!> @param [in] prev Path crossed by the ray until the last intersection plane
!> @param [in] isel_old Index specifying the direction of the perpendicular to the last intersection plane (0=x, 1=y, 2=z).
!> @param [in] dplane Distance from the ray origin to the last intersection plane along the perpendicular to that plane.
!> @param [in] src_lum Source luminosity associated with the ray
!> @param [in] cc ID number of the last intersected cell
!> @param [in[ ray_type Ray type of the blocked ray
subroutine add_to_high_ray_list(theta,phi,nside,prev,isel_old,dplane,src_lum,cc,ray_type)
IMPLICIT NONE
REAL(KIND=REAL64) :: theta, phi, prev, dplane,src_lum(0:lnum-1)
integer :: nside,isel_old,cc
!integer, parameter :: inum=10000
logical :: arr_exist 
integer :: label
integer :: ray_type

! check if high res list already exist 
arr_exist=allocated(ray_high_list)

if (.not.(arr_exist)) then  ! if not, create ! 

ihigh=0
allocate(ray_high_list(0:inum-1))
dim_high=inum

endif 

if (ihigh > dim_high-1) then
label=ray_type_high 
!print *, 'ihigh too high'
!read(*,*)
call expand_ray_list(label)
!stop
endif 


! assign value 
ray_high_list(ihigh)%theta=theta   
ray_high_list(ihigh)%phi=phi
ray_high_list(ihigh)%nside=nside
if (isel_old /= 0) then ! isel=0 happens if ray density increase at emitting cell
   !ray_high_list(ihigh)%dplane=prev/cproj(isel_old)
   ray_high_list(ihigh)%dplane=dplane
   ray_high_list(ihigh)%prev=prev
else
   ray_high_list(ihigh)%dplane=0.
   ray_high_list(ihigh)%prev=0.
endif
   
!!$if ( cproj(isel_old) == 0 ) then
!!$print *, 'ma guarda un po'
!!$print *, cproj
!!$print *, prev,isel_old
!!$print *, ray_high_list(ihigh)%dplane
!!$read(*,*)
!!$endif 
ray_high_list(ihigh)%isel=isel_old
allocate(ray_high_list(ihigh)%src_lum(0:lnum-1))
ray_high_list(ihigh)%src_lum=src_lum
!ray_high_list(ihigh)%ray_type='high'
ray_high_list(ihigh)%ray_type=ray_type
ray_high_list(ihigh)%cc_old=cc

ihigh=ihigh+1

end subroutine add_to_high_ray_list

!> This routine adds the blocked ray to the low ray list
!> @param theta Theta angle of the blocked ray direction
!> @param phi Phi angle of the blocked ray direction
!> @param nside HEALPix nside parameter associated with the blocked ray
!> @param prev Path length crossed by the ray
!> @param isel_old Index specifying the direction of the perpendicular to the last intersection plane (0=x, 1=y, 2=z).
!> @param cproj Projection cosines of the ray direction
!> @param src_lum Source luminosity
!> @param cc ID number of the last intersected cell 
subroutine add_to_low_ray_list(theta,phi,nside,prev,isel_old,cproj,src_lum,cc)
IMPLICIT NONE
REAL(KIND=REAL64) :: theta, phi, prev, cproj(3),src_lum(0:lnum-1)
integer :: nside,isel_old,cc
!integer, parameter :: inum=10000
logical :: arr_exist 
integer :: label

! check if low res list already exist 
arr_exist=allocated(ray_low_list)

if (.not.(arr_exist)) then  ! if not, create ! 
 
ilow=0
allocate(ray_low_list(0:inum-1))
dim_low=inum

endif 

if (ilow > dim_low-1) then
label=ray_type_low 
!print *, 'ilow too high',ilow
!read(*,*)
call expand_ray_list(label)

endif 

! assign value 
!write(35,*) 'isel old add=',isel_old

ray_low_list(ilow)%theta=theta  
ray_low_list(ilow)%phi=phi
ray_low_list(ilow)%nside=nside
ray_low_list(ilow)%dplane=prev/cproj(isel_old)
ray_low_list(ilow)%prev=prev
ray_low_list(ilow)%isel=isel_old
allocate(ray_low_list(ilow)%src_lum(0:lnum-1))
ray_low_list(ilow)%src_lum=src_lum
ray_low_list(ilow)%ray_type=ray_type_low
ray_low_list(ilow)%cc_old=cc

!write(35,*) 'new ray low dplane=',ray_low_list(ilow)%dplane 

ilow=ilow+1



end subroutine add_to_low_ray_list




END MODULE
