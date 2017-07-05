MODULE visual_routines
  use smooth_grid_routines
  use rt_routines
  use sed_routines
  IMPLICIT NONE 
  !> @param map_arr(x,y,lambda) External observer surface brightness maps at the stellar emission or dust emission wavelengths, for a single line-of-sight, which are stored in the local i_obs_arr() array.
  real(kind=real64), allocatable :: map_arr(:,:,:)
  !> @param map_in_arr(line-of-sight ID,lambda) Internal observer surface brightness maps at the stellar emission or dust emission wavelengths, for a single line-of-sight, which are stored in the local i_obs_in_arr() array.
  real(kind=real64), allocatable :: map_in_arr(:,:)
  
  !> @param npixel_maps_hd Number of pixels per side for the surface brightness maps stored in map_arr_hd() within map_projection().
  integer :: npixel_maps_hd
  
  !> @param pixel_size External observer map pixel size in the same units as csize() (only pc allowed for the moment).
  real(kind=real64) :: pixel_size
  !> @param pixel_size_hd High resolution external observer map pixel size in the same units as csize() (only pc allowed for the moment).
  real(kind=real64) :: pixel_size_hd
  !> @param area_pixel Pixel area for the observer map
  real(kind=real64) :: area_pixel
  !> @param area_pixel_hd Pixel area for the high resolution observer map
  real(kind=real64) :: area_pixel_hd
  type texture
     !> @param prof2d Projected cell normalized brightness profile. 
     real(kind=real64), allocatable :: prof2d(:,:)
  end type texture
  !> @param texture_arr Contains the normalized texture profiles for each cube subsidivion level. Needed to project a cube surface brightness profile on an arbitrary inclined map. Note that the first index of this array is 1 not zero.  
  type(texture), allocatable :: texture_arr(:)

  !> @param im0 "Start" indeces for the normalized brightness profile array. This is used so only the subsection of the prof2d() array is used when adding contributions to the surface brigthness map. Each index correspond to the different cell subdivision levels. 
  integer, allocatable :: im0(:)

  !> @param im1 "End" indeces for the normalized brightness profile array. This is used so only the subsection of the prof2d() array is used when adding contributions to the surface brigthness map. Each index correspond to the different cell subdivision levels. 
  integer, allocatable :: im1(:) 

  !> @param nt_arr Number of pixels in each cell texture profile within texture_arr.
  integer, allocatable :: nt_arr(:)

  !> @param hd_xfactor Ratio between npixel_maps_hd() and npixel_map().
  integer, parameter :: hd_xfactor = 4 

  !> @param iq_px List of subscripts corresponding to the HEALPix pixel number in the NESTED scheme. Used in the map calculations for the internal observer. 
  integer, allocatable :: iq_px(:)

  !> @param pix_list List of HEALPix pixels that are used when calculating the ray intersection lengths with the grid cells. This calculation is needed when calculating the projection of the cells on the observer maps. 
  integer, allocatable :: pix_list(:)

  !> @param neighbour_list List of HEALPix pixel neighbours to pixels for which the intersection length with a certain grid cell has already been calculated. Neighbour pixels are always considered until all their intersection lengths are zero (see calc_cube_texture_sphere() ). 
  integer, allocatable :: neighbour_list(:)

  !> @param mproj HEALPix map containing the normalized projection of the grid cells onto the observer sky.  
  real(kind=real64), allocatable :: mproj(:)

  !$OMP THREADPRIVATE (iq_px, pix_list, neighbour_list, mproj)
  
CONTAINS

!> Calls routines to calculate surface brightness maps for the external observer. In MPI mode, every node calculates the maps corresponding to the wavelengths stored in the local i_obs_arr(). Then reduction of the map_arr_out() array is performed. If no communication mode is used, for simplicity only the master node calculates the maps. 
subroutine make_maps
integer :: idir, il, k 
integer :: kall, knode
integer :: i0, i1 
real(kind=real64) :: theta, phi
real(kind=real64) :: obs_vec(0:2), xi_dir(0:2), yi_dir(0:2)

call set_i_opacity_arrays(i0,i1)

!skip if dust RT and cnflag_dust /= TRUE
if (rt_algorithm_ID == rta_dust .or. rt_algorithm_ID == rta_dust2d) then
   if (.not. cnflag_dust) return
endif

if (.not. print_maps .and. .not. print_maps_in) return 

! allocate map arrays (single direction and multiple direction) 
call create_map_arrays

if (print_maps .and. tot_ndir > 0) then 

   if (main_prc) print *, 'calculating external observer maps....'

   ! create cube texture arrays
   call create_texture_array

   ! loop on directions 
   do idir = 0, tot_ndir-1 

      if (lnum_node_maps == 0) exit   ! no locally stored i_obs for which surface brightness maps have to be calculated 
      if (no_communications .and. .not. main_prc) exit ! in no communication mode, only main_prc calculates the maps (this can be done better but it is the simplest solution for the moment)
   
      if (main_prc) print *, 'idir =', idir 
   ! define line-of-sight angles 
      theta = dir_i_out(idir,1)
      phi = dir_i_out(idir, 2)
   
      ! create texture for the input line of sight direction 
      call calc_texture_array(theta, phi, obs_vec, xi_dir, yi_dir)

      ! make projection
      call map_projection(idir, obs_vec, xi_dir, yi_dir)
   
      ! store in output multiwavelength/ multidirectional array
      call find_starting_kall(kall) ! counter for all wavelengths stored in map_arr_out          
      knode = 0 ! counter for wavelength stored in map_arr
      do il = 0, lnum-1
         if (ind_out_maps(kall) == il+i0) then
            if (iq_sca_node(il)) then             
               map_arr_out(:,:,kall,idir) = map_arr(:,:,knode)
               knode = knode + 1 
            endif
            kall = kall + 1 
            if (kall == lnum_maps) exit
         endif
      end do

   enddo

   ! reduce map_arr_out 
   call reduce_map_arr_out

   ! deallocate map_arr (not map_arr_out that has to be output)
   deallocate(map_arr)

   call print_done

endif

if (print_maps_in .and. tot_ndir_in > 0) then

   if (main_prc) print *, 'calculating internal observer maps....'

   ! loop on directions 
   do idir = 0, tot_ndir_in-1 

      if (lnum_node_maps == 0) exit   ! no locally stored i_obs for which maps have to be calculated 
      if (no_communications .and. .not. main_prc) exit ! in no communication mode, only main_prc calculates the maps (this can be done better but it is the simplest solution for the moment)
   
      if (main_prc) print *, 'idir =', idir 

      ! make projection
      call map_in_projection(idir)
   
      ! store in output multiwavelength/ multidirectional array
      call find_starting_kall(kall) ! counter for all wavelengths stored in map_in_arr_out 
      knode = 0 ! counter for wavelength stored in map_in_arr
      do il = 0, lnum-1
         if (ind_out_maps(kall) == il+i0) then
            if (iq_sca_node(il)) then            
               map_in_arr_out(:,kall,idir) = map_in_arr(:,knode)
               knode = knode + 1 
            endif
            kall = kall + 1 
            if (kall == lnum_maps) exit
         endif
      end do

   enddo

   ! reduce map_arr_out 
   call reduce_map_in_arr_out

   ! deallocate map_in_arr (not map_in_arr_out that has to be output)
   deallocate(map_in_arr)
  
   call print_done

endif

! convert maps to units of MJy/sr 
if (main_prc) then
   print *, 'convert map units to MJy/sr...'
   call convert_maps_to_MJy_sr
endif

call print_done

end subroutine make_maps

!> Covert the map units to MJy/sr. 
subroutine convert_maps_to_MJy_sr 
  integer :: il , i0, i1,k
  integer :: kall

  if ((rt_algorithm_ID == rta_projection) .and. (param_to_project == 'optical_depth')) return 

  call set_i_opacity_arrays(i0,i1)
  
select case (units_i_obs)
       
case ('erg/s/Hz/pc^2/sr')

   if (print_maps) then 
      map_arr_out=map_arr_out*10.0_real64**(-7) ! this is to convert erg -> Joule
      map_arr_out=map_arr_out*10.0_real64**(20)/parsec**2   ! MJy/sr
   endif
   if (print_maps_in) then
      map_in_arr_out=map_in_arr_out*10.0_real64**(-7) ! this is to convert erg -> Joule
      map_in_arr_out=map_in_arr_out*10.0_real64**(20)/parsec**2   ! MJy/sr
   endif
      
case ('W/Hz/pc^2/sr')

   if (print_maps) then 
      map_arr_out=map_arr_out*10.0_real64**(20)/parsec**2  !  MJy/sr
   endif
   if (print_maps_in) then 
      map_in_arr_out=map_in_arr_out*10.0_real64**(20)/parsec**2  !  MJy/sr
   endif
      
case ('W/m/pc^2/sr')

   if (print_maps) then 
      call find_starting_kall(kall) ! counter for all wavelengths stored in map_arr_out                 
      do il = 0, lnum -1
         if (ind_out_maps(kall) == il+i0) then                                   
               map_arr_out(:,:,kall,:) = map_arr_out(:,:,kall,:)*lambda_arr_SI(il+i0)**2/cspeed
               map_arr_out(:,:,kall,:) = map_arr_out(:,:,kall,:)*10.0_real64**(20)/parsec**2
            kall = kall + 1 
            if (kall == lnum_maps) exit
         endif
      enddo
   endif

   if (print_maps_in) then 
      call find_starting_kall(kall) ! counter for all wavelengths stored in map_in_arr_out     
      do il = 0, lnum -1
        if (ind_out_maps(kall) == il+i0) then
           if (iq_sca_node(il)) then
              map_in_arr_out(:,kall,:) = map_in_arr_out(:,kall,:)*lambda_arr_SI(il+i0)**2/cspeed
              map_in_arr_out(:,kall,:) = map_in_arr_out(:,kall,:)*10.0_real64**(20)/parsec**2
           endif
           kall = kall + 1 
           if (kall == lnum_maps) exit
        endif
     enddo
  endif
      
case DEFAULT

   if (main_prc) print *, 'STOP(convert_maps_to_MJy_sr): which units for i_obs ? input units_i_obs!'
   stop
   
end select


end subroutine convert_maps_to_MJy_sr




!> Creates map_arr(), map_arr_hd() and map_out_arr(). It also sets size_map() and pixel_size(), pixel_size_hd().
subroutine create_map_arrays 

  if (print_maps) then 
  
     size_map = map_size_factor*modelsize   
     
     npixel_maps_hd = hd_xfactor*npixel_maps

     if (.not. allocated(map_arr)) allocate(map_arr(0:npixel_maps-1, 0:npixel_maps-1, 0:lnum_node_maps-1))
     if (.not. allocated(map_arr_out)) allocate(map_arr_out(0:npixel_maps-1, 0:npixel_maps-1, 0:lnum_maps-1, 0:tot_ndir-1)) ! note that the wavelength dimension is lnum_maps NOT lnum_node_maps for the out array
  
     map_arr = 0
     map_arr_out = 0 
  
     pixel_size = size_map/npixel_maps
     pixel_size_hd = size_map/npixel_maps_hd
     area_pixel = pixel_size**2
     area_pixel_hd = pixel_size_hd**2

  endif

  if (print_maps_in) then

     npix_maps = 12*(2**(2*kp_maps))
     if (.not. allocated(map_in_arr)) allocate(map_in_arr(0:npix_maps-1, 0:lnum_node_maps-1))
     if (.not. allocated(map_in_arr_out)) allocate(map_in_arr_out(0:npix_maps-1, 0:lnum_maps-1, 0:tot_ndir_in-1)) ! note that the wavelength dimension is lnum_maps NOT lnum_node_maps for the out array
     map_in_arr = 0
     map_in_arr_out = 0
     
  endif
  
  
end subroutine create_map_arrays

!> Creates texture_arr().   
subroutine create_texture_array
  integer :: i, nt  
  !> @param lside cell size 
  real(kind=real64) :: lside

  if (allocated(texture_arr)) return

  allocate(texture_arr(1:max_lvl), nt_arr(1:max_lvl))

  do i = 1, max_lvl 
     
     lside=csize(0)/(base(1)*base(2)**(i-1))
     
     nt=nint(npixel_maps_hd/size_map*lside*2.5)
     allocate(texture_arr(i)%prof2d(0:nt-1,0:nt-1))
     nt_arr(i) = nt 
     
  end do
  

  allocate(im0(1:max_lvl), im1(1:max_lvl))

end subroutine create_texture_array


!> Calculates texture_arr(), the normalized surface brightness profiles for the projections of the cells on the external observer maps. It also sets the direction of the axis on the observer maps. These directions are such that, when seeing the galaxy edge-on at a position along the 3D X axis, the Y' axis on the maps corresponds to the 3D Z axis while the X' axis on the map corresponds to the 3D Y axis.  
!> /todo maybe instead of 0.99 factor just subtract half pixel size 
subroutine calc_texture_array(theta, phi,obs_vec, xi_dir, yi_dir)
  integer :: i,j, ip 
  integer :: ix, iy
  !> @param cell coordinates 
  real(kind=real64) :: rc(0:2)
  !> @param ro observer coordinates 
  real(kind=real64) :: ro(0:2)
  !> @param theta0 phi0 lines of sight direction (it should not matter if in the reference frame of the observer or of the cell to be projected because the corresponding line is the same) 
  real(kind=real64) :: theta, phi
  !> @param obs_vec versor parallel to observer line of sight 
  real(kind=real64) :: obs_vec(0:2)
  !> @param z_dir Z-axis direction 
  real(kind=real64) :: z_dir(0:2) = (/0,0,1/)
  !> @param z_zi  distance of vertex z versor - projection plane
  real(kind=real64) :: z_zi
  !> @param yi_dir first reference axis on projection plane. When Z is not parallel to the line of sight, yi_dir is the projection of Z on the observer plane. Otherwise, yi_dir is chosen such to be the limit in the case of very small differences between Z and obs_vec.
  real(kind=real64) :: yi_dir(0:2)
  !> @param xi_dir second reference axis on projection plane. Always chosen as the cross product between obs_vec and yi_dir 
  real(kind=real64) :: xi_dir(0:2)
  
  real(kind=real64) :: dtheta 

  real(kind=real64), allocatable :: xpp(:), ypp(:)
  
  real(kind=real64) :: cellsize 
  !> @param length Length of the line-of-sight - cell intersection. 
  real(kind=real64) :: length

  ! initialize dtheta
  dtheta = 1E-3

  ! initialize texture_arr
  do i = 1, max_lvl 
     texture_arr(i)%prof2d(:,:) = 0. 
  enddo
  
  rc = (/0,0,0/) ! cell coordinates 

  obs_vec=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)] ! line-of-sight 
  
  z_zi=scal_prod(obs_vec,z_dir) ! distance of vertex z versor - projection plane

  if (abs(abs(z_zi)-1.) > 1E-6) then  ! when Z-axis not parallel to line-of-sight (theta /= 0 and theta /= 180)
     
    yi_dir=z_dir-z_zi*obs_vec    ! first reference axis on projection plane   

  else 
     
     if (abs(theta-pi) < 1E-6) dtheta= -dtheta
     obs_vec=[sin(theta+dtheta)*cos(phi),sin(theta+dtheta)*sin(phi),cos(theta+dtheta)]
     
     z_zi=scal_prod(obs_vec,z_dir)
     yi_dir=z_dir-z_zi*obs_vec   ! vector on the plane slightly inclined compared to the observer plane
     z_zi=scal_prod(yi_dir,z_dir)
     yi_dir=yi_dir-z_zi*z_dir  ! vector on the observer plane (projection of the previous yi_dir on XY plane)

     obs_vec=[sin(theta-dtheta)*cos(phi),sin(theta-dtheta)*sin(phi),cos(theta-dtheta)]

  endif
  
  call normalize(yi_dir)
  !call cross_prod(obs_vec, yi_dir, xi_dir) !previous version giving opposite X axis orientation. 
  call cross_prod(yi_dir,obs_vec, xi_dir)  ! see subroutine description. This line sets the direction X' on the observer maps. 

  do i = 1, max_lvl 

     allocate(xpp(0:nt_arr(i)-1), ypp(0:nt_arr(i)-1))
     
     if (mod(nt_arr(i),2) == 0) then ! there is a difference in the way the texture have to be calculated if the texture side contains an odd or an even number of pixels. this works although not completely sure why.   

        do ip = 0, nt_arr(i)-1
           ! these are the pixel centre coordinates
           xpp(ip) = ip*pixel_size_hd -nt_arr(i)*pixel_size_hd/2!+pixel_size_hd/2
        end do

     else 

        do ip = 0, nt_arr(i)-1
           ! these are the pixel centre coordinates
           xpp(ip) = ip*pixel_size_hd -nt_arr(i)*pixel_size_hd/2+pixel_size_hd/2
        end do

     endif

     ypp = xpp
     cellsize = csize_arr(i)*0.99 ! this factor is to avoid cell border overlapping
                 

     do ix = 0, nt_arr(i) -1 

        do iy = 0, nt_arr(i) -1 

           ro= rc+obs_vec*modelsize+xpp(ix)*xi_dir+ypp(iy)*yi_dir ! pixel position on the observer map in 3D
           
        
           ! determine whether ray intersect cube and intersection points 
           call find_ray_cell_intersections(rc,cellsize,ro,obs_vec, length)

           if (length > 0) then
              
              !length=(sum((p2-p1)**2.))**0.5
              
              texture_arr(i)%prof2d(ix,iy)=length

           endif

        end do

     end do

     texture_arr(i)%prof2d = texture_arr(i)%prof2d/sum(texture_arr(i)%prof2d)

     ! determine im0(i) and im1(i)
     im0(i) = 0
     im1(i) = nt_arr(i)-1
     
     do ix = 1, nt_arr(i) -1

        if (sum(texture_arr(i)%prof2d(ix,:)) > 0 .or. sum(texture_arr(i)%prof2d(:,ix)) > 0) then
           im0(i) = ix -1
           exit
        endif

     end do

     do ix = nt_arr(i) -1, 1, -1

        if (sum(texture_arr(i)%prof2d(ix,:)) > 0 .or. sum(texture_arr(i)%prof2d(:,ix)) > 0) then
           im1(i) = ix + 1
           exit
        endif

     end do
     
     deallocate(xpp, ypp)
     
  end do


end subroutine calc_texture_array

!> Finds the points of intersections between a line of sight and a cell. Used in the routines to calculate the cell normalized surface brightness profiles (both in the case of the external and internal observer).  
subroutine find_ray_cell_intersections(rc,cellsize,ro,obs_vec,length)
  integer :: j
  !> @param p1 p2 Coordinates of line-of-sight intersection points with a cell
  real(kind=real64) :: p1(0:2), p2(0:2)
  real(kind=real64) :: rc(0:2), cellsize, ro(0:2), theta,phi, obs_vec(0:2), length
  !> @param flag_ray Number of line-of-sight intersection points with a cell found  
  integer :: flag_ray 

  flag_ray = 0
  length = 0
  
  do j=0,5           ! loop on cell faces
     call calc_ray_inters(j,rc,cellsize,ro,obs_vec, flag_ray,p1,p2)
     if (flag_ray == 2) then
        length=(sum((p2-p1)**2.))**0.5          
        exit  ! as soon as you get 2 intersection, calculate length and get out from the loop
     endif
  enddo

end subroutine find_ray_cell_intersections


! this program calculates the intersection to a cell cube plane. 
subroutine calc_ray_inters(j,rc,cellsize,ro,obs_vec,flag_ray,p1,p2)
  integer :: i, j
  real(kind=real64) :: rc(0:2), cellsize, ro(0:2), theta,phi, p1(0:2), p2(0:2), po(0:2), obs_vec(0:2), rel_vec(0:2), pcoord(0:2)
  integer :: flag_ray
  integer :: isel_ray 
  integer :: n_vec(0:2)
  real(kind=real64) :: dx, num, den, dist 
  integer :: outcube(0:2)

! j meaning
! j =0  xc-csize plane
! j =1  xc+csize plane   
! j =2  yc+csize plane
! j =3  yc-csize plane
! j =4  zc+csize plane
! j =5  zc-csize plane

  
  select case(j)
  case(0,1)
     isel_ray=0
     n_vec=[1,0,0]  ! this is the vector perpendicular to the plane to be intersected 
  case(2,3)
     isel_ray=1
     n_vec=[0,1,0]
  case(4,5)
     isel_ray=2
     n_vec=[0,0,1]
  end select

! define point po on intersection plane
  if (j == 0 .or. j == 2 .or. j == 4) then 
     dx=cellsize/2.
  else 
     dx=-cellsize/2.
  endif
  
  po=rc+n_vec*dx   ! this is a point on the plane to be intersected 

! define point on the line (pixel position in the case of external observer or observer position in the case of internal observer)
!ro=ro   ! input variable           

! calculate distance from pixel centre/observer to plane -ray intersection 
  num = sum((po-ro)*n_vec)   ! these are scalar products   
  den = sum(obs_vec*n_vec)

  if (den == 0.) then 
     ! if den =0, then ray direction and intersection plane are parallel
     ! just return here
     return
  endif

  dist=num/den  ! this is the distance between the intersection point and the centre of the cell to be projected. 

  ! calculate coordinate intersection 
  pcoord=ro+dist*obs_vec  ! this is the point of intersection between the line-of-sight and the plane parallel to the cube face 

  ! calculate vector between cell centre and intersection point
  rel_vec=pcoord-rc   ! relative distance between cell centre and intersection point 

  ! check whether the intersection point is within cube face 
  outcube(:)=0
  do i=0,2 
     if ((abs(rel_vec(i)) > (1+1E-7)*cellsize/2.) .and. (i /= isel_ray)) then 
        outcube(i)=1
     endif

  enddo

  ! if the ray does not intersect the cube face, just return. If it does, add coordinates to either p1 or p2  

  if (sum(outcube) > 0) then 
     return
  else  
     if (flag_ray == 0) then  
        p1 = pcoord
        flag_ray =1
     else if (flag_ray == 1) then 
        p2 = pcoord
        flag_ray =2
     else 
        print *,'STOP(calc_ray_inters): weird, more than 2 intersections found...!'
        stop
     endif
 
  endif

end subroutine calc_ray_inters


!> Projects the 3D grid cells onto the observer map. The calculated maps are only those corresponding to the wavelengths stored in the local i_obs_arr() array. 
!> \todo make notation for Y axis on observer plane consistent
subroutine map_projection(idir, obs_vec, xi_dir, yi_dir)
  integer :: i, idir, ihost
  real(kind=real64) :: obs_vec(0:2), xi_dir(0:2), yi_dir(0:2), r(0:2), ri(0:2)
  real(kind=real64) :: vec_n, xi_p, yi_p, c0,c1
  integer :: ix_p, iy_p, ix0, ix1, iy0, iy1, delta_p, ix, iy
  integer :: ilvl, ip,il, iw
  real(kind=real64), allocatable ::xpp(:), ypp(:)
  real(kind=real64) :: tot_flux_new(0:lnum_node_maps-1), tot_flux_old(0:lnum_node_maps-1)
  integer, parameter :: ns = 4
  real(kind=real64) :: map_section(0:ns-1, 0:ns-1), map_arr_hd_old(0:npixel_maps_hd-1, 0:npixel_maps_hd-1, 0:lnum_node_maps-1)
  real(kind=real64), allocatable :: map_arr_hd_thread(:, :, :)
  !> @param map_arr_hd(x,y,lambda) External observer surface brightness maps at the stellar emission or dust emission wavelengths, for a single line-of-sight, which are stored in the local i_obs_arr() array. The difference between this array and map_arr() is that it contains maps with much higher resolution (that is, many more pixels). These maps are the ones used in the cell projections. Then they are rebinned to smaller size when passed to map_arr
  real(kind=real64) :: map_arr_hd(0:npixel_maps_hd-1, 0:npixel_maps_hd-1, 0:lnum_node_maps-1)
  real(kind=real64) :: map_arr_temp(0:npixel_maps-1, 0:npixel_maps-1, 0:lnum_node_maps-1)
  integer :: ilevel
  logical :: zero_value
  
  ! initialize maps
  map_arr = 0 
    
  ! create x,y coordinates array for high resolution map
  allocate(xpp(0:npixel_maps_hd-1), ypp(0:npixel_maps_hd-1))
  
  do ip = 0, npixel_maps_hd-1
     ! these are the pixel left/down coordinates 
     xpp(ip) = ip*pixel_size_hd -npixel_maps_hd*pixel_size_hd/2!+pixel_size_hd/2
  end do
  ypp = xpp

  call omp_set_num_threads(nproc)  

  do ilevel = 1, max_lvl

   !  if (ilevel /= max_lvl -1) cycle
     
     map_arr_hd = 0
     map_arr_temp = 0
     
     !$OMP PARALLEL DEFAULT(NONE), &
     !$OMP PRIVATE(i,r, ix_p, iy_p, ilvl, delta_p, ix, iy, ix0, ix1, iy0, iy1, il, iw, map_arr_hd_thread, map_section, zero_value), &
     !$OMP SHARED(tot_ncell, cchild, i_obs_arr,obs_vec, xi_dir, yi_dir, xpp,ypp, im0, im1, lnum_node_maps, texture_arr, csize, area_pixel_hd, idir, ccoord, lvl, map_arr_hd, npixel_maps_hd, tot_flux_old, map_arr_hd_old, npixel_maps, map_arr,map_arr_temp, iq_maps_id, ind_out_maps, ilevel) 

     allocate(map_arr_hd_thread(0:npixel_maps_hd-1, 0:npixel_maps_hd-1, 0:lnum_node_maps-1))
     map_arr_hd_thread = 0 
  
     !$OMP DO SCHEDULE(DYNAMIC,chunk_size)
  
     do i = 0, tot_ncell -1

        if (cchild(i) /= -1 .or. sum(i_obs_arr(:,i, idir)) == 0 .or. lvl(i) /= ilevel) cycle

        !print *, i
        r = ccoord(:,i)  ! cell coordinates  
        
        call find_projected_point(r, obs_vec, xi_dir, yi_dir, xpp, ypp, ix_p, iy_p)
      
        ilvl = lvl(i)

        delta_p = im1(ilvl) - im0(ilvl)

        ix0 = ix_p - delta_p/2
        ix1 = ix_p + delta_p/2
        iy0 = iy_p - delta_p/2
        iy1 = iy_p + delta_p/2

        if (ix0 < 0 .or. ix1 > npixel_maps_hd-1 .or. iy0 < 0 .or. iy1 > npixel_maps_hd-1) cycle

        do il = 0, lnum_node_maps-1
           iw = iq_maps_id(il)
           map_arr_hd_thread(ix0:ix1, iy0:iy1, il) = map_arr_hd_thread(ix0:ix1, iy0:iy1, il) + i_obs_arr(iw, i, idir)*texture_arr(ilvl)%prof2d(im0(ilvl):im1(ilvl), im0(ilvl):im1(ilvl))*(csize(i)**2)/area_pixel_hd           
        end do
        
     end do

     !$OMP END DO NOWAIT

     !$OMP CRITICAL
     map_arr_hd = map_arr_hd + map_arr_hd_thread

     !$OMP END CRITICAL

     deallocate(map_arr_hd_thread)
  
     !$OMP BARRIER

     !$OMP MASTER
  
     ! total flux before median smoothing/ rebinning 
     do il = 0, lnum_node_maps-1
        tot_flux_old(il) = sum(map_arr_hd(:,:,il))*area_pixel_hd       
     enddo

     ! median smoothing
     map_arr_hd_old = map_arr_hd
     map_arr_hd = 0 
     
     !$OMP END MASTER
     !$OMP BARRIER

     !$OMP DO SCHEDULE(DYNAMIC,npixel_maps_hd/nproc)
     do ix = 0, npixel_maps_hd-1
        do iy = 0, npixel_maps_hd -1
       
           ! if (sum(map_arr_hd_old(ix, iy, :)) == 0) cycle ! this was used to avoid using median filtering in the map regions outside the model size. however, it can leave artifacts when the model is seen aligned along its major axis. By removing it, the median filtering works well everywhere except for the pixels immediately aside of the model extent. They should not be considered as physical!!!!
           zero_value = .FALSE.
           if (sum(map_arr_hd_old(ix, iy, :)) == 0) zero_value = .TRUE.

           ix0 = ix - ns/2 
           ix1 = ix + ns/2
           iy0 = iy - ns/2
           iy1 = iy + ns/2

           ! don't perform median smoothing close to the map edge. this is necessary because the median section function does not accept arrays with dimensions different from ns x ns. 
           if (ix0 < 0 .or. iy0 < 0 .or. ix1 > npixel_maps_hd-1 .or. iy1 > npixel_maps_hd-1) then 
              map_arr_hd(ix, iy, :) = map_arr_hd_old(ix, iy, :) 
              cycle
           endif

           do il = 0, lnum_node_maps -1 
              map_section = map_arr_hd_old(ix0:ix1, iy0:iy1, il) 
              map_arr_hd(ix, iy, il) = median_section(map_section,ns, zero_value)          
           end do
        
        end do
     end do
     !$OMP END DO NOWAIT

     !$OMP BARRIER
  
     ! rebinning

     !$OMP DO SCHEDULE(DYNAMIC,npixel_maps_hd/nproc)
     do ix = 0, npixel_maps -1

        ix0 = ix*hd_xfactor
        ix1= (ix+1)*hd_xfactor-1
        
        if (sum(map_arr_hd(ix0:ix1, :,:)) == 0) cycle 
        
        do iy = 0, npixel_maps-1
           
           iy0 = iy*hd_xfactor
           iy1= (iy+1)*hd_xfactor-1 
           
           do il = 0, lnum_node_maps-1 
              map_arr_temp(ix,iy,il) = sum(map_arr_hd(ix0:ix1, iy0:iy1,il))/hd_xfactor**2
           enddo
           
        end do

     end do
     !$OMP END DO NOWAIT

     !$OMP END PARALLEL
  
     ! total flux after smoothing and rebinning
     do il = 0, lnum_node_maps-1
        tot_flux_new(il) = sum(map_arr_temp(:,:,il))*area_pixel
     enddo

     ! renormalization
     do il = 0, lnum_node_maps-1
        if (tot_flux_new(il) > 0 ) then 
           map_arr(:,:,il) = map_arr(:,:,il) + map_arr_temp(:,:,il)*tot_flux_old(il)/tot_flux_new(il)
        endif
     enddo

  enddo
     
  deallocate(xpp, ypp)

  ! project point sources
  if (tot_p_src > 0) then

     ! create x,y map coordinates for output map
     allocate(xpp(0:npixel_maps-1), ypp(0:npixel_maps-1))
  
     do ip = 0, npixel_maps-1
        ! these are the pixel left/down coordinates 
        xpp(ip) = ip*pixel_size -npixel_maps*pixel_size/2!+pixel_size/2
     end do
     ypp = xpp   

     do i = 0, tot_p_src-1

        ihost = cell_src(i)  ! index host cell 
      
        r = ccoord_p_src(:,i)  ! cell coordinates 
        
        call find_projected_point(r, obs_vec, xi_dir, yi_dir, xpp, ypp, ix_p, iy_p)
        do il = 0, lnum_node_maps -1
           iw = iq_maps_id(il)
           map_arr(ix_p,iy_p,il)=map_arr(ix_p,iy_p,il)+i_obs_arr(iw, i+tot_ncell, idir)*(csize(ihost)**2)/area_pixel

        end do

     end do

     deallocate(xpp, ypp)
     
  end if
  
end subroutine map_projection

!> Returns the average of the values in a map section. Used for average smoothing. Note that zero elements are not considered in the average smoothing. 
real(kind=real64) function average_section(map_section, ns)
  real(kind=real64) :: map_section(0:, 0:), tot_flux
  integer :: ns
  integer :: ix, iy, itot

  tot_flux = 0
  itot = 0 
  
  do ix = 0, ns-1
     do iy = 0, ns -1
        if (map_section(ix,iy) == 0) cycle 
        tot_flux = tot_flux + map_section(ix,iy)
        itot = itot + 1
        
     enddo
  enddo

  if (itot > 0) then 
     average_section = tot_flux/itot
  else
     average_section = 0
  endif

end function average_section

!> calculates median of the input 2D array
real(kind=real64) function median_section(map_section, ns, zero_value)
  real(kind=real64) :: map_section(0:, 0:)
  integer :: ns
  integer :: ntot 
  real(kind=real64) :: lin_arr(0:ns**2-1)
  integer ::  lin_ord(0:ns**2-1)
  integer :: ix, iy, i, iel   
  logical :: zero_value

  ntot = ns**2
  
  ! convert to linear array
  do ix = 0, ns -1
     lin_arr(ix*ns:(ix+1)*ns-1) = map_section(:, ix)      
  end do
  
  ! sorting
  call quick_sort(lin_arr, lin_ord)

  do i = 0, ntot-1     ! this is to avoid considering zero elements 

     if (lin_arr(i) /= 0) exit 

  end do

  ! take median value. The condition containing zero_value is there to avoid that pixels outside the model area get non zero value because of the median smoothing. For those pixels the majority of surrounding pixels is zero and they have zero value. 
  if (i < ntot .and. .not. (zero_value .and. i > (ntot-1)/2)) then  
     iel = (ntot - i)/2 + i     
  else
     iel = 0
  endif
  
  median_section = lin_arr(iel)
  

end function median_section


subroutine find_projected_point(r, obs_vec, xi_dir, yi_dir, xpp, ypp, ix_p, iy_p)
  real(kind=real64) :: r(0:2), obs_vec(0:2), vec_n, ri(0:2), xi_p, yi_p, xi_dir(0:2), yi_dir(0:2)
  real(kind=real64) :: xpp(0:), ypp(0:), c0, c1
  integer :: ix_p, iy_p, ix1, iy1

  
  vec_n = scal_prod(obs_vec, r)  ! distance from projection plane
  
  ri = r -vec_n*obs_vec  ! projection vector in xyz grid 
     
  xi_p=scal_prod(ri, xi_dir)    ! coordinates on projection plane reference frame
  yi_p=scal_prod(ri, yi_dir)   

  call value_locate(xi_p, xpp, ix_p, ix1, .FALSE.)
  c0=abs(xpp(ix_p)-xi_p)
  c1=abs(xpp(ix1)-xi_p)
  !if (c1 < c0) ix_p=ix1

  call value_locate(yi_p, ypp, iy_p, iy1, .FALSE.)
  c0=abs(ypp(iy_p)-yi_p)
  c1=abs(ypp(iy1)-yi_p)
  !if (c1 < c0) iy_p=iy1


end subroutine find_projected_point


!> Calculates dot product between two 1-dimensional vectors.
real(kind=real64) function scal_prod(vector1, vector2)
  real(kind=real64) :: vector1(0:), vector2(0:)
  integer :: info1, info2
  integer :: i

  scal_prod = 0.0          ! dummy value
  info1 = size (vector1)
  info2 = size (vector2)

  if (info2 /= info1 ) then 
     print *, 'STOP(scal_prod): vector 1 and 2 do not have same no. of elements!!!' 
     stop
  else 
     do i = 0, info1-1 
        scal_prod = scal_prod + vector1 (i) * vector2 (i)
     enddo
  endif
  
end function scal_prod

!> Normalises input vector. 
subroutine normalize(vector)
  real(kind=real64) :: vector(0:)
  
  vector = vector / sqrt(sum(vector**2))
  
end subroutine normalize


!> Calculates cross-product between two input arrays. Note that the all arrays need to have the same size. 
subroutine cross_prod(vector1, vector2, out_arr)
  real(kind=real64) :: vector1(0:), vector2(0:),out_arr(0:) 
  integer :: info1, info2,info3

  info1 = size(vector1)
  info2 = size(vector2)
  info3 = size(out_arr)

  if (info1 /= info2 .or. info2 /= info3) then 
     print *, 'STOP(Cross_prod): Input arrays do not have the same size!'
     STOP
  else 
     out_arr (0) = vector1 (1) * vector2 (2) - vector1 (2) * vector2 (1)
     out_arr (1) = - (vector1 (0) * vector2 (2) - vector1 (2) * vector2 (0))
     out_arr (2) = vector1 (0) * vector2 (1) - vector1 (1) * vector2 (0)
  endif

end subroutine cross_prod

!> Calculates surface brightness maps for internal observer by projecting emission of from each cell onto HealPix sphere centred at the observer position. 
subroutine map_in_projection(idir)
  integer :: i, idir, il,iw
  !real(kind=real64), allocatable :: mproj(:)
  real(kind=real64), allocatable :: map_in_arr_temp(:,:)
  real(kind=real64) :: rc(0:2), ro(0:2), rel_vec(0:2)
  real(kind=real64) :: omega_hp, dist2, cellsize, omega_cell, omega_ratio
  real(kind=real64) :: theta, phi
  integer :: outcube(0:2), nside_map, ip
  !integer, allocatable :: iq_px(:)
  integer ::  n_px
  integer :: ihost
  
  map_in_arr = 0 
  
  ro = ccoord_obs(:, idir)

  nside_map = 2**kp_maps
  omega_hp = 4*pi/npix_maps

  call omp_set_num_threads(nproc) 

  !$OMP PARALLEL DEFAULT(NONE), &
  !$OMP PRIVATE(i, il, iw,map_in_arr_temp, rc, rel_vec, dist2, cellsize, omega_cell, omega_ratio, theta, phi, ip, outcube, n_px ), &
  !$OMP SHARED(tot_ncell, npix_maps, lnum_node_maps, cchild, i_obs_in_arr, idir, ccoord, csize, omega_hp,  nside_map, map_in_arr, ro, id_mpi, iq_maps_id, ind_out_maps) 
  
  allocate(mproj(0:npix_maps-1), map_in_arr_temp(0:npix_maps-1, 0:lnum_node_maps-1))
  allocate(iq_px(0:npix_maps-1))
  allocate(pix_list(0:npix_maps-1), neighbour_list(0:npix_maps-1))
  mproj = 0
  map_in_arr_temp = 0
  iq_px = 0
  pix_list = 0
  neighbour_list = 0 

  !$OMP DO SCHEDULE(DYNAMIC,chunk_size)
  
  do i = 0, tot_ncell-1

     if (cchild(i) /= -1 .or. sum(i_obs_in_arr(:,i, idir)) == 0) cycle

    ! print *, i,id_mpi
     
     rc = ccoord(:,i)  ! cell coordinates  

     rel_vec = rc - ro ! vector from observer to cell centre

     dist2 = sum(rel_vec**2)  ! distance squared 

     cellsize = csize(i)  

     if (dist2 /= 0) then 
        omega_cell = cellsize**2/dist2  ! solid angle subtending the cell to be projected
     else
        omega_cell = 4*pi
     endif

     omega_ratio = omega_cell/omega_hp  ! if omega_ratio < 1 then the cell solid angle is smaller than a spherical pixel solid angle. In this case there is only one pixel over which the cell brightness is projected.

     ! three cases to be considered
     ! 1) the cell solid angle is smaller than a pixel cell -> determine the only pixel containing the entire flux coming from that cell
     ! 2) the cell contains the observer --> use uniform surface brightness approximation
     ! 3) all other cases --> calculate surface brightness profile and add to total map
     
     ! case (1)
     if (omega_ratio < 1) then

        ! find out projection pixel
        call find_theta_phi_obs_in(rc,ro, theta,phi) ! see comments for same call below
!!$        theta = pi - theta  ! calculate opposite angle because find_theta_phi finds the angles for the vector from the cell to the observer. 
!!$        phi = phi + pi
!!$        if (phi > 2*pi) phi = phi - 2*pic

        call ang2pix_nest(nside_map, theta, phi,ip)

        do il = 0, lnum_node_maps-1
           iw = iq_maps_id(il)
           map_in_arr_temp(ip, il) =  map_in_arr_temp(ip,il) + i_obs_in_arr(iw,i,idir)*omega_ratio
        enddo
        
        cycle 
     endif

     call calc_outcube(rel_vec, cellsize, outcube)

     ! case(2)
     if (sum(outcube) == 0) then

        ! use uniform surface brightness approximation
        do il = 0, lnum_node_maps-1
           iw = iq_maps_id(il)
           map_in_arr_temp(:, il) =  map_in_arr_temp(:,il) + i_obs_in_arr(iw,i,idir)/2. 
        enddo
        
        cycle
     endif

     ! case(3)
     call calc_cube_texture_sphere(cellsize,rc,ro,omega_cell,nside_map, n_px)
     
     do il = 0, lnum_node_maps-1 
        iw = iq_maps_id(il) 
        map_in_arr_temp(iq_px(0:n_px-1), il) =  map_in_arr_temp(iq_px(0:n_px-1),il) + i_obs_in_arr(iw,i,idir)*mproj(iq_px(0:n_px-1))*omega_ratio 
     enddo

     mproj(iq_px(0:n_px-1)) = 0
     
     
  end do

  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  map_in_arr = map_in_arr + map_in_arr_temp

  !$OMP END CRITICAL

  deallocate(mproj, map_in_arr_temp, iq_px, pix_list, neighbour_list)
  
  !$OMP END PARALLEL

  ! project point sources 
  if (tot_p_src > 0) then

     do i = 0, tot_p_src-1 
     
        ihost = cell_src(i)  ! index host cell

        rc = ccoord_p_src(:,i)  ! point source position
        
        rel_vec = rc - ro ! vector from observer to cell centre

        dist2 = sum(rel_vec**2)  ! distance squared 

        cellsize = csize(ihost)  

        if (dist2 /= 0) then 
           omega_cell = cellsize**2/dist2  ! solid angle subtending the cell to be projected. Remember that also point sources the i_obs values are calculated assuming the source area coincides with the cell area. 
        else
           !omega_cell = 4*pi
           cycle   ! if point source and observer coincide, projection is not performed
        endif

        omega_ratio = omega_cell/omega_hp  
        
        call find_theta_phi_obs_in(rc,ro, theta,phi)

        call ang2pix_nest(nside_map, theta, phi,ip)
        
        do il = 0, lnum_node_maps-1
           iw = iq_maps_id(il)
           map_in_arr(ip, il) =  map_in_arr(ip,il) + i_obs_in_arr(iw,i+tot_ncell,idir)*omega_ratio  
                      
        enddo

     enddo

  end if
  

end subroutine map_in_projection


!> Calculates normalized surface brightness profile of a projected cell over the observer sphere. This is done by calculating the intersection lengths of the line-of-sights corresponding to the HEALPix spherical pixels (and the observer position). The code determines the spherical region over which to calculate these intersections (calculating them over the entire sphere is inefficient and not necessary because cells cover only a limited region of the sky). This is done by calculating the intersections for neighbour pixels until all the neighbour pixels have zero intersection length. 
subroutine calc_cube_texture_sphere(cellsize,rc,ro,omega_cell,nside_map, n_px)
  integer :: idir,i, ipix, i0
  real(kind=real64) :: cellsize, rc(0:2), ro(0:2), omega_cell, obs_vec(0:2)
  integer :: nside_map
 ! real(kind=real64) :: mproj(0:npix_maps-1)
 ! integer :: iq_px(0:npix_maps-1)
  integer :: n_px
 ! integer :: pix_list(0:npix_maps-1)
  integer :: n_pix_list
 ! integer :: neighbour_list(0:npix_maps-1):q
  integer :: n_neigh_list
  real(kind=real64) :: theta, phi, length
  integer, parameter :: nneigh_max = 8 ! maximum number of pixel neighbour
  integer :: neigh_arr(0:nneigh_max), nneigh
  integer :: ip
  integer :: flag_ray
  
  ! initialize variables
 ! mproj = 0
 ! iq_px = 0
  n_px = 0 
  !pix_list = 0
  n_pix_list = 0
  !neighbour_list = 0
  n_neigh_list = 0 

 
  
  ! find starting pixel and its intersection length
  call find_theta_phi_obs_in(rc,ro, theta,phi) ! note that here I inverted rc and ro with respect to the call in rt_loop_iobs. This is because the angle is in the reference frame of the observer in this case. 
!!$  theta = pi - theta  ! calculate opposite angle because find_theta_phi finds the angles for the vector from the cell to the observer. 
!!$  phi = phi + pi
!!$  if (phi > 2*pi) phi = phi - 2*pi
  call ang2pix_nest(nside_map, theta, phi,ip)

  obs_vec=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)] 
  
  call find_ray_cell_intersections(rc,cellsize,ro,obs_vec, length)

  n_px = 1
  iq_px(0) = ip
 
  if (length > 0) then
     mproj(ip) = length
  else 
     mproj(ip) = 1   ! if no intersection despite omega_ratio > 1, then just return central pixel     
     return 
  endif

  ! add first neighbour pixels to neighbour list 
  call neighbours_nest(nside_map, ip, neigh_arr, nneigh)
  neighbour_list(0:nneigh-1) = neigh_arr(0:nneigh-1)
  n_neigh_list = nneigh

  
  
  ! Find intersections
  do

     if (n_neigh_list == 0) exit
     
     pix_list(0:n_pix_list-1) = 0
     pix_list(0:n_neigh_list-1) = neighbour_list(0:n_neigh_list-1)
     n_pix_list = n_neigh_list
     n_neigh_list = 0

     do i = 0, n_pix_list -1

        ipix = pix_list(i)
        if (mproj(ipix) > 0) cycle ! intersection already calculated
        
        call pix2ang_nest(nside_map, ipix, theta, phi)   

        flag_ray = 0 
        obs_vec=[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)] 
  
        call find_ray_cell_intersections(rc,cellsize,ro,obs_vec, length)

        if (length > 0 ) then

           mproj(ipix) = length
           n_px = n_px +1 
           iq_px(n_px-1) = ipix 
           
           call neighbours_nest(nside_map, ipix, neigh_arr, nneigh)
           i0 = n_neigh_list
           neighbour_list(0+i0:nneigh-1+i0) = neigh_arr(0:nneigh-1)
           n_neigh_list = n_neigh_list + nneigh

        endif 

     end do

  end do

  ! sort n_px

  ! normalize mproj
  mproj(iq_px(0:n_px-1)) = mproj(iq_px(0:n_px-1))/sum(mproj(iq_px(0:n_px-1)))
  
  

end subroutine calc_cube_texture_sphere

!> Finds the index of the first element of ind_out_maps() whose maps are being calculated. In the stellar RT this index is zero. In the dust RT it can be different. 
subroutine find_starting_kall(kall)
  integer :: kall 
  integer :: i0, i1 
  
  call set_i_opacity_arrays(i0,i1)

  do kall = 0, lnum_maps -1 

     if (ind_out_maps(kall) >= i0) exit 

  end do


end subroutine find_starting_kall



END MODULE visual_routines
