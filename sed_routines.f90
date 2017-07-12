MODULE sed_routines
  use smooth_grid_routines
  use io_routines
  use iso_fortran_env
  IMPLICIT NONE

  !> @param t0 Minimum dust temperature allowed.
  real(kind=real64), parameter :: t0=0
  !> @param t1 Maximum dust temperature allowed
  real(kind=real64), parameter :: t1=1E4
  !> @param tol Tolerance parameter used in zbrent function.
  real(kind=real64), parameter :: tol=1E-8
  !> @param conv_ufield_ifield Factor to be applied to the stellar emission radiation field energy density in order to obtain the radiation field intensity integrated over the entire solid angle.
  real(kind=real64) :: conv_ufield_ifield
  !> @param phot_energy Photon energy for each wavelength stored in lambda_arr_SI().
  real(kind=real64), allocatable :: phot_energy(:)
  !> @param abs_int_rad_stars Absorbed energy flux for a given radiation field produced by stellar radiation. It is in W/m/m^2 when calculating the equilibrium dust temperature but it can be transformed to photons/s/m/m^2 when calculating the moments of the dosage function (the latter is needed only for the dust stochastically heated calculation)
  real(kind=real64), allocatable :: abs_int_rad_stars(:)
  !> @param abs_int_rad_dust Same as abs_int_rad_stars() but for the radiation field produced by dust emission. 
  real(kind=real64), allocatable :: abs_int_rad_dust(:)
  !> @param Rd_arr Moments of dosage function. The zero moment is the total number of photons per unit time heating a dust grain. The first moment is the total heating rate in [W] and the second is the sum of f_lambda*E_lambda^2 where f_lambda is the total number of photons absorbed by a grain at wavelength lambda and E_lambda is the energy of a photon with wavelength lambda. 
  real(kind=real64) :: Rd_arr(0:2)
  !> @param large_grain_energy TRUE when large grain energy limit has to be used to derive the temperature distribution in the stochastically heated dust emission calculation. See Voit 1991.
  logical :: large_grain_energy
  !> @param Emin Minimum grain enthalpy
  real(kind=real64) :: Emin 
  !> @param Emax Maximum grain enthalpy
  real(kind=real64) :: Emax 
  !> @param Tmin Minimum grain temperature
  real(kind=real64) :: Tmin 
  !> @param Tmin_prev Minimum grain temperature in previous temperature iteration. 
  real(kind=real64) :: Tmin_prev 
  !> @param Tmax Maximum grain temperature
  real(kind=real64) :: Tmax 
  !> @param Tmax_prev Maximum grain temperature in previous temperature iteration. 
  real(kind=real64) :: Tmax_prev 
  !> @param E_arr Grain enthalpy array
  real(kind=real64), allocatable :: E_arr(:)
  !> @param T_arr Grain temperature array
  real(kind=real64), allocatable :: T_arr(:)
  !> @param Pt Grain temperature probability distribution. Note that each point is already multiplied by the temperature bin size. So, this is P(T)d(T) = f(E)dE.
  real(kind=real64), allocatable :: Pt(:)
  !> @param Qp_arr Planck-averaged Qabs array at the temperatures T_arr(). 
  real(kind=real64), allocatable :: Qp_arr(:)

  !> @param delta_E_arr_bin Bin sizes for array E_arr()
  real(kind=real64), allocatable :: delta_E_arr_bin(:)
  !> @param n_temp_pt Number of sampling point for the temperature probability distribution Pt(). This has to be > 100 to avoid problems with the calculation of Re0(), Re1(), Re2(). Value 250 should give good accuracy in most cases. 
  integer, parameter :: n_temp_pt = 300
  !> @param Rd_integrated Array of integrals of the photon absorption rate over wavelengths. Corresponding to integrals of "dosage function" in Voit 1991, but the order is inverted because the integration is over wavelength not energy. 
  real(kind=real64), allocatable :: Rd_integrated(:)
  !> @param aa The elements of this matrix are the photon absorption rate per unit energy interval. Do not confuse this with A matrix in Gahathakurta & Draine. This term is only used to calculate source term SeE (equ 53 Voit 1991). 
  real(kind=real64), allocatable :: aa(:,:)
  !> @param bb Energy level transition matrix as defined in Gahathakurta & Draine 1989. BB(i,j) gives the total transition rate from level j to all levels >=i. 
  real(kind=real64), allocatable :: bb(:,:)
  !> @param Re0 Zero moment of the dosage funtion integrated up to epsilon (see Voit 1991)
  real(kind=real64), allocatable :: Re0(:)
  !> @param Re1 First moment of the dosage funtion integrated up to epsilon (see Voit 1991)
  real(kind=real64), allocatable :: Re1(:)
  !> @param Re2 Second moment of the dosage funtion integrated up to epsilon (see Voit 1991)
  real(kind=real64), allocatable :: Re2(:)
  !> @param Edot_arr Cooling rates at each temperature defined in the interval where to calculate the temperature probability distribution
  real(kind=real64), allocatable :: Edot_arr(:)
  !> @param tot_dust_em Dust emission spectra after integration over the grain size distributions. 
  real(kind=real64), allocatable :: tot_dust_em(:)
  !> @param dust_em_arr_fa Dust emission spectra for each grain size of the grain size distribution. First index is the grain size index; second index is the wavelength.  
  real(kind=real64), allocatable :: dust_em_arr_fa(:,:)
  !> @param int_rf_uv Wavelength integrated radiation field intensity in the UV up to 4430 A for each cell. 
  real(kind=real64), allocatable :: int_rf_uv(:)
  !> @param int_rf_opt Wavelength integrated radiation field intensity in the optical from 4430 A to longer wavelnegths for each cell. Note that the radiation field produced by the dust emission is not included. 
  real(kind=real64), allocatable :: int_rf_opt(:)
  !> @param rf_uv_arr list of wavelength integrated UV radiation field intensities defining the bins used in the adaptive SED approach
  real(kind=real64), allocatable :: rf_uv_arr(:)
  !> @param rf_opt_arr list of wavelength integrated optical radiation field intensities defining the bins used in the adaptive SED approach
  real(kind=real64), allocatable :: rf_opt_arr(:)
  !> @param u_av_uv_opt Average radiation field spectra of stellar emission in each bin. First index is wavelength. Second index is UV radiation field index corresponding to bin in rf_uv_arr(). Third index is optical radiation field index corresponding to bin in rf_opt_arr().
  real(kind=real64), allocatable :: u_av_uv_opt(:,:,:)
  !> @param u_av_dust Average radiation field spectra of dust emission in each bin. First index is wavelength. Second index is UV radiation field index corresponding to bin in rf_uv_arr(). Third index is optical radiation field index corresponding to bin in rf_opt_arr(). Note that, although the radiation field due to dust emission is included in the average radiation field, its integrated intensity is not considered when binning the radiation field.  
  real(kind=real64), allocatable :: u_av_dust(:,:,:)
  !> @param count_spectra_uv_opt Counter of spectra contained in each UV and optical radiation field bin defined by rf_uv_arr() and rf_opt_arr().
  real(kind=real64), allocatable :: count_spectra_uv_opt(:,:)
  

  !> \todo why is iso_fortran_env needed here ? 

  !$OMP THREADPRIVATE(abs_int_rad_stars,abs_int_rad_dust, large_grain_energy, Emin, Emax, Tmin, Tmax, E_arr, T_arr, Pt, Tmin_prev, Tmax_prev, delta_E_arr_bin,  Qp_arr, Rd_integrated,aa, bb, Re0, Re1, Re2, rd_arr, Edot_arr,tot_dust_em, dust_em_arr_fa) 
 
CONTAINS

  !> Calculates stellar or dust emission SED from the i_obs() arrays. If keyword print_sed is not set TRUE, then it just returns.
  !> \todo DONE_TO_TEST This is not going to work in the no_communication mode.... Fix this. NEW: Done. See if it works.  
  subroutine calc_sed
    integer :: i, k, i0, i1
    real(kind=real64), allocatable :: out_sed_arr(:,:), in_sed_arr(:,:)
    integer :: ierr
    
    if (.not. print_sed) return

    !skip if dust RT and cnflag_dust /= TRUE
    if (rt_algorithm_ID == rta_dust .or. rt_algorithm_ID == rta_dust2d) then
       if (.not. cnflag_dust) return
    endif
    
    if (main_prc) print *, 'calculating integrated emission SED...'

    ! allocate temporary arrays 
    allocate(out_sed_arr(0:lnum_tot-1,0:tot_ndir-1),in_sed_arr(0:lnum_tot-1,0:tot_ndir-1))
    out_sed_arr = 0
    in_sed_arr = 0

    ! derive sed

    call set_i_opacity_arrays(i0,i1)
    
    do i = 0, lnum-1
       if (no_communications .and. .not. main_prc) exit !if no communications, only main MPI process makes output files.   
       if (iq_sca_node(i)) then
          call set_wavelength_index(i,k) 
          i_obs = i_obs_arr(k,:,:)
          call sum_i_obs(i+i0,out_sed_arr)          
       endif

    end do

    ! reduce sed_arr_dir and sed_arr calculated by the different MPI processes (each process has updated only a part of the arrays corresponding to the wavelengths of the i_obs arrays stored locally.
   
    call mpi_allreduce(out_sed_arr, in_sed_arr, lnum_tot*tot_ndir, MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
     
    ! assign to the right array
    select case(rt_type)
    case(rtt_output_part2, rtt_read_i_obs_part2)
       sed_arr_dir = in_sed_arr
    case(rtt_i_obs, rtt_read_i_obs)
       sed_arr = in_sed_arr
    case default
       if (main_prc) then
          print *, 'STOP(calc_sed): rt_type not recognized!'
          print *, 'rt_type = ', rt_type
       endif
       call stop_prc

    end select

    deallocate(in_sed_arr, out_sed_arr)

    call print_done

  end subroutine calc_sed

  
  !> Allocates sed arrays sed_arr(), sed_arr_dir(), sed_arr_sca().
  subroutine create_sed_arr

    if (.not. allocated(sed_arr)) then 
       allocate(sed_arr(0:lnum_tot-1,0:tot_ndir-1),sed_arr_dir(0:lnum_tot-1,0:tot_ndir-1))
    endif
    sed_arr=0
    sed_arr_dir=0
       
  end subroutine create_sed_arr

!> ! Sums up the i_obs at different wavelengths to obtain the total output sed in the input directions. Used only for the "sed" algorithm. 
  subroutine calc_sed_arrays 
    
    integer :: i
    logical :: file_exists
    integer :: ierr
    character(len=lcar) :: filename
    
 !!! check needed i_obs files exist
    if (main_prc) then 
       do i = 0, lnum -1

          if (print_output_part2) then 
             filename = trim(adjustl(dir_runs))//file_i_obs_part2_arr(i)
             inquire(file=filename, exist=file_exists)
             if (.not. file_exists) then
                print *, 'STOP: ', trim(filename), ' not found!'
                STOP
             endif
          endif

          filename = trim(adjustl(dir_runs))//file_i_obs_arr(i)
          inquire(file=filename, exist=file_exists)
          if (.not. file_exists) then
             print *, 'STOP: ', trim(filename), ' not found!'
             STOP
          endif

       end do

    endif
    call mpi_barrier(MPI_COMM_WORLD, ierr)

!!! read i_obs files and calculate sed_arr and sed_arr_dir    
    if (print_output_part2) then 
       rt_type=rtt_read_i_obs_part2

       call read_output

       call calc_sed
    endif
    
    rt_type=rtt_read_i_obs

    call read_output

    call calc_sed

  end subroutine calc_sed_arrays

  !> Sums i_obs arrays to obtain emission SED.
  subroutine sum_i_obs(il,sed_arr)
    real(kind=real64) :: sed_arr(0:,0:) !!! NOTE: this is not the sed_arr defined above 
    integer :: il,j
    
    select case (units_i_obs)
       
       case ('erg/s/Hz/pc^2/sr')
    
          i_obs=i_obs*10.0_real64**(-7) ! this is to convert erg -> Joule
          i_obs=i_obs*10.0_real64**(26)/parsec**2.  ! Jy/sr HERE NOT MJy/sr

       case ('W/Hz/pc^2/sr')
    
          i_obs=i_obs*10.0_real64**(26)/parsec**2.  ! Jy/sr HERE NOT MJy/sr
          
       case ('W/m/pc^2/sr')

          i_obs = i_obs*lambda_arr_SI(il)**2/cspeed
          i_obs= i_obs*10.0_real64**(26)/parsec**2.
       
       case DEFAULT

       if (main_prc) print *, 'which units for i_obs ? input units_i_obs!'
       stop

    end select

    do j=0,tot_ndir-1 
       sed_arr(il,j)=sum(i_obs(0:tot_ncell-1,j)*(csize(:)/dist_obs)**2)
       if (tot_p_src > 0) then
          sed_arr(il,j)= sed_arr(il,j)+sum(i_obs(tot_ncell:tot_ncell+tot_p_src-1,j)*(csize(cell_src)/dist_obs)**2)
       endif

    end do

  end subroutine sum_i_obs

!> Calculates the dust emission luminosity density and stores it into dens_stars_arr(). 
  subroutine set_dust_emission

    call omp_set_num_threads(nproc) 
    
    select case(dust_heating_type_ID)
    case(dt_eff)
       call calc_dens_dustem
    case(dt_equ)
       call calc_dens_dustem_equ
    case(dt_sto)
       call calc_dens_dustem_sto
    case(dt_sto_lib)
       call calc_dens_dustem_sto_lib  
    case default
       if (main_prc) then
          print *, 'STOP(set_dust_emission): dust_heating_type not recognized!'
          print *, 'dust_heating_type = ', dust_heating_type
       endif
       call stop_prc
    end select
 
    call reduce_dens_stars_arr

    call assign_parent_dens_arr  

    call check_dens_stars_arr
    
  end subroutine set_dust_emission

  !> Compare the newly calculated dens_stars_arr() with dens_stars_arr_prev(). If the relative difference is smaller than conv_en_lim(), it assigns TRUE to cnflag_dust(), thus stopping the dust heating iterations. It not, it subtracts from dens_stars_arr() the luminosity already processed in the previous dust heating iterations and it updates dens_stars_arr_prev(). 
!> \todo maybe you should add remove_negative here when subtracting. 
  subroutine check_dens_stars_arr
    real(kind=real64), allocatable :: temp_dens_stars_arr(:,:)
    integer :: n_no_conv_el ! number of no convergence elements in dens_stars_arr 
    integer :: i

    if (rt_algorithm_ID == rta_i_obs_dust) return 

    if (main_prc) print *, 'checking convergence check_dens_star_arr....'
    
    ! check how many dens_stars_arr elements have not dust emission luminosity converged

    n_no_conv_el = 0

    do i = 0, tot_ncell-1

       if (cchild(i) /= -1 .or. (sum(dens_stars_arr(:,i)) == 0)) cycle 

       n_no_conv_el = n_no_conv_el + count((abs(dens_stars_arr(:,i)-dens_stars_arr_prev(:,i))/dens_stars_arr(:,i)) > conv_en_lim)
       
    end do
    
    if (main_prc) print *, 'Number of no-convergence elements in dens_stars_arr: ', n_no_conv_el

    if (n_no_conv_el > 0 .or. iterations_dustem == 0) then ! at least one dust heating run (even when no dust is present, so needed arrays are allocated for SED output) 
    
       allocate(temp_dens_stars_arr(0:lnum_dust-1,0:tot_ncell-1)) 

       temp_dens_stars_arr = dens_stars_arr

       dens_stars_arr = dens_stars_arr - dens_stars_arr_prev

       dens_stars_arr_prev = temp_dens_stars_arr
    
       deallocate(temp_dens_stars_arr)

       iterations = 0 
       iterations_dustem = iterations_dustem + 1 

       if (main_prc) then 
          print *, 'New dust heating iteration is going to start soon...'
          print *, 'Dust heating iteration number = ', iterations_dustem
       endif
          
    else

       if (main_prc) then
          print *, 'Dust heating convergence criteria fulfilled! Stop to dust heating iterations.'
       endif
       cnflag_dust = .TRUE.
       
    endif
       
    call print_done

  end subroutine check_dens_stars_arr
    
  
!> Calculates the constant factor conv_ufield_ifield(). Then transforms the stellar emission radiation field energy density u_final_uv_opt() into the integrated radiation field intensity in units of W/m^2/m    
  subroutine convert_ufield_ifield
    
    integer :: il 

    if (main_prc) print *, 'converting stellar emission ufield into ifield ...'
    
    select case (units_ufield)

    case ('erg/Hz/pc^3')

       conv_ufield_ifield = parsec**(-3)*cspeed*10.0_real64**(-7) ! W/Hz/m^2
 
!!$       ufield_arr=ufield_arr/parsec**3*cspeed   ! erg/s/Hz/m^2   XXX NOTE: this is now the integrated intensity of the radiation field NOT the energy density
!!$       ufield_arr=ufield_arr*10.0_real64**(-7)         ! W/Hz/m^2

    case ('J/Hz/pc^3')
       
       conv_ufield_ifield = parsec**(-3)*cspeed ! W/Hz/m^2

    case DEFAULT 

       if (main_prc) then 
          print *, 'STOP: which units for ufield ? Something wrong in calc_conv_ufield_ifield!'
          print *, 'units_ufield = ', trim(units_ufield)
       endif
       call stop_prc

    end select

    u_final_uv_opt = u_final_uv_opt*conv_ufield_ifield  ! W/Hz/m^2

    do il=0,lnum_stars-1  ! convert into W/m/m^2
       u_final_uv_opt(il,:)=u_final_uv_opt(il,:)*cspeed/((lambda_arr_SI(il)**2))   ! W/m/m^2  
    enddo

    call print_done

  end subroutine convert_ufield_ifield
  
  !> Calculates the dust emission luminosity density in the case of the effective grain mode.
  subroutine calc_dens_dustem 
      
    integer :: i,j, i0,im
    real(kind=real64) :: t_dust,nu,abs_en

    if (main_prc) print *, 'calculating effective grain dust emission...'
    
    i0 = i_lambda_dust(0)

    dens_stars_arr = 0

    !$OMP PARALLEL DEFAULT(NONE), &
    !$OMP PRIVATE(im, i, abs_en, t_dust, j) &
    !$OMP SHARED(lnum, lnum_stars, lnum_dust,np_mpi, id_mpi, kabs_arr, u_final_uv_opt, iterations_dustem, delta_lambda_bin_stars, delta_lambda_bin_dust,dens_stars_arr, kext_ref, dens_ref, lambda_arr_SI, lnum_tot,i0, cchild, tot_ncell, u_final_arr, main_prc)

    allocate(abs_int_rad_stars(0:lnum_stars-1), abs_int_rad_dust(0:lnum_dust-1))
    abs_int_rad_stars = 0
    abs_int_rad_dust =0

    !$OMP DO SCHEDULE(DYNAMIC,chunk_size)
   
    do im=0, tot_ncell/np_mpi  !!! loop on cells NOTE: no -1 here     
       
       i = im*np_mpi+id_mpi

       if (i > tot_ncell -1) cycle 
       if ((cchild(i) /= -1).or.(dens_ref(i) == 0)) cycle

       ! NOTE: u_final_uv_opt has to be in units per wavelength not per frequency here. kabs_arr has to be in m^2 (effective grain cross section)

       call calc_t_dust_equil(kabs_arr,0._real64,u_final_uv_opt(:,i),u_final_arr(:,i)/parsec**3*cspeed,t_dust,abs_en) 

       do j=0,lnum-1
            
          dens_stars_arr(j,i)=4*pi*kabs_arr(j+i0)/kext_ref*dens_ref(i)*bplanck(t_dust, lambda_arr_SI(j+i0))  ! W/m^2/pc/m  Remember that dens(i)*csize(i)= tau_cell. Dens_ref is in units of 1/pc usually. It does not matter if it is not, multiplication by csize later which has to be in the corresponding units. Note also 4pi factor.
   
       end do
       
    end do

    !$OMP END DO NOWAIT

    deallocate(abs_int_rad_stars, abs_int_rad_dust)

    !$OMP END PARALLEL    

    call set_units_dens_stars_arr
    
    call print_done

  end subroutine calc_dens_dustem
  

  !>  Calculates the equilibrium emission for the grain mixture case. 
  subroutine calc_dens_dustem_equ
    
   integer :: i,j,ic,ig,imax,im,i0
   real(kind=real64) :: t_dust,nu,abs_en


   if (main_prc) print *, 'calculating grain mixture equilibrium dust emission...'  
   ! set dust wavelength grid starting index
   i0 = i_lambda_dust(0)

   dens_stars_arr = 0

   !$OMP PARALLEL DEFAULT(NONE), &
   !$OMP PRIVATE(i,j,ic,ig,abs_en,t_dust,imax,im), &
   !$OMP SHARED(tot_ncell,cchild,n_dust_comp,qabs_arr_fa,n_dust_size,dust_fa,delta_dust_size,tau_nh_ref,lambda_arr_SI,dust_size_fa,n_dust_maxsize_fa,lnum,np_mpi, id_mpi, dens_stars_arr, lnum_tot, i0,delta_lambda_bin_stars, delta_lambda_bin_dust, iterations_dustem, u_final_uv_opt, dens_ref, lnum_stars, lnum_dust,u_final_arr)

   allocate(dust_em_arr_fa(0:n_dust_maxsize_fa-1,0:lnum-1),tot_dust_em(0:lnum-1))
   dust_em_arr_fa = 0 
   tot_dust_em = 0
   allocate(abs_int_rad_stars(0:lnum_stars-1), abs_int_rad_dust(0:lnum_dust-1))
   abs_int_rad_stars = 0
   abs_int_rad_dust = 0 

!$OMP DO SCHEDULE(DYNAMIC,chunk_size)

   do im=0, tot_ncell/np_mpi  !!! loop on cells NOTE: no -1 here   

      i = im*np_mpi+id_mpi
      
      if (i > tot_ncell -1) cycle 
      if ((cchild(i) /= -1).or.(dens_ref(i) == 0)) cycle

      if (mod(i,1000) == 0) print *, i
       
      do ic=0,n_dust_comp-1 

         do ig=0,n_dust_size(ic)-1 
          
            call calc_t_dust_equil(qabs_arr_fa(ic,ig,:),dust_size_fa(ic,ig),u_final_uv_opt(:,i),u_final_arr(:,i)/parsec**3*cspeed,t_dust,abs_en)  
            
            do j=0,lnum-1
                dust_em_arr_fa(ig,j)=pi*(dust_size_fa(ic,ig)**2)*qabs_arr_fa(ic,ig,j+i0)*bplanck(t_dust, lambda_arr_SI(j+i0))  ! W/m/sr
                dust_em_arr_fa(ig,j)=dust_em_arr_fa(ig,j)*dust_fa(ic,ig) !W/m/sr*(m*H)^-1

             end do          
             
          end do         

          imax=n_dust_size(ic)

          do j=0, lnum-1 
        
             tot_dust_em(j)=4*pi*sum(dust_em_arr_fa(0:imax-1,j)*delta_dust_size(ic,0:imax-1))          ! this is an integration on the grain size and solid angle [W/m/H]

          end do

          dens_stars_arr(:,i)=dens_stars_arr(:,i)+tot_dust_em/tau_nh_ref*dens_ref(i) ! W/m/m^2/pc Same comment as above for calc_dens_dustem.  

       end do
       
    end do

    !$OMP END DO NOWAIT

    deallocate(dust_em_arr_fa,tot_dust_em)
    deallocate(abs_int_rad_stars, abs_int_rad_dust)

    !$OMP END PARALLEL

    call set_units_dens_stars_arr

    call print_done

  end subroutine calc_dens_dustem_equ

!> Allocates arrays used within the stochastically heated dust emission routines calc_dens_dustem_sto() and calc_dens_dustem_sto_lib(). 
  subroutine allocate_sto_arrays

    allocate(dust_em_arr_fa(0:n_dust_maxsize_fa-1,0:lnum-1),tot_dust_em(0:lnum-1))
    dust_em_arr_fa = 0 
    tot_dust_em = 0
    allocate(abs_int_rad_stars(0:lnum_stars-1), abs_int_rad_dust(0:lnum_dust-1))
    abs_int_rad_stars = 0
    abs_int_rad_dust =0
    allocate(E_arr(0:n_temp_pt-1), Pt(0:n_temp_pt-1), T_arr(0:n_temp_pt-1), delta_E_arr_bin(0:n_temp_pt-1), qp_arr(0:n_temp_pt-1), edot_arr(0:n_temp_pt-1))
    E_arr = 0 
    Pt = 0
    T_arr = 0
    delta_E_arr_bin = 0
    Qp_arr = 0
    Edot_arr = 0
    tmin_prev = -1 
    tmax_prev = -1 
    allocate(Rd_integrated(0:lnum_tot-1))
    rd_integrated = 0 
    allocate(aa(0:n_temp_pt-1,0:n_temp_pt-1), bb(0:n_temp_pt-1,0:n_temp_pt-1), Re0(0:n_temp_pt-1), Re1(0:n_temp_pt-1), Re2(0:n_temp_pt-1))
    aa = 0
    bb = 0 
    Re0 = 0
    Re1 = 0 
    Re2 = 0 

  end subroutine allocate_sto_arrays

!> Deallocates arrays used within the stochastically heated dust emission routines calc_dens_dustem_sto() and calc_dens_dustem_sto_lib().
  subroutine deallocate_sto_arrays

    deallocate(dust_em_arr_fa,tot_dust_em)
    deallocate(abs_int_rad_stars, abs_int_rad_dust)
    deallocate(E_arr, Pt, T_arr, delta_E_arr_bin,  qp_arr, Edot_arr)
    deallocate(Rd_integrated)
    deallocate(aa, bb, Re0, Re1, Re2)
    
  end subroutine deallocate_sto_arrays


  !>  Calculates the stochastically heated emission at each position.  
  subroutine calc_dens_dustem_sto
    
   integer :: i,j,ic,ig,imax,im,i0,it
   real(kind=real64) :: t_dust,nu,abs_en
   
   if (main_prc) print *, 'calculating grain mixture stochastically heated dust emission...'  
   ! set dust wavelength grid starting index
   i0 = i_lambda_dust(0)
   
   dens_stars_arr = 0

   !$OMP PARALLEL DEFAULT(NONE), &
   !$OMP PRIVATE(i,j,it,ic,ig,abs_en,t_dust,imax,im), &
   !$OMP SHARED(tot_ncell,cchild,n_dust_comp,qabs_arr_fa,n_dust_size,dust_fa,delta_dust_size,tau_nh_ref,lambda_arr_SI,dust_size_fa,n_dust_maxsize_fa,lnum,np_mpi, id_mpi, dens_stars_arr, lnum_tot, i0,delta_lambda_bin_stars, delta_lambda_bin_dust, iterations_dustem, u_final_uv_opt, dens_ref, lnum_stars, lnum_dust,u_final_arr)
   
!!$   allocate(dust_em_arr_fa(0:n_dust_maxsize_fa-1,0:lnum-1),tot_dust_em(0:lnum-1))
!!$   dust_em_arr_fa = 0 
!!$   tot_dust_em = 0
!!$   allocate(abs_int_rad_stars(0:lnum_stars-1), abs_int_rad_dust(0:lnum_dust-1))
!!$   abs_int_rad_stars = 0
!!$   abs_int_rad_dust =0
!!$   allocate(E_arr(0:n_temp_pt-1), Pt(0:n_temp_pt-1), T_arr(0:n_temp_pt-1), delta_E_arr_bin(0:n_temp_pt-1), qp_arr(0:n_temp_pt-1), edot_arr(0:n_temp_pt-1))
!!$   E_arr = 0 
!!$   Pt = 0
!!$   T_arr = 0
!!$   delta_E_arr_bin = 0
!!$   Qp_arr = 0
!!$   Edot_arr = 0
!!$   tmin_prev = -1 
!!$   tmax_prev = -1 
!!$   allocate(Rd_integrated(0:lnum_tot-1))
!!$   rd_integrated = 0 
!!$   allocate(aa(0:n_temp_pt-1,0:n_temp_pt-1), bb(0:n_temp_pt-1,0:n_temp_pt-1), Re0(0:n_temp_pt-1), Re1(0:n_temp_pt-1), Re2(0:n_temp_pt-1))
!!$   aa = 0
!!$   bb = 0 
!!$   Re0 = 0
!!$   Re1 = 0 
!!$   Re2 = 0 

   call allocate_sto_arrays 

   !$OMP DO SCHEDULE(DYNAMIC,chunk_size)
   do im=0, tot_ncell/np_mpi  !!! loop on cells NOTE: no -1 here   

      i = im*np_mpi+id_mpi
      
      if (i > tot_ncell -1) cycle 
      if ((cchild(i) /= -1).or.(dens_ref(i) == 0)) cycle

      if (mod(i,1000) == 0) print *, i
       
      do ic=0,n_dust_comp-1 

         large_grain_energy = .TRUE.
         dust_em_arr_fa = 0
         tmin_prev = -1 
         tmax_prev = -1 

         do ig=n_dust_size(ic)-1,0,-1 ! note: this is reversed compared to above. The reason for this is that larger grains stay closer to equilibrium than smaller grains. So, for those grains the Gaussian approximation for f(E)dE is used (see Voit 1991).  
          
            ! equilibrium temperature
            call calc_t_dust_equil(qabs_arr_fa(ic,ig,:),dust_size_fa(ic,ig),u_final_uv_opt(:,i),u_final_arr(:,i)/parsec**3*cspeed,t_dust,abs_en)  ! abs_en in [W/m^2] 

            if (abs_en == 0) cycle 
            
            ! moments of dosage function 
            call calc_rd_arr(dust_size_fa(ic,ig))

            if (large_grain_energy) then 

               ! calculate gaussian approximation 
               call calc_gaussian_fE(t_dust,ic,ig)

            endif

            if (.not. large_grain_energy) then ! IMPORTANT: do not join this IF statement with that before. Large_grain_energy can be assigned to FALSE in the previous IF statement. 

               call calc_full_fE(ic,ig)

            endif

                     
            do it = 0, n_temp_pt-1
               do j=0,lnum-1

                  dust_em_arr_fa(ig,j)=dust_em_arr_fa(ig,j)+pi*(dust_size_fa(ic,ig)**2)*qabs_arr_fa(ic,ig,j+i0)*pt(it)*bplanck(t_arr(it), lambda_arr_SI(j+i0))*dust_fa(ic,ig)  !W/m/sr*(m*H)^-1

               end do
            end do
             
         end do         

          imax=n_dust_size(ic)

          do j=0, lnum-1 
        
             tot_dust_em(j)=4*pi*sum(dust_em_arr_fa(0:imax-1,j)*delta_dust_size(ic,0:imax-1))          ! this is an integration on the grain size and solid angle [W/m/H]

          end do

          dens_stars_arr(:,i)=dens_stars_arr(:,i)+tot_dust_em/tau_nh_ref*dens_ref(i) ! W/m/m^2/pc Same comment as above for calc_dens_dustem.  

       end do

    end do

    !$OMP END DO NOWAIT

    call deallocate_sto_arrays

!!$    deallocate(dust_em_arr_fa,tot_dust_em)
!!$    deallocate(abs_int_rad_stars, abs_int_rad_dust)
!!$    deallocate(E_arr, Pt, T_arr, delta_E_arr_bin,  qp_arr, Edot_arr)
!!$    deallocate(Rd_integrated)
!!$    deallocate(aa, bb, Re0, Re1, Re2)

    !$OMP END PARALLEL

    call set_units_dens_stars_arr

    call print_done

  end subroutine calc_dens_dustem_sto

!>  Calculates the stochastically heated emission at each position using the SED adaptive library approach of Natale et al.(2015).  
  subroutine calc_dens_dustem_sto_lib
    
   integer :: i,j,ic,ig,imax,im,i0,it
   integer :: iuv, iopt, iuv0, iuv1, iopt0, iopt1
   real(kind=real64) :: t_dust,nu,abs_en
   integer :: tot_nbins 
   
   if (main_prc) print *, 'calculating grain mixture stochastically heated dust emission...'  
   ! set dust wavelength grid starting index
   i0 = i_lambda_dust(0)
   
   dens_stars_arr = 0

   ! binning of stellar emission radiation field 
   call bin_rad_field

   ! set total number of UV/optical bins 
   tot_nbins = n_int_rf_bins**2

   ! allocate/initialize dust emission SED library array
   if (.not. allocated(tot_dust_em_sed)) allocate(tot_dust_em_sed(0:lnum-1, 0:n_int_rf_bins-1, 0:n_int_rf_bins-1))
   tot_dust_em_sed = 0 
   
   ! calculate dust emission SED library 

   !$OMP PARALLEL DEFAULT(NONE), &
   !$OMP PRIVATE(i,j,it,ic,ig,abs_en,t_dust,imax,im,iuv, iopt), &
   !$OMP SHARED(tot_ncell,cchild,n_dust_comp,qabs_arr_fa,n_dust_size,dust_fa,delta_dust_size,tau_nh_ref,lambda_arr_SI,dust_size_fa,n_dust_maxsize_fa,lnum,np_mpi, id_mpi, dens_stars_arr, lnum_tot, i0,delta_lambda_bin_stars, delta_lambda_bin_dust, iterations_dustem, u_final_uv_opt, dens_ref, lnum_stars, lnum_dust,u_final_arr,tot_nbins, n_int_rf_bins, u_av_uv_opt, u_av_dust, tot_dust_em_sed, count_spectra_uv_opt)

!!$   allocate(dust_em_arr_fa(0:n_dust_maxsize_fa-1,0:lnum-1),tot_dust_em(0:lnum-1))
!!$   dust_em_arr_fa = 0 
!!$   tot_dust_em = 0
!!$   allocate(abs_int_rad_stars(0:lnum_stars-1), abs_int_rad_dust(0:lnum_dust-1))
!!$   abs_int_rad_stars = 0
!!$   abs_int_rad_dust =0
!!$   allocate(E_arr(0:n_temp_pt-1), Pt(0:n_temp_pt-1), T_arr(0:n_temp_pt-1), delta_E_arr_bin(0:n_temp_pt-1), qp_arr(0:n_temp_pt-1), edot_arr(0:n_temp_pt-1))
!!$   E_arr = 0 
!!$   Pt = 0
!!$   T_arr = 0
!!$   delta_E_arr_bin = 0
!!$   Qp_arr = 0
!!$   Edot_arr = 0
!!$   tmin_prev = -1 
!!$   tmax_prev = -1 
!!$   allocate(Rd_integrated(0:lnum_tot-1))
!!$   rd_integrated = 0 
!!$   allocate(aa(0:n_temp_pt-1,0:n_temp_pt-1), bb(0:n_temp_pt-1,0:n_temp_pt-1), Re0(0:n_temp_pt-1), Re1(0:n_temp_pt-1), Re2(0:n_temp_pt-1))
!!$   aa = 0
!!$   bb = 0 
!!$   Re0 = 0
!!$   Re1 = 0 
!!$   Re2 = 0 

   call allocate_sto_arrays

   !$OMP DO SCHEDULE(DYNAMIC,10)
   do im=0, tot_nbins/np_mpi  !!! loop on cells NOTE: no -1 here   

      i = im*np_mpi+id_mpi
      
      if (i > tot_nbins -1) cycle

      iuv = i/n_int_rf_bins 
      iopt = mod(i, n_int_rf_bins)   

      if (count_spectra_uv_opt(iuv,iopt) == 0) cycle 

      !print *, iuv, iopt

      if (mod(i,100) == 0) print *, i
       
      do ic=0,n_dust_comp-1 

         large_grain_energy = .TRUE.
         dust_em_arr_fa = 0
         tmin_prev = -1 
         tmax_prev = -1 

         do ig=n_dust_size(ic)-1,0,-1 ! note: this is reversed compared to above. The reason for this is that larger grains stay closer to equilibrium than smaller grains. So, for those grains the Gaussian approximation for f(E)dE is used (see Voit 1991).  
          
            ! equilibrium temperature
            call calc_t_dust_equil(qabs_arr_fa(ic,ig,:),dust_size_fa(ic,ig),u_av_uv_opt(:,iuv,iopt),u_av_dust(:,iuv,iopt)/parsec**3*cspeed,t_dust,abs_en)  ! abs_en in [W/m^2] 

            if (abs_en == 0) cycle 
            
            ! moments of dosage function 
            call calc_rd_arr(dust_size_fa(ic,ig))

            if (large_grain_energy) then 

               ! calculate gaussian approximation 
               call calc_gaussian_fE(t_dust,ic,ig)

            endif

            if (.not. large_grain_energy) then ! IMPORTANT: do not join this IF statement with that before. Large_grain_energy can be assigned to FALSE in the previous IF statement. 

               call calc_full_fE(ic,ig)

            endif

                     
            do it = 0, n_temp_pt-1
               do j=0,lnum-1

                  dust_em_arr_fa(ig,j)=dust_em_arr_fa(ig,j)+pi*(dust_size_fa(ic,ig)**2)*qabs_arr_fa(ic,ig,j+i0)*pt(it)*bplanck(t_arr(it), lambda_arr_SI(j+i0))*dust_fa(ic,ig)  !W/m/sr*(m*H)^-1

               end do
            end do
             
         end do

         imax=n_dust_size(ic)

         do j=0, lnum-1 
        
            tot_dust_em(j)=4*pi*sum(dust_em_arr_fa(0:imax-1,j)*delta_dust_size(ic,0:imax-1))          ! this is an integration on the grain size and solid angle [W/m/H]

         end do

         tot_dust_em_sed(:, iuv,iopt) = tot_dust_em_sed(:, iuv,iopt) + tot_dust_em
         ! dens_stars_arr(:,i)=dens_stars_arr(:,i)+tot_dust_em/tau_nh_ref*dens_ref(i) ! W/m/m^2/pc Same comment as above for calc_dens_dustem.  

      end do

   end do

    !$OMP END DO NOWAIT

!!$    deallocate(dust_em_arr_fa,tot_dust_em)
!!$    deallocate(abs_int_rad_stars, abs_int_rad_dust)
!!$    deallocate(E_arr, Pt, T_arr, delta_E_arr_bin,  qp_arr, Edot_arr)
!!$    deallocate(Rd_integrated)
!!$    deallocate(aa, bb, Re0, Re1, Re2)

    call deallocate_sto_arrays

    !$OMP END PARALLEL

    ! reduce dust emission SED library 
    call reduce_tot_dust_em_sed

    ! assign spectra to each dusty cell according to MPI redistribution so that reduction of dens_stars_arr works in the same way as for all calc_dens_dustem routines 

    do im = 0, tot_ncell/np_mpi 
       
       i = im*np_mpi+id_mpi
      
       if (i > tot_ncell -1) cycle 
       if ((cchild(i) /= -1).or.(dens_ref(i) == 0)) cycle

       call value_locate(int_rf_uv(i), rf_uv_arr, iuv0, iuv1, .FALSE.)
       call value_locate(int_rf_opt(i), rf_opt_arr, iopt0, iopt1, .FALSE.)

       iuv = iuv0
       iopt = iopt0

       dens_stars_arr(:,i)=tot_dust_em_sed(:,iuv,iopt)/tau_nh_ref*dens_ref(i) ! W/m/m^2/pc Same comment as above for calc_dens_dustem.  


    end do

    call set_units_dens_stars_arr

    call print_done

  end subroutine calc_dens_dustem_sto_lib





!> Calculates the equilibrium dust temperature given the absorption Q coefficient, the radiation field intensity of the stellar emission and of the dust emission in [W/m/m^2]. The input qabs_arr coefficient can also have units [e.g. m^2], as for the effective grain emission calculation, without changing the output equilibrium temperature. The output absorbed energy abs_en is in units of W/m^2 times the units of the input qabs_arr.  
!> \todo remove dust_size here, not needed ! 
subroutine calc_t_dust_equil(qabs_arr, dust_size, rf_stars, rf_dust, t_dust, abs_en)
real(kind=real64) :: qabs_arr(0:), rf_stars(0:), rf_dust(0:), t_dust, dust_size
!real(kind=real64) :: abs_int_rad_stars(0:lnum_stars-1) ! absorbed stellar emission energy / wavelength
!real(kind=real64) :: abs_int(0:lnum-1),abs_en
real(kind=real64) :: abs_en
integer :: i0, i 

i0 = i_lambda_dust(0)

abs_int_rad_stars=qabs_arr(0:lnum_stars-1)*rf_stars  ! W/m/m^2

if (iterations_dustem > 0) then ! the >0 is because iterations_dustem is added + 1 before the following dust heating iteration.  
   abs_int_rad_dust=qabs_arr(i0:lnum_tot-1)*rf_dust  ! W/m/m^2  ! absorbed dust emission energy / wavelength
else 
   abs_int_rad_dust = 0 
endif

abs_en=sum(abs_int_rad_stars*delta_lambda_bin_stars)+sum(abs_int_rad_dust*delta_lambda_bin_dust)

if (abs_en < 0) then 
   print *, 'STOP! something weird with the equilibrium dust temperature calculation:'
   print *, 'abs_en=',abs_en
   print *, 'id_mpi =', id_mpi
   stop
endif

t_dust=zbrent_tdust(qabs_arr(i0:lnum_tot-1),abs_en) 

end subroutine calc_t_dust_equil

!> Calculates the moments of the dosage function (see Voit 1991, ApJ, 379, 122, needed for the dust stochastically heated calculation). 
subroutine calc_rd_arr(dust_size)
integer :: i0, i 
real(kind=real64) ::  dust_size


i0 = i_lambda_dust(0)

! Find moments of dosage function (Voit 1991, note that integration is on wavelength not energy as in the paper). 

abs_int_rad_stars = abs_int_rad_stars/phot_energy(0:lnum_stars-1)  ! photons/ m/m^2/s

abs_int_rad_dust=abs_int_rad_dust/phot_energy(i0:lnum_tot-1)

rd_arr(0) = sum(abs_int_rad_stars*delta_lambda_bin_stars) + sum(abs_int_rad_dust*delta_lambda_bin_dust)

rd_arr(1) = sum(abs_int_rad_stars*phot_energy(0:lnum_stars-1)*delta_lambda_bin_stars) + sum(abs_int_rad_dust*phot_energy(i0:lnum_tot-1)*delta_lambda_bin_dust)

rd_arr(2) = sum(abs_int_rad_stars*(phot_energy(0:lnum_stars-1)**2)*delta_lambda_bin_stars) + sum(abs_int_rad_dust*(phot_energy(i0:lnum_tot-1)**2)*delta_lambda_bin_dust)

rd_arr = rd_arr*pi*dust_size**2

!print *, rd_arr(0), rd_arr(1)/(1.6*1E-19), rd_arr(2)/(1.6*1E-19)**2. ! results in eV units

end subroutine calc_rd_arr


  !> Sets the units of dens_stars_arr() in the dust RT algorithms. These units are W/m/pc^3 where meters refers to the wavelength.
  !> \todo DONE Maybe should output u_final_arr and i_obs_arr in the same units for stellar and dust emission. Units of i_obs_arr for dust emission always in W/m/pc^2/sr.
  subroutine set_units_dens_stars_arr

    ! set right units
    select case(units_csize)
    case('pc')
       dens_stars_arr = dens_stars_arr*parsec**2  ! W/m/pc^3
    case default
       if (main_prc) print *, 'STOP(set_units_dens_star_arr): units_csize not recognized! Add new case here.!'
    end select

  end subroutine set_units_dens_stars_arr
  
  !> Calculates the difference between absorbed and emitted energy given the temperature T_dust.  
  real(kind=real64) function abs_en_diff(t_dust,kabs_arr_planck,abs_en,em_int)

    real(kind=real64) :: abs_en,em_en, t_dust, kabs_arr_planck(0:)
    integer :: i,ig, i0 
    real(kind=real64) :: em_int(0:)

    !calculate emitted luminosity

    i0 = i_lambda_dust(0)  ! starting index dust wavelength grid
    
    do i=0,lnum-1
       em_int(i)=(kabs_arr_planck(i)*bplanck(T_dust, lambda_arr_SI(i+i0))) ! W/m/sr
    end do
    
    em_en=4*pi*sum(em_int*delta_lambda_bin_dust)  ! note the 4pi coefficient here!
    
    abs_en_diff=em_en-abs_en
    
   
  end function abs_en_diff

!> Returns black body specific intensity in W/m^2/m/sr at wavelength la [m] and for temperature T_source [K]  (SI units)  
real(kind=real64) function bplanck(T_source, la)

  real(kind=real64) :: t_source, la,a1,a2

  a1=2*hplanck*cspeed**2./la**5.
  a2=exp(hplanck*cspeed/(la*kboltz*T_source))-1

  bplanck=a1*(a2**(-1.))  ! W/m^3/sr

end function bplanck

!>  Sets the lambda array in the mks units ( lambda_arr_SI()). It also sets delta_lambda_bin(), delta_lambda_bin_stars() and delta_lambda_bin_dust(). 
  subroutine set_lambda_arr_si
    integer :: i0

  allocate(lambda_arr_SI(0:lnum_tot-1),lambda_arr_SI_bin(0:lnum_tot-2),delta_lambda_bin(0:lnum_tot-1))
  allocate(delta_lambda_bin_stars(0:lnum_stars-1), delta_lambda_bin_dust(0:lnum_dust-1))

  select case (units_lambda)

  case ('um')
     lambda_arr_SI = lambda_arr*1E-6
     lambda_ref_SI = lambda_ref*1E-6
     
  case DEFAULT
     
     print *, 'which units for lambda ?'

     stop

  end select

  ! return if no dust RT will be performed (for example if there are not enough wavelengths in the input). 
  if (no_dust_rt .or. grid_creation) return 

  lambda_arr_SI_bin=10_real64**((log10(lambda_arr_SI(1:lnum_tot-1))+log10(lambda_arr_SI(0:lnum_tot-2)))/2.) ! if you change this with below line, update documentation for delta_lambda_bin as well. 
 ! lambda_arr_SI_bin=(lambda_arr_SI(1:lnum_tot-1)+lambda_arr_SI(0:lnum_tot-2))/2.
  delta_lambda_bin(1:lnum_tot-2)=lambda_arr_SI_bin(1:lnum_tot-2)-lambda_arr_SI_bin(0:lnum_tot-3)
  delta_lambda_bin(0)=lambda_arr_SI_bin(0)-lambda_arr_SI(0)
  delta_lambda_bin(lnum_tot-1)=lambda_arr_SI(lnum_tot-1)-lambda_arr_SI_bin(lnum_tot-2)

  delta_lambda_bin_stars(1:lnum_stars-2)=lambda_arr_SI_bin(1:lnum_stars-2)-lambda_arr_SI_bin(0:lnum_stars-3)
  delta_lambda_bin_stars(0)=lambda_arr_SI_bin(0)-lambda_arr_SI(0)
  delta_lambda_bin_stars(lnum_stars-1)=lambda_arr_SI(lnum_stars-1)-lambda_arr_SI_bin(lnum_stars-2)

  i0 = i_lambda_dust(0)
  
  delta_lambda_bin_dust(1:lnum_dust-2)=lambda_arr_SI_bin(1+i0:lnum_tot-2)-lambda_arr_SI_bin(0+i0:lnum_tot-3)
  delta_lambda_bin_dust(0)=lambda_arr_SI_bin(0+i0)-lambda_arr_SI(0+i0)
  delta_lambda_bin_dust(lnum_dust-1)=lambda_arr_SI(lnum_tot-1)-lambda_arr_SI_bin(lnum_tot-2)

  deallocate(lambda_arr_SI_bin)
  
end subroutine set_lambda_arr_si

!> Finds equilibrium dust temperature using Brent's method (modified version of the Van Wijngaarden–Dekker–Brent Method, see Numerical recipes in Fortran 77, Press et al.). 
real(kind=real64) FUNCTION ZBRENT_tdust(kabs_arr_planck,abs_en)

  IMPLICIT real (kind=real64) (A-H,O-Z)
   
  integer :: iter,itmax
  real(kind=real64) :: kabs_arr_planck(0:lnum-1),abs_en,em_int(0:lnum-1)
  PARAMETER (ITMAX=100,EPS=3.E-16)      
  A=t0
  B=t1
  FA=abs_en_diff(A,kabs_arr_planck,abs_en,em_int)
  FB=abs_en_diff(B,kabs_arr_planck,abs_en,em_int)
  !     error: Root not bracketed 
  IF(FB*FA.GT.0.) THEN
     WRITE(*,*) 'STOP(zbrent): root of FUNC not bracketed.'
     WRITE(*,*) 'lower bracket A = ',A,' FUNC(A) = ',FA 
     WRITE(*,*) 'upper bracket B = ',B,' FUNC(B) = ',FB 
     STOP
  END IF
  FC=FB
  DO ITER=1,ITMAX
     !print *, iter
     IF(FB*FC.GT.0.) THEN
        C=A
        FC=FA
        D=B-A
        E=D
     ENDIF
     IF(ABS(FC).LT.ABS(FB)) THEN
        A=B
        B=C
        C=A
        FA=FB
        FB=FC
        FC=FA
     ENDIF
     TOL1=2.*EPS*ABS(B)+0.5*TOL
     XM=.5*(C-B)
     IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
        ZBRENT_tdust=B
        RETURN
     ENDIF
     IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
        S=FB/FA
        IF(A.EQ.C) THEN
           P=2.*XM*S
           Q=1.-S
        ELSE
           Q=FA/FC
           R=FB/FC
           P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
           Q=(Q-1.)*(R-1.)*(S-1.)
        ENDIF
        IF(P.GT.0.) Q=-Q
        P=ABS(P)
        IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
           E=D
           D=P/Q
        ELSE
           D=XM
           E=D
        ENDIF
     ELSE
        D=XM
        E=D
     ENDIF
     A=B
     FA=FB
     IF(ABS(D) .GT. TOL1) THEN
        B=B+D
     ELSE
        B=B+SIGN(TOL1,XM)
     ENDIF
     FB=abs_en_diff(B,kabs_arr_planck,abs_en,em_int)
     
  end do
  !     error message for too many iterations
  WRITE(*,*) 'STOP(zbrent_tdust): number of iterations required to'
  WRITE(*,*) 'converge on the root of FUNC exceeded maximum'
  WRITE(*,*) 'allowed.  This can be remedied by increasing'
  WRITE(*,*) 'parameter ITMAX or by relaxing'
  WRITE(*,*) 'the convergence tolerance TOL'
  STOP      
 
END function zbrent_tdust

!> Reads the input dust model files.  Derive the dust model optical properties used in the calculation, interpolated at the same wavelengths as in lambda_arr() and same grain sizes as in the input grain size distributions. 
subroutine prepare_dust_model

  ! load tabulated opacity parameters 
  call load_opacity_param 

  ! load grain size distributions
  call load_fa_arr

  ! interpolate opacity values to same grain size as those for the grain size distribution
  call interpolate_q_grain_fa

  ! derive average opacities kabs, ksca, kext and scattering coefficient gsca
  call calc_total_opacity

  ! read input table of kabs, ksca, kext, gsca if required in the input file. keyword input_av_opacities(). 
  call read_av_opacities

  ! set albedo
  ksca_arr_norm = ksca_arr/kext_arr

  if (dust_heating_type_ID /= dt_sto .and. dust_heating_type_ID /= dt_sto_lib) return ! the following is needed only for the stochastically heated dust emission calculation

  ! derive Planck-averaged Qabs values for a set of wavelength.  
  call calc_planck_av_qabs

  ! derive photon energy as a function of wavelength 
  call calc_photon_energy

  ! derive specific enthalphy / specific heat capacity for the grains 
  call load_cT_hT_tables 

  
end subroutine prepare_dust_model


!> Loads the Qabs, Qsca, Qext and gsca factors from the tables in the dust opacity directory (from the TRUST benchmark project) and interpolate them to the input wavelength grid. Note that a further interpolation is needed to match the grain sizes to the values of the input tables. This is done in interpolate_q_grain_fa().
!> \todo Insert acknowledgements for using the opacity tables from TRUST and those from Draine and Li 2006. Where are the TRUST tables from ? Insert reference in output message when checking dust_model and dust_opacity_tables. 

subroutine load_opacity_param

  character(LEN=lcar) :: file_comp_arr(0:max_n_dust_comp-1)  
  integer :: i,u,j,k,id
  integer, parameter :: nhead = 15 
  real(kind=real64) :: arr6(0:5)
  real(kind=real64), allocatable :: tlambda_arr_in(:)
  real(kind=real64), allocatable :: tqabs_arr_in(:,:,:)
  real(kind=real64), allocatable :: tqsca_arr_in(:,:,:)
  real(kind=real64), allocatable :: tqext_arr_in(:,:,:)
  real(kind=real64), allocatable :: tgsca_arr_in(:,:,:)
  real(kind=real64), allocatable :: marr(:),carr(:)
 
  
  if (main_prc) print *, 'loading single grain opacity tables...'

  select case(dust_opacity_tables)
  case('TRUST')
     ! file names 
     file_q_gra='./DUST_OPACITY/TRUST/Gra_121_1201.dat'
     file_q_sil='./DUST_OPACITY/TRUST/suvSil_121_1201.dat'
     file_q_pah_neu='./DUST_OPACITY/TRUST/PAH_28_1201_neu.dat'

     ! initialise arrays
     n_dust_comp = 3 
     !allocate(n_dust_size_qabs(0:n_dust_comp-1))
     n_dust_size_qabs(0)=121   ! Graphite 
     n_dust_size_qabs(1)=121   ! Silicates 
     n_dust_size_qabs(2)=28    ! neutral PAH
     n_dust_wave_qabs=1201 ! number of wavelengths for which Qabs is tabulated 

  case('DraineLi06')
     ! file names 
     file_q_gra='./DUST_OPACITY/DraineLi06/Gra01'
     file_q_sil='./DUST_OPACITY/DraineLi06/Si01'
     file_q_pah_neu='./DUST_OPACITY/DraineLi06/PAHneu06'
     file_q_pah_ion='./DUST_OPACITY/DraineLi06/PAHion06'

     ! initialise arrays
     n_dust_comp = 4     
     !allocate(n_dust_size_qabs(0:n_dust_comp-1))
     !iq_dust_model = [1,1,1,1]
     n_dust_size_qabs(0)=81   ! Graphite 
     n_dust_size_qabs(1)=81   ! Silicates 
     n_dust_size_qabs(2)=30    ! neutral PAH
     n_dust_size_qabs(3)=30    ! ionized PAH
     n_dust_wave_qabs=1201 ! number of wavelengths for which Qabs is tabulated     
         
  case('user') ! all values but n_dust_comp
     n_dust_comp = count(iq_dust_model)
     
  case default
     if (main_prc) print *, 'STOP(load_opacity_param): dust_model not recognized!'
     call stop_prc
     
  end select

  file_comp_arr=(/file_q_gra, file_q_sil,file_q_pah_neu, file_q_pah_ion/)
  n_dust_maxsize_qabs=maxval(n_dust_size_qabs)   
  
  allocate(tqabs_arr_in(0:n_dust_comp-1,0:n_dust_maxsize_qabs-1,0:n_dust_wave_qabs-1))
  allocate(tlambda_arr_in(0:n_dust_wave_qabs-1))
  allocate(dust_size_qabs(0:n_dust_comp-1,0:n_dust_maxsize_qabs-1))

  allocate(qabs_arr_in(0:n_dust_comp-1,0:n_dust_maxsize_qabs-1,0:lnum_tot-1)) 
  allocate(qsca_arr_in(0:n_dust_comp-1,0:n_dust_maxsize_qabs-1,0:lnum_tot-1))
  allocate(qext_arr_in(0:n_dust_comp-1,0:n_dust_maxsize_qabs-1,0:lnum_tot-1))
  allocate(gsca_arr_in(0:n_dust_comp-1,0:n_dust_maxsize_qabs-1,0:lnum_tot-1))

 ! if (.not. use_lambda_grid) then 
     allocate(qabs_ref_in(0:n_dust_comp-1,0:n_dust_maxsize_qabs-1)) 
     allocate(qsca_ref_in(0:n_dust_comp-1,0:n_dust_maxsize_qabs-1))
     allocate(qext_ref_in(0:n_dust_comp-1,0:n_dust_maxsize_qabs-1))
     allocate(gsca_ref_in(0:n_dust_comp-1,0:n_dust_maxsize_qabs-1))
 ! endif
     
  allocate(tqsca_arr_in(0:n_dust_comp-1,0:n_dust_maxsize_qabs-1,0:n_dust_wave_qabs-1))
  allocate(tqext_arr_in(0:n_dust_comp-1,0:n_dust_maxsize_qabs-1,0:n_dust_wave_qabs-1))
  allocate(tgsca_arr_in(0:n_dust_comp-1,0:n_dust_maxsize_qabs-1,0:n_dust_wave_qabs-1))
  tqabs_arr_in=0
  tqsca_arr_in=0
  tqext_arr_in=0
  tgsca_arr_in=0
  tlambda_arr_in=0
  dust_size_qabs=0
  qabs_arr_in=0
  qsca_arr_in=0
  qext_arr_in=0
  gsca_arr_in=0

 ! if (.not. use_lambda_grid) then
     qabs_ref_in=0
     qsca_ref_in=0
     qext_ref_in=0
     gsca_ref_in=0
 ! endif
  
  ! read opacity coefficients from files
  i = 0
  do id = 0, max_n_dust_comp-1 

     if (.not. iq_dust_model(id)) cycle
     
     open(newunit=u,file=file_comp_arr(i),status='old')
     
     do j=0, nhead-1 ! skip header  
        read(u,*)
     end do

     do j=0, n_dust_size_qabs(id)-1 
        read(u,*) dust_size_qabs(i,j) ! grain size [um]
        read(u,*)
        do k=0,n_dust_wave_qabs-1 
           read(u,*) arr6 !x(=2pi a/wave) Wave (micron)    Q_abs     Q_sca    Q_ext    g
           
           tqabs_arr_in(i,j,k)=arr6(2)   ! Qabs
           tqsca_arr_in(i,j,k)=arr6(3)   ! Qsca
           tqext_arr_in(i,j,k)=arr6(4)   ! Qext
           tgsca_arr_in(i,j,k)=arr6(5)   ! g
           
           if ((i == 0).and.(j == 0)) then  ! reads lambda for i=0 
              tlambda_arr_in(k)=arr6(1)
              if (k > 0) then  ! check ascending order for wavelengths 
                 if (tlambda_arr_in(k) < tlambda_arr_in(k-1)) then 
                    if (main_prc) print *, 'STOP(load_opacity_param): wavelengths should be in ascending order in input dust opacity tables!'
                    call stop_prc
                 endif
              endif
           
           else 
              
              if (tlambda_arr_in(k) /= arr6(1)) then ! checks lambda is the same for i /= 0 
                 if (main_prc) print *, 'STOP(load_opacity_param): something wrong with qabs file. Not the same wavelength array in all opacity tables'                 
                 print *, i,j,k
                 call stop_prc
              endif
              
           endif
           
        end do ! nwave
        read(u,*)       
     end do ! nsize
     
     close(u)
     
     i = i + 1
  end do ! composition

!!! transform to SI units 
  tlambda_arr_in=tlambda_arr_in*1E-6  !!! um -> m 
  dust_size_qabs=dust_size_qabs*1E-6 !! um -> m 

  !! interpolate to find qabs_arr_in and the other array values for the RT wavelength grid in lambda_arr 

  allocate(marr(0:n_dust_comp-1),carr(0:n_dust_comp-1))
  marr=0
  carr=0

  k=0

  do i=0, n_dust_wave_qabs-1 

     if (k > lnum_tot-1) exit 

     if (tlambda_arr_in(i) < lambda_arr_SI(k)) cycle

     do j=0, n_dust_maxsize_qabs-1 

        marr=(tqabs_arr_in(:,j,i)-tqabs_arr_in(:,j,i-1))/(tlambda_arr_in(i)-tlambda_arr_in(i-1))
        carr=tqabs_arr_in(:,j,i)-marr*tlambda_arr_in(i)
        qabs_arr_in(:,j,k)=marr*lambda_arr_SI(k)+carr  !! !this Qabs array is evaluated at the same wavelenghts as the RT wavelength grid. However, it is not evaluated at the same grain sizes of the size distribution functions. Same for other quantities below. 
        marr=(tqsca_arr_in(:,j,i)-tqsca_arr_in(:,j,i-1))/(tlambda_arr_in(i)-tlambda_arr_in(i-1))
        carr=tqsca_arr_in(:,j,i)-marr*tlambda_arr_in(i)
        qsca_arr_in(:,j,k)=marr*lambda_arr_SI(k)+carr 

        marr=(tqext_arr_in(:,j,i)-tqext_arr_in(:,j,i-1))/(tlambda_arr_in(i)-tlambda_arr_in(i-1))
        carr=tqext_arr_in(:,j,i)-marr*tlambda_arr_in(i)
        qext_arr_in(:,j,k)=marr*lambda_arr_SI(k)+carr

        marr=(tgsca_arr_in(:,j,i)-tgsca_arr_in(:,j,i-1))/(tlambda_arr_in(i)-tlambda_arr_in(i-1))
        carr=tgsca_arr_in(:,j,i)-marr*tlambda_arr_in(i)
        gsca_arr_in(:,j,k)=marr*lambda_arr_SI(k)+carr 
      
     end do

     k=k+1

  end do

! interpolate to the value of lambda_ref 

!  if (.not. use_lambda_grid) then 

     do i=0, n_dust_wave_qabs-1 

        if (tlambda_arr_in(i) < lambda_ref_SI) cycle

        do j=0, n_dust_maxsize_qabs-1 

           marr=(tqabs_arr_in(:,j,i)-tqabs_arr_in(:,j,i-1))/(tlambda_arr_in(i)-tlambda_arr_in(i-1))
           carr=tqabs_arr_in(:,j,i)-marr*tlambda_arr_in(i)
           qabs_ref_in(:,j)=marr*lambda_ref_SI+carr  !! !this Qabs array is evaluated at the reference wavelength. However, it is not evaluated at the same grain sizes of the size distribution functions. Same for other quantities below. 
           marr=(tqsca_arr_in(:,j,i)-tqsca_arr_in(:,j,i-1))/(tlambda_arr_in(i)-tlambda_arr_in(i-1))
           carr=tqsca_arr_in(:,j,i)-marr*tlambda_arr_in(i)
           qsca_ref_in(:,j)=marr*lambda_ref_SI+carr 

           marr=(tqext_arr_in(:,j,i)-tqext_arr_in(:,j,i-1))/(tlambda_arr_in(i)-tlambda_arr_in(i-1))
           carr=tqext_arr_in(:,j,i)-marr*tlambda_arr_in(i)
           qext_ref_in(:,j)=marr*lambda_ref_SI+carr

           marr=(tgsca_arr_in(:,j,i)-tgsca_arr_in(:,j,i-1))/(tlambda_arr_in(i)-tlambda_arr_in(i-1))
           carr=tgsca_arr_in(:,j,i)-marr*tlambda_arr_in(i)
           gsca_ref_in(:,j)=marr*lambda_ref_SI+carr 
      
        end do

        exit  ! only one wavelength to calculate. See also "cycle" above. 
        
     end do

!  endif
  
  deallocate(tqabs_arr_in,tqsca_arr_in,tqext_arr_in,tgsca_arr_in,tlambda_arr_in,marr,carr)

  call print_done
  
end subroutine load_opacity_param

!>  Loads the grain size distribution from the standard TRUST tables or from user provided files.
!> \todo How to check for blank spaces at the end of ASCII files ? This is valid for all of them e.g. file_dir_out 
subroutine load_fa_arr

  character(LEN=lcar) :: file_comp_arr(0:max_n_dust_comp-1)
  integer :: u,i,j,nlines, id
  integer, parameter :: nhead = 4 ! file header lines 


  if (main_prc) print *, 'load grain size distributions...'
  
  ! Assign grain size distribution filenames for pre-defined dust models. For user defined dust models they are selected in the input file. 
  select case (dust_model)
  case('TRUST') 
     file_gra_fa='./DUST_OPACITY/TRUST/ZDA_BARE_GR_S_SzDist_Gra.dat' ! Graphite
     file_sil_fa='./DUST_OPACITY/TRUST/ZDA_BARE_GR_S_SzDist_Sil.dat' ! Silicates
     file_pah_neu_fa='./DUST_OPACITY/TRUST/ZDA_BARE_GR_S_SzDist_PAH.dat' ! Neutral PAH molecules
  case('DraineLi06')
     file_gra_fa='./DUST_OPACITY/DraineLi06/DraineLi06_SzDist_Gra01.dat' ! Graphite
     file_sil_fa='./DUST_OPACITY/DraineLi06/DraineLi06_SzDist_Si01.dat' ! Silicates
     file_pah_neu_fa='./DUST_OPACITY/DraineLi06/DraineLi06_SzDist_PAHneu06.dat' ! Neutral PAH molecules
     file_pah_ion_fa='./DUST_OPACITY/DraineLi06/DraineLi06_SzDist_PAHion06.dat' ! Neutral PAH molecules
  case('user') ! input file names checked in check_input
     
  case default
     if (main_prc) print *, 'STOP(load_fa_arr): dust_model not recognized!'
     call stop_prc
  end select


allocate(n_dust_size(0:n_dust_comp-1))
n_dust_size = 0 
file_comp_arr=(/file_gra_fa, file_sil_fa,file_pah_neu_fa, file_pah_ion_fa/)

!!! determine number of sizes in the fa tables 

i = 0 
do id=0,max_n_dust_comp-1 

   if (.not. iq_dust_model(id)) cycle
   
   open(newunit=u,file=file_comp_arr(i),status='old')
   call count_lines(u, nlines)
   n_dust_size(i)=nlines-nhead
   close(u)

   i = i + 1
   
end do

!!! create fa_arr and read tables 

n_dust_maxsize_fa=maxval(n_dust_size)

allocate(dust_size_fa(0:n_dust_comp-1,0:n_dust_maxsize_fa-1),dust_fa(0:n_dust_comp-1,0:n_dust_maxsize_fa-1))
dust_size_fa=0
dust_fa=0

i = 0 
do id=0,n_dust_comp-1 

   if (.not. iq_dust_model(id)) cycle
   
   open(newunit=u,file=file_comp_arr(i),status='old')

   do j=0, nhead-1 ! skip header  
      read(u,*)
   end do

   do j=0, n_dust_size(i)-1
      read(u,*) dust_size_fa(i,j), dust_fa(i,j)  ! a [um]     f(a) [cm^-1 H^-1]   

     if (j > 0) then 
         if (dust_size_fa(i,j) < dust_size_fa(i,j-1)) then 
            if (main_prc) print *, 'STOP(load_fa_arr): grain size should be in ascending order in the input tables!'
            call stop_prc
         endif
      endif

  end do ! nsize

   close(u)

   i = i + 1
   
end do

!!! transform size array to SI system
dust_size_fa=dust_size_fa*1E-6 ! um -> m
dust_fa=dust_fa*1E2    ! cm^-1 H^-1 ->  m^-1 H^-1

call print_done

end subroutine load_fa_arr

!> Interpolate the qabs, qsca, qext and gsca arrays to the same grain sizes of the input grain size distribution tables
subroutine interpolate_q_grain_fa
  
  real(kind=real64), allocatable :: marr(:),carr(:),dust_size_fa_bin(:,:)
  real(kind=real64) :: m,c
  integer :: k,i,ic,ig0,ig1,imax,id

  if (main_prc) print *, 'interpolate single grain opacity coefficients...'
  
  !!! allocate qabs, qsca, qext and gsca arrays that will be used during the calculation
  allocate(qabs_arr_fa(0:n_dust_comp-1,0:n_dust_maxsize_fa-1,0:lnum_tot-1))
  allocate(qsca_arr_fa(0:n_dust_comp-1,0:n_dust_maxsize_fa-1,0:lnum_tot-1))
  allocate(qext_arr_fa(0:n_dust_comp-1,0:n_dust_maxsize_fa-1,0:lnum_tot-1))
  allocate(gsca_arr_fa(0:n_dust_comp-1,0:n_dust_maxsize_fa-1,0:lnum_tot-1))
  
  qabs_arr_fa=0
  qsca_arr_fa=0
  qext_arr_fa=0
  gsca_arr_fa=0

  !if (.not. use_lambda_grid) then
     allocate(qabs_ref_fa(0:n_dust_comp-1,0:n_dust_maxsize_fa-1))
     allocate(qsca_ref_fa(0:n_dust_comp-1,0:n_dust_maxsize_fa-1))
     allocate(qext_ref_fa(0:n_dust_comp-1,0:n_dust_maxsize_fa-1))
     allocate(gsca_ref_fa(0:n_dust_comp-1,0:n_dust_maxsize_fa-1))

     qabs_ref_fa=0
     qsca_ref_fa=0
     qext_ref_fa=0
     gsca_ref_fa=0

  !endif
  
  allocate(marr(0:lnum_tot-1),carr(0:lnum_tot-1))
  marr=0
  carr=0

 
! interpolate on grain sizes 

  ic = 0   
  do id=0,max_n_dust_comp-1

     if (.not. iq_dust_model(id)) cycle
     
     do i=0, n_dust_size(ic)-1 

        call value_locate(dust_size_fa(ic,i), dust_size_qabs(ic,0:n_dust_size_qabs(id)-1),ig0,ig1, .FALSE.)  ! note id index (not ic) for n_dust_size_qabs  

        marr=(qabs_arr_in(ic,ig1,:)-qabs_arr_in(ic,ig0,:))/(dust_size_qabs(ic,ig1)-dust_size_qabs(ic,ig0))
        carr=qabs_arr_in(ic,ig1,:)-marr*dust_size_qabs(ic,ig1)

        qabs_arr_fa(ic,i,:)=marr*dust_size_fa(ic,i)+carr  !! !the Qabs_arr_in array is evaluated at the same wavelenghts as the RT wavelength grid. However, it is not evaluated at the same grain sizes of the size distribution functions. Here qabs_arr_fa is calculated which contains values interpolated to those grain sizes.

        marr=(qsca_arr_in(ic,ig1,:)-qsca_arr_in(ic,ig0,:))/(dust_size_qabs(ic,ig1)-dust_size_qabs(ic,ig0))
        carr=qsca_arr_in(ic,ig1,:)-marr*dust_size_qabs(ic,ig1)

        qsca_arr_fa(ic,i,:)=marr*dust_size_fa(ic,i)+carr 


        marr=(qext_arr_in(ic,ig1,:)-qext_arr_in(ic,ig0,:))/(dust_size_qabs(ic,ig1)-dust_size_qabs(ic,ig0))
        carr=qext_arr_in(ic,ig1,:)-marr*dust_size_qabs(ic,ig1)

        qext_arr_fa(ic,i,:)=marr*dust_size_fa(ic,i)+carr 

        
        marr=(gsca_arr_in(ic,ig1,:)-gsca_arr_in(ic,ig0,:))/(dust_size_qabs(ic,ig1)-dust_size_qabs(ic,ig0))
        carr=gsca_arr_in(ic,ig1,:)-marr*dust_size_qabs(ic,ig1)

        gsca_arr_fa(ic,i,:)=marr*dust_size_fa(ic,i)+carr

        
      !  if (.not. use_lambda_grid) then

           m=(qabs_ref_in(ic,ig1)-qabs_ref_in(ic,ig0))/(dust_size_qabs(ic,ig1)-dust_size_qabs(ic,ig0))
           c=qabs_ref_in(ic,ig1)-m*dust_size_qabs(ic,ig1)

           qabs_ref_fa(ic,i)=m*dust_size_fa(ic,i)+c  !! !the Qabs_arr_in array is evaluated at the same wavelenghts as the RT wavelength grid. However, it is not evaluated at the same grain sizes of the size distribution functions. Here qabs_arr_fa is calculated which contains values interpolated to those grain sizes.

           m=(qsca_ref_in(ic,ig1)-qsca_ref_in(ic,ig0))/(dust_size_qabs(ic,ig1)-dust_size_qabs(ic,ig0))
           c=qsca_ref_in(ic,ig1)-m*dust_size_qabs(ic,ig1)

           qsca_ref_fa(ic,i)=m*dust_size_fa(ic,i)+c 

           
           m=(qext_ref_in(ic,ig1)-qext_ref_in(ic,ig0))/(dust_size_qabs(ic,ig1)-dust_size_qabs(ic,ig0))
           c=qext_ref_in(ic,ig1)-m*dust_size_qabs(ic,ig1)

           qext_ref_fa(ic,i)=m*dust_size_fa(ic,i)+c 

        
           m=(gsca_ref_in(ic,ig1)-gsca_ref_in(ic,ig0))/(dust_size_qabs(ic,ig1)-dust_size_qabs(ic,ig0))
           c=gsca_ref_in(ic,ig1)-m*dust_size_qabs(ic,ig1)

           gsca_ref_fa(ic,i)=m*dust_size_fa(ic,i)+c

       ! endif
        
      
        end do

        ic = ic + 1
        
     end do
     

!!!! create delta_dust_size

allocate(delta_dust_size(0:n_dust_comp-1, 0:n_dust_maxsize_fa-1),dust_size_fa_bin(0:n_dust_comp-1, 0:n_dust_maxsize_fa-2))
delta_dust_size=0
dust_size_fa_bin=0

do ic=0,n_dust_comp-1 

  imax=n_dust_size(ic) 

  !dust_size_fa_bin(ic,0:imax-2)=10_real64**((log10(dust_size_fa(ic,1:imax-1))+log10(dust_size_fa(ic,0:imax-2)))/2.)
  dust_size_fa_bin(ic,0:imax-2)=(dust_size_fa(ic,1:imax-1)+dust_size_fa(ic,0:imax-2))/2.
  
  delta_dust_size(ic,1:imax-2)=dust_size_fa_bin(ic,1:imax-2)-dust_size_fa_bin(ic,0:imax-3)
  delta_dust_size(ic,0)=dust_size_fa_bin(ic,0)-dust_size_fa(ic,0)
  delta_dust_size(ic,imax-1)=dust_size_fa(ic,imax-1)-dust_size_fa_bin(ic,imax-2)

end do

deallocate(marr,carr,dust_size_fa_bin)

call print_done

end subroutine interpolate_q_grain_fa

!> Loads specific enthalphy and specific heat capacity as a function of grain temperature from the input tables. In the standard mode, the tables are those provided in the TRUST RT benchmark project. Alternatively input tables can also be used, provided that they are in the same form as the TRUST tables. 
subroutine load_cT_hT_tables 
integer :: i, j, u,id 
integer :: nhead  
character(LEN=lcar) :: file_calorimetry_arr(0:max_n_dust_cal_type-1)
integer :: max_n_dust_temp_cal

! check if standard values should be used
  select case(dust_opacity_tables)
  case('TRUST', 'DraineLi06')
     file_calorimetry_Gra = './DUST_OPACITY/TRUST/Graphitic_Calorimetry_1000.dat'
     file_calorimetry_Sil = './DUST_OPACITY/TRUST/Silicate_Calorimetry_1000.dat'
     n_dust_temp_cal(0) = 1000 ! Number of tabulated temperatures for Graphite
     n_dust_temp_cal(1) = 1000 ! and for Silicates
  case('user')
  case default
     print *, 'STOP(load_cT_hT_tables): dust_opacity_tables not recognized!'
     stop
  end select

max_n_dust_temp_cal = maxval(n_dust_temp_cal)

! set indeces to call correct Ct table for each grain composition
allocate(iq_ct_table(0:n_dust_comp-1))
i = 0
do id=0,max_n_dust_comp-1 

   if (.not. iq_dust_model(id)) cycle
   select case(id)
   case(0,2,3)
      iq_ct_table(i) = 0
   case(1)
      iq_ct_table(i) = 1
   end select

   i = i + 1
   
end do

! allocate grain enthalphy arrays 
allocate(cal_temp(0:max_n_dust_temp_cal-1, 0:max_n_dust_cal_type-1))
allocate(grain_enthalpy(0:max_n_dust_temp_cal-1, 0:max_n_dust_cal_type-1))
allocate(grain_heat_capacity(0:max_n_dust_temp_cal-1, 0:max_n_dust_cal_type-1))
cal_temp = 0 
grain_enthalpy = 0
grain_heat_capacity = 0

! read files 

file_calorimetry_arr = (/file_calorimetry_Gra, file_calorimetry_Sil/)

nhead = 3
grain_density_arr = 0

do i= 0, max_n_dust_cal_type -1 

   ! if no graphite/PAH cycle
   if (i == 0 .and. .not. (iq_dust_model(0) .or. any(iq_dust_model(2:3)))) cycle
   ! if no Silicates cycle
   if (i == 1 .and. .not. (iq_dust_model(1))) cycle 

   open(newunit=u,file=file_calorimetry_arr(i),status='old')

   do j=0, nhead-1 ! skip header  
      read(u,*)
   end do

   read(u,*) grain_density_arr(i)  ! gm/cm^3
   
   do j = 0, n_dust_temp_cal(i)-1 
       ! K   erg/gm     erg/gm/K
      read(u,*) cal_temp(j,i), grain_enthalpy(j,i), grain_heat_capacity(j,i)
      !print *, cal_temp(j,i), grain_enthalpy(j,i), grain_heat_capacity(j,i)
      
   end do

   grain_enthalpy(:,i) = grain_enthalpy(:,i)*grain_density_arr(i) ! erg/cm^3
   grain_heat_capacity(:,i) = grain_heat_capacity(:,i)*grain_density_arr(i) !erg/cm^3/K
   
   close(u)

end do 

! change to units per volume and to SI units 
grain_enthalpy = grain_enthalpy*1E-1  ! erg/cm^3 -> J/m^3 
grain_heat_capacity = grain_heat_capacity*1E-1 ! erg/cm^3/K -> J/m^3/K


end subroutine load_cT_hT_tables


!> Finds indeces of an array corresponding to the array value bin containing the input value val(). Note that the array should normally be sorted in ascending order. In case val is lower than array(0), it returns the indeces of the first two elements.  If val is higher than array(size(array)-1), it returns indeces of last two elements. The input array can be in the reversing order if the reverse_order keyword is set TRUE.    
subroutine value_locate(val, array,ig0,ig1, reverse_order) 
  IMPLICIT NONE
  real(kind=real64) :: val, array(0:)
  integer :: ig0, ig1, num,i
  logical :: reverse_order

  num=size(array)
  
if (.not. reverse_order) then 

   if (val <= array(0)) then 
  
      ig0=0
      ig1=1
   
   else if (val >= array(num-1)) then

      ig0 = num-2
      ig1 = num-1

   else
   
      do i=0, num-1       
         if ((val > array(i)).and.(val <= array(i+1))) then 
            
            ig0=i
            ig1=i+1
            exit
         endif
      end do ! i 
        
   end if

else 

   if (val >= array(0)) then 
  
      ig0=0
      ig1=1
   
   else if (val <= array(num-1)) then

      ig0 = num-2
      ig1 = num-1

   else
   
      do i=0, num-1       
         if ((val <= array(i)).and.(val > array(i+1))) then 
            
            ig0=i
            ig1=i+1
            exit
         endif
      end do ! i 
        
   end if

endif

end subroutine value_locate


!> Interpolates linearly at a position x within the bin [x0,x1] given the boundary values y0=f(x0) and y1=f(x1). y0 and y1 are defined as arrays, so the subroutine can be used for multiple values provided that the x coordinates x0, x1 and x are the same.   
!> \todo You decided to use this routine only for single pairs of numbers, not arrays. You can remove the (0:) here as well as in the variables given in the input call. 
subroutine lin_interpolate(y0,y1,x0,x1,x,y)

  real(kind=real64) :: x0,x1,x, y0(0:),y1(0:), y(0:) 
 ! real(kind=real64), allocatable :: m(:), c(:)
  real(kind=real64) :: m(0:0), c(0:0)
  integer :: num

  num = size(y0)
  if (num /= size(y1) .or. num /= size(y)) then ! stop if not same sizes
     print *, 'STOP(lin_interpolate): sizes input y arrays are not the same!'
     print *, 'sizes =', size(y0), size(y1), size(y)
     stop
  end if

  if (x1-x0 == 0) then
     print *, 'STOP(lin_interpolate): input x values are equal!'
     print *, 'x0 =', x0
     print *, 'x1 =', x1
     stop
  endif

 ! allocate(m(0:num-1), c(0:num-1))

  m = (y1-y0)/(x1-x0)
  c = y1 - m*x1

  y = m*x + c 

 ! deallocate(m,c)

end subroutine lin_interpolate



!> Calculates integrated opacity coefficients kabs_arr(), ksca_arr(), kext_arr() and scattering phase function parameter gsca_arr() from the tabulated opacity parameters and the input grain size distributions. The resulting values for the integrated coefficients depends slightly on the specific interpolation scheme. For this reason, we add the option to upload these values from an input file as well using the keyword load_av_opacities() and file_int_opacities(). In this case, the code checks that the values obtained by this routine and those uploaded through the input file file_av_opacities() do not vary more than 5%. This is necessary to guarantee consistency between the dust emission and the stellar emission calculations within a few percents. 
!> \todo DONE_TO_TEST about lambda_ref when lambda_grids are used... make it such that lambda_ref has to coincide with one of the wavelengths of the lambda grid and scale dust density from there. this would allow the use of the same main grid for different models.    
subroutine calc_total_opacity
  integer :: i, ic, ig 
  integer :: imax

  if (main_prc) print *, 'calculate integrated/averaged opacity coefficients...'

  ! allocate opacity arrays for integrated properties 
  allocate(kext_arr(0:lnum_tot-1),kabs_arr(0:lnum_tot-1), ksca_arr(0:lnum_tot-1), gsca_arr(0:lnum_tot -1), ksca_arr_norm(0:lnum_tot-1))
  kext_arr = 0
  kabs_arr = 0
  ksca_arr = 0
  gsca_arr = 0
  tau_nh_ref = 0
  kext_ref = 0 
  tot_n_dust = 0 ! total number of grain particles 
  
  ! calculate total number of dust particles per H and
  ! integrate over species and size distribution at each wavelength. This gives tau / N_H = [m^2/H]

  do ic = 0, n_dust_comp-1
  
     imax = n_dust_size(ic)

     tot_n_dust = tot_n_dust+sum(delta_dust_size(ic,0:imax-1)*dust_fa(ic, 0:imax-1))
     
     do i = 0, lnum_tot-1
        
        kext_arr(i) = kext_arr(i) + pi*sum((dust_size_fa(ic,0:imax-1)**2)*qext_arr_fa(ic,0:imax-1,i)*delta_dust_size(ic,0:imax-1)*dust_fa(ic, 0:imax-1)) ! m^2/H
        
        kabs_arr(i) = kabs_arr(i) + pi*sum((dust_size_fa(ic,0:imax-1)**2)*qabs_arr_fa(ic,0:imax-1,i)*delta_dust_size(ic,0:imax-1)*dust_fa(ic, 0:imax-1))

        ksca_arr(i) = ksca_arr(i) + pi*sum((dust_size_fa(ic,0:imax-1)**2)*qsca_arr_fa(ic,0:imax-1,i)*delta_dust_size(ic,0:imax-1)*dust_fa(ic, 0:imax-1))

        gsca_arr(i) = gsca_arr(i) + pi*sum((dust_size_fa(ic,0:imax-1)**2)*qsca_arr_fa(ic,0:imax-1,i)*gsca_arr_fa(ic,0:imax-1,i)*delta_dust_size(ic,0:imax-1)*dust_fa(ic, 0:imax-1))
        
     end do

     tau_nh_ref = tau_nh_ref + pi*sum((dust_size_fa(ic,0:imax-1)**2)*qext_ref_fa(ic,0:imax-1)*delta_dust_size(ic,0:imax-1)*dust_fa(ic, 0:imax-1))
     
  end do

  ! convert to average values for gsca_arr and cross section per grain for k coefficients 
  gsca_arr = gsca_arr/ksca_arr
  kext_arr = kext_arr/tot_n_dust ! m^2
  kabs_arr= kabs_arr/tot_n_dust  ! m^2
  ksca_arr = ksca_arr/tot_n_dust ! m^2
  kext_ref = tau_nh_ref/tot_n_dust 

!!$  do i = 0, lnum_tot-1
!!$
!!$     print *, lambda_arr(i), kabs_arr(i)*1E4, ksca_arr(i)*1E4, gsca_arr(i)
!!$     
!!$  end do  
!!$  print *, lambda_ref, kext_ref, tau_nh_ref*1E4
  
  
call print_done
  
end subroutine calc_total_opacity

!> If keyword input_av_opacities is TRUE, it reads input table of integrated/averaged opacity coefficients kext(), kabs(), kext() and scattering phase function gsca. Interpolates factors to input wavelength grid lambda_arr() and to reference wavelength lambda_ref(). Checks values are within 5% of those derived directly by integration/average of grain opacities coefficients.
subroutine read_av_opacities
  real (Kind=real64), allocatable :: tlambda(:), tkabs(:), tksca(:), tkext(:), tgsca(:),ttau_nh(:)
  real(kind=real64), allocatable :: kext_arr_tmp(:), kabs_arr_tmp(:), ksca_arr_tmp(:), gsca_arr_tmp(:) 
  real(kind=real64) :: kext_ref_tmp, tau_nh_ref_tmp, diff(0:3)
  
  real (Kind=real64) :: a1,a2
  integer :: i, j, u, i0,i1, opnum, nlines, nhead,k 
  real (kind=real64) :: m, c 

  ! return if file_av_opacities does not have to be read
  
  if (.not. input_av_opacities) return

  if (main_prc) print *, 'reading integrated/averaged opacity coefficients from input file...'
  
  ! read input opacity file with integrated/averaged coefficients

  open(newunit=u,file=file_av_opacities, status='old')

  call count_lines(u, nlines)

  nhead = 4 
  opnum = nlines - nhead ! the first four lines are comments 
  
  allocate(tlambda(0:opnum-1), tkabs(0:opnum-1), tksca(0:opnum-1), tkext(0:opnum-1), tgsca(0:opnum-1),ttau_nh(0:opnum-1))

  do i=1,nhead    ! skip file header 
     read (u,*)
  end do

  do i=0,opnum-1

     ! um    cm^2/Ndust 
     read(u,*) tlambda(i), tkabs(i), tksca(i), ttau_nh(i),a2, tgsca(i)

     if (i > 0) then
        if (tlambda(i) < tlambda(i-1)) then  ! check lambda in ascending order !!!
           if (main_prc) print *, 'STOP: lambda not in ascending order in ', trim(file_av_opacities), '!'
           STOP
        endif
     endif
        
  enddo

  close(u)

  tkext=tkabs+tksca
  
  ! interpolate to lambda_arr and assign to temporary arrays 

  allocate(kext_arr_tmp(0:lnum_tot-1), kabs_arr_tmp(0:lnum_tot-1), ksca_arr_tmp(0:lnum_tot-1), gsca_arr_tmp(0:lnum_tot-1)) 
  
    k=0

   ! do i=0, n_dust_wave-1 
    do i = 0, opnum-1

       if (k > lnum_tot-1) exit 
       
       if (tlambda(i) < lambda_arr(k)) cycle
       
       m = (tkext(i) - tkext(i-1))/(tlambda(i)-tlambda(i-1))
       c = tkext(i) - m*tlambda(i)
       kext_arr_tmp(k) = m*lambda_arr(k)+c

       m = (tkabs(i) - tkabs(i-1))/(tlambda(i)-tlambda(i-1))
       c = tkabs(i) - m*tlambda(i)
       kabs_arr_tmp(k) = m*lambda_arr(k)+c
       
       m = (tksca(i) - tksca(i-1))/(tlambda(i)-tlambda(i-1))
       c = tksca(i) - m*tlambda(i)
       ksca_arr_tmp(k) = m*lambda_arr(k)+c

       m = (tgsca(i) - tgsca(i-1))/(tlambda(i)-tlambda(i-1))
       c = tgsca(i) - m*tlambda(i)
       gsca_arr_tmp(k) = m*lambda_arr(k)+c
       
       k=k+1

    end do

    ! interpolate to lambda_ref

    call value_locate(lambda_ref, tlambda, i0,i1, .FALSE.)
    m = (tkext(i1) - tkext(i0))/(tlambda(i1)-tlambda(i0))
    c = tkext(i1) - m*tlambda(i1)
    kext_ref_tmp = m*lambda_ref+c

 
    
    ! convert to mks units
    kext_arr_tmp = kext_arr_tmp*1E-4  ! [cm^2] - > [m^2]
    kabs_arr_tmp = kabs_arr_tmp*1E-4 
    ksca_arr_tmp = ksca_arr_tmp*1E-4
    
    kext_ref_tmp = kext_ref_tmp*1E-4
    tau_nh_ref_tmp = kext_ref_tmp*tot_n_dust ! note: this is not derived by interpolating input grid.
    
    
!!$    do i = 0, lnum_tot-1
!!$
!!$       print *, lambda_arr(i), kabs_arr_tmp(i), ksca_arr_tmp(i), kext_arr_tmp(i), gsca_arr_tmp(i)
!!$       print *, lambda_arr(i), kabs_arr(i), ksca_arr(i), kext_arr(i), gsca_arr(i)
!!$       read(*,*)
!!$    enddo

    ! check that values are consistent with those previously calculated within 5 %
    do i = 0, lnum_tot -1

       diff(0) = abs(kext_arr_tmp(i)-kext_arr(i))/(kext_arr_tmp(i)+kext_arr(i))*2

       diff(1) = abs(kabs_arr_tmp(i)-kabs_arr(i))/(kabs_arr_tmp(i)+kabs_arr(i))*2
       
       diff(2) = abs(ksca_arr_tmp(i)-ksca_arr(i))/(ksca_arr_tmp(i)+ksca_arr(i))*2
       
       diff(3) = abs(gsca_arr_tmp(i)-gsca_arr(i))/abs(gsca_arr_tmp(i)+gsca_arr(i))*2
       if (abs(gsca_arr(i)) < 1E-5 .and. gsca_arr_tmp(i) < 1E-5) diff(3) = 0 !ACHTUNG gsca not checked if very low

       if (any(diff(0:3) > 0.05) ) then    
          if (main_prc) then 
             print *, 'STOP: some values in ', trim(file_av_opacities), 'do not match values derived from the single grain opacities and the grain size distributions within 5%'
             print *, 'lambda = ', lambda_arr(i)
             print *, 'kext (file) =', kext_arr_tmp(i)
             print *, 'kext (calculated) =', kext_arr(i)
             print *, 'diff = ', diff(0)*100, '%'
             print *, ''
             print *, 'kabs (file) =', kabs_arr_tmp(i)
             print *, 'kabs (calculated) =', kabs_arr(i)
             print *, 'diff = ', diff(1)*100, '%'
             print *, ''
             print *, 'ksca (file) =', ksca_arr_tmp(i)
             print *, 'ksca (calculated) =', ksca_arr(i)
             print *, 'diff = ', diff(2)*100, '%'
             print *, ''
             print *, 'gsca (file) =', gsca_arr_tmp(i)
             print *, 'gsca (calculated) =', gsca_arr(i)
             print *, 'diff = ', diff(3)*100, '%'
             print *, ''             
             STOP
          endif
       endif
       
    end do

    
    ! check reference values are also consistent
    !kext_ref 
    diff(0) = abs(kext_ref_tmp-kext_ref)/kext_ref
    diff(1) = abs(tau_nh_ref_tmp-tau_nh_ref)/tau_nh_ref
    diff(2) = 0. ; diff(3) = 0 
    
    if (any(diff > 0.05)) then
       if (main_prc) then 
          print *, 'STOP: something wrong. Opacity reference values derived from input_av_opacities do not match those calculated previously from single grain opaticities and grain size distributions within 5%!' 
          print *, 'kext_ref (file) = ', kext_ref_tmp
          print *, 'kext_ref (calculated) =', kext_ref
          print *, ''
          print *, 'tau_nh_ref (file) = ', tau_nh_ref_tmp
          print *, 'tau_nh_ref (calculated) =', tau_nh_ref
          print *, ''
          STOP
       endif
    endif

    ! replace values with those derived from the input opacity file
    kext_arr = kext_arr_tmp
    kabs_arr = kabs_arr_tmp
    ksca_arr = ksca_arr_tmp
    gsca_arr = gsca_arr_tmp
    kext_ref = kext_ref_tmp
    tau_nh_ref = tau_nh_ref_tmp

deallocate(kext_arr_tmp, kabs_arr_tmp, ksca_arr_tmp, gsca_arr_tmp)

call print_done

end subroutine read_av_opacities



!> Sets lambda_arr_SI_HD_dust() values
!> \todo Should I really use this routine. More interpolations needed.... 
subroutine set_lambda_arr_SI_HD_dust 
  integer, parameter :: num = 100 ! number of wavelengths 
  real(kind=real64) :: l0, l1 ! minimum and maximum wavelength [m]
  real(kind=real64) :: step_size
  integer :: i 
  
  allocate(lambda_arr_SI_HD_dust(0:num-1), lambda_arr_SI_HD_dust_bin(0:num-2), delta_lambda_bin_HD_dust(0:num-1)) ! all in [m]

  l0 = 1E-6  ! 1 um 
  l1 = 1000E-6  ! 1000 um 

  step_size = (log10(l1)-log10(l0))/(num-1)

  do i = 0, num-1 

    lambda_arr_SI_HD_dust(i) = 10._real64**(i*step_size + log10(l0))
    
  end do

  lambda_arr_SI_HD_dust_bin=10_real64**((log10(lambda_arr_SI_HD_dust(1:num-1))+log10(lambda_arr_SI_HD_dust(0:num-2)))/2.)
  
  delta_lambda_bin_HD_dust(1:num-2) = lambda_arr_SI_HD_dust_bin(1:num-2)-lambda_arr_SI_HD_dust_bin(0:num-3)
  delta_lambda_bin_HD_dust(0)=lambda_arr_SI_HD_dust_bin(0)-lambda_arr_SI_HD_dust(0)
  delta_lambda_bin_HD_dust(num-1)=lambda_arr_SI_HD_dust(num-1)-lambda_arr_SI_HD_dust_bin(num-2)

  deallocate(lambda_arr_SI_HD_dust_bin) ! this array is not needed anymore 

end subroutine set_lambda_arr_SI_HD_dust

!> Calculates planck averaged Qabs for the grains in the input grain size distributions and for a set of temperatures.
subroutine  calc_planck_av_qabs
  real(kind=real64),parameter :: t_min = 1, t_max = 2500
  real(kind=real64) :: step_size
  integer :: i, ic, ig, it  
  real(kind=real64) :: num, den
  real(kind=real64), allocatable :: bplanck_arr(:,:)

  if (main_prc) print *, 'calculating Planck averaged Qabs coefficients...'

  allocate(qabs_arr_planck(0:n_dust_comp-1,0:n_dust_maxsize_fa-1,0:n_temp_planck-1), T_arr_planck(0:n_temp_planck-1))
  allocate(bplanck_arr(0:n_temp_planck-1,0:lnum_tot-1))
  qabs_arr_planck = 0
  T_arr_planck = 0
  bplanck_arr = 0
  
  step_size = (log10(t_max)-log10(t_min))/(n_temp_planck-1)

  ! derive temperatures 
  do i = 0, n_temp_planck-1
     t_arr_planck(i) = 10._real64**(i*step_size+log10(t_min))
  end do

  ! evaluate planck function at the wavelength grid point 
  do it = 0, n_temp_planck-1 
     do i= 0, lnum_tot-1               
        bplanck_arr(it,i) = bplanck(t_arr_planck(it), lambda_arr_SI(i))
     end do
  end do

  ! derive planck averaged Qabs coefficients 
  do ic=0,n_dust_comp-1   
     do ig=0, n_dust_size(ic)-1    
        do it = 0, n_temp_planck-1           

           num = sum(qabs_arr_fa(ic,ig,:)*bplanck_arr(it,:)*delta_lambda_bin)
           den = sum(bplanck_arr(it,:)*delta_lambda_bin)           
           qabs_arr_planck(ic,ig,it) = num/den
           
        end do
     end do
  end do

  call print_done 
  
end subroutine calc_planck_av_qabs
  
!> Calculates energy of photons as a function of wavelength for the set of wavelength stored in lambda_arr_SI. The energy is in [J].
subroutine calc_photon_energy
 
allocate(phot_energy(0:lnum_tot-1))

phot_energy = hplanck*cspeed/lambda_arr_SI

end subroutine calc_photon_energy

!> Calculates pt() in the large grain energy approximation (see section 4.4 Voit 1991). In case 1) the mean photon energy is very small compared to the grain equilibrium energy and 2) the width of the Gaussian approximation is more than 10% the equilibrium dust temperature, the code assign large_grain_energy = .FALSE. So the stochastically heated emission is calculated numerically from the current grain size to smaller sizes. 
subroutine calc_gaussian_fE(t_dust, ic,ig) 
!> @param T_dust Equilibrium dust temperature 
real(kind=real64) :: t_dust
!> @param ct heat capacity C_T at equilibrium dust temperature 
real(kind=real64) :: ct  
!> @param Em Grain enthalpy at equilibrium dust temperature 
real(kind=real64) :: Em  
!> @param Qp Planck-averaged Qabs() at equilibrium dust temperature.
real(kind=real64) :: Qp 
!> @param dqp_dt dQp/dT at equilibrium temperature. Used to calculate dEdot/dE
real(kind=real64) :: dqp_dt 
!> @param dEdot_dE dEdot/dE at equilibrium temperature. This is the derivative of the cooling rate with respect to the grain enthalpy.
real(kind=real64) :: dEdot_dE
!> @param sig Sigma of the gaussian derived for f(E) (the grain enthalpy probability distribution).
real(kind=real64) :: sig
!> @param sigT Sigma of the gaussian derived for f(T) (the grain temperature probability distribution).
real(kind=real64) :: sigT
!> @param vol Grain volume
real(kind=real64) :: vol 

integer :: i0, i1, ic,ig   
real(kind=real64) :: x0,x1,x
real(kind=real64) :: y0(0:0), y1(0:0), y(0:0)


! interpolate grain heat capacity and grain enthalpy at equilibrium temperature
call value_locate(t_dust, cal_temp(:, iq_ct_table(ic)), i0,i1, .FALSE.)  
x0 = cal_temp(i0, iq_ct_table(ic))
x1 = cal_temp(i1, iq_ct_table(ic))
y0 = grain_heat_capacity(i0, iq_ct_table(ic))
y1= grain_heat_capacity(i1, iq_ct_table(ic))
x = t_dust
call lin_interpolate(y0,y1,x0,x1,t_dust,y)
ct = y(0) ! heat capacity
         
y0 = grain_enthalpy(i0, iq_ct_table(ic))
y1=  grain_enthalpy(i1, iq_ct_table(ic))
call lin_interpolate(y0,y1,x0,x1,t_dust,y)
Em = y(0)*4./3.*pi*dust_size_fa(ic,ig)**3.  ! grain enthalpy
         
! interpolate planck averaged emissivity at equilibrium temperature
call value_locate(t_dust, t_arr_planck , i0,i1, .FALSE.) 
x0 = t_arr_planck(i0)
x1 = t_arr_planck(i1)
y0 = qabs_arr_planck(ic,ig,i0)
y1 = qabs_arr_planck(ic,ig,i1)
call lin_interpolate(y0,y1,x0,x1,t_dust,y)
qp = y(0) ! planck averaged emissivity      
         
! calculate dQabs/dT          
dqp_dt = (y1(0) - y0(0))/(x1-x0)
         
! calculate dEdot_dE
!dEdot_dE = 4*pi*dust_size_fa(ic,ig)**2*kboltz*t_dust**3*(4*qp+t_dust*dqp_dt) ! original formula. below without some constant factors which are removed later
dEdot_dE = 4*sigmaSB*t_dust**3*(4*qp+t_dust*dqp_dt)   
         
! calculate sigma equilibrium distribution for energy 
sig = 0.5 * Rd_arr(2) * dust_size_fa(ic,ig)*ct* 4./3./dEdot_dE 
sig = sqrt(sig)
         
! calculate sigma equilibrium distribution for temperature
vol = 4./3.*pi*dust_size_fa(ic,ig)**3
sigT = sig/ct/vol
         
! decide if using Gaussian analytical solution for large grain energies or starting numerical calculations of stochastic heating following Guhathakurta & Draine 1989 plus the analytical solutions of Voit 1991 (see section 4.5 
if (rd_arr(1)/rd_arr(0)/Em < 0.01 .and. 2*sigT/t_dust < 0.1) then
            
   ! define energy grid around equilibrium energy
   Emin = Em -8*sig   
   if (Emin < 0.01*Em) Emin = 0.01*Em
   Emax = Em + 8*sig
                        
   call make_log_array(Emin, Emax, E_arr)

   call make_delta_array(E_arr, delta_E_arr_bin)

   ! set energy probability distribution 
   pt = exp(-0.5*((E_arr-Em)/sig)**2.)  ! Gaussian approximation (not normalized yet)
            
   pt = pt*delta_E_arr_bin/sum(pt*delta_E_arr_bin)  ! Note: this is f(E)dE = P(T)dT
            
   ! set temperature array corresponding to energy array             
   call convert_E_arr_to_T_arr(E_arr, T_arr, ic,ig)        
            
else

   large_grain_energy = .FALSE.
   Tmin = t_dust - 5*sigT   ! this is the initial temperature interval where to look for in the stochastic heating algorithm below. Some iterations might be needed if range too small.  
   if (Tmin < t_dust/2.) Tmin = T_dust/2.
   Tmax = t_dust + 5*sigT
                    
endif


end subroutine calc_gaussian_fE

!> Calculates the temperature probability distribution pt() using the numerical method of Gahathakurta & Draine 1989 with the modification of Voit 1991 (section 4.5). 
subroutine calc_full_fE(ic,ig)
integer ::i,  ic, ig,is, i0
real(kind=real64) :: a1(0:n_temp_pt-1), b1(0:n_temp_pt-1), c1(0:n_temp_pt-1), l1(0:n_temp_pt-1), l2(0:n_temp_pt-1)
!> @param SeE source function defined in Equ. 53 of Voit 1991
real(kind=real64) :: SeE
!> @param n_temp_sub Number of points for the integration of Fe1 (see Equ. 50 voit 1991).  
integer :: n_temp_sub 
!> @param E_arr_sub Energy array for the sampling of the energy bin [E-e, E]
real(kind=real64), allocatable :: E_arr_sub(:)
!> @param pt_sub Temperature probability distribution on the grid defined by E_arr_sub
real(kind=real64), allocatable :: pt_sub(:)
!> @param cc1 Argument of first exponential in the analytical solution for pt_sub().
real(kind=real64), allocatable :: cc1(:)
!> @param cc2 Argument of second exponential in the analytical solution for pt_sub().
real(kind=real64), allocatable :: cc2(:)
!> @param delta_E_arr_sub Bin sizes for E_arr_sub. 
real(kind=real64), allocatable :: delta_E_arr_sub(:)
!> @param f0 first boundary condition for solving the pt_sub equation 
real(kind=real64) :: f0
!> @param f0 second boundary condition for solving the pt_sub() equation (see Equ. 54 Voit 1991) 
real(kind=real64) :: ff
!> @param c2 coefficient differential equation pt_sub().
real(kind=real64) :: c2
!> @param d2 coefficient differential equation pt_sub().
real(kind=real64) :: d2
!> @param a2 coefficient differential equation pt_sub().
real(kind=real64) :: a2
!> @param b2 coefficient differential equation pt_sub().
real(kind=real64) :: b2
!> @param rd_interpol Interpolated values of rd_integrated().
real(kind=real64), allocatable :: rd_interpol(:)
!> @param percT percentage factor used to vary temperature range. 
real(kind=real32), parameter :: percT = 0.3  

n_temp_sub = 50

a1=0
b1=0
c1=0
l1=0
l2=0
Re0 = 0
Re1 = 0
Re2 = 0

allocate(E_arr_sub(0:n_temp_sub-1), pt_sub(0:n_temp_sub-1), cc1(0:n_temp_sub-1), cc2(0:n_temp_sub-1), delta_E_arr_sub(0:n_temp_sub-1),rd_interpol(0:n_temp_sub-1) )
E_arr_sub = 0
pt_sub = 0
cc1 =0
cc2 = 0 
delta_E_arr_sub =0 
rd_interpol =0 

!!$i0 = i_lambda_dust(0)
!!$
!!$abs_int_rad_stars=pi*dust_size_fa(ic,ig)**2*abs_int_rad_stars  ! 1/(m s) absorbed stellar emission photon rate per unit wavelength interval
!!$
!!$if (iterations_dustem > 1) then 
!!$   abs_int_rad_dust=pi*dust_size_fa(ic,ig)**2*abs_int_rad_dust ! 1/(m s)  ! absorbed dust emission photon rate per unit wavelength interval
!!$else 
!!$   abs_int_rad_dust = 0 
!!$endif

! Calculate arrays of integrals of photon absorption rate 
call calc_integrals_photon_abs_rate(dust_size_fa(ic,ig))

do ! iterate on temperature range if necessary  
            
   if (tmin /= tmin_prev .or. tmax /= tmax_prev) then 

      call make_log_array(Tmin, Tmax, T_arr)  ! temperature array in the range [Tmin, Tmax] 

   endif

   call convert_T_arr_to_E_arr(T_arr, E_arr, ic,ig)    ! corresponding Enthalpy array
         
   call make_delta_array(E_arr, delta_E_arr_bin)    ! bin size array for integration over enthalpy
   call interpolate_qabs_arr_planck(T_arr, qp_arr, ic,ig) 

   

   ! Calculate transition matrices AA and BB
   call calc_transition_matrices

   ! Calculate cooling rate at the temperatures in Tarr and assign probability rate for cooling to next lower level 
   call calc_Edot_arr(dust_size_fa(ic,ig), T_arr) 
         
   ! Calculate dosage function Rd0 and its moments Re0, Re1 and Re2
   call calc_dosage_function_moment_integrals
           
   ! Calculate constants for analytic solution in the last integration bin [E-e, E] (see Equ. 50 Voit 1991)
   do i = 1, n_temp_pt-2

      if (Re2(i)/rd_arr(2) > 1E-6) then
         a1(i) =  (Edot_arr(i+1)-Re1(i))/(Re2(i)/2.)  ! coefficient of df/dE
         b1(i) = -(Rd_arr(0)-Re0(i))/(Re2(i)/2.)      ! coefficient of f 
         c1(i) =  2./Re2(i)                           ! coefficient of SeE
         l1(i) = -0.5*a1(i)+SQRT(0.25*a1(i)**2-b1(i))  ! exponential coefficients exp(l*E)
         l2(i) = -0.5*a1(i)-SQRT(0.25*a1(i)**2-b1(i))

      else
         a1(i) = Edot_arr(i+1)/rd_arr(0) 
         l1(i) = Rd_arr(0)/Edot_arr(i+1)
      endif
    
   end do

  
        
   ! Loop on energy levels to calculate P values 

   pt = 0
   is = 1
   do ! make sure startin bb(is-1, is) > 0 
      if (bb(is-1,is) > 0) exit
      is = is+1
      
   end do

   pt(is-1) = 1E-10  ! initial arbitrary value for first element
        
   do i = is, n_temp_pt-2

      if (pt(i-1) < 0) pt(i-1) = 0
      ! FIRST APPROXIMATION (Guhatakurta & Draine 1989 equivalent to Equ 49 Voit 1991)
      pt(i) = - sum(BB(i-1,0:i-1)*Pt(0:i-1)) / BB(i-1,i) ! probability within temperature/enthalpy bin

      ! analytical approximation (Voit 1991)

      ! source function Se(E)
      SeE = sum(Pt(0:i-1)*AA(i-1,0:i-1))  ! remember that Pt is the total probability within the temperature/enthalpy bin. This is why there is no delta_T or delta_E term within the integration. Compared to Equ. 53, the sum is performed from the highest energies to the lowest (results are the same, only order of the integral is different).     

      ! create subgrid between E and E+dE for integrating Fe1(E) (equ 50 Voit 1991)
      Emin = E_arr(i-1)
      Emax = E_arr(i)
          
      call make_linear_array(Emin, Emax, E_arr_sub)

      if (Re2(i)/rd_arr(2) > 1E-6) then

         ff = 2./Re2(i)*(pt(i)*BB(i-1,i)+(Emax - Emin)*SeE+(Re1(i)-Edot_arr(i))*pt(i-1)/(Emax-Emin))  ! second boundary condition, equ. 54 Voit 1991
      
         C2 = -c1(i)/(b1(i))*SeE 
         f0 = pt(i-1)/(Emax-Emin)
         D2 = l2(i)-l1(i)
               
         A2 = 0.
         B2 = 0.
         IF (D2 /= 0) THEN  ! in this case A2 = B2 
            A2 = (l2(i)*(f0-C2)-ff)/D2   
            B2 = (l1(i)*(f0-C2)-ff)/D2
         ENDIF
              
         cc1 = (l1(i)*(E_arr_sub-E_arr(i-1)))
         cc2 = (l2(i)*(E_arr_sub-E_arr(i-1)))
               
         pt_sub = (A2*EXP(cc1)-B2*EXP(cc2))+C2  ! solution Eq. 52 Voit 1991
               
      else 

         f0 = pt(i-1)/(E_arr(i)-E_arr(i-1))
         pt_sub = (f0-SeE/Rd_arr(0)) * EXP(l1(i)*(E_arr_sub-E_arr(i-1)))+SeE/Rd_arr(0) ! solution Eq. 52 Voit 1991 when all Ren terms are negligible. 

      endif
            

      ! define energy array with transition from each energy level to last energy level
           
      !delta_E_arr_sub = (Emax+delta_E_arr_sub_bin(n_temp_sub-1)*2) - E_arr_sub   ! OLD ONE 
      delta_E_arr_sub = E_arr(i+1) - E_arr_sub
             
      ! interpolate Rd_integrated array at these energies      
      call interpolate_Rd_integrated(n_temp_sub-1, delta_E_arr_sub, rd_interpol)
      
                  
      ! integrate to calculate Fe1 and add contribution to pt array
      pt(i) = pt(i)  + sum(rd_interpol(1:n_temp_sub-1)*pt_sub(1:n_temp_sub-1)*(E_arr_sub(1:n_temp_sub-1)-E_arr_sub(0:n_temp_sub-2)))/bb(i-1,i) ! bb(i-1,i) is here because is equal to Edot/delta_E. Edot is in formula 50 Voit 1991. delta_E is also needed so we get the probability contribution to the entire energy bin of this term 
         
   end do

   where(pt < 0)   ! there can be very small negative values at the tail of the distribution. 
      pt =0
   end where

   pt = pt/sum(pt)        
         
   tmin_prev = tmin
   tmax_prev = tmax 

   ! check probability at the boundaries of temperature range
   if (pt(0) > 1E-20) then
      tmin = tmin - percT*tmin
   endif

   if (pt(n_temp_pt-2) > 1E-20) then
      tmax = tmax + percT*tmax
   endif

!!$   print *, 'first tries ' 
!!$     print *, 'p(T)   T'  ! slightly different from IDL code but check effect of the different wavelength grid. 
!!$     do i = 0, n_temp_pt-1 
!!$        print *, pt(i), t_arr(i)
!!$     enddo
!!$     read(*,*)


   ! exit if temperature limits are unchanged
   if (tmin == tmin_prev .and. tmax == tmax_prev) exit 
         
end do

deallocate(E_arr_sub, pt_sub, cc1, cc2, delta_E_arr_sub, rd_interpol)
 
if (count(pt /= pt) > 0) then ! if there are NANs set to zero. This can happen for few grains in some weird cases. To investigate origin. For now, just set to zero. 
   print *, 'WARNING: NAN in PT dust. Set to 0!'
   pt = 0
endif

!imax_arr = maxloc(pt)-1 
!print *, 'max T= ', t_arr(imax_arr(0))
!read(*,*)

end subroutine calc_full_fE

!> Makes array of values between xmin and xmax, equally spaced in logarithmic space. 
subroutine make_log_array(xmin, xmax, xarr)
real(kind=real64) :: xmin, xmax, xarr(0:)

integer :: num, i 
real(kind=real64) :: step_size 

! check input array size 
num = size(xarr)

! define step size and assign array values 
step_size = (log10(xmax) - log10(xmin))/(num-1)

do i = 0, num-1 
   xarr(i) = 10._real64**(i*step_size + log10(xmin))
end do

end subroutine make_log_array

!> Makes array of bin sizes that can be used for integration. 
!> \todo add check that the array is in ascending order ? 
subroutine make_delta_array(xarr, delta_xarr_bin)
  real(kind=real64) :: xarr(0:), delta_xarr_bin(0:)
  integer :: num 
  real(kind=real64), allocatable :: xarr_bin(:)

! check input array size 
num = size(xarr)
if (size(delta_xarr_bin) /= num ) then 
   print *, 'STOP(make_delta_array):: input delta_xarr_bin array does not have the right size'
   stop
endif

! define array bin size

allocate(xarr_bin(0:num-2))

xarr_bin = 10_real64**((log10(xarr(1:num-1))+log10(xarr(0:num-2)))/2.)
delta_xarr_bin(1:num-2)=xarr_bin(1:num-2)-xarr_bin(0:num-3)
delta_xarr_bin(0)=xarr_bin(0)-xarr(0)
delta_xarr_bin(num-1)=xarr(num-1)-xarr_bin(num-2)

deallocate(xarr_bin)

end subroutine make_delta_array

!> Converts grain enthalpy array into temperature array. 
subroutine convert_E_arr_to_T_arr(E_arr, T_arr, ic,ig) 
integer :: i, num , ic, ig 
real(kind=real64) :: E_arr(0:), T_arr(0:), vol
integer :: i0, i1 
real(kind=real64) :: x0, x1, y0(0:0), y1(0:0),x, y(0:0)

! check array size  
num =size(E_arr)
if (size(T_arr) /= num) then 
   print *, 'STOP(convert_E_arr_to_T_arr):: input T_arr array does not have the right size'
   stop
endif 

vol = 4./3.*pi*dust_size_fa(ic,ig)**3

do i = 0, num-1  ! this can be expensive and useless when you are close to equilibrium and all values are within the same energy bin. 

   call value_locate(E_arr(i), grain_enthalpy(:, iq_ct_table(ic))*vol, i0,i1, .FALSE.)
   x0 = grain_enthalpy(i0, iq_ct_table(ic))*vol
   x1=  grain_enthalpy(i1, iq_ct_table(ic))*vol
   y0 = cal_temp(i0, iq_ct_table(ic))
   y1 = cal_temp(i1, iq_ct_table(ic))
   x = E_arr(i)
   call lin_interpolate(y0,y1,x0,x1,x,y)
   T_arr(i) = y(0)   ! this is needed when calculating grain emission in the range of temperatures

end do

end subroutine convert_E_arr_to_T_arr


!> Interpolates Planck averaged emissivity table at the set of input temperatures.
subroutine interpolate_qabs_arr_planck(T_arr, qp_arr, ic,ig) 
real(kind=real64) :: T_arr(0:), qp_arr(0:)
integer :: ic, ig, i, num, i0, i1  
real(kind=real64) :: x0, x1, x, y0(0:0), y1(0:0), y(0:0)

! check array size  
num =size(T_arr)
if (size(qp_arr) /= num) then 
   print *, 'STOP(intepolate_qabs_arr_planck):: input qp_arr array does not have the right size'
   stop
endif 

do i = 0, num-1
   call value_locate(t_arr(i), t_arr_planck , i0,i1, .FALSE.) 
   x0 = t_arr_planck(i0)
   x1 = t_arr_planck(i1)
   y0 = qabs_arr_planck(ic,ig,i0)
   y1 = qabs_arr_planck(ic,ig,i1)
   x = t_arr(i)
   call lin_interpolate(y0,y1,x0,x1,x,y)
   qp_arr(i) = y(0)

end do

end subroutine interpolate_qabs_arr_planck

!> Calculates integrals of the photon absorbtion rate from a minimum lambda up to all the other possible values. 
subroutine calc_integrals_photon_abs_rate(dust_size)
integer :: i0, i  
real(kind=real64) :: dust_size
real(kind=real64) :: dG

i0 = i_lambda_dust(0)

abs_int_rad_stars=pi*dust_size**2*abs_int_rad_stars  ! 1/(m s) absorbed stellar emission photon rate per unit wavelength interval

if (iterations_dustem > 0) then ! see comment for iterations_dustem in calc_t_dust_equil.   
   abs_int_rad_dust=pi*dust_size**2*abs_int_rad_dust ! 1/(m s)  ! absorbed dust emission photon rate per unit wavelength interval
else 
   abs_int_rad_dust = 0 
endif

Rd_integrated(0) = abs_int_rad_stars(0)*delta_lambda_bin_stars(0)
if (i0 == 0) then 
   Rd_integrated(0) = Rd_integrated(0)+ abs_int_rad_dust(0)*delta_lambda_bin_dust(0) 
endif

do i = 1, lnum_tot -1 

   dG = 0 

   if (i <= lnum_stars-1) then 
      dG = dG + abs_int_rad_stars(i)*delta_lambda_bin_stars(i) 
   endif

   if (i >= i0) then 
      dG = dG + abs_int_rad_dust(i-i0)*delta_lambda_bin_dust(i-i0)
   endif

   Rd_integrated(i) = Rd_integrated(i-1) + dG
   
end do



end subroutine calc_integrals_photon_abs_rate

!> Calculates transitions matrices AA() and BB(). 
subroutine calc_transition_matrices
integer :: i 
real(kind=real64) :: delta_E_arr(0:n_temp_pt-1), Rd_Interpol(0:n_temp_pt-1)

! Rd_integrated is defined on the energy grid phot_energy (descending order). 
! Need to interpolate to get the values of the absorbed photon rate at the energy defined in delta_E_arr (which corresponds to the transitions between the levels).

bb = 0
aa = 0

do i = 1, n_temp_pt-2

   delta_E_arr(0:i) = E_arr(i+1) - E_arr(0:i)  ! array of energy differences between level i+1 and all lower levels. 

   call interpolate_Rd_integrated(i, delta_E_arr, rd_interpol)

   aa(i-1, 0:i-1) = (rd_interpol(1:i)-rd_interpol(0:i-1))/(E_arr(1:i)-E_arr(0:i-1))  ! using Earr difference instead of delta_E_arr_bin because the latter is for integration only. This is the dosage function evaluated at the energy levels up to i. This term is needed when calculating the source term SeE (see equ 53 in Voit 53). Note the order of the rd_interpol terms due to the fact that rd_interpol is in ascending order. The energy difference at the denominator is the one that has to be divided to pt later in order to perform the integral on Ed. 

   bb(i-1, 0:i-1) = -rd_interpol(1:i)  ! the minus sign here is correct. formulas are slightly changed when calculating pt and ff to take this sign into account
   

end do

end subroutine calc_transition_matrices

!> Calculates cooling rate at all temperatures in Tarr. It also calculates continuous cooling term BB(i,i+1).
subroutine calc_Edot_arr(dust_size, T_arr)
real(kind=real64) :: dust_size
real(kind=real64) :: T_arr(0:)
integer :: i 

Edot_arr = 4*pi*dust_size**2*sigmaSB*qp_arr*T_arr**4

do i = 0, n_temp_pt-2 

   bb(i, i+1) = Edot_arr(i+1)/(E_arr(i+1)-E_arr(i))
   
end do

end subroutine calc_Edot_arr


!> Calculates dosage function by differentiating the integrals of the photon absorption rate. These integrals are stored in the Rd_integrated array for the energies corresponding to the wavelength grid. The dosage function has to be calculated at the energies corresponding to the enthalpy difference for the levels defined by the grain temperature grid. So, first an interpolation of the Rd_integrated array is necessary. Then differentiation. Then calculation of moments of dosage function is performed up to all energies corresponding to the energy bins.  
subroutine calc_dosage_function_moment_integrals
  integer :: i
  integer, parameter :: n_delta_e = 100
  real(kind=real64) :: rd_interpol_e(0:n_delta_e-1), delta_E_arr_e(0:n_delta_e-1)
  real(kind=real64) :: max_delta_E_arr
  real(kind=real64) :: Re_e(0:n_delta_e-1), dosage_function_e(0:n_delta_e-1) 
  real(kind=real64) :: delta_E_arr(0:n_temp_pt-1)

delta_E_arr = 0  ! transition bin energy for the original energy grid
delta_E_arr_e = 0 ! transition bin energy for smoother grid from 0 to max delta_E_arr  
rd_interpol_e = 0
Re0 = 0
Re1= 0
Re2 = 0 
Re_e = 0

! define delta energy array (this is the epsilon in Voit 1991). 
delta_E_arr(0) = 0 
delta_E_arr(1:n_temp_pt-1) = (E_arr(1:n_temp_pt-1) - E_arr(0:n_temp_pt-2)) 

! ACHTUNG: TO derive the dosage function accurately within the range of energies covered by the energy bins, you need finer interpolation first. Then integrate to get the moments accurately as a function of bin energy (epsilon) and finally interpolate the moments at the exact values of the energy bins.... 
max_delta_E_arr = maxval(delta_E_arr)
call make_log_array(0.01_real64*1E-19, max_delta_E_arr, delta_E_arr_e)
delta_E_arr_e(0) = 0 

! interpolate Rd_integrated array at the values defined in delta_E_arr
call interpolate_Rd_integrated(n_delta_e-1, delta_E_arr_e, rd_interpol_e)

! calculate dosage function through difference (note here I do not divide by an energy bin size because this is then multiplied later when calculating the moments). Note the order of the rd_interpol terms due to the fact that rd_interpol is in DESCENDING order (not like in calc_transition_matrices). 

dosage_function_e(1:n_delta_e-1) = (rd_interpol_e(0:n_delta_e-2)-rd_interpol_e(1:n_delta_e-1))  

Re_e(0) = 0

! calculate zero moment dosage function up to all delta_E_arr values
do i = 1, n_delta_e-1

      ! Equ. 51 Voit 1991 
   Re_e(i) = Re_e(i-1) + dosage_function_e(i)
   
end do

call interpolate_array(n_temp_pt-1, delta_E_arr_e, Re_e, delta_E_arr, Re0)

Re_e = 0
! calculate first moment dosage function
do i = 1, n_delta_e-1
   
   Re_e(i) = Re_e(i-1) +  dosage_function_e(i)*delta_E_arr_e(i-1)
   
end do

call interpolate_array(n_temp_pt-1, delta_E_arr_e, Re_e, delta_E_arr, Re1)

Re_e = 0 

! calculate second moment dosage function
do i = 1, n_delta_e-1

   Re_e(i) = Re_e(i-1) +  dosage_function_e(i)*delta_E_arr_e(i-1)**2
   
end do

call interpolate_array(n_temp_pt-1, delta_E_arr_e, Re_e, delta_E_arr, Re2)

end subroutine calc_dosage_function_moment_integrals

!> Performs linear interpolation for all the values in the arrays xarr_out and yarr_out given the values in the arrays xarr_in and yarr_in. 
subroutine interpolate_array(num, xarr_in, yarr_in, xarr_out, yarr_out)
integer :: i, num, i0,i1 
real(kind=real64) :: x0, x1, x, y0(0:0), y1(0:0), y(0:0)
real(kind=real64) :: xarr_in(0:), yarr_in(0:), xarr_out(0:), yarr_out(0:)

do i = 0, num 
   
   call value_locate(xarr_out(i), xarr_in, i0,i1, .FALSE.) 
   x0 = xarr_in(i0)
   x1 = xarr_in(i1)
   y0 = yarr_in(i0)
   y1 = yarr_in(i1)
   x = xarr_out(i)
   call lin_interpolate(y0,y1,x0,x1,x,y)   
   yarr_out(i) = y(0)   

end do

end subroutine interpolate_array

!> Interpolates Rd_integrated() array at the energies corresponding to the transitions between grain enthalphy levels 
subroutine interpolate_Rd_integrated(num, delta_E_arr_in,rd_interpol_in)
integer :: i, num, i0,i1 
real(kind=real64) :: x0, x1, x, y0(0:0), y1(0:0), y(0:0), delta_E_arr_in(0:), rd_interpol_in(0:)

do i = 0, num
   
   call value_locate(delta_E_arr_in(i), phot_energy, i0,i1, .TRUE.) ! note the reverse order
   x0 = phot_energy(i0)
   x1 = phot_energy(i1)
   y0 = Rd_integrated(i0)
   y1 = Rd_integrated(i1)
   x = delta_E_arr_in(i)
   call lin_interpolate(y0,y1,x0,x1,x,y)   
   rd_interpol_in(i) = y(0)   

end do

end subroutine interpolate_Rd_integrated

!> Makes array of values between xmin and xmax, equally spaced in linear space. 
subroutine make_linear_array(xmin, xmax, xarr)
real(kind=real64) :: xmin, xmax, xarr(0:)

integer :: num, i 
real(kind=real64) :: step_size 

! check input array size 
num = size(xarr)

! define step size and assign array values 
!step_size = (log10(xmax) - log10(xmin))/(num-1)
step_size = (xmax - xmin)/(num-1)

do i = 0, num-1 
  ! xarr(i) = 10._real64**(i*step_size + log10(xmin))
   xarr(i) = i*step_size + xmin
end do

end subroutine make_linear_array

!> Converts temperature array into grain enthalpy array (reverse of convert_E_arr_to_T_arr).
!> \todo maybe insert minimum temperature achievable (something like in IDL code for the case when small grains reach very low temperatures) 
subroutine convert_T_arr_to_E_arr(T_arr, E_arr, ic,ig) 
integer :: i, num , ic, ig 
real(kind=real64) :: E_arr(0:), T_arr(0:), vol, tc_min
integer :: i0, i1 
real(kind=real64) :: x0, x1, y0(0:0), y1(0:0),x, y(0:0)

! check array size  
num =size(T_arr)
if (size(E_arr) /= num) then 
   print *, 'STOP(convert_T_arr_to_E_arr):: input E_arr array does not have the right size'
   stop
endif 

vol = 4./3.*pi*dust_size_fa(ic,ig)**3

tc_min = cal_temp(0, iq_ct_table(ic))

do i = 0, num-1  ! this can be expensive and useless when you are close to equilibrium and all values are within the same energy bin. 

   if (T_arr(i) >= tc_min ) then  
   
      call value_locate(T_arr(i), cal_temp(:, iq_ct_table(ic)), i0,i1, .FALSE.)
      x0 = cal_temp(i0, iq_ct_table(ic))
      x1 = cal_temp(i1, iq_ct_table(ic))
      y0 = grain_enthalpy(i0, iq_ct_table(ic))*vol
      y1=  grain_enthalpy(i1, iq_ct_table(ic))*vol  
      x = T_arr(i)
      call lin_interpolate(y0,y1,x0,x1,x,y)
      E_arr(i) = y(0)   ! this is needed when calculating grain emission in the range of temperatures
    
    else

       x0 = cal_temp(0, iq_ct_table(ic))
       y0 = grain_enthalpy(0, iq_ct_table(ic))*vol
       x = T_arr(i)
       
       E_arr(i) = (x/x0)**3*y0(0)   !!! interpolation at low temperatures 

    endif
      
end do

end subroutine convert_T_arr_to_E_arr

!> Reads and interpolates stellar emission library defined in the input (see file_stellar_library() and stellar_library()). 
subroutine prepare_stellar_emission_library

  call read_stellar_library
  call interpolate_lum_to_mass_lib


end subroutine prepare_stellar_emission_library


 
!> Interpolates the values of the stellar luminosity-to-mass ratio in lum_to_mass_lib() to the wavelengths in the wavelength grid lambda_arr(). It stores the values in lum_to_mass_int() and deallocate lum_to_mass_lib().  
subroutine interpolate_lum_to_mass_lib
  integer :: il, ia, iz 
  integer :: i0, i1
  real(kind=real64) :: x, x0, x1, y(0:0), y0(0:0), y1(0:0)

  allocate(lum_to_mass_int(0:lnum_tot-1, 0:nage_lib-1, 0:nmet_lib-1))

  do il = 0, lnum_tot-1 

     call value_locate(lambda_arr(il), lambda_lib, i0, i1, .FALSE.)
     
     do ia = 0, nage_lib -1 

        do iz = 0, nmet_lib -1 

           x0 = lambda_lib(i0)
           x1 = lambda_lib(i1)
           y0 = lum_to_mass_lib(i0, ia, iz)
           y1 = lum_to_mass_lib(i1, ia,iz)
           x = lambda_arr(il)
           call lin_interpolate(y0,y1,x0,x1,x,y)
           lum_to_mass_int(il,ia,iz) = y(0)            

        end do
        
     end do

  end do

  deallocate(lum_to_mass_lib)

end subroutine interpolate_lum_to_mass_lib

!> Performs bilinear interpolation given the function farr, evaluated at the points (x1,y1), (x2,y1), (x1, y2), (x2,y2), and the output value coordinates (x,y). 
subroutine bilinear_interpolation(xarr, yarr, farr, x, y, fout) 
  real(kind=real64) :: xarr(2), yarr(2), farr(2,2), x,y,fout
  real(kind=real64) :: dx_dy(2,2), den

  dx_dy(1,1) = (xarr(2)-x)*(yarr(2)-y)
  dx_dy(2,1) = (x-xarr(1))*(yarr(2)-y)
  dx_dy(1,2) = (xarr(2)-x)*(y-yarr(1))
  dx_dy(2,2) = (x-xarr(1))*(y-yarr(1))
  den = (xarr(2)-xarr(1))*(yarr(2)-yarr(1))

  fout = sum(farr*dx_dy)/den

end subroutine bilinear_interpolation

!> Sets stellar particle luminosities star_lum() for an input wavelength lambda_in. 
subroutine set_star_particle_luminosity(lambda_in)
  real(kind=real64) :: lambda_in 
  integer :: ia0, ia1
  integer :: iz0, iz1
  integer :: i
  real(kind=real64) :: z_part, t_part, lm_ratio_part
  integer ::  il
  real(kind=real64) :: xarr(2), yarr(2), farr(2,2)

  if (main_prc) print *, 'deriving stellar particle luminosity for wavelength [um]:', lambda_in

! check z_sun value 

if (z_sun < 0 .or. z_sun > 0.03) then
   if (main_prc) then 
      print *, 'STOP(set_star_particle_luminosity): z_sun value not allowed!'
      print *, 'z_sun =', z_sun
      print *, 'allowed range =  [0,0.03]'
   endif
   call stop_prc
endif

! --------------------------------------
! check lambda_in within the lambda_arr

call find_lambda_index(lambda_in, il)

!-----------------------------------------------------------------
! interpolate stellar mass to luminosity ratio at the particle age and metallicity

do i=0,tot_star_particles-1  

   t_part= agestar(i)*1.E9    !!! particle age in yr 

   call value_locate(t_part, age_lib, ia0, ia1, .FALSE.)

   z_part = z_sun*10.**fehstar(i)   ! particle metallicity 

   call value_locate(z_part, met_lib, iz0, iz1, .FALSE.)

   xarr = (/met_lib(iz0),met_lib(iz1)/)   ! x1, x2
   yarr = (/age_lib(ia0), age_lib(ia1)/)  ! y1, y2
   

   farr(1,1) = log10(lum_to_mass_int(il, ia0, iz0))
   farr(1,2) = log10(lum_to_mass_int(il, ia1, iz0))
   farr(2,1) = log10(lum_to_mass_int(il, ia0, iz1))
   farr(2,2) = log10(lum_to_mass_int(il, ia1, iz1))

   call bilinear_interpolation(xarr, yarr, farr, z_part, t_part, lm_ratio_part) 

   lm_ratio_part = 10.**lm_ratio_part

   star_lum(i)=mstar(i)*lm_ratio_part 

   if (star_lum(i) < 0) then 
      print *, 'STOP(set_star_particle_luminosity): star lum negative! Check interpolation stellar emission library!'
      print *, star_lum(i)
      print *, farr
      stop
   endif

end do

call print_done


end subroutine set_star_particle_luminosity

!> Finds index in lambda_arr corresponding to the input wavelength lambda_in. 
subroutine find_lambda_index(lambda_in, il)
  real(kind=real64) :: lambda_in
  integer :: il,iq(0:0)

  iq = minloc(abs(lambda_arr-lambda_in))
  il = iq(0)-1

  if (abs(lambda_arr(il)-lambda_in)/lambda_in > 1E-3) then 
     print *, 'STOP(find_lambda_index): lambda_in not found within lambda_arr!'
     print *, 'lambda_in = ', lambda_in
     stop
  endif

end subroutine find_lambda_index

!> Calculates integrals of the stellar emission radiation field spectra in the UV (<4430 A) and optical (> 4430) and assign each cell to a 2-dim bin of radiation field total intensity. Used in the SED library approach for the stochastically heated dust emission.   
subroutine bin_rad_field
  real(kind=real64) :: lambda_sep = 0.443  ! wavelength between UV and optical range 
  integer :: i0, i1, nw_uv, nw_opt 
  integer :: i, iuv0, iuv1, iopt0, iopt1, iuv, iopt
  real(kind=real64), allocatable :: lambda_uv(:), lambda_opt(:)
  real(kind=real64) :: max_int_uv, min_int_uv
  real(kind=real64) :: max_int_opt, min_int_opt
  
! Allocate arrays for the integrated intensity radiation field in UV and optical (including only stellar emission) 
if (.not. allocated(int_rf_uv)) then 
   allocate(int_rf_uv(0:tot_ncell-1), int_rf_opt(0:tot_ncell-1))
   int_rf_uv = 0 
   int_rf_opt = 0 
else
   return  ! return if binning of the radiation field already performed. This is fine because at the moment you are including only the stellar emission which does not change during the dust heating iterations. 
endif

! allocate UV and optical wavelength arrays 
call value_locate(lambda_sep, lambda_arr,i0,i1, .FALSE.)
nw_uv = i1 ! This is at least 1 
nw_opt = lnum_stars-i0 
allocate(lambda_uv(0:nw_uv-1), lambda_opt(0:nw_opt-1))
lambda_uv = lambda_arr(0:nw_uv-1)*1E-6  ! um -> m 
lambda_opt = lambda_arr(i0: i0+nw_opt-1)*1E-6  

! Calculate integrated radiation fields for dusty cells 
do i = 0, tot_ncell -1 

   if (cchild(i) /= -1 .or. dens_ref(i) == 0) cycle
  
   ! make integration radiation field here ! !!! 
   int_rf_uv(i) = sum((u_final_uv_opt(0:nw_uv-2,i)+u_final_uv_opt(1:nw_uv-1,i))/2./(lambda_uv(1:nw_uv-1)-lambda_uv(0:nw_uv-2)))

   int_rf_opt(i) = sum((u_final_uv_opt(i0:i0+nw_opt-2,i)+u_final_uv_opt(1+i0:i0+nw_opt-1,i))/2./(lambda_opt(1:nw_opt-1)-lambda_opt(0:nw_opt-2)))
  
end do 

! Find max and minimum intensities 
max_int_uv = maxval(int_rf_uv)
min_int_uv = minval(int_rf_uv, mask=int_rf_uv > 0)

max_int_opt = maxval(int_rf_opt)
min_int_opt = minval(int_rf_opt, mask=int_rf_opt > 0)

! Make arrays of integrated intensities values defining the bins 
allocate(rf_uv_arr(0:n_int_rf_bins), rf_opt_arr(0:n_int_rf_bins))  ! without -1 so number of bins is n_int_rf_bins
call make_log_array(min_int_uv, max_int_uv,rf_uv_arr)
call make_log_array(min_int_opt, max_int_opt,rf_opt_arr)

! allocate average radiation field spectra array 
allocate(u_av_uv_opt(0:lnum_stars-1, 0:n_int_rf_bins-1, 0:n_int_rf_bins-1), u_av_dust(0:lnum-1, 0:n_int_rf_bins-1, 0:n_int_rf_bins-1),count_spectra_uv_opt(0:n_int_rf_bins-1, 0:n_int_rf_bins-1))
u_av_uv_opt = 0
u_av_dust = 0  
count_spectra_uv_opt = 0

! Calculate average RFED within each bin (for both stellar and dust emission) 

do i = 0, tot_ncell-1 

   if (cchild(i) /= -1 .or. dens_ref(i) == 0) cycle

   call value_locate(int_rf_uv(i), rf_uv_arr, iuv0, iuv1, .FALSE.)
   call value_locate(int_rf_opt(i), rf_opt_arr, iopt0, iopt1, .FALSE.)

   u_av_uv_opt(:,iuv0,iopt0) = u_av_uv_opt(:,iuv0,iopt0) + u_final_uv_opt(:,i)
   u_av_dust(:,iuv0,iopt0) = u_av_dust(:,iuv0,iopt0) + u_final_arr(:,i)  ! note that this is in different units compared to u_final_uv_opt. see convert_ufield_ifield and calc_dens_dustem subroutines 
   count_spectra_uv_opt(iuv0,iopt0) = count_spectra_uv_opt(iuv0,iopt0) + 1

end do 

do iuv = 0, n_int_rf_bins-1 
   do iopt = 0, n_int_rf_bins-1 

      u_av_uv_opt(:,iuv,iopt) = u_av_uv_opt(:,iuv,iopt)/count_spectra_uv_opt(iuv,iopt)
      u_av_dust(:,iuv,iopt) = u_av_dust(:,iuv,iopt)/count_spectra_uv_opt(iuv,iopt)

   enddo
end do 

deallocate(lambda_uv, lambda_opt)

end subroutine bin_rad_field



END MODULE sed_routines
