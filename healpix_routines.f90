MODULE healpix_routines
  use iso_fortran_env
  
  integer(KIND=int32), private, save, dimension(0:1023) :: pix2x=0, pix2y=0
  integer(KIND=int32), private, save, dimension(0:127) :: x2pix1=0,y2pix1=0
  character(len=3000)           :: string
  integer, parameter, private :: LCH=48
  INTEGER(KIND=int32), private, PARAMETER :: ns_max4=8192
  INTEGER(KIND=int32), private, PARAMETER :: ns_max=ns_max4
  integer(int32), parameter :: oddbits=89478485   ! 2^0 + 2^2 + 2^4+..+2^26
  integer(int32), parameter :: evenbits=178956970 ! 2^1 + 2^3 + 2^4+..+2^27 
  REAL(KIND=real64),PARAMETER, private :: halfpi=asin(1.), pi=2*halfpi, twopi= 2*pi, twothird=2./3.
CONTAINS

!=======================================================================
!     pix2ang_nest
!
!     renders theta and phi coordinates of the nominal pixel center
!     for the pixel number ipix (NESTED scheme)
!     given the map resolution parameter nside
!=======================================================================
!#ifdef DOI8B
!  subroutine pix2ang_nest_8(nside, ipix, theta, phi)
!    integer(int32), parameter :: MKD = I8B
!#else
  subroutine pix2ang_nest  (nside, ipix, theta, phi)
    !integer(int32), parameter :: MKD = INT32
!#endif
    INTEGER(KIND=INT32), INTENT(IN)  :: nside
    INTEGER(KIND=INT32), INTENT(IN)  :: ipix
    REAL(KIND=REAL64),     INTENT(OUT) :: theta, phi

    INTEGER(KIND=INT32) :: npix, npface, ipf
    INTEGER(KIND=INT32) :: ip_low, ip_trunc, ip_med, ip_hi, &
         &     jrt, jr, nr, jpt, jp, kshift, nl4, scale, i, ismax
    INTEGER(KIND=INT32) :: ix, iy, face_num
    REAL(KIND=REAL64)     :: z, fn, fact1, fact2

    ! coordinate of the lowest corner of each face
    INTEGER(KIND=INT32), dimension(1:12) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
    INTEGER(KIND=INT32), dimension(1:12) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2
    !-----------------------------------------------------------------------
!#ifdef DOI8B
!    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
!       call pix2ang_nest(nside, int(ipix,kind=int32), theta, phi)
!       return
!    endif
!#else
!    if (nside > ns_max4) call fatal_error("nside out of range")    
     if (nside > ns_max4) STOP "nside out of range"
!#endif
    npix = nside2npix(nside)       ! total number of points    
!    if (ipix <0 .or. ipix>npix-1) call fatal_error ("ipix out of range")
    if (ipix <0 .or. ipix>npix-1) STOP "ipix out of range"
    !     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) call mk_pix2xy()

    npface = nside * int(nside, kind=INT32)
    nl4    = 4*nside

    !     finds the face, and the number in the face
    face_num = ipix/npface  ! face number in {0,11}
    ipf = MODULO(ipix,npface)  ! pixel number in the face {0,npface-1}

    fn = real(nside, kind=real64)
    fact1 = 1.0_real64/(3.0_real64*fn*fn)
    fact2 = 2.0_real64/(3.0_real64*fn)

    !     finds the x,y on the face (starting from the lowest corner)
    !     from the pixel number
    if (nside <= ns_max4) then
       ip_low = iand(ipf,1023)       ! content of the last 10 bits
       ip_trunc =    ipf/1024        ! truncation of the last 10 bits
       ip_med = iand(ip_trunc,1023)  ! content of the next 10 bits
       ip_hi  =      ip_trunc/1024   ! content of the high weight 10 bits

       ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
       iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
    else
       ix = 0
       iy = 0
       scale = 1
       ismax = 4
       do i=0, ismax
          ip_low = iand(ipf,1023)
          ix = ix + scale * pix2x(ip_low)
          iy = iy + scale * pix2y(ip_low)
          scale = scale * 32
          ipf   = ipf/1024
       enddo
       ix = ix + scale * pix2x(ipf)
       iy = iy + scale * pix2y(ipf)
    endif

    !     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy  ! 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy  ! 'horizontal' in {-nside+1,nside-1}

    !     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1   ! ring number in {1,4*nside-1}

    if (jr < nside) then     ! north pole region
       nr = jr
       z = 1.0_real64 - nr * fact1 * nr
       kshift = 0

    else if (jr <= 3*nside) then ! equatorial region
       nr = nside
       z  = (2*nside-jr)*fact2
       kshift = iand(jr - nside, 1)

    else if (jr > 3*nside) then ! south pole region
       nr = nl4 - jr
       z = - 1.0_real64 + nr * fact1 * nr
       kshift = 0
    endif

    theta = ACOS(z)

    !     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1_INT32 + kshift)/2  ! 'phi' number in the ring in {1,4*nr}
    if (jp > nl4) jp = jp - nl4
    if (jp < 1)   jp = jp + nl4
    
    phi = (jp - (kshift+1)*0.5_real64) * (halfpi / nr)
   
    return

!#ifdef DOI8B
!  end subroutine pix2ang_nest_8
!#else
  end subroutine pix2ang_nest
!#endif

  subroutine mk_pix2xy()
    !=======================================================================
    !     constructs the array giving x and y in the face from pixel number
    !     for the nested (quad-cube like) ordering of pixels
    !
    !     the bits corresponding to x and y are interleaved in the pixel number
    !     one breaks up the pixel number by even and odd bits
    !=======================================================================
    INTEGER(KIND=INT32) ::  kpix, jpix, ix, iy, ip, id
    
    !cc cf block data      data      pix2x(1023) /0/
    !-----------------------------------------------------------------------
    !      print *, 'initiate pix2xy'
    do kpix=0,1023          ! pixel number
       jpix = kpix
       IX = 0
       IY = 0
       IP = 1               ! bit position (in x and y)
!        do while (jpix/=0) ! go through all the bits
       do
          if (jpix == 0) exit ! go through all the bits
          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in ix
          jpix = jpix/2
          IX = ID*IP+IX

          ID = MODULO(jpix,2)  ! bit value (in kpix), goes in iy
          jpix = jpix/2
          IY = ID*IP+IY

          IP = 2*IP         ! next bit (in x and y)
       enddo
       pix2x(kpix) = IX     ! in 0,31
       pix2y(kpix) = IY     ! in 0,31       
    enddo

    return
  end subroutine mk_pix2xy

  function nside2npix(nside) result(npix_result)
    !=======================================================================
    ! given nside, returns npix such that npix = 12*nside^2
    !  nside should be a power of 2 smaller than ns_max
    !  if not, -1 is returned
    ! EH, Feb-2000
    ! 2009-03-04: returns i8b result, faster
    !=======================================================================
    INTEGER(KIND=int64)             :: npix_result
    INTEGER(KIND=INT32), INTENT(IN) :: nside

    INTEGER(KIND=int64) :: npix
    CHARACTER(LEN=*), PARAMETER :: code = "nside2npix"
    !=======================================================================

    npix = (12_int64*nside)*nside   
    if (nside < 1 .or. nside > ns_max .or. iand(nside-1,nside) /= 0) then
       print*,code//": Nside="//trim(string_i(nside))//" is not a power of 2."
       npix = -1
    endif
    npix_result = npix

    return
  end function nside2npix

  !========================================
  ! function string(arg, format)
  ! accepts: logical, i8, I4, SP and DP
  !=======================================

  !--------------------------------
  function string_i(arg, format) result(str)
    integer(int32) :: arg
    character(len=*),   optional :: format
    character(len=LCH)           :: str
    if (present(format)) then
       write(str,format) arg
    else
       write(str,*) arg
    endif
    return
  end function string_i

!=======================================================================
!     ang2pix_nest
!
!     renders the pixel number ipix (NESTED scheme) for a pixel which contains
!     a point on a sphere at coordinates theta and phi, given the map
!     resolution parameter nside
!
! 2009-03-09: calculations done directly at nside rather than ns_max
!             with ifort, for x and y integers>0, 
!             iand(x,y-1) is slightly faster than modulo(x,y)
!=======================================================================
!#ifdef DOI8B
!  subroutine ang2pix_nest_8(nside, theta, phi, ipix)
!    integer(int32), parameter :: INT32 = i8b
!#else
  subroutine ang2pix_nest  (nside, theta, phi, ipix)
    !integer(kind=int32), parameter :: INT32 = int32
!#endif
    INTEGER(KIND=int32), INTENT(IN)  :: nside
    REAL(KIND=real64),     INTENT(IN)  :: theta, phi
    INTEGER(KIND=int32), INTENT(OUT) :: ipix

    integer(kind=int32) :: ipf, scale, scale_factor
    REAL(KIND=real64)     ::  z, za, tt, tp, tmp
    INTEGER(KIND=int32) :: jp, jm, ifp, ifm, face_num, &
         &     ix, iy, ix_low, iy_low, ntt, i, ismax, ipix4
    !-----------------------------------------------------------------------
!#ifdef DOI8B
!    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
!       call ang2pix_nest(nside, theta, phi, ipix4)
!       ipix = int(ipix4, kind=INT32)
!       return
!    endif
!    if (nside <1 .or. nside > ns_max ) call fatal_error("nside out of range")
!#else
!    if (nside <1 .or. nside > ns_max4) call fatal_error("nside out of range")
    if (nside <1 .or. nside > ns_max4) STOP "nside out of range"
!#endif
    
    if (theta<0.0_real64 .or. theta>pi)  then
       print*,"ANG2PIX_NEST: theta : ",theta," is out of range [0,Pi]"
!       call fatal_error
       STOP
    endif
    if (x2pix1(127) <= 0) call mk_xy2pix1()
    
    z  = COS(theta)
    za = ABS(z)
    tt = MODULO(phi, twopi) / halfpi  ! in [0,4[

    if (za <= twothird) then ! equatorial region

       !        (the index of edge lines increase when the longitude=phi goes up)
       jp = INT(nside*(0.5_real64 + tt - z*0.75_real64)) !  ascending edge line index
       jm = INT(nside*(0.5_real64 + tt + z*0.75_real64)) ! descending edge line index

       !        finds the face
       ifp = jp / nside  ! in {0,4}
       ifm = jm / nside
       if (ifp == ifm) then          ! faces 4 to 7
          face_num = iand(ifp,3) + 4
       else if (ifp < ifm) then     ! (half-)faces 0 to 3
          face_num = iand(ifp,3)
       else                            ! (half-)faces 8 to 11
          face_num = iand(ifm,3) + 8
       endif

       ix =         iand(jm, nside-1)
       iy = nside - iand(jp, nside-1) - 1

    else ! polar region, za > 2/3

       ntt = INT(tt)
       if (ntt >= 4) ntt = 3
       tp = tt - ntt
       tmp = SQRT( 3.0_real64*(1.0_real64 - za) )  ! in ]0,1]

       !        (the index of edge lines increase when distance from the closest pole goes up)
       jp = INT( nside * tp          * tmp ) ! line going toward the pole as phi increases
       jm = INT( nside * (1.0_real64 - tp) * tmp ) ! that one goes away of the closest pole
       jp = MIN(nside-1, jp) ! for points too close to the boundary
       jm = MIN(nside-1, jm)

       !        finds the face and pixel's (x,y)
       if (z >= 0) then
          face_num = ntt  ! in {0,3}
          ix = nside - jm - 1
          iy = nside - jp - 1
       else
          face_num = ntt + 8 ! in {8,11}
          ix =  jp
          iy =  jm
       endif

    endif
   
    if (nside <= ns_max4) then 
       ix_low = iand(ix, 127)
       iy_low = iand(iy, 127)
       ipf =     x2pix1(ix_low) + y2pix1(iy_low) &
            & + (x2pix1(ix/128) + y2pix1(iy/128)) * 16384      
    else
       scale = 1_INT32
       scale_factor = 16384_INT32 ! 128*128
       ipf = 0_INT32
       ismax = 1 ! for nside in [2^14, 2^20]
       if (nside >  1048576 ) ismax = 3
       do i=0, ismax
          ix_low = iand(ix, 127) ! last 7 bits
          iy_low = iand(iy, 127) ! last 7 bits
          ipf = ipf + (x2pix1(ix_low)+y2pix1(iy_low)) * scale
          scale = scale * scale_factor
          ix  =     ix / 128 ! truncate out last 7 bits
          iy  =     iy / 128
       enddo
       ipf =  ipf + (x2pix1(ix)+y2pix1(iy)) * scale       
    endif    
    
    ipix = ipf + face_num* int(nside,INT32) * nside    ! in {0, 12*nside**2 - 1}
    
    return

!#ifdef DOI8B
!  end subroutine ang2pix_nest_8
!#else
  end subroutine ang2pix_nest
!#endif

  subroutine mk_xy2pix1()
    !=======================================================================
    !     sets the array giving the number of the pixel lying in (x,y)
    !     x and y are in {1,128}
    !     the pixel number is in {0,128**2-1}
    !
    !     if  i-1 = sum_p=0  b_p * 2^p
    !     then ix = sum_p=0  b_p * 4^p
    !          iy = 2*ix
    !     ix + iy in {0, 128**2 -1}
    !=======================================================================
    INTEGER(KIND=INT32):: k,ip,i,j,id
    !=======================================================================

    do i = 0,127           !for converting x,y into
       j  = i           !pixel numbers
       k  = 0
       ip = 1

       do
          if (j==0) then
             x2pix1(i) = k
             y2pix1(i) = 2*k
             exit
          else
             id = MODULO(J,2)
             j  = j/2
             k  = ip*id+k
             ip = ip*4
          endif
       enddo

    enddo

    return
  end subroutine mk_xy2pix1


  !====================================================================
! The following is a routine which finds the 7 or 8 neighbours of
! any pixel in the nested scheme of the HEALPIX pixelisation.
!====================================================================
!  neighbours_nest
!
!   Returns list n(8) of neighbours of pixel ipix (in NESTED scheme)
!   the neighbours are ordered in the following way:
!   First pixel is the one to the south (the one west of the south
! direction is taken
! for the pixels which don't have a southern neighbour). From
! then on the neighbours are ordered in the clockwise direction
! about the pixel with number ipix.
!
!   nneigh is the number of neighbours (mostly 8, 8 pixels have 7 neighbours)
!
!   Benjamin D. Wandelt October 1997
!   Added to pix_tools in March 1999
!   added 'return' for case nside=1, EH, Oct 2005
!   corrected bugs in case nside=1 and ipix=7, 9 or 11, EH, June 2006
!   2009-06-16: deals with Nside > 8192
!====================================================================
!#ifdef DOI8B
!  subroutine neighbours_nest_8(nside, ipix, n, nneigh)
!    use bit_manipulation
!    integer(kind=int32), parameter  ::   INT32 = I8B
!#else
  subroutine neighbours_nest(nside, ipix, n, nneigh)
  !  use bit_manipulation
   ! integer(kind=int32), parameter  ::   INT32 = INT32
!#endif
    !====================================================================
    integer(kind=int32), intent(in)::  nside
    integer(kind=INT32), intent(in)::  ipix
    integer(kind=INT32), intent(out), dimension(1:):: n
    integer(kind=int32), intent(out):: nneigh

    integer(kind=int32) :: ix,ixm,ixp,iy,iym,iyp,ixo,iyo
    integer(kind=int32) :: face_num,other_face, ipix4
    integer(kind=int32) :: ia,ib,ibp,ibm,ib2,icase
    integer(kind=INT32) :: npix,ipf,ipo
    integer(kind=INT32) :: local_magic1,local_magic2,nsidesq
    integer(kind=int32), dimension(1:8) :: n4
    character(len=*), parameter :: code = "neighbours_nest"

!     integer(kind=int32), intrinsic :: IAND

    !--------------------------------------------------------------------
!#ifdef DOI8B
!    if (nside <= ns_max4) then ! use faster 32-bit routine whenever possible
!       ipix4 = int(ipix, kind=int32)
!       call neighbours_nest(nside, ipix4, n4, nneigh)
!       n(1:nneigh) = int(n4(1:nneigh), kind=INT32)
!       return
!    endif
!    if (nside <1 .or. nside > ns_max ) call fatal_error(code//"> nside out of range")
!#else
    if (nside <1 .or. nside > ns_max4) then

    ! call fatal_error(code//"> nside out of range")
       print *, 'neighbours_nest:  nside out of range'
       stop
       endif
!#endif
    npix = nside2npix(nside) ! total number of points
    nsidesq = npix / 12
    if (ipix <0 .or. ipix>npix-1) then
       
       !call fatal_error(code//"> ipix out of range")
       print *, 'neighbours_nest: ipix out of range'
       stop
    endif

    ! quick and dirty hack for Nside=1

    if (nside == 1) then
       nneigh = 6
       if (ipix==0 ) n(1:6) = (/ 8, 4, 3, 2, 1, 5 /)
       if (ipix==1 ) n(1:6) = (/ 9, 5, 0, 3, 2, 6 /)
       if (ipix==2 ) n(1:6) = (/10, 6, 1, 0, 3, 7 /)
       if (ipix==3 ) n(1:6) = (/11, 7, 2, 1, 0, 4 /)
       if (ipix==4 ) n(1:6) = (/11, 7, 3, 0, 5, 8 /)
       if (ipix==5 ) n(1:6) = (/ 8, 4, 0, 1, 6, 9 /)
       if (ipix==6 ) n(1:6) = (/ 9, 5, 1, 2, 7,10 /)
       if (ipix==7 ) n(1:6) = (/10, 6, 2, 3, 4,11 /)
       if (ipix==8 ) n(1:6) = (/10,11, 4, 0, 5, 9 /)
       if (ipix==9 ) n(1:6) = (/11, 8, 5, 1, 6,10 /)
       if (ipix==10) n(1:6) = (/ 8, 9, 6, 2, 7,11 /)
       if (ipix==11) n(1:6) = (/ 9,10, 7, 3, 4, 8 /)
       return
    endif

    !     initiates array for (x,y)-> pixel number -> (x,y) mapping
    if (x2pix1(127) <= 0) call mk_xy2pix1()

    local_magic1=(nsidesq-1)/3
    local_magic2=2*local_magic1
    face_num=ipix/nsidesq

    ipf=modulo(ipix,nsidesq)   !Pixel number in face

    call pix2xy_nest(nside,ipf,ix,iy)
    ixm=ix-1
    ixp=ix+1
    iym=iy-1
    iyp=iy+1

    nneigh=8                  !Except in special cases below

    !     Exclude corners
    if(ipf==local_magic2)     then !WestCorner
       icase=5
       goto 100
    endif
    if(ipf==(nsidesq-1)) then !NorthCorner
       icase=6
       goto 100
    endif
    if(ipf==0)           then !SouthCorner
       icase=7
       goto 100
    endif
    if(ipf==local_magic1)     then !EastCorner
       icase=8
       goto 100
    endif

    !     Detect edges
    if(IAND(ipf,local_magic1)==local_magic1) then !NorthEast
       icase=1
       goto 100
    endif
    if(IAND(ipf,local_magic1)==0)      then !SouthWest
       icase=2
       goto 100
    endif
    if(IAND(ipf,local_magic2)==local_magic2) then !NorthWest
       icase=3
       goto 100
    endif
    if(IAND(ipf,local_magic2)==0)      then !SouthEast
       icase=4
       goto 100
    endif

    !     Inside a face
    call xy2pix_nest(nside, ixm, iym, face_num, n(1))
    call xy2pix_nest(nside, ixm, iy , face_num, n(2))
    call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
    call xy2pix_nest(nside, ix , iyp, face_num, n(4))
    call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
    call xy2pix_nest(nside, ixp, iy , face_num, n(6))
    call xy2pix_nest(nside, ixp, iym, face_num, n(7))
    call xy2pix_nest(nside, ix , iym, face_num, n(8))
    return

100 continue

    ia= face_num/4            !in {0,2}
    ib= modulo(face_num,4)       !in {0,3}
    ibp=modulo(ib+1,4)
    ibm=modulo(ib+4-1,4)
    ib2=modulo(ib+2,4)

    if(ia==0) then          !North Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ibp
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          ipo=modulo(swapLSBMSB(ipf),nsidesq)    !East-West flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1 , iyo, other_face, n(5))
          n(6)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(7))
       case(2)              !SouthWest edge
          other_face=4+ib
          ipo=modulo(invLSB(ipf),nsidesq)        !SW-NE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
          n(2)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(3))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       case(3)              !NorthWest edge
          other_face=0+ibm
          ipo=modulo(swapLSBMSB(ipf),nsidesq)    !East-West flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(3))
          n(4)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(5))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       case(4)              !SouthEast edge
          other_face=4+ibp
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(7))
          n(8)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
       case(5)              !West corner
          nneigh=7
          other_face=4+ib
          n(2)=other_face*nsidesq+nsidesq-1
          n(1)=n(2)-2
          other_face=0+ibm
          n(3)=other_face*nsidesq+local_magic1
          n(4)=n(3)+2
          n(5)=ipix+1
          n(6)=ipix-1
          n(7)=ipix-2
       case(6)              !North corner
          n(1)=ipix-3
          n(2)=ipix-1
          n(8)=ipix-2
          other_face=0+ibm
          n(4)=other_face*nsidesq+nsidesq-1
          n(3)=n(4)-2
          other_face=0+ib2
          n(5)=other_face*nsidesq+nsidesq-1
          other_face=0+ibp
          n(6)=other_face*nsidesq+nsidesq-1
          n(7)=n(6)-1
       case(7)              !South corner
          other_face=8+ib
          n(1)=other_face*nsidesq+nsidesq-1
          other_face=4+ib
          n(2)=other_face*nsidesq+local_magic1
          n(3)=n(2)+2
          n(4)=ipix+2
          n(5)=ipix+3
          n(6)=ipix+1
          other_face=4+ibp
          n(8)=other_face*nsidesq+local_magic2
          n(7)=n(8)+1
       case(8)              !East corner
          nneigh=7
          n(2)=ipix-1
          n(3)=ipix+1
          n(4)=ipix+2
          other_face=0+ibp
          n(6)=other_face*nsidesq+local_magic2
          n(5)=n(6)+1
          other_face=4+ibp
          n(7)=other_face*nsidesq+nsidesq-1
          n(1)=n(7)-1
       end select ! north

    elseif(ia==1) then      !Equatorial region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=0+ib
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo , iyo+1, other_face, n(5))
          n(6)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(7))
       case(2)              !SouthWest edge
          other_face=8+ibm
          ipo=modulo(invLSB(ipf),nsidesq)        !SW-NE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
          n(2)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(3))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       case(3)              !NorthWest edge
          other_face=0+ibm
          ipo=modulo(invMSB(ipf),nsidesq)    !NW-SE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(3))
          n(4)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(5))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       case(4)              !SouthEast edge
          other_face=8+ib
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          ipo=modulo(invMSB(ipf),nsidesq) !SE-NW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(7))
          n(8)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
       case(5)              !West corner
          other_face=8+ibm
          n(2)=other_face*nsidesq+nsidesq-1
          n(1)=n(2)-2
          other_face=4+ibm
          n(3)=other_face*nsidesq+local_magic1
          other_face=0+ibm
          n(4)=other_face*nsidesq
          n(5)=n(4)+1
          n(6)=ipix+1
          n(7)=ipix-1
          n(8)=ipix-2
       case(6)              !North corner
          nneigh=7
          n(1)=ipix-3
          n(2)=ipix-1
          other_face=0+ibm
          n(4)=other_face*nsidesq+local_magic1
          n(3)=n(4)-1
          other_face=0+ib
          n(5)=other_face*nsidesq+local_magic2
          n(6)=n(5)-2
          n(7)=ipix-2
       case(7)              !South corner
          nneigh=7
          other_face=8+ibm
          n(1)=other_face*nsidesq+local_magic1
          n(2)=n(1)+2
          n(3)=ipix+2
          n(4)=ipix+3
          n(5)=ipix+1
          other_face=8+ib
          n(7)=other_face*nsidesq+local_magic2
          n(6)=n(7)+1
       case(8)              !East corner
          other_face=8+ib
          n(8)=other_face*nsidesq+nsidesq-1
          n(1)=n(8)-1
          n(2)=ipix-1
          n(3)=ipix+1
          n(4)=ipix+2
          other_face=0+ib
          n(6)=other_face*nsidesq
          n(5)=n(6)+2
          other_face=4+ibp
          n(7)=other_face*nsidesq+local_magic2
       end select ! equator
    else                    !South Pole region
       select case(icase)
       case(1)              !NorthEast edge
          other_face=4+ibp
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          ipo=modulo(invLSB(ipf),nsidesq)    !NE-SW flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo , iyo+1, other_face, n(5))
          n(6)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(7))
       case(2)              !SouthWest edge
          other_face=8+ibm
          ipo=modulo(swapLSBMSB(ipf),nsidesq)        !W-E flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(1))
          n(2)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(3))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
       case(3)              !NorthWest edge
          other_face=4+ib
          ipo=modulo(invMSB(ipf),nsidesq)    !NW-SE flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo-1, iyo, other_face, n(3))
          n(4)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo+1, iyo, other_face, n(5))
          call xy2pix_nest(nside, ixm, iym, face_num, n(1))
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ix , iym, face_num, n(8))
          call xy2pix_nest(nside, ixp, iym, face_num, n(7))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
       case(4)              !SouthEast edge
          other_face=8+ibp
          call xy2pix_nest(nside, ixm, iy , face_num, n(2))
          call xy2pix_nest(nside, ixm, iyp, face_num, n(3))
          call xy2pix_nest(nside, ix , iyp, face_num, n(4))
          call xy2pix_nest(nside, ixp, iyp, face_num, n(5))
          call xy2pix_nest(nside, ixp, iy , face_num, n(6))
          ipo=modulo(swapLSBMSB(ipf),nsidesq) !E-W flip
          call pix2xy_nest(nside,ipo,ixo,iyo)
          call xy2pix_nest(nside, ixo, iyo+1, other_face, n(7))
          n(8)=other_face*nsidesq+ipo
          call xy2pix_nest(nside, ixo, iyo-1, other_face, n(1))
       case(5)              !West corner
          nneigh=7
          other_face=8+ibm
          n(2)=other_face*nsidesq+local_magic1
          n(1)=n(2)-1
          other_face=4+ib
          n(3)=other_face*nsidesq
          n(4)=n(3)+1
          n(5)=ipix+1
          n(6)=ipix-1
          n(7)=ipix-2
       case(6)              !North corner
          n(1)=ipix-3
          n(2)=ipix-1
          other_face=4+ib
          n(4)=other_face*nsidesq+local_magic1
          n(3)=n(4)-1
          other_face=0+ib
          n(5)=other_face*nsidesq
          other_face=4+ibp
          n(6)=other_face*nsidesq+local_magic2
          n(7)=n(6)-2
          n(8)=ipix-2
       case(7)              !South corner
          other_face=8+ib2
          n(1)=other_face*nsidesq
          other_face=8+ibm
          n(2)=other_face*nsidesq
          n(3)=n(2)+1
          n(4)=ipix+2
          n(5)=ipix+3
          n(6)=ipix+1
          other_face=8+ibp
          n(8)=other_face*nsidesq
          n(7)=n(8)+2
       case(8)              !East corner
          nneigh=7
          other_face=8+ibp
          n(7)=other_face*nsidesq+local_magic2
          n(1)=n(7)-2
          n(2)=ipix-1
          n(3)=ipix+1
          n(4)=ipix+2
          other_face=4+ibp
          n(6)=other_face*nsidesq
          n(5)=n(6)+2
       end select ! south
    endif

    return

!#ifdef DOI8B
!  end subroutine neighbours_nest_8
!#else
  end subroutine neighbours_nest
!#endif
  
!=======================================================================
!     gives the pixel number ipix (NESTED)
!     corresponding to ix, iy and face_num
!
!     Benjamin D. Wandelt 13/10/97
!     using code from HEALPIX toolkit by K.Gorski and E. Hivon
!     2009-06-15: deals with Nside > 8192
!     2012-03-02: test validity of ix_in and iy_in instead of undefined ix and iy
!=======================================================================
!#ifdef DOI8B
!  subroutine xy2pix_nest_8(nside, ix_in, iy_in, face_num, ipix)
!    integer(kind=int32), parameter  ::   INT32 = I8B
!#else
  subroutine xy2pix_nest(nside, ix_in, iy_in, face_num, ipix)
    !integer(kind=int32), parameter  ::   INT32 = INT32
!#endif
    !=======================================================================
    INTEGER(KIND=INT32), INTENT(IN) ::  nside, ix_in, iy_in, face_num
    INTEGER(KIND=INT32), INTENT(OUT) :: ipix
    INTEGER(KIND=INT32) ::  ix, iy, ix_low, ix_hi, iy_low, iy_hi, i, ismax
    integer(kind=INT32) :: ipf, scale, scale_factor
    character(len=*), parameter :: code = "xy2pix_nest"

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) then 
       !  call fatal_error(code//"> nside out of range")
       print *, 'xy2pix_nest : nside out of range'
       stop
    endif
    
    if (ix_in<0 .or. ix_in>(nside-1)) then
       print *, 'xy2pix_nest : ix out of range' 
       ! call fatal_error(code//"> ix out of range")
       stop
    endif
    
    if (iy_in<0 .or. iy_in>(nside-1)) then 
       !  call fatal_error(code//"> iy out of range")
       print *, 'xy2pix_nest : iy out of range'
       stop
       endif 
    if (x2pix1(127) <= 0) call mk_xy2pix1()

    ix = ix_in
    iy = iy_in
    if (nside <= ns_max4) then 
       ix_low = iand(ix, 127)
       iy_low = iand(iy, 127)
       ipf =     x2pix1(ix_low) + y2pix1(iy_low) &
            & + (x2pix1(ix/128) + y2pix1(iy/128)) * 16384
    else
       scale = 1_INT32
       scale_factor = 16384_INT32 ! 128*128
       ipf = 0_INT32
       ismax = 1 ! for nside in [2^14, 2^20]
       if (nside >  1048576 ) ismax = 3
       do i=0, ismax
          ix_low = iand(ix, 127) ! last 7 bits
          iy_low = iand(iy, 127) ! last 7 bits
          ipf = ipf + (x2pix1(ix_low)+y2pix1(iy_low)) * scale
          scale = scale * scale_factor
          ix  =     ix / 128 ! truncate out last 7 bits
          iy  =     iy / 128
       enddo
       ipf =  ipf + (x2pix1(ix)+y2pix1(iy)) * scale
    endif
    ipix = ipf + face_num* int(nside,INT32) * nside    ! in {0, 12*nside**2 - 1}
    return

!#ifdef DOI8B
!  end subroutine xy2pix_nest_8
!#else
  end subroutine xy2pix_nest
!#endif


!=======================================================================
!  pix2xy_nest
!     gives the x, y coords in a face from pixel number within the face (NESTED)
!
!     Benjamin D. Wandelt 13/10/97
!
!     using code from HEALPIX toolkit by K.Gorski and E. Hivon
!     2009-06-15: deals with Nside > 8192
!     2012-03-02: test validity of ipf_in instead of undefined ipf
!                 define ipf as INT32
!     2012-08-27:  corrected bug on (ix,iy) for Nside > 8192 (MARK)
!=======================================================================
!#ifdef DOI8B
!  subroutine pix2xy_nest_8(nside, ipf_in, ix, iy)
!    integer(kind=int32), parameter  ::   INT32 = I8B
!#else
  subroutine pix2xy_nest  (nside, ipf_in, ix, iy)
   ! integer(kind=int32), parameter  ::   INT32 = INT32
!#endif
    INTEGER(KIND=INT32), INTENT(IN)  :: nside
    INTEGER(KIND=INT32), INTENT(IN)  :: ipf_in
    INTEGER(KIND=INT32), INTENT(OUT) :: ix, iy

    integer(kind=INT32) :: ipf
    INTEGER(KIND=INT32) ::  ip_low, ip_trunc, ip_med, ip_hi, scale, i, ismax
    character(len=*), parameter :: code = "pix2xy_nest"

    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) then 
    !call fatal_error(code//"> nside out of range")
       print *, 'pix2xy_nest: nside out of range'
       stop
       endif 
       
       if (ipf_in<0 .or. ipf_in>nside*nside-1) then 
          !  call fatal_error(code//"> ipix out of range")
          print *, 'pix2xy_nest: ipix out of range'
          stop
          endif
          
    if (pix2x(1023) <= 0) call mk_pix2xy()

    ipf = ipf_in
    if (nside <= ns_max4) then
       ip_low = iand(ipf,1023_INT32)   ! content of the last 10 bits
       ip_trunc =    ipf/1024        ! truncation of the last 10 bits
       ip_med = iand(ip_trunc,1023)  ! content of the next 10 bits
       ip_hi  =      ip_trunc/1024   ! content of the high weight 10 bits

       ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
       iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
    else
       ix = 0
       iy = 0
       scale = 1
       ismax = 4
       do i=0, ismax
          ip_low = iand(ipf,1023_INT32)
          ix = ix + scale * pix2x(ip_low)
          iy = iy + scale * pix2y(ip_low) ! corrected 2012-08-27
          scale = scale * 32
          ipf   = ipf/1024
       enddo
       ix = ix + scale * pix2x(ipf)
       iy = iy + scale * pix2y(ipf) ! corrected 2012-08-27
    endif

    return

!#ifdef DOI8B
!  end subroutine pix2xy_nest_8
!#else
  end subroutine pix2xy_nest
!#endif  


!! Returns i with even and odd bit positions interchanged.
function swapLSBMSB(i)
  integer(int32) :: swapLSBMSB
  integer(int32), intent(in) :: i

  swapLSBMSB = IAND(i,evenbits)/2 + IAND(i,oddbits)*2
end function swapLSBMSB  


!! Returns i with even (0,2,4,...) bits inverted.
function invMSB(i)
  integer(int32) :: invMSB
  integer(int32), intent(in) :: i

  invMSB = IEOR(i,evenbits)
end function invMSB

!! Returns i with odd (1,3,5,...) bits inverted.
function invLSB(i)
  integer(int32) :: invLSB
  integer(int32), intent(in) :: i

  invLSB = IEOR(i,oddbits)
end function invLSB

  
  
END MODULE healpix_routines
