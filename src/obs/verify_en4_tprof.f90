PROGRAM verify_en4_tprof
  USE common
  USE params_model
  USE vars_model
  USE common_oceanmodel
  USE params_obs,                ONLY: nobs, id_t_obs, id_s_obs, id_u_obs, id_v_obs, id_eta_obs
  USE params_obs,                ONLY: DO_INSITU_to_POTTEMP, DO_POTTEMP_to_INSITU
  USE params_obs,                ONLY: DO_REMOVE_65N
  USE vars_obs
  USE common_obs_oceanmodel
  USE gsw_pot_to_insitu,         ONLY: t_from_pt, p_from_z, sa_from_sp, pt_from_t
  USE read_en4,                 ONLY: argo_data, read_en4_nc, ocn_profile
  ! For dynamic instantiation and use of namelists:
  USE input_nml_oceanmodel, ONLY: read_input_namelist

  IMPLICIT NONE

  CHARACTER(slen) :: obsinfile='obsin.nc'     !IN (default)
  CHARACTER(slen) :: guesfile='gues'          !IN (default) i.e. prefix to '.ocean_temp_salt.res.nc'
  CHARACTER(slen) :: obsoutfile='obsout.dat'  !OUT(default)

  TYPE(ocn_profile) :: prof
  REAL(r_size), ALLOCATABLE :: depth1d(:) ! nlev+1
  REAL(r_size), ALLOCATABLE :: sim1d(:)   ! nlev+1
  REAL(r_size), ALLOCATABLE :: oerr(:)    ! prof%nprof
  REAL(r_size), ALLOCATABLE :: olon(:)    ! prof%nprof
  INTEGER :: ish, ieh
  
  REAL(r_size), ALLOCATABLE :: v3d(:,:,:,:)
  REAL(r_size), ALLOCATABLE :: v2d(:,:,:)
  REAL(r_size), ALLOCATABLE :: depth2d(:,:,:) ! (nlon,nlat,nlev)

  REAL(r_size) :: dk,tg,qg
  REAL(r_size) :: ri,rj,rk
  INTEGER :: i,j,k,n

  ! For observation error modification:
  LOGICAL :: remap_obs_coords = .true.  ! Remap to mom4p1's -280 to 80 longitude range

  !-----------------------------------------------------------------------------
  ! Debugging parameters
  !-----------------------------------------------------------------------------
  INTEGER :: bdyobs=2                 !STEVE: use of boundary obs.
                                      !       1 := less restrictive, remove obs inside boundary
                                      !       2 := remove all observations touching a boundary
  LOGICAL :: debug_obsfilter = .false.
  LOGICAL :: debug_hdxf_0 = .false.   ! This error occured because there was not a model representation 
                                      ! of the observed value (i.e. SST obs with no SST model field).
                                      ! Solution was to populate a SST model field (v2d) with 
                                      ! surface temp data from the model (v3d(:,:,1)).
  !STEVE: to adjust writing to output file
  LOGICAL :: verbose = .false.
  LOGICAL :: dodebug1 = .false.

  INTEGER :: cnt_obs_u=0, cnt_obs_v=0, cnt_obs_t=0, cnt_obs_s=0
  INTEGER :: cnt_obs_ssh=0, cnt_obs_sst=0, cnt_obs_sss=0, cnt_obs_eta=0
  INTEGER :: cnt_obs_x=0, cnt_obs_y=0, cnt_obs_z=0
  INTEGER :: cnt_thin=0
  INTEGER, DIMENSION(nv3d+nv2d), SAVE :: cnt_obs = 0

  !STEVE: for debugging observation culling:
  INTEGER :: cnt_yout=0, cnt_xout=0, cnt_zout=0, cnt_triout=0
  INTEGER :: cnt_rigtnlon=0, cnt_nearland=0, cnt_oerlt0=0, cnt_altlatrng=0
  INTEGER :: cnt_nl1=0, cnt_nl2=0, cnt_nl3=0

  !-----------------------------------------------------------------------------
  ! Instantiations specific to this observation type:
  !-----------------------------------------------------------------------------
  ! For potential temperature conversion to in situ:
  REAL(r_size) :: p,pt,sp
  ! For in situ temperature conversion to potential temperature:
  REAL(r_size) :: sa,t,p_ref
  INTEGER :: nobs_s ! number of salinity obs (should match number of temperature obs)
  LOGICAL :: DO_REMOVE_BLACKSEA=.true.
  INTEGER :: cnt_blacksea=0
  
  !----------------------------------------------------------------------------
  ! Read in namelist parameters
  !----------------------------------------------------------------------------
  CALL read_input_namelist !STEVE: this has been moved to input_nml_{oceanmodel}.f90 since it needed to be slightly different for each model

  !-----------------------------------------------------------------------------
  ! Initialize the common_oceanmodel module, and process command line options
  !-----------------------------------------------------------------------------
  CALL set_common_oceanmodel
  CALL process_command_line !(get: -obsin <obsinfile> -gues <guesfile> -obsout <obsoutfile>)

  !-----------------------------------------------------------------------------
  ! Read observations from file
  !-----------------------------------------------------------------------------
  CALL prof%read_nc(trim(obsinfile))
  CALL prof%inspect()

  !-----------------------------------------------------------------------------
  ! Update the coordinate to match the model grid
  ! (Stolen from common_obs_mom4.f90::read_obs
  !-----------------------------------------------------------------------------
  allocate(oerr(prof%nprof),olon(prof%nprof))
  olon(:) = prof%lon(:)
  oerr = 1.d0
  if (remap_obs_coords) then
    CALL center_obs_coords(olon(:),oerr,prof%nprof)
  endif

  !-----------------------------------------------------------------------------
  ! Read model forecast for this member
  !-----------------------------------------------------------------------------
  ALLOCATE( v3d(nlon,nlat,nlev,nv3d) )
  ALLOCATE( v2d(nlon,nlat,nv2d) )
  WRITE(6,*) "Reading background/forecast data on model grid..."
  CALL read_diag(guesfile,v3d,v2d)
  WRITE(6,*) '****************'

!  !-----------------------------------------------------------------------------
!  ! Cycle through all observations
!  !-----------------------------------------------------------------------------
!  WRITE(6,*) "Cycling through observations..."
!  ohx=0.0d0
!  oqc=0
!  DO n=1,nobs
!    !---------------------------------------------------------------------------
!    ! Count bad obs errors associated with observations, and skip the ob
!    !---------------------------------------------------------------------------
!    if (oerr(n) <= 0) then
!      !STEVE: this occurred for a few synthetic obs, perhaps due to Dave's code generating obs errors
!      cnt_oerlt0 = cnt_oerlt0 + 1
!      CYCLE
!    endif
!
!    !---------------------------------------------------------------------------
!    ! Convert the physical coordinate to model grid coordinate (note: real, not integer)
!    !---------------------------------------------------------------------------
!    CALL phys2ijk(elem(n),rlon(n),rlat(n),rlev(n),ri,rj,rk) !(OCEAN)
!
!    !---------------------------------------------------------------------------
!    ! Filter in the tripolar region until localization is examined in the arctic !(ISSUE)
!    !---------------------------------------------------------------------------
!    if (DO_REMOVE_65N .and. rlat(n) > 65) then
!      if (verbose) WRITE(6,'(A)') "Latitude above 65.0N, in tripolar region. Removing observation..."
!      cnt_triout = cnt_triout + 1
!      CYCLE
!    endif
!    
!    !---------------------------------------------------------------------------
!    ! Filter out observations in the black sea, which may cause problems
!    !---------------------------------------------------------------------------
!    if (DO_REMOVE_BLACKSEA) then
!      if (40 < rlat(n) .and. rlat(n) < 50 .and. &
!          25 < rlon(n) .and. rlon(n) < 45) then
!        cnt_blacksea = cnt_blacksea + 1
!        CYCLE
!      endif
!    endif
!
!    !---------------------------------------------------------------------------
!    ! Filter out observations that are out of range for the grid
!    !---------------------------------------------------------------------------
!    if (CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) then
!      if (verbose) WRITE(6,'(A)') '* X-coordinate out of range'
!      if (verbose) WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', olon=', rlon(n)
!      cnt_xout = cnt_xout + 1
!      CYCLE
!    endif
!    if (CEILING(rj) < 2 .OR. nlat < CEILING(rj)) then
!      if (verbose) WRITE(6,'(A)') '* Y-coordinate out of range'
!      if (verbose) WRITE(6,'(A,F6.2,A,F6.2)') '*   rj=',rj,', olat=',rlat(n)
!      cnt_yout = cnt_yout + 1
!      CYCLE
!    endif
!    !STEVE: check against kmt, not nlev (OCEAN)
!    if (CEILING(rk) > nlev) then
!      CALL itpl_2d(kmt0,ri,rj,dk)
!      WRITE(6,'(A)') '* Z-coordinate out of range: rk > nlev'
!      cnt_zout = cnt_zout + 1
!      CYCLE
!    endif
!    if (CEILING(rk) < 2 .AND. rk < 1.00001d0) then   !(OCEAN)
!      rk = 1.00001d0                                 !(OCEAN)
!    endif                                            !(OCEAN)
!
!    !---------------------------------------------------------------------------
!    ! Check the observation against boundaries
!    !---------------------------------------------------------------------------
!    !STEVE: Check to make sure it's in the ocean, as determined       (OCEAN)
!    !       by mom4's topography map.
!    ! (note: must do it after coordinate checks, or the coordinate
!    !        could be outside of the range of the kmt array)
!    boundary_points : if (ri > nlon) then
!      !STEVE: I have to check what it does for this case...
!      !       but it causes an error in the next line if ri > nlon
!      if (verbose) WRITE(6,'(A)') '* coordinate is not on mom4 model grid: ri > nlon'
!      cnt_rigtnlon = cnt_rigtnlon + 1
!      if (cnt_rigtnlon <= 1) then
!        WRITE(6,*) "STEVE: ri > nlon (cnt_rigtnlon)"
!        WRITE(6,*) "ri = ", ri
!        WRITE(6,*) "nlon = ", nlon
!        WRITE(6,*) "rj = ", rj
!        WRITE(6,*) "elem(n) = ", elem(n)
!        WRITE(6,*) "rlon(n) = ", rlon(n)
!        WRITE(6,*) "rlat(n) = ", rlat(n)
!        WRITE(6,*) "rlev(n) = ", rlev(n)
!        WRITE(6,*) "rk = ", rk
!      endif
!      CYCLE
!    else
!      !STEVE: check this, case 1 allows more observations, case 2 is more restrictive
!      select case (bdyobs)
!      case(1)
!        if (kmt(NINT(ri),NINT(rj)) .lt. 1) then
!          if (debug_obsfilter) then
!            WRITE(6,'(A)') '* coordinate is on or too close to land, according to kmt'
!            WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj
!            WRITE(6,*) "kmt cell = ", kmt(NINT(ri),NINT(rj))
!          endif
!          cnt_nearland = cnt_nearland + 1
!          cnt_nl1 = cnt_nl1+1
!          CYCLE
!        elseif (kmt(NINT(ri),NINT(rj)) .lt. rk) then
!          if (debug_obsfilter) then
!            WRITE(6,'(A)') '* coordinate is on or too close to underwater topography, according to kmt'
!            WRITE(6,'(A,F6.2,A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj, ', rk=',rk
!            WRITE(6,*) "kmt cell = ", kmt(NINT(ri),NINT(rj))
!          endif
!          cnt_nearland = cnt_nearland + 1
!          cnt_nl2 = cnt_nl2+1
!          CYCLE
!        endif
!      case(2)
!        if(kmt(CEILING(ri),CEILING(rj)) .lt. 1 .or. &
!             kmt(CEILING(ri),FLOOR(rj)) .lt. 1 .or. &
!             kmt(FLOOR(ri),CEILING(rj)) .lt. 1 .or. &
!             kmt(FLOOR(ri),FLOOR(rj)) .lt. 1) THEN
!
!          if (debug_obsfilter) then
!            WRITE(6,'(A)') '* coordinate is too close to land, according to kmt'
!            WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj
!            WRITE(6,*) "kmt cell = ", kmt(FLOOR(ri):CEILING(ri),FLOOR(rj):CEILING(rj))
!          endif
!          cnt_nearland = cnt_nearland + 1
!          cnt_nl1 = cnt_nl1+1
!          CYCLE
!        elseif(kmt(CEILING(ri),CEILING(rj)) .lt. rk .or. &
!                  kmt(CEILING(ri),FLOOR(rj)) .lt. rk .or. &
!                  kmt(FLOOR(ri),CEILING(rj)) .lt. rk .or. &
!                  kmt(FLOOR(ri),FLOOR(rj)) .lt. rk) THEN
!
!          if (debug_obsfilter) then
!            WRITE(6,'(A)') '* coordinate is too close to underwater topography, according to kmt'
!            WRITE(6,'(A,F6.2,A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj, ', rk=',rk
!            WRITE(6,*) "kmt cell = ", kmt(FLOOR(ri):CEILING(ri),FLOOR(rj):CEILING(rj))
!          endif
!          cnt_nearland = cnt_nearland + 1
!          cnt_nl2 = cnt_nl2+1
!          CYCLE
!        endif
!      end select
!    endif boundary_points
!
!    !---------------------------------------------------------------------------
!    ! observation operator (computes H(x)) for specified member
!    !---------------------------------------------------------------------------
!    CALL Trans_XtoY(elem(n),ri,rj,rk,v3d,v2d,ohx(n))
!
!    endif
!
!    oqc(n) = 1
!  enddo !1:nobs

  !-----------------------------------------------------------------------------
  ! Print out the counts of observations removed for various reasons
  !-----------------------------------------------------------------------------
  WRITE(6,*) "In obsop_tprof.f90:: observations removed for:"
  WRITE(6,*) "cnt_oerlt0 = ", cnt_oerlt0
  WRITE(6,*) "cnt_xout = ", cnt_xout
  WRITE(6,*) "cnt_yout = ", cnt_yout
  WRITE(6,*) "cnt_zout = ", cnt_zout
  WRITE(6,*) "cnt_triout = ", cnt_triout
  WRITE(6,*) "cnt_rigtnlon = ", cnt_rigtnlon
  WRITE(6,*) "cnt_nearland = ", cnt_nearland
  WRITE(6,*) "cnt_nl1 = ", cnt_nl1
  WRITE(6,*) "cnt_nl2 = ", cnt_nl2
  WRITE(6,*) "cnt_altlatrng = ", cnt_altlatrng
  WRITE(6,*) "cnt_thin = ", cnt_thin
  if (DO_REMOVE_BLACKSEA) WRITE(6,*) "cnt_blacksea = ", cnt_blacksea

  !-----------------------------------------------------------------------------
  ! Write the observations and their associated innovations to file
  !-----------------------------------------------------------------------------
  !CALL write_obs2(obsoutfile,nobs,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr)

  !DEALLOCATE( elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr,v3d,v2d )

CONTAINS

SUBROUTINE process_command_line
!===============================================================================
! Process command line arguments 
!===============================================================================
IMPLICIT NONE
INTEGER, PARAMETER :: slen2=1024
CHARACTER(slen2) :: arg1,arg2
INTEGER :: i, ierr
INTEGER, DIMENSION(3) :: values

! STEVE: add input error handling!
! inputs are in the format "-x xxx"
do i=1,COMMAND_ARGUMENT_COUNT(),2
  CALL GET_COMMAND_ARGUMENT(i,arg1)
  PRINT *, "In obsop_tprof.f90::"
  PRINT *, "Argument ", i, " = ",TRIM(arg1)

  select case (arg1)
    case('-obsin')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      obsinfile = arg2
    case('-gues')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      guesfile = arg2
    case('-obsout')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      obsoutfile = arg2
    case('-rm65N')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) DO_REMOVE_65N
    case('-rmBlackSea')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) DO_REMOVE_BLACKSEA
    case('-debug')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) dodebug1
    case('-remap')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) remap_obs_coords
    case default
      PRINT *, "ERROR: option is not supported: ", arg1
      PRINT *, "(with value : ", trim(arg2), " )"
      stop 1
  end select
enddo

END SUBROUTINE process_command_line

END PROGRAM verify_en4_tprof
