PROGRAM verify_en4_tprof
  USE common,       ONLY : slen, r_size
  USE params_model, ONLY : nlon, nlat, nlev, nv3d, nv2d, iv3d_t, iv3d_s, iv3d_h
  USE vars_model,         ONLY: lon2d, lat2d, kmt, wet
  USE common_oceanmodel,  ONLY: set_common_oceanmodel, read_diag
  USE params_obs,                ONLY: nobs, id_t_obs, id_s_obs
  USE params_obs,                ONLY: DO_REMOVE_65N
  USE common_obs_oceanmodel,     ONLY: phys2ij, phys2ij_nearest, center_obs_coords
  USE read_en4,                 ONLY: ocn_profile, OCN_PROFILE_VALID, OCN_PROFILE_NOT_VALID
  USE m_interp,                 ONLY: spline_1d, interp_1d, t_interp
  ! For dynamic instantiation and use of namelists:
  USE input_nml_oceanmodel, ONLY: read_input_namelist

  IMPLICIT NONE

  CHARACTER(slen) :: obsinfile='obsin.nc'     !IN (default)
  CHARACTER(slen) :: guesfile='gues.nc'          !IN (default) i.e. prefix to '.ocean_temp_salt.res.nc'
  CHARACTER(slen) :: ovar='potm'
  CHARACTER(slen) :: obsoutfile='obsout_obsbkgd.nc'  !OUT(default)
  INTEGER         :: intp_opt = 1 ! 1-spline; 2-linear

  TYPE(ocn_profile) :: prof
  REAL(r_size), ALLOCATABLE :: depth1d(:) ! nlev+1
  REAL(r_size), ALLOCATABLE :: sim1d(:)   ! nlev+1
  REAL(r_size), ALLOCATABLE :: oerr(:)    ! prof%nprof
  REAL(r_size), ALLOCATABLE :: olon(:)    ! prof%nprof
  TYPE(t_interp) :: ip
  INTEGER :: ivar
  INTEGER :: ios, ioe, ims, ime
  INTEGER :: cnt_coldt
  
  REAL(r_size), ALLOCATABLE :: v3d(:,:,:,:)
  REAL(r_size), ALLOCATABLE :: v2d(:,:,:)
  REAL(r_size), ALLOCATABLE :: depth3d(:,:,:) ! (nlon,nlat,nlev)
  REAL(r_size), ALLOCATABLE :: maxdepth2d(:,:) ! (nlon,nlat)
  INTEGER, ALLOCATABLE :: maxnlev(:,:) ! (nlon,nlat)
  REAL(r_size), PARAMETER :: MAX_VALID_H = 10000.d0

  REAL(r_size) :: ri,rj,rk
  REAL(r_size) :: rii, rjj ! idx for nearest search
  INTEGER :: i,j,k,n
  INTEGER :: ii, jj  !idx for nearest search

  ! For observation error modification:
  LOGICAL :: remap_obs_coords = .true.  ! Remap to mom4p1's -280 to 80 longitude range

  !-----------------------------------------------------------------------------
  ! Debugging parameters
  !-----------------------------------------------------------------------------
  LOGICAL :: debug_obsfilter = .false.
  !LOGICAL :: debug_prof = .false.
  LOGICAL :: debug_prof = .true.
  !STEVE: to adjust writing to output file
  LOGICAL :: verbose = .false.

  !STEVE: for debugging observation culling:
  INTEGER :: cnt_yout=0, cnt_xout=0, cnt_triout=0
  INTEGER :: cnt_rigtnlon=0, cnt_nearland=0
  INTEGER :: cnt_nl1=0, cnt_nl2=0, cnt_nl3=0

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

  if (trim(ovar)=='potm') then
     ivar = iv3d_t
  elseif (trim(ovar)=='psal') then
     ivar = iv3d_s
  else
     write(6,*) "[err] verify_en4: unrecognized ovar (",trim(ovar),")."
     stop (10)
  endif

  !-----------------------------------------------------------------------------
  ! Read observations from file
  !-----------------------------------------------------------------------------
  CALL prof%read_nc(trim(obsinfile))
  CALL prof%inspect()

  !-----------------------------------------------------------------------------
  ! Update the coordinate to match the model grid
  ! (Stolen from common_obs_mom4.f90::read_obs
  !-----------------------------------------------------------------------------
  ALLOCATE(olon(prof%nprof))
  olon(:) = prof%lon(:)
  if (remap_obs_coords) then
    ALLOCATE(oerr(prof%nprof))
    oerr = 1.d0 ! placeholder, not used
    CALL center_obs_coords(olon(:),oerr,prof%nprof)
    DEALLOCATE(oerr)
  endif

  !-----------------------------------------------------------------------------
  ! Read model forecast for this member
  !-----------------------------------------------------------------------------
  ALLOCATE( v3d(nlon,nlat,nlev,nv3d) )
  ALLOCATE( v2d(nlon,nlat,nv2d) )
  WRITE(6,*) "Reading background/forecast data on model grid..."
  CALL read_diag(guesfile,v3d,v2d)

  !-----------------------------------------------------------------------------
  ! Generate model detph & related info
  !-----------------------------------------------------------------------------
  ! create 3d depth (note depth for some layers can have missing values)
  ALLOCATE( depth3d(nlon,nlat,nlev) )
  depth3d(:,:,1) = v3d(:,:,1,iv3d_h)
  DO k = 2, nlev
     !print*, "h:, k, min, max=", k, minval(v3d(:,:,k,iv3d_h)), maxval(v3d(:,:,k,iv3d_h),mask=v3d(:,:,k,iv3d_h)<10000.)
     depth3d(:,:,k) = depth3d(:,:,k-1) + v3d(:,:,k,iv3d_h) 
  ENDDO
  write(6,*) "max depth3d:, min, max=", minval(depth3d(:,:,nlev)), maxval(depth3d(:,:,nlev))

  ! create max valid model nlev based on h
  ALLOCATE( maxnlev(nlon,nlat) )
  ALLOCATE( maxdepth2d(nlon,nlat))
  maxnlev(:,:) = 0; maxdepth2d = 0.d0
  DO j = 1, nlat
     DO i = 1, nlon
        h_search: DO k = nlev, 1, -1
           if (v3d(i,j,k,iv3d_h) > 0.0011d0 .and. v3d(i,j,k,iv3d_h) < MAX_VALID_H  ) then
              maxnlev(i,j) = k
              exit h_search
           endif
        ENDDO h_search
        if (maxnlev(i,j)>=1) then
           maxdepth2d(i,j) = depth3d(i,j,maxnlev(i,j))
        else
           maxdepth2d(i,j) = 0.d0
        endif
     ENDDO
  ENDDO
  write(6,*) "maxnlev:, min, max=", minval(maxnlev), maxval(maxnlev)
  write(6,*) "maxdepth2d:, min, max=", minval(maxdepth2d), maxval(maxdepth2d)
  
  !-----------------------------------------------------------------------------
  ! Cycle through all profiles
  !-----------------------------------------------------------------------------
  cnt_coldt = 0
  ALLOCATE(depth1d(nlev+1))
  ALLOCATE(sim1d(nlev+1))
  WRITE(6,*) "Cycling through observations..."
  DO n=1, prof%nprof

    ! let simulation QC all be bad at first
    prof%simValid(:,n) = OCN_PROFILE_NOT_VALID

    ! Convert the physical coordinate to model grid coordinate (note: real, not integer)
    CALL phys2ij(olon(n),prof%lat(n),ri,rj) !(OCEAN)

    ! Filter in the tripolar region until localization is examined in the arctic !(ISSUE)
    if (DO_REMOVE_65N .and. prof%lat(n) > 65) then
      if (verbose) WRITE(6,'(A)') "Latitude above 65.0N, in tripolar region. Removing observation..."
      cnt_triout = cnt_triout + 1
      CYCLE
    endif
    
    ! Filter out observations in the black sea, which may cause problems
    if (DO_REMOVE_BLACKSEA) then
      if (40 < prof%lat(n) .and. prof%lat(n) < 50 .and. &
          25 < olon(n) .and. olon(n) < 45) then
        cnt_blacksea = cnt_blacksea + 1
        CYCLE
      endif
    endif

    ! Filter out observations that are out of range for the grid
    if (CEILING(ri) < 2 .OR. nlon+1 < CEILING(ri)) then
      if (verbose) WRITE(6,'(A)') '* X-coordinate out of range'
      if (verbose) WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', olon=', olon(n)
      cnt_xout = cnt_xout + 1
      CYCLE
    endif
    if (CEILING(rj) < 2 .OR. nlat < CEILING(rj)) then
      if (verbose) WRITE(6,'(A)') '* Y-coordinate out of range'
      if (verbose) WRITE(6,'(A,F6.2,A,F6.2)') '*   rj=',rj,', olat=',prof%lat(n)
      cnt_yout = cnt_yout + 1
      CYCLE
    endif

    ! Check the observation against boundaries
    !STEVE: Check to make sure it's in the ocean, as determined       (OCEAN)
    !       by mom4's topography map.
    ! (note: must do it after coordinate checks, or the coordinate
    !        could be outside of the range of the kmt array)
    boundary_points : if (ri > nlon) then
      !STEVE: I have to check what it does for this case...
      !       but it causes an error in the next line if ri > nlon
      if (verbose) WRITE(6,'(A)') '* coordinate is not on mom4 model grid: ri > nlon'
      cnt_rigtnlon = cnt_rigtnlon + 1
      if (cnt_rigtnlon <= 1) then
        WRITE(6,*) "STEVE: ri > nlon (cnt_rigtnlon)"
        WRITE(6,*) "ri = ", ri
        WRITE(6,*) "nlon = ", nlon
        WRITE(6,*) "rj = ", rj
        WRITE(6,*) "rlon(n) = ", olon(n)
        WRITE(6,*) "rlat(n) = ", prof%lat(n)
      endif
      CYCLE
    else
        ! ensure kmt >=1
        if(kmt(CEILING(ri),CEILING(rj)) .lt. 1 .or. &
             kmt(CEILING(ri),FLOOR(rj)) .lt. 1 .or. &
             kmt(FLOOR(ri),CEILING(rj)) .lt. 1 .or. &
             kmt(FLOOR(ri),FLOOR(rj)) .lt. 1) THEN

          if (debug_obsfilter) then
            WRITE(6,'(A)') '* coordinate is too close to land, according to kmt'
            WRITE(6,'(A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj
            WRITE(6,*) "kmt cell = ", kmt(FLOOR(ri):CEILING(ri),FLOOR(rj):CEILING(rj))
          endif
          cnt_nearland = cnt_nearland + 1
          cnt_nl1 = cnt_nl1+1
          CYCLE
        endif
        
        ! ensure all nearby grid points surrounding the obs are over ocean
        if(wet(CEILING(ri),CEILING(rj)) < 0.5 .or. &
           wet(CEILING(ri),FLOOR(rj))   < 0.5 .or. &
           wet(FLOOR(ri),CEILING(rj))   < 0.5 .or. &
           wet(FLOOR(ri),FLOOR(rj))     < 0.5 ) THEN ! ocn-1, land-0

          if (debug_obsfilter) then
            WRITE(6,'(A)') '* coordinate is too close to underwater topography, according to kmt'
            WRITE(6,'(A,F6.2,A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj, ', rk=',rk
            WRITE(6,*) "kmt cell = ", kmt(FLOOR(ri):CEILING(ri),FLOOR(rj):CEILING(rj))
          endif
          cnt_nearland = cnt_nearland + 1
          cnt_nl2 = cnt_nl2+1
          CYCLE
        endif

        ! ensure all nearby grid points have valid ocean depth
        if(maxnlev(CEILING(ri),CEILING(rj)) < 1  .or. &
           maxnlev(CEILING(ri),FLOOR(rj))   < 1  .or. &
           maxnlev(FLOOR(ri),CEILING(rj))   < 1  .or. &
           maxnlev(FLOOR(ri),FLOOR(rj))     < 1 ) THEN ! ocn-1, land-0

          if (debug_obsfilter) then
            WRITE(6,'(A)') '* coordinate is too close to underwater topography, according to kmt'
            WRITE(6,'(A,F6.2,A,F6.2,A,F6.2)') '*   ri=',ri,', rj=',rj, ', rk=',rk
            WRITE(6,*) "kmt cell = ", kmt(FLOOR(ri):CEILING(ri),FLOOR(rj):CEILING(rj))
          endif
          cnt_nearland = cnt_nearland + 1
          cnt_nl2 = cnt_nl2+1
          CYCLE
        endif
 
    endif boundary_points

   ! find the nearest grid point to the profile
   CALL phys2ij_nearest(olon(n),prof%lat(n),rii,rjj)
   ii = NINT(rii); jj = NINT(rjj)

   ! might be redundant with bydcheck but let's make sure
   IF (wet(ii,jj) < 0.5 .or. maxnlev(ii,jj) < 1)  CYCLE

   ! generate the model profile at this grid (add surface as additional 1 sea surface level)
   depth1d(2:nlev+1) = depth3d(ii,jj,1:nlev)
   sim1d(2:nlev+1)   = v3d(ii,jj,1:nlev,ivar)
   depth1d(1) = 0.d0; sim1d(1) = sim1d(2)

   ! model level used for interp (+1 due to additional sea surface)
   ims = 1; ime = maxnlev(ii,jj) + 1

   ! search the starting and ending depth index of the observed profile that model can simulate
   ios = 1; ioe = 1
   search_olev: do k = prof%nlevValid(n), 1, -1
      if (prof%depth(k,n) < maxdepth2d(ii,jj)) then
         ioe = k
         exit search_olev
      endif
   enddo search_olev

   ! debug
   if (debug_prof) then
      WRITE(6,*) "DEBUG profile idx: ", n
      WRITE(6,*) "olat,olon=", prof%lat(n), prof%lon(n), olon(n)
      WRITE(6,*) "mlat,mlon=", lat2d(ii,jj), lon2d(ii,jj)
      do k = 1, nlev+1
         print*, "model: k, dep1d, imelevel=", k, depth1d(k), sim1d(k), k==ime
      enddo
      print*, "ime-1, maxmlev, maxdepth, dep1d=", ime, maxnlev(ii,jj), maxdepth2d(ii,jj), depth1d(ime)
   
      print*, "ios, ioe=", ios, ioe
      do k = 1, prof%nlevValid(n)
         print*, "obs: k, depth=", k, prof%depth(k,n),prof%obsVal(k,n),ioe == k
      enddo
   endif

   ! interp for those calculable levels & assign QC flags
   prof%simValid(:,n) = OCN_PROFILE_NOT_VALID
   if (intp_opt == 1 ) then
   call spline_1d(    depth1d(ims:ime), &
                        sim1d(ims:ime),  &
                   prof%depth(ios:ioe,n), &
                  prof%simVal(ios:ioe,n) )
   elseif (intp_opt == 2) then
       call interp_1d(    depth1d(ims:ime), &
                            sim1d(ims:ime),  &
                       prof%depth(ios:ioe,n), &
                      prof%simVal(ios:ioe,n), &
                               ip )
   else
      WRITE(6,*) "[err] verify_en4.f90: unrecognized intp_opt=",intp_opt
      stop (13)
   endif
   do k = ios, ioe
      if (prof%simVal(k,n) > -4.d0 .and. prof%simVal(k,n) <=45.0) then
          ! both temp and sal should be within this valid range
          prof%simValid(k,n) = OCN_PROFILE_VALID
      endif
   enddo


   if (debug_prof) then
      do k = 1, prof%nlevValid(n)
         print*, "obs: n, k, depth, obs, oqc, model, mqc =", n, k, prof%depth(k,n),prof%obsVal(k,n), prof%obsValid(k,n), &
                 prof%simVal(k,n), prof%simValid(k,n)
      enddo
   endif

  if (debug_prof .and. ANY(prof%simVal(ios:ioe,n)<-2.0) ) then
     WRITE(6,*) "iterate to the next profile"

      do k = 1, prof%nlevValid(n)
         if (prof%simValid(k,n) == OCN_PROFILE_VALID.and.prof%simVal(k,n)<-2.d0) then
            cnt_coldt = cnt_coldt+1
         endif
      enddo

      do k = ims,ime
         write(6,*) "interp prof: n, k, h, sim=",n, k,depth1d(k),sim1d(k)
      enddo

     !if (n==228) pause 0
     !READ(*,*)
  endif
  enddo ! [n =1, prof%nprof]

  WRITE(*,*) "======n_cold =", cnt_coldt

  CALL prof%inspect()

  !-----------------------------------------------------------------------------
  ! Print out the counts of observations removed for various reasons
  !-----------------------------------------------------------------------------
  WRITE(6,*) "In obsop_tprof.f90:: observations removed for:"
  WRITE(6,*) "cnt_xout = ", cnt_xout
  WRITE(6,*) "cnt_yout = ", cnt_yout
  WRITE(6,*) "cnt_triout = ", cnt_triout
  WRITE(6,*) "cnt_rigtnlon = ", cnt_rigtnlon
  WRITE(6,*) "cnt_nearland = ", cnt_nearland
  WRITE(6,*) "cnt_nl1 = ", cnt_nl1
  WRITE(6,*) "cnt_nl2 = ", cnt_nl2
  if (DO_REMOVE_BLACKSEA) WRITE(6,*) "cnt_blacksea = ", cnt_blacksea

  !-----------------------------------------------------------------------------
  ! Write the observations and their associated innovations to file
  !-----------------------------------------------------------------------------

  CALL prof%write_nc(trim(obsoutfile))
  

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
    case('-ovar')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      ovar = arg2
    case('-intp')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      write(arg2,"(I10)") intp_opt
    case('-rm65N')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) DO_REMOVE_65N
    case('-rmBlackSea')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      PRINT *, "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) DO_REMOVE_BLACKSEA
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
