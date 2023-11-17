MODULE read_en4
!===============================================================================
! This program reads netcdf to convert profiles of temperature and salinity to 
! a format readable by letkf. The obsop.f90 program uses this module to read these
! in directly.
!
! Observation errors should be computed externally using Dave Behringer's technique based
! on the temperature gradients, or read from an observation error file
! (e.g. when doing adaptive obs error)
!
! Either the in situ temperatures will have to be converted to potential temperature,
! or the obs operator will have to transform the model potential temperature to
! in situ equivalent. The latter will be easier because there is no guarantee
! that temperature and salinity are observed simultaneously.
!
!===============================================================================
  
  USE common,                     ONLY: r_sngl, r_size, r_dble, slen
  USE params_obs,                 ONLY: id_t_obs, id_s_obs
  USE compute_profile_error,      ONLY: cmpTz
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: read_en4_nc, argo_data, inspect_obs_data
  PUBLIC :: ocn_profile
  PUBLIC :: select_prof_en4_nc
  

  INTEGER,PARAMETER :: OCN_PROFILE_VALID = 1
  INTEGER,PARAMETER :: OCN_PROFILE_NOT_VALID = 4

  TYPE ocn_profile
    INTEGER :: nprof = 0
    INTEGER :: nlev  = 0
    REAL(r_size),ALLOCATABLE :: obsVal(:,:) ! nlev * nprof
    REAL(r_size),ALLOCATABLE :: depth(:,:) ! nlev*nprof
    INTEGER,     ALLOCATABLE :: obsValid(:,:) ! nlev*nprof: 1-good, -9: bad
    REAL(r_size),ALLOCATABLE :: simVal(:,:)      ! nlev*nprof
    INTEGER,     ALLOCATABLE :: simValid(:,:) ! nlev*nprof:
    REAL(r_size),ALLOCATABLE :: lon(:) ! nprof
    REAL(r_size),ALLOCATABLE :: lat(:) ! nprof
    REAL(r_size),ALLOCATABLE :: time(:) ! nprof
    INTEGER,     ALLOCATABLE :: nlevValid(:) ! nprof

    CONTAINS
      PROCEDURE,PUBLIC,PASS(this) :: ALLOCATE => ALLOCATE_ocn_profile
      PROCEDURE,PUBLIC,PASS(this) :: DEALLOCATE => DEALLOCATE_ocn_profile
      PROCEDURE,PUBLIC,PASS(this) :: write_nc => write_nc_ocn_profile
      PROCEDURE,PUBLIC,PASS(this) :: read_nc  => read_nc_ocn_profile
      PROCEDURE,PUBLIC,PASS(this) :: inspect  => inspect_ocn_profile
  ENDTYPE
 
  TYPE argo_data
    REAL(r_size) :: x_grd(3)  ! longitude, latitude, and z depth (m)
    REAL(r_size) :: value     ! actual physical value of the parameter measured at this grid point
    REAL(r_size) :: lev       ! grid z level
    REAL(r_size) :: oerr      ! observation standard error
    REAL(r_size) :: hour      ! Hour of observation
    CHARACTER(9) :: plat      ! Platform
    CHARACTER(3) :: ptyp      ! Profile type
    CHARACTER(3) :: sid       ! Source id
    CHARACTER(1) :: qkey      ! Quality key
    INTEGER :: typ    ! observation variable type (e.g., PRES_TYPE)
    INTEGER :: nlevs  ! number of levels with data, counting from the top, including levels with missing data that have obs below them.
    INTEGER :: id     ! id number used in observation files to identify the observation
    INTEGER :: rid    ! id of the record, in order that it is read in
    INTEGER :: lid    ! id of the level for each record (upon skipping missing data for some levels)
    LOGICAL :: kept   ! tells letkf whether this obs is kept for assimilation
  END TYPE argo_data
  !TYPE(argo_data), ALLOCATABLE, DIMENSION(:) :: obs_data
  
  !! Read temp data
  !infile = trim(indir1)//'/'//YYYY//MM//DD//'tmp.nc' !STEVE: sample
  !typ = id_t_obs
  !CALL read_argo_nc(infile,typ,obs_data,nobs)
  !print *, "temp nobs = ", nobs
  
  !! Write letkf file
  !do i=1,nobs
  !!STEVE: the following are required for miyoshi's letkf observation input:
  !!1 = obelm
  !!2 = lon
  !!3 = lat
  !!4 = lev
  !!5 = value
  !!6 = oberr
  ! wk(1) = obs_data(i)%typ
  ! wk(2) = obs_data(i)%x_grd(1)
  ! wk(3) = obs_data(i)%x_grd(2)
  ! wk(4) = obs_data(i)%x_grd(3)
  ! wk(5) = obs_data(i)%value
  ! wk(6) = obs_data(i)%oerr
  ! WRITE(fid) wk
  !enddo
  
CONTAINS

SUBROUTINE allocate_ocn_profile(this,nlev,nprof)
  IMPLICIT NONE
  CLASS(ocn_profile),INTENT(INOUT) :: this
  INTEGER,INTENT(IN) :: nlev, nprof

  this % nlev = nlev
  this % nprof = nprof

  ALLOCATE( this % lon(nprof) )
  ALLOCATE( this % lat(nprof) )

  ALLOCATE( this % time(nprof) )

  ALLOCATE( this % nlevValid(nprof) )
  ALLOCATE( this % obsValid(nlev,nprof) )
  ALLOCATE( this % simValid(nlev,nprof) )

  ALLOCATE( this % depth(nlev,nprof) )
  ALLOCATE( this % obsVal(nlev,nprof) )
  ALLOCATE( this % simVal(nlev,nprof) )

  this % nlevValid(:) = 0
  this % obsValid(:,:)  = OCN_PROFILE_NOT_VALID
  this % simValid(:,:)  = OCN_PROFILE_NOT_VALID

ENDSUBROUTINE


SUBROUTINE deallocate_ocn_profile(this)
  IMPLICIT NONE
  CLASS(ocn_profile),INTENT(INOUT) :: this

  this % nlev = 0 
  this % nprof = 0

  if (ALLOCATED(this%lon)) DEALLOCATE( this%lon)
  if (ALLOCATED(this%lat)) DEALLOCATE( this%lat)
  if (ALLOCATED(this%time)) DEALLOCATE( this%time)
  if (ALLOCATED(this%nlevValid)) DEALLOCATE( this%nlevValid)
  if (ALLOCATED(this%obsVal)) DEALLOCATE( this%obsVal)
  if (ALLOCATED(this%obsValid)) DEALLOCATE( this%obsValid)
  if (ALLOCATED(this%simVal)) DEALLOCATE( this%simVal)
  if (ALLOCATED(this%simValid)) DEALLOCATE( this%simValid)
  if (ALLOCATED(this%depth)) DEALLOCATE( this%depth)

ENDSUBROUTINE

SUBROUTINE read_nc_ocn_profile(this,fnin)
  USE m_ncio,  ONLY: nc_get_fid, nc_close_fid, nc_rdvar1d, nc_rdvar2d, &
                     nc_rddim
  IMPLICIT NONE
  
  CLASS(ocn_profile),INTENT(INOUT) :: this
  CHARACTER(*),INTENT(IN) :: fnin

  INTEGER :: nlev, nprof
  INTEGER :: fid
  INTEGER,ALLOCATABLE :: i4buf1d(:), i4buf2d(:,:)
  REAL(8),ALLOCATABLE :: r8buf1d(:), r8buf2d(:,:)

!----------------------------------------------------------
! open the file & read in dim
  CALL nc_get_fid(trim(fnin),fid)
  CALL nc_rddim(fid,"nlev",nlev)
  CALL nc_rddim(fid,"nprof",nprof)
  CALL this%deallocate()
  CALL this%allocate(nprof=nprof,nlev=nlev)

!----------------------------------------------------------
! read in the vars
  ALLOCATE(i4buf1d(nprof),r8buf1d(nprof))
  ALLOCATE(i4buf2d(nlev,nprof),r8buf2d(nlev,nprof))
  CALL nc_rdvar1d(fid,"lon",this%lon)
  CALL nc_rdvar1d(fid,"lat",this%lat)
  CALL nc_rdvar1d(fid,"time",this%time)
  CALL nc_rdvar1d(fid,"nlevValid",this%nlevValid)
  CALL nc_rdvar2d(fid,"obsValid",this%obsValid)
  CALL nc_rdvar2d(fid,"simValid",this%simValid)
  CALL nc_rdvar2d(fid,"depth",   this%depth)
  CALL nc_rdvar2d(fid,"obsVal",  this%obsVal)
  CALL nc_rdvar2d(fid,"simVal",  this%simVal)

!----------------------------------------------------------
! close the file
  CALL nc_close_fid(fid)

  DEALLOCATE(i4buf1d,i4buf2d,r8buf1d,r8buf2d)

  CALL this%inspect()

ENDSUBROUTINE

SUBROUTINE inspect_ocn_profile(this)
  IMPLICIT NONE

  CLASS(ocn_profile),INTENT(IN) :: this
  
  write(6,*) "inspect ocn_profile:"
  write(6,*) "nlev, nprof = ", this%nlev, this%nprof
  if (allocated(this%lon)) then
     write(6,*) "lon: min, max=", minval(this%lon),maxval(this%lon)
  else
     write(6,*) "lon: not allocated"
  endif

  if (allocated(this%lat)) then
     write(6,*) "lat: min, max=", minval(this%lat),maxval(this%lat)
  else
     write(6,*) "lat: not allocated"
  endif

  if (allocated(this%time)) then
     write(6,*) "time: min, max=", minval(this%time),maxval(this%time)
  else
     write(6,*) "time: not allocated"
  endif

  if (allocated(this%nlevValid)) then
     write(6,*) "nlevValid: min, max=", minval(this%nlevValid),maxval(this%nlevValid)
  else
     write(6,*) "nlevValid: not allocated"
  endif

  if (allocated(this%obsValid)) then
     write(6,*) "obsValid: min, max=", minval(this%obsValid),maxval(this%obsValid)
  else
     write(6,*) "obsValid: not allocated"
  endif

  if (allocated(this%simValid)) then
     write(6,*) "simValid: min, max=", minval(this%simValid),maxval(this%simValid)
  else
     write(6,*) "simValid: not allocated"
  endif

  if (allocated(this%depth)) then
     write(6,*) "depth: min, max=", minval(this%depth,mask=this%obsValid==OCN_PROFILE_VALID), &
               maxval(this%depth,mask=this%obsValid==OCN_PROFILE_VALID)
  else
     write(6,*) "depth: not allocated"
  endif

  if (allocated(this%obsVal)) then
     write(6,*) "obsVal: min, max=", minval(this%obsVal,mask=this%obsValid==OCN_PROFILE_VALID), &
                maxval(this%obsVal,mask=this%obsValid==OCN_PROFILE_VALID)
  else
     write(6,*) "obsVal: not allocated"
  endif

  if (allocated(this%simVal)) then
     write(6,*) "simVal: min, max=", minval(this%simVal,mask=this%simValid==OCN_PROFILE_VALID),maxval(this%simVal,mask=this%simValid==OCN_PROFILE_VALID)
  else
     write(6,*) "simVal: not allocated"
  endif

ENDSUBROUTINE

SUBROUTINE write_nc_ocn_profile(this,fnout)
  USE m_ncio,    ONLY: nc_get_fid, nc_close_fid, &
                       nc_wrtvar1d, nc_wrtvar2d
  USE mod_f90gionc, ONLY: NC_createFile, NC_createDims, NC_createVars, &
                          NF90_INT, NF90_DOUBLE
  IMPLICIT NONE

  CLASS(ocn_profile),INTENT(IN) :: this
  CHARACTER(*),INTENT(IN) :: fnout
  ! vars for nc
  INTEGER,PARAMETER :: ndims = 2
  CHARACTER(80) :: dimname(ndims), dimlev, dimprof
  LOGICAL       :: dimunlimited(ndims)
  INTEGER       :: dimlen(ndims)

  INTEGER,PARAMETER :: nvars = 9
  CHARACTER(80) :: varname(nvars)
  INTEGER       :: varmaxdim(nvars) = 0
  CHARACTER(80) :: vardimname(4,nvars)
  INTEGER       :: vartype(nvars)

  INTEGER :: fid
  INTEGER :: ivar

  
  ! predefined
  dimname = ["nlev ", "nprof"] 
  dimlev  = dimname(1)
  dimprof = dimname(2)
  dimlen  = [ this%nlev, this%nprof ]
  dimunlimited = [.false., .true.]

  varname = ["lon      ", &  ! (nprof)
             "lat      ", &  ! (nprof)
             "time     ", &  ! (nprof)
             "nlevValid", &  ! (nprof)
             "obsValid ", &  ! (nlev,nprof)
             "simValid ", &  ! (nlev,nprof)
             "depth    ", &  ! (nlev,nprof)
             "obsVal   ", &  ! (nlev,nprof)
             "simVal   " ]   ! (nlev,nprof)
  vartype = [NF90_DOUBLE, & ! lon
             NF90_DOUBLE, & ! lat
             NF90_DOUBLE, & ! time
             NF90_INT, & ! nlevValid
             NF90_INT, & ! obsValid
             NF90_INT, & ! simValid
             NF90_DOUBLE, & ! depth
             NF90_DOUBLE, & ! obsVal
             NF90_DOUBLE ]  ! simVal
  varmaxdim(1:4) = 1
  varmaxdim(5:nvars) = 2
  do ivar = 1, 4
    vardimname(1:1,ivar) = [ dimprof ]
  enddo
  do ivar = 5, nvars
    vardimname(1:2,ivar) = [ dimlev,dimprof ]
  enddo

  ! nc create files
  CALL NC_CreateFile(trim(fnout))
  CALL NC_CreateDims(trim(fnout),ndims,dimname,dimunlimited,dimlen)
  call NC_CreateVars(trim(fnout),nvars,varname, vartype, varmaxdim, vardimname)

  ! write out profile
  CALL nc_get_fid(trim(fnout),fid,writemode=.true.)
  CALL nc_wrtvar1d(fid,varname(1),real(this%lon, 8)) ! lon
  CALL nc_wrtvar1d(fid,varname(2),real(this%lat, 8)) ! lat
  CALL nc_wrtvar1d(fid,varname(3),real(this%time,8)) ! time
  CALL nc_wrtvar1d(fid,varname(4),this%nlevValid) ! nlevValid
  CALL nc_wrtvar2d(fid,varname(5),this%obsValid) ! obsValid
  CALL nc_wrtvar2d(fid,varname(6),this%simValid) ! simValid
  CALL nc_wrtvar2d(fid,varname(7),real(this%depth, 8)) ! depth
  CALL nc_wrtvar2d(fid,varname(8),real(this%obsVal,8)) ! obsVal
  CALL nc_wrtvar2d(fid,varname(9),real(this%simVal,8)) ! simVal
  CALL nc_close_fid(fid)


END SUBROUTINE

SUBROUTINE select_prof_en4_nc(infile,typ,prof,Syyyymmddhh,delta_seconds)
  USE m_ncio,    ONLY: nc_get_fid, nc_close_fid, nc_rddim, nc_rdvar2d, nc_rdvar1d, &
                       nc_rdatt, nc_fndvar
  USE m_datetime, ONLY: t_datetime, t_timedelta, datetime, timedelta
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) :: infile
  INTEGER, INTENT(IN) :: typ
  CLASS(ocn_profile),INTENT(INOUT) :: prof
  CHARACTER(10),INTENT(IN), OPTIONAL :: Syyyymmddhh
                                        !1234567890
  INTEGER,      INTENT(IN), OPTIONAL :: delta_seconds  ! qc'ed if abs(obstime-Syyyymmddhh)>delta_seconds

! ncvars
  INTEGER :: fid
  INTEGER :: nprof, nlv
  INTEGER :: nprof_valid
  CHARACTER(len=14),dimension(1) :: srefdt = "19500101000000"
  CHARACTER(len=14),dimension(1),parameter :: srefdt_predefined = "19500101000000"
  CHARACTER(len=:), ALLOCATABLE :: qc_pos(:)     ! slen=nprof, *1
  CHARACTER(len=:), ALLOCATABLE :: qc_prof(:) ! slen=nprof, *1
  CHARACTER(len=:), ALLOCATABLE :: qc_var(:)     ! slen=nlv, *nprof
  CHARACTER(64) :: varname, varqcname, profqcname
  INTEGER,parameter :: days_between_19500101_19780101 = 10227
  type(t_datetime) :: datetime_19780101, qc_datetime, min_datetime, max_datetime
  type(t_datetime),ALLOCATABLE :: obs_datetime(:) ! nprof
  INTEGER :: min_total_seconds(1), max_total_seconds(1)
  real(r_dble),ALLOCATABLE :: rdbuf1d(:)   ! nprof
  real(r_sngl),ALLOCATABLE :: depth2d(:,:) ! (nlv,nprof)
  real(r_sngl),ALLOCATABLE :: var2d(:,:) ! (nlv,nprof)
  real(r_dble),ALLOCATABLE :: lat(:) ! nprof
  real(r_dble),ALLOCATABLE :: lon(:) ! nprof
  logical,ALLOCATABLE :: valid2d(:,:) ! overall qc (nlv,nprof)
  logical,ALLOCATABLE :: valid(:)     ! overall profile qc (nprof)
  real(r_dble) :: fillValue, valid_min, valid_max

  REAL(r_size) :: se0, seF

  INTEGER :: i,j,k,n
  REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: stde
  LOGICAL :: dodebug=.false.
  REAL(r_sngl) :: max_value=999
  REAL(r_sngl) :: max_depth=99999
  CHARACTER(*),parameter :: myname = "select_prof_en4_nc"

 
  
  !-------------------------------------------------------------------------------
  ! Open netcdf file & read dimension
  !-------------------------------------------------------------------------------
  CALL nc_get_fid(trim(infile),fid)
  CALL nc_rddim(fid, "N_PROF", nprof)
  CALL nc_rddim(fid, "N_LEVELS", nlv)
  write(6,*) "[msg] "//trim(myname)//"::fn, nprof, nlevels=",trim(infile),nprof,nlv

  ALLOCATE(valid(nprof),valid2d(nlv,nprof))
  ALLOCATE(lon(nprof),lat(nprof))
  ALLOCATE(depth2d(nlv,nprof))
  ALLOCATE(var2d(nlv,nprof))
  ALLOCATE(character(len=nprof)::qc_pos(1))
  ALLOCATE(character(len=nprof)::qc_prof(1))
  ALLOCATE(character(len=nlv)::qc_var(nprof))
  valid(:) = .true.
  valid2d(:,:) = .true.

  !-------------------------------------------------------------------------------
  ! read vars
  !-------------------------------------------------------------------------------
  ! spatial info
  CALL nc_rdvar1d(fid,"LONGITUDE", lon)
  CALL nc_rdatt(fid,"LONGITUDE", "valid_min",valid_min)
  CALL nc_rdatt(fid,"LONGITUDE", "valid_max",valid_max)
  CALL nc_rdatt(fid,"LONGITUDE", "_fillvalue",fillValue)
  print*, "valid_min, valid_max, fillValue=",valid_min, valid_max, fillValue
  where (lon<valid_min .or. lon>valid_max .or. nint(lon)==nint(fillValue))
       valid = .false.
  endwhere

  CALL nc_rdvar1d(fid,"LATITUDE",  lat)
  CALL nc_rdatt(fid,"LATITUDE", "valid_min",valid_min)
  CALL nc_rdatt(fid,"LATITUDE", "valid_max",valid_max)
  CALL nc_rdatt(fid,"LATITUDE", "_fillvalue",fillValue)
  print*, "valid_min, valid_max, fillValue=",valid_min, valid_max, fillValue
  where (lat<valid_min .or. lat>valid_max .or. nint(lat)==nint(fillValue))
       valid = .false.
  endwhere

  CALL nc_rdvar1d(fid,"POSITION_QC",nprof,qc_pos)
  do i = 1, nprof
     if (qc_pos(1)(i:i) /= '1') valid(i) = .false.
  enddo
  write(6,*) "lon: min, max=", minval(lon,mask=valid),maxval(lon,mask=valid)
  write(6,*) "lat: min, max=", minval(lat,mask=valid),maxval(lat,mask=valid)

  nprof_valid = 0
  do i = 1, nprof
     nprof_valid = nprof_valid +1
  enddo
  print*, "[debug] nprof_valid =",nprof_valid
  
  ! time info
  CALL nc_rdvar1d(fid,"REFERENCE_DATE_TIME", 14,srefdt)
  write(6,*) "reference date in the file = ", trim(srefdt(1))
  if (srefdt(1) /= srefdt_predefined(1)) then
     write(6,*) "[err] "//trim(myname)//":: reference time in the file not the same as predefined value (", &
                trim(srefdt_predefined(1)), ")"
     stop (2)
  end if
  datetime_19780101 = datetime(1978,1,1,0,0,0) ! this is the earliest date supported by w3movdat_full.f

  ALLOCATE(rdbuf1d(nprof))
  ALLOCATE(obs_datetime(nprof))
  CALL nc_rdvar1d(fid,"JULD",rdbuf1d)
  do i = 1, nprof
     ! this is because the w3movdate_full.f does not support date before 19780101
     obs_datetime(i) = datetime_19780101 + timedelta(seconds=nint((rdbuf1d(i)-days_between_19500101_19780101)*24*3600))
  enddo
  ! print out min & max obs time from the file
  do i = 1, nprof
     rdbuf1d(i) = obs_datetime(i)%total_seconds_since19780101()
  enddo
  min_total_seconds = minval(rdbuf1d,mask=valid)
  max_total_seconds = maxval(rdbuf1d,mask=valid)
  min_datetime = datetime_19780101 + timedelta(seconds=min_total_seconds(1))
  max_datetime = datetime_19780101 + timedelta(seconds=max_total_seconds(1))
  write(6,*) "obs_start_time, obs_end_time=", min_datetime%string(), max_datetime%string()

  ! time filtering
  if (present(Syyyymmddhh) .and. present(delta_seconds)) then
     read(Syyyymmddhh,"(I4,I2,I2,I2)") qc_datetime%year, qc_datetime%month, &
                                       qc_datetime%day,  qc_datetime%hour
     qc_datetime%minute = 0; qc_datetime%second = 0
     write(6,*) "qc_datetime, delta_seconds=", qc_datetime%string(), delta_seconds
     nprof_valid = 0
     do i = 1, nprof
        rdbuf1d(i) = obs_datetime(i)%total_seconds_since19780101()
        if (abs(qc_datetime%total_seconds_since19780101() - &
                obs_datetime(i)%total_seconds_since19780101()) > delta_seconds) then
           valid(i) = .false.  ! (nlv,nprof)
        end if
        if (valid(i)) nprof_valid = nprof_valid + 1
     enddo
     min_total_seconds = minval(rdbuf1d,mask=valid)
     max_total_seconds = maxval(rdbuf1d,mask=valid)
     min_datetime = datetime_19780101 + timedelta(seconds=min_total_seconds(1))
     max_datetime = datetime_19780101 + timedelta(seconds=max_total_seconds(1))
     write(6,*) "min_time, max_time=", min_datetime%string(), max_datetime%string()
     write(6,*) "nprof_valid=",nprof_valid
  endif

  ! read measured quantity
  if (typ == id_t_obs) then
     ! read in potential temperature and qc flags for profiles & all levels
     varname       = "POTM_CORRECTED"
     varqcname     = "POTM_CORRECTED_QC"
     profqcname    = "PROFILE_POTM_QC"
  elseif (typ == id_s_obs) then
     varname       = "PSAL_CORRECTED"
     varqcname     = "PSAL_CORRECTED_QC"
     profqcname    = "PROFILE_PSAL_QC"
  else
     write(6,*) "[err] "//trim(myname)//":: unrecognized typ=",typ, ", only support ",id_t_obs, "or", id_s_obs
     stop (8)
  endif
  write(6,*) "varname, varqcname, profqcname=",trim(varname),trim(varqcname),trim(profqcname)
  CALL nc_rdvar1d(fid,trim(profqcname),nprof,qc_prof)
  nprof_valid = 0
  do i = 1, nprof
     if (qc_prof(1)(i:i) /= '1') valid(i) = .false.
     if (valid(i)) nprof_valid = nprof_valid + 1
  enddo
  write(6,*) "nprof_valid = ", nprof_valid

  ! transfer profile QC to all levels
  valid2d(:,:) = .true.
  do i = 1, nprof
     if (.not.valid(i)) valid2d(:,i) = .false.
  enddo

  ! read depth
  CALL nc_rdvar2d(fid,"DEPH_CORRECTED",depth2d)
  CALL nc_rdatt(  fid,"DEPH_CORRECTED","valid_min",valid_min)
  CALL nc_rdatt(  fid,"DEPH_CORRECTED","valid_max",valid_max)
  CALL nc_rdatt(  fid,"DEPH_CORRECTED","_fillvalue",fillValue)
  write(6,*) "DEPH_CORRECTED: valid_min, valid_max, _fillvalue=", valid_min, valid_max, fillValue
  where (depth2d<valid_min .or. depth2d>valid_max .or. nint(depth2d)==nint(fillValue))
    valid2d = .false.
  endwhere
  if (dodebug) then
     write(6,*) "depth(1) = ", depth2d(:,1)
     write(6,*) "depth(2) = ", depth2d(:,2)
  endif

  ! read profile vars
  CALL nc_rdvar2d(fid,trim(varname),var2d)
  CALL nc_rdatt(  fid,trim(varname),"valid_min",valid_min)
  CALL nc_rdatt(  fid,trim(varname),"valid_max",valid_max)
  CALL nc_rdatt(  fid,trim(varname),"_fillvalue",fillValue)
  write(6,*) "var: valid_min, valid_max, _fillvalue=", valid_min, valid_max, fillValue
  where (var2d<valid_min .or. var2d>valid_max .or. nint(var2d)==nint(fillValue))
    valid2d = .false.
  endwhere
  write(6,*) "var: min, max=", minval(var2d,mask=valid2d), maxval(var2d,mask=valid2d)
  !print*, "var(2)=", var2d(:,2)

  ! read profile var qc
  CALL nc_rdvar1d(fid,trim(varqcname),nlv,qc_var)
  do i = 1, nprof
     do k = 1, nlv
        if (qc_var(i)(k:k) /= "1") valid2d(k,i) = .false.
     enddo
  enddo
  write(6,*) "var: min, max=", minval(var2d,mask=valid2d), maxval(var2d,mask=valid2d)


  !-------------------------------------------------------------------------------
  ! Close the netcdf file
  !-------------------------------------------------------------------------------
  CALL nc_close_fid(fid)


  !-------------------------------------------------------------------------------
  ! pack valid profile into ocn_profile structure
  !-------------------------------------------------------------------------------
  nprof_valid = 0
  do i = 1, nprof
     if (valid(i)) nprof_valid = nprof_valid + 1
  enddo
  write(6,*) "select nprof = ", nprof_valid
  call prof%allocate(nprof = nprof_valid, nlev = nlv)
 
  j = 0
  do i = 1, nprof
     if (.not.valid(i)) cycle
     j = j + 1  ! prof index in output
     prof % lon(j)       = lon(i)
     prof % lat(j)       = lat(i)
     prof % time(j)      = obs_datetime(i)%total_seconds_since19780101()/3600.d0
     prof % nlevValid(j) = 0
     n = 0
     do k = 1, nlv
        if (.not.valid2d(k,i)) cycle
        n  = n + 1
        prof % depth(n,j)    = depth2d(k,i)
        prof % obsVal(n,j)   = var2d(k,i)
        prof % obsValid(n,j) = OCN_PROFILE_VALID
     enddo
     prof % nlevValid(j) = n
     write(*,*) "finish writing prof: idx, nlev=", j, prof % nlevValid(j)
  enddo
  if (j/=nprof_valid) then
     write(*,*) "[err] "//trim(myname)//": error in prof indexing"
     stop (11)
  else
     ! check if depth for each profile monotonically increases
     do j = 1, prof % nprof
        if (prof%nlevValid(j) > 1) then
            do n = 2, prof % nlevValid(j)
               if (prof % depth(n,j) < prof % depth(n-1,j)) then
                  write(*,*) "[err] "//trim(myname)//": profile depth does not increase monotonically"
                  write(*,*) "prof: idx, depth(:)=", j, prof%depth(1:prof%nlevValid(j),j)
                  stop (12)
               endif
            enddo
        endif
     enddo
  endif


  call prof%inspect()

ENDSUBROUTINE
 



SUBROUTINE read_en4_nc(infile,typ,obs_data,nobs,Syyyymmddhh,delta_seconds)
!===============================================================================
! Read the argo profile data
!===============================================================================
  USE m_ncio,    ONLY: nc_get_fid, nc_close_fid, nc_rddim, nc_rdvar2d, nc_rdvar1d, &
                       nc_rdatt, nc_fndvar
  USE m_datetime, ONLY: t_datetime, t_timedelta, datetime, timedelta
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) :: infile
  INTEGER, INTENT(IN) :: typ
  TYPE(argo_data), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: obs_data
  INTEGER, INTENT(OUT) :: nobs
  CHARACTER(10),INTENT(IN), OPTIONAL :: Syyyymmddhh
                                        !1234567890
  INTEGER,      INTENT(IN), OPTIONAL :: delta_seconds  ! qc'ed if abs(obstime-Syyyymmddhh)>delta_seconds

! ncvars
  INTEGER :: fid
  INTEGER :: nprof, nlv
  INTEGER :: nprof_valid
  CHARACTER(len=14),dimension(1) :: srefdt = "19500101000000"
  CHARACTER(len=14),dimension(1),parameter :: srefdt_predefined = "19500101000000"
  CHARACTER(len=:), ALLOCATABLE :: qc_pos(:)     ! slen=nprof, *1
  CHARACTER(len=:), ALLOCATABLE :: qc_prof(:) ! slen=nprof, *1
  CHARACTER(len=:), ALLOCATABLE :: qc_var(:)     ! slen=nlv, *nprof
  CHARACTER(64) :: varname, varqcname, profqcname
  INTEGER,parameter :: days_between_19500101_19780101 = 10227
  type(t_datetime) :: datetime_19780101, qc_datetime, min_datetime, max_datetime
  type(t_datetime),ALLOCATABLE :: obs_datetime(:) ! nprof
  INTEGER :: min_total_seconds(1), max_total_seconds(1)
  real(r_dble),ALLOCATABLE :: rdbuf1d(:)   ! nprof
  real(r_sngl),ALLOCATABLE :: depth2d(:,:) ! (nlv,nprof)
  real(r_sngl),ALLOCATABLE :: var2d(:,:) ! (nlv,nprof)
  real(r_dble),ALLOCATABLE :: lat(:) ! nprof
  real(r_dble),ALLOCATABLE :: lon(:) ! nprof
  logical,ALLOCATABLE :: valid2d(:,:) ! overall qc (nlv,nprof)
  logical,ALLOCATABLE :: valid(:)     ! overall profile qc (nprof)
  real(r_dble) :: fillValue, valid_min, valid_max

  REAL(r_size) :: se0, seF

  INTEGER :: i,j,k,n
  REAL(r_size), ALLOCATABLE, DIMENSION(:,:) :: stde
  LOGICAL :: dodebug=.false.
  REAL(r_sngl) :: max_value=999
  REAL(r_sngl) :: max_depth=99999
  CHARACTER(*),parameter :: myname = "read_en4_nc"

 
  
  !-------------------------------------------------------------------------------
  ! Open netcdf file & read dimension
  !-------------------------------------------------------------------------------
  CALL nc_get_fid(trim(infile),fid)
  CALL nc_rddim(fid, "N_PROF", nprof)
  CALL nc_rddim(fid, "N_LEVELS", nlv)
  write(6,*) "[msg] "//trim(myname)//"::fn, nprof, nlevels=",trim(infile),nprof,nlv

  ALLOCATE(valid(nprof),valid2d(nlv,nprof))
  ALLOCATE(lon(nprof),lat(nprof))
  ALLOCATE(depth2d(nlv,nprof))
  ALLOCATE(var2d(nlv,nprof))
  ALLOCATE(character(len=nprof)::qc_pos(1))
  ALLOCATE(character(len=nprof)::qc_prof(1))
  ALLOCATE(character(len=nlv)::qc_var(nprof))
  valid(:) = .true.
  valid2d(:,:) = .true.

  !-------------------------------------------------------------------------------
  ! read vars
  !-------------------------------------------------------------------------------
  ! spatial info
  CALL nc_rdvar1d(fid,"LONGITUDE", lon)
  CALL nc_rdatt(fid,"LONGITUDE", "valid_min",valid_min)
  CALL nc_rdatt(fid,"LONGITUDE", "valid_max",valid_max)
  CALL nc_rdatt(fid,"LONGITUDE", "_fillvalue",fillValue)
  print*, "valid_min, valid_max, fillValue=",valid_min, valid_max, fillValue
  where (lon<valid_min .or. lon>valid_max .or. nint(lon)==nint(fillValue))
       valid = .false.
  endwhere

  CALL nc_rdvar1d(fid,"LATITUDE",  lat)
  CALL nc_rdatt(fid,"LATITUDE", "valid_min",valid_min)
  CALL nc_rdatt(fid,"LATITUDE", "valid_max",valid_max)
  CALL nc_rdatt(fid,"LATITUDE", "_fillvalue",fillValue)
  print*, "valid_min, valid_max, fillValue=",valid_min, valid_max, fillValue
  where (lat<valid_min .or. lat>valid_max .or. nint(lat)==nint(fillValue))
       valid = .false.
  endwhere

  CALL nc_rdvar1d(fid,"POSITION_QC",nprof,qc_pos)
  do i = 1, nprof
     if (qc_pos(1)(i:i) /= '1') valid(i) = .false.
  enddo
  write(6,*) "lon: min, max=", minval(lon,mask=valid),maxval(lon,mask=valid)
  write(6,*) "lat: min, max=", minval(lat,mask=valid),maxval(lat,mask=valid)

  nprof_valid = 0
  do i = 1, nprof
     nprof_valid = nprof_valid +1
  enddo
  print*, "[debug] nprof_valid =",nprof_valid
  
  ! time info
  CALL nc_rdvar1d(fid,"REFERENCE_DATE_TIME", 14,srefdt)
  write(6,*) "reference date in the file = ", trim(srefdt(1))
  if (srefdt(1) /= srefdt_predefined(1)) then
     write(6,*) "[err] "//trim(myname)//":: reference time in the file not the same as predefined value (", &
                trim(srefdt_predefined(1)), ")"
     stop (2)
  end if
  datetime_19780101 = datetime(1978,1,1,0,0,0) ! this is the earliest date supported by w3movdat_full.f

  ALLOCATE(rdbuf1d(nprof))
  ALLOCATE(obs_datetime(nprof))
  CALL nc_rdvar1d(fid,"JULD",rdbuf1d)
  do i = 1, nprof
     ! this is because the w3movdate_full.f does not support date before 19780101
     obs_datetime(i) = datetime_19780101 + timedelta(seconds=nint((rdbuf1d(i)-days_between_19500101_19780101)*24*3600))
  enddo
  ! print out min & max obs time from the file
  do i = 1, nprof
     rdbuf1d(i) = obs_datetime(i)%total_seconds_since19780101()
  enddo
  min_total_seconds = minval(rdbuf1d,mask=valid)
  max_total_seconds = maxval(rdbuf1d,mask=valid)
  min_datetime = datetime_19780101 + timedelta(seconds=min_total_seconds(1))
  max_datetime = datetime_19780101 + timedelta(seconds=max_total_seconds(1))
  write(6,*) "obs_start_time, obs_end_time=", min_datetime%string(), max_datetime%string()

  ! time filtering
  if (present(Syyyymmddhh) .and. present(delta_seconds)) then
     read(Syyyymmddhh,"(I4,I2,I2,I2)") qc_datetime%year, qc_datetime%month, &
                                       qc_datetime%day,  qc_datetime%hour
     qc_datetime%minute = 0; qc_datetime%second = 0
     write(6,*) "qc_datetime, delta_seconds=", qc_datetime%string(), delta_seconds
     nprof_valid = 0
     do i = 1, nprof
        rdbuf1d(i) = obs_datetime(i)%total_seconds_since19780101()
        if (abs(qc_datetime%total_seconds_since19780101() - &
                obs_datetime(i)%total_seconds_since19780101()) > delta_seconds) then
           valid(i) = .false.  ! (nlv,nprof)
        end if
        if (valid(i)) nprof_valid = nprof_valid + 1
     enddo
     min_total_seconds = minval(rdbuf1d,mask=valid)
     max_total_seconds = maxval(rdbuf1d,mask=valid)
     min_datetime = datetime_19780101 + timedelta(seconds=min_total_seconds(1))
     max_datetime = datetime_19780101 + timedelta(seconds=max_total_seconds(1))
     write(6,*) "min_time, max_time=", min_datetime%string(), max_datetime%string()
     write(6,*) "nprof_valid=",nprof_valid
  endif

  ! read measured quantity
  if (typ == id_t_obs) then
     ! read in potential temperature and qc flags for profiles & all levels
     varname       = "POTM_CORRECTED"
     varqcname     = "POTM_CORRECTED_QC"
     profqcname    = "PROFILE_POTM_QC"
  elseif (typ == id_s_obs) then
     varname       = "PSAL_CORRECTED"
     varqcname     = "PSAL_CORRECTED_QC"
     profqcname    = "PROFILE_PSAL_QC"
  else
     write(6,*) "[err] "//trim(myname)//":: unrecognized typ=",typ, ", only support ",id_t_obs, "or", id_s_obs
     stop (8)
  endif
  write(6,*) "varname, varqcname, profqcname=",trim(varname),trim(varqcname),trim(profqcname)
  CALL nc_rdvar1d(fid,trim(profqcname),nprof,qc_prof)
  nprof_valid = 0
  do i = 1, nprof
     if (qc_prof(1)(i:i) /= '1') valid(i) = .false.
     if (valid(i)) nprof_valid = nprof_valid + 1
  enddo
  write(6,*) "nprof_valid = ", nprof_valid

  ! transfer profile QC to all levels
  valid2d(:,:) = .true.
  do i = 1, nprof
     if (.not.valid(i)) valid2d(:,i) = .false.
  enddo

  ! read depth
  CALL nc_rdvar2d(fid,"DEPH_CORRECTED",depth2d)
  CALL nc_rdatt(  fid,"DEPH_CORRECTED","valid_min",valid_min)
  CALL nc_rdatt(  fid,"DEPH_CORRECTED","valid_max",valid_max)
  CALL nc_rdatt(  fid,"DEPH_CORRECTED","_fillvalue",fillValue)
  write(6,*) "DEPH_CORRECTED: valid_min, valid_max, _fillvalue=", valid_min, valid_max, fillValue
  where (depth2d<valid_min .or. depth2d>valid_max .or. nint(depth2d)==nint(fillValue))
    valid2d = .false.
  endwhere
  if (dodebug) then
     write(6,*) "depth(1) = ", depth2d(:,1)
     write(6,*) "depth(2) = ", depth2d(:,2)
  endif

  ! read profile vars
  CALL nc_rdvar2d(fid,trim(varname),var2d)
  CALL nc_rdatt(  fid,trim(varname),"valid_min",valid_min)
  CALL nc_rdatt(  fid,trim(varname),"valid_max",valid_max)
  CALL nc_rdatt(  fid,trim(varname),"_fillvalue",fillValue)
  write(6,*) "var: valid_min, valid_max, _fillvalue=", valid_min, valid_max, fillValue
  where (var2d<valid_min .or. var2d>valid_max .or. nint(var2d)==nint(fillValue))
    valid2d = .false.
  endwhere
  write(6,*) "var: min, max=", minval(var2d,mask=valid2d), maxval(var2d,mask=valid2d)
  !print*, "var(2)=", var2d(:,2)

  ! read profile var qc
  CALL nc_rdvar1d(fid,trim(varqcname),nlv,qc_var)
  do i = 1, nprof
     do k = 1, nlv
        if (qc_var(i)(k:k) /= "1") valid2d(k,i) = .false.
     enddo
  enddo
  write(6,*) "var: min, max=", minval(var2d,mask=valid2d), maxval(var2d,mask=valid2d)


  !-------------------------------------------------------------------------------
  ! Close the netcdf file
  !-------------------------------------------------------------------------------
  CALL nc_close_fid(fid)
 
!  
  !-------------------------------------------------------------------------------
  ! STEVE: need to compute the standard error as Dave did based on vertical temperature gradient
  !-------------------------------------------------------------------------------
  ALLOCATE(stde(nlv,nprof))
  stde=0.0 !STEVE: placeholder for now. Perhaps better to do this in the calling function (e.g. obsop_tprof.f90)
  if (typ .eq. id_t_obs) then
    se0 = 1.0
    seF = 1.5
  elseif (typ .eq. id_s_obs) then
    se0 = 0.05
    seF = 0.15
  endif
  do i = 1, nprof
     CALL cmpTz(stde(:,i),se0,seF,real(var2d(:,i),r_size),real(depth2d(:,i),r_size),nlv,fillValue)
     if (dodebug) write(6,*) "[msg] "//trim(myname)//":: iprof,stde(:,iprof) = ", stde(:,i)
  enddo

  !--------------------------------------------------------------------------------
  ! pack good obs into structure
  !--------------------------------------------------------------------------------
  nobs = 0
  nprof_valid = 0
  do i = 1, nprof
     if (valid(i)) nprof_valid = nprof_valid + 1
     do k = 1, nlv
        if (valid2d(k,i)) nobs = nobs + 1
     enddo
  enddo
  write(6,*) "[msg] "//trim(myname)//":: nprof_valid, nobs=", nprof_valid, nobs
  
  ALLOCATE(obs_data(nobs))
  n = 0
  do i=1,nprof
    do k=1,nlv
      !if (val < max_value .and. depth(k) < max_depth .and. abs(val - missing_value) > 1.0) then
      !if (var2d(k,i) < max_value .and. depth2d(k,i) < max_depth ) then
      if (valid2d(k,i)) then
        n = n+1
        obs_data(n)%typ      = typ
        obs_data(n)%x_grd(1) = lon(i)
        obs_data(n)%x_grd(2) = lat(i)
        obs_data(n)%x_grd(3) = depth2d(k,i)
        obs_data(n)%hour     = obs_datetime(i)%total_seconds_since19780101()/3600.d0
        obs_data(n)%value    = var2d(k,i)
        obs_data(n)%oerr     = stde(k,i)
        obs_data(n)%rid      = i         ! record id
        obs_data(n)%lid      = k         ! level id
        obs_data(n)%qkey     = "1"
      endif
    enddo
  enddo

  call inspect_obs_data(obs_data,trim(myname))
  
  ! Explicitly DEALLOCATE temporary arrays
  DEALLOCATE(stde,valid,valid2d,lon,lat,depth2d,var2d,obs_datetime,rdbuf1d)
  DEALLOCATE(qc_pos,qc_prof, qc_var)
  
END SUBROUTINE read_en4_nc

SUBROUTINE inspect_obs_data(obs_data, subname)
  IMPLICIT NONE
  TYPE(argo_data),INTENT(IN) :: obs_data(:)
  CHARACTER(*),INTENT(IN) :: subname

  WRITE(6,*) "[msg] "//trim(subname)//"::info"
  WRITE(6,*) "                nobs=", size(obs_data)
  WRITE(6,*) "       typ: min, max=", minval(obs_data(:)%typ), maxval(obs_data(:)%typ)
  WRITE(6,*) "  x_grd(1): min, max=", minval(obs_data(:)%x_grd(1)), maxval(obs_data(:)%x_grd(1))
  WRITE(6,*) "  x_grd(2): min, max=", minval(obs_data(:)%x_grd(2)), maxval(obs_data(:)%x_grd(2))
  WRITE(6,*) "  x_grd(3): min, max=", minval(obs_data(:)%x_grd(3)), maxval(obs_data(:)%x_grd(3))
  WRITE(6,*) "      hour: min, max=", minval(obs_data(:)%hour),  maxval(obs_data(:)%hour)
  WRITE(6,*) "     value: min, max=", minval(obs_data(:)%value), maxval(obs_data(:)%value)
  WRITE(6,*) "      oerr: min, max=", minval(obs_data(:)%oerr),  maxval(obs_data(:)%oerr)
  WRITE(6,*) "       rid: min, max=", minval(obs_data(:)%rid),   maxval(obs_data(:)%rid)
  WRITE(6,*) "       lid: min, max=", minval(obs_data(:)%lid),   maxval(obs_data(:)%lid)

END SUBROUTINE

END MODULE read_en4
