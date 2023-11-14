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
  
  PUBLIC :: read_en4_nc, argo_data
  PRIVATE :: inspect_obs_data
  
  
  INTEGER :: nobs
 
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
  CHARACTER(len=:), allocatable :: qc_pos(:)     ! slen=nprof, *1
  CHARACTER(len=:), allocatable :: qc_prof(:) ! slen=nprof, *1
  CHARACTER(len=:), allocatable :: qc_var(:)     ! slen=nlv, *nprof
  CHARACTER(64) :: varname, varqcname, profqcname
  integer,parameter :: days_between_19500101_19780101 = 10227
  type(t_datetime) :: datetime_19780101, qc_datetime, min_datetime, max_datetime
  type(t_datetime),allocatable :: obs_datetime(:) ! nprof
  integer :: min_total_seconds(1), max_total_seconds(1)
  real(r_dble),allocatable :: rdbuf1d(:)   ! nprof
  real(r_sngl),allocatable :: depth2d(:,:) ! (nlv,nprof)
  real(r_sngl),allocatable :: var2d(:,:) ! (nlv,nprof)
  real(r_dble),allocatable :: lat(:) ! nprof
  real(r_dble),allocatable :: lon(:) ! nprof
  logical,allocatable :: valid2d(:,:) ! overall qc (nlv,nprof)
  logical,allocatable :: valid(:)     ! overall profile qc (nprof)
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

  allocate(valid(nprof),valid2d(nlv,nprof))
  allocate(lon(nprof),lat(nprof))
  allocate(depth2d(nlv,nprof))
  allocate(var2d(nlv,nprof))
  allocate(character(len=nprof)::qc_pos(1))
  allocate(character(len=nprof)::qc_prof(1))
  allocate(character(len=nlv)::qc_var(nprof))
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

  allocate(rdbuf1d(nprof))
  allocate(obs_datetime(nprof))
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
  
  allocate(obs_data(nobs))
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
  
  ! Explicitly deallocate temporary arrays
  deallocate(stde,valid,valid2d,lon,lat,depth2d,var2d,obs_datetime,rdbuf1d)
  deallocate(qc_pos,qc_prof, qc_var)
  
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
