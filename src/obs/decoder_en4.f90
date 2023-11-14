PROGRAM decoder_en4
  USE common
  USE params_model
  USE params_obs, only : id_s_obs, id_t_obs
  USE vars_model
  USE common_oceanmodel
  !USE vars_obs
  USE common_obs_oceanmodel
  USE read_en4,      ONLY: argo_data, read_en4_nc

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! Command line inputs:
  !-----------------------------------------------------------------------------
  CHARACTER(slen) :: obsinfile='obsin.nc'     !IN (default)
  CHARACTER(slen) :: obsoutfile='obsout.dat'  !OUT(default)
  CHARACTER(slen) :: ovar = "potm" ! IN (default), or "psal"
  INTEGER      :: otyp = id_t_obs

  !-----------------------------------------------------------------------------
  ! Obs data arrays
  !-----------------------------------------------------------------------------
  TYPE(argo_data), ALLOCATABLE :: obs_data(:)
  INTEGER :: nobs
  CHARACTER(10) :: Syyyymmddhh = "YYYYMMDDHH"
  !                 1234567890
  INTEGER :: delta_seconds = -999 ! max difference from Syyymmddhh

  REAL(r_size), ALLOCATABLE :: elem(:)
  REAL(r_size), ALLOCATABLE :: rlon(:)
  REAL(r_size), ALLOCATABLE :: rlat(:)
  REAL(r_size), ALLOCATABLE :: rlev(:)
  REAL(r_size), ALLOCATABLE :: odat(:)
  REAL(r_size), ALLOCATABLE :: oerr(:)
  REAL(r_size), ALLOCATABLE :: ohx(:)
  REAL(r_size), ALLOCATABLE :: obhr(:)
  INTEGER     , ALLOCATABLE :: oqc(:)

  !-----------------------------------------------------------------------------
  ! Miscellaneous
  !-----------------------------------------------------------------------------
  INTEGER :: n, i

  !-----------------------------------------------------------------------------
  ! Debugging parameters
  !-----------------------------------------------------------------------------
  !STEVE: to adjust writing to output file
  LOGICAL :: print1st = .true.

  CHARACTER(*),PARAMETER :: myname = "decoder_en4"

  !-----------------------------------------------------------------------------
  ! Initialize the common_oceanmodel module, and process command line options
  !-----------------------------------------------------------------------------
  CALL process_command_line !(get: -obsin <obsinfile> -gues <guesfile> -obsout <obsoutfile>)


  !-----------------------------------------------------------------------------
  ! Read observations from file
  !-----------------------------------------------------------------------------
  SELECT CASE(trim(ovar)) 
    CASE("potm")
      otyp = id_t_obs
    CASE ("psal")
      otyp = id_s_obs
    CASE DEFAULT
      WRITE(6,*) "[err] "//trim(myname)//"::Unsupported ovar=", trim(ovar), &
                 ". ovar supported should either be ptemp or psal"
      STOP (10)
  END SELECT
 
  WRITE(6,*) "[msg] "//trim(myname)//"::ovar, otyp=", trim(ovar),otyp
  if (delta_seconds>0) then
     CALL read_en4_nc(trim(obsinfile),otyp,obs_data,nobs,Syyyymmddhh, delta_seconds)
  else
     CALL read_en4_nc(trim(obsinfile),otyp,obs_data,nobs)
  endif

  ALLOCATE( elem(nobs) )
  ALLOCATE( rlon(nobs) )
  ALLOCATE( rlat(nobs) )
  ALLOCATE( rlev(nobs) )
  ALLOCATE( odat(nobs) )
  ALLOCATE( oerr(nobs) )
  ALLOCATE( ohx(nobs) )
  ALLOCATE( oqc(nobs) )
  ALLOCATE( obhr(nobs) )

  print *, "[msg] "//trim(myname)//":: starting nobs = ", nobs
  do i=1,nobs
    elem(i) = obs_data(i)%typ
    rlon(i) = obs_data(i)%x_grd(1)
    rlat(i) = obs_data(i)%x_grd(2)
    rlev(i) = obs_data(i)%x_grd(3)
    odat(i) = obs_data(i)%value
    oerr(i) = obs_data(i)%oerr
    ohx(i)  = 0
    oqc(i)  = 1
    obhr(i) = obs_data(i)%hour
  enddo
  DEALLOCATE(obs_data)

  CALL inspect_letkf_obs()

  if (print1st) then  
    i=1
    print *, "i, elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr = ", i, elem(i),rlon(i),rlat(i),rlev(i),odat(i),oerr(i),ohx(i),oqc(i),obhr(i)
    i=nobs
    print *, "i, elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr = ", i, elem(i),rlon(i),rlat(i),rlev(i),odat(i),oerr(i),ohx(i),oqc(i),obhr(i)
  endif

  call write_obs3(trim(obsoutfile),nobs,elem,rlon,rlat,rlev, &
                  odat,oerr,obhr,oqc,qcflag_in=.true.)
 
  DEALLOCATE( elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr )

CONTAINS

SUBROUTINE inspect_letkf_obs()
  IMPLICIT NONE

  WRITE(6,*) "nobs = ", size(elem)
  WRITE(6,*) "elem: min, max=", minval(elem), maxval(elem)
  WRITE(6,*) "rlon: min, max=", minval(rlon), maxval(rlon)
  WRITE(6,*) "rlat: min, max=", minval(rlat), maxval(rlat)
  WRITE(6,*) "rlev: min, max=", minval(rlev), maxval(rlev)
  WRITE(6,*) "odat: min, max=", minval(odat), maxval(odat)
  WRITE(6,*) "oerr: min, max=", minval(oerr), maxval(oerr)
  WRITE(6,*) "ohx:  min, max=", minval(ohx),  maxval(ohx)
  WRITE(6,*) "oqc:  min, max=", minval(oqc),  maxval(oqc)
 
ENDSUBROUTINE 

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
  WRITE(6,*) "[msg] "//trim(myname)//":: "
  WRITE(6,*) "Argument ", i, " = ",TRIM(arg1)

  select case (arg1)
    case('-obsin')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      WRITE(6,*) "Argument ", i+1, " = ",TRIM(arg2)
      obsinfile = arg2
    case('-obsout')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      WRITE(6,*) "Argument ", i+1, " = ",TRIM(arg2)
      obsoutfile = arg2
    case('-qcyyyymmddhh')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      WRITE(6,*) "Argument ", i+1, " = ",TRIM(arg2)
      Syyyymmddhh = arg2
    case('-maxdt')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      WRITE(6,*) "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) delta_seconds
    case('-ovar')
      CALL GET_COMMAND_ARGUMENT(i+1,arg2)
      WRITE(6,*) "Argument ", i+1, " = ",TRIM(arg2)
      read (arg2,*) ovar
    case default
      WRITE(6,*) "ERROR: option is not supported: ", arg1
      WRITE(6,*) "(with value : ", trim(arg2), " )"
      stop 1
  end select
enddo

END SUBROUTINE process_command_line

END PROGRAM decoder_en4
