PROGRAM decoder_en4_nc
  USE common
  USE params_model
  USE params_obs, only : id_s_obs, id_t_obs
  USE vars_model
  USE common_oceanmodel
  !USE vars_obs
  USE common_obs_oceanmodel
  USE read_en4,      ONLY: ocn_profile, select_prof_en4_nc

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! Command line inputs:
  !-----------------------------------------------------------------------------
  CHARACTER(slen) :: obsinfile='obsin.nc'     !IN (default)
  CHARACTER(slen) :: obsoutfile='obsout.nc'  !OUT(default)
  CHARACTER(slen) :: ovar = "potm" ! IN (default), or "psal"
  INTEGER         :: otyp = id_t_obs
  CHARACTER(10) :: Syyyymmddhh = "YYYYMMDDHH"
  !                 1234567890
  INTEGER :: delta_seconds = -999 ! max difference from Syyymmddhh

  !-----------------------------------------------------------------------------
  ! Obs data arrays
  !-----------------------------------------------------------------------------
  TYPE(ocn_profile) :: prof


  CHARACTER(*),PARAMETER :: myname = "decoder_en4_nc4"

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
     CALL select_prof_en4_nc(trim(obsinfile),otyp,prof,Syyyymmddhh, delta_seconds)
  else
     CALL select_prof_en4_nc(trim(obsinfile),otyp,prof)
  endif


  call prof%write_nc(trim(obsoutfile))

  call prof%deallocate()

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

END PROGRAM decoder_en4_nc
