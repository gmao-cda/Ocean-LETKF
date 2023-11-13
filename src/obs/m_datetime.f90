! https://github.com/python/cpython/blob/main/Lib/_pydatetime.py
module m_datetime
  !use common,  only: r_size, r_sngl
  implicit none

  integer,parameter :: r_size = 8
  integer,parameter :: r_sngl = 4

  private
  public datetime, timedelta  ! constructor
  public t_datetime
  public t_timedelta
  !public movedate

  interface datetime
    procedure :: init_datetime
  end interface

  interface timedelta
    procedure :: init_timedelta
  end interface

  type :: t_datetime
    integer,public :: year
    integer,public :: month
    integer,public :: day
    integer,public :: hour
    integer,public :: minute
    integer,public :: second
  contains
    procedure,pass(this),public :: print => print_datetime
    procedure,pass(this),public :: total_seconds_since19780101
    procedure,pass(this),private :: copy_datetime
    procedure,pass(this),public :: string => string_datetime

    procedure,pass(dt),private :: add_dt_delta, add_delta_dt
    procedure,pass(dt),private :: minus_dt_delta, minus_dt_dt
    generic, public :: operator(+)   => add_dt_delta, add_delta_dt
    generic, public :: operator(-)   => minus_dt_delta, minus_dt_dt
    generic, public :: assignment(=) => copy_datetime
  end type t_datetime

  type :: t_timedelta
     integer,private :: days
     integer,private :: seconds
  contains
     procedure,pass(this),public :: total_seconds
     procedure,pass(this),public :: print => print_timedelta
     procedure,pass(this),private :: copy_timedelta
     generic, public :: assignment(=) => copy_timedelta
  end type

  integer,private :: lout_log = 6
  character(*),private,parameter :: modname = "m_datetime"

contains


!-------------------------------------------------------------------------------
! timedelta
!
subroutine print_timedelta(this,lout)
  implicit none
  class(t_timedelta),intent(in) :: this
  integer,intent(in),optional :: lout
  integer :: lout_ 

  lout_ = lout_log
  if (present(lout)) lout_ = lout
  write(lout_,"(A,I10,A,I10,A)") "t_timedelta(days=",this%days, ", seconds=",this%seconds,")"

end subroutine print_timedelta

elemental function init_timedelta(days,hours,minutes,seconds) result(delta)
  implicit none
  integer,intent(in),optional:: days, hours, minutes, seconds
  type(t_timedelta) :: delta

  character(*),parameter :: myname = modname//"::init_timedelta"
  integer :: days_, hours_, minutes_, seconds_
  integer :: total_seconds_

  days_ = 0; hours_ = 0; minutes_ = 0; seconds_ = 0
  if (present(days))    days_    = days
  if (present(hours))   hours_   = hours
  if (present(minutes)) minutes_ = minutes
  if (present(seconds)) seconds_ = seconds
  !call check_timedelta_fields_(days_,hours_,minutes_,seconds_)

  total_seconds_ = days_*24*3600 + hours_*3600 + minutes_*60 + seconds_
  delta%seconds = mod(total_seconds_,24*3600)
  delta%days    = (total_seconds_ - delta%seconds) / (24*3600)
  
end function init_timedelta

integer function total_seconds(this)
  implicit none
  class(t_timedelta),intent(in) :: this
  total_seconds = 24*3600*this%days + this%seconds
end function total_seconds

subroutine check_timedelta_fields_(days,hours,minutes,seconds)
  implicit none
  integer,intent(in) :: days, hours, minutes, seconds
  character(*),parameter :: myname = "check_delta_fields_"

  if (days<0) then
     write(lout_log,*) "[err] "//myname//": invalid range: days<0 for days=", days
     stop (30)
  endif

  if (hours<0) then
     write(lout_log,*) "[err] "//myname//": invalid range: hours<0 for hours=", hours
     stop (30)
  endif

  if (minutes<0) then
     write(lout_log,*) "[err] "//myname//": invalid range: minutes<0 for minutes=", minutes
     stop (30)
  endif

  if (seconds<0) then
     write(lout_log,*) "[err] "//myname//": invalid range: seconds<0 for seconds=", seconds
     stop (30)
  endif
end subroutine check_timedelta_fields_


subroutine copy_timedelta(this,from)
  implicit none
  class(t_timedelta),intent(inout) :: this
  class(t_timedelta),intent(in) :: from
  this%days = from%days; this%seconds = from%seconds
end subroutine copy_timedelta


!-------------------------------------------------------------------------------
! datatime
!

function string_datetime(this) result(string)
  implicit none
  class(t_datetime),intent(in) :: this
  character(len=29)            :: string  ! 2019/10/01 00:00:00

  write(string,"(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)") &
       this%year, "/", this%month, "/", this%day, " ", &
       this%hour, ":", this%minute, ":", this%second

end function string_datetime

subroutine print_datetime(this, lout)
  implicit none
  class(t_datetime),intent(in) :: this
  integer,intent(in),optional :: lout
  integer :: lout_ 
  
  lout_ = lout_log
  if (present(lout)) lout_ = lout

  write(lout_,"(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A)") &
        "t_datetime(", this%year, "/", this%month, "/", this%day, " ", &
        this%hour, ":", this%minute, ":", this%second,")"

end subroutine  print_datetime

function init_datetime(year,month,day,hour,minute,second) result (dt)
  implicit none

  integer,intent(in) :: year
  integer,intent(in) :: month
  integer,intent(in) :: day
  integer,intent(in),optional :: hour
  integer,intent(in),optional :: minute
  integer,intent(in),optional :: second
  type(t_datetime) :: dt

  integer :: year_, month_, day_, hour_, minute_, second_

  year_ = year; month_ = month; day_ = day
  hour_ = 0; minute_ = 0; second_ = 0
  if (present(hour))   hour_   = hour
  if (present(minute)) minute_ = minute
  if (present(second)) second_ = second
  call check_date_fields_(year_,month_,day_)
  call check_time_fields_(hour_,minute_,second_)

  dt%year = year_; dt%month = month_; dt%day = day_
  dt%hour = hour_; dt%minute = minute_; dt%second = second_

end function init_datetime


integer function days_in_month(year,month)
  implicit none
  integer,intent(in) :: year, month

  integer,parameter,dimension(12) :: DAYS_IN_MONTH_ = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
  if (month==2 .and. is_leap_year(year)) then
     days_in_month = DAYS_IN_MONTH_(month) + 1
  else
     days_in_month = DAYS_IN_MONTH_(month)
  endif

end function days_in_month

logical function is_leap_year(year)
  implicit none

  integer,intent(in) :: year

  if (mod(year,4) == 0 .and. (mod(year,100)/=0 .or. mod(year,400)==0) ) then
     is_leap_year = .true.
  else
     is_leap_year = .false.
  endif

end function is_leap_year


subroutine check_date_fields_(year,month,day)
  implicit none

  integer,intent(in) :: year, month, day
  
  integer :: maxday
  character(*),parameter :: myname = trim(modname)//"::check_date_fields_"
  if (year<1978 .or. year>9999) then 
  !if (year<1 .or. year>9999) then
     write(lout_log,*) "[err] "//trim(myname)//": invalid year range (year<1 .or. year>9999) for year:", year
     stop (20)
  endif
  if (month<1 .or. month>12) then
     write(lout_log,*) "[err] "//trim(myname)//": invalid year range (month<1 .or. month>12) for month:", month
     stop (21)
  endif
  maxday = days_in_month(year,month)
  if (day<1 .or. day > maxday) then
     write(lout_log,"(A,I2,A,I2)") "[err] "//trim(myname)//": invalid year range (day<1 .or. day>",maxday,") for day:", day
     stop (22)
  endif

end subroutine check_date_fields_


subroutine check_time_fields_(hour,minute,second)
  implicit none
  integer,intent(in) :: hour, minute, second
  character(*),parameter :: myname = modname//"::check_time_fields_"

  if (hour<0 .or. hour>23) then
     write(lout_log,*) "[err] "//trim(myname)//": invalid hour range (hour<1 .or. hour>9999) for hour:", hour
     stop (11)
  endif
  if (minute<0 .or. minute > 59) then
     write(lout_log,*) "[err] "//trim(myname)//": invalid minute range (minute<1 .or. minute>9999) for minute:", minute
     stop (12)
  endif
  if (second<0 .or. second>59) then
     write(lout_log,*) "[err] "//trim(myname)//": invalid second range (second<1 .or. second>9999) for second:", second
     stop (13)
  endif

end subroutine check_time_fields_

integer function total_seconds_since19780101(this) 
  implicit none

  class(t_datetime),intent(in) :: this

  integer :: sdate(5)
  integer :: total_minutes_since19780101

  sdate = [this%year, this%month, this%day, this%hour, this%minute] 
          ! YYYY/MON/DAY/TZONE/HR/MIN/SEC/MSEC 
  call W3FS21(sdate(1:5), total_minutes_since19780101)
  total_seconds_since19780101 = total_minutes_since19780101*60 + &
                                this%second

end function total_seconds_since19780101

elemental subroutine copy_datetime(this,from)
  implicit none
  class(t_datetime),intent(inout) :: this
  class(t_datetime),intent(in) :: from

  this%year = from%year; this%month = from%month; this%day   = from%day
  this%hour  = from%hour; this%minute = from%minute; this%second = from%second

end subroutine copy_datetime


function add_dt_delta(dt,delta) result(newdt)
  implicit none

  !type(t_datetime),intent(in) :: dt
  class(t_datetime),intent(in) :: dt
  type(t_timedelta),intent(in) :: delta

  type(t_datetime) :: newdt
  integer :: sdate(8), edate(8)
  real :: tdelta(5)

  sdate = [dt%year,dt%month,dt%day,0,dt%hour,dt%minute,dt%second,0] !YYYY/MON/DAY/TZONE/HR/MIN/SEC/MSEC
  tdelta = [real(delta%days),0.0,0.0,real(delta%seconds),0.0] ! DAY/HR/MIN/SEC/MSEC
  call w3movdat(tdelta, sdate, edate)
  newdt = datetime(year   = edate(1), &
                   month  = edate(2), &
                   day    = edate(3), &
                   hour   = edate(5), &
                   minute = edate(6), &
                   second = edate(7) )
  
end function add_dt_delta

function add_delta_dt(delta,dt) result(newdt)
implicit none

class(t_datetime),intent(in) :: dt
type(t_timedelta),intent(in) :: delta

type(t_datetime) :: newdt
integer :: sdate(8), edate(8)
real :: tdelta(5)

sdate = [dt%year,dt%month,dt%day,0,dt%hour,dt%minute,dt%second,0] !YYYY/MON/DAY/TZONE/HR/MIN/SEC/MSEC
tdelta = [real(delta%days),0.0,0.0,real(delta%seconds),0.0] ! DAY/HR/MIN/SEC/MSEC
call w3movdat(tdelta, sdate, edate)
newdt = datetime(year   = edate(1), &
                 month  = edate(2), &
                 day    = edate(3), &
                 hour   = edate(5), &
                 minute = edate(6), &
                 second = edate(7) )

end function add_delta_dt

function minus_dt_delta(dt,delta) result(newdt)
  implicit none

  !type(t_datetime),intent(in) :: dt
  class(t_datetime),intent(in) :: dt
  type(t_timedelta),intent(in) :: delta

  type(t_datetime) :: newdt
  integer :: sdate(8), edate(8)
  real :: tdelta(5)

  sdate = [dt%year,dt%month,dt%day,0,dt%hour,dt%minute,dt%second,0] !YYYY/MON/DAY/TZONE/HR/MIN/SEC/MSEC
  tdelta = [-real(delta%days),0.0,0.0,-real(delta%seconds),0.0] ! DAY/HR/MIN/SEC/MSEC
  call w3movdat(tdelta, sdate, edate)
  newdt = datetime(year   = edate(1), &
                   month  = edate(2), &
                   day    = edate(3), &
                   hour   = edate(5), &
                   minute = edate(6), &
                   second = edate(7) )
  
end function minus_dt_delta

function minus_dt_dt(dt,dt2) result(delta)
  implicit none

  class(t_datetime),intent(in) :: dt
  class(t_datetime),intent(in) :: dt2
  type(t_timedelta)  :: delta

  integer :: total_seconds_diff

  total_seconds_diff = dt%total_seconds_since19780101() - dt2%total_seconds_since19780101()
  delta = timedelta(seconds=total_seconds_diff)
  
end function minus_dt_dt


end module m_datetime


!program main
!  use m_datetime
!  implicit none
!
!  type(t_datetime) :: t(3)
!  type(t_timedelta) :: dt(2)
!  type(t_timedelta) :: arr_dt(2)
!
!  t(1) = datetime(2000,2,29,10,0,0)   
!  call t(1)%print()
!
!  dt(1) = timedelta(hours=-25)
!  call dt(1)%print()
!  print*, dt(1)%total_seconds()
!
!  !t(2) = movedate(t(1),dt(1))
!  !call t(2)%print()
!
!  t(2) = t(1) + dt(1)
!  call t(2)%print()
!
!  t(3) =  dt(1) + t(1)
!  call t(3)%print()
!
!  dt(2) = t(2)-t(1)
!  print*, dt(2)%total_seconds()
!
!  arr_dt = timedelta(hours=[1,2])
!  call arr_dt(1)%print()
!  call arr_dt(2)%print()
!
!
!end program
