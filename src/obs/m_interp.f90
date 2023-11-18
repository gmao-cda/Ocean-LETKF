module m_interp
!========================================================================================
!
! [PURPOSE:] interpolation
!
! [HISTORY:]
!   01/09/2014 Cheng Da     created as Module_NWP.f90 in CISS
!   10/06/2017 Cheng Da     imported from CISS, and renamed as m_interp.f90
!                           the elegant use of t_interp is copied from 
!                           Dr. Yongjun Zheng's interp_mod.f90
!
!
!----------------------------------------------------------------------------------------
  !use m_kind,     only : r_size
  use common,     only : r_size
  implicit none

  private

  integer,parameter,private :: lout_log = 6

  public :: t_interp
  public :: interp_1d
  public :: interp_2d
  public :: bisect_1d
  public :: spline_1d

  public :: nearest_2d


  private :: ll2dist

  real(r_size),parameter,private :: r_missing = -999.9d0


  type :: t_interp
    logical :: lgotip = .false.
    integer, allocatable :: kk(:)
    integer, allocatable :: ii(:,:)
    integer, allocatable :: jj(:,:)
    !real(r_size),allocatable :: w1(:,:)     ! w1(0:1,nx)  
    !real(r_size),allocatable :: w2(:,:,:,:) ! w2(0:1,0:1,ni,nj), to save time to compute weight
  endtype


  interface interp_1d
    module procedure interp_1d_110
    module procedure interp_1d_111
  endinterface

  interface bisect_1d
    module procedure bisect_1d_10
    module procedure bisect_1d_11
  endinterface


  contains

!-------------------------------------------------------------------------------
!
! (ii-1,jj)          (ii,ij)
!   O----------------O  
!   |                |     (1-wj)
!   |                |
!   |             X  +  rj
!   |                |     (wj)
!   O-------------+--O
! (ii-1,jj-1)    ri  (ii,ij-1)
!        wi      (1-wi) 
!
!
!   (ii-1,ij-1): (1-wi)*(1-wj)
!   (ii  ,ij-1):    wi *(1-wj)
!   (ii-1,jj  ): (1-wi)*   wj
!   (ii,  jj  ):    wi *   wj 
subroutine interp_2d(v2d,ri,rj,val)
  implicit none

  real(r_size),intent(in) :: v2d(:,:)  ! v2d(nlon,nlat)
  real(r_size),intent(in) :: ri, rj
  real(r_size),intent(out) :: val

  integer :: ii, ij, ni
  real(r_size) :: wi, wj

  ni=size(v2d,1)

  ii = ceiling(ri); wi = ri-(ii-1)
  ij = ceiling(rj); wj = rj-(ij-1)

  if (ii <= ni) then
     val  = v2d(ii-1,ij-1) * (1-wi) * (1-wj) + &
            v2d(ii  ,ij-1) *    wi  * (1-wj) + &
            v2d(ii-1,ij  ) * (1-wi) *    wj  + &
            v2d(ii  ,ij  ) *    wi  *    wj   
  else
     val  = v2d(ii-1,ij-1) * (1-wi) * (1-wj) + &
            v2d(1   ,ij-1) *    wi  * (1-wj) + &
            v2d(ii-1,ij  ) * (1-wi) *    wj  + &
            v2d(1   ,ij  ) *    wi  *    wj   
  endif 

endsubroutine



subroutine destroy_interp( ip )
  implicit none

  type(t_interp), intent(inout) :: ip

  if (allocated(ip%kk)) deallocate(ip%kk)
  if (allocated(ip%ii)) deallocate(ip%ii)
  if (allocated(ip%jj)) deallocate(ip%jj)

endsubroutine


subroutine interp_1d_111( x, f, x0, f0, ip )
  implicit none

  real(r_size),  intent(in   ) :: x(:)
  real(r_size),  intent(in   ) :: f(:)
  real(r_size),  intent(in   ) :: x0(:)
  real(r_size),  intent(  out) :: f0(:)
  type(t_interp),intent(inout) :: ip

  real(r_size) :: w
  integer :: k, n, m, i
  
  m=size(x)
  n=size(x0)
  if (.not.ip%lgotip) call bisect_1d_11( x, x0, ip ) 
  do i = 1, n
     k = ip%kk(i)
     if ( k>=1 .and. k<=m ) then
         w     = (x0(i)-x(k))/(x(k+1)-x(k))
         f0(i) = (1-w)*f(k) + w*f(k+1)
     else
         f0(i) = r_missing
     endif
  enddo
  if (.not.ip%lgotip) call destroy_interp( ip )

endsubroutine

!---------------------------------------------------------------------------------------
! linear interpolation:
!    o----o-------------o
!  x(i)   x0            x(i+1)
!  f(i)   f0            f(i+1)
!---------------------------------------------------------------------------------------
subroutine interp_1d_110( x, f, x0, f0 )
  implicit none

  real(r_size),intent(in)  :: x(:)
  real(r_size),intent(in)  :: f(:)
  real(r_size),intent(in)  :: x0
  real(r_size),intent(out) :: f0

  integer :: ierr 
  integer :: i, n
  real(r_size) :: w
 
  n=size(x)
  call bisect_1d( x, x0, i )
  if (i<1.or.i>n) then
     f0   = r_missing
  else
     w  = (x0-x(i))/(x(i+1)-x(i))
     f0 = (1-w)*f(i) + w*f(i+1)
  endif

endsubroutine

!---------------------------------------------------------------------------------------
! bisection search: assuming no repeated elements in the target array in ascending order
!---------------------------------------------------------------------------------------
subroutine bisect_1d_10( x, x0, ix0 )  
  implicit none

  real(r_size),intent(in) :: x(:)
  real(r_size),intent(in) :: x0
  integer,     intent(out) :: ix0

  integer :: n
  integer :: l, r, m

  ix0=0
  n=size(x)
  l=1; r=n
  if ( x0 > x(r) ) then
     ix0 = n+1
  elseif ( x0 < x(l) ) then
     ix0 = 0
  else
    do while ( l+1 < r ) 
        m = (l+r)/2
        if ( x(m)>x0 ) then
           r = m
        else
           l = m
        endif
    enddo
    ix0 = l
  endif   

endsubroutine


subroutine bisect_1d_11( x, x0, ip )
  
  real(r_size),intent(in) :: x(:)
  real(r_size),intent(in) :: x0(:)
  type(t_interp),intent(inout) :: ip

  integer :: nx0, i

  nx0=size(x0)
  if (allocated(ip%kk)) deallocate(ip%kk)
  allocate( ip%kk(nx0) )
  do i = 1, nx0
     call bisect_1d_10( x, x0(i), ip%kk(i) )
  enddo

endsubroutine


!---------------------------------------------------------------------------------------
!  return the model lat & lon index closed to the obs. For the longitude, the obs 
!  between lon(end) & lon(1) are also taken care.
!---------------------------------------------------------------------------------------
subroutine nearest_2d( olon, olat, mlon, mlat, ix, iy )
  implicit none

  real(r_size),intent(in) :: olon, olat
  real(r_size),intent(in) :: mlon(:), mlat(:)
  integer,     intent(inout) :: ix, iy

  real(r_size) :: dis, mindis
  real(r_size) :: mlon2(0:1), mlat2(0:1)
  integer :: i2(0:1), j2(0:1)
  integer :: i, j

  mindis = 9.999e30
  if ( iy < size(mlat) ) then
     mlat2(0:1) = mlat(iy:iy+1)
        j2(0:1) = (/iy,iy+1/)
  else
     mlat2(0:1) = mlat(iy-1:iy)
        j2(0:1) = (/iy-1,iy/)
  endif
  if ( ix < size(mlon) ) then
     mlon2(0:1) = mlon(ix:ix+1)
        i2(0:1) = (/ix,ix+1/)
  else
     mlon2(0:1) = (/mlon(size(mlon)), mlon(1)/)
        i2(0:1) = (/size(mlon), 1/)
  endif
  do i = 0, 1
     do j = 0, 1
        dis = ll2dist( olat, olon, mlat2(j), mlon2(i) )
        !print*, "obs: lat, lon, mlat, mlon=", olat, olon, mlat2(j),mlon2(i), dis
        if ( dis<mindis ) then
           ix = i2(i); iy = j2(j); mindis = dis
        endif
     enddo
  enddo
  !print*, "----> ix, iy=", ix, iy
  !print*, "----> lon(ix), lat(iy)=", mlon(ix), mlat(iy)

endsubroutine
        
!----------------------------------------------------------------------------------------
!    Spline interpolation. The boundary condition are the natural boundary conditions, 
!    i.e., S"(x)=0 at end points.
!    The known points are xi(:) & fi(:), and the value fo(:) at the desired output 
!    location xo(:) are output.
!   
!    directly translated from my C version of spline.
!----------------------------------------------------------------------------------------
subroutine spline_1d( xi, fi, xo, fo )
 implicit none
! passed args
  real(r_size),intent(in)  :: xi(:), fi(:)
  real(r_size),intent(in)  :: xo(:)
  real(r_size),intent(out) :: fo(:)
! local vars
  integer :: n, ns
  real(r_size) :: h(size(xi)),a(size(xi)),b(size(xi)),c(size(xi)),Q(size(xi)),m(size(xi))
  real(r_size) :: c0,c1,c2,c3
  integer :: i, k
  character(*),parameter :: subname="spline_1d"

  n=size(xi); ns=size(xo)
  if (size(xi)/=size(fi) .or. size(xo)/=size(fo)) then
     write(lout_log,*) "[error] ", trim(subname), ": dimension mismatches:"
     write(lout_log,*) "   size: xi, fi=", size(xi),size(fi)
     write(lout_log,*) "   size: xo, fo=", size(xo),size(fo)
     stop 14
  endif

  h(:)=0.0_r_size
  h(1:n-1)=xi(2:n)-xi(1:n-1)
  a(:)=1.0_r_size
  b(:)=0.0_r_size;c(:)=0.0_r_size;Q(:)=0.0_r_size
  do i=2,n-1
     Q(i)=6*((fi(i+1)-fi(i))/h(i)-(fi(i)-fi(i-1))/h(i-1));
     a(i)=2*(h(i-1)+h(i));
     b(i)=h(i);
     c(i)=h(i-1);
  enddo
  call solve_band_matrix(n,a,b,c,Q,m)
  do k=1,ns
     i=1
     do i=1,n-1
        if ( xo(k)>=xi(i) .and. xo(k)<xi(i+1) ) exit
     enddo
     if ( i == n ) then
        i=i-1
     endif
     c0=fi(i)
     c1=(fi(i+1)-fi(i))/h(i)-h(i)*m(i)/2-h(i)*(m(i+1)-m(i))/6
     c2=m(i)/2
     c3=(m(i+1)-m(i))/h(i)/6
     fo(k)=c0+c1*(xo(k)-xi(i))+c2*(xo(k)-xi(i))**2+c3*(xo(k)-xi(i))**3
  enddo

endsubroutine

!----------------------------------------------------------------------------------------
!    solve band matrix AX=Q using Thomas method.
!    a(n) are the diagonal elements of the A
!    b(n) are the elements above a(n), from 1->n-1
!    c(n) are the elements below c(n), from 2->n
!----------------------------------------------------------------------------------------
subroutine solve_band_matrix(n, a, b, c, Q, x)

  implicit none
! passed args
  integer,intent(in)::n
  real(r_size),intent(in) :: a(n),b(n),c(n),Q(n)
  real(r_size),intent(out) :: x(n)
! local vars
  real(r_size) :: as(n), Qs(n)
  integer :: i

  x(:)=0.0_r_size
  as(1)=a(1)
  do i = 2,n
     as(i)=a(i)-c(i)*b(i-1)/as(i-1)
  enddo
  Qs(1)=Q(1)
  do i = 2,n
     Qs(i)=Q(i)-c(i)*Qs(i-1)/as(i-1)
  enddo
  x(n)=Qs(n)/as(n)
  do i = n-1,1,-1
     x(i)=(Qs(i)-b(i)*x(i+1))/as(i)
  enddo

endsubroutine

!-------------------------------------------------------------------------------
!  calculate great circle distance between two points
!-------------------------------------------------------------------------------
elemental function ll2dist( lat1, lon1, lat2, lon2 ) result ( dis )
 implicit none
! Passed args
  real(r_size),intent(in    ) :: lat1
  real(r_size),intent(in    ) :: lon1
  real(r_size),intent(in    ) :: lat2
  real(r_size),intent(in    ) :: lon2
  real(r_size) :: dis
! Local vars
  real(r_size) :: temp
! Constant 
  real(r_size),parameter :: earth_Rad = 6372.8d0 ! from GSI m_distance
  real(r_size),parameter :: pi = acos(-1.d0)

    
  temp =sin((lat1-lat2)*pi/180/2)**2+cos(lat1*pi/180)*&
             cos(lat2*pi/180)*(sin((lon1-lon2)*pi/180/2)**2)
  temp = sqrt(temp)
  dis = 2*asin(temp)
  dis = earth_Rad*dis

endfunction


endmodule

