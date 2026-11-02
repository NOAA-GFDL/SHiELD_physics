subroutine splat_r8(jmax, slat, wlat)
  use machine, only : kind_io8
  implicit none
  integer, intent(in) :: jmax
  real(kind_io8), intent(out) :: slat(jmax), wlat(jmax)
  real(kind_io8), parameter :: pi = 3.14159265358979323846_kind_io8
  real(kind_io8) :: x(jmax/2), pk(jmax/2), pkm1(jmax/2), pkm2(jmax/2)
  real(kind_io8) :: sp, spmax, eps
  integer :: j, n, jh, jhe

  eps = 10.0_kind_io8*epsilon(1.0_kind_io8)
  jh  = jmax/2
  jhe = (jmax+1)/2
  do j = 1, jh
    x(j) = cos(pi*(real(j,kind_io8)-0.25_kind_io8)/(real(jmax,kind_io8)+0.5_kind_io8))
  enddo
  spmax = 1.0_kind_io8
  do while (spmax > eps)
    spmax = 0.0_kind_io8
    pkm1 = 1.0_kind_io8
    pk   = x
    do n = 2, jmax
      pkm2 = pkm1
      pkm1 = pk
      pk   = ((2*n-1)*x*pkm1 - (n-1)*pkm2)/n
    enddo
    do j = 1, jh
      sp   = pk(j)*(1.0_kind_io8-x(j)**2)/(jmax*(pkm1(j)-x(j)*pk(j)))
      x(j) = x(j) - sp
      spmax = max(spmax, abs(sp))
    enddo
  enddo
  do j = 1, jh
    slat(j) = x(j)
    wlat(j) = (2.0_kind_io8*(1.0_kind_io8-x(j)**2))/(jmax*pkm1(j))**2
    slat(jmax+1-j) = -slat(j)
    wlat(jmax+1-j) =  wlat(j)
  enddo
  if (jhe > jh) then
    slat(jhe) = 0.0_kind_io8
    wlat(jhe) = 2.0_kind_io8/real(jmax,kind_io8)**2
    do n = 2, jmax, 2
      wlat(jhe) = wlat(jhe)*real(n,kind_io8)**2/real(n-1,kind_io8)**2
    enddo
  endif
end subroutine splat_r8
