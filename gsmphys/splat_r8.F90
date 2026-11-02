subroutine splat_r8(jmax, slat, wlat)
  use machine, only : kind_io8
  implicit none

  integer, intent(in) :: jmax
  real(kind_io8), intent(out) :: slat(jmax), wlat(jmax)

  integer, parameter :: jz = 50
  real(kind_io8), parameter :: pi = 3.14159265358979_kind_io8
  real(kind_io8), parameter :: c = &
      (1.0_kind_io8 - (2.0_kind_io8/pi)**2)*0.25_kind_io8
  real(kind_io8), parameter :: eps = &
      10.0_kind_io8*epsilon(1.0_kind_io8)

  real(kind_io8), parameter :: bz(jz) = [ &
      2.4048255577_kind_io8,   5.5200781103_kind_io8, &
      8.6537279129_kind_io8,  11.7915344391_kind_io8, &
     14.9309177086_kind_io8,  18.0710639679_kind_io8, &
     21.2116366299_kind_io8,  24.3524715308_kind_io8, &
     27.4934791320_kind_io8,  30.6346064684_kind_io8, &
     33.7758202136_kind_io8,  36.9170983537_kind_io8, &
     40.0584257646_kind_io8,  43.1997917132_kind_io8, &
     46.3411883717_kind_io8,  49.4826098974_kind_io8, &
     52.6240518411_kind_io8,  55.7655107550_kind_io8, &
     58.9069839261_kind_io8,  62.0484691902_kind_io8, &
     65.1899648002_kind_io8,  68.3314693299_kind_io8, &
     71.4729816036_kind_io8,  74.6145006437_kind_io8, &
     77.7560256304_kind_io8,  80.8975558711_kind_io8, &
     84.0390907769_kind_io8,  87.1806298436_kind_io8, &
     90.3221726372_kind_io8,  93.4637187819_kind_io8, &
     96.6052679510_kind_io8,  99.7468198587_kind_io8, &
    102.888374254_kind_io8,  106.029930916_kind_io8, &
    109.171489649_kind_io8,  112.313050280_kind_io8, &
    115.454612653_kind_io8,  118.596176630_kind_io8, &
    121.737742088_kind_io8,  124.879308913_kind_io8, &
    128.020877005_kind_io8,  131.162446275_kind_io8, &
    134.304016638_kind_io8,  137.445588020_kind_io8, &
    140.587160352_kind_io8,  143.728733573_kind_io8, &
    146.870307625_kind_io8,  150.011882457_kind_io8, &
    153.153458019_kind_io8,  156.295034268_kind_io8 ]

  real(kind_io8) :: pk(jmax/2), pkm1(jmax/2), pkm2(jmax/2)
  real(kind_io8) :: slatd(jmax/2)
  real(kind_io8) :: r, sp, spmax
  integer :: j, n, jh, jhe

  ! Gaussian latitudes: original SPLAT IDRT == 4 branch.
  jh  = jmax/2
  jhe = (jmax+1)/2

  r = 1.0_kind_io8 / &
      sqrt((real(jmax,kind_io8)+0.5_kind_io8)**2+c)

  do j = 1, min(jh,jz)
    slatd(j) = cos(bz(j)*r)
  enddo

  do j = jz+1, jh
    slatd(j) = cos((bz(jz)+real(j-jz,kind_io8)*pi)*r)
  enddo

  spmax = 1.0_kind_io8
  do while (spmax > eps)
    spmax = 0.0_kind_io8

    do j = 1, jh
      pkm1(j) = 1.0_kind_io8
      pk(j)   = slatd(j)
    enddo

    do n = 2, jmax
      do j = 1, jh
        pkm2(j) = pkm1(j)
        pkm1(j) = pk(j)
        pk(j) = ((2*n-1)*slatd(j)*pkm1(j) - &
                 (n-1)*pkm2(j))/n
      enddo
    enddo

    do j = 1, jh
      sp = pk(j)*(1.0_kind_io8-slatd(j)**2) / &
           (jmax*(pkm1(j)-slatd(j)*pk(j)))
      slatd(j) = slatd(j)-sp
      spmax = max(spmax,abs(sp))
    enddo
  enddo

  ! Weights and hemispheric symmetry.
  do j = 1, jh
    slat(j) = slatd(j)
    wlat(j) = 2.0_kind_io8*(1.0_kind_io8-slatd(j)**2) / &
              (jmax*pkm1(j))**2

    slat(jmax+1-j) = -slat(j)
    wlat(jmax+1-j) =  wlat(j)
  enddo

  ! Equatorial point for an odd number of latitudes.
  if (jhe > jh) then
    slat(jhe) = 0.0_kind_io8
    wlat(jhe) = 2.0_kind_io8/real(jmax,kind_io8)**2

    do n = 2, jmax, 2
      wlat(jhe) = wlat(jhe)*real(n,kind_io8)**2 / &
                            real(n-1,kind_io8)**2
    enddo
  endif

end subroutine splat_r8
