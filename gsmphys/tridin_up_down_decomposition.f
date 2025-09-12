!!   Original subroutine 'tridin' taken from gsmphys/moninedmf.f
!!   and decomposed here into two routines:
!!   tridin_down: for the downward sweep, elminting the lower diagonal
!!   tridin_up: for the upward sweep, eliminating the upper diagonal
!!   The main motivation here is for implicit land coupling where
!!   the temperature&moisture at the end of the _down sweep should be
!!   passed to the land model, updated, then passed back for the upward
!!   sweep.
!!   Joseph.mouallem@noaa.gov
!!   Aug, 2025


!!!   Original Matrix:
!!!
!!!   │ cm1  cu1   0    0  │ │ x1 │   │ r1 │
!!!   │ cl2  cm2  cu2   0  │ │ x2 │ = │ r2 │
!!!   │  0   cl3  cm3  cu3 │ │ x3 │   │ r3 │
!!!   │  0    0   cl3  cm3 │ │ x4 │   │ r4 │


!!!   Apply tridin_down

!!!   │  1   au1    0     0  │ │ x1 │   │ a11 │
!!!   │  0    1    au2    0  │ │ x2 │ = │ a12 │
!!!   │  0    0     1    au3 │ │ x3 │   │ a13 │
!!!   │  0    0     0     1  │ │ x4 │   │ a14 │


!!!   Apply tridin_up

!!!   │  1    0     0     0  │ │ x1 │   │ a11 │
!!!   │  0    1     0     0  │ │ x2 │ = │ a12 │
!!!   │  0    0     1     0  │ │ x3 │   │ a13 │
!!!   │  0    0     0     1  │ │ x4 │   │ a14 │


      subroutine tridin_down(im,km,ntrac1,cl,cm,cu,r1,r2,au,a1,a2)

!!    Downward sweep for heat and moisture system
!!    Based on tridin but focused on heat (r1,a1) and moisture+tracers (r2,a2)
!!    cl, cm, cu, r1, r2 are intputs
!!    au, a1, a2 are outputs

      use machine     , only : kind_phys
      implicit none
      integer             is,k,kk,km,ntrac1,im,i
      real(kind=kind_phys) fk(im)
      real(kind=kind_phys) cl(im,2:km), cm(im,km), cu(im,km-1),         &
     &                     r1(im,km),   r2(im,km*ntrac1),                &
     &                     au(im,km-1), a1(im,km), a2(im,km*ntrac1),     &
     &                     fkk(im,2:km-1)

! for the first row
      do i=1,im
        fk(i)   = 1./cm(i,1)                 ! Normalize the diag in first row
        au(i,1) = fk(i)*cu(i,1)              ! Modified upper diagonal
        a1(i,1) = fk(i)*r1(i,1)              ! Update RHS
      enddo
      
! Same but for other tracers
      do k = 1, ntrac1 
        is = (k-1) * km
        do i = 1, im 
          a2(i,1+is) = fk(i) * r2(i,1+is)       ! Moisture + other tracers
        enddo 
      enddo   
      
! iterate from 2 to km-1 levels
      do k=2,km-1
        do i=1,im
          fkk(i,k) = 1./(cm(i,k)-cl(i,k)*au(i,k-1))  ! New normalization
          au(i,k)  = fkk(i,k)*cu(i,k)                ! Upper diag 
          a1(i,k)  = fkk(i,k)*(r1(i,k)-cl(i,k)*a1(i,k-1))  ! RHS 
        enddo
      enddo
      
! Same for tracers
      do kk = 1, ntrac1
        is = (kk-1) * km
        do k=2,km-1
          do i=1,im
            a2(i,k+is) = fkk(i,k)*(r2(i,k+is)-cl(i,k)*a2(i,k+is-1))
          enddo
        enddo
      enddo
      
! Solve for the lowest level
      do i=1,im
        fk(i)   = 1./(cm(i,km)-cl(i,km)*au(i,km-1))
        a1(i,km) = fk(i)*(r1(i,km)-cl(i,km)*a1(i,km-1))  ! TEMP !
      enddo
      
! Same for tracers
      do k = 1, ntrac1
        is = (k-1) * km
        do i = 1, im
          a2(i,km+is) = fk(i)*(r2(i,km+is)-cl(i,km)*a2(i,km+is-1))
        enddo
      enddo
c-----------------------------------------------------------------------
      return
      end




      subroutine tridin_up(im,km,ntrac1,au,a1,a2)

!!    Upward sweep
!!    Now solving from bot to top
!!    au, a1, a2 are inputs
!!    a1, a2 are outputs

      use machine     , only : kind_phys
      implicit none
      integer             is,k,kk,km,ntrac1,im,i
      real(kind=kind_phys) au(im,km-1), a1(im,km), a2(im,km*ntrac1)

! For temp
      do k=km-1,1,-1
        do i=1,im
          a1(i,k) = a1(i,k) - au(i,k)*a1(i,k+1)
        enddo
      enddo
      
! for moisture and tracers
      do kk = 1, ntrac1
        is = (kk-1) * km
        do k=km-1,1,-1
          do i=1,im
            a2(i,k+is) = a2(i,k+is) - au(i,k)*a2(i,k+is+1)
          enddo
        enddo
      enddo
c-----------------------------------------------------------------------
      return
      end
