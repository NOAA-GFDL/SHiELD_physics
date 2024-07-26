
  if (extrap) then
    do kp = 1, np
      logp = log(plev(kp))
      do j = 1, jm
        do i = 1, im
          if (logp <= pm(1,i,j)) then
              ap(i,j,kp) = a(i,j,1) + (a(i,j,3)-a(i,j,1)) * (logp-pm(1,i,j))/(pm(3,i,j)-pm(1,i,j))
          else if (logp >= pm(km,i,j)) then
              ap(i,j,kp) = a(i,j,km) + (a(i,j,km)-a(i,j,km-2)) * (logp-pm(km,i,j))/(pm(km,i,j)-pm(km-2,i,j))
          else 
            do k = 1, km-1
              if (logp >= pm(k,i,j) .and. logp <= pm(k+1,i,j)) then
                ap(i,j,kp) = a(i,j,k) + (a(i,j,k+1)-a(i,j,k)) * (logp-pm(k,i,j))/(pm(k+1,i,j)-pm(k,i,j))
                goto 1000
              endif
            enddo
          endif
1000    continue
        enddo
      enddo
    enddo
  else
    do kp = 1, np
      logp = log(plev(kp))
      do j = 1, jm
        do i = 1, im
          if (logp <= pm(1,i,j)) then
              ap(i,j,kp) = a(i,j,1)
          else if (logp >= pm(km,i,j)) then
              ap(i,j,kp) = a(i,j,km)
          else 
            do k = 1, km-1
              if (logp >= pm(k,i,j) .and. logp <= pm(k+1,i,j)) then
                ap(i,j,kp) = a(i,j,k) + (a(i,j,k+1)-a(i,j,k)) * (logp-pm(k,i,j))/(pm(k+1,i,j)-pm(k,i,j))
                goto 2000
              endif
            enddo
          endif
2000    continue
        enddo
      enddo
    enddo
  endif
