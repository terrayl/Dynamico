subroutine distancef(vart, varf, ntt, ntf, nxy, dist)
! Estimate the Euclidean distance between the given atmospheric circulation 
! variable, say X, daily or monthly mean and all other daily or monthly means within the training period
! 
! 
! INPUT->  vart(nxy,ntt)      : X array for all months/days of the training period
!          varf(nxy,ntf)      : X array for the month/day under scrutiny for the analogue search
!          ntt                : number of months/days in the training period 
!          ntf                : number of months/days under scrutiny for the analogue search 
!          nxy                : number of gridpoints
!
! OUTPUT-> dist(ntf,ntt)      : distance array
! 
! 
  implicit none
   integer                     ::ntt, ntf, nxy
   real                        vart(nxy, ntt), varf(nxy, ntf), dist(ntf, ntt)
   integer                     ::n, m

   !$OMP PARALLEL DO PRIVATE(m)
   do n=1,ntt

      do m=1,ntf
         dist(m,n)              = sqrt(sum(( vart(:,n) - varf(:,m) )**2))
      end do

   end do
   !$OMP END PARALLEL DO

end subroutine distancef
